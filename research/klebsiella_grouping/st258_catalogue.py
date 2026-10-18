"""Audit a historical ST258 export against a profile-fingerprint DNA bridge.

The original scope remains historical and unversioned. Numeric labels are
accepted only after every sample's provider cgmlstSt fingerprint is reproduced
from the frozen oriented reference SHA1s and literal novel allele hashes.
"""

import argparse
from collections import Counter
from dataclasses import asdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3

from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

from chronoclade.cgmlst.dna import CATALOGUE_SCHEMA
from chronoclade.cgmlst.locus_alignment import file_hash, save_json
from chronoclade.datasets import load_dataset


def write_csv(path, rows, fields=None):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest_object(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def audit(dataset_path, public_snapshot, reference_fastas, index, caller_source, output, supplementary_catalogues=()):
    dataset_path, public_snapshot, reference_fastas, index, caller_source, output = (
        Path(p).resolve() for p in (dataset_path, public_snapshot, reference_fastas, index, caller_source, output))
    dataset = load_dataset(dataset_path)
    if len(dataset.profiles) != 1:
        raise ValueError("Expected exactly one historical profile matrix")
    matrix = dataset.profiles[0]
    sample_info = {r["sample_id"]: r for r in dataset.samples}
    if any(r["species"] != "Klebsiella pneumoniae" for r in dataset.samples):
        raise ValueError("Expected a single Klebsiella pneumoniae taxon scope")
    details = json.loads(public_snapshot.read_text())["raw_details"]
    with gzip.open(index / "metadata.json.gz", "rt") as stream:
        genes = json.load(stream)["genes"]
    if set(matrix.catalogue.loci) - set(genes):
        raise ValueError("Historical loci lie outside the explicit current reference panel")
    inventory, observed_keys = {}, set()
    role_masks = {role: np.array([sample_info[s]["role"] == role for s in matrix.sample_ids])
                  for role in ("input", "context")}
    for j, locus in enumerate(matrix.catalogue.loci):
        codes, counts = np.unique(matrix.codes[:, j], return_counts=True)
        for code, count in zip(codes, counts, strict=True):
            if not code:
                continue
            allele = matrix.categories[int(code) - 1]
            kind = "numeric" if allele.isdigit() else "sha1" if re.fullmatch(r"[0-9a-f]{40}", allele) else "unsupported"
            if kind == "unsupported":
                raise ValueError(f"Unsupported observed allele at {locus}")
            key = locus, allele
            observed_keys.add(key)
            inventory[key] = {"locus": locus, "allele": allele, "kind": kind, "called_cells": int(count),
                              **{role + "_cells": int(np.count_nonzero((matrix.codes[:, j] == code) & mask))
                                 for role, mask in role_masks.items()}}
    numeric_hashes = {}
    db = index / "allele.db"
    with sqlite3.connect(f"file:{db}?mode=ro&immutable=1", uri=True) as connection:
        for dna_hash, locus, allele in connection.execute("SELECT hash,gene,st FROM alleles WHERE reverse=0"):
            key = locus, str(allele)
            if key not in observed_keys:
                continue
            if key in numeric_hashes and numeric_hashes[key] != dna_hash:
                raise ValueError("Conflicting numeric-to-oriented-SHA1 mapping")
            numeric_hashes[key] = dna_hash
    absent_numeric = [key for key in observed_keys if key[1].isdigit() and key not in numeric_hashes]
    if absent_numeric:
        raise ValueError(f"Cannot verify profile bridge: {len(absent_numeric)} numeric definitions absent")
    fingerprints = []
    # This is the pinned caller's hasha(checksums,{algorithm:'sha1'}): array
    # elements update SHA1 consecutively, without a separator; absent calls=''.
    for i, ident in enumerate(matrix.sample_ids):
        profile = {locus: matrix.categories[int(code) - 1] for locus, code in
                   zip(matrix.catalogue.loci, matrix.codes[i], strict=True) if code}
        source_id = sample_info[ident]["source_genome_id"]
        provider = details.get(source_id, {})
        expected = provider.get("cgmlstSt")
        if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{40}", expected):
            raise ValueError(f"No historical provider DNA-profile fingerprint for {ident}")
        hashes = [numeric_hashes[locus, profile[locus]] if profile.get(locus, "").isdigit()
                  else profile.get(locus, "") for locus in genes]
        actual = hashlib.sha1("".join(hashes).encode()).hexdigest()
        if actual != expected:
            raise ValueError(f"Historical/current profile DNA fingerprint differs for {ident}")
        fingerprints.append({"sample_id": ident, "source_genome_id": source_id,
                             "provider_cgmlstSt": expected, "reconstructed_cgmlstSt": actual,
                             "verified": True})
    del details
    sequences, reference_files = {}, []
    for locus in matrix.catalogue.loci:
        path = reference_fastas / f"{locus}.fa.gz"
        reference_files.append({"path": str(path), "sha256": file_hash(path)})
        requested = {allele for loc, allele in observed_keys if loc == locus}
        requested_hashes = {allele for allele in requested if not allele.isdigit()}
        with gzip.open(path, "rt") as stream:
            for record in SeqIO.parse(stream, "fasta"):
                dna = str(record.seq).upper()
                forward = hashlib.sha1(dna.lower().encode()).hexdigest()
                reverse = hashlib.sha1(str(Seq(dna).reverse_complement()).lower().encode()).hexdigest()
                matches = []
                if record.id in requested:
                    if forward != numeric_hashes[locus, record.id]:
                        raise ValueError(f"Full reference FASTA/SQLite DNA mismatch at {locus}:{record.id}")
                    matches.append((record.id, "profile-fingerprint-confirmed-numeric"))
                if forward in requested_hashes:
                    matches.append((forward, "exact-forward-sequence-sha1"))
                if reverse in requested_hashes and reverse != forward:
                    matches.append((reverse, "exact-reverse-complement-sha1-canonicalized"))
                for allele, method in matches:
                    key = locus, allele
                    if set(dna) - set("ACGT"):
                        inventory[key]["resolution"] = "reference_sequence_ambiguous"
                        continue
                    if key in sequences and sequences[key] != dna:
                        raise ValueError("One source allele resolves to conflicting canonical DNA")
                    sequences[key] = dna
                    inventory[key].update(resolution=method, reference_allele=record.id,
                                          length=len(dna), canonical_sequence_sha1=forward)
    supplementary_sources = []
    for catalogue in supplementary_catalogues:
        catalogue = Path(catalogue).resolve()
        manifest = json.loads(catalogue.read_text())
        if manifest.get("schema") != CATALOGUE_SCHEMA or manifest.get("status") != "complete":
            raise ValueError("Supplementary DNA must be a completed frozen catalogue")
        dna_csv = (catalogue.parent / manifest["artifact"]["path"]).resolve()
        if not dna_csv.is_relative_to(catalogue.parent) or file_hash(dna_csv) != manifest["artifact"]["sha256"]:
            raise ValueError("Supplementary catalogue DNA checksum or local path is invalid")
        supplementary_sources.extend({"path": str(p), "sha256": file_hash(p)} for p in (catalogue, dna_csv))
        provenance = catalogue.parent / "provenance.json"
        if provenance.exists():
            supplementary_sources.append({"path": str(provenance), "sha256": file_hash(provenance)})
        with dna_csv.open() as stream:
            for row in csv.DictReader(stream):
                dna = row["dna"].upper()
                if not dna or set(dna) - set("ACGT"):
                    raise ValueError("Supplementary catalogue contains unsupported DNA")
                dna_hash = hashlib.sha1(dna.lower().encode()).hexdigest()
                key = row["locus"], dna_hash
                if key not in observed_keys or key in sequences:
                    continue
                sequences[key] = dna
                inventory[key].update(resolution="exact-sha1-verified-native-call-catalogue",
                                      reference_allele=f"{catalogue}:{row['allele']}",
                                      length=len(dna), canonical_sequence_sha1=dna_hash)
    output.mkdir(parents=True, exist_ok=True)
    for key, row in inventory.items():
        row.setdefault("resolution", "unknown_hash_sequence_unavailable")
        for field in ("reference_allele", "length", "canonical_sequence_sha1"):
            row.setdefault(field, "")
        row["dna_available"] = key in sequences
    write_csv(output / "allele-audit.csv", list(inventory.values()))
    write_csv(output / "fingerprint-verification.csv", fingerprints)
    write_csv(output / "alleles.csv", [{"locus": locus, "allele": allele, "dna": dna}
                                      for (locus, allele), dna in sorted(sequences.items())])
    unresolved_by_sample, samples = [], []
    for i, ident in enumerate(matrix.sample_ids):
        unresolved = [{"sample_id": ident, "locus": locus, "allele": matrix.categories[int(code) - 1]}
                      for locus, code in zip(matrix.catalogue.loci, matrix.codes[i], strict=True)
                      if code and (locus, matrix.categories[int(code) - 1]) not in sequences]
        unresolved_by_sample.extend(unresolved)
        called = int(np.count_nonzero(matrix.codes[i]))
        samples.append({"sample_id": ident, "role": sample_info[ident]["role"], "called_loci": called,
                        "missing_calls": len(matrix.catalogue.loci) - called,
                        "unresolved_dna_calls": len(unresolved), "resolved_loci": called - len(unresolved)})
    write_csv(output / "sample-audit.csv", samples)
    write_csv(output / "unresolved-cells.csv", unresolved_by_sample, ["sample_id", "locus", "allele"])
    bridge = {"schema": "chronoclade.research.provider-dna-profile-bridge", "schema_version": 1,
              "original_dataset_id": dataset.dataset_id, "original_scope": asdict(matrix.catalogue),
              "policy": "Accept numeric mappings only after every historical provider cgmlstSt fingerprint matches reconstructed ordered oriented allele SHA1s",
              "historical_database_version": "unknown; no equality with current reference asserted",
              "historical_fingerprints_verified": len(fingerprints),
              "fingerprint_order": genes, "comparison_loci": list(matrix.catalogue.loci),
              "comparison_universe_complete": False,
              "source_files": [{"path": str(p), "sha256": file_hash(p)}
                               for p in (dataset_path, public_snapshot, db, index / "metadata.json.gz", caller_source,
                                         caller_source.parent.parent / "node_modules/hasha/index.js")],
              "reference_fasta_files": reference_files,
              "supplementary_sequence_sources": supplementary_sources,
              "allele_audit_sha256": file_hash(output / "allele-audit.csv"),
              "fingerprint_verification_sha256": file_hash(output / "fingerprint-verification.csv"),
              "dna_csv_sha256": file_hash(output / "alleles.csv")}
    bridge_hash = digest_object(bridge)
    save_json(output / "bridge.json", bridge)
    scope = {"species": "Klebsiella pneumoniae", "scheme_id": "research:st258-provider-profile-dna-bridge",
             "scheme_version": "observed-panel-sha1-bridge-v1", "database_version": None,
             "database_sha256": bridge_hash}
    save_json(output / "catalogue.json", {"schema": CATALOGUE_SCHEMA, "schema_version": 1,
        "status": "complete", "scope": scope, "artifact": {"path": "alleles.csv", "sha256": bridge["dna_csv_sha256"]},
        "bridge": {"path": "bridge.json", "canonical_content_sha256": bridge_hash},
        "coverage": "Resolved observed DNA only; original 628-locus panel is explicitly incomplete"})
    missing_keys = [key for key in observed_keys if key not in sequences]
    affected = sum(row["unresolved_dna_calls"] > 0 for row in samples)
    summary = {"samples": len(samples), "input_samples": int(role_masks["input"].sum()),
               "context_samples": int(role_masks["context"].sum()), "observed_loci": len(matrix.catalogue.loci),
               "complete_locus_universe": False, "provider_fingerprints_verified": len(fingerprints),
               "observed_alleles": len(observed_keys), "observed_kind_counts": dict(Counter(r["kind"] for r in inventory.values())),
               "resolved_alleles": len(sequences), "unresolved_alleles": len(missing_keys),
               "resolution_counts": dict(Counter(r["resolution"] for r in inventory.values())),
               "called_cells": int(np.count_nonzero(matrix.codes)), "original_missing_cells": int((matrix.codes == 0).sum()),
               "unresolved_dna_cells": len(unresolved_by_sample), "samples_affected_by_unresolved_dna": affected,
               "unresolved_input_cells": sum(r["unresolved_dna_calls"] for r in samples if r["role"] == "input"),
               "resolved_loci_min": min(r["resolved_loci"] for r in samples),
               "resolved_loci_median": float(np.median([r["resolved_loci"] for r in samples])),
               "resolved_loci_max": max(r["resolved_loci"] for r in samples),
               "bridge_scope": scope, "bridge_content_sha256": bridge_hash,
               "scope_limitation": "This is a newly fingerprint-verified observed reference bridge, not identification of the historical provider database version",
               "missing_policy": "Unresolved allele DNA remains missing evidence; never a match or a zero substitution call"}
    save_json(output / "summary.json", summary)
    (output / "AUDIT.txt").write_text(
        f"Historical ST258 provider export: {len(samples)} samples, {len(matrix.catalogue.loci)} observed loci (incomplete universe).\n"
        f"All {len(fingerprints)} provider DNA-profile fingerprints reproduced exactly; numeric mapping is independently confirmed.\n"
        f"Observed alleles {len(observed_keys)}; resolved DNA {len(sequences)}; unresolved {len(missing_keys)}.\n"
        f"Unresolved called cells {len(unresolved_by_sample)} across {affected} samples; original missing cells {summary['original_missing_cells']}.\n"
        "Historical database version remains unknown. The bridge is a distinct frozen comparison representation.\n"
        "Missing DNA must be N/unavailable, never a match. Per-sample and per-cell exclusions are retained.\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default="tmp/st258-native-e2e/prepared/dataset.json")
    parser.add_argument("--public-snapshot", default="tmp/greek-st258-cgmlst/run_context/sources/provider/public_0.json")
    parser.add_argument("--reference-fastas", default="tmp/typing-setup-smoke/databases/mlst-schemes/cgmlst_schemes/pubmlst_klebsiella_18")
    parser.add_argument("--index", default="tmp/typing-setup-smoke/databases/mlst-index/klebsiella_1")
    parser.add_argument("--caller-source", default="tmp/typing-setup-smoke/sources/mlst/src/mlst.js")
    parser.add_argument("--output", default="tmp/science-st258/catalogue")
    parser.add_argument("--supplementary-catalogue", action="append",
                        default=["tmp/science-zhou2017/locus-dna/catalogue.json"])
    args = parser.parse_args()
    print(json.dumps(audit(args.dataset, args.public_snapshot, args.reference_fastas, args.index, args.caller_source,
                           args.output, args.supplementary_catalogue), indent=2))


if __name__ == "__main__":
    main()
