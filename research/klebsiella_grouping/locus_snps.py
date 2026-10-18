"""Resolve observed native cgMLST DNA and compare the frozen labelled panel.

Run from the repository root with ``pixi run python -m
research.klebsiella_grouping.locus_snps --help``. No network or caller rerun.
"""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sqlite3

from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from scipy.stats import spearmanr

from chronoclade.cgmlst.dna import CATALOGUE_SCHEMA, DNAComparator
from chronoclade.context_partitions import profile_records
from chronoclade.context_refinement import _allele, _profile
from chronoclade.datasets import load_dataset
from chronoclade.query_typing import _database_hash


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_csv(path, rows, fields=None):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def called_dna(hit, contigs):
    """Match the pinned caller's one-based, inclusive slicing and orientation."""
    dna = contigs[hit["contig"]][hit["contigStart"] - 1:hit["contigEnd"]]
    if hit["reverse"]:
        dna = str(Seq(dna).reverse_complement())
    if hashlib.sha1(dna.lower().encode()).hexdigest() != hit["checksum"]:
        raise ValueError("Extracted call DNA differs from native SHA1 checksum")
    return dna.upper()


def native_calls(root, sample_ids):
    result = {}
    for path in sorted((root / "sources/provider/native_typing_attempts").glob("*/query_typing.json")):
        assignments = json.loads(path.read_text())["assignments"]
        outputs = list(path.parent.glob("*/cgmlst.json"))
        if len(assignments) != 1 or len(outputs) != 1:
            raise ValueError("Expected one preserved native call per typing attempt")
        ident = assignments[0]["sample_id"]
        if ident in sample_ids:
            if ident in result:
                raise ValueError(f"Duplicate native output for {ident}")
            result[ident] = assignments[0], path, outputs[0]
    if set(result) != sample_ids:
        raise ValueError("Native call sample identities differ from labelled panel")
    return result


def build_catalogue(dataset_path, labels_path, output, reference_fastas):
    dataset_path = Path(dataset_path).resolve()
    dataset = load_dataset(dataset_path)
    rows = profile_records(dataset)
    with Path(labels_path).open() as stream:
        labels = list(csv.DictReader(stream))
    if len(labels) != len({row["sample_id"] for row in labels}):
        raise ValueError("Duplicate label sample identity")
    ids = {row["sample_id"] for row in labels}
    if {row["sample_id"] for row in rows} != ids:
        raise ValueError("Prepared identities must exactly match labelled panel")
    scope = {"species": rows[0]["species"]}
    for field in ("scheme", "scheme_version", "database_version", "database_sha256"):
        scope["scheme_id" if field == "scheme" else field] = rows[0]["cgmlst_" + field]
    reference = dataset.parameters["native_typing_reference"]
    config = Path(reference["config_path"])
    if sha256(config) != reference["config_sha256"]:
        raise ValueError("Native reference config checksum changed")
    settings = json.loads(config.read_text())["organisms"][scope["species"]]["cgmlst"]
    index = Path(settings["index_dir"])
    if _database_hash(index) != scope["database_sha256"]:
        raise ValueError("Full caller index fingerprint differs from prepared scope")
    evidence = native_calls(dataset_path.parent, ids)
    source_paths = {dataset_path, Path(labels_path).resolve(), config}
    source_paths.update(path for path in index.rglob("*") if path.is_file())
    sequences, audit, unresolved = {}, [], []
    known = set()
    with sqlite3.connect(f"file:{index / 'allele.db'}?mode=ro", uri=True) as connection:
        for row in rows:
            ident = row["sample_id"]
            assignment, query_path, calls_path = evidence[ident]
            source_paths.update((query_path, calls_path))
            assembly = dataset_path.parent / row["assembly_reference"]
            source_paths.add(assembly)
            if sha256(assembly) != assignment["assembly_sha256"]:
                raise ValueError(f"Assembly checksum changed: {ident}")
            contigs = {record.description: str(record.seq) for record in SeqIO.parse(assembly, "fasta")}
            calls = json.loads(calls_path.read_text())["alleles"]
            profile = _profile(row)
            for locus in row["cgmlst_loci"]:
                allele = _allele(profile.get(locus))
                if allele is None:
                    unresolved.append({"sample_id": ident, "locus": locus, "allele": "", "reason": "missing_call"})
                    continue
                hit = calls.get(locus)
                if hit is None or str(hit["id"]) != allele or hit["gene"] != locus:
                    raise ValueError(f"Prepared/native locus allele mapping differs: {ident} {locus}")
                dna = called_dna(hit, contigs)
                if set(dna) - set("ACGT"):
                    unresolved.append({"sample_id": ident, "locus": locus, "allele": allele, "reason": "ambiguous_sequence"})
                    continue
                key = locus, allele
                if key in sequences and sequences[key] != dna:
                    raise ValueError(f"One locus/allele has conflicting observed DNA: {key}")
                if hit["exact"]:
                    match = connection.execute(
                        "SELECT 1 FROM alleles WHERE hash=? AND gene=? AND st=? AND reverse=0",
                        (hit["checksum"], locus, int(allele)),
                    ).fetchone()
                    if not match:
                        raise ValueError(f"Native exact call lacks oriented SQLite mapping: {key}")
                    known.add(key)
                elif allele != hit["checksum"]:
                    raise ValueError("Novel call identifier must equal its observed DNA SHA1")
                sequences[key] = dna
                audit.append({"sample_id": ident, "locus": locus, "allele": allele,
                              "exact": hit["exact"], "length": len(dna),
                              "sequence_sha1_lowercase": hit["checksum"],
                              "sequence_sha256_uppercase": hashlib.sha256(dna.encode()).hexdigest(),
                              "native_call_path": str(calls_path), "assembly_path": str(assembly),
                              "contig": hit["contig"], "start": hit["contigStart"],
                              "end": hit["contigEnd"], "reverse": bool(hit["reverse"])})
    # Independently verify numeric identifiers against the actual frozen reference DNA.
    for locus in sorted({key[0] for key in known}):
        # Index FASTAs contain selected BLAST representatives, not every allele.
        fasta = Path(reference_fastas).resolve() / f"{locus}.fa.gz"
        source_paths.add(fasta)
        with gzip.open(fasta, "rt") as stream:
            reference_dna = {r.id: str(r.seq).upper() for r in SeqIO.parse(stream, "fasta")}
        for key in sorted(key for key in known if key[0] == locus):
            if reference_dna.get(key[1]) != sequences[key]:
                raise ValueError(f"Reference FASTA differs from observed exact call: {key}")
    output.mkdir(parents=True, exist_ok=True)
    dna_csv = output / "alleles.csv"
    write_csv(dna_csv, [{"locus": locus, "allele": allele, "dna": dna}
                        for (locus, allele), dna in sorted(sequences.items())])
    catalogue = output / "catalogue.json"
    write_json(catalogue, {"schema": CATALOGUE_SCHEMA, "schema_version": 1, "status": "complete",
                          "scope": scope, "artifact": {"path": dna_csv.name, "sha256": sha256(dna_csv)},
                          "coverage": "observed called alleles only; not a full scheme reference"})
    write_csv(output / "call-sequence-audit.csv", audit)
    write_csv(output / "unresolved.csv", unresolved, ["sample_id", "locus", "allele", "reason"])
    provenance = {"schema": "chronoclade.research.observed-locus-dna", "schema_version": 1,
                  "dataset_id": dataset.dataset_id, "scope": scope,
                  "mapping": "Actual assembly coordinates, caller orientation and lowercase SHA1; exact calls additionally checked against SQLite and reference FASTA",
                  "sources": [{"path": str(p), "sha256": sha256(p)} for p in sorted(source_paths)],
                  "observed_calls": len(audit), "catalogue_records": len(sequences),
                  "numeric_reference_alleles": len(known), "novel_observed_alleles": len(sequences) - len(known),
                  "unresolved": dict(Counter(row["reason"] for row in unresolved))}
    write_json(output / "provenance.json", provenance)
    return rows, catalogue, provenance


def distribution(values):
    return {"min": int(min(values)), "median": float(np.median(values)), "max": int(max(values))}


def correlation(rows, left, right):
    a, b = ([float(row[key]) for row in rows] for key in (left, right))
    return float(spearmanr(a, b).statistic) if len(set(a)) > 1 and len(set(b)) > 1 else None


def verify_concatenated(records, comparator, pairs, by_locus, output):
    """Independent vectorised SNP oracle on every observed fixed-length locus.

    These concatenated positions have the same positional assumption as the
    production comparator. This is an arithmetic check, not an alignment audit.
    Missing calls receive '?' and joint callable columns alone contribute.
    """
    lengths = {locus: {len(dna) for _, dna in alleles} for locus, alleles in by_locus.items()}
    fixed = {locus: next(iter(sizes)) for locus, sizes in lengths.items() if len(sizes) == 1}
    variable = set(lengths) - set(fixed)
    vectors = {}
    for ident, row in records.items():
        profile = _profile(row)
        joined = "".join(comparator.sequences.get((locus, _allele(profile.get(locus))), "?" * size)
                         for locus, size in sorted(fixed.items()))
        vectors[ident] = np.frombuffer(joined.encode(), dtype=np.uint8)
    audit = []
    for pair in pairs:
        a, b = (vectors[pair[key]] for key in ("sample_a", "sample_b"))
        mask = (a != ord("?")) & (b != ord("?"))
        difference, sites = int(np.count_nonzero((a != b) & mask)), int(mask.sum())
        profiles = [_profile(records[pair[key]]) for key in ("sample_a", "sample_b")]
        extra_difference = extra_sites = 0
        for locus in variable:
            x, y = [comparator.sequences.get((locus, _allele(profile.get(locus)))) for profile in profiles]
            if x is not None and y is not None and len(x) == len(y):
                extra_difference += sum(a != b for a, b in zip(x, y, strict=True))
                extra_sites += len(x)
        if difference + extra_difference != pair["locus_snps"] or sites + extra_sites != pair["compared_sites"]:
            raise ValueError("Concatenated fixed-locus arithmetic differs from DNAComparator")
        audit.append({"sample_a": pair["sample_a"], "sample_b": pair["sample_b"],
                      "concatenated_fixed_locus_snps": difference, "concatenated_callable_sites": sites,
                      "additional_variable_length_locus_snps": extra_difference,
                      "additional_variable_length_locus_sites": extra_sites,
                      "matches_comparator_after_documented_variable_locus_addition": True})
    write_csv(output / "concatenation-verification.csv", audit)
    return {"verified_pairs": len(audit), "fixed_length_loci": len(fixed),
            "fixed_length_concatenation_sites": sum(fixed.values()),
            "length_variable_loci": sorted(variable),
            "status": "all numerator and denominator checks passed",
            "limitation": "Arithmetic verification only; equal length does not exclude compensating indels or prove positional homology"}


def compare(dataset, labels, primary_pairs, snp_pairs, output, reference_fastas, min_loci=599):
    output = Path(output).resolve()
    records, catalogue, provenance = build_catalogue(dataset, labels, output, reference_fastas)
    records = {row["sample_id"]: row for row in records}
    with Path(primary_pairs).open() as stream:
        primary = list(csv.DictReader(stream))
    with Path(snp_pairs).open() as stream:
        snps = list(csv.DictReader(stream))
    def lookup(rows):
        indexed = {frozenset((p["sample_a"], p["sample_b"])): p for p in rows}
        if len(indexed) != len(rows) or any(len(k) != 2 for k in indexed):
            raise ValueError("Duplicate or self pair")
        return indexed
    primary_lookup, snp_lookup = lookup(primary), lookup(snps)
    if set(primary_lookup) != set(snp_lookup) or len(primary) != len(records) * (len(records) - 1) // 2:
        raise ValueError("Expected all matched pairs on the identical labelled panel")
    comparator = DNAComparator(catalogue, list(records.values()))
    pairs = []
    for pair in primary:
        left, right = (records[pair[key]] for key in ("sample_a", "sample_b"))
        scheme_size = len(left["cgmlst_loci"])
        if not 1 <= min_loci <= scheme_size:
            raise ValueError("Minimum sequence loci must lie within scheme size")
        # Obtain all evidence with overlap=0, then apply the declared absolute minimum.
        evidence, reason = comparator.pair(left, right, 0)
        if evidence is None:
            raise ValueError(f"Unexpected unavailable pair with zero overlap requirement: {reason}")
        snp = snp_lookup[frozenset((left["sample_id"], right["sample_id"]))]
        pairs.append({"sample_a": left["sample_id"], "sample_b": right["sample_id"],
                      "epi_label": pair["epi_label"], "allele_mismatches": int(pair["allele_mismatches"]),
                      "allele_shared_loci": int(pair["shared_loci"]),
                      "locus_snps": evidence["snp_differences"],
                      "compared_loci": evidence["shared_called_loci"], "compared_sites": evidence["shared_called_sites"],
                      "eligible": evidence["shared_called_loci"] >= min_loci,
                      "locus_mismatch_fraction": evidence["distance"],
                      "allele_mismatches_on_sequence_compared_loci": evidence["allele_differences"],
                      "missing_call_loci": len(evidence["missing_call_loci"]),
                      "missing_sequence_loci": len(evidence["missing_sequence_loci"]),
                      "unequal_length_loci": len(evidence["unequal_length_loci"]),
                      "clonal_snps": int(snp["clonal_snps"]), "core_callable_sites": int(snp["callable_sites"])})
    write_csv(output / "pairs.csv", pairs)
    comparator.save_audit(output)
    by_locus = defaultdict(list)
    for (locus, allele), dna in comparator.sequences.items():
        by_locus[locus].append((allele, dna))
    verification = verify_concatenated(records, comparator, pairs, by_locus, output)
    variable = [{"locus": locus, "observed_alleles": len(alleles),
                 "observed_lengths": sorted({len(dna) for _, dna in alleles}),
                 "distinct_sequences": len({dna for _, dna in alleles}),
                 "equal_length_allele_pair_substitutions": [
                     sum(x != y for x, y in zip(a, b, strict=True))
                     for (_, a), (_, b) in combinations(alleles, 2) if len(a) == len(b)]}
                for locus, alleles in sorted(by_locus.items()) if len(alleles) > 1]
    write_json(output / "variable-loci.json", variable)
    zero = [p for p in pairs if p["allele_mismatches"] == 0]
    summary = {"samples": len(records), "pairs": len(pairs), "eligible_pairs": sum(p["eligible"] for p in pairs),
               "minimum_compared_loci": min_loci, "scheme_loci": len(next(iter(records.values()))["cgmlst_loci"]),
               "catalogue_records": provenance["catalogue_records"], "observed_calls": provenance["observed_calls"],
               "unresolved": provenance["unresolved"], "variable_loci": len(variable),
               "length_variable_loci": [v for v in variable if len(v["observed_lengths"]) > 1],
               "distributions": [], "zero_allele_pairs": len(zero),
               "zero_allele_locus_snps": distribution([p["locus_snps"] for p in zero]) if zero else None,
               "zero_allele_core_snps": distribution([p["clonal_snps"] for p in zero]) if zero else None,
               "zero_allele_pairs_with_nonzero_core_snps": sum(p["clonal_snps"] > 0 for p in zero),
               "zero_locus_pairs_with_nonzero_original_alleles": sum(p["locus_snps"] == 0 and p["allele_mismatches"] > 0 for p in pairs),
               "pairs_with_length_exclusions": sum(p["unequal_length_loci"] > 0 for p in pairs),
               "unequal_length_pair_locus_exclusions": sum(p["unequal_length_loci"] for p in pairs),
               "all_retained_differing_alleles_have_one_substitution": all(
                   difference == 1 for locus in variable
                   for difference in locus["equal_length_allele_pair_substitutions"]),
               "independent_concatenation_verification": verification,
               "alignment_audit": "Not performed: MAFFT unavailable in current environment; equal length does not exclude compensating indels",
               "limitations": comparator.limitations,
               "source_files": [{"path": str(Path(p).resolve()), "sha256": sha256(p)} for p in (primary_pairs, snp_pairs, __file__)]}
    for label in ("all", "positive", "negative", "unknown"):
        subset = pairs if label == "all" else [p for p in pairs if p["epi_label"] == label]
        summary["distributions"].append({"label": label, "pairs": len(subset),
            **{key: distribution([p[key] for p in subset]) for key in ("allele_mismatches", "locus_snps", "clonal_snps", "compared_loci", "compared_sites", "unequal_length_loci")},
            "allele_locus_spearman": correlation(subset, "allele_mismatches", "locus_snps"),
            "allele_core_spearman": correlation(subset, "allele_mismatches", "clonal_snps"),
            "locus_core_spearman": correlation(subset, "locus_snps", "clonal_snps")})
    write_json(output / "summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default="tmp/science-zhou2017/prepared/dataset.json")
    parser.add_argument("--labels", default="validation/science/zhou2017/labels.csv")
    parser.add_argument("--primary-pairs", default="tmp/science-zhou2017/evaluation-primary/pairs.csv")
    parser.add_argument("--snp-pairs", default="tmp/science-zhou2017/evaluation-primary/allele-snp-pairs.csv")
    parser.add_argument("--output", default="tmp/science-zhou2017/locus-dna")
    parser.add_argument("--reference-fastas", default="tmp/typing-setup-smoke/databases/mlst-schemes/cgmlst_schemes/pubmlst_klebsiella_18")
    parser.add_argument("--min-loci", type=int, default=599)
    args = parser.parse_args()
    print(json.dumps(compare(args.dataset, args.labels, args.primary_pairs, args.snp_pairs, args.output, args.reference_fastas, args.min_loci), indent=2))


if __name__ == "__main__":
    main()
