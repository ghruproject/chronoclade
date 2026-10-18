"""Adapt frozen EnteroBase allele files to the pinned Pathogenwatch caller.

The public allele catalogue does not include a validated HierCC background.
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

from Bio import SeqIO

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.enterobase_schemes import install_scheme
from chronoclade.progress import progress


_NORMALISATION = "numeric-id-acgt-v1"
_ALLELE_ID = re.compile(r"(?:.+[_-])?([0-9]+)")


def _normalise(source: Path, destination: Path) -> dict:
    """Preserve numbered ACGT alleles; record exclusions as upstream does."""
    seen = set()
    skipped = 0
    temporary = destination.with_suffix(destination.suffix + ".part")
    try:
        with gzip.open(source, "rt") as reader, gzip.open(temporary, "wt") as writer:
            for record in SeqIO.parse(reader, "fasta"):
                match = _ALLELE_ID.fullmatch(record.id)
                sequence = str(record.seq).upper()
                if match is None or not sequence or set(sequence) - set("ACGT"):
                    skipped += 1
                    continue
                allele = str(int(match.group(1)))
                if allele == "0" or allele in seen:
                    raise ValueError(f"Duplicate or missing-value allele ID in {source.name}")
                seen.add(allele)
                writer.write(f">{allele}\n{sequence}\n")
        if not seen:
            raise ValueError(f"No callable alleles in {source.name}")
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return {"called_alleles": len(seen), "excluded_alleles": skipped,
            "sha256": file_sha256(destination)}


def prepare_caller_scheme(root: Path, scheme_id: str, shortname: str) -> dict:
    """Download/resume checked definitions and create a separately checked view."""
    source_dir = root / "databases" / "enterobase" / scheme_id
    reference = install_scheme(scheme_id, source_dir)
    version = "snapshot:" + reference["reference_sha256"]
    scheme_dir = root / "databases" / "mlst-schemes" / "enterobase" / shortname
    scheme_dir.mkdir(parents=True, exist_ok=True)
    receipt_path = scheme_dir / "normalisation.json"
    previous = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
    if previous and (previous.get("reference_sha256") != reference["reference_sha256"]
                     or previous.get("method") != _NORMALISATION):
        raise ValueError("Existing caller allele view belongs to a different reference or method")
    old_files = {item["locus"]: item for item in previous.get("files", [])}
    normalized = []
    expanded = sum(item.get("uncompressed_bytes", 0) for item in reference["files"])
    progress("setup", f"{scheme_id}: normalizing {len(reference['loci']):,} loci "
             f"({expanded:,} expanded FASTA bytes; source retained compressed)")
    for entry in reference["files"]:
        locus = entry["locus"]
        source = source_dir / entry["path"]
        if file_sha256(source) != entry["sha256"]:
            raise ValueError(f"Frozen allele file changed: {locus}")
        target = scheme_dir / f"{locus}.fa.gz"
        old = old_files.get(locus)
        if old and target.is_file() and file_sha256(target) == old.get("sha256"):
            record = old
        else:
            record = {"locus": locus, "source_sha256": entry["sha256"],
                      **_normalise(source, target)}
        normalized.append(record)
        if len(normalized) % 100 == 0 or len(normalized) == len(reference["files"]):
            progress("setup", f"{scheme_id}: {len(normalized)}/{len(reference['files'])} caller loci normalized")
    write_json(receipt_path, {"schema": "chronoclade.enterobase-caller-view", "schema_version": 1,
                             "method": _NORMALISATION, "reference_sha256": reference["reference_sha256"],
                             "files": normalized})
    (scheme_dir / "README.md").write_text(
        f"# {scheme_id} caller allele view\n\n"
        f"Loci: {len(reference['loci']):,}. Expanded source FASTA bytes: {expanded:,}.\n"
        "The original compressed files remain under databases/enterobase/.\n"
        "This directory contains numbered ACGT caller alleles; normalization counts and hashes "
        "are in normalisation.json. Excluded records are not silently converted to valid calls.\n\n"
        f"Reference version: {version}.\n"
        "This is a public allele catalogue, not a validated HierCC background. "
        "Full typing readiness is recorded separately by setup.\n"
    )
    write_json(scheme_dir / "metadata.json", {"genes": reference["loci"], "last_updated": version})
    return {"shortname": shortname, "host": "enterobase", "host_path": scheme_id,
            "scheme_id": scheme_id, "type": "cgmlst", "cite": [],
            "name": scheme_id, "targets": [{"name": "Escherichia coli", "taxid": 562}],
            "db_path": str(scheme_dir), "last_updated": version,
            "reference_sha256": reference["reference_sha256"],
            "reference_manifest": str(source_dir / "reference.json"),
            "normalisation_sha256": file_sha256(receipt_path)}
