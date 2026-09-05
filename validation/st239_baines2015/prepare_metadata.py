#!/usr/bin/env python3
"""Create ChronoClade metadata for the recorded ST239 accession set."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def find_assembly(directory: Path, accession: str) -> Path:
    matches = sorted(directory.glob(f"{accession}.*"))
    fasta = [path for path in matches if ".fa" in path.name or ".fasta" in path.name]
    if len(fasta) != 1:
        raise ValueError(f"Expected one assembly for {accession}; found {len(fasta)}")
    return fasta[0].resolve()


def write_metadata(manifest: Path, assemblies: Path, output: Path) -> None:
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    fieldnames = [
        "sample_id",
        "assembly",
        "collection_date",
        "location",
        "species",
        "lineage",
        "origin",
        "is_reference",
        "patient_id",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "sample_id": row["sample_id"],
                    "assembly": find_assembly(assemblies, row["sample_id"]),
                    "collection_date": row["collection_date"],
                    "location": row["country"],
                    "species": "Staphylococcus_aureus",
                    "lineage": "ST239",
                    "origin": row["origin"],
                    "is_reference": "false",
                    "patient_id": "",
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("validation/st239_baines2015/input_manifest.csv"),
    )
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_metadata(args.manifest, args.assemblies, args.output)


if __name__ == "__main__":
    main()
