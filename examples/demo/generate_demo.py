#!/usr/bin/env python3
"""Generate a deterministic two-lineage beyondMLST demonstration dataset."""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

GENOME_LENGTH = 500_000
SAMPLES_PER_LINEAGE = 12
MUTATIONS_PER_STEP = 1
BASES = "ACGT"


def write_fasta(path: Path, sample_id: str, sequence: str) -> None:
    lines = [f">{sample_id}"]
    lines.extend(sequence[index : index + 80] for index in range(0, len(sequence), 80))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def mutate(sequence: list[str], positions: list[int]) -> None:
    alternatives = {"A": "C", "C": "G", "G": "T", "T": "A"}
    for position in positions:
        sequence[position] = alternatives[sequence[position]]


def generate(output: Path) -> Path:
    output = output.expanduser().resolve()
    assemblies = output / "assemblies"
    assemblies.mkdir(parents=True, exist_ok=True)

    rng = random.Random(20260818)
    reference = [rng.choice(BASES) for _ in range(GENOME_LENGTH)]
    positions = list(range(100, GENOME_LENGTH - 100))
    rng.shuffle(positions)
    positions_per_lineage = sum(range(1, SAMPLES_PER_LINEAGE + 1)) * MUTATIONS_PER_STEP
    signal_positions = positions[:positions_per_lineage]
    null_positions = positions[positions_per_lineage : 2 * positions_per_lineage]

    chronological_dates = list(range(2012, 2024))
    shuffled_dates = [2018, 2012, 2022, 2015, 2020, 2013, 2023, 2017, 2014, 2021, 2016, 2019]
    lineages = [
        ("ST_SIGNAL", chronological_dates, signal_positions),
        ("ST_NO_SIGNAL", shuffled_dates, null_positions),
    ]

    rows: list[dict[str, str]] = []
    for lineage, dates, mutation_positions in lineages:
        position_cursor = 0
        for index in range(SAMPLES_PER_LINEAGE):
            # Give each isolate its own branch from a common ancestor. Branch
            # length grows with sample index, producing a positive clock only
            # when those indices retain chronological dates.
            sequence = reference.copy()
            mutation_count = (index + 1) * MUTATIONS_PER_STEP
            mutate(
                sequence,
                mutation_positions[position_cursor : position_cursor + mutation_count],
            )
            position_cursor += mutation_count
            sample_id = f"{lineage}_{index + 1:02d}"
            assembly = assemblies / f"{sample_id}.fasta"
            write_fasta(assembly, sample_id, "".join(sequence))
            local = index < 8
            rows.append(
                {
                    "sample_id": sample_id,
                    "assembly": str(assembly.relative_to(output)),
                    "collection_date": f"{dates[index]}-06-15",
                    "location": "KIMS" if local else "Regional_context",
                    "species": "Escherichia_coli",
                    "lineage": lineage,
                    "origin": "local" if local else "context",
                    "is_reference": "true" if index == 0 else "false",
                    "patient_id": f"DEMO_{lineage}_{index + 1:02d}",
                }
            )

    metadata = output / "metadata.csv"
    with metadata.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("demo_run/input"),
        help="Directory for generated assemblies and metadata",
    )
    args = parser.parse_args()
    metadata = generate(args.output)
    print(f"Generated demonstration metadata: {metadata}")


if __name__ == "__main__":
    main()
