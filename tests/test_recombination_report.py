from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from chronoclade.recombination_report import (
    RecombinationReportError,
    read_importation_intervals,
    write_recombination_evidence,
)


def test_recombination_evidence_separates_imports_from_incomplete_columns(
    tmp_path: Path,
) -> None:
    alignment = tmp_path / "alignment.fasta"
    alignment.write_text(">A\nACGTNNACGT\n>B\nACGTAAACGT\n", encoding="utf-8")
    filtered = tmp_path / "filtered.fasta"
    # Positions 3-4 are imported; 5-6 are incomplete in at least one sequence.
    filtered.write_text(">A\nACGTAC\n>B\nACGTAC\n", encoding="utf-8")
    importations = tmp_path / "imports.txt"
    importations.write_text("Node\tBeg\tEnd\nNODE_1\t3\t4\n", encoding="utf-8")

    reference = tmp_path / "reference.fasta"
    reference.write_text(">chr1\nACGTA\n>plasmid\nAACGT\n", encoding="utf-8")
    summary = write_recombination_evidence(
        alignment=alignment,
        filtered_alignment=filtered,
        importations=importations,
        output_directory=tmp_path,
        reference=reference,
        bins=5,
    )

    assert summary["inferred_importation_intervals"] == 1
    assert summary["branches_with_inferred_importation"] == 1
    assert summary["recombination_removed_sites"] == 2
    assert summary["incomplete_only_removed_sites"] == 2
    assert summary["retained_clonal_sites"] == 6
    assert summary["longest_inferred_intervals"] == [
        {
            "node": "NODE_1",
            "alignment_start": 3,
            "alignment_end": 4,
            "length": 2,
            "reference_segments": ["chr1:3-4"],
            "crosses_reference_record_boundary": False,
        }
    ]
    assert summary["boundary_crossing_intervals"] == 0
    assert (tmp_path / "recombination_map.svg").is_file()
    assert (tmp_path / "recombination_map.png").is_file()
    assert json.loads((tmp_path / "recombination_summary.json").read_text())["sequences"] == 2
    with (tmp_path / "recombination_genome_profile.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 5
    assert sum(int(row["recombination_removed_sites"]) for row in rows) == 2
    assert sum(int(row["incomplete_only_removed_sites"]) for row in rows) == 2


def test_cross_record_importation_is_flagged(tmp_path: Path) -> None:
    alignment = tmp_path / "alignment.fasta"
    alignment.write_text(">A\nACGTACGT\n>B\nACGTACGT\n", encoding="utf-8")
    filtered = tmp_path / "filtered.fasta"
    filtered.write_text(">A\nACGT\n>B\nACGT\n", encoding="utf-8")
    importations = tmp_path / "imports.txt"
    importations.write_text("Node\tBeg\tEnd\nNODE_1\t3\t6\n", encoding="utf-8")
    reference = tmp_path / "reference.fasta"
    reference.write_text(">one\nACGT\n>two\nACGT\n", encoding="utf-8")

    summary = write_recombination_evidence(
        alignment=alignment,
        filtered_alignment=filtered,
        importations=importations,
        output_directory=tmp_path,
        reference=reference,
    )

    assert summary["boundary_crossing_intervals"] == 1
    assert summary["longest_inferred_intervals"][0]["reference_segments"] == [
        "one:3-4",
        "two:1-2",
    ]


def test_importation_coordinates_are_validated(tmp_path: Path) -> None:
    path = tmp_path / "imports.txt"
    path.write_text("Node\tBeg\tEnd\nNODE_1\t0\t4\n", encoding="utf-8")

    with pytest.raises(RecombinationReportError, match="Invalid importation interval"):
        read_importation_intervals(path)


def test_filtered_alignment_must_match_the_documented_filter(tmp_path: Path) -> None:
    alignment = tmp_path / "alignment.fasta"
    alignment.write_text(">A\nACGT\n>B\nACGT\n", encoding="utf-8")
    filtered = tmp_path / "filtered.fasta"
    filtered.write_text(">A\nACG\n>B\nACG\n", encoding="utf-8")
    importations = tmp_path / "imports.txt"
    importations.write_text("Node\tBeg\tEnd\n", encoding="utf-8")

    with pytest.raises(RecombinationReportError, match="filtered-alignment length"):
        write_recombination_evidence(
            alignment=alignment,
            filtered_alignment=filtered,
            importations=importations,
            output_directory=tmp_path,
        )
