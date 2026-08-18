from pathlib import Path

import pytest

from beyondmlst.metadata import (
    MetadataError,
    Sample,
    assembly_n50,
    read_metadata,
    select_reference,
)


def write_fasta(path: Path, lengths: list[int]) -> None:
    records = []
    for index, length in enumerate(lengths):
        records.append(f">contig_{index}\n{'A' * length}\n")
    path.write_text("".join(records), encoding="utf-8")


def test_reads_relative_assemblies_and_partial_dates(tmp_path: Path) -> None:
    assembly = tmp_path / "sample.fasta"
    write_fasta(assembly, [10])
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "sample_id,assembly,collection_date,location,species,lineage,origin\n"
        "S1,sample.fasta,2024-03,KIMS,E_coli,ST131,local\n",
        encoding="utf-8",
    )

    samples = read_metadata(metadata)

    assert samples[0].assembly == assembly
    assert samples[0].collection_date == "2024-03"


def test_rejects_duplicate_identifiers(tmp_path: Path) -> None:
    assembly = tmp_path / "sample.fasta"
    write_fasta(assembly, [10])
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "sample_id,assembly,collection_date,location,species,lineage,origin\n"
        "S1,sample.fasta,2024,KIMS,E_coli,ST131,local\n"
        "S1,sample.fasta,2025,KIMS,E_coli,ST131,local\n",
        encoding="utf-8",
    )

    with pytest.raises(MetadataError, match="duplicate sample_id"):
        read_metadata(metadata)


def test_n50_and_reference_selection(tmp_path: Path) -> None:
    fragmented = tmp_path / "fragmented.fasta"
    contiguous = tmp_path / "contiguous.fasta"
    write_fasta(fragmented, [5, 5, 5, 5])
    write_fasta(contiguous, [12, 8])
    samples = [
        Sample("A", fragmented, "2020", "X", "E_coli", "ST1", "local"),
        Sample("B", contiguous, "2021", "X", "E_coli", "ST1", "local"),
    ]

    assert assembly_n50(fragmented) == 5
    assert assembly_n50(contiguous) == 12
    assert select_reference(samples).sample_id == "B"


def test_explicit_reference_wins(tmp_path: Path) -> None:
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"
    write_fasta(first, [5])
    write_fasta(second, [20])
    samples = [
        Sample("A", first, "2020", "X", "E_coli", "ST1", "local", True),
        Sample("B", second, "2021", "X", "E_coli", "ST1", "local"),
    ]

    assert select_reference(samples).sample_id == "A"
