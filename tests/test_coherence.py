import csv
from pathlib import Path

from chronoclade.coherence import screen_alignment
from chronoclade.metadata import Sample


def _samples(tmp_path: Path, count: int) -> list[Sample]:
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig\nAAAA\n", encoding="utf-8")
    return [
        Sample(f"S{index}", assembly, "2020", "X", "E_coli", "ST1", "local")
        for index in range(count)
    ]


def test_coherence_screen_flags_one_extreme_outlier(tmp_path: Path) -> None:
    samples = _samples(tmp_path, 10)
    alignment = tmp_path / "alignment.fasta"
    records = []
    for index, sample in enumerate(samples):
        sequence = list("A" * 2000)
        sequence[index] = "C"
        if sample.sample_id == "S9":
            sequence[:400] = "C" * 400
        records.append(f">{sample.sample_id}\n{''.join(sequence)}\n")
    alignment.write_text("".join(records), encoding="utf-8")

    result = screen_alignment(alignment, samples, output=tmp_path / "coherence.tsv")

    assert result["status"] == "failed"
    assert result["flagged_samples"] == ["S9"]
    with (tmp_path / "coherence.tsv").open(newline="", encoding="utf-8") as handle:
        rows = {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["S9"]["flagged"] == "true"
    assert float(rows["S9"]["ratio_to_cohort_median"]) > 100


def test_coherence_screen_does_not_reject_two_balanced_groups(tmp_path: Path) -> None:
    samples = _samples(tmp_path, 10)
    alignment = tmp_path / "alignment.fasta"
    records = []
    for index, sample in enumerate(samples):
        background = "A" if index < 5 else "C"
        records.append(f">{sample.sample_id}\n{background * 2000}\n")
    alignment.write_text("".join(records), encoding="utf-8")

    result = screen_alignment(alignment, samples, output=tmp_path / "coherence.tsv")

    assert result["status"] == "passed"
    assert result["flagged_samples"] == []
