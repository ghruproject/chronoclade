from pathlib import Path
from types import SimpleNamespace

import pytest

from chronoclade import recombination
from chronoclade.evidence import read_alignment


def test_phipack_profile_respects_reference_contig_boundaries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference = tmp_path / "reference.fasta"
    reference.write_text(">one\n" + "A" * 1_200 + "\n>two\n" + "A" * 1_100 + "\n")
    alignment = tmp_path / "alignment.fasta"
    alignment.write_text(">A\n" + "A" * 2_300 + "\n>B\n" + "A" * 2_300 + "\n")

    observed_lengths: list[int] = []

    def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        run_dir = Path(str(kwargs["cwd"]))
        observed_lengths.append(int(command[command.index("-n") + 1]))
        (run_dir / "Profile.csv").write_text("50,0.005\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(recombination.subprocess, "run", fake_run)
    summary = recombination.run_phipack_filter(
        alignment=alignment,
        reference=reference,
        output_alignment=tmp_path / "filtered.fasta",
        profile_output=tmp_path / "profile.tsv",
        regions_output=tmp_path / "regions.tsv",
        summary_output=tmp_path / "summary.json",
        threads=2,
        force=False,
    )

    assert sorted(observed_lengths) == [1_100, 1_200]
    assert summary["recombination_detected"] is True
    assert summary["recombination_regions"] == 2
    assert summary["tested_core_blocks"] == 2
    assert summary["masked_alignment_sites"] == 200
    filtered = read_alignment(tmp_path / "filtered.fasta")
    assert filtered["A"][:100] == "N" * 100
    assert filtered["A"][1_200:1_300] == "N" * 100
    regions = (tmp_path / "regions.tsv").read_text(encoding="utf-8")
    assert "1\t0\t100\t0\t100" in regions
    assert "2\t0\t100\t1200\t1300" in regions


def test_phipack_profile_does_not_treat_missing_data_as_a_segment_join(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference = tmp_path / "reference.fasta"
    reference.write_text(">one\n" + "A" * 2_500 + "\n", encoding="utf-8")
    alignment = tmp_path / "alignment.fasta"
    alignment.write_text(
        ">A\n" + "A" * 2_500 + "\n>B\n" + "A" * 1_200 + "N" * 100 + "A" * 1_200 + "\n",
        encoding="utf-8",
    )

    observed_lengths: list[int] = []

    def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        run_dir = Path(str(kwargs["cwd"]))
        observed_lengths.append(int(command[command.index("-n") + 1]))
        (run_dir / "Profile.csv").write_text("50,0.5\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(recombination.subprocess, "run", fake_run)
    summary = recombination.run_phipack_filter(
        alignment=alignment,
        reference=reference,
        output_alignment=tmp_path / "filtered.fasta",
        profile_output=tmp_path / "profile.tsv",
        regions_output=tmp_path / "regions.tsv",
        summary_output=tmp_path / "summary.json",
        threads=2,
        force=False,
    )

    assert observed_lengths == [2_500]
    assert summary["tested_core_blocks"] == 1
    assert summary["recombination_detected"] is False


def test_profile_blocks_are_bounded_without_crossing_contigs() -> None:
    chunks = recombination._profile_chunks([600_000, 300_000], {"A": ""})

    assert [(chunk.contig, chunk.window_start, chunk.window_end) for chunk in chunks] == [
        (1, 0, 250_000),
        (1, 250_000, 500_000),
        (1, 500_000, 600_000),
        (2, 0, 250_000),
        (2, 250_000, 300_000),
    ]
