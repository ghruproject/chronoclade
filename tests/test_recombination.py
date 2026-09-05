from pathlib import Path
from types import SimpleNamespace

import pytest

from chronoclade import recombination


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
    summary = recombination.run_phipack_screen(
        alignment=alignment,
        reference=reference,
        profile_output=tmp_path / "profile.tsv",
        significant_blocks_output=tmp_path / "significant_blocks.tsv",
        summary_output=tmp_path / "summary.json",
        threads=2,
        force=False,
    )

    assert sorted(observed_lengths) == [1_100, 1_200]
    assert summary["recombination_detected"] is True
    assert summary["significant_blocks"] == 2
    assert summary["tested_core_blocks"] == 2
    assert "no alignment sites were masked" in summary["localisation_limit"]
    assert "without a multiple-testing correction" in summary["multiple_testing_limit"]
    blocks = (tmp_path / "significant_blocks.tsv").read_text(encoding="utf-8")
    assert "1\t0\t1200\t0.005\t1" in blocks
    assert "2\t0\t1100\t0.005\t1" in blocks


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
    summary = recombination.run_phipack_screen(
        alignment=alignment,
        reference=reference,
        profile_output=tmp_path / "profile.tsv",
        significant_blocks_output=tmp_path / "significant_blocks.tsv",
        summary_output=tmp_path / "summary.json",
        threads=2,
        force=False,
    )

    assert observed_lengths == [2_500]
    assert summary["tested_core_blocks"] == 1
    assert summary["recombination_detected"] is False


def test_profile_blocks_are_bounded_without_crossing_contigs() -> None:
    chunks = recombination._profile_chunks([600_000, 300_000])

    assert [(chunk.contig, chunk.window_start, chunk.window_end) for chunk in chunks] == [
        (1, 0, 250_000),
        (1, 250_000, 500_000),
        (1, 500_000, 600_000),
        (2, 0, 250_000),
        (2, 250_000, 300_000),
    ]
