import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from beyondmlst.metadata import Sample
from beyondmlst import temporal
from beyondmlst.temporal import TemporalError, parse_clock


def test_parse_treetime_clock(tmp_path: Path) -> None:
    output = tmp_path / "molecular_clock.txt"
    output.write_text(
        "Root-Tip-Regression:\n --rate:\t1.234e-06\n --r^2:  \t0.67\n",
        encoding="utf-8",
    )

    assert parse_clock(output) == {"rate": 1.234e-06, "r_squared": 0.67}


def test_rejects_unrecognised_clock_output(tmp_path: Path) -> None:
    output = tmp_path / "molecular_clock.txt"
    output.write_text("not a clock\n", encoding="utf-8")

    with pytest.raises(TemporalError):
        parse_clock(output)


def test_date_randomisations_run_with_bounded_concurrency(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    observed = tmp_path / "molecular_clock.txt"
    observed.write_text("--rate: 1e-6\n--r^2: 0.5\n", encoding="utf-8")
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig\nAAAA\n", encoding="utf-8")
    samples = [
        Sample(f"S{index}", assembly, str(2020 + index), "KIMS", "E_coli", "ST1", "local")
        for index in range(4)
    ]
    active = 0
    maximum_active = 0
    lock = threading.Lock()

    def fake_randomised_clock(**_: object) -> dict[str, float]:
        nonlocal active, maximum_active
        with lock:
            active += 1
            maximum_active = max(maximum_active, active)
        time.sleep(0.02)
        with lock:
            active -= 1
        return {"rate": 1e-6, "r_squared": 0.2}

    monkeypatch.setattr(temporal, "_run_randomised_clock", fake_randomised_clock)
    result = temporal.run_date_randomisation(
        tree=tmp_path / "tree.nwk",
        sequence_length=4,
        samples=samples,
        observed_clock=observed,
        randomisations=6,
        randomisation_jobs=3,
        seed=42,
        output=tmp_path / "temporal_signal.json",
    )

    assert result["successful_randomisations"] == 6
    assert 1 < maximum_active <= 3


def test_randomised_clock_allows_negative_null_rates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: list[str] = []

    def fake_run(command: list[str], **_: object) -> SimpleNamespace:
        captured.extend(command)
        run_dir = Path(command[command.index("--outdir") + 1])
        run_dir.mkdir()
        (run_dir / "molecular_clock.txt").write_text(
            "--rate: -1e-7\n--r^2: 0.1\n", encoding="utf-8"
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(temporal.subprocess, "run", fake_run)
    result = temporal._run_randomised_clock(
        tree=tmp_path / "tree.nwk",
        sequence_length=100,
        date_path=tmp_path / "dates.csv",
        run_dir=tmp_path / "randomised",
    )

    assert "--allow-negative-rate" in captured
    assert captured[captured.index("--clock-filter") + 1] == "0"
    assert result == {"rate": -1e-7, "r_squared": 0.1}
