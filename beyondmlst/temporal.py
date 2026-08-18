"""Temporal-signal helpers."""

from __future__ import annotations

import json
import random
import re
import subprocess
import tempfile
from pathlib import Path

from beyondmlst.metadata import Sample

RATE = re.compile(r"--rate:\s*([+-]?[0-9.]+(?:e[+-]?\d+)?)", re.IGNORECASE)
R_SQUARED = re.compile(r"--r\^2:\s*([+-]?[0-9.]+(?:e[+-]?\d+)?)", re.IGNORECASE)


class TemporalError(RuntimeError):
    """Raised when TreeTime output cannot be interpreted."""


def parse_clock(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8")
    rate = RATE.search(text)
    r_squared = R_SQUARED.search(text)
    if not rate or not r_squared:
        raise TemporalError(f"Could not parse TreeTime clock output: {path}")
    return {"rate": float(rate.group(1)), "r_squared": float(r_squared.group(1))}


def _write_dates(path: Path, samples: list[Sample], dates: list[str] | None = None) -> None:
    values = dates if dates is not None else [sample.collection_date for sample in samples]
    rows = ["sample_id,collection_date"]
    rows.extend(f"{sample.sample_id},{value}" for sample, value in zip(samples, values, strict=True))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def run_date_randomisation(
    *,
    tree: Path,
    sequence_length: int,
    samples: list[Sample],
    observed_clock: Path,
    randomisations: int,
    seed: int,
    output: Path,
) -> dict[str, object]:
    """Compare observed root-to-tip fit with date-permuted TreeTime fits."""

    observed = parse_clock(observed_clock)
    random_metrics: list[dict[str, float]] = []
    rng = random.Random(seed)
    original_dates = [sample.collection_date for sample in samples]

    if randomisations > 0:
        with tempfile.TemporaryDirectory(prefix="date-randomisation-", dir=output.parent) as tmp:
            temporary = Path(tmp)
            for index in range(randomisations):
                shuffled = original_dates.copy()
                rng.shuffle(shuffled)
                date_path = temporary / f"dates_{index:04d}.csv"
                run_dir = temporary / f"run_{index:04d}"
                _write_dates(date_path, samples, shuffled)
                command = [
                    "treetime",
                    "clock",
                    "--tree",
                    str(tree),
                    "--dates",
                    str(date_path),
                    "--name-column",
                    "sample_id",
                    "--date-column",
                    "collection_date",
                    "--sequence-length",
                    str(sequence_length),
                    "--reroot",
                    "least-squares",
                    "--outdir",
                    str(run_dir),
                    "--verbose",
                    "0",
                ]
                completed = subprocess.run(command, capture_output=True, text=True, check=False)
                clock_path = run_dir / "molecular_clock.txt"
                if completed.returncode == 0 and clock_path.is_file():
                    random_metrics.append(parse_clock(clock_path))

    exceedances = sum(
        metric["r_squared"] >= observed["r_squared"] for metric in random_metrics
    )
    p_value = (exceedances + 1) / (len(random_metrics) + 1)
    result: dict[str, object] = {
        "observed": observed,
        "requested_randomisations": randomisations,
        "successful_randomisations": len(random_metrics),
        "p_value_r_squared": p_value,
        "randomised": random_metrics,
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result
