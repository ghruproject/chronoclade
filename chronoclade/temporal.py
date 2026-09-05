"""Temporal-signal helpers."""

from __future__ import annotations

import json
import hashlib
import random
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from chronoclade.metadata import Sample

RATE = re.compile(r"--rate:\s*([+-]?[0-9.]+(?:e[+-]?\d+)?)", re.IGNORECASE)
RATE_STD = re.compile(
    r"--rate:\s*[+-]?[0-9.]+(?:e[+-]?\d+)?\s*\+/-\s*([0-9.]+(?:e[+-]?\d+)?)",
    re.IGNORECASE,
)
R_SQUARED = re.compile(r"--r\^2:\s*([+-]?[0-9.]+(?:e[+-]?\d+)?)", re.IGNORECASE)


class TemporalError(RuntimeError):
    """Raised when TreeTime output cannot be interpreted."""


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_clock(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8")
    rate = RATE.search(text)
    r_squared = R_SQUARED.search(text)
    if not rate or not r_squared:
        raise TemporalError(f"Could not parse TreeTime clock output: {path}")
    return {"rate": float(rate.group(1)), "r_squared": float(r_squared.group(1))}


def parse_clock_with_uncertainty(path: Path) -> dict[str, float]:
    """Parse TreeTime's fitted rate and its reported one-standard-deviation error."""

    result = parse_clock(path)
    rate_std = RATE_STD.search(path.read_text(encoding="utf-8"))
    if not rate_std:
        raise TemporalError(f"TreeTime did not report rate uncertainty: {path}")
    result["rate_std"] = float(rate_std.group(1))
    result["rate_lower_95"] = result["rate"] - 1.96 * result["rate_std"]
    result["rate_upper_95"] = result["rate"] + 1.96 * result["rate_std"]
    return result


def _write_dates(path: Path, samples: list[Sample], dates: list[str] | None = None) -> None:
    values = dates if dates is not None else [sample.collection_date for sample in samples]
    rows = ["sample_id,collection_date"]
    rows.extend(
        f"{sample.sample_id},{value}" for sample, value in zip(samples, values, strict=True)
    )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _run_randomised_clock(
    *,
    tree: Path,
    sequence_length: int,
    date_path: Path,
    run_dir: Path,
) -> dict[str, float] | None:
    """Run one independent TreeTime clock permutation."""

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
        "--allow-negative-rate",
        "--clock-filter",
        "0",
        "--outdir",
        str(run_dir),
        "--verbose",
        "0",
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    clock_path = run_dir / "molecular_clock.txt"
    if completed.returncode != 0 or not clock_path.is_file():
        return None
    return parse_clock(clock_path)


def run_date_randomisation(
    *,
    tree: Path,
    sequence_length: int,
    samples: list[Sample],
    observed_clock: Path,
    randomisations: int,
    randomisation_jobs: int,
    seed: int,
    output: Path,
) -> dict[str, object]:
    """Compare observed root-to-tip fit with date-permuted TreeTime fits."""

    if randomisation_jobs < 1:
        raise ValueError("randomisation_jobs must be at least 1")

    observed = parse_clock(observed_clock)
    random_metrics: list[dict[str, float]] = []
    rng = random.Random(seed)
    original_dates = [sample.collection_date for sample in samples]

    if randomisations > 0:
        with tempfile.TemporaryDirectory(prefix="date-randomisation-", dir=output.parent) as tmp:
            temporary = Path(tmp)
            jobs: list[tuple[Path, Path]] = []
            for index in range(randomisations):
                shuffled = original_dates.copy()
                rng.shuffle(shuffled)
                date_path = temporary / f"dates_{index:04d}.csv"
                run_dir = temporary / f"run_{index:04d}"
                _write_dates(date_path, samples, shuffled)
                jobs.append((date_path, run_dir))

            def run_job(job: tuple[Path, Path]) -> dict[str, float] | None:
                date_path, run_dir = job
                return _run_randomised_clock(
                    tree=tree,
                    sequence_length=sequence_length,
                    date_path=date_path,
                    run_dir=run_dir,
                )

            workers = min(randomisation_jobs, randomisations)
            with ThreadPoolExecutor(max_workers=workers) as executor:
                metrics = executor.map(run_job, jobs)
                random_metrics.extend(metric for metric in metrics if metric is not None)

    exceedances = sum(metric["r_squared"] >= observed["r_squared"] for metric in random_metrics)
    p_value = (exceedances + 1) / (len(random_metrics) + 1)
    result: dict[str, object] = {
        "method": "root_to_tip",
        "test_label": "Root-to-tip date-permutation screen",
        "seed": seed,
        "sequence_length": sequence_length,
        "tree_sha256": file_sha256(tree),
        "observed": observed,
        "requested_randomisations": randomisations,
        "successful_randomisations": len(random_metrics),
        "p_value_r_squared": p_value,
        "randomised": random_metrics,
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def _run_randomised_full_tree(
    *,
    tree: Path,
    sequence_length: int,
    date_path: Path,
    run_dir: Path,
) -> dict[str, float] | None:
    """Refit the complete TreeTime dating model for one date permutation."""

    command = [
        "treetime",
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
        "--time-marginal",
        "only-final",
        "--covariation",
        "--clock-filter",
        "0",
        "--outdir",
        str(run_dir),
        "--verbose",
        "0",
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    clock_path = run_dir / "molecular_clock.txt"
    if completed.returncode != 0 or not clock_path.is_file():
        return None
    try:
        return parse_clock_with_uncertainty(clock_path)
    except TemporalError:
        return None


def run_full_tree_date_randomisation(
    *,
    tree: Path,
    sequence_length: int,
    samples: list[Sample],
    observed_clock: Path,
    randomisations: int,
    randomisation_jobs: int,
    seed: int,
    output: Path,
) -> dict[str, object]:
    """Compare a full TreeTime fit with complete refits using permuted tip dates."""

    if randomisation_jobs < 1:
        raise ValueError("randomisation_jobs must be at least 1")
    observed = parse_clock_with_uncertainty(observed_clock)
    random_metrics: list[dict[str, float]] = []
    rng = random.Random(seed)
    original_dates = [sample.collection_date for sample in samples]

    if randomisations > 0:
        with tempfile.TemporaryDirectory(prefix="full-date-randomisation-", dir=output.parent) as tmp:
            temporary = Path(tmp)
            jobs: list[tuple[Path, Path]] = []
            for index in range(randomisations):
                shuffled = original_dates.copy()
                rng.shuffle(shuffled)
                date_path = temporary / f"dates_{index:04d}.csv"
                run_dir = temporary / f"run_{index:04d}"
                _write_dates(date_path, samples, shuffled)
                jobs.append((date_path, run_dir))

            def run_job(job: tuple[Path, Path]) -> dict[str, float] | None:
                date_path, run_dir = job
                return _run_randomised_full_tree(
                    tree=tree,
                    sequence_length=sequence_length,
                    date_path=date_path,
                    run_dir=run_dir,
                )

            workers = min(randomisation_jobs, randomisations)
            with ThreadPoolExecutor(max_workers=workers) as executor:
                metrics = executor.map(run_job, jobs)
                random_metrics.extend(metric for metric in metrics if metric is not None)

    observed_rate = observed["rate"]
    observed_lower = observed["rate_lower_95"]
    observed_upper = observed["rate_upper_95"]
    cr1_overlaps = [
        metric
        for metric in random_metrics
        if metric["rate_lower_95"] <= observed_rate <= metric["rate_upper_95"]
    ]
    cr2_overlaps = [
        metric
        for metric in random_metrics
        if max(observed_lower, metric["rate_lower_95"])
        <= min(observed_upper, metric["rate_upper_95"])
    ]
    result: dict[str, object] = {
        "method": "full_tree",
        "test_label": "Full TreeTime tip-date randomisation",
        "seed": seed,
        "sequence_length": sequence_length,
        "tree_sha256": file_sha256(tree),
        "criterion": "cr2",
        "observed": observed,
        "requested_randomisations": randomisations,
        "successful_randomisations": len(random_metrics),
        "cr1_passed": not cr1_overlaps,
        "cr2_passed": not cr2_overlaps,
        "cr1_overlapping_randomisations": len(cr1_overlaps),
        "cr2_overlapping_randomisations": len(cr2_overlaps),
        "randomised": random_metrics,
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result
