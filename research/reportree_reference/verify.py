#!/usr/bin/env python3
"""Compare a fresh ReporTree Neisseria run with the checked-in example.

Run from any directory with the task's Python environment:
    python research/reportree_reference/verify.py

The script writes tmp/reportree-reference/verification.json. It does not run
ReporTree; it waits for execution.json to record a successful run.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def read_tsv(path: Path, *, comment_lines: bool = False) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        if comment_lines:
            lines = (line for line in handle if not line.startswith("#"))
            reader = csv.DictReader(lines, delimiter="\t")
        else:
            reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def stable_regions(path: Path) -> dict[str, tuple[str, str, int]]:
    # The data header itself starts with #block_id, after two prose comments.
    with path.open(encoding="utf-8", newline="") as handle:
        lines = iter(handle)
        next(lines, None)
        next(lines, None)
        reader = csv.DictReader(lines, delimiter="\t")
        rows = list(reader)
    regions: dict[str, tuple[str, str, int]] = {}
    for row in rows:
        block_id = row.get("#block_id", row.get("block_id", ""))
        if not block_id:
            continue
        regions[block_id] = (
            row["first_partition"],
            row["last_partition"],
            int(row["len_block"]),
        )
    return regions


def metrics(path: Path) -> dict[tuple[str, str], float]:
    _, rows = read_tsv(path)
    return {
        (row["col2"], row["col1"]): float(row["AdjWallace"])
        for row in rows
    }


def membership_map(path: Path) -> tuple[list[str], dict[str, dict[str, tuple[str, ...]]]]:
    _, rows = read_tsv(path)
    if not rows:
        raise ValueError(f"No rows in {path}")
    id_column = "sequence" if "sequence" in rows[0] else "#ID"
    partition_columns = [
        name for name in rows[0]
        if re.fullmatch(r"MST-\d+x\d+(?:\.\d+)?", name)
    ]
    result: dict[str, dict[str, tuple[str, ...]]] = {}
    for column in partition_columns:
        groups: dict[str, list[str]] = {}
        for row in rows:
            groups.setdefault(row[column], []).append(row[id_column])
        result[column] = {
            label: tuple(sorted(samples)) for label, samples in groups.items()
        }
    return partition_columns, result


def partition_summary(path: Path) -> dict[tuple[str, tuple[str, ...]], dict[str, str]]:
    _, rows = read_tsv(path)
    result: dict[tuple[str, tuple[str, ...]], dict[str, str]] = {}
    for row in rows:
        if not row.get("partition"):
            continue
        members = tuple(sorted(filter(None, row.get("samples", "").split(","))))
        result[(row["partition"], members)] = row
    return result


def normalized_country_summary(value: str) -> tuple[tuple[str, str], ...]:
    """Compare country percentages independent of category ordering."""
    entries = re.findall(r"(.+?)\s+\((\d+(?:\.\d+)?)%\)", value)
    return tuple(sorted((country.strip().lstrip(",").strip(), percentage) for country, percentage in entries))


def group_summary_diffs(reference: Path, fresh: Path) -> dict[str, Any]:
    ref = partition_summary(reference)
    new = partition_summary(fresh)
    identity_union = set(ref) | set(new)
    missing = sorted(identity_union - set(new))
    extra = sorted(identity_union - set(ref))
    fields = (
        "cluster_length", "first_seq_date", "last_seq_date", "timespan_days",
        "country", "n_country",
    )
    differences: list[dict[str, Any]] = []
    for key in sorted(set(ref) & set(new)):
        unequal = {
            field: {"reference": ref[key].get(field, ""), "fresh": new[key].get(field, "")}
            for field in fields
            if (
                normalized_country_summary(ref[key].get(field, ""))
                != normalized_country_summary(new[key].get(field, ""))
                if field == "country"
                else ref[key].get(field, "") != new[key].get(field, "")
            )
        }
        if unequal:
            differences.append({"partition": key[0], "samples": list(key[1]), "fields": unequal})
    return {
        "reference_group_count": len(ref),
        "fresh_group_count": len(new),
        "missing_groups_by_sample_identity": len(missing),
        "extra_groups_by_sample_identity": len(extra),
        "summary_field_mismatch_count": len(differences),
        "summary_field_mismatches": differences[:100],
        "missing_group_examples": [
            {"partition": p, "samples": list(samples)} for p, samples in missing[:20]
        ],
        "extra_group_examples": [
            {"partition": p, "samples": list(samples)} for p, samples in extra[:20]
        ],
        "fields_compared": list(fields),
    }


def execution_return_code(data: dict[str, Any]) -> int | None:
    for key in ("returncode", "return_code", "exit_code", "exitCode"):
        value = data.get(key)
        if isinstance(value, int):
            return value
    return None


def version_from_log(path: Path) -> str | None:
    if not path.exists():
        return None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = re.search(r"\bversion\s+([^\s]+)\s+last updated on\s+([^\s]+)", line)
            if match:
                return f"{match.group(1)} (last updated {match.group(2)})"
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--task-dir", type=Path,
        default=Path(__file__).resolve().parents[2] / "tmp" / "reportree-reference",
        help="Directory containing upstream/ and reproduced/ (default: this checkout's tmp/reportree-reference)",
    )
    args = parser.parse_args()
    task_dir = args.task_dir.resolve()
    reproduced = task_dir / "reproduced"
    reference = task_dir / "upstream" / "examples" / "Neisseria" / "output"
    execution_path = reproduced / "execution.json"
    report_path = task_dir / "verification.json"

    result: dict[str, Any] = {
        "schema_version": 1,
        "reference": str(reference),
        "reproduced": str(reproduced),
        "execution": {"path": str(execution_path), "successful": False},
        "source_versions": {
            "reference_example": version_from_log(reference / "NG_822.log"),
            "fresh_source": version_from_log(reproduced / "NG_822.log"),
        },
        "checks": {},
        "overall": "NOT_RUN",
    }

    if not execution_path.exists():
        result["execution"]["reason"] = "execution.json is absent; fresh run completion is unconfirmed"
        report_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(f"Verification not run: missing {execution_path}")
        return 2

    execution = read_json(execution_path)
    code = execution_return_code(execution)
    expected_outputs_present = execution.get("expected_outputs_present") is True
    required_outputs = (
        "NG_822.log", "NG_822_stableRegions.tsv", "NG_822_metrics.tsv",
        "NG_822_metadata_w_partitions.tsv", "NG_822_partitions_summary.tsv",
    )
    missing_outputs = [name for name in required_outputs if not (reproduced / name).is_file()]
    result["execution"].update({
        "return_code": code,
        "expected_outputs_present": expected_outputs_present,
        "missing_required_outputs": missing_outputs,
        "record": execution,
    })
    if code != 0 or not expected_outputs_present or missing_outputs:
        result["execution"]["reason"] = (
            "run is not verified complete: require return code 0, "
            "expected_outputs_present true, and all comparison inputs present"
        )
        report_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(
            "Verification not run: "
            f"return_code={code!r}, expected_outputs_present={expected_outputs_present}, "
            f"missing_outputs={missing_outputs}"
        )
        return 2
    result["execution"]["successful"] = True

    # Stability regions: compare endpoints and run length, independent of block IDs.
    ref_regions = stable_regions(reference / "NG_822_stableRegions.tsv")
    new_regions = stable_regions(reproduced / "NG_822_stableRegions.tsv")
    ref_region_values = Counter(ref_regions.values())
    new_region_values = Counter(new_regions.values())
    region_match = ref_region_values == new_region_values
    result["checks"]["stability_regions"] = {
        "passed": region_match,
        "reference_count": sum(ref_region_values.values()),
        "fresh_count": sum(new_region_values.values()),
        "matching_endpoint_length_records": sum((ref_region_values & new_region_values).values()),
        "missing": [list(value) for value in (ref_region_values - new_region_values).elements()],
        "extra": [list(value) for value in (new_region_values - ref_region_values).elements()],
    }

    # Compare each adjacent-threshold Adjusted Wallace point estimate.
    ref_metrics = metrics(reference / "NG_822_metrics.tsv")
    new_metrics = metrics(reproduced / "NG_822_metrics.tsv")
    tolerance = 1e-8
    score_diffs = []
    for pair in sorted(set(ref_metrics) | set(new_metrics)):
        if pair not in ref_metrics or pair not in new_metrics:
            score_diffs.append({"col2": pair[0], "col1": pair[1], "reason": "missing score"})
            continue
        difference = new_metrics[pair] - ref_metrics[pair]
        if not math.isclose(new_metrics[pair], ref_metrics[pair], rel_tol=0.0, abs_tol=tolerance):
            score_diffs.append({
                "col2": pair[0], "col1": pair[1],
                "reference": ref_metrics[pair], "fresh": new_metrics[pair],
                "difference": difference,
            })
    result["checks"]["adjusted_wallace"] = {
        "passed": not score_diffs and len(ref_metrics) == 472 and len(new_metrics) == 472,
        "expected_score_count": 472,
        "expected_count_basis": (
            "partitioning log range 0..473 is end-exclusive: 473 partitions yield 472 adjacent transitions"
        ),
        "reference_count_matches_expected": len(ref_metrics) == 472,
        "reference_score_count": len(ref_metrics),
        "fresh_score_count": len(new_metrics),
        "absolute_tolerance": tolerance,
        "mismatch_count": len(score_diffs),
        "mismatches": score_diffs[:100],
    }

    # Compare partitions by isolate membership sets, never generated cluster names.
    ref_columns, ref_memberships = membership_map(reference / "NG_822_metadata_w_partitions.tsv")
    new_columns, new_memberships = membership_map(reproduced / "NG_822_metadata_w_partitions.tsv")
    ref_names = set(ref_columns)
    new_names = set(new_columns)
    membership_diffs = []
    for name in sorted(ref_names | new_names):
        if name not in ref_memberships or name not in new_memberships:
            membership_diffs.append({"partition": name, "reason": "partition column missing"})
            continue
        ref_groups = Counter(ref_memberships[name].values())
        new_groups = Counter(new_memberships[name].values())
        if ref_groups != new_groups:
            membership_diffs.append({
                "partition": name,
                "reference_cluster_count": len(ref_groups),
                "fresh_cluster_count": len(new_groups),
                "matching_membership_sets": sum((ref_groups & new_groups).values()),
                "missing_membership_sets": sum((ref_groups - new_groups).values()),
                "extra_membership_sets": sum((new_groups - ref_groups).values()),
            })
    result["checks"]["selected_partition_memberships"] = {
        "passed": not membership_diffs,
        "reference_partition_count": len(ref_names),
        "fresh_partition_count": len(new_names),
        "partitions_compared": sorted(ref_names & new_names),
        "membership_mismatch_count": len(membership_diffs),
        "mismatches": membership_diffs,
        "identity_basis": "sample-ID sets per cluster; cluster labels ignored",
    }

    # Compare cluster counts, geographic composition, and date summaries after
    # resolving each group by its actual isolate membership.
    summary = group_summary_diffs(
        reference / "NG_822_partitions_summary.tsv",
        reproduced / "NG_822_partitions_summary.tsv",
    )
    summary["passed"] = all(
        summary[key] == 0 for key in (
            "missing_groups_by_sample_identity", "extra_groups_by_sample_identity",
            "summary_field_mismatch_count",
        )
    )
    result["checks"]["group_country_date_summaries"] = summary

    result["overall"] = "PASS" if all(check["passed"] for check in result["checks"].values()) else "FAIL"
    report_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Verification {result['overall']}: {report_path}")
    for name, check in result["checks"].items():
        print(f"  {name}: {'PASS' if check['passed'] else 'FAIL'}")
    return 0 if result["overall"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
