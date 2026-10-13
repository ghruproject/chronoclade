"""Compare immutable independent selections on explicitly shared anchors."""

import json
from itertools import combinations
from pathlib import Path
from statistics import median


from chronoclade.artifacts import write_json
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import load_selection_ensemble, _resolve, _read
from chronoclade.stage_artifacts import (
    load_result,
    load_tree_result,
    reference,
    header,
)
from chronoclade.report_components.forest import date_forest
from chronoclade.time_stage.comparison_coverage import selection_coverage
from chronoclade.time_stage.comparison_evidence import (
    anchor_date_evidence,
    tree_audit,
    shared_anchor_splits,
    split_occurrence,
)
from chronoclade.time_stage.comparison_report import write_comparison_report


def _spread(values):
    values = [float(value) for value in values if value is not None]
    return {
        "assessable_runs": len(values),
        "minimum": min(values) if values else None,
        "median": median(values) if values else None,
        "maximum": max(values) if values else None,
        "range": max(values) - min(values) if values else None,
    }


def _clusters(tree_path, tree, ids, cutoff):
    spec = tree["artifacts"].get("snp_distances.json")
    if not spec:
        raise WorkflowError("Tree has no corrected SNP-distance evidence")
    distances = json.loads(_resolve(tree_path.parent, spec, "snp_distances").read_text())
    parent = {ident: ident for ident in tree["selected_sample_ids"]}

    def find(ident):
        while parent[ident] != ident:
            ident = parent[ident]
        return ident

    for row in distances:
        if row["callable_sites"] > 0 and row["clonal_snps"] <= cutoff:
            parent[find(row["sample_2"])] = find(row["sample_1"])
    callability = {
        frozenset((row["sample_1"], row["sample_2"])): row["callable_sites"] for row in distances
    }
    return {
        "|".join((left, right)): (find(left) == find(right))
        if callability.get(frozenset((left, right)), 0) > 0
        else None
        for left, right in combinations(ids, 2)
    }


def compare_time_runs(
    ensemble_manifest, run_manifests, *, output, cluster_snp_cutoff=10, missing_runs=None
):
    if cluster_snp_cutoff < 0:
        raise WorkflowError("Cluster SNP cutoff must be nonnegative")
    ensemble_path = Path(ensemble_manifest).resolve()
    ensemble = load_selection_ensemble(ensemble_path)
    output = Path(output).resolve()
    if (output / "comparison.json").exists():
        raise WorkflowError("Comparison output already exists; use an independent directory")
    output.mkdir(parents=True, exist_ok=True)
    available = {}
    for raw_path in run_manifests:
        path = Path(raw_path).resolve()
        raw = _read(path)
        if raw.get("stage") == "time" and raw.get("status") in {"failed", "pending", "running"}:
            result = raw
            if not isinstance(result.get("selection_id"), str) or not isinstance(
                result.get("source"), dict
            ):
                raise WorkflowError(
                    "Failed/pending time job requires source tree and selection identity"
                )
        else:
            result = load_result(path, "time")
        ident = result["selection_id"]
        if ident in available:
            raise WorkflowError("Duplicate time runs for one selection")
        available[ident] = (path, result)
    expected = {entry["selection_id"] for entry in ensemble["selections"]}
    missing_runs = missing_runs or {}
    if set(missing_runs) - expected or set(missing_runs) & set(available):
        raise WorkflowError("Missing-run outcomes must describe only absent declared selections")
    if any(
        not isinstance(row, dict)
        or row.get("status") not in {"pending", "running", "failed", "unavailable"}
        for row in missing_runs.values()
    ):
        raise WorkflowError("Invalid missing-run outcome")
    if set(available) - expected:
        raise WorkflowError("Time run does not belong to declared selection ensemble")
    _, coverage = selection_coverage(ensemble_path, ensemble, output)
    coverage_by_id = {row["selection_id"]: row for row in coverage["runs"]}
    runs, splits_by_run, tree_union = [], {}, set()
    anchors, targets = ensemble["common_anchor_ids"], ensemble["target_ids"]
    for selection_entry in ensemble["selections"]:
        ident = selection_entry["selection_id"]
        planned = coverage_by_id[ident]
        base = {
            key: planned[key]
            for key in (
                "selection_id",
                "selected_sample_count",
                "selected_input_count",
                "selected_context_count",
                "selection_manifest",
            )
        }
        if ident not in available:
            runs.append(
                {
                    **base,
                    "status": missing_runs.get(ident, {}).get("status", "unavailable"),
                    "reason": missing_runs.get(ident, {}).get(
                        "reason", "No complete time manifest supplied; run may be pending or failed"
                    ),
                }
            )
            continue
        path, result = available[ident]
        tree_path = _resolve(path.parent, result["source"]["tree"], "source.tree", external=True)
        tree = load_tree_result(tree_path)
        if tree["source"]["selection"]["sha256"] != selection_entry["sha256"]:
            raise WorkflowError("Time tree references a different selection checksum")
        base.update(tree_audit(tree_path, tree, output))
        base["tree"] = reference(output, tree_path)
        base["target_clusters"] = _clusters(tree_path, tree, targets, cluster_snp_cutoff)
        splits_by_run[ident] = shared_anchor_splits(tree_path, tree, anchors)
        tree_union.update(tree["selected_sample_ids"])
        if result.get("stage") == "time":
            runs.append(
                {
                    **base,
                    "status": result["status"],
                    "time_job": reference(output, path),
                    "reason": result.get("reason") or "Time job has not completed",
                }
            )
            continue
        confidence = result.get("clock_confidence") or {}
        ancestor = anchor_date_evidence(path, result, anchors, output=output)
        runs.append(
            {
                **base,
                "status": result["dating_status"],
                "time": reference(output, path),
                "reason": result["temporal_assessment"]["reason"],
                "rate": confidence.get("rate"),
                "root_date": confidence.get("root_numeric_date"),
                "common_anchor_mrca_date": ancestor["estimate"],
                "common_anchor_mrca": ancestor,
                "dating_override": result["dating_override"],
            }
        )
    assessable = [run for run in runs if run.get("status") == "dated"]
    pair_labels = sorted({pair for run in runs for pair in run.get("target_clusters", {})})
    agreements = [
        {
            "target_pair": pair,
            "assessable_runs": sum(
                run.get("target_clusters", {}).get(pair) is not None for run in runs
            ),
            "co_clustered_runs": sum(
                run.get("target_clusters", {}).get(pair) is True for run in runs
            ),
        }
        for pair in pair_labels
    ]
    coverage.update(
        available_tree_union_ids=sorted(tree_union),
        available_tree_union_count=len(tree_union),
        available_tree_union_fraction=len(tree_union) / coverage["pool_sample_count"],
    )
    result = {
        **header("time-comparison"),
        "source": {"ensemble": reference(output, ensemble_path)},
        "common_anchor_ids": anchors,
        "target_ids": targets,
        "requested_runs": len(runs),
        "dated_runs": len(assessable),
        "runs": runs,
        "selection_coverage": coverage,
        "shared_anchor_splits": split_occurrence(splits_by_run, len(runs), anchors),
        "sensitivity": {
            field: _spread([run.get(field) for run in assessable])
            for field in ("rate", "root_date", "common_anchor_mrca_date")
        },
        "cluster_policy": {
            "method": "single-linkage corrected SNP distances",
            "snp_cutoff": cluster_snp_cutoff,
            "calibration": "user-declared descriptive sensitivity threshold; no universal biological cutoff",
        },
        "target_pair_agreement": agreements,
        "interpretation": "Between-selection spread measures selection sensitivity, not confidence intervals. Shared-anchor MRCA dates are comparable; whole-selection roots may represent different ancestors. Cluster comparison uses exact target pairs, never run-local cluster numbers.",
    }
    forest = date_forest(
        [{"label": run["selection_id"], **(run.get("common_anchor_mrca") or {})} for run in runs],
        result["sensitivity"]["common_anchor_mrca_date"],
    )
    if forest is not None:
        forest_path = output / "shared-anchor-dates.svg"
        forest_path.write_text(forest, encoding="utf-8")
        result["shared_anchor_date_forest"] = reference(output, forest_path)
    else:
        result["shared_anchor_date_forest"] = None
    write_json(output / "comparison.json", result)
    write_comparison_report(output, result)
    return result
