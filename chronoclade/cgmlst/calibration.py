"""Nested LIN-group descriptive rates and transductive held-out-date checks."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Annotated

import numpy as np
import typer
from Bio import Phylo

from chronoclade.artifacts import write_json
from chronoclade.clock_rooting import optimise_clock_root
from chronoclade.context_partitions import partition_records, profile_records
from chronoclade.datasets import load_dataset
from chronoclade.errors import WorkflowError
from chronoclade.matrix_distances import BinaryDistanceEvidence
from chronoclade.metadata_dates import sample_date_interval
from chronoclade.report_components.root_comparison import _nonnegative_diagnostic_copy
from chronoclade.sample_labels import sample_label
from chronoclade.selection_manifest import _read, _resolve, load_selection_ensemble
from chronoclade.stage_artifacts import header, reference
from .calibration_report import NOTE, _plots, _report
from .clock_comparison import _fit, exact_panel_nj, matched_dates
from .scale_tree import IndexedTree, rapidnj_tree

UNITS = "estimated allele-distance units/year"


def rolling_date_folds(rows, *, folds=3):
    """Hold out complete latest calendar-year blocks with strictly earlier training."""
    intervals = {row["sample_id"]: sample_date_interval(row) for row in rows}
    years = sorted(
        {int(value["year_min"]) for value in intervals.values() if value["reason"] is None}
    )
    results = []
    for year in years[-folds:]:
        train = [
            row
            for row in rows
            if intervals[row["sample_id"]]["reason"] is None
            and intervals[row["sample_id"]]["year_max"] < year
        ]
        test = [
            row
            for row in rows
            if intervals[row["sample_id"]]["reason"] is None
            and int(intervals[row["sample_id"]]["year_min"]) == year
            and intervals[row["sample_id"]]["year_max"] < year + 1
        ]
        used = {row["sample_id"] for row in train + test}
        excluded = [row["sample_id"] for row in rows if row["sample_id"] not in used]
        results.append(
            {"held_out_year": year, "train": train, "test": test, "unused_sample_ids": excluded}
        )
    return results


def heldout_date_checks(tree, rows, loci, *, catalogue_complete, directory=None, folds=3):
    """Train root ONLY on earlier dates; forward distance MAE uses same-root baselines."""
    results = []
    for block in rolling_date_folds(rows, folds=folds):
        train, test = block["train"], block["test"]
        item = {
            "held_out_year": block["held_out_year"],
            "training_sample_ids": [r["sample_id"] for r in train],
            "test_sample_ids": [r["sample_id"] for r in test],
            "unused_sample_ids": block["unused_sample_ids"],
            "n_training": len(train),
            "n_test": len(test),
            "status": "insufficient",
            "reason": "fewer_than_three_training_dates_or_empty_test",
            "training_rate": None,
            "training_r_squared": None,
            "model_date_mae_years": None,
            "baseline_date_mae_years": None,
            "persistence_date_mae_years": None,
            "training_mean_midpoint_year": None,
            "predictions": [],
            "model_distance_mae": None,
            "baseline_distance_mae": None,
            "persistence_distance_mae": None,
            "validation_scope": "transductive held-out collection dates; all group genetics retained",
        }
        train_years = [sample_date_interval(row)["midpoint_year"] for row in train]
        if train_years and test:
            mean_date = float(np.mean(train_years))
            item["training_mean_midpoint_year"] = mean_date
            item["training_latest_midpoint_year"] = max(train_years)
            item["persistence_date_mae_years"] = float(
                np.mean(
                    [
                        abs(sample_date_interval(row)["midpoint_year"] - max(train_years))
                        for row in test
                    ]
                )
            )
            item["baseline_date_mae_years"] = float(
                np.mean(
                    [abs(sample_date_interval(row)["midpoint_year"] - mean_date) for row in test]
                )
            )
        if len(set(train_years)) < 3 or not test:
            results.append(item)
            continue
        rooted, audit = optimise_clock_root(tree, train)
        item["root_audit"] = audit
        if rooted is None:
            item["reason"] = audit["reason"]
            results.append(item)
            continue
        rooted = IndexedTree.from_tree(rooted)
        fit = _fit(
            rooted, train, "cgMLST mismatch fraction", loci, catalogue_complete=catalogue_complete
        )
        item.update(
            training_fit=fit,
            training_rate=fit["slope"],
            training_r_squared=fit["r_squared"],
            rate_units=UNITS,
            normalization_loci=loci,
        )
        if directory is not None:
            Phylo.write(
                rooted,
                Path(directory) / f"fold-{block['held_out_year']}-root.nwk",
                "newick",
                format_branch_length="%.17g",
            )
        rate = fit["slope"]
        if rate is None or not np.isfinite(rate):
            item.update(status="unavailable_fit", reason=fit["reason"])
            results.append(item)
            continue
        train_distances = [p["root_to_tip"] for p in fit["points"]]
        mean_distance = float(np.mean(train_distances))
        latest_training_year = max(int(sample_date_interval(r)["year_min"]) for r in train)
        latest_points = [
            p["root_to_tip"] for p in fit["points"] if int(p["year_min"]) == latest_training_year
        ]
        persistence_distance = float(np.mean(latest_points))
        item.update(
            training_mean_root_distance=mean_distance,
            training_latest_year_mean_root_distance=persistence_distance,
            distance_units=fit["distance_units"],
            primary_prediction="heldout date to same-root allele distance",
        )
        predictions = []
        for row in test:
            interval = sample_date_interval(row)
            distance = rooted.distance(row["sample_id"]) * loci
            predicted_distance = fit["centered_intercept"] + rate * (
                interval["midpoint_year"] - fit["reference_year"]
            )
            point = {
                "sample_id": row["sample_id"],
                "label": sample_label(row),
                "observed_midpoint_year": interval["midpoint_year"],
                "year_min": interval["year_min"],
                "year_max": interval["year_max"],
                "observed_root_to_tip_allele_distance": distance,
                "predicted_root_to_tip_allele_distance": float(predicted_distance),
                "absolute_distance_error": abs(predicted_distance - distance),
                "baseline_absolute_distance_error": abs(mean_distance - distance),
                "persistence_absolute_distance_error": abs(persistence_distance - distance),
            }
            if rate > 0:
                predicted_date = (
                    fit["reference_year"] + (distance - fit["centered_intercept"]) / rate
                )
                point.update(
                    predicted_midpoint_year=float(predicted_date),
                    absolute_error_years=abs(predicted_date - interval["midpoint_year"]),
                    baseline_absolute_error_years=abs(
                        item["training_mean_midpoint_year"] - interval["midpoint_year"]
                    ),
                )
            predictions.append(point)
        item.update(
            status="evaluated" if rate > 0 else "nonpositive_rate",
            reason=None
            if rate > 0
            else "Nonpositive training slope: no positive accumulation rate or inverse date diagnostic",
            predictions=predictions,
            model_distance_mae=float(np.mean([p["absolute_distance_error"] for p in predictions])),
            baseline_distance_mae=float(
                np.mean([p["baseline_absolute_distance_error"] for p in predictions])
            ),
            persistence_distance_mae=float(
                np.mean([p["persistence_absolute_distance_error"] for p in predictions])
            ),
        )
        item["distance_mae_improvement_over_baseline"] = (
            item["baseline_distance_mae"] - item["model_distance_mae"]
        )
        if rate > 0:
            item["model_date_mae_years"] = float(
                np.mean([p["absolute_error_years"] for p in predictions])
            )
            item["mae_improvement_over_baseline_years"] = (
                item["baseline_date_mae_years"] - item["model_date_mae_years"]
            )
        results.append(item)
    return results


def _group_tree(evidence, ids, directory):
    if len(ids) <= 1000:
        tree, loci = exact_panel_nj(evidence, ids)
        tree.root_with_outgroup(ids[0])
        return tree, loci, {"backend": "Biopython NJ", "exact_selected_submatrix": True}
    positions = {name: i for i, name in enumerate(evidence.sample_ids)}
    if not set(ids).issubset(positions):
        raise WorkflowError("Full LIN group has genomes absent from frozen distance evidence")
    index = [positions[name] for name in ids]
    matrix = np.asarray(evidence.arrays["distance"][np.ix_(index, index)])
    if not np.all(np.isfinite(matrix)):
        raise WorkflowError(
            "Full LIN group has unavailable distance pairs; no genomes are silently dropped"
        )
    loci = np.unique(evidence.arrays["scheme_loci"][np.ix_(index, index)])
    if len(loci) != 1 or loci[0] < 1:
        raise WorkflowError("Full LIN group needs one recorded locus universe")
    tree, backend = rapidnj_tree(matrix, ids, directory / "original.nwk")
    return tree, int(loci[0]), backend


def calibrate_lin_groups(selection_ensemble, *, output, folds=3):
    source = Path(selection_ensemble).expanduser().resolve()
    output = Path(output).expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise WorkflowError("Use an empty calibration output to preserve saved evidence")
    ensemble = load_selection_ensemble(source)
    dataset_path = _resolve(source.parent, ensemble["source"]["dataset"], "dataset", external=True)
    evidence_path = _resolve(
        source.parent, ensemble["source"]["distance_evidence"], "distances", external=True
    )
    dataset = load_dataset(dataset_path)
    original = {r["sample_id"]: dict(r) for r in dataset.samples}
    evidence = BinaryDistanceEvidence(evidence_path)
    first_selection = _read(_resolve(source.parent, ensemble["selections"][0], "selection"))
    declared_pool = set(first_selection["query_ids"] + first_selection["available_context_ids"])
    if declared_pool != set(evidence.sample_ids):
        raise WorkflowError("Distance evidence must bind the exact declared selection pool")
    typed = [row for row in profile_records(dataset) if row["sample_id"] in declared_pool]
    original = {name: row for name, row in original.items() if name in declared_pool}
    output.mkdir(parents=True, exist_ok=True)
    result = {
        **header("lin-temporal-calibration"),
        "status": "running",
        "note": NOTE,
        "source": {
            "ensemble": reference(output, source),
            "dataset": reference(output, dataset_path),
            "distance_evidence": reference(output, evidence_path),
        },
        "parameters": {
            "lin_levels": [6, 7],
            "latest_year_folds": folds,
            "minimum_collection_year": 1900,
            "rate_units": UNITS,
            "confidence_intervals": "not estimated: tips are nonindependent",
            "date_randomisations": 0,
            "scenario_policy": "2*r*H assumes shared ancestor H years ago; both descendant branches span H years since their shared ancestor, not collection-date separation; not an empirical threshold",
        },
        "groups": [],
        "partition_audits": {},
    }
    for level in (6, 7):
        blocks, audit = partition_records(typed, lin_level=level)
        result["partition_audits"][str(level)] = audit
        for block in blocks:
            ids = [r["sample_id"] for r in block["records"]]
            rows = [original[name] for name in ids]
            eligible, date_audit = matched_dates(rows, ids)
            directory = output / f"LIN{level}-{block['block_id']}"
            directory.mkdir()
            group = {
                "group_id": directory.name,
                "level": level,
                "prefix": block["prefix"],
                "selected_sample_ids": ids,
                "n_genomes": len(ids),
                "input_accessions": [sample_label(r) for r in rows if r["role"] == "input"],
                "n_inputs": sum(r["role"] == "input" for r in rows),
                "n_public": sum(r["role"] != "input" for r in rows),
                "n_dated": len(eligible),
                "n_distinct_dates": date_audit["n_distinct_dates"],
                "date_filter": date_audit["date_filter"],
                "date_exclusions": date_audit["exclusions"],
                "status": "insufficient_dates",
                "reason": "fewer than three distinct eligible date midpoints",
                "full_date_fit": None,
                "heldout_dates": [],
                "conditional_pair_distance_scenarios": [],
            }
            if len(ids) >= 3 and date_audit["n_distinct_dates"] >= 3:
                sys.setrecursionlimit(max(sys.getrecursionlimit(), 4 * len(ids) + 100))
                tree, loci, backend = _group_tree(evidence, ids, directory)
                catalogs = [
                    m.catalogue
                    for m in dataset.profiles
                    if set(ids).issubset(m.sample_ids) and len(m.catalogue.loci) == loci
                ]
                if len(catalogs) != 1:
                    raise WorkflowError("LIN group must bind one prepared locus catalogue")
                diagnostic, policy = _nonnegative_diagnostic_copy(tree)
                diagnostic = IndexedTree.from_tree(diagnostic)
                rooted, root_audit = optimise_clock_root(diagnostic, eligible)
                group.update(
                    tree_backend=backend,
                    normalization_loci=loci,
                    catalogue_complete=catalogs[0].complete,
                    branch_policy=policy,
                    full_root_audit=root_audit,
                )
                Phylo.write(
                    tree, directory / "original.nwk", "newick", format_branch_length="%.17g"
                )
                Phylo.write(
                    diagnostic,
                    directory / "nonnegative.nwk",
                    "newick",
                    format_branch_length="%.17g",
                )
                if rooted is None:
                    group.update(status="root_unavailable", reason=root_audit["reason"])
                if rooted is not None:
                    fit = _fit(
                        rooted,
                        eligible,
                        "cgMLST mismatch fraction",
                        loci,
                        catalogue_complete=catalogs[0].complete,
                    )
                    checks = heldout_date_checks(
                        diagnostic,
                        eligible,
                        loci,
                        catalogue_complete=catalogs[0].complete,
                        directory=directory,
                        folds=folds,
                    )
                    rates = [f["training_rate"] for f in checks if f["training_rate"] is not None]
                    group.update(
                        status="descriptive_fit",
                        reason=None,
                        full_date_fit=fit,
                        heldout_dates=checks,
                        dated_midpoint_year_range=[
                            min(p["midpoint_year"] for p in fit["points"]),
                            max(p["midpoint_year"] for p in fit["points"]),
                        ],
                        partial_date_count=sum(p["precision"] != "day" for p in fit["points"]),
                        training_fold_rate_range=[min(rates), max(rates)] if rates else None,
                        training_fold_distinct_root_edge_count=len(
                            {
                                f["root_audit"].get("optimal_edge_id")
                                for f in checks
                                if f.get("root_audit")
                            }
                        ),
                        uncertainty_note="Collection intervals use midpoints. Training-fold rate range is sensitivity, not a confidence interval; no independent-tip CI is computed.",
                    )
                    Phylo.write(
                        rooted,
                        directory / "full-date-clock-root.nwk",
                        "newick",
                        format_branch_length="%.17g",
                    )
                    if fit["slope"] is not None and fit["slope"] > 0:
                        group["conditional_pair_distance_scenarios"] = [
                            {
                                "shared_ancestor_horizon_years": h,
                                "expected_pair_distance_units": 2 * fit["slope"] * h,
                                "formula": "2*r*H",
                                "interpretation": "Conditional scenario, not a validated cutoff, allele-event count, or confidence interval",
                            }
                            for h in (5, 10)
                        ]
                    _plots(group, directory / "calibration.svg")
            result["groups"].append(group)
            write_json(directory / "calibration.json", group)
            write_json(output / "calibration.json", result)
    result["status"] = "complete"
    write_json(output / "calibration.json", result)
    _report(result, output)
    return result


def refresh_calibration_report(result_path):
    """Refresh accession display labels from bound metadata without refitting roots."""
    path = Path(result_path).resolve()
    result = _read(path)
    if (
        result.get("schema") != "chronoclade.lin-temporal-calibration"
        or result.get("status") != "complete"
    ):
        raise WorkflowError("Expected a completed LIN temporal calibration result")
    for key, spec in result["source"].items():
        _resolve(path.parent, spec, "source." + key, external=True)
    dataset = load_dataset(
        _resolve(path.parent, result["source"]["dataset"], "dataset", external=True)
    )
    rows = {row["sample_id"]: dict(row) for row in dataset.samples}
    for group in result["groups"]:
        ids = group["selected_sample_ids"]
        if len(ids) != len(set(ids)) or not set(ids).issubset(rows):
            raise WorkflowError("Calibration identifiers must match bound prepared metadata")
        group["input_accessions"] = [
            sample_label(rows[name]) for name in ids if rows[name]["role"] == "input"
        ]
        fits = [group.get("full_date_fit")] + [
            fold.get("training_fit") for fold in group["heldout_dates"]
        ]
        for fit in fits:
            for point in (fit or {}).get("points", []):
                point["label"] = sample_label(rows[point["sample_id"]])
        for fold in group["heldout_dates"]:
            for point in fold["predictions"]:
                point["label"] = sample_label(rows[point["sample_id"]])
        write_json(path.parent / group["group_id"] / "calibration.json", group)
    write_json(path, result)
    _report(result, path.parent)
    return result


def calibrate_lin_command(
    selection_ensemble: Annotated[Path, typer.Argument(help="Frozen full-pool selection ensemble")],
    output: Annotated[Path, typer.Option("--output", "-o")] = Path("lin_temporal_calibration"),
    folds: Annotated[int, typer.Option("--folds", min=1, max=10)] = 3,
):
    """Trial nested LIN6/LIN7 rates with training-only roots and held-out date blocks."""
    try:
        calibrate_lin_groups(selection_ensemble, output=output, folds=folds)
    except (WorkflowError, ValueError, OSError, RuntimeError, KeyError) as error:
        from chronoclade.cli import _fail

        _fail(error)
    typer.echo(f"LIN calibration trial: {output.resolve() / 'report.html'}")
