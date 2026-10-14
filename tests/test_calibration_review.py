"""Independent leakage and frozen-root predictive-metric checks."""

import copy
from io import StringIO

import numpy as np
import pytest
from Bio import Phylo

import chronoclade.cgmlst.calibration as calibration
from chronoclade.metadata_dates import sample_date_interval


def _rows():
    return [
        {"sample_id": chr(97 + i), "collection_date": f"{2000 + i}-01-01", "role": "input"}
        for i in range(6)
    ]


def _tree():
    return Phylo.read(StringIO("(a:.01,b:.02,c:.03,d:.04,e:.05,f:.06);"), "newick")


def test_heldout_dates_do_not_choose_training_root_or_fit(monkeypatch):
    source, records = _tree(), _rows()
    before = copy.deepcopy(records)
    changes = copy.deepcopy(records)
    changes[-1]["collection_date"] = "2005-12-31"
    optimise = calibration.optimise_clock_root
    calls = []

    def capture(tree, rows):
        calls.append([(row["sample_id"], row["collection_date"]) for row in rows])
        return optimise(tree, rows)

    monkeypatch.setattr(calibration, "optimise_clock_root", capture)
    original = calibration.heldout_date_checks(
        source, records, 100, catalogue_complete=False, folds=1
    )[0]
    changed = calibration.heldout_date_checks(
        source, changes, 100, catalogue_complete=False, folds=1
    )[0]

    assert (
        calls[0] == calls[1] == [(row["sample_id"], row["collection_date"]) for row in records[:-1]]
    )
    assert original["root_audit"] == changed["root_audit"]
    assert original["root_audit"]["dated_tip_ids"] == list("abcde")
    assert original["training_fit"] == changed["training_fit"]
    assert original["training_rate"] == changed["training_rate"]
    assert original["training_sample_ids"] == changed["training_sample_ids"] == list("abcde")
    assert original["test_sample_ids"] == changed["test_sample_ids"] == ["f"]
    old, new = original["predictions"][0], changed["predictions"][0]
    assert (
        old["observed_root_to_tip_allele_distance"] == new["observed_root_to_tip_allele_distance"]
    )
    date_difference = (
        sample_date_interval(changes[-1])["midpoint_year"]
        - sample_date_interval(records[-1])["midpoint_year"]
    )
    assert new["predicted_root_to_tip_allele_distance"] - old[
        "predicted_root_to_tip_allele_distance"
    ] == pytest.approx(original["training_rate"] * date_difference)
    assert records == before


def test_forward_ols_and_constant_baselines_use_same_frozen_root_and_units(monkeypatch):
    source = Phylo.read(StringIO("(a:.011,b:.021,c:.032,d:.038,e:.053,f:.059);"), "newick")
    records = _rows()
    # Freeze the root independently so the checks cannot share a mistake in
    # root-to-tip computation with the optimizer's audit.
    monkeypatch.setattr(calibration, "optimise_clock_root", lambda tree, train: (tree, {}))
    output = calibration.heldout_date_checks(source, records, 628, catalogue_complete=False)
    for fold in output:
        train_ids, test_ids = fold["training_sample_ids"], fold["test_sample_ids"]
        dates = np.array(
            [
                sample_date_interval(row)["midpoint_year"]
                for row in records
                if row["sample_id"] in train_ids
            ]
        )
        distances = np.array([source.distance(name) * 628 for name in train_ids])
        centred = dates - dates.mean()
        design = np.column_stack((np.ones(len(dates)), centred))
        intercept, slope = np.linalg.lstsq(design, distances, rcond=None)[0]
        fitted = intercept + slope * centred
        assert fold["training_rate"] == pytest.approx(slope)
        assert fold["training_fit"]["fraction_slope"] * 628 == pytest.approx(slope)
        assert fold["training_fit"]["residual_mse"] == pytest.approx(
            np.mean((distances - fitted) ** 2)
        )
        assert fold["training_fit"]["residual_mse"] == pytest.approx(
            np.mean((distances / 628 - fitted / 628) ** 2) * 628**2
        )
        assert "observed export loci" in fold["distance_units"]
        test_dates = np.array(
            [
                sample_date_interval(row)["midpoint_year"]
                for row in records
                if row["sample_id"] in test_ids
            ]
        )
        observed = np.array([source.distance(name) * 628 for name in test_ids])
        predicted = intercept + slope * (test_dates - dates.mean())
        latest_year = max(int(value) for value in dates)
        latest_distance = distances[np.floor(dates) == latest_year].mean()
        assert fold["training_mean_root_distance"] == pytest.approx(distances.mean())
        assert fold["training_latest_year_mean_root_distance"] == pytest.approx(latest_distance)
        assert fold["model_distance_mae"] == pytest.approx(np.mean(np.abs(observed - predicted)))
        assert fold["baseline_distance_mae"] == pytest.approx(
            np.mean(np.abs(observed - distances.mean()))
        )
        assert fold["persistence_distance_mae"] == pytest.approx(
            np.mean(np.abs(observed - latest_distance))
        )
        assert fold["model_distance_mae"] == pytest.approx(
            np.mean([p["absolute_distance_error"] for p in fold["predictions"]])
        )


def test_date_interval_overlapping_test_year_is_neither_training_nor_test():
    records = _rows()
    records.append(
        {
            "sample_id": "overlap",
            "date_start": "2004-12-31",
            "date_end": "2005-01-02",
            "date_precision": "interval",
        }
    )
    fold = calibration.rolling_date_folds(records, folds=1)[0]
    assert fold["held_out_year"] == 2005
    assert [row["sample_id"] for row in fold["train"]] == list("abcde")
    assert [row["sample_id"] for row in fold["test"]] == ["f"]
    assert fold["unused_sample_ids"] == ["overlap"]


def test_negative_rate_retains_forward_metrics_and_blocks_secondary_inversion(monkeypatch):
    source = Phylo.read(StringIO("(a:.06,b:.05,c:.04,d:.03,e:.02,f:.01);"), "newick")
    monkeypatch.setattr(calibration, "optimise_clock_root", lambda tree, train: (tree, {}))
    folds = calibration.heldout_date_checks(source, _rows(), 100, catalogue_complete=False)
    for fold in folds:
        assert fold["training_rate"] == pytest.approx(-1)
        assert fold["status"] == "nonpositive_rate"
        assert fold["model_distance_mae"] == pytest.approx(0, abs=1e-12)
        assert fold["baseline_distance_mae"] > 0
        assert fold["model_date_mae_years"] is None
        assert all(point.get("predicted_midpoint_year") is None for point in fold["predictions"])


def test_native_calibration_binds_full_binary_pool_and_reports_singleton(tmp_path):
    import json

    from chronoclade.cgmlst.scale_distance import write_distances
    from chronoclade.datasets import LocusCatalogue, from_profile_records, write_dataset
    from chronoclade.selection_manifest import load_selection_manifest, write_selection_ensemble

    loci = tuple(f"locus{i}" for i in range(8))
    records = []
    for index, name in enumerate(("q", "c1", "c2", "c3", "c4", "c5", "q_single", "outside_pool")):
        singleton = name == "q_single"
        differences = 3 if singleton else index
        records.append(
            {
                "sample_id": name,
                "origin": "local" if name in {"q", "q_single"} else "context",
                "species": "Klebsiella pneumoniae",
                "mlst_st": "258",
                "country": "Greece",
                "biosample": "SAMN0002" if singleton else "",
                "collection_date": f"{2000 + index}-01-01",
                "cgmlst_scheme": "test",
                "cgmlst_scheme_version": "1",
                "cgmlst_loci": list(loci),
                "cgmlst_profile": {
                    locus: "2" if i < differences else "1" for i, locus in enumerate(loci)
                },
                "cgmlst_status": "resolved",
                "cglin_scheme": "scgMLST629_S",
                "cglin_scheme_version": "1",
                "cglin_raw": "0,0,1,0,1,0," + ("98" if singleton else "0") + ",0,0,0",
                "cglin_status": "complete",
            }
        )
    prepared = write_dataset(
        from_profile_records(records, [LocusCatalogue("test", "1", loci, "test")]),
        tmp_path / "prepared",
    )
    pool = records[:-1]
    evidence, _, _ = write_distances(pool, tmp_path / "binary", min_overlap=0.9)
    queries = [row for row in pool if row["origin"] == "local"]
    context = [row for row in pool if row["origin"] == "context"]
    ensemble = write_selection_ensemble(
        queries,
        context,
        {"distance_evidence": evidence},
        dataset_manifest=prepared,
        output=tmp_path / "selection",
        size=1,
        nearest_per_query=1,
    )
    entry = json.loads(ensemble.read_text())["selections"][0]
    selection = load_selection_manifest(ensemble.parent / entry["path"])
    assert len(selection["selected_sample_ids"]) < len(pool)
    output = tmp_path / "calibration"
    result = calibration.calibrate_lin_groups(ensemble, output=output, folds=1)
    assert result["status"] == "complete"
    assert len(result["groups"]) == 3
    full = next(group for group in result["groups"] if group["level"] == 6)
    assert full["n_genomes"] == 7
    assert set(full["selected_sample_ids"]) == {row["sample_id"] for row in pool}
    assert "outside_pool" not in full["selected_sample_ids"]
    assert full["n_inputs"] == 2 and full["n_public"] == 5
    assert full["status"] == "descriptive_fit"
    assert full["normalization_loci"] == 8
    assert full["full_date_fit"]["n_points"] == 7
    assert full["full_date_fit"]["slope"] is not None
    singleton = next(group for group in result["groups"] if group["prefix"][-1] == 98)
    assert singleton["n_genomes"] == singleton["n_dated"] == 1
    assert singleton["status"] == "insufficient_dates"
    assert singleton["full_date_fit"] is None
    assert singleton["conditional_pair_distance_scenarios"] == []
    assert singleton["input_accessions"] == ["SAMN0002"]
    stored = json.loads((output / "calibration.json").read_text())
    assert stored["status"] == "complete"
    report = (output / "report.html").read_text()
    assert "SAMN0002" in report
    assert "8 loci in the frozen export" in report
    assert "No biological cutoff is recommended" in report
    assert "outside_pool" not in report
