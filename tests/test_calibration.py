"""Blocked dates, scaling, and insufficient-data calibration contracts."""

from io import StringIO
import numpy as np
import pytest
from Bio import Phylo
from chronoclade.cgmlst.calibration import rolling_date_folds, heldout_date_checks


def rows():
    return [
        {"sample_id": chr(97 + i), "collection_date": f"{2000 + i}-01-01", "role": "input"}
        for i in range(6)
    ]


def test_latest_year_blocks_train_strictly_before_test():
    folds = rolling_date_folds(rows())
    assert [f["held_out_year"] for f in folds] == [2003, 2004, 2005]
    for fold in folds:
        assert all(int(r["collection_date"][:4]) < fold["held_out_year"] for r in fold["train"])
        assert all(int(r["collection_date"][:4]) == fold["held_out_year"] for r in fold["test"])
        assert not {r["sample_id"] for r in fold["train"]} & {r["sample_id"] for r in fold["test"]}


def test_training_only_root_and_scaled_prediction(monkeypatch):
    tree = Phylo.read(StringIO("(a:.01,b:.02,c:.03,d:.04,e:.05,f:.06);"), "newick")
    calls = []

    def frozen(current, train):
        calls.append([r["sample_id"] for r in train])
        return current, {"dated_tip_ids": calls[-1]}

    monkeypatch.setattr("chronoclade.cgmlst.calibration.optimise_clock_root", frozen)
    checks = heldout_date_checks(tree, rows(), 100, catalogue_complete=False)
    assert calls == [["a", "b", "c"], ["a", "b", "c", "d"], ["a", "b", "c", "d", "e"]]
    for check in checks:
        assert check["training_rate"] == pytest.approx(1)
        assert check["model_date_mae_years"] == pytest.approx(0, abs=1e-8)
        assert check["model_distance_mae"] == pytest.approx(0, abs=1e-8)
        assert check["baseline_distance_mae"] > 0
        assert check["persistence_distance_mae"] == pytest.approx(1)
        assert check["baseline_date_mae_years"] > 0
        assert check["persistence_date_mae_years"] == pytest.approx(1)
        fit = check["training_fit"]
        assert fit["residual_mse"] == pytest.approx(
            np.mean([p["residual"] ** 2 for p in fit["points"]])
        )
        assert "observed export loci" in fit["slope_units"]


def test_insufficient_training_never_reuses_full_date_root(monkeypatch):
    tree = Phylo.read(StringIO("(a:.1,b:.2,c:.3);"), "newick")

    def forbidden(*args):
        raise AssertionError("Root optimizer must not run for insufficient training")

    monkeypatch.setattr("chronoclade.cgmlst.calibration.optimise_clock_root", forbidden)
    checks = heldout_date_checks(tree, rows()[:3], 100, catalogue_complete=False)
    assert all(c["status"] == "insufficient" for c in checks)
    assert all(c["model_date_mae_years"] is None for c in checks)


def test_negative_rate_blocks_date_inversion(monkeypatch):
    tree = Phylo.read(StringIO("(a:.06,b:.05,c:.04,d:.03,e:.02,f:.01);"), "newick")
    monkeypatch.setattr(
        "chronoclade.cgmlst.calibration.optimise_clock_root", lambda tree, train: (tree, {})
    )
    checks = heldout_date_checks(tree, rows(), 100, catalogue_complete=False)
    assert all(c["status"] == "nonpositive_rate" for c in checks)
    assert all(c["model_date_mae_years"] is None for c in checks)
