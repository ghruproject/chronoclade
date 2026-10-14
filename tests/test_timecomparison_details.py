"""Saved ancestor uncertainty, descriptive topology and coverage audit contracts."""

import json
from io import StringIO

import pytest
from Bio import Phylo

from chronoclade.errors import WorkflowError
from chronoclade.report_components.forest import date_forest
from chronoclade.selection_manifest import load_selection_ensemble
from chronoclade.stage_artifacts import reference
from chronoclade.time_stage.comparison import compare_time_runs, _spread
from chronoclade.time_stage.comparison_coverage import selection_coverage
from chronoclade.time_stage.comparison_evidence import (
    anchor_date_evidence,
    date_summary,
    read_node_dates,
    shared_anchor_splits,
    split_occurrence,
)
from chronoclade.time_stage.workflow import run_time
from chronoclade.tree_stage.workflow import run_tree
from test_tree_stage import selection_fixture, native_stubs as _tree_native_stubs


@pytest.fixture
def native_stubs(monkeypatch):
    return _tree_native_stubs.__wrapped__(monkeypatch)


def dated_evidence(tmp_path, *, csv_text=None):
    tree = Phylo.read(StringIO("((a:1,b:1,c:1)ANCESTOR:1,(d:1,e:1)OTHER:1)ROOT;"), "newick")
    dated = tmp_path / "dated.nexus"
    Phylo.write(tree, dated, "nexus")
    dates = tmp_path / "node_dates.csv"
    dates.write_text(
        csv_text
        or "node,numeric_date,lower_bound,upper_bound\nROOT,2000,1999,2001\nANCESTOR,2012,2010,2014\n"
    )
    return {
        "dated_tree": reference(tmp_path, dated),
        "artifacts": {"node_dates.csv": reference(tmp_path, dates)},
    }


def test_matches_internal_anchor_ancestor_and_reports_optional_descendants(tmp_path):
    result = dated_evidence(tmp_path)
    evidence = anchor_date_evidence(tmp_path / "time.json", result, ["a", "b"])
    assert evidence["node"] == "ANCESTOR"
    assert evidence["estimate"] == 2012 and evidence["lower_bound"] == 2010
    assert evidence["upper_bound"] == 2014
    assert evidence["descendant_sample_ids"] == ["a", "b", "c"]
    assert evidence["anchor_count"] == 2 and evidence["descendant_count"] == 3
    assert evidence["interval_level"] == "90% max-posterior region"


def test_missing_ancestor_or_anchor_is_an_evidence_error(tmp_path):
    result = dated_evidence(
        tmp_path, csv_text="node,numeric_date,lower_bound,upper_bound\nROOT,2000,1999,2001\n"
    )
    with pytest.raises(WorkflowError, match="identify the shared-anchor ancestor"):
        anchor_date_evidence(tmp_path / "time.json", result, ["a", "b"])
    with pytest.raises(WorkflowError, match="occur in the dated"):
        anchor_date_evidence(tmp_path / "time.json", result, ["a", "missing"])


def test_unsupported_interval_level_cannot_be_labelled_as_ninety_percent(tmp_path):
    result = dated_evidence(tmp_path)
    result["clock_confidence"] = {"interval_level": "95% confidence interval"}
    with pytest.raises(WorkflowError, match="Unsupported saved dating interval"):
        anchor_date_evidence(tmp_path / "time.json", result, ["a", "b"])


@pytest.mark.parametrize(
    "text",
    [
        "node,numeric_date\nx,2000\n",
        "node,numeric_date,lower_bound,upper_bound\nx,nan,1999,2001\n",
        "node,numeric_date,lower_bound,upper_bound\nx,2000,1999,inf\n",
        "node,numeric_date,lower_bound,upper_bound\nx,2000,2001,2002\n",
        "node,numeric_date,lower_bound,upper_bound\nx,2000,,2001\n",
        "node,numeric_date,lower_bound,upper_bound\nx,2000,1999,2001\nx,2000,1999,2001\n",
        "node,numeric_date,lower_bound,upper_bound\nx,2000,1999\n",
        "node,numeric_date,lower_bound,upper_bound\nx,2000,1999,2001,extra\n",
    ],
)
def test_rejects_malformed_node_date_csv(tmp_path, text):
    path = tmp_path / "dates.csv"
    path.write_text(text)
    with pytest.raises(WorkflowError, match="Node-date CSV"):
        read_node_dates(path)


def test_missing_bounds_stay_unavailable_and_read_failures_are_explained(tmp_path):
    path = tmp_path / "dates.csv"
    path.write_text("node,numeric_date,lower_bound,upper_bound\nx,2000,,\n")
    assert read_node_dates(path)["x"] == {
        "estimate": 2000,
        "lower_bound": None,
        "upper_bound": None,
    }
    with pytest.raises(WorkflowError, match="Cannot read node-date CSV"):
        read_node_dates(tmp_path / "missing.csv")


def test_original_partial_dates_and_saved_eligibility_are_preserved():
    rows = [
        {"sample_id": "year", "collection_date": "2018"},
        {"sample_id": "month", "collection_date": "2020-02"},
        {"sample_id": "interval", "collection_date": "2019/2021"},
        {"sample_id": "missing", "collection_date": ""},
    ]
    evidence = date_summary(rows, ["year", "month", "interval"])
    assert evidence["dated_sample_count"] == 3 and evidence["undated_sample_count"] == 1
    assert evidence["date_interval"] == {"start": "2018-01-01", "end": "2021-12-31"}
    assert evidence["date_precision_counts"] == {"year": 1, "month": 1, "interval": 1}


def test_unrooted_splits_are_reroot_invariant_and_project_optional_tips(tmp_path):
    biological = Phylo.read(StringIO("((a:1,b:1,x:1):1,(c:1,d:1):1);"), "newick")
    path = tmp_path / "corrected.nwk"
    Phylo.write(biological, path, "newick")
    before = shared_anchor_splits(
        tmp_path / "tree.json", {"corrected_tree": reference(tmp_path, path)}, list("abcd")
    )
    biological.root_with_outgroup("a")
    Phylo.write(biological, path, "newick")
    after = shared_anchor_splits(
        tmp_path / "tree.json", {"corrected_tree": reference(tmp_path, path)}, list("abcd")
    )
    assert before == after == {(("a", "b"), ("c", "d"))}
    evidence = split_occurrence({"one": before, "two": after, "unresolved": set()}, 4, list("abcd"))
    assert evidence["assessable_runs"] == 3 and evidence["requested_runs"] == 4
    assert evidence["splits"][0]["frequency"] == 2 / 3
    assert evidence["statistical_branch_support"]["status"] == "unavailable"
    assert "not bootstrap support" in evidence["interpretation"]
    assert split_occurrence({"small": set()}, 1, list("abc"))["assessable_runs"] == 0


def test_forest_separates_regions_from_median_range_and_keeps_absent_run():
    spread = _spread([2012, 2014, None])
    assert spread["median"] == 2013 and spread["range"] == 2
    svg = date_forest(
        [
            {"label": "first", "estimate": 2012, "lower_bound": 2010, "upper_bound": 2015},
            {"label": "second", "estimate": 2014, "lower_bound": None, "upper_bound": None},
            {"label": "failed", "estimate": None},
        ],
        spread,
    )
    assert "90% max-posterior" in svg and "not a confidence interval" in svg
    assert "failed" in svg and "bounds unavailable" in svg and "Median 2013" in svg
    assert date_forest([{"label": "failed"}], _spread([])) is None


def test_comparison_keeps_failures_coverage_and_saved_tree_audit(tmp_path, native_stubs):
    selection, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "tree")
    run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    result = compare_time_runs(
        ensemble,
        [tmp_path / "time/time.json"],
        output=tmp_path / "comparison",
        missing_runs={
            "selection-002": {"status": "failed", "reason": "Assembly acquisition failed"}
        },
    )
    coverage = result["selection_coverage"]
    assert len(coverage["runs"]) == 2 and coverage["pool_sample_count"] == 3
    assert coverage["requested_selection_union_count"] == 3
    assert coverage["available_tree_union_count"] == 3
    assert coverage["runs"][0]["representative_coverage"]["origin"] == "saved_selection_artifact"
    assert result["runs"][0]["dated_sample_count"] == 0
    assert result["runs"][0]["complete_alignment_sites"] == 2000
    assert result["runs"][0]["reference"]["sha256"]
    assert result["runs"][1]["status"] == "failed"
    assert "dated_sample_count" not in result["runs"][1]
    html = (tmp_path / "comparison/report.html").read_text()
    assert "Assembly acquisition failed" in html and "Coverage CSV" in html
    assert "Shared genomes" in html and "Unavailable" in html
    assert "Statistical branch support is unavailable" in html
    assert result["shared_anchor_date_forest"] is None
    assert result["sample_labels"] == {"q": "q", "c": "c", "d": "d"}


def test_comparison_renders_accessions_without_changing_pair_identity(tmp_path, native_stubs):
    from chronoclade.time_stage.comparison_report import write_comparison_report

    selection, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "tree")
    run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    output = tmp_path / "comparison"
    result = compare_time_runs(ensemble, [tmp_path / "time/time.json"], output=output)
    pair = {"target_pair": "PW_one|PW_two", "assessable_runs": 5, "co_clustered_runs": 5}
    result["target_pair_agreement"] = [pair]
    result["sample_labels"] = {
        "PW_one": "ERR8772835 (PW_one)", "PW_two": "ERR8772768"
    }
    write_comparison_report(output, result)
    html = (output / "report.html").read_text()
    assert "ERR8772835 · ERR8772768" in html
    assert "PW_one" not in html and "PW_two" not in html
    assert pair["target_pair"] == "PW_one|PW_two"


def test_corrupt_saved_coverage_csv_is_rejected_before_report(tmp_path):
    selection, ensemble = selection_fixture(tmp_path)
    data = json.loads(selection.read_text())
    csv_path = selection.parent / data["representative_coverage"]["artifact"]["path"]
    csv_path.write_text("corrupted")
    with pytest.raises(WorkflowError, match="sha256 does not match"):
        compare_time_runs(ensemble, [], output=tmp_path / "comparison")
    assert not (tmp_path / "comparison/comparison.json").exists()


def test_older_selection_coverage_reuses_shared_engine_without_mutating_sources(tmp_path):
    _, ensemble = selection_fixture(tmp_path)
    manifest = json.loads(ensemble.read_text())
    for entry in manifest["selections"]:
        path = ensemble.parent / entry["path"]
        data = json.loads(path.read_text())
        data.pop("representative_coverage")
        path.write_text(json.dumps(data))
        entry.update(reference(ensemble.parent, path))
    ensemble.write_text(json.dumps(manifest))
    frozen = {
        path: path.read_bytes()
        for path in (
            ensemble,
            *(ensemble.parent / entry["path"] for entry in manifest["selections"]),
        )
    }
    output = tmp_path / "comparison"
    output.mkdir()
    _, coverage = selection_coverage(ensemble, load_selection_ensemble(ensemble), output)
    assert (
        coverage["runs"][0]["representative_coverage"]["origin"]
        == "computed_from_saved_distance_evidence"
    )
    assert (output / "selection-001-coverage.csv").exists()
    assert all(path.read_bytes() == original for path, original in frozen.items())
