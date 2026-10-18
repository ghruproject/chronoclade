"""Ensemble execution retains exact bindings, failed runs and checked resume."""

import json

import pytest

from chronoclade.errors import WorkflowError
from chronoclade.execution.ensemble import (
    load_tree_ensemble, register_tree_ensemble, run_tree_ensemble, run_time_ensemble,
)
from chronoclade.tree_stage.workflow import run_tree
from test_tree_stage import selection_fixture, native_stubs as _native_stubs


@pytest.fixture
def native_stubs(monkeypatch):
    return _native_stubs.__wrapped__(monkeypatch)


def test_full_ensemble_independent_resume_and_automatic_comparison(tmp_path, native_stubs):
    _, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    output = tmp_path / "trees"
    result = run_tree_ensemble(ensemble, output=output, ensemble_jobs=2)
    assert result["completed_runs"] == 2 and result["failed_runs"] == 0
    assert 'chronoclade time ' in (output / 'README.md').read_text()
    assert len(native_stubs) == 2
    checked, _ = load_tree_ensemble(output / "tree_ensemble.json")
    assert checked == result
    assert run_tree_ensemble(ensemble, output=output, ensemble_jobs=2) == result
    assert len(native_stubs) == 2
    dated = run_time_ensemble(output / "tree_ensemble.json", output=tmp_path / "dates")
    assert dated["completed_runs"] == 2
    assert all(row["dating_status"] == "unsupported" for row in dated["runs"])
    comparison_path = tmp_path / "dates" / dated["comparison"]["path"]
    comparison = json.loads(comparison_path.read_text())
    assert comparison["requested_runs"] == 2 and comparison["dated_runs"] == 0
    assert [r["status"] for r in comparison["runs"]] == ["unsupported", "unsupported"]
    assert run_time_ensemble(output / "tree_ensemble.json", output=tmp_path / "dates") == dated


def test_failed_selection_retained_and_secret_never_published(tmp_path, native_stubs, monkeypatch):
    _, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    from chronoclade.tree_stage import workflow
    original = workflow.run_tree

    def one_failure(path, **kwargs):
        if json.loads(path.read_text())["selection_id"] == "selection-002":
            raise WorkflowError("failed with private-secret")
        return original(path, **kwargs)

    monkeypatch.setattr(workflow, "run_tree", one_failure)
    output = tmp_path / "trees"
    index = run_tree_ensemble(ensemble, output=output, api_key="private-secret")
    assert [r["status"] for r in index["runs"]] == ["complete", "failed"]
    assert "private-secret" not in (output / "tree_ensemble.json").read_text()
    dated = run_time_ensemble(output / "tree_ensemble.json", output=tmp_path / "dates")
    assert [r["status"] for r in dated["runs"]] == ["complete", "failed"]
    comparison = json.loads((tmp_path / "dates" / dated["comparison"]["path"]).read_text())
    assert [r["status"] for r in comparison["runs"]] == ["unsupported", "failed"]
    assert "[redacted]" in comparison["runs"][1]["reason"]


def test_changed_parameters_require_explicit_replacement(tmp_path, native_stubs):
    _, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree_ensemble(ensemble, output=tmp_path / "trees", seed=1)
    with pytest.raises(WorkflowError, match="parameters changed"):
        run_tree_ensemble(ensemble, output=tmp_path / "trees", seed=2)


def test_no_completed_trees_guide_requires_recovery_before_time(tmp_path, monkeypatch):
    _, ensemble = selection_fixture(tmp_path, dates=("", "", ""))

    def fail(*args, **kwargs):
        raise WorkflowError('Missing assemblies; see assembly_audit.json')

    monkeypatch.setattr('chronoclade.tree_stage.workflow.run_tree', fail)
    output = tmp_path / 'trees'
    result = run_tree_ensemble(ensemble, output=output)
    assert result['completed_runs'] == 0
    assert result['failed_runs'] == 2
    guide = (output / 'README.md').read_text()
    assert 'Fix the failures before running time' in guide
    assert 'tree_ensemble.json' in guide
    assert 'chronoclade time ' not in guide


def test_registration_checks_sources_and_retains_missing_tree(tmp_path, native_stubs):
    selection, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "single")
    index = register_tree_ensemble(ensemble, [tmp_path / "single/tree.json"], output=tmp_path / "index")
    assert [r["status"] for r in index["runs"]] == ["complete", "unavailable"]
    load_tree_ensemble(tmp_path / "index/tree_ensemble.json")
    with pytest.raises(WorkflowError, match="Duplicate"):
        register_tree_ensemble(ensemble, [tmp_path / "single/tree.json"] * 2, output=tmp_path / "duplicate")
    # Registration is immutable even when a stage's latest pointer is replaced.
    (tmp_path / "single/tree.json").write_text("{}")
    load_tree_ensemble(tmp_path / "index/tree_ensemble.json")
    frozen = (tmp_path / "index" / index["runs"][0]["result"]["path"]).resolve()
    frozen.write_text("{}")
    with pytest.raises(WorkflowError):
        load_tree_ensemble(tmp_path / "index/tree_ensemble.json")


def test_registration_rejects_changed_trees_and_preserves_prior_index(tmp_path, native_stubs):
    selection, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "single")
    register_tree_ensemble(ensemble, [tmp_path / "single/tree.json"], output=tmp_path / "index")
    previous = (tmp_path / "index/tree_ensemble.json").read_bytes()
    run_tree(selection, output=tmp_path / "single", force=True)
    with pytest.raises(WorkflowError, match="parameters changed"):
        register_tree_ensemble(ensemble, [tmp_path / "single/tree.json"], output=tmp_path / "index")
    assert (tmp_path / "index/tree_ensemble.json").read_bytes() == previous
    load_tree_ensemble(tmp_path / "index/tree_ensemble.json")


def test_invalid_parallelism_does_not_start_jobs(tmp_path):
    with pytest.raises(WorkflowError, match="positive"):
        run_tree_ensemble(tmp_path / "absent.json", output=tmp_path / "trees", ensemble_jobs=0)


def test_checked_time_reuse_avoids_execution_and_validates_gate(tmp_path, native_stubs, monkeypatch):
    _, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree_ensemble(ensemble, output=tmp_path / "trees")
    original = run_time_ensemble(tmp_path / "trees/tree_ensemble.json", output=tmp_path / "dates")
    paths = [tmp_path / "dates" / row["result"]["path"] for row in original["runs"]]
    monkeypatch.setattr("chronoclade.time_stage.workflow.run_time",
                        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("no rerun")))
    result = run_time_ensemble(tmp_path / "trees/tree_ensemble.json", output=tmp_path / "reused",
                               reuse_results=paths)
    assert result["completed_runs"] == 2
    assert all(row["execution"] == "reused_checksum_verified_result" for row in result["runs"])
    with pytest.raises(WorkflowError, match="gate parameters"):
        run_time_ensemble(tmp_path / "trees/tree_ensemble.json", output=tmp_path / "wrong",
                           reuse_results=paths, allow_unsupported=True)
    with pytest.raises(WorkflowError, match="uniquely"):
        run_time_ensemble(tmp_path / "trees/tree_ensemble.json", output=tmp_path / "duplicate",
                           reuse_results=paths + paths[:1])


@pytest.mark.parametrize("change", ["legacy", "different_seed", "resolve_polytomies"])
def test_time_reuse_rejects_different_inference_policy(tmp_path, native_stubs, change):
    from chronoclade.time_stage.workflow import run_time
    selection, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "tree")
    register_tree_ensemble(ensemble, [tmp_path / "tree/tree.json"], output=tmp_path / "trees")
    run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    path = tmp_path / "time/time.json"
    altered = json.loads(path.read_text())
    if change == "legacy":
        altered["parameters"] = {"allow_unsupported": False}
    elif change == "different_seed":
        altered["parameters"]["dating_seed"] += 1
    else:
        altered["parameters"]["keep_polytomies"] = False
    path.write_text(json.dumps(altered))
    with pytest.raises(WorkflowError, match="native dating inference policy or seed"):
        run_time_ensemble(tmp_path / "trees/tree_ensemble.json", output=tmp_path / "reused", reuse_results=[path])
    assert not (tmp_path / "reused/time_ensemble.json").exists()
