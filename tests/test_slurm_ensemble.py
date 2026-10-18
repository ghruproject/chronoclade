"""Exercise array generation, independent workers, partial reduction and retries."""

import json
from pathlib import Path
import subprocess

import pytest
from typer.testing import CliRunner

from chronoclade.artifacts import file_sha256
from chronoclade.cli import app
from chronoclade.errors import WorkflowError
from chronoclade.execution.ensemble import load_tree_ensemble
from chronoclade.execution.jobs import SlurmResources
from chronoclade.execution.slurm import collect, run_task, write_slurm_ensemble
from test_tree_stage import selection_fixture, native_stubs as _native_stubs


@pytest.fixture(autouse=True)
def isolated_allocation(monkeypatch):
    # These worker tests declare their own synthetic allocation, including when
    # pytest itself is running inside a smaller SLURM test-shard allocation.
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)


@pytest.fixture
def native_stubs(monkeypatch):
    return _native_stubs.__wrapped__(monkeypatch)


def prepare(tmp_path, stage="tree", source=None, **kwargs):
    if source is None:
        _, source = selection_fixture(tmp_path, dates=("", "", ""))
    return write_slurm_ensemble(
        stage, source, output=tmp_path / (stage + " jobs"),
        results=tmp_path / (stage + " results"),
        resources=SlurmResources(cpus=2, account="team", partition="short"),
        max_parallel=2, **kwargs,
    )


def test_array_scripts_and_dependency_submission(tmp_path, monkeypatch):
    plan = prepare(tmp_path)
    for name in ["array.slurm", "collect.slurm", "submit.sh"]:
        subprocess.run(["bash", "-n", str(plan.parent / name)], check=True)
    script = (plan.parent / "array.slurm").read_text()
    assert "#SBATCH --array=0-1%2" in script
    assert "#SBATCH --cpus-per-task=2" in script
    assert "OPENBLAS_NUM_THREADS=1" in script
    assert "export PATH=" in script
    assert "SLURM_ARRAY_TASK_ID" in script
    # Fake scheduler verifies actual quoting and dependency wiring without jobs.
    scheduler = tmp_path / "bin"
    scheduler.mkdir()
    (scheduler / "sbatch").write_text(
        '#!/usr/bin/env python3\nimport json,sys\n'
        f'with open({str(tmp_path / "calls")!r}, "a") as f: f.write(json.dumps(sys.argv[1:])+"\\n")\n'
        'print("123;cluster")\n'
    )
    (scheduler / "sbatch").chmod(0o700)
    import os
    monkeypatch.setenv("PATH", str(scheduler) + os.pathsep + os.environ["PATH"])
    subprocess.run(["bash", str(plan.parent / "submit.sh")], check=True)
    calls = [json.loads(line) for line in (tmp_path / "calls").read_text().splitlines()]
    assert (plan.parent / "submissions.tsv").read_text().splitlines() == [
        f"array\t123\t{file_sha256(plan)}", "collector\t123",
    ]
    assert calls == [["--parsable", str(plan.parent / "array.slurm")],
                     ["--parsable", "--dependency=afterany:123", str(plan.parent / "collect.slurm")]]


def test_workers_collect_tree_and_time_without_rerunning(tmp_path, native_stubs, monkeypatch):
    plan = prepare(tmp_path)
    run_task(plan, 1)
    run_task(plan, 0)
    assert collect(plan) == 0
    trees = tmp_path / "tree results/tree_ensemble.json"
    index, _ = load_tree_ensemble(trees)
    assert [r["selection_id"] for r in index["runs"]] == ["selection-001", "selection-002"]
    time_plan = prepare(tmp_path, "time", trees)
    run_task(time_plan, 0)
    run_task(time_plan, 1)
    monkeypatch.setattr("chronoclade.time_stage.workflow.run_time",
                        lambda *a, **k: pytest.fail("Collector reran dating"))
    assert collect(time_plan) == 0
    dated = json.loads((tmp_path / "time results/time_ensemble.json").read_text())
    assert dated["completed_runs"] == 2
    assert all(r["dating_status"] == "unsupported" for r in dated["runs"])
    assert (tmp_path / "time results" / dated["report"]["path"]).is_file()
    assert len(native_stubs) == 2


def test_missing_tasks_remain_visible_and_collector_does_not_run_them(tmp_path, native_stubs):
    plan = prepare(tmp_path)
    run_task(plan, 0)
    assert collect(plan) == 2
    trees = tmp_path / "tree results/tree_ensemble.json"
    index, _ = load_tree_ensemble(trees)
    assert [r["status"] for r in index["runs"]] == ["complete", "unavailable"]
    time_plan = prepare(tmp_path, "time", trees)
    assert collect(time_plan) == 2
    dated = json.loads((tmp_path / "time results/time_ensemble.json").read_text())
    assert [r["status"] for r in dated["runs"]] == ["failed", "unavailable"]
    assert len(native_stubs) == 1


def test_recollection_after_retry_preserves_previous_index(tmp_path, native_stubs):
    plan = prepare(tmp_path)
    run_task(plan, 0)
    assert collect(plan) == 2
    index = tmp_path / "tree results/tree_ensemble.json"
    original = index.read_bytes()
    run_task(plan, 1)
    assert collect(plan) == 0
    assert any(p.read_bytes() == original for p in index.parent.glob("tree_ensemble-*.json"))
    result, _ = load_tree_ensemble(index)
    assert result["completed_runs"] == 2


def test_bad_inputs_resources_and_changed_source_rejected(tmp_path, monkeypatch):
    _, source = selection_fixture(tmp_path)
    with pytest.raises(ValueError, match="fit cpus"):
        prepare(tmp_path, source=source, options={"threads": 3})
    with pytest.raises(ValueError, match="Unsupported"):
        prepare(tmp_path, source=source, options={"api_key": "secret"})
    plan = prepare(tmp_path, source=source)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    with pytest.raises(WorkflowError, match="allocation"):
        run_task(plan, 0)
    source.write_text(source.read_text() + "\n")
    with pytest.raises(WorkflowError, match="sha256"):
        collect(plan)


def test_failed_retry_removes_previous_success_receipt(tmp_path, native_stubs, monkeypatch):
    plan = prepare(tmp_path)
    run_task(plan, 0)
    assert (plan.parent / "task-0.json").exists()
    def fail(*args, **kwargs):
        raise WorkflowError("worker failed")
    monkeypatch.setattr("chronoclade.tree_stage.workflow.run_tree", fail)
    with pytest.raises(WorkflowError, match="worker failed"):
        run_task(plan, 0)
    assert not (plan.parent / "task-0.json").exists()
    assert collect(plan) == 2


def test_cli_prepares_but_does_not_submit(tmp_path, monkeypatch):
    _, source = selection_fixture(tmp_path)
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: pytest.fail("Unexpected submission"))
    result = CliRunner().invoke(app, ["slurm-ensemble", "tree", str(source),
                                     "--out", str(tmp_path / "jobs"),
                                     "--results", str(tmp_path / "results"), "--cpus", "2"])
    assert result.exit_code == 0, result.output
    assert "Submit:" in result.output
    assert Path(tmp_path / "jobs/array.json").is_file()


def test_failure_status_is_redacted_and_collector_explains_retry(tmp_path, monkeypatch):
    plan = prepare(tmp_path)
    monkeypatch.setattr("chronoclade.pathogenwatch.load_api_key", lambda: "private-test-key")
    def fail(*args, **kwargs):
        started = json.loads((plan.parent / "task-0.status.json").read_text())
        assert started["status"] == "running"
        assert started["started_at"]
        raise WorkflowError("request failed private-test-key")
    monkeypatch.setattr("chronoclade.tree_stage.workflow.run_tree", fail)
    with pytest.raises(WorkflowError, match="request failed ") as raised:
        run_task(plan, 0)
    assert "private-test-key" not in str(raised.value)
    status = json.loads((plan.parent / "task-0.status.json").read_text())
    assert status["status"] == "failed"
    assert status["reason"] == "request failed [redacted]"
    assert status["finished_at"]
    assert collect(plan) == 2
    index = json.loads((tmp_path / "tree results/tree_ensemble.json").read_text())
    assert "request failed [redacted]" in index["runs"][0]["reason"]
    assert "retry array index 0" in index["runs"][0]["reason"]
    assert "No task start" in index["runs"][1]["reason"]


def test_cli_preserves_tree_scientific_options(tmp_path):
    _, source = selection_fixture(tmp_path)
    result = CliRunner().invoke(app, [
        "slurm-ensemble", "tree", str(source), "--out", str(tmp_path / "jobs"),
        "--results", str(tmp_path / "results"), "--coherence-policy", "diagnostic",
        "--temporal-p-value", "0.01", "--no-assess-temporal",
    ])
    assert result.exit_code == 0, result.output
    options = json.loads((tmp_path / "jobs/array.json").read_text())["options"]
    assert options["coherence_policy"] == "diagnostic"
    assert options["temporal_p_value"] == 0.01
    assert options["assess_temporal"] is False


def test_module_entrypoint_enables_progress_and_summary(tmp_path, monkeypatch, capsys):
    from chronoclade.execution import slurm
    plan = prepare(tmp_path)
    monkeypatch.setattr("sys.argv", ["slurm", str(plan), "collect"])
    assert slurm.main() == 2
    captured = capsys.readouterr()
    assert "Collected tree: 0/2 complete" in captured.err
    assert "retry array index 0" in captured.err


def test_array_submission_id_survives_collector_submission_failure(tmp_path, monkeypatch):
    import os
    plan = prepare(tmp_path)
    scheduler = tmp_path / "bin"
    scheduler.mkdir()
    (scheduler / "sbatch").write_text(
        '#!/usr/bin/env bash\n'
        'if [[ "$*" == *--dependency=* ]]; then exit 1; fi\n'
        'echo 321\n'
    )
    (scheduler / "sbatch").chmod(0o700)
    monkeypatch.setenv("PATH", str(scheduler) + os.pathsep + os.environ["PATH"])
    submitted = subprocess.run(["bash", str(plan.parent / "submit.sh")], check=False)
    assert submitted.returncode != 0
    assert (plan.parent / "submissions.tsv").read_text() == f"array\t321\t{file_sha256(plan)}\n"


def test_missing_time_scheduler_reason_persisted_in_index_and_comparison(tmp_path, native_stubs, monkeypatch):
    tree_plan = prepare(tmp_path)
    run_task(tree_plan, 0)
    run_task(tree_plan, 1)
    assert collect(tree_plan) == 0
    plan = prepare(tmp_path, stage='time', source=tmp_path / 'tree results/tree_ensemble.json')
    run_task(plan, 0)
    monkeypatch.setattr('chronoclade.execution.accounting.task_accounting', lambda *a, **k:
        {1: {'availability': 'available', 'job_id': '321_1', 'state': 'OUT_OF_MEMORY', 'exit_code': '0:125'}})
    assert collect(plan) == 2
    index = json.loads((tmp_path / 'time results/time_ensemble.json').read_text())
    row = index['runs'][1]
    assert row['scheduler']['state'] == 'OUT_OF_MEMORY'
    assert 'Scheduler OUT_OF_MEMORY' in row['reason']
    comparison = tmp_path / 'time results' / index['comparison']['path']
    assert 'Scheduler OUT_OF_MEMORY' in comparison.read_text()



def test_missing_time_diagnostics_cannot_override_completeness(tmp_path, native_stubs):
    from chronoclade.execution.ensemble import run_time_ensemble
    tree_plan = prepare(tmp_path)
    run_task(tree_plan, 0)
    run_task(tree_plan, 1)
    assert collect(tree_plan) == 0
    with pytest.raises(WorkflowError, match="only attach reason and scheduler"):
        run_time_ensemble(tmp_path / 'tree results/tree_ensemble.json',
            output=tmp_path / 'malicious time', collect_only=True,
            missing_outcomes={'selection-001': {'status': 'complete', 'result': {'path': 'fake'}}})
    assert not (tmp_path / 'malicious time/time_ensemble.json').exists()
