"""Scheduler diagnosis stays bounded, sanitized and separate from completion."""

import subprocess
from types import SimpleNamespace

import pytest

from chronoclade.artifacts import file_sha256
from chronoclade.execution.accounting import task_accounting


@pytest.mark.parametrize(
    "state,exit_code",
    [("CANCELLED by 123", "0:15"), ("TIMEOUT", "0:15"), ("OUT_OF_MEMORY", "0:125")],
)
def test_bound_task_accounting_categories(tmp_path, monkeypatch, state, exit_code):
    path = tmp_path / "plan.json"
    path.write_text("{}")
    plan = {"tasks": [{"sample_id": "a"}]}
    status = {0: {"sample_id": "a", "plan_sha256": file_sha256(path), "slurm_job_id": "123"}}
    calls = []

    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        return SimpleNamespace(
            returncode=0,
            stdout=f"123|120_0|{state}|{exit_code}|\n123.batch|123.batch|FAILED|1:0|\n",
        )

    monkeypatch.setattr(subprocess, "run", run)
    outcome = task_accounting(path, plan, status, identity_key="sample_id")[0]
    assert outcome["state"] == state.split()[0] and outcome["exit_code"] == exit_code
    assert len(calls) == 1 and calls[0][1]["timeout"] == 5


@pytest.mark.parametrize("mode", ["missing", "timeout", "lag", "unknown", "wrong-plan"])
def test_unavailable_and_untrusted_accounting(tmp_path, monkeypatch, mode):
    path = tmp_path / "plan.json"
    path.write_text("{}")
    plan = {"tasks": [{"sample_id": "a"}]}
    status = {0: {"sample_id": "a", "plan_sha256": file_sha256(path), "slurm_job_id": "123"}}
    if mode == "wrong-plan":
        status[0]["plan_sha256"] = "untrusted"

    def run(*args, **kwargs):
        if mode == "wrong-plan":
            pytest.fail("Unbound receipt must not query scheduler")
        if mode == "missing":
            raise FileNotFoundError()
        if mode == "timeout":
            raise subprocess.TimeoutExpired("sacct", 5)
        return SimpleNamespace(
            returncode=0, stdout="" if mode == "lag" else "123|123|PRIVATE UNTRUSTED TEXT|0:0|"
        )

    monkeypatch.setattr(subprocess, "run", run)
    outcome = task_accounting(path, plan, status, identity_key="sample_id")[0]
    assert outcome["availability"] == "unavailable" and "state" not in outcome


def test_only_plan_bound_submission_and_exact_array_row(tmp_path, monkeypatch):
    path = tmp_path / "plan.json"
    path.write_text("{}")
    plan = {"tasks": [{"sample_id": "a"}]}
    ledger = tmp_path / "submissions.tsv"
    ledger.write_text("array\t999\narray\t888\twrongplan\n")
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: pytest.fail("No bound job"))
    assert "job_id" not in task_accounting(path, plan, {0: {}}, identity_key="sample_id")[0]
    ledger.write_text(f"array\t120\t{file_sha256(path)}\n")
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *a, **k: SimpleNamespace(
            returncode=0, stdout="123|120_0|CANCELLED by 99|0:15|\n124|120_1|OUT_OF_MEMORY|0:9|"
        ),
    )
    outcome = task_accounting(path, plan, {0: {}}, identity_key="sample_id")[0]
    assert outcome["job_id"] == "120_0" and outcome["state"] == "CANCELLED"
