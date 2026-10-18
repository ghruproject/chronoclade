"""Ownership and CLI-surface checks for resumable prepare work."""

from pathlib import Path

import pytest

import chronoclade
from chronoclade.datasets import DatasetError
from chronoclade.prepare_recovery import DirectoryLock, RecoveryError, begin_job, owned_directory
from chronoclade.prepare_provider import _resume_checkpoint_ids
from chronoclade.query_typing import QueryTypingError
from typer.testing import CliRunner


def test_imported_package_origin_is_this_staging_copy():
    staging_root = Path(__file__).resolve().parents[1]
    assert Path(chronoclade.__file__).resolve().is_relative_to(staging_root)


def test_job_manifest_is_bound_and_lock_is_process_owned(tmp_path):
    work = tmp_path / "target.work"
    binding = {"target": str(tmp_path / "target"), "source": "a" * 64}
    first = begin_job(work, binding, resume=False)
    assert first["job_id"]
    assert begin_job(work, binding, resume=True)["job_id"] == first["job_id"]
    with pytest.raises(RecoveryError, match="identity-mismatched"):
        begin_job(work, {**binding, "source": "b" * 64}, resume=True)
    with DirectoryLock(tmp_path / "job.lock"):
        with pytest.raises(RecoveryError, match="Another process"):
            with DirectoryLock(tmp_path / "job.lock"):
                pass


def test_lock_refuses_symlink_without_touching_target(tmp_path):
    target = tmp_path / "unowned.txt"
    target.write_text("preserve this evidence\n")
    lock = tmp_path / "job.lock"
    lock.symlink_to(target)
    with pytest.raises(RecoveryError, match="symlink"):
        with DirectoryLock(lock):
            pass
    assert target.read_text() == "preserve this evidence\n"


def test_resume_accepts_absent_cache_but_rejects_incomplete_or_linked_cache(tmp_path):
    cache = tmp_path / "native_typing_cache"
    assert _resume_checkpoint_ids(cache) == []
    cache.mkdir()
    with pytest.raises(RecoveryError, match="Cannot read typing recovery manifest"):
        _resume_checkpoint_ids(cache)
    cache.rmdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    cache.symlink_to(outside, target_is_directory=True)
    with pytest.raises(DatasetError, match="symlink"):
        _resume_checkpoint_ids(cache)


def test_resume_accepts_only_valid_uninitialized_typing_manifest(tmp_path):
    import json

    from chronoclade.prepare_recovery import canonical_digest

    cache = tmp_path / "native_typing_cache"
    cache.mkdir()
    identity = {"schema": "chronoclade.query-typing-recovery-identity", "samples": []}
    manifest = {
        "schema": "chronoclade.query-typing-recovery",
        "schema_version": 1,
        "status": "incomplete",
        "identity": identity,
        "identity_sha256": canonical_digest(identity),
        "checkpoint_initialized": False,
    }
    (cache / "typing-recovery.json").write_text(json.dumps(manifest))
    assert _resume_checkpoint_ids(cache) == []

    manifest["identity_sha256"] = "0" * 64
    (cache / "typing-recovery.json").write_text(json.dumps(manifest))
    with pytest.raises(QueryTypingError, match="Missing typing checkpoint"):
        _resume_checkpoint_ids(cache)


def test_owned_attempt_paths_reject_symlink_children(tmp_path):
    work = tmp_path / "target.work"
    begin_job(work, {"target": str(tmp_path / "target")}, resume=False)
    outside = tmp_path / "unowned"
    outside.mkdir()
    (outside / "marker").write_text("untouched")
    (work / "attempts").symlink_to(outside, target_is_directory=True)
    with pytest.raises(RecoveryError, match="symlink"):
        owned_directory(work, "attempts/next-attempt", create=True)
    assert not (outside / "next-attempt").exists()
    assert (outside / "marker").read_text() == "untouched"


def test_resume_rejects_symlinked_attempt_root_before_adopting_manifest(tmp_path):
    work = tmp_path / "target.work"
    manifest = begin_job(work, {"target": str(tmp_path / "target")}, resume=False)
    manifest["attempts"] = [{"attempt_id": "old", "evidence_directory": "attempts/old"}]
    (work / "prepare-recovery.json").write_text(__import__("json").dumps(manifest))
    outside = tmp_path / "unowned"
    outside.mkdir()
    (work / "attempts").symlink_to(outside, target_is_directory=True)
    with pytest.raises(RecoveryError, match="symlink"):
        begin_job(work, {"target": str(tmp_path / "target")}, resume=True)
    assert list(outside.iterdir()) == []


def test_prepare_and_type_queries_help_expose_explicit_resume():
    from click import unstyle
    from chronoclade.cli import app

    runner = CliRunner()
    prepare = runner.invoke(app, ["prepare", "--help"])
    typing = runner.invoke(app, ["type-queries", "--help"])
    assert prepare.exit_code == 0, prepare.output
    assert typing.exit_code == 0, typing.output
    assert "--resume" in unstyle(prepare.output)
    assert "--resume" in unstyle(typing.output)


def test_resume_rejects_frozen_input_and_requires_native_config(tmp_path):
    from chronoclade.prepare_stage import run_prepare

    path = tmp_path / "frozen.json"
    path.write_text("[]\n")
    with pytest.raises(DatasetError, match="only for live"):
        run_prepare(path, tmp_path / "out", input_kind="profiles", resume=True)
    (tmp_path / "assemblies.csv").write_text(
        "sample_id,assembly,species\nq,missing.fa,Klebsiella pneumoniae\n"
    )
    with pytest.raises(DatasetError, match="requires --typing-config"):
        run_prepare(tmp_path / "assemblies.csv", tmp_path / "out2",
                    input_kind="assemblies", resume=True)
