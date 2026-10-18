"""Synthetic checkpoint tests for the native query-typing recovery boundary."""

import hashlib

import pytest

from chronoclade.metadata import Sample
from chronoclade.query_typing import QueryTypingError, type_query_assemblies


def _sample(tmp_path, name):
    assembly = tmp_path / f"{name}.fa"
    assembly.write_text(f">{name}\nACGT\n")
    return Sample(name, assembly, "2020", "UK", "Klebsiella pneumoniae", "ST1", "local")


def test_completed_sample_checkpoint_survives_interruption_and_resume(tmp_path, monkeypatch):
    samples = [_sample(tmp_path, name) for name in ("first", "second")]
    config = tmp_path / "config.json"
    config.write_text("{}\n")
    provenance = {"execution": {"identity_sha256": "fixture"}}

    def identity(current, _config, profiles_only):
        assert not profiles_only
        return {
            "schema": "chronoclade.query-typing-recovery-identity",
            "schema_version": 1,
            "mode": "full",
            "config_sha256": "a" * 64,
            "execution": provenance["execution"],
            "samples": [
                {
                    "sample_id": sample.sample_id,
                    "species": "klebsiella pneumoniae",
                    "assembly_sha256": hashlib.sha256(sample.assembly.read_bytes()).hexdigest(),
                    "typing_provenance": provenance,
                }
                for sample in current
            ],
        }

    calls = []
    fail_second_attempts = 2

    def run_one(samples_one, _config, _output, *, existing_assignments=None,
                profiles_only=False, progress_index=1, progress_total=None, **_private):
        nonlocal fail_second_attempts
        sample = samples_one[0]
        calls.append(sample.sample_id)
        if sample.sample_id == "second" and fail_second_attempts:
            fail_second_attempts -= 1
            raise RuntimeError("synthetic native process failure")
        return [{
            "sample_id": sample.sample_id,
            "species": sample.species,
            "assembly_sha256": hashlib.sha256(sample.assembly.read_bytes()).hexdigest(),
            "typing_mode": "full",
            "typing_status": "typed",
            "typing_provenance": provenance,
            "cgmlst_scheme": "klebsiella_1",
            "cgmlst_scheme_version": "fixture-v1",
            "cgmlst_profile": {"locus-a": "1"},
            "cgmlst_loci": ["locus-a"],
        }]

    monkeypatch.setattr("chronoclade.query_typing._typing_recovery_identity", identity)
    monkeypatch.setattr("chronoclade.query_typing._type_query_assemblies_once", run_one)
    monkeypatch.setattr("chronoclade.query_typing.load_typing_config", lambda *_a, **_k: {})
    output = tmp_path / "typing"
    with pytest.raises(RuntimeError, match="synthetic native"):
        type_query_assemblies(samples, config, output)

    partial = output / "query_typing.partial.json"
    assert [row["sample_id"] for row in __import__("json").loads(partial.read_text())["assignments"]] == ["first"]
    saved_raw = output / "attempts"
    assert saved_raw.is_dir()
    assert calls == ["first", "second"]

    import json
    manifest_path = output / "typing-recovery.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["completed_sample_ids"] = []  # simulate SIGINT before advisory manifest update
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(RuntimeError, match="synthetic native"):
        type_query_assemblies(samples, config, output, resume=True)
    partial_rows = json.loads(partial.read_text())["assignments"]
    assert [row["sample_id"] for row in partial_rows] == ["first"]

    result = type_query_assemblies(samples, config, output, resume=True)
    assert [row["sample_id"] for row in result] == ["first", "second"]
    assert result[0]["typing_reused"] is True
    assert calls == ["first", "second", "second", "second"]
    assert (output / "query_typing.json").is_file()
    assert __import__("json").loads((output / "typing-recovery.json").read_text())["status"] == "complete"


def test_changed_input_and_corrupt_checkpoint_fail_before_native_calls(tmp_path, monkeypatch):
    sample = _sample(tmp_path, "query")
    config = tmp_path / "config.json"
    config.write_text("{}\n")
    provenance = {"execution": {"identity_sha256": "fixture"}}

    def identity(current, _config, _profiles_only):
        return {
            "schema": "chronoclade.query-typing-recovery-identity",
            "schema_version": 1,
            "mode": "full",
            "config_sha256": "a" * 64,
            "execution": provenance["execution"],
            "samples": [{
                "sample_id": item.sample_id,
                "species": "klebsiella pneumoniae",
                "assembly_sha256": hashlib.sha256(item.assembly.read_bytes()).hexdigest(),
                "typing_provenance": provenance,
            } for item in current],
        }

    calls = []
    def run_one(samples_one, *_args, **_kwargs):
        calls.append(samples_one[0].sample_id)
        return [{
            "sample_id": samples_one[0].sample_id,
            "species": samples_one[0].species,
            "assembly_sha256": hashlib.sha256(samples_one[0].assembly.read_bytes()).hexdigest(),
            "typing_mode": "full", "typing_status": "typed", "typing_provenance": provenance,
        }]

    monkeypatch.setattr("chronoclade.query_typing._typing_recovery_identity", identity)
    monkeypatch.setattr("chronoclade.query_typing._type_query_assemblies_once", run_one)
    monkeypatch.setattr("chronoclade.query_typing.load_typing_config", lambda *_a, **_k: {})
    output = tmp_path / "typing"
    type_query_assemblies([sample], config, output)
    original_calls = list(calls)
    sample.assembly.write_text(">query\nTGCA\n")
    with pytest.raises(QueryTypingError, match="identity differs"):
        type_query_assemblies([sample], config, output, resume=True)
    assert calls == original_calls
    sample.assembly.write_text(">query\nACGT\n")
    partial = output / "query_typing.partial.json"
    import json
    envelope = json.loads(partial.read_text())
    envelope["assignments"][0]["typing_status"] = "tampered"
    partial.write_text(json.dumps(envelope))
    with pytest.raises(QueryTypingError, match="integrity"):
        type_query_assemblies([sample], config, output, resume=True)
    assert calls == original_calls


def test_resume_accepts_atomic_zero_row_checkpoint(tmp_path, monkeypatch):
    sample = _sample(tmp_path, "first")
    config = tmp_path / "config.json"
    config.write_text("{}\n")
    provenance = {"execution": {"identity_sha256": "fixture"}}

    def identity(current, _config, _profiles_only):
        return {
            "schema": "chronoclade.query-typing-recovery-identity", "schema_version": 1,
            "mode": "full", "config_sha256": "a" * 64,
            "execution": provenance["execution"],
            "samples": [{"sample_id": item.sample_id, "species": "klebsiella pneumoniae",
                         "assembly_sha256": hashlib.sha256(item.assembly.read_bytes()).hexdigest(),
                         "typing_provenance": provenance} for item in current],
        }

    calls = []
    def run_one(samples, *_args, **_kwargs):
        calls.append(samples[0].sample_id)
        if len(calls) == 1:
            raise RuntimeError("synthetic failure before first checkpoint row")
        return [{"sample_id": samples[0].sample_id, "species": samples[0].species,
                 "assembly_sha256": hashlib.sha256(samples[0].assembly.read_bytes()).hexdigest(),
                 "typing_mode": "full", "typing_status": "typed",
                 "typing_provenance": provenance}]

    monkeypatch.setattr("chronoclade.query_typing._typing_recovery_identity", identity)
    monkeypatch.setattr("chronoclade.query_typing._type_query_assemblies_once", run_one)
    monkeypatch.setattr("chronoclade.query_typing.load_typing_config", lambda *_a, **_k: {})
    output = tmp_path / "typing"
    with pytest.raises(RuntimeError, match="before first"):
        type_query_assemblies([sample], config, output)
    assert __import__("json").loads((output / "query_typing.partial.json").read_text())["assignments"] == []
    # Simulate interruption between manifest creation and the atomic zero-row
    # checkpoint write. Only this explicit manifest state may resume without it.
    manifest_path = output / "typing-recovery.json"
    manifest = __import__("json").loads(manifest_path.read_text())
    manifest["checkpoint_initialized"] = False
    manifest_path.write_text(__import__("json").dumps(manifest))
    (output / "query_typing.partial.json").unlink()
    assert type_query_assemblies([sample], config, output, resume=True)[0]["sample_id"] == "first"
    assert calls == ["first", "first"]


def test_refuses_legacy_output_and_symlink_without_mutation(tmp_path):
    legacy = tmp_path / "legacy"
    legacy.mkdir()
    (legacy / "query_typing.json").write_text('{"assignments": []}\n')
    before = (legacy / "query_typing.json").read_bytes()
    with pytest.raises(QueryTypingError, match="already exists"):
        type_query_assemblies([], tmp_path / "unused.json", legacy)
    assert (legacy / "query_typing.json").read_bytes() == before

    destination = tmp_path / "owned"
    destination.mkdir()
    marker = destination / "marker"
    marker.write_text("untouched")
    link = tmp_path / "link"
    link.symlink_to(destination, target_is_directory=True)
    with pytest.raises(QueryTypingError, match="symlink"):
        type_query_assemblies([], tmp_path / "unused.json", link)
    assert marker.read_text() == "untouched"


def test_standalone_typing_rejects_symlinked_attempt_root(tmp_path, monkeypatch):
    sample = _sample(tmp_path, "query")
    config = tmp_path / "config.json"
    config.write_text("{}\n")
    provenance = {"execution": {"identity_sha256": "fixture"}}
    monkeypatch.setattr("chronoclade.query_typing.load_typing_config", lambda *_a, **_k: {})
    monkeypatch.setattr(
        "chronoclade.query_typing._typing_recovery_identity",
        lambda current, *_args: {
            "mode": "full", "execution": provenance["execution"],
            "samples": [{
                "sample_id": current[0].sample_id,
                "species": "klebsiella pneumoniae",
                "assembly_sha256": hashlib.sha256(current[0].assembly.read_bytes()).hexdigest(),
                "typing_provenance": provenance,
            }],
        },
    )
    output = tmp_path / "typing"
    output.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    marker = outside / "marker"
    marker.write_text("preserve")
    (output / "attempts").symlink_to(outside, target_is_directory=True)
    with pytest.raises(QueryTypingError, match="symlink"):
        type_query_assemblies([sample], config, output)
    assert marker.read_text() == "preserve"
    assert list(outside.iterdir()) == [marker]
