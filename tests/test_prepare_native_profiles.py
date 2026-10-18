"""Ordinary prepare can call frozen cgMLST profiles without lineage references."""

import json
import sys

import pytest

from chronoclade import prepare_provider, query_typing
from chronoclade.pathogenwatch import PathogenwatchClient
from test_query_typing import indexed_scheme, profile


def inputs(tmp_path, *, ready=False):
    version = "snapshot:" + "a" * 64
    indexed_scheme(tmp_path / "index", version=version, genes=("a", "b", "c", "d"))
    config = tmp_path / "typing.json"
    config.write_text(json.dumps({
        "ready": ready,
        "organisms": {"Escherichia coli": {"cgmlst": {
            "command": [sys.executable], "scheme": "ecoli_1", "scheme_version": version,
            "tool_version": "fixture-1", "index_dir": "index",
        }}},
    }))
    (tmp_path / "query.fa").write_text(">query\nACGT\n")
    csv = tmp_path / "queries.csv"
    csv.write_text("sample_id,assembly,species,mlst_st\nquery,query.fa,Escherichia coli,95\n")
    return csv, config


def resolve(csv, config, work, **kwargs):
    return prepare_provider.resolve_queries(
        csv, work, input_kind="assemblies", typing_config=config,
        client=PathogenwatchClient(api_key=""), enrich_metadata=False, **kwargs,
    )


def caller(monkeypatch):
    calls = []
    def run(command, data, destination):
        calls.append(command)
        assert "--scheme=ecoli_1" in command
        assert b"ACGT" in data
        return profile()
    monkeypatch.setattr(query_typing, "_run", run)
    monkeypatch.setattr(prepare_provider, "assign_existing_profiles",
                        lambda *a, **k: pytest.fail("Optional assigner must not run"))
    return calls


@pytest.mark.parametrize("ready", [False, True])
def test_prepare_automatically_uses_valid_caller_without_assignment(tmp_path, monkeypatch, ready):
    csv, config = inputs(tmp_path, ready=ready)
    calls = caller(monkeypatch)
    rows, _, _, ledger, _, operations = resolve(csv, config, tmp_path / "work")
    assert len(calls) == 1
    assert rows[0]["cgmlst_profile"] == {"a": "1"}
    assert rows[0]["cgmlst_loci"] == ["a", "b", "c", "d"]
    assert rows[0]["cgmlst_locus_universe_complete"] is True
    assert rows[0]["typing_mode"] == "profiles-only"
    assert "hiercc_codes" not in rows[0]
    assert operations["reference_database"] == "ready_validated_configuration"
    assert operations["native_typing_mode"] == "profiles-only"
    assert operations["typing"] == "run_missing_profile_inputs"
    assert operations["lineage_assignment"] == "not_run_optional_assignment_unavailable"
    assert operations["lineage_assignment_unavailable_reason"]
    assert ledger[0]["status"] == "ready"
    assert ledger[0]["lineage_status"] == "unavailable_optional"
    exported = json.loads((tmp_path / "work/native_typing/query_typing.json").read_text())
    assert exported["typing_mode"] == "profiles-only"


def test_caller_only_prepare_resumes_its_native_checkpoint(tmp_path, monkeypatch):
    csv, config = inputs(tmp_path)
    calls = caller(monkeypatch)
    cache = tmp_path / "checkpoint"
    first = resolve(csv, config, tmp_path / "first", native_checkpoint=cache,
                    native_attempt=tmp_path / "first-attempt")
    identity = json.loads((cache / "typing-recovery.json").read_text())["identity"]
    assert identity["mode"] == "profiles-only"
    assert "assignment" not in identity["samples"][0]["typing_provenance"]
    second = resolve(csv, config, tmp_path / "second", native_checkpoint=cache,
                     native_attempt=tmp_path / "second-attempt", resume_typing=True)
    assert len(calls) == 1
    assert second[0][0]["typing_reused"] is True
    assert second[0][0]["cgmlst_profile"] == first[0][0]["cgmlst_profile"]
    portable = json.loads((tmp_path / "second/native_typing/query_typing.json").read_text())
    assert portable["typing_mode"] == "profiles-only"


def test_existing_profile_does_not_trigger_assembly_retyping_for_optional_lineage(tmp_path, monkeypatch):
    csv, config = inputs(tmp_path)
    calls = caller(monkeypatch)
    first = resolve(csv, config, tmp_path / "first")
    exported = tmp_path / "first/native_typing/query_typing.json"
    second = resolve(csv, config, tmp_path / "second", query_typing=exported)
    assert len(calls) == 1
    assert second[0][0]["cgmlst_profile"] == first[0][0]["cgmlst_profile"]
    assert second[-1]["typing"] == "not_run"
    assert second[3][0]["status"] == "ready"


def test_invalid_caller_never_claims_profile_only_readiness(tmp_path, monkeypatch):
    csv, config = inputs(tmp_path)
    calls = caller(monkeypatch)
    (tmp_path / "index/allele.db").unlink()
    rows, _, _, ledger, _, operations = resolve(csv, config, tmp_path / "work")
    assert not calls
    assert not rows[0].get("cgmlst_profile")
    assert operations["reference_database"] == "unavailable_configured_database_or_tool_not_ready"
    assert operations["reference_database_unavailable_reason"]
    assert ledger[0]["status"] == "incomplete"


def test_explicit_full_mode_does_not_fall_back(tmp_path, monkeypatch):
    csv, config = inputs(tmp_path)
    calls = caller(monkeypatch)
    rows, _, _, ledger, _, operations = resolve(csv, config, tmp_path / "work", profiles_only=False)
    assert not calls
    assert operations["typing"] == "not_run"
    assert ledger[0]["status"] == "incomplete"
    assert not rows[0].get("cgmlst_profile")
