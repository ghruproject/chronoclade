"""Preparation operation labels describe this invocation, not assembly history."""

import json

import pytest

from chronoclade import prepare_provider
from chronoclade.pathogenwatch import PathogenwatchClient
from chronoclade.prepare_stage import _report, run_prepare


@pytest.mark.parametrize("enrich", [False, True])
def test_local_prepare_does_not_claim_provider_retrieval(tmp_path, monkeypatch, enrich):
    assembly = tmp_path / "query.fa"
    assembly.write_text(">query\nACGT\n")
    metadata = tmp_path / "assemblies.csv"
    metadata.write_text(
        "sample_id,assembly,species,run_accessions\nquery,query.fa,Klebsiella pneumoniae,ERR123\n"
    )
    client = PathogenwatchClient(api_key="")

    def no_provider_request(*args, **kwargs):
        pytest.fail("Local assembly preparation must not request provider data")

    monkeypatch.setattr(client, "request_json", no_provider_request)
    ena_calls = []

    def ena(rows):
        ena_calls.extend(row["sample_id"] for row in rows)
        return (
            [],
            [],
            [{"sample_id": "query", "operation": "ena_metadata", "status": "unavailable"}],
        )

    monkeypatch.setattr(prepare_provider, "enrich_exact_ena", ena)
    result = run_prepare(
        metadata,
        tmp_path / "prepared",
        input_kind="assemblies",
        client=client,
        enrich_metadata=enrich,
    )
    audit = json.loads(result.audit_path.read_text())
    assert audit["operations"]["provider_retrieval"] == "not_run_local_assemblies"
    assert audit["operations"]["typing"] == "not_run"
    assert ena_calls == (["query"] if enrich else [])
    ledger = [
        json.loads(line)
        for line in (result.dataset_manifest.parent / "retrieval.jsonl").read_text().splitlines()
    ]
    assert any(row["operation"] == "ena_metadata" for row in ledger) is enrich
    html = (result.dataset_manifest.parent / "report.html").read_text()
    assert "Provider retrieval: Not requested; local assemblies supplied" in html
    assert "run_exact_queries" not in html


def test_explicit_local_provider_ids_keep_exact_retrieval_attempt_status(tmp_path, monkeypatch):
    (tmp_path / "q.fa").write_text(">q\nACGT\n")
    metadata = tmp_path / "assemblies.csv"
    metadata.write_text(
        "sample_id,assembly,species,source_genome_id,numeric_source_id\n"
        "q,q.fa,Klebsiella pneumoniae,exact-public-id,123\n"
    )
    attempts = []

    def exports(rows, client, output, provider, **kwargs):
        attempts.append(output.name)
        return rows

    monkeypatch.setattr(prepare_provider, "_grouped_analysis_exports", exports)
    result = prepare_provider.resolve_queries(
        metadata,
        tmp_path / "work",
        input_kind="assemblies",
        client=PathogenwatchClient(api_key="fixture-only"),
        enrich_metadata=False,
    )
    assert attempts == ["query_exports_grouped", "query_lineages", "query_hiercc"]
    assert result[-1]["provider_retrieval"] == "run_exact_queries"


@pytest.mark.parametrize(
    "status,label",
    [
        ("run_exact_queries", "Exact query retrieval attempted"),
        ("not_run_offline", "Not requested; frozen evidence imported"),
    ],
)
def test_existing_operation_statuses_have_readable_labels(status, label):
    html = _report(
        {
            "samples": [],
            "operations": {
                "provider_retrieval": status,
                "typing": "not_run",
                "reference_database": "unavailable",
            },
        }
    )
    assert f"Provider retrieval: {label}" in html


@pytest.mark.parametrize(
    "typing,reference,expected",
    [
        (
            "run_missing_profile_or_lineage_inputs",
            "ready_validated_configuration",
            (
                "Typing: Native typing run where needed",
                "Reference database: Configured and validated",
            ),
        ),
        (
            "not_run_frozen_import",
            "not_checked_not_required_for_import",
            (
                "Typing: Not run; frozen typing imported",
                "Reference database: Not checked; not required for imported evidence",
            ),
        ),
        (
            "not_run",
            "unavailable_not_configured",
            ("Typing: Not run", "Reference database: Not configured"),
        ),
        (
            "not_run",
            "unavailable_configured_database_or_tool_not_ready",
            (
                "Typing: Not run",
                "Reference database: Unavailable; configured database or tool is not ready",
            ),
        ),
    ],
)
def test_existing_typing_and_reference_statuses_readably_preserve_evidence(
    typing, reference, expected
):
    audit = {
        "samples": [],
        "operations": {
            "provider_retrieval": "not_run_local_assemblies",
            "typing": typing,
            "reference_database": reference,
        },
    }
    before = json.dumps(audit, sort_keys=True)
    html = _report(audit)
    assert all(phrase in html for phrase in expected)
    assert json.dumps(audit, sort_keys=True) == before


def test_unknown_operation_statuses_are_visible_and_escaped():
    html = _report(
        {
            "samples": [],
            "operations": {
                "provider_retrieval": '<future-provider>&"',
                "typing": "<future-typing>",
                "reference_database": "<future-reference>",
            },
        }
    )
    assert "Provider retrieval: &lt;future-provider&gt;&amp;&quot;" in html
    assert "Typing: &lt;future-typing&gt;" in html
    assert "Reference database: &lt;future-reference&gt;" in html
    assert "<future-" not in html
