"""Context dispatch uses shared HClink import; unknown scope is not fabricated."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from chronoclade import profile_context_provider as provider, profile_inputs
from chronoclade.context_partitions import partition_records
from chronoclade.datasets import (
    DatasetError,
    LocusCatalogue,
    from_profile_records,
    load_dataset,
    write_dataset,
)
from chronoclade.pathogenwatch_hiercc import parse_hiercc_export
from chronoclade.typing_scopes import typing_scope

ID = "iLGx7naVnHXTvwPF8DEgAo"
RAW = Path(__file__).parent / "fixtures/pathogenwatch-hclink-562-4-id108048.csv.gz"
RAW_SHA = "c5ccab9dabf1f1ab4a831ea385c2819bac8db9ba5d187ab63556c9a29a5fdd76"


def query():
    return {
        "sample_id": "synthetic-input",
        "role": "input",
        "species": "Escherichia coli",
        "mlst_st": "131",
        "cgmlst_scheme": "fixture-profile",
        "cgmlst_scheme_version": "fixture-v1",
        "cgmlst_profile": {"a": "1"},
        "cgmlst_loci": ["a"],
        "cgmlst_locus_universe_complete": True,
        "hiercc_scheme": "pathogenwatch:HierCC",
        "hiercc_scheme_version": "explicit-synthetic-reference",
        "hiercc_codes": {"HC50": "7"},
        "hiercc_status": "resolved",
    }


class Client:
    api_key = "fixture-only-not-a-credential"
    base_url = "https://fixture.invalid"

    def __init__(self, candidate):
        self.candidate = candidate
        self.requests = []

    def supported_organisms(self):
        return [{"fullName": "Escherichia coli", "organismId": "562"}]

    def freeze_catalogue(self, path, **kwargs):
        assert kwargs == {"organism_id": "562", "st": "131"}
        value = {
            "rows": [deepcopy(self.candidate)],
            "provenance": {"complete": True, "expected_count": 1},
        }
        path.write_text(json.dumps(value))
        return value

    def request_json(self, method, path, **kwargs):
        self.requests.append((method, path, kwargs))
        assert method == "POST" and path == "/api/genomes/group"
        assert kwargs["body"] == {"ids": [108048]}
        return [
            {
                "ids": [108048],
                "organismId": "562",
                "downloads": [{"name": "hclink", "job": "hclink-562-4"}],
            }
        ]


def prepared(tmp_path):
    return write_dataset(
        from_profile_records(
            [query()],
            [
                LocusCatalogue(
                    "fixture-profile", "fixture-v1", ("a",), "synthetic test", complete=True
                )
            ],
        ),
        tmp_path / "prepared",
    )


def candidate():
    row = query()
    row.update(
        sample_id="PW_" + ID,
        source_genome_id=ID,
        numeric_source_id=108048,
        organism_id="562",
        role="context",
        cgmlst_profile={"a": "2"},
    )
    return row


def test_context_requests_real_hclink_batch_and_preserves_unknown_reference(tmp_path, monkeypatch):
    data = RAW.read_bytes()
    assert hashlib.sha256(data).hexdigest() == RAW_SHA
    row = candidate()
    for key in list(row):
        if key.startswith("hiercc_"):
            row.pop(key)
    client = Client(row)
    calls = []

    def download(route, **kwargs):
        calls.append((route, kwargs["body"]))
        assert route == "/api/downloads/hclink?job=hclink-562-4"
        assert kwargs["body"] == {"ids": "108048"}
        return data

    monkeypatch.setattr(profile_inputs, "request_download", download)
    with pytest.raises(DatasetError, match="does not mean there are no public relatives"):
        provider.discover_profile_context(
            prepared(tmp_path), tmp_path / "context", hiercc_level="HC50", client=client
        )
    assert calls == [("/api/downloads/hclink?job=hclink-562-4", {"ids": "108048"})]
    receipts = list((tmp_path / "context_retrieval").rglob("receipt.json"))
    assert len(receipts) == 1
    receipt = json.loads(receipts[0].read_text())["receipt"]
    assert receipt["response_sha256"] == RAW_SHA
    assert (
        receipt["returned_source_genome_ids"] == [ID] and receipt["missing_source_genome_ids"] == []
    )
    refusal = json.loads((tmp_path / "context_retrieval/context_compatibility.json").read_text())
    rejected = refusal["rejected_candidates"][0]
    assert rejected["candidate_scope"]["scheme_version"] is None
    assert rejected["candidate_scope"]["database_sha256"] is None
    assert not (tmp_path / "context/dataset.json").exists()
    parsed = parse_hiercc_export(data, [ID], job="hclink-562-4")[0]
    assert parsed["hiercc_codes"]["HC50"] == "7" and parsed["hiercc_blank_levels"] == [
        "HC0",
        "HC2",
        "HC5",
        "HC10",
        "HC20",
    ]
    assert parsed["hiercc_scheme_version"] is None and typing_scope(parsed, "hiercc") is None
    assert hashlib.sha256(RAW.read_bytes()).hexdigest() == RAW_SHA


def test_explicit_synthetic_known_scope_context_uses_hiercc_partition(tmp_path, monkeypatch):
    row = candidate()
    client = Client(row)
    calls = []
    original = provider._grouped_analysis_exports

    def exports(rows, client, output, provenance, **kwargs):
        calls.append(kwargs.get("download_names"))
        return original(rows, client, output, provenance, **kwargs)

    monkeypatch.setattr(provider, "_grouped_analysis_exports", exports)
    monkeypatch.setattr(
        profile_inputs,
        "request_download",
        lambda *a, **k: pytest.fail("already typed fixture requires no download"),
    )
    result = provider.discover_profile_context(
        prepared(tmp_path), tmp_path / "context", hiercc_level="HC50", client=client
    )
    assert calls[0] == {"hclink"} and {"klebsiella-lincodes"} not in calls
    assert client.requests == []
    dataset = load_dataset(result.dataset_manifest)
    blocks, audit = partition_records(provider.dataset_records(dataset), hiercc_level="HC50")
    assert len(blocks) == 1 and blocks[0]["kind"] == "hiercc" and blocks[0]["level"] == "HC50"
    assert sorted(r["role"] for r in blocks[0]["records"]) == ["context", "input"]
    assert audit["unresolved_inputs"] == []
    assert {r["sample_id"] for r in dataset.samples} == {"synthetic-input", "PW_" + ID}
    assert (
        json.loads(result.audit_path.read_text())["provider_context"]["context_assemblies"]
        == "not_requested"
    )


def test_actual_export_without_reference_scope_still_cannot_partition():
    row = query()
    row.update(parse_hiercc_export(RAW.read_bytes(), [ID], job="hclink-562-4")[0])
    blocks, audit = partition_records([row], hiercc_level="HC50")
    assert blocks == [] and audit["unresolved_inputs"][0]["sample_id"] == "synthetic-input"


def test_shared_batch_identity_guard_rejects_unrequested_real_response(tmp_path, monkeypatch):
    row = candidate()
    row["source_genome_id"] = "different-requested-id"
    for key in list(row):
        if key.startswith("hiercc_"):
            row.pop(key)
    client = Client(row)
    monkeypatch.setattr(profile_inputs, "request_download", lambda *a, **k: RAW.read_bytes())
    monkeypatch.setattr(profile_inputs.time, "sleep", lambda *_: None)
    provenance = {}
    output = provider._grouped_analysis_exports(
        [row], client, tmp_path / "lineages", provenance, download_names={"hclink"}
    )
    assert not output[0].get("hiercc_codes")
    assert provenance["analysis_exports"][0]["status"] == "unavailable"
    rejected = provenance["analysis_exports"][0]["batch_receipts"]
    assert rejected and all(r["status"] == "rejected_identity" for r in rejected)
    assert all(r["returned_source_genome_ids"] is None for r in rejected)
    assert provenance["analysis_export_identity_retries"]
    assert not list((tmp_path / "lineages").rglob("receipt.json"))
    assert hashlib.sha256(RAW.read_bytes()).hexdigest() == RAW_SHA
