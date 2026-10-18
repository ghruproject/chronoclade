"""Authentic exported labels have a descriptive scope independent of native typing."""

import hashlib
import json
from pathlib import Path
import base64

import pytest

from chronoclade.context_partitions import partition_records
from chronoclade.pathogenwatch import content_hash
from chronoclade.pathogenwatch_hiercc import parse_hiercc_export
from chronoclade.profile_context_provider import _matching
from chronoclade.provider_hiercc import (
    DISCLOSURE,
    bind_exports,
    grouping_scope,
    normalized_server,
    reported_scope,
    source_snapshot,
)
from chronoclade.report_components.relatedness import scope_summary
from chronoclade.typing_scopes import typing_scope

FIXTURE = Path(__file__).parent / "fixtures/pathogenwatch-hiercc-provider-proof"
ID = "iLGx7naVnHXTvwPF8DEgAo"
PROVENANCE = json.loads((FIXTURE / "provenance.json").read_text())
RECEIPT = next(e for e in PROVENANCE["analysis_exports"] if e["download"] == "hclink")[
    "batch_receipts"
][0]
DATA = base64.b64decode((FIXTURE / "response.b64").read_text())
JOB = "hclink-562-4"


def proof(tmp_path):
    """Captured response; newly recorded source envelope is explicitly synthetic."""
    row = {
        "sample_id": "query",
        "source_genome_id": ID,
        "numeric_source_id": "108048",
        "species": "Escherichia coli",
        "organism_id": "562",
        "mlst_st": "131",
        "role": "input",
    }
    row.update(parse_hiercc_export(DATA, [ID], job=JOB)[0])
    cap = {
        "groups": [
            {"organismId": "562", "ids": [108048], "downloads": [{"name": "hclink", "job": JOB}]}
        ],
        "requested_numeric_ids": [108048],
        "retrieved_at": "2026-10-18T00:00:00Z",
        "request_contract": {
            "schema": "chronoclade.export-capability-request",
            "schema_version": 1,
            "base_url": "https://pathogen.watch",
            "requested_numeric_ids": [108048],
            "requested_source_ids": [ID],
            "row_scope": [{"species": "Escherichia coli", "organism_id": "562"}],
        },
    }
    cap["sha256"] = content_hash(cap)
    path = tmp_path / "export_capabilities.json"
    path.write_text(json.dumps(cap))
    (tmp_path / "response.bin").write_bytes(DATA)
    receipt = {
        "download": "hclink",
        "job": JOB,
        "base_url": "https://pathogen.watch",
        "endpoint": "/api/downloads/hclink",
        "route": "/api/downloads/hclink?job=" + JOB,
        "source_genome_ids": [ID],
        "numeric_request_ids": ["108048"],
        "returned_source_genome_ids": [ID],
        "missing_source_genome_ids": [],
        "status": "accepted",
        "retry_attempt": 1,
        "request_evidence": "recorded_exact_http_request",
        "response_sha256": hashlib.sha256(DATA).hexdigest(),
        "response_bytes": len(DATA),
        "saved_response_path": "response.bin",
    }
    provenance = {
        "analysis_exports": [
            {
                "download": "hclink",
                "job": JOB,
                "status": "downloaded",
                "requested_count": 1,
                "batch_receipts": [receipt],
            }
        ]
    }
    return row, path, provenance


def test_captured_labels_with_synthetic_complete_source_allow_descriptive_partition(tmp_path):
    row, path, prov = proof(tmp_path)
    bound = bind_exports([row], prov, path)[0]
    assert typing_scope(bound, "hiercc") is None
    assert bound["hiercc_database_version"] is None and bound["hiercc_database_sha256"] is None
    assert grouping_scope(bound) is not None
    context = dict(bound, sample_id="other", role="context")
    assert _matching(bound, context, "hiercc", "HC50")
    blocks, audit = partition_records([bound, context], hiercc_level="HC50")
    assert len(blocks) == 1 and len(blocks[0]["records"]) == 2 and not audit["unresolved_inputs"]
    assert blocks[0]["scheme_version"] is None
    assert blocks[0]["grouping_policy"] == "provider_reported_hiercc_v1"
    report = scope_summary({"hiercc_level": "HC50"}, {}, records=[bound, context])
    assert DISCLOSURE in report and "not assigned locally" in report
    assert source_snapshot([bound])["policy"] == "provider_reported_hiercc_v1"
    # Same namespace across separately observed accepted envelopes, snapshot differs.
    cap = json.loads(path.read_text())
    cap["retrieved_at"] = "2026-10-18T01:00:00Z"
    cap["sha256"] = content_hash({k: v for k, v in cap.items() if k != "sha256"})
    path.write_text(json.dumps(cap))
    later = bind_exports([row], prov, path)[0]
    assert grouping_scope(bound) == grouping_scope(later)
    assert source_snapshot([bound]) != source_snapshot([later])


@pytest.mark.parametrize(
    "mutation",
    ["raw", "receipt", "uuid", "numeric", "job", "capability", "native", "known_reference"],
)
def test_tampered_or_foreign_proof_fails_closed(tmp_path, mutation):
    row, path, prov = proof(tmp_path)
    receipt = prov["analysis_exports"][0]["batch_receipts"][0]
    if mutation == "raw":
        (tmp_path / "response.bin").write_bytes(b"wrong")
    elif mutation == "receipt":
        receipt["request_evidence"] = "derived_legacy_request"
    elif mutation == "uuid":
        row["source_genome_id"] = "unrequested"
    elif mutation == "numeric":
        row["numeric_source_id"] = "999"
    elif mutation == "job":
        receipt["job"] = "hclink-562-other"
    elif mutation == "capability":
        path.write_text(path.read_text().replace("108048", "108049"))
    elif mutation == "native":
        row["hiercc_source"] = "native_hclink"
    else:
        row["hiercc_database_sha256"] = "a" * 64
    if mutation == "uuid":
        assert grouping_scope(bind_exports([row], prov, path)[0]) is None
    else:
        with pytest.raises(ValueError):
            bind_exports([row], prov, path)


def test_arbitrary_rows_cannot_copy_a_serialized_scope_and_conflicts_do_not_match(tmp_path):
    row, path, prov = proof(tmp_path)
    valid = bind_exports([row], prov, path)[0]
    copied = json.loads(json.dumps(valid))
    assert reported_scope(copied) is None
    assert not _matching(valid, copied, "hiercc", "HC50")
    other = dict(valid, hiercc_database_sha256="a" * 64)
    assert not _matching(valid, other, "hiercc", "HC50")
    for field, value in [("hiercc_codes", {"HC50": "8"}), ("hiercc_status", "unassigned")]:
        assert not _matching(valid, dict(valid, **{field: value}), "hiercc", "HC50")
    assert normalized_server("https://PATHOGEN.WATCH:443/") == "https://pathogen.watch"
    assert normalized_server("https://next.pathogen.watch") != "https://pathogen.watch"


def test_real_new_exact_query_proof_preserves_unknown_release(tmp_path):
    import shutil

    source = tmp_path / "exact"
    shutil.copytree(FIXTURE, source)
    raw = source / RECEIPT["saved_response_path"]
    raw.parent.mkdir(parents=True)
    raw.write_bytes(DATA)
    row = json.loads((FIXTURE / "row.json").read_text())
    result = bind_exports([row], PROVENANCE, source / "export_capabilities.json")[0]
    assert grouping_scope(result) is not None
    assert typing_scope(result, "hiercc") is None
    assert result["hiercc_database_version"] is None
    assert result["hiercc_database_sha256"] is None
    assert result["source_genome_id"] == ID
    assert result["hiercc_codes"]["HC50"] == "7"


def test_old_actual_capability_does_not_invent_missing_recorded_request(tmp_path):
    row = json.loads((FIXTURE / "row.json").read_text())
    with pytest.raises(ValueError, match="recorded exact capability"):
        bind_exports([row], PROVENANCE, FIXTURE / "legacy_export_capabilities.json")


@pytest.mark.parametrize("publication_interrupt", [False, True])
@pytest.mark.parametrize("extra_context_locus", [False, True])
def test_normal_fake_provider_prepare_context_and_portable_proof(
    tmp_path, monkeypatch, publication_interrupt, extra_context_locus
):
    """Ordinary APIs; all responses synthetic and every request is intercepted."""
    import gzip
    import shutil
    from chronoclade import profile_inputs
    from chronoclade.datasets import load_dataset
    from chronoclade.pathogenwatch import PathogenwatchClient
    from chronoclade.prepare_stage import run_prepare
    from chronoclade.profile_context_provider import discover_profile_context, dataset_records

    other = "bLGx7naVnHXTvwPF8DEgAo"
    details = {
        ID: {
            "uuid": ID,
            "numeric_source_id": "108048",
            "species": "Escherichia coli",
            "organismId": "562",
            "mlst_st": "131",
            "projectAccess": "PUBLIC",
            "country": "Vietnam",
        },
        other: {
            "uuid": other,
            "numeric_source_id": "108049",
            "species": "Escherichia coli",
            "organismId": "562",
            "mlst_st": "131",
            "projectAccess": "PUBLIC",
            "country": "UK",
        },
    }

    def transport(method, endpoint, *, body, params, headers):
        if endpoint == "/api/organisms/supported":
            return [{"fullName": "Escherichia coli", "organismId": "562"}]
        if endpoint == "/api/genomes/details":
            return details[params["id"]]
        assert method == "POST" and endpoint == "/api/genomes/group"
        return [
            {
                "organismId": "562",
                "ids": body["ids"],
                "downloads": [
                    {"name": "hclink", "job": JOB},
                    {"name": "cgmlst", "job": "cgmlst-562-1"},
                ],
            }
        ]

    download_calls = []

    def download(route, **kwargs):
        download_calls.append((route, dict(kwargs.get("body", {}))))
        numeric = kwargs["body"]["ids"]
        ident = ID if numeric == "108048" else other
        assert numeric in {"108048", "108049"}
        if route.startswith("/api/downloads/hclink?"):
            return gzip.compress(gzip.decompress(DATA).replace(ID.encode(), ident.encode()))
        assert route == "/api/downloads/cgmlst?job=cgmlst-562-1"
        calls = {"a": "1", "b": "2"}
        if extra_context_locus:
            calls.update({f"gene{i}": "1" for i in range(8)})
            if ident == other:
                calls["observed_context_extra"] = "3"
        return (
            "Genome ID,Gene,Allele ID\n"
            + "".join(f"{ident},{gene},{allele}\n" for gene, allele in calls.items())
        ).encode()

    client = PathogenwatchClient(api_key="synthetic-fixture", transport=transport)

    def catalogue(path, **kwargs):
        assert kwargs == {"organism_id": "562", "st": "131"}
        value = {
            "rows": [
                {
                    **details[other],
                    "source_genome_id": other,
                    "sample_id": "PW_" + other,
                    "organism_id": "562",
                    "role": "context",
                    "aliases": [other],
                }
            ],
            "provenance": {"complete": True, "expected_count": 1},
        }
        path.write_text(json.dumps(value))
        return value

    client.freeze_catalogue = catalogue
    monkeypatch.setattr(profile_inputs, "request_download", download)
    query_path = tmp_path / "input.txt"
    query_path.write_text(ID + "\n")
    prepared = run_prepare(
        query_path,
        tmp_path / "prepared",
        input_kind="accessions",
        species="Escherichia coli",
        client=client,
        enrich_metadata=False,
    )
    readiness = json.loads((prepared.dataset_manifest.parent / "readiness.json").read_text())
    assert readiness["readiness"] == "incomplete"
    expected_audit = None
    if publication_interrupt:
        from chronoclade import prepare_stage

        original_prepare = prepare_stage.run_prepare

        def interrupted(*args, **kwargs):
            if Path(args[1]) == tmp_path / "context":
                raise OSError("synthetic publication interrupted")
            return original_prepare(*args, **kwargs)

        monkeypatch.setattr(prepare_stage, "run_prepare", interrupted)
        with pytest.raises(OSError, match="publication interrupted"):
            discover_profile_context(
                prepared.dataset_manifest, tmp_path / "context", hiercc_level="HC50", client=client
            )
        expected_document = json.loads((tmp_path / "context_retrieval/profiles.json").read_text())
        expected_audit = expected_document["parameters"]["profile_context"][
            "provider_reported_hiercc"
        ]
        assert expected_audit["reported_samples"] == 2
        exported_before_resume = list(download_calls)
        monkeypatch.setattr(prepare_stage, "run_prepare", original_prepare)
    result = discover_profile_context(
        prepared.dataset_manifest,
        tmp_path / "context",
        hiercc_level="HC50",
        client=client,
        resume=publication_interrupt,
    )
    if publication_interrupt:
        assert download_calls == exported_before_resume
        resumed_document = json.loads((tmp_path / "context_retrieval/profiles.json").read_text())
        assert resumed_document["records"] == expected_document["records"]
    dataset = load_dataset(result.dataset_manifest)
    rows = dataset_records(dataset, source_root=result.dataset_manifest)
    assert len(rows) == 2 and all(grouping_scope(row) for row in rows)
    if extra_context_locus:
        assert dataset.profiles[0].catalogue.complete is False
        assert len(dataset.profiles[0].catalogue.loci) == 11
        assert dataset.profiles[0].profile("PW_" + ID)["observed_context_extra"] is None
        assert dataset.profiles[0].profile("PW_" + other)["observed_context_extra"] == "3"
    blocks, audit = partition_records(rows, hiercc_level="HC50")
    assert len(blocks) == 1 and not audit["unresolved_inputs"]
    assert {r["source_genome_id"] for r in rows} == {ID, other}
    assert all(r["hiercc_database_version"] is None for r in rows)
    assert all(not r.get("assembly_reference") for r in rows)
    assert _matching(rows[0], rows[1], "hiercc", "HC50")
    context = json.loads(result.audit_path.read_text())["provider_context"]
    assert context["provider_reported_hiercc"]["reported_samples"] == 2
    if expected_audit is not None:
        assert context["provider_reported_hiercc"] == expected_audit
    assert context["context_assemblies"] == "not_requested"
    moved = tmp_path / "moved"
    shutil.copytree(result.dataset_manifest.parent, moved)
    assert source_snapshot(
        dataset_records(load_dataset(moved / "dataset.json"), source_root=moved)
    ) == source_snapshot(rows)

    from chronoclade.cgmlst.workflow import run_cgmlst

    run_cgmlst(
        result.dataset_manifest,
        tmp_path / "analysis",
        hiercc_level="HC50",
        context_size=1,
        replicates=1,
        bootstrap_replicates=0,
    )
    numerical = json.loads(
        next((tmp_path / "analysis").glob("*/numerical_analysis.json")).read_text()
    )
    nearest = numerical["analysis"]["nearest_neighbours"]
    assert len(nearest) == 1
    assert nearest[0]["query_id"] == "PW_" + ID
    assert nearest[0]["context_id"] == "PW_" + other
    assert nearest[0]["allele_differences"] == 0
    assert nearest[0]["shared_called_loci"] == (10 if extra_context_locus else 2)
    assert nearest[0]["tied_neighbours"] == 1
    stage_report = (tmp_path / "analysis/index.html").read_text()
    assert DISCLOSURE in stage_report
    block_reports = list((tmp_path / "analysis").glob("*/profile_report.html"))
    assert len(block_reports) == 1 and DISCLOSURE in block_reports[0].read_text()
    saved = json.loads((tmp_path / "analysis/cgmlst.json").read_text())
    assert (
        saved["parameters"]["provider_reported_hiercc"]["source_snapshot_sha256"]
        == source_snapshot(rows)["source_snapshot_sha256"]
    )
    # Unchanged source reaches the ordinary completed-run boundary, rather than
    # falsely failing JSON tuple/list equality in the resume fingerprint.
    with pytest.raises(ValueError, match="already complete"):
        run_cgmlst(
            result.dataset_manifest,
            tmp_path / "analysis",
            hiercc_level="HC50",
            context_size=1,
            replicates=1,
            bootstrap_replicates=0,
            resume=True,
        )
    # A changed composite source snapshot cannot silently reuse numerical work.
    checkpoint = tmp_path / "analysis/resume.json"
    old = json.loads(checkpoint.read_text())
    old["parameters"]["provider_reported_hiercc"]["source_snapshot_sha256"] = "0" * 64
    checkpoint.write_text(json.dumps(old))
    with pytest.raises(ValueError, match="source or parameters differ"):
        run_cgmlst(
            result.dataset_manifest,
            tmp_path / "analysis",
            hiercc_level="HC50",
            context_size=1,
            replicates=1,
            bootstrap_replicates=0,
            resume=True,
        )


@pytest.mark.parametrize("code,assigned", [("0", True), ("", False)])
def test_blank_is_not_zero_and_no_code_is_inferred(tmp_path, code, assigned):
    row, path, prov = proof(tmp_path)
    raw = f"Genome ID,Genome Name,Closest ST,Distance to closest ST,HC50\n{ID},fixture,1,0,{code}\n".encode()
    (tmp_path / "response.bin").write_bytes(raw)
    receipt = prov["analysis_exports"][0]["batch_receipts"][0]
    receipt.update(response_sha256=hashlib.sha256(raw).hexdigest(), response_bytes=len(raw))
    row.update(parse_hiercc_export(raw, [ID], job=JOB)[0])
    bound = bind_exports([row], prov, path)[0]
    blocks, _ = partition_records([bound], hiercc_level="HC50")
    assert bool(blocks) == assigned
    assert bound["hiercc_codes"] == ({"HC50": "0"} if assigned else {})


def test_other_provider_is_not_promoted_into_pathogenwatch_namespace(tmp_path):
    row, path, prov = proof(tmp_path)
    cap = json.loads(path.read_text())
    cap["request_contract"]["base_url"] = "https://other.invalid"
    cap["sha256"] = content_hash({k: v for k, v in cap.items() if k != "sha256"})
    path.write_text(json.dumps(cap))
    prov["analysis_exports"][0]["batch_receipts"][0]["base_url"] = "https://other.invalid"
    result = bind_exports([row], prov, path)[0]
    assert grouping_scope(result) is None
    assert "outside" in result["hiercc_reported_grouping_limitation"]


def test_scope_checkpoint_revalidates_raw_proof_not_serialized_token(tmp_path):
    from chronoclade.provider_hiercc import rebind_scope_exports

    row, path, provenance = proof(tmp_path)
    original = bind_exports([row], provenance, path)
    serialized = json.loads(json.dumps(original))
    assert reported_scope(serialized[0]) is None
    restored = rebind_scope_exports(serialized, provenance, path)
    assert source_snapshot(restored) == source_snapshot(original)
    (tmp_path / "response.bin").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="hash/size"):
        rebind_scope_exports(serialized, provenance, path)
