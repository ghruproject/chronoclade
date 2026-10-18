"""Exact provider joins, complete profile pools and readiness are independently audited."""

import json
import shutil
from dataclasses import asdict
from pathlib import Path

import pytest

from chronoclade.datasets import DatasetError, LocusCatalogue, load_dataset
from chronoclade.pathogenwatch import PathogenwatchClient
from chronoclade.prepare_provider import (
    _resolve_accessions,
    catalogues_for_records,
    enrich_exact_ena,
)
from chronoclade.prepare_stage import run_prepare
from chronoclade.profile_context_provider import discover_profile_context

SPECIES = "Klebsiella pneumoniae"
COLLECTION = "jX5cwsoyJ1KDquMqssUzAD"


def record(ident="query", **kwargs):
    return dict(
        sample_id=ident,
        species=SPECIES,
        role="input",
        mlst_st="39",
        cgmlst_scheme="test",
        cgmlst_scheme_version="v1",
        cgmlst_profile={"a": "1"},
        cglin_scheme="scgMLST629_S",
        cglin_scheme_version="lin-v1",
        cglin_status="complete",
        cglin_raw="0,0,1,2,3,4,5,6,7,8",
        **kwargs,
    )


def test_observed_catalogue_is_never_complete():
    row = record()
    catalogue = catalogues_for_records([row])[0]
    assert catalogue.complete is False
    assert catalogue.database_version is None
    row.update(cgmlst_loci=["a", "b"], cgmlst_locus_universe_complete=True)
    assert catalogues_for_records([row])[0].complete is True


def test_exact_ena_preserves_override_and_rejects_unrelated_identity():
    rows = [
        {
            "sample_id": "q",
            "run_accessions": ["ERR123"],
            "country": "Greece",
            "host": "",
            "collection_date": "",
        }
    ]
    provenance, conflicts, ledger = enrich_exact_ena(
        rows,
        fetch=lambda _: [
            {
                "run_accession": "ERR123",
                "sample_accession": "SAMEA1",
                "country": "UK",
                "host": "human",
                "collection_date": "2019",
            }
        ],
    )
    assert rows[0]["country"] == "Greece"
    assert rows[0]["host"] == "human" and rows[0]["collection_date"] == "2019"
    assert conflicts[0]["value"] == "UK"
    assert any(p["selected"] and p["field"] == "host" for p in provenance)
    rows = [{"sample_id": "q", "run_accessions": ["ERR123"]}]
    _, _, ledger = enrich_exact_ena(
        rows, fetch=lambda _: [{"run_accession": "ERR124", "country": "UK"}]
    )
    assert "country" not in rows[0] and ledger[0]["status"] == "failed"


def test_accession_resolution_saves_exact_decoded_search_and_detail_evidence(tmp_path):
    import hashlib

    class Client:
        def supported_organisms(self):
            return [{"fullName": SPECIES, "organismId": 1280}]

        def request_json(self, method, endpoint, *, body=None, params=None):
            if endpoint == "/api/search/genomes":
                assert method == "POST"
                assert body == {"organismId": "1280", "searchText": "SRR123", "qc": [True, False]}
                assert params == {"limit": 100, "sort": "id"}
                return {"meta": {"count": 1}, "genomes": [{"uuid": "source-A"}]}
            assert (method, endpoint, params) == (
                "GET",
                "/api/genomes/details",
                {"id": "source-A"},
            )
            return {
                "uuid": "source-A",
                "species": SPECIES,
                "organismId": "1280",
                "runAccession": "SRR123",
                "projectAccess": "PUBLIC",
            }

    query = tmp_path / "accessions.txt"
    query.write_text("SRR123\n")
    work = tmp_path / "work"
    rows = _resolve_accessions(query, SPECIES, Client(), work=work)
    assert rows[0]["source_genome_id"] == "source-A"
    receipt_path = work / "accession_resolution" / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    assert receipt["hash_definition"].startswith("sha256 of UTF-8 canonical decoded JSON")
    item = receipt["records"][0]
    assert (item["reported_count"], item["returned_count"], item["result_ids"]) == (
        1,
        1,
        ["source-A"],
    )
    assert item["request"]["path"] == "0000-search-request.json"
    assert item["response"]["path"] == "0000-search-response.json"
    assert item["details"][0]["source_genome_id"] == "source-A"
    for evidence in (item["request"], item["response"], item["details"][0]):
        saved = receipt_path.parent / evidence["path"]
        assert saved.is_file()
        assert evidence["bytes"] == saved.stat().st_size
    response = json.loads((receipt_path.parent / item["response"]["path"]).read_text())
    canonical = json.dumps(
        response, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")
    assert item["response"]["canonical_json_sha256"] == hashlib.sha256(canonical).hexdigest()


def test_normal_prepare_publishes_successful_accession_receipt_and_artifact_hashes(tmp_path):
    import hashlib

    from chronoclade.prepare_stage import run_prepare

    source_id = "ooy8fFUy3fXqa7WaWX15eS"

    def transport(method, endpoint, *, body, params, headers):
        assert "X-API-Key" not in headers
        if endpoint == "/api/organisms/supported":
            return [{"fullName": SPECIES, "organismId": 1280}]
        if endpoint == "/api/search/genomes":
            return {"meta": {"count": 1}, "genomes": [{"uuid": source_id}]}
        if endpoint == "/api/genomes/details":
            return {
                "uuid": source_id,
                "species": SPECIES,
                "organismId": "1280",
                "runAccession": "SRR123",
                "projectAccess": "PUBLIC",
            }
        raise AssertionError((method, endpoint, body, params))

    query = tmp_path / "accessions.txt"
    query.write_text("SRR123\n")
    result = run_prepare(
        query,
        tmp_path / "published",
        input_kind="accessions",
        species=SPECIES,
        client=PathogenwatchClient(api_key="", transport=transport),
        enrich_metadata=False,
    )
    bundle = result.dataset_manifest.parent
    receipt_rel = "sources/provider/accession_resolution/receipt.json"
    receipt_path = bundle / receipt_rel
    assert receipt_path.is_file()
    provider = load_dataset(result.dataset_manifest).parameters["provider_prepare"]
    assert provider["accession_resolution_evidence"]["path"] == "accession_resolution/receipt.json"
    assert provider["accession_resolution_evidence"]["sha256"] == hashlib.sha256(
        receipt_path.read_bytes()
    ).hexdigest()
    prepare = json.loads((bundle / "prepare.json").read_text())
    artifact = next(item for item in prepare["artifacts"] if item["path"] == receipt_rel)
    assert artifact["sha256"] == hashlib.sha256(receipt_path.read_bytes()).hexdigest()
    assert artifact["bytes"] == receipt_path.stat().st_size


@pytest.mark.parametrize(
    "response, message",
    [
        ({"meta": {"count": 2}, "genomes": [{"uuid": "A"}]}, "Accession search is incomplete"),
        (
            {"meta": {"count": 1}, "genomes": [{"uuid": "A"}, {"uuid": "B"}]},
            "Accession search is incomplete",
        ),
    ],
)
def test_accession_search_rejection_keeps_captured_response_in_work(
    tmp_path, response, message
):
    class Client:
        def supported_organisms(self):
            return [{"fullName": SPECIES, "organismId": 1280}]

        def request_json(self, *_args, **_kwargs):
            return response

    query = tmp_path / "accessions.txt"
    query.write_text("SRR123\n")
    work = tmp_path / "work"
    with pytest.raises(DatasetError, match=message):
        _resolve_accessions(query, SPECIES, Client(), work=work)
    receipt = json.loads((work / "accession_resolution" / "receipt.json").read_text())
    assert receipt["records"][0]["status"] == "rejected_incomplete_search"
    response_path = work / "accession_resolution" / receipt["records"][0]["response"]["path"]
    assert json.loads(response_path.read_text()) == response


def test_accession_ambiguity_preserves_all_exact_detail_receipts(tmp_path):
    from chronoclade.profile_inputs import ProfileInputError

    identifiers = ["ooy8fFUy3fXqa7WaWX15eS", "jX5cwsoyJ1KDquMqssUzAD"]

    class Client:
        def supported_organisms(self):
            return [{"fullName": SPECIES, "organismId": 1280}]

        def request_json(self, method, endpoint, *, body=None, params=None):
            if endpoint == "/api/search/genomes":
                return {
                    "meta": {"count": 2},
                    "genomes": [{"uuid": ident} for ident in identifiers],
                }
            ident = params["id"]
            return {
                "uuid": ident,
                "species": SPECIES,
                "organismId": "1280",
                "runAccession": "SRR123",
                "projectAccess": "PUBLIC",
            }

    query = tmp_path / "accessions.txt"
    query.write_text("SRR123\n")
    work = tmp_path / "work"
    with pytest.raises(ProfileInputError, match="2 exact matches"):
        _resolve_accessions(query, SPECIES, Client(), work=work)
    receipt = json.loads((work / "accession_resolution" / "receipt.json").read_text())
    record = receipt["records"][0]
    assert [detail["source_genome_id"] for detail in record["details"]] == identifiers
    assert all(
        (work / "accession_resolution" / detail["path"]).is_file()
        for detail in record["details"]
    )


def test_collection_live_prepare_without_typing_never_fetches_context(tmp_path):
    calls = []

    def transport(method, path, *, body, params, headers):
        calls.append(path)
        if path == "/api/collections/details":
            return {"uuid": COLLECTION, "size": 2, "organismName": SPECIES, "downloads": []}
        if path == "/api/collections/genomes":
            index = params["page"]
            return {
                "genomes": [
                    {"uuid": "source" + str(index), "name": "q" + str(index), "mlst": "39"}
                ],
                "hasMore": index == 1,
            }
        raise AssertionError(path)

    client = PathogenwatchClient(api_key="", transport=transport)
    result = run_prepare(
        "https://pathogen.watch/collections/" + COLLECTION + "-slug",
        tmp_path / "out",
        input_kind="collection",
        client=client,
        enrich_metadata=False,
    )
    dataset = load_dataset(result.dataset_manifest)
    assert len(dataset.samples) == 2 and not dataset.profiles
    audit = json.loads(result.audit_path.read_text())
    assert audit["operations"]["reference_database"] == "unavailable_not_configured"
    assert audit["operations"]["context_discovery"] == "not_run"
    assert calls.count("/api/collections/genomes") == 2
    assert all(r["role"] == "input" for r in dataset.samples)
    bundle = result.dataset_manifest.parent
    provider = dataset.parameters["provider_prepare"]
    snapshot_reference = provider["collection_snapshot"]
    assert not Path(snapshot_reference).is_absolute()
    assert (bundle / snapshot_reference).is_file()
    assert ".prepare-provider-" not in snapshot_reference
    snapshot = json.loads((bundle / snapshot_reference).read_text())
    assert snapshot["snapshot_sha256"] == provider["collection_snapshot_sha256"]
    saved_input = json.loads((bundle / "sources/input.json").read_text())
    assert saved_input["parameters"]["provider_prepare"]["collection_snapshot"] == snapshot_reference
    moved = tmp_path / "moved-collection"
    shutil.move(bundle, moved)
    assert (moved / snapshot_reference).is_file()
    relocated = load_dataset(moved / "dataset.json")
    assert relocated.parameters["provider_prepare"]["collection_snapshot"] == snapshot_reference
    assert json.loads((moved / snapshot_reference).read_text())["snapshot_sha256"] == provider[
        "collection_snapshot_sha256"
    ]


@pytest.mark.parametrize("mlst_st, expected_ready", [(None, False), ("39", True)])
def test_ordinary_assembly_prepare_contract_with_native_typing_stub(
    tmp_path, monkeypatch, mlst_st, expected_ready
):
    from chronoclade import pasteur_profile_normalization, prepare_provider, query_typing

    assembly = tmp_path / "query.fa"
    assembly.write_text(">query\nACGT\n")
    fields = ["sample_id", "assembly", "species"]
    values = ["query", "query.fa", SPECIES]
    if mlst_st is not None:
        fields.append("mlst_st")
        values.append(mlst_st)
    metadata = tmp_path / "assemblies.csv"
    metadata.write_text(",".join(fields) + "\n" + ",".join(values) + "\n")
    config = tmp_path / "typing-config.json"
    config.write_text("{}\n")
    calls = []

    monkeypatch.setattr(query_typing, "load_typing_config", lambda *_a, **_k: {"organisms": {}})
    monkeypatch.setattr(prepare_provider, "assign_existing_profiles", lambda *_a, **_k: set())
    # This dispatch/readiness control stubs native calls, not reference validity.
    # Import the binding module before mocking the caller loader so its module-
    # level import cannot retain a patched loader after this test ends.
    reference_stub = {"fixture": "native-reference-binding-not-tested"}
    monkeypatch.setattr(
        pasteur_profile_normalization, "reference_descriptor", lambda *_a, **_k: reference_stub
    )

    def stub_native_typing(samples, _config, _output, **_recovery_options):
        calls.extend(sample.sample_id for sample in samples)
        return [
            {
                "sample_id": sample.sample_id,
                "species": sample.species,
                "typing_status": "typed",
                "cgmlst_scheme": "test-cgmlst",
                "cgmlst_scheme_version": "v1",
                "cgmlst_profile": {"locus-a": "1", "locus-b": "2"},
                "cgmlst_loci": ["locus-a", "locus-b"],
                "cgmlst_locus_universe_complete": True,
                "cglin_scheme": "scgMLST629_S",
                "cglin_scheme_version": "v1",
                "cglin_status": "complete",
                "cglin_raw": "1.2.3.4.5.6.7.8.9.10",
            }
            for sample in samples
        ]

    monkeypatch.setattr(query_typing, "type_query_assemblies", stub_native_typing)
    result = run_prepare(
        metadata,
        tmp_path / "prepared",
        input_kind="assemblies",
        typing_config=config,
        enrich_metadata=False,
        client=PathogenwatchClient(api_key=""),
    )
    audit = json.loads(result.audit_path.read_text())
    sample = audit["samples"][0]
    assert load_dataset(result.dataset_manifest).parameters["native_typing_reference"] == reference_stub
    assert calls == ["query"]
    assert sample["profile_ready"] is True
    assert sample["lineage_ready"] is True
    assert sample["mlst_st_ready"] is expected_ready
    assert sample["readiness"] == ("ready" if expected_ready else "incomplete")
    if not expected_ready:
        assert "verified seven-gene MLST ST" in " ".join(sample["reasons"])


def test_live_local_assemblies_unready_database_retains_input(tmp_path):
    (tmp_path / "q.fa").write_text(">q\nACGT\n")
    (tmp_path / "metadata.csv").write_text(
        "sample_id,assembly,species,country\nq,q.fa,Klebsiella pneumoniae,Greece\n"
    )
    (tmp_path / "config.json").write_text('{"ready":false}')
    result = run_prepare(
        tmp_path / "metadata.csv",
        tmp_path / "out",
        input_kind="assemblies",
        typing_config=tmp_path / "config.json",
        enrich_metadata=False,
        client=PathogenwatchClient(api_key=""),
    )
    ds = load_dataset(result.dataset_manifest)
    assert ds.samples[0]["country"] == "Greece"
    assert (result.dataset_manifest.parent / ds.samples[0]["assembly_reference"]).is_file()
    audit = json.loads(result.audit_path.read_text())
    assert audit["samples"][0]["profile_ready"] is False
    assert audit["operations"]["reference_database"].startswith("unavailable_configured")


def test_context_retrieves_all_matching_profiles_without_assemblies(tmp_path, monkeypatch):
    import chronoclade.profile_context_provider as provider

    query = record()
    source = tmp_path / "profiles.json"
    source.write_text(json.dumps([query]))
    cats = tmp_path / "catalogues.json"
    cats.write_text(json.dumps([asdict(LocusCatalogue("test", "v1", ("a",), "fixture"))]))
    prepared = run_prepare(source, tmp_path / "prepared", catalogues=cats)
    public = []
    for index in range(120):
        row = record("PW_source" + str(index))
        row.update(
            source_genome_id="source" + str(index),
            origin="context",
            aliases=["source" + str(index)],
            numeric_source_id=index + 1,
            cgmlst_profile={},
        )
        public.append(row)
    # Two separately retrieved genomes belong to the same BioSample unit.
    public[0]["biosample_accessions"] = ["SAMEA123"]
    public[1]["biosample"] = "SAMEA123"
    unrelated = record("other")
    unrelated.update(source_genome_id="other", cglin_raw="0,0,1,2,9,4,5,6,7,8", cgmlst_profile={})
    public.append(unrelated)

    class Client:
        api_key = "fake"

        def supported_organisms(self):
            return [{"fullName": SPECIES, "organismId": "573"}]

        def freeze_catalogue(self, path, **kwargs):
            value = {
                "rows": public,
                "provenance": {"complete": True, "expected_count": 121, "search_page_count": 2},
            }
            path.write_text(json.dumps(value))
            return value

    requested = []

    def exports(rows, client, output, provenance, download_names=None, resume=False):
        if download_names is None:
            requested.extend(r["source_genome_id"] for r in rows)
            return [dict(r, cgmlst_profile={"a": "1"}) for r in rows]
        return rows

    monkeypatch.setattr(provider, "_grouped_analysis_exports", exports)
    result = discover_profile_context(
        prepared.dataset_manifest, tmp_path / "context", client=Client(), lin_level=5
    )
    ds = load_dataset(result.dataset_manifest)
    assert len(ds.samples) == 120 and len(requested) == 120
    assert ds.samples[0] == load_dataset(prepared.dataset_manifest).samples[0]
    audit = json.loads(result.audit_path.read_text())["provider_context"]
    assert audit["usable"] == 119 and audit["context_candidates"] == 120
    assert audit["context_deduplication"]["raw_records"] == 120
    assert audit["context_deduplication"]["sample_units"] == 119
    duplicate = next(
        d for d in audit["context_deduplication"]["decisions"]
        if d["sample_unit_id"] == "SAMEA123"
    )
    assert duplicate["source_genome_ids"] == ["source0", "source1"]
    assert audit["pools"][0]["requested"] == 120
    assert audit["pools"][0]["level_counts"] == {str(level): 120 for level in range(5, 11)}
    assert audit["context_assemblies"] == "not_requested"
    assert all(not r.get("assembly_reference") for r in ds.samples[1:])


def test_context_resume_reuses_completed_scope_and_frozen_catalogues(tmp_path, monkeypatch):
    import chronoclade.profile_context_provider as provider

    identity = {"value": "source-v1"}
    monkeypatch.setattr(provider, "execution_provenance", lambda *a: {"identity_sha256": identity["value"]})
    q1 = record("q1")
    q2 = record("q2")
    q2["mlst_st"] = "40"
    source = tmp_path / "profiles.json"
    source.write_text(json.dumps([q1, q2]))
    cats = tmp_path / "catalogues.json"
    cats.write_text(json.dumps([asdict(LocusCatalogue("test", "v1", ("a",), "fixture"))]))
    prepared = run_prepare(source, tmp_path / "prepared", catalogues=cats)

    class Client:
        api_key = "fake"
        base_url = "https://example.test"

        def __init__(self):
            self.freezes = []

        def supported_organisms(self):
            return [{"fullName": SPECIES, "organismId": "573"}]

        def freeze_catalogue(self, path, *, organism_id, st, resume=False):
            self.freezes.append(st)
            row = record("PW_source" + st)
            row.update(source_genome_id="source" + st, origin="context", aliases=["source" + st])
            envelope = {
                "rows": [row],
                "provenance": {
                    "complete": True,
                    "expected_count": 1,
                    "query": {"organismId": organism_id, "mlst": [st]},
                    "organism_capabilities": {"fullName": SPECIES, "organismId": organism_id},
                    "base_url": self.base_url,
                },
            }
            path.write_text(json.dumps(envelope))
            return envelope

    client = Client()
    actual_loader = provider.load_catalogue
    monkeypatch.setattr(
        provider,
        "load_catalogue",
        lambda path: json.loads(path.read_text()),
    )
    export_calls = []
    fail_once = {"value": True}

    def exports(rows, _client, _output, _provenance, download_names=None, resume=False):
        ids = [r["source_genome_id"] for r in rows]
        export_calls.extend(ids)
        if download_names is None and "source40" in ids and fail_once["value"]:
            fail_once["value"] = False
            raise RuntimeError("simulated interrupted profile export")
        if download_names is None:
            return [dict(r, cgmlst_profile={"a": "1"}) for r in rows]
        return rows

    monkeypatch.setattr(provider, "_grouped_analysis_exports", exports)
    output = tmp_path / "resumed"
    with pytest.raises(RuntimeError, match="interrupted"):
        discover_profile_context(prepared.dataset_manifest, output, client=client)

    retrieval = tmp_path / "resumed_retrieval"
    assert sorted(client.freezes) == ["39", "40"]
    assert (retrieval / "public_0.json").is_file()
    assert (retrieval / "public_1.json").is_file()
    assert (retrieval / "matched_0.json").is_file()
    first_scope_calls = export_calls.count("source39")

    with pytest.raises(DatasetError, match="inputs changed"):
        discover_profile_context(
            prepared.dataset_manifest, output, client=client, resume=True, lin_level=6
        )
    assert len(client.freezes) == 2

    identity["value"] = "source-v2"
    with pytest.raises(DatasetError, match="inputs changed"):
        discover_profile_context(prepared.dataset_manifest, output, client=client, resume=True)
    identity["value"] = "source-v1"

    # This mocked loader stands in for load_catalogue's full snapshot validation.
    monkeypatch.setattr(provider, "load_catalogue", lambda path: json.loads(path.read_text()))
    result = discover_profile_context(prepared.dataset_manifest, output, client=client, resume=True)
    assert len(client.freezes) == 2
    assert export_calls.count("source39") == first_scope_calls == 2
    assert len(export_calls) == 6
    assert result.dataset_manifest.is_file()
    call_count = len(export_calls)
    reused = discover_profile_context(
        prepared.dataset_manifest, output, client=client, resume=True
    )
    assert reused == result
    assert len(export_calls) == call_count and len(client.freezes) == 2
    with pytest.raises(DatasetError, match="already exists"):
        discover_profile_context(prepared.dataset_manifest, output, client=client)
    identity["value"] = "source-v2"
    with pytest.raises(DatasetError, match="does not match resume inputs"):
        discover_profile_context(prepared.dataset_manifest, output, client=client, resume=True)
    assert len(export_calls) == call_count and len(client.freezes) == 2
    monkeypatch.setattr(provider, "load_catalogue", actual_loader)


def test_context_resume_reuses_source_copy_after_packaging_failure(tmp_path, monkeypatch):
    import chronoclade.profile_context_provider as provider
    import chronoclade.prepare_stage as prepare_stage

    source = tmp_path / "profiles.json"
    source.write_text(json.dumps([record()]))
    cats = tmp_path / "catalogues.json"
    cats.write_text(json.dumps([asdict(LocusCatalogue("test", "v1", ("a",), "fixture"))]))
    prepared = run_prepare(source, tmp_path / "prepared", catalogues=cats)

    class Client:
        api_key = "fake"
        base_url = "https://example.test"

        def __init__(self):
            self.freezes = 0

        def supported_organisms(self):
            return [{"fullName": SPECIES, "organismId": "573"}]

        def freeze_catalogue(self, path, *, organism_id, st, resume=False):
            self.freezes += 1
            row = record("PW_source")
            row.update(source_genome_id="source", origin="context", aliases=["source"])
            envelope = {
                "rows": [row],
                "provenance": {
                    "complete": True,
                    "expected_count": 1,
                    "query": {"organismId": organism_id, "mlst": [st]},
                    "organism_capabilities": {"fullName": SPECIES, "organismId": organism_id},
                    "base_url": self.base_url,
                },
            }
            path.write_text(json.dumps(envelope))
            return envelope

    client = Client()
    export_calls = []

    def exports(rows, _client, _output, _provenance, download_names=None, resume=False):
        export_calls.extend(r["source_genome_id"] for r in rows)
        if download_names is None:
            return [dict(row, cgmlst_profile={"a": "1"}) for row in rows]
        return rows

    monkeypatch.setattr(provider, "_grouped_analysis_exports", exports)
    original_report = prepare_stage._report
    fail_once = {"value": True}

    def report_with_interruption(audit):
        if fail_once["value"]:
            fail_once["value"] = False
            raise RuntimeError("simulated packaging failure")
        return original_report(audit)

    monkeypatch.setattr(prepare_stage, "_report", report_with_interruption)
    output = tmp_path / "packaging"
    with pytest.raises(RuntimeError, match="packaging failure"):
        discover_profile_context(prepared.dataset_manifest, output, client=client)
    source_copy = tmp_path / "packaging_retrieval" / "source_dataset"
    assert (source_copy / "dataset.json").is_file()

    result = discover_profile_context(
        prepared.dataset_manifest, output, client=client, resume=True
    )
    assert result.dataset_manifest.is_file()
    assert client.freezes == 1
    assert export_calls == ["source", "source"]


def test_context_ecoli_requires_explicit_hiercc_level(tmp_path):
    source = tmp_path / "p.json"
    source.write_text(json.dumps([{"sample_id": "q", "species": "Escherichia coli"}]))
    cats = tmp_path / "c.json"
    cats.write_text(json.dumps([asdict(LocusCatalogue("test", "v1", ("a",), "fixture"))]))
    prepared = run_prepare(source, tmp_path / "prepared", catalogues=cats)
    with pytest.raises(DatasetError, match="explicit"):
        discover_profile_context(prepared.dataset_manifest, tmp_path / "context")


def test_existing_profile_runs_assigner_without_assembly(tmp_path, monkeypatch):
    from chronoclade import query_typing
    from chronoclade.prepare_provider import assign_existing_profiles

    row = record()
    row.pop("cglin_raw")
    row["cglin_status"] = "unassigned"
    calls = []
    monkeypatch.setattr(
        query_typing,
        "_provenance",
        lambda tool, fields: {"databases": {"index_dir": {"sha256": "a" * 64}}},
    )

    def run(command, data, destination):
        calls.append((command, json.loads(data)))
        return {"LINcode": [0, 0, 1, 2, 3, 4, 5, 6, 7, 8]}

    monkeypatch.setattr(query_typing, "_run", run)
    config = {
        "organisms": {
            SPECIES: {
                "cgmlst": {"scheme": "test", "scheme_version": "v1", "index_dir": tmp_path},
                "assignment": {
                    "kind": "plincer",
                    "scheme": "scgMLST629_S",
                    "scheme_version": "lin-v1",
                    "canonical_loci": ["a", "b"],
                    "command": ["plincer"],
                    "scheme_file": tmp_path / "scheme.toml",
                    "profiles_file": tmp_path / "profiles.json",
                    "alleles_db": tmp_path / "alleles.sqlite",
                },
            }
        }
    }
    assert assign_existing_profiles([row], config, tmp_path / "assignment") == {"query"}
    assert calls[0][0][1:3] == ["classify", "-"]
    assert calls[0][1]["code"] == "1_"
    assert row["cglin_status"] == "complete"
    assert "assembly" not in row
    row.pop("cglin_raw")
    row["cgmlst_scheme_version"] = "different"
    assert assign_existing_profiles([row], config, tmp_path / "other") == set()
    assert len(calls) == 1


def test_lineage_and_crosswalk_cannot_overwrite_authoritative_identity(tmp_path):
    from chronoclade.profile_context_provider import dataset_records

    source = tmp_path / "records.json"
    row = record()
    source.write_text(
        json.dumps(
            {
                "records": [row],
                "lineages": [
                    {
                        "sample_id": "query",
                        "kind": "cglin",
                        "scheme_id": "scgMLST629_S",
                        "scheme_version": "lin-v1",
                        "database_version": "canonical",
                        "assignment": row["cglin_raw"],
                        "resolution_status": "complete",
                        "assignment_method": "fixture",
                        "evidence": {
                            "sample_id": "intruder",
                            "role": "context",
                            "species": "intruder",
                            "country": "intruder",
                            "cgmlst_profile": {"a": "99"},
                        },
                    }
                ],
                "crosswalk": [
                    {
                        "sample_id": "query",
                        "role": "context",
                        "country": "intruder",
                        "source_genome_id": "source",
                        "numeric_source_id": 83,
                        "source_checksum": "a" * 40,
                        "source_length": 12345,
                        "qc_pass": True,
                    }
                ],
            }
        )
    )
    cats = tmp_path / "catalogues.json"
    cats.write_text(json.dumps([asdict(LocusCatalogue("test", "v1", ("a",), "fixture"))]))
    result = run_prepare(source, tmp_path / "prepared", catalogues=cats)
    adapted = dataset_records(load_dataset(result.dataset_manifest))[0]
    assert adapted["sample_id"] == "query" and adapted["role"] == "input"
    assert adapted["species"] == SPECIES and adapted["country"] is None
    assert adapted["cgmlst_profile"] == {"a": "1"}
    assert adapted["cglin_database_version"] == "canonical"
    assert adapted["source_genome_id"] == "source"
    assert adapted["numeric_source_id"] == 83
    assert adapted["source_checksum"] == "a" * 40
    assert adapted["source_length"] == 12345 and adapted["qc_pass"] is True


def test_finer_context_retrieves_only_compatible_prefix_before_profiles(tmp_path, monkeypatch):
    import chronoclade.profile_context_provider as provider

    source = tmp_path / "profiles.json"
    source.write_text(json.dumps([record()]))
    cats = tmp_path / "catalogues.json"
    cats.write_text(json.dumps([asdict(LocusCatalogue("test", "v1", ("a",), "fixture"))]))
    prepared = run_prepare(source, tmp_path / "prepared", catalogues=cats)
    matching = record("PW_match")
    matching.update(source_genome_id="match", role="context", cgmlst_profile={})
    excluded = dict(matching, sample_id="PW_outside", source_genome_id="outside",
                    cglin_raw="0,0,1,2,3,4,5,99,7,8")
    incompatible = dict(matching, sample_id="PW_version", source_genome_id="version",
                        cglin_scheme_version="another-scheme-version")

    class Client:
        def supported_organisms(self):
            return [{"fullName": SPECIES, "organismId": "573"}]

        def freeze_catalogue(self, path, **kwargs):
            envelope = {"rows": [matching, excluded, incompatible],
                        "provenance": {"complete": True, "expected_count": 3}}
            path.write_text(json.dumps(envelope))
            return envelope

    requested = []

    def exports(rows, client, output, provenance, download_names=None, resume=False):
        if download_names is not None:
            return rows
        requested.extend(row["source_genome_id"] for row in rows)
        return [dict(row, cgmlst_profile={"a": "1"}) for row in rows]

    monkeypatch.setattr(provider, "_grouped_analysis_exports", exports)
    result = discover_profile_context(prepared.dataset_manifest, tmp_path / "context",
                                      lin_level=8, client=Client())
    assert requested == ["match"]
    dataset = load_dataset(result.dataset_manifest)
    assert len(dataset.samples) == 2
    audit = json.loads(result.audit_path.read_text())["provider_context"]
    pool = audit["pools"][0]
    assert pool["discovered"] == 3 and pool["requested"] == pool["retrieved"] == 1
    assert pool["level_counts"]["7"] == 2 and pool["level_counts"]["8"] == 1
    assert audit["policy"]["lin_level"] == 8
    assert audit["context_assemblies"] == "not_requested"


@pytest.mark.parametrize("failure_kind", ["cglin", "cgmlst", None])
def test_public_context_distinguishes_incompatible_scope_from_empty_pool(
    tmp_path, monkeypatch, failure_kind
):
    import chronoclade.profile_context_provider as provider

    source = tmp_path / "profiles.json"
    source.write_text(json.dumps([record()]))
    cats = tmp_path / "catalogues.json"
    cats.write_text(json.dumps([asdict(LocusCatalogue("test", "v1", ("a",), "fixture"))]))
    prepared = run_prepare(source, tmp_path / "prepared", catalogues=cats)
    candidate = record("PW_public")
    candidate.update(source_genome_id="public", role="context", origin="context")
    if failure_kind:
        candidate[failure_kind + "_scheme_version"] = "other-reference"
    else:
        # Compatible scope, genuinely outside the requested prefix.
        candidate["cglin_raw"] = "0,0,1,2,99,4,5,6,7,8"

    class Client:
        base_url = "https://example.test"

        def supported_organisms(self):
            return [{"fullName": SPECIES, "organismId": "573"}]

        def freeze_catalogue(self, path, **kwargs):
            envelope = {"rows": [candidate], "provenance": {
                "complete": True, "expected_count": 1,
            }}
            path.write_text(json.dumps(envelope))
            return envelope

    monkeypatch.setattr(provider, "_grouped_analysis_exports", lambda rows, *args, **kwargs: rows)
    output = tmp_path / "context"
    if failure_kind:
        with pytest.raises(DatasetError, match="does not mean there are no public relatives"):
            discover_profile_context(prepared.dataset_manifest, output, client=Client())
        assert not output.exists()
        audit = json.loads((tmp_path / "context_retrieval/context_compatibility.json").read_text())
        assert audit["status"] == "failed_incompatible_typing_scope"
        assert audit["kind"] == failure_kind
        assert audit["rejected_candidates"][0]["source_genome_id"] == "public"
        assert audit["context_assemblies"] == "not_requested"
    else:
        result = discover_profile_context(prepared.dataset_manifest, output, client=Client())
        audit = json.loads(result.audit_path.read_text())["provider_context"]
        assert audit["usable"] == 0
        assert audit["pools"][0]["requested"] == 0
        assert audit["pools"][0]["incompatible_numeric_lineage_candidates"] == []
        assert not (tmp_path / "context_retrieval/context_compatibility.json").exists()


def test_native_context_reassigns_all_st_profiles_before_finer_lin_scope(tmp_path, monkeypatch):
    import chronoclade.profile_context_provider as provider
    import chronoclade.pasteur_profile_normalization as normalization

    source, cats = tmp_path / "profiles.json", tmp_path / "catalogues.json"
    source.write_text(json.dumps([record()]))
    cats.write_text(json.dumps([asdict(LocusCatalogue("test", "v1", ("a",), "fixture"))]))
    prepared = run_prepare(source, tmp_path / "prepared", catalogues=cats)
    reference = tmp_path / "native.json"
    reference.write_text("{}")
    monkeypatch.setattr(normalization, "reference_descriptor", lambda path: {"config_sha256": "a" * 64})
    requested = []

    class Normalizer:
        def __init__(self, *_):
            self.audits = []

        def normalize(self, rows):
            # A profile outside the OLD reported LIN8 prefix acquires the
            # query's actual CURRENT assignment after real profile analysis.
            return [dict(r, cglin_raw="0,0,1,2,3,4,5,6,7,8",
                         cglin_scheme_version="current-reference") for r in rows]

    monkeypatch.setattr(normalization, "PasteurProfileNormalizer", Normalizer)
    candidate = record("PW_public")
    candidate.update(source_genome_id="public", origin="context", role="context",
                     cgmlst_profile={}, cglin_raw="0,0,1,2,3,4,5,99,7,8",
                     cglin_scheme_version="old-job")

    class Client:
        base_url = "https://example.test"

        def supported_organisms(self):
            return [{"fullName": SPECIES, "organismId": "573"}]

        def freeze_catalogue(self, path, **kwargs):
            value = {"rows": [candidate], "provenance": {"complete": True, "expected_count": 1}}
            path.write_text(json.dumps(value))
            return value

    def exports(rows, *args, download_names=None, **kwargs):
        if download_names is not None:
            return rows
        requested.extend(r["source_genome_id"] for r in rows)
        return [dict(r, cgmlst_profile={"a": "1"}) for r in rows]

    monkeypatch.setattr(provider, "_grouped_analysis_exports", exports)
    result = discover_profile_context(prepared.dataset_manifest, tmp_path / "context",
                                      client=Client(), lin_level=8, typing_config=reference)
    assert requested == ["public"]
    audit = json.loads(result.audit_path.read_text())["provider_context"]
    assert audit["usable"] == 1
    assert audit["pools"][0]["current_reference_reassignment"]["prefilter"] == "all_ST_profiles_before_current_LIN_scope"
    assert audit["pools"][0]["level_counts"]["8"] == 1


def test_real_batch_checkpoints_survive_interrupted_normalizer_and_reference_guard(
    tmp_path, monkeypatch
):
    import csv
    import hashlib
    import io
    import chronoclade.profile_inputs as inputs
    import chronoclade.profile_context_provider as provider
    import chronoclade.pasteur_profile_normalization as normalization

    source, cats = tmp_path / "input.json", tmp_path / "catalogues.json"
    source.write_text(json.dumps([record()]))
    cats.write_text(json.dumps([asdict(LocusCatalogue("test", "v1", ("a",), "fixture"))]))
    prepared = run_prepare(source, tmp_path / "prepared", catalogues=cats)
    reference = tmp_path / "native.json"
    reference.write_text("reference-A")
    monkeypatch.setattr(
        normalization,
        "reference_descriptor",
        lambda path: {"config_sha256": hashlib.sha256(path.read_bytes()).hexdigest()},
    )
    normalization_counts = []
    interrupt = [True]

    class Normalizer:
        def __init__(self, *args):
            self.audits = []

        def normalize(self, values):
            if any(r["role"] == "context" for r in values):
                normalization_counts.append(len(values))
                if interrupt[0]:
                    interrupt[0] = False
                    raise RuntimeError("classifier interrupted")
            return [
                dict(
                    r,
                    cglin_raw="0,0,1,2,3,4,5,6,7,8",
                    cglin_scheme_version="current",
                    cgmlst_scheme="test",
                    cgmlst_scheme_version="v1",
                    cgmlst_profile={"a": "1"},
                )
                for r in values
            ]

    monkeypatch.setattr(normalization, "PasteurProfileNormalizer", Normalizer)
    public = []
    for i in range(201):
        row = record(f"PW_id-{i}")
        row.update(
            source_genome_id=f"id-{i}",
            numeric_source_id=i + 1,
            organism_id="573",
            role="context",
            origin="context",
            cgmlst_profile={},
        )
        public.append(row)
    events = []

    class Client:
        base_url = "https://pathogen.watch"
        api_key = "not-recorded"

        def supported_organisms(self):
            return [{"fullName": SPECIES, "organismId": "573"}]

        def freeze_catalogue(self, path, **kwargs):
            events.append("catalogue")
            envelope = {
                "rows": public,
                "provenance": {
                    "complete": True,
                    "expected_count": len(public),
                    "query": {"organismId": "573", "mlst": ["39"]},
                    "organism_capabilities": {"fullName": SPECIES, "organismId": "573"},
                    "base_url": self.base_url,
                    "species": SPECIES,
                    "retrieved_at": "2026-10-18T00:00:00Z",
                },
            }
            path.write_text(json.dumps(envelope))
            return envelope

        def request_json(self, method, path, **kwargs):
            if path == "/api/search/genomes":
                events.append("preview")
                return {"meta": {"count": len(public)}, "genomes": []}
            events.append("capability")
            return [
                {
                    "ids": kwargs["body"]["ids"],
                    "organismId": "573",
                    "downloads": [{"name": "cgmlst", "job": "cgmlst-573-2"}],
                }
            ]

    def download(route, **kwargs):
        ids = kwargs["body"]["ids"].split(",")
        events.append(tuple(ids))
        stream = io.StringIO()
        writer = csv.writer(stream)
        writer.writerow(["Genome ID", "Gene", "Allele ID"])
        for number in ids:
            writer.writerow([f"id-{int(number) - 1}", "g1", "1"])
        return stream.getvalue().encode()

    monkeypatch.setattr(inputs, "request_download", download)
    monkeypatch.setattr(provider, "load_catalogue", lambda path: json.loads(path.read_text()))
    monkeypatch.setattr("chronoclade.context_preflight.load_catalogue", lambda path: json.loads(path.read_text()))
    output = tmp_path / "context"
    client = Client()
    with pytest.raises(RuntimeError, match="classifier interrupted"):
        discover_profile_context(
            prepared.dataset_manifest, output, client=client, typing_config=reference, lin_level=8
        )
    checkpoint = (
        tmp_path
        / "context_retrieval/scope_0/profiles/organism_0/batch_cgmlst_cgmlst-573-2_0/receipt.json"
    )
    assert checkpoint.is_file()
    before = list(events)
    reference.write_text("reference-B")
    with pytest.raises(DatasetError, match="resume inputs changed"):
        discover_profile_context(
            prepared.dataset_manifest,
            output,
            client=client,
            typing_config=reference,
            lin_level=8,
            resume=True,
        )
    assert events == before
    reference.write_text("reference-A")
    result = discover_profile_context(
        prepared.dataset_manifest,
        output,
        client=client,
        typing_config=reference,
        lin_level=8,
        resume=True,
    )
    assert result.dataset_manifest.is_file()
    assert events == before
    assert normalization_counts == [201, 201]


def test_live_local_prepare_resume_reuses_completed_checkpoint_and_publishes_portably(
    tmp_path, monkeypatch
):
    from chronoclade import pasteur_profile_normalization, prepare_provider, query_typing
    from chronoclade.artifacts import file_sha256
    from chronoclade.datasets import load_dataset
    from chronoclade.prepare_stage import run_prepare

    assemblies = []
    for ident in ("first", "second"):
        assembly = tmp_path / f"{ident}.fa"
        assembly.write_text(f">{ident}\nACGT\n")
        assemblies.append(assembly)
    metadata = tmp_path / "queries.csv"
    metadata.write_text(
        "sample_id,assembly,species,lineage,mlst_st,country,collection_date\n"
        "first,first.fa,Klebsiella pneumoniae,ST1,1,UK,2020\n"
        "second,second.fa,Klebsiella pneumoniae,ST1,1,UK,2020\n"
    )
    config = tmp_path / "typing.json"
    config.write_text("{}\n")
    monkeypatch.setattr(
        query_typing,
        "load_typing_config",
        lambda *_a, **_k: {
            "organisms": {
                "Klebsiella pneumoniae": {
                    "cgmlst": {"index_loci": ["locus-a"], "scheme": "klebsiella_1"},
                    "assignment": {"kind": "plincer"},
                }
            }
        },
    )
    monkeypatch.setattr(
        query_typing,
        "_provenance",
        lambda *_a, **_k: {"databases": {"index_dir": {"sha256": "a" * 64}}},
    )
    monkeypatch.setattr(prepare_provider, "assign_existing_profiles", lambda *_a, **_k: set())
    monkeypatch.setattr(
        pasteur_profile_normalization,
        "reference_descriptor",
        lambda *_a, **_k: {"fixture": "native-reference-binding-not-tested"},
    )

    calls = []
    fail_second = True

    def fake_native(samples, _config, output, *, existing_assignments=None,
                    profiles_only=False, progress_index=1, progress_total=None, **_private):
        nonlocal fail_second
        sample = samples[0]
        calls.append(sample.sample_id)
        output.mkdir(parents=True, exist_ok=True)
        (output / f"{sample.sample_id}.native.json").write_text('{"raw":"kept"}\n')
        if sample.sample_id == "second" and fail_second:
            fail_second = False
            raise RuntimeError("synthetic native failure")
        identity = query_typing._typing_recovery_identity([sample], _config, False)
        row = identity["samples"][0]
        return [{
            "sample_id": sample.sample_id,
            "species": sample.species,
            "assembly_sha256": file_sha256(sample.assembly),
            "typing_mode": "full",
            "typing_status": "typed",
            "typing_provenance": row["typing_provenance"],
            "cgmlst_scheme": "klebsiella_1",
            "cgmlst_scheme_version": "fixture-v1",
            "cgmlst_profile": {"locus-a": "1"},
            "cgmlst_loci": ["locus-a"],
            "cgmlst_locus_universe_complete": True,
        }]

    monkeypatch.setattr(query_typing, "_type_query_assemblies_once", fake_native)
    target = tmp_path / "prepared"
    with pytest.raises(DatasetError, match="Resume with:.*--resume") as error:
        run_prepare(
            metadata,
            target,
            input_kind="assemblies",
            typing_config=config,
            enrich_metadata=False,
        )
    assert "--typing-config" in str(error.value)
    assert "chronoclade prepare" in str(error.value)
    work = target.with_name("prepared.work")
    assert work.is_dir()
    first_attempt = next((work / "attempts").iterdir())
    assert list((first_attempt / "native_typing_attempts").rglob("first.native.json"))
    assert list((first_attempt / "native_typing_attempts").rglob("second.native.json"))
    assert calls == ["first", "second"]

    original_csv = metadata.read_text()
    metadata.write_text(original_csv + "\n")
    with pytest.raises(DatasetError, match="identity-mismatched"):
        run_prepare(metadata, target, input_kind="assemblies", typing_config=config,
                    enrich_metadata=False, resume=True)
    assert calls == ["first", "second"]
    metadata.write_text(original_csv)

    original_assembly = assemblies[1].read_text()
    assemblies[1].write_text(">second\nTGCA\n")
    with pytest.raises(DatasetError, match="identity differs"):
        run_prepare(metadata, target, input_kind="assemblies", typing_config=config,
                    enrich_metadata=False, resume=True)
    assert calls == ["first", "second"]
    assemblies[1].write_text(original_assembly)

    result = run_prepare(
        metadata,
        target,
        input_kind="assemblies",
        typing_config=config,
        enrich_metadata=False,
        resume=True,
    )
    assert calls == ["first", "second", "second"]
    assert not work.exists()
    preserved = target / "sources" / "provider" / "prior_attempts" / first_attempt.name
    assert list((preserved / "native_typing_attempts").rglob("first.native.json"))
    dataset = load_dataset(result.dataset_manifest)
    assert dataset.sample_ids == ("first", "second")
    for row in dataset.samples:
        reference = row["assembly_reference"]
        assert (target / reference).is_file()
        assert not Path(reference).is_absolute()
    assert dataset.parameters["prepare"]["assemblies"]
    assert json.loads(result.stage_manifest.read_text())["status"] == "complete"
    with pytest.raises(DatasetError, match="already exists"):
        run_prepare(metadata, target, input_kind="assemblies", typing_config=config,
                    enrich_metadata=False, resume=True)


def test_live_local_prepare_resume_before_native_cache_initialization(tmp_path, monkeypatch):
    from chronoclade import pasteur_profile_normalization, prepare_provider, query_typing
    from chronoclade.artifacts import file_sha256
    from chronoclade.datasets import load_dataset
    from chronoclade.prepare_stage import run_prepare

    assembly = tmp_path / "query.fa"
    assembly.write_text(">query\nACGT\n")
    metadata = tmp_path / "queries.csv"
    metadata.write_text(
        "sample_id,assembly,species,lineage,mlst_st,country,collection_date\n"
        "query,query.fa,Klebsiella pneumoniae,ST1,1,UK,2020\n"
    )
    config = tmp_path / "typing.json"
    config.write_text("{}\n")
    monkeypatch.setattr(
        query_typing,
        "load_typing_config",
        lambda *_a, **_k: {"organisms": {"Klebsiella pneumoniae": {
            "cgmlst": {"index_loci": ["locus-a"], "scheme": "klebsiella_1"},
            "assignment": {"kind": "plincer"},
        }}},
    )
    monkeypatch.setattr(
        query_typing, "_provenance",
        lambda *_a, **_k: {"databases": {"index_dir": {"sha256": "a" * 64}}},
    )
    monkeypatch.setattr(prepare_provider, "assign_existing_profiles", lambda *_a, **_k: set())
    monkeypatch.setattr(
        pasteur_profile_normalization, "reference_descriptor",
        lambda *_a, **_k: {"fixture": "native-reference-binding-not-tested"},
    )

    original_materialise = prepare_provider.materialise_assemblies
    materialise_calls = 0

    def fail_before_typing_once(*args, **kwargs):
        nonlocal materialise_calls
        materialise_calls += 1
        if materialise_calls == 1:
            raise RuntimeError("synthetic interruption before native cache initialization")
        return original_materialise(*args, **kwargs)

    monkeypatch.setattr(prepare_provider, "materialise_assemblies", fail_before_typing_once)
    native_calls = []

    def fake_native(samples, typing_config, output, **_kwargs):
        sample = samples[0]
        native_calls.append(sample.sample_id)
        identity = query_typing._typing_recovery_identity([sample], typing_config, False)
        return [{
            "sample_id": sample.sample_id,
            "species": sample.species,
            "assembly_sha256": file_sha256(sample.assembly),
            "typing_mode": "full",
            "typing_status": "typed",
            "typing_provenance": identity["samples"][0]["typing_provenance"],
            "cgmlst_scheme": "klebsiella_1",
            "cgmlst_scheme_version": "fixture-v1",
            "cgmlst_profile": {"locus-a": "1"},
            "cgmlst_loci": ["locus-a"],
            "cgmlst_locus_universe_complete": True,
        }]

    monkeypatch.setattr(query_typing, "_type_query_assemblies_once", fake_native)
    target = tmp_path / "prepared"
    work = target.with_name("prepared.work")
    with pytest.raises(DatasetError, match="before native cache initialization"):
        run_prepare(metadata, target, input_kind="assemblies", typing_config=config,
                    enrich_metadata=False)
    assert not (work / "native_typing_cache").exists()
    assert native_calls == []

    result = run_prepare(metadata, target, input_kind="assemblies", typing_config=config,
                         enrich_metadata=False, resume=True)
    assert native_calls == ["query"]
    assert not work.exists()
    dataset = load_dataset(result.dataset_manifest)
    assert dataset.sample_ids == ("query",)
    reference = dataset.samples[0]["assembly_reference"]
    assert not Path(reference).is_absolute()
    assert (target / reference).is_file()
    assert json.loads(result.stage_manifest.read_text())["status"] == "complete"
