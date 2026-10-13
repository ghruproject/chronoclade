"""Exact provider joins, complete profile pools and readiness are independently audited."""

import json
from dataclasses import asdict

import pytest

from chronoclade.datasets import DatasetError, LocusCatalogue, load_dataset
from chronoclade.pathogenwatch import PathogenwatchClient
from chronoclade.prepare_provider import (
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
    assert list(
        (result.dataset_manifest.parent / "sources/provider").rglob("collection_snapshot.json")
    )


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

    def exports(rows, client, output, provenance, download_names=None):
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
    assert audit["pools"][0]["level_counts"] == {"5": 120, "6": 120, "7": 120}
    assert audit["context_assemblies"] == "not_requested"
    assert all(not r.get("assembly_reference") for r in ds.samples[1:])


def test_context_resume_reuses_completed_scope_and_frozen_catalogues(tmp_path, monkeypatch):
    import chronoclade.profile_context_provider as provider

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

        def freeze_catalogue(self, path, *, organism_id, st):
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

    def exports(rows, _client, _output, _provenance, download_names=None):
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

        def freeze_catalogue(self, path, *, organism_id, st):
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

    def exports(rows, _client, _output, _provenance, download_names=None):
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
