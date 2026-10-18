"""Raw replay tests use only synthetic official-format responses, never credentials."""

from dataclasses import asdict
import json

import pytest
from typer.testing import CliRunner

from chronoclade.artifacts import write_json
from chronoclade.datasets import DatasetError, LocusCatalogue, load_dataset
from chronoclade.pathogenwatch import content_hash, normalize_record, load_catalogue
from chronoclade.provider_sources import ProviderSource, LEGACY_RECEIPT



def capture(root, count=3):
    root.mkdir()
    supported = [{"organismId": "573", "fullName": "Klebsiella pneumoniae", "typing": ["MLST"]}]
    searches, details = [], {}
    for index in range(count):
        ident = f"g{index:03d}"
        searches.append(
            dict(
                uuid=ident,
                projectAccess="PUBLIC",
                name=ident,
                organism="Klebsiella pneumoniae",
                location=None,
            )
        )
        details[ident] = dict(
            uuid=ident,
            id=index + 1,
            organismId="573",
            mlst="147",
            name=ident,
            qc=True,
            metadata=[["Country", "GB"], ["Date", "2019"]],
            sampleAccession=f"SAMN100{index}",
            runAccession=f"SRR100{index}",
        )
    rows = [normalize_record(s, details[s["uuid"]], organism_id="573", st="147") for s in searches]
    pages = [{"genomes": searches, "meta": {"count": count}}]
    query = {"organismId": "573", "mlst": ["147"], "qc": [True, False]}
    provenance = dict(
        source="pathogenwatch",
        base_url="https://pathogen.watch",
        query=query,
        query_sha256=content_hash(query),
        retrieved_at="2026-10-17T00:00:00Z",
        species="Klebsiella pneumoniae",
        organism_capabilities=supported[0],
        public_only=True,
        complete=True,
        expected_count=count,
        raw_record_count=count,
        raw_search_sha256=content_hash(pages),
        raw_details_sha256=content_hash(details),
        rows_sha256=content_hash(rows),
    )
    payload = dict(
        schema_version=1,
        provenance=provenance,
        rows=rows,
        raw_details=details,
        raw_search_pages=pages,
        raw_supported_organisms=supported,
    )
    payload["snapshot_sha256"] = content_hash(payload)
    write_json(root / "public_0.json", payload)
    directory = root / "scope_0_exact/profiles"
    directory.mkdir(parents=True)
    caps = dict(
        groups=[
            dict(
                ids=list(range(1, count + 1)),
                organismId="573",
                downloads=[dict(name="cgmlst", job="cgmlst-573-2")],
            )
        ],
        requested_numeric_ids=list(range(1, count + 1)),
    )
    caps["sha256"] = content_hash(caps)
    write_json(directory / "export_capabilities.json", caps)
    (directory / "organism_0").mkdir()
    for start in range(0, count, 100):
        (directory / "organism_0" / f"cgmlst_cgmlst-573-2_{start}.bin").write_text(
            "Genome ID,Gene,Allele ID\n"
            + "".join(f"g{i:03d},A,{i + 1}\n" for i in range(start, min(start + 100, count)))
        )
    return root


def rehash(path, key):
    d = json.loads(path.read_text())
    d[key] = content_hash({k: v for k, v in d.items() if k != key})
    write_json(path, d)


def test_legacy_complete_batches_and_portable_roundtrip(tmp_path):
    root = capture(tmp_path / "source", 101)
    source = ProviderSource(root)
    cap, rows, assignments = source.for_scope("Klebsiella pneumoniae", "147")
    assert len(rows) == len(assignments) == 101
    receipts = source.scopes[0]["scope"]["receipts"]
    assert len(receipts) == 2 and len(receipts[-1]["numeric_request_ids"]) == 1
    assert all(r["request_evidence"] == LEGACY_RECEIPT for r in receipts)
    snapshot = source.snapshot(tmp_path / "portable")
    restored = ProviderSource(snapshot)
    assert restored.for_scope("Klebsiella pneumoniae", "147")[2] == assignments
    assert restored.base_url == "https://pathogen.watch"
    second_snapshot = restored.snapshot(tmp_path / "portable-again")
    assert snapshot.read_bytes() == second_snapshot.read_bytes()
    with pytest.raises(DatasetError, match="exact requested"):
        restored.for_scope("Klebsiella pneumoniae", "39")


@pytest.mark.parametrize(
    "damage",
    [
        "missing",
        "unrequested",
        "retry",
        "duplicate",
        "subset",
        "organism",
        "job",
        "cap_hash",
        "catalogue_hash",
        "escape",
    ],
)
def test_rejects_unsafe_or_incomplete_raw_scope(tmp_path, damage):
    root = capture(tmp_path / "source")
    directory = root / "scope_0_exact/profiles"
    response = directory / "organism_0/cgmlst_cgmlst-573-2_0.bin"
    caps = directory / "export_capabilities.json"
    if damage == "missing":
        response.unlink()
    elif damage == "unrequested":
        response.write_text(response.read_text() + "alien,A,1\n")
    elif damage == "retry":
        response.with_name("cgmlst_cgmlst-573-2_0_retry1.bin").write_bytes(response.read_bytes())
    elif damage == "escape":
        outside = tmp_path / "outside.bin"
        outside.write_bytes(response.read_bytes())
        response.unlink()
        response.symlink_to(outside)
    elif damage.startswith("catalogue"):
        payload = json.loads((root / "public_0.json").read_text())
        payload["provenance"]["expected_count"] = 1
        write_json(root / "public_0.json", payload)
    else:
        d = json.loads(caps.read_text())
        if damage == "duplicate":
            d["groups"][0]["ids"] = [1, 1, 3]
        elif damage == "subset":
            d["requested_numeric_ids"] = [1, 2]
            d["groups"][0]["ids"] = [1, 2]
        elif damage == "organism":
            d["groups"][0]["organismId"] = "562"
        elif damage == "job":
            d["groups"][0]["downloads"] = []
        else:
            d["requested_numeric_ids"] = [999]
        write_json(caps, d)
        if damage != "cap_hash":
            rehash(caps, "sha256")
    with pytest.raises(DatasetError):
        ProviderSource(root)


def test_raw_replay_recomputes_metadata_without_weakening_normalized_loader(tmp_path, monkeypatch):
    import chronoclade.pathogenwatch as pw

    root = capture(tmp_path / "source")
    original = pw.normalize_record
    monkeypatch.setattr(
        pw,
        "normalize_record",
        lambda *a, **k: dict(original(*a, **k), country="Changed current method"),
    )
    with pytest.raises(RuntimeError, match="normalised rows"):
        load_catalogue(root / "public_0.json")
    source = ProviderSource(root)
    assert (
        source.for_scope("Klebsiella pneumoniae", "147")[1][0]["country"]
        == "Changed current method"
    )


def test_missing_profiles_explicit_and_raw_file_tamper_rejected_after_seal(tmp_path):
    root = capture(tmp_path / "source")
    f = root / "scope_0_exact/profiles/organism_0/cgmlst_cgmlst-573-2_0.bin"
    f.write_text("Genome ID,Gene,Allele ID\ng000,A,1\n")
    source = ProviderSource(root)
    assert source.scopes[0]["scope"]["receipts"][0]["missing_source_genome_ids"] == ["g001", "g002"]
    manifest = source.snapshot(tmp_path / "sealed")
    d = json.loads(manifest.read_text())
    asset = manifest.parent / d["scopes"][0]["receipts"][0]["saved_response_path"]
    asset.write_text("Genome ID,Gene,Allele ID\ng000,A,2\n")
    with pytest.raises(DatasetError, match="receipt"):
        ProviderSource(manifest)


def test_recorded_request_receipt_never_relabelled_legacy(tmp_path):
    root = capture(tmp_path / "source")
    manifest = ProviderSource(root).snapshot(tmp_path / "sealed")
    d = json.loads(manifest.read_text())
    d["scopes"][0]["receipts"][0].update(
        request_evidence="recorded_exact_http_request", retry_attempt=2
    )
    write_json(manifest, d)
    rehash(manifest, "snapshot_sha256")
    source = ProviderSource(manifest)
    assert (
        source.scopes[0]["scope"]["receipts"][0]["request_evidence"]
        == "recorded_exact_http_request"
    )


def test_ordinary_cli_offline_replay_invokes_current_normalizer(tmp_path, monkeypatch):
    from chronoclade.cli import app
    from chronoclade.prepare_stage import run_prepare
    import chronoclade.pasteur_profile_normalization as norm
    import chronoclade.profile_context_provider as provider

    root = capture(tmp_path / "source")
    source = tmp_path / "query.json"
    cats = tmp_path / "catalogues.json"
    row = dict(
        sample_id="q",
        species="Klebsiella pneumoniae",
        origin="local",
        mlst_st="147",
        collection_date="2018",
        country="Greece",
        cgmlst_scheme="native",
        cgmlst_scheme_version="v1",
        cgmlst_loci=["A", "B"],
        cgmlst_locus_universe_complete=True,
        cgmlst_profile={"A": "1", "B": "1"},
        cgmlst_status="resolved",
        cglin_scheme="scgMLST629_S",
        cglin_scheme_version="current",
        cglin_raw="0,0,1,2,3,4,5",
        cglin_status="complete",
    )
    write_json(source, [row])
    write_json(cats, [asdict(LocusCatalogue("native", "v1", ("A", "B"), "fixture"))])
    prepared = run_prepare(source, tmp_path / "prepared", catalogues=cats)
    calls = []

    class Normalizer:
        def __init__(self, *args):
            self.audits = []

        def normalize(self, rows):
            calls.append([r["sample_id"] for r in rows])
            return [
                dict(
                    r,
                    cgmlst_scheme="native",
                    cgmlst_scheme_version="v1",
                    cgmlst_loci=["A", "B"],
                    cgmlst_locus_universe_complete=True,
                    cglin_scheme="scgMLST629_S",
                    cglin_scheme_version="current",
                    cglin_raw="0,0,1,2,3,4,5",
                    cglin_status="complete",
                )
                for r in rows
            ]

    monkeypatch.setattr(norm, "resolve_reference", lambda *a: tmp_path / "config.json")
    monkeypatch.setattr(norm, "reference_descriptor", lambda *a: {"config_sha256": "a" * 64})
    monkeypatch.setattr(norm, "PasteurProfileNormalizer", Normalizer)

    def forbidden(*args, **kwargs):
        raise AssertionError("Network/client path forbidden in frozen source replay")

    import socket

    monkeypatch.setattr(socket, "create_connection", forbidden)
    monkeypatch.setattr(provider, "PathogenwatchClient", forbidden)
    monkeypatch.setattr(provider, "_grouped_analysis_exports", forbidden)
    result = CliRunner().invoke(
        app,
        [
            "cgmlst",
            str(prepared.dataset_manifest),
            "--partition-by", "lineage",
            "--out",
            str(tmp_path / "run"),
            "--context-source",
            str(root),
            "--bootstrap-replicates",
            "0",
        ],
    )
    assert result.exit_code == 0, result.output + str(result.exception)
    assert calls == [["q"], ["PW_g000", "PW_g001", "PW_g002"]]
    context = load_dataset(tmp_path / "run_context/dataset.json")
    assert len(context.samples) == 4
    assert (
        context.parameters["profile_context"]["raw_provider_source"]["schema"]
        == "chronoclade.provider-source"
    )
    assert (
        context.parameters["profile_context"]["pools"][0]["exports"]["as_of"]
        == "2026-10-17T00:00:00Z"
    )
    assert (tmp_path / "run_context/sources/provider/raw_source/provider-source.json").is_file()
    assert context.parameters["profile_context"]["source_mode"] == "frozen_raw_replay"
    assert context.parameters["profile_context"]["source_as_of"][0]["st"] == "147"
    original_execution = provider.execution_provenance
    monkeypatch.setattr(
        provider,
        "execution_provenance",
        lambda *a, **k: dict(original_execution(*a, **k), test_method_revision="changed"),
    )
    second = CliRunner().invoke(
        app,
        [
            "cgmlst",
            str(prepared.dataset_manifest),
            "--partition-by", "lineage",
            "--out",
            str(tmp_path / "second"),
            "--context-source",
            str(root),
            "--bootstrap-replicates",
            "0",
        ],
    )
    assert second.exit_code == 0, second.output + str(second.exception)
    assert calls == [["q"], ["PW_g000", "PW_g001", "PW_g002"]] * 2
    fresh = load_dataset(tmp_path / "second_context/dataset.json")
    assert (
        fresh.parameters["profile_context"]["raw_provider_source"]
        == context.parameters["profile_context"]["raw_provider_source"]
    )
    assert (
        fresh.parameters["profile_context"]["resume_fingerprint"]["execution"]
        != context.parameters["profile_context"]["resume_fingerprint"]["execution"]
    )


@pytest.mark.parametrize(
    "field,value",
    [
        ("base_url", "https://different.example"),
        ("job", "unadvertised-job"),
        ("source_genome_ids", ["alien"]),
        ("numeric_request_ids", ["3", "2", "1"]),
        ("route", "/api/downloads/cgmlst?job=wrong"),
    ],
)
def test_sealed_receipts_bind_scope_job_order_and_provider(tmp_path, field, value):
    root = capture(tmp_path / "source")
    manifest = ProviderSource(root).snapshot(tmp_path / "sealed")
    data = json.loads(manifest.read_text())
    data["scopes"][0]["receipts"][0][field] = value
    write_json(manifest, data)
    rehash(manifest, "snapshot_sha256")
    with pytest.raises(DatasetError, match="receipt"):
        ProviderSource(manifest)


def test_partial_batch_manifest_and_malformed_scope_rejected(tmp_path):
    root = capture(tmp_path / "source", 101)
    manifest = ProviderSource(root).snapshot(tmp_path / "sealed")
    data = json.loads(manifest.read_text())
    data["scopes"][0]["receipts"].pop()
    write_json(manifest, data)
    rehash(manifest, "snapshot_sha256")
    with pytest.raises(DatasetError, match="incomplete"):
        ProviderSource(manifest)
    data["scopes"] = [None]
    write_json(manifest, data)
    rehash(manifest, "snapshot_sha256")
    with pytest.raises(DatasetError, match="Malformed"):
        ProviderSource(manifest)


def test_self_consistent_subset_cannot_override_raw_search_total(tmp_path):
    root = capture(tmp_path / "source")
    path = root / "public_0.json"
    value = json.loads(path.read_text())
    value["rows"] = value["rows"][:2]
    value["raw_details"].pop("g002")
    value["raw_search_pages"][0]["genomes"] = value["raw_search_pages"][0]["genomes"][:2]
    value["provenance"].update(
        expected_count=2,
        raw_record_count=2,
        rows_sha256=content_hash(value["rows"]),
        raw_details_sha256=content_hash(value["raw_details"]),
        raw_search_sha256=content_hash(value["raw_search_pages"]),
    )
    write_json(path, value)
    rehash(path, "snapshot_sha256")
    with pytest.raises(DatasetError, match="page counts"):
        ProviderSource(root)


def test_future_atomic_batches_preserve_recorded_request_receipts_in_source_roundtrip(
    tmp_path, monkeypatch
):
    import shutil
    from types import SimpleNamespace
    from chronoclade.profile_inputs import _grouped_analysis_exports
    import chronoclade.profile_inputs as inputs

    root = capture(tmp_path / "source", 3)
    # Replace only this test-generated legacy export with a genuine shared-helper capture.
    shutil.rmtree(root / "scope_0_exact")
    catalogue = load_catalogue(root / "public_0.json")
    rows = catalogue["rows"]

    def group(*args, **kwargs):
        return [
            {
                "ids": kwargs["body"]["ids"],
                "organismId": "573",
                "downloads": [{"name": "cgmlst", "job": "cgmlst-573-2"}],
            }
        ]

    client = SimpleNamespace(
        api_key="not-persisted", base_url="https://pathogen.watch", request_json=group
    )

    def download(route, **kwargs):
        ids = kwargs["body"]["ids"].split(",")
        return (
            "Genome ID,Gene,Allele ID\n" + "".join(f"g{int(i) - 1:03d},A,{i}\n" for i in ids)
        ).encode()

    monkeypatch.setattr(inputs, "request_download", download)
    values = _grouped_analysis_exports(rows, client, root / "scope_0/profiles", {})
    assert len(values) == 3
    source = ProviderSource(root)
    receipts = source.scopes[0]["scope"]["receipts"]
    assert all(r["request_evidence"] == "recorded_exact_http_request" for r in receipts)
    restored = ProviderSource(source.snapshot(tmp_path / "portable"))
    assert (
        restored.for_scope("Klebsiella pneumoniae", "147")[2]
        == source.for_scope("Klebsiella pneumoniae", "147")[2]
    )
    assert (
        restored.scopes[0]["scope"]["receipts"][0]["request_evidence"]
        == "recorded_exact_http_request"
    )
