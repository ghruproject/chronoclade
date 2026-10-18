import json
import csv
import gzip
import io
from pathlib import Path

from chronoclade import profile_inputs
from chronoclade.pathogenwatch import PathogenwatchClient
from chronoclade.datasets import load_dataset
from chronoclade.prepare_provider import _has_lineage
from chronoclade.prepare_stage import run_prepare

FIXTURE = Path(__file__).parent / "fixtures" / "pathogenwatch-hclink-562-4-id108048.csv.gz"
IDENT = "iLGx7naVnHXTvwPF8DEgAo"


def cgmlst_bytes():
    stream = io.StringIO(newline="")
    writer = csv.writer(stream)
    writer.writerow(["Genome ID", "Gene", "Allele ID"])
    writer.writerow([IDENT, "locus1", "1"])
    return gzip.compress(stream.getvalue().encode())


def test_live_prepare_resolution_fetches_both_exports_and_preserves_unknown_reference_scope(
    tmp_path, monkeypatch
):
    calls = []
    downloads = [
        {"name": "cgmlst", "job": "cgmlst-562-1"},
        {"name": "hclink", "job": "hclink-562-4"},
    ]

    def transport(method, path, *, body, params, headers):
        calls.append((method, path, body, params))
        if path == "/api/organisms/supported":
            return [{"fullName": "Escherichia coli", "organismId": 562}]
        if path == "/api/search/genomes":
            return {"meta": {"count": 1}, "genomes": [{"uuid": IDENT}]}
        if path == "/api/genomes/details":
            return {
                "uuid": IDENT,
                "id": 108048,
                "name": "SAMD00110466",
                "organismId": 562,
                "projectAccess": "PUBLIC",
                "runAccession": "DRR120871",
                "sampleAccession": "SAMD00110466",
                "mlst": "131",
            }
        if path == "/api/genomes/group":
            return [{"organismId": 562, "ids": [108048], "downloads": downloads}]
        raise AssertionError((method, path))

    def request_download(route, *, api_key, body, base_url):
        assert api_key == "fixture-only-secret"
        assert body == {"ids": "108048"}
        if route == "/api/downloads/cgmlst?job=cgmlst-562-1":
            return cgmlst_bytes()
        if route == "/api/downloads/hclink?job=hclink-562-4":
            return FIXTURE.read_bytes()
        raise AssertionError(route)

    monkeypatch.setattr(profile_inputs, "request_download", request_download)
    client = PathogenwatchClient(api_key="fixture-only-secret", transport=transport)
    # Keep all inputs confined to this test's output.
    accessions = tmp_path / "accessions.csv"
    accessions.write_text("accession\nDRR120871\n")
    prepared = run_prepare(
        accessions,
        tmp_path / "prepared",
        input_kind="accessions",
        species="Escherichia coli",
        client=client,
        enrich_metadata=False,
    )
    dataset = load_dataset(prepared.dataset_manifest)
    assert len(dataset.profiles) == 1
    assert dataset.profiles[0].sample_ids == ("PW_" + IDENT,)
    assert dataset.profiles[0].codes.tolist() == [[1]]
    bundle = prepared.dataset_manifest.parent
    provider_snapshot = json.loads(
        (bundle / "sources/provider/live_profiles.json").read_text()
    )
    record = provider_snapshot["records"][0]
    assert record["source_genome_id"] == IDENT
    assert record["cgmlst_profile"] == {"locus1": "1"}
    assert record["hiercc_codes"]["HC50"] == "7"
    assert record["hiercc_codes"]["HC2350"] == "1"
    assert set(record["hiercc_blank_levels"]) == {"HC0", "HC2", "HC5", "HC10", "HC20"}
    assert record["hiercc_scheme_version"] is None
    assert record["hiercc_database_version"] is None
    assert record["hiercc_database_sha256"] is None
    assert not _has_lineage(record)
    hclink = next(row for row in dataset.lineages if row["kind"] == "hiercc")
    assert hclink["assignment"] == record["hiercc_codes"]
    assert hclink["database_version"] is None
    assert hclink["assignment_method"] == "pathogenwatch_hclink_analysis_export"
    assert hclink["evidence"]["hiercc_blank_levels"] == ["HC0", "HC2", "HC5", "HC10", "HC20"]
    audit = provider_snapshot["parameters"]["provider_prepare"]["analysis_exports"]
    assert [(item["download"], item["job"], item["status"]) for item in audit] == [
        ("cgmlst", "cgmlst-562-1", "downloaded"),
        ("hclink", "hclink-562-4", "downloaded"),
    ]
    assert "database release, version and fingerprint unavailable" in audit[1]["scope_provenance"]
    for item in audit:
        receipt = item["batch_receipts"][0]
        assert receipt["request_evidence"] == "recorded_exact_http_request"
        assert receipt["base_url"] == client.base_url
        assert receipt["endpoint"] == f"/api/downloads/{item['download']}"
        assert receipt["source_genome_ids"] == [IDENT]
        assert receipt["numeric_request_ids"] == ["108048"]
        assert receipt["retry_attempt"] == 1
        assert receipt["returned_source_genome_ids"] == [IDENT]
        assert receipt["missing_source_genome_ids"] == []
        assert len(receipt["response_sha256"]) == 64
        assert receipt["response_bytes"] > 0
        assert receipt["status"] == "accepted"
    hclink_files = list(
        (bundle / "sources/provider").rglob("batch_hclink_hclink-562-4_0/response.raw")
    )
    assert len(hclink_files) == 1
    hclink_file = hclink_files[0]
    assert hclink_file.with_name("receipt.json").is_file()
    assert hclink_file.read_bytes() == FIXTURE.read_bytes()
    assert sum(path == "/api/genomes/group" for _, path, _, _ in calls) == 1
    readiness = json.loads((bundle / "readiness.json").read_text())
    assert readiness["readiness"] == "incomplete"
    assert readiness["samples"][0]["lineage_ready"] is False
    assert any("known scope is unavailable" in reason for reason in readiness["samples"][0]["reasons"])


def test_grouped_hclink_is_independent_of_existing_cgmlst_profile(tmp_path, monkeypatch):
    row = {
        "sample_id": "query",
        "species": "Escherichia coli",
        "organism_id": "562",
        "numeric_source_id": 108048,
        "source_genome_id": IDENT,
        "cgmlst_profile": {"locus1": "1"},
    }
    client = PathogenwatchClient(
        api_key="fixture-only-secret",
        transport=lambda method, path, **kwargs: [
            {"organismId": 562, "ids": [108048], "downloads": [{"name": "hclink", "job": "hclink-562-4"}]}
        ] if path == "/api/genomes/group" else (_ for _ in ()).throw(AssertionError(path)),
    )
    monkeypatch.setattr(profile_inputs, "request_download", lambda *args, **kwargs: FIXTURE.read_bytes())
    result = profile_inputs._grouped_analysis_exports(
        [row], client, tmp_path / "exports", {}, download_names={"hclink"}
    )
    assert result[0]["cgmlst_profile"] == {"locus1": "1"}
    assert result[0]["hiercc_codes"]["HC50"] == "7"


def test_hclink_skips_non_ecoli_and_existing_hiercc(tmp_path, monkeypatch):
    called = []

    def transport(method, path, **kwargs):
        called.append(path)
        raise AssertionError("existing or non-E. coli assignments must not trigger capability retrieval")

    client = PathogenwatchClient(api_key="fixture-only-secret", transport=transport)
    rows = [
        {"species": "Klebsiella pneumoniae", "numeric_source_id": 1, "source_genome_id": "kp"},
        {"species": "Escherichia coli", "organism_id": "562", "numeric_source_id": 2,
         "source_genome_id": "ec", "hiercc_codes": {"HC10": "1"}},
        {"species": "Escherichia coli", "organism_id": "562", "numeric_source_id": 3,
         "source_genome_id": "ec-already", "hiercc_export_sha256": "a" * 64},
    ]
    result = profile_inputs._grouped_analysis_exports(
        rows, client, tmp_path / "exports", {}, download_names={"hclink"}
    )
    assert result == rows
    assert called == []


def test_klebsiella_lineage_route_still_filters_to_advertised_cglin(tmp_path, monkeypatch):
    row = {
        "sample_id": "kp",
        "species": "Klebsiella pneumoniae",
        "organism_id": "573",
        "numeric_source_id": 5730,
        "source_genome_id": "kp-source",
    }
    client = PathogenwatchClient(
        api_key="fixture-only-secret",
        transport=lambda method, path, **kwargs: [
            {
                "organismId": 573,
                "ids": [5730],
                "downloads": [
                    {"name": "cgmlst", "job": "cgmlst-573-1"},
                    {"name": "hclink", "job": "hclink-562-4"},
                    {"name": "klebsiella-lincodes", "job": "klebsiella-lincodes-1"},
                ],
            }
        ] if path == "/api/genomes/group" else (_ for _ in ()).throw(AssertionError(path)),
    )
    seen = []

    def fake_exports(members, advertised, _client, _out, _provenance, *, resume=False):
        assert resume is False
        seen.extend(item["name"] for item in advertised)
        return members

    monkeypatch.setattr(profile_inputs, "_analysis_exports", fake_exports)
    profile_inputs._grouped_analysis_exports(
        [row], client, tmp_path / "exports", {}, download_names={"klebsiella-lincodes"}
    )
    assert seen == ["klebsiella-lincodes"]
