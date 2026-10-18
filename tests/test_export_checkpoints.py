"""Portable tests for accepted raw checkpoints, distinct from normalized caches."""

import csv
import io
import json
from types import SimpleNamespace

import pytest

from chronoclade import profile_inputs as p
from chronoclade.export_checkpoints import ExportCheckpointError, checkpoint_path
from chronoclade.pathogenwatch import content_hash

JOB = "cgmlst-573-2"


def rows(n=201):
    return [
        {
            "source_genome_id": f"id-{i}",
            "numeric_source_id": i + 1,
            "species": "Klebsiella pneumoniae",
            "organism_id": "573",
            "mlst_st": "258",
        }
        for i in range(n)
    ]


def response(numeric):
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(["Genome ID", "Gene", "Allele ID"])
    for number in numeric:
        writer.writerow([f"id-{int(number) - 1}", "g1", "1"])
    return out.getvalue().encode()


def client(calls):
    def request(method, path, **kwargs):
        calls.append(("groups", kwargs["body"]["ids"]))
        return [
            {
                "ids": kwargs["body"]["ids"],
                "organismId": "573",
                "downloads": [{"name": "cgmlst", "job": JOB}],
            }
        ]

    return SimpleNamespace(
        base_url="https://pathogen.watch", api_key="NEVER-PERSIST-THIS", request_json=request
    )


def test_interrupt_after_accepted_batch_resumes_only_remaining_network(tmp_path, monkeypatch):
    calls = []
    c = client(calls)

    def initial(route, **kwargs):
        ids = kwargs["body"]["ids"].split(",")
        calls.append(("download", ids))
        if ids[0] == "101":
            raise KeyboardInterrupt("simulated process interruption")
        return response(ids)

    monkeypatch.setattr(p, "request_download", initial)
    with pytest.raises(KeyboardInterrupt):
        p._grouped_analysis_exports(rows(), c, tmp_path, {})
    assert checkpoint_path(tmp_path / "organism_0", "cgmlst", JOB, 0).is_file()
    assert not checkpoint_path(tmp_path / "organism_0", "cgmlst", JOB, 100).exists()

    def remaining(route, **kwargs):
        ids = kwargs["body"]["ids"].split(",")
        calls.append(("download", ids))
        return response(ids)

    monkeypatch.setattr(p, "request_download", remaining)
    c.request_json = lambda *a, **k: pytest.fail("frozen capabilities must resume without network")
    result = p._grouped_analysis_exports(rows(), c, tmp_path, {}, resume=True)
    assert len(result) == 201 and all(r["cgmlst_profile"] == {"g1": "1"} for r in result)
    assert sum(item[0] == "download" and item[1][0] == "1" for item in calls) == 1
    assert sum(item[0] == "groups" for item in calls) == 1
    for path in tmp_path.rglob("*.json"):
        assert "NEVER-PERSIST-THIS" not in path.read_text()


def test_all_exports_before_normalization_can_resume_without_download(tmp_path, monkeypatch):
    calls = []
    c = client(calls)
    monkeypatch.setattr(
        p, "request_download", lambda route, **kwargs: response(kwargs["body"]["ids"].split(","))
    )
    first = p._grouped_analysis_exports(rows(), c, tmp_path, {})
    capability_before = (tmp_path / "export_capabilities.json").read_bytes()
    capability_date = json.loads(capability_before)["retrieved_at"]
    # Classification is deliberately outside the raw checkpoint contract.
    with pytest.raises(RuntimeError, match="classifier"):
        raise RuntimeError("classifier interrupted after raw exports")
    monkeypatch.setattr(
        p, "request_download", lambda *a, **k: pytest.fail("no repeated raw export")
    )
    c.request_json = lambda *a, **k: pytest.fail("no repeated capability request")
    provenance = {}
    assert p._grouped_analysis_exports(rows(), c, tmp_path, provenance, resume=True) == first
    assert (tmp_path / "export_capabilities.json").read_bytes() == capability_before
    assert (
        json.loads((tmp_path / "export_capabilities.json").read_text())["retrieved_at"]
        == capability_date
    )
    assert len(provenance["analysis_export_raw_reuse"]) == 3


@pytest.mark.parametrize("change", ["bytes", "hash", "job", "order", "reference_scope", "schema", "route", "parse_identity"])
def test_bad_saved_batch_fails_closed(tmp_path, monkeypatch, change):
    c = client([])
    monkeypatch.setattr(
        p, "request_download", lambda route, **kwargs: response(kwargs["body"]["ids"].split(","))
    )
    p._grouped_analysis_exports(rows(2), c, tmp_path, {})
    cp = checkpoint_path(tmp_path / "organism_0", "cgmlst", JOB, 0)
    value = json.loads(cp.read_text())
    if change == "bytes":
        (tmp_path / value["receipt"]["saved_response_path"]).write_bytes(b"truncated")
    elif change == "hash":
        value["receipt"]["response_sha256"] = "0" * 64
    elif change == "job":
        value["contract"]["job"] = "cgmlst-foreign"
    elif change == "order":
        value["contract"]["source_genome_ids"].reverse()
    elif change == "reference_scope":
        value["contract"]["row_scope"][0]["cgmlst_scheme_version"] = "different"
    elif change == "schema":
        value["schema_version"] = 999
    elif change == "route":
        value["receipt"]["route"] = "/api/downloads/cgmlst?job=foreign"
    elif change == "parse_identity":
        import hashlib
        data=b"Genome ID,Gene,Allele ID\nforeign,g1,1\n"
        (tmp_path/value["receipt"]["saved_response_path"]).write_bytes(data)
        value["receipt"]["response_sha256"]=hashlib.sha256(data).hexdigest()
        value["receipt"]["response_bytes"]=len(data)
    value["sha256"] = content_hash({k: v for k, v in value.items() if k != "sha256"})
    cp.write_text(json.dumps(value))
    monkeypatch.setattr(
        p,
        "request_download",
        lambda *a, **k: pytest.fail("bad checkpoint cannot silently redownload"),
    )
    with pytest.raises(ExportCheckpointError):
        p._grouped_analysis_exports(rows(2), c, tmp_path, {}, resume=True)


def test_unknown_legacy_raw_file_never_adopted(tmp_path, monkeypatch):
    c = client([])
    out = tmp_path / "organism_0"
    out.mkdir()
    (out / f"cgmlst_{JOB}_0.bin").write_bytes(response(["1"]))
    monkeypatch.setattr(p, "request_download", lambda *a, **k: pytest.fail("unverified collision"))
    with pytest.raises(ExportCheckpointError, match="without accepted request receipt"):
        p._grouped_analysis_exports(rows(1), c, tmp_path, {})


def test_fresh_job_cannot_overwrite_same_output(tmp_path, monkeypatch):
    c = client([])
    monkeypatch.setattr(
        p, "request_download", lambda route, **kwargs: response(kwargs["body"]["ids"].split(","))
    )
    p._grouped_analysis_exports(rows(1), c, tmp_path, {})
    with pytest.raises(ExportCheckpointError, match="verified resume"):
        p._grouped_analysis_exports(rows(1), c, tmp_path, {})


@pytest.mark.parametrize("change", ["base_url", "order", "scope", "download_names"])
def test_capability_request_contract_mismatch_fails_closed(tmp_path, monkeypatch, change):
    c = client([])
    data = rows(2)
    monkeypatch.setattr(
        p, "request_download", lambda route, **kwargs: response(kwargs["body"]["ids"].split(","))
    )
    p._grouped_analysis_exports(data, c, tmp_path, {})
    kw = {}
    if change == "base_url":
        c.base_url = "https://different.example"
    elif change == "order":
        data.reverse()
    elif change == "scope":
        data[0]["mlst_st"] = "39"
    else:
        kw["download_names"] = {"hclink"}
        data[0].update(species="Escherichia coli", organism_id="562")
    with pytest.raises(ExportCheckpointError):
        p._grouped_analysis_exports(data, c, tmp_path, {}, resume=True, **kw)


def test_crash_between_raw_and_receipt_never_adopts_unpublished_bytes(tmp_path, monkeypatch):
    from chronoclade import export_checkpoints as checkpoints

    c = client([])
    downloads = []

    def fetch(route, **kwargs):
        downloads.append(kwargs["body"]["ids"])
        return response(kwargs["body"]["ids"].split(","))

    monkeypatch.setattr(p, "request_download", fetch)
    atomic = checkpoints.atomic_json

    def crash(path, value):
        assert (path.parent / "response.raw").is_file()
        raise KeyboardInterrupt("crash after raw stage, before receipt commit")

    monkeypatch.setattr(checkpoints, "atomic_json", crash)
    with pytest.raises(KeyboardInterrupt):
        p._grouped_analysis_exports(rows(1), c, tmp_path, {})
    checkpoint = checkpoint_path(tmp_path / "organism_0", "cgmlst", JOB, 0)
    assert not checkpoint.parent.exists()
    pending = list((tmp_path / "organism_0").glob(".pending-*/response.raw"))
    assert len(pending) == 1
    original = pending[0].read_bytes()
    monkeypatch.setattr(checkpoints, "atomic_json", atomic)
    c.request_json = lambda *a, **k: pytest.fail("capabilities already committed")
    result = p._grouped_analysis_exports(rows(1), c, tmp_path, {}, resume=True)
    assert result[0]["cgmlst_profile"] == {"g1": "1"}
    assert downloads == ["1", "1"]  # Safe refetch: unpublished raw is never adopted.
    assert pending[0].read_bytes() == original  # Preserve interruption evidence.
    assert checkpoint.is_file() and (checkpoint.parent / "response.raw").is_file()


def test_unknown_partial_published_batch_fails_closed(tmp_path, monkeypatch):
    c = client([])
    checkpoint = checkpoint_path(tmp_path / "organism_0", "cgmlst", JOB, 0)
    checkpoint.parent.mkdir(parents=True)
    (checkpoint.parent / "response.raw").write_bytes(response(["1"]))
    monkeypatch.setattr(
        p,
        "request_download",
        lambda *a, **k: pytest.fail("cannot overwrite unknown committed directory"),
    )
    with pytest.raises(ExportCheckpointError, match="Incomplete committed batch"):
        p._grouped_analysis_exports(rows(1), c, tmp_path, {})


def test_hclink_shared_raw_resume_retains_blank_and_known_levels(tmp_path, monkeypatch):
    import gzip

    ident = "ecoli-id"
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(
        ["Genome ID", "Genome Name", "Closest ST", "Distance to closest ST", "HC0", "HC50"]
    )
    writer.writerow([ident, "example", "123", "0", "", "7"])
    data = gzip.compress(out.getvalue().encode())
    calls = []

    def group(*args, **kwargs):
        calls.append("group")
        return [
            {
                "ids": [1],
                "organismId": "562",
                "downloads": [{"name": "hclink", "job": "hclink-562-4"}],
            }
        ]

    c = SimpleNamespace(
        base_url="https://pathogen.watch", api_key="not-persisted", request_json=group
    )
    row = {
        "source_genome_id": ident,
        "numeric_source_id": 1,
        "species": "Escherichia coli",
        "organism_id": "562",
    }
    monkeypatch.setattr(p, "request_download", lambda *a, **k: calls.append("download") or data)
    first = p._grouped_analysis_exports([row], c, tmp_path, {}, download_names={"hclink"})
    second = p._grouped_analysis_exports(
        [row], c, tmp_path, {}, download_names={"hclink"}, resume=True
    )
    assert first == second and calls == ["group", "download"]
    assert first[0]["hiercc_codes"] == {"HC50": "7"}
    assert first[0]["hiercc_blank_levels"] == ["HC0"]
    assert first[0]["hiercc_database_version"] is None
