import io
import json
import hashlib
import threading
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor

import pytest

from chronoclade.file_snapshots import copy_evidence_snapshot
from chronoclade.pathogenwatch_download import (
    DownloadError,
    download_assemblies,
    unpack_bulk,
    validate_fasta,
)


def test_download_resume_verifies_content_and_recovers_corruption(tmp_path):
    calls = []

    def fetch(source):
        calls.append(source)
        return b">contig\nACGTN\n"

    arguments = dict(
        rows=[{"source_genome_id": "UUID1"}],
        output=tmp_path / "out",
        cache_dir=tmp_path / "cache",
        api_key="secret",
        fetch=fetch,
    )
    paths, ledger = download_assemblies(**arguments)
    assert paths["UUID1"].read_bytes() == b">contig\nACGTN\n"
    assert ledger[0]["status"] == "downloaded"
    _, ledger = download_assemblies(**arguments)
    assert calls == ["UUID1"] and ledger[0]["status"] == "cached"
    meta = json.loads((tmp_path / "cache/UUID1.json").read_text())
    (tmp_path / "cache" / f"{meta['sha256']}.fasta").write_text("corrupt")
    _, ledger = download_assemblies(**arguments)
    assert calls == ["UUID1", "UUID1"] and ledger[0]["status"] == "downloaded"
    assert "secret" not in json.dumps(ledger)


def test_failures_preserve_successful_downloads(tmp_path):
    paths, ledger = download_assemblies(
        [{"source_genome_id": "a"}, {"source_genome_id": "b"}],
        output=tmp_path / "out",
        cache_dir=tmp_path / "cache",
        api_key="secret",
        fetch=lambda source: b">ok\nACGT\n" if source == "a" else b"<html>error</html>",
    )
    assert set(paths) == {"a"}
    assert [row["status"] for row in ledger] == ["downloaded", "failed"]
    assert not (tmp_path / "out/b.fasta").exists()


def test_bulk_matches_ids_not_archive_order_and_retains_partial_batch():
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w") as archive:
        archive.writestr("b.fasta", ">b\nACGT\n")
    assert set(unpack_bulk(stream.getvalue(), {"a": "UUIDa", "b": "UUIDb"})) == {"UUIDb"}
    with pytest.raises(DownloadError, match="joined"):
        unpack_bulk(stream.getvalue(), {"a": "UUIDa"})


@pytest.mark.parametrize("data", [b"", b">a\n", b"<html>", b">a\nACGTX\n", b">a\nACGT\n>b\n"])
def test_invalid_fasta_rejected(data):
    with pytest.raises(DownloadError):
        validate_fasta(data)


def test_identical_contents_download_concurrently(tmp_path):
    paths, ledger = download_assemblies(
        [{"source_genome_id": str(i)} for i in range(4)],
        output=tmp_path / "out",
        cache_dir=tmp_path / "cache",
        api_key="secret",
        workers=4,
        fetch=lambda source: b">contig\nACGT\n",
    )
    assert len(paths) == 4
    assert all(row["status"] == "downloaded" for row in ledger)
    assert len(list((tmp_path / "cache").glob("*.fasta"))) == 1


def test_provider_revision_invalidates_cached_assembly(tmp_path):
    calls = []

    def fetch(source):
        calls.append(source)
        return b">a\nACGT\n" if len(calls) == 1 else b">a\nTGCA\n"

    row = {"source_genome_id": "a", "source_checksum": "revision1"}
    kwargs = dict(
        output=tmp_path / "out", cache_dir=tmp_path / "cache", api_key="secret", fetch=fetch
    )
    download_assemblies([row], **kwargs)
    row["source_checksum"] = "revision2"
    paths, ledger = download_assemblies([row], **kwargs)
    assert len(calls) == 2 and ledger[0]["status"] == "downloaded"
    assert b"TGCA" in paths["a"].read_bytes()


def test_wrong_source_content_is_rejected(tmp_path):
    row = {
        "source_genome_id": "a",
        "source_checksum": hashlib.sha1(b">a\nACGT\n").hexdigest(),
        "source_length": 4,
    }
    paths, ledger = download_assemblies(
        [row],
        output=tmp_path / "out",
        cache_dir=tmp_path / "cache",
        api_key="secret",
        fetch=lambda source: b">a\nTGCA\n",
    )
    assert not paths and ledger[0]["reason"] == "source_checksum_mismatch"
    assert not list((tmp_path / "cache").glob("*.fasta"))


def _public_detail(source, checksum):
    return {
        "uuid": source,
        "projectAccess": "PUBLIC",
        "id": 417,
        "organismId": "573",
        "mlst": "147",
        "species": "Klebsiella pneumoniae",
        "checksum": checksum,
        "length": 5,
        "qc": True,
    }


def test_missing_checksum_requires_exact_public_identity_and_verifies_payload(
    tmp_path, monkeypatch
):
    content = b">contig\nACGTN\n"
    checksum = hashlib.sha1(content).hexdigest()
    identity_calls = []
    download_calls = []
    class PublicClient:
        def __init__(self, *, base_url, api_key):
            assert base_url == "https://pathogen.watch" and api_key == ""

        def request_json(self, method, route, *, params, authenticated):
            identity_calls.append((method, route, params, authenticated))
            detail = _public_detail("public-id", checksum)
            detail.pop("projectAccess")  # Details endpoint omits search visibility.
            return detail

    monkeypatch.setattr("chronoclade.pathogenwatch.PathogenwatchClient", PublicClient)
    row = {
        "source_genome_id": "public-id",
        "species": "Klebsiella pneumoniae",
        "organism_id": "573",
        "mlst_scheme": "mlst",
        "mlst_st": "147",
    }
    original_row = dict(row)
    paths, ledger = download_assemblies(
        [row],
        output=tmp_path / "out",
        cache_dir=tmp_path / "cache",
        api_key="secret",
        require_public_identity=True,
        fetch=lambda source: download_calls.append(source) or content,
    )
    assert identity_calls == [
        ("GET", "/api/genomes/details", {"id": "public-id"}, False)
    ]
    assert download_calls == ["public-id"]
    assert paths["public-id"].read_bytes() == content
    assert ledger[0]["source_checksum_verified"] is True
    assert ledger[0]["source_identity_verified"] is True
    assert ledger[0]["source_checksum"] == checksum and ledger[0]["source_length"] == 5
    assert ledger[0]["raw_details_sha256"]
    assert row == original_row


@pytest.mark.parametrize(
    "detail,reason",
    [
        ({"projectAccess": "PRIVATE"}, "public_identity_unverified"),
        ({"projectAccess": "PUBLIC", "uuid": "other"}, "public_identity_unverified"),
        (
            {
                "uuid": "public-id",
                "projectAccess": "PUBLIC",
                "organismId": "573",
                "mlst": "1",
                "species": "Klebsiella pneumoniae",
            },
            "public_identity_scope_mismatch",
        ),
    ],
)
def test_missing_checksum_fails_closed_on_unverified_public_identity(tmp_path, detail, reason):
    fetch_calls = []
    row = {
        "source_genome_id": "public-id",
        "species": "Klebsiella pneumoniae",
        "organism_id": "573",
        "mlst_scheme": "mlst",
        "mlst_st": "147",
    }
    paths, ledger = download_assemblies(
        [row],
        output=tmp_path / "out",
        cache_dir=tmp_path / "cache",
        api_key="secret",
        require_public_identity=True,
        identity_fetch=lambda _source: detail,
        fetch=lambda source: fetch_calls.append(source) or b">contig\nACGT\n",
    )
    assert not paths and fetch_calls == []
    assert ledger[0]["reason"] == reason


def test_new_public_checksum_reuses_and_upgrades_matching_old_cache(tmp_path):
    source = "public-id"
    content = b">contig\nACGTN\n"
    digest = hashlib.sha256(content).hexdigest()
    checksum = hashlib.sha1(content).hexdigest()
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / f"{digest}.fasta").write_bytes(content)
    (cache / f"{source}.json").write_text(
        json.dumps({"source_genome_id": source, "sha256": digest, "source_checksum": None})
    )
    download_calls = []
    row = {"source_genome_id": source, "species": "Klebsiella pneumoniae", "mlst_st": "147"}
    original_row = dict(row)
    paths, ledger = download_assemblies(
        [row],
        output=tmp_path / "out",
        cache_dir=cache,
        api_key="secret",
        require_public_identity=True,
        identity_fetch=lambda _source: _public_detail(source, checksum),
        fetch=lambda ident: download_calls.append(ident) or b"bad",
    )
    saved = json.loads((cache / f"{source}.json").read_text())
    assert paths[source].read_bytes() == content and not download_calls
    assert ledger[0]["status"] == "cached" and ledger[0]["source_checksum_verified"]
    assert ledger[0]["source_checksum"] == checksum and ledger[0]["source_length"] == 5
    assert row == original_row
    assert saved["source_checksum"] == checksum and saved["identity_verified"] is True


def test_cached_assemblies_use_independent_snapshot_copies(tmp_path, monkeypatch):
    content = b">contig\nACGTN\n"
    checksum = hashlib.sha1(content).hexdigest()
    fetch_calls = []
    copy_calls = []
    original_copy = copy_evidence_snapshot

    def tracked_copy(source, destination):
        method = original_copy(source, destination)
        copy_calls.append((source, destination, method))
        return method

    monkeypatch.setattr(
        "chronoclade.pathogenwatch_download.copy_evidence_snapshot", tracked_copy
    )
    row = {"source_genome_id": "source", "source_checksum": checksum}
    outputs = [tmp_path / "first", tmp_path / "second"]
    ledgers = []
    paths = []
    for output in outputs:
        result, ledger = download_assemblies(
            [row],
            output=output,
            cache_dir=tmp_path / "cache",
            api_key="secret",
            fetch=lambda _source: fetch_calls.append(_source) or content,
        )
        paths.append(result["source"])
        ledgers.append(ledger[0])

    source_blob = tmp_path / "cache" / f"{hashlib.sha256(content).hexdigest()}.fasta"
    assert fetch_calls == ["source"]
    assert [entry["status"] for entry in ledgers] == ["downloaded", "cached"]
    assert len(copy_calls) == 2
    assert all(entry["copy_method"] in {"copy_on_write", "copy"} for entry in ledgers)
    assert all(call[0] == source_blob and call[1].suffix == ".partial" for call in copy_calls)
    assert all(path.read_bytes() == content for path in paths)
    assert len({source_blob.stat().st_ino, *(path.stat().st_ino for path in paths)}) == 3


def test_overlapping_calls_share_one_verified_source_download(tmp_path):
    fetch_started = threading.Event()
    second_call_started = threading.Event()
    release_fetch = threading.Event()
    calls = []
    row = {"source_genome_id": "same-source"}

    def fetch(_source):
        calls.append(_source)
        fetch_started.set()
        assert release_fetch.wait(5)
        return b">contig\nACGTN\n"

    def run(output, started=None):
        if started:
            started.set()
        return download_assemblies(
            [row],
            output=output,
            cache_dir=tmp_path / "cache",
            api_key="secret",
            fetch=fetch,
        )

    with ThreadPoolExecutor(max_workers=2) as executor:
        first = executor.submit(run, tmp_path / "out1")
        assert fetch_started.wait(2)
        second = executor.submit(run, tmp_path / "out2", second_call_started)
        assert second_call_started.wait(2)
        time.sleep(0.1)
        assert calls == ["same-source"]
        release_fetch.set()
        first_paths, first_audit = first.result(timeout=5)
        second_paths, second_audit = second.result(timeout=5)

    assert first_paths["same-source"].read_bytes() == b">contig\nACGTN\n"
    assert second_paths["same-source"].read_bytes() == b">contig\nACGTN\n"
    assert sorted([first_audit[0]["status"], second_audit[0]["status"]]) == [
        "cached",
        "downloaded",
    ]
    assert calls == ["same-source"]


def test_failed_download_releases_source_cache_lock(tmp_path):
    calls = []

    def fetch(_source):
        calls.append(_source)
        if len(calls) == 1:
            return b"<html>error</html>"
        return b">contig\nACGTN\n"

    args = dict(
        rows=[{"source_genome_id": "retry-source"}],
        cache_dir=tmp_path / "cache",
        api_key="secret",
        fetch=fetch,
    )
    paths, ledger = download_assemblies(output=tmp_path / "failed", **args)
    assert not paths and ledger[0]["status"] == "failed"
    paths, ledger = download_assemblies(output=tmp_path / "retry", **args)
    assert paths["retry-source"].read_bytes() == b">contig\nACGTN\n"
    assert ledger[0]["status"] == "downloaded"
    assert calls == ["retry-source", "retry-source"]
