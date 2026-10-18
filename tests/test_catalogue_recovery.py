"""An interrupted capture must recover work without publishing a partial pool."""
import json
from pathlib import Path

import pytest

from chronoclade.pathogenwatch import PathogenwatchClient, PathogenwatchError, load_catalogue


class PublicService:
    def __init__(self, *, fail_page=False, fail_detail=False):
        self.calls = []
        self.fail_page = fail_page
        self.fail_detail = fail_detail

    def __call__(self, method, path, *, body, params, headers):
        assert "X-API-Key" not in headers
        self.calls.append((path, dict(params)))
        if path.endswith("supported"):
            return [{"organismId": "573", "typing": ["MLST"],
                     "fullName": "Klebsiella pneumoniae"}]
        if path.endswith("genomes"):
            if params.get("after") and self.fail_page:
                raise PathogenwatchError("Interrupted search")
            uuid = "b" if params.get("after") else "a"
            return {"meta": {"count": 2, "endCursor": uuid},
                    "genomes": [{"uuid": uuid, "projectAccess": "PUBLIC"}]}
        uuid = params["id"]
        if uuid == "b" and self.fail_detail:
            raise PathogenwatchError("Interrupted details")
        return {"uuid": uuid, "organismId": "573", "mlst": "147", "name": uuid,
                "projectAccess": "PUBLIC", "metadata": {"Date": "2019"}}


def client(service):
    return PathogenwatchClient(api_key="protected-value", transport=service, workers=1)


@pytest.mark.parametrize("failure", ["page", "detail"])
def test_resume_reuses_only_checked_responses_and_publishes_complete(tmp_path, failure):
    destination = tmp_path / "pool.json"
    service = PublicService(fail_page=failure == "page", fail_detail=failure == "detail")
    with pytest.raises(PathogenwatchError, match="Interrupted"):
        client(service).freeze_catalogue(destination, page_size=1)
    assert not destination.exists()
    checkpoint = tmp_path / "pool.json.retrieval"
    assert checkpoint.exists()
    service = PublicService()
    envelope = client(service).freeze_catalogue(destination, page_size=1, resume=True)
    assert load_catalogue(destination) == envelope
    assert len(envelope["raw_details"]) == 2
    assert envelope["provenance"]["complete"] is True
    assert envelope["provenance"]["reused_search_pages"] == (1 if failure == "page" else 2)
    assert envelope["provenance"]["reused_details"] == (0 if failure == "page" else 1)
    assert not any(path.endswith("genomes") and not params.get("after")
                   for path, params in service.calls)
    assert "protected-value" not in "".join(p.read_text() for p in tmp_path.rglob("*.json"))


def interrupted(tmp_path):
    destination = tmp_path / "pool.json"
    with pytest.raises(PathogenwatchError):
        client(PublicService(fail_detail=True)).freeze_catalogue(destination, page_size=1)
    return destination


def test_resume_rejects_query_or_page_size_changes_before_network(tmp_path):
    path = interrupted(tmp_path)
    service = PublicService()
    for kwargs in ({"st": "1", "page_size": 1}, {"page_size": 2}):
        with pytest.raises(PathogenwatchError, match="binding mismatch"):
            client(service).freeze_catalogue(path, resume=True, **kwargs)
    assert not service.calls
    assert not path.exists()


@pytest.mark.parametrize("kind", ["pages", "details"])
def test_tampered_saved_response_refuses_completion(tmp_path, kind):
    path = interrupted(tmp_path)
    receipt_path = next((tmp_path / "pool.json.retrieval" / kind).glob("*.json"))
    receipt = json.loads(receipt_path.read_text())
    receipt["response"]["tampered"] = True
    receipt_path.write_text(json.dumps(receipt))
    with pytest.raises(PathogenwatchError, match="checksum"):
        client(PublicService()).freeze_catalogue(path, page_size=1, resume=True)
    assert not path.exists()


def test_existing_capture_requires_explicit_resume_and_rejects_symlink(tmp_path):
    path = interrupted(tmp_path)
    with pytest.raises(PathogenwatchError, match="resume"):
        client(PublicService()).freeze_catalogue(path, page_size=1)
    receipt_path = next((tmp_path / "pool.json.retrieval" / "details").glob("*.json"))
    original = receipt_path.read_text()
    receipt_path.unlink()
    target = tmp_path / "external.json"
    target.write_text(original)
    receipt_path.symlink_to(target)
    with pytest.raises(PathogenwatchError, match="symlink"):
        client(PublicService()).freeze_catalogue(path, page_size=1, resume=True)
    assert not path.exists()


def test_simultaneous_writer_is_refused(tmp_path):
    from chronoclade.prepare_recovery import DirectoryLock

    path = tmp_path / "pool.json"
    with DirectoryLock(Path(str(path) + ".retrieval.lock")):
        with pytest.raises(PathogenwatchError, match="Another process"):
            client(PublicService()).freeze_catalogue(path, page_size=1)
    assert not path.exists()


def test_preflight_count_drift_refused_before_details(tmp_path):
    service = PublicService()
    destination = tmp_path / "pool.json"
    with pytest.raises(PathogenwatchError, match="changed since context preflight"):
        client(service).freeze_catalogue(destination, expected_count=1)
    assert not any(path.endswith("details") for path, _ in service.calls)
    assert not destination.exists()
