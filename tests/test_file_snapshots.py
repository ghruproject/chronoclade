"""COW snapshots preserve independent evidence and safely fall back."""

import pytest

from chronoclade.errors import WorkflowError

def test_evidence_snapshot_is_independent_and_fallback_cleans_partial_file(tmp_path, monkeypatch):
    import chronoclade.file_snapshots as snapshots

    source = tmp_path / "source.npy"
    source.write_bytes(b"original matrix data")
    destination = tmp_path / "snapshot.npy"
    method = snapshots.copy_evidence_snapshot(source, destination)
    assert method in {"copy_on_write", "copy"}
    assert source.stat().st_ino != destination.stat().st_ino
    source.write_bytes(b"changed")
    assert destination.read_bytes() == b"original matrix data"
    destination.write_bytes(b"snapshot changed")
    assert source.read_bytes() == b"changed"

    def failed_clone(original, cloned):
        cloned.write_bytes(b"partial")
        return False

    monkeypatch.setattr(snapshots, "_clone_file", failed_clone)
    fallback = tmp_path / "fallback.npy"
    assert snapshots.copy_evidence_snapshot(source, fallback) == "copy"
    assert fallback.read_bytes() == b"changed"
    with pytest.raises(WorkflowError, match="already exists"):
        snapshots.copy_evidence_snapshot(source, fallback)


def test_evidence_snapshot_rejects_dangling_destination_symlink(tmp_path):
    from chronoclade.file_snapshots import copy_evidence_snapshot

    source = tmp_path / "source"
    source.write_bytes(b"source")
    missing = tmp_path / "missing"
    destination = tmp_path / "snapshot"
    destination.symlink_to(missing)
    with pytest.raises(WorkflowError, match="already exists"):
        copy_evidence_snapshot(source, destination)
    assert destination.is_symlink()
    assert not missing.exists()


def test_selection_snapshot_alias_preserved():
    from chronoclade.file_snapshots import copy_evidence_snapshot
    from chronoclade.selection_manifest import copy_evidence_snapshot as alias

    assert alias is copy_evidence_snapshot
