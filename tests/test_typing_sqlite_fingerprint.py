"""SQLite reader locks must never split identical native typing references."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sqlite3
from threading import Barrier

import pytest

from chronoclade.query_typing import QueryTypingError, _database_hash


def database(tmp_path):
    root = tmp_path / 'index'
    root.mkdir()
    path = root / 'allele.db'
    with sqlite3.connect(path) as connection:
        connection.execute('PRAGMA journal_mode=WAL')
        connection.execute('CREATE TABLE alleles (value TEXT)')
        connection.execute("INSERT INTO alleles VALUES ('first')")
        connection.commit()
        connection.execute('PRAGMA wal_checkpoint(TRUNCATE)')
    return root, path


def test_real_concurrent_sqlite_readers_preserve_database_fingerprint(tmp_path):
    root, path = database(tmp_path)
    expected = _database_hash(root)
    barrier = Barrier(4)
    def reader(_):
        with sqlite3.connect(path.as_uri() + '?mode=ro', uri=True) as connection:
            assert connection.execute('SELECT COUNT(*) FROM alleles').fetchone()[0] == 1
            barrier.wait(timeout=10)
            return [_database_hash(root) for _ in range(5)]
    with ThreadPoolExecutor(max_workers=4) as pool:
        values = list(pool.map(reader, range(4)))
    assert all(value == expected for row in values for value in row)
    assert _database_hash(root) == expected


def test_harmless_shm_changes_and_empty_sqlite_sidecars_are_excluded(tmp_path):
    root, path = database(tmp_path)
    expected = _database_hash(root)
    for content in [b'initial reader locks', b'different concurrent locks']:
        path.with_name(path.name + '-shm').write_bytes(content)
        path.with_name(path.name + '-wal').touch()
        path.with_name(path.name + '-journal').touch()
        assert _database_hash(root) == expected
    path.with_name(path.name + '-shm').unlink()
    assert _database_hash(root) == expected


@pytest.mark.parametrize('suffix', ['-wal', '-journal'])
def test_pending_sqlite_data_is_rejected(tmp_path, suffix):
    root, path = database(tmp_path)
    path.with_name(path.name + suffix).write_bytes(b'uncheckpointed writes')
    with pytest.raises(QueryTypingError, match='uncheckpointed SQLite writes'):
        _database_hash(root)
    with pytest.raises(QueryTypingError, match='uncheckpointed SQLite writes'):
        _database_hash(path)


def test_actual_wal_writer_rejected_and_checkpointed_allele_change_detected(tmp_path):
    root, path = database(tmp_path)
    expected = _database_hash(root)
    with sqlite3.connect(path) as connection:
        connection.execute("INSERT INTO alleles VALUES ('second')")
        connection.commit()
        with pytest.raises(QueryTypingError, match='uncheckpointed SQLite writes'):
            _database_hash(root)
        connection.execute('PRAGMA wal_checkpoint(TRUNCATE)')
        assert _database_hash(root) != expected


@pytest.mark.parametrize('name', ['orphan-shm', 'plain-shm', 'plain-wal', 'plain-journal'])
def test_unrelated_or_orphan_sidecars_remain_content_bound(tmp_path, name):
    root, _ = database(tmp_path)
    (root / 'plain').write_text('Not a SQLite database')
    before = _database_hash(root)
    item = root / name
    item.write_text('first content')
    first = _database_hash(root)
    item.write_text('different content')
    assert first != before
    assert _database_hash(root) != first


def test_known_reader_sidecar_disappearing_after_scan_is_safe(tmp_path, monkeypatch):
    root, path = database(tmp_path)
    expected = _database_hash(root)
    sidecar = path.with_name(path.name + '-wal')
    sidecar.touch()
    original = Path.stat
    removed = [False]
    def disappearing(item, *args, **kwargs):
        if item == sidecar and not removed[0]:
            sidecar.unlink()
            removed[0] = True
        return original(item, *args, **kwargs)
    monkeypatch.setattr(Path, 'stat', disappearing)
    assert _database_hash(root) == expected
    assert removed[0]


def test_disappearing_stable_resource_is_not_silently_ignored(tmp_path, monkeypatch):
    root, _ = database(tmp_path)
    stable = root / 'metadata.json'
    stable.write_text('stable scientific data')
    original = Path.stat
    removed = [False]
    def disappearing(item, *args, **kwargs):
        if item == stable and not removed[0]:
            stable.unlink()
            removed[0] = True
        return original(item, *args, **kwargs)
    monkeypatch.setattr(Path, 'stat', disappearing)
    with pytest.raises(QueryTypingError, match='resource disappeared'):
        _database_hash(root)
