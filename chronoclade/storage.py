"""Lossless managed compression without changing scientific paths or checksums.

A logical file keeps its original name and SHA256 in scientific manifests. Its
adjacent storage receipt describes independently checked gzip/zstd bytes. Readers
stream these bytes; consumers needing a real file restore it atomically. CLI
sessions evict only unchanged temporary restorations when the command exits.
"""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
import fcntl
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import stat
import tempfile


class StorageError(ValueError):
    """Compressed evidence could not be checked or restored safely."""


_SUFFIX = ".cc-storage.json"
_CHUNK = 1024 * 1024
_SESSION = ContextVar("chronoclade_storage_session", default=None)
_ELIGIBLE = {".fasta", ".fa", ".fna", ".phy", ".skf", ".npy", ".jsonl"}
_CONTROLS = {"job.json", "dataset.json", "prepare.json", "cgmlst.json", "tree.json", "time.json", "resume.json"}


def receipt_path(path):
    path = Path(path)
    return path.with_name(path.name + _SUFFIX)


def _archive(path, codec):
    return Path(str(path) + (".cc.zst" if codec == "zstd" else ".cc.gz"))


def _regular(path):
    return path.is_file() and not path.is_symlink()


def available(path):
    path = Path(path)
    return path.is_file() or receipt_path(path).is_file()


def _raw_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def _record(path):
    receipt = receipt_path(path)
    try:
        if not _regular(receipt):
            raise StorageError(f"Storage receipt is missing or is a symlink: {receipt}")
        data = json.loads(receipt.read_text())
        if (not isinstance(data, dict) or data.get("schema") != "chronoclade.compressed-file"
                or data.get("schema_version") != 1 or data.get("logical_name") != Path(path).name
                or data.get("codec") not in {"gzip", "zstd"}):
            raise StorageError(f"Invalid storage receipt: {receipt}")
        for key in ("sha256", "stored_sha256"):
            if not isinstance(data.get(key), str) or not re.fullmatch(r"[0-9a-f]{64}", data[key]):
                raise StorageError(f"Invalid storage checksum: {receipt}")
        for key in ("bytes", "stored_bytes"):
            if type(data.get(key)) is not int or data[key] < 0:
                raise StorageError(f"Invalid storage size: {receipt}")
        if "mode" in data and (type(data["mode"]) is not int or not 0 <= data["mode"] <= 0o777):
            raise StorageError(f"Invalid storage permissions: {receipt}")
        archive = _archive(path, data["codec"])
        if not _regular(archive) or archive.stat().st_size != data["stored_bytes"]:
            raise StorageError(f"Compressed artifact missing or size mismatch: {archive}")
        return data, archive
    except (OSError, ValueError, TypeError) as error:
        if isinstance(error, StorageError):
            raise
        raise StorageError(f"Cannot read storage receipt {receipt}: {error}") from error


def _zstd():
    try:
        import zstandard
        return zstandard
    except ImportError as error:
        raise StorageError("Zstandard support requires the zstandard package; reinstall ChronoClade") from error


@contextmanager
def _decoded(archive, codec):
    try:
        with archive.open("rb") as raw:
            if codec == "gzip":
                with gzip.GzipFile(fileobj=raw) as stream:
                    yield stream
            else:
                # The C backend uses bytes here despite the 0.25 documentation
                # saying KiB. Check the frame explicitly to bound either backend.
                zstd = _zstd()
                if zstd.get_frame_parameters(raw.read(18)).window_size > 128 * 1024**2:
                    raise StorageError(f"Compressed artifact exceeds supported 128 MiB window: {archive}")
                raw.seek(0)
                with zstd.ZstdDecompressor(max_window_size=128 * 1024**2).stream_reader(raw) as stream:
                    yield stream
    except (OSError, EOFError, _zstd().ZstdError) as error:
        raise StorageError(f"Cannot decompress {archive}: {error}") from error


@contextmanager
def _checked(path):
    data, archive = _record(path)
    if _raw_hash(archive) != data["stored_sha256"]:
        raise StorageError(f"Compressed artifact checksum mismatch: {archive}")
    with _decoded(archive, data["codec"]) as source:
        # BufferedReader supports readline/TextIOWrapper for streaming FASTA/JSON.
        with io.BufferedReader(source) as stream:
            yield stream, data


@contextmanager
def open_binary(path):
    """Read exact logical bytes; packed files are verified before use."""
    path = Path(path)
    if path.is_file() or not receipt_path(path).exists():
        with path.open("rb") as stream:
            yield stream
        return
    with _checked(path) as (stream, data):
        # Validate the full logical object even when a caller reads only a prefix.
        digest, count = hashlib.sha256(), 0
        class Reader(io.RawIOBase):
            def readable(self):
                return True
            def readinto(self, buffer):
                nonlocal count
                block = stream.read(len(buffer))
                buffer[:len(block)] = block
                digest.update(block)
                count += len(block)
                if count > data["bytes"]:
                    raise StorageError(f"Decompressed artifact exceeds recorded size: {path}")
                return len(block)
        with io.BufferedReader(Reader()) as checked:
            yield checked
            for _ in iter(lambda: checked.read(_CHUNK), b""):
                pass
        if count != data["bytes"] or digest.hexdigest() != data["sha256"]:
            raise StorageError(f"Decompressed artifact checksum/size mismatch: {path}")


@contextmanager
def open_text(path, *, encoding="utf-8", newline=None):
    """Stream managed storage, or explicit user-supplied .gz text inputs."""
    path = Path(path)
    with open_binary(path) as binary:
        if path.suffix == ".gz":
            with gzip.GzipFile(fileobj=binary) as source:
                text = io.TextIOWrapper(source, encoding=encoding, newline=newline)
                try:
                    yield text
                finally:
                    text.detach()
        else:
            text = io.TextIOWrapper(binary, encoding=encoding, newline=newline)
            try:
                yield text
            finally:
                text.detach()


def read_text(path, *, encoding="utf-8"):
    with open_text(path, encoding=encoding) as stream:
        return stream.read()


def logical_size(path):
    path = Path(path)
    return path.stat().st_size if path.is_file() else _record(path)[0]["bytes"]


@contextmanager
def _lock(path):
    lock = Path(str(path) + ".cc-lock")
    descriptor = os.open(lock, os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0), 0o600)
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise StorageError(f"Storage lock is not a regular file: {lock}")
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        yield descriptor
    finally:
        os.close(descriptor)


def _temporary(parent):
    descriptor, name = tempfile.mkstemp(prefix=".cc-storage-", dir=parent)
    return os.fdopen(descriptor, "wb"), Path(name)


def _identity(path):
    info = path.stat()
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns


def materialize(path, *, permanent=False):
    """Restore exact bytes; a session holds read leases until tools finish."""
    path = Path(path)
    path = path.parent.resolve() / path.name
    session = _SESSION.get()
    if path.is_file() and not receipt_path(path).exists():
        return path
    if session is not None and path in session.leases:
        if permanent:
            session.restored.discard(path)
        return path
    with _lock(path) as descriptor:
        restored = not path.is_file()
        if path.is_symlink():
            raise StorageError(f"Refusing to restore over a symlink: {path}")
        if restored:
            handle, temporary = _temporary(path.parent)
            try:
                with handle, open_binary(path) as source:
                    for block in iter(lambda: source.read(_CHUNK), b""):
                        handle.write(block)
                    handle.flush()
                    os.fsync(handle.fileno())
                data, _ = _record(path)
                os.chmod(temporary, data.get("mode", 0o600))
                os.replace(temporary, path)
            finally:
                temporary.unlink(missing_ok=True)
        if session is not None and not permanent:
            # Duplicating the open file description retains this shared lock
            # after _lock closes its descriptor. Pack waits for all readers.
            fcntl.flock(descriptor, fcntl.LOCK_SH)
            session.leases[path] = os.dup(descriptor)
            if restored:
                session.restored.add(path)
    return path


def _write_record(path, data):
    handle, temporary = _temporary(path.parent)
    try:
        with handle:
            handle.write((json.dumps(data, sort_keys=True, indent=2) + "\n").encode())
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, data.get("mode", 0o600))
        os.replace(temporary, receipt_path(path))
    finally:
        temporary.unlink(missing_ok=True)


def pack_file(path, *, codec="zstd", min_saving=0.05):
    """Commit a verified archive before removing its unchanged logical file."""
    path = Path(path)
    path = path.parent.resolve() / path.name
    session = _SESSION.get()
    if session is not None and path in session.leases:
        raise StorageError("Cannot pack an active native/mmap input; close its StorageSession first")
    if codec not in {"gzip", "zstd"} or not 0 <= min_saving < 1:
        raise StorageError("Storage requires gzip/zstd and a saving fraction from 0 to less than 1")
    with _lock(path):
        if not path.exists():
            with open_binary(path) as stream:
                for _ in iter(lambda: stream.read(_CHUNK), b""):
                    pass
            return {"path": str(path), "status": "already_packed"}
        if not _regular(path):
            raise StorageError(f"Compression requires a regular, non-symlink file: {path}")
        identity = _identity(path)
        mode = stat.S_IMODE(path.stat().st_mode) & 0o777
        if receipt_path(path).exists():
            data, _ = _record(path)
            if identity[2] == data["bytes"] and _raw_hash(path) == data["sha256"]:
                with _checked(path) as (stream, _):
                    restored = hashlib.sha256()
                    for block in iter(lambda: stream.read(_CHUNK), b""):
                        restored.update(block)
                if restored.hexdigest() != data["sha256"]:
                    raise StorageError(f"Decompressed artifact checksum mismatch: {path}")
                if _identity(path) != identity:
                    raise StorageError(f"Source changed before storage commit: {path}")
                path.unlink()
                return {"path": str(path), "status": "packed", "bytes": data["bytes"],
                        "stored_bytes": data["stored_bytes"], "codec": data["codec"]}
        digest = hashlib.sha256()
        handle, temporary = _temporary(path.parent)
        try:
            with handle:
                if codec == "gzip":
                    compressor = gzip.GzipFile(filename="", mode="wb", fileobj=handle, compresslevel=1, mtime=0)
                else:
                    params = _zstd().ZstdCompressionParameters.from_level(
                        3, window_log=27, enable_ldm=True, write_checksum=True,
                    )
                    compressor = _zstd().ZstdCompressor(compression_params=params).stream_writer(
                        handle, size=identity[2], closefd=False,
                    )
                with compressor, path.open("rb") as source:
                    for block in iter(lambda: source.read(_CHUNK), b""):
                        digest.update(block)
                        compressor.write(block)
                handle.flush()
                os.fsync(handle.fileno())
            if _identity(path) != identity:
                raise StorageError(f"Source changed during compression: {path}")
            if temporary.stat().st_size >= identity[2] * (1 - min_saving):
                return {"path": str(path), "status": "not_smaller"}
            # Validate decompression before publishing anything or removing source.
            with _decoded(temporary, codec) as stream:
                restored, count = hashlib.sha256(), 0
                for block in iter(lambda: stream.read(_CHUNK), b""):
                    restored.update(block)
                    count += len(block)
            if count != identity[2] or restored.digest() != digest.digest():
                raise StorageError(f"Compression round-trip verification failed: {path}")
            data = {"schema": "chronoclade.compressed-file", "schema_version": 1,
                    "logical_name": path.name, "codec": codec, "mode": mode, "sha256": digest.hexdigest(),
                    "bytes": count, "stored_sha256": _raw_hash(temporary),
                    "stored_bytes": temporary.stat().st_size}
            archive = _archive(path, codec)
            if archive.is_symlink() or receipt_path(path).is_symlink():
                raise StorageError(f"Storage destination is a symlink: {path}")
            os.chmod(temporary, mode)
            os.replace(temporary, archive)
            _write_record(path, data)
            directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
            try:
                os.fsync(directory)
            finally:
                os.close(directory)
            if _identity(path) != identity:
                raise StorageError(f"Source changed before storage commit: {path}")
            path.unlink()
            return {"path": str(path), "status": "packed", "bytes": count,
                    "stored_bytes": data["stored_bytes"], "codec": codec}
        finally:
            temporary.unlink(missing_ok=True)


def _walk(root):
    root = Path(root).absolute()
    if root.is_symlink() or not root.is_dir():
        raise StorageError("Storage requires a real directory; directory symlinks are not traversed")
    for directory, dirs, files in os.walk(root, followlinks=False):
        dirs[:] = [name for name in dirs if not Path(directory, name).is_symlink() and name not in {".git", ".pixi", ".venv"}]
        for name in sorted(files):
            path = Path(directory, name)
            if not path.is_symlink():
                yield path


def reject_running(root):
    for path in _walk(root):
        if path.name == "job.json":
            try:
                state = json.loads(path.read_text())
                if not isinstance(state, dict):
                    raise StorageError(f"Cannot establish job state: {path}")
                if state.get("status") == "running":
                    raise StorageError(f"Refusing to compact a running job: {path}")
            except (OSError, json.JSONDecodeError) as error:
                raise StorageError(f"Cannot establish job state: {path}") from error


def pack(root, *, codec="zstd", min_bytes=_CHUNK, dry_run=False):
    if codec not in {"gzip", "zstd"} or not isinstance(min_bytes, int) or min_bytes < 0:
        raise StorageError("Storage requires gzip/zstd and a nonnegative minimum size")
    reject_running(root)
    candidates = [p for p in _walk(root) if p.is_file() and (
                  p.suffix in _ELIGIBLE or re.fullmatch(
                      r"(?:public_\d+|matched_\d+|input|collection_snapshot)\.json", p.name))
                  and p.name not in _CONTROLS and not p.name.endswith(_SUFFIX)
                  and not p.name.endswith(".fingerprint.json") and p.stat().st_size >= min_bytes]
    if dry_run:
        return [{"path": str(p), "bytes": p.stat().st_size, "status": "candidate"} for p in candidates]
    return [pack_file(path, codec=codec) for path in candidates]


def inspect(root, *, verify=False):
    rows = []
    for receipt in _walk(root):
        if not receipt.name.endswith(_SUFFIX):
            continue
        path = receipt.with_name(receipt.name[:-len(_SUFFIX)])
        data, _ = _record(path)
        if verify:
            # Always check archived bytes, even when a materialized copy exists.
            with _checked(path) as (stream, _):
                digest, count = hashlib.sha256(), 0
                for block in iter(lambda: stream.read(_CHUNK), b""):
                    digest.update(block)
                    count += len(block)
                if digest.hexdigest() != data["sha256"] or count != data["bytes"]:
                    raise StorageError(f"Decompressed artifact checksum/size mismatch: {path}")
            if path.exists() and (not _regular(path) or path.stat().st_size != data["bytes"]
                                  or _raw_hash(path) != data["sha256"]):
                raise StorageError(f"Materialized artifact checksum/size mismatch: {path}")
        rows.append({"path": str(path), **data, "materialized": path.is_file()})
    return rows


def restore(root):
    return [str(materialize(row["path"], permanent=True)) for row in inspect(root, verify=True)]


class StorageSession:
    """Restore native/mmap inputs temporarily and optionally compact new outputs."""
    def __init__(self, policy="plain"):
        if policy not in {"plain", "compact"}:
            raise StorageError("Storage policy must be plain or compact")
        self.policy, self.restored, self.outputs = policy, set(), set()
        self.leases = {}
        self.token = None

    def __enter__(self):
        self.token = _SESSION.set(self)
        return self

    def __exit__(self, exc_type, *_exc):
        # Release our shared leases before taking exclusive packing/eviction
        # locks; other sessions retain theirs until their native tools finish.
        for descriptor in self.leases.values():
            os.close(descriptor)
        self.leases.clear()
        try:
            if exc_type is None and self.policy == "compact":
                for output in sorted(self.outputs):
                    pack(output)
            for path in self.restored:
                with _lock(path):
                    if _regular(path) and receipt_path(path).is_file():
                        data, _ = _record(path)
                        if path.stat().st_size == data["bytes"] and _raw_hash(path) == data["sha256"]:
                            path.unlink()
        finally:
            _SESSION.reset(self.token)


def completed_output(path):
    session = _SESSION.get()
    if session is not None:
        session.outputs.add(Path(path).absolute())


def logical_files(root):
    """Enumerate payloads without exposing compression receipts as science data."""
    seen = set()
    for path in _walk(root):
        if path.name.endswith(_SUFFIX):
            logical = path.with_name(path.name[:-len(_SUFFIX)])
        elif path.name.endswith((".cc.zst", ".cc.gz", ".cc-lock")) or path.name.startswith(".cc-storage-"):
            continue
        else:
            logical = path
        if logical not in seen:
            seen.add(logical)
            yield logical


def current_policy():
    """Freeze the active CLI policy into generated job plans."""
    session = _SESSION.get()
    return session.policy if session is not None else "plain"
