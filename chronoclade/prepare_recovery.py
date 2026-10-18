"""Small fail-closed primitives for resumable local native typing work."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import stat
import uuid
from pathlib import Path

from chronoclade.artifacts import file_sha256, write_json


class RecoveryError(ValueError):
    """A saved recovery job cannot safely be adopted."""


class DirectoryLock:
    """Nonblocking process lock whose kernel ownership ends with the process."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.handle = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if self.path.is_symlink():
            raise RecoveryError(f"Recovery lock path is a symlink: {self.path}")
        flags = os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0)
        try:
            descriptor = os.open(self.path, flags, 0o600)
        except OSError as error:
            raise RecoveryError(f"Cannot safely open recovery lock {self.path}") from error
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            os.close(descriptor)
            raise RecoveryError(f"Recovery lock is not a regular file: {self.path}")
        self.handle = os.fdopen(descriptor, "r+")
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            self.handle.close()
            self.handle = None
            raise RecoveryError(f"Another process owns recovery lock {self.path}") from error
        self.handle.seek(0)
        self.handle.truncate()
        self.handle.write(f"pid={os.getpid()}\n")
        self.handle.flush()
        return self

    def __exit__(self, *_exc):
        if self.handle is not None:
            fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
            self.handle.close()
            self.handle = None


def canonical_digest(value) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def path_identity(value):
    if value is None:
        return None
    path = Path(value).expanduser().resolve()
    if path.is_file():
        return {"path": str(path), "sha256": file_sha256(path), "bytes": path.stat().st_size}
    if path.is_dir():
        files = sorted(item for item in path.rglob("*") if item.is_file())
        return {
            "path": str(path),
            "sha256": canonical_digest(
                [(item.relative_to(path).as_posix(), file_sha256(item)) for item in files]
            ),
            "files": len(files),
        }
    return {"path": str(path), "missing": True}


def load_object(path: Path, what: str) -> dict:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise RecoveryError(f"Cannot read {what} at {path}") from error
    if not isinstance(value, dict):
        raise RecoveryError(f"Invalid {what} at {path}")
    return value


def owned_directory(root: Path, relative: str | Path = ".", *, create: bool = False) -> Path:
    """Validate a directory beneath an owned root without following symlinks."""
    root = Path(root)
    if root.is_symlink():
        raise RecoveryError(f"Owned work root is a symlink: {root}")
    try:
        root_resolved = root.resolve(strict=True)
    except OSError as error:
        raise RecoveryError(f"Owned work root does not exist: {root}") from error
    relative = Path(relative)
    if relative.is_absolute() or ".." in relative.parts:
        raise RecoveryError(f"Recovery path escapes its owned work directory: {relative}")
    current = root
    for part in relative.parts:
        if part == ".":
            continue
        current = current / part
        if current.is_symlink():
            raise RecoveryError(f"Recovery path contains a symlink: {current}")
        if current.exists():
            if not current.is_dir():
                raise RecoveryError(f"Recovery path is not a directory: {current}")
        elif create:
            current.mkdir()
        resolved = current.resolve(strict=False)
        try:
            resolved.relative_to(root_resolved)
        except ValueError as error:
            raise RecoveryError(f"Recovery path escapes its owned work directory: {current}") from error
    return current


def begin_job(work: Path, binding: dict, *, resume: bool) -> dict:
    """Create a new job or validate the exact incomplete job requested for resume."""
    work = Path(work)
    manifest_path = work / "prepare-recovery.json"
    if work.is_symlink():
        raise RecoveryError(f"Recovery work path is a symlink and will not be adopted: {work}")
    if work.exists():
        if not resume:
            raise RecoveryError(
                f"Recovery work already exists at {work}; use --resume or choose a new --out"
            )
        if manifest_path.is_symlink():
            raise RecoveryError(f"Recovery manifest is a symlink and will not be adopted: {manifest_path}")
        manifest = load_object(manifest_path, "prepare recovery manifest")
        if (
            manifest.get("schema") != "chronoclade.prepare-recovery"
            or manifest.get("schema_version") != 1
            or manifest.get("binding") != binding
            or not isinstance(manifest.get("job_id"), str)
            or not manifest.get("job_id")
            or manifest.get("status") != "incomplete"
        ):
            raise RecoveryError(
                f"Recovery work at {work} is corrupt, unrelated, complete, or identity-mismatched; preserved"
            )
        owned_directory(work, "attempts")
        owned_directory(work, "native_typing_cache")
        for child in (work / "attempts", work / "native_typing_cache"):
            if child.exists():
                for descendant in child.rglob("*"):
                    if descendant.is_symlink():
                        raise RecoveryError(f"Recovery work contains a symlink: {descendant}")
        for attempt in manifest.get("attempts", []):
            evidence_directory = attempt.get("evidence_directory") if isinstance(attempt, dict) else None
            if not isinstance(evidence_directory, str):
                raise RecoveryError("Recovery manifest contains an invalid attempt path")
            owned_directory(work, evidence_directory)
        return manifest
    if resume:
        raise RecoveryError(f"No interrupted recovery work exists at {work}; omit --resume")
    work.mkdir(parents=True)
    manifest = {
        "schema": "chronoclade.prepare-recovery",
        "schema_version": 1,
        "status": "incomplete",
        "job_id": uuid.uuid4().hex,
        "binding": binding,
        "attempts": [],
    }
    write_json(manifest_path, manifest)
    return manifest



def bind_resolved_cohort(work: Path, manifest: dict, rows) -> None:
    """Bind the full ordered resolved cohort before any cached-profile adoption."""
    cohort = []
    for row in rows:
        assembly_hash = row.get("assembly_sha256")
        assembly = row.get("assembly")
        if not assembly_hash and isinstance(assembly, str) and Path(assembly).is_file():
            assembly_hash = file_sha256(Path(assembly))
        reference = row.get("assembly_reference") or assembly
        cohort.append({
            "sample_id": row.get("sample_id"),
            "species": " ".join(str(row.get("species") or "").replace("_", " ").casefold().split()),
            "source_genome_id": row.get("source_genome_id") or None,
            "assembly_sha256": assembly_hash,
            "assembly_reference_sha256": canonical_digest(reference) if reference else None,
        })
    cohort_digest = canonical_digest(cohort)
    previous = manifest.get("resolved_cohort")
    if previous is not None and (
        previous.get("sha256") != cohort_digest or previous.get("records") != cohort
    ):
        raise RecoveryError("Resolved cohort identity differs from the interrupted prepare job")
    if previous is None:
        manifest["resolved_cohort"] = {"sha256": cohort_digest, "records": cohort}
        save_job(work, manifest)


def save_job(work: Path, manifest: dict) -> None:
    write_json(Path(work) / "prepare-recovery.json", manifest)
