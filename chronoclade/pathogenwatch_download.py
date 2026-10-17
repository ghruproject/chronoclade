"""Bounded, authenticated assembly acquisition with verified content caches."""

from __future__ import annotations

import gzip
import fcntl
import hashlib
import io
import json
import re
import time
import uuid
import zipfile
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin, urlparse
from urllib.request import HTTPRedirectHandler, Request, build_opener

from chronoclade.file_snapshots import copy_evidence_snapshot
from chronoclade.pathogenwatch_identity import PublicIdentityError, resolve_public_source_identity


class DownloadError(RuntimeError):
    """A download failed validation or exhausted its retries."""


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


@contextmanager
def _source_cache_lock(cache_dir: Path, source: str):
    """Serialize cache lookup and production for one source across threads/processes."""
    lock_path = cache_dir / f".{source}.lock"
    with lock_path.open("a+b") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def request_download(
    route: str,
    *,
    api_key: str,
    body: dict | None = None,
    base_url: str = "https://pathogen.watch",
    attempts: int = 3,
    timeout: float = 120,
    stats: dict | None = None,
) -> bytes:
    """Follow the documented file redirect without forwarding authentication."""
    if stats is None:
        stats = {}
    for attempt in range(attempts):
        try:
            stats["requests"] = stats.get("requests", 0) + 1
            request = Request(
                base_url + route,
                data=json.dumps(body).encode() if body else None,
                headers={"X-API-Key": api_key, "Content-Type": "application/json"},
            )
            try:
                with build_opener(_NoRedirect()).open(request, timeout=timeout) as response:
                    return response.read()
            except HTTPError as error:
                if error.code not in {301, 302, 303, 307, 308}:
                    raise
                target = urljoin(base_url, error.headers.get("Location", ""))
                parsed = urlparse(target)
                if parsed.scheme != "https" or not parsed.hostname or parsed.username:
                    raise DownloadError("Download redirect is not a valid HTTPS file URL")
                # Signed URLs must never enter logs or manifests. No auth on file request.
                stats["requests"] = stats.get("requests", 0) + 1
                with build_opener().open(Request(target), timeout=timeout) as response:
                    return response.read()
        except HTTPError as error:
            if error.code not in {429, 500, 502, 503, 504}:
                raise DownloadError(f"Pathogenwatch download HTTP {error.code}") from None
        except (URLError, TimeoutError, OSError):
            pass
        if attempt + 1 < attempts:
            time.sleep(min(2**attempt, 8))
    raise DownloadError("Pathogenwatch download exhausted retries")


def validate_fasta(content: bytes) -> bytes:
    """Reject empty, malformed and non-nucleotide assembly files."""
    if content.startswith(b"\x1f\x8b"):
        try:
            content = gzip.decompress(content)
        except (OSError, EOFError) as error:
            raise DownloadError("Corrupt compressed FASTA") from error
    try:
        lines = content.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise DownloadError("FASTA is not ASCII") from error
    if not lines or not lines[0].startswith(">"):
        raise DownloadError("Downloaded content is not FASTA")
    length = 0
    current = 0
    seen_header = False
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seen_header and not current:
                raise DownloadError("FASTA contains an empty contig")
            current = 0
            seen_header = True
        elif re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", line):
            current += len(line)
            length += len(line)
        else:
            raise DownloadError("FASTA contains invalid sequence characters")
    if not length or not current:
        raise DownloadError("FASTA contains no sequence")
    return content


def unpack_bulk(content: bytes, id_to_source: dict[str, str]) -> dict[str, bytes]:
    """Map archive members only by exact IDs, never archive order or fuzzy names."""
    if content.startswith(b"\x1f\x8b"):
        content = gzip.decompress(content)
    try:
        archive = zipfile.ZipFile(io.BytesIO(content))
    except zipfile.BadZipFile as error:
        raise DownloadError("Bulk download is not a ZIP archive") from error
    result: dict[str, bytes] = {}
    for member in archive.infolist():
        if member.is_dir():
            continue
        name = Path(member.filename).name
        stem = re.sub(r"\.(?:fasta|fna|fa)(?:\.gz)?$", "", name, flags=re.I)
        source = id_to_source.get(stem)
        if source is None:
            raise DownloadError("Bulk archive member cannot be joined to a requested ID")
        if source in result:
            raise DownloadError("Bulk archive contains a duplicate genome ID")
        if member.file_size > 100_000_000:
            raise DownloadError("Bulk archive member exceeds the assembly size limit")
        result[source] = validate_fasta(archive.read(member))
    return result


def _atomic(path: Path, content: bytes) -> None:
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".partial")
    temporary.write_bytes(content)
    temporary.replace(path)


def download_assemblies(
    rows: list[dict],
    *,
    output: Path,
    cache_dir: Path,
    api_key: str,
    workers: int = 2,
    base_url: str = "https://pathogen.watch",
    fetch=None,
    identity_fetch=None,
    require_public_identity: bool = False,
) -> tuple[dict[str, Path], list[dict]]:
    """Use public-identity and content checks before accepting any cached assembly."""
    output.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    def one(row):
        # Identity lookups enrich only this worker's verification state. The
        # canonical prepared metadata remains exactly as supplied by its owner.
        verification_row = dict(row)
        source = str(row["source_genome_id"])
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", source):
            return (
                source,
                None,
                {"source_genome_id": source, "status": "failed", "reason": "unsafe_source_id"},
            )
        metadata = cache_dir / f"{source}.json"
        started = time.monotonic()
        cached = False
        requests = 0
        stats = {}
        try:
            with _source_cache_lock(cache_dir, source):
                expected = str(verification_row.get("source_checksum") or "")
                checksum_pattern = r"[a-fA-F0-9]{40}"
                identity_verified = bool(re.fullmatch(checksum_pattern, expected))
                identity_metadata = {}
                if require_public_identity and not identity_verified:
                    if identity_fetch is None:
                        from chronoclade.pathogenwatch import PathogenwatchClient

                        identity_client = PathogenwatchClient(base_url=base_url, api_key="")
                        def details(source_id):
                            return identity_client.request_json(
                                "GET",
                                "/api/genomes/details",
                                params={"id": source_id},
                                authenticated=False,
                            )
                    else:
                        details = identity_fetch
                    try:
                        detail = details(source)
                    except Exception as error:
                        raise PublicIdentityError("public_identity_unavailable") from error
                    identity_metadata = resolve_public_source_identity(
                        source, verification_row, detail
                    )
                    expected = str(identity_metadata.get("source_checksum") or "")
                    identity_verified = bool(re.fullmatch(checksum_pattern, expected))
                    requests += 1
                    if not identity_verified:
                        raise PublicIdentityError("public_identity_checksum_unavailable")
                data = None
                saved = None
                if metadata.is_file():
                    try:
                        saved = json.loads(metadata.read_text())
                        blob = cache_dir / (saved["sha256"] + ".fasta")
                        existing = blob.read_bytes()
                        if (
                            saved["source_genome_id"] == source
                            and hashlib.sha256(existing).hexdigest() == saved["sha256"]
                            and (
                                hashlib.sha1(existing).hexdigest() == expected.lower()
                                if identity_verified
                                else saved.get("source_checksum") == verification_row.get("source_checksum")
                            )
                        ):
                            data = validate_fasta(existing)
                            cached = True
                    except (ValueError, KeyError, OSError, DownloadError):
                        pass
                if data is None:
                    stats = {}
                    if fetch:
                        requests = 1
                        data = validate_fasta(fetch(source))
                    else:
                        data = validate_fasta(
                            request_download(
                                f"/api/genomes/download/{source}",
                                api_key=api_key,
                                base_url=base_url,
                                stats=stats,
                            )
                        )
                        requests = stats.get("requests", 0)
                if identity_verified and hashlib.sha1(data).hexdigest() != expected.lower():
                    raise DownloadError("source_checksum_mismatch")
                sequence_bases = sum(
                    len(line.strip()) for line in data.splitlines() if not line.startswith(b">")
                )
                source_length = verification_row.get("source_length")
                if source_length is None:
                    source_length = identity_metadata.get("source_length")
                if source_length is not None and sequence_bases != int(source_length):
                    raise DownloadError("source_length_mismatch")
                digest = hashlib.sha256(data).hexdigest()
                blob = cache_dir / f"{digest}.fasta"
                cache_record = {
                    "source_genome_id": source,
                    "sha256": digest,
                    "source_checksum": expected.lower() if identity_verified else verification_row.get("source_checksum"),
                    "identity_verified": identity_verified,
                    "source_length": source_length,
                    "raw_details_sha256": identity_metadata.get("raw_details_sha256")
                    or verification_row.get("raw_details_sha256"),
                }
                if not cached:
                    _atomic(blob, data)
                if not cached or (
                    identity_verified and saved.get("source_checksum") != expected.lower()
                ):
                    _atomic(metadata, json.dumps(cache_record).encode())
                destination = output / f"{source}.fasta"
                partial = destination.with_suffix(".fasta.partial")
                partial.unlink(missing_ok=True)
                copy_method = copy_evidence_snapshot(blob, partial)
                partial.replace(destination)
                return (
                    source,
                    destination.resolve(),
                    {
                        "source_genome_id": source,
                        "status": "cached" if cached else "downloaded",
                        "sha256": digest,
                        "bytes": len(data),
                        "sequence_bases": sequence_bases,
                        "source_checksum_verified": identity_verified,
                        "source_identity_verified": identity_verified,
                        "source_checksum": expected.lower() if identity_verified else None,
                        "source_length": source_length,
                        "raw_details_sha256": identity_metadata.get("raw_details_sha256"),
                        "copy_method": copy_method,
                        "requests": requests,
                        "seconds": round(time.monotonic() - started, 3),
                    },
                )
        except (DownloadError, PublicIdentityError, OSError, ValueError) as error:
            requests = stats.get("requests", requests)
            return (
                source,
                None,
                {
                    "source_genome_id": source,
                    "status": "failed",
                    "reason": str(error)
                    if str(error) in {
                        "source_checksum_mismatch",
                        "source_length_mismatch",
                        "public_identity_unverified",
                        "public_identity_scope_mismatch",
                        "public_identity_checksum_unavailable",
                        "public_identity_unavailable",
                    }
                    else "download_or_validation_failed",
                    "requests": requests,
                    "seconds": round(time.monotonic() - started, 3),
                },
            )

    with ThreadPoolExecutor(max_workers=max(1, min(workers, 4))) as executor:
        results = list(executor.map(one, rows))
    return (
        {source: path for source, path, _ in results if path is not None},
        [ledger for _, _, ledger in results],
    )
