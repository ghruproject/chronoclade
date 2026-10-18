"""Parse bounded Pathogenwatch HClink/HierCC analysis exports."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import math
import re
from typing import Iterable


class HierCCExportError(ValueError):
    """An HClink response cannot be safely joined to its requested genomes."""


class HierCCExportIdentityError(HierCCExportError):
    """An HClink response cannot be reconciled to its exact requested IDs."""


_MAX_COMPRESSED_BYTES = 25_000_000
_MAX_DECOMPRESSED_BYTES = 25_000_000
_MAX_ROWS = 1_000
_JOB = re.compile(r"[A-Za-z0-9_-]+")
_LEVEL = re.compile(r"HC([0-9]+)")


def _decode(data: bytes) -> bytes:
    if len(data) > _MAX_COMPRESSED_BYTES:
        raise HierCCExportError("HClink response exceeds compressed size limit")
    if data.startswith(b"\x1f\x8b"):
        try:
            with gzip.GzipFile(fileobj=io.BytesIO(data)) as stream:
                decoded = stream.read(_MAX_DECOMPRESSED_BYTES + 1)
        except (OSError, EOFError) as error:
            raise HierCCExportError("HClink response is corrupt gzip") from error
        if len(decoded) > _MAX_DECOMPRESSED_BYTES:
            raise HierCCExportError("HClink response exceeds decompressed size limit")
        return decoded
    if len(data) > _MAX_DECOMPRESSED_BYTES:
        raise HierCCExportError("HClink response exceeds decompressed size limit")
    return data


def _validated_rows(reader):
    try:
        yield from reader
    except csv.Error as error:
        raise HierCCExportError("HClink response contains malformed CSV") from error


def parse_hiercc_export(
    data: bytes, requested: Iterable[str], *, job: str
) -> list[dict]:
    """Normalize one exact HClink batch without inferring reference versions.

    Blank HierCC cells remain absent from ``hiercc_codes`` and are separately
    retained in ``hiercc_blank_levels``/``hiercc_raw_levels``. Numeric cluster
    zero is a valid assignment and remains ``"0"``.
    """
    if not isinstance(job, str) or not _JOB.fullmatch(job):
        raise HierCCExportError("HClink job identifier is malformed")
    requested_ids = list(requested)
    if not requested_ids or any(not isinstance(value, str) or not value for value in requested_ids):
        raise HierCCExportError("HClink export requires exact nonempty requested IDs")
    if len(requested_ids) != len(set(requested_ids)):
        raise HierCCExportError("HClink request repeats a source genome ID")

    decoded = _decode(data)
    try:
        reader = csv.DictReader(
            io.StringIO(decoded.decode("utf-8-sig"), newline=""), strict=True
        )
        headers = reader.fieldnames or []
    except UnicodeDecodeError as error:
        raise HierCCExportError("HClink response is not UTF-8 CSV") from error
    except csv.Error as error:
        raise HierCCExportError("HClink response contains malformed CSV header") from error
    if len(headers) != len(set(headers)):
        raise HierCCExportError("HClink export repeats a column header")
    required = {"Genome ID", "Genome Name", "Closest ST", "Distance to closest ST"}
    if not required.issubset(headers):
        raise HierCCExportError("HClink export is missing required identity/closest-ST columns")
    levels: dict[str, str] = {}
    for header in headers:
        match = _LEVEL.fullmatch(header)
        if match:
            level = "HC" + str(int(match[1]))
            if level in levels:
                raise HierCCExportError("HClink export repeats a normalized HierCC level")
            levels[level] = header
    if not levels:
        raise HierCCExportError("HClink export contains no recognized HierCC levels")

    expected = set(requested_ids)
    records: dict[str, dict] = {}
    for index, item in enumerate(_validated_rows(reader), start=1):
        if index > _MAX_ROWS:
            raise HierCCExportError("HClink export exceeds row limit")
        if None in item:
            raise HierCCExportError("HClink row contains extra unnamed CSV fields")
        if any(item.get(column) is None for column in required):
            raise HierCCExportError("HClink row is missing required column values")
        ident = (item.get("Genome ID") or "").strip()
        if ident not in expected:
            raise HierCCExportIdentityError("HClink export contains an unrequested genome ID")
        if ident in records:
            raise HierCCExportIdentityError("HClink export repeats a genome ID")
        raw_levels = {}
        codes = {}
        blank_levels = []
        for level, header in levels.items():
            raw = item.get(header)
            if raw is None:
                raise HierCCExportError("HClink row is missing a declared level cell")
            raw = raw.strip()
            raw_levels[level] = raw
            if not raw:
                blank_levels.append(level)
                continue
            if not raw.isascii() or not raw.isdigit():
                raise HierCCExportError("HClink cluster IDs must be nonnegative integers or blank")
            codes[level] = str(int(raw))

        closest_raw = (item.get("Closest ST") or "").strip()
        distance_raw = (item.get("Distance to closest ST") or "").strip()
        distance = None
        if distance_raw:
            try:
                distance = float(distance_raw)
            except ValueError as error:
                raise HierCCExportError("HClink closest-ST distance is not numeric") from error
            if not math.isfinite(distance) or distance < 0:
                raise HierCCExportError("HClink closest-ST distance must be finite and nonnegative")

        records[ident] = {
            "source_genome_id": ident,
            "hiercc_scheme": "pathogenwatch:HierCC",
            "hiercc_scheme_version": None,
            "hiercc_provider_job": job,
            "hiercc_database_version": None,
            "hiercc_database_sha256": None,
            "hiercc_database_version_status": "unavailable_not_returned_by_provider",
            "hiercc_source": "pathogenwatch_hclink_analysis_export",
            "hiercc_codes": codes,
            "hiercc_status": "resolved" if codes else "unassigned",
            "hiercc_genome_name": (item.get("Genome Name") or "").strip(),
            "hiercc_closest_st": closest_raw or None,
            "hiercc_distance_to_closest_st": distance,
            "hiercc_distance_to_closest_st_raw": distance_raw or None,
            "hiercc_raw_levels": raw_levels,
            "hiercc_blank_levels": blank_levels,
            "hiercc_export_sha256": hashlib.sha256(data).hexdigest(),
        }
    missing = expected - records.keys()
    if missing:
        raise HierCCExportIdentityError("HClink export omits requested genome IDs")
    return [records[ident] for ident in requested_ids]
