"""Resolve focal genomes and frozen typing before acquiring context assemblies.

Collection membership is a query dataset. Public discovery always omits credentials;
only an explicitly requested collection and analysis/assembly downloads may use them.
"""

from __future__ import annotations

from copy import deepcopy
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import time
from typing import Callable, Iterable
from urllib.parse import urlencode, urlsplit
import zipfile

from chronoclade.cglin import annotate_catalogue, download_cglin_export
from chronoclade.export_checkpoints import (
    ExportCheckpointError,
    atomic_json,
    batch_contract,
    checkpoint_path,
    load_checkpoint,
    persist_checkpoint,
    row_scope,
)
from chronoclade.metadata import SAFE_IDENTIFIER, Sample
from chronoclade.pathogenwatch import (
    PathogenwatchClient,
    PathogenwatchError,
    content_hash,
    normalize_country,
    normalize_dates,
)
from chronoclade.pathogenwatch_download import DownloadError, download_assemblies, request_download
from chronoclade.progress import progress_task


class ProfileInputError(ValueError):
    """Input identities or frozen profile data cannot be resolved safely."""


class CGMLSTExportIdentityError(ProfileInputError):
    """A response contains a genome outside the exact requested batch."""


def collection_id(value: str) -> str:
    """Accept a short UUID or a Pathogenwatch collection URL with an optional slug."""
    text = value.strip()
    if "://" in text:
        parsed = urlsplit(text)
        if parsed.scheme != "https" or parsed.hostname not in {
            "pathogen.watch",
            "next.pathogen.watch",
        }:
            raise ProfileInputError("Collection URL must use HTTPS on pathogen.watch")
        if parsed.username or parsed.password or parsed.port not in {None, 443}:
            raise ProfileInputError("Collection URL cannot contain credentials or a custom port")
        match = re.fullmatch(r"/collections/([^/]+)/?", parsed.path)
        if not match:
            raise ProfileInputError("Expected a Pathogenwatch /collections/UUID URL")
        text = match[1]
    ident = text.split("-", 1)[0]
    if not re.fullmatch(r"[1-9A-HJ-NP-Za-km-z]{20,22}", ident):
        raise ProfileInputError("Collection identifier must be a short UUID")
    return ident


def _csv(path: Path) -> list[dict]:
    try:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames:
                raise ProfileInputError("Metadata requires a CSV header")
            rows = [
                {key: (value or "").strip() for key, value in row.items() if key} for row in reader
            ]
    except OSError as error:
        raise ProfileInputError(f"Cannot read metadata: {path}") from error
    if not rows:
        raise ProfileInputError("Metadata has no samples")
    return rows


def _canonical(raw: dict, *, origin: str, species: str = "", base: Path | None = None) -> dict:
    row = deepcopy(raw)
    flat = {
        key.rsplit("/", 1)[-1].casefold(): value
        for key, value in raw.items()
        if key.startswith("metadata/")
    }
    metadata = raw.get("metadata") or {}
    if isinstance(metadata, list):
        metadata = dict(metadata)
    flat.update({str(key).casefold(): value for key, value in metadata.items()})
    source = str(raw.get("source_genome_id") or raw.get("uuid") or "")
    sample = str(raw.get("sample_id") or ("PW_" + source if source else raw.get("name") or ""))
    if not SAFE_IDENTIFIER.fullmatch(sample):
        raise ProfileInputError("Each sample requires a unique safe sample_id or source UUID")
    st = str(raw.get("mlst_st") or raw.get("mlst") or raw.get("typing/MLST/ST") or "")
    lineage = str(raw.get("lineage") or ("ST" + st if st and st != "-" else "Unassigned"))
    country_raw = (
        raw.get("country")
        or flat.get("country")
        or raw.get("location")
        or flat.get("location")
        or ""
    )
    region = raw.get("region") or flat.get("city/region") or flat.get("region") or ""
    city = raw.get("city") or flat.get("city") or ""
    if str(region).strip().casefold() in {"unknown", "none", "null", "na", "n/a", "not provided"}:
        region = ""
    if str(city).strip().casefold() in {"unknown", "none", "null", "na", "n/a", "not provided"}:
        city = ""
    date = (
        raw.get("collection_date")
        or flat.get("collection date")
        or flat.get("collection_date")
        or flat.get("date")
        or flat.get("year")
        or ""
    )
    row.update(
        sample_id=sample,
        origin=origin,
        species=str(raw.get("species") or species),
        lineage=lineage,
        mlst_st=st,
        source_genome_id=source,
        numeric_source_id=raw.get("numeric_source_id", raw.get("id")),
        country=normalize_country(country_raw),
        country_raw=raw.get("country_raw", country_raw),
        region=str(region),
        city=str(city),
        nuts2=str(raw.get("nuts2") or flat.get("nuts2") or ""),
        location=str(raw.get("location") or country_raw or "Unknown"),
        geography_raw=raw.get(
            "geography_raw", {"country": country_raw, "region": region, "city": city}
        ),
        aliases=list(raw.get("aliases") or []),
    )
    for target, raw_field in (
        ("source_checksum", "source_checksum"),
        ("source_length", "source_length"),
        ("source_n50", "source_n50"),
        ("source_contigs", "source_contigs"),
        ("qc_pass", "qc_pass"),
    ):
        value = raw.get(raw_field)
        if value is not None:
            row[target] = value
    if row.get("source_checksum") is None and raw.get("checksum") is not None:
        row["source_checksum"] = raw["checksum"]
    if row.get("source_length") is None and raw.get("length") is not None:
        row["source_length"] = raw["length"]
    if row.get("qc_pass") is None and raw.get("qc") is not None:
        row["qc_pass"] = raw["qc"]
    row.update(normalize_dates(raw, {"collection_date": date}, raw))
    # Empty or partial dates remain explicit, never replaced by upload dates.
    if raw.get("date_precision"):
        row.update(
            {
                key: raw[key]
                for key in (
                    "collection_date",
                    "date_precision",
                    "date_start",
                    "date_end",
                    "date_raw",
                )
                if key in raw
            }
        )
    for field in ("sampleAccession", "runAccession", "assemblyAccession", "biosample", "accession"):
        if raw.get(field):
            row["aliases"].append(str(raw[field]))
    for key, value in raw.items():
        if (
            key.startswith("metadata/")
            and any(
                x in key.lower()
                for x in ("run accession", "sample accession", "assembly accession")
            )
            and value
        ):
            row["aliases"].append(str(value))
    row["aliases"] += [sample, source, str(raw.get("name") or "")]
    row["biosample_accessions"] = sorted(
        set(row.get("biosample_accessions") or [])
        | {a for a in row["aliases"] if a.startswith("SAM")}
    )
    row["run_accessions"] = sorted(
        set(row.get("run_accessions") or [])
        | {a for a in row["aliases"] if re.fullmatch(r"[SED]RR[0-9]+", a)}
    )
    row["assembly_accessions"] = sorted(
        set(row.get("assembly_accessions") or [])
        | {a for a in row["aliases"] if a.startswith(("GCA_", "GCF_"))}
    )
    row["aliases"] = sorted(set(value for value in row["aliases"] if value))
    for source_field, collection_field in (
        ("source_length", "assemblyStats/Genome length"),
        ("source_n50", "assemblyStats/N50"),
        ("source_contigs", "assemblyStats/No. contigs"),
    ):
        if collection_field in raw:
            row[source_field] = raw[collection_field]
    assembly = raw.get("assembly")
    if (
        isinstance(assembly, str)
        and assembly
        and (base is not None or Path(assembly).is_absolute())
    ):
        path = Path(assembly).expanduser()
        path = path if path.is_absolute() else base / path
        if not path.is_file():
            raise ProfileInputError(f"Assembly does not exist for {sample}: {path}")
        row["assembly"] = str(path.resolve())
        row["assembly_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    else:
        row.pop("assembly", None)
    return row


def _collection(client: PathogenwatchClient, ident: str, output: Path) -> tuple[list[dict], dict]:
    details = client.request_json(
        "GET", "/api/collections/details", params={"uuid": ident}, authenticated=True
    )
    if not isinstance(details, dict) or details.get("uuid") != ident:
        raise ProfileInputError("Collection details do not match the requested UUID")
    size = details.get("size")
    if type(size) is not int or size < 1:
        raise ProfileInputError("Collection must advertise a positive genome count")
    rows, pages, seen = [], [], set()
    for page_number in range(1, size + 2):
        page = client.request_json(
            "GET",
            "/api/collections/genomes",
            params={"uuid": ident, "page": page_number},
            authenticated=True,
        )
        if (
            not isinstance(page, dict)
            or not isinstance(page.get("genomes"), list)
            or type(page.get("hasMore")) is not bool
        ):
            raise ProfileInputError("Unsupported collection pagination response")
        pages.append(page)
        if page["hasMore"] and not page["genomes"]:
            raise ProfileInputError("Collection pagination stalled")
        for member in page["genomes"]:
            uuid = member.get("uuid")
            if not uuid or uuid in seen:
                raise ProfileInputError("Collection repeats or omits a source UUID")
            seen.add(uuid)
            row = _canonical(member, origin="local", species=details.get("organismName", ""))
            row["organism_id"] = str(details.get("organismId") or "")
            rows.append(row)
        if not page["hasMore"]:
            break
    if len(rows) != size:
        raise ProfileInputError("Collection membership does not reconcile with its advertised size")
    snapshot = {"schema_version": 1, "details": details, "pages": pages}
    snapshot["snapshot_sha256"] = content_hash(snapshot)
    path = output / "collection_snapshot.json"
    path.write_text(json.dumps(snapshot, indent=2) + "\n")
    return rows, {
        "collection_uuid": ident,
        "collection_access": details.get("access", "unknown"),
        "collection_snapshot": str(path.resolve()),
        "collection_snapshot_sha256": snapshot["snapshot_sha256"],
        "downloads": details.get("downloads") or [],
        "member_count": len(rows),
    }


def parse_cgmlst_export(data: bytes, requested: Iterable[str], *, job: str) -> list[dict]:
    """Read the verified long-format allele export; never invent a full locus universe."""
    if data.startswith(b"\x1f\x8b"):
        data = gzip.decompress(data)
    if zipfile.is_zipfile(io.BytesIO(data)):
        archive = zipfile.ZipFile(io.BytesIO(data))
        files = [member for member in archive.infolist() if not member.is_dir()]
        if len(files) != 1 or files[0].file_size > 200_000_000:
            raise ProfileInputError("Expected one bounded cgMLST CSV export")
        data = archive.read(files[0])
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
    if not {"Genome ID", "Gene", "Allele ID"}.issubset(reader.fieldnames or []):
        raise ProfileInputError("Unsupported cgMLST export columns")
    expected, profiles, ambiguous = set(requested), {}, {}
    for item in reader:
        ident, locus, allele = item["Genome ID"], item["Gene"], item["Allele ID"].strip()
        if ident not in expected:
            raise CGMLSTExportIdentityError("cgMLST export contains an unrequested genome")
        if not locus:
            raise ProfileInputError("cgMLST export contains an empty locus")
        calls = profiles.setdefault(ident, {})
        missing = allele in {"", "0", "-", "?", "None", "null"}
        if not missing and not (
            allele.isdigit() and int(allele) > 0 or re.fullmatch(r"[a-fA-F0-9]{40}", allele)
        ):
            raise ProfileInputError("cgMLST export has an unsupported allele value")
        allele = "" if missing else str(int(allele)) if allele.isdigit() else allele.lower()
        if locus in calls and calls[locus] != allele:
            ambiguous.setdefault(ident, set()).add(locus)
        calls[locus] = allele
    result = []
    for ident, calls in profiles.items():
        conflicts = ambiguous.get(ident, set())
        result.append(
            {
                "source_genome_id": ident,
                "cgmlst_scheme": "pathogenwatch:" + job,
                "cgmlst_scheme_version": job,
                "cgmlst_loci": sorted(calls),
                "cgmlst_locus_universe_complete": False,
                "cgmlst_profile": {
                    key: value for key, value in calls.items() if value and key not in conflicts
                },
                "cgmlst_novel_alleles": {},
                "cgmlst_ambiguous_loci": sorted(conflicts),
                "cgmlst_source": "pathogenwatch_analysis_export",
                "cgmlst_status": "resolved" if any(calls.values()) else "unassigned",
                "cgmlst_export_sha256": hashlib.sha256(data).hexdigest(),
            }
        )
    return result


def _analysis_export_eligible(row: dict, name: str) -> bool:
    """Share exact-query eligibility rules between capability and download stages."""
    if not str(row.get("numeric_source_id", "")).isdigit():
        return False
    if name == "cgmlst":
        return not row.get("cgmlst_profile") and not row.get("cgmlst_novel_alleles")
    if name == "klebsiella-lincodes":
        return "klebsiella" in row.get("species", "").casefold() and not row.get("cglin_raw")
    if name == "hclink":
        return (
            row.get("species", "").casefold() == "escherichia coli"
            and str(row.get("organism_id") or row.get("organismId") or "") == "562"
            and not row.get("hiercc_codes")
            and not row.get("hiercc_export_sha256")
        )
    return False


def _download_parse_batches(
    rows: list[dict],
    *,
    name: str,
    job: str,
    output: Path,
    client: PathogenwatchClient,
    provenance: dict,
    parser: Callable[..., list[dict]],
    identity_error: type[Exception],
    receipts: list[dict],
    resume: bool = False,
) -> list[dict]:
    """Download/reconcile bounded exact-ID profile or HClink batches."""
    assignments = []
    endpoint = f"/api/downloads/{name}"
    route = endpoint + "?" + urlencode({"job": job})
    batch_count = (len(rows) + 99) // 100
    for batch_number, start in enumerate(range(0, len(rows), 100), start=1):
        batch = rows[start : start + 100]
        source_ids = [row["source_genome_id"] for row in batch]
        numeric_ids = [str(row["numeric_source_id"]) for row in batch]
        contract = batch_contract(batch, name=name, job=job, base_url=client.base_url)
        checkpoint = checkpoint_path(output, name, job, start)
        if checkpoint.exists():
            if not resume:
                raise ExportCheckpointError(
                    "Raw export checkpoint already exists; use resume or a fresh output scope"
                )
            accepted, parsed, _ = load_checkpoint(
                checkpoint, root=output.parent, parser=parser, expected=contract
            )
            receipts.append(accepted)
            provenance.setdefault("analysis_export_raw_reuse", []).append(
                {
                    "checkpoint": str(checkpoint.resolve()),
                    "job": job,
                    "response_sha256": accepted["response_sha256"],
                    "network_request": False,
                }
            )
            assignments.extend(parsed)
            continue
        if checkpoint.parent.exists():
            raise ExportCheckpointError(
                "Incomplete committed batch directory has no verified receipt; use a fresh output scope"
            )
        if list(output.glob(f"{'cgmlst' if name == 'cgmlst' else 'hclink'}_{job}_{start}*.bin")):
            raise ExportCheckpointError(
                "Raw export exists without accepted request receipt; use a fresh output scope"
            )
        for attempt in range(3):
            receipt = {
                "download": name,
                "request_evidence": "recorded_exact_http_request",
                "job": job,
                "base_url": client.base_url,
                "endpoint": endpoint,
                "route": route,
                "source_genome_ids": source_ids,
                "numeric_request_ids": numeric_ids,
                "retry_attempt": attempt + 1,
                "returned_source_genome_ids": None,
                "missing_source_genome_ids": None,
            }
            try:
                with progress_task(
                    "prepare",
                    f"{name} export {job} batch {batch_number}/{batch_count} "
                    f"({len(batch)} genomes; cumulative {start + len(batch)}/{len(rows)})",
                ):
                    data = request_download(
                        route,
                        api_key=client.api_key,
                        base_url=client.base_url,
                        body={"ids": ",".join(numeric_ids)},
                    )
            except DownloadError as error:
                receipt.update(status="download_error", failure_reason=str(error))
                receipts.append(receipt)
                raise
            suffix = f"_retry{attempt}" if attempt else ""
            prefix = "cgmlst" if name == "cgmlst" else "hclink"
            destination = checkpoint.parent / "response.raw"
            receipt.update(
                saved_response_path=destination.relative_to(output.parent).as_posix(),
                response_sha256=hashlib.sha256(data).hexdigest(),
                response_bytes=len(data),
            )
            try:
                parsed = parser(data, source_ids, job=job)
            except identity_error as error:
                rejected = output / ".rejected-exports"
                rejected.mkdir(exist_ok=True)
                destination = rejected / f"{prefix}_{job}_{start}{suffix}.raw"
                destination.write_bytes(data)
                receipt["saved_response_path"] = destination.relative_to(output.parent).as_posix()
                receipt.update(
                    status="rejected_identity",
                    failure_reason=str(error),
                    missing_source_genome_ids=source_ids,
                )
                receipts.append(receipt)
                provenance.setdefault("analysis_export_identity_retries", []).append(
                    {
                        "job": job,
                        "attempt": attempt + 1,
                        "requested_source_ids": source_ids,
                        "rejected_export": str(destination.resolve()),
                        "rejected_export_sha256": receipt["response_sha256"],
                    }
                )
                if attempt == 2:
                    raise
                time.sleep(2**attempt)
                continue
            except ValueError as error:
                rejected = output / ".rejected-exports"
                rejected.mkdir(exist_ok=True)
                destination = rejected / f"{prefix}_{job}_{start}{suffix}.raw"
                destination.write_bytes(data)
                receipt["saved_response_path"] = destination.relative_to(output.parent).as_posix()
                receipt.update(status="parse_error", failure_reason=str(error))
                receipts.append(receipt)
                raise
            returned = list(
                dict.fromkeys(
                    str(row.get("source_genome_id"))
                    for row in parsed
                    if row.get("source_genome_id") is not None
                )
            )
            returned_set = set(returned)
            receipt["returned_source_genome_ids"] = returned
            receipt["missing_source_genome_ids"] = [
                ident for ident in source_ids if ident not in returned_set
            ]
            receipt["status"] = "accepted"
            persist_checkpoint(checkpoint, contract, receipt, data)
            receipts.append(receipt)
            assignments.extend(parsed)
            break
    return assignments


def _analysis_exports(
    rows: list[dict],
    downloads: list[dict],
    client: PathogenwatchClient,
    output: Path,
    provenance: dict,
    *,
    resume: bool = False,
) -> list[dict]:
    if not client.api_key:
        provenance["analysis_export_status"] = "unavailable_no_credentials"
        return rows
    output.mkdir(parents=True, exist_ok=True)
    for advertised in downloads:
        name, job = advertised.get("name", ""), advertised.get("job", "")
        if name not in {"cgmlst", "klebsiella-lincodes", "hclink"} or not re.fullmatch(
            r"[A-Za-z0-9_-]+", job
        ):
            continue
        eligible = [row for row in rows if _analysis_export_eligible(row, name)]
        if not eligible:
            continue
        receipts: list[dict] = []
        try:
            if name == "klebsiella-lincodes":
                manifest = download_cglin_export(
                    eligible,
                    output / ("exports_" + job),
                    api_key=client.api_key,
                    job=job,
                    scheme_version=job,
                    base_url=client.base_url,
                )
                rows = annotate_catalogue(rows, manifest["assignments"])
                for row in rows:
                    if row.get("cglin_scheme_version") == job:
                        row["cglin_scope_source"] = "server_advertised_analysis_job"
                        row["cglin_database_version_status"] = "unavailable"
            elif name == "cgmlst":
                assignments = _download_parse_batches(
                    eligible,
                    name=name,
                    job=job,
                    output=output,
                    client=client,
                    provenance=provenance,
                    parser=parse_cgmlst_export,
                    identity_error=CGMLSTExportIdentityError,
                    receipts=receipts,
                    resume=resume,
                )
                # No complete universe/call fraction is inferred from the observed export.
                by_id = {row["source_genome_id"]: row for row in assignments}
                rows = [
                    dict(
                        row,
                        **{
                            key: value
                            for key, value in by_id.get(row["source_genome_id"], {}).items()
                            if key.startswith("cgmlst_")
                        },
                    )
                    for row in rows
                ]
            else:
                from chronoclade.pathogenwatch_hiercc import (
                    HierCCExportIdentityError,
                    parse_hiercc_export,
                )

                assignments = _download_parse_batches(
                    eligible,
                    name=name,
                    job=job,
                    output=output,
                    client=client,
                    provenance=provenance,
                    parser=parse_hiercc_export,
                    identity_error=HierCCExportIdentityError,
                    receipts=receipts,
                    resume=resume,
                )
                by_id = {row["source_genome_id"]: row for row in assignments}
                rows = [
                    dict(
                        row,
                        **{
                            key: value
                            for key, value in by_id.get(row["source_genome_id"], {}).items()
                            if key.startswith("hiercc_")
                        },
                    )
                    for row in rows
                ]
            provenance.setdefault("analysis_exports", []).append(
                {
                    "download": name,
                    "job": job,
                    "status": "downloaded",
                    "requested_count": len(eligible),
                    "batch_receipts": receipts,
                    "scope_provenance": (
                        "server-advertised HClink assignment job; database release, version and fingerprint unavailable"
                        if name == "hclink"
                        else "server-advertised analysis job; typing database version and fingerprint unavailable"
                    ),
                }
            )
        except (PathogenwatchError, DownloadError, ProfileInputError, ValueError) as error:
            provenance.setdefault("analysis_exports", []).append(
                {
                    "download": name,
                    "job": job,
                    "status": "unavailable",
                    "reason": str(error),
                    "batch_receipts": receipts,
                }
            )
    return rows


def _grouped_analysis_exports(
    rows: list[dict],
    client: PathogenwatchClient,
    output: Path,
    provenance: dict,
    *,
    download_names: set[str] | None = None,
    resume: bool = False,
) -> list[dict]:
    """Discover server-advertised jobs for exact already-resolved public IDs."""
    names = download_names or {"cgmlst"}
    first_export = len(provenance.get("analysis_exports", []))
    eligible = [row for row in rows if any(_analysis_export_eligible(row, name) for name in names)]
    if not eligible:
        return rows
    if not client.api_key:
        provenance["analysis_export_status"] = "unavailable_no_credentials"
        return rows
    output.mkdir(parents=True, exist_ok=True)
    requested = {int(row["numeric_source_id"]): row for row in eligible}
    try:
        capability_path = output / "export_capabilities.json"
        request_contract = {
            "schema": "chronoclade.export-capability-request",
            "schema_version": 1,
            "base_url": client.base_url,
            "requested_numeric_ids": sorted(requested),
            "requested_source_ids": [row["source_genome_id"] for row in eligible],
            "download_names": sorted(names),
            "row_scope": row_scope(eligible),
        }
        if capability_path.exists():
            if not resume:
                raise ExportCheckpointError(
                    "Existing capability evidence requires verified resume or a fresh scope"
                )
            try:
                snapshot = json.loads(capability_path.read_text())
                if (
                    snapshot.get("request_contract") != request_contract
                    or snapshot.get("requested_numeric_ids") != sorted(requested)
                    or snapshot.get("sha256")
                    != content_hash({k: v for k, v in snapshot.items() if k != "sha256"})
                ):
                    raise ValueError
                groups = snapshot["groups"]
            except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
                raise ExportCheckpointError(
                    "Saved capability request/hash differs; use a fresh scope"
                ) from error
        else:
            groups = client.request_json(
                "POST", "/api/genomes/group", body={"ids": sorted(requested)}, authenticated=True
            )
            snapshot = {
                "groups": groups,
                "requested_numeric_ids": sorted(requested),
                "retrieved_at": datetime.now(timezone.utc).isoformat(),
                "request_contract": request_contract,
            }
            snapshot["sha256"] = content_hash(snapshot)
            atomic_json(capability_path, snapshot)
        if not isinstance(groups, list):
            raise ProfileInputError("Unsupported genome-group export capability response")
        returned: set[int] = set()
        resolved: dict[str, dict] = {}
        for index, group in enumerate(groups):
            ids = group.get("ids")
            if not isinstance(ids, list) or any(
                type(ident) is not int or ident not in requested or ident in returned
                for ident in ids
            ):
                raise ProfileInputError(
                    "Genome-group capabilities contain unrequested/repeated IDs"
                )
            returned.update(ids)
            members = [requested[ident] for ident in ids]
            if any(
                row.get("organism_id") and str(row["organism_id"]) != str(group.get("organismId"))
                for row in members
            ):
                raise ProfileInputError(
                    "Genome-group export organism does not match public identity"
                )
            advertised = group.get("downloads") or []
            if download_names is not None:
                allowed = {
                    name
                    for name in download_names
                    if any(_analysis_export_eligible(row, name) for row in members)
                }
                advertised = [item for item in advertised if item.get("name") in allowed]
            if download_names is None and not any(
                item.get("name") == "cgmlst" for item in advertised
            ):
                provenance.setdefault("profile_export_unavailable", []).append(
                    {
                        "organism_id": group.get("organismId"),
                        "reason": "server_does_not_advertise_cgmlst_job",
                    }
                )
            exported = _analysis_exports(
                members, advertised, client, output / f"organism_{index}", provenance, resume=resume
            )
            resolved.update({row["source_genome_id"]: row for row in exported})
        if returned != set(requested):
            provenance.setdefault("profile_export_unavailable", []).append(
                {
                    "reason": "source_ids_omitted_by_export_capability_endpoint",
                    "numeric_ids": sorted(set(requested) - returned),
                }
            )
        result = [resolved.get(row.get("source_genome_id"), row) for row in rows]
        new_exports = provenance.get("analysis_exports", [])[first_export:]
        if any(item.get("download") == "hclink" and item.get("status") == "downloaded"
               for item in new_exports):
            from chronoclade.provider_hiercc import bind_exports
            result = bind_exports(result, {"analysis_exports": new_exports}, capability_path)
        return result
    except (PathogenwatchError, ProfileInputError) as error:
        provenance.setdefault("profile_export_unavailable", []).append({"reason": str(error)})
        return rows


def _accession_input(path: Path) -> tuple[list[str], list[dict]]:
    """Read a one-ID-per-line list or a CSV with an exact accession/source ID."""
    text = path.read_text(encoding="utf-8-sig")
    first = text.splitlines()[0] if text.splitlines() else ""
    overrides: list[dict] = []
    if "," in first or first.strip() in {"accession", "source_genome_id"}:
        overrides = _csv(path)
        identifiers = [
            row.get("accession") or row.get("source_genome_id") or "" for row in overrides
        ]
        if any(not ident for ident in identifiers):
            raise ProfileInputError(
                "Accession CSV requires accession or source_genome_id on every row"
            )
    else:
        identifiers = [
            line.strip()
            for line in text.splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    if not identifiers or len(identifiers) != len(set(identifiers)):
        raise ProfileInputError("Accession list must contain unique nonempty identifiers")
    return identifiers, overrides


def _join_accession(accession: str, candidates: list[dict]) -> dict:
    matches = [
        row
        for row in candidates
        if accession
        in set(row.get("aliases") or [])
        | set(row.get("run_accessions") or [])
        | set(row.get("assembly_accessions") or [])
        | set(row.get("biosample_accessions") or [])
        | {str(row.get("biosample") or ""), str(row.get("source_genome_id") or "")}
    ]
    if len(matches) != 1:
        raise ProfileInputError(
            f"Accession {accession} has {len(matches)} exact matches; supply an exact source_genome_id or frozen catalogue"
        )
    return deepcopy(matches[0])


def _apply_overrides(rows: list[dict], overrides: list[dict], base: Path) -> list[dict]:
    output = deepcopy(rows)
    for override in overrides:
        identifier = (
            override.get("source_genome_id")
            or override.get("accession")
            or override.get("sample_id")
        )
        matches = [
            row
            for row in output
            if identifier in row.get("aliases", []) or identifier == row["sample_id"]
        ]
        if len(matches) != 1:
            raise ProfileInputError(
                f"Metadata override {identifier} does not match exactly one query"
            )
        target = matches[0]
        sample_id = override.get("sample_id")
        if isinstance(sample_id, str) and sample_id.strip():
            target["sample_id"] = sample_id.strip()
        identity_fields = {"source_genome_id", "accession", "sample_id"}
        metadata_fields = {
            key: value
            for key, value in override.items()
            if key not in identity_fields
            and value is not None
            and (not isinstance(value, str) or value.strip())
        }
        if not metadata_fields:
            continue
        conflicts = target.setdefault("metadata_override_conflicts", [])
        for key, value in metadata_fields.items():
            if target.get(key) is not None and target.get(key) != "" and target.get(key) != value:
                conflicts.append(
                    {"field": key, "source_value": target.get(key), "manual_value": value}
                )
            target[key] = value
        if override.get("collection_date"):
            for key in ("date_precision", "date_raw", "date_start", "date_end"):
                target.pop(key, None)
        target["manual_metadata_override"] = True
        target.update(_canonical(target, origin="local", base=base))
    return output


def _sample(row: dict) -> Sample:
    return Sample(
        sample_id=row["sample_id"],
        assembly=Path(row["assembly"]),
        collection_date=row.get("collection_date", ""),
        location=row.get("location") or row.get("country") or "Unknown",
        species=row["species"],
        lineage=row.get("lineage") or "Unassigned",
        origin=row.get("origin", "local"),
        is_reference=str(row.get("is_reference", "")).lower() in {"true", "1", "yes"},
        patient_id=row.get("patient_id", ""),
        country=row.get("country", ""),
        region=row.get("region", ""),
        city=row.get("city", ""),
    )


def materialise_assemblies(
    records: Iterable[dict], *, output: Path, client: PathogenwatchClient | None = None
) -> list[Sample]:
    """Acquire only the records explicitly selected by the caller."""
    rows = deepcopy(list(records))
    missing = [row for row in rows if not row.get("assembly")]
    if missing:
        client = client or PathogenwatchClient()
        if not client.api_key:
            raise ProfileInputError(
                "Assembly acquisition requires a configured Pathogenwatch API key"
            )
        if any(not row.get("source_genome_id") for row in missing):
            raise ProfileInputError("Assembly acquisition requires exact source genome IDs")
        paths, audit = download_assemblies(
            missing,
            output=output / "assemblies",
            cache_dir=output / "assembly_cache",
            api_key=client.api_key,
            base_url=client.base_url,
        )
        (output / "assembly_acquisition.json").write_text(json.dumps(audit, indent=2) + "\n")
        if any(row["source_genome_id"] not in paths for row in missing):
            raise ProfileInputError(
                "Some selected assemblies could not be acquired; see assembly_acquisition.json"
            )
        for row in missing:
            row["assembly"] = str(paths[row["source_genome_id"]])
    return [_sample(row) for row in rows]


def _frozen_cglin(rows: list[dict], assignments: list[dict]) -> list[dict]:
    """Apply the explicit export as the sole cgLIN namespace, including missing IDs.

    Live/native typing can still supply allele profiles. Its cluster release must
    not replace the user's frozen cluster assignments or leak into unmatched IDs.
    """
    clean = [
        {key: value for key, value in row.items() if not key.startswith("cglin_") and key != "cgst"}
        for row in rows
    ]
    annotated = annotate_catalogue(clean, assignments)
    for row in annotated:
        row["cglin_scope_source"] = "user_frozen_export"
        row["cglin_database_version_status"] = (
            "unavailable" if row.get("cglin_scheme_version") == "unknown" else "user_declared"
        )
    return annotated
