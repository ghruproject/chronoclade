"""Validated descriptive scopes for already provider-reported HierCC labels.

These keys identify a provider namespace, never a native assignment database.
"""

from __future__ import annotations

from typing import NamedTuple
import hashlib
import json
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit

from chronoclade.pathogenwatch import content_hash
from chronoclade.pathogenwatch_hiercc import parse_hiercc_export

POLICY = "provider_reported_hiercc_v1"
DISCLOSURE = "Provider-reported HierCC; reference release not exposed."
_FIELD = "_provider_reported_hiercc"


class ReportedScope(NamedTuple):
    namespace: tuple[str, ...]
    snapshot_sha256: str
    response_sha256: str
    observed_at: str
    assignment_sha256: str


def _assignment_hash(row):
    fields = (
        "source_genome_id",
        "hiercc_source",
        "hiercc_provider_job",
        "hiercc_scheme",
        "hiercc_scheme_version",
        "hiercc_database_version",
        "hiercc_database_sha256",
        "hiercc_export_sha256",
        "hiercc_codes",
    )
    return content_hash({key: row.get(key) for key in fields})


def reported_metadata(row):
    scope = reported_scope(row)
    return (
        {"grouping_policy": POLICY, "provider_reported_namespace": scope.namespace} if scope else {}
    )


def normalized_server(value):
    parsed = urlsplit(value)
    if (
        parsed.scheme != "https"
        or not parsed.hostname
        or parsed.username
        or parsed.password
        or parsed.query
        or parsed.fragment
    ):
        raise ValueError("Reported HierCC requires a credential-free HTTPS provider URL")
    host = parsed.hostname.lower()
    if parsed.port not in (None, 443):
        host += f":{parsed.port}"
    return urlunsplit(("https", host, parsed.path.rstrip("/"), "", ""))


def reported_scope(row):
    """Only a validated in-memory proof object enables the descriptive policy."""
    value = row.get(_FIELD)
    return (
        value
        if isinstance(value, ReportedScope) and value.assignment_sha256 == _assignment_hash(row)
        else None
    )


def grouping_scope(row):
    scope = reported_scope(row)
    if scope is None:
        return None
    return row.get("hiercc_scheme"), scope.namespace


def _clear(rows):
    return without_internal_proof(rows)


def _read(path):
    return json.loads(Path(path).read_text())


def _validate_capability(path):
    value = _read(path)
    if value.get("sha256") != content_hash({k: v for k, v in value.items() if k != "sha256"}):
        raise ValueError("Reported HierCC capability snapshot hash differs")
    request = value.get("request_contract") or {}
    if (
        request.get("schema") != "chronoclade.export-capability-request"
        or request.get("schema_version") != 1
    ):
        raise ValueError("Reported HierCC requires a recorded exact capability request")
    server = normalized_server(request["base_url"])
    numeric = request["requested_numeric_ids"]
    ids = request["requested_source_ids"]
    scopes = request["row_scope"]
    if (
        not ids
        or len(ids) != len(set(ids))
        or len(ids) != len(numeric)
        or len(numeric) != len(set(numeric))
        or len(scopes) != len(ids)
        or numeric != value["requested_numeric_ids"]
    ):
        raise ValueError("Reported HierCC capability identity is incomplete")
    # The recorded request stores numeric IDs sorted independently of source order.
    # Group responses and batch receipts provide the exact paired crosswalk below.
    advertised = {}
    for group in value["groups"]:
        for ident in group["ids"]:
            if ident not in numeric or ident in advertised:
                raise ValueError("Reported HierCC capability IDs differ from request")
            advertised[ident] = group
    if set(advertised) != set(numeric):
        raise ValueError("Reported HierCC capability inventory is incomplete")
    return value, server, advertised


def bind_exports(rows, provenance, capability_path, *, artifact_hashes=None):
    """Reparse accepted raw exports; never infer source authority from row strings."""
    rows = [dict(row) for row in rows]
    path = Path(capability_path)
    capability, server, advertised = _validate_capability(path)
    if urlsplit(server).hostname not in {"pathogen.watch", "next.pathogen.watch"}:
        return [
            dict(
                row,
                hiercc_reported_grouping_limitation="Provider is outside the supported Pathogenwatch namespace",
            )
            for row in rows
        ]
    for export in provenance.get("analysis_exports", []):
        if export.get("download") != "hclink" or export.get("status") != "downloaded":
            continue
        receipts = export.get("batch_receipts", [])
        if export.get("requested_count") != sum(
            len(r.get("source_genome_ids", [])) for r in receipts
        ):
            raise ValueError("Reported HierCC complete export request inventory differs")
        for receipt in receipts:
            if normalized_server(receipt["base_url"]) != server:
                continue
            job = receipt["job"]
            ids, numeric = receipt["source_genome_ids"], receipt["numeric_request_ids"]
            if (
                receipt.get("status") != "accepted"
                or receipt.get("request_evidence") != "recorded_exact_http_request"
                or receipt.get("download") != "hclink"
                or job != export.get("job")
                or receipt.get("endpoint") != "/api/downloads/hclink"
                or receipt.get("route") != "/api/downloads/hclink?job=" + job
                or not ids
                or len(ids) != len(numeric)
                or len(ids) != len(set(ids))
                or len(numeric) != len(set(numeric))
            ):
                raise ValueError("Reported HierCC exact export receipt is invalid")
            if not set(ids).issubset(capability["request_contract"]["requested_source_ids"]):
                continue  # A different explicitly recorded export scope.
            if any(not str(n).isdigit() or str(int(n)) != str(n) or int(n) < 1 for n in numeric):
                raise ValueError("Reported HierCC numeric request IDs are not canonical")
            if (
                type(receipt.get("retry_attempt")) is not int
                or not 1 <= receipt["retry_attempt"] <= 3
            ):
                raise ValueError("Reported HierCC retry identity is invalid")
            for ident in ids:
                scope = capability["request_contract"]["row_scope"][
                    capability["request_contract"]["requested_source_ids"].index(ident)
                ]
                if (
                    str(scope.get("organism_id")) != "562"
                    or scope.get("species", "").casefold() != "escherichia coli"
                ):
                    raise ValueError("Reported HierCC capability organism scope differs")
            for ident in numeric:
                group = advertised.get(int(ident), {})
                if str(group.get("organismId")) != "562" or not any(
                    item.get("name") == "hclink" and item.get("job") == job
                    for item in group.get("downloads", [])
                ):
                    raise ValueError(
                        "Reported HierCC job/organism was not advertised for exact IDs"
                    )
            relative = Path(receipt["saved_response_path"])
            raw = (path.parent / relative).resolve()
            if relative.is_absolute() or not raw.is_relative_to(path.parent.resolve()):
                raise ValueError("Reported HierCC raw response escapes source scope")
            if (
                artifact_hashes is not None
                and artifact_hashes.get(raw) != receipt["response_sha256"]
            ):
                continue  # A receipt from a different inventoried export scope.
            data = raw.read_bytes()
            digest = hashlib.sha256(data).hexdigest()
            if digest != receipt["response_sha256"] or len(data) != receipt["response_bytes"]:
                raise ValueError("Reported HierCC raw response hash/size differs")
            if artifact_hashes is not None and artifact_hashes.get(raw) != digest:
                raise ValueError("Reported HierCC raw response is not bound by source inventory")
            parsed = parse_hiercc_export(data, ids, job=job)
            if (
                receipt.get("returned_source_genome_ids") != ids
                or receipt.get("missing_source_genome_ids") != []
            ):
                raise ValueError("Reported HierCC returned ID inventory differs")
            by_id = {
                item["source_genome_id"]: (item, str(number))
                for item, number in zip(parsed, numeric)
            }
            snapshot = content_hash(
                {"policy": POLICY, "capability": capability["sha256"], "receipt": receipt}
            )
            for row in rows:
                ident = row.get("source_genome_id")
                if ident not in by_id:
                    continue
                assignment, number = by_id[ident]
                if (
                    row.get("species", "").casefold() != "escherichia coli"
                    or str(row.get("numeric_source_id")) != number
                ):
                    raise ValueError("Reported HierCC sample identity/organism differs")
                if any(
                    row.get(key) != value
                    for key, value in assignment.items()
                    if key.startswith("hiercc_")
                    and key
                    in {
                        "hiercc_source",
                        "hiercc_provider_job",
                        "hiercc_scheme",
                        "hiercc_scheme_version",
                        "hiercc_database_version",
                        "hiercc_database_sha256",
                        "hiercc_export_sha256",
                        "hiercc_codes",
                    }
                ):
                    raise ValueError("Reported HierCC assignment differs from exact raw source")
                row[_FIELD] = ReportedScope(
                    (POLICY, server, "562", job, assignment["hiercc_scheme"]),
                    snapshot,
                    digest,
                    capability["retrieved_at"],
                    _assignment_hash(row),
                )
    return rows


def bind_dataset(rows, dataset, source_root):
    """Bind portable saved proof using its declared inventory, including source datasets."""
    rows = _clear(rows)
    if source_root is None or not any(
        r.get("hiercc_source") == "pathogenwatch_hclink_analysis_export" for r in rows
    ):
        return rows
    root = Path(source_root).resolve()
    if root.is_file():
        root = root.parent
    audit = root / "prepare.json"
    if not audit.exists():
        return rows
    inventory = {}
    for item in _read(audit)["artifacts"]:
        relative = Path(item["path"])
        path = (root / relative).resolve()
        if relative.is_absolute() or not path.is_relative_to(root):
            raise ValueError("Reported HierCC source inventory escapes dataset")
        inventory[path] = item["sha256"]
    # Enumerate only explicitly inventoried capability snapshots, never newest-file guesses.
    caps = [p for p in inventory if p.name == "export_capabilities.json"]
    provenances = [
        dataset.parameters.get("provider_prepare", {}),
        *(
            p.get("exports", {})
            for p in dataset.parameters.get("profile_context", {}).get("pools", [])
        ),
    ]
    for path in inventory:
        if path.name == "dataset.json" and "source_dataset" in path.parts:
            saved = _read(path).get("parameters", {})
            provenances.append(saved.get("provider_prepare", {}))
    for path in caps:
        if hashlib.sha256(path.read_bytes()).hexdigest() != inventory[path]:
            raise ValueError("Reported HierCC saved capability bytes changed")
        if "request_contract" not in _read(path):
            for row in rows:
                if row.get("hiercc_source") == "pathogenwatch_hclink_analysis_export":
                    row["hiercc_reported_grouping_limitation"] = (
                        "Exact capability request provenance was not recorded in this older snapshot"
                    )
            continue
        for provenance in provenances:
            rows = bind_exports(rows, provenance, path, artifact_hashes=inventory)
    return rows


def rebind_scope_exports(rows, provenance, capability_path):
    """A JSON checkpoint is not proof: revalidate its exact retained raw exports."""
    rows = without_internal_proof(rows)
    if not any(
        item.get("download") == "hclink" and item.get("status") == "downloaded"
        for item in provenance.get("analysis_exports", [])
    ):
        return rows
    rows = bind_exports(rows, provenance, capability_path)
    if any(
        row.get("hiercc_source") == "pathogenwatch_hclink_analysis_export"
        and row.get("hiercc_codes")
        and reported_scope(row) is None
        for row in rows
    ):
        raise ValueError("Cached provider-reported HierCC source proof could not be reconstructed")
    return rows


def source_snapshot(rows):
    scopes = [reported_scope(row) for row in rows if reported_scope(row)]
    if not scopes:
        return None
    return {
        "policy": POLICY,
        "source_snapshot_sha256": content_hash(sorted({scope.snapshot_sha256 for scope in scopes})),
        "source_receipt_sha256": sorted({scope.snapshot_sha256 for scope in scopes}),
        "reported_samples": len(scopes),
        "namespaces": [
            list(namespace) for namespace in sorted({scope.namespace for scope in scopes})
        ],
        "observed_at": sorted({scope.observed_at for scope in scopes}),
        "disclosure": DISCLOSURE,
    }


def without_internal_proof(rows):
    return [{key: value for key, value in row.items() if key != _FIELD} for row in rows]
