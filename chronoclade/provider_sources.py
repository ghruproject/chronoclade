"""Explicit raw Pathogenwatch source replay, never normalized-cache adoption."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import json
import re

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import DatasetError
from chronoclade.file_snapshots import copy_evidence_snapshot
from chronoclade.pathogenwatch import (
    PathogenwatchError,
    content_hash,
    load_raw_catalogue,
    replay_catalogue,
)
from chronoclade.profile_inputs import parse_cgmlst_export
from chronoclade.export_checkpoints import checkpoint_path, load_checkpoint, ExportCheckpointError

SCHEMA = "chronoclade.provider-source"
LEGACY_RECEIPT = "derived_from_saved_capability_order_and_legacy_batch100"


def _read(path):
    try:
        value = json.loads(Path(path).read_text())
        if not isinstance(value, dict):
            raise ValueError("Provider source JSON must be an object")
        return value
    except (OSError, ValueError) as error:
        raise DatasetError("Cannot read frozen provider source JSON") from error


def _inside(root, relative):
    if not isinstance(relative, str) or Path(relative).is_absolute():
        raise DatasetError("Provider source paths must be relative")
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file():
        raise DatasetError("Provider source path escapes its bundle or is missing")
    return path


def _ref(root, path):
    return {
        "path": path.relative_to(root).as_posix(),
        "sha256": file_sha256(path),
        "bytes": path.stat().st_size,
    }


def _verify_ref(root, value):
    if not isinstance(value, dict):
        raise DatasetError("Provider source requires checksummed file references")
    path = _inside(root, value.get("path"))
    if value.get("sha256") != file_sha256(path) or value.get("bytes") != path.stat().st_size:
        raise DatasetError("Provider source file hash/size mismatch")
    return path


def _catalogue(path):
    try:
        capture = load_raw_catalogue(path)
        rows = replay_catalogue(capture)
        if {r["source_genome_id"] for r in rows} != set(capture["raw_details"]):
            raise ValueError("Replay changed exact public source identities")
        provenance = capture["provenance"]
        query = provenance["query"]
        scheme = "mlst" if "mlst" in query else "mlst2"
        if not isinstance(query.get(scheme), list) or len(query[scheme]) != 1:
            raise ValueError("Exactly one MLST ST is required per source scope")
        organism = provenance["organism_capabilities"]
        if (
            str(organism.get("organismId")) != str(query["organismId"])
            or organism.get("fullName", "").casefold() != provenance["species"].casefold()
            or any(
                r.get("species", "").casefold() != provenance["species"].casefold() for r in rows
            )
            or not rows
            or not str(query[scheme][0]).isdigit()
            or ("mlst" in query and "mlst2" in query)
        ):
            raise ValueError("Raw catalogue organism/species/MLST scope is ambiguous")
        base_url = provenance.get("base_url")
        if not isinstance(base_url, str) or not re.fullmatch(
            r"https://[A-Za-z0-9.-]+", base_url.rstrip("/")
        ):
            raise ValueError("Provider source requires an explicit HTTPS base URL")
        if any(
            type(r.get("numeric_source_id")) not in {int, str}
            or not re.fullmatch(r"[0-9]+", str(r["numeric_source_id"]))
            for r in rows
        ):
            raise ValueError("Public source numeric IDs must be exact positive decimal identifiers")
        crosswalk = {int(r["numeric_source_id"]): r["source_genome_id"] for r in rows}
        if len(crosswalk) != len(rows) or any(i < 1 for i in crosswalk):
            raise ValueError("Public source UUID/numeric crosswalk is not one-to-one")
        return (
            capture,
            rows,
            crosswalk,
            {
                "base_url": base_url.rstrip("/"),
                "species": provenance["species"],
                "organism_id": str(query["organismId"]),
                "mlst_scheme": scheme,
                "st": str(query[scheme][0]),
                "query_sha256": content_hash(query),
                "retrieved_at": provenance["retrieved_at"],
                "selection_scope": "whole_st_catalogue",
                "source_ids_sha256": content_hash(sorted(crosswalk.values())),
                "numeric_ids_sha256": content_hash(sorted(crosswalk)),
                "expected_count": len(rows),
            },
        )
    except (KeyError, TypeError, ValueError, OSError, PathogenwatchError) as error:
        raise DatasetError(f"Invalid raw provider catalogue: {error}") from error


def _groups(path, crosswalk, organism):
    payload = _read(path)
    if payload.get("sha256") != content_hash({k: v for k, v in payload.items() if k != "sha256"}):
        raise DatasetError("Provider capability snapshot hash mismatch")
    requested = payload.get("requested_numeric_ids")
    if (
        not isinstance(requested, list)
        or any(type(i) is not int for i in requested)
        or len(set(requested)) != len(requested)
        or set(requested) != set(crosswalk)
    ):
        raise DatasetError("Whole-ST source requires exact full catalogue profile request coverage")
    seen, batches = set(), []
    groups = payload.get("groups")
    if not isinstance(groups, list) or not groups:
        raise DatasetError("Provider source has no capability groups")
    for index, group in enumerate(groups):
        if (
            not isinstance(group, dict)
            or not isinstance(group.get("downloads"), list)
            or any(not isinstance(d, dict) for d in group["downloads"])
        ):
            raise DatasetError("Malformed provider capability group")
        ids = group.get("ids")
        jobs = [d.get("job") for d in group.get("downloads", []) if d.get("name") == "cgmlst"]
        if (
            str(group.get("organismId")) != organism
            or not isinstance(ids, list)
            or not ids
            or any(type(i) is not int or i not in crosswalk or i in seen for i in ids)
            or len(set(ids)) != len(ids)
        ):
            raise DatasetError(
                "Provider capability group has wrong organism or repeated/unrequested IDs"
            )
        if (
            len(jobs) != 1
            or not isinstance(jobs[0], str)
            or not re.fullmatch(r"[A-Za-z0-9_-]+", jobs[0])
        ):
            raise DatasetError("Whole-ST source needs exactly one advertised cgMLST job per group")
        seen.update(ids)
        for start in range(0, len(ids), 100):
            batches.append((index, jobs[0], start, ids[start : start + 100]))
    if seen != set(crosswalk):
        raise DatasetError("Capability groups omit whole-ST catalogue IDs")
    return batches


def _receipt(root, response, job, numeric, crosswalk, base_url):
    uuids = [crosswalk[i] for i in numeric]
    try:
        parsed = parse_cgmlst_export(response.read_bytes(), uuids, job=job)
    except (ValueError, OSError) as error:
        raise DatasetError(f"Invalid frozen cgMLST response: {error}") from error
    returned = list(dict.fromkeys(r["source_genome_id"] for r in parsed))
    return {
        "download": "cgmlst",
        "job": job,
        "base_url": base_url,
        "endpoint": "/api/downloads/cgmlst",
        "route": "/api/downloads/cgmlst?job=" + job,
        "source_genome_ids": uuids,
        "numeric_request_ids": [str(i) for i in numeric],
        "retry_attempt": 1,
        "saved_response_path": response.relative_to(root).as_posix(),
        "response_sha256": file_sha256(response),
        "response_bytes": response.stat().st_size,
        "returned_source_genome_ids": returned,
        "missing_source_genome_ids": [ident for ident in uuids if ident not in set(returned)],
        "status": "accepted",
        "request_evidence": LEGACY_RECEIPT,
    }, parsed


class ProviderSource:
    """Validated raw responses with current parser output; no provider client."""

    def __init__(self, source):
        source = Path(source).expanduser().resolve()
        manifest = source / "provider-source.json" if source.is_dir() else source
        self.root = source if source.is_dir() else source.parent
        self.scopes = []
        if manifest.is_file():
            payload = _read(manifest)
            if (
                payload.get("schema") != SCHEMA
                or payload.get("schema_version") != 1
                or payload.get("snapshot_sha256")
                != content_hash({k: v for k, v in payload.items() if k != "snapshot_sha256"})
            ):
                raise DatasetError("Invalid provider source manifest/hash")
            scopes = payload.get("scopes")
            if not isinstance(scopes, list) or not scopes:
                raise DatasetError("Provider source has no complete scopes")
            for scope in scopes:
                if not isinstance(scope, dict):
                    raise DatasetError("Malformed provider source scope")
                self._load_scope(
                    _verify_ref(self.root, scope.get("catalogue")),
                    _verify_ref(self.root, scope.get("capabilities")),
                    scope,
                )
        else:
            if not source.is_dir():
                raise DatasetError("Provider source manifest does not exist")
            public = sorted(source.glob("public_*.json"))
            if not public:
                raise DatasetError("No frozen public catalogue exists in provider source")
            for catalogue in public:
                match = re.fullmatch(r"public_(\d+)\.json", catalogue.name)
                if not match:
                    raise DatasetError("Ambiguous legacy catalogue filename")
                capabilities = list(
                    source.glob(f"scope_{match[1]}_*/profiles/export_capabilities.json")
                )
                stable = source / f"scope_{match[1]}" / "profiles" / "export_capabilities.json"
                if stable.is_file():
                    capabilities.append(stable)
                if len(capabilities) != 1:
                    raise DatasetError(
                        "Legacy source requires exactly one capability capture per scope; no newest-file guessing"
                    )
                self._load_scope(catalogue, capabilities[0], None)
        keys = [
            (s["scope"]["scope"]["species"].casefold(), s["scope"]["scope"]["st"])
            for s in self.scopes
        ]
        if (
            len(set(keys)) != len(keys)
            or len({s["scope"]["scope"]["base_url"] for s in self.scopes}) != 1
        ):
            raise DatasetError("Provider source scopes repeat or mix base URLs")
        self.base_url = self.scopes[0]["scope"]["scope"]["base_url"]
        self.identity = {
            "schema": SCHEMA,
            "schema_version": 1,
            "scopes": [s["scope"] for s in self.scopes],
        }
        self.identity["snapshot_sha256"] = content_hash(self.identity)

    def _load_scope(self, catalogue, capabilities, saved):
        catalogue = _inside(self.root, catalogue.relative_to(self.root).as_posix())
        capabilities = _inside(self.root, capabilities.relative_to(self.root).as_posix())
        capture, rows, crosswalk, scope = _catalogue(catalogue)
        batches = _groups(capabilities, crosswalk, scope["organism_id"])
        receipts, assignments, expected_files = [], {}, set()
        if saved is not None:
            if (
                saved.get("scope") != scope
                or not isinstance(saved.get("receipts"), list)
                or any(not isinstance(r, dict) for r in saved.get("receipts", []))
            ):
                raise DatasetError("Provider source declared scope differs from raw catalogue")
            saved_receipts = saved["receipts"]
            if len(saved_receipts) != len(batches):
                raise DatasetError("Provider source batch list is incomplete")
        for position, (group, job, start, numeric) in enumerate(batches):
            checkpoint_receipt = None
            if saved is None:
                checkpoint = checkpoint_path(
                    capabilities.parent / f"organism_{group}", "cgmlst", job, start
                )
                if checkpoint.is_file():
                    try:
                        checkpoint_receipt, _, contract = load_checkpoint(
                            checkpoint, root=capabilities.parent, parser=parse_cgmlst_export
                        )
                    except ExportCheckpointError as error:
                        raise DatasetError(str(error)) from error
                    if (
                        contract.get("download") != "cgmlst"
                        or contract.get("job") != job
                        or contract.get("base_url") != scope["base_url"]
                        or contract.get("numeric_request_ids") != [str(i) for i in numeric]
                        or contract.get("source_genome_ids") != [crosswalk[i] for i in numeric]
                    ):
                        raise DatasetError(
                            "Recorded checkpoint request differs from whole-ST scope"
                        )
                    response = capabilities.parent / checkpoint_receipt["saved_response_path"]
                else:
                    response = (
                        capabilities.parent / f"organism_{group}" / f"cgmlst_{job}_{start}.bin"
                    )
                if not response.is_file() or not response.resolve().is_relative_to(self.root):
                    raise DatasetError("Provider source is missing an exact legacy batch response")
            else:
                response = _inside(self.root, saved_receipts[position].get("saved_response_path"))
            receipt, parsed = _receipt(
                self.root, response, job, numeric, crosswalk, scope["base_url"]
            )
            if checkpoint_receipt is not None:
                expected = dict(
                    checkpoint_receipt,
                    saved_response_path=response.relative_to(self.root).as_posix(),
                )
                receipt.update(
                    request_evidence=expected["request_evidence"],
                    retry_attempt=expected["retry_attempt"],
                )
                if receipt != expected:
                    raise DatasetError(
                        "Recorded raw batch receipt differs from verified source response"
                    )
            if saved is not None:
                prior = saved_receipts[position]
                marker = prior.get("request_evidence")
                attempt = prior.get("retry_attempt")
                if (
                    marker not in {LEGACY_RECEIPT, "recorded_exact_http_request"}
                    or type(attempt) is not int
                    or not 1 <= attempt <= 3
                ):
                    raise DatasetError("Invalid provider request evidence or retry attempt")
                receipt.update(request_evidence=marker, retry_attempt=attempt)
                if prior != receipt:
                    raise DatasetError(
                        "Provider batch receipt identity/hash differs from validated response"
                    )
            if response in expected_files:
                raise DatasetError("Provider response is repeated across batches")
            expected_files.add(response)
            receipts.append(receipt)
            for row in parsed:
                if row["source_genome_id"] in assignments:
                    raise DatasetError("Duplicate exported source across batches")
                assignments[row["source_genome_id"]] = row
        if saved is None and set(capabilities.parent.rglob("cgmlst_*.bin")) != {
            p for p in expected_files if p.name.startswith("cgmlst_") and p.suffix == ".bin"
        }:
            raise DatasetError(
                "Ambiguous/unrequested legacy response or retry; no filename guessing"
            )
        if saved is None:
            recorded = {
                p.parent / "response.raw"
                for p in capabilities.parent.glob("organism_*/batch_cgmlst_*/receipt.json")
            }
            if recorded != {p for p in expected_files if p.name == "response.raw"}:
                raise DatasetError(
                    "Unrequested or incomplete recorded batch checkpoint in provider source"
                )
        scope_record = {
            "scope": scope,
            "catalogue": _ref(self.root, catalogue),
            "capabilities": _ref(self.root, capabilities),
            "receipts": receipts,
            "complete": True,
        }
        if saved is not None and saved != scope_record:
            raise DatasetError("Provider source scope manifest differs from verified raw evidence")
        self.scopes.append(
            {"scope": scope_record, "capture": capture, "rows": rows, "assignments": assignments}
        )

    def for_scope(self, species, st):
        matches = [
            s
            for s in self.scopes
            if s["scope"]["scope"]["species"].casefold() == species.casefold()
            and s["scope"]["scope"]["st"] == str(st)
        ]
        if len(matches) != 1:
            raise DatasetError(
                "Frozen provider source does not contain the exact requested species/ST"
            )
        value = matches[0]
        return value["capture"], deepcopy(value["rows"]), deepcopy(value["assignments"])

    def snapshot(self, directory):
        """Publish independent raw evidence, with portable receipts and content identity."""
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        scopes = []
        for scope_index, value in enumerate(self.scopes):
            scope = deepcopy(value["scope"])
            paths = [scope["catalogue"]["path"], scope["capabilities"]["path"]] + [
                r["saved_response_path"] for r in scope["receipts"]
            ]
            copied = {}
            for file_index, relative in enumerate(paths):
                source = _inside(self.root, relative)
                digest = file_sha256(source)
                dest = directory / "assets" / f"{digest}-{scope_index}-{file_index}{source.suffix}"
                dest.parent.mkdir(exist_ok=True)
                if not dest.exists():
                    copy_evidence_snapshot(source, dest)
                if file_sha256(dest) != digest or file_sha256(source) != digest:
                    raise DatasetError("Provider source changed while snapshotting")
                copied[relative] = dest.relative_to(directory).as_posix()
            scope["catalogue"]["path"] = copied[scope["catalogue"]["path"]]
            scope["capabilities"]["path"] = copied[scope["capabilities"]["path"]]
            for receipt in scope["receipts"]:
                receipt["saved_response_path"] = copied[receipt["saved_response_path"]]
            scopes.append(scope)
        payload = {"schema": SCHEMA, "schema_version": 1, "scopes": scopes}
        payload["snapshot_sha256"] = content_hash(payload)
        path = directory / "provider-source.json"
        if path.exists() and _read(path) != payload:
            raise DatasetError("Existing provider source snapshot differs")
        write_json(path, payload)
        return path
