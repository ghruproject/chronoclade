"""Immutable request-bound raw exports; never cached normalized typing results."""

from pathlib import Path
import hashlib
import json
import tempfile
from urllib.parse import urlsplit, urlencode

from chronoclade.pathogenwatch import content_hash

SCHEMA = "chronoclade.analysis-export-batch"
VERSION = 1


class ExportCheckpointError(RuntimeError):
    """Saved raw evidence is incompatible or unverifiable; do not silently reuse."""


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def row_scope(rows):
    keys = (
        "species",
        "organism_id",
        "mlst_st",
        "cgmlst_scheme",
        "cgmlst_scheme_version",
        "cglin_scheme",
        "cglin_scheme_version",
    )
    return [{key: row.get(key) for key in keys} for row in rows]


def batch_contract(rows, *, name, job, base_url):
    parsed = urlsplit(base_url)
    if (
        parsed.scheme != "https"
        or not parsed.hostname
        or parsed.username
        or parsed.password
        or parsed.query
        or parsed.fragment
    ):
        raise ExportCheckpointError(
            "Raw export checkpoint requires a credential-free HTTPS base URL"
        )
    source_ids = [row["source_genome_id"] for row in rows]
    numeric_ids = [str(row["numeric_source_id"]) for row in rows]
    if (
        not rows
        or len(set(source_ids)) != len(rows)
        or len(set(numeric_ids)) != len(rows)
        or any(not isinstance(ident, str) or not ident for ident in source_ids)
        or any(
            not value.isdigit() or int(value) < 1 or str(int(value)) != value
            for value in numeric_ids
        )
    ):
        raise ExportCheckpointError(
            "Export batch requires unique exact source IDs and positive canonical numeric IDs"
        )
    return {
        "download": name,
        "job": job,
        "base_url": base_url,
        "endpoint": f"/api/downloads/{name}",
        "route": f"/api/downloads/{name}?" + urlencode({"job": job}),
        "source_genome_ids": [row["source_genome_id"] for row in rows],
        "numeric_request_ids": [str(row["numeric_source_id"]) for row in rows],
        "row_scope": row_scope(rows),
    }


def checkpoint_path(output, name, job, start):
    return Path(output) / f"batch_{name}_{job}_{start}" / "receipt.json"


def persist_checkpoint(path, contract, receipt, data):
    """Publish raw bytes and receipt together with one atomic directory rename."""
    path = Path(path)
    directory = path.parent
    if directory.exists():
        raise ExportCheckpointError("Existing committed batch directory cannot be overwritten")
    payload = {
        "schema": SCHEMA,
        "schema_version": VERSION,
        "contract": contract,
        "receipt": receipt,
    }
    payload["sha256"] = content_hash(payload)
    # Unpublished private stages are preserved after interruption but never reused.
    # response.raw is not a legacy export filename and has no accepted receipt.
    stage = Path(tempfile.mkdtemp(prefix=".pending-" + directory.name + "-", dir=directory.parent))
    (stage / "response.raw").write_bytes(data)
    atomic_json(stage / "receipt.json", payload)
    stage.rename(directory)


def load_checkpoint(path, *, root, parser, expected=None):
    """Reparse exact saved bytes and verify request/response identity every time."""
    try:
        payload = json.loads(Path(path).read_text())
        if (
            not isinstance(payload, dict)
            or payload.get("schema") != SCHEMA
            or type(payload.get("schema_version")) is not int
            or payload.get("schema_version") != VERSION
            or payload.get("sha256")
            != content_hash({k: v for k, v in payload.items() if k != "sha256"})
        ):
            raise ValueError("unsupported schema or invalid checkpoint hash")
        contract, receipt = payload["contract"], payload["receipt"]
        if not isinstance(contract, dict) or not isinstance(receipt, dict):
            raise ValueError("invalid request/receipt object")
        if expected is not None and contract != expected:
            raise ValueError("immutable export request changed")
        for field in (
            "download",
            "job",
            "base_url",
            "endpoint",
            "route",
            "source_genome_ids",
            "numeric_request_ids",
        ):
            if receipt.get(field) != contract.get(field):
                raise ValueError("receipt differs from exact request")
        if (
            receipt.get("status") != "accepted"
            or receipt.get("request_evidence") != "recorded_exact_http_request"
        ):
            raise ValueError("batch was not accepted with recorded request evidence")
        if type(receipt.get("retry_attempt")) is not int or not 1 <= receipt["retry_attempt"] <= 3:
            raise ValueError("invalid retry identity")
        relative = receipt["saved_response_path"]
        if not isinstance(relative, str) or Path(relative).is_absolute():
            raise ValueError("response path must be relative")
        root = Path(root).resolve()
        response = (root / relative).resolve()
        if not response.is_relative_to(root) or not response.is_file():
            raise ValueError("response path escapes evidence root or is missing")
        data = response.read_bytes()
        if len(data) != receipt.get("response_bytes") or hashlib.sha256(
            data
        ).hexdigest() != receipt.get("response_sha256"):
            raise ValueError("saved response hash/size differs")
        parsed = parser(data, contract["source_genome_ids"], job=contract["job"])
        returned = list(dict.fromkeys(str(row["source_genome_id"]) for row in parsed))
        if receipt.get("returned_source_genome_ids") != returned or receipt.get(
            "missing_source_genome_ids"
        ) != [ident for ident in contract["source_genome_ids"] if ident not in set(returned)]:
            raise ValueError("saved parsed identity accounting changed")
        return receipt, parsed, contract
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
        raise ExportCheckpointError(
            f"Unverified raw export checkpoint {Path(path).name}; use a fresh output scope"
        ) from error
