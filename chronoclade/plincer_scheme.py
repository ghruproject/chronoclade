"""Validate LIN percentage bins and preserve one pinned upstream correction."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import tomllib

ORIGINAL_SHA256 = "c52db6d411e8c5e677f4124c290464c3c2e46b9857861e15fc71daac8c35414c"
CORRECTION_ID = "plincer-kpsc-published-629-boundaries-v1"
SOURCE_URL = "https://raw.githubusercontent.com/pathogenwatch-oss/klebsiella-lincodes/c6d3d85e7484aee4d713429b78402b0590ac8fdc/scheme.toml"
THRESHOLD_SOURCE = "https://bigsdb.pasteur.fr/klebsiella/cgmlst-lincodes/"
OLD = b'{name = "6", min = 98.4102, max = 798.8871}'
NEW = b'{name = "6", min = 98.4102, max = 98.8871}'
NAMES = ["species", "subspecies", "sublineage", "clonalGroup", "5", "6", "7", "8", "9", "10"]


class SchemeValidationError(ValueError):
    """An assignment scheme is malformed or its correction is unverified."""


def _sha(data):
    return hashlib.sha256(data).hexdigest()


PUBLISHED_BOUNDS = [
    0,
    3.0207,
    6.9952,
    69.7933,
    93.1638,
    98.4102,
    98.8871,
    99.3641,
    99.6820,
    99.8410,
    100,
]


def validate_bytes(data: bytes, *, require_published=False) -> dict:
    try:
        scheme = tomllib.loads(data.decode("utf-8"))
    except (UnicodeError, ValueError) as error:
        raise SchemeValidationError("plincer scheme is not readable TOML") from error
    if (scheme.get("scheme_id"), scheme.get("host"), scheme.get("host_path")) != (
        18,
        "pasteur",
        "pubmlst_klebsiella_seqdef",
    ):
        raise SchemeValidationError("plincer scheme must identify Pasteur Klebsiella scheme 18")
    levels = scheme.get("levels")
    if not isinstance(levels, list) or len(levels) != 10:
        raise SchemeValidationError("plincer scheme must contain ten ordered LIN bins")
    previous = 0
    for name, level in zip(NAMES, levels):
        if not isinstance(level, dict) or level.get("name") != name:
            raise SchemeValidationError("plincer scheme has unexpected LIN bin names/order")
        lo, hi = level.get("min"), level.get("max")
        if (
            any(
                isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x)
                for x in (lo, hi)
            )
            or not 0 <= lo < hi <= 100
        ):
            raise SchemeValidationError(
                "plincer LIN bin bounds must be finite percentages in [0, 100]"
            )
        if lo != previous:
            raise SchemeValidationError("plincer LIN bins must be contiguous without gaps/overlap")
        previous = hi
    if previous != 100:
        raise SchemeValidationError("plincer LIN bins must cover the complete [0, 100] range")
    if (
        require_published
        and [levels[0]["min"], *[level["max"] for level in levels]] != PUBLISHED_BOUNDS
    ):
        raise SchemeValidationError(
            "plincer scheme differs from published KpSC 629-locus boundaries; custom schemes are unverified"
        )
    return scheme


def _correction(original, policy="level6-only"):
    if policy not in {"level6-only", "published-kpsc-629"}:
        raise SchemeValidationError("Unknown plincer correction policy")
    if _sha(original) != ORIGINAL_SHA256 or original.count(OLD) != 1:
        raise SchemeValidationError(
            "Unknown malformed plincer scheme; no automatic correction allowed"
        )
    corrected = original.replace(OLD, NEW, 1)
    transformations = [{"level": "6", "field": "max", "from": 798.8871, "to": 98.8871}]
    if policy == "published-kpsc-629":
        corrected = corrected.replace(b"99.3461", b"99.3641")
        transformations.extend(
            [
                {"level": "7", "field": "max", "from": 99.3461, "to": 99.3641},
                {"level": "8", "field": "min", "from": 99.3461, "to": 99.3641},
            ]
        )
    validate_bytes(corrected, require_published=policy == "published-kpsc-629")
    receipt = {
        "schema": "chronoclade.plincer-scheme-correction",
        "schema_version": 1,
        "correction_id": CORRECTION_ID,
        "policy": policy,
        "original_file": "scheme.upstream.toml",
        "original_sha256": _sha(original),
        "working_file": "scheme.published.toml",
        "working_sha256": _sha(corrected),
        "source_url": SOURCE_URL,
        "threshold_source": THRESHOLD_SOURCE,
        "transformations": transformations,
        "scientific_status": (
            "published_table_numeric_boundaries_match"
            if policy == "published-kpsc-629"
            else "published_boundary_discrepancy_remaining"
        ),
        "unresolved_boundaries": (
            []
            if policy == "published-kpsc-629"
            else [{"level7_max_and_level8_min": 99.3461, "published": 99.3641}]
        ),
        "basis": "round((629-7)/629*100, 4); published variant also round((629-4)/629*100, 4)",
        "assignment_policy_limitation": {
            "pinned_upstream_compare_max_missing_loci": 29,
            "published_max_missing_loci": 30,
            "status": "unchanged upstream comparison policy; numeric-boundary correction only",
        },
    }
    return corrected, receipt


def prepare_scheme(source: Path, destination: Path, *, policy="published-kpsc-629") -> dict | None:
    """Preserve upstream bytes and emit only an exact known correction."""
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve() or destination.name != "scheme.published.toml":
        raise SchemeValidationError(
            "Use a separate scheme.published.toml; historical source paths must remain unchanged"
        )
    original = source.read_bytes()
    if _sha(original) == ORIGINAL_SHA256:
        working, receipt = _correction(original, policy)
    else:
        validate_bytes(original, require_published=True)
        working, receipt = original, None
    destination.parent.mkdir(parents=True, exist_ok=True)
    upstream = destination.with_name("scheme.upstream.toml")
    if upstream.exists() and upstream.read_bytes() != original:
        raise SchemeValidationError(
            "Saved original plincer scheme differs; use a fresh reference directory"
        )
    receipt_path = destination.with_name("scheme.published.correction.json")
    if receipt is None and receipt_path.exists():
        raise SchemeValidationError("Unexpected correction receipt for an uncorrected scheme")
    if not upstream.exists():
        upstream.write_bytes(original)
    if destination.exists() and destination.read_bytes() != working:
        raise SchemeValidationError(
            "Existing working plincer scheme differs; use a fresh reference directory"
        )
    temporary = destination.with_name(destination.name + ".tmp")
    temporary.write_bytes(working)
    temporary.replace(destination)
    if receipt is not None:
        temporary_receipt = receipt_path.with_suffix(".json.tmp")
        temporary_receipt.write_text(json.dumps(receipt, indent=2) + "\n")
        temporary_receipt.replace(receipt_path)
    return receipt


def validate_scheme(path: Path) -> dict | None:
    """Reject malformed references; verify original/receipt for corrected bytes."""
    path = Path(path)
    working = path.read_bytes()
    validate_bytes(working, require_published=True)
    receipt_path = path.with_name("scheme.published.correction.json")
    known_working = {
        "cca91a9708ec4a0ad7cdc28249e2645bd949bbde638532c3d48a4ad3f4860eb3",
        "1f5d836a71d237e19be0e9c552fed293a62079ac685c717cc416e6ef484b088b",
    }
    if _sha(working) not in known_working and not receipt_path.exists():
        return None
    try:
        original = path.with_name("scheme.upstream.toml").read_bytes()
        actual_receipt = json.loads(receipt_path.read_text())
        expected, receipt = _correction(original, actual_receipt.get("policy"))
    except (OSError, ValueError) as error:
        raise SchemeValidationError(
            "Corrected plincer scheme lacks verified original/receipt"
        ) from error
    if path.name != receipt["working_file"] or working != expected or actual_receipt != receipt:
        raise SchemeValidationError("Corrected plincer scheme or transformation receipt changed")
    return receipt


def derive_config(source_config: Path, destination_config: Path) -> Path:
    """Derive a separate corrected descriptor; never build/rewrite old databases."""
    source_config, destination_config = (
        Path(source_config).resolve(),
        Path(destination_config).resolve(),
    )
    if source_config == destination_config:
        raise SchemeValidationError("Corrected configuration must have a new path")
    try:
        source_bytes = source_config.read_bytes()
        config = json.loads(source_bytes)
        organisms = config["organisms"]
        if not isinstance(organisms, dict):
            raise ValueError
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise SchemeValidationError("Cannot derive a corrected typing descriptor") from error
    changed = []
    # Preserve relative path meaning when the new descriptor is in another folder.
    for species, organism in organisms.items():
        for tool_name in ("cgmlst", "assignment"):
            tool = organism.get(tool_name, {})
            for field in (
                "index_dir",
                "scheme_file",
                "profiles_file",
                "alleles_db",
                "loci_file",
                "reference_db",
            ):
                if field in tool:
                    value = Path(tool[field]).expanduser()
                    tool[field] = str(
                        (source_config.parent / value).resolve()
                        if not value.is_absolute()
                        else value.resolve()
                    )
            if isinstance(tool.get("command"), list):
                tool["command"] = [
                    str((source_config.parent / arg).resolve())
                    if isinstance(arg, str) and (source_config.parent / arg).is_file()
                    else arg
                    for arg in tool["command"]
                ]
        assigner = organism.get("assignment", {})
        if assigner.get("kind") != "plincer":
            continue
        original = Path(assigner["scheme_file"])
        # A derivation is for the exact frozen original, not an arbitrary repair.
        if _sha(original.read_bytes()) != ORIGINAL_SHA256:
            raise SchemeValidationError("Derivation requires the exact pinned original scheme")
        published = original.with_name("scheme.published.toml")
        changed.append((assigner, original, published))
    if not changed:
        raise SchemeValidationError("Typing descriptor has no pinned plincer assignment to derive")
    # Validate target descriptor collision before creating any reference files.
    if destination_config.exists():
        raise SchemeValidationError(
            "Corrected configuration already exists; choose a fresh descriptor path"
        )
    for assigner, original, published in changed:
        prepare_scheme(original, published)
        validate_scheme(published)
        assigner["scheme_file"] = str(published)
    config["scheme_config_derivation"] = {
        "source_config": str(source_config),
        "source_config_sha256": _sha(source_bytes),
        "policy": "published-kpsc-629",
        "database_rebuild": False,
        "unchanged_resources": "caller index, profiles, allele databases and source configuration",
    }
    destination_config.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation protects existing descriptors, including historical runs.
    with destination_config.open("x") as handle:
        json.dump(config, handle, indent=2)
        handle.write("\n")
    return destination_config


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Derive a separate published-KpSC LIN descriptor without rebuilding databases"
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        output = derive_config(arguments.config, arguments.out)
    except (OSError, SchemeValidationError) as error:
        parser.exit(2, f"Scheme derivation failed: {error}\n")
    print(
        f"Derived {output}; existing scheme/config/databases unchanged; upstream missing-locus policy remains 29"
    )


if __name__ == "__main__":
    main()
