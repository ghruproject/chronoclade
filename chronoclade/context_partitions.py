"""Explicit, disjoint lineage blocks; no adaptive narrowing or provider access."""

import hashlib
import json

from chronoclade.cglin import resolved_prefix
from chronoclade.datasets import PreparedDataset
from chronoclade.typing_scopes import typing_scope
from chronoclade.provider_hiercc import reported_metadata

LIN_LEVELS = tuple(range(5, 11))


def validate_lin_level(level, *, label="LIN level"):
    """Require an explicit supported reported-prefix depth, excluding booleans."""
    if isinstance(level, bool) or not isinstance(level, int) or level not in LIN_LEVELS:
        raise ValueError(f"{label} must be an integer between 5 and 10")
    return level


def profile_records(dataset: PreparedDataset, *, source_root=None) -> list[dict]:
    """Adapt a validated bundle to the existing categorical scientific engine."""
    records = {
        row["sample_id"]: dict(row, origin="local" if row["role"] == "input" else "context")
        for row in dataset.samples
    }
    for matrix in dataset.profiles:
        cat = matrix.catalogue
        for ident in matrix.sample_ids:
            row = records[ident]
            if "cgmlst_profile" in row:
                raise ValueError("cgmlst stage requires one profile scheme per sample")
            row.update(
                cgmlst_profile=matrix.profile(ident),
                cgmlst_loci=list(cat.loci),
                cgmlst_scheme=cat.scheme_id,
                cgmlst_scheme_version=cat.scheme_version,
                cgmlst_database_version=cat.database_version,
                cgmlst_database_sha256=cat.database_sha256,
                cgmlst_locus_universe_complete=cat.complete,
                cgmlst_status="resolved",
            )
    seen_lineages = set()
    for assignment in dataset.lineages:
        kind = assignment.get("kind")
        if kind not in {"cglin", "hiercc"}:
            raise ValueError("Lineage evidence requires kind cglin or hiercc")
        evidence = assignment.get("evidence") or {}
        if not isinstance(evidence, dict):
            raise ValueError("Lineage evidence must be a mapping")
        identity = (assignment["sample_id"], kind)
        if identity in seen_lineages:
            raise ValueError(
                f"Repeated lineage evidence for sample/kind {identity!r}; resolve conflicts before cgmlst"
            )
        seen_lineages.add(identity)
        row = records[assignment["sample_id"]]
        kind = assignment["kind"]
        row.update({key: value for key, value in evidence.items() if key.startswith(kind + "_")})
        row.update(
            {
                kind + "_scheme": assignment["scheme_id"],
                kind + "_scheme_version": assignment["scheme_version"],
                kind + "_database_version": assignment.get("database_version"),
                kind + "_status": assignment["resolution_status"],
                kind + ("_raw" if kind == "cglin" else "_codes"): assignment["assignment"],
            }
        )
    from chronoclade.provider_hiercc import bind_dataset
    return bind_dataset(list(records.values()), dataset, source_root)


def _group_key(row, level, hiercc_level):
    if row.get("_ambiguous_mlst_namespace"):
        return None
    kind = "hiercc" if hiercc_level else "cglin"
    if row.get(kind + "_status") in {
        "failed",
        "missing",
        "unassigned",
        "unsupported",
        "conflict",
        "malformed",
    }:
        return None
    from chronoclade.provider_hiercc import grouping_scope
    scope = grouping_scope(row) if kind == "hiercc" else None
    scope = scope or typing_scope(row, kind)
    if scope is None or not row.get("mlst_st") or not row.get("species"):
        return None
    scheme, version = scope
    if hiercc_level:
        code = (row.get("hiercc_codes") or {}).get(hiercc_level)
        if not str(code or "").isdigit():
            return None
        prefix = (str(code),)
    else:
        prefix = resolved_prefix(row, level)
        if prefix is None:
            return None
    species = " ".join(row["species"].replace("_", " ").casefold().split())
    # A known fingerprint and an unknown fingerprint do not share a block.
    return (
        species,
        str(row["mlst_st"]),
        kind,
        scheme,
        version,
        row.get(kind + "_database_version"),
        row.get(kind + "_database_sha256"),
        row.get("_comparison_mlst_scheme", row.get("mlst_scheme")),
        prefix,
    )


def mlst_comparison_records(records):
    """Retain raw ST metadata while resolving only unambiguous comparison namespaces."""
    # Missing scheme metadata is not a new scheme. Use the only declared comparison
    # namespace when unambiguous, while retaining the sample's actual null metadata.
    known = {}
    for row in records:
        species = " ".join((row.get("species") or "").replace("_", " ").casefold().split())
        scope = (species, str(row.get("mlst_st") or ""))
        if row.get("mlst_scheme"):
            known.setdefault(scope, set()).add(row["mlst_scheme"])
    prepared = []
    for row in records:
        species = " ".join((row.get("species") or "").replace("_", " ").casefold().split())
        choices = known.get((species, str(row.get("mlst_st") or "")), set())
        scheme = row.get("mlst_scheme")
        comparison = scheme or (next(iter(choices)) if len(choices) == 1 else None)
        prepared.append(
            dict(
                row,
                _comparison_mlst_scheme=comparison,
                _ambiguous_mlst_namespace=not scheme and len(choices) > 1,
            )
        )
    return prepared


def partition_records(records, *, lin_level=5, context_lin_level=None, hiercc_level=None):
    """Match every frozen context record inside an explicit input lineage block."""
    validate_lin_level(lin_level)
    if context_lin_level is not None:
        validate_lin_level(context_lin_level, label="Context LIN level")
        if hiercc_level is not None:
            raise ValueError("Context LIN scope cannot be combined with HierCC")
    if hiercc_level is not None and (
        not hiercc_level.startswith("HC") or not hiercc_level[2:].isdigit()
    ):
        raise ValueError("HierCC level must be an explicit HC code, such as HC10")
    prepared = mlst_comparison_records(records)
    original_records = prepared
    scope_keys = {
        _group_key(row, context_lin_level, None)
        for row in prepared if row["role"] == "input"
    } - {None} if context_lin_level is not None else None
    records = [row for row in prepared if row["role"] == "input" or scope_keys is None
               or _group_key(row, context_lin_level, None) in scope_keys]
    unresolved_scope = [
        {"sample_id": row["sample_id"], "reason": "missing_or_unresolved_context_scope_typing"}
        for row in prepared if row["role"] == "input" and context_lin_level is not None
        and _group_key(row, context_lin_level, None) is None
    ]
    scoped_ids = {row["sample_id"] for row in records}
    scope_excluded = sorted(row["sample_id"] for row in prepared if row["sample_id"] not in scoped_ids)
    groups, unresolved = {}, []
    for row in records:
        if row["role"] != "input":
            continue
        key = _group_key(row, lin_level, hiercc_level)
        if key is None:
            unresolved.append(
                {"sample_id": row["sample_id"], "reason": "missing_or_unresolved_partition_typing"}
            )
        else:
            groups.setdefault(key, [])
    for row in records:
        key = _group_key(row, lin_level, hiercc_level)
        if key in groups:
            groups[key].append(row)
    blocks = []
    for key, rows in sorted(groups.items(), key=lambda item: repr(item[0])):
        digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:12]
        prefix = list(key[-1])
        reported = reported_metadata(rows[0]) if key[2] == "hiercc" else {}
        blocks.append(
            {
                "block_id": f"{key[2]}_{digest}",
                "species": rows[0]["species"].replace("_", " "),
                "mlst_st": key[1],
                "mlst_scheme": key[7],
                "kind": key[2],
                "scheme": key[3],
                "scheme_version": rows[0].get("hiercc_scheme_version") if reported else key[4],
                **reported,
                "database_version": key[5],
                "database_sha256": key[6],
                "prefix": prefix,
                "level": hiercc_level or lin_level,
                "records": sorted(rows, key=lambda r: r["sample_id"]),
            }
        )
    counts = []
    if not hiercc_level:
        for level in LIN_LEVELS:
            keys = {_group_key(row, level, None) for row in original_records if row["role"] == "input"} - {
                None
            }
            for key in sorted(keys, key=repr):
                matches = [row for row in original_records if _group_key(row, level, None) == key]
                counts.append(
                    {
                        "level": level,
                        "species": key[0],
                        "mlst_st": key[1],
                        "mlst_scheme": key[7],
                        "scheme": key[3],
                        "scheme_version": key[4],
                        "prefix": list(key[-1]),
                        "input_count": sum(row["role"] == "input" for row in matches),
                        "context_count": sum(row["role"] == "context" for row in matches),
                    }
                )
    used = {row["sample_id"] for block in blocks for row in block["records"]}
    return blocks, {
        "policy": "explicit-lineage-depth",
        "mlst_namespace_policy": "Declared schemes remain separate; missing scheme metadata is admitted only when at most one declared scheme exists for the species/ST. Original metadata is retained.",
        "ambiguous_mlst_namespace_ids": sorted(
            row["sample_id"] for row in records if row.get("_ambiguous_mlst_namespace")
        ),
        "lin_level": lin_level,
        "context_lin_level": context_lin_level,
        "context_scope_counts": {
            "prepared_context": sum(row["role"] == "context" for row in original_records),
            "scope_context": sum(row["role"] == "context" for row in records),
            "outside_scope_context": len(scope_excluded),
        },
        "outside_context_scope_ids": scope_excluded,
        "unresolved_context_scope_inputs": unresolved_scope,
        "hiercc_level": hiercc_level,
        "unresolved_inputs": unresolved,
        "available_level_counts": counts,
        "unmatched_context_ids": sorted(
            row["sample_id"]
            for row in original_records
            if row["role"] == "context" and row["sample_id"] not in used
        ),
        "context_scope": (
            f"Union of input-associated reported LIN{context_lin_level} prefixes in the frozen bundle; "
            f"analysed in LIN{lin_level} blocks"
            if context_lin_level is not None else
            "all matching context already present in the frozen bundle"
        ),
        "live_context_retrieval": "not_requested",
    }
