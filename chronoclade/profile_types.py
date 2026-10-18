"""Scope-aware lineage identifiers and categorical allele-profile access."""

from __future__ import annotations

import re
from typing import Any, Mapping

from chronoclade.typing_scopes import typing_scope as _scope
from chronoclade.cglin import DEFAULT_SCHEME, parse_code



def _text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _lineage(row: Mapping, kind: str, level: int | str) -> tuple | None:
    from chronoclade.provider_hiercc import grouping_scope
    scope = grouping_scope(row) if kind == "hiercc" else None
    scope = scope or _scope(row, kind)
    if scope is None or _text(row.get(f"{kind}_status")).casefold() in {
        "conflict",
        "malformed",
        "unsupported",
        "missing",
        "unassigned",
        "failed",
    }:
        return None
    if kind == "cglin":
        components, status = parse_code(row.get("cglin_raw"))
        if scope[0] != DEFAULT_SCHEME or status not in {"complete", "partial"}:
            return None
        return (*scope, components[:level]) if len(components) >= level else None
    codes = row.get("hiercc_codes")
    if not isinstance(codes, Mapping):
        return None
    code = _text(codes.get(level))
    # A provisional/unknown cluster is not a resolved HierCC label.
    return (*scope, code) if code.isdigit() else None


def _profile(row: Mapping) -> Mapping | None:
    if row.get("cgmlst_status") in {"conflict", "failed", "unassigned", "unsupported"}:
        return None
    profile = row.get("cgmlst_profile")
    novel = row.get("cgmlst_novel_alleles") or {}
    if not isinstance(profile, Mapping) or not isinstance(novel, Mapping):
        return None
    if set(profile).intersection(novel):
        return None
    combined = dict(profile, **novel)
    return combined or None


def _loci(row: Mapping, profile: Mapping) -> set[str] | None:
    declared = row.get("cgmlst_loci")
    if declared is None:
        # Imported dense profiles must enumerate their complete scheme.
        return set(profile)
    if not isinstance(declared, (list, tuple)) or not declared:
        return None
    loci = set(declared)
    if len(loci) != len(declared) or not set(profile).issubset(loci):
        return None
    return loci


def _allele(value: Any) -> str | None:
    value = _text(value)
    # Official allele identifiers are positive integers; zero is an absent call.
    if re.fullmatch(r"[0-9a-fA-F]{40}", value):
        return value.lower()
    return str(int(value)) if value.isdigit() and int(value) > 0 else None
