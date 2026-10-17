"""Fail-closed resolution of exact public Pathogenwatch source metadata."""

from __future__ import annotations

import re

from chronoclade.pathogenwatch import PathogenwatchError, normalize_record


class PublicIdentityError(ValueError):
    """Exact public source details could not be verified for an assembly."""


def resolve_public_source_identity(source_id: str, row: dict, detail: dict) -> dict:
    """Validate one details response and return checksum-bound source metadata."""
    if (
        not isinstance(detail, dict)
        or detail.get("uuid") != source_id
        or detail.get("projectAccess") not in (None, "PUBLIC")
    ):
        raise PublicIdentityError("public_identity_unverified")

    organism_id = str(detail.get("organismId") or "")
    expected_organism = str(row.get("organism_id") or "")
    species = str(row.get("species") or "")
    detail_species = str(detail.get("species") or "")
    if (
        not organism_id
        or (expected_organism and organism_id != expected_organism)
        or (species and detail_species and species.casefold() != detail_species.casefold())
    ):
        raise PublicIdentityError("public_identity_scope_mismatch")

    expected_st = str(row.get("mlst_st") or "")
    scheme = row.get("mlst_scheme")
    if scheme not in {"mlst", "mlst2"}:
        scheme = next(
            (candidate for candidate in ("mlst", "mlst2")
             if expected_st and str(detail.get(candidate) or "") == expected_st),
            "mlst" if detail.get("mlst") is not None else "mlst2",
        )
    st = expected_st or str(detail.get(scheme) or "")
    if not st or str(detail.get(scheme) or "") != st:
        raise PublicIdentityError("public_identity_scope_mismatch")

    search = {
        "uuid": source_id,
        "projectAccess": "PUBLIC",
        "organism": detail_species or species,
    }
    try:
        normalized = normalize_record(
            search,
            detail,
            organism_id=organism_id,
            st=st,
            mlst_scheme=scheme,
        )
    except (PathogenwatchError, KeyError, TypeError, ValueError) as error:
        raise PublicIdentityError("public_identity_scope_mismatch") from error

    checksum = str(normalized.get("source_checksum") or "")
    if not re.fullmatch(r"[a-fA-F0-9]{40}", checksum):
        raise PublicIdentityError("public_identity_checksum_unavailable")
    return {
        key: normalized.get(key)
        for key in (
            "source_checksum",
            "source_length",
            "source_n50",
            "source_contigs",
            "numeric_source_id",
            "qc_pass",
            "source_qc",
            "organism_id",
            "mlst_scheme",
            "mlst_st",
            "species",
            "raw_details_sha256",
        )
    }
