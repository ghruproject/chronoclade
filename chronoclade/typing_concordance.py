"""Compare native typing with frozen upstream assignments offline.

Both inputs use the ChronoClade assignment fields emitted by
``query_typing.json``. The utility requires exact sample and assembly identity,
scheme/version scope and locus universes before reporting biological
concordance. It never infers calls from nearby samples.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Any, Iterable, Mapping

from chronoclade.cglin import parse_code
from chronoclade.query_typing import QueryTypingError, normalise_typing_species


class TypingConcordanceError(ValueError):
    """Typing inputs cannot be compared safely."""


def _rows_by_id(rows: Iterable[Mapping[str, Any]], label: str) -> dict[str, Mapping[str, Any]]:
    result: dict[str, Mapping[str, Any]] = {}
    for row in rows:
        if not isinstance(row, Mapping):
            raise TypingConcordanceError(f"{label} assignments must be objects")
        ident = row.get("sample_id")
        if not isinstance(ident, str) or not ident or ident != ident.strip():
            raise TypingConcordanceError(f"{label} assignments require exact sample_id values")
        if ident in result:
            raise TypingConcordanceError(f"{label} assignments contain duplicate sample_id: {ident}")
        result[ident] = row
    return result


def _declared_version(value: Any) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return None
    version = value.strip()
    if version.casefold() in {"unknown", "unavailable", "none", "null", "n/a"}:
        return None
    return version


def _scope_status(
    native: Mapping[str, Any], frozen: Mapping[str, Any], kind: str
) -> tuple[str, list[str], str]:
    problems = []
    for field in ("scheme", "scheme_version"):
        key = f"{kind}_{field}"
        left, right = _declared_version(native.get(key)), _declared_version(frozen.get(key))
        if left is None or right is None:
            problems.append(f"missing_{field}")
        elif left != right:
            problems.append(f"different_{field}")
    versions = []
    for row, label in ((native, "native"), (frozen, "frozen")):
        values = {
            _declared_version(row[key])
            for key in (f"{kind}_database_version", f"{kind}_reference_version")
            if _declared_version(row.get(key)) is not None
        }
        if len(values) > 1:
            problems.append(f"conflicting_{label}_database_versions")
            versions.append(None)
        else:
            versions.append(next(iter(values), None))
    if (versions[0] is None) != (versions[1] is None):
        problems.append("missing_database_version_on_one_side")
    elif versions[0] is not None and versions[0] != versions[1]:
        problems.append("different_database_version")

    hash_key = f"{kind}_database_sha256"
    left_hash, right_hash = native.get(hash_key), frozen.get(hash_key)
    if left_hash is not None or right_hash is not None:
        def valid(value: Any) -> bool:
            return isinstance(value, str) and re.fullmatch(r"[0-9a-fA-F]{64}", value) is not None

        if (left_hash is None) != (right_hash is None):
            problems.append("missing_database_fingerprint_on_one_side")
        elif not valid(left_hash) or not valid(right_hash):
            problems.append("missing_or_invalid_database_fingerprint")
        elif left_hash.lower() != right_hash.lower():
            problems.append("different_database_fingerprint")
    reference_scope_problems = {
        "conflicting_native_database_versions",
        "conflicting_frozen_database_versions",
        "missing_database_version_on_one_side",
        "different_database_version",
        "missing_database_fingerprint_on_one_side",
        "missing_or_invalid_database_fingerprint",
        "different_database_fingerprint",
    }
    if reference_scope_problems.intersection(problems):
        reference_pinning = "incompatible"
    elif versions[0] is not None or (left_hash is not None and right_hash is not None):
        reference_pinning = "verified"
    else:
        reference_pinning = "unverified"
    return ("compatible" if not problems else "incompatible"), problems, reference_pinning


def _allele(value: Any, locus: str) -> str:
    if isinstance(value, int) and not isinstance(value, bool) and value > 0:
        return str(value)
    if isinstance(value, str):
        text = value.strip()
        if text.isdigit() and int(text) > 0:
            return str(int(text))
        if re.fullmatch(r"[0-9a-fA-F]{40}", text):
            return text.lower()
    raise TypingConcordanceError(f"Invalid cgMLST allele value at locus {locus}")


def _cgmlst_calls(row: Mapping[str, Any], label: str) -> tuple[list[str], dict[str, tuple[str, str | None]]]:
    loci = row.get("cgmlst_loci")
    profile = row.get("cgmlst_profile")
    novel = row.get("cgmlst_novel_alleles") or {}
    missing = row.get("cgmlst_missing_loci") or []
    ambiguous = row.get("cgmlst_ambiguous_loci") or []
    if (
        not isinstance(loci, list)
        or not loci
        or any(not isinstance(locus, str) or not locus for locus in loci)
        or len(set(loci)) != len(loci)
    ):
        raise TypingConcordanceError(f"{label} cgMLST requires a unique complete locus universe")
    if not isinstance(profile, Mapping) or not isinstance(novel, Mapping):
        raise TypingConcordanceError(f"{label} cgMLST calls must be mappings")
    if not isinstance(missing, list) or not isinstance(ambiguous, list):
        raise TypingConcordanceError(f"{label} missing and ambiguous loci must be lists")
    called = set(profile) | set(novel)
    if (
        set(profile) - set(loci)
        or set(novel) - set(loci)
        or set(missing) - set(loci)
        or set(ambiguous) - set(loci)
        or set(profile) & set(novel)
        or called & (set(missing) | set(ambiguous))
        or set(missing) & set(ambiguous)
    ):
        raise TypingConcordanceError(f"{label} cgMLST calls overlap or name undeclared loci")
    calls = {
        locus: ("called", _allele(allele, locus)) for locus, allele in profile.items()
    }
    calls.update({locus: ("called", _allele(allele, locus)) for locus, allele in novel.items()})
    calls.update({locus: ("missing", None) for locus in missing})
    calls.update({locus: ("ambiguous", None) for locus in ambiguous})
    # An omitted locus in a complete scheme universe has no usable call.
    for locus in loci:
        calls.setdefault(locus, ("missing", None))
    return loci, calls


def _compare_cgmlst(native: Mapping[str, Any], frozen: Mapping[str, Any]) -> dict[str, Any]:
    scope, problems, reference_pinning = _scope_status(native, frozen, "cgmlst")
    result: dict[str, Any] = {
        "scope_status": scope,
        "scope_problems": problems,
        "reference_pinning": reference_pinning,
    }
    if scope != "compatible":
        result["status"] = "not_comparable"
        return result
    native_loci, native_calls = _cgmlst_calls(native, "native")
    frozen_loci, frozen_calls = _cgmlst_calls(frozen, "frozen")
    if set(native_loci) != set(frozen_loci):
        result.update(
            status="not_comparable",
            locus_scope_status="incompatible",
            native_only_loci=sorted(set(native_loci) - set(frozen_loci)),
            frozen_only_loci=sorted(set(frozen_loci) - set(native_loci)),
        )
        return result

    disagreements = []
    allele_matches = allele_discordances = callability_discordances = 0
    status_matches = status_discordances = 0
    for locus in native_loci:
        left, right = native_calls[locus], frozen_calls[locus]
        if left == right:
            status_matches += 1
            if left[0] == "called":
                allele_matches += 1
            continue
        if (left[0] == "called") != (right[0] == "called"):
            callability_discordances += 1
        elif left[0] == "called":
            allele_discordances += 1
        else:
            status_discordances += 1
        disagreements.append(
            {
                "locus": locus,
                "native_status": left[0],
                "native_allele": left[1],
                "frozen_status": right[0],
                "frozen_allele": right[1],
            }
        )
    universe_complete = (
        native.get("cgmlst_locus_universe_complete") is True
        and frozen.get("cgmlst_locus_universe_complete") is True
    )
    verification_problems = []
    if reference_pinning != "verified":
        verification_problems.append("reference_pinning_unverified")
    if not universe_complete:
        verification_problems.append("canonical_locus_universe_unverified")
    result.update(
        status="compared",
        verification_status="verified" if not verification_problems else "unverified",
        verification_problems=verification_problems,
        full_scheme_equivalence_eligible=not verification_problems,
        locus_scope_status="compatible",
        locus_count=len(native_loci),
        native_called=sum(status == "called" for status, _ in native_calls.values()),
        frozen_called=sum(status == "called" for status, _ in frozen_calls.values()),
        allele_matches=allele_matches,
        allele_discordances=allele_discordances,
        callability_discordances=callability_discordances,
        missing_vs_ambiguous_discordances=status_discordances,
        disagreements=disagreements,
    )
    return result


def _cglin_value(row: Mapping[str, Any], label: str) -> tuple[tuple[int, ...], str, bool]:
    raw = row.get("cglin_raw", "")
    components, status = parse_code(raw)
    if status == "malformed" or status == "unsupported":
        raise TypingConcordanceError(f"{label} cgLIN code is malformed or unsupported")
    provisional = row.get("cglin_provisional", False)
    if not isinstance(provisional, bool):
        provisional = str(provisional).casefold() in {"true", "1", "yes"}
    return components, status, provisional


def _compare_cglin(native: Mapping[str, Any], frozen: Mapping[str, Any]) -> dict[str, Any]:
    if not any(key in native for key in ("cglin_raw", "cglin_scheme", "cglin_scheme_version")) and not any(
        key in frozen for key in ("cglin_raw", "cglin_scheme", "cglin_scheme_version")
    ):
        return {"status": "not_present"}
    scope, problems, reference_pinning = _scope_status(native, frozen, "cglin")
    result = {
        "scope_status": scope,
        "scope_problems": problems,
        "reference_pinning": reference_pinning,
    }
    if scope != "compatible":
        result["status"] = "not_comparable"
        return result
    native_value = _cglin_value(native, "native")
    frozen_value = _cglin_value(frozen, "frozen")
    result.update(
        status="match" if native_value == frozen_value else "discordant",
        verification_status="verified" if reference_pinning == "verified" else "unverified",
        full_scheme_equivalence_eligible=reference_pinning == "verified",
        native_components=list(native_value[0]),
        native_code_status=native_value[1],
        native_provisional=native_value[2],
        frozen_components=list(frozen_value[0]),
        frozen_code_status=frozen_value[1],
        frozen_provisional=frozen_value[2],
    )
    return result


def compare_typing_results(
    native_rows: Iterable[Mapping[str, Any]], frozen_rows: Iterable[Mapping[str, Any]]
) -> dict[str, Any]:
    """Compare two offline assignment sets using exact IDs and assembly hashes.

    Rows must carry ``sample_id``, ``species`` and ``assembly_sha256``. cgMLST
    rows also require a complete ``cgmlst_loci`` universe. Database fingerprints
    are checked when supplied by either side; scheme and version are always
    mandatory for each typing method being compared. Full cgMLST equivalence is
    verified only when both rows declare a complete canonical locus universe
    and a matching database version or fingerprint.
    """
    native = _rows_by_id(native_rows, "native")
    frozen = _rows_by_id(frozen_rows, "frozen")
    ids = sorted(set(native) | set(frozen))
    samples = []
    for ident in ids:
        left, right = native.get(ident), frozen.get(ident)
        if left is None or right is None:
            samples.append(
                {
                    "sample_id": ident,
                    "status": "missing_native" if left is None else "missing_frozen",
                    "cgmlst": {"status": "not_comparable"},
                    "cglin": {"status": "not_comparable"},
                }
            )
            continue
        sample_problems = []
        for field in ("assembly_sha256",):
            a, b = left.get(field), right.get(field)
            if not isinstance(a, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", a):
                sample_problems.append(f"invalid_native_{field}")
            if not isinstance(b, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", b):
                sample_problems.append(f"invalid_frozen_{field}")
            elif isinstance(a, str) and a.lower() != b.lower():
                sample_problems.append("different_assembly_sha256")
        try:
            native_species = normalise_typing_species(left.get("species"))
            frozen_species = normalise_typing_species(right.get("species"))
        except QueryTypingError:
            native_species = frozen_species = None
        if native_species is None or native_species != frozen_species:
            sample_problems.append("different_species")
        if sample_problems:
            samples.append(
                {
                    "sample_id": ident,
                    "status": "identity_mismatch",
                    "identity_problems": sample_problems,
                    "cgmlst": {"status": "not_comparable"},
                    "cglin": {"status": "not_comparable"},
                }
            )
            continue
        samples.append(
            {
                "sample_id": ident,
                "status": "compared",
                "assembly_sha256": left["assembly_sha256"].lower(),
                "cgmlst": _compare_cgmlst(left, right),
                "cglin": _compare_cglin(left, right),
            }
        )
    for sample in samples:
        if sample["status"] != "compared":
            sample["full_equivalence_eligible"] = False
            sample["full_concordance_verified"] = False
            continue
        cgmlst_row = sample["cgmlst"]
        cglin_row = sample["cglin"]
        lineage_eligible = cglin_row.get("status") == "not_present" or (
            cglin_row.get("full_scheme_equivalence_eligible") is True
        )
        eligible = (
            cgmlst_row.get("full_scheme_equivalence_eligible") is True and lineage_eligible
        )
        sample["full_equivalence_eligible"] = eligible
        sample["full_concordance_verified"] = bool(
            eligible
            and cgmlst_row["allele_discordances"] == 0
            and cgmlst_row["callability_discordances"] == 0
            and cgmlst_row["missing_vs_ambiguous_discordances"] == 0
            and cglin_row.get("status") in {"not_present", "match"}
        )
    summary = Counter(sample["status"] for sample in samples)
    cgmlst = [sample["cgmlst"] for sample in samples if sample["cgmlst"].get("status") == "compared"]
    cglin = [sample["cglin"] for sample in samples if sample["cglin"].get("status") in {"match", "discordant"}]
    return {
        "schema": "chronoclade.typing_concordance.v1",
        "sample_count": len(samples),
        "sample_status_counts": dict(sorted(summary.items())),
        "fully_verified_sample_count": sum(
            sample["full_concordance_verified"] for sample in samples
        ),
        "full_equivalence_eligible_sample_count": sum(
            sample["full_equivalence_eligible"] for sample in samples
        ),
        "cgmlst": {
            "compared_samples": len(cgmlst),
            "allele_matches": sum(row["allele_matches"] for row in cgmlst),
            "allele_discordances": sum(row["allele_discordances"] for row in cgmlst),
            "callability_discordances": sum(row["callability_discordances"] for row in cgmlst),
            "missing_vs_ambiguous_discordances": sum(
                row["missing_vs_ambiguous_discordances"] for row in cgmlst
            ),
        },
        "cglin": {
            "compared_samples": len(cglin),
            "matches": sum(row["status"] == "match" for row in cglin),
            "discordances": sum(row["status"] == "discordant" for row in cglin),
        },
        "samples": samples,
    }
