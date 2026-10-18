"""Strict offline comparison of native typing with frozen upstream output."""

import pytest

from typing_concordance_helpers import TypingConcordanceError, compare_typing_results


def assignment(**changes):
    row = {
        "sample_id": "sample-1",
        "species": "Klebsiella pneumoniae",
        "assembly_sha256": "a" * 64,
        "cgmlst_scheme": "klebsiella_1",
        "cgmlst_scheme_version": "2026-01-01",
        "cgmlst_database_version": "alleles-2026-01-01",
        "cgmlst_database_sha256": "b" * 64,
        "cgmlst_loci": ["locus-a", "locus-b", "locus-c", "locus-d"],
        "cgmlst_profile": {"locus-a": "1", "locus-d": "2"},
        "cgmlst_novel_alleles": {"locus-b": "c" * 40},
        "cgmlst_missing_loci": ["locus-c"],
        "cgmlst_ambiguous_loci": [],
        "cgmlst_locus_universe_complete": True,
        "cglin_scheme": "scgMLST629_S",
        "cglin_scheme_version": "2026-01-01",
        "cglin_database_version": "lin-2026-01-01",
        "cglin_database_sha256": "d" * 64,
        "cglin_raw": "1_2_3_4_5_6_7_8_9_10",
        "cglin_provisional": False,
    }
    row.update(changes)
    return row


def test_matching_frozen_typing_counts_known_novel_and_partial_lin_calls():
    result = compare_typing_results(
        [assignment()], [assignment(species="Klebsiella_pneumoniae")]
    )

    assert result["sample_status_counts"] == {"compared": 1}
    assert result["cgmlst"] == {
        "compared_samples": 1,
        "allele_matches": 3,
        "allele_discordances": 0,
        "callability_discordances": 0,
        "missing_vs_ambiguous_discordances": 0,
    }
    assert result["cglin"] == {"compared_samples": 1, "matches": 1, "discordances": 0}
    assert result["samples"][0]["cgmlst"]["native_called"] == 3
    assert result["samples"][0]["full_concordance_verified"]


def test_reports_locus_allele_and_callability_discordance_separately():
    native = assignment(
        cgmlst_profile={"locus-a": "9", "locus-d": "2"},
        cgmlst_loci=["locus-a", "locus-b", "locus-c", "locus-d", "locus-e"],
        cgmlst_novel_alleles={"locus-b": "c" * 40},
        cgmlst_missing_loci=["locus-c"],
        cgmlst_ambiguous_loci=["locus-e"],
    )
    frozen = assignment(
        cgmlst_profile={"locus-a": "1", "locus-c": "4", "locus-d": "2"},
        cgmlst_loci=["locus-a", "locus-b", "locus-c", "locus-d", "locus-e"],
        cgmlst_novel_alleles={"locus-b": "c" * 40},
        cgmlst_missing_loci=["locus-e"],
        cgmlst_ambiguous_loci=[],
        cglin_raw="1_2_3_4_*_*_*_*_*_*",
        cglin_provisional=True,
    )

    sample = compare_typing_results([native], [frozen])["samples"][0]
    cgmlst = sample["cgmlst"]
    assert cgmlst["allele_discordances"] == 1
    assert cgmlst["callability_discordances"] == 1
    assert cgmlst["missing_vs_ambiguous_discordances"] == 1
    assert [item["locus"] for item in cgmlst["disagreements"]] == ["locus-a", "locus-c", "locus-e"]
    assert sample["cglin"]["status"] == "discordant"


def test_does_not_compare_different_assemblies_or_reference_versions():
    wrong_assembly = assignment(assembly_sha256="e" * 64)
    result = compare_typing_results([assignment()], [wrong_assembly])
    assert result["samples"][0]["status"] == "identity_mismatch"
    assert "different_assembly_sha256" in result["samples"][0]["identity_problems"]

    wrong_version = assignment(cgmlst_scheme_version="2026-02-01")
    result = compare_typing_results([assignment()], [wrong_version])
    assert result["samples"][0]["cgmlst"]["status"] == "not_comparable"
    assert result["samples"][0]["cgmlst"]["scope_problems"] == ["different_scheme_version"]


def test_reference_fingerprint_mismatch_blocks_cgmlst_comparison():
    result = compare_typing_results(
        [assignment()], [assignment(cgmlst_database_sha256="f" * 64)]
    )
    assert result["samples"][0]["cgmlst"]["scope_status"] == "incompatible"
    assert result["samples"][0]["cgmlst"]["scope_problems"] == [
        "different_database_fingerprint"
    ]


def test_database_versions_must_be_present_on_both_sides_and_match():
    mismatch = compare_typing_results(
        [assignment()], [assignment(cgmlst_database_version="alleles-2026-02-01")]
    )["samples"][0]["cgmlst"]
    assert mismatch["status"] == "not_comparable"
    assert "different_database_version" in mismatch["scope_problems"]

    one_sided = assignment()
    one_sided.pop("cgmlst_database_version")
    missing = compare_typing_results([assignment()], [one_sided])["samples"][0]["cgmlst"]
    assert missing["status"] == "not_comparable"
    assert "missing_database_version_on_one_side" in missing["scope_problems"]


def test_unpinned_references_and_incomplete_locus_universes_cannot_verify_equivalence():
    native = assignment()
    frozen = assignment()
    for row in (native, frozen):
        row.pop("cgmlst_database_version")
        row.pop("cgmlst_database_sha256")
        row.pop("cglin_database_version")
        row.pop("cglin_database_sha256")
    unpinned = compare_typing_results([native], [frozen])["samples"][0]
    assert unpinned["cgmlst"]["status"] == "compared"
    assert unpinned["cgmlst"]["reference_pinning"] == "unverified"
    assert unpinned["full_equivalence_eligible"] is False
    assert unpinned["full_concordance_verified"] is False

    partial = assignment(cgmlst_locus_universe_complete=False)
    observed = compare_typing_results([partial], [assignment()])["samples"][0]
    assert observed["cgmlst"]["status"] == "compared"
    assert "canonical_locus_universe_unverified" in observed["cgmlst"]["verification_problems"]
    assert observed["full_concordance_verified"] is False


def test_reports_unpaired_samples_and_rejects_duplicate_ids():
    result = compare_typing_results([assignment()], [assignment(sample_id="sample-2")])
    assert result["sample_status_counts"] == {"missing_frozen": 1, "missing_native": 1}
    with pytest.raises(TypingConcordanceError, match="duplicate sample_id"):
        compare_typing_results([assignment(), assignment()], [])


def test_different_locus_universes_are_not_silently_intersected():
    frozen = assignment(
        cgmlst_loci=["locus-a", "locus-b", "locus-c"],
        cgmlst_profile={"locus-a": "1"},
        cgmlst_novel_alleles={"locus-b": "c" * 40},
        cgmlst_missing_loci=["locus-c"],
    )
    result = compare_typing_results([assignment()], [frozen])["samples"][0]["cgmlst"]
    assert result["status"] == "not_comparable"
    assert result["locus_scope_status"] == "incompatible"
    assert result["native_only_loci"] == ["locus-d"]
