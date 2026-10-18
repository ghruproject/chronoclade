"""Observed unions never assert a complete biological scheme or alter old calls."""

import pytest

from chronoclade.datasets import DatasetError, LocusCatalogue, from_profile_records
from chronoclade.prepare_provider import catalogues_for_records
from chronoclade.profile_analysis import _pair as compare_pair


def row(name, calls, *, version="v1"):
    return {
        "sample_id": name,
        "species": "Escherichia coli",
        "mlst_st": "131",
        "role": "input",
        "cgmlst_scheme": "synthetic-observed",
        "cgmlst_scheme_version": version,
        "cgmlst_profile": calls,
        "cgmlst_loci": list(calls),
        "cgmlst_locus_universe_complete": False,
    }


def test_2500_observed_loci_expand_without_mutating_old_calls_missingness_or_distances():
    loci = tuple(f"locus{i:04d}" for i in range(2500))
    calls = {locus: "1" for locus in loci}
    calls[loci[-1]] = None
    query = row("query", calls)
    prior = row("prior", {**calls, loci[1]: "2"})
    added = row("context", {**calls, "actually_observed_extra": "7"})
    catalogue = LocusCatalogue(
        "synthetic-observed", "v1", loci, "frozen query export", complete=False
    )
    baseline = from_profile_records([query, prior], [catalogue]).profiles[0]
    merged = catalogues_for_records([query, prior, added], [catalogue])
    assert merged[0].complete is False and merged[0].loci == (*loci, "actually_observed_extra")
    assert merged[0].database_version is None and merged[0].database_sha256 is None
    assert merged[0].source.startswith("Observed same-scope locus union; prior source:")
    assert catalogue.loci == loci and catalogue.source == "frozen query export"
    matrix = from_profile_records([query, prior, added], merged).profiles[0]
    for ident in ("query", "prior"):
        profile = matrix.profile(ident)
        assert {locus: profile[locus] for locus in loci} == baseline.profile(ident)
        assert profile["actually_observed_extra"] is None
    old_a = dict(query, cgmlst_profile=baseline.profile("query"))
    old_b = dict(prior, cgmlst_profile=baseline.profile("prior"))
    new_a = dict(query, cgmlst_profile=matrix.profile("query"), cgmlst_loci=list(merged[0].loci))
    new_b = dict(prior, cgmlst_profile=matrix.profile("prior"), cgmlst_loci=list(merged[0].loci))
    before, _ = compare_pair(old_a, old_b, 0.9)
    after, _ = compare_pair(new_a, new_b, 0.9)
    assert before["allele_differences"] == after["allele_differences"] == 1
    assert before["shared_called_loci"] == after["shared_called_loci"] == 2499
    assert before["distance"] == after["distance"]
    assert catalogues_for_records([query, prior, added], merged) == merged


def test_complete_declared_catalogue_stays_strict():
    catalogue = LocusCatalogue(
        "synthetic-observed", "v1", ("a",), "declared complete", complete=True
    )
    records = [row("query", {"a": "1"}), row("context", {"a": "1", "extra": "2"})]
    actual = catalogues_for_records(records, [catalogue])
    assert actual == (catalogue,)
    with pytest.raises(DatasetError, match="unsupported loci"):
        from_profile_records(records, actual)


@pytest.mark.parametrize(
    "field,value",
    [
        ("cgmlst_scheme_version", "v2"),
        ("cgmlst_database_version", "different-release"),
        ("cgmlst_database_sha256", "a" * 64),
    ],
)
def test_incompatible_scopes_never_extend_existing_observed_catalogue(field, value):
    catalogue = LocusCatalogue("synthetic-observed", "v1", ("a",), "observed", complete=False)
    first = row("query", {"a": "1"})
    other = dict(row("other", {"b": "2"}), **{field: value})
    result = catalogues_for_records([first, other], [catalogue])
    assert result[0] == catalogue and result[0].loci == ("a",)
    assert len(result) == 2 and result[1].loci == ("b",)
