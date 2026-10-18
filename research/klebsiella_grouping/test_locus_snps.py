"""Focused checks for native call extraction and independent sequence arithmetic."""

import hashlib
from types import SimpleNamespace

import pytest

from research.klebsiella_grouping.locus_snps import called_dna, verify_concatenated


@pytest.mark.parametrize(("reverse", "expected"), [(False, "AACG"), (True, "CGTT")])
def test_extract_inclusive_coordinates_and_verify_oriented_lowercase_checksum(reverse, expected):
    hit = {"contig": "contig full description", "contigStart": 3, "contigEnd": 6,
           "reverse": reverse, "checksum": hashlib.sha1(expected.lower().encode()).hexdigest()}
    assert called_dna(hit, {"contig full description": "TTAACGAA"}) == expected
    hit["checksum"] = "0" * 40
    with pytest.raises(ValueError, match="native SHA1"):
        called_dna(hit, {"contig full description": "TTAACGAA"})


def test_independent_concatenation_retains_missing_and_length_exclusions(tmp_path):
    records = {"a": {"cgmlst_profile": {"one": "1", "two": "1", "three": "1"}},
               "b": {"cgmlst_profile": {"one": "2", "two": "2"}}}
    sequences = {("one", "1"): "ACG", ("one", "2"): "ATG", ("two", "1"): "AA",
                 ("two", "2"): "A", ("three", "1"): "AT"}
    comparator = SimpleNamespace(sequences=sequences)
    by_locus = {"one": [("1", "ACG"), ("2", "ATG")],
                "two": [("1", "AA"), ("2", "A")], "three": [("1", "AT")]}
    pairs = [{"sample_a": "a", "sample_b": "b", "locus_snps": 1, "compared_sites": 3}]
    result = verify_concatenated(records, comparator, pairs, by_locus, tmp_path)
    assert result["verified_pairs"] == 1
    assert result["fixed_length_loci"] == 2
    assert result["length_variable_loci"] == ["two"]
    pairs[0]["locus_snps"] = 2
    with pytest.raises(ValueError, match="differs from DNAComparator"):
        verify_concatenated(records, comparator, pairs, by_locus, tmp_path)
