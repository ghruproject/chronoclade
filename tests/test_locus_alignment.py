"""Per-locus MSA reconstruction, evidence accounting and fail-closed checks."""

import csv
import json
import sys

import pytest

from chronoclade.cgmlst.dna import CATALOGUE_SCHEMA
from chronoclade.cgmlst.locus_alignment import (
    MAFFT_OPTIONS, align_catalogue, file_hash, read_fasta, validate_msa,
)


@pytest.fixture
def panel(tmp_path):
    scope = {"species": "Klebsiella pneumoniae", "scheme_id": "test",
             "scheme_version": "1", "database_version": None, "database_sha256": "a" * 64}
    rows = []
    profiles = {"a": {"shift": "1", "fixed": "1"},
                "b": {"shift": "2", "fixed": "2"}, "c": {"fixed": "1"}}
    for ident, profile in profiles.items():
        rows.append({"sample_id": ident, "species": scope["species"], "cgmlst_scheme": "test",
                     "cgmlst_scheme_version": "1", "cgmlst_database_version": None,
                     "cgmlst_database_sha256": "a" * 64,
                     "cgmlst_loci": ["shift", "fixed", "absent"], "cgmlst_profile": profile})
    dna = tmp_path / "alleles.csv"
    with dna.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerows([("locus", "allele", "dna"), ("shift", "1", "AACC"),
                         ("shift", "2", "ACCA"), ("fixed", "1", "GG"), ("fixed", "2", "GG")])
    catalogue = tmp_path / "catalogue.json"
    catalogue.write_text(json.dumps({"schema": CATALOGUE_SCHEMA, "schema_version": 1,
        "status": "complete", "scope": scope,
        "artifact": {"path": "alleles.csv", "sha256": file_hash(dna)}}))
    tool = tmp_path / "mafft"
    tool.write_text(f"#!{sys.executable}\nimport sys\n"
                    "if '--version' in sys.argv:\n print('fixture-aligner 1');sys.exit()\n"
                    "print('>dna_0001\\nAACC-\\n>dna_0002\\n-ACCA')\n")
    tool.chmod(0o755)
    return catalogue, rows, tool


def test_equal_length_alleles_are_aligned_and_indels_remain_separate(panel, tmp_path):
    catalogue, rows, tool = panel
    result = align_catalogue(catalogue, rows, tmp_path / "out", tool)
    assert result.samples == {"a": "AACC-GG", "b": "-ACCAGG", "c": "NNNNNGG"}
    assert result.aligned_alleles["shift", "1"] == "AACC-"
    assert result.aligned_alleles["fixed", "1"] == result.aligned_alleles["fixed", "2"]
    assert result.partitions[0]["method"] == "mafft"
    assert result.partitions[1]["method"] == "identical-sequence-exact"
    assert result.partitions[2]["status"] == "absent"
    assert result.partitions[2]["width"] == 0
    assert result.partitions[2]["start"] is result.partitions[2]["end"] is None
    assert result.provenance["mafft"]["options"] == list(MAFFT_OPTIONS)
    pairs = list(result.pairs(min_loci=2))
    assert len(pairs) == 3
    ab = pairs[0]
    assert (ab["locus_snps"], ab["compared_sites"], ab["gap_vs_base_columns"], ab["indel_tracts"]) == (0, 5, 2, 2)
    assert ab["compared_loci"] == 2 and ab["eligible"]
    ac = pairs[1]
    assert (ac["locus_snps"], ac["compared_sites"], ac["missing_columns"]) == (0, 2, 5)
    assert ac["compared_loci"] == 1 and not ac["eligible"]
    assert all(p["independent_arithmetic_verified"] for p in pairs)


@pytest.mark.parametrize("failure", ["exit", "change", "duplicate", "truncate"])
def test_alignment_failures_never_fall_back_to_positional_comparison(panel, tmp_path, failure):
    catalogue, rows, tool = panel
    body = {"exit": "sys.exit(3)", "change": "print('>dna_0001\\nAATC-\\n>dna_0002\\n-ACCA')",
            "duplicate": "print('>dna_0001\\nAACC-\\n>dna_0001\\n-ACCA')",
            "truncate": "print('>dna_0001\\nAACC-')"}[failure]
    tool.write_text(f"#!{sys.executable}\nimport sys\n"
                    "if '--version' in sys.argv:\n print('fixture-aligner 1');sys.exit()\n" + body + "\n")
    output = tmp_path / "failed"
    output.mkdir()
    (output / "alignment.json").write_text('{"status":"complete"}')
    with pytest.raises(ValueError):
        align_catalogue(catalogue, rows, output, tool)
    assert not (output / "alignment.json").exists()


def test_unequal_lengths_retain_joint_bases_and_report_indels(panel, tmp_path):
    catalogue, rows, tool = panel
    dna = catalogue.parent / "alleles.csv"
    dna.write_text(dna.read_text().replace("ACCA", "ACC"))
    manifest = json.loads(catalogue.read_text())
    manifest["artifact"]["sha256"] = file_hash(dna)
    catalogue.write_text(json.dumps(manifest))
    tool.write_text(tool.read_text().replace("-ACCA", "-ACC-"))
    result = align_catalogue(catalogue, rows, tmp_path / "unequal", tool)
    ab = next(result.pairs(min_loci=2))
    assert ab["compared_loci"] == 2
    assert ab["compared_sites"] == 5
    assert ab["gap_vs_base_columns"] == ab["indel_tracts"] == 1
    assert ab["both_gap_columns"] == 1
    assert ab["locus_snps"] == 0


def test_msa_validation_rejects_bad_width_and_unsupported_bases():
    with pytest.raises(ValueError, match="unequal widths"):
        validate_msa({"a": "AC", "b": "A-C"}, {"a": "AC", "b": "AC"})
    with pytest.raises(ValueError, match="unsupported"):
        validate_msa({"a": "AN"}, {"a": "AC"})
    with pytest.raises(ValueError, match="precedes"):
        read_fasta("AC\n>a\nAC\n")
