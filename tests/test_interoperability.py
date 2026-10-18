"""Downstream profile and tree exports preserve exact identities and metadata."""

import csv

import pytest
from Bio import Phylo

from chronoclade.interoperability import (
    write_grapetree_exports,
    write_microreact_exports,
)


def records():
    return [
        {
            "sample_id": "sample/a 1",
            "label": "Isolate A",
            "source_genome_id": "PW_001",
            "accession": "ERR001",
            "species": "Klebsiella pneumoniae",
            "collection_date": "2023-04",
            "country": "Greece",
            "region": "Attica",
            "city": "Athens",
            "host": "human",
            "hiercc_codes": {"HC10": "42", "HC5": "8"},
            "cgmlst_scheme": "scheme-x",
            "cgmlst_scheme_version": "v1",
            "cgmlst_profile": {"locus-a": "alpha", "locus-b": "7"},
        },
        {
            "sample_id": "sample:b-2",
            "label": "Isolate B",
            "species": "Klebsiella pneumoniae",
            "collection_date": None,
            "country": None,
            "host": "human",
            "cgmlst_scheme": "scheme-x",
            "cgmlst_scheme_version": "v1",
            "cgmlst_profile": {"locus-a": "beta", "locus-b": None},
        },
    ]


def test_grapetree_pair_has_shared_safe_ids_numeric_alleles_and_missing_zero(tmp_path):
    result = write_grapetree_exports(records(), tmp_path / "grapetree",
                                     loci=("locus-a", "locus-b", "locus-c"))
    with result["profiles"].open(encoding="utf-8", newline="") as handle:
        profile_rows = list(csv.reader(handle, delimiter="\t"))
    with result["metadata"].open(encoding="utf-8", newline="") as handle:
        metadata_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert profile_rows[0] == ["#Strain", "locus-a", "locus-b", "locus-c"]
    assert len(profile_rows) == 3
    assert all(row[0].startswith("CC_") for row in profile_rows[1:])
    assert [row[0] for row in profile_rows[1:]] == [row["ID"] for row in metadata_rows]
    assert all(call.isdigit() for row in profile_rows[1:] for call in row[1:])
    assert profile_rows[1][-1] == "0" and profile_rows[2][2:] == ["0", "0"]
    assert metadata_rows[0]["chronoclade_sample_id"] == "sample/a 1"
    assert metadata_rows[0]["sample_label"] == "ERR001"
    assert metadata_rows[0]["hiercc_codes"] == '{"HC10":"42","HC5":"8"}'
    assert metadata_rows[0]["city"] == "Athens"
    assert metadata_rows[1]["city"] == ""
    assert metadata_rows[1]["country"] == ""
    key_rows = list(csv.DictReader(result["allele_key"].open(encoding="utf-8", newline="")))
    assert {row["source_allele"] for row in key_rows} == {"alpha", "beta", "7"}
    assert "not establish complete scheme coverage" not in result["readme"].read_text()


def test_grapetree_inferred_loci_warn_about_incomplete_universe(tmp_path):
    result = write_grapetree_exports(records(), tmp_path / "inferred")
    assert "does not establish complete scheme coverage" in result["readme"].read_text()


def test_grapetree_rejects_mixed_scheme_and_repeated_identifiers(tmp_path):
    mixed = records()
    mixed[1]["cgmlst_scheme"] = "other"
    with pytest.raises(ValueError, match="one compatible"):
        write_grapetree_exports(mixed, tmp_path / "mixed")
    duplicate = records()
    duplicate[1]["sample_id"] = duplicate[0]["sample_id"]
    with pytest.raises(ValueError, match="duplicate sample_id"):
        write_grapetree_exports(duplicate, tmp_path / "duplicate")


def test_grapetree_rejects_different_database_and_locus_scope(tmp_path):
    rows = records()
    rows[0]["cgmlst_database_version"] = "db-1"
    rows[1]["cgmlst_database_version"] = "db-2"
    with pytest.raises(ValueError, match="scheme/database/locus scope"):
        write_grapetree_exports(rows, tmp_path / "version")
    rows[1]["cgmlst_database_version"] = "db-1"
    rows[1]["cgmlst_database_sha256"] = "a" * 64
    with pytest.raises(ValueError, match="scheme/database/locus scope"):
        write_grapetree_exports(rows, tmp_path / "fingerprint")


def test_grapetree_rejects_all_missing_calls_invalid_ids_and_metadata_fields(tmp_path):
    missing = records()
    missing[1]["cgmlst_profile"] = {"locus-a": "NA", "locus-b": "?"}
    with pytest.raises(ValueError, match=r"Records have no called cgMLST loci: \['sample:b-2'\]"):
        write_grapetree_exports(missing, tmp_path / "all-missing")
    assert not (tmp_path / "all-missing").exists()

    bad_id = records()
    bad_id[0]["sample_id"] = ""
    with pytest.raises(ValueError, match="nonempty sample_id"):
        write_grapetree_exports(bad_id, tmp_path / "bad-id")
    with pytest.raises(ValueError, match="Unsupported export metadata fields"):
        write_grapetree_exports(records(), tmp_path / "bad-field", metadata_fields=("patient_id",))

    occupied = tmp_path / "occupied"
    occupied.mkdir()
    sentinel = occupied / "keep.txt"
    sentinel.write_text("keep")
    with pytest.raises(FileExistsError, match="not empty"):
        write_grapetree_exports(records(), occupied)
    assert sentinel.read_text() == "keep"


def test_microreact_tree_and_metadata_share_safe_ids_and_keep_support(tmp_path):
    source_tree = tmp_path / "source.nwk"
    source_tree.write_text("(('sample/a 1':0.10000000000000002,'sample:b-2':0.2)95:0.3);")
    result = write_microreact_exports(source_tree, records(), tmp_path / "microreact")
    tree = Phylo.read(result["tree"], "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    with result["metadata"].open(encoding="utf-8", newline="") as handle:
        metadata = list(csv.DictReader(handle))
    assert all(ident.startswith("CC_") for ident in tips)
    assert tips == [row["ID"] for row in metadata]
    assert metadata[0]["chronoclade_sample_id"] == "sample/a 1"
    assert next(clade for clade in tree.find_clades() if clade.confidence is not None).confidence == 95
    assert tree.get_terminals()[0].branch_length == 0.10000000000000002
    assert set(row["sample_id"] for row in csv.DictReader(result["id_map"].open())) == {
        "sample/a 1", "sample:b-2"
    }


def test_microreact_rejects_missing_tree_metadata(tmp_path):
    tree = tmp_path / "source.nwk"
    tree.write_text("(unknown:0.1,sample:b-2:0.2);")
    with pytest.raises(ValueError, match="no matching metadata"):
        write_microreact_exports(tree, records(), tmp_path / "microreact")
