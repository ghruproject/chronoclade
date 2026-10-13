"""Direct DNA distances are sequence evidence, not arithmetic allele values."""

import csv
import json

import pytest
from typer.testing import CliRunner

from chronoclade.artifacts import file_sha256
from chronoclade.cgmlst.dna import CATALOGUE_SCHEMA, DNAComparator
from chronoclade.cgmlst.workflow import run_cgmlst
from chronoclade.cli import app
from chronoclade.context_partitions import profile_records
from chronoclade.datasets import LocusCatalogue, from_profile_records, write_dataset
from chronoclade.selection_manifest import load_selection_manifest


def fixture(tmp_path, sequences=None):
    records = []
    for ident, role, profile in [
        ("q", "local", {"a": "100", "b": "1"}),
        ("c1", "context", {"a": "101", "b": "1"}),
        ("c2", "context", {"a": "9000", "b": "2"}),
    ]:
        records.append({
            "sample_id": ident, "origin": role, "species": "Klebsiella pneumoniae",
            "mlst_st": "39", "collection_date": "2020", "country": "Greece",
            "cgmlst_scheme": "test", "cgmlst_scheme_version": "1",
            "cgmlst_database_version": "v1", "cgmlst_database_sha256": "a" * 64,
            "cgmlst_profile": profile, "cgmlst_status": "resolved",
            "cglin_scheme": "scgMLST629_S", "cglin_scheme_version": "1",
            "cglin_raw": "0,0,1,0,1,0,0,0,0,0", "cglin_status": "complete",
        })
    dataset = from_profile_records(records, [LocusCatalogue(
        "test", "1", ("a", "b"), "frozen fixture", "v1", database_sha256="a" * 64)])
    dataset_path = write_dataset(dataset, tmp_path / "prepared")
    dna_csv = tmp_path / "dna.csv"
    with dna_csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["locus", "allele", "dna"])
        writer.writerows(sequences or [
            ("a", "100", "AAAA"), ("a", "101", "CCCC"), ("a", "9000", "AAAT"),
            ("b", "1", "GGGG"), ("b", "2", "GGGT"),
        ])
    catalogue = tmp_path / "allele_catalogue.json"
    catalogue.write_text(json.dumps({
        "schema": CATALOGUE_SCHEMA, "schema_version": 1, "status": "complete",
        "scope": {"species": "Klebsiella pneumoniae", "scheme_id": "test",
                  "scheme_version": "1", "database_version": "v1", "database_sha256": "a" * 64},
        "artifact": {"path": dna_csv.name, "sha256": file_sha256(dna_csv)},
    }))
    return dataset, dataset_path, catalogue


def test_direct_dna_full_workflow_changes_neighbours_units_and_preserves_selection(tmp_path):
    _, dataset, catalogue = fixture(tmp_path)
    alleles = run_cgmlst(dataset, tmp_path / "alleles", bootstrap_replicates=2,
                        nearest_per_query=1, context_size=0)
    direct = run_cgmlst(dataset, tmp_path / "dna", distance_method="direct-dna", distance_threshold=.02,
                       allele_catalogue=catalogue, bootstrap_replicates=2,
                       nearest_per_query=1, context_size=0)
    def analysis(result):
        manifest = json.loads(result.manifest_path.read_text())
        directory = result.manifest_path.parent / manifest["blocks"][0]["block_id"]
        return json.loads((directory / "profile_analysis.json").read_text()), directory
    a, _ = analysis(alleles)
    d, directory = analysis(direct)
    assert a["distance_method"] == "alleles"
    assert a["nearest_neighbours"][0]["context_id"] == "c1"
    assert d["nearest_neighbours"][0]["context_id"] == "c2"
    assert d["nearest_neighbours"][0]["snp_differences"] == 2
    assert d["nearest_neighbours"][0]["shared_called_sites"] == 8
    assert d["nearest_neighbours"][0]["distance"] == .25
    assert d["root_to_tip"][0]["distance_units"] == "cgMLST nucleotide mismatch fraction"
    assert "allele_conversion" not in d["root_to_tip"][0]
    assert d["cohorts"][0]["valid_bootstrap_replicates"] == 2
    network = d["location_network"][0]
    assert network["distance_method"] == "direct-dna"
    assert network["nearest_edges"][0]["snp_differences_min"] == 2
    assert network["nearest_edges"][0]["shared_sites_max"] == 8
    with (directory / "pairwise_distances.csv").open() as stream:
        pairs = list(csv.DictReader(stream))
    assert all(row["shared_called_sites"] == "8" for row in pairs)
    selection = load_selection_manifest(next((directory / "selection_bundle/selections").glob("selection-*.json")))
    assert selection["source"]["method"] == "direct-dna"
    assert selection["selected_context_ids"] == ["c2"]
    assert selection["nearest_ranking"] == "raw_snp_differences_when_available"
    assert "cgMLST nucleotide mismatch fraction" in (directory / "cohort_1_nj.svg").read_text()
    assert direct.report_path.exists()


def test_missing_sequences_and_unequal_lengths_are_audited_not_matching(tmp_path):
    dataset, _, catalogue = fixture(tmp_path, sequences=[
        ("a", "100", "AAAA"), ("a", "101", "AA"), ("a", "9000", "AAAT"),
        ("b", "1", "GGGG"),
    ])
    rows = {row["sample_id"]: row for row in profile_records(dataset)}
    model = DNAComparator(catalogue, list(rows.values()))
    evidence, reason = model.pair(rows["q"], rows["c1"], .9)
    assert evidence is None and reason == "insufficient_sequence_overlap"
    evidence, reason = model.pair(rows["q"], rows["c1"], .5)
    assert reason is None and evidence["snp_differences"] == 0
    assert evidence["shared_called_sites"] == 4
    assert evidence["shared_called_loci"] == 1
    assert evidence["unequal_length_loci"] == ["a"]
    evidence, _ = model.pair(rows["q"], rows["c2"], .5)
    assert evidence["missing_sequence_loci"] == ["b"]
    assert evidence["snp_differences"] == 1
    audit = json.loads(model.save_audit(tmp_path).read_text())
    assert audit["exclusion_counts"]["unequal_length_loci"] == 2
    assert audit["provenance"]["dna_csv_sha256"]


@pytest.mark.parametrize("change,match", [
    ({"database_sha256": "b" * 64}, "incompatible"),
    ({"species": "Escherichia coli"}, "species"),
    ({"database_sha256": None}, "known database"),
])
def test_catalogue_scope_fails_closed_before_output(tmp_path, change, match):
    _, dataset, catalogue = fixture(tmp_path)
    payload = json.loads(catalogue.read_text())
    payload["scope"].update(change)
    catalogue.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match=match):
        run_cgmlst(dataset, tmp_path / "out", distance_method="direct-dna", distance_threshold=.02, allele_catalogue=catalogue)
    assert not (tmp_path / "out").exists()


def test_reject_ambiguous_bases_and_checksum_mutation(tmp_path):
    dataset, _, catalogue = fixture(tmp_path, sequences=[("a", "100", "AANA")])
    with pytest.raises(ValueError, match="ambiguous"):
        DNAComparator(catalogue, profile_records(dataset))
    (tmp_path / "dna.csv").write_text("locus,allele,dna\na,100,AAAA\n")
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        DNAComparator(catalogue, profile_records(dataset))


def test_cli_default_and_direct_dna_required_catalogue_errors(tmp_path):
    _, dataset, catalogue = fixture(tmp_path)
    runner = CliRunner()
    help_result = runner.invoke(app, ["cgmlst", "--help"])
    assert help_result.exit_code == 0
    assert "alleles" in help_result.output and "direct-dna" in help_result.output
    result = runner.invoke(app, ["cgmlst", str(dataset), "--out", str(tmp_path / "out"),
                                 "--distance-method", "direct-dna"])
    assert result.exit_code == 2
    assert "requires --allele-catalogue" in result.output
    assert not (tmp_path / "out").exists()
    result = runner.invoke(app, ["cgmlst", str(dataset), "--out", str(tmp_path / "out"),
                                 "--distance-method", "unknown"])
    assert result.exit_code == 2


def test_direct_snp_ties_retained_despite_distinct_allele_counts(tmp_path):
    _, dataset, catalogue = fixture(tmp_path, sequences=[
        ("a", "100", "AAAA"), ("a", "101", "AATT"), ("a", "9000", "AAAT"),
        ("b", "1", "GGGG"), ("b", "2", "GGGT"),
    ])
    result = run_cgmlst(dataset, tmp_path / "out", distance_method="direct-dna", distance_threshold=.02,
                        allele_catalogue=catalogue, bootstrap_replicates=0,
                        nearest_per_query=1, context_size=0)
    manifest = json.loads(result.manifest_path.read_text())
    directory = result.manifest_path.parent / manifest["blocks"][0]["block_id"]
    analysis = json.loads((directory / "profile_analysis.json").read_text())
    assert {row["context_id"] for row in analysis["nearest_neighbours"]} == {"c1", "c2"}
    assert {row["tied_neighbours"] for row in analysis["nearest_neighbours"]} == {2}
    selection = load_selection_manifest(next((directory / "selection_bundle/selections").glob("selection-*.json")))
    assert set(selection["selected_context_ids"]) == {"c1", "c2"}


def test_locus_bootstrap_weights_actual_sites_not_mean_gene_fractions(tmp_path):
    import numpy as np
    dataset, _, catalogue = fixture(tmp_path, sequences=[
        ("a", "100", "AAAA"), ("a", "101", "CCCC"), ("a", "9000", "AAAT"),
        ("b", "1", "GGGGGGGG"), ("b", "2", "GGGGGGGT"),
    ])
    rows = {row["sample_id"]: row for row in profile_records(dataset)}
    model = DNAComparator(catalogue, list(rows.values()))
    matrix = model.resample([rows["q"], rows["c2"]], ["a", "b"], np.array([0, 0, 1]), 1)
    assert matrix[0, 1] == 3 / 16
    assert np.array_equal(np.diag(matrix), [0, 0])


def test_direct_dna_requires_explicit_method_threshold(tmp_path):
    _, dataset, catalogue = fixture(tmp_path)
    with pytest.raises(ValueError, match="explicit --distance-threshold"):
        run_cgmlst(dataset, tmp_path / "out", distance_method="direct-dna", allele_catalogue=catalogue)
    assert not (tmp_path / "out").exists()
    result = CliRunner().invoke(app, ["cgmlst", str(dataset), "--out", str(tmp_path / "out"),
                                     "--distance-method", "direct-dna", "--allele-catalogue", str(catalogue)])
    assert result.exit_code == 2 and "explicit --distance-threshold" in result.output


def test_missing_input_dna_yields_partial_stage_and_retains_input_selection(tmp_path):
    _, dataset, catalogue = fixture(tmp_path, sequences=[
        ("a", "101", "CCCC"), ("a", "9000", "AAAT"),
        ("b", "1", "GGGG"), ("b", "2", "GGGT"),
    ])
    result = run_cgmlst(dataset, tmp_path / "out", distance_method="direct-dna",
                        distance_threshold=.02, allele_catalogue=catalogue, bootstrap_replicates=0)
    manifest = json.loads(result.manifest_path.read_text())
    assert manifest["status"] == "partial"
    block = manifest["blocks"][0]
    assert block["analysed_input_count"] == 0 and block["input_count"] == 1
    directory = result.manifest_path.parent / block["block_id"]
    analysis = json.loads((directory / "profile_analysis.json").read_text())
    assert analysis["exclusions"] == [{"sample_id": "q", "reason": "insufficient_sequence_overlap"}]
    selection = load_selection_manifest(next((directory / "selection_bundle/selections").glob("selection-*.json")))
    assert "q" in selection["selected_sample_ids"]


def test_newick_preserves_nucleotide_branch_precision(tmp_path):
    from Bio import Phylo
    _, dataset, catalogue = fixture(tmp_path, sequences=[
        ("a", "100", "AAAAAAA"), ("a", "101", "CCCCCCC"), ("a", "9000", "AAAAAAT"),
        ("b", "1", "GGGGGGG"), ("b", "2", "GGGGGGT"),
    ])
    result = run_cgmlst(dataset, tmp_path / "out", distance_method="direct-dna",
                        distance_threshold=.02, allele_catalogue=catalogue, bootstrap_replicates=0)
    manifest = json.loads(result.manifest_path.read_text())
    directory = result.manifest_path.parent / manifest["blocks"][0]["block_id"]
    tree = Phylo.read(directory / "cohort_1_nj.nwk", "newick")
    assert tree.distance("q", "c2") == pytest.approx(2 / 14, rel=1e-11)


def test_small_snp_fractions_still_flag_non_euclidean_ordination(tmp_path):
    from chronoclade.profile_analysis import analyse_profiles
    prefix = "A" * 120000
    dataset, _, catalogue = fixture(tmp_path, sequences=[
        ("a", "100", prefix + "AAAA"), ("a", "101", prefix + "AATT"),
        ("a", "9000", prefix + "TTAA"), ("a", "105", prefix + "TTTT"),
        ("b", "1", "G" * 120000),
    ])
    rows = profile_records(dataset)
    rows[-1]["cgmlst_profile"]["b"] = "1"
    rows.append(dict(rows[1], sample_id="c3", cgmlst_profile={"a": "105", "b": "1"}))
    model = DNAComparator(catalogue, rows)
    analysis = analyse_profiles(rows, output=tmp_path / "ordination", distance_model=model,
                                bootstrap_replicates=0, min_overlap=1, distance_threshold=.000001)
    cohort = analysis["cohorts"][0]
    assert any("negative PCoA eigenvalues" in item for item in cohort["warnings"])
    assert min(cohort["pcoa_eigenvalues"]) > -1e-10
