"""Mandatory identity, interchangeable representatives and portable tree inputs."""

import csv
import json
import shutil

import pytest

from chronoclade.datasets import PreparedDataset, samples_from_records, write_dataset
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import (
    load_selection_ensemble,
    load_selection_manifest,
    write_selection_ensemble,
)
from chronoclade.selections import (
    DistanceEvidence,
    read_distance_evidence,
    representative_coverage,
    select_assembly_context,
    select_context_ensemble,
    selected_context_records,
)


def row(identifier, **kwargs):
    return {
        "sample_id": identifier,
        "origin": "context",
        "country": "UK",
        "collection_date": "2020",
        **kwargs,
    }


def analysis(tmp_path, records):
    path = tmp_path / "distances.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "sample_id_1",
                "sample_id_2",
                "distance",
                "allele_differences",
                "shared_called_loci",
            ],
        )
        writer.writeheader()
        writer.writerows(records)
    return {"paths": {"pairwise_distances": str(path)}}


def pair(a, b, raw, normalized):
    return {
        "sample_id_1": a,
        "sample_id_2": b,
        "allele_differences": raw,
        "distance": normalized,
        "shared_called_loci": 100,
    }


@pytest.mark.parametrize("distance_method", ["cgmlst", "alleles"])
def test_raw_allele_boundary_ties_union_with_pins_overrun(tmp_path, distance_method):
    queries = [row("q1"), row("q2")]
    context = [row(i) for i in ("a", "b", "c", "d", "pin", "far")]
    distances = analysis(
        tmp_path,
        [
            pair("q1", "a", 2, 0.02),
            pair("q1", "b", 2, 0.20),
            pair("q1", "c", 2, 0.01),
            pair("q1", "far", 3, 0.005),
            pair("q2", "d", 1, 0.01),
        ],
    )
    distances["distance_method"] = distance_method
    selections, ensemble = select_context_ensemble(
        queries,
        context,
        distances,
        size=1,
        nearest_per_query=1,
        include=["pin"],
        replicates=3,
    )
    first = selections[0]
    assert first["selected_context_ids"] == ["pin", "a", "d", "b", "c"]
    assert first["budget_overrun"] == 4
    assert first["mandatory_sample_ids"] == first["selected_sample_ids"]
    assert first["nearest_neighbours"][0]["boundary_tie_ids"] == ["a", "b", "c"]
    assert first["nearest_neighbours"][0]["tie_expanded"]
    assert ensemble["distinct_selections"] == 1
    assert selections[1]["duplicate_of"] == "selection-001"
    assert selections[1]["alternative_status"] == "duplicate"
    assert selections[1]["variation_status"] == "no_nonmandatory_choices"
    assert "far" not in first["selected_context_ids"]


def test_zero_budget_keeps_inputs_nearest_and_pins(tmp_path):
    distances = analysis(tmp_path, [pair("q", "a", 1, 0.01)])
    contexts, audit = select_assembly_context(
        [row("q")],
        [row("a"), row("pin")],
        distances,
        size=0,
        include=["pin"],
    )
    assert {r["sample_id"] for r in contexts} == {"a", "pin"}
    assert audit["budget_overrun"] == 2


def test_pinned_nearest_preserves_both_reasons(tmp_path):
    distances = analysis(tmp_path, [pair("q", "a", 1, 0.01)])
    selections, _ = select_context_ensemble([row("q")], [row("a")], distances, include=["a"])
    assert selections[0]["decisions"][0]["reasons"] == [
        "user_requested",
        "cgmlst_neighbour:q:rank=1",
    ]


def test_alternatives_change_only_optional_ids_and_freeze_quotas():
    context = [row(f"a{i}") for i in range(12)]
    distances = {"distance_evidence": DistanceEvidence((pair("q", "a0", 0, 0),))}
    selections, ensemble = select_context_ensemble(
        [row("q")],
        context,
        distances,
        size=4,
        replicates=5,
        seed=73,
    )
    repeated = select_context_ensemble(
        [row("q")], list(reversed(context)), distances, size=4, replicates=5, seed=73
    )
    assert (selections, ensemble) == repeated
    assert ensemble["distinct_selections"] == 5
    for selection in selections:
        assert selection["mandatory_sample_ids"] == ["q", "a0"]
        assert selection["quotas"] == selections[0]["quotas"]
        assert selection["strata_coverage"] == [
            {"stratum": ["unassigned", "2020", "UK", "Unknown"], "available": 12, "selected": 4}
        ]
    overlap = ensemble["pairwise_overlap"][0]
    left, right = (set(s["selected_sample_ids"]) for s in selections[:2])
    assert overlap["shared_ids"] == sorted(left & right)
    assert overlap["jaccard"] == len(left & right) / len(left | right)


def test_missing_dates_are_retained_and_reported():
    selections, _ = select_context_ensemble([], [row("x", collection_date=None)], {}, size=1)
    assert selections[0]["selected_context_ids"] == ["x"]
    assert selections[0]["missing_metadata"]["collection_date"] == ["x"]
    assert selections[0]["strata_coverage"][0]["stratum"][1] == "Unknown"


def test_distance_adapter_labels_snp_without_cgmlst_policy():
    evidence = DistanceEvidence(
        ({"sample_id_1": "q", "sample_id_2": "x", "distance": 0.5},),
        method="snp",
        distance_definition="SNP differences",
    )
    selections, _ = select_context_ensemble([row("q")], [row("x")], {"distance_evidence": evidence})
    assert selections[0]["source_method"] == "snp"
    assert selections[0]["nearest_ranking"] == "supplied_distance"
    assert selections[0]["decisions"][0]["reason"] == "snp_neighbour:q:rank=1"


def test_direct_dna_selection_preserves_raw_snp_ties_and_callable_sites(tmp_path):
    records = [
        dict(pair("q", "a", 1, 0.0001), snp_differences=8, shared_called_sites=80000),
        dict(pair("q", "b", 3, 0.005), snp_differences=5, shared_called_sites=1000),
        dict(pair("q", "c", 2, 0.0005), snp_differences=5, shared_called_sites=10000),
        dict(pair("q", "unavailable", 0, 0.00001), snp_differences=None,
             shared_called_sites=0),
    ]
    distance_path = tmp_path / "direct.csv"
    with distance_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    analysis_data = {
        "distance_method": "direct-dna",
        "distance_definition": "SNP differences / jointly called nucleotide sites",
        "paths": {"pairwise_distances": str(distance_path)},
    }
    queries, context = [row("q", role="input")], [row(i) for i in ("a", "b", "c", "unavailable")]
    prepared = write_dataset(
        PreparedDataset(samples=samples_from_records(queries + context)), tmp_path / "prepared"
    )
    ensemble_path = write_selection_ensemble(
        queries, context, analysis_data, dataset_manifest=prepared,
        output=tmp_path / "selections", size=0, nearest_per_query=1,
    )
    ensemble = load_selection_ensemble(ensemble_path)
    selection_path = ensemble_path.parent / ensemble["selections"][0]["path"]
    selection = load_selection_manifest(selection_path)
    assert selection["source_method"] == "direct-dna"
    assert selection["source"]["method"] == "direct-dna"
    assert selection["nearest_ranking"] == "raw_snp_differences_when_available"
    assert selection["required_nearest_context_ids"] == ["b", "c"]
    assert selection["nearest_neighbours"][0]["boundary_distance"] == 5
    assert selection["nearest_neighbours"][0]["tie_expanded"]
    evidence = read_distance_evidence(analysis_data)
    _, coverage = representative_coverage(queries, context, ["b", "c"], evidence)
    observed = next(item for item in coverage if item["sample_id"] == "q")
    assert observed["ranking_metric"] == "snp_differences"
    assert observed["nearest_value"] == 5
    assert observed["nearest_selected_ids"] == ["b", "c"]
    assert [item["shared_called_sites"] for item in observed["nearest_tie_evidence"]] == [1000, 10000]
    assert [item["allele_differences"] for item in observed["nearest_tie_evidence"]] == [3, 2]


def test_incomplete_raw_evidence_never_mixes_raw_and_normalized_ranks(tmp_path):
    distances = analysis(tmp_path, [pair("q", "a", 2, 0.2), pair("q", "b", None, 0.001)])
    selections, _ = select_context_ensemble(
        [row("q")], [row("a"), row("b")], distances, size=0, nearest_per_query=1
    )
    assert selections[0]["required_nearest_context_ids"] == ["a"]


def test_selection_pool_and_evidence_fail_closed():
    with pytest.raises(WorkflowError, match="disjoint"):
        select_context_ensemble([row("a")], [row("a")], {})
    with pytest.raises(WorkflowError, match="uniquely"):
        select_context_ensemble(
            [], [row("a", accession="alias"), row("b", accession="alias")], {}, include=["alias"]
        )
    with pytest.raises(WorkflowError, match="Conflicting"):
        select_context_ensemble(
            [row("q")],
            [row("a")],
            {
                "distance_evidence": DistanceEvidence(
                    (pair("q", "a", 1, 0.01), pair("a", "q", 2, 0.02))
                )
            },
        )


@pytest.fixture
def bundle(tmp_path):
    queries, context = [row("q", origin="local")], [row("a"), row("b")]
    dataset = PreparedDataset(samples=samples_from_records(queries + context))
    dataset_manifest = write_dataset(dataset, tmp_path / "prepared")
    distances = analysis(tmp_path, [pair("q", "a", 1, 0.01), pair("q", "b", 2, 0.02)])
    manifest = write_selection_ensemble(
        queries,
        context,
        distances,
        dataset_manifest=dataset_manifest,
        output=tmp_path / "runs",
        replicates=2,
        size=1,
        nearest_per_query=1,
    )
    return manifest, queries, context, distances


def test_manifest_roundtrip_relative_portability_and_exact_tree_ids(bundle, tmp_path):
    manifest, queries, context, _ = bundle
    ensemble = load_selection_ensemble(manifest)
    selection_path = manifest.parent / ensemble["selections"][0]["path"]
    selection = load_selection_manifest(selection_path)
    assert selection["selected_sample_ids"] == ["q", "a"]
    assert selected_context_records(selection, queries, context)[0]["sample_id"] == "a"
    moved = tmp_path / "moved"
    moved.mkdir()
    shutil.copytree(tmp_path / "prepared", moved / "prepared")
    shutil.copytree(tmp_path / "runs", moved / "runs")
    shutil.rmtree(tmp_path / "prepared")
    shutil.rmtree(tmp_path / "runs")
    assert load_selection_ensemble(moved / "runs" / "ensemble.json") == ensemble
    assert selection["source"]["distance_evidence"]["sha256"]
    assert selection["source"]["dataset"]["sha256"]


def test_manifest_tampered_distances_and_dataset_rejected(bundle):
    manifest, _, _, _ = bundle
    evidence = manifest.parent / "distance_evidence.json"
    evidence.write_text("{}")
    with pytest.raises(WorkflowError, match="source.distance_evidence.sha256"):
        load_selection_ensemble(manifest)


def test_manifest_mandatory_and_order_validation(bundle):
    manifest, _, _, _ = bundle
    ensemble = load_selection_ensemble(manifest)
    selection_path = manifest.parent / ensemble["selections"][0]["path"]
    selection = load_selection_manifest(selection_path)
    selection["selected_sample_ids"] = ["a", "q"]
    selection_path.write_text(json.dumps(selection))
    with pytest.raises(WorkflowError, match="selected_sample_ids"):
        load_selection_manifest(selection_path)
    with pytest.raises(WorkflowError, match="selections\\[0\\].sha256"):
        load_selection_ensemble(manifest)


def test_selection_writer_does_not_overwrite_or_publish_unknown_pool(bundle):
    manifest, queries, context, distances = bundle
    dataset_manifest = manifest.parent.parent / "prepared" / "dataset.json"
    with pytest.raises(WorkflowError, match="already exists"):
        write_selection_ensemble(
            queries, context, distances, dataset_manifest=dataset_manifest, output=manifest.parent
        )
    output = manifest.parent.parent / "invalid"
    with pytest.raises(WorkflowError, match="outside source.dataset"):
        write_selection_ensemble(
            queries,
            context + [row("unknown")],
            distances,
            dataset_manifest=dataset_manifest,
            output=output,
        )
    assert not output.exists()


def test_partition_exact_block_identity_enforced(bundle):
    from chronoclade.artifacts import file_sha256, write_json

    manifest, queries, context, distances = bundle
    dataset_manifest = manifest.parent.parent / "prepared" / "dataset.json"
    partition_path = manifest.parent.parent / "partition.json"
    partition = {
        "schema": "chronoclade.cgmlst.partitions",
        "schema_version": 1,
        "dataset_id": json.loads(dataset_manifest.read_text())["dataset_id"],
        "dataset_sha256": file_sha256(dataset_manifest),
        "blocks": [{"block_id": "block-1", "sample_ids": ["q", "a", "b"]}],
    }
    write_json(partition_path, partition)
    output = manifest.parent.parent / "partitioned_selection"
    ensemble_path = write_selection_ensemble(
        queries,
        context,
        dict(distances, partition={"block_id": "block-1"}),
        dataset_manifest=dataset_manifest,
        partition_manifest=partition_path,
        output=output,
    )
    ensemble = load_selection_ensemble(ensemble_path)
    assert ensemble["source"]["partition_block_id"] == "block-1"
    subset_output = manifest.parent.parent / "wrong_subset"
    with pytest.raises(WorkflowError, match="exact source.partition"):
        write_selection_ensemble(
            queries,
            context[:1],
            dict(distances, partition={"block_id": "block-1"}),
            dataset_manifest=dataset_manifest,
            partition_manifest=partition_path,
            output=subset_output,
        )
    assert not subset_output.exists()


def test_distance_ids_validated_before_any_publication(bundle):
    manifest, queries, context, _ = bundle
    evidence = DistanceEvidence((pair("q", "outside", 1, 0.01),))
    output = manifest.parent.parent / "unknown_evidence"
    with pytest.raises(WorkflowError, match="Distance evidence.*outside"):
        write_selection_ensemble(
            queries,
            context,
            {"distance_evidence": evidence},
            dataset_manifest=manifest.parent.parent / "prepared" / "dataset.json",
            output=output,
        )
    assert not output.exists()


def test_overlap_evidence_and_dataset_hash_are_validated(bundle):
    manifest, _, _, _ = bundle
    ensemble = load_selection_ensemble(manifest)
    ensemble["pairwise_overlap"][0]["intersection"] = 0
    manifest.write_text(json.dumps(ensemble))
    with pytest.raises(WorkflowError, match="pairwise_overlap"):
        load_selection_ensemble(manifest)
    dataset_manifest = manifest.parent.parent / "prepared" / "dataset.json"
    dataset_manifest.write_text(dataset_manifest.read_text() + "\n")
    with pytest.raises(WorkflowError, match="source.dataset.sha256"):
        load_selection_ensemble(manifest)


def test_selection_date_strata_use_validated_intervals_not_raw_prefixes():
    context = [row('canonical',collection_date=None,date_start='2018-01-01',date_end='2020-12-31',date_precision='interval'),
               row('invalid',collection_date='09/01/2019'), row('future',collection_date='2099')]
    selected,_ = select_context_ensemble([],context,{},size=3)
    audit=selected[0]
    assert {tuple(row['stratum']) for row in audit['strata_coverage']} == {
        ('unassigned','2018..2020','UK','Unknown'), ('unassigned','Unknown','UK','Unknown')}
    assert audit['missing_metadata']['collection_date'] == ['future','invalid']


def test_full_pool_coverage_records_raw_distance_ties_and_noncomparable_samples(tmp_path):
    queries = [row("q", role="input")]
    context = [row(ident) for ident in ("x", "y", "z", "orphan")]
    records = (
        pair("q", "x", 1, 0.1) | {"shared_called_loci": 10, "scheme_loci": 100},
        pair("q", "y", 1, 0.01) | {"shared_called_loci": 100, "scheme_loci": 100},
        pair("q", "z", 3, 0.03) | {"shared_called_loci": 100, "scheme_loci": 100},
        pair("x", "z", 2, 0.02) | {"shared_called_loci": 100, "scheme_loci": 100},
        pair("y", "z", 2, 0.2) | {"shared_called_loci": 10, "scheme_loci": 100},
    )
    dataset = write_dataset(
        PreparedDataset(samples=samples_from_records(queries + context)), tmp_path / "prepared"
    )
    ensemble_path = write_selection_ensemble(
        queries,
        context,
        {"distance_evidence": DistanceEvidence(records)},
        dataset_manifest=dataset,
        output=tmp_path / "selection_bundle",
        size=0,
        nearest_per_query=1,
    )
    ensemble = load_selection_ensemble(ensemble_path)
    selection_path = ensemble_path.parent / ensemble["selections"][0]["path"]
    selection = load_selection_manifest(selection_path)
    coverage = selection["representative_coverage"]
    assert coverage["scope"] == "full_available_selection_pool"
    assert coverage["sample_count"] == 5
    assert coverage["artifact"]["sha256"]

    with (selection_path.parent / coverage["artifact"]["path"]).open(newline="") as stream:
        coverage_rows = {item["sample_id"]: item for item in csv.DictReader(stream)}
    assert set(coverage_rows) == {"q", "x", "y", "z", "orphan"}
    z = coverage_rows["z"]
    assert z["selected"] == "False"
    assert z["ranking_metric"] == "allele_differences"
    assert z["nearest_value"] == "2"
    assert json.loads(z["nearest_selected_ids"]) == ["x", "y"]
    ties = json.loads(z["nearest_tie_evidence"])
    assert [(item["sample_id"], item["allele_differences"], item["shared_called_loci"])
            for item in ties] == [("x", 2, 100), ("y", 2, 10)]
    assert coverage_rows["orphan"]["coverage_status"] == "no_comparable_selected_representative"

    coverage_path = selection_path.parent / coverage["artifact"]["path"]
    coverage_path.write_text(coverage_path.read_text() + "tampered\n")
    with pytest.raises(WorkflowError, match="representative_coverage.artifact.sha256"):
        load_selection_manifest(selection_path)


def test_full_pool_coverage_normalizes_csv_distance_fields(tmp_path):
    path = tmp_path / "string_distances.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "sample_id_1", "sample_id_2", "distance", "allele_differences",
                "shared_called_loci", "call_overlap",
            ],
        )
        writer.writeheader()
        writer.writerows([
            {"sample_id_1": "q", "sample_id_2": "r1", "distance": "0.1",
             "allele_differences": "2", "shared_called_loci": "20", "call_overlap": "0.2"},
            {"sample_id_1": "q", "sample_id_2": "r2", "distance": "0.05",
             "allele_differences": "2", "shared_called_loci": "10", "call_overlap": "0.1"},
        ])
    evidence = read_distance_evidence({"paths": {"pairwise_distances": path}})
    summary, rows = representative_coverage(
        [row("q")], [row("r1"), row("r2")], ["r1", "r2"], evidence,
        ranking_metric="raw_allele_differences_when_available",
    )
    q = next(item for item in rows if item["sample_id"] == "q")
    assert q["nearest_value"] == 2
    assert [item["shared_called_loci"] for item in q["nearest_tie_evidence"]] == [20, 10]
    assert [item["scheme_loci"] for item in q["nearest_tie_evidence"]] == [100, 100]
    assert [item["distance"] for item in q["nearest_tie_evidence"]] == [0.1, 0.05]
    assert summary["nearest_value_quantiles"]["median"] == 0


def test_binary_coverage_uses_raw_ties_and_keeps_unavailable_rows(tmp_path):
    from chronoclade.cgmlst.scale_distance import write_distances

    loci = ["a", "b", "c", "d"]
    definitions = [
        ("q", "input", [1, 1, 1, 1], "db-a"),
        ("x", "context", [2, 1, 1, 1], "db-a"),
        ("y", "context", [1, 2, 1, 1], "db-a"),
        ("unavailable", "context", [1, 1, 1, 2], "db-b"),
    ]
    records = [
        {
            "sample_id": ident,
            "role": role,
            "origin": "local" if role == "input" else "context",
            "country": "UK",
            "species": "Klebsiella pneumoniae",
            "cgmlst_scheme": "test",
            "cgmlst_scheme_version": "1",
            "cgmlst_database_sha256": database,
            "cgmlst_loci": loci,
            "cgmlst_profile": dict(zip(loci, calls)),
        }
        for ident, role, calls, database in definitions
    ]
    evidence, _, _ = write_distances(records, tmp_path / "binary", min_overlap=0.9)
    queries, context = records[:1], records[1:]
    summary, rows = representative_coverage(queries, context, ["x", "y"], evidence)
    by_id = {record["sample_id"]: record for record in rows}
    assert summary["ranking_metric"] == "allele_differences"
    assert by_id["q"]["nearest_value"] == 1
    assert by_id["q"]["nearest_selected_ids"] == ["x", "y"]
    assert [row["shared_called_loci"] for row in by_id["q"]["nearest_tie_evidence"]] == [4, 4]
    assert by_id["unavailable"]["coverage_status"] == "no_comparable_selected_representative"

    dataset = write_dataset(PreparedDataset(samples=samples_from_records(records)), tmp_path / "prepared")
    ensemble = write_selection_ensemble(
        queries,
        context,
        {"distance_evidence": evidence},
        dataset_manifest=dataset,
        output=tmp_path / "binary_selection",
        size=0,
        nearest_per_query=1,
    )
    selection_path = ensemble.parent / load_selection_ensemble(ensemble)["selections"][0]["path"]
    selection = load_selection_manifest(selection_path)
    assert selection["representative_coverage"]["sample_count"] == 4
    coverage_path = selection_path.parent / selection["representative_coverage"]["artifact"]["path"]
    with coverage_path.open(newline="") as stream:
        coverage_rows = {row["sample_id"]: row for row in csv.DictReader(stream)}
    assert {ident for ident, row in coverage_rows.items() if row["selected"] == "True"} == {"q", "x", "y"}
    assert coverage_rows["unavailable"]["coverage_status"] == "no_comparable_selected_representative"
    with selection_path.open() as stream:
        persisted = json.load(stream)
    assert persisted["representative_coverage"]["artifact"]["sha256"]


def test_binary_full_pool_coverage_keeps_only_linear_output_memory():
    from chronoclade.matrix_distances import BinaryDistanceEvidence

    count = 15_000
    ids = [f"s{index:05d}" for index in range(count)]
    evidence = object.__new__(BinaryDistanceEvidence)
    evidence.method = "cgmlst"
    evidence.distance_definition = "fraction of mismatching jointly called cgMLST loci"
    evidence.sample_ids = tuple(ids)
    lookups = []

    def get(left, right):
        lookups.append((left, right))
        return {
            "distance": 0.1,
            "allele_differences": 1,
            "shared_called_loci": 10,
            "scheme_loci": 10,
            "call_overlap": 1.0,
        }

    evidence.get = get
    rows = [dict(sample_id=ident) for ident in ids]
    summary, coverage = representative_coverage(rows[:1], rows[1:], [ids[0]], evidence)
    assert len(coverage) == count
    assert len(lookups) == count - 1
    assert summary["sample_count"] == count
    assert coverage[0]["coverage_status"] == "selected_identity_only"


def test_explicit_reserve_keeps_all_nearest_ties_and_optional_variation(tmp_path):
    queries = [row("q", role="input")]
    context = [row(i) for i in ("a", "b", "c", "far1", "far2", "far3")]
    evidence = {"distance_evidence": DistanceEvidence(tuple(pair("q", i, 1, 0.01) for i in "abc"))}
    selections, ensemble = select_context_ensemble(
        queries, context, evidence, size=1, nearest_per_query=1,
        representative_reserve=1, replicates=3,
    )
    assert ensemble["distinct_selections"] == 3
    for selection in selections:
        assert selection["required_nearest_context_ids"] == ["a", "b", "c"]
        assert selection["context_budget"] == 1
        assert selection["budget_overrun"] == 2
        assert selection["effective_context_budget"] == 4
        assert len(selection["selected_context_ids"]) == 4
    dataset = write_dataset(PreparedDataset(samples=samples_from_records(queries + context)), tmp_path / "prepared")
    path = write_selection_ensemble(queries, context, evidence, dataset_manifest=dataset,
                                   output=tmp_path / "reserved", size=1, nearest_per_query=1,
                                   representative_reserve=1)
    ensemble = load_selection_ensemble(path)
    selection_path = path.parent / ensemble["selections"][0]["path"]
    saved = load_selection_manifest(selection_path)
    saved["representative_reserve"] = 0
    selection_path.write_text(json.dumps(saved))
    with pytest.raises(WorkflowError, match="effective_context_budget"):
        load_selection_manifest(selection_path)


@pytest.mark.parametrize("reserve", [-1, True, 1.5])
def test_invalid_representative_reserve_rejected(reserve):
    with pytest.raises(WorkflowError, match="Representative reserve"):
        select_context_ensemble([], [], {}, representative_reserve=reserve)


def test_evidence_snapshot_is_independent_and_fallback_cleans_partial_file(tmp_path, monkeypatch):
    import chronoclade.selection_manifest as manifests

    source = tmp_path / "source.npy"
    source.write_bytes(b"original matrix data")
    destination = tmp_path / "snapshot.npy"
    method = manifests.copy_evidence_snapshot(source, destination)
    assert method in {"copy_on_write", "copy"}
    assert source.stat().st_ino != destination.stat().st_ino
    source.write_bytes(b"changed")
    assert destination.read_bytes() == b"original matrix data"

    def failed_clone(original, cloned):
        cloned.write_bytes(b"partial")
        return False

    monkeypatch.setattr(manifests, "_clone_file", failed_clone)
    fallback = tmp_path / "fallback.npy"
    assert manifests.copy_evidence_snapshot(source, fallback) == "copy"
    assert fallback.read_bytes() == b"changed"
    with pytest.raises(WorkflowError, match="already exists"):
        manifests.copy_evidence_snapshot(source, fallback)
