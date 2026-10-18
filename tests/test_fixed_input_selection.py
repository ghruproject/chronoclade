"""Independent all-input inclusion must not masquerade as distance selection."""

import json

import pytest

from chronoclade.artifacts import file_sha256
from chronoclade.datasets import PreparedDataset, samples_from_records, write_dataset
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import load_selection_manifest, write_fixed_input_selection
from chronoclade.tree_stage.workflow import selection_rows


def fixture(tmp_path, ids=("a", "b", "c")):
    rows = [
        dict(
            sample_id=ident,
            origin="local",
            species="Klebsiella pneumoniae",
            mlst_st=st,
            collection_date="2010",
        )
        for ident, st in zip(ids, ("258", "512", "1199"), strict=True)
    ]
    rows.append(
        dict(sample_id="outside", origin="context", species="Klebsiella pneumoniae", mlst_st="258")
    )
    dataset = write_dataset(
        PreparedDataset(samples=samples_from_records(rows)), tmp_path / "prepared"
    )
    metadata = tmp_path / "published.csv"
    metadata.write_text("sample_id\na\nb\nc\n")
    inclusion = tmp_path / "inclusion.json"
    inclusion.write_text(
        json.dumps(
            {
                "schema": "chronoclade.fixed-input-inclusion",
                "schema_version": 1,
                "selection_policy": "fixed_all_inputs",
                "sample_ids": list(ids),
                "reason": "Independently frozen complete published study panel",
                "source_metadata": {"path": metadata.name, "sha256": file_sha256(metadata)},
            }
        )
    )
    return dataset, inclusion, metadata


def selection(tmp_path):
    dataset, inclusion, metadata = fixture(tmp_path)
    path = write_fixed_input_selection(dataset, inclusion, output=tmp_path / "selection")
    return path, dataset, inclusion, metadata


def mutate(path, change):
    value = json.loads(path.read_text())
    change(value)
    path.write_text(json.dumps(value))


def test_exact_all_input_mixed_st_policy_and_tree_handoff(tmp_path):
    path, dataset, _, _ = selection(tmp_path)
    value = load_selection_manifest(path)
    assert value["selected_sample_ids"] == ["a", "b", "c"]
    assert value["selected_context_ids"] == value["nearest_neighbours"] == []
    assert value["source"]["distance_evidence"]["status"] == "not_required"
    assert "distance_definition" not in value["source"]
    actual, rows = selection_rows(path, value)
    assert actual == dataset
    assert [row["mlst_st"] for row in rows] == ["258", "512", "1199"]


@pytest.mark.parametrize("ids", [["a", "b"], ["b", "a", "c"], ["a", "b", "c", "outside"]])
def test_missing_reordered_or_context_inclusion_rejected(tmp_path, ids):
    dataset, inclusion, _ = fixture(tmp_path)
    mutate(inclusion, lambda value: value.update(sample_ids=ids))
    with pytest.raises(WorkflowError, match="all prepared inputs"):
        write_fixed_input_selection(dataset, inclusion, output=tmp_path / "selection")
    assert not (tmp_path / "selection/selection.json").exists()


@pytest.mark.parametrize("field", ["selected_sample_ids", "mandatory_sample_ids", "query_ids"])
def test_missing_or_added_selected_input_rejected(tmp_path, field):
    path, *_ = selection(tmp_path)
    mutate(path, lambda value: value.update({field: ["a", "b"]}))
    with pytest.raises(WorkflowError, match="all-and-only"):
        load_selection_manifest(path)


@pytest.mark.parametrize(
    "field",
    [
        "selected_context_ids",
        "available_context_ids",
        "pinned_context_ids",
        "required_nearest_context_ids",
    ],
)
def test_no_context_ranking_can_be_smuggled_into_fixed_policy(tmp_path, field):
    path, *_ = selection(tmp_path)
    mutate(path, lambda value: value.update({field: ["outside"]}))
    with pytest.raises(WorkflowError, match="cannot select or rank"):
        load_selection_manifest(path)


def test_nearest_claim_rejected(tmp_path):
    path, *_ = selection(tmp_path)
    mutate(
        path, lambda value: value.update(nearest_neighbours=[{"query_id": "a", "required_ids": []}])
    )
    with pytest.raises(WorkflowError, match="nearest-neighbour"):
        load_selection_manifest(path)


@pytest.mark.parametrize("policy", ["unsupported", None])
def test_not_required_evidence_cannot_bypass_other_policies(tmp_path, policy):
    path, *_ = selection(tmp_path)
    mutate(path, lambda value: value.update(selection_policy=policy))
    with pytest.raises(WorkflowError):
        load_selection_manifest(path)


def test_independent_metadata_and_inclusion_checksums_are_binding(tmp_path):
    path, _, inclusion, metadata = selection(tmp_path)
    saved = metadata.read_text()
    metadata.write_text(saved + "extra\n")
    with pytest.raises(WorkflowError, match="sha256"):
        load_selection_manifest(path)
    metadata.write_text(saved)
    inclusion.write_text(inclusion.read_text() + "\n")
    with pytest.raises(WorkflowError, match="sha256"):
        load_selection_manifest(path)


def test_prepared_dataset_artifact_changes_rejected(tmp_path):
    path, dataset, *_ = selection(tmp_path)
    samples = dataset.parent / "samples.csv"
    samples.write_text(samples.read_text().replace("2010", "2011"))
    with pytest.raises(WorkflowError, match="source.dataset is invalid"):
        load_selection_manifest(path)


def test_existing_output_evidence_never_overwritten(tmp_path):
    path, dataset, inclusion, _ = selection(tmp_path)
    before = path.read_bytes()
    with pytest.raises(WorkflowError, match="immutable"):
        write_fixed_input_selection(dataset, inclusion, output=path.parent)
    assert path.read_bytes() == before


@pytest.mark.parametrize(
    "field", ["representative_coverage", "nearest_ranking", "nearest_per_query"]
)
def test_ranked_coverage_claims_rejected(tmp_path, field):
    path, *_ = selection(tmp_path)
    mutate(path, lambda value: value.update({field: {"status": "complete"}}))
    with pytest.raises(WorkflowError, match="distance ranking"):
        load_selection_manifest(path)
