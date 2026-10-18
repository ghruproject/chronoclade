"""Bind an independently declared complete input set without distance selection."""

from __future__ import annotations

import os
from pathlib import Path
import re

from chronoclade import __version__
from chronoclade.artifacts import write_json
from chronoclade.errors import WorkflowError

POLICY = "fixed_all_inputs"
INCLUSION_SCHEMA = "chronoclade.fixed-input-inclusion"
DISTANCE_REASON = (
    "All prepared inputs are retained; no distance-based context selection is performed."
)


def validate_source(manifest, base, dataset):
    """Validate independent inclusion evidence and absence of distance claims."""
    from chronoclade.selection_manifest import _ids, _read, _resolve

    source = manifest["source"]
    if source.get("method") != POLICY or source.get("partition") is not None:
        raise WorkflowError(
            "Fixed input selection requires its explicit policy without a lineage partition"
        )
    if source.get("distance_evidence") != {"status": "not_required", "reason": DISTANCE_REASON}:
        raise WorkflowError(
            "Fixed input distance evidence must be explicitly not_required with its reason"
        )
    path = _resolve(
        base, source.get("inclusion_evidence"), "source.inclusion_evidence", external=True
    )
    inclusion = _read(path)
    if (
        inclusion.get("schema") != INCLUSION_SCHEMA
        or type(inclusion.get("schema_version")) is not int
        or inclusion["schema_version"] != 1
        or inclusion.get("selection_policy") != POLICY
    ):
        raise WorkflowError("Expected versioned fixed-input inclusion evidence")
    if not isinstance(inclusion.get("reason"), str) or not inclusion["reason"].strip():
        raise WorkflowError("Independent fixed-input inclusion reason is required")
    _resolve(
        path.parent, inclusion.get("source_metadata"), "inclusion.source_metadata", external=True
    )
    expected = [row["sample_id"] for row in dataset.samples if row["role"] == "input"]
    if not expected or _ids(inclusion, "sample_ids") != expected:
        raise WorkflowError(
            "Inclusion evidence must bind all prepared inputs exactly in dataset order"
        )


def validate_selection(manifest, dataset):
    """All-and-only inputs: no context, ranking, neighbour or subset claims."""
    from chronoclade.selection_manifest import _ids

    expected = [row["sample_id"] for row in dataset.samples if row["role"] == "input"]
    if any(
        field in manifest
        for field in ("representative_coverage", "nearest_ranking", "nearest_per_query")
    ):
        raise WorkflowError(
            "Fixed input selection cannot claim distance ranking or representative coverage"
        )
    for field in ("query_ids", "selected_sample_ids", "mandatory_sample_ids"):
        if _ids(manifest, field) != expected:
            raise WorkflowError(
                "Fixed input selection must retain all-and-only inputs in dataset order"
            )
    for field in (
        "selected_context_ids",
        "available_context_ids",
        "pinned_context_ids",
        "required_nearest_context_ids",
    ):
        if _ids(manifest, field):
            raise WorkflowError("Fixed input selection cannot select or rank context genomes")
    if manifest.get("nearest_neighbours") != [] or manifest.get("decisions") != []:
        raise WorkflowError(
            "Fixed input selection cannot claim nearest-neighbour or context decisions"
        )
    if manifest.get("input_decisions") != [
        {"sample_id": ident, "mandatory": True, "reasons": ["input"]} for ident in expected
    ]:
        raise WorkflowError("Fixed input decisions must cover every input without ranking")
    if any(
        type(manifest.get(field)) is not int or manifest[field] != 0
        for field in (
            "context_budget",
            "effective_context_budget",
            "budget_overrun",
            "selected_contexts",
            "available_contexts",
        )
    ):
        raise WorkflowError("Fixed input selection must have zero context budgets/counts")
    if not isinstance(manifest.get("selection_id"), str) or not re.fullmatch(
        r"[A-Za-z0-9][A-Za-z0-9_.-]*", manifest["selection_id"]
    ):
        raise WorkflowError("Fixed input selection requires a safe selection_id")
    if type(manifest.get("seed")) is not int:
        raise WorkflowError("Fixed input selection seed must be an integer")


def write_fixed_input_selection(
    dataset_manifest, inclusion_evidence, *, output, selection_id="all-inputs", seed=42
):
    """Write a normal tree input manifest; inclusion is independent external evidence.

    This does not calculate profiles, distances, rankings or a phylogeny. Context
    rows may exist in the prepared bundle but are never selected by this policy.
    """
    from chronoclade.datasets import load_dataset
    from chronoclade.selection_manifest import _reference, load_selection_manifest

    dataset_path = Path(dataset_manifest).expanduser().resolve()
    inclusion_path = Path(inclusion_evidence).expanduser().resolve()
    dataset = load_dataset(dataset_path)
    directory = Path(output).expanduser().resolve()
    if directory.exists() and any(directory.iterdir()):
        raise WorkflowError("Fixed input output must be empty; existing evidence is immutable")
    directory.mkdir(parents=True, exist_ok=True)
    ids = [row["sample_id"] for row in dataset.samples if row["role"] == "input"]
    manifest = {
        "schema": "chronoclade.selection",
        "schema_version": 1,
        "status": "complete",
        "software": {"name": "chronoclade", "version": __version__},
        "selection_policy": POLICY,
        "selection_id": selection_id,
        "seed": seed,
        "source": {
            "dataset": _reference(directory, dataset_path),
            "dataset_id": dataset.dataset_id,
            "method": POLICY,
            "inclusion_evidence": _reference(directory, inclusion_path),
            "distance_evidence": {"status": "not_required", "reason": DISTANCE_REASON},
        },
        "query_ids": ids,
        "selected_sample_ids": ids,
        "mandatory_sample_ids": ids,
        "selected_context_ids": [],
        "available_context_ids": [],
        "pinned_context_ids": [],
        "required_nearest_context_ids": [],
        "nearest_neighbours": [],
        "decisions": [],
        "input_decisions": [
            {"sample_id": ident, "mandatory": True, "reasons": ["input"]} for ident in ids
        ],
        "context_budget": 0,
        "effective_context_budget": 0,
        "budget_overrun": 0,
        "selected_contexts": 0,
        "available_contexts": 0,
    }
    temporary = directory / "selection.pending.json"
    destination = directory / "selection.json"
    try:
        write_json(temporary, manifest)
        load_selection_manifest(temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination
