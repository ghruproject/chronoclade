"""Reuse selection coverage evidence, including older immutable selections."""

import csv
import json

from chronoclade.datasets import load_dataset
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import _resolve, _read
from chronoclade.selections import DistanceEvidence, representative_coverage
from chronoclade.stage_artifacts import reference


def selection_coverage(ensemble_path, ensemble, output):
    """Every requested selection contributes coverage, even if its tree failed."""
    dataset = load_dataset(
        _resolve(
            ensemble_path.parent, ensemble["source"]["dataset"], "source.dataset", external=True
        )
    )
    by_id = {row["sample_id"]: row for row in dataset.samples}
    selections = {}
    coverage = []
    evidence = None
    union = set()
    cumulative = []
    expected_pool = None
    for index, entry in enumerate(ensemble["selections"], 1):
        path = _resolve(ensemble_path.parent, entry, "selection")
        selection = _read(path)  # Already checked by load_selection_ensemble.
        pool_ids = set(selection["query_ids"] + selection["available_context_ids"])
        if expected_pool is not None and pool_ids != expected_pool:
            raise WorkflowError("Comparison selections must describe the same available pool")
        expected_pool = pool_ids
        ident = selection["selection_id"]
        selections[ident] = selection
        union.update(selection["selected_sample_ids"])
        cumulative.append({"selection_id": ident, "cumulative_unique_genomes": len(union)})
        saved = selection.get("representative_coverage")
        if saved is not None:
            artifact = _resolve(path.parent, saved["artifact"], "representative_coverage.artifact")
            summary = {
                **saved,
                "artifact": reference(output, artifact),
                "origin": "saved_selection_artifact",
            }
        else:
            # Use the same coverage engine as selection publication; do not resample.
            if evidence is None:
                evidence_path = _resolve(
                    path.parent,
                    selection["source"]["distance_evidence"],
                    "source.distance_evidence",
                    external=True,
                )
                data = _read(evidence_path)
                if data.get("schema") == "chronoclade.binary_distances":
                    from chronoclade.matrix_distances import BinaryDistanceEvidence

                    evidence = BinaryDistanceEvidence(evidence_path)
                else:
                    evidence = DistanceEvidence(
                        tuple(data["records"]), data["method"], data["distance_definition"]
                    )
            summary, rows = representative_coverage(
                [by_id[i] for i in selection["query_ids"]],
                [by_id[i] for i in selection["available_context_ids"]],
                selection["selected_sample_ids"],
                evidence,
                ranking_metric=selection.get("nearest_ranking"),
            )
            artifact = output / f"selection-{index:03d}-coverage.csv"
            with artifact.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                for row in rows:
                    writer.writerow(
                        {
                            key: json.dumps(value, separators=(",", ":"))
                            if isinstance(value, (list, dict))
                            else value
                            for key, value in row.items()
                        }
                    )
            summary.update(
                artifact=reference(output, artifact), origin="computed_from_saved_distance_evidence"
            )
        coverage.append(
            {
                "selection_id": ident,
                "representative_coverage": summary,
                "strata_coverage": selection.get("strata_coverage", []),
                "strata_scope": selection.get("strata_coverage_scope", "context_pool"),
                "selected_sample_count": len(selection["selected_sample_ids"]),
                "selected_input_count": len(selection["query_ids"]),
                "selected_context_count": len(selection["selected_context_ids"]),
                "selection_manifest": reference(output, path),
            }
        )
    first = next(iter(selections.values()))
    pool = sorted(first["query_ids"] + first["available_context_ids"])
    return selections, {
        "scope": "full_available_selection_pool",
        "partition_block_id": ensemble["source"].get("partition_block_id"),
        "pool_sample_count": len(pool),
        "input_count": len(first["query_ids"]),
        "context_pool_count": len(first["available_context_ids"]),
        "requested_selection_union_ids": sorted(union),
        "requested_selection_union_count": len(union),
        "requested_selection_union_fraction": len(union) / len(pool),
        "cumulative_unique_genomes": cumulative,
        "pairwise_overlap": ensemble["pairwise_overlap"],
        "runs": coverage,
        "interpretation": "Coverage and overlap describe every requested saved selection, including failed/missing analyses. Genetic/time/location/host strata counts cover the context pool; every input is retained separately. This available block/pool is not every genome in the prepared dataset.",
    }
