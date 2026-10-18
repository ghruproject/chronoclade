"""Saved biological-tree evidence for independent-selection comparisons."""

from collections import Counter
import csv
import math
from pathlib import Path

from Bio import Phylo

from chronoclade.errors import WorkflowError
from chronoclade.metadata_dates import sample_date_interval
from chronoclade.selection_manifest import _resolve
from chronoclade.stage_artifacts import read_dated_tree, reference, fingerprint


def anchor_date_evidence(path, result, anchors, *, output=None):
    """Match the actual dated-tree ancestor to its saved TreeTime node record."""
    level = (result.get("clock_confidence") or {}).get("interval_level")
    if level is not None and level != "90% max-posterior region":
        raise WorkflowError("Unsupported saved dating interval level")
    unavailable = {
        "status": "unavailable",
        "estimate": None,
        "lower_bound": None,
        "upper_bound": None,
        "interval_level": "90% max-posterior region",
    }
    if result["dated_tree"] is None:
        return {**unavailable, "reason": "No dated tree"}
    biological = read_dated_tree(_resolve(path.parent, result["dated_tree"], "dated_tree"))
    if not anchors or not set(anchors).issubset(t.name for t in biological.get_terminals()):
        raise WorkflowError("Shared anchors must occur in the dated biological tree")
    ancestor = biological.common_ancestor(anchors)
    descendants = sorted(t.name for t in ancestor.get_terminals())
    evidence = {
        **unavailable,
        "node": ancestor.name,
        "descendant_sample_ids": descendants,
        "descendant_count": len(descendants),
        "anchor_count": len(anchors),
    }
    spec = result.get("artifacts", {}).get("node_dates.csv")
    if spec is None:
        return {**evidence, "reason": "No saved node-date table"}
    source = _resolve(path.parent, spec, "node_dates")
    by_node = read_node_dates(source)
    node = by_node.get(str(ancestor.name))
    if node is None:
        raise WorkflowError("Node-date table does not identify the shared-anchor ancestor")
    return {
        **evidence,
        **node,
        "status": "complete",
        "interval_available": node["lower_bound"] is not None,
        "source": reference(output or path.parent, source),
    }


def input_target_date_evidence(path, result, targets, anchors, common_evidence, *, output=None):
    """Reuse saved dates for the declared input target, separately from shared contexts."""
    same = set(targets) == set(anchors)
    if not targets:
        evidence = {
            "status": "unavailable", "estimate": None,
            "lower_bound": None, "upper_bound": None,
            "interval_level": "90% max-posterior region",
            "reason": "No input targets declared",
        }
    elif same:
        evidence = common_evidence
    else:
        evidence = anchor_date_evidence(path, result, targets, output=output)
    return {
        **evidence,
        "target_sample_ids": list(targets),
        "target_count": len(targets),
        "extra_descendant_sample_ids": sorted(
            set(evidence.get("descendant_sample_ids", [])) - set(targets)
        ),
        "same_target_as_common_anchors": same,
    }


def read_node_dates(path):
    """Fail on malformed saved evidence rather than silently inventing intervals."""
    try:
        with Path(path).open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            required = {"node", "numeric_date", "lower_bound", "upper_bound"}
            if not required.issubset(reader.fieldnames or []):
                raise WorkflowError("Node-date CSV requires node, numeric_date and interval bounds")
            result = {}
            for row in reader:
                ident = row.get("node")
                if not ident or ident in result or None in row:
                    raise WorkflowError("Node-date CSV requires unique, complete node records")
                point = _finite(row.get("numeric_date"), "numeric_date")
                lower = _finite(row.get("lower_bound"), "lower_bound", optional=True)
                upper = _finite(row.get("upper_bound"), "upper_bound", optional=True)
                if (
                    (lower is None) != (upper is None)
                    or lower is not None
                    and not lower <= point <= upper
                ):
                    raise WorkflowError("Node-date CSV interval does not contain its estimate")
                result[ident] = {"estimate": point, "lower_bound": lower, "upper_bound": upper}
            return result
    except (OSError, UnicodeError, csv.Error) as error:
        raise WorkflowError(f"Cannot read node-date CSV: {error}") from error


def _finite(value, field, *, optional=False):
    if optional and value == "":
        return None
    try:
        number = float(value)
    except (ValueError, TypeError) as error:
        raise WorkflowError(f"Node-date CSV has invalid {field}") from error
    if not math.isfinite(number):
        raise WorkflowError(f"Node-date CSV has nonfinite {field}")
    return number


def date_summary(rows, eligible_ids):
    eligible = set(eligible_ids)
    dates = [sample_date_interval(row) for row in rows if row["sample_id"] in eligible]
    intervals = [date for date in dates if date.get("start") and date.get("end")]
    return {
        "dated_sample_count": len(eligible),
        "undated_sample_count": len(rows) - len(eligible),
        "date_interval": {
            "start": min((d["start"] for d in intervals), default=None),
            "end": max((d["end"] for d in intervals), default=None),
        },
        "date_precision_counts": dict(Counter(d.get("precision", "unknown") for d in dates)),
        "policy": "Original metadata interval bounds; clock eligibility comes from the saved tree assessment.",
    }


def tree_audit(path, tree, output):
    rows = _read_rows(_resolve(path.parent, tree["metadata"], "metadata"))
    labels = {
        row["sample_id"]: row.get("label") or row.get("accession") or row["sample_id"]
        for row in rows
    }
    mask_files = {
        name: reference(output, _resolve(path.parent, spec, name))
        for name, spec in tree["artifacts"].items()
        if name
        in {"clonalframeml.importation_status.txt", "clonalframeml.position_cross_reference.txt"}
    }
    return {
        **date_summary(rows, tree["dated_sample_ids"]),
        "temporal_status": tree["temporal_assessment"]["code"],
        "reference_sample_id": tree["reference_sample_id"],
        "reference_label": labels[tree["reference_sample_id"]],
        "reference": reference(output, _resolve(path.parent, tree["reference"], "reference")),
        "alignment": reference(output, _resolve(path.parent, tree["alignment"], "alignment")),
        "alignment_columns": tree["alignment_columns"],
        "complete_alignment_sites": tree["complete_alignment_sites"],
        "recombination_mask": {"summary": tree["recombination"], "artifacts": mask_files},
        "recombination_removed_sites": tree["recombination"].get("recombination_removed_sites"),
        "incomplete_only_removed_sites": tree["recombination"].get("incomplete_only_removed_sites"),
        "retained_clonal_percent": tree["recombination"].get("retained_clonal_percent"),
        "analysis_parameters": tree["parameters"],
        "root_policy": tree["root_policy"],
    }


def _read_rows(path):
    import json

    return json.loads(path.read_text())


def shared_anchor_splits(path, tree, anchors):
    """Project edge bipartitions onto fixed anchors, deduplicating suppressed edges."""
    biological = Phylo.read(
        _resolve(path.parent, tree["corrected_tree"], "corrected_tree"), "newick"
    )
    anchor_set = set(anchors)
    if not anchor_set.issubset(t.name for t in biological.get_terminals()):
        raise WorkflowError("Shared anchors are absent from corrected tree")
    splits = set()
    for node in biological.find_clades():
        left = anchor_set & {tip.name for tip in node.get_terminals()}
        right = anchor_set - left
        if min(len(left), len(right)) >= 2:
            sides = sorted((tuple(sorted(left)), tuple(sorted(right))))
            splits.add(tuple(sides))
    return splits


def split_occurrence(by_run, requested_runs, anchors):
    counts = Counter(split for splits in by_run.values() for split in splits)
    return {
        "status": "complete" if len(anchors) >= 4 and by_run else "unavailable",
        "requested_runs": requested_runs,
        "assessable_runs": len(by_run) if len(anchors) >= 4 else 0,
        "statistical_branch_support": {
            "status": "unavailable",
            "reason": "No calibrated branch-support analysis supplied",
        },
        "interpretation": "Unrooted split occurrence across selections is descriptive frequency, not bootstrap support, posterior clade probability or statistical branch support. Missing trees do not enter its denominator; unresolved trees count as absence of an observed split.",
        "splits": [
            {
                "split_id": fingerprint(split),
                "side_a": list(split[0]),
                "side_b": list(split[1]),
                "observed_runs": count,
                "assessable_runs": len(by_run),
                "frequency": count / len(by_run),
            }
            for split, count in sorted(counts.items())
        ],
    }
