"""Frozen descriptive retrieval benchmark; no threshold fitting or pairwise CIs."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from datetime import date
import hashlib
import json
from pathlib import Path

import numpy as np

from chronoclade.cglin import resolved_prefix
from chronoclade.datasets import load_dataset

THRESHOLDS = (0, 1, 2, 4, 7, 10)
DEPTHS = tuple(range(5, 11))
REQUIRED = {"sample_id", "study_id", "event_id", "epi_status", "control_for_event", "source"}


def read_labels(path, sample_ids):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if not REQUIRED.issubset(reader.fieldnames or ()):
            raise ValueError(f"Labels require columns {sorted(REQUIRED)}")
        rows = list(reader)
    labels = {}
    for row in rows:
        row = {key: (value or "").strip() for key, value in row.items()}
        ident = row["sample_id"]
        if ident in labels or ident not in sample_ids:
            raise ValueError(f"Duplicate or unknown labelled sample: {ident}")
        if not row["study_id"] or not row["source"]:
            raise ValueError(f"Study and independent source required: {ident}")
        if row["epi_status"] not in {"outbreak", "control", "unknown"}:
            raise ValueError(f"Invalid epi_status: {ident}")
        if row["epi_status"] == "outbreak" and not row["event_id"]:
            raise ValueError(f"Outbreak requires documented event_id: {ident}")
        if row["epi_status"] != "control" and row["control_for_event"]:
            raise ValueError(f"Only controls may have control_for_event: {ident}")
        labels[ident] = row
    if set(labels) != set(sample_ids):
        raise ValueError("Labels must cover exactly the prepared sample IDs")
    return labels


def epi_label(left, right):
    if left["study_id"] != right["study_id"]:
        return "unknown"
    if (
        left["epi_status"] == right["epi_status"] == "outbreak"
        and left["event_id"] == right["event_id"]
    ):
        return "positive"
    for outbreak, control in ((left, right), (right, left)):
        if (
            outbreak["epi_status"] == "outbreak"
            and control["epi_status"] == "control"
            and control["control_for_event"] == outbreak["event_id"]
        ):
            return "negative"
    return "unknown"


def time_bin(left, right):
    # Exact ISO day only. Partial dates remain unknown, never imputed.
    try:
        if len(left) != 10 or len(right) != 10:
            return "unknown"
        days = abs((date.fromisoformat(left) - date.fromisoformat(right)).days)
    except (ValueError, TypeError):
        return "unknown"
    for limit, name in ((30, "0-30"), (90, "31-90")):
        if days <= limit:
            return name
    return ">90"


def lineage_records(dataset):
    result = {}
    for row in dataset.lineages:
        if row.get("kind") != "cglin":
            continue
        ident = row["sample_id"]
        if ident in result:
            raise ValueError(f"Multiple cgLIN assignments for {ident}; freeze one snapshot")
        evidence = dict(row.get("evidence") or {})
        evidence["cglin_raw"] = row["assignment"]
        evidence["cglin_status"] = row["resolution_status"]
        # Unknown scope is not evidence of a common assignment reference.
        scope = (
            row.get("scheme_id"),
            row.get("scheme_version"),
            evidence.get("cglin_database_sha256")
            or row.get("database_version")
            or evidence.get("cglin_reference_sha256")
            or evidence.get("cglin_frozen_export_sha256"),
        )
        prefixes = {depth: resolved_prefix(evidence, depth) for depth in DEPTHS}
        if (
            not all(scope)
            or any(str(x).lower() == "unknown" for x in scope)
            or row["resolution_status"] not in {"resolved", "complete", "partial", "provisional"}
        ):
            prefixes = dict.fromkeys(DEPTHS)
        result[ident] = (scope, prefixes)
    return result


def evaluate(dataset, labels, *, min_shared_loci=599):
    if type(min_shared_loci) is not int or not 1 <= min_shared_loci <= 629:
        raise ValueError("min_shared_loci must be an integer from 1 through 629")
    if len(dataset.profiles) != 1:
        raise ValueError("Require exactly one frozen allele matrix")
    matrix = dataset.profiles[0]
    catalogue = matrix.catalogue
    if not catalogue.complete or len(catalogue.loci) != 629:
        raise ValueError("Primary analysis requires a complete verified 629-locus catalogue")
    if not catalogue.database_sha256:
        raise ValueError("Primary analysis requires the frozen allele reference database SHA256")
    if set(matrix.sample_ids) != set(dataset.sample_ids):
        raise ValueError("Matrix must contain every prepared sample (missing calls may be zero)")
    lineages = lineage_records(dataset)
    metadata = {row["sample_id"]: row for row in dataset.samples}
    ids = matrix.sample_ids
    samples, pairs = [], []
    for i, ident in enumerate(ids):
        called = int(np.count_nonzero(matrix.codes[i]))
        row = {
            "sample_id": ident,
            "called_loci": called,
            "missing_loci": 629 - called,
            "genome_eligible": called >= min_shared_loci,
            **labels[ident],
        }
        for depth in DEPTHS:
            row[f"LIN{depth}_assigned"] = bool(lineages.get(ident, (None, {}))[1].get(depth))
        samples.append(row)
    for i, left in enumerate(ids):
        for j in range(i + 1, len(ids)):
            right = ids[j]
            a, b = matrix.codes[i], matrix.codes[j]
            shared = (a != 0) & (b != 0)
            n = int(np.count_nonzero(shared))
            mismatch = int(np.count_nonzero(a[shared] != b[shared])) if n else None
            la, lb = labels[left], labels[right]
            reason = ""
            if la["study_id"] != lb["study_id"]:
                reason = "different_study"
            elif la.get("patient_id") and la.get("patient_id") == lb.get("patient_id"):
                reason = "same_patient"
            elif n < min_shared_loci:
                reason = "insufficient_joint_calls"
            row = {
                "sample_a": left,
                "sample_b": right,
                "study_id": la["study_id"],
                "epi_label": epi_label(la, lb),
                "shared_loci": n,
                "allele_mismatches": mismatch,
                "eligible": not reason,
                "exclusion_reason": reason,
                "time_bin": time_bin(
                    metadata[left].get("collection_date", ""),
                    metadata[right].get("collection_date", ""),
                ),
            }
            for threshold in THRESHOLDS:
                row[f"alleles_le_{threshold}"] = mismatch <= threshold if not reason else None
            for depth in DEPTHS:
                sa, pa = lineages.get(left, (None, {}))
                sb, pb = lineages.get(right, (None, {}))
                row[f"LIN{depth}"] = (
                    pa[depth] == pb[depth] if sa == sb and pa.get(depth) and pb.get(depth) else None
                )
            pairs.append(row)
    candidates = [f"alleles_le_{t}" for t in THRESHOLDS] + [f"LIN{d}" for d in DEPTHS]
    queries = []
    for ident in ids:
        eligible = [p for p in pairs if p["eligible"] and ident in (p["sample_a"], p["sample_b"])]
        excluded = [
            p
            for p in pairs
            if not p["eligible"]
            and ident in (p["sample_a"], p["sample_b"])
            and p["epi_label"] == "positive"
        ]
        total_positive = sum(p["epi_label"] == "positive" for p in eligible)
        minimum = min((p["allele_mismatches"] for p in eligible), default=None)
        nearest = [p for p in eligible if p["allele_mismatches"] == minimum]
        for candidate in [*candidates, "nearest_all_ties"]:
            retrieved = (
                nearest
                if candidate == "nearest_all_ties"
                else [p for p in eligible if p[candidate] is True]
            )
            recovered = sum(p["epi_label"] == "positive" for p in retrieved)
            queries.append(
                {
                    "sample_id": ident,
                    "study_id": labels[ident]["study_id"],
                    "event_id": labels[ident]["event_id"],
                    "candidate": candidate,
                    "eligible_background": len(eligible),
                    "eligible_links": total_positive,
                    "excluded_links": len(excluded),
                    "recovered_links": recovered,
                    "recovery": recovered / total_positive if total_positive else None,
                    "burden": len(retrieved),
                    "minimum_allele_mismatches": minimum,
                    "nearest_ties_available": len(nearest),
                    "nearest_ties_retained": sum(p in retrieved for p in nearest),
                    "retrieved_ids": json.dumps(
                        [
                            p["sample_b"] if p["sample_a"] == ident else p["sample_a"]
                            for p in retrieved
                        ]
                    ),
                    "negative_retrieved_ids": json.dumps(
                        [
                            p["sample_b"] if p["sample_a"] == ident else p["sample_a"]
                            for p in retrieved
                            if p["epi_label"] == "negative"
                        ]
                    ),
                    "negative_retrieved": sum(p["epi_label"] == "negative" for p in retrieved),
                    "unknown_retrieved": sum(p["epi_label"] == "unknown" for p in retrieved),
                    "abstaining_pairs": (
                        0
                        if candidate == "nearest_all_ties"
                        else sum(p[candidate] is None for p in eligible)
                    ),
                }
            )
    grouped = defaultdict(list)
    for row in queries:
        grouped[(row["study_id"], row["event_id"], row["candidate"])].append(row)
    events = []
    for (study, event, candidate), rows in sorted(grouped.items()):
        recoveries = [r["recovery"] for r in rows if r["recovery"] is not None]
        events.append(
            {
                "study_id": study,
                "event_id": event,
                "candidate": candidate,
                "individuals": len(rows),
                "linked_queries": len(recoveries),
                "mean_query_recovery": float(np.mean(recoveries)) if recoveries else None,
                "median_burden": float(np.median([r["burden"] for r in rows])),
                "max_burden": max(r["burden"] for r in rows),
            }
        )
    candidate_counts = []
    for candidate in candidates:
        eligible = [p for p in pairs if p["eligible"]]
        retrieved = [p for p in eligible if p[candidate] is True]
        candidate_counts.append(
            {
                "candidate": candidate,
                "eligible_pairs": len(eligible),
                "abstaining_pairs": sum(p[candidate] is None for p in eligible),
                "retrieved_pairs": len(retrieved),
                "positive_retrieved": sum(p["epi_label"] == "positive" for p in retrieved),
                "negative_retrieved": sum(p["epi_label"] == "negative" for p in retrieved),
                "unknown_retrieved": sum(p["epi_label"] == "unknown" for p in retrieved),
            }
        )
    summary = {
        "schema": "chronoclade.research.klebsiella-grouping.v1",
        "dataset_id": dataset.dataset_id,
        "individuals": len(ids),
        "coverage_policy": {
            "min_shared_loci": min_shared_loci,
            "canonical_loci": 629,
            "minimum_shared_fraction": min_shared_loci / 629,
        },
        "studies": len({r["study_id"] for r in labels.values()}),
        "documented_events": len(
            {
                (r["study_id"], r["event_id"])
                for r in labels.values()
                if r["epi_status"] == "outbreak"
            }
        ),
        "unordered_pair_counts": dict(Counter(p["epi_label"] for p in pairs)),
        "eligible_pair_counts": dict(Counter(p["epi_label"] for p in pairs if p["eligible"])),
        "excluded_pair_counts": dict(
            Counter(p["exclusion_reason"] for p in pairs if not p["eligible"])
        ),
        "uncertainty": "Descriptive counts only; pairs are dependent. No independent-pair confidence intervals.",
        "catalogue": {
            "scheme_id": catalogue.scheme_id,
            "scheme_version": catalogue.scheme_version,
            "database_sha256": catalogue.database_sha256,
        },
    }
    return {
        "samples": samples,
        "pairs": pairs,
        "queries": queries,
        "events": events,
        "candidate_counts": candidate_counts,
        "summary": summary,
    }


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        if rows:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--min-shared-loci",
        default=599,
        type=int,
        help="Primary 599; protocol sensitivities 567, 598, 629",
    )
    args = parser.parse_args()
    dataset = load_dataset(args.dataset)
    labels = read_labels(args.labels, dataset.sample_ids)
    result = evaluate(dataset, labels, min_shared_loci=args.min_shared_loci)
    args.output.mkdir(parents=True, exist_ok=False)
    for name in ("samples", "pairs", "queries", "events", "candidate_counts"):
        write_csv(args.output / f"{name}.csv", result[name])
    result["summary"]["inputs_sha256"] = {
        name: hashlib.sha256(path.read_bytes()).hexdigest()
        for name, path in (("dataset_manifest", args.dataset), ("labels", args.labels))
    }
    (args.output / "summary.json").write_text(json.dumps(result["summary"], indent=2) + "\n")


if __name__ == "__main__":
    main()
