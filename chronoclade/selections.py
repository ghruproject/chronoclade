"""Mandatory context policy and reproducible alternatives over supplied distances.

No distances or clock fits are computed here. Raw allele or SNP differences rank
neighbours; the supplied normalized distance is retained for diversity and audit.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from chronoclade.errors import WorkflowError
from chronoclade.metadata_dates import sample_date_interval


@dataclass(frozen=True)
class DistanceEvidence:
    """Stage-neutral adapter: frozen pair records and an explicit distance basis."""

    records: tuple[dict, ...]
    method: str = "cgmlst"
    distance_definition: str = "allele differences / shared called loci"

    def __post_init__(self):
        if not self.method or not self.distance_definition:
            raise WorkflowError("Distance evidence requires method and distance definition")


def read_distance_evidence(analysis: dict) -> DistanceEvidence:
    from chronoclade.matrix_distances import BinaryDistanceEvidence
    evidence = analysis.get("distance_evidence")
    if isinstance(evidence, (DistanceEvidence, BinaryDistanceEvidence)):
        return evidence
    if isinstance(evidence, dict) and evidence.get("schema") == "chronoclade.binary_distances":
        binary_path = analysis.get("paths", {}).get("distance_evidence")
        if not binary_path:
            raise WorkflowError("Binary distance descriptor requires its manifest path")
        matrix = BinaryDistanceEvidence(binary_path)
        if matrix.descriptor != evidence:
            raise WorkflowError("Binary distance descriptor disagrees with its manifest")
        return matrix
    if evidence is not None:
        raise WorkflowError("distance_evidence must be a validated distance adapter")
    binary = analysis.get("paths", {}).get("distance_evidence")
    if binary:
        return BinaryDistanceEvidence(binary)
    path = analysis.get("paths", {}).get("pairwise_distances")
    if path:
        if not Path(path).is_file():
            raise WorkflowError(f"Distance evidence does not exist: {path}")
        with Path(path).open(newline="", encoding="utf-8") as handle:
            return DistanceEvidence(
                tuple(dict(row) for row in csv.DictReader(handle)),
                method=analysis.get("distance_method", "cgmlst"),
                distance_definition=analysis.get(
                    "distance_definition", "allele differences / shared called loci"
                ),
            )
    return DistanceEvidence(
        (),
        method=analysis.get("distance_method", "cgmlst"),
        distance_definition=analysis.get(
            "distance_definition", "allele differences / shared called loci"
        ),
    )


def _number(value):
    try:
        result = float(value)
    except (ValueError, TypeError):
        return None
    return result if math.isfinite(result) and result >= 0 else None


def _raw_difference_field(evidence):
    return {"cgmlst": "allele_differences", "alleles": "allele_differences",
            "direct-dna": "snp_differences"}.get(evidence.method)


def _distance_maps(evidence):
    from chronoclade.matrix_distances import BinaryDistanceEvidence
    if isinstance(evidence, BinaryDistanceEvidence):
        return evidence.value_map(), evidence.value_map(raw=True)
    normalized, nearest = {}, {}
    raw_field = _raw_difference_field(evidence)
    raw_available = raw_field is not None and any(
        _number(row.get(raw_field)) is not None for row in evidence.records
    )
    for row in evidence.records:
        a, b = row.get("sample_id_1"), row.get("sample_id_2")
        if not isinstance(a, str) or not a or not isinstance(b, str) or not b:
            raise WorkflowError("Distance evidence requires nonempty sample_id_1/sample_id_2")
        distance = _number(row.get("distance"))
        raw = _number(row.get(raw_field)) if raw_field else None
        rank = raw if raw_available else distance
        for values, value in ((normalized, distance), (nearest, rank)):
            if value is None:
                continue
            if (a, b) in values and values[a, b] != value:
                raise WorkflowError(f"Conflicting distance evidence for {a!r}, {b!r}")
            values[a, b] = values[b, a] = value
    return normalized, nearest


def selected_context_records(
    selection: dict, queries: list[dict], context: list[dict]
) -> list[dict]:
    """Validate a saved selection against its exact pool before acquiring assemblies."""
    by_id = {row["sample_id"]: row for row in context}
    context_ids = selection.get("selected_context_ids", [])
    sample_ids = selection.get("selected_sample_ids", [])
    query_ids = [row["sample_id"] for row in queries]
    if (
        len(by_id) != len(context)
        or selection.get("available_context_ids") != sorted(by_id)
        or selection.get("query_ids") != query_ids
        or len(context_ids) != len(set(context_ids))
        or not set(context_ids).issubset(by_id)
        or sample_ids != query_ids + context_ids
        or len(sample_ids) != len(set(sample_ids))
    ):
        raise WorkflowError(
            "Shared context selection does not match the resolved query lineage pool"
        )
    decisions = selection.get("decisions", [])
    if [decision.get("sample_id") for decision in decisions] != context_ids or any(
        not isinstance(decision.get("reason"), str) or not decision["reason"]
        for decision in decisions
    ):
        raise WorkflowError("Shared context selection decisions do not match selected identifiers")
    return [dict(by_id[d["sample_id"]], selection_reason=d["reason"]) for d in decisions]


def _pool(queries, context, include):
    query_ids = [row["sample_id"] for row in queries]
    by_id = {row["sample_id"]: row for row in context}
    if len(by_id) != len(context) or len(query_ids) != len(set(query_ids)):
        raise WorkflowError("Duplicate sample identifiers in selection pool")
    if set(query_ids).intersection(by_id):
        raise WorkflowError("Input and contextual sample identifiers must be disjoint")
    if any(
        not isinstance(identifier, str) or not identifier for identifier in query_ids + list(by_id)
    ):
        raise WorkflowError("Selection requires nonempty string sample identifiers")
    pinned = []
    for identifier in include or []:
        matches = [
            r["sample_id"]
            for r in context
            if identifier in {r["sample_id"], r.get("source_genome_id"), r.get("accession")}
        ]
        if len(matches) != 1:
            raise WorkflowError(
                f"Requested context genome {identifier!r} is not uniquely in the pool"
            )
        if matches[0] not in pinned:
            pinned.append(matches[0])
    return query_ids, by_id, pinned


def _mandatory(query_ids, by_id, pinned, nearest, count, method):
    reasons = {identifier: ["user_requested"] for identifier in pinned}
    evidence, queues = [], {}
    for query in sorted(query_ids):
        candidates = sorted(
            (identifier for identifier in by_id if (query, identifier) in nearest),
            key=lambda identifier: (nearest[query, identifier], identifier),
        )
        cutoff = nearest[query, candidates[min(count, len(candidates)) - 1]] if candidates else None
        retained = [identifier for identifier in candidates if nearest[query, identifier] <= cutoff]
        queues[query] = retained
        evidence.append(
            {
                "query_id": query,
                "requested_neighbours": count,
                "boundary_distance": cutoff,
                "required_ids": retained,
                "boundary_tie_ids": [i for i in retained if nearest[query, i] == cutoff],
                "tie_expanded": len(retained) > count,
            }
        )
    # Round-robin ordering keeps the first required neighbour of every input visible.
    selected = list(pinned)
    for rank in range(max((len(q) for q in queues.values()), default=0)):
        for query, candidates in queues.items():
            if rank >= len(candidates):
                continue
            identifier = candidates[rank]
            if identifier not in reasons:
                selected.append(identifier)
                reasons[identifier] = []
            reasons[identifier].append(f"{method}_neighbour:{query}:rank={rank + 1}")
    return selected, reasons, evidence


def _strata(by_id, analysis):
    groups = {
        identifier: str(group["group_id"])
        for group in analysis.get("genetic_groups", [])
        for identifier in group.get("sample_ids", [])
    }
    strata = defaultdict(list)
    for identifier, row in sorted(by_id.items()):
        interval = sample_date_interval(row)
        years = (interval["start"][:4], interval["end"][:4]) if interval["status"] == "valid" else None
        date_bin = (years[0] if years[0] == years[1] else "..".join(years)) if years else "Unknown"
        cell = (
            groups.get(identifier, "unassigned"),
            date_bin,
            str(row.get("nuts2") or row.get("region") or row.get("country") or "Unknown"),
            str(row.get("host") or "Unknown"),
        )
        strata[cell].append(identifier)
    return dict(strata)


def _quotas(strata, mandatory, size, seed):
    """Freeze metadata quotas once; only representative identity varies per run."""
    cells = [
        cell for cell, identifiers in sorted(strata.items()) if set(identifiers) - set(mandatory)
    ]
    random.Random(seed).shuffle(cells)
    slots = min(max(0, size - len(mandatory)), math.ceil(size / 4), len(cells))
    return {cell: 1 for cell in cells[:slots]}


def _select(by_id, mandatory, required_reasons, distances, strata, quotas, size, seed, method):
    selected = list(mandatory)
    reasons = {identifier: list(values) for identifier, values in required_reasons.items()}
    rng = random.Random(seed)
    for cell in quotas:
        candidates = [i for i in strata[cell] if i not in reasons]
        identifier = rng.choice(candidates)
        selected.append(identifier)
        reasons[identifier] = ["genetic_group_time_region_representative"]
    while len(selected) < min(size, len(by_id)):
        remaining = sorted(set(by_id) - set(selected))
        scores = {
            identifier: min(
                distances[identifier, other]
                for other in selected
                if (identifier, other) in distances
            )
            for identifier in remaining
            if any((identifier, other) in distances for other in selected)
        }
        if scores:
            best = max(scores.values())
            identifier = rng.choice([i for i in remaining if scores.get(i) == best])
            reason = f"{method}_diversity_representative"
        else:
            identifier = rng.choice(remaining)
            reason = "metadata_background_no_comparable_profile"
        selected.append(identifier)
        reasons[identifier] = [reason]
    return selected, reasons


def _coverage(strata, selected):
    counts = Counter(selected)
    return [
        {"stratum": list(cell), "available": len(ids), "selected": sum(counts[i] for i in ids)}
        for cell, ids in sorted(strata.items())
    ]


def representative_coverage(queries, context, selected_ids, evidence, *, ranking_metric=None):
    """Measure every available pool member against its nearest selected genome.

    Binary evidence is queried one scalar pair at a time, keeping temporary memory
    proportional to the pool size rather than copying its square matrices. Dense
    evidence is scanned once and only per-genome best rows and ties are retained.
    cgMLST uses raw allele differences and direct DNA uses raw SNP differences
    for ranking. Other methods use their saved distance. Selected samples have a
    zero-distance identity entry in addition to
    any comparable external representatives.
    """
    from chronoclade.matrix_distances import BinaryDistanceEvidence

    pool = sorted([row["sample_id"] for row in queries + context])
    pool_set = set(pool)
    representatives = list(selected_ids)
    representative_set = set(representatives)
    if len(representative_set) != len(representatives) or not representative_set.issubset(pool_set):
        raise WorkflowError("Coverage representatives must be unique members of the full pool")

    raw_field = _raw_difference_field(evidence)
    raw_ranking = raw_field is not None and (
        isinstance(evidence, BinaryDistanceEvidence)
        or ranking_metric == f"raw_{raw_field}_when_available"
        or ranking_metric is None
        and any(_number(row.get(raw_field)) is not None for row in evidence.records)
    )
    metric = raw_field if raw_ranking else "distance"
    best = {ident: None for ident in pool}
    comparable_counts = Counter()

    def consider(sample_id, representative_id, record, *, identity=False):
        def numeric(field, *, integer=False):
            value = _number(record.get(field))
            if value is not None and integer and float(value).is_integer():
                return int(value)
            return value

        shared = numeric("shared_called_loci", integer=True)
        call_overlap = numeric("call_overlap")
        scheme = numeric("scheme_loci", integer=True)
        if scheme is None and shared is not None and call_overlap:
            scheme = round(shared / call_overlap)
        if call_overlap is None and shared is not None and scheme:
            call_overlap = shared / scheme
        candidate = {
            "sample_id": representative_id,
            "comparison": "same_sample_identity" if identity else "pairwise_distance",
            "distance": numeric("distance"),
            "allele_differences": numeric("allele_differences", integer=True),
            "shared_called_loci": shared,
            "scheme_loci": scheme,
            "call_overlap": call_overlap,
            "panel_loci": numeric("panel_loci", integer=True),
        }
        if evidence.method == "direct-dna":
            candidate["snp_differences"] = numeric("snp_differences", integer=True)
            candidate["shared_called_sites"] = numeric("shared_called_sites", integer=True)
        value = 0 if identity and raw_ranking else (0.0 if identity else _number(record.get(metric)))
        if raw_ranking and value is not None and float(value).is_integer():
            value = int(value)
        if value is None:
            return
        comparable_counts[sample_id] += 0 if identity else 1
        candidate["ranking_value"] = value
        current = best[sample_id]
        if current is None or value < current["value"] - (0 if raw_ranking else 1e-12):
            best[sample_id] = {"value": value, "ties": [candidate]}
        elif (
            value == current["value"]
            if raw_ranking
            else math.isclose(value, current["value"], rel_tol=0, abs_tol=1e-12)
        ):
            current["ties"].append(candidate)

    if isinstance(evidence, BinaryDistanceEvidence):
        for sample_id in pool:
            for representative_id in representatives:
                if sample_id == representative_id:
                    continue
                record = evidence.get(sample_id, representative_id)
                if record is not None:
                    consider(sample_id, representative_id, record)
    else:
        # One pass over the existing pair rows avoids building another O(n²) map.
        for row in evidence.records:
            left, right = row.get("sample_id_1"), row.get("sample_id_2")
            if left not in pool_set or right not in pool_set or left == right:
                continue
            if right in representative_set:
                consider(left, right, row)
            if left in representative_set:
                consider(right, left, row)

    # An included genome covers itself exactly, while the record keeps that
    # identity separate from callable pair evidence and its denominator.
    for ident in representatives:
        consider(ident, ident, {"distance": 0.0, "allele_differences": 0, "snp_differences": 0}, identity=True)

    rows = []
    for ident in pool:
        result = best[ident]
        ties = sorted(result["ties"], key=lambda item: item["sample_id"]) if result else []
        identity_only = bool(ties) and all(
            tie["comparison"] == "same_sample_identity" for tie in ties
        )
        rows.append({
            "sample_id": ident,
            "selected": ident in representative_set,
            "candidate_representatives": max(0, len(representatives) - int(ident in representative_set)),
            "comparable_representatives": comparable_counts[ident],
            "coverage_status": (
                "selected_identity_only" if identity_only else
                "covered" if ties else "no_comparable_selected_representative"
            ),
            "ranking_metric": metric,
            "nearest_value": result["value"] if result else None,
            "nearest_tie_count": len(ties),
            "nearest_tie_evidence": ties,
            "nearest_selected_ids": [tie["sample_id"] for tie in ties],
            "nearest_distance": min(
                (value for tie in ties if (value := _number(tie["distance"])) is not None),
                default=None,
            ),
        })

    values = sorted(row["nearest_value"] for row in rows if row["nearest_value"] is not None)

    def quantile(fraction):
        if not values:
            return None
        position = (len(values) - 1) * fraction
        lower, upper = math.floor(position), math.ceil(position)
        if lower == upper:
            return values[lower]
        return values[lower] + (values[upper] - values[lower]) * (position - lower)

    summary = {
        "schema": "chronoclade.selection-representative-coverage",
        "schema_version": 1,
        "status": "complete",
        "scope": "full_available_selection_pool",
        "sample_count": len(pool),
        "scope_sample_ids_sha256": hashlib.sha256(
            json.dumps(pool, separators=(",", ":")).encode()
        ).hexdigest(),
        "selected_sample_ids": representatives,
        "selected_count": len(representatives),
        "ranking_metric": metric,
        "distance_definition": evidence.distance_definition,
        "coverage_status_counts": dict(sorted(Counter(row["coverage_status"] for row in rows).items())),
        "comparable_sample_count": len(values),
        "noncomparable_sample_count": len(pool) - len(values),
        "nearest_value_quantiles": {
            "minimum": values[0] if values else None,
            "q25": quantile(0.25),
            "median": quantile(0.5),
            "q75": quantile(0.75),
            "maximum": values[-1] if values else None,
        },
    }
    return summary, rows


def select_context_ensemble(
    queries: list[dict],
    context: list[dict],
    analysis: dict,
    *,
    replicates: int = 1,
    size: int = 50,
    nearest_per_query: int = 3,
    representative_reserve: int = 0,
    include: list[str] | None = None,
    seed: int = 42,
) -> tuple[list[dict], dict]:
    """Keep mandatory IDs fixed, varying eligible representatives with recorded seeds.

    ``size`` remains the requested context budget. ``representative_reserve`` is
    an explicit minimum number of optional slots above the mandatory union;
    it can extend that budget, but never removes required neighbours or ties.
    The available pool may contain fewer optional genomes than requested.
    """
    if (
        type(size) is not int
        or size < 0
        or type(nearest_per_query) is not int
        or nearest_per_query < 1
    ):
        raise WorkflowError("Context size must be non-negative and nearest-per-query positive")
    if type(replicates) is not int or replicates < 1 or type(seed) is not int:
        raise WorkflowError("Replicates must be positive and seed must be an integer")
    if type(representative_reserve) is not int or representative_reserve < 0:
        raise WorkflowError("Representative reserve must be a non-negative integer")
    query_ids, by_id, pinned = _pool(queries, context, include)
    all_records = {row["sample_id"]: row for row in queries} | by_id
    supplied = read_distance_evidence(analysis)
    distances, nearest = _distance_maps(supplied)
    mandatory, reasons, neighbour_evidence = _mandatory(
        query_ids,
        by_id,
        pinned,
        nearest,
        nearest_per_query,
        supplied.method,
    )
    # Explicit optional room preserves every required nearest boundary tie.
    effective_size = max(size, len(mandatory) + representative_reserve)
    strata = _strata(by_id, analysis)
    quotas = _quotas(strata, mandatory, effective_size, seed)
    selections, seen = [], {}
    no_choices = len(mandatory) >= effective_size or len(by_id) <= effective_size
    for index in range(replicates):
        attempted_seeds = []
        for attempt in range(1 if no_choices or index == 0 else 32):
            run_seed = seed + index + attempt * replicates
            attempted_seeds.append(run_seed)
            selected, decisions = _select(
                by_id,
                mandatory,
                reasons,
                distances,
                strata,
                quotas,
                effective_size,
                run_seed,
                supplied.method,
            )
            key = tuple(sorted(selected))
            if key not in seen:
                break
        duplicate_of = seen.get(key)
        selection_id = f"selection-{index + 1:03d}"
        seen.setdefault(key, selection_id)
        audit = {
            "selection_id": selection_id,
            "query_ids": query_ids,
            "available_context_ids": sorted(by_id),
            "selected_context_ids": selected,
            "selected_sample_ids": query_ids + selected,
            "mandatory_sample_ids": query_ids + mandatory,
            "mandatory_context_ids": mandatory,
            "pinned_context_ids": pinned,
            "required_nearest_context_ids": [
                i
                for i in mandatory
                if any(reason.startswith(f"{supplied.method}_neighbour:") for reason in reasons[i])
            ],
            "context_budget": size,
            "representative_reserve": representative_reserve,
            "effective_context_budget": effective_size,
            "available_contexts": len(by_id),
            "selected_contexts": len(selected),
            "budget_overrun": max(0, len(mandatory) - size),
            "nearest_per_query": nearest_per_query,
            "seed": run_seed,
            "attempted_seeds": attempted_seeds,
            "source_method": supplied.method,
            "distance_definition": supplied.distance_definition,
            "nearest_ranking": f"raw_{_raw_difference_field(supplied)}_when_available"
            if _raw_difference_field(supplied)
            else "supplied_distance",
            "nearest_neighbours": neighbour_evidence,
            "decisions": [
                {
                    "sample_id": i,
                    "reason": decisions[i][0],
                    "reasons": decisions[i],
                    "mandatory": i in mandatory,
                }
                for i in selected
            ],
            "input_decisions": [
                {"sample_id": i, "reasons": ["input"], "mandatory": True} for i in query_ids
            ],
            "quotas": [
                {"stratum": list(cell), "optional_slots": value} for cell, value in quotas.items()
            ],
            "strata_coverage": _coverage(strata, selected),
            "strata_coverage_scope": "context_pool; every input is retained separately",
            "date_strata_policy": "validated canonical year or year range; missing/invalid/future dates are Unknown",
            "missing_metadata": {
                field: [i for i in sorted(query_ids + selected)
                        if (sample_date_interval(all_records[i])["status"] != "valid"
                            if field == "collection_date" else not all_records[i].get(field))]
                for field in ("collection_date", "country", "host")
            },
            "queries_without_selected_comparable_neighbour": [
                q for q in sorted(query_ids) if not any((q, i) in nearest for i in selected)
            ],
            "alternative_status": "duplicate" if duplicate_of else "distinct",
            "duplicate_of": duplicate_of,
            "duplicate_reason": (
                "no_nonmandatory_choices"
                if no_choices
                else "no_distinct_set_after_recorded_attempts"
            )
            if duplicate_of
            else None,
            "variation_status": "no_nonmandatory_choices"
            if no_choices
            else "no_distinct_alternative_found"
            if duplicate_of
            else "variable_representatives",
            "scope": "Selection within the supplied distance catalogue; no clock-fit optimisation.",
        }
        selections.append(audit)
    overlaps = []
    for index, left in enumerate(selections):
        a = set(left["selected_sample_ids"])
        for right in selections[index + 1 :]:
            b = set(right["selected_sample_ids"])
            overlaps.append(
                {
                    "left": left["selection_id"],
                    "right": right["selection_id"],
                    "shared_ids": sorted(a & b),
                    "intersection": len(a & b),
                    "union": len(a | b),
                    "jaccard": len(a & b) / len(a | b) if a | b else 1.0,
                }
            )
    return selections, {
        "seed": seed,
        "requested_replicates": replicates,
        "common_anchor_ids": query_ids + mandatory,
        "target_ids": query_ids,
        "pairwise_overlap": overlaps,
        "cumulative_unique_context_ids": sorted(
            {i for s in selections for i in s["selected_context_ids"]}
        ),
        "distinct_selections": len(seen),
        "quotas": selections[0]["quotas"],
    }


def select_assembly_context(
    queries: list[dict],
    context: list[dict],
    analysis: dict,
    *,
    size: int = 50,
    nearest_per_query: int = 3,
    representative_reserve: int = 0,
    include: list[str] | None = None,
    seed: int = 42,
) -> tuple[list[dict], dict]:
    """Legacy single-set API using the mandatory/ensemble policy."""
    selections, _ = select_context_ensemble(
        queries,
        context,
        analysis,
        size=size,
        nearest_per_query=nearest_per_query,
        representative_reserve=representative_reserve,
        include=include,
        seed=seed,
    )
    audit = selections[0]
    by_id = {row["sample_id"]: row for row in context}
    records = [
        dict(by_id[d["sample_id"]], selection_reason=d["reason"]) for d in audit["decisions"]
    ]
    # Older consumers store pool identity separately. Preserve input-order independence
    # of the single-selection audit while the manifest API retains input order exactly.
    legacy_audit = {
        k: v
        for k, v in audit.items()
        if k
        not in {
            "query_ids",
            "selected_sample_ids",
            "mandatory_sample_ids",
            "input_decisions",
        }
    }
    return records, legacy_audit
