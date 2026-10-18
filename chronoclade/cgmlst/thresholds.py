"""Descriptive single-linkage sweeps over supplied raw pairwise distances.

Missing pairs are absent edges, never zero distances. This module does not choose a
biological cutoff, infer distance units, or assess temporal stability.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from itertools import groupby
from math import isfinite
from numbers import Integral, Real


def _distance(value: Real, label: str) -> int | float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{label} must be a finite nonnegative number")
    if value < 0 or (not isinstance(value, Integral) and not isfinite(value)):
        raise ValueError(f"{label} must be a finite nonnegative number")
    return int(value) if isinstance(value, Integral) else float(value)


class _Components:
    def __init__(self, count: int):
        self.parent = list(range(count))
        self.size = [1] * count
        self.node = list(range(count))

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, first: int, second: int) -> None:
        first, second = self.find(first), self.find(second)
        if first == second:
            return
        if (self.size[first], -first) < (self.size[second], -second):
            first, second = second, first
        self.parent[second] = first
        self.size[first] += self.size[second]

    def clusters(self, sample_ids: list[str]) -> list[list[str]]:
        groups = defaultdict(list)
        for index, sample_id in enumerate(sample_ids):
            groups[self.find(index)].append(sample_id)
        return sorted(groups.values())


def threshold_sweep(
    sample_ids: Iterable[str],
    edges: Iterable[tuple[str, str, Real]],
    thresholds: Iterable[Real] | None = None,
    *,
    distance_unit: str | None = None,
) -> dict:
    """Return JSON-compatible partitions for edges with distance <= threshold.

    IDs must be unique nonempty strings. Undirected duplicate edges with identical
    distances are accepted once; conflicting duplicates and self-edges are errors.
    Thresholds are deduplicated and sorted. By default, use distinct observed edge
    weights (so an edgeless graph has no requested partitions). ``initial_clusters``
    always records every sample, including isolated samples and edgeless inputs.

    Sort edges once and incrementally union connected components. Equal-weight
    edges act simultaneously: the compact merge history uses multifurcating nodes,
    avoiding an arbitrary binary order within ties. Leaves 0..n-1 index the sorted
    returned ``sample_ids``; subsequent node IDs refer to ``merge_history`` entries.
    History covers all observed edges, even beyond the requested thresholds.

    ``stable_ranges`` records runs of unchanged membership across consecutive
    requested thresholds only, including one-threshold runs. Endpoints are observed
    thresholds, not inferred continuous intervals or temporal/biological stability.
    Units are a caller-supplied label and never change or normalize the distances.
    Cost is O(E log E + E alpha(N) + N K), plus sorting the returned memberships,
    for E unique edges, N samples and K requested thresholds.
    """
    ids = list(sample_ids)
    if any(not isinstance(sample_id, str) or not sample_id for sample_id in ids):
        raise ValueError("sample IDs must be nonempty strings")
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate sample IDs")
    ids.sort()
    if distance_unit is not None and (
        not isinstance(distance_unit, str) or not distance_unit.strip()
    ):
        raise ValueError("distance_unit must be a nonempty string or None")
    index = {sample_id: position for position, sample_id in enumerate(ids)}
    unique_edges = {}
    for edge in edges:
        try:
            first, second, value = edge
        except (TypeError, ValueError) as error:
            raise ValueError("each edge must contain two sample IDs and a distance") from error
        if (
            not isinstance(first, str) or not isinstance(second, str)
            or first not in index or second not in index
        ):
            raise ValueError(f"edge contains unknown sample ID: {first!r}, {second!r}")
        if first == second:
            raise ValueError(f"self-edge is not allowed: {first!r}")
        weight = _distance(value, "edge distance")
        pair = tuple(sorted((index[first], index[second])))
        if pair in unique_edges and unique_edges[pair] != weight:
            raise ValueError(f"conflicting duplicate edge: {first!r}, {second!r}")
        unique_edges[pair] = weight
    ordered_edges = sorted((weight, *pair) for pair, weight in unique_edges.items())
    requested = sorted(
        set(unique_edges.values()) if thresholds is None else {
            _distance(value, "threshold") for value in thresholds
        }
    )
    components = _Components(len(ids))
    partitions = []
    stable_ranges = []
    history = []
    next_threshold = 0

    def snapshot(threshold: int | float) -> None:
        clusters = components.clusters(ids)
        partition_index = len(partitions)
        if partitions and clusters == partitions[-1]["clusters"]:
            stable_ranges[-1]["end_threshold"] = threshold
            stable_ranges[-1]["threshold_count"] += 1
        else:
            stable_ranges.append({
                "start_threshold": threshold,
                "end_threshold": threshold,
                "threshold_count": 1,
                "partition_index": partition_index,
            })
        partitions.append({
            "threshold": threshold,
            "clusters": clusters,
            "n_clusters": len(clusters),
            "n_non_singleton_clusters": sum(len(group) > 1 for group in clusters),
            "largest_cluster_size": max(map(len, clusters), default=0),
        })

    for weight, batch in groupby(ordered_edges, key=lambda edge: edge[0]):
        while next_threshold < len(requested) and requested[next_threshold] < weight:
            snapshot(requested[next_threshold])
            next_threshold += 1
        tied_edges = list(batch)
        prior_roots = {
            components.find(item) for _, first, second in tied_edges for item in (first, second)
        }
        for _, first, second in tied_edges:
            components.union(first, second)
        merged = defaultdict(list)
        for root in prior_roots:
            merged[components.find(root)].append(components.node[root])
        for root, children in sorted(merged.items(), key=lambda item: min(item[1])):
            if len(children) > 1:
                node_id = len(ids) + len(history)
                history.append({
                    "distance": weight,
                    "node_id": node_id,
                    "children": sorted(children),
                    "size": components.size[root],
                })
                components.node[root] = node_id
    while next_threshold < len(requested):
        snapshot(requested[next_threshold])
        next_threshold += 1
    return {
        "method": "single_linkage",
        "distance_unit": distance_unit,
        "sample_ids": ids,
        "n_samples": len(ids),
        "n_edges": len(unique_edges),
        "initial_clusters": [[sample_id] for sample_id in ids],
        "thresholds": requested,
        "partitions": partitions,
        "merge_history": history,
        "stable_ranges": stable_ranges,
        "stable_range_scope": "unchanged membership across consecutive requested thresholds only",
    }
