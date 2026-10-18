"""Exact single linkage from dense distances with linear auxiliary memory.

A minimum spanning forest preserves the connected components at every inclusive
threshold. Its edges are evidence-preserving compression, not inferred pairwise
measurements. An induced sample subset requires a freshly computed forest.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np

from .thresholds import _distance, threshold_sweep


def _sample_indexes(sample_ids, distances, eligible, sample_indices):
    ids = list(sample_ids)
    if any(not isinstance(sample_id, str) or not sample_id for sample_id in ids):
        raise ValueError("sample IDs must be nonempty strings")
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate sample IDs")
    if not isinstance(distances, np.ndarray) or distances.shape != (len(ids), len(ids)):
        raise ValueError("distances must be a square numpy array matching sample IDs")
    if distances.dtype.kind not in "iuf":
        raise ValueError("distances must contain real numeric values")
    if eligible is not None and (
        not isinstance(eligible, np.ndarray)
        or eligible.shape != distances.shape or eligible.dtype.kind != "b"
    ):
        raise ValueError("eligible must be a boolean numpy array matching distances")
    indexes = list(range(len(ids))) if sample_indices is None else list(sample_indices)
    if any(
        isinstance(index, (bool, np.bool_)) or not isinstance(index, (int, np.integer))
        or not 0 <= index < len(ids) for index in indexes
    ):
        raise ValueError("sample_indices must contain in-range integer indexes")
    if len(set(indexes)) != len(indexes):
        raise ValueError("duplicate sample_indices")
    indexes.sort(key=lambda index: ids[index])
    return ids, np.asarray(indexes, dtype=np.intp)


def minimum_spanning_forest(
    sample_ids: Sequence[str],
    distances: np.ndarray,
    *,
    eligible: np.ndarray | None = None,
    sample_indices: Iterable[int] | None = None,
) -> list[tuple[str, str, int | float]]:
    """Return deterministic original edges preserving all threshold partitions.

    Matrix positions follow ``sample_ids``. NaN and positive infinity mean absent
    pairs. A false eligibility mask also removes a pair, including a stored raw
    zero. Without a mask, every finite off-diagonal zero is a genuine edge. This
    matters for raw count arrays that retain placeholder counts for unavailable
    pairs: supply their separately recorded eligibility, never infer it from zero.

    Validate exact symmetry and finite nonnegative available distances row by row;
    diagonal entries and values masked unavailable are ignored. Only the induced
    ``sample_indices`` graph is examined. A subset is recomputed from its original
    distances, never obtained by pruning a forest for the full graph.

    Vectorized Prim scans the square matrix in O(N²) time and O(N) auxiliary
    memory, plus <=N-1 returned edges. Arrays may be memory maps. It does not copy
    the square matrix or build Python objects for all pairwise edges. Integer
    weights retain their native precision; zero edges require no offset trick.
    Lexical sample order breaks all ties, independent of input matrix order.
    """
    ids, indexes = _sample_indexes(sample_ids, distances, eligible, sample_indices)
    names = [ids[index] for index in indexes]
    n = len(indexes)
    visited = np.zeros(n, dtype=bool)
    known = np.zeros(n, dtype=bool)
    best = np.zeros(n, dtype=distances.dtype)
    parent = np.full(n, -1, dtype=np.intp)
    edges = []
    for _ in range(n):
        candidates = np.flatnonzero(~visited & known)
        node = int(candidates[np.argmin(best[candidates])]) if len(candidates) else int(
            np.flatnonzero(~visited)[0]
        )
        visited[node] = True
        if parent[node] >= 0:
            first, second = sorted((names[node], names[parent[node]]))
            edges.append((first, second, best[node].item()))
        source = indexes[node]
        row = distances[source, indexes]
        column = distances[indexes, source]
        mask = np.ones(n, dtype=bool) if eligible is None else eligible[source, indexes].copy()
        reverse_mask = mask.copy() if eligible is None else eligible[indexes, source].copy()
        mask[node] = reverse_mask[node] = False
        if not np.array_equal(mask, reverse_mask):
            raise ValueError("eligibility must be symmetric on the selected samples")
        if np.any(mask & ((row < 0) | (column < 0))):
            raise ValueError("available distances must be nonnegative")
        available = mask & np.isfinite(row)
        reverse_available = reverse_mask & np.isfinite(column)
        if not np.array_equal(available, reverse_available) or np.any(
            row[available] != column[available]
        ):
            raise ValueError("available distances must be symmetric on the selected samples")
        update = available & ~visited & (~known | (row < best))
        tied = available & ~visited & known & (row == best) & (node < parent)
        best[update] = row[update]
        parent[update | tied] = node
        known[update] = True
    return sorted(edges, key=lambda edge: (edge[2], edge[0], edge[1]))


def matrix_threshold_sweep(
    sample_ids: Sequence[str],
    distances: np.ndarray,
    thresholds=None,
    *,
    eligible: np.ndarray | None = None,
    sample_indices: Iterable[int] | None = None,
    distance_unit: str | None = None,
) -> dict:
    """Reuse the shared threshold engine on an exact matrix spanning forest.

    Default snapshots are at forest-edge weights, sufficient to include every
    membership change. Use ``thresholds=[]`` for linear-size merge history without
    materializing any partitions; query it with ``summarise_threshold_history``.
    Requested snapshots otherwise contain O(N K) member IDs. Counts of forest
    edges are not counts of available original pairs or direct-distance evidence.
    """
    # Materialize a potentially one-shot subset iterator once for both consumers.
    indexes = None if sample_indices is None else list(sample_indices)
    edges = minimum_spanning_forest(sample_ids, distances, eligible=eligible, sample_indices=indexes)
    selected = list(sample_ids) if indexes is None else [sample_ids[index] for index in indexes]
    result = threshold_sweep(selected, edges, thresholds, distance_unit=distance_unit)
    result["representation"] = "minimum_spanning_forest"
    result["n_forest_edges"] = result.pop("n_edges")
    result["threshold_basis"] = "requested cuts or distinct forest-edge weights"
    return result


def summarise_threshold_history(
    sweep: dict, thresholds: Iterable, *, focal_sample_ids: Iterable[str] = ()
) -> list[dict]:
    """Query counts and focal neighbourhood references without copying partitions.

    ``sweep`` comes from the shared engine or ``matrix_threshold_sweep``. For each
    focal sample, ``node_id`` names its exact component in that sweep's merge
    history (leaf IDs index sorted sample IDs). ``focal_samples`` counts selected
    focal IDs in that component. No memberships are copied into the output.
    History processing is linear, followed by a path-compressed root lookup for
    each queried focal sample. Memory is O(N + M + K Q), for M merge nodes, K cuts
    and Q focal samples. Thresholds are sorted and deduplicated.
    """
    ids = sweep["sample_ids"]
    focal = list(focal_sample_ids)
    if len(set(focal)) != len(focal) or not set(focal) <= set(ids):
        raise ValueError("focal_sample_ids must be distinct known sample IDs")
    index = {sample_id: i for i, sample_id in enumerate(ids)}
    focal.sort()
    requested = sorted({_distance(value, "threshold") for value in thresholds})
    history = sweep["merge_history"]
    parents = list(range(len(ids) + len(history)))
    sizes = [1] * len(ids) + [merge["size"] for merge in history]
    focal_set = set(focal)
    focal_counts = [int(sample_id in focal_set) for sample_id in ids] + [0] * len(history)
    n_clusters, n_non_singleton, largest = len(ids), 0, int(bool(ids))
    cursor = 0
    results = []

    def root(node):
        while parents[node] != node:
            parents[node] = parents[parents[node]]
            node = parents[node]
        return node

    for threshold in requested:
        while cursor < len(history) and history[cursor]["distance"] <= threshold:
            merge = history[cursor]
            node, children = merge["node_id"], merge["children"]
            for child in children:
                parents[child] = node
            focal_counts[node] = sum(focal_counts[child] for child in children)
            n_clusters -= len(children) - 1
            n_non_singleton += 1 - sum(sizes[child] > 1 for child in children)
            largest = max(largest, merge["size"])
            cursor += 1
        neighbourhoods = {}
        for sample_id in focal:
            node = root(index[sample_id])
            neighbourhoods[sample_id] = {
                "node_id": node, "size": sizes[node], "focal_samples": focal_counts[node],
            }
        results.append({
            "threshold": threshold, "n_clusters": n_clusters,
            "n_non_singleton_clusters": n_non_singleton, "largest_cluster_size": largest,
            "focal_neighbourhoods": neighbourhoods,
        })
    return results
