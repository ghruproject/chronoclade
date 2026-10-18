"""Dense matrix compression preserves original single-linkage connectivity."""

import json

import numpy as np
import pytest

from chronoclade.cgmlst.matrix_thresholds import (
    matrix_threshold_sweep,
    minimum_spanning_forest,
    summarise_threshold_history,
)


def naive_partition(ids, matrix, cutoff, eligible=None):
    adjacency = {sample_id: set() for sample_id in ids}
    for i, first in enumerate(ids):
        for j, second in enumerate(ids[:i]):
            if (eligible is None or eligible[i, j]) and np.isfinite(matrix[i, j]):
                if matrix[i, j] <= cutoff:
                    adjacency[first].add(second)
                    adjacency[second].add(first)
    remaining = set(ids)
    groups = []
    while remaining:
        pending = [min(remaining)]
        component = set()
        while pending:
            sample_id = pending.pop()
            if sample_id not in component:
                component.add(sample_id)
                pending.extend(adjacency[sample_id] - component)
        remaining -= component
        groups.append(sorted(component))
    return sorted(groups)


def test_forest_matches_independent_dense_graph_zeros_missing_ties_and_memmap(tmp_path):
    rng = np.random.default_rng(414)
    ids = [f"sample-{i:02}" for i in range(24)]
    matrix = rng.choice([0., 0., 1., 2., 7., np.nan, np.inf], size=(24, 24))
    matrix = np.triu(matrix, 1) + np.triu(matrix, 1).T
    matrix[-1, :] = matrix[:, -1] = np.nan
    np.fill_diagonal(matrix, 0)
    mask = rng.random((24, 24)) > .3
    mask = np.triu(mask, 1) | np.triu(mask, 1).T
    path = tmp_path / "matrix.npy"
    np.save(path, matrix)
    mapped = np.load(path, mmap_mode="r")
    forest = minimum_spanning_forest(ids, mapped, eligible=mask)
    assert len(forest) <= len(ids) - 2
    for first, second, weight in forest:
        i, j = ids.index(first), ids.index(second)
        assert mask[i, j] and weight == matrix[i, j]
    result = matrix_threshold_sweep(ids, mapped, [0, .5, 1, 2, 7, 99], eligible=mask)
    for partition in result["partitions"]:
        assert partition["clusters"] == naive_partition(ids, matrix, partition["threshold"], mask)
    assert "n_edges" not in result  # Forest size must not masquerade as original pair coverage.
    assert result["n_forest_edges"] == len(forest)
    assert json.loads(json.dumps(result, allow_nan=False)) == result
    order = rng.permutation(len(ids))
    permuted_ids = [ids[i] for i in order]
    assert minimum_spanning_forest(
        permuted_ids, matrix[np.ix_(order, order)], eligible=mask[np.ix_(order, order)]
    ) == forest


def test_induced_sampling_recomputes_forest_instead_of_pruning():
    ids = ["c", "b", "a", "isolated"]
    matrix = np.array([[0, 1, 2, np.nan], [1, 0, 1, np.nan],
                       [2, 1, 0, np.nan], [np.nan, np.nan, np.nan, 0.]])
    assert minimum_spanning_forest(ids, matrix) == [("a", "b", 1.), ("b", "c", 1.)]
    subset = minimum_spanning_forest(ids, matrix, sample_indices=[0, 2, 3])
    assert subset == [("a", "c", 2.)]
    result = matrix_threshold_sweep(ids, matrix, [0, 1, 2], sample_indices=iter([0, 2, 3]))
    assert result["sample_ids"] == ["a", "c", "isolated"]
    assert result["partitions"][1]["clusters"] == [["a"], ["c"], ["isolated"]]
    assert result["partitions"][2]["clusters"] == [["a", "c"], ["isolated"]]
    # Masked zeros are absent edges, and do not silently join disconnected samples.
    zeros = np.zeros((3, 3), dtype=np.uint32)
    mask = np.zeros((3, 3), dtype=bool)
    mask[0, 1] = mask[1, 0] = True
    assert minimum_spanning_forest(["a", "b", "c"], zeros, eligible=mask) == [("a", "b", 0)]
    assert matrix_threshold_sweep([], np.zeros((0, 0)), [0])["partitions"][0]["clusters"] == []


def test_history_queries_match_memberships_without_materializing_partitions():
    ids = ["d", "b", "a", "c", "isolated"]
    matrix = np.array([[0, 3, 3, 3, np.nan], [3, 0, 0, 1, np.nan],
                       [3, 0, 0, 1, np.nan], [3, 1, 1, 0, np.nan],
                       [np.nan, np.nan, np.nan, np.nan, 0.]])
    compact = matrix_threshold_sweep(ids, matrix, [])
    assert compact["partitions"] == []
    thresholds = [0, .5, 1, 2, 3, 4]
    focal = ["d", "a", "isolated"]
    summaries = summarise_threshold_history(compact, reversed(thresholds), focal_sample_ids=focal)
    nodes = {i: {sample_id} for i, sample_id in enumerate(compact["sample_ids"])}
    for merge in compact["merge_history"]:
        nodes[merge["node_id"]] = set().union(*(nodes[child] for child in merge["children"]))
    for summary in summaries:
        expected = naive_partition(ids, matrix, summary["threshold"])
        assert summary["n_clusters"] == len(expected)
        assert summary["n_non_singleton_clusters"] == sum(len(group) > 1 for group in expected)
        assert summary["largest_cluster_size"] == max(map(len, expected))
        for sample_id, query in summary["focal_neighbourhoods"].items():
            group = next(set(group) for group in expected if sample_id in group)
            assert nodes[query["node_id"]] == group
            assert query["size"] == len(group)
            assert query["focal_samples"] == len(set(focal) & group)
    # Exact uint64 weights must not pass through imprecise floating-point best distances.
    high = 2 ** 63
    large = np.array([[0, high + 1, high], [high + 1, 0, high + 2],
                      [high, high + 2, 0]], dtype=np.uint64)
    exact = matrix_threshold_sweep(["a", "b", "c"], large, [high, high + 1])
    assert exact["partitions"][0]["clusters"] == [["a", "c"], ["b"]]
    assert exact["partitions"][1]["clusters"] == [["a", "b", "c"]]


def test_invalid_matrix_masks_ids_and_subset_fail_closed():
    base = np.array([[0., 1], [1, 0]])
    cases = [
        (["a", "a"], base, {}, "duplicate sample"),
        (["a", ""], base, {}, "nonempty strings"),
        (["a"], base, {}, "square numpy"),
        (["a", "b"], base.astype(complex), {}, "real numeric"),
        (["a", "b"], base, {"eligible": np.ones((2, 2), dtype=int)}, "boolean"),
        (["a", "b"], base, {"sample_indices": [0, 0]}, "duplicate sample_indices"),
        (["a", "b"], base, {"sample_indices": [2]}, "in-range"),
        (["a", "b"], base, {"sample_indices": [True]}, "in-range"),
        (["a", "b"], np.array([[0., -1], [-1, 0]]), {}, "nonnegative"),
        (["a", "b"], np.array([[0., 1], [2, 0]]), {}, "symmetric"),
        (["a", "b"], np.array([[0., np.nan], [1, 0]]), {}, "symmetric"),
        (["a", "b"], base, {"eligible": np.array([[True, False], [True, True]])}, "symmetric"),
    ]
    for ids, matrix, kwargs, message in cases:
        with pytest.raises(ValueError, match=message):
            minimum_spanning_forest(ids, matrix, **kwargs)
    compact = matrix_threshold_sweep(["a", "b"], base, [])
    for focal in [["missing"], ["a", "a"]]:
        with pytest.raises(ValueError, match="distinct known"):
            summarise_threshold_history(compact, [0], focal_sample_ids=focal)
    with pytest.raises(ValueError, match="threshold"):
        summarise_threshold_history(compact, [-1])
