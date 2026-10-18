"""Single-linkage sweeps preserve raw boundaries and incomplete evidence."""

import itertools
import json
import random

import pytest

from chronoclade.cgmlst.thresholds import threshold_sweep


def naive_components(sample_ids, edges, threshold):
    """Independent breadth-first graph traversal, rebuilt at each threshold."""
    neighbours = {sample_id: set() for sample_id in sample_ids}
    for first, second, distance in edges:
        if distance <= threshold:
            neighbours[first].add(second)
            neighbours[second].add(first)
    remaining = set(sample_ids)
    clusters = []
    while remaining:
        pending = [min(remaining)]
        component = set()
        while pending:
            sample_id = pending.pop()
            if sample_id not in component:
                component.add(sample_id)
                pending.extend(neighbours[sample_id] - component)
        remaining -= component
        clusters.append(sorted(component))
    return sorted(clusters)


def test_sweep_matches_independent_graph_at_raw_boundaries():
    rng = random.Random(58)
    sample_ids = [f"sample-{index:02}" for index in range(15)]
    edges = [
        (first, second, rng.choice([0, 1, 2.5, 7, 11]))
        for first, second in itertools.combinations(sample_ids, 2)
        if rng.random() < 0.14
    ]
    thresholds = [11, 0, 0.5, 1, 2, 2.5, 3, 7, 8, 11, 20]
    result = threshold_sweep(iter(reversed(sample_ids)), iter(edges), iter(thresholds),
                             distance_unit="raw nucleotide differences")
    assert result["thresholds"] == sorted(set(thresholds))
    assert result["distance_unit"] == "raw nucleotide differences"
    for partition in result["partitions"]:
        expected = naive_components(sample_ids, edges, partition["threshold"])
        assert partition["clusters"] == expected
        assert partition["n_clusters"] == len(expected)
        assert partition["n_non_singleton_clusters"] == sum(len(c) > 1 for c in expected)
        assert partition["largest_cluster_size"] == max(map(len, expected))
    assert json.loads(json.dumps(result, allow_nan=False)) == result
    # Raw integer counts retain their exact <= boundary without float rounding.
    huge = 10 ** 400
    exact = threshold_sweep(["a", "b"], [("a", "b", huge + 1)], [huge, huge + 1])
    assert exact["partitions"][0]["clusters"] == [["a"], ["b"]]
    assert exact["partitions"][1]["clusters"] == [["a", "b"]]


def test_ties_have_simultaneous_deterministic_merge_history():
    sample_ids = ["d", "b", "a", "c", "isolated"]
    edges = [("a", "b", 1), ("b", "c", 1), ("a", "c", 1), ("c", "d", 3)]
    expected = threshold_sweep(sample_ids, edges, [0, 1, 2, 3, 4])
    for permutation in itertools.permutations(edges):
        reversed_edges = [(second, first, weight) for first, second, weight in permutation]
        assert threshold_sweep(reversed(sample_ids), reversed_edges, [4, 3, 2, 1, 0]) == expected
    assert expected["merge_history"] == [
        {"distance": 1, "node_id": 5, "children": [0, 1, 2], "size": 3},
        {"distance": 3, "node_id": 6, "children": [3, 5], "size": 4},
    ]
    assert expected["partitions"][1]["clusters"] == [["a", "b", "c"], ["d"], ["isolated"]]
    assert expected["stable_ranges"] == [
        {"start_threshold": 0, "end_threshold": 0, "threshold_count": 1, "partition_index": 0},
        {"start_threshold": 1, "end_threshold": 2, "threshold_count": 2, "partition_index": 1},
        {"start_threshold": 3, "end_threshold": 4, "threshold_count": 2, "partition_index": 3},
    ]
    # Histories reference earlier nodes, and describe the same components as BFS.
    nodes = {index: {sample_id} for index, sample_id in enumerate(expected["sample_ids"])}
    active = set(nodes)
    for merge in expected["merge_history"]:
        assert set(merge["children"]) <= active
        nodes[merge["node_id"]] = set().union(*(nodes[i] for i in merge["children"]))
        assert len(nodes[merge["node_id"]]) == merge["size"]
        active.difference_update(merge["children"])
        active.add(merge["node_id"])
        assert sorted(sorted(nodes[i]) for i in active) == naive_components(
            sample_ids, edges, merge["distance"]
        )


def test_missing_pairs_singletons_empty_graphs_and_default_weights():
    result = threshold_sweep(["a", "b", "c", "d"], [("a", "b", 2), ("b", "c", 5)])
    assert result["thresholds"] == [2, 5]
    assert result["partitions"][-1]["clusters"] == [["a", "b", "c"], ["d"]]
    assert result["n_edges"] == 2
    # No a-c edge is required: a single-linkage chain joins them at 5.
    identical = threshold_sweep(["b", "a"], [("a", "b", 0), ("b", "a", 0)], [0])
    assert identical["n_edges"] == 1
    assert identical["partitions"][0]["clusters"] == [["a", "b"]]
    empty = threshold_sweep(["b", "a"], [])
    assert empty["initial_clusters"] == [["a"], ["b"]]
    assert empty["partitions"] == empty["merge_history"] == empty["stable_ranges"] == []
    assert threshold_sweep(["b", "a"], [], [0, 99])["partitions"][-1]["clusters"] == [
        ["a"], ["b"]
    ]
    assert threshold_sweep([], [], [0])["partitions"] == [{
        "threshold": 0, "clusters": [], "n_clusters": 0,
        "n_non_singleton_clusters": 0, "largest_cluster_size": 0,
    }]
    # Requesting no snapshots still gives the complete observed merge history.
    assert threshold_sweep(["a", "b"], [("a", "b", 5)], [])["merge_history"] == [
        {"distance": 5, "node_id": 2, "children": [0, 1], "size": 2}
    ]


def test_invalid_ids_edges_distances_thresholds_and_units_fail_closed():
    invalid = [
        (["a", "a"], [], None, "duplicate sample"),
        ([""], [], None, "nonempty strings"),
        ([1], [], None, "nonempty strings"),
        (["a", "b"], [("a", "missing", 1)], None, "unknown sample"),
        (["a", "b"], [([], "b", 1)], None, "unknown sample"),
        (["a"], [("a", "a", 0)], None, "self-edge"),
        (["a", "b"], [("a", "b", 1), ("b", "a", 2)], None, "conflicting duplicate"),
        (["a", "b"], [("a", "b")], None, "each edge"),
    ]
    for value in [-1, float("nan"), float("inf"), -float("inf"), True, "1", None]:
        invalid.append((["a", "b"], [("a", "b", value)], None, "edge distance"))
        invalid.append((["a", "b"], [], [value], "threshold"))
    for sample_ids, edges, thresholds, message in invalid:
        with pytest.raises(ValueError, match=message):
            threshold_sweep(sample_ids, edges, thresholds)
    with pytest.raises(ValueError, match="distance_unit"):
        threshold_sweep([], [], distance_unit=" ")
