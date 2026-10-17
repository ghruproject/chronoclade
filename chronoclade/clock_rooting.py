"""Continuous-edge root-to-tip OLS rooting with linear sufficient-statistic DP.

The objective follows the residual-MSE rooting option in Clockor2
(Featherstone et al., 2024, doi:10.1093/sysbio/syae003):
https://github.com/clockor2/clockor2/blob/662a997970f1326cc584012330369546de0f92da/src/features/engine/bestFittingRoot.ts
The tree DP and continuous minimisation below are independently derived.
"""

from __future__ import annotations

import copy
import math
from collections.abc import Iterable, Mapping

import numpy as np
from Bio.Phylo.BaseTree import Clade, Tree

from chronoclade.metadata_dates import sample_date_interval


def optimise_clock_root(tree: Tree, samples: Iterable[Mapping]) -> tuple[Tree | None, dict]:
    """Minimise dated root-to-tip OLS residual MSE over eligible tree edges.

    Dates use the midpoint from :func:`sample_date_interval`. All nonnegative
    edges are searched continuously; signed negative lengths are retained but
    cannot host the root. A bifurcating root is an artificial split of one edge;
    its joined edge is eligible only when both original halves are nonnegative.

    A postorder pass accumulates dated-tip count, centred-date sum, distance
    sum, distance-square sum and date-distance product for each subtree. A
    preorder reroot pass obtains the three whole-tree distance moments at every
    node. Those moments define an OLS residual quadratic along each edge, so
    searching, copying, rerooting and checking the unrooted topology/lengths
    take O(V) time and memory. Sorting date IDs for the audit adds O(T log T).
    IndexedTree inputs additionally rebuild their O(V log V) LCA index.
    Unnamed clades have compact deterministic preorder IDs instead of expanded
    descendant labels. No per-node descendant sets or pairwise distance matrix
    are used.

    Returns a rooted independent copy and an audit. A skipped fit returns None
    and an explicit reason. The input tree and metadata are never modified.
    """
    nodes = _preorder(tree.root)
    terminals = [node for node in nodes if not node.clades]
    tip_names = [tip.name for tip in terminals]
    dated, excluded = _dated_tip_years(list(samples), set(tip_names))
    ids = {
        node: str(node.name) if node.name else f"node:{index}" for index, node in enumerate(nodes)
    }
    audit = {
        "method": "fixed_topology_continuous_clock_root_ols_residual_mse",
        "status": "skipped",
        "reason": None,
        "root_applied": False,
        "dated_tip_count": len(dated),
        "dated_tip_ids": sorted(dated),
        "excluded_dates": excluded,
        "excluded_negative_edges": [],
        "nonfinite_branch_edges": [],
        "candidate_edge_count": 0,
        "search_complexity": "O(V) rooting time and memory; O(T log T) audit sorting; "
        "O(V log V) when rebuilding an indexed-tree LCA cache",
        "optimal_edge_id": None,
        "optimal_edge_parent_id": None,
        "optimal_edge_child_id": None,
        "root_position_from_parent": None,
        "edge_length": None,
        "slope": None,
        "r_squared": None,
        "mse": None,
        "pairwise_distances_preserved": None,
    }
    for parent in nodes:
        for child in parent.clades:
            length = _branch_length(child)
            edge_id = f"{ids[parent]}->{ids[child]}"
            if not math.isfinite(length) or length < 0:
                audit["excluded_negative_edges"].append(
                    {
                        "edge_id": edge_id,
                        "length": length,
                        "reason": "negative_or_nonfinite_length",
                    }
                )
            if not math.isfinite(length):
                audit["nonfinite_branch_edges"].append({"edge_id": edge_id, "length": length})
    if audit["nonfinite_branch_edges"]:
        audit["reason"] = "nonfinite_branch_lengths"
        return None, audit
    if len(tip_names) < 3:
        audit["reason"] = "fewer_than_three_tree_tips"
        return None, audit
    if len(dated) < 3:
        audit["reason"] = "fewer_than_three_dated_tips"
        return None, audit

    # Centre years and use stable centred moment merges. On platforms where
    # longdouble is wider than float64, its extra precision is also useful.
    sample_ids = sorted(dated)
    years = np.asarray([dated[name] for name in sample_ids], dtype=np.longdouble)
    year_centered = years - years.mean()
    year_ss = np.dot(year_centered, year_centered)
    if not np.isfinite(year_ss) or year_ss <= 0:
        audit["reason"] = "no_collection_date_variation"
        return None, audit
    if len(set(years.tolist())) < 3:
        audit["reason"] = "fewer_than_three_distinct_dates"
        return None, audit
    date_values = dict(zip(sample_ids, year_centered))
    subtree, totals, outside = _distance_moments(nodes, date_values)
    count = len(sample_ids)

    candidates = []
    root_children = tree.root.clades
    collapsed_root = len(root_children) == 2
    for parent in nodes:
        for child in parent.clades:
            length = _branch_length(child)
            if length >= 0 and not (collapsed_root and parent is tree.root):
                candidates.append((parent, child, length, f"{ids[parent]}->{ids[child]}"))
    if collapsed_root:
        left, right = root_children
        left_length, right_length = _branch_length(left), _branch_length(right)
        edge_id = f"{ids[left]}|{ids[right]}"
        if left_length >= 0 and right_length >= 0:
            candidates.append((right, left, left_length + right_length, edge_id))
        else:
            audit["excluded_negative_edges"].append(
                {
                    "edge_id": edge_id,
                    "length": left_length + right_length,
                    "reason": "bifurcating_root_contains_negative_or_nonfinite_edge",
                }
            )
    audit["candidate_edge_count"] = len(candidates)

    best = None
    for parent, child, length, edge_id in candidates:
        inside = subtree[child]
        other = outside[child]
        n, x, y, xx, yy, xy = inside
        m, ox, oy, oxx, oyy, oxy = other
        # Moving z towards the child shifts its distances by -z and all
        # others by +z. Within-side variances/covariances stay constant.
        # Let w=n*m/N, dx=mean(x_inside)-mean(x_other), and dy the
        # corresponding distance difference at the parent. Decomposing OLS
        # into within-side and between-side sums gives:
        #   z* = (dy - (C_within/X_within)*dx)/2
        #   SSE(z) = Y_within - C_within**2/X_within + q*(z-z*)**2
        #   q = 4*w*X_within / X_total.
        # This is the exact residual quadratic, with no subtraction of large
        # raw square sums and no cancellation when q is almost zero.
        if n and m:
            within_x, within_y, within_c = xx + oxx, yy + oyy, xy + oxy
            dx, dy = x - ox, y - oy + 2 * length
            optimum = (dy - within_c * dx / within_x) / 2
            position = np.clip(optimum, 0, length)
            q = (np.longdouble(4 * n * m) / count) * within_x / year_ss
            minimum = within_y - within_c**2 / within_x
            sse = max(np.longdouble(0), minimum + q * (position - optimum) ** 2)
            scale = max(
                abs(within_y),
                abs(within_c**2 / within_x),
                abs(q * (position - optimum) ** 2),
                np.finfo(np.longdouble).tiny,
            )
        else:
            # All dated tips lie on one side: root movement only changes the
            # fitted intercept, so choose the parent endpoint deterministically.
            _, _, _, _, distance_variance, date_distance = totals[parent]
            position = np.longdouble(0)
            sse = max(np.longdouble(0), distance_variance - date_distance**2 / year_ss)
            scale = max(
                abs(distance_variance),
                abs(date_distance**2 / year_ss),
                np.finfo(np.longdouble).tiny,
            )
        # Numerical ties are resolved by stable compact IDs, including endpoint
        # roots shared by several edges. The tolerance scales with the moments.
        error = 128 * np.finfo(np.longdouble).eps * scale
        if (
            best is None
            or sse < best[0] - max(error, best[1])
            or (abs(sse - best[0]) <= max(error, best[1]) and edge_id < best[5])
        ):
            best = (sse, error, parent, child, length, edge_id, position)
    if best is None:
        audit["reason"] = "no_nonnegative_edges"
        return None, audit
    _, _, parent, child, length, edge_id, position = best

    # Evaluate the winning edge once with explicit residual vectors, avoiding
    # cancellation in the reported SSE and refining its stationary point.
    adjacency, _ = _tree_graph(tree)
    base = _distances_from_node(adjacency, parent, terminals)
    child_names = {node.name for node in _preorder(child) if not node.clades}
    signs = np.asarray(
        [-1 if name in child_names else 1 for name in sample_ids], dtype=np.longdouble
    )
    y0 = np.asarray([base[name] for name in sample_ids], dtype=np.longdouble)
    r0, rs = _residualise(y0, year_centered), _residualise(signs, year_centered)
    quadratic = np.dot(rs, rs)
    position = (
        np.longdouble(0) if quadratic <= 0 else np.clip(-np.dot(r0, rs) / quadratic, 0, length)
    )
    residual = r0 + position * rs
    sse = np.dot(residual, residual)
    values = y0 + position * signs
    total_ss = np.dot(values - values.mean(), values - values.mean())
    slope = float(np.dot(year_centered, values - values.mean()) / year_ss)
    r_squared = float(1 - sse / total_ss) if total_ss > 0 else None
    position = float(position)

    result_tree, clones = _root_copy(tree, nodes, adjacency, parent, child, length, position)
    distances_preserved = _pairwise_distances_match(tree, result_tree, clones)
    audit.update(
        {
            "status": "rooted",
            "root_applied": True,
            "optimal_edge_id": edge_id,
            "optimal_edge_parent_id": ids[parent],
            "optimal_edge_child_id": ids[child],
            "root_position_from_parent": position,
            "edge_length": length,
            "slope": slope,
            "r_squared": r_squared,
            "mse": float(sse / count),
            "pairwise_distances_preserved": distances_preserved,
        }
    )
    if not distances_preserved:
        audit.update(
            {"status": "skipped", "reason": "pairwise_distance_check_failed", "root_applied": False}
        )
        return None, audit
    return result_tree, audit


def _merge_moments(left, right):
    """Stable merge of (count, date mean, distance mean, date M2, distance M2, C).

    M2 is the unnormalised centred sum of squares; C is the centred cross sum.
    These are equivalent to n, sum(x), sum(y), sum(x*x), sum(y*y), sum(x*y),
    but do not subtract large calendar-year or root-distance offsets.
    """
    ln, lx, ly, lxx, lyy, lxy = left
    rn, rx, ry, rxx, ryy, rxy = right
    if not ln:
        return right
    if not rn:
        return left
    n = ln + rn
    weight = np.longdouble(ln * rn) / n
    dx, dy = rx - lx, ry - ly
    return (
        n,
        lx + dx * rn / n,
        ly + dy * rn / n,
        lxx + rxx + weight * dx**2,
        lyy + ryy + weight * dy**2,
        lxy + rxy + weight * dx * dy,
    )


def _shift_moments(moments, length):
    n, x, y, xx, yy, xy = moments
    return (n, x, y + length, xx, yy, xy) if n else moments


def _distance_moments(nodes, dates):
    """Compute subtree and whole-tree distance moments in two linear passes.

    A branch adds its signed length to the side's mean distance and leaves its
    centred moments unchanged. Postorder merges give each downward side.
    Preorder prefix/suffix merges give each complementary side, shifted across
    the branch to the child. Each side is merged a constant number of times,
    including at polytomies. No unstable subtraction of whole-tree raw square
    sums or repeated enumeration of a node's descendants is needed.
    """
    zero = np.longdouble(0)
    empty = (0, zero, zero, zero, zero, zero)
    subtree = {}
    for node in reversed(nodes):
        moments = (
            (1, dates[node.name], zero, zero, zero, zero)
            if (not node.clades and node.name in dates)
            else empty
        )
        for child in node.clades:
            moments = _merge_moments(
                moments, _shift_moments(subtree[child], np.longdouble(_branch_length(child)))
            )
        subtree[node] = moments
    outside = {nodes[0]: empty}
    totals = {}
    for parent in nodes:
        totals[parent] = _merge_moments(subtree[parent], outside[parent])
        sides = [
            _shift_moments(subtree[child], np.longdouble(_branch_length(child)))
            for child in parent.clades
        ]
        prefix = [outside[parent]]
        for side in sides:
            prefix.append(_merge_moments(prefix[-1], side))
        suffix = [empty] * (len(sides) + 1)
        for index in range(len(sides) - 1, -1, -1):
            suffix[index] = _merge_moments(sides[index], suffix[index + 1])
        for index, child in enumerate(parent.clades):
            complement = _merge_moments(prefix[index], suffix[index + 1])
            outside[child] = _shift_moments(complement, np.longdouble(_branch_length(child)))
    return subtree, totals, outside


def _dated_tip_years(
    samples: Iterable[Mapping], tree_tip_names: set[str]
) -> tuple[dict, list[dict]]:
    by_id = {}
    excluded = []
    counts = {}
    for sample in samples:
        sample_id = sample.get("sample_id")
        if sample_id:
            sample_id = str(sample_id)
            counts[sample_id] = counts.get(sample_id, 0) + 1
    for sample in samples:
        sample_id = sample.get("sample_id")
        if not sample_id:
            excluded.append({"sample_id": None, "reason": "missing_sample_id"})
            continue
        sample_id = str(sample_id)
        if sample_id not in tree_tip_names:
            excluded.append({"sample_id": sample_id, "reason": "sample_not_in_tree"})
            continue
        if counts[sample_id] > 1:
            excluded.append({"sample_id": sample_id, "reason": "duplicate_sample_id"})
            continue
        interval = sample_date_interval(sample)
        if interval["status"] != "valid":
            excluded.append(
                {"sample_id": sample_id, "reason": interval["reason"] or interval["status"]}
            )
            continue
        by_id[sample_id] = float(interval["midpoint_year"])
    represented = {str(sample.get("sample_id")) for sample in samples if sample.get("sample_id")}
    excluded.extend(
        {"sample_id": tip, "reason": "missing_metadata"}
        for tip in sorted(tree_tip_names - represented)
    )
    return by_id, excluded


def _preorder(root: Clade) -> list[Clade]:
    """Iterative traversal also supports trees deeper than Python's stack."""
    nodes, stack = [], [root]
    while stack:
        node = stack.pop()
        nodes.append(node)
        stack.extend(reversed(node.clades))
    return nodes


def _tree_graph(tree: Tree) -> tuple[dict[Clade, list[tuple[Clade, float]]], dict[Clade, float]]:
    nodes = _preorder(tree.root)
    adjacency = {node: [] for node in nodes}
    root_distances = {tree.root: 0.0}
    for parent in nodes:
        for child in parent.clades:
            length = _branch_length(child)
            adjacency[parent].append((child, length))
            adjacency[child].append((parent, length))
            root_distances[child] = root_distances[parent] + length
    return adjacency, root_distances


def _distances_from_node(adjacency, start: Clade, tips: list[Clade]) -> dict[str, float]:
    wanted = {tip.name for tip in tips}
    distances = {}
    stack = [(start, None, np.longdouble(0))]
    while stack:
        node, previous, distance = stack.pop()
        if not node.clades and node.name in wanted:
            distances[node.name] = distance
        for neighbour, length in adjacency[node]:
            if neighbour is not previous:
                stack.append((neighbour, node, distance + length))
    return distances


def _residualise(values: np.ndarray, year_centered: np.ndarray) -> np.ndarray:
    centered = values - values.mean()
    denominator = np.dot(year_centered, year_centered)
    return centered - year_centered * (np.dot(year_centered, centered) / denominator)


def _suppress_root(adjacency, root):
    """Remove only an artificial degree-two root, joining its signed halves."""
    if len(adjacency[root]) != 2:
        return
    (left, left_length), (right, right_length) = adjacency.pop(root)
    adjacency[left] = [(node, length) for node, length in adjacency[left] if node is not root]
    adjacency[right] = [(node, length) for node, length in adjacency[right] if node is not root]
    adjacency[left].append((right, left_length + right_length))
    adjacency[right].append((left, left_length + right_length))


def _root_copy(tree, nodes, adjacency, parent, child, length, position):
    """Clone and orient the unrooted graph without recursive Bio.Phylo APIs."""
    # Exclude clades/root from deepcopy: recursively copying a comb would hit
    # Python's recursion limit even though the numerical search is linear.
    clones = {}
    memo = {}
    for node in nodes:
        clone = copy.copy(node)
        clones[node] = clone
        memo[id(node)] = clone
    for node in nodes:
        clone = clones[node]
        clone.__dict__ = copy.deepcopy(
            {key: value for key, value in node.__dict__.items() if key != "clades"}, memo
        )
        clone.clades = []
    result = copy.copy(tree)
    # IndexedTree's topology-dependent LCA cache must be rebuilt, not copied.
    # Its cached node keys would otherwise describe the old rooting.
    index_fields = {"_nodes", "_node_index", "_name_index", "_depths", "_levels", "_up"}
    excluded = {"root"} | (index_fields if callable(getattr(tree, "reindex", None)) else set())
    result.__dict__ = copy.deepcopy(
        {key: value for key, value in tree.__dict__.items() if key not in excluded}, memo
    )
    graph = {node: list(neighbours) for node, neighbours in adjacency.items()}
    _suppress_root(graph, tree.root)
    new_root = Clade(branch_length=tree.root.branch_length)
    new_root.clades = [clones[parent], clones[child]]
    result.root, result.rooted = new_root, True
    stack = [(parent, child, position), (child, parent, length - position)]
    while stack:
        node, previous, incoming_length = stack.pop()
        clone = clones[node]
        clone.branch_length = incoming_length
        for neighbour, neighbour_length in graph[node]:
            if neighbour is not previous:
                clone.clades.append(clones[neighbour])
                stack.append((neighbour, node, neighbour_length))
    if callable(getattr(result, "reindex", None)):
        result.reindex()
    return result, clones


def _pairwise_distances_match(first: Tree, second: Tree, clones: dict[Clade, Clade]) -> bool:
    """Prove all pairwise distances agree by checking the unrooted edges.

    The copy's explicit node correspondence permits an exact O(V) topology
    comparison. Suppress the artificial roots, then compare every original
    node, tip label and signed edge length. Identical paths imply identical
    pairwise sums, without enumerating O(T squared) tip pairs.
    """
    first_graph, _ = _tree_graph(first)
    second_graph, _ = _tree_graph(second)
    if {clones[node] for node in first_graph if not node.clades} != {
        node for node in second_graph if not node.clades
    }:
        return False
    _suppress_root(first_graph, first.root)
    _suppress_root(second_graph, second.root)
    reverse = {clone: original for original, clone in clones.items()}
    if len(first_graph) != len(second_graph):
        return False
    if any(node not in reverse for node in second_graph):
        return False
    if {reverse[node] for node in second_graph} != first_graph.keys():
        return False
    # Build neighbour dictionaries once. Repeated linear scans of a star's
    # adjacency list would itself be quadratic.
    for original, neighbours in first_graph.items():
        clone = clones[original]
        if clone.name != original.name:
            return False
        actual = {reverse[node]: length for node, length in second_graph[clone]}
        if actual.keys() != {node for node, _ in neighbours}:
            return False
        if any(
            not math.isclose(length, actual[node], rel_tol=1e-9, abs_tol=1e-9)
            for node, length in neighbours
        ):
            return False
    return True


def _branch_length(clade: Clade) -> float:
    return float(clade.branch_length or 0.0)
