"""Choose a clock-compatible root on a fixed phylogeny by residual MSE."""

from __future__ import annotations

import copy
import math
from collections.abc import Iterable, Mapping

import numpy as np
from Bio.Phylo.BaseTree import Clade, Tree

from chronoclade.metadata_dates import sample_date_interval


def optimise_clock_root(tree: Tree, samples: Iterable[Mapping]) -> tuple[Tree | None, dict]:
    """Minimise dated root-to-tip OLS residual MSE over eligible tree edges.

    Branch lengths and topology are held fixed. Partial collection dates use the
    midpoint returned by :func:`sample_date_interval`. Every nonnegative edge is
    searched continuously; negative branches remain in the tree but cannot host
    the root. A bifurcating root is treated as an artificial split of one edge.

    ``samples`` contains mappings with ``sample_id`` and date fields accepted by
    ``sample_date_interval``. Returns a rooted deep copy and an audit dictionary.
    A skipped fit returns None and an explicit reason.

    Runtime is O(E * V) for E candidate edges and V tree nodes, with O(V)
    working vectors. Descendant sets and edge labels can use O(V squared)
    memory on highly unbalanced trees. No dense distance matrix is built.
    """
    result_tree = copy.deepcopy(tree)
    samples = list(samples)
    terminals = list(tree.get_terminals())
    tip_names = [tip.name for tip in terminals]

    dated, excluded = _dated_tip_years(samples, set(tip_names))
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
        "search_complexity": "O(E*V) time; O(V) working vectors, up to O(V squared) descendant bookkeeping",
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
    audit["excluded_negative_edges"] = [
        {"edge_id": _edge_id(parent, child), "length": _branch_length(child),
         "reason": "negative_or_nonfinite_length"}
        for parent in tree.find_clades(order="preorder") for child in parent.clades
        if not math.isfinite(_branch_length(child)) or _branch_length(child) < 0
    ]
    audit["nonfinite_branch_edges"] = [
        {"edge_id": _edge_id(parent, child), "length": _branch_length(child)}
        for parent in tree.find_clades(order="preorder") for child in parent.clades
        if not math.isfinite(_branch_length(child))
    ]

    if audit["nonfinite_branch_edges"]:
        audit["reason"] = "nonfinite_branch_lengths"
        return None, audit

    if len(tip_names) < 3:
        audit["reason"] = "fewer_than_three_tree_tips"
        return None, audit
    if len(dated) < 3:
        audit["reason"] = "fewer_than_three_dated_tips"
        return None, audit

    sample_ids = sorted(dated)
    years = np.asarray([dated[name] for name in sample_ids], dtype=float)
    year_centered = years - years.mean()
    year_ss = float(np.dot(year_centered, year_centered))
    if not math.isfinite(year_ss) or year_ss <= 0:
        audit["reason"] = "no_collection_date_variation"
        return None, audit
    if len(set(years.tolist())) < 3:
        audit["reason"] = "fewer_than_three_distinct_dates"
        return None, audit

    adjacency, _ = _tree_graph(tree)
    sample_tip_names = set(sample_ids)
    descendants = {
        clade: {tip.name for tip in clade.get_terminals() if tip.name in sample_tip_names}
        for clade in tree.find_clades(order="postorder")
    }

    candidates = []
    root_children = list(tree.root.clades)
    collapsed_root = len(root_children) == 2
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            length = _branch_length(child)
            edge_id = _edge_id(parent, child)
            if not math.isfinite(length) or length < 0:
                continue
            if collapsed_root and parent is tree.root:
                continue
            candidates.append((parent, child, length, edge_id, False))

    if collapsed_root:
        left, right = root_children
        left_length, right_length = _branch_length(left), _branch_length(right)
        if (math.isfinite(left_length) and math.isfinite(right_length)
                and left_length >= 0 and right_length >= 0):
            length = left_length + right_length
            candidates.append((right, left, length,
                               f"{_clade_id(left)}|{_clade_id(right)}", True))
        else:
            audit["excluded_negative_edges"].append({
                "edge_id": f"{_clade_id(left)}|{_clade_id(right)}",
                "length": left_length + right_length,
                "reason": "bifurcating_root_contains_negative_or_nonfinite_edge",
            })

    audit["candidate_edge_count"] = len(candidates)
    best = None
    for parent, child, length, edge_id, collapsed in candidates:
        base = _distances_from_node(adjacency, parent, terminals)
        child_side = descendants[child]
        signs = np.asarray([-1.0 if name in child_side else 1.0 for name in sample_ids])
        y0 = np.asarray([base[name] for name in sample_ids], dtype=float)
        r0 = _residualise(y0, year_centered)
        rs = _residualise(signs, year_centered)
        quadratic = float(np.dot(rs, rs))
        position = 0.0 if quadratic <= 0 else -float(np.dot(r0, rs)) / quadratic
        position = min(max(position, 0.0), length)
        residual = r0 + position * rs
        mse = float(np.dot(residual, residual) / len(sample_ids))
        slope = float(np.dot(year_centered, y0 + position * signs) / year_ss)
        centered_y = y0 + position * signs
        total_ss = float(np.dot(centered_y - centered_y.mean(), centered_y - centered_y.mean()))
        r_squared = 1.0 - float(np.dot(residual, residual)) / total_ss if total_ss > 0 else None
        key = (mse, edge_id)
        if best is None or key < best[0]:
            best = (key, parent, child, length, edge_id, collapsed, position,
                    slope, r_squared, mse)

    if best is None:
        audit["reason"] = "no_nonnegative_edges"
        return None, audit

    (_, parent, child, length, edge_id, collapsed, position,
     slope, r_squared, mse) = best
    outgroup_length = length - position
    if collapsed:
        # Bio.Phylo drops the old bifurcating root and initially leaves the
        # opposite side with the full joined edge. Set its complementary share.
        result_left = _find_clade(tree, result_tree, child)
        result_tree.root_with_outgroup(result_left, outgroup_branch_length=outgroup_length)
        opposite = next(
            clade for clade in result_tree.root.clades if clade is not result_left
        )
        opposite.branch_length = position
    else:
        result_child = _find_clade(tree, result_tree, child)
        old_root = result_tree.root
        result_tree.root_with_outgroup(result_child, outgroup_branch_length=outgroup_length)
        if parent is tree.root:
            # Bio.Phylo preserves the whole old basal edge on the old root when
            # the selected edge directly borders a polytomous root. Split it.
            old_root.branch_length = position

    distances_preserved = _pairwise_distances_match(tree, result_tree)
    audit.update({
        "status": "rooted",
        "reason": None,
        "root_applied": True,
        "optimal_edge_id": edge_id,
        "optimal_edge_parent_id": _clade_id(parent),
        "optimal_edge_child_id": _clade_id(child),
        "root_position_from_parent": position,
        "edge_length": length,
        "slope": slope,
        "r_squared": r_squared,
        "mse": mse,
        "pairwise_distances_preserved": distances_preserved,
    })
    if not distances_preserved:
        audit.update({
            "status": "skipped",
            "reason": "pairwise_distance_check_failed",
            "root_applied": False,
        })
        return None, audit
    return result_tree, audit


def _dated_tip_years(samples: Iterable[Mapping], tree_tip_names: set[str]) -> tuple[dict, list[dict]]:
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
            excluded.append({"sample_id": sample_id, "reason": interval["reason"] or interval["status"]})
            continue
        by_id[sample_id] = float(interval["midpoint_year"])
    represented = {str(sample.get("sample_id")) for sample in samples if sample.get("sample_id")}
    excluded.extend({"sample_id": tip, "reason": "missing_metadata"}
                    for tip in sorted(tree_tip_names - represented))
    return by_id, excluded


def _tree_graph(tree: Tree) -> tuple[dict[Clade, list[tuple[Clade, float]]], dict[Clade, float]]:
    adjacency = {clade: [] for clade in tree.find_clades()}
    root_distances = {tree.root: 0.0}
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            length = _branch_length(child)
            adjacency[parent].append((child, length))
            adjacency[child].append((parent, length))
            root_distances[child] = root_distances[parent] + length
    return adjacency, root_distances


def _distances_from_node(adjacency, start: Clade, tips: list[Clade]) -> dict[str, float]:
    wanted = {tip.name for tip in tips}
    distances = {}
    stack = [(start, None, 0.0)]
    while stack:
        node, previous, distance = stack.pop()
        if node.is_terminal() and node.name in wanted:
            distances[node.name] = distance
        for neighbour, length in adjacency[node]:
            if neighbour is not previous:
                stack.append((neighbour, node, distance + length))
    return distances


def _residualise(values: np.ndarray, year_centered: np.ndarray) -> np.ndarray:
    centered = values - values.mean()
    denominator = float(np.dot(year_centered, year_centered))
    return centered - year_centered * (float(np.dot(year_centered, centered)) / denominator)


def _pairwise_distances_match(first: Tree, second: Tree) -> bool:
    first_tips = {tip.name: tip for tip in first.get_terminals()}
    second_tips = {tip.name: tip for tip in second.get_terminals()}
    if first_tips.keys() != second_tips.keys():
        return False
    first_graph, _ = _tree_graph(first)
    second_graph, _ = _tree_graph(second)
    for name, start in first_tips.items():
        one = _distances_from_node(first_graph, start, list(first_tips.values()))
        two = _distances_from_node(second_graph, second_tips[name], list(second_tips.values()))
        if any(not math.isclose(one[other], two[other], rel_tol=1e-9, abs_tol=1e-9)
               for other in first_tips):
            return False
    return True


def _find_clade(source: Tree, target_tree: Tree, target: Clade) -> Clade:
    path = _clade_path(source.root, target)
    clade = target_tree.root
    for index in path:
        clade = clade.clades[index]
    return clade


def _clade_path(root: Clade, target: Clade) -> list[int]:
    def search(current: Clade, path: list[int]) -> list[int] | None:
        if current is target:
            return path
        for index, child in enumerate(current.clades):
            found = search(child, [*path, index])
            if found is not None:
                return found
        return None

    found = search(root, [])
    if found is None:
        raise ValueError("Chosen outgroup clade is not in the input tree")
    return found


def _clade_id(clade: Clade) -> str:
    if clade.name:
        return str(clade.name)
    names = sorted(str(tip.name) for tip in clade.get_terminals())
    return "{" + ",".join(names) + "}"


def _edge_id(parent: Clade, child: Clade) -> str:
    return f"{_clade_id(parent)}->{_clade_id(child)}"


def _branch_length(clade: Clade) -> float:
    return float(clade.branch_length or 0.0)
