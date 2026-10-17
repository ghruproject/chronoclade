import io

import numpy as np
import pytest
from Bio import Phylo

import chronoclade.clock_rooting as clock_rooting
from chronoclade.clock_rooting import optimise_clock_root
from chronoclade.metadata_dates import sample_date_interval


def tree(text):
    return Phylo.read(io.StringIO(text), "newick")


def samples(dates):
    return [{"sample_id": name, "collection_date": date} for name, date in dates.items()]


def _mse(dates, distances):
    years = np.asarray(
        [
            sample_date_interval({"collection_date": date})["midpoint_year"]
            for date in dates.values()
        ]
    )
    values = np.asarray([distances[name] for name in dates])
    design = np.column_stack((np.ones(len(years)), years))
    residual = values - design @ np.linalg.lstsq(design, values, rcond=None)[0]
    return float(np.mean(residual**2))


def _brute_grid_mse(source, dates, steps=2000):
    best = float("inf")
    for parent in source.find_clades():
        for child in parent.clades:
            length = float(child.branch_length or 0)
            if length < 0:
                continue
            descendants = {tip.name for tip in child.get_terminals()}
            for position in np.linspace(0, length, steps + 1):
                distances = {}
                for sample_id in dates:
                    base = source.distance(parent, sample_id)
                    distances[sample_id] = base + (
                        -position if sample_id in descendants else position
                    )
                best = min(best, _mse(dates, distances))
    return best


def test_analytic_optimum_matches_dense_brute_force_grid():
    original = tree("(a:0.2,b:0.5,(c:1,d:0.7):0.8,e:0.3);")
    dates = {"a": "2001", "b": "2007-05", "c": "2012-03-04", "d": "2017", "e": "2022-09-15"}
    expected_grid_mse = _brute_grid_mse(original, dates)

    rooted, audit = optimise_clock_root(original, samples(dates))

    assert rooted is not None
    assert audit["status"] == "rooted"
    assert audit["mse"] <= expected_grid_mse + 1e-8
    assert audit["mse"] == pytest.approx(expected_grid_mse, abs=1e-6)
    assert audit["optimal_edge_id"]
    assert 0 <= audit["root_position_from_parent"] <= audit["edge_length"]
    assert audit["slope"] is not None
    assert audit["r_squared"] is not None


def test_rooting_returns_deepcopy_and_preserves_all_pairwise_distances():
    original = tree("(a:0.2,b:0.5,(c:1,d:0.7):0.8,e:0.3);")
    before = {(a, b): original.distance(a, b) for a in "abcde" for b in "abcde"}
    dates = {"a": "2001", "b": "2007", "c": "2012", "d": "2017", "e": "2022"}

    rooted, audit = optimise_clock_root(original, samples(dates))

    assert rooted is not original
    assert audit["pairwise_distances_preserved"] is True
    for pair, distance in before.items():
        assert rooted.distance(*pair) == pytest.approx(distance)
        assert original.distance(*pair) == pytest.approx(distance)


def test_bifurcating_root_artificial_edge_is_searched_and_preserved():
    original = tree("((a:1,b:1):2,(c:1,d:1):2);")
    before = {(a, b): original.distance(a, b) for a in "abcd" for b in "abcd"}
    dates = {"a": "2001", "b": "2012", "c": "2015", "d": "2024"}

    rooted, audit = optimise_clock_root(original, samples(dates))

    assert rooted is not None
    assert audit["status"] == "rooted"
    assert "|" in audit["optimal_edge_id"]
    assert audit["edge_length"] == pytest.approx(4)
    assert audit["pairwise_distances_preserved"] is True
    for pair, distance in before.items():
        assert rooted.distance(*pair) == pytest.approx(distance)


def test_missing_invalid_and_unmatched_dates_are_audited_and_skip_with_fewer_than_three():
    original = tree("(a:1,b:1,c:1,d:1);")
    metadata = [
        {"sample_id": "a", "collection_date": "2018"},
        {"sample_id": "b", "collection_date": "2020-02"},
        {"sample_id": "c", "collection_date": "not-a-date"},
        {"sample_id": "outside", "collection_date": "2021"},
    ]

    rooted, audit = optimise_clock_root(original, metadata)

    assert rooted is None
    assert audit["reason"] == "fewer_than_three_dated_tips"
    assert {item["reason"] for item in audit["excluded_dates"]} >= {
        "invalid_collection_date",
        "sample_not_in_tree",
        "missing_metadata",
    }


def test_no_date_variation_skips_and_negative_edges_are_excluded_but_retained():
    original = tree("(a:-1,b:0,c:1,d:2);")
    dates = {name: "2020-06-15" for name in "abcd"}

    rooted, audit = optimise_clock_root(original, samples(dates))

    assert rooted is None
    assert audit["reason"] == "no_collection_date_variation"
    assert len(audit["excluded_negative_edges"]) == 1

    dates = {"a": "2010", "b": "2015", "c": "2020", "d": "2025"}
    rooted, audit = optimise_clock_root(original, samples(dates))
    assert rooted is not None
    assert len(audit["excluded_negative_edges"]) == 1
    assert next(tip for tip in rooted.get_terminals() if tip.name == "a").branch_length == -1


def test_two_distinct_midpoints_skip_even_when_three_tips_are_dated():
    original = tree("(a:1,b:1,c:1);")
    dates = {"a": "2020-01-01", "b": "2020-01-01", "c": "2021-01-01"}

    rooted, audit = optimise_clock_root(original, samples(dates))

    assert rooted is None
    assert audit["dated_tip_count"] == 3
    assert audit["reason"] == "fewer_than_three_distinct_dates"


def test_nonfinite_branch_lengths_skip_before_they_can_contaminate_fit():
    original = tree("(a:1,b:1,c:1,d:1);")
    original.find_any(name="a").branch_length = float("nan")
    dates = {"a": "2010", "b": "2015", "c": "2020", "d": "2025"}

    rooted, audit = optimise_clock_root(original, samples(dates))

    assert rooted is None
    assert audit["reason"] == "nonfinite_branch_lengths"
    assert len(audit["nonfinite_branch_edges"]) == 1


def test_pairwise_distance_check_failure_does_not_return_a_root(monkeypatch):
    original = tree("(a:1,b:1,c:1,d:1);")
    dates = {"a": "2010", "b": "2015", "c": "2020", "d": "2025"}
    monkeypatch.setattr(clock_rooting, "_pairwise_distances_match", lambda *_: False)

    rooted, audit = optimise_clock_root(original, samples(dates))

    assert rooted is None
    assert audit["status"] == "skipped"
    assert audit["reason"] == "pairwise_distance_check_failed"
    assert audit["pairwise_distances_preserved"] is False


def _continuous_reference(source, dates):
    """Independent explicit-distance OLS + bounded SciPy minimisation."""
    minimize_scalar = pytest.importorskip("scipy.optimize").minimize_scalar
    years = np.asarray(
        [
            sample_date_interval({"collection_date": value})["midpoint_year"]
            for value in dates.values()
        ]
    )
    years -= years.mean()
    design = np.column_stack((np.ones(len(years)), years))
    candidates = [
        (parent, child, float(child.branch_length or 0))
        for parent in source.find_clades()
        for child in parent.clades
        if (child.branch_length or 0) >= 0
        and not (parent is source.root and len(source.root.clades) == 2)
    ]
    if len(source.root.clades) == 2:
        left, right = source.root.clades
        if (left.branch_length or 0) >= 0 and (right.branch_length or 0) >= 0:
            candidates.append((right, left, (left.branch_length or 0) + (right.branch_length or 0)))
    best = float("inf")
    for parent, child, length in candidates:
        child_names = {tip.name for tip in child.get_terminals()}
        base = np.asarray([source.distance(parent, name) for name in dates])
        signs = np.asarray([-1 if name in child_names else 1 for name in dates])

        def objective(position):
            values = base + position * signs
            residual = values - design @ np.linalg.lstsq(design, values, rcond=None)[0]
            return float(np.mean(residual**2))

        best = min(best, objective(0), objective(length))
        if length:
            result = minimize_scalar(
                objective, bounds=(0, length), method="bounded", options={"xatol": 1e-13}
            )
            best = min(best, result.fun)
    return best


@pytest.mark.parametrize(
    "newick",
    [
        "(((a:0.2,b:0):0.3,c:0.7):0.1,(d:0.3,e:0.8):0.4);",
        "(a:0,b:0,(c:0,d:0,e:0):0);",
        "(a:-0.1,b:0.5,(c:0.3,d:-0.2,e:0.4):0.7);",
        "(((a:0.1,b:0.2):-0.05,c:0.3):0.1,(d:0.7,e:0):0.9);",
        # One original basal half is negative; its joined edge is ineligible.
        "((a:0.1,b:0.2):-0.05,(c:0.3,d:0.7,e:0):0.9);",
        # Entirely undated side: its eligible edges have a flat OLS objective.
        "(((a:0.1,b:0.2):0.3,c:0.7):0.1,(d:1,e:0.1):0.9);",
    ],
)
@pytest.mark.parametrize("dated_names", ["abcde", "abc"])
def test_stable_dp_matches_independent_continuous_ols_and_signed_distances(newick, dated_names):
    original = tree(newick)
    dates = {name: f"{2000 + 3 * index}-02-11" for index, name in enumerate(dated_names)}
    expected = _continuous_reference(original, dates)
    before = {(a, b): original.distance(a, b) for a in "abcde" for b in "abcde"}

    rooted, audit = optimise_clock_root(original, samples(dates))

    assert rooted is not None
    assert audit["mse"] == pytest.approx(expected, abs=1e-10)
    assert _mse(dates, {name: rooted.distance(name) for name in dates}) == pytest.approx(
        audit["mse"], abs=1e-10
    )
    for pair, distance in before.items():
        assert rooted.distance(*pair) == pytest.approx(distance, abs=1e-12)
    assert audit["dated_tip_count"] == len(dated_names)
    if len(dated_names) < 5:
        assert {
            item["sample_id"]
            for item in audit["excluded_dates"]
            if item["reason"] == "missing_metadata"
        } == {"d", "e"}


def test_dp_messages_match_explicit_distance_moments_at_every_node():
    original = tree("(((a:0.1,b:-0.2):0.3,c:0):0.7,d:0.9,e:0.4);")
    dates = {"a": -2.0, "c": 0.5, "e": 1.5}
    nodes = clock_rooting._preorder(original.root)
    subtree, totals, outside = clock_rooting._distance_moments(nodes, dates)
    for node in nodes:
        inside_names = {tip.name for tip in node.get_terminals()}
        for moment, names in [
            (subtree[node], set(dates) & inside_names),
            (outside[node], set(dates) - inside_names),
            (totals[node], set(dates)),
        ]:
            n, x, y, xx, yy, xy = moment
            assert n == len(names)
            if not names:
                continue
            date_values = np.asarray([dates[name] for name in sorted(names)])
            distances = np.asarray([original.distance(node, name) for name in sorted(names)])
            assert x == pytest.approx(date_values.mean())
            assert y == pytest.approx(distances.mean())
            assert xx == pytest.approx(np.sum((date_values - x) ** 2), abs=1e-12)
            assert yy == pytest.approx(np.sum((distances - y) ** 2), abs=1e-12)
            assert xy == pytest.approx(np.dot(date_values - x, distances - y), abs=1e-12)


def test_rooting_ignores_root_stem_and_keeps_independent_clade_metadata():
    original = tree("(a:0.2,b:0.5,(c:1,d:0.7):0.8,e:0.3):999;")
    original.root.confidence = 88
    original.find_any(name="c").comment = "example"
    original.find_any(name="c").annotations = {"country": ["GB"]}
    dates = {"a": "2001", "b": "2007", "c": "2012", "d": "2017", "e": "2022"}
    no_stem = tree("(a:0.2,b:0.5,(c:1,d:0.7):0.8,e:0.3);")
    rooted, audit = optimise_clock_root(original, samples(dates))
    _, stemless_audit = optimise_clock_root(no_stem, samples(dates))
    assert rooted is not None
    assert audit["mse"] == stemless_audit["mse"]
    assert audit["root_position_from_parent"] == stemless_audit["root_position_from_parent"]
    tip = rooted.find_any(name="c")
    assert tip.comment == "example"
    tip.annotations["country"].append("FR")
    assert original.find_any(name="c").annotations == {"country": ["GB"]}


def test_compact_ids_and_zero_edge_ties_are_repeatable():
    original = tree("(a:0,b:0,(c:0,d:0):0);")
    dates = {"a": "2001", "b": "2007", "c": "2012", "d": "2017"}
    first, one = optimise_clock_root(original, samples(dates))
    second, two = optimise_clock_root(original, list(reversed(samples(dates))))
    assert first is not None and second is not None
    assert one["optimal_edge_id"] == two["optimal_edge_id"] == "node:0->a"
    assert one["mse"] == two["mse"] == 0
    assert one["root_position_from_parent"] == 0
    assert all("{" not in item["edge_id"] for item in one["excluded_negative_edges"])


def _comb_tree(size):
    from Bio.Phylo.BaseTree import Clade, Tree

    root = Clade(
        clades=[Clade(name="tip0", branch_length=0.1), Clade(name="tip1", branch_length=0.2)]
    )
    for index in range(2, size):
        root.branch_length = 0.01
        root = Clade(clades=[root, Clade(name=f"tip{index}", branch_length=0.1)])
    return Tree(root=root)


def test_full_7815_tip_deep_tree_uses_linear_work_without_recursive_or_pairwise_apis(monkeypatch):
    from Bio.Phylo.BaseTree import TreeMixin

    original = _comb_tree(7815)
    metadata = [
        {"sample_id": f"tip{i}", "collection_date": str(2000 + i % 25)} for i in range(7815)
    ]
    merges = 0
    merge = clock_rooting._merge_moments

    def counted(left, right):
        nonlocal merges
        merges += 1
        return merge(left, right)

    def forbidden(*args, **kwargs):
        raise AssertionError("quadratic/recursive tree API called")

    monkeypatch.setattr(clock_rooting, "_merge_moments", counted)
    monkeypatch.setattr(TreeMixin, "get_terminals", forbidden)
    monkeypatch.setattr(TreeMixin, "find_clades", forbidden)
    monkeypatch.setattr(TreeMixin, "distance", forbidden)
    rooted, audit = optimise_clock_root(original, metadata)

    assert rooted is not None
    node_count = 2 * 7815 - 1
    assert merges <= 7 * node_count
    assert audit["candidate_edge_count"] == node_count - 2
    assert audit["dated_tip_count"] == 7815
    assert audit["pairwise_distances_preserved"] is True
    assert len([node for node in clock_rooting._preorder(rooted.root) if not node.clades]) == 7815


def test_indexed_tree_returns_fresh_distance_index_after_rooting():
    from chronoclade.cgmlst.scale_tree import IndexedTree

    original = IndexedTree.from_tree(tree("(a:0.2,b:0.5,(c:1,d:0.7):0.8,e:0.3);"))
    dates = {"a": "2001", "b": "2007", "c": "2012", "d": "2017", "e": "2022"}
    rooted, audit = optimise_clock_root(original, samples(dates))
    assert isinstance(rooted, IndexedTree)
    assert rooted.root is rooted._nodes[0]
    assert _mse(dates, {name: rooted.distance(name) for name in dates}) == pytest.approx(
        audit["mse"]
    )
    for a in dates:
        for b in dates:
            assert rooted.distance(a, b) == pytest.approx(original.distance(a, b))


def test_linear_topology_proof_rejects_mutated_edge(monkeypatch):
    original = tree("(a:0.2,b:0.5,(c:1,d:0.7):0.8,e:0.3);")
    dates = {"a": "2001", "b": "2007", "c": "2012", "d": "2017", "e": "2022"}
    root_copy = clock_rooting._root_copy

    def corrupt(*args):
        result, clones = root_copy(*args)
        clones[original.find_any(name="c")].branch_length += 0.01
        return result, clones

    monkeypatch.setattr(clock_rooting, "_root_copy", corrupt)
    rooted, audit = optimise_clock_root(original, samples(dates))
    assert rooted is None
    assert audit["reason"] == "pairwise_distance_check_failed"


@pytest.mark.parametrize("seed", range(12))
def test_random_signed_polytomies_match_continuous_reference(seed):
    from Bio.Phylo.BaseTree import Clade, Tree

    rng = np.random.default_rng(seed)
    groups = [Clade(name=f"tip{i}", branch_length=float(rng.uniform(-0.2, 1.2))) for i in range(9)]
    while len(groups) > 3:
        arity = min(int(rng.integers(2, 4)), len(groups) - 1)
        chosen = [groups.pop(int(rng.integers(len(groups)))) for _ in range(arity)]
        groups.append(Clade(branch_length=float(rng.uniform(-0.2, 1.2)), clades=chosen))
    original = Tree(root=Clade(clades=groups))
    dates = {f"tip{i}": f"{2000 + 2 * i}-03-15" for i in range(7)}
    rooted, audit = optimise_clock_root(original, samples(dates))
    assert rooted is not None
    assert audit["mse"] == pytest.approx(_continuous_reference(original, dates), abs=1e-10)


def test_all_negative_edges_skip_without_altering_source():
    original = tree("(a:-1,b:-2,c:-3);")
    rooted, audit = optimise_clock_root(original, samples({"a": "2000", "b": "2001", "c": "2002"}))
    assert rooted is None
    assert audit["reason"] == "no_nonnegative_edges"
    assert audit["candidate_edge_count"] == 0
    assert [tip.branch_length for tip in original.get_terminals()] == [-1, -2, -3]


def test_extra_terminal_introduced_by_bad_reroot_is_rejected(monkeypatch):
    from Bio.Phylo.BaseTree import Clade

    original = tree("(a:1,b:2,c:3);")
    root_copy = clock_rooting._root_copy

    def corrupt(*args):
        result, clones = root_copy(*args)
        result.root.clades.append(Clade(name="extra", branch_length=0))
        return result, clones

    monkeypatch.setattr(clock_rooting, "_root_copy", corrupt)
    rooted, audit = optimise_clock_root(original, samples({"a": "2000", "b": "2001", "c": "2002"}))
    assert rooted is None
    assert audit["reason"] == "pairwise_distance_check_failed"
