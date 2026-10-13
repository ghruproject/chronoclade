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
    years = np.asarray([sample_date_interval({"collection_date": date})["midpoint_year"]
                        for date in dates.values()])
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
                    distances[sample_id] = base + (-position if sample_id in descendants else position)
                best = min(best, _mse(dates, distances))
    return best


def test_analytic_optimum_matches_dense_brute_force_grid():
    original = tree("(a:0.2,b:0.5,(c:1,d:0.7):0.8,e:0.3);")
    dates = {"a": "2001", "b": "2007-05", "c": "2012-03-04",
             "d": "2017", "e": "2022-09-15"}
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
    before = {(a, b): original.distance(a, b)
              for a in "abcde" for b in "abcde"}
    dates = {"a": "2001", "b": "2007", "c": "2012", "d": "2017", "e": "2022"}

    rooted, audit = optimise_clock_root(original, samples(dates))

    assert rooted is not original
    assert audit["pairwise_distances_preserved"] is True
    for pair, distance in before.items():
        assert rooted.distance(*pair) == pytest.approx(distance)
        assert original.distance(*pair) == pytest.approx(distance)


def test_bifurcating_root_artificial_edge_is_searched_and_preserved():
    original = tree("((a:1,b:1):2,(c:1,d:1):2);")
    before = {(a, b): original.distance(a, b)
              for a in "abcd" for b in "abcd"}
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
        "invalid_collection_date", "sample_not_in_tree", "missing_metadata"
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
