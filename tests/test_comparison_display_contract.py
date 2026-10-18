from chronoclade.cgmlst.clock_comparison import _figure
from chronoclade.time_stage import comparison_report


def test_matched_clock_svg_uses_short_allele_display_label(tmp_path):
    fit = {
        "slope": 0.1, "centered_intercept": 2.0, "reference_year": 2000.0,
        "r_squared": 0.8, "n_points": 3,
        "distance_units": "scheme-equivalent allele-distance units on 629 declared loci",
        "display_distance_units": "Estimated allele differences",
        "points": [{"midpoint_year": 1999.0, "root_to_tip": 1.0},
                   {"midpoint_year": 2001.0, "root_to_tip": 3.0}],
    }
    panel = {"nonnegative_original_root": fit, "clock_compatible": fit}
    path = tmp_path / "comparison.svg"
    _figure({"cgmlst": panel, "corrected_snp": panel}, path)
    svg = path.read_text()
    assert "Estimated allele differences" in svg
    assert "scheme-equivalent allele-distance units on 629 declared loci" not in svg


def test_shared_target_caption_explains_single_linkage_cutoff(monkeypatch, tmp_path):
    captured = {}
    monkeypatch.setattr(comparison_report, "write_report", lambda _path, **kwargs: captured.update(kwargs))
    coverage = {
        "interpretation": "coverage", "cumulative_unique_genomes": [],
        "pairwise_overlap": [], "runs": [],
    }
    result = {
        "runs": [], "selection_coverage": coverage, "common_anchor_ids": [],
        "input_target_ids": [], "target_ids": [], "sensitivity": {},
        "target_pair_agreement": [], "shared_anchor_splits": {"interpretation": "splits", "splits": []},
        "interpretation": "comparison", "dated_runs": 0, "sample_labels": {},
        "cluster_policy": {"method": "single-linkage corrected SNP distances", "snp_cutoff": 10},
    }
    comparison_report.write_comparison_report(tmp_path, result)
    html = captured["extra"]
    assert "configured cutoff of 10 SNPs" in html
    assert "intermediate genomes" in html
    assert "not every reported pair is directly within the cutoff" in html
    assert "not a validated biological or outbreak cutoff" in html
    assert "distinct from published cgLIN groups" in html

    assert html.index("<h3>Shared-target clustering</h3>") < html.index("Co-clustering uses")
