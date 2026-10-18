"""Readable evidence tables for independent dating comparisons."""

from html import escape

from chronoclade.report_components.stage import write_report
from chronoclade.sample_labels import display_label


def _display(value):
    if value is None:
        return "Unavailable"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, float):
        return format(value, ".7g")
    return str(value)


def _table(columns, items):
    heads = "".join("<th>" + escape(label) + "</th>" for _, label in columns)
    body = "".join(
        "<tr>"
        + "".join("<td>" + escape(_display(item.get(key))) + "</td>" for key, _ in columns)
        + "</tr>"
        for item in items
    )
    return (
        '<div class="table-scroll"><table><thead><tr>'
        + heads
        + "</tr></thead><tbody>"
        + body
        + "</tbody></table></div>"
    )


def _artifact(spec, label):
    if not spec:
        return "Unavailable"
    return (
        '<a href="'
        + escape(spec["path"], quote=True)
        + '">'
        + escape(label)
        + '</a> <code title="'
        + escape(spec["sha256"], quote=True)
        + '">'
        + spec["sha256"][:12]
        + "…</code>"
    )


def write_comparison_report(output, result):
    runs, coverage = result["runs"], result["selection_coverage"]
    labels = result.get("sample_labels", {})
    anchor_count = len(result.get("common_anchor_ids", []))
    input_ids = result.get("input_target_ids", result.get("target_ids", []))
    same_targets = set(input_ids) == set(result.get("common_anchor_ids", []))
    extra = (
        f"<p>The shared-anchor ancestor uses {anchor_count} genomes present in every "
        "selection, including inputs and any shared required public comparisons. "
        f"The declared input target contains {len(input_ids)} genomes."
        + (" These are the same target.</p>" if same_targets else "</p>")
        + "<h3>Independent runs</h3>"
    ) + _table(
        [
            ("selection_id", "Selection"),
            ("status", "Dating status"),
            ("rate", "Rate (substitutions/site/year)"),
            ("common_anchor_mrca_date", "Shared-anchor MRCA year"),
            ("root_date", "Whole-selection root year"),
            ("dating_override", "Override"),
            ("reason", "Assessment / failure"),
        ],
        runs,
    )
    extra += "<h3>Within-run dating uncertainty</h3><p>Saved TreeTime 90% max-posterior regions belong to each run's actual shared-anchor ancestor. Its descendants can include optional context genomes. Missing bounds remain unavailable.</p>"
    if result.get("shared_anchor_date_forest"):
        extra += (
            '<img style="width:100%;height:auto" src="'
            + escape(result["shared_anchor_date_forest"]["path"], quote=True)
            + '" alt="Within-run MRCA intervals and separately labelled selection sensitivity"/>'
        )
    intervals = []
    for run in runs:
        row = run.get("common_anchor_mrca") or {}
        intervals.append(
            {**row, "selection_id": run["selection_id"], "descendants": row.get("descendant_count")}
        )
    extra += _table(
        [
            ("selection_id", "Selection"),
            ("node", "Saved ancestor"),
            ("descendants", "Descendant genomes"),
            ("estimate", "MRCA year"),
            ("lower_bound", "90% region lower"),
            ("upper_bound", "90% region upper"),
        ],
        intervals,
    )
    input_intervals = []
    if not same_targets:
        for run in runs:
            row = run.get("input_target_mrca") or {}
            if run.get("status") == "dated" and row.get("estimate") is not None:
                input_intervals.append({
                    **row, "selection_id": run["selection_id"],
                    "extra_descendants": len(row.get("extra_descendant_sample_ids", [])),
                })
    if input_intervals:
        extra += (
            f"<h3>Input-genome ancestor ({len(input_ids)} genomes)</h3>"
            "<p>Each estimate is the ancestor of the same declared input genomes. "
            "Additional descendants can change between selections. Bounds are the saved "
            "within-run TreeTime 90% max-posterior region; missing bounds remain unavailable.</p>"
            + _table([
                ("selection_id", "Selection"), ("estimate", "Ancestor year"),
                ("lower_bound", "90% region lower"), ("upper_bound", "90% region upper"),
                ("descendant_count", "Full descendants"),
                ("extra_descendants", "Additional descendants"),
            ], input_intervals)
        )
    extra += (
        "<h3>Between-selection sensitivity</h3><p>The minimum, median, maximum and range of point estimates describe selection sensitivity; they are not confidence intervals.</p>"
        + _table(
            [
                ("metric", "Metric"),
                ("assessable_runs", "Assessable runs"),
                ("minimum", "Minimum"),
                ("median", "Median"),
                ("maximum", "Maximum"),
                ("range", "Range"),
            ],
            [
                {"metric": "Input-genome ancestor year" if field == "input_target_mrca_date" else field,
                 **spread}
                for field, spread in result["sensitivity"].items()
                if field != "input_target_mrca_date" or input_intervals
            ],
        )
    )
    extra += "<h3>Dating and alignment audit</h3>" + _table(
        [
            ("selection_id", "Selection"),
            ("selected_sample_count", "Selected"),
            ("dated_sample_count", "Dated"),
            ("undated_sample_count", "Undated"),
            ("interval", "Metadata date interval"),
            ("temporal_status", "Temporal gate"),
            ("reference_label", "Reference genome"),
            ("alignment_columns", "Alignment columns"),
            ("complete_alignment_sites", "Complete alignment sites"),
        ],
        [
            {
                **run,
                "interval": " to ".join(str(v) for v in run.get("date_interval", {}).values())
                if run.get("date_interval")
                else None,
            }
            for run in runs
        ],
    )
    extra += "<h4>Independent recombination filtering</h4>" + _table(
        [
            ("selection_id", "Selection"),
            ("recombination_removed_sites", "Imported columns removed"),
            ("incomplete_only_removed_sites", "Additional incomplete columns removed"),
            ("retained_clonal_percent", "Retained clonal columns (%)"),
        ],
        runs,
    )
    for run in runs:
        if "alignment" not in run:
            continue
        extra += (
            "<details><summary>"
            + escape(run["selection_id"])
            + " reference, alignment and recombination mask checksums</summary><p>"
            + _artifact(run["reference"], "Reference FASTA")
            + " · "
            + _artifact(run["alignment"], "Saved alignment")
            + "</p>"
        )
        extra += (
            "<ul>"
            + "".join(
                "<li>" + _artifact(spec, name) + "</li>"
                for name, spec in run["recombination_mask"]["artifacts"].items()
            )
            + "</ul><p>Full filtering counts, parameters and checksums are saved in comparison.json.</p></details>"
        )
    extra += (
        "<h3>Full available block / pool coverage</h3><p>"
        + escape(coverage["interpretation"])
        + "</p>"
        + _table(
            [
                ("pool_sample_count", "Available genomes"),
                ("input_count", "Inputs"),
                ("context_pool_count", "Available context"),
                ("requested_selection_union_count", "Unique genomes across requested selections"),
                ("requested_selection_union_fraction", "Requested-selection union fraction"),
                ("available_tree_union_count", "Unique genomes in available trees"),
                ("available_tree_union_fraction", "Available-tree union fraction"),
            ],
            [coverage],
        )
    )
    extra += _table(
        [("selection_id", "Selection"), ("cumulative_unique_genomes", "Cumulative unique genomes")],
        coverage["cumulative_unique_genomes"],
    )
    overlap = [
        {**row, "selections": row["left"] + " / " + row["right"]}
        for row in coverage["pairwise_overlap"]
    ]
    extra += "<h4>Selection overlap</h4>" + _table(
        [
            ("selections", "Selections"),
            ("intersection", "Shared genomes"),
            ("jaccard", "Jaccard overlap"),
        ],
        overlap,
    )
    for row in coverage["runs"]:
        representative = row["representative_coverage"]
        extra += (
            "<details><summary>"
            + escape(row["selection_id"])
            + " representative-distance and strata coverage</summary><p>"
            + _artifact(representative["artifact"], "Coverage CSV")
            + " · "
            + escape(representative["ranking_metric"])
            + ": "
            + escape(representative["distance_definition"])
            + "</p>"
        )
        extra += _table(
            [
                ("sample_count", "Pool genomes"),
                ("comparable_sample_count", "Comparable"),
                ("noncomparable_sample_count", "Non-comparable"),
                ("minimum", "Minimum nearest value"),
                ("q25", "Q25"),
                ("median", "Median"),
                ("q75", "Q75"),
                ("maximum", "Maximum"),
            ],
            [{**representative, **representative["nearest_value_quantiles"]}],
        )
        extra += (
            "<p>Genetic, date, location and host strata cover context genomes. All inputs are retained separately.</p>"
            + _table(
                [
                    ("genetic", "Genetic group"),
                    ("date", "Date stratum"),
                    ("location", "Location"),
                    ("host", "Host"),
                    ("available", "Available context"),
                    ("selected", "Selected context"),
                ],
                [
                    {
                        **group,
                        **dict(zip(("genetic", "date", "location", "host"), group["stratum"])),
                    }
                    for group in row["strata_coverage"]
                ],
            )
            + "</details>"
        )
    extra += "<h3>Shared-target clustering</h3>"
    cluster_policy = result.get("cluster_policy", {})
    cutoff = cluster_policy.get("snp_cutoff")
    method = cluster_policy.get("method")
    if method == "single-linkage corrected SNP distances" and cutoff is not None:
        extra += (
            "<p>Co-clustering uses connected components from single-linkage on corrected "
            f"clonal-SNP distances, with a configured cutoff of {escape(_display(cutoff))} SNPs "
            "for pairs with callable sites. A pair can connect through intermediate genomes, "
            "so not every reported pair is directly within the cutoff. This descriptive threshold "
            "is not a validated biological or outbreak cutoff; the pairwise results are distinct "
            "from published cgLIN groups.</p>"
        )
    extra += _table(
        [
            ("target_pair", "Input genome pair"),
            ("assessable_runs", "Assessable runs"),
            ("co_clustered_runs", "Co-clustered runs"),
        ],
        [
            {
                **row,
                "target_pair": " · ".join(
                    display_label(labels.get(ident, ident))
                    for ident in row["target_pair"].split("|")
                ),
            }
            for row in result["target_pair_agreement"]
        ],
    )
    splits = result["shared_anchor_splits"]
    extra += (
        "<h3>Descriptive unrooted split occurrence</h3><p>"
        + escape(splits["interpretation"])
        + " Statistical branch support is unavailable.</p>"
        + _table(
            [
                ("split_id", "Split identifier"),
                ("observed_runs", "Observed runs"),
                ("assessable_runs", "Available trees"),
                ("frequency", "Descriptive frequency"),
            ],
            [{**row, "split_id": row["split_id"][:12]} for row in splits["splits"]],
        )
    )
    extra += "<p>" + escape(result["interpretation"]) + "</p>"
    write_report(
        output / "report.html",
        title="Selection sensitivity",
        report_label="TIME TREE COMPARISON REPORT",
        overview=f"{len(result['common_anchor_ids'])} shared anchor genomes · {len(runs)} requested selections",
        result={
            "selected_sample_ids": result["common_anchor_ids"],
            "temporal_assessment": {
                "code": "sensitivity",
                "supported": False,
                "reason": f"{result['dated_runs']} of {len(runs)} requested runs produced dates",
            },
        },
        extra=extra,
    )
