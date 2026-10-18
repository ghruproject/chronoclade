"""Plot all frozen candidate results; no ranking, fitted threshold or default change."""

import argparse
import csv
import json
from collections import Counter
from html import escape
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from chronoclade.report_components.styles import report_styles
from chronoclade.report_components.tables import write_report_html


def rows(path):
    return list(csv.DictReader(path.open()))


def table(headers, records):
    return (
        '<div class="table-scroll"><table><thead><tr>'
        + "".join("<th>" + escape(h) + "</th>" for h in headers)
        + "</tr></thead><tbody>"
        + "".join(
            "<tr>" + "".join("<td>" + escape(str(v)) + "</td>" for v in row) + "</tr>"
            for row in records
        )
        + "</tbody></table></div>"
    )


def render(directory, labels_path):
    directory = Path(directory)
    summary = json.loads((directory / "summary.json").read_text())
    query = rows(directory / "queries.csv")
    pairs = rows(directory / "pairs.csv")
    counts = rows(directory / "candidate_counts.csv")
    labels = {row["sample_id"]: row for row in rows(Path(labels_path))}
    methods = [row["candidate"] for row in counts] + ["nearest_all_ties"]
    metrics = []
    for method in methods:
        selected = [
            q
            for q in query
            if q["candidate"] == method and labels[q["sample_id"]]["epi_status"] == "outbreak"
        ]
        recoveries = [float(q["recovery"]) for q in selected if q["recovery"]]
        metrics.append(
            dict(
                method=method,
                recovery=float(np.mean(recoveries)) if recoveries else None,
                controls=float(np.mean([int(q["negative_retrieved"]) for q in selected])),
                max_controls=max(int(q["negative_retrieved"]) for q in selected),
                size=float(np.mean([int(q["burden"]) for q in selected])),
                abstaining=sum(int(q["abstaining_pairs"]) for q in selected),
                nearest=sum(int(q["nearest_ties_retained"]) for q in selected),
                nearest_total=sum(int(q["nearest_ties_available"]) for q in selected),
            )
        )
    names = [
        m["method"].replace("alleles_le_", "≤").replace("nearest_all_ties", "Nearest ties")
        for m in metrics
    ]
    colours = [
        "#3c7383"
        if m["method"].startswith("alleles")
        else "#754b8e"
        if m["method"].startswith("LIN")
        else "#59616c"
        for m in metrics
    ]
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, layout="constrained")
    axes[0].bar(
        names,
        [100 * m["recovery"] if m["recovery"] is not None else 0 for m in metrics],
        color=colours,
    )
    axes[0].set_ylabel("Known event peers recovered (%)")
    axes[0].set_ylim(0, 105)
    axes[1].bar(names, [m["controls"] for m in metrics], color=colours)
    axes[1].set_ylabel("Sporadic controls admitted\n(mean per outbreak query)")
    axes[1].set_ylim(0, 5.2)
    axes[1].tick_params(axis="x", rotation=45)
    for ax in axes:
        ax.set_axisbelow(True)
        ax.grid(axis="y", alpha=0.15)
    fig.suptitle("One published ST11 event: recovery and review burden", fontsize=13)
    fig.savefig(directory / "candidate-comparison.svg")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 4), layout="constrained")
    sets = []
    for label, title, colour in [
        ("positive", "Same documented event", "#754b8e"),
        ("negative", "Event versus sporadic control", "#3c7383"),
        ("unknown", "Control versus control (unknown)", "#7f8792"),
    ]:
        values = [
            int(p["allele_mismatches"])
            for p in pairs
            if p["epi_label"] == label and p["eligible"] == "True"
        ]
        sets.append((label, values))
        if values:
            ax.hist(
                values,
                bins=np.arange(
                    -0.5,
                    max(int(p["allele_mismatches"]) for p in pairs if p["allele_mismatches"]) + 1.5,
                ),
                weights=np.ones(len(values)) / len(values),
                histtype="step",
                linewidth=1.8,
                label=f"{title} (n={len(values)} pairs)",
                color=colour,
            )
    ax.set_xlabel("Observed cgMLST allele differences")
    ax.set_ylabel("Fraction of pairs in each label")
    ax.legend(frameon=False, fontsize=9)
    fig.savefig(directory / "distance-distributions.svg")
    plt.close(fig)
    coverage = [int(p["shared_loci"]) for p in pairs]
    pair_summary = [
        [
            {
                "positive": "Same event",
                "negative": "Event versus sporadic",
                "unknown": "Sporadic versus sporadic (unknown)",
            }[label],
            len(values),
            min(values) if values else "Unavailable",
            float(np.median(values)) if values else "Unavailable",
            max(values) if values else "Unavailable",
        ]
        for label, values in sets
    ]
    metric_rows = [
        [
            m["method"].replace("alleles_le_", "≤ ").replace("nearest_all_ties", "Nearest ties"),
            f"{100 * m['recovery']:.1f}%" if m["recovery"] is not None else "Unavailable",
            f"{m['controls']:.2f}",
            m["max_controls"],
            f"{m['size']:.2f}",
            f"{m['nearest']}/{m['nearest_total']}",
            m["abstaining"],
        ]
        for m in metrics
    ]
    body = (
        '<section class="stage"><div class="stage-body"><h2>Dataset and limits</h2><p>31 ST11 assemblies: 26 author-labelled outbreak isolates and five explicitly sporadic controls. This is one development event, not held-out validation. Labels used epidemiology plus PFGE/MLST; they are not patient-level transmission truth.</p><p>All genomes were retyped with one complete 629-locus catalogue. Shared calls across pairs: '
        + str(min(coverage))
        + "–"
        + str(max(coverage))
        + ". Primary minimum: "
        + str(summary["coverage_policy"]["min_shared_loci"])
        + ". Missing calls are excluded, not treated as agreement.</p><p>Candidate cutoffs and LIN levels were fixed before typing. Every result is shown; no winner or new default is selected.</p></div></section>"
    )
    body += (
        '<section class="stage"><div class="stage-body"><h2>Recovery versus control admission</h2><p>Each outbreak query has 25 other documented event peers and five controls. Bars average the 26 outbreak queries. Pairs share isolates; these are descriptive counts, with no independent-pair confidence intervals. LIN abstentions are shown in the table.</p>'
        + (directory / "candidate-comparison.svg").read_text()
        + table(
            [
                "Rule",
                "Event-peer recovery",
                "Mean controls admitted",
                "Maximum controls admitted",
                "Mean neighbourhood size",
                "Nearest ties retained",
                "LIN pair abstentions",
            ],
            metric_rows,
        )
        + "</div></section>"
    )
    body += (
        '<section class="stage"><div class="stage-body"><h2>Observed distance distributions</h2><p>Different labels can overlap genetically. A low distance does not establish direct transmission; a shared LIN prefix is not a maximum pairwise-distance guarantee.</p>'
        + (directory / "distance-distributions.svg").read_text()
        + table(
            [
                "Comparison",
                "Eligible pairs",
                "Minimum alleles",
                "Median alleles",
                "Maximum alleles",
            ],
            pair_summary,
        )
        + "</div></section>"
    )
    time_counts = Counter((p["epi_label"], p["time_bin"]) for p in pairs if p["eligible"] == "True")
    body += (
        '<section class="stage"><div class="stage-body"><h2>Collection-time coverage</h2><p>Dates are reported observations. These bins describe available comparisons; they do not fit a time-dependent cutoff.</p>'
        + table(
            ["Pair label", "Date separation (days)", "Pairs"],
            [[a, b, n] for (a, b), n in sorted(time_counts.items())],
        )
        + "</div></section>"
    )
    snp_file = directory / "snp-comparison.json"
    if snp_file.exists():
        snp = json.loads(snp_file.read_text())
        descriptions = {
            "positive": "Same event",
            "negative": "Event versus sporadic",
            "unknown": "Sporadic versus sporadic (unknown)",
        }
        snp_rows = [
            [descriptions[r["label"]], r["pairs"], r["snp_min"], r["snp_median"], r["snp_max"]]
            for r in snp["label_distributions"]
        ]
        zero = snp["zero_allele_pairs"]
        body += (
            '<section class="stage"><div class="stage-body"><h2>Matched SNP comparison</h2>'
            "<p>The original strict run stopped after flagging the five known controls. "
            "This separate diagnostic run retains all 31 genomes and all flags. "
            "It is exploratory; temporal assessment and dating were disabled.</p>"
            + (directory / "allele-snp-comparison.svg").read_text()
            + table(
                ["Relationship", "Pairs", "Minimum SNPs", "Median SNPs", "Maximum SNPs"], snp_rows
            )
            + "<p>Pair-specific callable sites: "
            + f"{snp['callable_sites_min']:,}–{snp['callable_sites_max']:,}. "
            + "The corrected tree contains an edge separating all 26 event isolates from the five controls. "
            "No bootstrap support was estimated.</p>"
            + f"<p>Among {zero['pairs']} zero-allele-difference pairs, corrected SNP distances "
            + f"range from {zero['snp_min']} to {zero['snp_max']} (median {zero['snp_median']:g}). "
            + "SNPs measure variation across the filtered core alignment; allele types cover a different, smaller set of genes.</p>"
            + '<p><a href="'
            + escape(snp["tree_report_relative_path"], quote=True)
            + '">Open the diagnostic SNP-tree report</a> · '
            + '<a href="allele-snp-pairs.csv">Matched allele/SNP pairs</a> · '
            + '<a href="snp-comparison.json">SNP comparison evidence</a></p>'
            + "</div></section>"
        )
    body += '<section class="stage"><div class="stage-body"><h2>Evidence and next step</h2><p><a href="pairs.csv">All pair observations</a> · <a href="queries.csv">Every query neighbourhood</a> · <a href="samples.csv">Typing coverage</a> · <a href="summary.json">Reference and input hashes</a></p><p>Next evaluate the unchanged candidate rules on an independent event. Neither this single study nor its source paper validates a universal 629-locus cutoff.</p></div></section>'
    html = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Klebsiella relatedness pilot</title><style>'
        + report_styles()
        + 'svg{max-width:100%;height:auto}</style></head><body><main class="shell"><header class="identity"><div class="identity-mark">ChronoClade · RESEARCH PILOT</div><div class="identity-copy"><h1>Klebsiella ST11: grouping comparison</h1><p>Zhou2017 · 31 ST11 genomes · one development study</p></div></header>'
        + body
        + "</main></body></html>"
    )
    write_report_html(directory / "report.html", html, encoding="utf-8")
    (directory / "plot-summary.json").write_text(
        json.dumps({"metrics": metrics, "pair_summaries": pair_summary}, indent=2) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    args = parser.parse_args()
    render(args.results, args.labels)
