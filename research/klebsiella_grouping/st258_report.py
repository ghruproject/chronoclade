"""Plain-language, offline report for the ST258 distance and sampling comparison."""

import argparse
from collections import Counter, defaultdict
import csv
from html import escape
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from chronoclade.report_components.styles import report_styles
from chronoclade.report_components.tables import write_report_html
from research.klebsiella_grouping.report import table

MAIN = ("Allele profiles", "Aligned locus SNPs")
COLORS = ("#71527f", "#36828c")


def read_csv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def section(title, content, ident):
    return f'<section class="stage" id="{ident}"><div class="stage-body"><h2>{escape(title)}</h2>{content}</div></section>'


def fmt(value):
    return f"{value:,}" if isinstance(value, int) else str(value)


def short_labels(labels):
    return ", ".join(labels[:4]) + (f" +{len(labels) - 4} more" if len(labels) > 4 else "")


def render(directory):
    directory = Path(directory)
    data = json.loads((directory / "dataset-summary.json").read_text())
    nearest = json.loads((directory / "nearest-comparison.json").read_text())
    groups = json.loads((directory / "groups.json").read_text())
    thresholds = read_csv(directory / "threshold-summary.csv")
    samples = read_csv(directory / "sampling.csv")
    catalogue = json.loads((directory.parent / "catalogue/summary.json").read_text())
    shared = sum(row["allele_snp_nearest_overlap"] > 0 for row in nearest)
    exact = sum(row["same_nearest_set"] for row in nearest)
    matched = sum(row["matched_allele_snp_same_set"] for row in nearest)
    total = data["total_pairs"]
    eligible = data["eligible_pairs"]
    chain_example = next(r for r in thresholds if r["method"] == MAIN[0] and int(r["cutoff"]) == 2)
    body = section(
        "Start here",
        f"<p>This compares <strong>{data['inputs']} Greek input genomes</strong> with "
        f"<strong>{data['contexts']:,} public comparison genomes</strong>, all from the cached ST258 dataset. "
        "No new assemblies were downloaded.</p>"
        "<p><strong>Allele distance</strong> counts genes with different allele calls. "
        "<strong>Locus SNP distance</strong> counts base changes after aligning the available DNA sequences "
        "separately for each gene, then adds those counts across the cgMLST genes. "
        "Indels and missing calls are recorded separately.</p>"
        f"<p><strong>Closest public relatives:</strong> {shared}/{data['inputs']} inputs share at least one "
        f"closest comparison between the methods; {exact}/{data['inputs']} have exactly the same set of closest comparisons. "
        "Ties are retained, so a closest set may contain many genomes. "
        f"Using only loci with DNA available to both methods raises exact agreement to {matched}/{data['inputs']}; "
        "missing DNA explains part of the difference.</p>"
        "<p><strong>Genetic neighbourhoods:</strong> a cutoff joins genomes through chains of short distances. "
        "Their endpoints can be farther apart than that cutoff. These are exploratory genetic groups; "
        "they do not establish transmission or an outbreak.</p>"
        f"<p>For example, at two allele differences, {int(chain_example['direct_contexts']):,} public genomes "
        f"are directly within that distance of an input, but connecting paths bring in {int(chain_example['linked_contexts']):,}. "
        "A small joining cutoff can therefore produce a very broad group.</p>"
        "<p>Read downward: dataset → closest relatives → neighbourhood size → sampling sensitivity → collection time.</p>",
        "start",
    )
    country_rows = [
        [r["country"], r["inputs"], r["contexts"], r["inputs"] + r["contexts"]]
        for r in data["countries"]
    ]
    unavailable = catalogue["unresolved_dna_cells"] / catalogue["called_cells"] * 100
    body += section(
        "What data went into the comparison?",
        table(["Country", "Your inputs", "Public comparisons", "Total"], country_rows)
        + "<p>We compare the <strong>628 cgMLST genes represented in the cached export</strong>. Unreported loci remain missing. "
        f"DNA could not be recovered for {catalogue['unresolved_dna_cells']:,} individual allele calls ({unavailable:.3f}% of calls). "
        "Those calls stay missing; unavailable sequences are concentrated in particular alleles rather than randomly distributed.</p>"
        + f"<p>Both methods use the same {eligible:,}/{total:,} eligible genome pairs, requiring at least "
        "597 jointly available loci (95% of this observed panel) in both representations. "
        "An unavailable comparison is not a zero distance. "
        + (
            f"All {data['inputs']} inputs have enough callable loci.</p>"
            if not data["unusable_inputs"]
            else f"{len(data['unusable_inputs'])} inputs have insufficient callable loci.</p>"
        )
        + (
            "<details><summary>Reported host metadata</summary>"
            + table(
                ["Reported host", "Genomes"],
                [
                    [host, n]
                    for host, n in sorted(
                        data["host_counts"].items(), key=lambda pair: (-pair[1], pair[0])
                    )
                ],
            )
            + "<p>Original host labels are retained apart from case and whitespace. Missing host data remain Unknown.</p></details>"
            if data["host_reported"]
            else ""
        ),
        "dataset",
    )
    nearest_rows = []
    for row in nearest:
        a, b = row[MAIN[0]], row[MAIN[1]]
        nearest_rows.append(
            [
                row["input"],
                a["distance"],
                a["ties"],
                short_labels(a["labels"]),
                b["distance"],
                b["ties"],
                short_labels(b["labels"]),
                row["allele_snp_nearest_overlap"],
            ]
        )
    body += section(
        "Do the methods find the same closest relatives?",
        "<p>Each input is compared with the public genomes. The table shows the smallest observed raw distance "
        "and every tie count; four accession labels are previewed. Distances have different units and "
        "equal numbers do not imply equivalent genetic boundaries.</p>"
        + table(
            [
                "Input",
                "Fewest allele differences",
                "Allele ties",
                "Allele-nearest comparisons",
                "Fewest locus SNPs",
                "SNP ties",
                "SNP-nearest comparisons",
                "Comparisons shared",
            ],
            nearest_rows,
        )
        + '<p><a href="nearest-neighbours-audit.csv">Download all nearest comparisons, dates, countries and coverage</a>. '
        "SNP counts omit aligned gap columns, which are separately included in that evidence file. "
        "Compared loci and nucleotide sites can vary between pairs, even though both methods use the same eligible pairs. "
        "Zero means no differences in the jointly compared evidence; it does not establish whole-genome identity.</p>"
        + f"<details><summary>Could missing DNA explain some differences?</summary><p>We also repeated allele comparisons "
        "using only loci with jointly available aligned nucleotide evidence. With that matched-locus check, "
        f"{matched}/{data['inputs']} inputs have exactly the same closest set as the SNP method "
        f"(versus {exact}/{data['inputs']} using the original allele calls). "
        "The full matched-locus threshold and sampling results are included in the downloads. "
        "This check separates differences in available evidence from differences in how that evidence is measured.</p></details>",
        "nearest",
    )
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for axis, method, color in zip(axes, MAIN, COLORS, strict=True):
        subset = [row for row in thresholds if row["method"] == method]
        axis.plot(
            [int(r["cutoff"]) for r in subset],
            [int(r["linked_contexts"]) for r in subset],
            marker="o",
            color=color,
            label="Connected through paths",
        )
        axis.plot(
            [int(r["cutoff"]) for r in subset],
            [int(r["direct_contexts"]) for r in subset],
            linestyle="--",
            color=color,
            label="Directly within cutoff",
        )
        axis.legend(frameon=False, fontsize=8)
        axis.set(
            title=method,
            xlabel="Allele differences" if method == MAIN[0] else "Aligned locus SNPs",
            ylabel="Public genomes linked to any input",
        )
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.15)
        axis.set_xlim(-0.5, 20.5)
    fig.tight_layout()
    fig.savefig(directory / "neighbourhood-size.svg")
    plt.close(fig)
    main_thresholds = [
        [
            r["method"],
            r["cutoff"],
            r["input_neighbourhoods"],
            fmt(int(r["direct_contexts"])),
            fmt(int(r["linked_contexts"])),
            max(
                row["largest_available_internal_distance"]
                for row in groups
                if row["method"] == r["method"] and row["cutoff"] == int(r["cutoff"])
            ),
        ]
        for r in sorted(thresholds, key=lambda r: (int(r["cutoff"]), r["method"]))
        if r["method"] in MAIN
    ]
    body += section(
        "How wide is the genetic neighbourhood?",
        "<p>These plots show how much public background joins your inputs as the cutoff increases. "
        "The horizontal axes have different units. Each cutoff is a resolution to inspect, not a recommended outbreak threshold. "
        "<strong>Single linkage permits chains:</strong> a very large component can contain pairs much farther apart than its joining cutoff.</p>"
        "<p>Directly within the cutoff means a genome meets that distance to at least one of your inputs. "
        "Connected through paths also includes genomes reached through other comparison genomes. "
        "The plots focus on cutoffs up to 20 in each method's units; the table includes all tested cutoffs.</p>"
        + (directory / "neighbourhood-size.svg").read_text()
        + table(
            [
                "Method",
                "Raw cutoff",
                "Separate input neighbourhoods",
                "Public comparisons directly within cutoff",
                "Public comparisons linked through paths",
                "Largest available distance within an input neighbourhood",
            ],
            main_thresholds,
        ),
        "neighbourhoods",
    )
    aggregation = defaultdict(list)
    for row in samples:
        aggregation[row["method"], int(row["cutoff"]), float(row["context_fraction"])].append(row)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.4))
    robust_rows = []
    for axis, method in zip(axes, MAIN, strict=True):
        cuts = sorted({key[1] for key in aggregation if key[0] == method})
        image = np.full((len(cuts), 3), np.nan)
        for i, cut in enumerate(cuts):
            for j, fraction in enumerate([0.25, 0.5, 0.75]):
                rows = aggregation[method, cut, fraction]
                values = [float(r["peer_recovery_mean"]) for r in rows if r["peer_recovery_mean"]]
                if values:
                    image[i, j] = np.mean(values) * 100
                splits = [int(r["inputs_split"]) for r in rows]
                robust_rows.append(
                    [
                        method,
                        cut,
                        f"{fraction:.0%}",
                        f"{min(values) * 100:.1f}–{max(values) * 100:.1f}%"
                        if values
                        else "No peers",
                        f"{min(splits)}–{max(splits)} / 41",
                    ]
                )
        axis.imshow(image, aspect="auto", vmin=0, vmax=100, cmap="YlGnBu")
        axis.set(
            xticks=range(3),
            xticklabels=["25%", "50%", "75%"],
            yticks=range(len(cuts)),
            yticklabels=cuts,
            title=method,
            xlabel="Public comparisons retained",
            ylabel="Raw distance cutoff",
        )
        for i in range(len(cuts)):
            for j in range(3):
                value = image[i, j]
                axis.text(
                    j,
                    i,
                    "—" if np.isnan(value) else f"{value:.0f}%",
                    ha="center",
                    va="center",
                    color="white" if value > 65 else "#16191e",
                    fontsize=9,
                )
    fig.suptitle("Group members still connected after removing public context", fontsize=12)
    fig.tight_layout()
    fig.savefig(directory / "sampling-robustness.svg")
    plt.close(fig)
    body += section(
        "Do those neighbourhoods depend on which public genomes were sampled?",
        "<p>All 41 inputs stay in every run. We retained 25%, 50% or 75% of public comparisons, "
        "with five fixed-seed alternatives at each size, and rebuilt the groups each time. "
        "The locus alignments and pairwise distances stay fixed; this tests which genomes were included, not alignment uncertainty. "
        "<strong>The score only assesses genomes retained in that run:</strong> it does not penalise deliberately removed genomes.</p>"
        "<p>For each input, we ask: of its original group members that are still present, what percentage "
        "remain connected to it? The figure averages those percentages across inputs and then across the five runs. "
        "An input is counted as split if at least one retained member is no longer connected.</p>"
        "<p>A low percentage means a group split when bridging public genomes were removed. "
        "A high percentage means the retained members stayed connected under this test. "
        "Random removal tests dataset sensitivity; it does not reproduce uneven real-world surveillance.</p>"
        + (directory / "sampling-robustness.svg").read_text()
        + "<details><summary>Ranges across the five alternative samples</summary>"
        + table(
            [
                "Method",
                "Cutoff",
                "Contexts retained",
                "Retained group members still connected",
                "Inputs whose neighbourhood split",
            ],
            robust_rows,
        )
        + "</details><p>Ranges are observed replicate minima and maxima, not confidence intervals. "
        "The displayed percentage averages per-input recovery, then averages the five runs. "
        "Cases with no retained genetic peers have no peer-recovery score.</p>",
        "sampling",
    )
    counts = Counter((row["year"], row["role"]) for row in data["years"])
    years = sorted({row["year"] for row in data["years"]})
    fig, axis = plt.subplots(figsize=(9, 3.6))
    for role, color, offset, label in [
        ("input", COLORS[0], -0.2, "Your inputs"),
        ("context", COLORS[1], 0.2, "Public comparisons"),
    ]:
        denominator = sum(n for (year, r), n in counts.items() if r == role)
        axis.bar(
            np.arange(len(years)) + offset,
            [100 * counts[y, role] / denominator if denominator else 0 for y in years],
            width=0.38,
            color=color,
            label=f"{label} (n={denominator:,})",
        )
    axis.set(
        xticks=range(len(years)),
        xticklabels=years,
        ylabel="% of dated genomes within each dataset",
        xlabel="Collection year",
    )
    axis.tick_params(axis="x", rotation=60)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(directory / "collection-years.svg")
    plt.close(fig)
    group_rows = []
    for row in groups:
        if row["method"] not in MAIN:
            continue
        countries = ", ".join(
            f"{name} ({n})"
            for name, n in sorted(row["countries"].items(), key=lambda p: (-p[1], p[0]))[:4]
        )
        group_rows.append(
            [
                row["method"],
                row["cutoff"],
                row["group"],
                row["input_accessions"],
                row["inputs"],
                row["contexts"],
                f"{row['dated']}/{row['genomes']}",
                row["earliest"] or "Unknown",
                row["latest"] or "Unknown",
                countries,
            ]
        )
    quality = table(
        ["Dataset", "Date status", "Genomes"],
        [
            [
                {"input": "Your inputs", "context": "Public comparisons"}.get(r["role"], r["role"]),
                {
                    "eligible": "Usable date",
                    "before1900": "Before 1900",
                    "invalid_collection_date": "Invalid date",
                    "missing_collection_date": "Date not provided",
                    "future_collection_date": "Future date",
                }.get(r["status"], r["status"]),
                r["count"],
            ]
            for r in data["date_quality"]
        ],
    )
    body += section(
        "Where and when were the members sampled?",
        "<p>Dates and countries describe the genetic groups; they do not determine membership. "
        "The plot compares percentages within each dated dataset because the public collection is much larger. "
        "Year/month precision is preserved in evidence, and the tables show outer date bounds. "
        "Missing or invalid dates and dates before 1900 are excluded from time summaries. "
        "Genomes remain in the genetic analysis.</p>"
        + (directory / "collection-years.svg").read_text()
        + "<details><summary>Time and country detail for every input neighbourhood</summary>"
        + table(
            [
                "Method",
                "Cutoff",
                "Group at this cutoff",
                "Input accessions",
                "Inputs",
                "Public comparisons",
                "Dated",
                "Earliest possible",
                "Latest possible",
                "Largest country counts",
            ],
            group_rows,
        )
        + "</details><details><summary>Date coverage and exclusions</summary>"
        + quality
        + "</details>"
        "<p>These are sampled date distributions, not an estimate of an allele clock, incidence, geographic origin or transmission direction. "
        "A country with more public genomes may simply have sequenced or deposited more genomes.</p>",
        "time",
    )
    body += section(
        "Methods and downloads",
        "<p>Allele profiles retain their original categorical calls. SNPs come from separately aligned cgMLST loci, "
        "with MAFFT provenance and locus coordinates saved. The same pair-coverage mask is used for all comparisons. "
        "An exact minimum spanning forest represents the same single-linkage groups as the complete available graph. "
        "Each resampled graph is rebuilt from its retained distances. Locus SNPs have not been recombination-filtered.</p>"
        "<p>This is an unlabelled public-context comparison. It cannot establish a universally valid epidemiological cutoff. "
        "The original provider reference version remains unknown; the sequence bridge was checked against preserved historical profile fingerprints.</p>"
        '<ul><li><a href="threshold-summary.csv">All cutoffs and neighbourhood counts</a></li>'
        '<li><a href="sampling.csv">Every sampling run and score</a></li>'
        '<li><a href="groups.csv">Group metadata and input accessions</a></li>'
        '<li><a href="nearest-neighbours-audit.csv">Nearest-relative evidence, including ties and indels</a></li>'
        '<li><a href="receipt.json">Source identity, seeds and execution receipt</a></li>'
        '<li><a href="../catalogue/AUDIT.txt">DNA availability and historical-reference bridge audit</a></li>'
        '<li><a href="../alignment/alignment.json">Alignment provenance and locus files</a></li></ul>'
        + (
            '<details><summary>Independent numerical checks</summary><ul>'
            + ''.join(f'<li><a href="audits/{p.name}">{escape(p.stem.replace("-", " "))}</a></li>' for p in sorted((directory / 'audits').glob('*.json')))
            + '</ul><p>These receipts state exactly what was independently checked. The sampling audit covers 288 of 480 score rows.</p></details>'
            if (directory / 'audits').is_dir() else ''
        ),
        "methods",
    )
    nav = (
        '<nav class="contents">'
        + "".join(
            f'<a href="#{ident}">{name}</a>'
            for ident, name in [
                ("start", "Summary"),
                ("dataset", "Dataset"),
                ("nearest", "Closest relatives"),
                ("neighbourhoods", "Neighbourhoods"),
                ("sampling", "Sampling"),
                ("time", "Time and place"),
                ("methods", "Methods"),
            ]
        )
        + "</nav>"
    )
    html = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ST258: allele and locus-SNP comparison</title><style>'
        + report_styles()
        + 'svg{max-width:100%;height:auto}.contents{display:flex;gap:18px;flex-wrap:wrap;padding:18px}details{margin:18px 0}summary{cursor:pointer;font-weight:600}</style></head><body><main class="shell"><header class="identity"><div class="identity-mark">ChronoClade · METHOD COMPARISON</div><div class="identity-copy"><h1>ST258: do allele and SNP views agree?</h1><p>41 Greek inputs · 7,774 public comparisons · sampling and collection-time checks</p></div></header>'
        + nav
        + body
        + "</main></body></html>"
    )
    write_report_html(directory / "report.html", html, encoding="utf-8")
    (directory / "README.md").write_text(
        "# ST258 method comparison\n\nStart with report.html; figures and table controls are embedded, so it opens offline.\nKeep this directory and its parent for evidence downloads.\n\n- threshold-summary.csv: genetic neighbourhood sizes at each descriptive cutoff.\n- sampling.csv: all15 context samples per method, with fixed seeds and induced-graph scores.\n- groups.csv: date/country descriptions for input-containing groups.\n- nearest-neighbours-audit.csv: all closest public relatives, ties, denominators and indels.\n- receipt.json: source references, hashes, seeds and analysis policy.\n- ../catalogue/: historical fingerprint bridge and unresolved DNA audit.\n- ../alignment/: per-locus MSA files, concatenation, partitions and provenance.\n- ../nucleotide-matrices/: memory-mapped evidence arrays used for the calculations.\n\nMissing evidence is never zero distance. Raw allele and SNP cutoffs are different units.\nSingle-linkage chains can join endpoints farther apart than the joining cutoff.\nThe report describes genetic comparisons, not demonstrated transmission or validated outbreak boundaries.\n"
    )
    print(directory / "report.html")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    render(args.directory)
