"""Compare descriptive single-linkage sweeps without selecting a biological cutoff."""

import argparse
from collections import Counter
import csv
from datetime import date
import hashlib
from html import escape
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from chronoclade.cgmlst.thresholds import threshold_sweep
from chronoclade.metadata_dates import sample_date_interval
from chronoclade.report_components.styles import report_styles
from chronoclade.report_components.tables import write_report_html
from research.klebsiella_grouping.report import table
from research.klebsiella_grouping.evaluate import epi_label


def read_rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def metadata_summary(members, metadata):
    """Calendar envelope preserves uncertainty; no dates define membership."""
    intervals = [sample_date_interval(metadata[ident]) for ident in members]
    valid = [value for value in intervals if value["status"] == "valid"]
    states = {}
    for field in ("country", "region", "city", "host"):
        counts = Counter(metadata[ident].get(field) or "Unknown" for ident in members)
        known = {
            k: n for k, n in counts.items() if k.lower() not in {"unknown", "none", "na", "n/a"}
        }
        if known:
            states[field] = dict(sorted(counts.items()))
    return {
        "dated": len(valid),
        "date_precision": dict(sorted(Counter(value["precision"] for value in intervals).items())),
        "earliest_possible": min((value["start"] for value in valid), default=None),
        "latest_possible": max((value["end"] for value in valid), default=None),
        "states": states,
    }


def group_observation(members, metadata, labels, lookup):
    known = []
    absent = 0
    positive = negative = 0
    for i, left in enumerate(members):
        for right in members[:i]:
            relationship = epi_label(labels[left], labels[right])
            positive += relationship == "positive"
            negative += relationship == "negative"
            row = lookup.get(frozenset((left, right)))
            if row is None:
                absent += 1
                continue
            known.append(row["distance"])
    member_set = set(members)
    outside = [row["distance"] for pair, row in lookup.items() if len(pair & member_set) == 1]
    return {
        "members": members,
        "size": len(members),
        "event_isolates": sum(labels[ident]["epi_status"] == "outbreak" for ident in members),
        "controls": sum(labels[ident]["epi_status"] == "control" for ident in members),
        "observed_diameter": max(known, default=0) if not absent else None,
        "maximum_available_pair_distance": max(known, default=None),
        "unavailable_internal_pairs": absent,
        "closest_available_outside_distance": min(outside, default=None),
        "positive_links": positive,
        "negative_links": negative,
        **metadata_summary(members, metadata),
    }


def analyse(mode, unit, samples, rows, field, metadata, labels):
    edges, lookup = [], {}
    for row in rows:
        if row.get("eligible", "True") not in {True, "True", "true", "1"}:
            continue
        value = float(row[field])
        if not value.is_integer():
            raise ValueError("This pilot requires integer raw mismatch counts")
        item = {"distance": int(value), "epi_label": row["epi_label"]}
        key = frozenset((row["sample_a"], row["sample_b"]))
        if key in lookup:
            raise ValueError("Duplicate pair in pilot input")
        lookup[key] = item
        edges.append((row["sample_a"], row["sample_b"], item["distance"]))
    maximum = int(max((edge[2] for edge in edges), default=0))
    sweep = threshold_sweep(samples, edges, range(maximum + 1), distance_unit=unit)
    groups, curves = {}, []
    for partition in sweep["partitions"]:
        threshold = partition["threshold"]
        total_positive = total_negative = 0
        for members in partition["clusters"]:
            key = hashlib.sha256("\n".join(members).encode()).hexdigest()[:12]
            if key not in groups:
                groups[key] = {
                    "group_id": key,
                    "first_cutoff": threshold,
                    "last_cutoff": threshold,
                    **group_observation(members, metadata, labels, lookup),
                }
            groups[key]["last_cutoff"] = threshold
            total_positive += groups[key]["positive_links"]
            total_negative += groups[key]["negative_links"]
        curves.append(
            {
                "mode": mode,
                "cutoff": threshold,
                "clusters": len(partition["clusters"]),
                "largest": max(map(len, partition["clusters"]), default=0),
                "event_pairs_connected": total_positive,
                "event_control_pairs_connected": total_negative,
            }
        )
    denominator_fields = {
        "Allele profiles": ["shared_loci"],
        "Locus nucleotides": ["compared_loci", "compared_sites"],
        "Corrected core SNPs": ["callable_sites"],
    }[mode]
    coverage = {}
    for denominator in denominator_fields:
        values = [
            int(row[denominator])
            for row in rows
            if frozenset((row["sample_a"], row["sample_b"])) in lookup
        ]
        coverage[denominator] = {"min": min(values, default=None), "max": max(values, default=None)}
    return {
        "mode": mode,
        "unit": unit,
        "eligible_pairs": len(edges),
        "coverage": coverage,
        "zero_nucleotide_pairs_with_allele_differences": sum(
            int(row[field]) == 0 and int(row["allele_mismatches"]) > 0
            for row in rows
            if row.get("eligible") in {True, "True"}
        )
        if mode == "Locus nucleotides"
        else None,
        "sweep": sweep,
        "groups": list(groups.values()),
        "curve": curves,
    }


def render(output, results, metadata):
    memberships = {
        group["group_id"]: group["members"] for result in results for group in result["groups"]
    }
    names = {
        ident: f"Group {index:02d}"
        for index, ident in enumerate(
            sorted(memberships, key=lambda ident: (-len(memberships[ident]), ident)), 1
        )
    }
    fig, axes = plt.subplots(1, len(results), figsize=(5 * len(results), 3.8), squeeze=False)
    colours = {"event_pairs_connected": "#71527f", "event_control_pairs_connected": "#36828c"}
    for axis, result in zip(axes[0], results, strict=True):
        for field, name in (
            ("event_pairs_connected", "Event pairs connected / 325"),
            ("event_control_pairs_connected", "Event–control pairs connected / 130"),
        ):
            axis.step(
                [row["cutoff"] for row in result["curve"]],
                [row[field] for row in result["curve"]],
                where="post",
                label=name,
                color=colours[field],
                linewidth=2,
            )
        axis.set(title=result["mode"], xlabel=result["unit"], ylabel="Labelled pairs connected")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.18)
        axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "threshold-curves.svg")
    plt.close(fig)
    parts = [
        '<section class="stage"><div class="stage-body"><h2>Genetic resolution across cutoffs</h2>'
        "<p>These are three different distance scales on the same 31 genomes. "
        "No cutoff is selected. Curves count labelled pairs connected through a path, "
        "not necessarily pairs whose direct distance meets the cutoff. "
        "The 325 event pairs and 130 event–control pairs share isolates.</p>",
        (output / "threshold-curves.svg").read_text(),
        "</div></section>",
    ]
    all_summaries = []
    for result in results:
        interesting = sorted(
            (group for group in result["groups"] if group["size"] > 1),
            key=lambda group: (group["first_cutoff"], -group["size"], group["group_id"]),
        )
        summary = []
        for group in interesting:
            states = "; ".join(
                f"{field}: " + ", ".join(f"{k} ({v})" for k, v in counts.items())
                for field, counts in group["states"].items()
            )
            dates = f"{group['earliest_possible'] or 'Unknown'} – {group['latest_possible'] or 'Unknown'}"
            summary.append(
                [
                    names[group["group_id"]],
                    f"{group['first_cutoff']:g}–{group['last_cutoff']:g}",
                    group["size"],
                    group["event_isolates"],
                    group["controls"],
                    group["observed_diameter"]
                    if group["observed_diameter"] is not None
                    else "Unavailable",
                    group["closest_available_outside_distance"]
                    if group["closest_available_outside_distance"] is not None
                    else "None available",
                    f"{group['dated']}/{group['size']}",
                    dates,
                    states or "Unknown",
                ]
            )
        for group in result["groups"]:
            all_summaries.append(
                {
                    "mode": result["mode"],
                    "display_group": names[group["group_id"]],
                    **{
                        k: json.dumps(v) if isinstance(v, (dict, list)) else v
                        for k, v in group.items()
                    },
                }
            )
        parts.extend(
            [
                '<section class="stage"><div class="stage-body"><h2>'
                + escape(result["mode"])
                + "</h2><p>Matched eligible pairs: "
                + str(result["eligible_pairs"])
                + ". Coverage: "
                + escape(
                    "; ".join(
                        f"{key}: {value['min']}–{value['max']}"
                        for key, value in result["coverage"].items()
                    )
                )
                + ". Numerators and coverage are specific to this distance definition.</p>",
                (
                    "<p>Alleles were aligned per locus with MAFFT, reconstructed by genome and concatenated. "
                    "SNPs count jointly A/C/G/T differences; indels remain explicit gap evidence. "
                    "Zero SNPs does not mean an identical gene panel. Length-changing loci are retained.</p>"
                    if result.get("locus_alignment") == "aligned" else
                    "<p><strong>Excluded differences:</strong> "
                    + str(result["zero_nucleotide_pairs_with_allele_differences"])
                    + " pairs have zero compared nucleotide differences despite differing allele profiles, "
                    "because their length-changing allele pairs were excluded. Zero here does not mean "
                    "identical gene panels. Every retained differing allele contributes one substitution "
                    "in this panel, so this run adds no substitution resolution over those retained allele calls.</p>"
                    if result["mode"] == "Locus nucleotides"
                    else ""
                ),
                "<p>Each row is one distinct membership set. Cutoff range means membership "
                "is unchanged at these integer cuts in this panel; it does not establish an outbreak. "
                "Internal diameter can exceed the joining cutoff because single linkage permits chains.</p>",
                "<p>Group labels identify the same membership set across all three views in this report; "
                "they are not LIN codes or permanent nomenclature. Full membership is in the CSV download.</p>",
                table(
                    [
                        "Group",
                        "Cutoffs",
                        "Genomes",
                        "Event",
                        "Controls",
                        "Internal diameter",
                        "Closest outside",
                        "Dated",
                        "Sampled date envelope",
                        "Reported states",
                    ],
                    summary,
                ),
                "</div></section>",
            ]
        )
        (output / (result["mode"].lower().replace(" ", "-") + ".json")).write_text(
            json.dumps(result, indent=2) + "\n"
        )
    write_rows(output / "groups.csv", all_summaries)
    write_rows(output / "thresholds.csv", [row for result in results for row in result["curve"]])
    fig, axis = plt.subplots(figsize=(9, 3.8))
    for i, (ident, row) in enumerate(
        sorted(metadata.items(), key=lambda pair: (pair[1]["collection_date"], pair[0]))
    ):
        interval = sample_date_interval(row)
        if interval["status"] != "valid":
            continue
        start, end = (date.fromisoformat(interval[key]) for key in ("start", "end"))
        axis.plot([start, end], [i, i], color="#36828c")
        axis.scatter([start], [i], s=16, color="#71527f")
    axis.set(
        xlabel="Observed collection date (intervals preserve date precision)",
        ylabel="All 31 genomes",
    )
    axis.set_yticks([])
    axis.spines[["top", "right"]].set_visible(False)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(output / "sampling-dates.svg")
    plt.close(fig)
    parts.extend(
        [
            '<section class="stage"><div class="stage-body"><h2>Time as a separate layer</h2>'
            "<p>Dates never change genetic membership. Group tables show observed dates and "
            "the outer calendar envelope, including date uncertainty. Missing dates remain missing. "
            "The short outbreak collection cannot calibrate a general alleles-per-year rule. "
            "Country is reported for the 26 documented Beijing event isolates; control locations remain unknown.</p>",
            (output / "sampling-dates.svg").read_text(),
            "</div></section>",
            '<section class="stage"><div class="stage-body"><h2>Methods and evidence</h2>'
            "<p>Locus DNA and corrected core SNPs use different genomic targets. "
            "The locus alignment has not been recombination-filtered. Corrected core SNPs come from "
            "the separately flagged diagnostic run, retaining all 31 genomes. Threshold ranges describe "
            "this sample collection only. No confidence intervals assume pairwise independence.</p>"
            '<p><a href="groups.csv">Every membership set</a> · <a href="thresholds.csv">All threshold counts</a> · '
            '<a href="receipt.json">Input hashes</a> · <a href="../evaluation-primary/report.html">Original labelled comparison</a>'
            ' · <a href="../locus-aligned/summary.json">Aligned-locus comparison audit</a></p>'
            '<p>Normal ChronoClade reports: <a href="../cgmlst/st_f99ec4b312b5/profile_report.html">Allele profiles</a> · '
            '<a href="../cgmlst-locus-aligned-final/st_f99ec4b312b5/profile_report.html">Locus nucleotides</a>. '
            "The direct-DNA run uses the documented descriptive fraction 0.00001 to exercise the workflow; "
            "it is not a calibrated epidemiological cutoff. This sweep does not use that selected fraction.</p></div></section>",
        ]
    )
    html = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ChronoClade distance comparison</title><style>'
        + report_styles()
        + 'svg{max-width:100%;height:auto}</style></head><body><main class="shell"><header class="identity"><div class="identity-mark">ChronoClade · RESEARCH PILOT</div><div class="identity-copy"><h1>ST11: allele and SNP neighbourhoods</h1><p>31 genomes · three distance definitions · exploratory threshold analysis</p></div></header>'
        + "".join(parts)
        + "</main></body></html>"
    )
    write_report_html(output / "report.html", html, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--locus-pairs", type=Path, required=True)
    parser.add_argument("--locus-mode", choices=("aligned", "positional"), default="aligned")
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    metadata = {row["sample_id"]: row for row in read_rows(args.metadata)}
    labels = {row["sample_id"]: row for row in read_rows(args.labels)}
    if set(metadata) != set(labels) or len(metadata) != 31:
        raise ValueError("Expected the frozen 31-genome panel")
    samples = sorted(metadata)
    primary = read_rows(args.primary / "allele-snp-pairs.csv")
    locus = read_rows(args.locus_pairs)
    for rows in (primary, locus):
        identities = {row[key] for row in rows for key in ("sample_a", "sample_b")}
        keys = {frozenset((row["sample_a"], row["sample_b"])) for row in rows}
        if (
            identities != set(samples)
            or len(rows) != 465
            or len(keys) != 465
            or any(len(key) != 2 for key in keys)
        ):
            raise ValueError("Pair inputs must cover the same 31 genomes exactly")
        if any(
            row["epi_label"] != epi_label(labels[row["sample_a"]], labels[row["sample_b"]])
            for row in rows
        ):
            raise ValueError("Pair labels differ from the frozen independent labels")
    available = [
        {
            frozenset((row["sample_a"], row["sample_b"]))
            for row in rows
            if row.get("eligible", "True") in {True, "True", "true", "1"}
        }
        for rows in (primary, locus)
    ]
    common = set.intersection(*available)
    for rows in (primary, locus):
        for row in rows:
            row["eligible"] = frozenset((row["sample_a"], row["sample_b"])) in common
    results = [
        analyse(
            "Allele profiles",
            "Allele mismatches",
            samples,
            primary,
            "allele_mismatches",
            metadata,
            labels,
        ),
        analyse(
            "Locus nucleotides",
            "Nucleotide differences in cgMLST loci",
            samples,
            locus,
            "locus_snps",
            metadata,
            labels,
        ),
        analyse(
            "Corrected core SNPs",
            "Recombination-filtered core SNPs",
            samples,
            primary,
            "clonal_snps",
            metadata,
            labels,
        ),
    ]
    results[1]["locus_alignment"] = args.locus_mode
    args.output.mkdir(parents=True, exist_ok=False)
    render(args.output, results, metadata)
    sources = [args.primary / "allele-snp-pairs.csv", args.locus_pairs, args.metadata, args.labels]
    (args.output / "receipt.json").write_text(
        json.dumps(
            {
                "sources": {
                    str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources
                },
                "method": "inclusive single linkage, all integer raw-distance cutoffs, no biological selection",
                "time_policy": "descriptive sample-date intervals, no fitted rate or time filter",
                "code": {
                    str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in (Path(__file__), Path("chronoclade/cgmlst/thresholds.py"))
                },
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
