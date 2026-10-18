"""A short, offline explanation of the actual ReporTree reproduction."""

from collections import Counter
import csv
from html import escape
import json
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from chronoclade.report_components.styles import report_styles
from chronoclade.report_components.tables import write_report_html
from research.klebsiella_grouping.report import table
from research.klebsiella_grouping.st258_report import section


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def cut(name):
    return int(float(re.search(r"-(\d+(?:\.\d+)?)x", name)[1]))


def regions(path):
    result = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        block, first, last, length = line.split("\t")
        # A transition 40->39 through54->53 covers partition levels39..54.
        result.append(
            {
                "block": block,
                "lower": cut(first.split("->")[1]),
                "upper": cut(last.split("->")[0]),
                "transitions": int(length),
            }
        )
    return result


def main():
    root = Path(__file__).resolve().parents[2]
    base = root / "tmp/reportree-reference"
    fresh = base / "reproduced"
    verified = json.loads((base / "verification.json").read_text())
    execution = json.loads((fresh / "execution.json").read_text())
    blocks = regions(fresh / "NG_822_stableRegions.tsv")
    metrics = rows(fresh / "NG_822_metrics.tsv")
    metadata = rows(fresh / "NG_822_metadata_w_partitions.tsv")
    partitions = rows(fresh / "NG_822_partitions.tsv")
    result = base / "report"
    result.mkdir(exist_ok=True)
    # Fresh output, never the expected files, drives every displayed result.
    fig, axes = plt.subplots(2, 1, figsize=(10, 5.7), sharex=True)
    x = [cut(r["col2"]) for r in metrics]
    axes[0].plot(x, [float(r["AdjWallace"]) for r in metrics], color="#71527f", lw=1.2)
    axes[0].axhline(0.99, color="#555", linestyle="--", lw=0.8)
    axes[0].set(ylabel="Membership agreement\n(Adjusted Wallace)", ylim=(-0.04, 1.04))
    names = list(partitions[0])[1:]
    axes[1].plot(
        [cut(n) for n in names], [len({r[n] for r in partitions}) for n in names], color="#36828c"
    )
    axes[1].set(
        ylabel="Groups, including singletons",
        xlabel="Allele-distance threshold on the MST",
        yscale="log",
    )
    for axis in axes:
        for block in blocks:
            axis.axvspan(block["lower"], block["upper"], alpha=0.13, color="#36828c")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.12)
    fig.tight_layout()
    fig.savefig(result / "reference-stability.svg")
    plt.close(fig)
    body = section(
        "The example we followed",
        "<p>We ran <strong>ReporTree itself</strong> on the authors' "
        "Neisseria gonorrhoeae example: <strong>3,791 genomes, 822 loci and 21 countries</strong>. "
        "It builds an MST from allele profiles, groups genomes at successive thresholds, then "
        "summarises those groups using the supplied metadata.</p>"
        "<p><strong>The question is simple:</strong> are there ranges of genetic resolution where "
        "group membership changes little? This is different from removing public genomes.</p>"
        f"<p>The fresh run found <strong>{len(blocks)} stable ranges</strong>. The independent comparison "
        f"result is <strong>{escape(str(verified['overall']))}</strong>: stable ranges, membership scores, "
        "reported group memberships and country/date summaries match the supplied reference. "
        f"The run took {execution['elapsed_seconds'] / 60:.1f} minutes.</p>",
        "summary",
    )
    body += section(
        "What does ReporTree call stable?",
        "<p>At each step, the threshold increases by one allele. ReporTree compares the new grouping "
        "with the previous one using Adjusted Wallace. Following the published example, a range "
        "qualifies when that score stays at least <strong>0.99 for nine consecutive steps</strong>. "
        "These settings come from the example; they are not validated Klebsiella cutoffs.</p>"
        "<p>Green shading marks qualifying ranges. The lower panel shows how many groups remain. "
        "Stability means similar membership across thresholds; it does not mean persistence over years "
        "or prove an outbreak. Identical zero-threshold groups are represented once when scoring, "
        "following the tool's default redundancy handling.</p>"
        + (result / "reference-stability.svg").read_text()
        + table(
            ["Stable partition levels", "Consecutive comparisons"],
            [[f"{b['lower']}–{b['upper']}", b["transitions"]] for b in blocks],
        ),
        "stability",
    )
    chosen = [blocks[0]["lower"], blocks[1]["lower"]]
    summary = []
    detail = []
    for level in chosen:
        name = f"MST-{level}x1.0"
        groups = {}
        for row in metadata:
            groups.setdefault(row[name], []).append(row)
        summary.append(
            [
                level,
                len(groups),
                sum(len(v) >= 2 for v in groups.values()),
                sum(len(v) == 1 for v in groups.values()),
                max(map(len, groups.values())),
            ]
        )
        for label, members in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0])):
            if len(members) < 2:
                continue
            countries = Counter(r["country"] for r in members)
            dates = [
                r["date"] for r in members if r.get("date") and r["date"] not in ("nan", "EMPTY")
            ]
            detail.append(
                [
                    level,
                    label,
                    len(members),
                    len(countries),
                    "; ".join(f"{k}: {n}" for k, n in countries.most_common()),
                    min(dates) if dates else "Unknown",
                    max(dates) if dates else "Unknown",
                ]
            )
    body += section(
        "Two resolutions of the same dataset",
        f"<p>The first two stable ranges begin at {chosen[0]} and {chosen[1]} alleles. "
        "These are example resolutions to inspect, not two competing answers. The broader level "
        "merges some of the finer groups. The authors' existing group names are preserved where applicable.</p>"
        + table(
            ["Threshold", "All groups", "Groups with ≥2 genomes", "Singletons", "Largest group"],
            summary,
        )
        + "<details><summary>Countries and supplied dates for the groups at these two resolutions</summary>"
        + table(
            [
                "Threshold",
                "Group",
                "Genomes",
                "Countries",
                "Country counts",
                "Earliest supplied date",
                "Latest supplied date",
            ],
            detail,
        )
        + "</details><p>This is the useful report pattern: select a genetic resolution, see the groups, "
        "then describe where and when their members were sampled.</p>",
        "groups",
    )
    # Apply the same upstream stability calculation to our existing exact partitions.
    fig, axes = plt.subplots(2, 2, figsize=(10, 5.8), sharex="col")
    strows = []
    for j, (slug, title) in enumerate(
        (("allele-profiles", "ST258 · alleles"), ("aligned-locus-snps", "ST258 · locus SNPs"))
    ):
        directory = base / "st258" / slug
        m = rows(directory / "ST258_metrics.tsv")
        b = regions(directory / "ST258_stableRegions.tsv")
        p = rows(directory / "partitions.tsv")
        x = [cut(r["col2"]) for r in m]
        axes[0, j].plot(x, [float(r["AdjWallace"]) for r in m], color="#71527f", lw=1.3)
        axes[0, j].axhline(0.99, color="#555", linestyle="--", lw=0.8)
        axes[0, j].set(title=title, ylim=(-0.04, 1.04))
        cols = list(p[0])[1:]
        axes[1, j].plot(
            [cut(c) for c in cols], [len({r[c] for r in p}) for c in cols], color="#36828c"
        )
        axes[1, j].set(xlabel="Allele cutoff" if j == 0 else "Locus-SNP cutoff", yscale="log")
        if not b:
            strows.append([title, "None", "—", "—"])
        for block in b:
            sizes = Counter(r[f"single-{block['lower']}x1.0"] for r in p)
            strows.append(
                [
                    title,
                    f"{block['lower']}–{block['upper']}",
                    len(sizes),
                    ", ".join(map(str, sorted(sizes.values(), reverse=True))),
                ]
            )
        for axis in axes[:, j]:
            for block in b:
                axis.axvspan(block["lower"], block["upper"], alpha=0.13, color="#36828c")
            axis.spines[["top", "right"]].set_visible(False)
    axes[0, 0].set_ylabel("Membership agreement")
    axes[1, 0].set_ylabel("Groups, including singletons")
    fig.tight_layout()
    fig.savefig(result / "st258-stability.svg")
    plt.close(fig)
    body += section(
        "What happens when we apply that test to ST258?",
        "<p>We supplied our existing full ST258 single-linkage partitions to <strong>ReporTree's "
        "same stability analysis</strong>, at every integer threshold up to the final merge. "
        "The 7,815 genomes and shared coverage rules are unchanged. This applies its stability test; "
        "it does not rerun GrapeTree or claim its missing-data distances are identical to ours.</p>"
        + (result / "st258-stability.svg").read_text()
        + table(
            [
                "Representation",
                "Stable partition levels",
                "Groups at start",
                "Group sizes at start",
            ],
            strows,
        )
        + "<p><strong>Check the sizes before interpreting a plateau.</strong> A stable result consisting "
        "of almost the whole dataset plus one outlier is not a useful set of fine-scale clusters. "
        "No region under these settings means the test did not select one; it does not prove "
        "that every biologically relevant subgroup is absent.</p>",
        "st258",
    )
    body += section(
        "What this establishes—and what comes next",
        "<p>This anchors ChronoClade to a concrete existing workflow: a genetic hierarchy, "
        "membership stability across thresholds, and country/date descriptions. It does not "
        "automatically supply an epidemiologically meaningful cutoff for ST258. We should retain "
        "the hierarchy and make local groups around the inputs easy to inspect, rather than invent "
        "a single winning threshold.</p><p>Alleles remain the default; aligned locus SNPs remain the "
        "finer sequence comparison. There is no conversion between their thresholds.</p>"
        '<ul><li><a href="../verification.json">Independent reproduction checks</a></li>'
        '<li><a href="../reproduced/execution.json">Exact command, versions and input hashes</a></li>'
        '<li><a href="../reproduced/NG_822_metadata_w_partitions.tsv">Fresh Neisseria metadata and groups</a></li>'
        '<li><a href="../reproduced/NG_822.nwk">Fresh GrapeTree MST as Newick</a></li>'
        '<li><a href="../reproduced/NG_822_flt_samples_matrix.tsv">Allele matrix used by GrapeTree</a></li>'
        '<li><a href="../st258/execution.json">ST258 commands and settings</a></li>'
        '<li><a href="https://github.com/insapathogenomics/ReporTree/wiki/5.-Examples#large-scale-genetic-clustering-and-linkage-to-antibiotic-resistance-data-eg-neisseria-gonorrhoeae">Authors’ worked example, with GrapeTree figures</a></li></ul>'
        "<p>Fresh run: ReporTree 2.7.0. Supplied reference: 2.0.0. No upstream source code was changed. "
        "Source and dependency revisions are saved alongside the run.</p>",
        "methods",
    )
    nav = (
        '<nav class="contents">'
        + "".join(
            f'<a href="#{i}">{n}</a>'
            for i, n in (
                ("summary", "Summary"),
                ("stability", "Stability"),
                ("groups", "Example groups"),
                ("st258", "ST258"),
                ("methods", "Evidence"),
            )
        )
        + "</nav>"
    )
    html = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ChronoClade · ReporTree worked example</title><style>'
        + report_styles()
        + 'svg{max-width:100%;height:auto}.contents{display:flex;gap:18px;flex-wrap:wrap;padding:18px}details{margin:18px 0}summary{cursor:pointer;font-weight:600}</style></head><body><main class="shell"><header class="identity"><div class="identity-mark">ChronoClade · WORKED EXAMPLE</div><div class="identity-copy"><h1>Following ReporTree, step by step</h1><p>A published allele-profile example, reproduced locally, then applied to ST258</p></div></header>'
        + nav
        + body
        + "</main></body></html>"
    )
    write_report_html(result / "report.html", html, encoding="utf-8")
    with (result / "groups-at-two-resolutions.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "threshold",
                "group",
                "genomes",
                "countries",
                "country_counts",
                "first_date",
                "last_date",
            ]
        )
        writer.writerows(detail)
    print(result / "report.html")


if __name__ == "__main__":
    main()
