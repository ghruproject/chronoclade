"""Join the frozen primary pairs to the separately flagged diagnostic SNP run."""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path

from Bio import Phylo
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata


def compare(primary, tree_manifest, labels_path):
    primary = Path(primary)
    tree_manifest = Path(tree_manifest).resolve()
    result = json.loads(tree_manifest.read_text())
    if result["parameters"]["coherence_policy"] != "diagnostic":
        raise ValueError("This amended comparison requires the separately declared diagnostic run")
    with (primary / "pairs.csv").open() as stream:
        pairs = list(csv.DictReader(stream))
    with Path(labels_path).open() as stream:
        labels = {r["sample_id"]: r for r in csv.DictReader(stream)}
    if set(result["selected_sample_ids"]) != set(labels):
        raise ValueError("SNP identities must exactly match all labelled genomes")
    corrected = tree_manifest.parent / result["corrected_tree"]["path"]
    source = corrected.parent / "snp_distances.json"
    snps = json.loads(source.read_text())
    lookup = {frozenset((p["sample_1"], p["sample_2"])): p for p in snps}
    if len(lookup) != len(pairs) or len(lookup) != len(snps):
        raise ValueError("Missing or duplicate SNP pairs")
    joined = []
    for pair in pairs:
        match = lookup[frozenset((pair["sample_a"], pair["sample_b"]))]
        joined.append(
            {**pair, "clonal_snps": match["clonal_snps"], "callable_sites": match["callable_sites"]}
        )
    with (primary / "allele-snp-pairs.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(joined[0]))
        writer.writeheader()
        writer.writerows(joined)
    distributions = []
    for label in ["positive", "negative", "unknown"]:
        subset = [p for p in joined if p["epi_label"] == label]
        snp_values = [p["clonal_snps"] for p in subset]
        allele_values = [int(p["allele_mismatches"]) for p in subset]
        x, y = rankdata(allele_values), rankdata(snp_values)
        distributions.append(
            {
                "label": label,
                "pairs": len(subset),
                "snp_min": min(snp_values),
                "snp_median": float(np.median(snp_values)),
                "snp_max": max(snp_values),
                "allele_min": min(allele_values),
                "allele_max": max(allele_values),
                "spearman_point": float(np.corrcoef(x, y)[0, 1])
                if np.std(x) > 0 and np.std(y) > 0
                else None,
            }
        )
    case_ids = {ident for ident, row in labels.items() if row["epi_status"] == "outbreak"}
    phylogeny = Phylo.read(corrected, "newick")
    tips = {tip.name for tip in phylogeny.get_terminals()}
    if tips != set(labels):
        raise ValueError("Corrected tree tips differ from frozen labels")
    event_split = any(
        {t.name for t in clade.get_terminals()} in (case_ids, tips - case_ids)
        for clade in phylogeny.find_clades()
    )
    nearest = []
    snp_nearest_ids = {}
    for ident in sorted(case_ids):
        candidate = [p for p in joined if ident in {p["sample_a"], p["sample_b"]}]
        amin = min(int(p["allele_mismatches"]) for p in candidate)
        smin = min(p["clonal_snps"] for p in candidate)
        aid = {
            p["sample_b"] if p["sample_a"] == ident else p["sample_a"]
            for p in candidate
            if int(p["allele_mismatches"]) == amin
        }
        sid = {
            p["sample_b"] if p["sample_a"] == ident else p["sample_a"]
            for p in candidate
            if p["clonal_snps"] == smin
        }
        snp_nearest_ids[ident] = sid
        nearest.append(
            {
                "query": ident,
                "allele_ties": len(aid),
                "snp_ties": len(sid),
                "snp_nearest_retained_in_allele_ties": len(aid & sid),
                "snp_nearest_total": len(sid),
            }
        )
    candidate_retention = defaultdict(lambda: {"retained": 0, "available": 0})
    with (primary / "queries.csv").open() as stream:
        for row in csv.DictReader(stream):
            if row["sample_id"] not in snp_nearest_ids:
                continue
            wanted = snp_nearest_ids[row["sample_id"]]
            retrieved = set(json.loads(row["retrieved_ids"]))
            candidate_retention[row["candidate"]]["retained"] += len(wanted & retrieved)
            candidate_retention[row["candidate"]]["available"] += len(wanted)
    zero = [p["clonal_snps"] for p in joined if int(p["allele_mismatches"]) == 0]
    fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
    for label, title, colour in [
        ("positive", "Same event", "#754b8e"),
        ("negative", "Event versus sporadic", "#3c7383"),
        ("unknown", "Sporadic versus sporadic (unknown)", "#7f8792"),
    ]:
        points = defaultdict(int)
        for p in joined:
            if p["epi_label"] == label:
                points[int(p["allele_mismatches"]), p["clonal_snps"]] += 1
        ax.scatter(
            [p[0] for p in points],
            [p[1] for p in points],
            s=[15 + 8 * np.sqrt(n) for n in points.values()],
            alpha=0.65,
            label=title,
            color=colour,
            edgecolors="white",
            linewidths=0.4,
        )
    ax.set_xlabel("Observed cgMLST allele differences")
    ax.set_ylabel("Recombination-filtered core SNP differences")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title(
        "Same 31 genomes · diagnostic cohort policy · marker area reflects tied pairs", fontsize=11
    )
    fig.savefig(primary / "allele-snp-comparison.svg")
    plt.close(fig)
    summary = {
        "scope": "Secondary exploratory diagnostic; strict failure retained, all31 genomes retained, no temporal assessment",
        "candidate_snp_nearest_retention": dict(candidate_retention),
        "identities": len(tips),
        "pairs": len(joined),
        "callable_sites_min": min(p["callable_sites"] for p in joined),
        "callable_sites_max": max(p["callable_sites"] for p in joined),
        "event_split_present": event_split,
        "label_distributions": distributions,
        "zero_allele_pairs": {
            "pairs": len(zero),
            "snp_min": min(zero),
            "snp_median": float(np.median(zero)),
            "snp_max": max(zero),
        },
        "outbreak_queries": nearest,
        "snp_nearest_retained_in_allele_ties": sum(
            r["snp_nearest_retained_in_allele_ties"] for r in nearest
        ),
        "snp_nearest_total": sum(r["snp_nearest_total"] for r in nearest),
        "corrected_coherence_flags": result["coherence_screening"]["corrected"]["flagged_samples"],
        "temporal_code": result["temporal_assessment"]["code"],
        "tree_report_relative_path": os.path.relpath(
            tree_manifest.parent / result["report"]["path"], primary.resolve()
        ),
        "sources_sha256": {
            "tree_manifest": hashlib.sha256(tree_manifest.read_bytes()).hexdigest(),
            "snp_pairs": hashlib.sha256(source.read_bytes()).hexdigest(),
            "corrected_tree": hashlib.sha256(corrected.read_bytes()).hexdigest(),
        },
    }
    (primary / "snp-comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    args = parser.parse_args()
    compare(args.primary, args.tree, args.labels)
