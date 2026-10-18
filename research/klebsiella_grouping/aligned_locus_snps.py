"""Run the frozen labelled-panel comparison after per-locus MAFFT alignment."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

from chronoclade.cgmlst.locus_alignment import align_catalogue, file_hash, save_json
from chronoclade.context_partitions import profile_records
from chronoclade.datasets import load_dataset
from research.klebsiella_grouping.locus_snps import correlation, distribution, write_csv


def read_pairs(path):
    with Path(path).open() as stream:
        rows = list(csv.DictReader(stream))
    result = {frozenset((row["sample_a"], row["sample_b"])): row for row in rows}
    if len(result) != len(rows) or any(len(key) != 2 for key in result):
        raise ValueError("Duplicate or self pair in comparison input")
    return result


def run(dataset_path, catalogue, primary_pairs, core_pairs, positional_pairs, output, mafft, min_loci=599):
    dataset = load_dataset(dataset_path)
    records = profile_records(dataset)
    primary, core, positional = (read_pairs(path) for path in (primary_pairs, core_pairs, positional_pairs))
    if not set(primary) == set(core) == set(positional):
        raise ValueError("Comparisons must retain the identical frozen pair universe")
    n = len(records)
    if len(primary) != n * (n - 1) // 2:
        raise ValueError("Comparison requires every labelled sample pair")
    alignment = align_catalogue(catalogue, records, output, mafft)
    rows = []
    for evidence in alignment.pairs(min_loci):
        key = frozenset((evidence["sample_a"], evidence["sample_b"]))
        p, c, previous = primary[key], core[key], positional[key]
        rows.append({**evidence, "epi_label": p["epi_label"],
                     "allele_mismatches": int(p["allele_mismatches"]),
                     "shared_loci": int(p["shared_loci"]),
                     "clonal_snps": int(c["clonal_snps"]), "callable_sites": int(c["callable_sites"]),
                     "positional_locus_snps": int(previous["locus_snps"]),
                     "positional_compared_loci": int(previous["compared_loci"]),
                     "positional_compared_sites": int(previous["compared_sites"])})
    output = Path(output)
    write_csv(output / "pairs.csv", rows)
    partitions = alignment.partitions
    zero = [row for row in rows if row["allele_mismatches"] == 0]
    missing = json.loads((output / "missing-calls.json").read_text())
    summary = {"schema": "chronoclade.research.aligned-locus-nucleotides", "schema_version": 1,
               "dataset_id": dataset.dataset_id, "samples": n, "pairs": len(rows),
               "eligible_pairs": sum(row["eligible"] for row in rows), "minimum_compared_loci": min_loci,
               "canonical_loci": len(partitions), "represented_loci": sum(p["width"] > 0 for p in partitions),
               "absent_loci": [p["locus"] for p in partitions if not p["width"]],
               "mafft_loci": [p["locus"] for p in partitions if p["method"] == "mafft"],
               "alignment_sites": alignment.provenance["alignment_sites"],
               "missing": dict(Counter(r["reason"] for r in missing)),
               "independent_arithmetic_verified_pairs": sum(row["independent_arithmetic_verified"] for row in rows),
               "pairs_changed_from_positional": sum(row["locus_snps"] != row["positional_locus_snps"] for row in rows),
               "zero_allele_pairs": len(zero),
               "zero_allele_aligned_snps": distribution([r["locus_snps"] for r in zero]) if zero else None,
               "zero_allele_core_snps": distribution([r["clonal_snps"] for r in zero]) if zero else None,
               "zero_aligned_snps_nonzero_alleles": sum(r["locus_snps"] == 0 and r["allele_mismatches"] > 0 for r in rows),
               "zero_aligned_snps_with_indels": sum(r["locus_snps"] == 0 and r["gap_vs_base_columns"] > 0 for r in rows),
               "distributions": [], "alignment_provenance": "alignment.json",
               "sources": [{"path": str(Path(p).resolve()), "sha256": file_hash(p)}
                           for p in (dataset_path, catalogue, primary_pairs, core_pairs, positional_pairs, __file__)],
               "limitations": ["Per-locus MSA depends on the frozen observed allele set and MAFFT parameters.",
                               "Joint ACGT substitution counts exclude aligned gaps and missing calls; indel evidence is reported separately.",
                               "No recombination correction or core-genome equivalence is implied."]}
    for label in ("all", "positive", "negative", "unknown"):
        subset = rows if label == "all" else [row for row in rows if row["epi_label"] == label]
        summary["distributions"].append({"label": label, "pairs": len(subset),
            **{field: distribution([r[field] for r in subset]) for field in
               ("locus_snps", "compared_loci", "compared_sites", "gap_vs_base_columns", "indel_tracts",
                "both_gap_columns", "missing_columns", "allele_mismatches", "clonal_snps")},
            "allele_locus_spearman": correlation(subset, "allele_mismatches", "locus_snps"),
            "locus_core_spearman": correlation(subset, "locus_snps", "clonal_snps")})
    save_json(output / "summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default="tmp/science-zhou2017/prepared/dataset.json")
    parser.add_argument("--catalogue", default="tmp/science-zhou2017/locus-dna/catalogue.json")
    parser.add_argument("--primary-pairs", default="tmp/science-zhou2017/evaluation-primary/pairs.csv")
    parser.add_argument("--core-pairs", default="tmp/science-zhou2017/evaluation-primary/allele-snp-pairs.csv")
    parser.add_argument("--positional-pairs", default="tmp/science-zhou2017/locus-dna/pairs.csv")
    parser.add_argument("--output", default="tmp/science-zhou2017/locus-aligned")
    parser.add_argument("--min-loci", type=int, default=599)
    parser.add_argument("--mafft", required=True)
    args = parser.parse_args()
    summary = run(args.dataset, args.catalogue, args.primary_pairs, args.core_pairs,
                  args.positional_pairs, args.output, args.mafft, args.min_loci)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
