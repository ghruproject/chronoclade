"""Matched ST258 neighbourhoods and repeated context-removal checks.

Descriptive research analysis: no epidemiological cutoff or clock is fitted.
"""

import argparse
from collections import Counter
import csv
from datetime import date
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from chronoclade.cgmlst.matrix_thresholds import minimum_spanning_forest
from chronoclade.cgmlst.thresholds import threshold_sweep
from chronoclade.datasets import load_dataset
from chronoclade.metadata_dates import sample_date_interval
from chronoclade.sample_labels import sample_label
from chronoclade.cgmlst.locus_alignment import file_hash
from chronoclade.location_network.reconstruction import observed_host

CUTS = {
    "Allele profiles": [0, 1, 2, 4, 7, 10, 15, 20, 30, 50],
    "Aligned locus SNPs": [0, 1, 2, 4, 7, 10, 15, 20, 30, 50, 100, 200],
    "Matched-locus alleles": [0, 1, 2, 4, 7, 10, 15, 20, 30, 50],
}
SEEDS = (17, 31, 43, 59, 71)
FRACTIONS = (0.25, 0.5, 0.75)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    if rows:
        with Path(path).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def display_names(samples):
    names = {row["sample_id"]: sample_label(row) for row in samples}
    missing = 0
    for ident in names:
        if names[ident].startswith("PW_"):
            missing += 1
            names[ident] = f"Genome without accession {missing}"
    return names


def dated(row):
    value = sample_date_interval(row)
    if value["status"] == "valid" and date.fromisoformat(value["start"]).year >= 1900:
        return value
    return None


def describe_group(members, metadata, inputs):
    rows = [metadata[ident] for ident in members]
    dates = [value for row in rows if (value := dated(row))]
    return {
        "genomes": len(members),
        "inputs": sum(ident in inputs for ident in members),
        "contexts": sum(ident not in inputs for ident in members),
        "dated": len(dates),
        "earliest": min((v["start"] for v in dates), default=None),
        "latest": max((v["end"] for v in dates), default=None),
        "countries": dict(Counter(row.get("country") or "Unknown" for row in rows)),
        "hosts": dict(Counter(observed_host(row.get("host")) or "Unknown" for row in rows)),
        "date_precisions": dict(Counter(v["precision"] for v in dates)),
    }


def partition_maps(ids, forest, cuts):
    sweep = threshold_sweep(ids, forest, cuts)
    result = {}
    for partition in sweep["partitions"]:
        result[partition["threshold"]] = [frozenset(members) for members in partition["clusters"]]
    return result, sweep


def internal_distances(members, index, array, mask):
    positions = np.array([index[ident] for ident in members], dtype=int)
    largest = 0
    available = 0
    for start in range(0, len(positions), 128):
        rows = positions[start : start + 128]
        values = array[np.ix_(rows, positions)]
        valid = mask[np.ix_(rows, positions)]
        if np.any(valid):
            largest = max(largest, int(values[valid].max()))
        available += int(valid.sum())
    diagonal = sum(bool(mask[i, i]) for i in positions)
    expected = len(positions) * (len(positions) - 1) // 2
    return {
        "largest_available_internal_distance": largest,
        "unavailable_internal_pairs": expected - (available - diagonal) // 2,
    }


def nearest_sets(array, mask, inputs, contexts):
    result = {}
    for index in inputs:
        eligible = contexts[mask[index, contexts]]
        if not len(eligible):
            result[int(index)] = {"distance": None, "ids": []}
            continue
        values = array[index, eligible]
        minimum = int(values.min())
        result[int(index)] = {
            "distance": minimum,
            "ids": [int(i) for i in eligible[values == minimum]],
        }
    return result


def reordered_array(source, indexes, output):
    if output.exists():
        raise ValueError(f"Output already exists: {output}")
    result = np.lib.format.open_memmap(output, mode="w+", dtype=source.dtype, shape=source.shape)
    for start in range(0, len(indexes), 128):
        stop = min(start + 128, len(indexes))
        result[start:stop] = source[np.ix_(indexes[start:stop], indexes)]
    result.flush()
    return result


def run(root, dataset_path, cache_path, output):
    started = time.monotonic()
    dataset = load_dataset(dataset_path)
    ids = list(dataset.profiles[0].sample_ids)
    metadata = {row["sample_id"]: row for row in dataset.samples}
    names = display_names(dataset.samples)
    inputs = {ident for ident, row in metadata.items() if row.get("role") == "input"}
    if len(ids) != 7815 or len(inputs) != 41:
        raise ValueError("Expected the declared7815-genome41-input panel")
    input_indices = np.array([i for i, ident in enumerate(ids) if ident in inputs])
    context_indices = np.array([i for i, ident in enumerate(ids) if ident not in inputs])
    descriptor = json.loads(cache_path.read_text())
    for field in ("allele_differences", "shared_called_loci"):
        source = cache_path.parent / descriptor["paths"][field]["path"]
        if file_hash(source) != descriptor["paths"][field]["sha256"]:
            raise ValueError(f"Native allele evidence checksum changed: {field}")
    cache_index = {ident: i for i, ident in enumerate(descriptor["sample_ids"])}
    if set(cache_index) != set(ids):
        raise ValueError("Allele cache has different identities")
    source_raw = np.load(
        cache_path.parent / descriptor["paths"]["allele_differences"]["path"], mmap_mode="r"
    )
    source_calls = np.load(
        cache_path.parent / descriptor["paths"]["shared_called_loci"]["path"], mmap_mode="r"
    )
    output.mkdir(parents=True, exist_ok=False)
    indexes = np.array([cache_index[ident] for ident in ids])
    allele = reordered_array(source_raw, indexes, output / "allele-counts.npy")
    calls = reordered_array(source_calls, indexes, output / "allele-joint-loci.npy")
    directory = root / "nucleotide-matrices"
    descriptor_dna = json.loads((directory / "matrix.json").read_text())
    if descriptor_dna["sample_ids"] != ids:
        raise ValueError("DNA matrix order differs")
    for field, spec in descriptor_dna["arrays"].items():
        if file_hash(directory / spec["path"]) != spec["sha256"]:
            raise ValueError(f"Aligned nucleotide evidence checksum changed: {field}")
    dna = np.load(directory / "snps.npy", mmap_mode="r")
    sites = np.load(directory / "sites.npy", mmap_mode="r")
    loci = np.load(directory / "joint_loci.npy", mmap_mode="r")
    gap = np.load(directory / "gap_columns.npy", mmap_mode="r")
    matched = np.load(directory / "resolved_alleles.npy", mmap_mode="r")
    mask = np.lib.format.open_memmap(
        output / "matched-eligible.npy", mode="w+", dtype="bool", shape=allele.shape
    )
    for start in range(0, len(ids), 128):
        end = min(start + 128, len(ids))
        mask[start:end] = (
            (calls[start:end] >= 597) & (loci[start:end] >= 597) & (sites[start:end] > 0)
        )
    mask.flush()
    print(f"Matched eligible pairs: {(int(mask.sum()) - int(np.trace(mask))) // 2}", flush=True)
    methods = {
        "Allele profiles": allele,
        "Aligned locus SNPs": dna,
        "Matched-locus alleles": matched,
    }
    baselines = {}
    summaries = []
    groups = []
    all_nearest = {}
    nearest_rows = []
    index_by_id = {ident: i for i, ident in enumerate(ids)}
    for method, array in methods.items():
        diameter_cache = {}
        forest = minimum_spanning_forest(ids, array, eligible=mask)
        partitions, sweep = partition_maps(ids, forest, CUTS[method])
        baselines[method] = partitions
        save(output / (method.lower().replace(" ", "-") + "-sweep.json"), sweep)
        all_nearest[method] = nearest_sets(array, mask, input_indices, context_indices)
        focal_values = array[np.ix_(input_indices, context_indices)]
        focal_mask = mask[np.ix_(input_indices, context_indices)]
        for cutoff, components in partitions.items():
            focal = [members for members in components if members & inputs]
            linked = set().union(*focal) - inputs
            summaries.append(
                {
                    "method": method,
                    "cutoff": cutoff,
                    "all_components": len(components),
                    "input_neighbourhoods": len(focal),
                    "linked_contexts": len(linked),
                    "direct_contexts": int(
                        np.count_nonzero(np.any((focal_values <= cutoff) & focal_mask, axis=0))
                    ),
                    "largest_component": max(map(len, components)),
                    "inputs_without_context": sum(
                        len(members - inputs) == 0 for members in focal for _ in members & inputs
                    ),
                }
            )
            for number, members in enumerate(
                sorted(focal, key=lambda group: (-len(group), sorted(group))), 1
            ):
                values = describe_group(members, metadata, inputs)
                if members not in diameter_cache:
                    diameter_cache[members] = internal_distances(members, index_by_id, array, mask)
                groups.append(
                    {
                        "method": method,
                        "cutoff": cutoff,
                        "group": number,
                        "input_accessions": "; ".join(
                            names[ident] for ident in sorted(members & inputs)
                        ),
                        **values,
                        **diameter_cache[members],
                    }
                )
        print(f"Baseline complete: {method}; elapsed {time.monotonic() - started:.1f}s", flush=True)
    nearest_comparison = []
    for index in input_indices:
        ident = ids[index]
        record = {"input": names[ident], "input_sample_id": ident}
        for method, values in all_nearest.items():
            item = values[int(index)]
            record[method] = {
                "distance": item["distance"],
                "ties": len(item["ids"]),
                "labels": [names[ids[j]] for j in item["ids"]],
                "countries": dict(
                    Counter(metadata[ids[j]].get("country") or "Unknown" for j in item["ids"])
                ),
            }
            for j in item["ids"]:
                nearest_rows.append(
                    {
                        "input": names[ident],
                        "method": method,
                        "comparison": names[ids[j]],
                        "distance": item["distance"],
                        "country": metadata[ids[j]].get("country") or "Unknown",
                        "collection_date": metadata[ids[j]].get("collection_date") or "",
                        "joint_DNA_loci": int(loci[index, j]),
                        "joint_allele_loci": int(calls[index, j]),
                        "joint_DNA_sites": int(sites[index, j]),
                        "gap_vs_base_columns": int(gap[index, j]),
                        "input_sample_id": ident,
                        "context_sample_id": ids[j],
                    }
                )
        a = set(all_nearest["Allele profiles"][int(index)]["ids"])
        b = set(all_nearest["Aligned locus SNPs"][int(index)]["ids"])
        c = set(all_nearest["Matched-locus alleles"][int(index)]["ids"])
        record["allele_snp_nearest_overlap"] = len(a & b)
        record["allele_snp_nearest_union"] = len(a | b)
        record["same_nearest_set"] = bool(a and b) and a == b
        record["matched_allele_snp_same_set"] = bool(c and b) and c == b
        nearest_comparison.append(record)
    save(output / "nearest-comparison.json", nearest_comparison)
    write_csv(output / "nearest-neighbours-audit.csv", nearest_rows)
    sampling = []
    sampling_members = []
    for fraction in FRACTIONS:
        for seed in SEEDS:
            chosen = np.random.default_rng(seed).choice(
                context_indices, round(len(context_indices) * fraction), replace=False
            )
            retained = np.sort(np.r_[input_indices, chosen])
            retained_ids = {ids[i] for i in retained}
            sampling_members.append(
                {"fraction": fraction, "seed": seed, "retained_indices": [int(i) for i in retained]}
            )
            for method, array in methods.items():
                forest = minimum_spanning_forest(ids, array, eligible=mask, sample_indices=retained)
                partitions, _ = partition_maps([ids[i] for i in retained], forest, CUTS[method])
                for cutoff, components in partitions.items():
                    sampled = {
                        ident: members
                        for members in components
                        if members & inputs
                        for ident in members & inputs
                    }
                    full = {
                        ident: members
                        for members in baselines[method][cutoff]
                        if members & inputs
                        for ident in members & inputs
                    }
                    recoveries = []
                    split = 0
                    positive_input_pairs = retained_input_pairs = 0
                    for ident in sorted(inputs):
                        expected = (full[ident] & retained_ids) - {ident}
                        actual = sampled[ident] - {ident}
                        if not actual <= expected:
                            raise ValueError(
                                "Induced-subset clustering created an impossible merge"
                            )
                        if expected:
                            recoveries.append(len(actual) / len(expected))
                        split += actual != expected
                        positive_input_pairs += len((full[ident] & inputs) - {ident})
                        retained_input_pairs += len((sampled[ident] & inputs) - {ident})
                    sampling.append(
                        {
                            "method": method,
                            "cutoff": cutoff,
                            "context_fraction": fraction,
                            "seed": seed,
                            "contexts": len(chosen),
                            "inputs": len(input_indices),
                            "inputs_split": split,
                            "peer_recovery_mean": float(np.mean(recoveries))
                            if recoveries
                            else None,
                            "peer_recovery_min": min(recoveries, default=None),
                            "baseline_input_pairs": positive_input_pairs // 2,
                            "retained_input_pairs": retained_input_pairs // 2,
                        }
                    )
            print(
                f"Sampling {fraction:.0%}, seed{seed} complete; elapsed {time.monotonic() - started:.1f}s",
                flush=True,
            )
    save(output / "sampling-members-audit.json", {"sample_ids": ids, "samples": sampling_members})
    write_csv(output / "sampling.csv", sampling)
    write_csv(output / "threshold-summary.csv", summaries)
    save(output / "groups.json", groups)
    write_csv(
        output / "groups.csv",
        [
            {
                **row,
                "countries": json.dumps(row["countries"]),
                "hosts": json.dumps(row["hosts"]),
                "date_precisions": json.dumps(row["date_precisions"]),
            }
            for row in groups
        ],
    )
    country = Counter((row.get("country") or "Unknown", row.get("role")) for row in dataset.samples)
    countries = [
        {"country": name, "inputs": country[name, "input"], "contexts": country[name, "context"]}
        for name in sorted(
            {key[0] for key in country},
            key=lambda n: (-country[n, "input"] - country[n, "context"], n),
        )
    ]
    years = []
    date_quality = Counter()
    for row in dataset.samples:
        interval = sample_date_interval(row)
        valid = dated(row)
        date_quality[
            (row.get("role"), "eligible" if valid else interval.get("reason") or "before1900")
        ] += 1
        if valid:
            first, last = (
                date.fromisoformat(valid["start"]).year,
                date.fromisoformat(valid["end"]).year,
            )
            if first == last:
                years.append(
                    {
                        "year": first,
                        "role": row.get("role"),
                        "country": row.get("country") or "Unknown",
                    }
                )
    save(
        output / "dataset-summary.json",
        {
            "samples": len(ids),
            "inputs": len(inputs),
            "contexts": len(context_indices),
            "countries": countries,
            "date_quality": [
                {"role": role, "status": status, "count": count}
                for (role, status), count in date_quality.items()
            ],
            "years": years,
            "host_reported": sum(
                observed_host(row.get("host")) is not None for row in dataset.samples
            ),
            "host_counts": dict(
                Counter(observed_host(row.get("host")) or "Unknown" for row in dataset.samples)
            ),
            "eligible_pairs": (int(mask.sum()) - int(np.trace(mask))) // 2,
            "total_pairs": len(ids) * (len(ids) - 1) // 2,
            "minimum_joint_loci": 597,
            "coverage_universe": 628,
            "unusable_inputs": [names[ids[i]] for i in input_indices if not mask[i, i]],
        },
    )
    save(
        output / "receipt.json",
        {
            "dataset_id": dataset.dataset_id,
            "dataset_manifest": str(dataset_path),
            "source_allele_evidence": str(cache_path),
            "nucleotide_evidence": str(directory / "matrix.json"),
            "source_sha256": {
                str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in [dataset_path, cache_path, directory / "matrix.json", Path(__file__)]
            },
            "seeds": SEEDS,
            "context_fractions": FRACTIONS,
            "cutoffs": CUTS,
            "elapsed_seconds": time.monotonic() - started,
            "policy": "same common pair eligibility; induced-subset forests recomputed; no biological cutoff selected",
        },
    )
    print(f"Analysis complete in {time.monotonic() - started:.1f}s", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("tmp/science-st258"))
    parser.add_argument(
        "--dataset", type=Path, default=Path("tmp/st258-native-e2e/prepared/dataset.json")
    )
    parser.add_argument(
        "--allele-cache",
        type=Path,
        default=Path(
            "tmp/st258-native-e2e/cgmlst/cglin_37d98c4c5fe0/binary_distances/distance_evidence.json"
        ),
    )
    parser.add_argument("--output", type=Path, default=Path("tmp/science-st258/comparison"))
    args = parser.parse_args()
    run(args.root, args.dataset, args.allele_cache, args.output)
