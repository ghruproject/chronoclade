"""Cheap lineage-coherence checks before phylogenetic optimisation."""

from __future__ import annotations

import csv
import math
import statistics
from pathlib import Path

from chronoclade.evidence import EvidenceError, pairwise_distances, read_alignment
from chronoclade.metadata import Sample


MIN_ROBUST_Z = 10.0
MIN_MEDIAN_RATIO = 5.0


def screen_alignment(
    alignment: Path,
    samples: list[Sample],
    *,
    output: Path,
) -> dict[str, object]:
    """Identify isolated, extreme raw-distance outliers in a mapped alignment.

    This is deliberately a relative cohort check rather than a bacterial SNP
    cutoff. It is intended to catch accession, species or lineage mistakes,
    not to define transmission clusters.
    """

    sequences = read_alignment(alignment)
    pairwise = pairwise_distances(sequences, samples)
    by_sample: dict[str, list[tuple[int, float]]] = {sample.sample_id: [] for sample in samples}
    for record in pairwise:
        proportion = record["snp_proportion"]
        if proportion is None:
            continue
        value = (int(record["clonal_snps"]), float(proportion))
        by_sample[str(record["sample_1"])].append(value)
        by_sample[str(record["sample_2"])].append(value)

    medians = {
        sample_id: statistics.median(value[1] for value in values)
        for sample_id, values in by_sample.items()
        if values
    }
    if len(medians) != len(samples):
        missing = sorted(set(by_sample) - set(medians))
        raise EvidenceError(
            "Lineage-coherence screen found no callable comparisons for: " + ", ".join(missing)
        )

    centre = statistics.median(medians.values())
    mad = statistics.median(abs(value - centre) for value in medians.values())
    rows: list[dict[str, object]] = []
    flagged: list[str] = []
    for sample_id in sorted(medians):
        value = medians[sample_id]
        ratio = value / centre if centre > 0 else (math.inf if value > 0 else 1.0)
        robust_z = (
            0.6745 * (value - centre) / mad if mad > 0 else (math.inf if value > centre else 0.0)
        )
        is_outlier = robust_z >= MIN_ROBUST_Z and ratio >= MIN_MEDIAN_RATIO
        if is_outlier:
            flagged.append(sample_id)
        rows.append(
            {
                "sample_id": sample_id,
                "median_raw_snps": statistics.median(pair[0] for pair in by_sample[sample_id]),
                "median_snp_proportion": f"{value:.8g}",
                "cohort_median_snp_proportion": f"{centre:.8g}",
                "ratio_to_cohort_median": f"{ratio:.4g}",
                "robust_z": f"{robust_z:.4g}",
                "flagged": str(is_outlier).lower(),
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return {
        "status": "failed" if flagged else "passed",
        "sample_count": len(samples),
        "flagged_samples": flagged,
        "cohort_median_snp_proportion": centre,
        "median_absolute_deviation": mad,
        "output": str(output),
    }
