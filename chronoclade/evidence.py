"""Read alignments and compare recombination-filtered SNP evidence."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from chronoclade.metadata import Sample


class EvidenceError(RuntimeError):
    """Raised when actionable evidence cannot be derived safely."""


@dataclass(frozen=True)
class _EncodedSequence:
    valid: int
    low_bit: int
    high_bit: int


def read_alignment(path: Path) -> dict[str, str]:
    """Read one equal-length FASTA alignment into memory."""

    sequences: dict[str, list[str]] = {}
    current = ""
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current or current in sequences:
                    raise EvidenceError(f"Invalid or duplicate FASTA identifier in {path}")
                sequences[current] = []
            elif current:
                sequences[current].append(line.upper())
            else:
                raise EvidenceError(f"Alignment sequence appears before its header: {path}")
    result = {name: "".join(parts) for name, parts in sequences.items()}
    lengths = {len(sequence) for sequence in result.values()}
    if not result or len(lengths) != 1:
        raise EvidenceError(f"Alignment is empty or contains unequal sequence lengths: {path}")
    return result


def _is_context(sample: Sample) -> bool:
    return sample.origin.strip().casefold() == "context"


def _year(value: str) -> str:
    return value[:4] if len(value) >= 4 and value[:4].isdigit() else ""


def _month(value: str) -> str:
    if len(value) < 7 or value[4] != "-" or not value[:4].isdigit():
        return ""
    month = value[5:7]
    return value[:7] if month.isdigit() and 1 <= int(month) <= 12 else ""


def _period_comparison(first: Sample, second: Sample) -> str:
    first_month, second_month = _month(first.collection_date), _month(second.collection_date)
    if first_month and first_month == second_month:
        return "same_month"
    first_year, second_year = _year(first.collection_date), _year(second.collection_date)
    if first_year and first_year == second_year:
        return "same_year"
    if first_year and second_year:
        return "between_years"
    return "date_unavailable"


def _patient_comparison(first: Sample, second: Sample) -> str:
    if not first.patient_id or not second.patient_id:
        return "patient_unavailable"
    return "same_patient" if first.patient_id == second.patient_id else "different_patients"


def pairwise_distances(sequences: dict[str, str], samples: list[Sample]) -> list[dict[str, object]]:
    """Count clonal SNPs and pair-specific callable sites for every pair."""

    missing = [sample.sample_id for sample in samples if sample.sample_id not in sequences]
    if missing:
        raise EvidenceError("Filtered alignment is missing samples: " + ", ".join(missing))
    encoded = _encode_sequences(sequences)
    records: list[dict[str, object]] = []
    for index, first in enumerate(samples):
        first_sequence = encoded[first.sample_id]
        for second in samples[index + 1 :]:
            second_sequence = encoded[second.sample_id]
            callable_mask = first_sequence.valid & second_sequence.valid
            difference_mask = callable_mask & (
                (first_sequence.low_bit ^ second_sequence.low_bit)
                | (first_sequence.high_bit ^ second_sequence.high_bit)
            )
            callable_sites = callable_mask.bit_count()
            snps = difference_mask.bit_count()
            first_context, second_context = _is_context(first), _is_context(second)
            if first_context and second_context:
                comparison = "context_context"
            elif first_context or second_context:
                comparison = "focal_context"
            else:
                comparison = "focal_focal"
            records.append(
                {
                    "sample_1": first.sample_id,
                    "sample_2": second.sample_id,
                    "origin_1": first.origin,
                    "origin_2": second.origin,
                    "date_1": first.collection_date,
                    "date_2": second.collection_date,
                    "patient_1": first.patient_id,
                    "patient_2": second.patient_id,
                    "comparison": comparison,
                    "period_comparison": _period_comparison(first, second),
                    "patient_comparison": _patient_comparison(first, second),
                    "clonal_snps": snps,
                    "callable_sites": callable_sites,
                    "snp_proportion": snps / callable_sites if callable_sites else None,
                }
            )
    return records


def _encode_sequences(sequences: dict[str, str]) -> dict[str, _EncodedSequence]:
    """Pack A/C/G/T validity and two-bit base codes into Python integers."""

    lookup = np.full(256, 4, dtype=np.uint8)
    for base, code in ((ord("A"), 0), (ord("C"), 1), (ord("G"), 2), (ord("T"), 3)):
        lookup[base] = code
    result: dict[str, _EncodedSequence] = {}
    for sample_id, sequence in sequences.items():
        raw = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
        codes = lookup[raw]
        valid = codes < 4
        low_bit = valid & ((codes & 1) == 1)
        high_bit = valid & ((codes & 2) == 2)
        result[sample_id] = _EncodedSequence(
            valid=int.from_bytes(np.packbits(valid, bitorder="little").tobytes(), "little"),
            low_bit=int.from_bytes(np.packbits(low_bit, bitorder="little").tobytes(), "little"),
            high_bit=int.from_bytes(np.packbits(high_bit, bitorder="little").tobytes(), "little"),
        )
    return result


def nearest_contexts(
    samples: list[Sample], records: list[dict[str, object]]
) -> list[dict[str, object]]:
    """Return each focal sample's nearest context by final clonal SNP distance."""

    focal_ids = {sample.sample_id for sample in samples if not _is_context(sample)}
    context_ids = {sample.sample_id for sample in samples if _is_context(sample)}
    result: list[dict[str, object]] = []
    for focal_id in sorted(focal_ids):
        comparisons = []
        for record in records:
            pair = {str(record["sample_1"]), str(record["sample_2"])}
            if focal_id not in pair:
                continue
            context = next((sample_id for sample_id in pair if sample_id in context_ids), None)
            if context:
                comparisons.append(
                    (
                        int(record["clonal_snps"]),
                        -int(record["callable_sites"]),
                        context,
                        record,
                    )
                )
        if comparisons:
            snps, _, context_id, record = min(comparisons)
            result.append(
                {
                    "focal_sample": focal_id,
                    "context_sample": context_id,
                    "clonal_snps": snps,
                    "callable_sites": int(record["callable_sites"]),
                    "snp_proportion": record["snp_proportion"],
                }
            )
    return result
