"""Derive transparent public-health evidence from recombination-filtered outputs."""

from __future__ import annotations

import csv
import json
import statistics
from dataclasses import dataclass
from html import escape
from pathlib import Path

import numpy as np

from beyondmlst.metadata import Sample


class EvidenceError(RuntimeError):
    """Raised when actionable evidence cannot be derived safely."""


@dataclass(frozen=True)
class _EncodedSequence:
    valid: int
    low_bit: int
    high_bit: int


@dataclass
class _TreeNode:
    name: str = ""
    children: tuple[_TreeNode, ...] = ()

    @property
    def tips(self) -> set[str]:
        if not self.children:
            return {self.name}
        return set().union(*(child.tips for child in self.children))


class _NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.position = 0

    def parse(self) -> _TreeNode:
        node = self._subtree()
        self._skip_space()
        if self.position < len(self.text) and self.text[self.position] == ";":
            self.position += 1
        self._skip_space()
        if self.position != len(self.text):
            raise EvidenceError("Could not parse the recombination-corrected Newick tree")
        return node

    def _subtree(self) -> _TreeNode:
        self._skip_space()
        children: list[_TreeNode] = []
        if self._peek() == "(":
            self.position += 1
            while True:
                children.append(self._subtree())
                self._skip_space()
                token = self._peek()
                if token == ",":
                    self.position += 1
                    continue
                if token == ")":
                    self.position += 1
                    break
                raise EvidenceError("Could not parse the recombination-corrected Newick tree")
        name = self._label()
        self._skip_branch_length()
        return _TreeNode(name=name, children=tuple(children))

    def _label(self) -> str:
        self._skip_space()
        if self._peek() in {"'", '"'}:
            quote = self._peek()
            self.position += 1
            start = self.position
            while self.position < len(self.text) and self.text[self.position] != quote:
                self.position += 1
            value = self.text[start : self.position]
            if self.position < len(self.text):
                self.position += 1
            return value
        start = self.position
        while self.position < len(self.text) and self.text[self.position] not in ":,();":
            self.position += 1
        return self.text[start : self.position].strip()

    def _skip_branch_length(self) -> None:
        self._skip_space()
        if self._peek() != ":":
            return
        self.position += 1
        while self.position < len(self.text) and self.text[self.position] not in ",();":
            self.position += 1

    def _skip_space(self) -> None:
        while self.position < len(self.text) and self.text[self.position].isspace():
            self.position += 1

    def _peek(self) -> str:
        return self.text[self.position] if self.position < len(self.text) else ""


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


def _summary(values: list[int]) -> dict[str, int | float | None]:
    if not values:
        return {"comparisons": 0, "minimum": None, "median": None, "maximum": None}
    return {
        "comparisons": len(values),
        "minimum": min(values),
        "median": statistics.median(values),
        "maximum": max(values),
    }


def summarise_distances(records: list[dict[str, object]]) -> dict[str, object]:
    """Summarise the exact comparisons used by the public-health report."""

    categories: dict[str, dict[str, int | float | None]] = {}
    for field, names in (
        ("comparison", ("focal_focal", "focal_context", "context_context")),
        ("period_comparison", ("same_month", "same_year", "between_years")),
        ("patient_comparison", ("same_patient", "different_patients")),
    ):
        for name in names:
            values = [
                int(record["clonal_snps"])
                for record in records
                if record[field] == name
                and (field == "comparison" or record["comparison"] == "focal_focal")
            ]
            categories[name] = _summary(values)
    callable_values = [int(record["callable_sites"]) for record in records]
    return {"categories": categories, "callable_sites": _summary(callable_values)}


def _write_rows(path: Path, rows: list[dict[str, object]]) -> Path:
    if not rows:
        path.write_text("", encoding="utf-8")
        return path
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def _matrix(
    sample_ids: list[str],
    records: list[dict[str, object]],
    field: str,
    *,
    diagonal: dict[str, int] | None = None,
) -> list[list[int]]:
    lookup: dict[frozenset[str], int] = {
        frozenset((str(record["sample_1"]), str(record["sample_2"]))): int(record[field])
        for record in records
    }
    result: list[list[int]] = []
    for first in sample_ids:
        row = []
        for second in sample_ids:
            row.append(
                (diagonal or {}).get(first, 0)
                if first == second
                else lookup[frozenset((first, second))]
            )
        result.append(row)
    return result


def _write_matrix(path: Path, sample_ids: list[str], matrix: list[list[int]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", *sample_ids])
        for sample_id, values in zip(sample_ids, matrix, strict=True):
            writer.writerow([sample_id, *values])
    return path


def write_distance_heatmap(
    path: Path,
    sample_ids: list[str],
    matrix: list[list[int]],
    origins: dict[str, str],
) -> Path:
    """Write a dependency-free SVG clonal-SNP heatmap with exact tooltips."""

    count = len(sample_ids)
    cell = max(4.0, min(18.0, 720.0 / max(count, 1)))
    margin = 210.0 if count <= 60 else 70.0
    plot = count * cell
    width = max(620.0, margin + plot + 90)
    height = margin + plot + 85
    maximum = max((value for row in matrix for value in row), default=0) or 1
    cells: list[str] = []
    for row_index, row in enumerate(matrix):
        for column_index, value in enumerate(row):
            fraction = value / maximum
            red = round(248 - 128 * fraction)
            green = round(250 - 206 * fraction)
            blue = round(252 - 190 * fraction)
            colour = "#e7ecef" if row_index == column_index else f"rgb({red},{green},{blue})"
            first, second = sample_ids[row_index], sample_ids[column_index]
            cells.append(
                f'<rect x="{margin + column_index * cell:.2f}" '
                f'y="{margin + row_index * cell:.2f}" width="{cell:.2f}" '
                f'height="{cell:.2f}" fill="{colour}"><title>'
                f"{escape(first)} vs {escape(second)}: {value} clonal SNPs"
                "</title></rect>"
            )
    labels: list[str] = []
    if count <= 60:
        font_size = max(7.0, min(11.0, cell * 0.72))
        for index, sample_id in enumerate(sample_ids):
            colour = "#8b3f12" if origins.get(sample_id, "").casefold() == "context" else "#17394d"
            coordinate = margin + index * cell + cell * 0.7
            labels.append(
                f'<text x="{margin - 7}" y="{coordinate:.2f}" text-anchor="end" '
                f'font-size="{font_size:.2f}" fill="{colour}">{escape(sample_id)}</text>'
            )
            labels.append(
                f'<text transform="translate({coordinate:.2f},{margin - 7}) rotate(-90)" '
                f'text-anchor="start" font-size="{font_size:.2f}" fill="{colour}">'
                f"{escape(sample_id)}</text>"
            )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}" role="img" aria-labelledby="title description">
  <title id="title">Recombination-filtered pairwise SNP heatmap</title>
  <desc id="description">Pairwise clonal SNP counts. Darker cells contain more differences. Exact values are available in the downloadable matrix.</desc>
  <rect width="100%" height="100%" fill="white"/>
  <text x="20" y="30" font-family="system-ui,sans-serif" font-size="18" font-weight="700" fill="#17394d">Recombination-filtered pairwise SNPs</text>
  <text x="20" y="52" font-family="system-ui,sans-serif" font-size="12" fill="#5b6b76">Focal labels are blue; contextual labels are brown. Scale: 0–{maximum} SNPs.</text>
  {"".join(cells)}
  {"".join(labels)}
</svg>
"""
    path.write_text(svg, encoding="utf-8")
    return path


def _date_key(sample: Sample) -> tuple[str, str]:
    return (_year(sample.collection_date) or "9999", sample.collection_date)


def patient_sensitivity(
    samples: list[Sample], records: list[dict[str, object]]
) -> dict[str, object]:
    """Describe a deterministic one-isolate-per-patient distance sensitivity set."""

    focal = [sample for sample in samples if not _is_context(sample)]
    by_patient: dict[str, list[Sample]] = {}
    retained = [sample for sample in focal if not sample.patient_id]
    for sample in focal:
        if sample.patient_id:
            by_patient.setdefault(sample.patient_id, []).append(sample)
    for patient_samples in by_patient.values():
        retained.append(
            min(patient_samples, key=lambda sample: (*_date_key(sample), sample.sample_id))
        )
    retained_ids = {sample.sample_id for sample in retained}
    excluded = sorted(sample.sample_id for sample in focal if sample.sample_id not in retained_ids)
    full_values = [
        int(record["clonal_snps"]) for record in records if record["comparison"] == "focal_focal"
    ]
    retained_values = [
        int(record["clonal_snps"])
        for record in records
        if record["sample_1"] in retained_ids
        and record["sample_2"] in retained_ids
        and record["comparison"] == "focal_focal"
    ]
    coded = sum(bool(sample.patient_id) for sample in focal)
    return {
        "focal_samples": len(focal),
        "samples_with_patient_id": coded,
        "patient_metadata_complete": bool(focal) and coded == len(focal),
        "representative_samples": sorted(retained_ids),
        "excluded_repeated_patient_samples": excluded,
        "full_focal_distances": _summary(full_values),
        "one_per_patient_distances": _summary(retained_values),
        "interpretation": (
            "Repeated patient isolates were removed deterministically by retaining the earliest "
            "dated isolate per coded patient."
            if excluded
            else "No repeated coded patient isolates were removed."
        ),
    }


def _candidate_groups(node: _TreeNode, focal_ids: set[str]) -> list[set[str]]:
    tips = node.tips
    focal_tips = tips & focal_ids
    if not focal_tips:
        return []
    if tips <= focal_ids:
        return [focal_tips]
    if not node.children:
        return [{node.name}] if node.name in focal_ids else []
    groups: list[set[str]] = []
    for child in node.children:
        groups.extend(_candidate_groups(child, focal_ids))
    return groups


def _record_lookup(records: list[dict[str, object]]) -> dict[frozenset[str], dict[str, object]]:
    return {
        frozenset((str(record["sample_1"]), str(record["sample_2"]))): record for record in records
    }


def topology_groups(
    tree: Path, samples: list[Sample], records: list[dict[str, object]]
) -> dict[str, object]:
    """Describe maximal focal-only groups in the rooted corrected topology."""

    root = _NewickParser(tree.read_text(encoding="utf-8")).parse()
    by_id = {sample.sample_id: sample for sample in samples}
    focal_ids = {sample.sample_id for sample in samples if not _is_context(sample)}
    context_ids = {sample.sample_id for sample in samples if _is_context(sample)}
    if root.tips != set(by_id):
        missing = sorted(set(by_id) - root.tips)
        extra = sorted(root.tips - set(by_id))
        raise EvidenceError(f"Tree/metadata samples differ; missing={missing}, extra={extra}")
    raw_groups = _candidate_groups(root, focal_ids)
    lookup = _record_lookup(records)
    groups: list[dict[str, object]] = []
    for index, sample_ids in enumerate(
        sorted(raw_groups, key=lambda value: sorted(value)), start=1
    ):
        members = [by_id[sample_id] for sample_id in sorted(sample_ids)]
        within = [
            int(lookup[frozenset((first, second))]["clonal_snps"])
            for member_index, first in enumerate(sorted(sample_ids))
            for second in sorted(sample_ids)[member_index + 1 :]
        ]
        contextual = [
            (
                int(lookup[frozenset((focal_id, context_id))]["clonal_snps"]),
                focal_id,
                context_id,
            )
            for focal_id in sample_ids
            for context_id in context_ids
        ]
        nearest_context = min(contextual) if contextual else None
        dates = sorted({member.collection_date for member in members}, key=str)
        groups.append(
            {
                "group_id": f"LG{index}",
                "sample_ids": sorted(sample_ids),
                "sample_count": len(sample_ids),
                "patient_count": len(
                    {member.patient_id for member in members if member.patient_id}
                ),
                "locations": sorted({member.location for member in members}),
                "first_collection_date": dates[0] if dates else "",
                "last_collection_date": dates[-1] if dates else "",
                "distinct_collection_dates": len(dates),
                "spans_sampling_periods": len(dates) >= 2,
                "within_group_clonal_snps": _summary(within),
                "nearest_context": (
                    {
                        "sample_id": nearest_context[2],
                        "nearest_focal": nearest_context[1],
                        "clonal_snps": nearest_context[0],
                    }
                    if nearest_context
                    else None
                ),
            }
        )
    return {
        "rooting": "TreeTime least-squares root",
        "branch_support_available": False,
        "focal_monophyletic": len(groups) == 1,
        "candidate_local_group_count": len(groups),
        "groups": groups,
        "interpretation": (
            "Candidate groups are maximal focal-only clades in the rooted recombination-corrected "
            "tree. They are review aids, not confirmed transmission clusters."
        ),
    }


def add_group_distance_summaries(
    summary: dict[str, object],
    records: list[dict[str, object]],
    topology: dict[str, object],
) -> None:
    """Add within/between candidate-group clonal distances to an existing summary."""

    groups = topology.get("groups", [])
    if not isinstance(groups, list):
        raise EvidenceError("Topology groups must be a list")
    membership: dict[str, str] = {}
    for group in groups:
        if not isinstance(group, dict):
            continue
        group_id = str(group.get("group_id", ""))
        sample_ids = group.get("sample_ids", [])
        if isinstance(sample_ids, list):
            membership.update({str(sample_id): group_id for sample_id in sample_ids})
    within: list[int] = []
    between: list[int] = []
    for record in records:
        if record["comparison"] != "focal_focal":
            continue
        first_group = membership.get(str(record["sample_1"]))
        second_group = membership.get(str(record["sample_2"]))
        if not first_group or not second_group:
            continue
        target = within if first_group == second_group else between
        target.append(int(record["clonal_snps"]))
    categories = summary.get("categories")
    if not isinstance(categories, dict):
        raise EvidenceError("Distance summary categories must be a mapping")
    categories["within_candidate_groups"] = _summary(within)
    categories["between_candidate_groups"] = _summary(between)


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


SCENARIO_LABELS = {
    "persistent_local_lineage": "Consistent with a persistent local lineage",
    "multiple_introductions": "Consistent with multiple introductions",
    "mixed": "Consistent with local persistence plus additional introductions",
    "indeterminate": "Indeterminate",
}


def _plural(count: int, singular: str, plural: str | None = None) -> str:
    return singular if count == 1 else (plural or singular + "s")


def synthesise_scenario(
    samples: list[Sample],
    topology: dict[str, object],
    distance_summary: dict[str, object],
    sensitivity: dict[str, object],
    temporal_assessment: dict[str, object],
) -> dict[str, object]:
    """Create a cautious, auditable four-scenario working interpretation."""

    focal = [sample for sample in samples if not _is_context(sample)]
    contexts = [sample for sample in samples if _is_context(sample)]
    groups = topology["groups"]
    if not isinstance(groups, list):
        raise EvidenceError("Topology groups must be a list")
    spanning = [
        group
        for group in groups
        if isinstance(group, dict)
        and bool(group.get("spans_sampling_periods"))
        and int(group.get("sample_count", 0)) >= 2
    ]
    if len(focal) < 2 or not contexts or len({sample.collection_date for sample in focal}) < 2:
        code = "indeterminate"
    elif len(groups) == 1 and spanning:
        code = "persistent_local_lineage"
    elif len(groups) >= 2 and spanning:
        code = "mixed"
    elif len(groups) >= 2:
        code = "multiple_introductions"
    else:
        code = "indeterminate"

    evidence = [
        {
            "id": "clonal_distances",
            "status": "supported",
            "finding": (
                f"Recombination-filtered SNP and callable-site counts were calculated for "
                f"{len(samples) * (len(samples) - 1) // 2} genome "
                f"{_plural(len(samples) * (len(samples) - 1) // 2, 'pair')}."
            ),
        },
        {
            "id": "longitudinal_span",
            "status": "supported"
            if len({sample.collection_date for sample in focal}) >= 2
            else "unavailable",
            "finding": (
                f"{len(focal)} focal {_plural(len(focal), 'genome')} "
                f"{'covers' if len(focal) == 1 else 'cover'} "
                f"{len({sample.collection_date for sample in focal})} distinct collection "
                f"{_plural(len({sample.collection_date for sample in focal}), 'date')}."
            ),
        },
        {
            "id": "context",
            "status": "suggestive" if contexts else "unavailable",
            "finding": (
                f"{len(contexts)} public context {_plural(len(contexts), 'genome')} "
                f"{'was' if len(contexts) == 1 else 'were'} analysed; neighbour selection was "
                "bounded and cannot guarantee the globally nearest public genome."
                if contexts
                else "No public context genomes were analysed."
            ),
        },
        {
            "id": "topology",
            "status": "suggestive",
            "finding": (
                f"The rooted corrected tree contains {len(groups)} maximal focal-only candidate "
                f"{_plural(len(groups), 'group')}; branch support is not yet propagated through "
                "ClonalFrameML."
            ),
        },
        {
            "id": "patient_sensitivity",
            "status": (
                "supported" if bool(sensitivity["patient_metadata_complete"]) else "unavailable"
            ),
            "finding": (
                str(sensitivity["interpretation"])
                if bool(sensitivity["patient_metadata_complete"])
                else "Patient identifiers are incomplete, so patient-level sensitivity is limited."
            ),
        },
        {
            "id": "temporal_signal",
            "status": "supported" if bool(temporal_assessment["supported"]) else "unavailable",
            "finding": str(temporal_assessment["reason"]),
        },
    ]
    reasons = [
        f"The local isolates form {len(groups)} separate {_plural(len(groups), 'group')} in "
        "the recombination-corrected tree.",
        f"{len(spanning)} local {_plural(len(spanning), 'group')} "
        f"{'contains' if len(spanning) == 1 else 'contain'} isolates collected on more than "
        "one date.",
        f"The comparison includes {len(contexts)} public context "
        f"{_plural(len(contexts), 'genome')}.",
    ]
    categories = distance_summary.get("categories", {})
    if isinstance(categories, dict):
        within_groups = categories.get("within_candidate_groups", {})
        between_groups = categories.get("between_candidate_groups", {})
        if (
            isinstance(within_groups, dict)
            and isinstance(between_groups, dict)
            and within_groups.get("comparisons")
            and between_groups.get("comparisons")
        ):
            reasons.append(
                f"Median clonal distance was {float(within_groups['median']):.3g} SNPs within "
                f"candidate groups and {float(between_groups['median']):.3g} SNPs between them."
            )
    actions: dict[str, list[str]] = {
        "persistent_local_lineage": [
            "Review whether isolates in the persistent group overlap by patient, ward, facility or exposure.",
            "Continue targeted sampling to determine whether the lineage remains locally detectable.",
        ],
        "multiple_introductions": [
            "Review referral, travel, community and source histories separately for each candidate local group.",
            "Expand public or regional context around each group before estimating introduction pathways.",
        ],
        "mixed": [
            "Investigate the longitudinal candidate group for local persistence and the separated groups for external acquisition pathways.",
            "Prioritise missing patient, facility and exposure metadata that could distinguish these explanations.",
        ],
        "indeterminate": [
            "Add longitudinal focal isolates and public context, or improve their metadata, before assigning a circulation scenario.",
            "Review the exact clonal distances and corrected topology rather than applying a universal SNP cutoff.",
        ],
    }
    rules = {
        "persistent_local_lineage": (
            "One focal-only topology group spans at least two sampling dates and public context "
            "is present."
        ),
        "multiple_introductions": (
            "At least two separated focal-only topology groups are present, with no group "
            "showing longitudinal persistence in the sampled dates."
        ),
        "mixed": (
            "At least two separated focal-only topology groups are present and at least one "
            "contains multiple focal genomes spanning sampling dates."
        ),
        "indeterminate": (
            "The minimum focal, longitudinal or contextual evidence gate was not met."
        ),
    }
    if not bool(sensitivity["patient_metadata_complete"]):
        actions[code].append(
            "Complete coded patient identifiers to assess repeated-patient sampling."
        )
    return {
        "code": code,
        "label": SCENARIO_LABELS[code],
        "confidence": "moderate" if code != "indeterminate" and len(contexts) >= 3 else "low",
        "provisional": True,
        "decision_rule": rules[code],
        "reasons": reasons,
        "recommended_follow_up": actions[code],
        "evidence": evidence,
        "guardrail": (
            "This is a working genomic interpretation, not a count of introduction events or "
            "evidence of direct transmission. Epidemiological review is required."
        ),
        "distance_summary": distance_summary,
    }


def build_public_health_evidence(
    *,
    filtered_alignment: Path,
    rooted_tree: Path,
    samples: list[Sample],
    output: Path,
    temporal_assessment: dict[str, object],
) -> dict[str, object]:
    """Build all machine-readable and visual evidence for one lineage report."""

    sequences = read_alignment(filtered_alignment)
    records = pairwise_distances(sequences, samples)
    sample_order = [
        sample.sample_id
        for sample in sorted(
            samples, key=lambda item: (_is_context(item), _date_key(item), item.sample_id)
        )
    ]
    snp_matrix = _matrix(sample_order, records, "clonal_snps")
    callable_matrix = _matrix(
        sample_order,
        records,
        "callable_sites",
        diagonal={
            sample_id: sum(sequences[sample_id].count(base) for base in "ACGT")
            for sample_id in sample_order
        },
    )
    pairwise_path = _write_rows(output / "clonal_pairwise_distances.tsv", records)
    snp_path = _write_matrix(output / "clonal_snp_matrix.tsv", sample_order, snp_matrix)
    callable_path = _write_matrix(
        output / "pairwise_callable_sites.tsv", sample_order, callable_matrix
    )
    heatmap_path = write_distance_heatmap(
        output / "clonal_snp_heatmap.svg",
        sample_order,
        snp_matrix,
        {sample.sample_id: sample.origin for sample in samples},
    )
    summary = summarise_distances(records)
    sensitivity = patient_sensitivity(samples, records)
    topology = topology_groups(rooted_tree, samples, records)
    add_group_distance_summaries(summary, records, topology)
    nearest = nearest_contexts(samples, records)
    scenario = synthesise_scenario(
        samples,
        topology,
        summary,
        sensitivity,
        temporal_assessment,
    )
    result = {
        "scenario": scenario,
        "distance_summary": summary,
        "topology": topology,
        "patient_sensitivity": sensitivity,
        "nearest_contexts_by_clonal_distance": nearest,
        "outputs": {
            "pairwise_distances": str(pairwise_path),
            "snp_matrix": str(snp_path),
            "callable_sites_matrix": str(callable_path),
            "heatmap": str(heatmap_path),
        },
    }
    (output / "public_health_evidence.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result
