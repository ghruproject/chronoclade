"""Reader-facing evidence for ClonalFrameML recombination filtering."""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


class RecombinationReportError(ValueError):
    """Raised when ClonalFrameML outputs cannot be reconciled."""


@dataclass(frozen=True)
class ImportationInterval:
    """One 1-based, closed ClonalFrameML importation interval on a tree branch."""

    node: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def read_importation_intervals(path: Path) -> list[ImportationInterval]:
    """Read the standard ``Node Beg End`` ClonalFrameML interval file."""

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["Node", "Beg", "End"]:
            raise RecombinationReportError(
                f"Expected ClonalFrameML columns Node, Beg, End in {path}"
            )
        intervals: list[ImportationInterval] = []
        for row_number, row in enumerate(reader, start=2):
            try:
                interval = ImportationInterval(
                    node=str(row["Node"]).strip(),
                    start=int(row["Beg"]),
                    end=int(row["End"]),
                )
            except (KeyError, TypeError, ValueError) as error:
                raise RecombinationReportError(
                    f"Could not parse importation interval on row {row_number} of {path}"
                ) from error
            if not interval.node or interval.start < 1 or interval.end < interval.start:
                raise RecombinationReportError(
                    f"Invalid importation interval on row {row_number} of {path}"
                )
            intervals.append(interval)
    return intervals


def _alignment_incomplete_mask(path: Path) -> tuple[np.ndarray, int]:
    """Return sites with an ambiguous base in any sequence and the sequence count."""

    mask: np.ndarray | None = None
    position = 0
    sequences = 0
    with path.open("rb") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(b">"):
                if mask is not None and position != mask.size:
                    raise RecombinationReportError(f"Unequal sequence lengths in {path}")
                sequences += 1
                position = 0
                continue
            if not sequences:
                raise RecombinationReportError(f"Sequence data precede the FASTA header in {path}")
            bases = np.frombuffer(line.upper(), dtype=np.uint8)
            if mask is None:
                mask = np.zeros(len(bases), dtype=np.bool_)
            elif position + len(bases) > mask.size:
                if sequences == 1:
                    mask.resize(position + len(bases), refcheck=False)
                else:
                    raise RecombinationReportError(f"Unequal sequence lengths in {path}")
            mask[position : position + len(bases)] |= ~np.isin(bases, (65, 67, 71, 84))
            position += len(bases)
    if mask is None or not sequences:
        raise RecombinationReportError(f"Alignment contains no sequences: {path}")
    if position != mask.size:
        raise RecombinationReportError(f"Unequal sequence lengths in {path}")
    return mask, sequences


def _fasta_alignment_length(path: Path) -> int:
    length = 0
    seen = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen:
                    break
                seen = True
            elif seen:
                length += len(line)
    if not seen:
        raise RecombinationReportError(f"Alignment contains no sequences: {path}")
    return length


def _reference_records(path: Path | None, alignment_length: int) -> list[tuple[str, int, int]]:
    """Return reference-record names and 1-based cumulative alignment bounds."""

    if path is None:
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[tuple[str, int]] = []
    name = ""
    length = 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records.append((name, length))
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name:
        records.append((name, length))
    if not records or sum(length for _, length in records) != alignment_length:
        return []
    mapped: list[tuple[str, int, int]] = []
    start = 1
    for record_name, record_length in records:
        end = start + record_length - 1
        mapped.append((record_name, start, end))
        start = end + 1
    return mapped


def _reference_segments(
    interval: ImportationInterval, records: list[tuple[str, int, int]]
) -> list[str]:
    segments = []
    for name, record_start, record_end in records:
        overlap_start = max(interval.start, record_start)
        overlap_end = min(interval.end, record_end)
        if overlap_start <= overlap_end:
            segments.append(
                f"{name}:{overlap_start - record_start + 1}-{overlap_end - record_start + 1}"
            )
    return segments


def _write_intervals(
    path: Path,
    intervals: list[ImportationInterval],
    records: list[tuple[str, int, int]],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "node",
                "alignment_start",
                "alignment_end",
                "length",
                "reference_segments",
                "crosses_reference_record_boundary",
            ]
        )
        for item in intervals:
            segments = _reference_segments(item, records)
            writer.writerow(
                [
                    item.node,
                    item.start,
                    item.end,
                    item.length,
                    ";".join(segments),
                    str(len(segments) > 1).lower(),
                ]
            )


def _write_profile(
    path: Path,
    *,
    branch_depth: np.ndarray,
    recombinant: np.ndarray,
    incomplete_only: np.ndarray,
    bins: int,
) -> list[dict[str, int | float]]:
    length = recombinant.size
    edges = np.linspace(0, length, bins + 1, dtype=int)
    rows: list[dict[str, int | float]] = []
    for index, (left, right) in enumerate(zip(edges[:-1], edges[1:], strict=True), start=1):
        width = right - left
        recombinant_sites = int(recombinant[left:right].sum())
        incomplete_sites = int(incomplete_only[left:right].sum())
        rows.append(
            {
                "bin": index,
                "alignment_start": left + 1,
                "alignment_end": right,
                "branch_importation_count": int(branch_depth[index - 1]),
                "recombination_removed_sites": recombinant_sites,
                "incomplete_only_removed_sites": incomplete_sites,
                "retained_sites": width - recombinant_sites - incomplete_sites,
                "recombination_removed_percent": 100.0 * recombinant_sites / width,
                "incomplete_only_removed_percent": 100.0 * incomplete_sites / width,
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _write_figure(rows: list[dict[str, int | float]], output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams["svg.hashsalt"] = "chronoclade-recombination"
    from matplotlib import pyplot

    x = np.asarray([(int(row["alignment_start"]) + int(row["alignment_end"])) / 2 for row in rows])
    depth = np.asarray([int(row["branch_importation_count"]) for row in rows])
    recombinant = np.asarray([float(row["recombination_removed_percent"]) for row in rows])
    incomplete = np.asarray([float(row["incomplete_only_removed_percent"]) for row in rows])

    figure, axes = pyplot.subplots(2, 1, figsize=(11.2, 6.0), sharex=True, constrained_layout=True)
    for axis in axes:
        axis.set_facecolor("#ffffff")
        axis.spines[["top", "right"]].set_visible(False)
        axis.spines[["left", "bottom"]].set_color("#343842")
        axis.tick_params(colors="#343842", labelsize=9)
        axis.grid(axis="y", color="#d7d9df", linewidth=0.7, alpha=0.75)
    axes[0].fill_between(x, depth, step="mid", color="#2855a6", alpha=0.86)
    axes[0].set_ylabel("Branches with\nan inferred import")
    axes[0].set_title(
        "Where ClonalFrameML inferred recombination",
        loc="left",
        fontsize=15,
        fontweight="bold",
        color="#17191f",
    )
    axes[1].stackplot(
        x,
        recombinant,
        incomplete,
        colors=["#d63c2f", "#9ba1ad"],
        labels=["Recombination", "Incomplete sites only"],
        step="mid",
    )
    axes[1].set_ylim(0, 100)
    axes[1].set_ylabel("Columns removed\nwithin each bin (%)")
    axes[1].set_xlabel("Reference-ordered alignment position")
    axes[1].legend(frameon=False, ncols=2, loc="upper left", fontsize=9)
    axes[1].ticklabel_format(axis="x", style="sci", scilimits=(6, 6), useMathText=True)
    for suffix, dpi in ((".svg", None), (".png", 220)):
        figure.savefig(
            output.with_suffix(suffix),
            dpi=dpi,
            facecolor="#ffffff",
            metadata={"Title": "ClonalFrameML recombination and alignment filtering"},
        )
    pyplot.close(figure)


def write_recombination_evidence(
    *,
    alignment: Path,
    filtered_alignment: Path,
    importations: Path,
    output_directory: Path,
    reference: Path | None = None,
    bins: int = 500,
) -> dict[str, object]:
    """Write exact interval tables, a genome profile, figures and a compact summary."""

    intervals = read_importation_intervals(importations)
    incomplete, sequences = _alignment_incomplete_mask(alignment)
    length = incomplete.size
    records = _reference_records(reference, length)
    recombinant = np.zeros(length, dtype=np.bool_)
    bin_count = min(max(1, bins), length)
    edges = np.linspace(0, length, bin_count + 1, dtype=int)
    branch_sets: list[set[str]] = [set() for _ in range(bin_count)]
    for interval in intervals:
        if interval.end > length:
            raise RecombinationReportError(
                f"Importation interval ends beyond the alignment ({interval.end} > {length})"
            )
        recombinant[interval.start - 1 : interval.end] = True
        first_bin = min(
            bin_count - 1, int(np.searchsorted(edges, interval.start - 1, side="right") - 1)
        )
        last_bin = min(
            bin_count - 1, int(np.searchsorted(edges, interval.end - 1, side="right") - 1)
        )
        for index in range(max(0, first_bin), last_bin + 1):
            branch_sets[index].add(interval.node)

    incomplete_only = incomplete & ~recombinant
    retained = ~(recombinant | incomplete)
    filtered_length = _fasta_alignment_length(filtered_alignment)
    retained_sites = int(retained.sum())
    if filtered_length != retained_sites:
        raise RecombinationReportError(
            "ClonalFrameML filtered-alignment length does not equal complete, non-imported "
            f"columns ({filtered_length} != {retained_sites})"
        )

    interval_path = output_directory / "recombination_intervals.tsv"
    profile_path = output_directory / "recombination_genome_profile.csv"
    _write_intervals(interval_path, intervals, records)
    rows = _write_profile(
        profile_path,
        branch_depth=np.asarray([len(nodes) for nodes in branch_sets]),
        recombinant=recombinant,
        incomplete_only=incomplete_only,
        bins=bin_count,
    )
    _write_figure(rows, output_directory / "recombination_map.svg")

    recombinant_sites = int(recombinant.sum())
    incomplete_only_sites = int(incomplete_only.sum())
    longest = sorted(intervals, key=lambda item: item.length, reverse=True)[:10]
    boundary_crossing = sum(len(_reference_segments(item, records)) > 1 for item in intervals)
    summary: dict[str, object] = {
        "coordinate_system": "1-based closed positions in the reference-ordered alignment",
        "alignment_sites": int(length),
        "sequences": sequences,
        "inferred_importation_intervals": len(intervals),
        "branches_with_inferred_importation": len({item.node for item in intervals}),
        "recombination_removed_sites": recombinant_sites,
        "recombination_removed_percent": 100.0 * recombinant_sites / length,
        "incomplete_only_removed_sites": incomplete_only_sites,
        "incomplete_only_removed_percent": 100.0 * incomplete_only_sites / length,
        "retained_clonal_sites": retained_sites,
        "retained_clonal_percent": 100.0 * retained_sites / length,
        "filtered_alignment_sites": filtered_length,
        "reference_records": len(records),
        "boundary_crossing_intervals": boundary_crossing,
        "longest_inferred_intervals": [
            {
                "node": item.node,
                "alignment_start": item.start,
                "alignment_end": item.end,
                "length": item.length,
                "reference_segments": _reference_segments(item, records),
                "crosses_reference_record_boundary": len(_reference_segments(item, records)) > 1,
            }
            for item in longest
        ],
        "filtering_rule": (
            "ClonalFrameML removed every column inferred as imported on any branch and every "
            "column containing an ambiguous base in any sequence."
        ),
    }
    (output_directory / "recombination_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary
