"""Fast Parsnp-style recombination profiling with PhiPack."""

from __future__ import annotations

import csv
import json
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from chronoclade.evidence import read_alignment
from chronoclade.errors import WorkflowError
from chronoclade.metadata import fasta_lengths


PROFILE_WINDOW = 100
PROFILE_STEP = 100
PROFILE_THRESHOLD = 0.01
PROFILE_BLOCK_SIZE = 250_000


@dataclass(frozen=True)
class ProfileChunk:
    """One bounded PhiPack Profile job within a reference contig."""

    index: int
    contig: int
    contig_offset: int
    window_start: int
    window_end: int
    keep_start: int
    keep_end: int


def _profile_chunks(
    contig_lengths: list[int], sequences: dict[str, str]
) -> list[ProfileChunk]:
    """Return bounded analysis blocks that never cross a reference-record join."""

    del sequences  # Missing calls remain within their reference segment; PhiPack handles them.

    chunks: list[ProfileChunk] = []
    offset = 0
    for contig, length in enumerate(contig_lengths, start=1):
        for block_start in range(0, length, PROFILE_BLOCK_SIZE):
            block_end = min(length, block_start + PROFILE_BLOCK_SIZE)
            if block_end - block_start < 1_000:
                continue
            chunks.append(
                ProfileChunk(
                    index=len(chunks),
                    contig=contig,
                    contig_offset=offset,
                    window_start=block_start,
                    window_end=block_end,
                    keep_start=block_start,
                    keep_end=block_end,
                )
            )
        offset += length
    return chunks


def _write_fasta(path: Path, sequences: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in sequences.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def _mask_sequence(sequence: str, regions: list[dict[str, object]]) -> str:
    parts: list[str] = []
    cursor = 0
    for region in regions:
        start = int(region["alignment_start"])
        end = int(region["alignment_end"])
        parts.extend((sequence[cursor:start], "N" * (end - start)))
        cursor = end
    parts.append(sequence[cursor:])
    return "".join(parts)


def _run_profile_chunk(
    chunk: ProfileChunk,
    sequences: dict[str, str],
    temporary: Path,
) -> list[tuple[int, int, float]]:
    run_dir = temporary / f"chunk_{chunk.index:04d}"
    run_dir.mkdir()
    alignment = run_dir / "alignment.fasta"
    absolute_start = chunk.contig_offset + chunk.window_start
    absolute_end = chunk.contig_offset + chunk.window_end
    _write_fasta(
        alignment,
        {name: sequence[absolute_start:absolute_end] for name, sequence in sequences.items()},
    )
    command = [
        "Profile",
        "-o",
        "-v",
        "-n",
        str(chunk.window_end - chunk.window_start),
        "-w",
        str(PROFILE_WINDOW),
        "-m",
        str(PROFILE_STEP),
        "-f",
        alignment.name,
    ]
    completed = subprocess.run(
        command,
        cwd=run_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    (run_dir / "Profile.log").write_text(
        "COMMAND\n" + " ".join(command) + "\n\nOUTPUT\n" + completed.stdout + completed.stderr,
        encoding="utf-8",
    )
    profile = run_dir / "Profile.csv"
    if completed.returncode != 0 or not profile.is_file():
        raise WorkflowError(
            f"PhiPack Profile failed for reference contig {chunk.contig}; "
            f"see {run_dir / 'Profile.log'}"
        )

    rows: list[tuple[int, int, float]] = []
    for raw_line in profile.read_text(encoding="utf-8").splitlines():
        try:
            position_text, p_value_text = raw_line.split(",", maxsplit=1)
            local_position = int(position_text)
            p_value = float(p_value_text)
        except ValueError:
            continue
        contig_position = chunk.window_start + local_position
        if chunk.keep_start <= contig_position < chunk.keep_end:
            rows.append((chunk.contig, contig_position, p_value))
    return rows


def _merged_mask(profile: list[tuple[int, int, float]], contig_lengths: list[int]) -> list[dict[str, object]]:
    intervals: list[tuple[int, int, int, float]] = []
    for contig, position, p_value in profile:
        if p_value < 0 or p_value >= PROFILE_THRESHOLD:
            continue
        start = max(0, position - PROFILE_WINDOW // 2)
        end = min(contig_lengths[contig - 1], position + PROFILE_WINDOW // 2)
        intervals.append((contig, start, end, p_value))
    intervals.sort()

    merged: list[dict[str, object]] = []
    offsets = [0]
    for length in contig_lengths[:-1]:
        offsets.append(offsets[-1] + length)
    for contig, start, end, p_value in intervals:
        if merged and int(merged[-1]["contig"]) == contig and start <= int(merged[-1]["end"]):
            merged[-1]["end"] = max(int(merged[-1]["end"]), end)
            merged[-1]["minimum_p_value"] = min(float(merged[-1]["minimum_p_value"]), p_value)
            merged[-1]["profile_hits"] = int(merged[-1]["profile_hits"]) + 1
            merged[-1]["alignment_end"] = offsets[contig - 1] + int(merged[-1]["end"])
            continue
        merged.append(
            {
                "contig": contig,
                "start": start,
                "end": end,
                "alignment_start": offsets[contig - 1] + start,
                "alignment_end": offsets[contig - 1] + end,
                "minimum_p_value": p_value,
                "profile_hits": 1,
            }
        )
    return merged


def run_phipack_filter(
    *,
    alignment: Path,
    reference: Path,
    output_alignment: Path,
    profile_output: Path,
    regions_output: Path,
    summary_output: Path,
    threads: int,
    force: bool,
) -> dict[str, object]:
    """Mask Parsnp-style PhiPack Profile regions without crossing contig joins."""

    alignment = alignment.expanduser().resolve()
    reference = reference.expanduser().resolve()
    output_alignment = output_alignment.expanduser().resolve()
    profile_output = profile_output.expanduser().resolve()
    regions_output = regions_output.expanduser().resolve()
    summary_output = summary_output.expanduser().resolve()

    if not force and all(
        path.is_file() for path in (output_alignment, profile_output, regions_output, summary_output)
    ):
        return json.loads(summary_output.read_text(encoding="utf-8"))

    sequences = read_alignment(alignment)
    alignment_size = len(next(iter(sequences.values())))
    contig_lengths = fasta_lengths(reference)
    if sum(contig_lengths) != alignment_size:
        raise WorkflowError(
            "SKA alignment length does not equal the concatenated reference-contig length; "
            "PhiPack windows cannot be mapped safely"
        )
    chunks = _profile_chunks(contig_lengths, sequences)
    if not chunks:
        raise WorkflowError("Reference contains no contig of at least 1,000 bases for PhiPack")

    with tempfile.TemporaryDirectory(prefix="chronoclade-phipack-", dir=alignment.parent) as tmp:
        temporary = Path(tmp)
        with ThreadPoolExecutor(max_workers=min(threads, len(chunks))) as executor:
            results = executor.map(
                lambda chunk: _run_profile_chunk(chunk, sequences, temporary),
                chunks,
            )
            profile = [row for rows in results for row in rows]

    profile.sort()
    with profile_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["reference_contig", "position", "p_value"])
        writer.writerows(profile)

    regions = _merged_mask(profile, contig_lengths)
    with regions_output.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "contig",
            "start",
            "end",
            "alignment_start",
            "alignment_end",
            "minimum_p_value",
            "profile_hits",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(regions)

    masked = bytearray(alignment_size)
    for region in regions:
        masked[int(region["alignment_start"]) : int(region["alignment_end"])] = b"\x01" * (
            int(region["alignment_end"]) - int(region["alignment_start"])
        )
    filtered = {name: _mask_sequence(sequence, regions) for name, sequence in sequences.items()}
    _write_fasta(output_alignment, filtered)

    valid_p_values = [p_value for _, _, p_value in profile if p_value >= 0]
    summary: dict[str, object] = {
        "method": "parsnp_style_phipack_profile",
        "window_size": PROFILE_WINDOW,
        "step_size": PROFILE_STEP,
        "maximum_profile_block_size": PROFILE_BLOCK_SIZE,
        "p_value_threshold": PROFILE_THRESHOLD,
        "reference_contigs": len(contig_lengths),
        "tested_core_blocks": len(chunks),
        "profile_positions": len(profile),
        "minimum_p_value": min(valid_p_values) if valid_p_values else None,
        "recombination_regions": len(regions),
        "masked_alignment_sites": sum(masked),
        "recombination_detected": bool(regions),
        "boundary_rule": (
            "Reference FASTA records were divided into bounded analysis blocks; no block crossed "
            "a reference-contig join. Ambiguous and missing calls remained within their reference "
            "segment and were handled by PhiPack."
        ),
    }
    summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary
