"""Fast recombination screening with PhiPack Profile."""

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


def _profile_chunks(contig_lengths: list[int]) -> list[ProfileChunk]:
    """Return bounded analysis blocks that never cross a reference-record join."""

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


def _run_profile_chunk(
    chunk: ProfileChunk,
    sequences: dict[str, str],
    temporary: Path,
) -> list[tuple[int, int, int, int, float]]:
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

    rows: list[tuple[int, int, int, int, float]] = []
    for raw_line in profile.read_text(encoding="utf-8").splitlines():
        try:
            position_text, p_value_text = raw_line.split(",", maxsplit=1)
            local_position = int(position_text)
            p_value = float(p_value_text)
        except ValueError:
            continue
        contig_position = chunk.window_start + local_position
        if chunk.window_start <= contig_position < chunk.window_end:
            rows.append(
                (
                    chunk.contig,
                    chunk.window_start,
                    chunk.window_end,
                    contig_position,
                    p_value,
                )
            )
    return rows


def run_phipack_screen(
    *,
    alignment: Path,
    reference: Path,
    profile_output: Path,
    significant_blocks_output: Path,
    summary_output: Path,
    threads: int,
    force: bool,
) -> dict[str, object]:
    """Detect PHI-positive analysis blocks without claiming tract localisation."""

    alignment = alignment.expanduser().resolve()
    reference = reference.expanduser().resolve()
    profile_output = profile_output.expanduser().resolve()
    significant_blocks_output = significant_blocks_output.expanduser().resolve()
    summary_output = summary_output.expanduser().resolve()

    if not force and all(
        path.is_file() for path in (profile_output, significant_blocks_output, summary_output)
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
    chunks = _profile_chunks(contig_lengths)
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
        writer.writerow(
            ["reference_contig", "block_start", "block_end", "profile_position", "p_value"]
        )
        writer.writerows(profile)

    significant_results = [row for row in profile if 0 <= row[4] < PROFILE_THRESHOLD]
    block_hits: dict[tuple[int, int, int], list[float]] = {}
    for contig, block_start, block_end, _position, p_value in significant_results:
        block_hits.setdefault((contig, block_start, block_end), []).append(p_value)
    significant_blocks = [
        (*block, min(p_values), len(p_values))
        for block, p_values in sorted(block_hits.items())
    ]
    with significant_blocks_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "reference_contig",
                "block_start",
                "block_end",
                "minimum_p_value",
                "significant_profile_positions",
            ]
        )
        writer.writerows(significant_blocks)

    valid_p_values = [row[4] for row in profile if row[4] >= 0]
    summary: dict[str, object] = {
        "method": "phipack_profile_screen",
        "profile_settings_source": "Parsnp PhiPack Profile settings",
        "window_size": PROFILE_WINDOW,
        "step_size": PROFILE_STEP,
        "maximum_profile_block_size": PROFILE_BLOCK_SIZE,
        "p_value_threshold": PROFILE_THRESHOLD,
        "reference_contigs": len(contig_lengths),
        "tested_core_blocks": len(chunks),
        "profile_positions": len(profile),
        "significant_profile_positions": len(significant_results),
        "minimum_p_value": min(valid_p_values) if valid_p_values else None,
        "significant_blocks": len(significant_blocks),
        "recombination_detected": bool(significant_blocks),
        "interpretation": (
            "One or more blocks were PHI-positive at the unadjusted screening threshold; "
            "run the full ClonalFrameML workflow before interpreting topology or dates."
            if significant_blocks
            else "No block was PHI-positive at the configured screening threshold."
        ),
        "multiple_testing_limit": (
            "The p < 0.01 threshold is applied to each computational block without a "
            "multiple-testing correction. This favours sensitivity for triage; a positive "
            "screen is an escalation signal, not a recombination tract call."
        ),
        "boundary_rule": (
            "Reference FASTA records were divided into bounded analysis blocks; no block crossed "
            "a reference-contig join. Ambiguous and missing calls remained within their reference "
            "segment and were handled by PhiPack."
        ),
        "localisation_limit": (
            "The fixed blocks are computational screening units. A significant block indicates "
            "PHI evidence of incompatibility somewhere within it; it does not localise a "
            "recombinant tract and no alignment sites were masked."
        ),
    }
    summary_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary
