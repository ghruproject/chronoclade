"""Shared native SNP-tree and dating operations for the tree and time stages."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from chronoclade.coherence import screen_alignment
from chronoclade.location_network.marginal import UNKNOWN
from chronoclade.errors import WorkflowError
from chronoclade.metadata import Sample


@dataclass(frozen=True)
class LineageFiles:
    """Canonical paths for one lineage analysis."""

    directory: Path
    inputs: Path
    metadata: Path
    states: Path
    alignment: Path
    coherence_screen: Path
    clonal_coherence_screen: Path
    ska_prefix: Path
    ska_file: Path
    iqtree_prefix: Path
    starting_tree: Path
    clonalframe_prefix: Path
    tree: Path
    importations: Path
    filtered_alignment: Path
    clock_dir: Path
    clock_file: Path
    time_tree: Path

    @classmethod
    def in_directory(cls, directory: Path) -> LineageFiles:
        return cls(
            directory=directory,
            inputs=directory / "inputs.tsv",
            metadata=directory / "metadata.csv",
            states=directory / "states.csv",
            alignment=directory / "core_alignment.fasta",
            coherence_screen=directory / "lineage_coherence.tsv",
            clonal_coherence_screen=directory / "clonal_lineage_coherence.tsv",
            ska_prefix=directory / "ska",
            ska_file=directory / "ska.skf",
            iqtree_prefix=directory / "iqtree",
            starting_tree=directory / "iqtree.treefile",
            clonalframe_prefix=directory / "clonalframeml",
            tree=directory / "clonalframeml.labelled_tree.newick",
            importations=directory / "clonalframeml.importation_status.txt",
            filtered_alignment=directory / "clonalframeml.filtered.fasta",
            clock_dir=directory / "clock",
            clock_file=directory / "clock" / "molecular_clock.txt",
            time_tree=directory / "timetree" / "timetree.nexus",
        )


def _write_csv(
    path: Path,
    fieldnames: list[str],
    rows: Iterable[dict[str, object]],
    *,
    delimiter: str = ",",
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def _write_lineage_inputs(directory: Path, samples: list[Sample]) -> tuple[Path, Path, Path]:
    inputs = directory / "inputs.tsv"
    inputs.write_text(
        "".join(f"{sample.sample_id}\t{sample.assembly}\n" for sample in samples),
        encoding="utf-8",
    )
    metadata = directory / "metadata.csv"
    fields = [
        "sample_id",
        "assembly",
        "collection_date",
        "location",
        "country",
        "region",
        "city",
        "species",
        "lineage",
        "origin",
        "is_reference",
        "patient_id",
    ]
    _write_csv(
        metadata,
        fields,
        (
            {
                **asdict(sample),
                "assembly": str(sample.assembly),
            }
            for sample in samples
        ),
    )
    states = directory / "states.csv"
    _write_csv(
        states,
        ["sample_id", "location", "origin"],
        (
            {
                "sample_id": sample.sample_id,
                "location": (
                    "?" if sample.location.strip().casefold() in UNKNOWN else sample.location
                ),
                "origin": sample.origin,
            }
            for sample in samples
        ),
    )
    return inputs, metadata, states


def _run_command(
    command: list[str],
    *,
    log: Path,
    expected: Path | tuple[Path, ...],
    force: bool,
    cwd: Path | None = None,
    inputs: tuple[Path, ...] = (),
) -> None:
    expected_paths = (expected,) if isinstance(expected, Path) else expected
    fingerprint_path = log.with_suffix(".fingerprint.json")
    fingerprint = {
        "command": command,
        "inputs": [_file_fingerprint(path) for path in inputs],
    }
    previous = None
    if fingerprint_path.is_file():
        try:
            previous = json.loads(fingerprint_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = None
    if all(path.exists() for path in expected_paths) and not force and previous == fingerprint:
        from chronoclade.progress import progress
        progress("native", f"Reusing checked {Path(command[0]).name} output")
        return
    for path in expected_paths:
        if path.is_file():
            path.unlink()
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND\n" + " ".join(command) + "\n\nOUTPUT\n")
        from chronoclade.progress import progress_task
        with progress_task("native", Path(command[0]).name, details=f"log: {log}"):
            completed = subprocess.run(
                command,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=False,
                cwd=cwd,
            )
            missing = [path for path in expected_paths if not path.exists()]
            if completed.returncode != 0 or missing:
                missing_text = ", ".join(str(path) for path in missing)
                raise WorkflowError(
                    f"Command failed or did not create expected output(s): {missing_text}. See log: {log}"
                )
    temporary_fingerprint = fingerprint_path.with_suffix(".json.tmp")
    temporary_fingerprint.write_text(json.dumps(fingerprint, indent=2) + "\n", encoding="utf-8")
    temporary_fingerprint.replace(fingerprint_path)


def _file_fingerprint(path: Path) -> dict[str, object]:
    """Hash one direct stage input so stale outputs are not silently reused."""

    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise WorkflowError(f"Stage input does not exist: {resolved}")
    digest = hashlib.sha256()
    with resolved.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def alignment_length(path: Path) -> int:
    """Return the number of columns in an equal-length FASTA alignment."""

    lengths: list[int] = []
    current = 0
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                    current = 0
            else:
                current += len(line)
    if current:
        lengths.append(current)
    if not lengths:
        raise WorkflowError(f"Alignment contains no sequences: {path}")
    if len(set(lengths)) != 1:
        raise WorkflowError(f"Alignment sequences are not the same length: {path}")
    return lengths[0]


def complete_alignment_sites(path: Path) -> int:
    """Count columns containing only unambiguous A/C/G/T bases."""

    length = alignment_length(path)
    incomplete = bytearray(length)
    position = 0
    seen_header = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header and position != length:
                    raise WorkflowError(f"Alignment sequences are not the same length: {path}")
                seen_header = True
                position = 0
                continue
            if not seen_header:
                raise WorkflowError(f"Alignment sequence appears before its FASTA header: {path}")
            for base in line.upper():
                if position >= length:
                    raise WorkflowError(f"Alignment sequences are not the same length: {path}")
                if base not in "ACGT":
                    incomplete[position] = 1
                position += 1
    complete = length - sum(incomplete)
    if complete == 0:
        raise WorkflowError(f"Alignment contains no complete A/C/G/T sites: {path}")
    return complete


def _run_core_phylogeny(
    files: LineageFiles,
    members: list[Sample],
    reference: Sample,
    threads: int,
    force: bool,
    coherence_policy: str = "strict",
    seed: int | None = None,
) -> tuple[Path, Path, dict[str, object] | None]:
    if coherence_policy not in {"strict", "diagnostic"}:
        raise WorkflowError("Coherence policy must be strict or diagnostic")
    mapping_stages = [
        (
            [
                "ska",
                "build",
                "-f",
                str(files.inputs),
                "-o",
                str(files.ska_prefix),
                "--threads",
                str(threads),
            ],
            "ska_build.log",
            files.ska_file,
            (files.inputs, *(sample.assembly for sample in members)),
        ),
        (
            [
                "ska",
                "map",
                str(reference.assembly),
                str(files.ska_file),
                "-o",
                str(files.alignment),
                "--ambig-mask",
                "--repeat-mask",
                "--threads",
                str(threads),
            ],
            "ska_map.log",
            files.alignment,
            (reference.assembly, files.ska_file),
        ),
    ]
    for command, log_name, expected, inputs in mapping_stages:
        _run_command(
            command,
            log=files.directory / "logs" / log_name,
            expected=expected,
            force=force,
            cwd=files.directory,
            inputs=inputs,
        )

    coherence = screen_alignment(
        files.alignment,
        members,
        output=files.coherence_screen,
    )
    if coherence["flagged_samples"]:
        from chronoclade.progress import progress

        progress(
            "tree",
            "Raw-distance review: retaining all selected genomes until the "
            "recombination-filtered screen; raw flags can reflect recombination. "
            f"Flagged genomes: {', '.join(coherence['flagged_samples'])}. "
            f"Evidence: {files.coherence_screen}",
        )

    phylogeny_stages = [
        (
            [
                "iqtree",
                *(["-seed", str(seed)] if seed is not None else []),
                "-s",
                str(files.alignment),
                "-m",
                "GTR+G",
                "-T",
                str(threads),
                "-pre",
                str(files.iqtree_prefix),
                "-redo",
            ],
            "iqtree.log",
            files.starting_tree,
            (files.alignment,),
        ),
        (
            [
                "ClonalFrameML",
                str(files.starting_tree),
                str(files.alignment),
                str(files.clonalframe_prefix),
                "-ignore_incomplete_sites",
                "true",
                "-output_filtered",
                "true",
                "-num_threads",
                str(threads),
            ],
            "clonalframeml.log",
            (files.tree, files.importations, files.filtered_alignment),
            (files.starting_tree, files.alignment),
        ),
    ]
    for command, log_name, expected, inputs in phylogeny_stages:
        _run_command(
            command,
            log=files.directory / "logs" / log_name,
            expected=expected,
            force=force,
            cwd=files.directory,
            inputs=inputs,
        )

    clonal_coherence = screen_alignment(
        files.filtered_alignment,
        members,
        output=files.clonal_coherence_screen,
    )
    if clonal_coherence["flagged_samples"] and coherence_policy == "strict":
        names = ", ".join(str(name) for name in clonal_coherence["flagged_samples"])
        raise WorkflowError(
            "Lineage-coherence screen found extreme recombination-filtered distance "
            f"outlier(s): {names}. Relative-distance flags can reflect legitimate cohort "
            f"diversity and do not establish incorrect genome identity. Review the cohort "
            f"or explicitly choose diagnostic policy. Evidence: {files.clonal_coherence_screen}"
        )
    return files.tree, files.filtered_alignment, None


def _clock_dates_path(files: LineageFiles) -> Path:
    """Prepare precise external clock inputs without changing saved metadata."""
    from chronoclade.metadata_dates import treetime_date

    source = files.directory / "clock_dates.csv"
    if not source.is_file():
        source = files.metadata
    with source.open(newline="", encoding="utf-8") as stream:
        rows = [
            {"sample_id": row["sample_id"],
             "collection_date": treetime_date(row.get("collection_date"))}
            for row in csv.DictReader(stream)
        ]
    target = files.directory / "treetime_dates.csv"
    _write_csv(target, ["sample_id", "collection_date"], rows)
    return target


def _run_observed_clock(files: LineageFiles, tree: Path, sequence_length: int, force: bool) -> None:
    _run_command(
        [
            "treetime",
            "clock",
            "--tree",
            str(tree),
            "--dates",
            str(_clock_dates_path(files)),
            "--name-column",
            "sample_id",
            "--date-column",
            "collection_date",
            "--sequence-length",
            str(sequence_length),
            "--reroot",
            "least-squares",
            "--allow-negative-rate",
            "--clock-filter",
            "0",
            "--plot-rtt",
            "root_to_tip_regression.svg",
            "--outdir",
            str(files.clock_dir),
        ],
        log=files.directory / "logs" / "clock.log",
        expected=files.clock_file,
        force=force,
        inputs=(tree, _clock_dates_path(files)),
    )


def _dated_members(members: list[Sample]) -> list[Sample]:
    """Keep genuinely usable date values without inventing dates for other tips."""
    from chronoclade.metadata import MetadataError, _validate_date

    from datetime import date

    result = []
    for sample in members:
        try:
            _validate_date(sample.collection_date, row_number=0)
            # Partial dates are represented by their earliest possible day.
            parts = [int(part) for part in sample.collection_date.split("-")]
            earliest = date(parts[0], parts[1] if len(parts) > 1 else 1,
                            parts[2] if len(parts) > 2 else 1)
            if earliest > date.today():
                continue
        except MetadataError:
            continue
        result.append(sample)
    return result


def _run_dated_tree(
    files: LineageFiles,
    sequence_length: int,
    force: bool,
    tree: Path | None = None,
    *,
    reroot: bool = False,
    seed: int | None = None,
    keep_polytomies: bool = False,
) -> None:
    rooted_tree = files.clock_dir / "rerooted.newick"
    input_tree = tree if tree is not None else files.tree
    command = [
        "treetime",
        "--tree",
        str(input_tree if reroot else rooted_tree),
        "--dates",
        str(_clock_dates_path(files)),
        "--name-column",
        "sample_id",
        "--date-column",
        "collection_date",
        "--sequence-length",
        str(sequence_length),
        "--confidence",
        "--time-marginal",
        "only-final",
        "--covariation",
        "--clock-filter",
        "0",
        "--plot-tree",
        "timetree.svg",
        "--plot-rtt",
        "root_to_tip_regression.svg",
        "--outdir",
        str(files.directory / "timetree"),
    ]
    if seed is not None:
        command.extend(["--rng-seed", str(seed)])
    if keep_polytomies:
        command.append("--keep-polytomies")
    if reroot:
        command.extend(["--reroot", "least-squares"])
    else:
        command.append("--keep-root")
    _run_command(
        command,
        log=files.directory / "logs" / "timetree.log",
        expected=files.time_tree,
        force=force,
        inputs=((input_tree if reroot else rooted_tree), _clock_dates_path(files)),
    )
