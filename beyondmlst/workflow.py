"""Execution of the beyondMLST workflow."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from beyondmlst.evidence import EvidenceError, build_public_health_evidence
from beyondmlst.metadata import (
    Sample,
    group_samples,
    select_reference,
    slugify_lineage,
)
from beyondmlst.report import (
    assess_temporal_signal,
    write_lineage_report,
    write_summary_report,
)
from beyondmlst.temporal import run_date_randomisation

REQUIRED_TOOLS = (
    "ska",
    "iqtree",
    "ClonalFrameML",
    "treetime",
)


class WorkflowError(RuntimeError):
    """Raised when an external workflow stage fails."""


@dataclass(frozen=True)
class ResourcePlan:
    """Bounded concurrency derived from the user's global CPU budget."""

    total_threads: int
    lineage_workers: int
    threads_per_lineage: int
    randomisation_workers_per_lineage: int


def allocate_resources(
    *,
    total_threads: int,
    lineage_jobs: int,
    randomisation_jobs: int,
    ready_lineages: int,
) -> ResourcePlan:
    """Allocate lineage and permutation workers without oversubscribing CPUs."""

    if min(total_threads, lineage_jobs, randomisation_jobs) < 1:
        raise ValueError("thread and job counts must be at least 1")
    if ready_lineages < 0:
        raise ValueError("ready_lineages cannot be negative")

    lineage_workers = min(lineage_jobs, ready_lineages, total_threads)
    if lineage_workers == 0:
        return ResourcePlan(total_threads, 0, total_threads, min(randomisation_jobs, total_threads))

    threads_per_lineage = max(1, total_threads // lineage_workers)
    return ResourcePlan(
        total_threads=total_threads,
        lineage_workers=lineage_workers,
        threads_per_lineage=threads_per_lineage,
        randomisation_workers_per_lineage=min(randomisation_jobs, threads_per_lineage),
    )


def tool_status() -> dict[str, str | None]:
    return {tool: shutil.which(tool) for tool in REQUIRED_TOOLS}


def native_platform_supported() -> bool:
    """Return whether the locked Pixi tools support this operating system."""

    return sys.platform.startswith(("linux", "darwin"))


def check_tools() -> None:
    if not native_platform_supported():
        raise WorkflowError(
            "Native execution is supported on macOS and Linux. This platform is not in the "
            "locked Pixi environment."
        )
    missing = [tool for tool, path in tool_status().items() if path is None]
    if missing:
        raise WorkflowError(
            "Missing workflow tools: "
            + ", ".join(missing)
            + ". Run this command through `pixi run`."
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
                "location": sample.location,
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
) -> None:
    expected_paths = (expected,) if isinstance(expected, Path) else expected
    if all(path.exists() for path in expected_paths) and not force:
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND\n" + " ".join(command) + "\n\nOUTPUT\n")
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


def plan(samples: list[Sample], *, min_samples: int) -> list[dict[str, object]]:
    """Return the deterministic lineage plan without running external tools."""

    result: list[dict[str, object]] = []
    for (species, lineage), members in sorted(group_samples(samples).items()):
        reference = select_reference(members)
        result.append(
            {
                "species": species,
                "lineage": lineage,
                "slug": slugify_lineage(species, lineage),
                "sample_count": len(members),
                "distinct_dates": len({sample.collection_date for sample in members}),
                "locations": sorted({sample.location for sample in members}),
                "reference": reference.sample_id,
                "status": "ready" if len(members) >= min_samples else "too_few_samples",
            }
        )
    return result


def read_context_manifest(path: Path | None) -> list[dict[str, str]]:
    """Read a context-selection manifest for report provenance."""

    if path is None:
        return []
    path = path.expanduser().resolve()
    if not path.is_file():
        raise WorkflowError(f"Context manifest does not exist: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "species", "lineage"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise WorkflowError(
                "Context manifest must be tab-separated and contain sample_id, species and lineage"
            )
        return [dict(row) for row in reader]


def context_evidence(
    members: list[Sample], manifest_rows: list[dict[str, str]], *, directory: Path
) -> dict[str, object]:
    """Summarise contextual sampling and preserve the lineage manifest subset."""

    directory.mkdir(parents=True, exist_ok=True)
    origins: dict[str, int] = {}
    for member in members:
        key = member.origin.strip().lower()
        origins[key] = origins.get(key, 0) + 1
    contexts = [member for member in members if member.origin.strip().lower() == "context"]
    local = [member for member in members if member.origin.strip().lower() == "local"]
    relevant_rows = [
        row
        for row in manifest_rows
        if row.get("species") == members[0].species and row.get("lineage") == members[0].lineage
    ]
    nearest: list[dict[str, object]] = []
    for row in relevant_rows:
        raw_distance = row.get("min_ska_distance", "")
        try:
            distance: float | None = float(raw_distance) if raw_distance else None
        except ValueError:
            distance = None
        nearest.append(
            {
                "sample_id": row.get("sample_id", ""),
                "country": row.get("country", ""),
                "collection_date": row.get("collection_date", ""),
                "nearest_focal": row.get("nearest_focal", ""),
                "min_ska_distance": distance,
                "selection_reason": row.get("selection_reason", ""),
            }
        )
    nearest.sort(
        key=lambda row: (
            row["min_ska_distance"] is None,
            row["min_ska_distance"] if row["min_ska_distance"] is not None else float("inf"),
            str(row["sample_id"]),
        )
    )

    subset_path: Path | None = None
    if relevant_rows:
        subset_path = directory / "context_manifest.tsv"
        fieldnames = list(relevant_rows[0])
        _write_csv(subset_path, fieldnames, relevant_rows, delimiter="\t")
    return {
        "local_samples": len(local),
        "context_samples": len(contexts),
        "retrospective_samples": origins.get("retrospective", 0),
        "origin_counts": origins,
        "context_locations": sorted({member.location for member in contexts}),
        "manifest_available": bool(relevant_rows),
        "manifest_path": str(subset_path) if subset_path else None,
        "nearest_screening_contexts": nearest[:10],
        "interpretation": (
            "Contextual genomes are available for topology and relatedness assessment."
            if contexts
            else "No contextual genomes were supplied; introductions cannot be distinguished "
            "from local circulation."
        ),
    }


def _run_lineage(
    item: dict[str, object],
    members: list[Sample],
    *,
    output: Path,
    threads: int,
    randomisation_jobs: int,
    randomisations: int,
    temporal_p_value: float,
    seed: int,
    force: bool,
    context_manifest_rows: list[dict[str, str]],
) -> dict[str, object]:
    """Run one ready lineage; independent lineages may call this concurrently."""

    lineage_dir = output / str(item["slug"])
    lineage_dir.mkdir(parents=True, exist_ok=True)
    reference = select_reference(members)
    inputs, metadata, states = _write_lineage_inputs(lineage_dir, members)
    alignment = lineage_dir / "core_alignment.fasta"
    ska_prefix = lineage_dir / "ska"
    ska_file = lineage_dir / "ska.skf"
    iqtree_prefix = lineage_dir / "iqtree"
    starting_tree = lineage_dir / "iqtree.treefile"
    clonalframe_prefix = lineage_dir / "clonalframeml"
    tree = lineage_dir / "clonalframeml.labelled_tree.newick"
    importations = lineage_dir / "clonalframeml.importation_status.txt"
    filtered_alignment = lineage_dir / "clonalframeml.filtered.fasta"
    context_summary = context_evidence(members, context_manifest_rows, directory=lineage_dir)

    _run_command(
        [
            "ska",
            "build",
            "-f",
            str(inputs),
            "-o",
            str(ska_prefix),
            "--threads",
            str(threads),
        ],
        log=lineage_dir / "logs" / "ska_build.log",
        expected=ska_file,
        force=force,
        cwd=lineage_dir,
    )
    _run_command(
        [
            "ska",
            "map",
            str(reference.assembly),
            str(ska_file),
            "-o",
            str(alignment),
            "--ambig-mask",
            "--repeat-mask",
            "--threads",
            str(threads),
        ],
        log=lineage_dir / "logs" / "ska_map.log",
        expected=alignment,
        force=force,
        cwd=lineage_dir,
    )
    _run_command(
        [
            "iqtree",
            "-s",
            str(alignment),
            "-m",
            "GTR+G",
            "-T",
            str(threads),
            "-pre",
            str(iqtree_prefix),
            "-redo",
        ],
        log=lineage_dir / "logs" / "iqtree.log",
        expected=starting_tree,
        force=force,
        cwd=lineage_dir,
    )
    _run_command(
        [
            "ClonalFrameML",
            str(starting_tree),
            str(alignment),
            str(clonalframe_prefix),
            "-ignore_incomplete_sites",
            "true",
            "-output_filtered",
            "true",
            "-num_threads",
            str(threads),
        ],
        log=lineage_dir / "logs" / "clonalframeml.log",
        expected=(tree, importations, filtered_alignment),
        force=force,
        cwd=lineage_dir,
    )

    sequence_length = complete_alignment_sites(alignment)
    clock_dir = lineage_dir / "clock"
    clock_file = clock_dir / "molecular_clock.txt"
    _run_command(
        [
            "treetime",
            "clock",
            "--tree",
            str(tree),
            "--dates",
            str(metadata),
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
            str(clock_dir),
        ],
        log=lineage_dir / "logs" / "clock.log",
        expected=clock_file,
        force=force,
    )

    temporal_path = lineage_dir / "temporal_signal.json"
    if force or not temporal_path.exists():
        temporal = run_date_randomisation(
            tree=tree,
            sequence_length=sequence_length,
            samples=members,
            observed_clock=clock_file,
            randomisations=randomisations,
            randomisation_jobs=randomisation_jobs,
            seed=seed,
            output=temporal_path,
        )
    else:
        temporal = json.loads(temporal_path.read_text(encoding="utf-8"))

    assessment = assess_temporal_signal(
        temporal,
        p_value_threshold=temporal_p_value,
    )
    temporal_supported = bool(assessment["supported"])
    rooted_tree = clock_dir / "rerooted.newick"
    public_health = build_public_health_evidence(
        filtered_alignment=filtered_alignment,
        rooted_tree=rooted_tree if rooted_tree.is_file() else tree,
        samples=members,
        output=lineage_dir,
        temporal_assessment=assessment,
    )

    time_tree = lineage_dir / "timetree" / "timetree.nexus"
    if temporal_supported:
        _run_command(
            [
                "treetime",
                "--tree",
                str(tree),
                "--dates",
                str(metadata),
                "--name-column",
                "sample_id",
                "--date-column",
                "collection_date",
                "--sequence-length",
                str(sequence_length),
                "--confidence",
                "--time-marginal",
                "only-final",
                "--reroot",
                "least-squares",
                "--covariation",
                "--clock-filter",
                "0",
                "--plot-tree",
                "timetree.svg",
                "--plot-rtt",
                "root_to_tip_regression.svg",
                "--outdir",
                str(lineage_dir / "timetree"),
            ],
            log=lineage_dir / "logs" / "timetree.log",
            expected=time_tree,
            force=force,
        )

    time_tree_available = temporal_supported and time_tree.exists()
    location_tree = time_tree if time_tree_available else tree
    location_output = lineage_dir / "location" / "annotated_tree.nexus"
    _run_command(
        [
            "treetime",
            "mugration",
            "--tree",
            str(location_tree),
            "--states",
            str(states),
            "--name-column",
            "sample_id",
            "--attribute",
            "location",
            "--confidence",
            "--outdir",
            str(lineage_dir / "location"),
        ],
        log=lineage_dir / "logs" / "location.log",
        expected=location_output,
        force=force,
    )

    report = {
        **item,
        "reference_path": str(reference.assembly),
        "alignment_length": alignment_length(alignment),
        "complete_alignment_sites": sequence_length,
        "resources": {
            "threads": threads,
            "date_randomisation_jobs": randomisation_jobs,
        },
        "temporal_status": assessment["code"],
        "temporal_signal_reason": assessment["reason"],
        "temporal_signal_supported": temporal_supported,
        "temporal_signal": temporal,
        "context": context_summary,
        "public_health_status": public_health["scenario"]["code"],
        "public_health_label": public_health["scenario"]["label"],
        "public_health_confidence": public_health["scenario"]["confidence"],
        "public_health": public_health,
        "outputs": {
            "alignment": str(alignment),
            "starting_tree": str(starting_tree),
            "tree": str(tree),
            "recombination_importations": str(importations),
            "filtered_alignment": str(filtered_alignment),
            "root_to_tip_plot": str(clock_dir / "root_to_tip_regression.svg"),
            "date_randomisation_plot": str(lineage_dir / "date_randomisation.svg"),
            "timetree": str(time_tree) if time_tree_available else None,
            "timetree_plot": (
                str(lineage_dir / "timetree" / "timetree.svg") if time_tree_available else None
            ),
            "location_tree": str(location_output),
            "public_health_evidence": str(lineage_dir / "public_health_evidence.json"),
            "clonal_pairwise_distances": str(lineage_dir / "clonal_pairwise_distances.tsv"),
            "clonal_snp_matrix": str(lineage_dir / "clonal_snp_matrix.tsv"),
            "pairwise_callable_sites": str(lineage_dir / "pairwise_callable_sites.tsv"),
            "clonal_snp_heatmap": str(lineage_dir / "clonal_snp_heatmap.svg"),
            "html_report": str(lineage_dir / "report.html"),
        },
    }
    write_lineage_report(report, directory=lineage_dir, p_value_threshold=temporal_p_value)
    (lineage_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def run_workflow(
    samples: list[Sample],
    *,
    output: Path,
    threads: int,
    lineage_jobs: int,
    randomisation_jobs: int,
    randomisations: int,
    temporal_p_value: float,
    min_samples: int,
    seed: int,
    force: bool,
    context_manifest: Path | None = None,
) -> dict[str, object]:
    """Run every lineage through alignment, recombination and dating stages."""

    check_tools()
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    lineage_plan = plan(samples, min_samples=min_samples)
    groups = group_samples(samples)
    manifest_rows = read_context_manifest(context_manifest)
    ready_items = [item for item in lineage_plan if item["status"] == "ready"]
    resources = allocate_resources(
        total_threads=threads,
        lineage_jobs=lineage_jobs,
        randomisation_jobs=randomisation_jobs,
        ready_lineages=len(ready_items),
    )

    def run_ready_item(indexed_item: tuple[int, dict[str, object]]) -> dict[str, object]:
        index, item = indexed_item
        key = (str(item["species"]), str(item["lineage"]))
        try:
            return _run_lineage(
                item,
                groups[key],
                output=output,
                threads=resources.threads_per_lineage,
                randomisation_jobs=resources.randomisation_workers_per_lineage,
                randomisations=randomisations,
                temporal_p_value=temporal_p_value,
                seed=seed + index,
                force=force,
                context_manifest_rows=manifest_rows,
            )
        except (EvidenceError, OSError, WorkflowError) as error:
            raise WorkflowError(f"Lineage {key[0]} / {key[1]} failed: {error}") from error

    completed_by_slug: dict[str, dict[str, object]] = {}
    if ready_items:
        indexed_items = list(enumerate(ready_items))
        with ThreadPoolExecutor(max_workers=resources.lineage_workers) as executor:
            for report in executor.map(run_ready_item, indexed_items):
                completed_by_slug[str(report["slug"])] = report

    reports = [
        completed_by_slug[str(item["slug"])] if item["status"] == "ready" else item
        for item in lineage_plan
    ]

    summary = {
        "workflow": "beyondmlst",
        "resources": asdict(resources),
        "lineages": reports,
        "context_manifest": str(context_manifest.expanduser().resolve())
        if context_manifest
        else None,
        "guardrail": (
            "Location-state reconstructions are exploratory and do not by themselves establish "
            "direct transmission or a definitive number of introductions."
        ),
    }
    summary["report"] = str(write_summary_report(summary, output=output))
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary
