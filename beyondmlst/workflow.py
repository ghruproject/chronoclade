"""Execution of the beyondMLST workflow."""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Iterable

from beyondmlst.metadata import (
    Sample,
    assembly_length,
    group_samples,
    select_reference,
    slugify_lineage,
)
from beyondmlst.temporal import run_date_randomisation

REQUIRED_TOOLS = (
    "generate_ska_alignment.py",
    "run_gubbins.py",
    "mask_gubbins_aln.py",
    "treetime",
)


class WorkflowError(RuntimeError):
    """Raised when an external workflow stage fails."""


def tool_status() -> dict[str, str | None]:
    return {tool: shutil.which(tool) for tool in REQUIRED_TOOLS}


def native_platform_supported() -> bool:
    """Gubbins' locked Bioconda build is validated natively on Linux."""

    return sys.platform.startswith("linux")


def check_tools() -> None:
    if not native_platform_supported() and os.environ.get("BEYONDMLST_ALLOW_UNSUPPORTED") != "1":
        raise WorkflowError(
            "Native execution is currently supported on Linux. On macOS or Windows, use the "
            "beyondMLST Docker image; set BEYONDMLST_ALLOW_UNSUPPORTED=1 only for development."
        )
    missing = [tool for tool, path in tool_status().items() if path is None]
    if missing:
        raise WorkflowError(
            "Missing workflow tools: "
            + ", ".join(missing)
            + ". Run this command through `pixi run`."
        )


def _write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
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


def _run_command(command: list[str], *, log: Path, expected: Path, force: bool) -> None:
    if expected.exists() and not force:
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND\n" + " ".join(command) + "\n\nOUTPUT\n")
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if completed.returncode != 0 or not expected.exists():
        raise WorkflowError(
            f"Command failed or did not create {expected}. See log: {log}"
        )


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


def run_workflow(
    samples: list[Sample],
    *,
    output: Path,
    threads: int,
    randomisations: int,
    temporal_p_value: float,
    min_samples: int,
    seed: int,
    force: bool,
) -> dict[str, object]:
    """Run every lineage through alignment, recombination and dating stages."""

    check_tools()
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    reports: list[dict[str, object]] = []

    for item in plan(samples, min_samples=min_samples):
        species = str(item["species"])
        lineage = str(item["lineage"])
        members = group_samples(samples)[(species, lineage)]
        lineage_dir = output / str(item["slug"])
        lineage_dir.mkdir(parents=True, exist_ok=True)
        if item["status"] != "ready":
            reports.append(item)
            continue

        reference = select_reference(members)
        inputs, metadata, states = _write_lineage_inputs(lineage_dir, members)
        alignment = lineage_dir / "core_alignment.fasta"
        gubbins_prefix = lineage_dir / "gubbins"
        tree = lineage_dir / "gubbins.final_tree.tre"
        gff = lineage_dir / "gubbins.recombination_predictions.gff"
        masked = lineage_dir / "core_alignment.recombination_masked.fasta"

        _run_command(
            [
                "generate_ska_alignment.py",
                "--reference",
                str(reference.assembly),
                "--input",
                str(inputs),
                "--out",
                str(alignment),
                "--threads",
                str(threads),
            ],
            log=lineage_dir / "logs" / "alignment.log",
            expected=alignment,
            force=force,
        )
        _run_command(
            [
                "run_gubbins.py",
                "--prefix",
                str(gubbins_prefix),
                "--threads",
                str(threads),
                "--tree-builder",
                "iqtree-fast",
                "--first-tree-builder",
                "rapidnj",
                "--first-model",
                "JC",
                "--model",
                "GTR",
                str(alignment),
            ],
            log=lineage_dir / "logs" / "gubbins.log",
            expected=tree,
            force=force,
        )
        _run_command(
            [
                "mask_gubbins_aln.py",
                "--aln",
                str(alignment),
                "--gff",
                str(gff),
                "--out",
                str(masked),
            ],
            log=lineage_dir / "logs" / "mask_recombination.log",
            expected=masked,
            force=force,
        )

        sequence_length = assembly_length(reference.assembly)
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
                seed=seed,
                output=temporal_path,
            )
        else:
            temporal = json.loads(temporal_path.read_text(encoding="utf-8"))

        observed = temporal["observed"]
        temporal_supported = bool(
            observed["rate"] > 0
            and temporal["successful_randomisations"] == randomisations
            and temporal["p_value_r_squared"] <= temporal_p_value
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
                    "--outdir",
                    str(lineage_dir / "timetree"),
                ],
                log=lineage_dir / "logs" / "timetree.log",
                expected=time_tree,
                force=force,
            )

        location_tree = time_tree if time_tree.exists() else tree
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
            "sequence_length": sequence_length,
            "temporal_signal_supported": temporal_supported,
            "temporal_signal": temporal,
            "outputs": {
                "alignment": str(alignment),
                "tree": str(tree),
                "recombination_gff": str(gff),
                "masked_alignment": str(masked),
                "timetree": str(time_tree) if time_tree.exists() else None,
                "location_tree": str(location_output),
            },
        }
        (lineage_dir / "report.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        reports.append(report)

    summary = {
        "workflow": "beyondmlst",
        "lineages": reports,
        "guardrail": (
            "Location-state reconstructions are exploratory and do not by themselves establish "
            "direct transmission or a definitive number of introductions."
        ),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary
