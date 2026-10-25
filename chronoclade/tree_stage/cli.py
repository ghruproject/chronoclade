"""Lazy dispatch for the independent assembly-tree command."""

from pathlib import Path
from typing import Annotated
import os
import typer
from rich.console import Console
from enum import Enum


class CoherencePolicy(str, Enum):
    strict = "strict"
    diagnostic = "diagnostic"


class RecombinationMethod(str, Enum):
    phipack = "phipack"
    clonalframeml = "clonalframeml"
    verticall = "verticall"


def register(app):
    app.command("tree")(tree_command)


def tree_command(
    selection: Annotated[Path, typer.Argument(help="Validated selection or ensemble manifest")],
    output: Annotated[Path, typer.Option("--out", "--output", "-o")] = Path("tree_results"),
    threads: Annotated[int, typer.Option(min=1)] = 1,
    reference_sample_id: Annotated[str | None, typer.Option(help="Exact selected input genome to use as reference; default highest N50 input")] = None,
    recombination_method: Annotated[RecombinationMethod, typer.Option(help="phipack (fast default), clonalframeml, or verticall reference-relative masking")] = RecombinationMethod.phipack,
    phipack_window: Annotated[int, typer.Option(min=200, help="PHIPack scanning window in mapped bases")] = 10000,
    phipack_step: Annotated[int, typer.Option(min=1, help="PHIPack scanning step in mapped bases")] = 5000,
    phipack_threshold: Annotated[float, typer.Option(min=0.000001, max=1, help="PHIPack Benjamini-Hochberg FDR threshold")] = 0.05,
    support_replicates: Annotated[int, typer.Option(min=0, help="Optional filtered-alignment UFBoot2 replicates: 0 disables; otherwise at least 1000. Separate supported ML tree; dating topology unchanged.")] = 0,
    randomisations: Annotated[int, typer.Option(min=0, help="Clustered date permutations; default 100. Smaller quick screens have coarser empirical p-values.")] = 100,
    randomisation_jobs: Annotated[int, typer.Option(min=1)] = 1,
    ensemble_jobs: Annotated[int, typer.Option(min=1)] = 1,
    temporal_p_value: Annotated[float, typer.Option(min=0.000001, max=1)] = 0.05,
    coherence_policy: Annotated[CoherencePolicy, typer.Option(help="Raw-distance flags trigger review; strict aborts on recombination-filtered flags; diagnostic retains all with audits")] = CoherencePolicy.strict,
    seed: int = 42,
    assess_temporal: bool = True,
    force: bool = False,
    api_key_env: str = "PATHOGENWATCH_API_KEY",
    base_url: str = "https://pathogen.watch",
    assembly_cache: Annotated[Path | None, typer.Option(help="Shared verified assembly cache")] = None,
):
    """Build the exact selected assembly phylogeny and clustered temporal screen."""
    from chronoclade.tree_stage.workflow import run_tree
    from chronoclade.errors import WorkflowError
    from chronoclade.selection_manifest import _read

    try:
        runner = run_tree
        options = {}
        if reference_sample_id is not None:
            options["reference_sample_id"] = reference_sample_id
        if assembly_cache is not None:
            options["assembly_cache"] = assembly_cache
        if _read(selection).get("schema") == "chronoclade.selection-ensemble":
            from chronoclade.execution.ensemble import run_tree_ensemble
            runner = run_tree_ensemble
            options["ensemble_jobs"] = ensemble_jobs
        result = runner(
            selection,
            output=output,
            threads=threads,
            support_replicates=support_replicates,
            randomisations=randomisations,
            randomisation_jobs=randomisation_jobs,
            temporal_p_value=temporal_p_value,
            seed=seed,
            assess_temporal=assess_temporal,
            force=force,
            coherence_policy=coherence_policy.value,
            recombination_method=recombination_method.value,
            phipack_window=phipack_window,
            phipack_step=phipack_step,
            phipack_threshold=phipack_threshold,
            api_key=_api_key(api_key_env),
            base_url=base_url,
            **options,
        )
    except (ValueError, OSError, WorkflowError) as error:
        Console().print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(2) from error
    if result.get("schema") == "chronoclade.tree-ensemble":
        Console().print(f"Tree ensemble: {output.resolve() / 'tree_ensemble.json'}")
        Console().print(f"Completed: {result['completed_runs']}; failed: {result['failed_runs']}")
        if result["failed_runs"]:
            for run in result.get("runs", []):
                if run.get("status") == "failed":
                    Console().print(
                        f"Failed {run['selection_id']}: {run.get('reason', 'See ensemble manifest')}",
                        style="bold red", markup=False,
                    )
                    break
            raise typer.Exit(2)
    else:
        Console().print(f"Tree manifest: {output.resolve() / 'tree.json'}")
        Console().print(f"Temporal assessment: {result['temporal_assessment']['code']}")
        if (result.get("coherence_screening") or {}).get("warning"):
            Console().print(result["coherence_screening"]["warning"], style="bold yellow", markup=False)

    from chronoclade.storage import completed_output
    completed_output(output)


def _api_key(environment):
    from chronoclade.pathogenwatch import load_api_key

    return os.environ.get(environment) or load_api_key()
