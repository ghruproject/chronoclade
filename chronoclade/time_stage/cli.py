"""Lazy independent dating and ensemble-comparison dispatch."""

from pathlib import Path
from typing import Annotated
import typer
from rich.console import Console


def register(app):
    app.command("time")(time_command)
    app.command("time-compare")(comparison_command)


def time_command(
    tree: Annotated[Path, typer.Argument(help="Validated tree.json or tree_ensemble.json")],
    output: Annotated[Path, typer.Option("--out", "--output", "-o")] = Path("time_results"),
    force: bool = False,
    allow_unsupported: bool = False,
    ensemble_jobs: Annotated[int, typer.Option(min=1)] = 1,
    cluster_snp_cutoff: Annotated[int, typer.Option(min=0)] = 10,
    reuse_runs: Annotated[list[Path] | None, typer.Option("--reuse-run", help="Checked time.json for this ensemble; repeat")] = None,
):
    """Date the saved tree when its temporal assessment supports dating."""
    from chronoclade.time_stage.workflow import run_time
    from chronoclade.errors import WorkflowError
    from chronoclade.selection_manifest import _read

    try:
        runner = run_time
        options = {}
        if _read(tree).get("schema") == "chronoclade.tree-ensemble":
            from chronoclade.execution.ensemble import run_time_ensemble
            runner = run_time_ensemble
            options.update(ensemble_jobs=ensemble_jobs, cluster_snp_cutoff=cluster_snp_cutoff,
                           reuse_results=reuse_runs or ())
        elif reuse_runs:
            raise WorkflowError("--reuse-run requires a tree ensemble")
        result = runner(tree, output=output, force=force, allow_unsupported=allow_unsupported, **options)
    except (ValueError, OSError, WorkflowError) as error:
        Console().print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(2) from error
    if result.get("schema") == "chronoclade.time-ensemble":
        Console().print(f"Time ensemble: {output.resolve() / 'time_ensemble.json'}")
        Console().print(f"Completed: {result['completed_runs']}; failed: {result['failed_runs']}")
        Console().print(f"Comparison: {output.resolve() / result['report']['path']}")
    else:
        Console().print(f"Time manifest: {output.resolve() / 'time.json'}")
        Console().print(f"Dating status: {result['dating_status']}")


def comparison_command(
    ensemble: Annotated[Path, typer.Argument(help="Selection ensemble.json")],
    runs: Annotated[
        list[Path], typer.Option("--run", help="Complete time.json; repeat for each run")
    ],
    output: Annotated[Path, typer.Option("--out", "-o")] = Path("time_comparison"),
    cluster_snp_cutoff: Annotated[int, typer.Option(min=0)] = 10,
):
    """Compare independent selection dates/rates and shared-target clustering."""
    from chronoclade.time_stage.comparison import compare_time_runs
    from chronoclade.errors import WorkflowError

    try:
        compare_time_runs(ensemble, runs, output=output, cluster_snp_cutoff=cluster_snp_cutoff)
    except (ValueError, OSError, WorkflowError) as error:
        Console().print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(2) from error
    Console().print(f"Comparison: {output.resolve() / 'report.html'}")
