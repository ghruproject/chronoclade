"""Command-line interface for beyondMLST."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from beyondmlst import __version__
from beyondmlst.metadata import MetadataError, read_metadata
from beyondmlst.workflow import (
    WorkflowError,
    native_platform_supported,
    plan,
    run_workflow,
    tool_status,
)

app = typer.Typer(
    name="beyondmlst",
    help="Recombination-aware temporal and contextual bacterial phylogenetics.",
    no_args_is_help=True,
)
console = Console()


def _fail(error: Exception) -> None:
    console.print(f"[bold red]Error:[/bold red] {error}")
    raise typer.Exit(code=2)


@app.command()
def version() -> None:
    """Print the installed version."""

    console.print(__version__)


@app.command()
def preflight() -> None:
    """Check that external workflow tools are available."""

    table = Table(title="beyondMLST preflight")
    table.add_column("Tool")
    table.add_column("Status")
    ok = True
    for tool, path in tool_status().items():
        if path:
            table.add_row(tool, f"[green]{path}[/green]")
        else:
            table.add_row(tool, "[red]missing[/red]")
            ok = False
    if native_platform_supported():
        table.add_row("native platform", "[green]supported[/green]")
    else:
        table.add_row("native platform", "[red]use Docker (Linux container)[/red]")
        ok = False
    console.print(table)
    if not ok:
        raise typer.Exit(code=2)


@app.command()
def validate(
    metadata: Annotated[Path, typer.Argument(help="Input metadata CSV")],
    min_samples: Annotated[
        int, typer.Option("--min-samples", min=3, help="Minimum genomes per lineage")
    ] = 10,
) -> None:
    """Validate files and show the lineage analysis plan."""

    try:
        samples = read_metadata(metadata)
        items = plan(samples, min_samples=min_samples)
    except (MetadataError, WorkflowError) as error:
        _fail(error)
    console.print(f"[green]Validated {len(samples)} samples across {len(items)} lineages.[/green]")
    _print_plan(items)


def _print_plan(items: list[dict[str, object]]) -> None:
    table = Table(title="Lineage plan")
    for heading in ("Species", "Lineage", "Samples", "Dates", "Reference", "Status"):
        table.add_column(heading)
    for item in items:
        table.add_row(
            str(item["species"]),
            str(item["lineage"]),
            str(item["sample_count"]),
            str(item["distinct_dates"]),
            str(item["reference"]),
            str(item["status"]),
        )
    console.print(table)


@app.command()
def run(
    metadata: Annotated[Path, typer.Argument(help="Input metadata CSV")],
    output: Annotated[
        Path, typer.Option("--output", "-o", help="Output directory")
    ] = Path("beyondmlst_results"),
    threads: Annotated[int, typer.Option("--threads", "-t", min=1)] = 4,
    date_randomisations: Annotated[
        int,
        typer.Option(
            "--date-randomisations",
            min=0,
            help="Number of tip-date permutations; use at least 100 for analysis",
        ),
    ] = 100,
    temporal_p_value: Annotated[
        float, typer.Option("--temporal-p-value", min=0.0, max=1.0)
    ] = 0.05,
    min_samples: Annotated[int, typer.Option("--min-samples", min=3)] = 10,
    seed: Annotated[int, typer.Option("--seed")] = 20260818,
    force: Annotated[bool, typer.Option("--force", help="Rerun completed stages")] = False,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Validate and print the plan only")
    ] = False,
) -> None:
    """Run lineage-specific recombination, temporal and location analyses."""

    try:
        samples = read_metadata(metadata)
        items = plan(samples, min_samples=min_samples)
        if dry_run:
            _print_plan(items)
            console.print(json.dumps(items, indent=2))
            return
        summary = run_workflow(
            samples,
            output=output,
            threads=threads,
            randomisations=date_randomisations,
            temporal_p_value=temporal_p_value,
            min_samples=min_samples,
            seed=seed,
            force=force,
        )
    except (MetadataError, WorkflowError, OSError) as error:
        _fail(error)
    supported = sum(
        bool(item.get("temporal_signal_supported")) for item in summary["lineages"]
    )
    console.print(
        f"[green]Completed {len(summary['lineages'])} lineage records; "
        f"{supported} passed the temporal-signal gate.[/green]"
    )
    console.print(f"Summary: {(output.expanduser().resolve() / 'summary.json')}")


if __name__ == "__main__":
    app()
