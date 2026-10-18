"""Command-line interface for ChronoClade."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from chronoclade import __version__
from chronoclade.prepare_cli import prepare_command
from chronoclade.cgmlst.cli import cgmlst_command
from chronoclade.cgmlst.refresh import refresh_command
from chronoclade.cgmlst.clock_comparison import compare_clocks_command
from chronoclade.cgmlst.calibration import calibrate_lin_command
from chronoclade.execution.cli import job_command
from chronoclade.export_cli import export_command
from chronoclade.scheme_cli import list_schemes_command, setup_scheme_command
from chronoclade.pasteur_authorization import authorize_pasteur_command
from chronoclade.tree_stage.cli import register as register_tree
from chronoclade.time_stage.cli import register as register_time


app = typer.Typer(
    name="chronoclade",
    help="Recombination-aware temporal and contextual bacterial phylogenetics.",
    no_args_is_help=True,
)
console = Console()

app.command("prepare")(prepare_command)
app.command("cgmlst")(cgmlst_command)
app.command("refresh-cgmlst-report")(refresh_command)
app.command("job")(job_command)
app.command("export")(export_command)
app.command("list-schemes")(list_schemes_command)
app.command("setup-scheme")(setup_scheme_command)
app.command("authorize-pasteur")(authorize_pasteur_command)
register_tree(app)
register_time(app)

app.command("compare-clocks")(compare_clocks_command)
app.command("calibrate-lin")(calibrate_lin_command)


@app.callback()
def command_options(
    quiet: Annotated[bool, typer.Option(help="Hide progress messages; retain summaries and errors")] = False,
):
    from chronoclade.progress import configure_progress
    configure_progress(quiet=quiet)


def _fail(error: Exception) -> None:
    console.print(f"[bold red]Error:[/bold red] {error}")
    raise typer.Exit(code=2)


@app.command()
def version() -> None:
    """Print the installed version."""

    console.print(__version__)


@app.command()
def preflight(
    stage: Annotated[str, typer.Option(help="Software checks for all, prepare, cgmlst, tree or time")] = "all",
) -> None:
    """Check stage software; references, credentials and inputs are checked separately."""

    from chronoclade.native_tools import native_platform_supported, tool_status

    requirements = {"all": None, "prepare": (), "cgmlst": (),
                    "tree": ("ska", "iqtree", "ClonalFrameML", "treetime"),
                    "time": ("treetime",)}
    if stage not in requirements:
        _fail(ValueError("Stage must be all, prepare, cgmlst, tree or time"))
    table = Table(title=f"ChronoClade software preflight · {stage}")
    table.add_column("Tool")
    table.add_column("Status")
    ok = True
    for tool, path in tool_status().items():
        if requirements[stage] is not None and tool not in requirements[stage]:
            continue
        if path:
            table.add_row(tool, f"[green]{path}[/green]")
        else:
            table.add_row(tool, "[red]missing[/red]")
            ok = False
    if stage in {"prepare", "cgmlst"}:
        table.add_row("Python profile workflow", "[green]available[/green]")
    elif native_platform_supported():
        table.add_row("native platform", "[green]supported[/green]")
    else:
        table.add_row("native platform", "[red]unsupported by the Pixi environment[/red]")
        ok = False
    console.print(table)
    console.print("Software check only. Typing references, credentials and input readiness are assessed separately.")
    if stage == "cgmlst":
        console.print("Large profile analyses additionally require RapidNJ; see the installation guide.")
    if not ok:
        raise typer.Exit(code=2)


@app.command("setup-typing")
def setup_typing_command(
    destination: Annotated[
        Path,
        typer.Option("--destination", help="Pinned native typing tools and reference databases"),
    ] = Path("~/.local/share/chronoclade/typing"),
    species: Annotated[
        list[str] | None,
        typer.Option("--species", help="klebsiella or ecoli; repeat to install both"),
    ] = None,
    prepare_db: Annotated[
        bool,
        typer.Option(
            "--prepare-db",
            help="Also download and index publicly available cgMLST alleles; may take hours",
        ),
    ] = False,
    pasteur_secrets: Annotated[
        Path | None,
        typer.Option(
            "--pasteur-secrets",
            help=(
                "Protected Pasteur client credentials for Klebsiella LIN builds and, "
                "with a cached authorization token, authenticated cgMLST downloads"
            ),
        ),
    ] = None,
    enterobase_key_file: Annotated[
        Path | None,
        typer.Option(
            "--enterobase-key-file",
            help="Protected EnteroBase key file for an E. coli HierCC database build",
        ),
    ] = None,
) -> None:
    """Install Pathogenwatch typing tools via uv and report database readiness."""
    from chronoclade.typing_setup import TypingSetupError, setup_typing

    try:
        path, report = setup_typing(
            destination,
            species=species or ["klebsiella", "ecoli"],
            prepare_db=prepare_db,
            pasteur_secrets=pasteur_secrets,
            enterobase_key_file=enterobase_key_file,
        )
    except (TypingSetupError, OSError) as error:
        _fail(error)
    console.print(f"Typing setup: {path}")
    console.print(json.dumps(report["databases"], indent=2))


@app.command("type-queries")
def type_queries_command(
    metadata: Annotated[Path, typer.Argument(help="Query genome metadata CSV")],
    typing_config: Annotated[
        Path,
        typer.Option(
            "--typing-config", help="Pinned tools and reference databases from setup-typing"
        ),
    ],
    output: Annotated[Path, typer.Option("--output", "-o")] = Path("chronoclade_typing"),
    profiles_only: Annotated[bool, typer.Option(help="Call cgMLST profiles only; do not assign LIN/HierCC or require its reference")] = False,
    resume: Annotated[bool, typer.Option("--resume", help="Resume an interrupted exact-identity typing run")] = False,
) -> None:
    """Call cgMLST and Klebsiella LIN or E. coli HierCC locally."""
    from chronoclade.metadata import MetadataError, read_metadata
    from chronoclade.query_typing import QueryTypingError, type_query_assemblies

    try:
        results = type_query_assemblies(
            read_metadata(metadata), typing_config, output.expanduser().resolve(),
            profiles_only=profiles_only,
            resume=resume,
        )
    except (QueryTypingError, MetadataError, OSError, ValueError) as error:
        checkpoint = output.expanduser().resolve()
        resume_command = f"chronoclade type-queries {metadata} --typing-config {typing_config} --output {checkpoint} --resume"
        _fail(f"{error}\nCheckpoint: {checkpoint}\nResume with: {resume_command}")
    console.print(
        f"Typed {len(results)} queries; results: {output.expanduser().resolve() / 'query_typing.json'}"
    )


@app.command("setup-scale")
def setup_scale(output: Annotated[Path | None, typer.Option("--out")] = None,
                jobs: Annotated[int, typer.Option(min=1)] = 2):
    """Install the pinned RapidNJ backend used for complete large profile trees."""
    from chronoclade.cgmlst.scale_setup import install_rapidnj
    try:
        executable = install_rapidnj(output, jobs=jobs)
    except (ValueError, OSError, RuntimeError) as error:
        _fail(error)
    console.print(f"RapidNJ backend: {executable}", markup=False)


if __name__ == "__main__":
    app()
