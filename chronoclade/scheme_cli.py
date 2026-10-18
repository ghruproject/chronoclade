"""Public allele reference discovery and frozen installation commands."""
from pathlib import Path
from typing import Annotated

from rich.console import Console
from rich.table import Table
import typer

from chronoclade.enterobase_schemes import SCHEMES, inspect_scheme, install_scheme

console = Console()


def list_schemes_command(
    inspect: Annotated[bool, typer.Option('--inspect', help='Check public indexes now; does not download allele files')] = False,
):
    """List registered EnteroBase references; availability is separate from tested support."""
    table = Table(title='EnteroBase cgMLST allele references')
    for name in ('Scheme', 'Organism group', 'Public index', 'Loci'):
        table.add_column(name)
    failed = False
    for scheme, organism in SCHEMES.items():
        status, loci = 'not checked', '—'
        if inspect:
            try:
                snapshot = inspect_scheme(scheme)
                status, loci = 'available', str(len(snapshot['loci']))
            except (ValueError, OSError) as error:
                status, failed = f'unavailable: {error}', True
        table.add_row(scheme, organism, status, loci)
    console.print(table)
    console.print('Allele references do not include verified HierCC assignment or species acceptance.')
    if failed:
        raise typer.Exit(2)


def setup_scheme_command(
    scheme: Annotated[str, typer.Argument(help='Exact registered EnteroBase scheme ID')],
    destination: Annotated[Path, typer.Option('--destination', help='Dedicated reference directory')],
    workers: Annotated[int, typer.Option(min=1, max=8)] = 4,
):
    """Download and verify one public allele catalogue; safely resume the same snapshot."""
    try:
        manifest = install_scheme(scheme, destination, workers=workers)
    except (ValueError, OSError) as error:
        console.print(f'Error: {error}', markup=False)
        raise typer.Exit(2) from error
    console.print(f"Verified {len(manifest['loci'])} loci: {destination.expanduser().resolve()/'reference.json'}", markup=False)
    console.print('Alleles ready. Native caller indexing and HierCC assignment are separate setup steps.')
