"""Inspect, compact and restore reusable analysis files."""
from pathlib import Path
from typing import Annotated
import json
import typer
from rich.console import Console
from chronoclade import storage

app = typer.Typer(help="Lossless storage management for idle analysis directories.")


def _run(action, root, **kwargs):
    try:
        rows = action(root.expanduser(), **kwargs)
    except (OSError, ValueError) as error:
        Console().print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(2) from error
    return rows


@app.command("pack")
def pack_command(
    directory: Path,
    codec: Annotated[str, typer.Option(help="zstd (large-window, default) or gzip")] = "zstd",
    min_size_mb: Annotated[float, typer.Option(min=0, help="Minimum logical file size in MiB")] = 1,
    dry_run: bool = False,
):
    """Compact large payloads after verifying an exact round trip."""
    rows = _run(storage.pack, directory, codec=codec, min_bytes=int(min_size_mb * 1024**2), dry_run=dry_run)
    if dry_run:
        for row in rows:
            typer.echo(f"{row['bytes'] / 1024**2:.2f} MiB\t{row['path']}")
        typer.echo(f"{len(rows)} candidates; {sum(r['bytes'] for r in rows) / 1024**2:.1f} MiB logical (preview only)")
        return
    saved = sum(row.get("bytes", 0) - row.get("stored_bytes", row.get("bytes", 0)) for row in rows)
    Console().print(f"{len(rows)} files; {saved / 1024**2:.1f} MiB saved" + (" (preview only)" if dry_run else ""))


@app.command("status")
def status_command(directory: Path, as_json: Annotated[bool, typer.Option("--json")] = False):
    """Show packed sizes and whether native/mmap copies are present."""
    rows = _run(storage.inspect, directory)
    if as_json:
        typer.echo(json.dumps(rows, indent=2))
    else:
        logical = sum(row["bytes"] for row in rows)
        stored = sum(row["stored_bytes"] for row in rows)
        materialized = sum(row["bytes"] for row in rows if row["materialized"])
        Console().print(f"{len(rows)} packed files: {logical / 1024**2:.1f} MiB logical; "
                        f"{stored / 1024**2:.1f} MiB compressed; {materialized / 1024**2:.1f} MiB restored")


@app.command("verify")
def verify_command(directory: Path):
    """Verify both compressed bytes and decompressed scientific content."""
    rows = _run(storage.inspect, directory, verify=True)
    Console().print(f"Verified {len(rows)} compressed files and their logical checksums")


@app.command("restore")
def restore_command(directory: Path):
    """Persist exact original files for other tools, retaining checked archives."""
    rows = _run(storage.restore, directory)
    Console().print(f"Restored {len(rows)} files")
