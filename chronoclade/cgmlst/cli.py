"""Arguments for the independent categorical-profile stage."""

from enum import Enum
from pathlib import Path
from typing import Annotated
import typer
from rich.console import Console


class AnalysisBackend(str, Enum):
    auto = "auto"
    dense = "dense"
    scale = "scale"


class PartitionBy(str, Enum):
    lineage = "lineage"
    st = "st"


class DistanceMethod(str, Enum):
    alleles = "alleles"
    direct_dna = "direct-dna"


def cgmlst_command(
    dataset: Annotated[Path, typer.Argument(help="Prepared dataset.json")],
    output: Annotated[Path, typer.Option("--out", "--output", "-o")] = Path("cgmlst_results"),
    partition_by: Annotated[PartitionBy, typer.Option(help="Lineage blocks (default), or supplied-input-only species/ST/profile blocks without lineage assignment or context retrieval")] = PartitionBy.lineage,
    lin_level: Annotated[int, typer.Option(min=5, max=10, help="Reported LIN prefix depth for analysis blocks and default context retrieval; not a distance cutoff")] = 5,
    context_lin_level: Annotated[int | None, typer.Option(min=5, max=10, help="Optional context-pool union of input-associated LIN prefixes, independent of analysis block depth")] = None,
    display_lin_level: Annotated[int, typer.Option(min=5, max=10, help="Independent reported LIN prefix depth for tree/ordination colours; does not redefine nearest neighbours")] = 7,
    hiercc_level: Annotated[str | None, typer.Option(help="Explicit E. coli HC level")] = None,
    backend: Annotated[AnalysisBackend, typer.Option(help="auto preserves size-based choice; scale explicitly uses compiled distances and RapidNJ")] = AnalysisBackend.auto,
    subsamples: Annotated[int, typer.Option(min=1)] = 1,
    context_size: Annotated[int, typer.Option(min=0)] = 50,
    representative_reserve: Annotated[int, typer.Option(min=0, help="Minimum optional representative room beyond mandatory neighbours")] = 0,
    nearest_per_query: Annotated[int, typer.Option(min=1)] = 3,
    include: Annotated[
        list[str] | None, typer.Option(help="Pin a context ID; repeat for several")
    ] = None,
    resume: Annotated[bool, typer.Option(help="Reuse saved numerical analysis after packaging failure with identical source and parameters")] = False,
    fetch_context: Annotated[bool, typer.Option(help="Retrieve matching public profiles before analysis")] = False,
    context_source: Annotated[Path | None, typer.Option(help="Replay a validated frozen whole-ST raw provider source without network; currently requires a Klebsiella native reference")] = None,
    public_typing: Annotated[Path | None, typer.Option(help="Public profile export")] = None,
    cglin_export: Annotated[Path | None, typer.Option(help="Compatible lineage export")] = None,
    catalogues: Annotated[Path | None, typer.Option(help="Explicit canonical locus catalogues")] = None,
    typing_config: Annotated[Path | None, typer.Option(help="Native reference configuration for canonical public-profile normalization; defaults to the verified prepare reference")] = None,
    seed: int = 42,
    bootstrap_replicates: Annotated[int | None, typer.Option(min=0, max=200, help="Default: 30 for dense blocks, 0 for large blocks; recorded per block")] = None,
    min_overlap: Annotated[float, typer.Option(min=0, max=1)] = 0.9,
    distance_threshold: Annotated[float | None, typer.Option(min=0, max=1, help="Cutoff in the selected mismatch fraction; default .02 for alleles, required for direct-dna")] = None,
    distance_method: Annotated[DistanceMethod, typer.Option(help="Allele profiles (default) or direct nucleotide comparisons within cgMLST genes")] = DistanceMethod.alleles,
    allele_catalogue: Annotated[Path | None, typer.Option(help="Scope-checked frozen allele DNA catalogue required for direct-dna")] = None,
):
    """Analyse every matching profile in explicit frozen lineage blocks."""
    from chronoclade.cgmlst.workflow import run_cgmlst
    from chronoclade.errors import WorkflowError

    console = Console()
    try:
        if partition_by == PartitionBy.st and (fetch_context or context_source or public_typing or cglin_export or typing_config or context_lin_level is not None or hiercc_level is not None):
            raise ValueError("ST profile-only partition cannot fetch context or use lineage/context options")
        source = dataset.expanduser()
        if fetch_context or context_source:
            from chronoclade.profile_context_provider import discover_profile_context
            prepared = discover_profile_context(
                source, output.expanduser().parent / (output.name + "_context"),
                lin_level=context_lin_level if context_lin_level is not None else lin_level, hiercc_level=hiercc_level,
                public_typing=public_typing, cglin_export=cglin_export, catalogues=catalogues,
                typing_config=typing_config, context_source=context_source,
                resume=resume,
            )
            source = prepared.dataset_manifest
        result = run_cgmlst(
            source,
            output.expanduser(),
            lin_level=lin_level,
            partition_by=partition_by.value,
            display_lin_level=display_lin_level,
            context_lin_level=context_lin_level,
            backend=backend.value,
            hiercc_level=hiercc_level,
            resume=resume,
            representative_reserve=representative_reserve,
            replicates=subsamples,
            context_size=context_size,
            nearest_per_query=nearest_per_query,
            include=include,
            seed=seed,
            bootstrap_replicates=bootstrap_replicates,
            min_overlap=min_overlap,
            distance_threshold=distance_threshold,
            distance_method=distance_method.value,
            allele_catalogue=allele_catalogue.expanduser() if allele_catalogue else None,
        )
    except (ValueError, OSError, WorkflowError) as error:
        console.print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(code=2) from error
    console.print(f"cgMLST manifest: {result.manifest_path}")
    console.print(f"cgMLST report: {result.report_path}")
