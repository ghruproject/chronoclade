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
    profiles = "profiles"
    lineage = "lineage"
    st = "st"


class DistanceMethod(str, Enum):
    alleles = "alleles"
    direct_dna = "direct-dna"


def cgmlst_command(
    dataset: Annotated[Path, typer.Argument(help="Prepared dataset.json")],
    output: Annotated[Path, typer.Option("--out", "--output", "-o")] = Path("cgmlst_results"),
    partition_by: Annotated[PartitionBy, typer.Option(help="Compatible species/profile blocks (default); explicit lineage blocks or legacy input-only ST blocks")] = PartitionBy.profiles,
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
    context_max_candidates: Annotated[int | None, typer.Option(min=1, help="Refuse a public same-ST pool above this size before detailed retrieval; never truncate it")] = None,
    context_preflight: Annotated[bool, typer.Option(help="Preview all public species/ST counts and cost without details, profiles, assemblies or analysis; write a plan to --out")] = False,
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
    """Analyse compatible cgMLST profiles with optional lineage annotations."""
    from chronoclade.cgmlst.workflow import run_cgmlst
    from chronoclade.errors import WorkflowError
    from chronoclade.pathogenwatch import PathogenwatchError

    console = Console()
    try:
        if partition_by == PartitionBy.st and (fetch_context or context_preflight or context_source or public_typing or cglin_export or typing_config or context_lin_level is not None or hiercc_level is not None):
            raise ValueError("ST profile-only partition cannot fetch context or use lineage/context options")
        if partition_by == PartitionBy.profiles and (context_lin_level is not None or hiercc_level is not None):
            raise ValueError("Lineage-scope options require --partition-by lineage")
        if partition_by == PartitionBy.profiles and lin_level != 5:
            raise ValueError("--lin-level filtering requires --partition-by lineage")
        if context_max_candidates is not None and not (fetch_context or context_preflight):
            raise ValueError("--context-max-candidates requires --fetch-context or --context-preflight")
        source = dataset.expanduser()
        if context_preflight:
            if context_source or public_typing or cglin_export:
                raise ValueError("Count-only context preflight cannot use frozen context or typing imports")
            from chronoclade.context_preflight import plan_context, write_context_plan, require_allowed_pools
            plan = plan_context(source, max_candidates=context_max_candidates)
            plan_path = write_context_plan(plan, output)
            console.print(f"Context preflight: {plan_path}")
            cost = plan['cost']
            console.print(f"{cost['public_candidate_records']} public candidate records; "
                          f"up to {cost['cgmlst_export_batches_per_analysis_upper_bound']} export batches per cgMLST analysis; "
                          f"one dense distance matrix upper bound: {cost['one_float64_distance_matrix_bytes_upper_bound']:,} bytes")
            if not plan['download_credentials_available']:
                console.print("Protected profile-download credentials are unavailable; this public count preview does not require them.")
            require_allowed_pools(plan)
            return
        if fetch_context or context_source:
            from chronoclade.profile_context_provider import discover_profile_context
            prepared = discover_profile_context(
                source, output.expanduser().parent / (output.name + "_context"),
                lin_level=context_lin_level if context_lin_level is not None else lin_level, hiercc_level=hiercc_level,
                public_typing=public_typing, cglin_export=cglin_export, catalogues=catalogues,
                typing_config=typing_config, context_source=context_source,
                resume=resume,
                lineage_filter=partition_by == PartitionBy.lineage,
                max_candidates=context_max_candidates,
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
    except (ValueError, OSError, WorkflowError, PathogenwatchError) as error:
        console.print(f"Error: {error}", style="bold red", markup=False)
        if fetch_context and isinstance(error, PathogenwatchError):
            console.print("Public context cache retained. For interrupted requests, repeat the same "
                          "command with --resume. Resolve count or integrity failures before retrying.",
                          markup=False)
        raise typer.Exit(code=2) from error
    console.print(f"cgMLST manifest: {result.manifest_path}")
    console.print(f"cgMLST report: {result.report_path}")
    import json
    status = json.loads(result.manifest_path.read_text()).get("status")
    if status != "complete":
        console.print(f"Analysis {status}: inspect unresolved inputs and report before continuing", style="bold red")
        raise typer.Exit(code=2)
