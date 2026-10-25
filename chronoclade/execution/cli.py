"""Supporting job preparation; submission remains explicit."""
from pathlib import Path
from typing import Annotated
from enum import Enum
import typer


def job_command(stage: str, arguments: Annotated[list[str], typer.Argument()],
                output: Annotated[Path, typer.Option('--out', '-o')] = Path('stage_job'),
                slurm: bool = False, cpus: int = 4, memory_gb: int = 16,
                walltime: str = '04:00:00', account: str | None = None,
                partition: str | None = None, gpus: int = 0,
                executable: str = 'chronoclade'):
    """Write an exact local/SLURM command; use -- before stage flags."""
    from .jobs import JobSpec, SlurmResources, write_job
    try:
        from chronoclade.storage import current_policy
        spec = JobSpec(stage, tuple(arguments), Path.cwd(), executable, current_policy())
        resources = SlurmResources(cpus, memory_gb, walltime, account, partition, gpus) if slurm else None
        path = write_job(spec, output, resources=resources)
    except (ValueError, OSError) as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(2) from exc
    typer.echo(f'Prepared job: {path}')


class EnsembleStage(str, Enum):
    tree = 'tree'
    time = 'time'


def slurm_ensemble_command(
    stage: EnsembleStage,
    source: Annotated[Path, typer.Argument(help='Selection ensemble.json or tree_ensemble.json')],
    output: Annotated[Path, typer.Option('--out', '-o', help='Fresh job directory')] = Path('ensemble_job'),
    results: Annotated[Path, typer.Option(help='Fresh scientific results directory')] = Path('ensemble_results'),
    max_parallel: Annotated[int, typer.Option(min=1)] = 4,
    cpus: Annotated[int, typer.Option(min=1, help='CPUs per array task')] = 4,
    memory_gb: Annotated[int, typer.Option(min=1, help='Memory per array task and collector')] = 16,
    walltime: str = '04:00:00',
    account: str | None = None,
    partition: str | None = None,
    threads: Annotated[int | None, typer.Option(min=1, help='Tree threads; defaults to cpus')] = None,
    randomisation_jobs: Annotated[int | None, typer.Option(min=1, help='Tree permutation workers; defaults to cpus')] = None,
    randomisations: Annotated[int, typer.Option(min=0)] = 100,
    reference_sample_id: Annotated[str | None, typer.Option(help='Tree reference: exact selected input sample ID')] = None,
    recombination_method: Annotated[str, typer.Option(help='Tree filtering: phipack, clonalframeml or verticall')] = 'phipack',
    phipack_window: Annotated[int, typer.Option(min=200)] = 10000,
    phipack_step: Annotated[int, typer.Option(min=1)] = 5000,
    phipack_threshold: Annotated[float, typer.Option(min=0.000001, max=1)] = 0.05,
    support_replicates: Annotated[int, typer.Option(min=0)] = 0,
    seed: int = 42,
    temporal_p_value: Annotated[float, typer.Option(min=0.000001, max=1)] = 0.05,
    coherence_policy: Annotated[str, typer.Option(help='strict or diagnostic; same policy as tree')] = 'strict',
    assess_temporal: bool = True,
    allow_unsupported: bool = False,
    cluster_snp_cutoff: Annotated[int, typer.Option(min=0)] = 10,
):
    """Prepare parallel SLURM selections and a dependent ensemble collector."""
    from .jobs import SlurmResources
    from .slurm import write_slurm_ensemble
    from chronoclade.errors import WorkflowError
    options = (dict(threads=threads or cpus, randomisation_jobs=randomisation_jobs or cpus,
                    randomisations=randomisations, support_replicates=support_replicates, seed=seed,
                    temporal_p_value=temporal_p_value, coherence_policy=coherence_policy,
                    assess_temporal=assess_temporal, recombination_method=recombination_method,
                    phipack_window=phipack_window, phipack_step=phipack_step,
                    phipack_threshold=phipack_threshold)
               if stage == EnsembleStage.tree else dict(cluster_snp_cutoff=cluster_snp_cutoff, allow_unsupported=allow_unsupported))
    if reference_sample_id is not None:
        if stage != EnsembleStage.tree:
            raise typer.BadParameter('--reference-sample-id is only valid for tree arrays')
        options['reference_sample_id'] = reference_sample_id
    try:
        path = write_slurm_ensemble(
            stage.value, source, output=output, results=results,
            resources=SlurmResources(cpus, memory_gb, walltime, account, partition),
            max_parallel=max_parallel, options=options,
        )
    except (ValueError, OSError, WorkflowError) as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(2) from exc
    typer.echo(f'Prepared array: {path}')
    typer.echo(f'Submit: bash {output.resolve() / "submit.sh"}')
    typer.echo(f'Up to {max_parallel} tasks, each requesting {cpus} CPUs and {memory_gb} GiB')


def slurm_prepare_command(
    source: Annotated[Path, typer.Argument(help='Assembly CSV; other prepare input kinds are not supported')],
    typing_config: Annotated[Path, typer.Option(help='Validated native caller configuration')],
    output: Annotated[Path, typer.Option('--out', '-o', help='Fresh job directory')] = Path('prepare_jobs'),
    results: Annotated[Path, typer.Option(help='Fresh final prepared dataset directory')] = Path('prepared'),
    max_parallel: Annotated[int, typer.Option(min=1)] = 4,
    cpus: Annotated[int, typer.Option(min=1, help='CPUs per assembly; typing callers use one thread')] = 1,
    memory_gb: Annotated[int, typer.Option(min=1)] = 4,
    walltime: str = '01:00:00',
    account: str | None = None,
    partition: str | None = None,
    species: str | None = None,
    enrich_metadata: bool = False,
):
    """Prepare local assemblies concurrently, then collect the complete original cohort."""
    from .jobs import SlurmResources
    from .prepare_slurm import write_slurm_prepare
    try:
        plan = write_slurm_prepare(
            source, typing_config=typing_config, output=output, results=results,
            resources=SlurmResources(cpus, memory_gb, walltime, account, partition),
            max_parallel=max_parallel, species=species, enrich_metadata=enrich_metadata)
    except (OSError, ValueError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(2) from error
    typer.echo(f'Prepared parallel preparation: {plan}')
    typer.echo(f'Submit: bash {plan.parent / "submit.sh"}')
    typer.echo(f'Up to {max_parallel} assemblies concurrently; final dataset: {results.resolve()}')
