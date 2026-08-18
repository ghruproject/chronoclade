# Design and scope

## Why a Python CLI with Pixi?

The target users already have a successful precedent in BactScout. A small
Typer-based Python application gives a clear user interface and makes metadata
validation, reporting and workflow decisions testable. Pixi supplies the
compiled dependencies from Bioconda without asking users to install each tool
manually.

A full workflow engine is not required for the first pilot. The Python runner
uses explicit, logged subprocesses and stage-level resume checks. If datasets
later require multi-node scheduling, the same stages can be moved behind a
Nextflow implementation without changing the metadata contract or CLI.

## MVP boundary

The MVP starts from assembled genomes. Read QC and assembly belong upstream in
BactScout and GHRU-assembly. It does not automatically search NCBI because the
choice and balance of contextual genomes materially affect conclusions about
introductions.

The pipeline performs:

1. Strict metadata and file validation.
2. Separation by confirmed species and lineage.
3. Explicit reference selection, or an N50-based fallback.
4. Whole-genome alignment using Gubbins' SKA2 helper.
5. Recombination inference and a clonal-frame tree using Gubbins.
6. Root-to-tip analysis and date randomisation using TreeTime.
7. Time-tree estimation only if the temporal-signal gate passes.
8. Exploratory discrete location-state reconstruction.

## Planned extensions

1. Validate the full run against one dominant *K. pneumoniae* lineage and one
   *E. coli* lineage from the GHRU data.
2. Add a curated contextual-genome manifest and reproducible NCBI download
   command, keeping selection decisions reviewable.
3. Add BactDating as a confirmatory dating backend in a pinned container.
4. Export a combined Microreact project with metadata, tree and timeline.
5. Add Nextflow execution for Slurm and larger multi-lineage studies.
6. Add sensitivity analyses for same-patient longitudinal isolates, contextual
   subsampling and alternative references.

## Interpretation boundary

The workflow can identify patterns consistent with locally persistent clades
or several phylogenetically distinct local clusters embedded among contextual
genomes. The number of introduction events is not directly observed. It is a
model- and sampling-dependent inference and must be described with uncertainty.
