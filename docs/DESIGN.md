# Design and scope

The public-health question set and interpretation rules are recorded in
[PUBLIC_HEALTH_QUESTIONS.md](PUBLIC_HEALTH_QUESTIONS.md).

## Why a Python CLI with Pixi?

The target users already have a successful precedent in BactScout. A small
Typer-based Python application gives a clear user interface and makes metadata
validation, reporting and workflow decisions testable. Pixi supplies the
compiled dependencies from Bioconda without asking users to install each tool
manually.

A full workflow engine is not required. The Python runner uses explicit, logged
subprocesses and stage-level resume checks. It runs independent lineages and
TreeTime date permutations concurrently within one global CPU budget.
Multi-node scheduling can be reconsidered if real datasets demonstrate a need
for it, without changing the metadata contract or CLI.

## MVP boundary

The MVP starts from assembled genomes. Read QC and assembly belong upstream in
BactScout and GHRU-assembly. Contextual genomes may be acquired with
[`atbfetcher`](https://github.com/happykhan/atbfetcher), but beyondMLST does not
select them as an unreviewed random sample. The choice and balance of context
materially affect conclusions about introductions.

The pipeline performs:

1. Strict metadata and file validation.
2. Separation by confirmed species and lineage.
3. Explicit reference selection, or an N50-based fallback.
4. Reference-ordered whole-genome alignment using SKA2 directly.
5. An initial maximum-likelihood tree using IQ-TREE.
6. Recombination inference and corrected branch lengths using ClonalFrameML.
7. Pairwise clonal SNP and callable-site evidence from the filtered alignment.
8. Topology-defined candidate local groups, longitudinal summaries and a
   one-isolate-per-patient distance sensitivity view.
9. Root-to-tip analysis and date randomisation using TreeTime.
10. Time-tree estimation only if the temporal-signal gate passes.
11. Exploratory discrete location-state reconstruction.
12. A decision-first HTML report with a transparent scenario ledger,
    recommended follow-up, exact distance outputs and temporal diagnostics.

## Contextual-genome acquisition

The repository includes a compact, all-species AllTheBacteria metadata snapshot
for discovery and filtering. `atbfetcher` remains the assembly
acquisition layer. This avoids making every user install the full ATB metadata
database. For an analysis of a focal ST, beyondMLST will:

1. Obtain the complete high-quality same-ST candidate list from the versioned
   compact Parquet snapshot.
2. Intersect it with reviewable geography, collection-period, host and source
   strata stored in that snapshot.
3. Fetch the candidate assemblies by accession.
4. Use a fast SKA distance screen to retain close neighbours of every focal
   isolate within the bounded downloaded pool, plus a reproducible stratified
   background across location and time.
5. Freeze the selected accessions in a context manifest before tree inference.

The manifest records the `atbfetcher` version, ATB metadata snapshot, MLST
scheme and ST, filters, random seed, accession, provenance, collection-date
precision, quality fields and the reason each genome was retained. A user may
provide an updated compact table with `--metadata-table` or provide the final
context manifest manually.

The bundled 2025-05 snapshot contains all high-quality, downloadable ATB
bacterial assemblies with a perfect assigned ST, across every represented
species and MLST scheme. Its generator, source URLs, row counts, licence and
SHA-256 digest are versioned beside it. The full 27 GB SQLite database is needed
only by a maintainer when generating a new snapshot, not by analysts running
beyondMLST.

The bounded pool does not guarantee retrieval of the globally nearest public
genomes when an ST contains thousands of records. Reports state this
explicitly. A whole-database sketch search seeded by each focal isolate is the
preferred extension before provisional introduction-oriented conclusions are
promoted beyond moderate confidence. The current automated scenario is
therefore explicitly provisional.

## Planned extensions

1. Validate the full run against one dominant *K. pneumoniae* lineage and one
   *E. coli* lineage from the GHRU data.
2. Validate the `atbfetcher` acquisition and down-selection defaults against
   large public *K. pneumoniae* and *E. coli* ST collections.
3. Export a combined Microreact project with metadata, tree and timeline.
4. Extend the current one-isolate-per-patient distance view to full topology
   reruns, contextual subsampling and alternative-reference sensitivity.

## Import-detection models

DetectImports is not part of the core workflow. It is designed for settings
with dense sampling in a population where local transmission is common and
imports are exceptional. Those assumptions are unlikely to hold for routine
GHRU AMR surveillance. beyondMLST therefore uses contextual placement,
recombination-corrected distances and topology to describe evidence consistent
with persistence or multiple introductions, without labelling individual
samples as imports.

## Interpretation boundary

The workflow can identify patterns consistent with locally persistent clades
or several phylogenetically distinct local clusters embedded among contextual
genomes. The number of introduction events is not directly observed. It is a
model- and sampling-dependent inference and must be described with uncertainty.
