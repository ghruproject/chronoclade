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
7. Root-to-tip analysis and date randomisation using TreeTime.
8. Time-tree estimation only if the temporal-signal gate passes.
9. Exploratory discrete location-state reconstruction.
10. A reader-facing HTML report with TreeTime and date-randomisation visuals,
   explicit verdicts and links to the underlying evidence.

## Contextual-genome acquisition

`atbfetcher` is the acquisition layer for AllTheBacteria and RefSeq assemblies.
For an analysis of a focal ST, beyondMLST will:

1. Obtain the complete high-quality same-ST candidate list with
   `atbfetcher mlst-query`.
2. Intersect it with reviewable geography, collection-period, host and source
   strata from the same ATB SQLite snapshot used by `atbfetcher query`.
3. Fetch the candidate assemblies by accession.
4. Use a fast SKA distance screen to retain close neighbours of every focal
   isolate within the bounded downloaded pool, plus a reproducible stratified
   background across location and time.
5. Freeze the selected accessions in a context manifest before tree inference.

The manifest records the `atbfetcher` version, ATB metadata snapshot, MLST
scheme and ST, filters, random seed, accession, provenance, collection-date
precision, quality fields and the reason each genome was retained. A user may
provide the same manifest manually, so downloading the large ATB SQLite
database is optional rather than a requirement for every beyondMLST run.

The pinned `atbfetcher` release currently has a retired R2 URL for its MLST
table. beyondMLST first uses the normal command; if that cache is absent and the
query fails, it obtains the current official ATB `mlst.parquet` from OSF at the
cache path expected by `atbfetcher`, retries the command, and records the
fallback in the audit JSON.

The bounded pool does not guarantee retrieval of the globally nearest public
genomes when an ST contains thousands of records. Reports state this
explicitly. A whole-database sketch search seeded by each focal isolate is the
preferred extension before introduction-oriented conclusions are automated.

## Planned extensions

1. Validate the full run against one dominant *K. pneumoniae* lineage and one
   *E. coli* lineage from the GHRU data.
2. Validate the `atbfetcher` acquisition and down-selection defaults against
   large public *K. pneumoniae* and *E. coli* ST collections.
3. Export a combined Microreact project with metadata, tree and timeline.
4. Add sensitivity analyses for same-patient longitudinal isolates, contextual
   subsampling and alternative references.

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
