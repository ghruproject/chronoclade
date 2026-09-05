# ChronoClade

`chronoclade` is a reproducible bacterial genomic epidemiology workflow for
asking what happens **beyond an MLST assignment**: are closely related isolates
consistent with sustained local circulation, or with repeated introductions?

The command-line tool takes assemblies and sample metadata, analyses each
species/lineage separately, removes recombination, tests whether the data have
temporal signal, and only creates a dated phylogeny when that test is passed.
Local, retrospective, and public contextual genomes can be supplied in the
same run. A decision-first HTML report presents a provisional public-health
scenario and recommended follow-up first. Detailed distances, trees, temporal
diagnostics and audit files remain available in a collapsed technical section.

> [!IMPORTANT]
> A phylogeny is not a transmission tree. Location-state reconstruction is
> exploratory and is sensitive to uneven or incomplete contextual sampling.
> `chronoclade` reports the evidence; it does not automatically label individual
> transmission or introduction events.

## Workflow

```text
assemblies + metadata
        |
        v
validate and split by species/lineage
        |
        v
select/confirm a lineage reference
        |
        v
SKA2 reference-ordered whole-genome alignment
        |
        v
IQ-TREE starting maximum-likelihood phylogeny
        |
        v
ClonalFrameML recombination inference and corrected tree
        |
        +--> clonal SNP + pairwise callable-site matrices
        |         |
        |         +--> longitudinal/patient summaries + heatmap
        |
        +--> topology-defined candidate local groups
        |         |
        |         +--> transparent scenario evidence ledger
        |
        +--> TreeTime root-to-tip analysis + date randomisation
        |         |
        |         +--> dated tree only when temporal signal passes
        |
        +--> exploratory location-state reconstruction
```

## Current status

This repository contains an early working MVP. It validates inputs and
orchestrates installed SKA2, IQ-TREE, ClonalFrameML and TreeTime tools. It needs
validation on the first GHRU *Klebsiella* and *E. coli* datasets before a stable
release.

## Installation

[Pixi](https://pixi.sh/) installs the Python CLI and all compiled bioinformatics
tools together. The lock file covers Apple-silicon and Intel macOS plus x86_64
and arm64 Linux, so Docker is not required.

```bash
git clone https://github.com/ghruproject/chronoclade.git
cd chronoclade
pixi install
pixi run chronoclade preflight
```

The first install may take several minutes while SKA2, IQ-TREE, ClonalFrameML
and TreeTime are downloaded. Windows users can run the same Pixi commands under
WSL2; native Windows is not part of the MVP support matrix.

## Input metadata

One CSV row is required per assembly:

```csv
sample_id,assembly,collection_date,location,species,lineage,origin,is_reference,patient_id
KPN001,assemblies/KPN001.fasta,2023-03,KIMS,Klebsiella_pneumoniae,ST15,local,true,P001
KPN002,assemblies/KPN002.fasta,2024-01-17,KIMS,Klebsiella_pneumoniae,ST15,local,false,P002
KPN_REF,context/KPN_REF.fasta,2021,Philippines,Klebsiella_pneumoniae,ST15,context,false,
```

Required columns:

- `sample_id`: unique, tree-safe identifier containing letters, numbers, `.`,
  `_` or `-`
- `assembly`: FASTA path, resolved relative to the metadata file
- `collection_date`: `YYYY`, `YYYY-MM`, `YYYY-MM-DD`, `YYYY-XX-XX`, or a
  TreeTime numeric/range date
- `location`: the geographic or institutional state to reconstruct
- `species`: confirmed species or species-complex label
- `lineage`: lineage analysed as a unit, normally an ST or genomic cluster
- `origin`: `local`, `retrospective`, or `context`

Optional columns:

- `is_reference`: mark exactly one preferred reference per lineage with
  `true`; otherwise the assembly with the highest N50 is selected
- `patient_id`: coded identifier used for within/between-patient summaries and
  a deterministic one-isolate-per-patient sensitivity view; it is never used
  to infer direct transmission

Use de-identified metadata only. Do not commit patient-level metadata or
sequence data to this repository.

## Commands

Validate metadata without running external tools:

```bash
pixi run chronoclade validate metadata.csv
```

Preview the planned lineage analyses:

```bash
pixi run chronoclade run metadata.csv --dry-run
```

Run the workflow:

```bash
pixi run chronoclade run metadata.csv \
  --output chronoclade_results \
  --threads 8 \
  --lineage-jobs 2 \
  --randomisation-jobs 4 \
  --date-randomisations 100
```

### Prepare public context genomes

[`atbfetcher`](https://github.com/happykhan/atbfetcher) is installed in the
same Pixi environment. ChronoClade includes a compact AllTheBacteria 2025-05
metadata snapshot covering every bacterial species and MLST scheme with a
high-quality, downloadable genome and a perfect ST assignment. Prepare a
bounded same-ST context set directly:

The current snapshot is 32 MiB and contains 2,047,053 genomes across 942
species and 146 MLST schemes.

```bash
pixi run chronoclade prepare-context focal_metadata.csv \
  --scheme ecoli_achtman_4 \
  --st 131 \
  --output context/ST131 \
  --candidate-pool 500 \
  --max-context 150 \
  --nearest-per-focal 3 \
  --threads 8
```

No 27 GB SQLite database is required. Use `--dry-run` to freeze and review the
accession pool before downloading assemblies. Geography, year, host and
isolation-source filters are optional; repeat `--country` to retain several
country prefixes. `--metadata-table` can supply a newer or locally generated
snapshot without changing the workflow.

The command balances the same-ST candidate pool across country and year before
download, screens candidates against all focal isolates with `ska distance`,
and then retains nearest neighbours plus a stratified background. Continue the
analysis using the exact generated inputs:

```bash
pixi run chronoclade run context/ST131/combined_metadata.csv \
  --context-manifest context/ST131/context_manifest.tsv \
  --output chronoclade_results
```

The SKA distance is a fast candidate-selection measurement. It is not reported
as a recombination-corrected transmission threshold; final relatedness is
assessed from the full clonal analysis. “Nearest” means nearest among the
bounded downloaded screening pool, not necessarily nearest among every public
genome. Whole-database sketch-based neighbour retrieval is a post-MVP priority.

`--threads` is the total CPU budget. Independent lineages run concurrently up
to `--lineage-jobs`; the budget is divided between them. Within each lineage,
independent TreeTime permutations run concurrently up to
`--randomisation-jobs`. The effective values are capped so the workflow does
not request more CPUs than the global budget.

The date-randomisation test compares the observed root-to-tip fit with fits
obtained after permuting collection dates. A dated tree is generated only when
the observed rate is positive and the randomisation p-value is at or below the
configured threshold (default `0.05`). Automatic clock-outlier filtering is
disabled during this gate so that genomes are not silently removed based on
their temporal fit.

The public-health summary is separate from the temporal-signal gate. It uses
the recombination-filtered distances, corrected rooted topology, longitudinal
span, patient sensitivity and contextual placement to assign one cautious
working interpretation: persistent local lineage, multiple introductions,
mixed, or indeterminate. The evidence ledger exposes every input to that
interpretation. No universal SNP threshold is applied, and confidence is
capped at moderate until branch support is propagated through the corrected
tree and public-neighbour retrieval is exhaustive.

## Main outputs

Each lineage directory contains:

- `inputs.tsv` and `metadata.csv`: exact inputs used
- `core_alignment.fasta`: reference-ordered SKA2 alignment
- `iqtree.treefile`: starting maximum-likelihood phylogeny
- `clonalframeml.labelled_tree.newick`: recombination-corrected phylogeny
- `clonalframeml.importation_status.txt`: inferred recombination intervals by branch
- `clonalframeml.filtered.fasta`: alignment containing non-recombinant sites
- `clonal_pairwise_distances.tsv`: pairwise clonal SNPs, callable sites and
  longitudinal/patient comparison classes
- `clonal_snp_matrix.tsv` and `pairwise_callable_sites.tsv`: exact square
  matrices for audit and reuse
- `clonal_snp_heatmap.svg` and `.png`: corrected-distance heatmap for editing
  or direct use
- `public_health_evidence.json`: scenario, evidence ledger, candidate local
  groups, final contextual neighbours and patient sensitivity
- `clock/`: observed root-to-tip analysis
- `clock/root_to_tip_regression.svg`: TreeTime root-to-tip visual
- `root_to_tip.png`: high-resolution root-to-tip visual
- `temporal_signal.json`: observed and randomised temporal-signal statistics
- `date_randomisation.csv`: tidy observed and per-permutation statistics
- `date_randomisation.svg` and `.png`: observed R² and rate against permuted dates
- `timetree/`: dated tree and uncertainty, when supported
- `timetree_with_confidence.svg`: the dated phylogeny on a calendar axis with
  90% confidence intervals for inferred internal-node dates
- `timetree.png`: high-resolution raster companion to the confidence-interval
  phylogeny
- `timetree_confidence.csv` and `node_dates.csv`: clock-rate uncertainty,
  inferred root date, interval-width summaries and per-node 90% intervals
- `location/`: exploratory ancestral location-state reconstruction
- `context_manifest.tsv`: selected public genomes and their acquisition and
  screening provenance, when supplied
- `report.json`: machine-readable lineage report
- `report.html`: a stage-by-stage walkthrough from root-to-tip exploration to
  date randomisation, conditional time scaling and public-health interpretation
- `supporting_results.zip`: one portable bundle of the CSV/TSV, JSON, SVG, PNG
  and tree files used by the report

The top-level `index.html` links all lineage reports. `summary.json` records
every lineage, skipped analysis, command and output.

`prepare-context` additionally writes:

- `same_st_accessions.txt`: every high-quality same-ST accession discovered
- `candidate_pool.tsv`: the metadata-balanced pre-download pool
- `context_distances.tsv`: complete focal/candidate SKA screening distances
- `screened_candidates.tsv`: all successfully screened candidates
- `context_manifest.tsv`: the frozen selected context set and inclusion reasons
- `combined_metadata.csv`: focal plus selected context samples, ready for `run`
- `context_selection.json`: metadata snapshot, filters, versions and attrition

The bundled table and its source manifest are in `chronoclade/data/`. Maintainers
can regenerate it with `scripts/build_atb_context_snapshot.py` when a new ATB
release is adopted.

## Demonstration dataset

Generate two synthetic lineages, one with temporal signal and one with dates
deliberately shuffled, then run the complete native Pixi workflow:

```bash
pixi run python examples/demo/generate_demo.py --output demo_run/input
pixi run chronoclade run demo_run/input/metadata.csv \
  --output demo_run/results \
  --threads 8 \
  --lineage-jobs 2 \
  --randomisation-jobs 4 \
  --date-randomisations 20
```

Twenty permutations keep the demonstration quick; use at least 100 for real
analyses. Open `demo_run/results/index.html` for the combined report.

### Controlled public-health scenarios

Generate four lightweight reports that exercise the intended public-health
interpretations without running the external phylogenetic tools:

```bash
pixi run python examples/scenarios/generate_scenarios.py \
  --output scenario_reports
```

Open `scenario_reports/index.html` to compare persistent local lineage,
multiple introductions, mixed and indeterminate outputs. These fixtures have
supplied synthetic trees and alignments. They test scenario logic and report
presentation; they are not biological validation datasets and must not be
presented as real outbreaks.

## Scientific guardrails

- Analyse *Klebsiella* and *E. coli* separately and stratify them into
  sufficiently close lineages before alignment.
- Treat root-to-tip regression as exploratory; date randomisation is the
  workflow gate.
- Do not interpret a local-only tree as evidence against introductions.
- Balance contextual sampling across geography and time where possible.
- Treat context selection as part of the analysis: retain its manifest and
  report candidate attrition and focal-neighbour coverage.
- Review longitudinal isolates from the same patient explicitly; genomic
  proximity alone is not proof of direct transmission.
- Report TreeTime dates and rates only when the date-randomisation gate passes,
  and retain the complete clock diagnostics for review.

## Development

```bash
pixi run test
pixi run lint
```

See [docs/DESIGN.md](docs/DESIGN.md) for scope,
[docs/PUBLIC_HEALTH_QUESTIONS.md](docs/PUBLIC_HEALTH_QUESTIONS.md) for the
interpretation framework, and [docs/ROADMAP.md](docs/ROADMAP.md) for the path
from the working pipeline to an actionable public-health report.
