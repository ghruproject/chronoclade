# beyondMLST

`beyondmlst` is a reproducible bacterial genomic epidemiology workflow for
asking what happens **beyond an MLST assignment**: are closely related isolates
consistent with sustained local circulation, or with repeated introductions?

The command-line tool takes assemblies and sample metadata, analyses each
species/lineage separately, removes recombination, tests whether the data have
temporal signal, and only creates a dated phylogeny when that test is passed.
Local, retrospective, and public contextual genomes can be supplied in the
same run.

> [!IMPORTANT]
> A phylogeny is not a transmission tree. Location-state reconstruction is
> exploratory and is sensitive to uneven or incomplete contextual sampling.
> `beyondmlst` reports the evidence; it does not automatically label individual
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
Gubbins recombination inference and clonal-frame tree
        |
        +--> TreeTime root-to-tip analysis + date randomisation
        |         |
        |         +--> dated tree only when temporal signal passes
        |
        +--> exploratory location-state reconstruction
```

## Current status

This repository contains an early working MVP. It validates inputs and
orchestrates installed SKA2/Gubbins/TreeTime tools. It needs validation on the
first GHRU *Klebsiella* and *E. coli* datasets before a stable release.

## Installation

[Pixi](https://pixi.sh/) is the recommended native installation method on
Linux and HPC systems because it installs the Python CLI and compiled
bioinformatics tools together.

```bash
git clone https://github.com/ghruproject/beyondmlst.git
cd beyondmlst
pixi install
pixi run beyondmlst preflight
```

The first install may take several minutes while Gubbins, SKA2, IQ-TREE and
TreeTime are downloaded.

### macOS and Windows

Use the Linux container on macOS and Windows. The current Bioconda Gubbins
build can install through Rosetta on Apple silicon but its compiled
recombination step is not reliable there.

```bash
git clone https://github.com/ghruproject/beyondmlst.git
cd beyondmlst
docker build -t beyondmlst .
docker run --rm -v "$PWD:/data" -w /data \
  beyondmlst \
  run metadata.csv --output beyondmlst_results --threads 8
```

All paths in `metadata.csv` must resolve inside the mounted `/data` directory.
The CI-built `ghcr.io/ghruproject/beyondmlst:main` image is also available to
authenticated organisation members. The organisation currently disables
public package visibility, so anonymous pulls are not available.

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
- `patient_id`: coded identifier retained in copied metadata for downstream
  interpretation; it is not used to infer transmission

Use de-identified metadata only. Do not commit patient-level metadata or
sequence data to this repository.

## Commands

Validate metadata without running external tools:

```bash
pixi run beyondmlst validate metadata.csv
```

Preview the planned lineage analyses:

```bash
pixi run beyondmlst run metadata.csv --dry-run
```

Run the workflow:

```bash
pixi run beyondmlst run metadata.csv \
  --output beyondmlst_results \
  --threads 8 \
  --date-randomisations 100
```

The date-randomisation test compares the observed root-to-tip fit with fits
obtained after permuting collection dates. A dated tree is generated only when
the observed rate is positive and the randomisation p-value is at or below the
configured threshold (default `0.05`).

## Main outputs

Each lineage directory contains:

- `inputs.tsv` and `metadata.csv`: exact inputs used
- `core_alignment.fasta`: reference-ordered SKA2 alignment
- `gubbins.final_tree.tre`: recombination-filtered phylogeny
- `core_alignment.recombination_masked.fasta`: masked whole-genome alignment
- `clock/`: observed root-to-tip analysis
- `temporal_signal.json`: observed and randomised temporal-signal statistics
- `timetree/`: dated tree and uncertainty, when supported
- `location/`: exploratory ancestral location-state reconstruction
- `report.json`: machine-readable lineage report

`summary.json` records every lineage, skipped analysis, command and output.

## Scientific guardrails

- Analyse *Klebsiella* and *E. coli* separately and stratify them into
  sufficiently close lineages before alignment.
- Treat root-to-tip regression as exploratory; date randomisation is the
  workflow gate.
- Do not interpret a local-only tree as evidence against introductions.
- Balance contextual sampling across geography and time where possible.
- Review longitudinal isolates from the same patient explicitly; genomic
  proximity alone is not proof of direct transmission.
- Confirm important dates/rates with BactDating or BEAST before publication.

## Development

```bash
pixi run test
pixi run lint
```

See [docs/DESIGN.md](docs/DESIGN.md) for scope and planned extensions.
