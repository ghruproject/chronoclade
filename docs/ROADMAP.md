# Roadmap to an actionable public-health report

## Working foundation

The current MVP now provides the reproducible analytical foundation:

- one Pixi environment for the Python CLI and compiled tools;
- species/ST-separated SKA2, IQ-TREE and ClonalFrameML analysis;
- root-to-tip graphics, date randomisation and a TreeTime dated tree only when
  the temporal-signal gate passes;
- same-ST discovery for all represented bacterial species and schemes from a
  bundled compact AllTheBacteria snapshot, with assembly fetching through
  `atbfetcher`;
- deterministic country/year-balanced context pooling, SKA screening and a
  frozen context manifest;
- recombination-filtered pairwise SNP and callable-site matrices, longitudinal
  distance summaries and an SVG heatmap;
- topology-defined candidate local groups, final clonal contextual neighbours
  and one-isolate-per-patient distance sensitivity;
- decision-first HTML and JSON reports with a transparent four-scenario evidence
  ledger, recommended follow-up and explicit interpretation limits; and
- Linux/macOS CI that executes the native time-tree and location-tree paths.

This is sufficient to run and audit a provisional public-health scenario. Its
confidence is deliberately capped until branch support, exhaustive neighbour
retrieval and GHRU validation are complete.

## Next release: decision evidence

### 1. Clonal distance calibration

The workflow now calculates pairwise SNP differences from the
recombination-filtered alignment, with callable sites for every comparison. The
next step is to calibrate these outputs on real datasets and add:

- within-candidate-group and between-candidate-group distributions; and
- an optional, externally justified analysis-specific threshold.

No universal SNP threshold should be built in. An externally justified,
analysis-specific threshold may be supplied and must be shown in the report.

### 2. Supported topology and cluster summaries

Add branch support to the starting phylogeny and propagate it onto the corrected
topology-defined groups. The existing groups and scenario must remain
provisional until this support is available. Add a metadata-coloured tree that
explicitly distinguishes monophyletic local structure from local groups
interspersed among contextual genomes.

### 3. Better public neighbour retrieval

The present bounded, stratified context pool cannot guarantee that it contains
the nearest genome among all public records. Add a whole-database sketch query
for each focal isolate or preliminary focal cluster, then combine those hits
with the geography/time-stratified background before the SKA screen. Retain the
query database version, score, rank and attrition in the manifest.

### 4. Longitudinal epidemiology

Extend the optional metadata contract with coded facility, ward or sampling
site fields. Produce a patient/location timeline and rerun topology inference
using one representative isolate per patient. Compare the full and deduplicated
topologies; the current report already compares their corrected-distance
summaries.

### 5. Scenario synthesis

Calibrate the implemented transparent evidence table for four summaries:

1. consistent with a persistent local lineage;
2. consistent with multiple introductions;
3. local persistence plus additional introductions; or
4. indeterminate.

The summary should cite the exact topology, distance, context and temporal
results that support or contradict it. It must remain indeterminate when
context coverage, branch support, callable sites or temporal resolution are
insufficient. It must not claim direct transmission or count introductions as
observed events.

### 6. Validation and calibration

Before a stable release, run blinded review on at least:

- one dominant *Klebsiella pneumoniae* ST and one *Escherichia coli* ST from
  the GHRU longitudinal data;
- a dataset expected to show long-term local persistence;
- a dataset expected to contain several divergent local groups;
- a mixed or deliberately under-contextualised dataset; and
- one-isolate-per-patient and alternative-context sensitivity runs.

For each dataset, record whether two independent public-health reviewers reach
the same interpretation from the report and which additional evidence they
request.

## Decisions to make with the first GHRU datasets

- Define “local” for each analysis: facility, city, region or country.
- Confirm the species/ST combinations and whether any ST must be split into a
  narrower genomic lineage before analysis.
- Decide which contextual locations and collection periods are essential,
  while retaining a global background.
- Confirm the minimum metadata that can be used safely: month precision,
  coded patient, facility and ward/site.
- Choose the number of context genomes from runtime and sensitivity evidence,
  not an arbitrary fixed default.
- Agree who signs off the four-scenario interpretation and how an
  “indeterminate” result is communicated.

## Workflow-engine decision

Keep the Python runner for the next release. It already bounds parallel
lineages and date permutations within one CPU budget. Reconsider Nextflow only
if real use requires multi-node scheduling, a read-to-report pipeline, or
large cohort retries that cannot be handled reliably by the current
stage-resume model.
