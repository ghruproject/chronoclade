# Reading the report

Each lineage report follows the order in which an analyst should inspect the
evidence. The temporal test and the public-health interpretation are separate.
A lineage can fail the clock test and still have useful clonal distances and
topology.

## 1. Check the analysis basis

The report begins with sample count, sampling span and the number of complete
clonal sites. It records that the raw-distance screen, ClonalFrameML correction
and post-recombination distance screen completed before temporal testing. These
checks detect grossly divergent genomes; they do not prove every accession or
metadata value is correct.

## 2. Explore the clock relationship

The root-to-tip plot places collection date on the horizontal axis and genetic
distance from the fitted root on the vertical axis. The slope estimates a rate;
R² describes the fit of the regression.

Use this plot to inspect the direction of the slope, scatter, outliers and
obvious clusters. Do not decide that temporal signal is present from R² alone.

## 3. Test against random dates

ChronoClade permutes collection dates among the tips and repeats the TreeTime
clock fit. The report plots the null distributions for R² and clock rate beside
the observed values.

The empirical p-value is:

```text
(permuted R² values at least as large as observed R² + 1)
---------------------------------------------------------
            (successful permutations + 1)
```

The default gate requires all requested permutations to complete, a positive
observed rate, and p <= 0.05. ChronoClade stops time scaling when any condition
fails.

## 4. Inspect the dated tree

Passing datasets receive a TreeTime phylogeny on a calendar axis. Horizontal
intervals show the uncertainty in internal-node dates. The report also gives
the clock rate and its standard deviation, root estimate, and summaries of
node-interval width.

## 5. Review topology, distance and context

The working interpretation uses the recombination-filtered topology, clonal SNP
distances, longitudinal span and placement of public context genomes. The
evidence table states which observations support the interpretation and which
inputs are missing.

The report uses four labels:

| Label | Sampled pattern |
| --- | --- |
| Persistent local lineage | One focal-only group spans more than one sampling date |
| Multiple introductions | Focal isolates form at least two separated groups and none spans dates |
| Mixed | Several focal groups are present and at least one spans dates |
| Indeterminate | Focal, longitudinal or contextual evidence is insufficient |

These are working genomic descriptions. An analyst should compare them with
patient movement, ward, referral, travel and sampling information.

## Downloadable evidence

Every chart has a vector SVG and, where useful, a high-resolution PNG. CSV and
TSV files hold the numerical values. Newick or Nexus files preserve the trees;
JSON files preserve the report inputs and decisions. `supporting_results.zip`
packages the reader-facing evidence without the large intermediate alignments.

## Fast-screen report

Fast mode has three deliberately limited stages: PhiPack recombination
detection, root-to-tip inspection on an uncorrected tree, and an unclustered
date-permutation screen. It never presents a dated tree or a
circulation/introduction interpretation. A PHI-positive block at the unadjusted
screening threshold is an escalation signal for the full ClonalFrameML analysis,
not a recombinant tract.
