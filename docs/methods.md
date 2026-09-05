# Methods

## Whole-genome phylogeny

ChronoClade analyses each `species` and `lineage` group independently. SKA2
builds a split-k-mer index and maps samples to a selected reference to produce a
reference-ordered alignment. Before tree inference, ChronoClade compares each
genome's median raw SNP proportion with the cohort distribution. A genome is
flagged only when it is both at least five times the cohort median and an
extreme robust outlier (robust z-score at least 10). The workflow stops and
writes `lineage_coherence.tsv` when this screen fails, avoiding an expensive
tree and recombination run on an obvious accession, species or lineage error.
After ClonalFrameML, the same relative check is repeated on the filtered
alignment and written to `clonal_lineage_coherence.tsv` before clock analysis.
This relative check is an input safeguard, not a universal bacterial SNP
threshold or a transmission definition.

IQ-TREE estimates the starting maximum-likelihood phylogeny with a GTR+G model.

ClonalFrameML estimates recombination on that tree and writes a corrected tree
plus an alignment of sites not assigned to imported regions. Pairwise clonal SNP
counts use A, C, G and T sites callable in both members of each pair. The report
always gives the callable-site count beside the SNP count. These corrected
distances provide a second opportunity to identify unusual genomes after
recombination has been modelled.

## Root-to-tip analysis

TreeTime reroots the recombination-corrected phylogeny by least squares and fits
root-to-tip distance against sampling date. ChronoClade disables automatic
clock-outlier removal at this gate so that samples are not discarded because
they weaken the temporal fit.

Root-to-tip regression is useful for visual inspection. It is not treated as a
formal temporal-signal test because population structure can associate genetic
distance with sampling date even when the dataset cannot estimate evolutionary
times reliably.[^murray]

## Date-randomisation test

ChronoClade permutes sampling dates among tips while retaining the corrected
tree and sequence length. TreeTime repeats the same least-squares clock analysis
for each permutation. The empirical one-sided p-value compares observed R² with
the permutation distribution and includes a one-count correction.

A dataset passes when:

1. every requested randomisation completes;
2. the observed rate is positive; and
3. the corrected empirical p-value is at or below the configured threshold.

The default is 100 permutations and p <= 0.05. This is a pragmatic gate for the
workflow, not a universal definition of temporal signal. Structured datasets
may require clustered permutations or sensitivity analyses that preserve known
population groups.[^duche]

## Time scaling

Passing datasets are analysed with TreeTime using marginal time inference,
covariation and 90% confidence intervals. The output tree has branch positions
on a calendar-time axis. ChronoClade exports the root estimate, clock-rate
uncertainty and each internal node's date interval.

The node intervals quantify uncertainty within the fitted TreeTime model. They
do not account for uncertainty caused by incomplete sampling, metadata error or
an inappropriate clock model, so those limitations must remain visible in the
interpretation.

## Location states and public context

TreeTime's mugration model reconstructs the supplied `location` state on the
dated tree, or on the corrected genetic tree when temporal signal is absent.
This reconstruction is exploratory. Its result depends on how locations were
defined and sampled.

The optional context workflow queries a compact AllTheBacteria metadata
snapshot for high-quality same-ST assemblies. It balances the candidate pool
across country and year, screens candidates against every focal sample with SKA
distance, then retains nearby genomes and a stratified background. The manifest
records the bounded search. "Nearest" means nearest within the downloaded pool,
not nearest among all public bacterial genomes.

## References

- Sagulenko P, Puller V, Neher RA. 2018. [TreeTime: Maximum-likelihood phylodynamic analysis](https://doi.org/10.1093/ve/vex042). *Virus Evolution* 4:vex042.
- Didelot X, Croucher NJ, Bentley SD, Harris SR, Wilson DJ. 2018. [Bayesian inference of ancestral dates on bacterial phylogenetic trees](https://doi.org/10.1093/nar/gky783). *Nucleic Acids Research* 46:e134.
- Didelot X, Wilson DJ. 2015. [ClonalFrameML: efficient inference of recombination in whole bacterial genomes](https://doi.org/10.1371/journal.pcbi.1004041). *PLoS Computational Biology* 11:e1004041.
- Murray GGR et al. 2016. [The effect of genetic structure on molecular dating and tests for temporal signal](https://doi.org/10.1111/2041-210X.12500). *Methods in Ecology and Evolution* 7:80-89.

[^murray]: Murray et al. showed that temporal and genetic structure can inflate root-to-tip regressions and date-randomisation tests.
[^duche]: A clustered permutation should be considered when collection date is confounded with a well-supported genetic group.
