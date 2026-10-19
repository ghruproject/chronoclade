# Methods

## Stages and profile analysis

ChronoClade separates preparation, profile analysis and selection, corrected
SNP-tree inference, and dating into independent modules. Preparation resolves
inputs and typing; cgMLST analyses compatible profiles and selects representatives;
tree acquires the selected assemblies, infers a recombination-adjusted phylogeny
and assesses temporal signal; time dates a saved tree only when its evidence
supports dating. Exploratory profile-date regressions and formal tree-stage
temporal assessment answer different questions; the former never passes the
temporal gate.

Each profile comparison requires compatible species, cgMLST scheme/version,
database provenance where available, and a compatible locus set. Allele
identifiers are categorical; numeric allele IDs are not evolutionary
distances. Pairwise distances are the number of mismatches divided by the
number of jointly called loci. Missing calls are excluded from the numerator
and the report gives the shared-call count and overlap. Comparisons with
incompatible schemes or locus sets are not combined into one cohort.

Some exports list only observed loci and do not declare the full canonical
scheme universe. For such compatible records, ChronoClade uses the union of
observed loci as the overlap denominator and flags that canonical completeness
is unknown. This does not estimate the proportion of the full scheme covered.
Complete-comparability cohorts are deterministic partitions for ordination and
neighbour-joining views; they are not themselves genetic groups.

Genetic groups use deterministic complete linkage at the configured
`--group-distance`, expressed as a mismatch fraction. Locus-bootstrap
co-assignment is descriptive support for those candidate groups. If no valid
bootstrap replicates are available, support is unknown, not evidence against a
group. A singleton has no within-group pair support to estimate.

Genetic group stability under locus resampling and concentration in time/place
are reported separately.

The machine-readable collection-year summary lists the observed years in the group; gaps do not
demonstrate uninterrupted presence. Concentration summarizes year/location
cells among records with both date and place metadata; it does not test for an
outbreak or establish transmission. Both summaries describe the sample
submitted, not population prevalence. The existing `temporal_persistence` JSON
key retains this descriptive collection-year meaning for compatibility.

The profile-stage nearest relative is selected from the available compatible
public profiles and reports shared called loci. Country figures include
profile-available members of a comparable cohort, combine query and public
context records, and exclude records without profiles. They do not split
country counts by origin or provide a separate region-count view in the current
profile summary. The separate query/context coverage counts include records
without profiles.

## Whole-genome phylogeny

The following assembly methods apply to the `tree` module. Clock methods apply
to temporal assessment in `tree` and dating in `time`.

ChronoClade analyses each saved compatible cgMLST block independently. Optional
lineage partitioning can narrow those blocks. SKA2
builds a split-k-mer index and maps samples to a selected reference to produce a
reference-ordered alignment. Before tree inference, ChronoClade compares each
genome's median raw SNP proportion with the cohort distribution. A genome is
flagged only when it is both at least five times the cohort median and an
extreme robust outlier (robust z-score at least 10). Raw flags are written to
`lineage_coherence.tsv` and retained for review; recombination can inflate raw
distances. After the selected recombination filter, the same relative check is repeated on the filtered
alignment and written to `clonal_lineage_coherence.tsv` before clock analysis.
The default strict policy stops on filtered-alignment flags. Explicit diagnostic
policy retains all selected genomes and carries a visible warning into the reports.
Neither policy prunes samples. A flag can reflect legitimate cohort diversity;
it does not establish incorrect identity. This relative check is an input
safeguard, not a universal bacterial SNP threshold or a transmission definition.

IQ-TREE estimates the starting maximum-likelihood phylogeny with a GTR+G model.

The default recombination method is PHIPack Profile window filtering. It scans
10,000 mapped bases at a time in 5,000-base steps with a 100-base PHI neighbourhood,
without crossing reference-record boundaries. Benjamini–Hochberg adjustment is
applied across all finite window tests; windows with q ≤ 0.05 are removed in full
from every genome, along with columns not callable as A/C/G/T in every genome.
IQ-TREE then estimates the filtered tree anew. Untestable windows are retained
and explicitly counted; they are not evidence of no recombination. Window masks
are exploratory regions, not branch-specific importation calls. See the
[PHIPack manual](https://gensoft.pasteur.fr/docs/PhiPack/1.0/phimanual.pdf) for the
profile test and its unadjusted native output. ChronoClade checks native test
status as well as numeric CSV output because an unavailable approximation can
leave an uninitialised CSV value. The method and settings are recorded in tree
and time results. PHIPack results do not imply equivalent filtering or ages to
the earlier ClonalFrameML benchmark.

The optional `--recombination-method verticall` runs native Verticall
reference-to-query comparisons and sample-specific horizontal/unaligned masking
on each reference sequence independently. Safe called-segment assembly
snapshots meet its ambiguity requirements without changing original assemblies.
An exact auxiliary reference row supplies coordinates and is excluded after
masking; every selected sample is retained. Native masked alignments are saved,
and the downstream common alignment uses only columns called in every selected
sample before fresh IQ-TREE inference. This additional complete-site policy is
explicit, rather than implying that the native masked alignment has no missing
data. Secondary results use the first native result, not a result selected for
preferred clock behavior. Region calls and complete-column accounting are
retained. See the [Verticall alignment workflow](https://github.com/rrwick/Verticall/wiki/Alignment-tree-workflow).

With `--recombination-method clonalframeml`, ClonalFrameML estimates recombination on that tree and reports 1-based, closed
alignment intervals for each affected branch. With `-output_filtered true` and
`-ignore_incomplete_sites true`, it constructs one shared clonal alignment by
removing every column assigned to an import on any branch and every column with
an ambiguous base in any sequence. ChronoClade verifies that these two masks
exactly account for the filtered-alignment length. It exports the original calls,
a normalised interval table, a binned genome profile and SVG/PNG figures.
Reference-record coordinates are recovered from the ordered reference FASTA.
Intervals that cross a record boundary are retained for audit but flagged because
the adjacency was introduced by concatenation rather than by the chromosome.

A branch-level interval is evidence that ClonalFrameML assigned an import to
that lineage of the tree; it does not mean that every sampled genome acquired
the segment. Pairwise clonal SNP counts use A, C, G and T sites callable in both
members of each pair. The report always gives the callable-site count beside the
SNP count. These corrected distances provide a second opportunity to identify
unusual genomes after recombination has been modelled.

## Tree-stage root-to-tip analysis

TreeTime reroots the recombination-corrected phylogeny by least squares and fits
root-to-tip distance against sampling date. ChronoClade disables automatic
clock-outlier removal at this gate so that samples are not discarded because
they weaken the temporal fit.

Original year/month precision stays in the metadata; clock inputs serialize
years as `YYYY-XX-XX` and months as `YYYY-MM-XX` so TreeTime receives the
corresponding intervals, including in date permutations. This conversion does
not add support for censored dates or other interval formats.

Native clock eligibility also accepts explicit finite numeric years and
TreeTime's `[start:end]` numeric ranges, preserving their bounds. Informal ranges
such as `2006-2009` and censored dates such as `before 2009` are not converted
automatically: those tips remain in the selected tree without a calibration
date until a justified date model is supplied. No missing bound is invented.

Root-to-tip regression is useful for visual inspection. It is not treated as a
formal temporal-signal test because population structure can associate genetic
distance with sampling date even when the dataset cannot estimate evolutionary
times reliably.[^murray]

## Date-randomisation test

The tree stage's default `root-to-tip` method permutes sampling dates among
maximal monophyletic clusters sharing the same original date string while
retaining the corrected tree and sequence length. TreeTime repeats the same
least-squares clock analysis for each permutation. The empirical one-sided
p-value compares observed R² with the permutation distribution and includes a
one-count correction. R² drives this screen; the accompanying rate distribution
is descriptive.

Internal compatibility helpers also implement full TreeTime refits for each
permutation and weak CR1/strict CR2 interval criteria. They are not selectable
through the normal `tree` CLI. Those helpers use approximate 95% rate intervals
from TreeTime's reported standard error, not Bayesian credible intervals; they
do not replace the normal clustered root-to-tip gate described here.

A dataset passes when:

1. every requested randomisation completes;
2. the observed rate is positive; and
3. the corrected empirical p-value is at or below the configured threshold.

The default is 100 permutations and p <= 0.05. This is a pragmatic gate for the
workflow, not a universal definition of temporal signal. Structured datasets
may require additional sensitivity analyses that preserve known population
groups beyond the same-date clusters used here.[^duche]

## Stage boundaries

The cgMLST module analyses compatible profiles and saves representative
selections; it does not build a whole-genome tree. The tree module performs the
corrected genomic analysis and assesses temporal signal. When fewer than three
usable distinct collection dates are available, temporal analysis is recorded as
not assessed; it is not represented as a failed permutation result. The time
module reads the saved tree and assessment rather than rebuilding the phylogeny.

## Time scaling

When temporal evidence supports dating, the time module uses TreeTime with
marginal time inference, covariation and 90% model intervals. The output tree has branch positions
on a calendar-time axis. ChronoClade exports the root estimate, clock-rate
uncertainty and each internal node's date interval.

The node intervals quantify uncertainty within the fitted TreeTime model. They
do not account for uncertainty caused by incomplete sampling, metadata error or
an inappropriate clock model, so those limitations must remain visible in the
interpretation.

## Location states and public context

cgMLST reports can reconstruct country, host, region or city changes on the
saved tree. Missing metadata stays unknown, and a field with only one observed
state cannot support a useful network. Reconstructed changes depend on the tree,
root, metadata and sampling frame; they do not prove transmission or migration.
Geographic summaries describe available records, not prevalence, incidence or
country of acquisition. See [cgLIN](cglin.md).

The current public context workflow retrieves compatible profiles and metadata
from the known species/ST discovery pools. The default analysis partitions by
species and compatible profile scope; optional lineage filtering must be requested
explicitly. Missing ST does not block local analysis but limits public discovery.
cgMLST analyses those profiles
and selects representatives; tree acquires assemblies only for the selected IDs.
A nearest relative is nearest within the retrieved, compatible pool, not among
all public bacterial genomes. Provider coverage and incomplete annotations limit
interpretation.

Native typing uses Pathogenwatch's cgMLST caller and `plincer` or exact `hclink`.
Missing, ambiguous and novel alleles, incomplete assignments and database
provenance remain explicit. HClink is a reference-linking inference rather than
an official new EnteroBase cluster designation. See [query typing](query-typing.md).

## References

- Sagulenko P, Puller V, Neher RA. 2018. [TreeTime: Maximum-likelihood phylodynamic analysis](https://doi.org/10.1093/ve/vex042). *Virus Evolution* 4:vex042.
- [StrainHub project](https://github.com/abschneider/StrainHub). ChronoClade's country view is stylistically inspired by this tool but uses TreeTime's location reconstruction.
- Didelot X, Croucher NJ, Bentley SD, Harris SR, Wilson DJ. 2018. [Bayesian inference of ancestral dates on bacterial phylogenetic trees](https://doi.org/10.1093/nar/gky783). *Nucleic Acids Research* 46:e134.
- Didelot X, Wilson DJ. 2015. [ClonalFrameML: efficient inference of recombination in whole bacterial genomes](https://doi.org/10.1371/journal.pcbi.1004041). *PLoS Computational Biology* 11:e1004041.
- Bruen TC, Philippe H, Bryant D. 2006. [A simple and robust statistical test for detecting the presence of recombination](https://doi.org/10.1534/genetics.105.048975). *Genetics* 172:2665-2681.
- Murray GGR et al. 2016. [The effect of genetic structure on molecular dating and tests for temporal signal](https://doi.org/10.1111/2041-210X.12466). *Methods in Ecology and Evolution* 7:80-89.

[^murray]: Murray et al. showed that temporal and genetic structure can inflate root-to-tip regressions and date-randomisation tests.
[^duche]: A clustered permutation should be considered when collection date is confounded with a well-supported genetic group.
