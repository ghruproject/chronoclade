# First relatedness pilot: Zhou2017 ST11

**Method correction:** the positional locus-nucleotide results below are retained
as superseded evidence. The supported runtime now aligns observed alleles per
locus with MAFFT, reconstructs each genome's aligned loci and concatenates them.
See the final section for the corrected comparison.

18 October 2026. Development-study findings, not an external validation or a
recommended new default. Protocol and source labels were frozen before native
typing and distance evaluation. The independent numerical audit reproduced all
465 genome pairs and 403 query-rule rows in each coverage analysis.

## Data and workflow

The [source paper](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.00371/full)
reports 26 outbreak isolates and five sporadic controls, all ST11. Its event
assignment used epidemiology plus PFGE/MLST. It is independent of our future
cgMLST calls, but is not independently adjudicated transmission truth. Controls
are explicitly unrelated to the event; control–control pairs remain unknown.

All 31 versioned assemblies were resolved across multiple archive projects,
downloaded and checked against the frozen metadata/WGS identifiers. The ordinary
ChronoClade native prepare command called profiles and LIN assignments using one
existing reference snapshot. All records have 626–628 of the canonical 629 loci
called. Pairwise joint coverage is 624–628; all 465 pairs pass the primary
599-locus minimum. Missing calls were excluded from comparisons.

The source used 1,143 targets. These results retype its genomes under 629 loci;
they do not reproduce or transfer its numerical cutoff.

## Measured separation

| Published relationship | Genome pairs | Allele differences: minimum / median / maximum |
| --- | ---: | --- |
| Same documented event | 325 | 0 / 0 / 6 |
| Event versus explicit sporadic controls | 130 | 10 / 13 / 20 |
| Sporadic versus sporadic (primary label unknown) | 10 | 4 / 7.5 / 13 |

The panel has a separation between within-event and event–control distances.
This is one event with five controls, not 130 independent negative observations.
No confidence interval assumes pairwise independence.

## Predeclared neighbourhood rules

All methods search the same 31-genome panel, excluding the query itself. For each
outbreak query, 25 other event isolates and five controls are available.

| Rule | Event pairs retained / 325 | Event–control links admitted / 130 |
| --- | ---: | ---: |
| Same LIN5 | 325 | 130 |
| Same LIN6 | 325 | 104 |
| Same LIN7 | 325 | 0 |
| Same LIN8 | 276 | 0 |
| Same LIN9 | 231 | 0 |
| Same LIN10 | 231 | 0 |
| At most 0 allele differences | 210 | 0 |
| At most 1 | 232 | 0 |
| At most 2 | 273 | 0 |
| At most 4 | 320 | 0 |
| At most 7 | 325 | 0 |
| At most 10 | 325 | 2 |

LIN7 retains every event peer without admitting the five controls in this panel.
LIN8 retains 84.9% of event-pair links; deeper levels retain 71.1%. Broader LIN5
and LIN6 admit the controls. This shows the recovery/burden trade-off and why
finer is not automatically better. It does not make LIN7 or seven alleles a
universal outbreak rule. Earlier ST258 broad-scope concerns remain separate.

Nearest-tie reporting remains distinct from event recovery. LIN8 retains 463 of
505 directed nearest ties across the 26 outbreak queries; a prefix-only search
would miss some closest available comparisons even within this panel.

## Code identity and coverage

One pair shares the complete assigned code `0_0_105_2_0_0_2_8_0_0` but differs at
two jointly called loci. Both are native complete assignments from the same
snapshot; the observed mismatches include categorical novel-allele hashes.
The numerical audit preserves the differing loci, values and upstream matched
cgSTs. An assigned full code must not be described as proof of profile identity.
The mechanism of reference assignment is a separate upstream-method question;
this experiment does not change codes or calls to force agreement.

The 90% and 95% shared-call sensitivities retain all pairs and reproduce the
primary results. Requiring 629/629 joint calls retains no pairs: every genome
lacks a call at `KP1_0515_S`. That sensitivity is uninformative, not evidence of
zero genetic differences. The canonical locus was preserved as missing rather
than silently removed or replaced by an allele.

## Time and SNP comparison

The event's positive pairs occupy only the 0–30 and 31–90 day bins. Some negative
comparisons share those date spans. This collection cannot estimate a general
time-dependent cutoff or an alleles-per-year rule.

The standard strict SNP run retained all 31 genomes through alignment and
recombination correction. Its corrected-distance screen flagged exactly the five
published sporadic controls and stopped. The failed run is preserved; it is not
reported as successful default-policy acceptance. A separately frozen
[diagnostic amendment](AMENDMENT_01.md) retains all 31 for an explicitly flagged
SNP comparison, without temporal assessment or calendar dating.

That diagnostic run completed with all 31 genomes retained. Independent counting
from the recombination-filtered alignment reproduced all 465 SNP pairs and their
callable-site denominators (5,016,446–5,016,447).

| Published relationship | Corrected SNPs: minimum / median / maximum |
| --- | --- |
| Same documented event | 0 / 6 / 16 |
| Event versus explicit sporadic controls | 71 / 88 / 123 |
| Sporadic versus sporadic (primary label unknown) | 36 / 65.5 / 93 |

The corrected tree has an edge separating the 26 event isolates from the five
controls. No bootstrap support was estimated. Among 210 pairs with zero allele
differences, corrected SNP distances range from 0 to 14 (median 3). The filtered
core alignment and 629-gene typing scheme measure different feature sets.

Across the 26 outbreak queries, LIN7 retains all 51 directed corrected-SNP
nearest-neighbour ties; LIN8 retains 48/51 and LIN9/10 retain 45/51. Selecting
only the closest allele-distance ties retains 48/51. These are neighbourhood
retrieval checks, not direct-transmission predictions.

The diagnostic run stayed within its declared laptop limits: peak monitored
process-family memory was approximately 5.94 GiB and output approximately
581 MiB. Its completion does not turn the preserved strict-policy failure into
a successful strict run.

## Local evidence

- Protocol: `docs/research/klebsiella-grouping/PROTOCOL.md` and its hash freeze.
- Labels, versioned accessions, source/sequence checksums, independent input and
  numerical audits: `validation/science/zhou2017/`.
- Prepared dataset and ordinary cgMLST report: `tmp/science-zhou2017/`.
- Candidate tables and figure report: `tmp/science-zhou2017/evaluation-primary/`.
- Original strict SNP failure: `tmp/science-zhou2017/tree/`.
- Separate diagnostic run: `tmp/science-zhou2017/tree-diagnostic/`.

Large inputs remain local. No runtime defaults were changed. The next scientific
step is an independent event/lineage with reliable accession and label mapping,
then explicit assessment of sampling-density and time-span effects.

## Locus-nucleotide and threshold extension

Amendment 02 records that the earlier allele/core results were already known.
The new comparison uses 667 observed allele sequences extracted from verified
native calls, with 19,458 called sequences checked against their hashes and
known numeric alleles checked against the frozen reference. Missing calls remain
missing; no observed called sequence was unresolved.

All 465 pairs pass the 599-locus sequence-overlap minimum. Compared sequence
coverage varies from 618–628 loci and 504,706–511,735 sites. Nine loci have
different observed allele lengths: the current positional comparator excludes
397 pair-locus comparisons across 186 genome pairs. Equal-length sequences were
not independently aligned, so compensating indels remain an additional limitation.

| Published relationship | Locus nucleotide differences: minimum / median / maximum |
| --- | --- |
| Same documented event | 0 / 0 / 5 |
| Event versus explicit sporadic controls | 7 / 10 / 19 |
| Sporadic versus sporadic (primary label unknown) | 1 / 6.5 / 13 |

Every retained differing allele comparison has exactly one positional nucleotide
difference. This panel therefore provides no extra substitution resolution from
the locus-nucleotide route. All 210 zero-allele pairs also have zero differences
at the compared gene positions; 205 of these have nonzero corrected core SNPs.
Conversely, 21 pairs with two allele differences become zero nucleotide
differences because both unequal-length allele pairs are excluded. Those zero
values must not be described as identical cgMLST gene panels. Proper alignment
and length-variation handling are required before promoting this route.

The ordinary `chronoclade cgmlst --distance-method direct-dna` command completed
on all 31 genomes with the checked catalogue. Its explicit fraction 0.00001
exercises the documented workflow, not a fitted epidemiological boundary.

The shared incremental single-linkage engine sweeps every integer raw-distance
cutoff. Independent graph reconstruction reproduced all 165 requested partitions
and all 73 distinct membership summaries across the three distance definitions.
The 26 event isolates remain one separate membership set over these observed
cutoffs:

| Distance | Single-linkage cutoff range | Maximum internal distance | Closest outside distance |
| --- | --- | ---: | ---: |
| Allele mismatches | 3–9 | 6 | 10 |
| Positional locus nucleotide differences | 3–6 | 5 | 7 |
| Corrected core SNPs | 7–70 | 16 | 71 |

These are descriptive results for this panel, not recommended cutoff ranges.
Single-linkage joining cutoffs differ from all-pairs diameters. The comparisons
use the same genome pairs, but different genomic denominators and exclusions.
The sampled event envelope is 8 October–23 December 2011. Date intervals and
reported location metadata accompany every membership set without affecting the
genetic clustering; no rate or time-dependent threshold is fitted.

Updated local comparison: `tmp/science-zhou2017/threshold-comparison-final/report.html`.
Normal nucleotide report: `tmp/science-zhou2017/cgmlst-locus-dna-final/index.html`.
Audit: `validation/science/zhou2017/threshold-audit.json`; locus sequence and
independent concatenation evidence: `tmp/science-zhou2017/locus-dna/`.
The shared engine is ready for integration; this new threshold view is a research
report, not yet part of the normal product reports. See `NEXT_PHASE.md` for the
remaining scientific and integration gates.

## Correction: align each locus before concatenating

The user identified the inadequate positional shortcut. The normal direct-DNA
workflow now uses per-locus MAFFT alignments, including variable equal-length
alleles, reconstructs sample-specific sequences and concatenates in canonical
locus order. Missing calls become Ns; alignment gaps remain gaps and are reported
separately from substitutions. The globally absent locus remains a zero-width
partition with explicit absent status, not an invented column.

The alignment contains 511,735 columns, with 39 variable loci aligned by MAFFT,
589 invariant loci reconstructed exactly and one absent locus. All 465 pairs
remain eligible: 624–628 compared loci and 507,211–511,735 jointly A/C/G/T sites.
The 30 equal-length variable locus comparisons still each contain one substitution.

| Published relationship | Aligned locus SNPs: minimum / median / maximum |
| --- | --- |
| Same documented event | 0 / 0 / 8 |
| Event versus explicit sporadic controls | 7 / 10 / 21 |
| Sporadic versus sporadic (primary label unknown) | 1 / 6.5 / 13 |

Alignment restores four substitutions in `glnL_S` plus 156 gap-versus-base columns
in one indel tract. Those four substitutions were lost by the positional method
in all 30 pairs involving GCA_001613665.1. The other eight unequal-length variable
loci have indel evidence without additional substitutions in this panel.

The 21 pairs with two allele differences and zero locus SNPs now retain those
loci and explicitly record their indels. Zero substitutions still does not mean
sequence identity, but it is no longer caused by dropping the entire locus.
Per-locus scalar arithmetic and independent vectorised concatenation arithmetic
agree for all 465 pairs. Alignment uncertainty remains; the gene panel has not
been recombination-filtered.

Corrected report: `tmp/science-zhou2017/threshold-aligned/report.html`.
The normal direct-DNA report and selection bundle are generated by the standard
CLI under `tmp/science-zhou2017/cgmlst-locus-aligned-final/`. The bundle retains
all 31 mandatory inputs. No biological cutoff was selected or changed.
