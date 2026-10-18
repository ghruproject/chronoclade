# Relatedness and grouping method review

Date: 17 October 2026. Scope: bounded primary-literature and upstream-method
review for `RES-01`. This review recommends a research deferral for automatic
outbreak thresholds. It does not change the approved beta scope, validate a
threshold, or claim that a published method is validated for ChronoClade.

## Decision

Keep reported cgLIN levels and raw allele-distance neighbours as separate
descriptive evidence. The local level-8 pilot is a plausible way to show a
smaller ST258 neighbourhood beside the wider background, but it is not an
outbreak cluster or a transmission inference. Do not choose an automatic
distance cutoff or a time-varying cutoff from the current ST258 example.

For `RES-01`, defer outcome claims and threshold selection until an independent,
scheme-compatible set of epidemiologically linked and unlinked isolates is
frozen with adequate provenance. The observed root-to-tip slope is not a
pairwise relatedness threshold. A time-aware outbreak model is a candidate only
for a study with the outbreak mechanism, dates and epidemiological information
its assumptions require.

## Evidence from published and upstream methods

| Approach | What it estimates or assigns | What its evidence supports | Limits for ChronoClade |
| --- | --- | --- | --- |
| Stable, nested cgLIN nomenclature | The *K. pneumoniae* species-complex system uses 10 ordered bins on the 629-locus cgMLST scheme. Its first four levels capture population structure and clonal groups; downstream bins encode increasingly close cgMLST similarity. A new cgST inherits a prefix from its closest already encoded profile. | Hennart et al. built the taxonomy from more than 7,000 genomes and defined levels using phylogenetic/clustering discontinuities and explicit similarity bins. The Pasteur upstream specification documents how the identifier is assigned, missing-locus handling, provisional assignments and reference dependence. This is a strong baseline nomenclature for describing KpSC diversity. [Hennart et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9254007/); [Pasteur cgLIN specification](https://bigsdb.pasteur.fr/klebsiella/cgmlst-lincodes/) | A prefix is a stable nomenclature value assigned through a nearest-previously-encoded profile procedure, not a newly fitted complete-linkage cluster of the present cohort. It must not be presented as a guarantee that every pair in the prefix is within the bin's nominal mismatch bound. The upstream scheme is exactly 629 loci and permits at most 30 missing loci; incomplete/provisional codes and source reference snapshots matter. ChronoClade's frozen ST258 profile export declares 628 observed loci, `complete=False`, and no database version/fingerprint, despite imported codes appearing complete. The local audit therefore cannot validate absolute prefix thresholds against the source scheme. |
| Fixed pairwise allele cutoff with epidemiological controls | A threshold groups samples when their pairwise allele distance is within a selected bound. | A Kp study evaluated a 1,143-target cgMLST scheme on a hospital ST11 outbreak: 26 outbreak isolates had a maximum distance of 10 alleles, while unrelated ST11 genomes provided comparison controls. The same paper found pairs of unrelated ST11 genomes at fewer than 10 differences, and explicitly cautioned that the ≤10 boundary needed more study. [Zhou et al. 2017](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.00371/full) | The value 10 is not portable to the Pasteur 629-locus scheme or to ST258: it came from a different 1,143-locus scheme, one short hospital outbreak and a limited comparison panel. Close allele distance can include epidemiologically unrelated isolates; a fixed number alone does not define an outbreak. This paper is a useful evaluation design, not a recommended ChronoClade cutoff. |
| Time-window subclustering plus epidemiological review | Within a broad cluster, a stricter distance and date window can create reviewable subclusters; subsequent typing and epidemiological investigation refine them. | In a New York *Salmonella Enteritidis* investigation, a cluster of 265 isolates over 2.5 years had grown by single-link chaining from a five-allele/60-day rule to a 14-allele span. A retrospective zero-allele/60-day cgMLST refinement made 34 groups; wgMLST, hqSNP and review of exposures, geography and dates further refined most groups. The authors say epidemiological links were essential and laboratory results alone were insufficient. [Baker et al. 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10282972/) | The study reports roughly 3,000 cgMLST loci and 25,000 wgMLST loci. It demonstrates an investigation workflow and the danger of single-link chaining, not a portable zero-allele cutoff: the work concerned a highly clonal serovar, a different scheme/database and epidemiologists reviewing candidate links. |
| Outbreak-specific model-based threshold (`sameStrain`) | A forward model simulates expected distances for a single-strain, point-source food/environmental outbreak from mutation rate, contamination duration and observed sample dates; it uses the 99th percentile of simulated pairwise distances as an outbreak-specific threshold. | Duval et al. evaluated the framework on simulated outbreaks and 16 published foodborne outbreaks. The classification agreed with original outbreak/sporadic labels for 14 of 16 datasets; two *Listeria monocytogenes* datasets were discordant. Their simulations found specificity can be poor with low mutation rates, and fewer than 10 isolates can bias duration and rate estimates downwards. [Duval et al. 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10156608/) | It targets a particular common-source mechanism and is not a general population-neighbourhood or human-to-human transmission rule. It needs plausible outbreak duration or mutation-rate information, sample dates, a suitable distance model and enough outbreak isolates. It is not feasible as a routine default for the current global ST258 context search without a defined source hypothesis and independent Kp data. |

The cgLIN source states that bins 5–10 are based on 10 down to 0 mismatches in
the 629-locus scheme. Similarity calculations account for missing alleles, and
profiles with more than 30 missing loci are not assigned a cgLIN. Importantly,
the code is assigned to cgST profiles and each new profile is encoded relative
to the closest previously encoded profile. The local ST258 audit shows why the
actual pairwise distance matrix still matters: observed shared-prefix members
can have larger endpoint distances than the bin's nominal neighbour mismatch
bound, and 37 nearest input–public pairs cross an input's own level-8 prefix.
See [ST258_GROUPING_AUDIT.md](ST258_GROUPING_AUDIT.md) for the frozen counts,
coverage and source-provenance findings.

### A neighbourhood is not an outbreak or transmission claim

A neighbourhood view answers “which records are genetically close to these
inputs in this searched dataset, under this scheme and coverage rule?” It is a
query and display scope. An outbreak threshold makes a stronger classification
claim about membership in a common epidemiological event. A transmission claim
is stronger again: it concerns direct or indirect spread and requires suitable
epidemiological and temporal evidence. Allele similarity, a shared prefix, a
tree clade or a root-to-tip association alone cannot establish those claims.

A root-to-tip regression asks whether genetic divergence from a selected root
is associated with sampling time across tips. It does not estimate the
distribution of pairwise distances for linked versus unlinked samples and does
not specify an outbreak process. The [TreeTime documentation](https://treetime.readthedocs.io/en/latest/tutorials/clock.html)
describes its clock regression as a way to assess clock-like behavior and
estimate substitution rate from root-to-tip distances. The conclusion that
this cannot, by itself, classify an outbreak pair is a methodological inference
from the different estimands. `sameStrain` is time-aware for a different
reason: it simulates a stated single-source process using mutation rate,
outbreak duration and sampling dates. A dynamic threshold therefore requires a
specified biological model and independent labels or external validation; a
slope-derived multiplier or extrapolation is insufficient.

## Preregistered candidate evaluation

This plan is a prerequisite for evaluating thresholds; it does not select one.

1. **Freeze an independent labelled benchmark before fitting.** Assemble
   scheme-compatible *K. pneumoniae* datasets that include high-confidence
   epidemiological links (for example, documented shared source or a
   well-defined hospital outbreak) and credible negative comparisons. Define
   label rules with an epidemiologist who does not see ChronoClade distance
   results. “No link recorded” is unknown, not a negative label. Include more
   than one outbreak, site and time period. Split by whole study/site or event
   so related isolates and their pairwise combinations cannot leak between
   training and held-out evaluation.
2. **Lock identity, scheme and callability.** Record accession and assembly
   identity/checksum, cgMLST locus catalogue and fingerprint/version, caller and
   settings, cgLIN assignment source/version, dates and label provenance.
   Require the exact same compatible locus universe for a pair. Preserve raw
   mismatch count, jointly called loci and missingness for every distance;
   preregister exclusions/coverage strata. Do not translate allele cutoffs
   between different scheme sizes or reference versions.
3. **Name the task before comparing methods.** For an exploratory neighbourhood,
   compare the reported cgLIN prefixes at levels 5–10 with all actual nearest
   ties ranked by raw mismatches, showing the jointly called-locus denominator
   and global-nearest exceptions. This evaluates useful display scope and
   retention, not outbreak sensitivity. For an outbreak-classification task,
   compare the fixed cgLIN and distance baselines on the locked 629-locus scheme.
   If that scheme and its reference fingerprint are verified, preregister a
   pairwise cutoff sensitivity grid of 0, 1, 2, 4, 7 and 10 mismatches. These
   values are the fine-bin boundaries documented for KpSC cgLIN, used here only
   as candidate pairwise cutoffs to test independently; they do not inherit
   cgLIN prefix semantics. Exclude thresholds imported from other schemes as
   direct candidates unless their underlying profiles can be recomputed under
   this scheme. Test a time-aware model only on events meeting its single-source
   assumptions and only if its inputs and implementation can be independently
   reproduced.
4. **Evaluate once on held-out events.** Report pair-level and event-level
   sensitivity, specificity, precision/recall with uncertainty, false merges,
   missed linked isolates, and whether chaining merges separate labelled
   events. Include performance by elapsed-time span, lineage/ST, profile
   completeness, shared-locus coverage, and sampling density. Repeat after
   controlled downsampling and missing-call perturbations. Report nearest-tie
   retention and distances after any proposed scope reduction. Do not choose a
   cutoff on the held-out data or use that same set as proof of its performance.
5. **Set acceptance criteria before revealing results.** An epidemiology/public
   health reviewer must specify which false merges and missed links are
   tolerable for the intended alert/review use, and how those criteria trade
   off. A candidate advances only if it improves on or clearly complements the
   reported-lineage and nearest-distance baselines on more than one independent
   held-out event set, with no unexplained scheme/coverage failures. Otherwise
   keep the existing grouping, show descriptive neighbours and document the
   negative result. Never label a method a “pass” because it matches the
   threshold used to construct its own labels.

## Feasibility and beta disposition

The existing frozen ST258 profiles are sufficient for descriptive prefix and
nearest-neighbour comparisons without constructing a second dense matrix; the
local audit already reproduces the stored nearest pairs and measures
cross-prefix exceptions. They are not sufficient to estimate classifier
performance because the source snapshot is incomplete and no independently
adjudicated linked/unlinked labels are attached. Recovering the exact upstream
reference snapshot and obtaining independent epidemiological labels are the
concrete prerequisites.

No additional outbreak classifier is justified for the core beta from this
review. Keep the level-8 report pilot as a candidate descriptive neighbourhood
view, retain the wider background and all true nearest ties, and show actual
pair distances and coverage. Defer fixed or dynamic outbreak cutoffs to the
preregistered evaluation above. This is a research deferral for `RES-01`; it
does not block the core workflow and does not turn the existing ST258 pilot into
a validated production grouping policy.
