# Klebsiella grouping benchmark: independent design review

Date: 18 October 2026. Status: proposed methods, to freeze before examining new
benchmark results. This document reviews the design independently of new
performance results. Its inputs are `docs/beta/GROUPING_METHOD_REVIEW.md` and
`docs/beta/ST258_GROUPING_AUDIT.md`. It changes no runtime behaviour and makes
no performance or outbreak-validation claim.

## Question and estimand

The primary question is whether a method gives a useful, manageable collection
of close-relative candidates around a query genome in a defined searched
collection. Independent epidemiological links test whether known relevant
records are recovered; they do not exhaust the set of genetically close
relatives. Direct transmission, direction of spread and universal outbreak
membership are outside the primary estimand.

Keep three objects distinct: a query's retrieved neighbourhood, an undirected
pairwise relationship, and a connected component formed by chaining such
relationships. A prefix neighbourhood and a raw-distance neighbourhood can be
compared as retrieval rules without claiming that their members form the same
kind of biological cluster. A good result can be a useful complementary view,
several incomparable rules, or insufficient evidence for any change.

## Minimum experiment and freeze

A feasible first experiment uses existing downloadable, scheme-compatible
collections and already available SNP evidence. Target at least three
independent studies: one development study and two wholly held-out studies
with independently documented events, preferably different sites and periods.
This is a feasibility benchmark, not an adequately powered universal validation.
If fewer than two held-out studies have informative positive and negative
comparisons, report a descriptive pilot and stop short of method promotion.
An additional ST beyond ST258 is a scope requirement for any claim across
lineages; otherwise scope conclusions to the sampled lineage.

Before producing method outputs, save a dated, checksummed protocol and manifests:

- Study/event eligibility, search dates and sources; explicit reasons for
  inclusion/exclusion; all eligible genomes rather than only the published
  genetically tight subset. Document when original ascertainment used WGS.
- Accession, assembly/read identity and checksum, patient/episode identity where
  available, dates with precision, site/ward, lineage, and sampling frame.
- Labels with supporting source passages/pages, adjudicator, confidence and
  unresolved disagreements; a separate blinded label file.
- Exact locus catalogue and checksum, allele definitions/reference snapshot,
  caller/version/settings, cgLIN source and assignment snapshot, and any SNP
  pipeline, reference genome, callable sites, masks and tree provenance.
- Development/test study assignment, rules below, frozen query/background sets,
  metric definitions, random seeds and planned deviations policy.

Split by the highest shared recruitment/investigation unit: publications reusing
one cohort, multi-site outbreak studies and repeated patient records remain in
one partition. Check duplicate accessions, assemblies and patients across
studies before freezing. Test queries search only their held-out study's frozen
background for the primary analysis. An external broad-context search is a
separate descriptive experiment with unknown labels, not additional negatives.
Do not use the existing ST258 pilot to select a winner or to certify its chosen
resolution. It is a development/engineering example if reused at all.

## Independent epidemiological labels

Two reviewers, ideally including an epidemiologist, apply rules without seeing
cgLIN, allele distances, SNP distances or tree positions. Preserve disagreements
as unknown unless resolved from independent records. Genomic evidence in a
paper may identify the collection, but cannot supply its benchmark truth.
An event name inferred solely from a genomic cutoff is not an independent label.

Use pair-level labels with reasons:

- **Supported link:** independently documented shared exposure/source or
  epidemiological connection within a defined investigation, with compatible
  dates. Retain direct documented contact, common-source exposure and broad
  event membership as separate evidence strata; the latter does not establish
  a direct pairwise transmission link.
- **Supported separate event:** affirmative evidence of different exposure
  chains/source events, or a prospectively defined and investigated comparison
  cohort with explicit evidence against the specified event association.
  Same ST controls are especially useful, but matching ST is not a label.
- **Unknown:** absent contact records, geography/date separation alone,
  “sporadic” without its ascertainment definition, genomic-only assignments,
  ambiguous patient identity or discordant independent accounts.

Never treat no recorded link as a negative. Supported separation concerns the
specified epidemiological association; such a pair may still be genetically
close. Record how comprehensive the investigation was. Run a high-confidence
label analysis first and a separately labelled broader event-membership
sensitivity analysis. If every purported negative depends on genomic evidence,
there is no valid specificity estimate.

Use one preselected isolate per patient episode for primary metrics, selected by
earliest date and accession tie-break before genomic inspection. Retain serial
isolates in a secondary analysis. Report when patient/episode identity is
unavailable, since this limits independence and the de-duplication check.

## Fixed candidate rules

All rules use the same eligible query/background identities wherever possible.
Exclude the query itself and same-episode duplicate from its retrieval set.

| Candidate | Frozen rule | Interpretation |
| --- | --- | --- |
| Reported LIN5–LIN10 | Retrieve background records sharing the query's complete reported prefix at each of levels 5, 6, 7, 8, 9 and 10, with compatible namespace/snapshot. | Six nomenclature retrieval rules; no nominal diameter guarantee. Unresolved prefixes abstain at that depth. |
| Raw cgMLST thresholds | Retrieve records with observed mismatch count at most 0, 1, 2, 4, 7 or 10 over jointly called compatible loci. | Six candidate pairwise rules using published fine-bin boundaries as a grid, not established outbreak cutoffs or equivalences to prefixes. |
| Raw nearest baseline | Retrieve all minimum-count ties, plus separately the top three ranks including all boundary ties. | A relative retrieval baseline; nearest available can still be distant. Report denominator and tie burden. |
| Existing SNP distances | Rank neighbours using an existing reproducible, within-study SNP distance matrix and retain the same nearest/top-three tie rules. | Comparator on the common eligible subset. No new SNP cutoff is selected from held-out data. |
| Existing SNP topology | Describe whether recovered sets occupy a supported clade and document mixed/missed event tips. | Secondary concordance evidence, not truth labels or a fitted clustering rule. |

Do not choose a SNP clade or support threshold after viewing epidemiological
labels. Record tree support values as supplied. If no reproducible SNP matrix or
suitable tree exists, mark this comparator unavailable; do not expand this first
experiment into a new alignment/inference project. An allele-derived tree is
not an independent SNP comparator. If published SNP cutoffs are explored, freeze
them from the original study's stated pipeline before test evaluation, label
that analysis secondary and make no cross-reference numerical equivalence.

For all threshold and prefix rules the primary neighbourhood is direct query
membership. As a secondary safety diagnostic, form connected components from
these relations and measure chaining: endpoints exceeding the tested threshold,
supported separate events combined, and supported links split. Do not silently
substitute components for direct neighbourhoods.

## Missing calls and reference compatibility

The primary allele analysis requires the canonical 629-locus universe to be
verified. Require each profile to have at least 599 called loci and each pair at
least 599 jointly called loci (at most 30 unavailable in the joint comparison).
The joint-pair requirement is an additional conservative design choice, not an
upstream cgLIN assignment rule. Preserve mismatch count, jointly called count,
individual missingness and reasons for every excluded/abstaining comparison.

A mismatch count of zero means no observed differences at jointly called loci;
it does not mean identical profiles or genomes. Never count missing calls as
matches, normalize counts to 629 to reassign a threshold, or translate a cutoff
between schemes. If canonical identity cannot be recovered, retain the collection
as a separately named observed-locus descriptive dataset and exclude it from the
confirmatory candidate-grid comparison. This applies even if imported codes
look complete.

Report a complete-joint-call sensitivity subset (629 jointly called loci).
Also bound missingness effects without imputation: with observed mismatches m
and U unavailable canonical loci, full-profile distance is between m and m+U.
For threshold t, classify a pair as certainly included if m+U <= t, certainly
excluded if m > t, and otherwise indeterminate; report retrieval bounds.
These bounds supplement the operational observed-count rule. They are likely
uninformative at small thresholds, which is itself an important limitation.

Use one verified cgLIN assignment snapshot per analysis, or stratify snapshots;
never pool incompatible labels as if components matched. A common caller name
or advertised job is insufficient provenance. Recompute incompatible allele
profiles only if the required inputs are already available and that step is
frozen before test unblinding. SNP distances are comparable within a common
reference, mask and callable-site definition; report separate reference strata,
use matched eligible subsets, and never concatenate numbers from different
references. Show full-cohort coverage alongside every matched-subset result.

## Outcomes and analysis units

The primary outcome pair is **recovery of supported linked candidates** versus
**review burden**. For each query with at least one eligible supported link,
calculate the fraction of its eligible supported links retrieved. Average first
within event, then within study; present each study and the equal-study macro
average. Report eligible linked queries and unrecoverable links excluded for
coverage separately. For all eligible queries report median, quartiles and
maximum neighbourhood size, unknown-label fraction and the number of supported
separate-event records retrieved. Show the recovery/burden frontier across all
fixed candidates; large retrieval sets do not receive credit for size alone.

Secondary outcomes are pair sensitivity on supported links, supported-separate
pair exclusion rate, missed-link and false-merge counts, nearest-tie retention,
abstention/coverage, and component fragmentation/chaining. Count each unordered
pair once for pair summaries. Report patient/event/study counts as well as pair
counts. “Specificity” means exclusion of supported separate-event pairs under
these label rules, not proof that retrieved records are close or distant.

Precision among known-labelled pairs may be tabulated only with its labelled
coverage and ascertainment caveat. Overall precision/PPV is not identified when
unknowns are substantial or controls were selected artificially. Give bounds on
retrieved supported-link proportion by assigning unknowns in turn to each label;
do not replace unknowns with negatives. Results compare epidemiological recovery
and review burden, not a gold-standard genetic-relatedness classifier.

## Dependence and uncertainty with small samples

Pair counts are not independent sample sizes. Use study as the highest resampling
unit, preserving all patients, events and pairs within each sampled study. Use
paired resampling for candidate differences and give study-specific differences.
With at least ten independent held-out studies, a prespecified 10,000-replicate
study bootstrap can provide exploratory percentile intervals for equal-study
macro metrics. This ten-study rule is pragmatic, not a guarantee of coverage.

For the minimum two held-out studies, do not produce apparently precise pairwise
binomial intervals or rely on an unstable cluster bootstrap. Show both study
estimates, their range, event-level estimates and leave-one-study-out influence;
state that population uncertainty is poorly identified. If an event-level binary
endpoint such as “all supported links recovered” is reported, an exact binomial
interval is appropriate only for genuinely independent event units; events
sharing a study/site are not automatically independent. With just two studies,
a study-level binary interval is valid under an explicit independent-study
sampling assumption but necessarily very wide. No significance testing across
the 12 candidate rules and no ranking by tiny point-estimate differences.

## Prespecified sensitivity analyses

- **Time:** use absolute sampling-date gap bins 0–30, 31–180, 181–365 and >365
  days, with unknown/imprecise dates separate. For year/month intervals, assign
  a bin only when the entire possible gap lies inside it. Report sparse cells
  rather than combine bins after seeing results. Compare recovery/burden over
  these strata; do not fit a slope-derived threshold. A chronological subset is
  secondary and cannot replace whole-study separation.
- **Density:** retain the frozen queries and downsample their background to
  50% and 25% within study/event 100 times using frozen seeds. Include the full
  sample result. Report variation in nearest identity, ties, recovery and burden;
  this variation is a perturbation distribution, not a population confidence
  interval. Measure only links whose partners remain available in each draw and
  separately report the loss of available linked partners.
- **Callability:** report the complete-joint-call subset, missingness strata and
  the bounds above. If controlled masking is run, use only complete profiles,
  mask 1%, 3% and 5% of calls with frozen seeds, and report inclusion changes and
  eligibility losses. Random masking does not model systematic locus failure.
- **Lineage and evidence:** show ST, epidemiological evidence strength,
  repeated-isolate inclusion and SNP-compatible-subset results. Treat strata as
  descriptive when counts are small.

## Decision rule and stopping point

Freeze practical tolerance for missed independently linked records, supported
separate-event inclusions and review workload with the intended reviewer before
unblinding. Record numerical limits and their operational rationale in the
protocol; without agreed tolerances, this experiment cannot authorize an
automatic production default. Neither a distance boundary nor a prefix level is
accepted merely because it reproduces the label construction.

For this minimum experiment, call a view **promising for further evaluation**
only if its recovery/burden improvement or complementary value is visible in
both held-out studies, all agreed tolerances hold in each, and scheme,
callability and chaining failures are explained. Report ties and tradeoffs
explicitly. A candidate that retrieves fewer genomes while dropping relevant
links is not automatically better. Absence of a clear improvement leaves the
existing descriptive options in place; uncertainty is an acceptable result.
Automatic default selection and outbreak claims require a larger independent
validation and an explicit operational decision beyond this pilot.

The deliverable is one frozen identity/label manifest, an eligibility flow,
a per-query/per-pair results table with denominators, per-study recovery/burden
plots and sensitivity tables, and a deviations log. Keep method specification,
observed outcomes and implementation changes in separate documents. No new
benchmark outcomes were inspected for this design review.
