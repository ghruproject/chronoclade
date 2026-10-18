# Klebsiella relatedness pilot: protocol 1

Declared 18 October 2026 before computing the new Zhou-panel profiles or distance
results. Existing ST258 engineering and Bowers capsule-label results are known;
this is not a blinded confirmatory study. Save the protocol SHA256 beside the
input freeze before running the evaluator. Amendments must precede their affected
analysis and retain the original record.

## Primary question

How do reported LIN prefixes and simple measured allele-distance neighbourhoods
recover author-defined outbreak isolates while admitting explicitly unrelated
same-ST controls? Measure neighbourhood recovery and review burden. Do not claim
direct transmission, geographical origin, an optimal universal cutoff, or that
absence of an epidemiological record proves unrelatedness.

## Development dataset and labels

Use all 31 strains in Table 1 of Zhou et al. (2017), DOI
10.3389/fmicb.2017.00371: 26 outbreak-labelled strains and five author-defined
sporadic controls, all ST11. The authors used epidemiological information plus
PFGE/MLST in defining the outbreak; those labels are not independently adjudicated
patient-contact truth. Preserve this ascertainment limitation.

Freeze Table 1 metadata, exact strain/BioSample/assembly-version joins, source
URLs and bytes before typing. Case-only aliases require an explicit crosswalk;
no fuzzy name matches or replacement genomes. Dates come from Table 1, with their
precision retained. Unknown patient IDs/hosts stay unknown. Hospital letters are
site labels, not cities. Do not infer country from an author's affiliation.

For the primary pair labels:

- Two outbreak-labelled isolates: same documented event, positive for event
  recovery. This is not proof of direct patient-to-patient transmission.
- Outbreak versus an explicitly sporadic isolate: author-defined event-negative
  comparison. Five controls represent five isolates, not 130 independent negatives.
- Two sporadic controls: unknown for the primary analysis; do not make extra
  negative observations from their combinations.

This is one development event with five controls, not a held-out validation set.
NIH2011 is a candidate independent positive-event panel, pending exact sample and
patient mapping. Bowers capsule labels are genetic proxies only. No method is
promoted to a production default from this pilot.

## Typing and eligibility

Run all 31 current versioned assemblies through ChronoClade native preparation,
using one fixed, ready Pasteur 629-locus reference/caller/LIN snapshot. Record
assembly SHA256, runtime/tool versions, the complete ordered locus catalogue,
reference configuration and versions/fingerprints in output evidence. Reuse the
existing installation; do not call this a fresh reference installation.

The paper used 1,143 targets. We recompute its genomes under the 629-locus scheme;
we are comparing methods on its labelled collection, not reproducing its numerical
allele cutoff or translating thresholds between schemes.

Require the same complete 629-locus catalogue across genomes. Primary pair
comparisons require at least 599 jointly called loci. This is a conservative
analysis coverage policy, not an official cgLIN pairwise rule. Report all exclusions
and per-genome/per-pair callability. Missing calls never count as agreement or
mismatch. Unresolved/provisional LIN positions stay unavailable; never pad codes.

## Fixed candidate rules

1. Shared valid reported LIN prefix at levels 5, 6, 7, 8, 9 and 10.
2. Raw categorical mismatch distance at most 0, 1, 2, 4, 7 or 10, subject to the
   shared-locus rule. The grid is fixed in advance from nomenclature boundaries;
   it does not inherit the meaning of a LIN prefix or validate those cutoffs.
3. All minimum-distance ties as a descriptive nearest-neighbour baseline.

Evaluate direct query neighbourhoods first. Connected-component chaining, if
computed, is a secondary diagnostic and must report its maximum endpoint distance.
Do not treat prefix groups as complete-linkage clusters.

## Outcomes

Report every rule, not only the apparent best:

- Recovered same-event records per outbreak query / its 25 eligible peers.
- Explicit control isolates admitted per outbreak query, with their identities.
- Neighbourhood size and retained nearest ties, including cross-prefix exceptions.
- Pair counts and raw mismatch/shared-locus distributions for positive, negative
  and unknown labels; unassigned LIN and missing-call exclusions separately.
- Specificity/recall summaries only with their actual single-event/five-control
  denominator and no misleading independent-pair confidence intervals. Do not
  tune a winner on this collection.
- Descriptive stratification by known collection-date separation (0–30, 31–90,
  over 90 days), preserving intervals when dates lack day precision. No
  time-dependent cutoff is estimated from a root-to-tip slope.

Coverage sensitivity at 90%, 95% and 100% jointly called loci is secondary, with
all denominators explicit. A sensitivity result is not a new preferred rule.

## SNP comparison and limits

If resources permit, use the standard strict ChronoClade tree route on the exact
31-genome cohort and preserve all identities. Compare corrected SNP differences,
callable sites and tree relationships with the profile neighbourhoods. If strict
coherence or acquisition fails, retain the failure; do not remove a control to
obtain a favourable tree. A predeclared separate 26-outbreak-only tree may be
reported as a different scope and cannot assess control separation.

SNP relationships are complementary genetic evidence, not independent
 epidemiological truth. The existing Bowers 101-genome corrected tree may support
a separate descriptive method comparison, explicitly reusing prior results.

## Decision

A useful pilot can expose missed event peers, admitted sporadic controls,
chaining, or disagreement between LIN and actual distances. It may identify
candidate views for external evaluation. It cannot establish a universal cutoff,
time-aware clustering, general sensitivity/specificity, or transmission direction.
Keep current explicit scopes/defaults unchanged. A held-out study and additional
lineages are required before proposing a default; a negative result is acceptable.
