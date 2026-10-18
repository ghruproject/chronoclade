# ST258 comparison and sampling robustness

User-authorised objective, 18 October 2026: compare allele profiles and correctly
aligned cgMLST-locus nucleotide sequences, test sensitivity to context sampling,
attach observed collection-time summaries, and explain the findings plainly.
All work stays local. This objective is tracked by an active app goal.

## Dataset and evidence

Use the cached 7,815-genome dataset: 41 Greek input genomes and 7,774 public
comparisons. Preserve the original provider schema and unknown historical
database fingerprint. Resolve allele DNA only with an auditable identity bridge;
do not assert that the current Pasteur catalogue is the historical provider
snapshot. The 628 observed exported loci are an explicitly incomplete locus
universe, not a newly claimed complete 629-locus scheme.

Primary comparisons use the same genome identities and pairs with at least
597 jointly available loci in each representation (ceiling of 95% of 628).
Record sequence-unresolved calls and pair exclusions. Missing DNA is unavailable
evidence, never a zero distance. Align observed allele sequences separately per
locus, reconstruct and concatenate genomes; count joint A/C/G/T substitutions
and retain separate indel evidence. No recombination filtering is implied.

## Analysis

- Use exact single-linkage partitions. An exact minimum spanning forest can
  compactly represent all cuts; it must preserve zero and unavailable distances.
  Report internal distances separately from joining thresholds.
- Display descriptive raw-distance cuts 0, 1, 2, 4, 7, 10, 15, 20, 30 and 50,
  with 100 and 200 additionally available for locus SNPs. Equal numbers are not
  equivalent biological boundaries between the two measurements. Select no
  automatic outbreak cutoff or winning method.
- Distinguish context genomes directly within a cutoff of an input from those
  linked through chains. This additional descriptive count was specified before
  inspecting the new SNP matrices; native allele connectivity was already known.
- Keep all 41 inputs in every sampling run. Retain 25%, 50% and 75% of contexts
  in five alternative simple random samples each, using fixed seeds 17, 31, 43,
  59 and 71. Recompute the induced-subset forest, rather than pruning the full
  forest. This tests removal of sampled background, not real-world sampling bias.
- Compare each retained input neighbourhood with its full-data neighbourhood
  restricted to the same retained genomes. Report splitting and retained-peer
  recovery, not a Jaccard penalty merely for genomes deliberately removed.
- Rank input-to-context neighbours by actual raw distances, retaining ties and
  displaying coverage. Compare overlap of nearest sets without fitting cutoffs.
- Attach location and observed date intervals to genetic groups. Use the existing
  >=1900 date policy for time summaries, reporting invalid/missing/excluded dates.
  Do not use dates to define genetic groups or fit an allele clock.

## Acceptance and presentation

Verify native allele-array identities; independently check aligned SNP counts on
a declared subset, all focal-to-context pairs where practical, and zero/indel
cases. Check exact forest partitions against independent graph traversal on
smaller subsets. Keep numerical receipts, seeds, source hashes and method notes.
Measure runtime, memory and disk growth on this laptop; do not substitute a small
panel while calling it the full dataset.

The report must answer, in this order: what was analysed; what agrees and differs
between allele and nucleotide views; which neighbourhoods depend on sampling;
where and when their members were sampled; what this does and does not establish.
Use the existing ChronoClade style, compact figures and searchable tables opening
at ten rows. Put detailed methods and machine-readable evidence behind downloads.
Use accessions or sample labels, with internal keys restricted to audit files.
Explicitly state that this unlabelled public context cannot establish transmission,
the origin of a lineage, or a universally valid epidemiological cutoff.

## Completed run

The full-panel comparison and independent checks are complete. See
[results and verification](ST258_COMPARISON_RESULTS.md) and
`tmp/science-st258/comparison/acceptance.json`. No universal cutoff was selected.
