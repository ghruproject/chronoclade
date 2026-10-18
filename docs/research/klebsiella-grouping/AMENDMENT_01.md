# Amendment 1: retain controls for a diagnostic SNP comparison

18 October 2026. Declared after the strict tree run stopped and before inspecting
SNP pair counts. The primary allele/LIN protocol and results remain unchanged.

The standard strict run on the exact 31 genomes completed alignment, initial tree
and recombination correction. Its corrected-distance coherence screen flagged all
five author-defined sporadic controls and none of the 26 outbreak isolates. It
stopped rather than publishing a successful tree result. Preserve this run and
its failure; do not reinterpret it as successful strict acceptance.

For the scientific comparison, run the same fixed 31-genome selection in a fresh
output directory with the supported `--coherence-policy diagnostic` option.
Retain every control, every input, the same seed, references and four threads.
No assemblies are replaced, no dates are filtered and no native result is copied
into a new completed manifest. The normal command recomputes its artifacts.

This secondary analysis is exploratory and flagged for cohort diversity. Compare
corrected SNP counts, callable sites and topology with allele/LIN neighbourhoods.
It is not evidence of model adequacy, valid transmission reconstruction, an
outbreak cutoff or passing default quality checks. Do not assess temporal signal
or fit calendar dates for this heterogeneous diagnostic cohort; use
`--no-assess-temporal --randomisations 0`. Collection-date strata remain the
predeclared descriptive analysis.

Run bounds: 20 minutes, four native threads, 8 GiB process-family RSS and 2 GiB
new output for this diagnostic run. If it fails or exceeds these practical
bounds, retain the failure and report SNP comparison unavailable; do not prune
controls to obtain a preferred result.
