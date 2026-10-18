# Neighbourhood analysis: implementation and acceptance

18 October 2026. Local work; no scientific defaults or GitHub publication.

## Current delivery

- Shared raw-distance single-linkage threshold engine, independent of allele or
  nucleotide representation. Incremental unions, deterministic tied edges,
  explicit unavailable pairs and singleton retention.
- Matched 31-genome ST11 comparison of allele counts, cgMLST-locus nucleotide
  counts and corrected core SNPs. Sequence identities and coverage must be
  auditable; positional gene comparisons must retain their limitations.
- Research report in the existing ChronoClade style: threshold behaviour,
  membership sets, internal diameters, source labels and observed metadata.
  Time is a separate descriptive layer and does not determine genetic groups.

Acceptance: all pair identities match the frozen panel; independent calculations
reproduce threshold partitions; nucleotide numerators and denominators are
verified; missing/unequal-length loci remain explicit; generated figures and
evidence links work in the browser. No automatic best cutoff or clock inference.

## Next scientific gates

1. Freeze a second study with defensible event/control labels before applying the
   unchanged methods. Report failures and overlap, not just clean separations.
2. Evaluate added background and sample removal: does a neighbourhood retain its
   membership and separation when sampling changes? Keep threshold insensitivity
   separate from robustness to sampling.
3. Assess the locus-nucleotide method against a properly aligned gene panel,
   including length-changing and compensating indels. Decide its supported
   alignment policy before promoting it as a general SNP workflow.
4. Validate representative alternative subsamples against full-panel diversity,
   closest SNP neighbours and downstream topology/date-estimate variation.
5. Integrate the reviewed neighbourhood view into the normal cgMLST reports for
   both distance routes. Reuse the engine and metadata summaries; do not retain
   separate implementations for allele and nucleotide reports.

Dense large-panel performance needs its own measured gate. The initial engine
stores available edges and requested partitions; a 31-genome pilot does not
demonstrate feasibility for 10,000 dense genomes. Compact merge histories and
exact sparse/MST representations are candidates for that separate optimisation.

Time-aware enrichment and evolutionary cutoff models remain research questions.
They require suitable longitudinal/background data and validation, not an assumed
alleles-per-year conversion from this short outbreak.
