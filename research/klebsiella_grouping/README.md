# Descriptive Klebsiella retrieval evaluator

This standalone research utility fits no thresholds and changes no runtime
behaviour or defaults. Do not run a real benchmark before its protocol, sample
identity manifest and independently adjudicated labels are frozen. A single
study is a descriptive within-study pilot, not external validation.

From the repository root:

```sh
pixi run python -m research.klebsiella_grouping.evaluate \
  --dataset /absolute/path/prepared/manifest.json \
  --labels /absolute/path/frozen-labels.csv \
  --output /absolute/path/new-results-directory
```

`--min-shared-loci` defaults to 599. Protocol-declared sensitivity runs use
567 (ceiling of 90% of 629), 598 (ceiling of 95%) and 629 (100%), each with a
new output directory. `summary.json.coverage_policy` records the integer cutoff,
canonical denominator 629 and actual fraction. No result selects a new policy.

The output directory must not already exist. Inputs are read only. `--dataset`
is a normal `chronoclade.prepared-dataset` v1 manifest accepted by core
`load_dataset`; its checked artifacts are loaded through that API. Exactly one
`AlleleMatrix` must cover exactly all sample IDs. Its `LocusCatalogue` must have
`complete=true`, exactly 629 distinct ordered locus names and a non-null
`database_sha256`. The investigator must establish that these are the intended
canonical loci and frozen reference: counting 629 names cannot authenticate a
catalogue. Mixed schemes/reference snapshots must be evaluated separately.
Zero matrix code means missing; it never contributes an allele match.

Required labels CSV columns (exact spelling):

| Field | Contract |
| --- | --- |
| `sample_id` | Exactly one row per prepared sample ID; no extra or omitted IDs. |
| `study_id` | Nonempty independent study/cohort identifier. |
| `event_id` | Required for `outbreak`; same documented event within a study establishes the positive relationship. Blank for unassigned controls. |
| `epi_status` | One of `outbreak`, `control`, `unknown`. |
| `control_for_event` | Only permitted for `control`: exact event identifier for which an independent source explicitly supports the negative comparison; otherwise blank. A generic “sporadic” name is insufficient without that evidence. One event per row is supported. |
| `source` | Nonempty independent evidence citation/location. Investigator checks independence from genomic cutoffs; the utility cannot adjudicate prose. |
| `patient_id` | Optional; equal nonempty values within a study exclude the pair from retrieval. Prefer one frozen isolate per patient/episode. Unknown IDs do not assert independence. |

Blank `control_for_event` leaves outbreak–control pairs unknown. Control–control,
different-event outbreak pairs and other unsupported relationships are unknown.
Cross-study pairs are reported but excluded from primary retrieval. Known
positive pairs excluded for missingness or patient identity remain visible as
`excluded_links` for each query. `source` is a provenance requirement, not
machine verification of epidemiological independence.

LIN inputs are prepared `lineages` rows with `kind="cglin"`, `assignment` (the
upstream code syntax accepted by core `resolved_prefix`), `resolution_status`,
`scheme_id`, `scheme_version`, and `database_version`. When database version is
absent, the preserved `evidence.cglin_database_sha256` (the native plincer
snapshot key), `evidence.cglin_reference_sha256` or
`evidence.cglin_frozen_export_sha256` can identify a common frozen snapshot.
Multiple cgLIN assignments per sample fail rather than being arbitrarily
chosen. Missing/unknown scope abstains. Status must be `resolved`, `complete`,
`partial` or `provisional`; each depth must be independently resolved by the
core prefix API. Provisional assignments require the core native prefix-consensus
evidence; no missing components are filled. Scope identity must match for a pair.
Unassigned LIN is reported separately and never assigned an artificial prefix.

The fixed methods are direct query retrieval using shared LIN5–LIN10 prefixes,
and raw pairwise observed mismatch thresholds 0, 1, 2, 4, 7 and 10, plus all
minimum-distance ties as a directed per-query nearest baseline. These are
candidate boundaries, not validated outbreak cutoffs. Pair eligibility requires
at least the chosen shared-locus cutoff (primary 599) and the same study, excluding known same-patient
pairs. No missing-call scaling, imputation or connected-component chaining is
performed. All methods use that same eligible background. LIN abstentions count
as unrecovered links against the common eligible-link denominator, and their
counts are shown. With no jointly called loci, mismatch count is blank (`None`
in memory), never a false zero; low-coverage pairs retain their observed raw
counts but cannot be retrieved.

Outputs:

- `samples.csv`: per-record coverage, eligibility and LIN assignment availability,
  retaining label provenance.
- `pairs.csv`: each unordered pair once; observed count, joint denominator,
  exclusion reason, independent label and each candidate relationship. LIN
  equality may be recorded for an excluded pair but it is never retrieved.
- `queries.csv`: candidate recovery of eligible supported links, retrieval burden,
  negative/unknown retrievals, excluded links and abstention counts. Recovery is
  blank if no eligible supported links exist. All prepared records serve as
  queries; `role` does not silently remove controls from descriptive queries.
  Includes the directed `nearest_all_ties` baseline, retrieved/control identity
  JSON lists, minimum distance and available/retained nearest-tie counts.
- `events.csv`: study/event/candidate record counts, linked query counts,
  mean per-query recovery, median/max burden. Blank-event controls form a
  descriptive stratum, not an inferred event. These are sample-record counts;
  they are not asserted to be independent patients.
- `candidate_counts.csv`: unordered eligible/retrieved/abstaining pair counts,
  with retrieved positives, negatives and unknowns. No independent-pair CIs.
  The directed nearest baseline is deliberately excluded from this table.
- `summary.json`: dataset identity, input checksums, sample/study/event counts,
  pair eligibility/label totals and reference identity. These are descriptive
  counts; there is no statistical ranking, significance test or fitted winner.

`pairs.csv` gives descriptive date-gap bins 0–30, 31–90 and >90 days
only when both prepared `collection_date` fields are valid exact ISO dates.
Partial/unknown dates yield `unknown`. There is no slope or time-based threshold.
No precision/PPV or specificity is asserted from incomplete label ascertainment.
SNP comparisons, whole-study validation and uncertainty analysis remain outside
this small evaluator; see the preregisterable design review.

## Exploratory extension: three distances and threshold sweeps

`docs/research/klebsiella-grouping/AMENDMENT_02.md` declares the local comparison
of allele mismatches, nucleotide differences inside the same cgMLST gene panel,
and recombination-filtered core SNPs. These are distinct measurements, each with
its own denominator. Neither their thresholds nor their rates are interchangeable.

The shared `chronoclade.cgmlst.thresholds.threshold_sweep` function accepts IDs
and raw weighted undirected edges for either distance type. It sorts once and
incrementally joins components, treats tied edges simultaneously and preserves
isolated genomes. Missing edges are unavailable evidence, not zero distance.
Outputs include requested partitions, compact merge history and ranges of
unchanged membership across requested cutoffs. These ranges do not establish
epidemiological validity or persistence through time. Explicit partition output
uses O(samples × requested thresholds) space; edge storage is O(available pairs).
This implementation is not yet a benchmarked solution for dense 10,000-genome
distance matrices.

For the frozen 31-genome research panel:

```sh
pixi run python -m research.klebsiella_grouping.threshold_report \
  --primary tmp/science-zhou2017/evaluation-primary \
  --locus-pairs tmp/science-zhou2017/locus-dna/pairs.csv \
  --metadata tmp/science-zhou2017/prepared/samples.csv \
  --labels validation/science/zhou2017/labels.csv \
  --output tmp/science-zhou2017/threshold-comparison
```

The output includes a self-contained `report.html`, memberships and all integer
threshold counts as CSV, complete sweeps as JSON, and an input/code hash receipt.
Joining thresholds and internal pairwise diameters are reported separately:
single linkage permits chains whose endpoint distance exceeds the cutoff.
Collection-date intervals and reported country/region/city/host accompany genetic
groups; metadata does not change group membership. No threshold is automatically
selected, and no product report or default is replaced by this research runner.

Synthetic verification only:

```sh
pixi run python -m pytest -q research/klebsiella_grouping/test_evaluate.py
```

## Aligned locus method correction

The normal direct-DNA route now aligns observed alleles independently per locus
with MAFFT, reconstructs each genome and concatenates the aligned loci. Historical
`locus_snps.py` and positional reports are retained only as superseded evidence.
Use `aligned_locus_snps.py` for the corrected research comparison and pass its
`locus-aligned/pairs.csv` to `threshold_report --locus-mode aligned`.
The scalar per-locus and vectorised concatenation calculations must agree; gaps
and missing evidence are separate from SNPs. Saved partitions retain all canonical
loci, including explicit absent records for a globally missing gene.
