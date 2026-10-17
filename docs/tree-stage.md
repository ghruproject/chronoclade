# Independent assembly tree

`chronoclade tree <selection.json> --out tree/ --threads 4` consumes a validated
cgMLST selection. It checks every source checksum and builds precisely the
ordered selected IDs. It never retrieves other context or chooses replacements.
Linked local assemblies are validated and snapshotted independently where copy-on-write
is supported. Missing assemblies are acquired
by exact Pathogenwatch `source_genome_id` using the configured API key; credentials
are never saved. Selected public downloads verify exact source identity, upstream
SHA-1 and genome length when available; older caches are upgraded only after verification. A
per-source lock prevents duplicate downloads across concurrent jobs. Output copies
use the same safe copy-on-write snapshot helper as selection evidence, with a
regular-copy fallback and no hardlinks. Acquisition failure leaves an explicit `assembly_audit.json` in
the failed job and preserves any previous complete result.

SKA reference mapping, IQ-TREE GTR+G and ClonalFrameML reuse the existing scientific
engines. Native trees choose the shared input assembly with greatest N50, with a
stable ID tie-break, and pass the requested seed to IQ-TREE. Alignment failures
fail the exact selection; relative-distance coherence flags follow the explicit
strict or diagnostic policy below. Neither policy silently drops genomes. The saved result records raw and corrected trees,
filtered alignment and callable-site counts, executable hashes/package versions,
reference, importation/filtering evidence, SNP pair distances and the original
selection/dataset hashes. Undated samples remain in the phylogeny.

`tree.json` is published last; artifacts live in fingerprinted job directories.
An unchanged complete run validates and resumes without rebuilding. Parameter,
local assembly or executable changes invalidate that job. `--force` creates a
fresh job without deleting the previous complete result. A failed replacement
never publishes a replacement `tree.json`.

The offline report retains the full-report style, corrected biological tree,
country/host metadata including missing categories, exact location-history network
and temporal results. SNP neighbours have SNP/site evidence and do not inherit
cgMLST allele-distance labels. The corrected tree uses substitutions/site.

## Temporal gate

The default runs 20 date permutations with the least-squares root refitted for the
observed clock and every permutation. Maximal monophyletic clades sharing the
same exact collection date string are the exchangeability units, following
[Murray et al. 2016](https://doi.org/10.1111/2041-210X.12466). Tree topology is inferred
without dates; same-date clades are then identified on that frozen topology. One
date is permuted per cluster and expanded to its members. Year/month/day strings
are retained, permitting TreeTime's partial-date interpretation. Cluster date
multiplicities are preserved; individual-tip date multiplicities can change with
unequal cluster sizes. Tip weights remain unchanged in clock fitting.

The saved assessment includes every eligible ID, excluded undated IDs, actual
cluster memberships, seed, requested/successful refits, root policy and a
descriptive pairwise genetic-distance/date-gap correlation. That correlation is
not an independent-pair significance test. Positive rate and the existing
root-to-tip permutation threshold gate are required; incomplete permutations,
fewer than three distinct cluster dates or disabled assessment remain unsupported.
A passing clustered screen does not establish model adequacy or remove every
sampling bias. This is not the full-model CR2 uncertainty-overlap test.

`--no-assess-temporal` skips clock/permutation tools. Fewer than three usable
distinct dates writes a complete tree product with an unsupported assessment.
Tree never creates a time-scaled phylogeny; it saves `assessment.json` for time.

### Relative-distance coherence policy

Native `chronoclade tree` defaults to `--coherence-policy strict`, which stops on extreme relative-distance flags before or after recombination correction. These flags can reflect legitimate cohort diversity; they do not establish incorrect genome identity. For an intentionally broad selected cohort, explicitly use `--coherence-policy diagnostic` to retain every selected genome and continue with a visible report warning. Neither policy prunes samples.

Both screens retain TSV and JSON audits, including thresholds, cohort median, median absolute deviation, every sample statistic and flagged statistics. Completed `tree.json` embeds both audits and the policy. The time report carries forward any coherence warning from its source tree. The policy is bound into single-tree resume fingerprints and tree-ensemble parameters. Legacy `run` retains the strict default.

Native tree `--seed` controls both IQ-TREE inference and temporal permutations. The manifest records `phylogeny_seed` and `phylogeny_seed_policy: explicit_IQ_TREE`; these are bound into the resume fingerprint so older unseeded inference is not reused as seeded output. Legacy core callers without a supplied seed retain IQ-TREE default seed behavior.

Native assembly trees select the highest-N50 reference from the input genomes, with deterministic sample-ID tie breaking. The shared inputs remain present in every alternative selection, so varying public context cannot change the mapping reference or its sequence hash. This fixes the reference basis across ensemble alternatives, although callable columns can still differ with cohort composition. The recorded `reference_selection_policy: highest_N50_common_input` is bound into both tree and ensemble fingerprints. A selection without input genomes fails. Legacy lineage reference selection remains unchanged.
