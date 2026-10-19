# Independent assembly tree

`chronoclade tree <selection.json> --out tree/ --threads 4` consumes a validated
cgMLST selection. It checks every source checksum and builds precisely the
ordered selected IDs. It never retrieves other context or chooses replacements.

The mapping reference defaults to the highest-N50 input assembly. Set
`--reference-sample-id SAMPLE_ID` to use a specific selected input instead; for
an ensemble, that input must occur in every selection. This chooses the
biological sequence used for mapping, not a root or an outgroup. The reference
identity and policy are recorded and checked on reuse. The same option is
available for Slurm tree arrays.

Linked local assemblies are validated and snapshotted independently where copy-on-write
is supported. Missing assemblies are acquired
by exact Pathogenwatch `source_genome_id` using the configured API key; credentials
are never saved. Selected public downloads verify exact source identity, upstream
SHA-1 and genome length when available; older caches are upgraded only after verification. A
per-source lock prevents duplicate downloads across concurrent jobs. Output copies
use the same safe copy-on-write snapshot helper as selection evidence, with a
regular-copy fallback and no hardlinks. Acquisition failure leaves an explicit `assembly_audit.json` in
the failed job and preserves any previous complete result.

SKA reference mapping and IQ-TREE GTR+G are followed by the selected recombination
filter. The default is `--recombination-method phipack`, a fast exploratory
PHIPack Profile scan followed by window masking and a new IQ-TREE inference on
the filtered alignment. Use `--recombination-method clonalframeml` for the
branch-specific ClonalFrameML analysis used by earlier versions and benchmarks.
Use `--recombination-method verticall` for Verticall's reference-relative
pairwise and alignment-mask workflow, followed by a rebuilt IQ-TREE phylogeny.
Native trees choose the shared input assembly with greatest N50, with a
stable ID tie-break, and pass the requested seed to IQ-TREE. Alignment failures
fail the exact selection; relative-distance coherence flags follow the explicit
strict or diagnostic policy below. Neither policy silently drops genomes. The saved result records raw and corrected trees,
filtered alignment and callable-site counts, executable hashes/package versions,
reference, importation/filtering evidence, SNP pair distances and the original
selection/dataset hashes. Undated samples remain in the phylogeny.

PHIPack uses 10,000-base scanning windows, a 5,000-base step and a 100-base PHI
neighbourhood. `--phipack-window`, `--phipack-step` and `--phipack-threshold`
control the first two values and the Benjamini–Hochberg FDR threshold (default
0.05). Significant whole windows are removed from every sequence; incomplete
columns are also removed. Scans never cross reference-record boundaries, and
include each record's final partial step. Tests run concurrently up to `--threads`.
The saved `phipack/windows.csv` records mapped and per-record coordinates,
raw p-values, adjusted q-values and masking decisions; `phipack_summary.json`
records coverage and limitations. Windows with an unavailable PHI approximation
are reported as untested and retained. A lack of significant windows does not
establish absence of recombination. This is exploratory region filtering rather
than ClonalFrameML import inference; downstream dating retains that provenance.
The same options are available on `slurm-ensemble tree`, and all are bound into
the saved plan and result fingerprints.

Verticall runs each reference FASTA record independently because its alignment
mask expects a single reference sequence. Query assemblies are snapshotted with
contigs split at non-ACGT bases and unique temporary contig names; original
assemblies are unchanged. `verticall/assembly_repair.json` maps safe sample
aliases back to original IDs and records original/snapshot hashes and removed
ambiguity counts. The reference itself must be unambiguous; an unsuitable
reference is refused rather than silently replaced.

Each per-record pseudo-alignment includes an exact temporary reference row for
coordinate mapping. Only that auxiliary row is excluded from the native masked
output; all selected genomes, including the selected reference sample, remain.
Native first-result handling is explicit (`--multi first`), and omitted samples,
incomplete coordinate classifications or inconsistent masks stop the workflow.
The native sample-specific masked alignment is retained as
`verticall.masked.fasta`. For the existing tree/time callable-site policy, the
downstream alignment additionally removes columns not A/C/G/T in every selected
sample; this complete-site step is an explicit extension of the upstream
[alignment-tree workflow](https://github.com/rrwick/Verticall/wiki/Alignment-tree-workflow).
`verticall/regions.csv`, per-record TSVs/logs and `verticall_summary.json` retain
vertical/horizontal/unaligned calls, coordinate accounting and limitations.
This is reference-relative masking, not branch-specific import inference.

`tree.json` is published last; artifacts live in fingerprinted job directories.
An unchanged complete run validates and resumes without rebuilding. Parameter,
local assembly or executable changes invalidate that job. `--force` creates a
fresh job without deleting the previous complete result. A failed replacement
never publishes a replacement `tree.json`.

The offline report retains the full-report style, corrected biological tree,
country/host metadata including missing categories, exact location-history network
and temporal results. SNP neighbours have SNP/site evidence and do not inherit
cgMLST allele-distance labels. The corrected tree uses substitutions/site.
When usable host metadata is supplied, the state-network selector also offers
host changes reconstructed on that same corrected tree. Location remains the
default. Without host metadata the option is omitted; a single known host state
gets an explanatory note rather than invented links. Both views use the shared
interactive network and its export controls.

## Optional rapid branch support

Use `--support-replicates 1000` (or more) to run a separate IQ-TREE UFBoot2
inference with bootstrap NNI on the selected method's filtered alignment. The default
is `0`, disabled. The model is GTR+G; the requested seed and replicate count are
recorded and included in resume fingerprints for single and ensemble runs.

The supported Newick tree and IQ-TREE report are retained alongside a comparison
of unrooted splits against the corrected tree. Only identical splits have matching
support percentages; absent splits remain unassigned. This separate ML inference
does not replace the filtered tree or the topology passed to dating.
Bootstrap support, temporal signal and node-date uncertainty are different results.

The [IQ-TREE reference](https://github.com/iqtree/iqtree.github.io/blob/master/doc/Command-Reference.md)
specifies at least 1,000 UFBoot replicates and describes the additional bootstrap
NNI refinement. Inspect the saved `.iqtree` report for inference warnings.

## Temporal gate

The default runs 100 date permutations with the least-squares root refitted for the
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

With 100 successful permutations the minimum attainable empirical p-value is
1/101. `--randomisations 20` remains available for an explicitly smaller quick
screen, with a minimum of 1/21 and much coarser resolution. Both requested and
successful counts remain in the saved result; failed fits cannot improve support.

`--no-assess-temporal` skips clock/permutation tools. Fewer than three usable
distinct dates writes a complete tree product with an unsupported assessment.
Tree never creates a time-scaled phylogeny; it saves `assessment.json` for time.

### Relative-distance coherence policy

Native `chronoclade tree` defaults to `--coherence-policy strict`. Raw-alignment distance flags are audited and retained for review: recombination can inflate those distances. Strict policy stops only when extreme relative-distance flags remain after recombination correction. These flags can reflect legitimate cohort diversity; they do not establish incorrect genome identity. For an intentionally broad selected cohort, explicitly use `--coherence-policy diagnostic` to retain every selected genome and continue with a visible report warning. Neither policy prunes samples.

Both screens retain TSV and JSON audits, including thresholds, cohort median, median absolute deviation, every sample statistic and flagged statistics. Completed `tree.json` embeds both audits and the policy. The time report carries forward any coherence warning from its source tree. The policy is bound into single-tree resume fingerprints and tree-ensemble parameters. Legacy `run` retains the strict default.

Native tree `--seed` controls both IQ-TREE inference and temporal permutations. The manifest records `phylogeny_seed` and `phylogeny_seed_policy: explicit_IQ_TREE`; these are bound into the resume fingerprint so older unseeded inference is not reused as seeded output. Legacy core callers without a supplied seed retain IQ-TREE default seed behavior.

Native assembly trees select the highest-N50 reference from the input genomes, with deterministic sample-ID tie breaking. The shared inputs remain present in every alternative selection, so varying public context cannot change the mapping reference or its sequence hash. This fixes the reference basis across ensemble alternatives, although callable columns can still differ with cohort composition. The recorded `reference_selection_policy: highest_N50_common_input` is bound into both tree and ensemble fingerprints. A selection without input genomes fails. Legacy lineage reference selection remains unchanged.

## Progress and failure status

Normal CLI execution prints selected genome counts, the native tool currently
running, its log path and elapsed completion time. Long native commands report
that they are still running every 30 seconds. `chronoclade --quiet tree ...` hides
these progress messages while retaining the final result and errors. No command
arguments or credentials are copied into progress messages.

A tree or time ensemble with any computationally failed runs exits 2, retaining
its complete per-run manifest. A completed scientific result that does not support
dating exits 0 and remains explicitly unsupported. Automation should check both
the shell exit status and the saved per-run outcomes.
