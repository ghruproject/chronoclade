# Independent cgMLST analysis

`cgmlst` reads a validated prepared bundle. It partitions inputs at an explicit
cgLIN depth (5 by default; explicit depths 5–10 are supported), includes every matching context
record already in that bundle, and writes one existing-style report per block.
Species, supplied MLST ST and lineage scheme/version/database scope separate
blocks. Declared MLST schemes remain separate. A missing MLST scheme can enter
the sole declared comparison namespace for its species/ST while its metadata
remains unknown; multiple declared schemes make that missing namespace ambiguous
and excluded. It does not automatically narrow to a deeper prefix.

```bash
chronoclade prepare frozen-records.json --catalogues locus-catalogues.json --out prepared
chronoclade cgmlst prepared/dataset.json --out cgmlst --lin-level 5 \
  --subsamples 5 --context-size 50
```

For E. coli with existing compatible HierCC assignments, supply an explicit level:

```bash
chronoclade cgmlst prepared/dataset.json --out ecoli-cgmlst --hiercc-level HC10
```

Use `--fetch-context` to discover the public same-ST candidate pool and retrieve
compatible profiles and lineage metadata matching the requested LIN/HierCC block.
Context assemblies are acquired only by the later tree stage. Supplied public
profile/lineage exports are also supported. The partition index and `partitions.json`
record exact included IDs, unmatched context, unresolved inputs and available
LIN5–10 counts. Counts are scoped to the retrieved/prepared evidence; they do not
imply complete worldwide surveillance. Missing ST or lineage evidence is reported;
a neighbour's label is never copied to an input genome.

### Narrower context with independent display groups

For a documented finer scope, use an explicit prefix depth; no default is chosen
from a plot or a regression. For example:

```bash
chronoclade cgmlst prepared/dataset.json --fetch-context --context-lin-level 8 \
  --lin-level 5 --display-lin-level 8 --backend scale \
  --subsamples 5 --context-size 50 --bootstrap-replicates 0 --out cgmlst-lin8
```

These settings have different purposes:

| Setting | Purpose | Effect |
| --- | --- | --- |
| `--context-lin-level 8` | Public context scope | Retrieve/retain the union of input-associated LIN8 prefixes, including matches in another input’s prefix |
| `--lin-level 5` | Disjoint analysis blocks | Analyse that union within compatible LIN5 blocks, allowing cross-LIN8 nearest relatives |
| `--display-lin-level 8` | Reported group colours in NJ and MDS | Show reported LIN8 prefixes; coordinates and distances are unchanged |
| Nearest relatives | Distance-ranked comparisons within each block | Use actual allele differences, jointly called loci and ties; neither LIN depth is a distance cutoff |

Display depth defaults to 7 for compatibility. Representative-selection genetic
strata retain the existing reported LIN7 basis independently of display depth;
this is a reproducible selection convention, not a validated close-relative
threshold. A LIN8 block may therefore contain only one LIN7 selection stratum;
time, geography, distance and mandatory-neighbour rules still apply.

The ST258 LIN8 pilot motivated these explicit options, but LIN8 is not a universal
outbreak or close-relative definition. Narrowing context can exclude a nearest
relative outside every input-associated prefix. Neighbours are computed **within
each analysis block**. The example retains cross-LIN8 comparisons inside LIN5;
using `--lin-level 8` instead deliberately separates those prefixes. Omitting
`--context-lin-level` preserves the original behavior: public retrieval follows
`--lin-level`, and frozen contexts are restricted only by the analysis blocks.
The optional context scope filters both newly retrieved and already prepared
context records; inputs are retained and unresolved typing is audited. It cannot
be combined with HierCC. The index and `partitions.json` retain original prepared
pool counts, scope counts and excluded context IDs. Public
retrieval audits additionally distinguish same-ST discoveries, requested prefix
matches, returned profiles, compatibility exclusions and usable sample units;
these are accessible-provider counts rather than worldwide prevalence.

Every comparable profile remains in the NJ/ordination analysis and downloadable
tree. Selection does not prune this guide tree. The approved location network and
report presentation are reused. Ordination controls show dataset, country, host,
date and genetic-group colours on the same coordinates. Displayed genetic groups use compatible reported cgLIN prefixes. Optional distance-cut
diagnostics remain exploratory and do not replace lineage assignments.

`--backend auto` preserves the size-based backend choice. `--backend scale` uses
compiled allele distances, compiled complete-linkage and RapidNJ explicitly, even
for smaller blocks; it is useful for several hundred closely related profiles
where the dense Python clustering/NJ route is slow. `--backend dense` explicitly
uses the dense route and rejects blocks above 1,500 records. The requested and
actual backend are recorded in `cgmlst.json` and bound to resume parameters.
Scale supports allele profiles and requires RapidNJ (on macOS, build/install it
with `pixi run python scripts/install_rapidnj.py` as described in
[large profile analysis](cgmlst-scale.md)). Direct DNA uses auto/dense. The scale
backend reports its own clustering tie policy and ordination solver; choosing a
backend is explicit rather than a silent claim of identical numerical methods.
Use `--bootstrap-replicates 0` for the scale route; large/scale bootstrap support
is not implemented.

Blocks above 1,500 records use [compiled distances and RapidNJ](cgmlst-scale.md),
disk-backed checked distance evidence, leading-axis ordination and the full
interactive NJ tree. Default bootstrap counts are 30 for dense blocks and zero for large blocks,
recorded separately for each block. Large-block bootstrapping is not implemented;
requesting it fails explicitly. A 10,000 × 629 synthetic full analysis was tested.

## Report navigation

Open `index.html` for the block summary, then follow a block report. Each block
has a **back to all groups** link. Results are ordered as:

1. Combined input/public dataset: country totals, available region/city detail,
   host and date coverage.
2. Genetic relationships: nearest relatives, the complete NJ tree, NJ-distance
   MDS and reported cgLIN group summaries.
3. Temporal: exploratory distance versus collection date.
4. State network: interactive location or host reconstruction, with export controls.
5. Subsamples: selected counts, mandatory neighbours, alternative-set overlap
   and handoff manifests.
6. Methods and evidence files: coverage, optional cutoff diagnostics and downloads.

Searchable, sortable tables open with 10 rows. Region/city detail includes only
reported values and is omitted when neither field has data. Unknown countries
remain in country totals.

The NJ viewer opens with the full tree. Use **Zoom to inputs**, select a group
and **Zoom to selected group**, or find an accession and **Locate** it. **Read
labels**, vertical zoom, drag/pan and **Fit all** support closer inspection.
These controls change the viewport, not the saved topology or tip membership.

NJ and MDS selectors initially show the 20 largest groups plus every group
containing input genomes. Search covers all groups; **Show all groups** retains
a bounded scrollable list. Selecting a group highlights it in both figures.
No rare groups or ordination points are removed by the compact legend. The
MDS coordinate selector offers direct profile-distance PCoA for comparison;
colouring by country, host or date keeps the chosen coordinates fixed.

Root and block `README.md` files explain the handoffs. `files/index.html` groups
links into tables, trees, figures, networks, subsamples and methods. Evidence
stays at its recorded paths; the navigation folders do not duplicate matrices.
Reports embed figures and interactive views for opening without a server. Keep
the folder for evidence downloads. See [output files](outputs.md).

## Distance modes

Categorical allele-profile mismatch is the default. An allele number is a label,
not a numeric measurement: two different called alleles count as one locus
mismatch. Missing calls never count as matches.

For direct nucleotide differences within the cgMLST genes:

```bash
chronoclade cgmlst prepared/dataset.json --out cgmlst-dna \
  --distance-method direct-dna --allele-catalogue catalogue/catalogue.json \
  --distance-threshold 0.00001
```

Supply a complete `chronoclade.allele-sequence-catalogue` version 1 JSON manifest
referencing a checksummed CSV with exactly `locus,allele,dna` columns. Its species,
scheme ID/version and database version/fingerprint must match the prepared bundle.
The command validates the catalogue locally; it does not fetch assemblies or
infer sequences from allele numbers. Ambiguous DNA bases and duplicate catalogue
entries are rejected.

For each locus, observed allele sequences are aligned with MAFFT, including
equal-length variable alleles. Each genome's aligned loci are then reconstructed
and concatenated in canonical scheme order. Missing calls or missing sequences
become `N`; alignment gaps remain `-`. Entire loci are not dropped merely because
their alleles differ in length. Globally absent loci retain explicit zero-width
partition records.

SNP distance counts differences only at jointly A/C/G/T columns, divided by the
number of those columns. Gap-versus-base columns and indel tracts are reported
separately; zero SNPs does not establish sequence identity. Minimum overlap
applies to loci with jointly callable nucleotide sites. Nearest relatives rank
by raw SNP differences and retain ties. These are aligned gene-panel comparisons,
not genome-wide or recombination-adjusted SNP distances. Alignment uncertainty
remains a limitation.

The Pixi `cgmlst` and default environments include MAFFT. Each block's
`locus_alignment/` directory preserves the concatenated FASTA, locus partitions,
per-locus alignments and tool provenance; `sequence_comparison_audit.json` links
this evidence and records missing calls/sequences.

Direct DNA requires an explicit cutoff in nucleotide mismatch fraction. The
example cutoff is descriptive, not a validated epidemiological boundary. Allele
mode retains its default cutoff of 0.02. Reports reuse the same style and views,
with mode-specific labels, counts, coverage and temporal units. Bootstrap
replicates resample whole loci in both modes. Direct DNA currently supports up to
1,500 records per block; larger blocks fail explicitly and remain supported by
the allele-profile backend.

## Representative alternatives

`--subsamples N` requests alternative representative sets from each block.
All input genomes, explicit pins (`--include ID`) and declared nearest neighbours,
including ties at the requested boundary, are mandatory. The context-size budget
excludes input genomes. Mandatory contexts can exceed it: `budget_overrun`
records the excess rather than dropping required samples. To retain optional representatives
even when nearest-neighbour ties exhaust the budget, use `--representative-reserve 50`.
The effective context budget is the larger of `--context-size` and the mandatory
context count plus this reserve; the available pool can limit the result. This is
explicitly recorded in each selection rather than silently dropping ties.

Nonmandatory representatives vary across reproducible seeds while retaining the
selection strata and quota policy. Every selection records IDs, per-ID reasons,
coverage, missing metadata, source distance units and checksums. Pairwise overlap
and cumulative unique coverage are retained in the ensemble. A block with too few
choices can produce repeated sets; those are explicitly labelled as duplicates,
not presented as independent evidence of convergence.

Each block writes `selection_bundle/ensemble.json` and
`selection_bundle/selections/selection-001.json` (and subsequent alternatives).
The bundle copies its distance evidence and refers to its prepared dataset and
partition manifest using relative paths and hashes. Move the common parent of
these bundles together. The public loaders reject changed sources or members
that do not belong to the exact source partition.

These manifests feed the independent [tree](tree-stage.md) and [time](time-stage.md) commands directly. Ensemble comparison uses common input anchors and records missing or unsupported runs. Between-selection variation measures subsampling sensitivity, not a formal population confidence interval.

## Output ownership

`cgmlst.json` is published after all requested blocks finish. It references each
report and ensemble, the prepared dataset, partition evidence and exact analysis
settings. Existing nonempty output directories are rejected unless `--resume` is supplied.
A failed stage does not publish a completed manifest. Repeat the same command with
`--resume` to reuse hash-validated retrieval and numerical checkpoints. Changed
source data or analysis settings are rejected. A revised analysis needs a new output directory. The shared `job` command can execute the saved command locally or generate a SLURM script; it never submits automatically. Tree/time stages provide their own checked resume.


## Command-driven study workflow

Descriptive study metadata can be supplied during live preparation or added to an
existing frozen bundle with exact accession/sample IDs. Typing and membership are
preserved; changed descriptive values retain their previous value as a conflict
and their user-supplied provenance. Hospital, NUTS2 and specimen fields are retained.

```bash
chronoclade prepare accessions.csv --species "Klebsiella pneumoniae" \
  --metadata study-metadata.csv --out prepared
chronoclade cgmlst prepared/dataset.json --fetch-context --lin-level 5 \
  --subsamples 5 --context-size 50 --representative-reserve 50 --out cgmlst
# If interrupted, repeat the cgmlst command with --resume.
```

Public context is deduplicated after compatible profile retrieval using exact
BioSample/run/assembly-linked sample units. Every input remains present; duplicate
source IDs and representative decisions remain in the preparation audit. No
genetic-distance collapsing is performed. Public catalogue snapshots and completed
scope exports are retained in the sibling `_context_retrieval` directory; incomplete
exports are not reused. Accepted raw cgMLST/HClink batches have checked request
receipts and response hashes. With the same source and settings, `--resume`
can reuse those accepted downloads after interruption during native assignment;
the unfinished assignment is run again. This does not promise that every provider
request is eliminated on retry.

Native-reference reuse checks the actual caller/assignment resource contents,
not just the configuration path. Changing a database or scheme rejects an old
numerical resume. Prepared bundles with only a configuration fingerprint need an explicit
`--typing-config` and a new output directory for fresh analysis.

For an independently validated captured whole-ST raw source, use
`--context-source path/to/provider-source.json` instead of `--fetch-context`.
This currently requires a configured Klebsiella native reference. ChronoClade
checks the saved catalogue, export identities, receipts and file hashes, then
reassigns the profiles with that reference. The report records the captured
source date; replay does not imply a fresh public search or reuse old normalized
assignments.

Selection evidence uses independent filesystem copy-on-write snapshots where
supported, with a regular-copy fallback. It never uses mutable hardlinks. Retaining
many completed runs still consumes space; copy-on-write reduces duplicate matrix
storage within a run. Exact outgroup discovery is not automated by these commands;
the declared root remains recorded in the saved analysis evidence.

## Root-to-tip diagnostics

The Temporal section shows the exploratory NJ distance-versus-date regression.
The root-comparison table, methodological discussion and paired sensitivity
figure are not displayed in the report. Saved diagnostic evidence retains a
date-optimised root, found by searching continuous positions along tree edges
to minimise dated-tip regression residual mean squared error, using the same
principle as Clockor2. This is an exploratory allele-divergence diagnostic, not a
calibrated molecular clock.
Collection-date interval midpoints are used; dates before 1900 are excluded from
the regression and retained in the metadata and genetic analysis.

Negative NJ branches are set to zero only in a separate diagnostic tree. Both
root comparisons use that tree, with the number and magnitude of changed edges
saved in the diagnostic evidence. A signed-original baseline is retained as a sensitivity check. The NJ
tree displayed elsewhere and the location network keep their original branches
and root. Normalised cgMLST distances are scaled by the recorded locus count;
with variable missing-locus coverage these are allele-distance units rather than
literal accumulated allele changes. A date-optimised fit alone does not establish
temporal signal.

Refresh a completed report from its checked saved NJ and coordinate evidence,
without rebuilding the tree, ordination or selections. This updates report
controls and output-folder navigation as well as saved clock diagnostics. It
does not restore the removed root-comparison section. Alternatively, compare the
exact selections used by a native SNP-tree ensemble:

```bash
chronoclade refresh-cgmlst-report cgmlst/cgmlst.json
chronoclade compare-clocks tree/tree_ensemble.json --out matched-clock-comparison
```

The matched comparison constructs a fresh NJ tree from each exact selected
cgMLST distance submatrix. Both panels use the same date-eligible genomes and the
same continuous-edge residual-MSE rooting criterion. The report separates
allele-distance units/year from SNP-tree substitutions/site/year; their slopes
are not directly interchangeable. It compares descriptive fits and does not
perform another date-randomisation test or recalibrate the saved time trees.

## Experimental finer-group calibration

```bash
chronoclade calibrate-lin cgmlst/<block>/selection_bundle/ensemble.json \
  --folds 3 --output lin-calibration
```

This trial analyses the complete frozen profile pool, not the selected assembly
subsample. It builds fresh NJ trees for input-containing LIN6 and LIN7 groups;
singleton or insufficiently dated groups remain in the table without a fitted
rate. These nested groups are not independent replications.

For each of the latest calendar-year blocks, rooting and regression use only
earlier eligible dates. The held-out dates then predict root-to-tip allele
distance using that frozen root and fitted line. Mean absolute distance error is
compared with two predictions at the same root: the training-set mean and the
latest training-year mean. Held-out genomes still contribute genetic information
to the NJ tree, so this is a transductive date check, not validation on unseen
genomes. Fold counts and rate/root sensitivity are saved; tip independence and
formal confidence intervals are not assumed.

The trial does not publish biological cutoffs. `2*r*H` is saved only as a
conditional expected pair-distance scenario when both descendants are H years
from a shared ancestor under a constant rate. It does not follow from the
difference between their sampling dates, and it is neither an empirical cutoff
nor an uncertainty bound. Stronger predictive evidence and independently
labelled epidemiological comparisons would be needed to justify a cutoff.

## Supplied profiles without lineage assignments

For a relationship analysis of supplied inputs that have cgMLST profiles but no
verified LIN/HierCC assignment, explicitly select:

```bash
chronoclade cgmlst prepared/dataset.json --partition-by st --out input_profiles
```

This partitions by species, MLST namespace/ST and compatible profile scope.
It accepts input records only and rejects prepared public context, `--fetch-context`
and lineage/context-scope options. It does not assign lineage codes, retrieve
context or establish a close-relative cutoff. Missing lineage evidence stays
explicit. All supplied inputs remain required in the downstream selection; this
option does not subsample away the input set. The normal lineage-based workflow
remains the default.

The partition mode, exact memberships, exclusions and source bundle are recorded
in the analysis/partition manifests and selection provenance. Resume rejects a
changed mode. This input-only route is useful while a separate HierCC reference
is unavailable, but is not evidence that native HierCC assignment is supported.

## Reading the comparison scope

The stage summary and group report distinguish the saved search scope, the
input-associated context prefixes, the analysis block and the LIN level used
for colours. Source-record counts precede deduplication; public-comparison counts
refer to the retained genomes. Equal counts are combined to keep the summary
readable. A broad background block is not a close-relative cluster.

The detailed nearest-neighbour table adds one **LIN relationship** column for
cgLIN data. It says whether the reported prefixes match and how many consecutive
LIN levels they share. These are matching code components, not allele
differences. Distance rankings and ties remain based on the measured profile or
gene-SNP differences. Unresolved codes and unverified reference compatibility
remain explicit.

## Refresh a report

`chronoclade refresh-cgmlst-report cgmlst/cgmlst.json` repackages saved evidence
without refitting root/date diagnostics. To request those calculations explicitly,
add `--recompute-diagnostics`. Source trees, coordinates and selections are checked
and preserved in either case. Computation provenance remains attached to the
original results; rendering and any new diagnostics get separate provenance.
