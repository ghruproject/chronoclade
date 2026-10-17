# Independent cgMLST analysis

`cgmlst` reads a validated prepared bundle. It partitions inputs at an explicit
cgLIN depth (5 by default; 6 or 7 are options), includes every matching context
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
LIN5/6/7 counts. Counts are scoped to the retrieved/prepared evidence; they do not
imply complete worldwide surveillance. Missing ST or lineage evidence is reported;
a neighbour's label is never copied to an input genome.

Every comparable profile remains in the NJ/ordination analysis and downloadable
tree. Selection does not prune this guide tree. The approved location network and
report presentation are reused. Ordination controls show dataset, country, host,
date and genetic-group colours on the same coordinates. Displayed genetic groups use compatible reported cgLIN prefixes. Optional distance-cut
diagnostics remain exploratory and do not replace lineage assignments.

Blocks above 1,500 records use [compiled distances and RapidNJ](cgmlst-scale.md),
disk-backed checked distance evidence, leading-axis ordination and the full
interactive NJ tree. Default bootstrap counts are 30 for dense blocks and zero for large blocks,
recorded separately for each block. Large-block bootstrapping is not implemented;
requesting it fails explicitly. A 10,000 × 629 synthetic full analysis was tested.

## Report navigation

Open `index.html` for the block summary, then follow a block report. Each block
has a **back to all groups** link. Results are ordered as:

1. Combined input/public dataset: country totals, available region/NUTS2 detail,
   host and date coverage.
2. Genetic relationships: nearest relatives, the complete NJ tree, NJ-distance
   MDS and reported cgLIN group summaries.
3. Temporal: exploratory distance versus collection date.
4. State network: interactive location or host reconstruction, with export controls.
5. Subsamples: selected counts, mandatory neighbours, alternative-set overlap
   and handoff manifests.
6. Methods and evidence files: coverage, optional cutoff diagnostics and downloads.

Searchable, sortable tables open with 10 rows. Region/NUTS2 detail includes only
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

Equal-length homologous allele sequences are compared sitewise. Missing calls,
missing sequences and unequal-length loci are excluded and recorded in
`sequence_comparison_audit.json`. Minimum overlap applies to compared loci.
Distance is raw nucleotide differences divided by jointly compared nucleotide
sites; the nearest-relative ranking uses raw differences and retains ties.
These are gene-panel SNP comparisons, not genome-wide or recombination-adjusted
SNP distances. Positional comparison assumes the same locus coordinates;
compensating indels or misaligned alleles require separate alignment work.

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
exports are not reused.

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
