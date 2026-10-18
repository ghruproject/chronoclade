# Independent analysis modules

ChronoClade provides four independent modules: **prepare**, **cgmlst**, **tree** and
**time**. The cgMLST module uses allele profiles by default. You can choose direct
DNA comparisons over the same cgMLST genes when an allele sequence catalogue is
available. Both modes produce saved selections that pass to **tree**, followed
by **time**. Reports retain the established style and identify the input method
and distance units.

## Prepare inputs

```bash
pixi run -e prepare chronoclade prepare records.json --input-kind profiles \
  --catalogues catalogues.json --out prepared

pixi run -e prepare chronoclade prepare \
  https://next.pathogen.watch/collections/jX5cwsoyJ1KDquMqssUzAD-ellington-et-al-2015 \
  --input-kind collection --out collection-input
```

The collection becomes the input dataset. Preparation does not discover public
context. It records profile and lineage readiness, exact identity links, missing
metadata and provenance. Exact linked ENA metadata can fill missing host, country,
collection date and isolation source; conflicting values remain inspectable.

Assemblies require an explicit species and configured typing tools **and reference
databases**. Installed callers do not make an unavailable database ready. Existing
valid profiles and assignments can be imported while reference access is pending.
Credentials are read from runtime configuration; do not put them in input files.

## Explore the complete compatible profile block

```bash
pixi run -e cgmlst chronoclade cgmlst prepared/dataset.json --out cgmlst \
  --lin-level 5 --subsamples 5 --context-size 50
```

The default distance is categorical allele mismatch. To compare DNA sequences
for those called cgMLST alleles instead, provide a compatible sequence catalogue:

```bash
pixi run -e cgmlst chronoclade cgmlst prepared/dataset.json --out cgmlst-dna \
  --distance-method direct-dna --allele-catalogue catalogue/catalogue.json \
  --distance-threshold 0.00001
```

This mode aligns observed alleles separately for each locus with MAFFT, reconstructs
each genome's aligned loci and concatenates them. It counts differences only at
jointly A/C/G/T columns and divides by their count. Missing calls/sequences become
`N`; aligned indels remain gaps and are reported separately from SNPs. Length-changing
alleles are retained in the alignment. Catalogue sequences with non-ACGT bases are rejected;
the minimum overlap and coverage audit remain explicit. Nearest relatives rank
by raw SNP count, with ties retained. These are SNP comparisons within the
cgMLST genes, not genome-wide or recombination-corrected SNP distances.
A direct-DNA run requires an explicit distance cutoff in nucleotide mismatch
fraction; the example is an exploratory setting, not a validated cluster boundary.

The cgMLST report uses the chosen distance consistently and adjusts the labels
and evidence tables to match the input. Allele mode reports raw allele
mismatches, shared called loci and mismatch fractions. Direct DNA mode reports
raw SNP counts, jointly compared nucleotide sites and normalized DNA distances.

Use `--fetch-context` to discover and retrieve compatible public profiles and
metadata first. The initial pool uses species/ST, followed by the explicitly
requested LIN depth (explicitly 5–10) or HierCC level.
An optional `--context-lin-level` restricts the union of input-associated public
prefixes independently of analysis block depth; `--display-lin-level` controls
report group colours. Neither defines a universal close-relative cutoff. Context assemblies are not
requested at this stage. The retrieval audit distinguishes discovered,
retrieved and usable records; an observed export is not a complete reference
locus catalogue.

The partition index links each disjoint block report and its selections.
Each block report states its share of the original inputs, leads with the saved
SNP-tree selection counts, and includes one-row-per-input closest-relative and
genetic-group summaries. The group table uses the same names as the ordination
legend and reports input/public counts, countries, hosts and collection years.
Detailed comparisons and selection evidence remain expandable.

Trees and ordination use the complete comparable profile block. Ordination defaults to
classical MDS (PCoA) of NJ patristic distances, summed along the branches between
tips. Negative branches are set to zero for this distance calculation only; the
saved NJ tree remains unchanged, and the count and total magnitude are audited.
A coordinate selector retains direct profile-distance PCoA for comparison.
The two-dimensional projection is approximate in both views.

Colour controls show dataset, country, region, host, cgLIN group or collection
date without changing coordinates within the chosen distance basis. Triangles
mark inputs and circles mark public comparisons. Hover, focus or select a point
to inspect accession and metadata, including overlapping genomes. Click a legend
entry to highlight its members. For cgLIN organisms, displayed groups use the
reported level-7 prefixes within each analysis block, using the same typing scheme.
Missing codes remain unassigned. Changing coordinates does not change membership.

Separate complete-linkage distance-cutoff diagnostics remain in the expandable
methods section. Those exploratory cuts do not define the displayed cgLIN groups.
Nearest neighbours report raw allele mismatches
or raw SNP counts, the matching shared-locus or callable-site denominator, and
all ties.

The report shows the cgLIN grouping level beside the ordination. Full prefixes
identify the groups, using the prepared typing assignments. The complete
NJ tree uses the same group colours; selecting a group in either the tree or the
ordination highlights its genomes in both. Groups need not be single tree clades.
The same reported groups supply selection strata. The tree opens with a full
overview. Zoom to inputs, a selected group or a
sample accession, then read tip labels without changing the saved topology or
signed branch lengths. NJ and MDS show the 20 largest groups plus all
input-containing groups initially; search and a scrollable full list retain
access to every group. Tables default to 10 rows. Region/NUTS2 detail is shown
only for reported values. The date regression has its own Temporal section;
root-sensitivity comparisons remain in saved evidence, outside the report.

The run and block folders include READMEs and categorized file browsers. These
link to original evidence without moving or duplicating large matrices. See
[the output layout](outputs.md#independent-module-outputs).

Multiple selections are alternative representative sets from the same block.
Every input, pinned context and required closest relative is retained, including
ties that exceed the requested context budget. Duplicate alternatives are labelled
when the pool offers insufficient choice.

Each selection also carries a checksummed coverage CSV for every available
genome in its block. It reports the nearest selected representative using the
chosen method's raw distance, retains shared-locus or callable-site denominators
and all boundary ties, and marks genomes with no comparable representative. The manifest
summarises the full-block nearest-distance distribution and records the exact
scope by sample count and ID checksum.

## Build a selected assembly tree and date it

```bash
pixi run -e tree chronoclade tree \
  cgmlst/BLOCK/selection_bundle/selections/selection-001.json --out tree-001

pixi run -e time chronoclade time tree-001/tree.json --out time-001
```

Tree acquires assemblies only for the exact selected IDs. It runs alignment,
phylogeny and recombination adjustment, then saves temporal evidence and its own
report. Failed acquisition does not choose replacements. Time reads the saved
tree, alignment/site-count and assessment evidence without rebuilding the tree.
Unsupported temporal evidence produces an explicit result rather than a dated
tree. Any supported override is explicit and recorded.

Repeat with the other saved selections, then compare:

```bash
pixi run -e time chronoclade time-compare \
  cgmlst/BLOCK/selection_bundle/ensemble.json \
  --run time-001/time.json --run time-002/time.json --out comparison
```

The comparison uses exact shared anchors and target pairs. It distinguishes the
shared-anchor ancestor date from whole-selection roots, which may describe
different ancestors. Between-selection ranges measure sampling sensitivity;
they are not confidence intervals or proof of convergence.

The earlier sequence-embedding work is parked in the
[research notes](research/parked-sequence-embeddings.md), with its evidence and
revisit criteria.

## Local and SLURM jobs

The same stage command can be saved as a local or SLURM job:

```bash
chronoclade job --slurm --account YOUR_ACCOUNT --partition YOUR_PARTITION \
  --cpus 8 --memory-gb 32 --out tree-job tree -- \
  cgmlst/BLOCK/selection_bundle/selections/selection-001.json --out tree-001 --threads 8
```

`job` writes a quoted script and a command/resource manifest. It does not submit.
Set the executable to an installed stage environment's `chronoclade` path when
necessary. Transfer the dataset, selections and their referenced files together,
and configure site resources and credentials independently. Script generation
and local testing do not establish successful execution on a real SLURM cluster.

## Execute alternative selections together

```bash
pixi run -e tree chronoclade tree \
  cgmlst/BLOCK/selection_bundle/ensemble.json --out tree-ensemble --ensemble-jobs 1

pixi run -e time chronoclade time \
  tree-ensemble/tree_ensemble.json --out time-ensemble --ensemble-jobs 1
```

Each selection has an independent run, fingerprint and saved status. A failed
selection does not erase successful alternatives. Increase `--ensemble-jobs`
only when resources support multiple native runs; each run also uses its own
thread setting. Tree ensembles share a validated assembly cache.

Time automatically writes the comparison report. Its shared-anchor forest shows
within-run dating regions separately from the between-selection range, with
actual descendant membership, coverage and descriptive topology agreement.
These alternatives measure sensitivity to selection, not statistical bootstrap
support or proof of convergence.

For already completed runs, repeat `--reuse-run PATH/time.json` on the time
ensemble command. Each result must pass artifact checks and match its selection,
source tree and dating-gate parameters. Reuse does not rerun dating; changed
inputs or parameters reject an existing output rather than mixing analyses.
