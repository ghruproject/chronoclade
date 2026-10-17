# Reading the reports

Each module has its own report: preparation, cgMLST, corrected SNP tree, and time
tree. Start with the module's saved manifest or report index. The cgMLST index
lists the input partitions and links to each block; **All groups** returns from a
block report to that summary.

Genome labels prefer a sequencing-run accession (SRR, ERR or DRR), then a
BioSample accession, then an assembly accession when available. User-supplied
input names are retained. Internal IDs remain in evidence and label mappings so
separate records stay traceable. See `sample_labels.csv` in a cgMLST block.

## cgMLST report order

### Dataset

The opening section summarizes input and public comparison genomes together,
before assembly subsampling. Country totals distinguish those roles. Region and
NUTS2 detail shows only reported values and disappears when there is no data;
unknown countries remain in the country denominator. Host/date coverage shows
what is available for later views.

Counts describe the retrieved, usable evidence. They do not measure population
prevalence or surveillance completeness. Country metadata does not establish
where infection was acquired.

### Genetic relationships

The nearest-relative table ranks comparable public genomes by raw allele
mismatches, or raw gene-panel SNP counts in direct-DNA mode. It retains tied
neighbours and the jointly compared loci or nucleotide sites.

The complete NJ tree uses every usable profile in its comparison cohort. It opens
as an overview. Use **Zoom to inputs**, **Zoom to selected group**, **Read labels**,
or find an accession and **Locate** it. Zoom and pan affect the viewport only.
Branches describe genetic distance; they do not represent calendar time.

MDS defaults to distances along that NJ tree. Negative branches are set to zero
for the ordination calculation only; the saved NJ tree keeps its signed branch
lengths. Direct profile-distance PCoA remains available through the coordinate
selector. Both are approximate two-dimensional projections, so points can overlap.
Hover, focus or select a point to inspect the accessions and metadata.

For Klebsiella blocks, displayed groups use compatible reported cgLIN level-7
prefixes. They are not inferred from the ordination or a 2% distance cutoff.
Group colours match between NJ and MDS; selecting a group highlights it in both.
The selectors initially show the 20 largest groups plus every input-containing
group. Search covers all groups, and **Show all groups** opens a bounded list.
Every group and genome remains available. Groups need not form a single tree clade.

Colour by dataset, country, region, host or date without changing the coordinates
within a chosen basis. Triangles mark inputs and circles mark public comparisons.

### Temporal

The exploratory date-versus-distance plot has its own section. Collection-date
intervals use their midpoints; pre-1900 dates are excluded from this regression
and retained elsewhere. A scaled allele-distance slope depends on locus coverage
and is not a direct mutation-rate estimate. Root-to-tip points share ancestry;
a fitted line alone does not establish temporal signal.

Root-comparison tables, sensitivity figures and their methodological discussion
are omitted from the report. The saved diagnostic evidence remains available.
Formal temporal assessment belongs to the corrected SNP-tree stage.

### State network

Choose location or host where available. The interactive view reconstructs
metadata state changes along the rooted tree. Dashed/alternative links expose
uncertainty across optimal assignments. Select nodes, change the centrality
measure, drag, zoom or pan, then save PNG or interactive HTML from the view.
The report does not repeat the network as a second static figure.

Arrows describe reconstructed ancestral-to-descendant state changes. They depend
on the tree, root, metadata and sampling frame, and do not prove transmission.

### Subsamples

This section states how many inputs and public comparisons enter each saved
assembly selection, how many neighbours are mandatory, and how much alternatives
overlap. All inputs, pinned records and required neighbours are retained. The
context budget can expand to preserve required ties and an explicitly requested
representative reserve. Optional representatives vary across reproducible seeds.

Download a selection or the ensemble manifest to run the separate `tree` module.
The full cgMLST NJ/MDS views remain based on the complete comparable profile pool.
Variation between selections measures sampling sensitivity; it is not a confidence
interval or proof of convergence.

### Methods and files

Tables default to 10 rows and support search, column sorting and page-size
selection. Detailed coverage and optional distance-cutoff diagnostics remain
expandable. Those optional cuts do not define the displayed cgLIN groups.

**Browse output files by type** opens categorized links; **Output folder guide**
opens the block README. Report figures and interactive views are embedded in the
HTML. Keep neighboring files for evidence downloads and stage handoffs. See the
[output layout](outputs.md#independent-module-outputs).

## SNP-tree and time reports

The tree module reports the selected assembly set, recombination adjustment,
corrected distances, relationships, state-change network and temporal assessment.
Raw and corrected lineage-coherence screens detect grossly divergent genomes;
they do not verify every metadata value. Review flagged accessions and lineage
assignments before overriding a strict screen.

Time reads the saved corrected tree and assessment. Supported dating produces a
calendar tree and model intervals. Unsupported evidence remains explicit. A
successful execution does not establish precise ancestral dates.

The ensemble comparison uses exact shared anchors. It distinguishes within-run
model uncertainty from variation between alternative selections; whole-selection
roots can describe different ancestors. See [tree](tree-stage.md) and
[time](time-stage.md) for the recorded checks and comparison policies.

## Older reports

The older `run --mode fast|full|finish` interface uses different stage paths and
may retain older report sections. See [the staged workflow](staged-workflow.md).
Preserved PhiPack results belong to the earlier assembly-screen workflow, not the
current independent cgMLST module.
