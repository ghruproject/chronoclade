# Output files

## Independent module outputs

Each command owns a separate output directory. Follow the report paths recorded
in its manifest; a single selection and an ensemble produce different manifests.

| Stage | Entry points | Handoff |
| --- | --- | --- |
| `prepare` | `report.html`, `readiness.json`, `prepare.json` | `dataset.json` and its referenced files |
| `cgmlst` | `index.html`, `README.md`, `files/index.html`, `cgmlst.json` | Each block's `selection_bundle/ensemble.json` or individual selection JSON |
| `tree` | Report linked from `tree.json` or `tree_ensemble.json` | The same tree manifest plus referenced evidence |
| `time` | Report linked from `time.json` or `time_ensemble.json` | Dating evidence and, for an ensemble, its comparison report |

A cgMLST run has this layout (`BLOCK` represents a generated partition directory):

```text
cgmlst/
  index.html                   # Summary of all blocks
  README.md                    # Start here and stage handoffs
  cgmlst.json                  # Run parameters, identities and block references
  partitions.json              # Partition counts, membership and exclusions
  files/index.html             # Navigation to each block
  BLOCK/
    profile_report.html        # Self-contained cgMLST report
    README.md                  # Block guide
    profile_analysis.json      # Report data and diagnostics
    numerical_analysis.json    # Frozen numerical checkpoint and checksums
    sample_labels.csv          # Display labels and stable internal IDs
    cohort_1_nj.nwk            # Complete cohort NJ tree
    cohort_1_tree_groups.html   # Standalone interactive NJ viewer
    cohort_1_ordination.html    # Standalone interactive MDS/PCoA viewer
    binary_distances/           # Full, potentially large distance evidence
    selection_bundle/
      ensemble.json            # Handoff to tree for all alternatives
      selections/
        selection-001.json     # Handoff for one alternative
    files/
      index.html
      tables/index.html
      trees/index.html
      figures/index.html
      networks/index.html
      subsamples/index.html
      methods/index.html
```

Cohort names and optional artifacts depend on the available data. The categorized
`files/` folders contain navigation pages linking to original artifacts. They do
not relocate evidence or make extra matrix copies. Saved root-sensitivity
calculations remain evidence files; that comparison is omitted from the cgMLST
report itself.

Report images and interactive views are embedded, so opening the HTML does not
require a server. Links to downloads, block indexes and file browsers still need
the corresponding output files. Share the report alone for viewing, or keep the
common parent directory with prepared data, selections and stage outputs for
reproducibility and subsequent stages.

Distance matrices grow quadratically. Copy-on-write selection snapshots reduce
duplication where supported, but repeated completed runs can still use several
gigabytes. Avoid editing or removing manifest-referenced evidence from runs you
intend to resume or pass downstream. See [large-block evidence](cgmlst-scale.md).

## Older run interface outputs

The paths below describe `run --mode fast|full|finish` and the older
`prepare-context` interface. They differ from the independent module layout above.

The top-level `index.html` links the report for each species/lineage group.
`summary.json` records analysed and skipped groups.

The `logs/` directory contains each external command and a JSON fingerprint of
its direct inputs. ChronoClade uses those fingerprints to decide whether a
completed stage can be reused.

## Public staged run outputs

The public run defaults to profile-first `fast` mode. It writes a root
`index.html`, `fast/profile_report.html` and profile-analysis data such as
`fast/profile_analysis.json`, profile distances, nearest neighbours, groups,
exclusions and the profile PCoA/NJ figures. Per-lineage profile reports link
only to figure and data files that were actually produced. `full` adds a
corrected assembly-based report and supporting results; `finish` adds temporal
evidence and a time tree only when the clock gate passes. Earlier stage reports
remain available after later stages complete. See
[`staged-workflow.md`](staged-workflow.md) for stage-level paths and semantics.

## Main lineage outputs

| File | Content |
| --- | --- |
| `report.html` | Reader-facing analysis report |
| `report.json` | Machine-readable report record |
| `supporting_results.zip` | Portable report evidence bundle |
| `core_alignment.fasta` | SKA2 reference-ordered alignment |
| `lineage_coherence.tsv` | Per-genome raw-distance screen run before IQ-TREE and ClonalFrameML |
| `clonal_lineage_coherence.tsv` | Per-genome recombination-filtered distance screen run before temporal analysis |
| `iqtree.treefile` | Starting maximum-likelihood tree |
| `clonalframeml.labelled_tree.newick` | Recombination-corrected tree |
| `clonalframeml.importation_status.txt` | Recombination intervals by branch |
| `clonalframeml.filtered.fasta` | Alignment after recombination filtering |
| `recombination_intervals.tsv` | Normalised branch-level calls with alignment and reference-record coordinates, including cross-record flags |
| `recombination_genome_profile.csv` | Binned importation depth and filtering values used in the map |
| `recombination_map.svg`, `.png` | Genome-wide recombination and alignment-filtering figure |
| `recombination_summary.json` | Counts of inferred intervals and columns removed or retained |
| `clonal_pairwise_distances.tsv` | Pairwise SNPs, callable sites and comparison classes |
| `clonal_snp_matrix.tsv` | Square clonal SNP matrix |
| `pairwise_callable_sites.tsv` | Square callable-site matrix |
| `public_health_evidence.json` | Scenario rules, groups, neighbours and sensitivity results |
| `nearest_neighbours.tsv`, `.csv`, `.json` | Per-focal nearest comparisons, exact ties, callable sites, separate tree-distance ranks and scope/provenance |
| `genetic_tree.svg`, `.png` | Genetic relationship tree before calendar dating; large trees use labelled overview pages |

## Temporal evidence

| File | Content |
| --- | --- |
| `clock/root_to_tip_regression.svg` | TreeTime root-to-tip plot |
| `clock/rtt.csv` | Root-to-tip values for every genome |
| `root_to_tip.png` | High-resolution regression plot |
| `temporal_signal.json` | Observed and permuted clock statistics |
| `date_randomisation.csv` | One row per observed or permuted fit |
| `date_randomisation.svg`, `.png` | Randomisation distributions |
| `timetree/timetree.nexus` | Time-scaled tree, when the gate passes |
| `timetree_with_confidence.svg` | Calendar tree with node intervals |
| `timetree.png` | High-resolution dated tree |
| `timetree_confidence.csv` | Clock and root uncertainty summary |
| `node_dates.csv` | Internal-node estimates and 90% intervals |

## Legacy per-lineage fast-screen evidence

The files below describe the earlier assembly-based fast screen. They may
appear in preserved results from that workflow, but the current public `fast`
mode is profile-first and does not run PhiPack, SKA alignment or IQ-TREE.

| File | Content |
| --- | --- |
| `phipack_profile.tsv` | PhiPack Profile positions and p-values, with reference-record coordinates |
| `phipack_significant_blocks.tsv` | PHI-positive profile results and their computational block coordinates |
| `phipack_summary.json` | Parameters, tested blocks, detection result and interpretation limits |
| `core_alignment.fasta` | Uncorrected SKA alignment supplied to the fast-mode IQ-TREE run |
| `iqtree_fast.treefile` | Uncorrected screening phylogeny used for root-to-tip analysis |

## Context preparation

`prepare-context` writes the complete same-ST accession list, the balanced
candidate pool, SKA screening distances, the selected context manifest,
combined metadata and a JSON audit of the selection process. Keep these files
with the analysis. A different public context sample can change the apparent
placement of focal isolates.

The Pathogenwatch route also writes the full frozen raw catalogue
`pathogenwatch_catalogue.json`, annotated `context_catalogue.json`,
`download_ledger.json` (successes, failures, checksums, bytes and timing), and
`context_geography/` with count/percentage SVG and PNG figures, source tables,
assignment coverage, deduplication audit and standalone HTML. The public panel
uses the full QC-eligible catalogue; selected-context and focal panels have
separate denominators. These outputs remain usable when dating is unsupported.
Manifests retain source IDs, accessions, cgLIN prefix keys and catalogue hashes.
Keep the referenced frozen catalogue together with the manifest. Lineage reports
copy the annotated catalogue, selection audit and figure tables to supporting
results archives.
