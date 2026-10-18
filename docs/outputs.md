# Output files

## Independent module outputs

Each command owns a separate output directory. Follow the report paths recorded
in its manifest; a single selection and an ensemble produce different manifests.

| Stage | Entry points | Handoff |
| --- | --- | --- |
| `prepare` | `README.md`, `report.html`, `readiness.json`, `prepare.json` | `dataset.json` and its referenced files |
| `cgmlst` | `index.html`, `README.md`, `files/index.html`, `cgmlst.json` | Each block's `selection_bundle/ensemble.json` or individual selection JSON |
| `tree` | `README.md`; report linked from `tree.json` or `tree_ensemble.json` | The same tree manifest plus referenced evidence |
| `time` | `README.md`; report linked from `time.json` or `time_ensemble.json` | Dating evidence and, for an ensemble, its comparison report |

New prepare, tree and time outputs include a root `README.md` explaining the
saved result, evidence folders, next command and recovery policy. Tree/time
ensemble guides distinguish individual `runs/` from comparison outputs. A
completed preparation may still retain missing typing: inspect `readiness.json`
before proceeding. A successful command does not establish biological validity.
Guides written in older outputs are not added retroactively until the stage is
run again through a code path that writes them.

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

Profile geography summaries use reported country, region and city metadata.
Missing region or city values stay absent in the report; NUTS2 remains in the
machine-readable evidence for compatibility but is not used as a substitute
for city or shown as a geography level. Host is separate sample metadata, not a
geographic field.

The state-change network can reconstruct **country, host, region or city** on
the same saved tree. Use **Reconstruct state** to switch fields. Optional fields
are offered only when reported for that tree; a field with only one observed
value shows a coverage message instead of an empty graph. Missing values stay
unknown. Region and city labels include known parent locations to distinguish
places with the same name. Each field has its own saved reconstruction audit.
The interactive views retain uncertainty controls and PNG/HTML exports. Arrows
represent inferred state changes, not confirmed transmission.

Report images and interactive views are embedded, so opening the HTML does not
require a server. Links to downloads, block indexes and file browsers still need
the corresponding output files. Share the report alone for viewing, or keep the
common parent directory with prepared data, selections and stage outputs for
reproducibility and subsequent stages.

Distance matrices grow quadratically. Copy-on-write selection snapshots reduce
duplication where supported, but repeated completed runs can still use several
gigabytes. Avoid editing or removing manifest-referenced evidence from runs you
intend to resume or pass downstream. See [large-block evidence](cgmlst-scale.md).

## Viewer export files

Use [`chronoclade export`](cli.md#viewer-exports) to create local GrapeTree or
Microreact files from a saved module result. The command creates a new directory;
it does not upload files or modify the source result.

| Viewer | Input files | Audit and instructions |
| --- | --- | --- |
| GrapeTree | `profiles.tsv`, `metadata.tsv` | `allele_key.csv`, `README.txt`, `export.json` |
| Microreact | `tree.newick`, `metadata.csv` | `id_map.csv`, `README.txt`, `export.json` |

The viewer files share deterministic safe `CC_…` identifiers. Metadata retains
`sample_label` for readable accession/strain labels and `chronoclade_sample_id`
for the original internal identity. In Microreact, select `ID` as the shared key
and use `sample_label` for display labels; `id_map.csv` also records the mapping.
Do not relabel just the tree or just the metadata independently.

GrapeTree profiles have a `#Strain` header and one column per exported locus.
Alleles are recoded into positive numerical categories **within each locus**;
`allele_key.csv` maps categories back to source allele identifiers or hashes.
`0` means missing, not an allele or an observed match. Recoded numbers are
categorical labels, not numerical measurements. The export records the supplied
locus universe; exporting an incomplete catalogue does not make it complete.

Microreact exports retain the selected saved tree topology, internal support
labels and branch lengths, with metadata restricted to its tips. A dated result
is converted to Newick for import. Branch units remain those of the saved stage;
exporting a tree does not turn branch length into a transmission estimate.

`export.json` records the source manifest checksum, selected block/cohort when
applicable, sample count, scope and output checksums. Export files are portable
viewer inputs; they do not replace ChronoClade stage manifests for resuming or and continuing the pipeline.
