# ChronoClade

ChronoClade compares bacterial input genomes with public context, builds
recombination-adjusted SNP trees, and assesses whether collection dates support
a time tree. The current workflow uses four independent commands with saved
inputs and separate reports.

| Stage | What it does | Main handoff |
| --- | --- | --- |
| `prepare` | Resolve collections, accessions or assemblies; obtain typing and metadata | `dataset.json` |
| `cgmlst` | Retrieve matching profiles, analyse the complete compatible block, select representatives | `selection_bundle/ensemble.json` |
| `tree` | Acquire selected assemblies, infer corrected SNP trees, assess temporal signal | `tree.json` or `tree_ensemble.json` |
| `time` | Date supported trees and compare alternative selections | `time.json` or `time_ensemble.json` |

[Install ChronoClade](installation.md){ .md-button .md-button--primary }
[Run the four-module workflow](independent-suite.md){ .md-button }

## Start with profiles

The cgMLST module uses categorical allele mismatches by default. Missing calls
are excluded from pairwise comparisons, with shared-locus coverage recorded.
Direct nucleotide comparisons within the cgMLST genes are an optional mode when
a compatible allele sequence catalogue is available.

Species and MLST ST define the public candidate pool. Compatible LIN or HierCC
assignments determine the requested context block. Profile analysis includes all
usable genomes in that block; subsampling prepares smaller assembly datasets for
the tree and time stages.

## Read the results

The cgMLST report starts with the combined input and comparison dataset, then
shows nearest relatives, a full NJ tree and NJ-distance MDS, exploratory date
regression, state-change networks and the saved selections. Tables default to
10 rows. Trees support zoom and accession lookup; group selectors are compact,
searchable and linked between NJ and MDS.

[Report guide](report.md) · [Output-folder guide](outputs.md) · [cgMLST details](cgmlst.md)

Each report embeds its figures and interactive views. Keep the output directory
for downloadable evidence and stage handoffs. cgMLST output includes READMEs and
categorized file navigation.

## Interpret within the evidence

Profile similarity and reconstructed location or host changes do not establish
direct transmission. Root-to-tip regression is exploratory. Tree-stage temporal
testing and supported time-tree estimation remain separate steps. Variation
between alternative representative sets measures selection sensitivity, not a
population confidence interval.

Typing availability, assembly access and metadata completeness depend on the
provider and reference configuration. Missing evidence remains explicit.
