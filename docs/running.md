# Running an analysis

ChronoClade uses four independent commands. Each stage saves a manifest and a
separate output directory; pass the saved handoff from one stage to the next.
See [the complete workflow](independent-suite.md) for options and interpretation,
or the [Apple Silicon Klebsiella walkthrough](mac-klebsiella-walkthrough.md)
for a recorded end-to-end example.

## Prepare inputs

Resolve a collection, accession list, assembly table or frozen profile bundle.
For assemblies, declare the species and configure a compatible typing caller and
reference database. Preparation records readiness and provenance; it does not
retrieve public context.

```bash
pixi run -e prepare chronoclade prepare records.json \
  --input-kind profiles --catalogues catalogues.json --out prepared
```

Inspect `prepared/readiness.json` and the report. Continue with
`prepared/dataset.json` only when the intended inputs and typing readiness are
suitable for the analysis.

## Analyse cgMLST and select representatives

```bash
pixi run -e cgmlst chronoclade cgmlst prepared/dataset.json \
  --fetch-context --subsamples 1 --context-size 10 --out cgmlst
```

The default profile mode needs compatible species/scheme/version/database scope,
not LIN or HierCC assignments. Missing ST does not block local profile analysis,
but prevents a same-ST public search for that input. All inputs and mandatory
nearest-neighbour ties remain selected even when they exceed `--context-size`.
The budget does not cap the number of public profiles retrieved or analysed.

The cgMLST stage retrieves compatible public profiles and metadata when requested,
analyses the full comparable block, and writes one or more saved selections. It
does not acquire context assemblies. Open `cgmlst/index.html`, review the block
report and selection evidence, then use the generated
`selection_bundle/ensemble.json` or an individual selection JSON as the tree
input. See [cgMLST options](cgmlst.md).

## Build the tree and assess temporal support

```bash
pixi run -e tree chronoclade tree \
  cgmlst/BLOCK/selection_bundle/ensemble.json --out tree

pixi run -e time chronoclade time tree/tree_ensemble.json --out time
```

Tree acquires assemblies for the exact selected IDs, infers a recombination-
adjusted SNP tree and records temporal evidence. Time reads that saved tree and
assessment; it does not rebuild the phylogeny. Unsupported temporal evidence is
reported explicitly without a dated tree. For one saved selection, use the
individual selection JSON, then pass `tree/tree.json` to `time`.

Keep stage directories together so manifests and linked evidence remain valid.
Each output README records the next handoff and recovery guidance. The same
commands can be run through a generated job script; see
[local and SLURM jobs](independent-suite.md#local-and-slurm-jobs).
