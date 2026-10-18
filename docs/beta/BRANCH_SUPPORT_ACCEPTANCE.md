# Optional branch support acceptance

17 October 2026. Local engineering evidence; biological validation remains open.

`chronoclade tree --support-replicates 1000` enables a separate IQ-TREE UFBoot2
run with bootstrap NNI on the ClonalFrameML-filtered alignment. The existing
corrected tree and time-stage input topology are unchanged. A split comparison
assigns values only when an identical unrooted split occurs in the supported ML
tree. A missing split is not labelled zero-support or given a borrowed value.

Method contract and source:
https://github.com/iqtree/iqtree.github.io/blob/master/doc/Command-Reference.md
At least 1,000 replicates are required; 0 disables. Seed, model, replicate count,
source alignment/corrected-tree hashes, inference artifacts and matching-split
counts are retained. Parameters invalidate single/ensemble resume when changed.

## Inspected evidence

- `tests/test_tree_support.py`: root-independent split matching, absent splits,
  exact selected IDs, invalid replicate counts and disabled path.
- `validation/beta/branch-support-tests.txt`: 30 passed across support and tree.
- `validation/beta/tree-support-integration-tests.txt`: 50 passed, 1 skipped across
  support/tree/time/export. Tests using native stubs remain engineering checks.
- `tmp/beta-branch-support/smoke_result.json`: actual installed IQ-TREE,
  synthetic eight-tip 4,000-site alignment,1,000 replicates,5.12 seconds; all five reference
  splits present and assigned support. This is a native execution check, not a
  biological benchmark or claim of calibrated support for Klebsiella.
- `validation/beta/branch-support-native-smoke.txt`:native result and artifact
  hashes. Saved IQ-TREE report and split CSV are beside the smoke result.

The source native smoke script will be frozen beside the evidence. A real
Klebsiella 12-genome preregistered comparison is running separately under
`validation/beta/bowers2015/`; do not count it as passed before inspecting results.
