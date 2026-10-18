# Independent dating and selection sensitivity

`chronoclade time tree/tree.json --out time/` validates all tree/source/artifact
checksums, exact raw/corrected/rooted tips, selected metadata and callable-site
counts. It does not reacquire assemblies, select genomes, rebuild alignment or
infer a genetic tree. The saved root is retained when the existing TreeTime
maximum-likelihood dating engine estimates rates and node dates.

The temporal assessment is copied exactly from tree. Unsupported evidence yields
a complete `time.json` and report with `dating_status: unsupported` and no dated
tree; it requires no dating executable. `--allow-unsupported` is an explicit,
saved override and still needs at least three eligible dated samples. `--force`
only controls execution and never overrides the scientific gate.

A dated product saves NEXUS, node dates, rate uncertainty, TreeTime confidence
outputs, time-tree figures and the established geographic network reconstructed
on the dated biological tree. Rates are substitutions/site/year; branches are
years. Node-date intervals use TreeTime's reported 90% max-posterior regions.
With usable host metadata, the same network selector also offers host-state
changes on the dated tree. Location remains the default; missing host data does
not create an empty host section, and one known host state cannot create links.
Native tools and site counts follow the established engine, rather than implying
that permutation spread supplies dating uncertainty.

Fingerprint-based independent jobs support resume; failed dating leaves existing
tree/cgMLST products and the previous valid time manifest intact.

## Multiple selections

```
chronoclade time-compare cgmlst/block/selection_bundle/ensemble.json \
  --run selection-001/time/time.json --run selection-002/time/time.json \
  --out sensitivity/ --cluster-snp-cutoff 10
```

Every run is bound to the ensemble's exact selection hash. The comparison retains
unsupported runs and missing requested runs with their denominators. A failed
time job may be supplied by its `job.json` so its reason is explicitly retained.
Shared-anchor MRCA dates identify the same sampled anchors across runs; rates and
whole-selection root dates are saved separately because whole-selection roots
may represent different ancestors. Target clustering compares exact target pairs,
using a declared single-linkage corrected-SNP threshold; run-local cluster numbers
are never matched across runs. A threshold of 10 is a descriptive setting, not a
universal biological cutoff.

The report provides readable run, sensitivity and target-pair tables with JSON
evidence downloads. Minimum, maximum and range quantify selection sensitivity,
not confidence intervals. Unsupported/failed runs never contribute manufactured
dates or rates. No minimum/maximal spread proves adequate representative sampling
or convergence; scientific tolerances and coverage remain a separate judgement.

SLURM uses the same commands through the supporting `chronoclade job` operation.
Independent job namespaces keep separate selections and parameters reviewable.
A written batch script is not evidence of a real cluster execution.

Native dating retains the saved corrected-tree root and polytomies: TreeTime receives `--keep-polytomies` and `--rng-seed` using the saved tree manifest’s `parameters.seed`. The time manifest records `dating_seed`, `keep_polytomies: true`, and `dating_inference_policy: seeded_saved_topology_keep_polytomies`. These settings are bound into its fingerprint so earlier dating results cannot be reused under the new inference policy. Legacy lineage dating retains its existing defaults unless callers explicitly supply the optional settings.

Time ensembles validate explicit `--reuse-run` products against the same native dating parameters derived from each saved tree. Results with an older policy, a different dating seed, or polytomy resolution enabled are rejected for reuse; their existing manifests remain readable. The per-selection dating parameters are also bound into the ensemble fingerprint.

Before publishing a dated result, native dating compares the saved rooted Newick tree with the dated tree using the exact sets of tips beneath each nontrivial internal clade, excluding the root. Node labels and branch lengths are ignored. Any lost or introduced rooted clade fails the job even when all tip IDs match. Dated results record the method, both clade counts, and `preserved: true` in `topology_preservation`; unsupported results record null.

`topology_validation_policy` is included in the shared native dating parameters and both single-run and ensemble fingerprints. Dating products created before this automatic check remain readable but cannot be reused as automatically validated products.
