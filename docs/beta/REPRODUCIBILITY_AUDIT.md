# Reproducibility and automatic-reuse audit

Audited locally on 17 October 2026 for REPRO-01. **Partial: two automatic-reuse
version gaps are reproduced. No production code or live-job fingerprints were
changed.** Existing prepared bundles, running native typing and tree jobs remain
untouched. These are cache/provenance tests, not biological validation.

## Observed contracts

| Stage/path | What binds reuse or import | Current conclusion |
| --- | --- | --- |
| `prepare` | A pre-existing output is refused. Import validates dataset and extra artifact checksums, sizes, schema and completeness. New result records the current software version. | Safe immutable publication; no automatic completed-result cache here. Historical imports being readable is legitimate and distinct from recomputation. |
| Native query typing reuse | Assembly SHA256, species, configured tool commands/files, tool version and reference database hashes. | Strong data/tool/reference binding. Local parsing/normalization method and ChronoClade source identity are not bound. |
| Public context `--resume` | Prepared-source hash, explicit LIN/HC scope, frozen export/catalogue hashes, native-reference descriptor where supplied, provider URL; completed bundle and artifacts are revalidated. | Strong request/data binding. Completed normalized context can be reused after a local normalization-method/version change; no source/method identity in its resume contract. |
| `cgmlst --resume` | Dataset and adapted-record hashes, all requested parameters, direct-DNA catalogue provenance and **ChronoClade version**. Numerical artifacts and saved selection inputs are checked. Completed `cgmlst.json` runs are refused. | Version changes are rejected for interrupted numerical reuse. Same-version source changes and scientific Python-library versions are not bound. |
| `tree` automatic reuse | Selection hash, parameters/policies/seeds, native-tool evidence and local assembly hashes. Exact IDs, source metadata, alignments and result artifacts are checked when loading. | **Current ChronoClade version/source/method identity is absent from the fingerprint.** A completed prior result is automatically returned after changing the current version. |
| `time` automatic reuse | Saved input-tree hash, dating seed/topology policies, override parameters and TreeTime tool evidence. Source snapshots and exact IDs/topology/assessment are checked. | **Current ChronoClade version/source/method identity is absent from the fingerprint.** The common early cache-return branch is independent of supported/unsupported dating status. |
| `refresh-clock-reports` | Validated saved dataset, numerical tree/coordinate checksums, preserved selection files. Recalculates clock-root diagnostics and regenerates figures/reports. | This is **not purely presentation refresh**. It computes new diagnostic results, but does not currently record a separate current diagnostic/runtime provenance or rendering provenance. |

`stage_artifacts.header()` records the software version, but `load_result()`
only requires a nonempty saved version; it does not compare it with the running
program. That is appropriate for reading historical results, but it does not
protect automatic reuse. A format `schema_version` describes the contract shape;
it is not a numerical-method revision.

Native-tool evidence records executable paths/hashes and Conda package metadata
when found. This is useful and demonstrably invalidates changed-tool cache keys.
However, an entry-point script's hash alone does not identify an updated Python
library behind that script, especially outside the supported pinned Conda route.
Relevant installed library versions/artifact provenance should also be recorded.

## Concrete reproducer evidence

Run without touching any real dataset or native binary:

```bash
pixi run python -m pytest validation/beta/reproducibility-reproducers.py -q -rx
```

Observed: **3 passed, 3 expected failures** in 6.33 seconds. Expected failures
are open defects, not successful acceptance gates:

| Reproducer | Observed result |
| --- | --- |
| Identical tree inputs/tools/parameters | Pass: completed result is reused; the stub core is called once. |
| Changed native-tool evidence | Pass: new fingerprint and a second core computation. |
| Interrupted cgMLST with changed current version | Pass: rejects resume before reusing numerical evidence. |
| Tree after current ChronoClade version changes | Expected failure: returns the old fingerprint/result; no new core computation. |
| Time after current ChronoClade version changes | Expected failure: returns the old fingerprint/result. The isolated reproducer uses an unsupported-date fixture; code inspection confirms this same cache branch precedes the dating-status branch. |
| Shared tree/time report with internal `label` plus public accession | Expected failure: visible Genome cell uses the internal label instead of the available accession. |

The version probes temporarily replace the current in-process header version in
an isolated test process. No installed file or live process is modified. They
prove the automatic version-blind reuse behavior; same-version source-change
risk is established by inspection of the fingerprint inputs, rather than claimed
as a separate real-package-upgrade experiment.

Evidence: `validation/beta/reproducibility-reproducers.py` and
`validation/beta/reproducibility-reproducers.txt`. Existing positive source-change,
metadata tamper, altered coordinate, changed tree and seed/topology-policy tests
remain applicable. Checksummed artifacts establish integrity of saved results;
they do not establish that those results were computed using the current method.

## Minimum proposed fix, before production edits

1. Add one shared execution-provenance helper, rather than separate ad hoc stage
   checks. Record package version, explicit numerical-method contract revision,
   relevant runtime-library versions and an installed ChronoClade source digest.
   Keep source identity separate from native tool/reference identities, which
   already have substantial coverage.
2. For the beta, conservatively hash all installed package Python source files.
   This also catches same-version development edits and transitive shared logic.
   It can cause extra recomputation after presentation-code changes, but avoids
   incorrectly excluding numerical code based on its directory name.
3. Bind automatic numerical reuse to that execution identity in cgMLST, tree,
   time, native normalization and context processing. Store it in checkpoints and
   published results. Historical artifacts remain readable; an old artifact
   without this identity is not silently upgraded into a newly verified cache.
4. Distinguish `computed_with` from `rendered_with`. An explicit packaging refresh
   may reuse frozen numerical artifacts with their original computation provenance
   and publish separately recorded rendering provenance. It must not relabel old
   numerical results as computed by the current program.
5. Split pure packaging refresh from optional **diagnostic recomputation**.
   Root/date recalculation needs its own numerical-method/runtime identity and
   must preserve the source numerical tree, coordinate and selection hashes.
6. Retain immutable job namespaces and prior valid results. On method/source
   mismatch, use a fresh namespace or a clear refusal requiring a new output;
   never rewrite an old fingerprint to make it appear compatible.

The shared helper can initially be a small module returning plain versioned
mappings. Avoid a plugin framework or a second job scheduler. A later optimisation
could use separate explicit numerical component revisions to allow more precise
reuse, once method ownership is clear and tests enforce revision discipline.

Directory-name classification would be unsafe today:
`report_components/tree_ordination.py` performs numerical MDS and negative-branch
handling; `report_components/lineage_groups.py` assigns groups used in selection;
`report_components/root_comparison.py` computes diagnostic fits. Merely excluding
`report_components/` from a numerical fingerprint would miss real method changes.

### Rollout and acceptance boundaries

Do not change these fingerprint contracts while the current large native typing
and bridge jobs are running. Introduce them after the live evidence is published,
then exercise fresh small interrupted/resumed jobs under the new contract. Preserve
old artifacts and explicitly label legacy reuse compatibility as unverified.

Minimum acceptance gates for the change:

- Same source/method/runtime/data/tools/parameters reuse safely.
- Changed scientific source, method revision or relevant dependency invalidates
  automatic reuse, including same-package-version source changes.
- Changed data, tool, reference and temporal policy continue to invalidate reuse.
- Old results still load read-only with their saved provenance.
- Pure report refresh preserves all numerical hashes and records new rendering
  provenance; diagnostic recalculation records its own method provenance.
- Reproducers above become ordinary passing regression tests, with no unexpected
  expected-failure suppression in the beta acceptance run.

This proposal is not yet implemented. REPRO-01 remains partial; successful existing
native runs and green source-integrity tests do not close the version/source gap.

## Identifier traceability: quick assessment

Analytical IDs remain stable across selection manifests, trees, profiles and
metadata. Prepared bundles retain public accession fields and crosswalk evidence;
interoperability exports maintain ID/crosswalk links. Those internal identifiers
should remain available in evidence files rather than be replaced by potentially
nonunique display names.

Reader-facing labels usually use `sample_labels`, which prefers run/BioSample/
assembly accessions and removes legacy accession-plus-PW suffixes. The quick
report reproducer nevertheless finds a gap: `report_components/stage.py` takes
an explicit `row['label']` before accession. A supplied or legacy internal label
therefore survives despite a public accession being present. Use the shared label
resolver for those visible tables/nearest-relative labels, retaining the exact ID
in the data contract and optional traceability attribute. No broad renaming of
scientific identifiers is proposed.


## Reader-label follow-up

The shared stage report now uses the existing accession-label helper for both
the selected-genome table and corrected-SNP nearest-neighbour cells. Its former
expected-failure probe is an ordinary assertion; `tests/test_stage_report_labels.py`
adds a permanent regression. The focused tree/time/report suites passed 35 tests
with 1 skip. This fixes only presentation; the two version-blind automatic-reuse
reproducers remain expected failures pending the shared provenance change.

A shared `execution_provenance.py` helper is now implemented and passes three
focused source/version/method/runtime tests. It deliberately hashes all installed
package Python sources, including numerical report components. It is not yet
wired into automatic stage reuse: that migration is deferred until the active
native runs finish. The version-blind reuse findings therefore remain open.

## Implemented integration and current acceptance

The shared execution identity is now bound into cgMLST numerical resume,
tree/time automatic reuse, native query typing and normalized public-context
resume. New results record `computed_with`; report packaging records
`rendered_with`. Historical artifacts remain readable with their original
provenance, but absent/different identity cannot silently become a current
numerical cache result.

`refresh-cgmlst-report` now defaults to packaging saved evidence without root/date
refits. `--recompute-diagnostics` is explicit and records a separate diagnostic
method/runtime identity. Both paths preserve the original computation identity,
source NJ tree, coordinate and selection evidence. The earlier descriptions of
version-blind reuse and mandatory diagnostic refresh are historical findings,
not the current implementation.

The independent read-only integration review passed 14 tests with zero expected
failures. Additional permanent stage tests cover same-version method changes,
version changes and historical read-only loading. A final full-suite source
freeze is underway; no final release pass is inferred from these focused tests.
