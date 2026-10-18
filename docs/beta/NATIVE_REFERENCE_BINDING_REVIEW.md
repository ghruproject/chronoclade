# Native reference content binding: applied follow-up

The version-2 full-ST run completed with exit 0 under its original frozen source.
The parent verified all 637 original reference resources before and after that
run (`validation/beta/source-freeze-v2-after-wholest.json`). This review does not claim
that the older configuration-only fingerprint detects arbitrary later mutations.
The reviewed fix was applied after releasing that freeze. Its isolated review
is retained in `validation/beta/reference-binding-sandbox/`. The applied
reference/export/progress controls passed 155 tests; a separate read-only check
bound the seven actual reference resources and confirmed that SQLite reads left
the descriptor unchanged (`validation/beta/reference-binding-applied-readonly.json`).
Real native interruption acceptance remains separate.

## Defect

Before this fix, `reference_descriptor` recorded configuration path/SHA256 and species.
Changing a scheme, profiles database or allele lookup at an unchanged configured
path can therefore leave the context resume fingerprint unchanged. Existing
normalized scope checkpoints or a completed context preparation could be reused
against different actual reference content. This is separate from raw download
checkpoint recovery.

## Small shared change

The same `reference_descriptor` API now creates an explicit
`chronoclade.native-typing-reference`, schema version 2, with a content checksum.
It uses existing `file_sha256`, `content_hash` and `_provenance` helpers:

- Assignment scheme, profiles, allele lookup and ordered-locus metadata file
  hashes, plus existing tool/reference provenance and verified scheme correction
  receipt. Prepared HClink reference directories use stable per-file content
  references when present; no prepared HClink database availability is asserted.
- Caller `allele.db`, caller scheme metadata and its referenced scheme-details
  file, plus ordered indexed/canonical panel hashes and existing tool provenance.
  These are the resources used for canonical profile translation. Assembly caller
  caches retain their separate full-index provenance; this descriptor does not
  pretend to checksum every BLAST index file.
- Configuration path/hash and species, with a final descriptor SHA256. The
  configuration is checked again after binding to reject a concurrent change.

There is no second custom hash algorithm or copied reference database. Nonempty
SQLite WAL/journal files are refused. Directory references exclude only sidecars
belonging to recognizable SQLite databases: SHM and empty WAL/journal state.
Unknown similarly named files remain part of the reference. Normal read-only
SQLite lock activity does not change a content binding.

`resolve_reference` verifies the saved schema/hash and actual current descriptor,
rather than checking only the config file. Legacy config-only descriptors remain
readable historical metadata, but cannot silently authorize automatic native
reference reuse. An explicit existing `--typing-config` can bind a **fresh**
analysis; the resulting context records a migration from the previous descriptor
hash to the current content hash. Historical descriptors and results are not
rewritten. A changed content hash still rejects an old numerical resume, even
when the same config path is explicitly supplied.

The provider calculates the descriptor once for its fingerprint/final metadata
and records explicit rebinding provenance. Native normalization, nearest-reference
matching and LIN thresholds are unchanged by this proposal.

## Evidence and limits

The isolated descriptor/normalizer/provider controls cover actual resource
mutation under unchanged configuration bytes, all assignment files, caller
allele lookup and metadata, legacy/unsupported descriptor refusal, explicit
rebinding, nonempty-WAL refusal, and recognized transient-sidecar handling.
Two provider tests mutate the actual profiles-file bytes at the same path: one
after interruption before `matched_0.json`, another after completed preparation.
Both reject resume before another catalogue/export operation.

A read-only control uses the real corrected configuration and seven actual
canonicalization/assignment resources, then opens both SQLite references read-only
and compares the complete descriptor again. It passed in approximately one second;
raw evidence is `real-reference-readonly-control.json`. No reference/config bytes,
API requests or native jobs were changed. The evidence includes the isolated
implementation SHA256 and active frozen package SHA256.

`regression.txt` records the exact isolated descriptor/normalizer/provider test
count. Existing normalizer fixtures were strengthened to include the metadata and
resource files now required by a real content binding; the checksum safeguards
were not bypassed. The separate combined sandbox applies both proposed patches
only to copied files and exercises their shared exporter/provider paths together.
These are technical cache-contract tests with simulated transports/classifiers,
not a substitute for the live biological whole-ST acceptance.

## Review/application after the active run

The raw batch and progress patches are now applied by the parent. The final
reference-binding patch was rebased onto that active source and passes
`git apply --check`; 71 focused controls passed in 10.87 seconds. After parent
review and permission to change the source, apply only this remaining patch:

```sh
git apply --check validation/beta/reference-binding-sandbox/focused.patch
git apply validation/beta/reference-binding-sandbox/focused.patch
.pixi/envs/default/bin/pytest -q tests/test_native_reference_binding.py \
  tests/test_export_checkpoints.py tests/test_pasteur_profile_normalization.py \
  tests/test_profile_inputs.py tests/test_prepare_provider.py tests/test_provider_sources.py
```

Run the remaining required suite and publish a fresh local source freeze before
native acceptance. Do not force existing v2 checkpoints into the new numerical
cache. Explicit `--context-source` still provides the normal route to reuse
independently validated historical raw sources with fresh reference bindings.

## E. coli profiles-only compatibility review

The full-assignment descriptor contract is preserved; no profiles-only descriptor
mode, new API or preparation inference was added. The ordinary valid route
(`type-queries --profiles-only`, then `prepare --query-typing ...` without a
`--typing-config`, then `cgmlst --partition-by st`) does not request this full
native descriptor. An extended existing caller/import regression makes
`reference_descriptor` raise if invoked and verifies that ordinary E. coli
profiles-only preparation still succeeds, retains its full locus universe, and
leaves HierCC explicitly not requested.

A caller-only configuration explicitly supplied where a full reference binding is
requested now gives a clear `DatasetError` requiring caller and lineage readiness
and directing profiles-only import to omit that full-reference option. There is
no automatic fallback and no invented lineage readiness.
