# Raw analysis-export resume: applied follow-up

This proposal is isolated in `validation/beta/context-batch-resume-sandbox/`.
The parent applied it only after the original frozen whole-ST native run
completed. It was not injected into that run. The parent reports 108 active
focused controls passed after applying raw-resume and the separate progress fix.
The scope is the shared cgMLST/HClink exporter and context retrieval workspace.
It adds no classifier, distance calculation or new analysis cache.

## Defect and minimum shared fix

Before this fix, the context provider created a new random scope directory on each
attempt. Its completed `matched_N.json` checkpoint is written only after native
normalization. Interrupting classification therefore loses ordinary access to
already downloaded raw profiles, despite their files remaining on disk.

The proposal uses deterministic `scope_N` directories inside the existing
fingerprint-bound retrieval workspace. The shared exporter captures capability
request, response, as-of timestamp and checksum together in one atomically
published `export_capabilities.json`, before exporting batches. Resume verifies
that snapshot against the same base URL, requested ordered source IDs, numeric
ID universe, download names and row scheme/organism scope. It uses the captured
job/group order, rather than discovering a newer job and mixing references.

Each accepted raw batch is published as one atomic directory:

```text
scope_0/profiles/
  export_capabilities.json
  organism_0/
    batch_cgmlst_cgmlst-573-2_0/
      receipt.json
      response.raw
```

The version-1 `chronoclade.analysis-export-batch` wrapper records exact base URL,
download, job, route, ordered numeric IDs and source UUIDs, row scope, accepted
receipt, returned/missing UUID accounting, response bytes and SHA256. It contains
no authentication value. The receipt fields remain compatible with the shared
`ProviderSource` representation. That source reader now discovers stable scope
directories and preserves recorded HTTP evidence when sealing future captures.
Explicit historical raw-source import retains its separately labelled legacy
capability-order evidence; ordinary resume never adopts legacy `.bin` files
without request receipts.

Every reuse checks wrapper schema/hash, immutable request equality, response
path containment, bytes/hash and parsing/identity reconciliation through the
existing parser. Tampered, truncated, unsupported, mismatched and unreceipted
published files fail closed. Fresh preparation cannot overwrite completed raw
checkpoints without an explicit resume request. Cached group/profile dates are
retained, not updated to imply a new download; audit provenance records each raw
batch reuse.

## Crash window

Raw bytes and their receipt are staged in a private `.pending-batch_*` directory.
Only after both writes complete is that directory atomically renamed to its
published name. An interruption between raw writing and receipt commit leaves
no accepted batch directory. Resume preserves the private unfinished evidence
and safely refetches that batch. It never adopts the unpublished raw file, and
requires no manual deletion. An unknown partial **published** directory is
corruption and fails closed; it is not treated as a normal private stage.
Rejected identity/parse-response bytes remain separate diagnostic artifacts and
are never used as accepted exports. Private pending artifacts currently have no
automatic garbage collection; existing output resource guards still count them.

## What is reused and what is recomputed

Only raw exports and captured export capability metadata are reused. After an
interruption in normalization, all available profiles are passed through native
normalization again. No old LIN assignment or normalized numerical result is
adopted by this helper. The current reference/config and package execution
fingerprints continue to govern context resume.

A focused integration test uses the actual shared downloader/parser/provider,
with simulated source transport and an interrupted classifier. All 201 raw
profiles are acquired once; resume makes no repeated catalogue, capability or
cgMLST download request, while classification is invoked again on all 201.
Changing the reference configuration rejects resume before those operations.

The specialized Klebsiella LIN exporter is unchanged and may still make its own
separate network requests on retry. This is **not** a guarantee of zero network
traffic throughout every context workflow. Its unchanged source job and as-of
provenance are retained; current-reference native reassignment is still rerun.

## Important separate reference-fingerprint gap

The existing `reference_descriptor` binds config bytes/path and species, but not
all reference resources beneath unchanged config paths. The integration check
therefore proves rejection of changed configuration, rather than in-place
mutation of a database with an unchanged configuration. Underlying resource
checksums should be added through the existing shared provenance helpers as a
separate reviewed fix. That implementation is now isolated in
`validation/beta/reference-binding-sandbox/`; see
`NATIVE_REFERENCE_BINDING_REVIEW.md`. No such mutation was made to active
reference files during this acceptance run.
The old frozen full-ST source replay remains unaffected by this isolated work.

## Acceptance and safe application

Synthetic transport tests cover interruption after one complete batch,
interruption after all exports but before classification, immutable request
mismatch, hash/truncation/parser identity failures, job/order/route/reference
scope mismatch, credential omission, unsupported schema, unknown legacy files,
fresh-output collision, the raw-before-receipt crash window and safe refetch.
The shared HClink control retains blank levels, known assignments and unknown
reference-release metadata. A future-source roundtrip preserves recorded exact
request evidence. Broader exporter/provider regressions use existing fixtures.
See `regression.txt` for the completed run's exact count and duration.

After all active native jobs are terminal and the parent releases the source
freeze, review and apply only the focused patch:

```sh
git apply --check validation/beta/context-batch-resume-sandbox/focused.patch
git apply validation/beta/context-batch-resume-sandbox/focused.patch
.pixi/envs/default/bin/pytest -q tests/test_export_checkpoints.py \
  tests/test_profile_inputs.py tests/test_prepare_provider.py tests/test_provider_sources.py
```

Then run the broader required suite and create a new source freeze. Existing
execution fingerprints intentionally reject an old package's numerical resume;
use the ordinary explicit sealed `--context-source` route for historical raw
captures. Do not copy old normalization checkpoints or rename old UUID scope
folders into the new workspace. Finally exercise an ordinary fresh small
`--fetch-context` run with an actual classifier interruption and `--resume` under
the same frozen source before claiming native end-to-end resume acceptance.
