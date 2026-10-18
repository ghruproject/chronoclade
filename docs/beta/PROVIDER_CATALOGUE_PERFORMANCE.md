# Public catalogue acquisition: first-run bottleneck

This is a read-only implementation and frozen-response audit alongside the bounded
ordinary whole-ST live run. No endpoint experiments, mock clients, source changes
or manual candidate filters were used.

## What the current code requests

`PathogenwatchClient.search` retrieves all search pages (100 records per page) and
reconciles their exact UUID set against the server count. `freeze_catalogue` then
requests `/api/genomes/details?id=<UUID>` for every result, with six workers by
default. The complete catalogue is written atomically only after all details are
received and normalized. A partial set of detail responses is not checkpointed.

The existing frozen ST258 snapshot contains 13,420 members, 135 search pages,
13,420 detail records and a 55,547,052-byte catalogue. Thus metadata acquisition
requires about 13,555 search/detail requests before the profile export work.
The bounded live run does not assume that today's server count is unchanged.

The requests are not redundant under the current input contract. Its actual
cached search rows contain:

- Numeric ID, UUID, name, organism label, typing string and public-project flag.
- Location and `metadataDate`, plus project/owner display fields.

They do **not** contain the full metadata dictionary, `startDate`/`endDate`
intervals, explicit organism ID and both MLST namespaces, QC result, assembly
checksum, sample/run/study accessions or host/isolation-source fields consumed
by `normalize_record`. Search `metadataDate` is insufficient to reconstruct the
existing date-precision/interval checks. Using search rows alone would discard
identity, geography, quality and missing-data evidence.

## Existing reuse paths

`load_catalogue` already validates snapshot/section hashes, public-only status,
query count and a complete replay from raw search/detail responses. The modern
profile-context provider reuses `public_0.json` through `--resume` **within its
own matching retrieval workspace**. Its execution fingerprint correctly prevents
reusing completed numerical context after a package/source change.

The supported `cgmlst --fetch-context` command does not expose a supplied
catalogue or shared metadata cache.
Its `--public-typing` accepts allele/lineage annotations, not an entire catalogue.
There is no normal flag to seed this run with the earlier complete frozen ST258
catalogue without modifying internal workspace files. The current acceptance
therefore uses the real client and live metadata acquisition.

The cached typed-profile envelope has 12,706 records from an older reported-LIN5
subset, not all 13,420 ST records. It cannot substitute for the whole-ST pool.
Its `{records: ...}` envelope is also not directly accepted by `--public-typing`,
which expects a list or `{assignments: ...}`. No hidden conversion/filter was used.

## Documented bulk metadata possibility

The already saved official OpenAPI document describes both
`POST /api/downloads/metadata` and `POST /api/downloads/summary`, taking a
comma-separated list of numeric genome IDs and returning a generated-file
redirect. These are documented bulk endpoints; the download capability list
otherwise advertises analysis jobs, not metadata fields.

The documentation does not specify export columns, and no saved real bulk
metadata/summary export was found in this audit. Therefore it is **not yet
verified** that either export contains every field required for public identity,
ST/scheme verification, collection-date intervals, QC, metadata conflicts and
checksum provenance. Do not replace validated details using a guessed schema.
A later bounded comparison of five exact public IDs against their existing
raw detail responses can establish coverage. Only absent fields would require
individual detail requests; exports must reconcile every requested ID, reject
unrequested/repeated identities and retain original bytes/checksums.

## Minimum shared improvement after the source freeze

1. Forward an explicit frozen `--catalogue` option into modern profile-context
   discovery, reusing `load_catalogue`. Verify organism/ST/MLST namespace, public
   status, exact expected counts and base URL; record source SHA256 and the
   snapshot's as-of time. No freshness claim should be inferred. Importing an old
   source snapshot is separate from reusing old numerical normalization.
2. Reuse the shared `progress_task` around search/details acquisition immediately,
   with periodic completed/expected counts. The current live log can otherwise
   remain silent after its initial input typing while thousands of requests run.
3. Add resumable source fragments: freeze search pages/UUID universe first, then
   atomically store verified detail responses by UUID/content hash. Reuse explicit
   frozen inputs across jobs; make refreshing metadata an explicit policy. Do not
   couple immutable raw-source availability to the numerical execution cache.
4. Once the real bulk schema has been validated, fetch bulk metadata/summary in
   bounded batches and use per-record details only for fields/identities the bulk
   exports cannot prove. Share this acquisition component across preparation and
   profile-context discovery rather than adding another independent client.

Acceptance for either path must show identical normalized rows, date intervals,
missing values, source identities and public-only scope against a frozen detail
fixture; interrupted acquisition must resume without repeating verified records.
No assembly download is needed for this work.

## Live evidence

Exact command, caps and measurements:
`validation/beta/native-wholest-live-run.json`.
Program log: `tmp/beta-native-wholest-live/cli.log`.
The run was capped at 30 minutes, 12 GiB sampled process-family RSS and
4 GiB newly allocated output. A budget stop or export failure is evidence of an
unmet end-to-end acceptance item; it must not be converted into a smaller pool
or a claimed complete result.

## Observed live acquisition phases

The ordinary live catalogue completed successfully with **13,421** records,
13,421 detail responses and 8,190 accession-deduplicated units. Its own elapsed
time was **195.893 seconds**. This is one additional source record compared with
the older frozen catalogue; an explicit cached import must retain its as-of scope.
Evidence is `validation/beta/native-wholest-catalogue-measurement.json`.

Profile acquisition then continued through the existing sequential 100-genome
export batches in `profile_inputs._analysis_exports`. All 135 saved profile exports arrived without identity-retry files. The resulting
worker request contains exactly 13,421 public rows, verified by full JSON parsing
after the run stopped. This proves acquisition and submission before LIN8
filtering, rather than completion of the assignment or final context preparation. The program log remained silent during these catalogue/export phases,
a usability gap independent of whether the bounded run succeeds.

After validated catalogue reuse, a second measured optimization is transport
batch sizing. The saved official application configuration advertises a download
limit of 5,000, but acceptance of larger cgMLST batches has not been established
here. Benchmark a few exact public requests at 100 versus 500 members before
changing a default. Preserve exact-ID reconciliation, redirect credential
isolation, retries and raw evidence at every size. Larger batches must never be
used to relax identity checks or silently omit unavailable genomes.

## Terminal bounded acceptance result

The unchanged guard stopped the command at **1,800.042 seconds**, exit **-15**.
The native all-ST normalization/reassignment subprocess was still running; its
last heartbeat was 660 seconds into that subprocess. Peak sampled process-family
RSS was **7,301,070,848 bytes (6.80 GiB)** and peak newly allocated output was
**329,289,728 bytes (314.0 MiB)**. There was no resource overrun: wall time was
the limiting acceptance criterion. These are concurrent local-run measurements,
not an isolated performance benchmark.

The frozen catalogue passed the real `load_catalogue` validator after stopping.
All 13,421 catalogue/detail records and all 135 profile export files were retained.
Full JSON parsing confirmed three input requests and 13,421 public requests in
the native reassignment worker inputs. No context assembly downloads occurred.
**The whole-ST normal CLI acceptance did not pass**: final current-reference LIN8
retained counts, unavailable-profile audit, prepared context and report were not
produced. Those denominators remain unknown, rather than zero. The evidence is
`validation/beta/native-wholest-terminal-acceptance.json` and the exact-command
resource receipt linked above. No smaller pool or manual filtering replaced it.

A resume limitation matters here. `public_0.json` is validated and reusable in
the matching retrieval workspace, but each new scope attempt creates a fresh
UUID directory before lineage/profile exports. The completed `matched_0.json`
checkpoint is written only after reassignment. Thus the interrupted raw exports
and in-progress assignment cannot be reused independently by the ordinary
provider. Preserve source checkpoints separately before optimizing arithmetic.

The parent's isolated exact-memoization proposal preserves upstream classification
functions and all tie handling; it needs complete old/new output equivalence
before a fresh bounded run. No optimization was applied to this frozen run.
Separately, the parent found upstream `scheme.toml` level-6 `max=798.8871` (outside
percentage range). The installed frozen reference was left unchanged. Its
provenance-preserving correction and boundary validation are a separate scientific
reference acceptance item; no complete typing-validation claim follows from this
throughput test.
