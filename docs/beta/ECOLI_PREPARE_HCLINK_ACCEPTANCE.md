# E. coli live-prepare HClink acceptance check

## Current result: provider import passes, reference scope remains unknown

After correcting the provider export route and parser, the ordinary `prepare` command was repeated for `DRR120871` / `SAMD00110466`. It exited successfully and imported both `cgmlst-562-1` and `hclink-562-4`. Exact requests, saved response hashes, returned identities and imported fields reconcile in `validation/beta/ecoli-hclink-final-live/verified-run/acceptance-evidence.json`.

The cgMLST export contains 2,500 observed loci. HClink assigns HC50, HC100, HC200, HC400, HC1100, HC1500 and HC2000 to 7, and HC2350 to 1. HC0, HC2, HC5, HC10 and HC20 are blank and remain unassigned. The reported closest profile and distance are retained as evidence rather than copied into the query profile identity.

Neither export supplies a database release or fingerprint. The locus catalogue remains explicitly incomplete and lineage readiness remains false because a comparable reference scope is unknown. This passes the single-record provider retrieval/import control; it does not validate native HierCC assignment for a new assembly or a biological benchmark. No assembly download, native typing or metadata enrichment ran.

## Earlier pre-fix run: retained regression evidence

This earlier single-query integration control predates the implementation described above; it is retained as before-state evidence, not the current result. On 17 October 2026, ordinary `chronoclade prepare` resolved the public Pathogenwatch identity for `DRR120871` / BioSample `SAMD00110466` and imported its server-provided *E. coli* cgMLST profile. That run did **not** retrieve or import the advertised HClink/HierCC assignment.

The exact invocation was:

```console
pixi run chronoclade prepare validation/beta/ecoli-hclink-prepare-one/input.csv --input-kind accessions --species 'Escherichia coli' --no-enrich-metadata --out validation/beta/ecoli-hclink-prepare-one/prepared
```

Exit status was 0. The command, one-accession CSV, stdout/stderr, status, frozen provider snapshots and their SHA-256 checksums are retained under `validation/beta/ecoli-hclink-prepare-one/`. `acceptance-evidence.json` is the compact machine-readable record; `prepared/prepare.json` checksums the published preparation bundle.

The live provider resolved the accession to public source genome `iLGx7naVnHXTvwPF8DEgAo` (sample ID `PW_iLGx7naVnHXTvwPF8DEgAo`), recorded ST131, Vietnam and collection year 2014, and downloaded the advertised `cgmlst-562-1` job. The imported profile has 2,500 observed loci and is marked `resolved`; the database version and fingerprint are unavailable, and its generated locus catalogue is explicitly incomplete. The provider snapshot records the job-level provenance and exact downloaded-export checksum rather than claiming a database release.

HClink was not requested in this preparation. The resulting provider snapshot contains no HierCC/HClink fields or HClink export entry, `lineages.jsonl` is empty, and readiness says the resolved HierCC assignment is unavailable. HC0–HC20 therefore remain unknown/unassigned in this result; no levels or database version were inferred. Readiness is `incomplete`, as expected when a complete locus catalogue and resolved HierCC assignment are required. This is a real provider-path gap for the requested acceptance, not evidence that the public HClink job itself is unavailable: a separate prior read-only API check for this same genome observed that job, but that response was deliberately not injected into this ordinary-prepare test.

No assembly was downloaded, no native caller or assigner ran, and ENA enrichment was disabled. This check neither validates biological concordance nor demonstrates an E. coli analysis report. The only successful claim is exact public identity resolution plus provider cgMLST profile import.
