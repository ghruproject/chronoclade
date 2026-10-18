# Live E. coli HClink import: completed implementation plan

## Status

The implementation described below is complete. Focused parser/provider regressions passed, and an ordinary live `prepare` run imported both advertised E. coli cgMLST and HClink exports. See [the current acceptance result](ECOLI_PREPARE_HCLINK_ACCEPTANCE.md) and `validation/beta/ecoli-hclink-final-live/verified-run/acceptance-evidence.json`. This verifies provider assignment import for one exact public genome; it does not establish a reusable/versioned HClink reference or new-assembly assignment.

## Initial failure, retained as implementation history

The earlier single-query live `prepare` check in `validation/beta/ecoli-hclink-prepare-one/` resolved `DRR120871` / BioSample `SAMD00110466` to public Pathogenwatch genome `iLGx7naVnHXTvwPF8DEgAo` (numeric ID 108048). Its capability snapshot advertised both `cgmlst-562-1` and `hclink-562-4`. Before the implementation, ordinary `prepare` imported the former but not the latter because the exact-query route only selected cgMLST and Klebsiella cgLIN jobs. The current result is documented above; the original run remains retained as before-state evidence.

## Existing contract and frozen response

The independent documented HClink request for this exact numeric ID is preserved as `validation/beta/ecoli-hclink-distribution/pathogenwatch-hclink-562-4-id108048.csv.gz` (173 bytes; SHA-256 `c5ccab9dabf1f1ab4a831ea385c2819bac8db9ba5d187ab63556c9a29a5fdd76`). Its decoded one-row record and request provenance are in `pathogenwatch-hclink-562-4-decoded.json` and `pathogenwatch-hclink-result.json` in the same directory.

The gzip CSV has exact identity columns (`Genome ID`, `Genome Name`), closest-profile fields (`Closest ST`, `Distance to closest ST`), and HierCC level columns. For this row, genome ID is `iLGx7naVnHXTvwPF8DEgAo`, closest ST is `246645`, distance is `35.7`, HC0/HC2/HC5/HC10/HC20 are blank, HC50 through HC2000 are `7`, and HC2350 is `1`. Blank cells mean no assignment at those levels; they must be omitted from the `hiercc_codes` mapping and retained as blank-level/raw export evidence, never converted to cluster `0`. The response gives the job `hclink-562-4` and display name HierCC, but no reference release/date, HierCC database version, or database fingerprint. Do not report the job ID as a database version. Preserve server/job and export-file provenance separately from `hiercc_scheme_version` and `hiercc_database_version`.

The existing `_record_lineages` adapter already emits one lineage evidence row for `hiercc_*` fields, requires nullable scheme/database metadata, and carries all `hiercc_*` fields into `evidence`. `public_typing._normalise` already validates a level-to-cluster map and accepts optional database SHA, but it requires explicit scheme and scheme version. Avoid forcing live server results through this frozen-import validator if it would require inventing a version. `typing_scope` requires a scheme and version for compatibility; retaining unknown scope must therefore keep this assignment out of comparable lineage partitions/readiness until a valid shared scope is known. This is a conservative limitation, not a reason to discard the server's raw assignment.

## Intended implementation

1. Add a focused HClink CSV decoder, preferably `chronoclade/pathogenwatch_hiercc.py`, parallel in responsibility to `parse_cgmlst_export`. It should decode gzip (and the observed UTF-8 BOM), validate required columns and exact requested `Genome ID` set, reject repeated/unrequested IDs, normalize recognized `HC`/`d` level headers, omit blank levels from codes, preserve the raw blank levels and closest-ST/distance values, and return assignment records plus raw export checksum. Do not use nearest ST or distance to derive any HC value.
2. Extend the authenticated server-advertised export path in `profile_inputs.py` for `name == "hclink"`. Use the advertised job verbatim, bounded ID batches, `POST /api/downloads/hclink?job=<job>` with numeric IDs (the existing `request_download` credential/redirect helper), per-batch exact identity checks, retry only according to the established identity retry contract, and preserve each compressed response in the query-export bundle. Keep the existing cgMLST parser and assignment behaviour unchanged.
3. Route exact E. coli query rows with missing HierCC evidence through this path even when cgMLST profiles are already present. Keep the existing Klebsiella cgLIN route intact; do not broaden provider downloads to all advertised typing jobs. Group capability requests by organism as today, filter to `hclink` only for the configured E. coli/organism-562 case, and do not make broad context requests.
4. Import actual nonblank codes as `hiercc_codes` with `hiercc_status="resolved"` when at least one level is assigned (otherwise `unassigned`). Record `hiercc_source`, provider job ID, export SHA-256, exact raw level values including blanks, and a database-version-status field such as `unavailable`; leave `hiercc_database_version` null and do not set a fabricated release. A recognized scheme label may be stored only if it comes from the returned capability metadata; leave `hiercc_scheme_version` null unless Pathogenwatch explicitly supplies a scheme/reference version. Emit an `analysis_exports` audit entry for success/unavailable with requested IDs, advertised job, response checksum, and the explicit statement that database release/version/fingerprint were not supplied.
5. Keep lineage evidence and report/readiness semantics conservative: the downloaded partial HierCC assignment is useful evidence, but unknown database scope must not become compatible with a native local HClink reference or trigger lineage partitioning. HC levels that are blank must render as unknown/unassigned, not an invented code.

## Planned checks, now exercised

- Parser test using the retained real one-row gzip export: assert checksum, exact source ID, preserved closest ST/distance, HC50..HC2000=`7`, HC2350=`1`, blank levels `{HC0, HC2, HC5, HC10, HC20}`, and null database release/version/fingerprint. Keep the small raw file as a frozen public fixture with its provenance JSON; do not regenerate its content.
- Parser rejection tests: wrong/unrequested Genome ID, repeated IDs, absent required columns, duplicate level headers, corrupt gzip, malformed nonblank cluster IDs, and response omission of a requested ID. Blank cells must be accepted and never mapped to zero.
- `_analysis_exports` tests: advertised HClink job is selected only for eligible E. coli records; batches contain numeric IDs; the parsed Genome ID must map back to the exact resolved source genome; a failed/mismatched response leaves existing cgMLST intact and records an unavailable reason.
- `resolve_queries`/normal preparation regression: with a fake PathogenWatch transport returning the existing capability shape and the retained HClink response, an accession query downloads cgMLST and HClink once each, writes both checksummed provider artifacts, imports both profile and lineage evidence, and records no invented HierCC database release. A complete existing HierCC assignment prevents a redundant HClink fetch; non-E. coli calls and the established Klebsiella cgLIN path retain current behaviour.
- Readiness/report test: live provider HC codes are retained as imported evidence, but null database scope is not advertised as a comparable/refined lineage partition. Verify HC0–HC20 are shown as unassigned/unknown in evidence and cannot become numeric cluster zero.
- The one-row real configured PathogenWatch `prepare` command was rerun in a fresh output directory. Both advertised jobs were imported, response receipts match the saved bytes and returned identity, HC blanks remain absent from codes, and the database version remains null. See the current acceptance evidence linked above.

## File ownership recommendation

Keep the implementation change narrow and assign one owner to these files after the source freeze lifts:

- `chronoclade/pathogenwatch_hiercc.py` (new): gzip/CSV parser and small normalized assignment result.
- `chronoclade/profile_inputs.py`: authenticated download and export audit integration; leave CGMLST and cgLIN parsing paths unchanged.
- `chronoclade/prepare_provider.py`: exact-query routing for missing E. coli HierCC, preserving existing cgLIN behavior.
- `tests/test_pathogenwatch_hiercc.py` (new): parser tests against the real retained response and rejection cases.
- `tests/test_prepare_provider.py` and/or `tests/test_profile_inputs.py`: narrow fake-transport routing and regression tests.
- `validation/beta/ecoli-hclink-prepare-one/` retains the pre-fix run; `validation/beta/ecoli-hclink-final-live/verified-run/` retains the post-fix live acceptance and response receipts.

The implementation and its live acceptance are complete. The distribution of a compatible indexed HierCC reference and native assignment of new assemblies remain separate open work.
