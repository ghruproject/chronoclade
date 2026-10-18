# Native typing beta audit

> Follow-up, 17 October 2026: Pasteur authorisation and both Klebsiella reference
> builds now succeed. The query configuration is ready. Eight exact assemblies
> have been typed; reported cgLIN codes agree 8/8, and all 14 changed allele
> identifiers resolve to identical DNA. See
> `validation/beta/native-typing-inputs/NATIVE_RESULTS.md` and the acceptance
> matrix. Findings below retain the initial pre-download audit context; E. coli
> and complete historical reference equivalence remain unverified.

Date: 17 October 2026. Scope: Pathogenwatch cgMLST/cgLIN readiness, exact
identity and reference handling, existing evidence, and an offline concordance
comparator. This audit does not establish biological equivalence.

## Finding

The native typing software is installed at the pinned upstream commits, but the
local setup is not ready for production typing. Neither species has a complete
ready configuration. The only checked-in exercise of the native executables is
a synthetic three-locus integration fixture. The repository does not contain a
paired set of frozen Pathogenwatch benchmark assemblies and upstream assignments
with matching checksums and reference versions. No native-versus-Pathogenwatch
biological concordance result can therefore be claimed.

## Current evidence and status

| Area | Status | Evidence and limit |
| --- | --- | --- |
| Pinned executables | Installed, verified by local manifest | `~/.local/share/chronoclade/typing/typing_setup.json` records mlst 8.0.0, plincer 7.0.0 and hclink 4.0.1 at the pinned commits. This confirms software installation only. |
| Klebsiella native database | Blocked | The same manifest reports `ready=false`, with no `klebsiella_1` caller index or scheme version and no complete plincer `scgMLST629_S` reference. The generated `query_config.json` also has `ready=false`. |
| Pasteur authorization | Blocked / external | No local Plincer access-token cache was present at the expected credential-cache location. The historical client-token request returned HTTP 200, but that does not prove the account's one-time approval or authorize a database build. No request or build was made in this audit. |
| E. coli/HierCC database | Blocked / external | The manifest reports `ready=false`; its hclink metadata, index and loci file are missing. The presence of other files alone does not establish a usable, version-compatible database. The configured API key was not read or used. |
| Native integration | Implemented and tested on synthetic fixtures | `validation/query_typing_native/README.md` and `results.json` record the actual pinned tools against controlled three-locus fixtures for both species. The document explicitly limits this to integration, not biological typing validation. |
| Frozen Pathogenwatch catalogue | Partial context evidence | `validation/st147_pathogenwatch/frozen_catalogue.json.gz` contains cgLIN assignments for 5,647 public rows, but no cgMLST profiles and no paired native assemblies. It cannot satisfy the typing benchmark. |
| Exact benchmark assemblies | Missing | The checked-in `validation/` files include manifests, ledgers and analysis outputs, but not a frozen assembly set paired to exported upstream cgMLST and cgLIN results. Existing exact accession mappings do not supply the assemblies and reference-version crosswalk required by the benchmark. |
| Import and retrieval safeguards | Implemented; tests pass | Public assignments join by exact Pathogenwatch source ID; duplicate conflicting profiles remain unresolved. Query typing reuse requires matching assembly SHA256, species and native tool/database provenance. Caller locus order is checked before assignment. These safeguards do not substitute for a concordance run. |

The existing uncommitted Pasteur authentication changes in `typing_setup.py` and
`tests/test_typing_setup.py` were preserved. This audit did not inspect or copy
credential values, start OAuth, call EnteroBase, download references, or mutate
an external account.

## Readiness defect found and fixed locally

`load_typing_config()` previously accepted an empty cgMLST `index_dir` as long
as the directory existed. An isolated temporary E. coli configuration marked
ready reproduced the false readiness. The loader now validates the caller
metadata/version, compressed scheme details, ordered locus list, allele FASTA
files, and allele lookup SQLite schema/content, then checks locus order against
the lineage assignment reference. The setup report uses the same caller-index
validator and checks nonempty, parseable lineage metadata and loci. Regression
tests cover both an empty runtime index and a nonempty directory containing
only an unrelated file; both are rejected. This closes a preflight defect, not
a typing-equivalence question.

The stricter E. coli setup readiness check also invalidated an earlier test
fixture that contained placeholder text instead of parseable metadata/loci.
That fixture now uses minimal structured metadata and a two-locus list. No
production reference database was added or repaired during this audit.

No retrieval defect was demonstrated in the existing source-ID or assignment
paths during this audit. The `typing_setup.py` readiness test covers missing
reference files and the public scheme builder's isolated manifest; broader
production database compatibility remains untested because the references are
not ready.

## Comparator added

`chronoclade/typing_concordance.py` provides `compare_typing_results()` for two
offline lists of normalized assignment rows. It requires unique exact sample
IDs, matching assembly SHA256, supported matching species, matching method
scheme/version and equal locus universes before comparing calls. A declared
database version or fingerprint must be present on both sides and match. If
neither is declared, observed calls remain visible but reference pinning is
explicitly unverified and full equivalence is ineligible. A full cgMLST result
also requires `cgmlst_locus_universe_complete: true` on both rows; a locus list
by itself does not prove that it covers the canonical scheme. cgMLST output
separates allele discordance, callability discordance and
missing-versus-ambiguous discordance, and lists each affected locus. cgLIN
compares the entire parsed prefix, complete or partial status, and provisional
flag. It reports one-sided samples and incompatible scopes rather than silently
intersecting, substituting a neighbor, or counting a version mismatch as a
biological disagreement.

The comparator is deliberately offline and is not wired into a CLI. Its tests
cover matching known and novel alleles, different alleles, called-versus-missing
calls, missing-versus-ambiguous calls, provisional/partial LIN differences,
assembly and reference mismatches, missing database pins, incomplete canonical
locus universes, unpaired samples and different locus universes. Those tests
validate comparison behaviour only; they are not native equivalence evidence.

## Minimal reproducible concordance plan

1. Freeze an input manifest containing exact Pathogenwatch genome IDs, sample IDs,
   species, assembly paths and SHA256 values. Include easy and difficult or
   incomplete assemblies. Keep each downloaded byte-for-byte assembly used for
   native typing.
2. Freeze the corresponding upstream cgMLST calls and full ordered locus
   universe and an explicit `cgmlst_locus_universe_complete` flag,
   missing/ambiguous call status where available, and cgLIN assignment.
   Record both scheme names and exact scheme/reference versions. Record tool
   versions/commits and reference database SHA256 values for the native run and
   the upstream export when available. Do not join by name, accession, order or
   nearest result.
3. Build the native databases only after the required upstream authorization is
   confirmed. Verify the setup manifest and `query_config.json` are ready, then
   run the same checksummed assemblies through `type-queries` with the pinned
   native tools.
4. Load both `assignments` arrays and call
   `compare_typing_results(native_rows, frozen_rows)`. Save its JSON report beside
   the immutable input manifest and upstream export. Inspect every identity
   mismatch, incompatible version, allele/callability locus and LIN discordance.
5. Acceptance requires all intended samples to have exact paired identity and
   compatible scheme/reference versions; every cgMLST call and callability
   difference must be explained from the raw assembly/tool evidence. Exact cgLIN
   equivalence is claimed only for rows using the same declared LIN reference
   version. Provisional and unresolved assignments remain separate outcomes.

The current hard blockers are the missing ready reference databases and the
missing matched biological benchmark inputs. The empty-index preflight issue is
an implementation defect to fix before relying on the readiness gate. E.
coli/HierCC evidence remains separate from Klebsiella cgLIN and cannot be
inferred from the Klebsiella synthetic fixture.
