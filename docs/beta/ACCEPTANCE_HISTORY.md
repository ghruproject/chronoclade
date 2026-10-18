# Beta audit and acceptance plan

## Current priority: Klebsiella on this Mac

User direction, 18 October 2026: finish and verify the local Apple Silicon
assembly-to-report workflow first. Linux x86_64 verification is deferred until
the user supplies an environment; it does not block this local acceptance, and
no x86_64 execution claim is made. E. coli/HierCC work is parked pending the
user's reference-sourcing work. Historical evidence below is retained.

The immediate closeout is recovery-fix integration, final software checks, a
real preparation interruption/resume, and one fresh-installed Klebsiella
`prepare → cgmlst → tree → time` run with inspected reports and normal command
handoffs. Reuse verified references explicitly; do not present reuse as a fresh
reference build. Existing positive dating benchmarks remain valid within their
recorded scopes; a small run returning unsupported is a valid outcome.
Human walkthroughs and historical-reference limitations stay visible, but do
not expand this closeout into new cluster-method research or platform testing.
All work remains local; no pushes, pull requests or releases.

## Report, supported dating and walkthrough closeout

The user-authorised steps 1–3 are complete locally. V16 source is
`7d71f3631420a610c33693ae9b73150e580d29ca0d58e9825334d84c27128dd3`
(134 Python files). Report changes preserve the existing style: a separate
caption band below tree axes, upstream report links and a responsive container
for the nearest-context table. Desktop and 390 px browser checks passed;
the final narrow page width is 390 px with no page-level horizontal overflow.

The supported 63-genome time-stage check ran on the current V16 build from a
saved corrected SNP tree, without constructing another SNP tree. It reproduced
the prior MRCA estimate 1995.6604 and 90% max-posterior region
1994.268418–1998.078415, retaining all 63 IDs and the same 61 rooted internal
clades. See `validation/beta/positive-subsampling-plan/time-closeout-selection001-v16/`.
The root independently validated the saved result contract and interval.
The new [Mac walkthrough](../mac-klebsiella-walkthrough.md) covers the ordinary
user route, with the recorded cached validation separated from it.

V15 full software checks passed: 1,249 tests, 2 skipped, Ruff and strict docs.
The final V16 delta only adds the nearest-table responsive container; its
focused stage/renderer checks passed (48 passed, 1 skipped), Ruff and strict docs
passed, and the positive native time-stage run completed with unchanged inputs.
No GitHub publication, x86 verification or human walkthrough is claimed.
Closeout evidence: `validation/beta/mac-polish-positive-v16/`.

## Current local verification

18 October 2026: V14 source (`75b9470654fb46e7dc25a6cc269ee98fd9bcca3fda21cab1a0929121a0534478`,
134 Python files) passed the complete suite: **1,249 passed, 2 skipped**, Ruff,
and the strict documentation build. See
`validation/beta/source-freeze-v14-checks.json`. The first fresh V13 run exposed
a macOS `/tmp` versus `/private/tmp` path-containment bug after both assemblies
typed successfully. V14 resolves the bundle root before checking artifact paths;
a regression confirms valid aliases pass and escaping artifact links fail.
The original failure is preserved in `execution-v13/science-result.json`.

Real two-assembly preparation recovery passed: interruption after the first
completed sample, exact resume typing only the second sample, unchanged first
sample raw results, equivalent final profiles/LIN calls, and refusal of changed
inputs or completed targets. The uninterrupted baseline bundle was independently
validated and adopted after a post-publication runner scan race; its CLI exit
code is unrecorded. Evidence:
`validation/beta/native-typing-recovery-native-plan/evidence/native-acceptance-completion.json`.

The fresh macOS arm64 Pixi installation passed source/runtime/tool verification.
Its new-assembly four-module scientific run completed successfully: prepare
35 s, cgMLST 881 s, tree 540 s, time 58 s (about 25 minutes total).
The profile block contains 2 inputs and 232 public comparisons. The exact
28-genome tree retains all 26 mandatory nearest-neighbour contexts; 100/100
clustered permutations yield p=0.0891, so the time stage correctly returns
unsupported without a dated tree or override. Independent review reconciles
all 80 tree artifact hashes and exact identities across raw, corrected and
rooted trees. See `validation/beta/mac-native-final-plan/execution-v14/root-independent-terminal-review.json`.
Typing references and
the configured typing runtime are explicitly reused; public context is replayed
from a complete cached provider snapshot. This is not a new live provider survey
or a fresh reference installation. Evidence:
`validation/beta/mac-native-final-plan/execution-v14/`.

The preparation and cgMLST reports were inspected in the browser, including
embedded images, accession labels, pagination, NJ zoom, linked LIN highlights,
host availability and interactive network export. Tree/time reports retain the
existing style and display the unsupported dating decision. Standalone file
visual review and two human walkthroughs remain separate. A small corrected-tree
figure footer overlap is recorded as a cosmetic follow-up.

## Previous integrated checkpoint

Date: 18 October 2026. Status: beta incomplete. The current V11 source is frozen at `e7d8011e8ad2dc2ff6ecebf60bdac000ab507c8ae644f3988cdd91378937142a` across 133 Python files. V11 full pytest passed (1,230 passed, 2 skipped), Ruff passed, and strict docs built successfully on the frozen source; `validation/beta/source-freeze-v11-checks.json` records the unchanged source identity. The prior V10 source passed 1,223 tests (2 skipped), Ruff and strict docs; see `validation/beta/source-freeze-v10-checks.json`. The earlier V9 full suite failed two brittle source-file-count assertions; both focused corrected files passed 32 tests against V10. Do not combine checks from different source snapshots. V7 macOS arm64 and Linux ARM64 runs passed a supplied pretyped eight-record fixture, not new-assembly typing. Five positive 63-tip selections have independent fixed-target/tree comparison. Native reference equivalence, full native assembly typing, Linux x86_64/HPC, controlled E. coli real-context/HierCC validation and two human sessions remain open.
Governing brief: [BETA_SPECIFICATION.md](BETA_SPECIFICATION.md).
No work in this plan authorises GitHub publication.

## Active delivery goal

The user has requested completion of this entire specification through audit and
acceptance, rather than stopping after the first fixes. Execute in this order:

1. Normal CLI context refinement, with separate search, analysis-block and display
   depths; preserve measured nearest neighbours across display-prefix boundaries.
2. Real small four-module execution, exact identities, resume/errors, supported
   or unsupported temporal outcomes and checked artifacts.
3. Independent preregistered biological comparisons and resource measurements.
4. Clean macOS/Linux installation, usable walkthroughs, actual-output visual
   demonstration and two independent human attempts.

Keep all unresolved matrix entries visible. Full completion requires the evidence
specified below, including the human checks; engineering tests cannot replace
them. All implementation and evidence remain local until publication is approved. V7 scoped installation and biological results remain pinned to V7. The current V11 source is frozen and full software checks passed (1,230 tests passed, 2 skipped, Ruff and strict docs); the terminal receipt records unchanged source identity. V10 checks passed on its separate frozen source. Separate human, biological-reference and platform gates remain open. Do not combine results from different software snapshots into one claimed run.

Automatic selection of search, context, analysis-block or display depth is not an explicit acceptance gate in the governing specification. Keep depth-threshold evaluation as optional research unless the specification is revised.

## Evidence rules

`verified` means a reproducible check with inspected artifacts. `partial` means
only some required evidence exists. `unverified` means implementation/documents
exist without a matching acceptance run. `blocked` names an external prerequisite;
`defect` names implementation work. A passing unit test is not biological validation.
Historical runs below are evidence inventories, not newly rerun results.

## Acceptance matrix

| ID | Requirement and acceptance condition | Current status | Evidence or next check |
| --- | --- | --- | --- |
| SPEC-01 | Approved brief saved with ST258 grouping amendment; local-only rule | verified | V7 audit maintains the same 12 sections and 33 acceptance rows. Governing specification remains unchanged; no GitHub publication is authorised. |
| INST-01 | Clean Apple Silicon install and small documented run, no developer repairs | partial | Fresh supported Pixi macOS arm64 completed version, preflight, prepare, cgMLST, strict tree (100 permutations) and time on the exact supplied eight-record pretyped fixture (3 inputs, 5 contexts); identity and source/input pins were independently checked. Time returned normal unsupported p=0.376 without override. This does not test native typing of new assemblies. `validation/beta/install-v7-pixi/terminal-review.json` |
| INST-02 | Clean Linux x86_64 install and same fixture result | partial; scoped Linux ARM64 pass, x86_64 remains open | Independent V7 Linux ARM64 review confirms install, version, preflight, prepare, cgMLST, strict tree (100 permutations) and time on the same exact eight-record pretyped fixture; all stages exited 0 and pins reconciled. Time returned normal unsupported p=0.376 without override. This does not establish clean-assembly typing, Linux x86_64, physical Linux, HPC or SLURM. `validation/beta/linux-aarch64-acceptance-v7-independent-review/terminal-review-v7.json` |
| INST-03 | Tool install, reference readiness and credentials distinguished | partial | Caller and lineage readiness remain distinct. Assembly-only context discovery needs verified numeric ST metadata before --fetch-context; native cgMLST does not infer ST. V7 native101 identity packaging is verified, but the 101-profile dataset is not declared lineage-ready: one context resolves only through LIN6. Provider fingerprint/equivalence and compatible E. coli HierCC reference remain open. |
| UX-01 | Stage/sample/progress/failure summary and next command understandable | partial | Targeted V6 navigation checks and V7 external-wheel asset audit remain scoped, not accessibility or human acceptance. The earlier E. coli report review identified a verbose allele-rate label and an inaccurate unassigned-group caption; the V8/V9 fixes are integrated and tested, and a saved-point preview confirms the concise units and neutral tree caption. Fresh report regeneration against the current main source remains separate. Full progress/readiness/next-command usability, verified-ST handoff and two human sessions remain open. validation/beta/report-navigation-accessibility-v6.json; validation/beta/install-v7/acceptance-evidence.json; validation/beta/ecoli-report-current-v7/browser-review.json |
| UX-02 | Two target users finish walkthrough independently | blocked | Two independent target users are still required; no human sessions are evidenced. Agent or browser checks do not fulfil this gate. |
| PREP-01 | Assemblies/accessions/collection URL resolve exact identities and conflicts | partial | The V7 repaired ordinary frozen101 packaging independently verified all 101 source UUID identity joins, unchanged native calls/lineage fields/catalogue and unchanged original hash/size/mtime pins; the prior failed wrapper record is preserved. This identity/import evidence reuses saved V6 native calls; it is not a new typing run or proof of complete lineage readiness: 100 records have full ten-component LIN, while SRR1561312 resolves only through LIN6. See validation/beta/positive-subsampling-plan/native101-packaging-v7/independent-terminal-review.json and validation/beta/bowers2015/full101/native-input-v6/independent-audit.md. |
| TYPE-01 | Exact assembly native cgMLST versus frozen Pathogenwatch calls, including missing loci | partial | One literal `nagE_S` discrepancy was resolved earlier. Follow-up current official Pasteur FASTAs found no exact match to six historical provider locus/digest targets. Of 29 native hit rows, 14 span the full best-allele length and 15 are shorter fragments; all remain non-exact under the native caller. Historical provider catalogue bytes/fingerprint are missing, so equivalence remains unresolved; this is not a demonstrated contradiction. `validation/beta/native-reference-followup-v8/sequence-comparison.json` |
| TYPE-02 | cgLIN native concordance under pinned compatible reference | partial | V7 packaging preserves native lineage fields and exact identities; prior 12-profile LIN strings agree after delimiter normalization. Provider reference database version/fingerprint is unavailable, so this is not full reference equivalence or a claim that 101 records are fully resolved; one context resolves only to LIN6. validation/beta/bowers2015/full101/native-input-v6/independent-audit.md; validation/beta/positive-subsampling-plan/native101-packaging-v7/independent-terminal-review.json |
| TYPE-03 | E. coli/HierCC native equivalence | partial | A metadata-only ST131 search inventoried 22,604 unique UUIDs across 227 pages, but is not a full-detail catalogue. The query’s provider-reported HClink row has HC50=7; the 22,604 metadata candidates have yielded 15,425 HC50-matched candidate records; profile batch 23/155 was active at 2026-10-18 05:08 UTC, so the terminal accepted profile pool remains unknown. The reference release is absent. The 2,500-locus query catalogue is incomplete; V10 preparation failed when a same-scope context introduced a locus absent from the query catalogue. V11 integrates union expansion only for an incomplete catalogue in the exact same scheme/version/database/hash scope; complete declared catalogues remain strict, existing query calls and missing values are preserved, and no database version or full locus universe is invented. Isolated staging converted the preflight to a 2,501-locus union, preserving original calls; the added query call is null and the synthetic context call is `1`, so this does not support biological inference. V11 full software checks passed; a full ST131-to-HC50 context retrieval is now in progress, but no result is available yet. Controlled real-context validation and the E. coli bulk/final acceptance gate remain open. Native HClink assignment, compatible reference setup, full ordinary context/report workflow and independent biological comparison also remain open. A separate prior-source ST95 pilot compared 28 published/native profile-distance pairs (2/28 exact counts; all 28 agree on <=40 membership), but underlying profiles are unavailable and it does not validate temporal performance, source attribution or HierCC assignment; see `validation/beta/ecoli-benchmark-gap-review/review.json`. The report review remains presentation-only. See `validation/beta/ecoli-provider-context-plan/inventory-v9/result.json`, `validation/beta/ecoli-provider-context-plan/input-audit.json`, `validation/beta/ecoli-provider-context-plan/observed-locus-union-preflight.json`, `validation/beta/ecoli-provider-context-plan/observed-union-staging/actual-query-preflight-fixed.json`, and `validation/beta/ecoli-report-current-v7/browser-review.json`. |
| CTX-01 | Species/ST discovery, explicit LIN/HC refinement, exact attrition/coverage | verified for captured ST258 snapshot | Previously accepted live ST258 context evidence remains scoped; it does not imply universal context completeness. Automatic context/search/analysis/display-depth selection is not an explicit beta gate in the governing specification. |
| DIST-01 | Default categorical allele and explicit direct gene-DNA distances, coverage and units agree with independent arithmetic | verified for frozen20 panel | Prior 190-pair frozen-panel arithmetic check remains scoped; no new distance acceptance is claimed. |
| REL-01 | Broad context groups never labelled as proof of close relatedness | partial | Broad context, finer LIN prefixes and distance-ranked neighbours remain separately described. No universal grouping or production cutoff is claimed. |
| REL-02 | Useful finer grouping retains actual distance-ranked nearest ties | verified for captured ST258 queries | Previously audited captured-ST258 minimum/top-three tie retention remains scoped to its provider/prefix search space; no global nearest-genome claim. |
| REL-03 | Missing loci and cross-prefix neighbours remain explicit | partial | Cross-prefix neighbours and missing-locus differences remain visible; compatible reference snapshot and assignment mechanism remain unresolved. |
| RES-01 | Dynamic thresholds evaluated with independent labels/held-out evidence | experimental | Held-out-date research remains transductive and does not validate a routine cutoff. Independent epidemiological labels and agreed error tradeoffs remain absent. RES-01 is research-nonblocking for the core beta. |
| SUB-01 | Inputs/pins/mandatory ties retained, deterministic alternatives and strata reasons | verified for frozen 97-profile selection fixture | Five distinct reproducible 63-tip sets retain all 12 inputs and mandatory ties. Independent raw categorical rankings across all 85 supplied contexts found no nearest-three omissions relative to the 83-context LIN5 block. The two LIN5-unmatched contexts remain in the supplied 97-profile pool; the normal block contains 95 profiles. `validation/beta/positive-subsampling-plan/independent-selection-review-v7.json` |
| SUB-02 | Shared-anchor sensitivity measured without calling it CI/convergence | verified for fixed-target and shared-topology comparison; biological limits remain | Independent comparison covers all five strict 63-tip runs. The fixed 13-ID unrooted RF distance is 0 for all ten pairs; nearest-tie ID sets are stable although alignment hashes differ. Fixed 12-target MRCA point estimates span 1.233204 years (1995.130927–1996.364131); this is a descriptive point range, not a confidence interval or proof of convergence. `validation/beta/positive-subsampling-plan/independent-five-comparison-v7/comparison.json` |
| TREE-01 | Exact selected identities, checksummed acquisition, no replacements | verified for named frozen cohorts and exact V7 package identities | The V7 frozen101 packaging independently verified all101 source UUID identity joins and unchanged native call tables. It reused saved V6 typing outputs; it is not a fresh V7 native101 run or new downstream analysis. Original wrapper failure remains preserved. |
| TREE-02 | Alignment/recombination/callable-site outputs agree with benchmark expectations | partial | Five strict trees completed on distinct 63-tip selections; the shared 13-tip unrooted topology has RF=0 pairwise. Alignment/mask hashes differ, so this is not full-tree or biological equivalence. Historical-mask and independent biological benchmark remain open. `validation/beta/positive-subsampling-plan/independent-five-comparison-v7/comparison.json` |
| TREE-03 | Optional rapid topology support, labelled source/policy | verified for named UFBoot run; broad calibration unclaimed | Optional support V9 ran on the exact eight-record fixture with 1,000 UFBoot--bnni replicates (491.1 s, 980 MB peak RSS, 235 MB output). All five corrected split identities were listed; four were present with reported split-support values of 90–100%, and one absent split retains null support; this is not a split-coverage metric. The independent review verified unchanged primary tree/alignment pins. Temporal result is not assessed; this scoped feature evaluation is not broad support calibration. See `validation/beta/optional-support-v9/retry-independent-review.json`. |
| TEMP-01 | Supported/unsupported/not-assessed/failed distinct; complete clustered permutation screen | verified for five fixed positive sets and scoped negative controls | All five positive V7 runs returned supported after 100 date permutations (p=1/101); the V7 negative platform fixture returned normal unsupported p=0.376 without dating override. These remain named-fixture results, not universal temporal validation. `validation/beta/positive-subsampling-plan/tree-five-terminal-review-v7.json`; `validation/beta/linux-aarch64-acceptance-v7-independent-review/terminal-review-v7.json` |
| TIME-01 | Saved topology retained, correct dates/units and 90% marginal node intervals | partial; fixed-target outputs independently reconciled | V8 time-compare review reconciles all 12 fixed input-target estimates/bounds and full 63-tip descendant identities across five V7 runs, plus unchanged 51-anchor fields. Matched-clock review checked all five exact panels and regression arithmetic; reported slopes are descriptive, not temporal-signal tests. Point estimates span 1.233204 years; the range is not a confidence interval. Per-run 90% max-posterior regions remain conditional within-run. `validation/beta/positive-subsampling-plan/time-comparison-execution-v8/independent-review.json`; `validation/beta/positive-subsampling-plan/matched-clocks-independent-review-v8/acceptance-receipt.json` |
| TIME-02 | No-signal dataset returns useful unsupported result without forced dates | verified for frozen negative controls | Earlier negative controls, including the supplied V6 ARM64 fixture, returned normal unsupported outcomes without forced dating; the V7 negative fixture complements, but does not replace, the five positive-set comparisons. |
| REP-01 | Existing style/network preserved; offline images/views and useful navigation | partial | V6 browser checks and V7 external-wheel asset checks remain scoped and do not establish WCAG or human acceptance. The earlier E. coli report review identified two caption defects; V8/V9 fixes are integrated and tested, and a saved-point preview confirms concise units and a neutral tree caption. Fresh regeneration against the current main source remains separate. Two independent human sessions remain open. validation/beta/report-navigation-accessibility-v6.json; validation/beta/install-v7/acceptance-evidence.json; validation/beta/ecoli-report-current-v7/browser-review.json |
| REP-02 | Every module has explained files, exports and README | partial | Prior stage guides and scoped wheel outputs remain evidence. The V7 supported Pixi terminal route passes four ordinary stages for the supplied pretyped profiles; two clean-user walkthroughs remain open. validation/beta/install-v7/acceptance-evidence.json |
| EXP-01 | GrapeTree profiles/metadata import without ID/missing-call corruption | verified | Actual published 95-profile native GrapeTree export imported in local GrapeTree 3.0.0. Independent reconciliation verified all 59,755 profile cells across 629 loci, 168 missing calls encoded as zero, and exact 95-tip tree/metadata/profile ID agreement. This is a target-app import for this dataset, not a Linux or human usability test; no equivalence with ChronoClade NJ is claimed. `validation/beta/published-exports-v7/identity-value-reconciliation.json` and `validation/beta/published-exports-v7/verify_published_exports.py`. |
| EXP-02 | Microreact metadata/tree import with exact shared IDs | verified | Actual published 63-tip dated tree and metadata imported in the hosted Microreact upload viewer. All 63 IDs mapped and topology plus branch lengths matched the ChronoClade dated tree; No Save Project action was performed. This verifies the imported upload only, not human usability or Linux support. `validation/beta/published-exports-v7/identity-value-reconciliation.json` and `validation/beta/published-exports-v7/verify_published_exports.py`. |
| PERF-01 | Date diagnostics do not repeat full profiles/provider payload per fitted point | verified | Prior 7,815-profile diagnostic payload benchmark remains scoped; V7 does not add a performance claim. |
| PERF-02 | Full large profile pool retained within measured resource envelope | verified for frozen full CLI output with separate read-only audit | Prior frozen V4 large-profile run and independent audit remain scoped; V7 source checks and selector run do not replace a broader performance gate. |
| REPRO-01 | Interrupt/resume and changed-source rejection validated per module | partial | Prior V4/V5 recovery evidence remains scoped. The V7 frozen101 packaging has a preserved original wrapper guard stop and a repaired ordinary import receipt; this does not establish scheduler/network recovery. |
| PUB-01 | Published biological targets/data/checksums declared before comparison | verified for frozen benchmark scopes | Previously frozen Bowers/E. coli targets remain the evidence scope; no publication was authorised. |
| DEMO-01 | Visual presentation uses actual outputs and labels unsupported/unverified claims | partial | Five-set tree/time comparison is independently reviewed within its fixed-target scope; published target-app exports are reconciled. The historical E. coli caption defects are fixed and tested in V8/V9, with a saved-point preview; fresh regeneration against the current main source is not established by that preview. V10 and V11 full software checks passed on their separately frozen sources. Two human sessions, E. coli real-context union validation/native HierCC, Linux x86_64/HPC (if claimed) and biological benchmark gates remain open. |

## Work sequence

1. Freeze this brief and audit the existing code/evidence without deleting old data.
2. Investigate ST258 useful grouping with finer LIN prefixes and actual neighbours.
3. Verify typing-reference readiness and build a reproducible concordance check.
4. Exercise a small four-module fixture, exports, errors/resume and unsupported dates.
5. Fix blockers: diagnostic payload bloat, interpretation/close-relative scope,
   missing native setup/handoffs and any failed invariants.
6. Run independent biological benchmarks and large-block resource measurements.
7. Produce a real-output demonstration and human usability package.
8. Report unresolved release blockers; do not silently mark external/human work done.

## Benchmark policy

- Small frozen fixture: deterministic profiles with known called/missing values,
  input/context roles, LIN assignments, partial dates and geographic/host metadata.
- Typing agreement: exact same assembly hashes and pinned upstream references,
  easy and difficult samples, locus concordance and assignment provenance.
- Kp biological benchmark: freeze chosen paper accessions, available date/location/
  host labels, expected clades and qualitative findings before analysis. Expected
  findings cannot be invented from the ChronoClade result itself.
- Positive and negative temporal datasets: frozen before fitting or filtering;
  no removal of dates because they weaken a desired fit.
- Large performance: retained7815-profile ST258 or synthetic10k; record count,
  wall time, peak RSS, disk and complete tip/coordinate identity. Large-SNP stages
  use explicitly selected subsets, not claimed whole-pool inference.
- Historical ST258 native execution used diagnostic coherence policy and20
  permutations, so it is not by itself a standard-policy beta pass or evidence
  of precise dates. These caveats must remain visible.

## Owners and artifact locations

- Parent: specification, acceptance matrix, tree/time audit, diagnostics and integration.
- Grouping agent: ST258_GROUPING_AUDIT.md and validation/beta/st258_grouping*.
- Typing agent: TYPING_AUDIT.md and standalone typing concordance utility/tests.
- Usability/export agent: USABILITY_EXPORT_AUDIT.md and interoperability utility/tests.
- All raw beta check logs and machine-readable results: validation/beta/.
- Frozen large evidence remains under tmp/st258-native-e2e; do not duplicate matrices.

## Wave 1 evidence: 17 October 2026

These are local engineering checks. They do not establish readiness for a beta
release, a validated outbreak cutoff or native typing concordance.

| Work | Observed result | Remaining acceptance work |
| --- | --- | --- |
| ST258 relatedness audit | Depth 5–10 counts; raw-nearest table reproduced exactly (842 pairs); broad public pool 7,774 versus LIN8 input-associated union 758 | Explain reference/code discordance; validate scope choices on independent labelled data |
| Explicit LIN8 report pilot | 41 inputs + 758 public genomes, 12 groups; all 842 nearest pairs and211 mandatory top-three neighbours retained; five distinct 302-genome selections;40.96s local run | Review usability and decide explicit scope/default integration; no SNP/time rerun or new biological cutoff claim |
| Shared display depth | Validated optional level 5–10 API, unchanged default 7; 14 lineage/tree/refresh checks passed | No production CLI/default change yet |
| Diagnostic evidence size | New separate root-comparison artifact 20.4 MB versus 783.6 MB,97.4% smaller; unchanged source tree and fit statistics | End-to-end regression/resources across the remaining modules and clean beta runs |
| Stage output guides | Shared prepare/tree/time root READMEs with handoff/recovery descriptions; 73 passed, 1 skipped | Two human walkthroughs, clean installs and user navigation |
| Viewer exports | Local export command and documentation, source/output checksums and ID maps | Published 95-profile GrapeTree and 63-tip Microreact uploads imported and reconciled; no Linux target-app or human usability acceptance. See [export audit](USABILITY_EXPORT_AUDIT.md) and `validation/beta/grapetree-import-validation.md`. |

Evidence:

- [Grouping audit and reproduction commands](ST258_GROUPING_AUDIT.md).
- `validation/beta/st258_grouping_results.json` and `st258_grouping_test.py`.
- `tmp/beta-st258-relatedness/scope_acceptance.json` and its report/derived bundle.
- `validation/beta/diagnostic-payload-benchmark.log` and `benchmark_diagnostic_payload.py`.
- `validation/beta/stage-guide-tests.txt`.

Historical source runs remain unchanged. No GitHub publication is authorised.

## First-wave integration verification

The stable source snapshot passed **853 tests, 2 skipped** in 288.58 seconds
(`validation/beta/stable-tests.txt`). Ruff and the strict documentation build
passed. Earlier logs retain failures observed while fixtures/code were being
updated; the stable rerun is the integration result. The bundled synthetic
example under `examples/beta/` also completed prepare → cgmlst → export using
the isolated installed wheel. This verifies engineering handoffs, not native
biological concordance or release readiness.

The original ST258 artifacts remain intact. All changes are local and unpushed.

## Pasteur reference and native typing follow-up

The authorised reference builds completed locally on 17 October 2026. Klebsiella
`query_config.json` reports ready; both references contain 629 loci in identical
order, version 2026-10-05. The LIN lookup contains 839,874 allele records and passes
SQLite quick-check. Caller/LIN builds took 323/365 seconds respectively. No
reference data or credentials were published.

`authorize-pasteur` now provides a reusable browser-consent/hidden-verifier flow;
setup accepts the resulting cached token for authenticated caller and LIN data.
The full software suite passed 873 tests with 2 skips; 59 focused setup/auth/query
tests passed after the final next-command error wording adjustment. Ruff and
strict documentation checks passed.

The eight-assembly native run and sequence-resolved observations are described
in `validation/beta/native-typing-inputs/NATIVE_RESULTS.md`. These are partial
biological concordance evidence, not complete historical reference validation.

## E. coli / EnteroBase expansion

Track the additional user request through the separate gates in
[ENTEROBASE_SUPPORT.md](ENTEROBASE_SUPPORT.md): public reference installation,
known-profile native concordance, HierCC assignment, normal public-context/report
workflow and preregistered published comparison. Other public schemes reuse the
reference installer but require organism-specific validation before advertising
full support. Frozen index inventory: `validation/beta/enterobase/`.

As of 17 October, catalogue parser/download integrity and CLI fixtures pass
12 tests; six public indexes are accessible, two return 403. E. coli native
caller indexing and published input freezing are complete. Ordinary prepare
imports a single provider-produced HClink assignment, but a compatible indexed
cgMLST-profile-to-HierCC mapping and native assignment for new assemblies remain
unverified. The probed `hierBG` files are inaccessible and are not assumed to be
required for a prepared HClink reference.
