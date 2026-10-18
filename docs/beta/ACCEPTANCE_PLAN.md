# Klebsiella beta acceptance

Updated 18 October 2026. This page records the current decision; earlier source
snapshots and their separate evidence remain in [acceptance history](ACCEPTANCE_HISTORY.md).

Latest installed-package check: [local beta closeout](LOCAL_BETA_CLOSEOUT.md).
The current wheel completed all four stages on eight supplied assemblies with
fresh Python dependencies, using the locked native tools. The optional per-locus
SNP route has its own bundled synthetic example and six verified pairwise results.
This adds current package evidence without replacing the separately scoped
native-typing, public-context and supported-dating benchmarks below.

## Decision and scope

The Apple Silicon Klebsiella workflow is ready for exploratory scientific use,
subject to the checks and limitations below. It is not a validated outbreak
classifier. The user authorised the final documentation, worked example and
independent agent walkthroughs, followed by a GitHub push if those checks pass.
A push does not constitute a versioned release or a claim of broad platform support.

Human recruitment is unavailable. Two independent agent walkthroughs replace the
human release gate for this beta: one follows the documented commands, the other
reads the reports without developer explanations. Record this as simulated
usability evidence. Real-user feedback becomes a beta follow-up.

## Current checklist

| Area | Status | Evidence and limits |
| --- | --- | --- |
| Mac installation and four-module handoffs | Passed for recorded small native run | Two Greek ST258 assemblies → 234 profiles → 28 selected assemblies → unsupported time result. Fresh Pixi environment; configured typing references and provider/assembly caches reused. About 25 minutes. |
| Native cgMLST/cgLIN typing | Passed within eight-assembly comparison | 4,998 direct allele-ID matches; all 14 changed identifiers resolve to the same DNA sequence; 12 jointly unusable calls. Eight reported cgLIN codes agree. Historical provider reference equivalence remains incomplete. |
| Context and close-relative reporting | Passed for recorded datasets | Species/ST discovery, explicit finer LIN context scope, separate analysis/display depths, measured distances and coverage. LIN prefixes are not validated epidemiological cutoffs. |
| Representative selection | Passed for fixed benchmark | Five distinct 63-tip sets retain all inputs and mandatory neighbours. Shared 13-tip topology is stable; fixed-target MRCA point estimates span 1.23 years. This range is not a confidence interval or proof of convergence. |
| SNP trees and temporal gate | Passed for named positive/negative examples | Strict recombination-adjusted trees; 100 clustered date permutations. Greek two-input run returns unsupported, p=0.0891; five positive selections pass their saved screens. |
| Supported time output | Passed for recorded 63-genome rerun | V16 preserved all IDs and rooted topology; MRCA 1995.6604, conditional 90% max-posterior region 1994.268418–1998.078415. Saved SNP tree reused. |
| Reports and exports | Passed for inspected examples | Existing style, embedded views, readable labels, table controls and upstream links. GrapeTree/Microreact imports reconciled on named datasets. Country/host/region/city networks share one engine; missing fields omitted. |
| Software checks | Passed before documentation closeout | 1,260 tests passed, 2 skipped; Ruff, strict documentation and browser network/export checks passed. Final closeout changes receive relevant checks only. |
| Simulated command walkthrough | Passed within simulated scope | Independent agent, README and documented fixture only. Fresh output directory, with install-cache reuse disclosed. |
| Simulated report walkthrough | Passed within simulated scope | Independent reader identifies dataset, neighbours, geography, selection and temporal verdict from reports. |
| Worked-example entry page | Passed | One page links the actual four-module negative run and separate supported dating example. Reuses existing evidence; does not rerun inference. |

## Evidence locations

These are local audit records unless included in the compact publication summary.
Do not combine different source snapshots into one purported run.

- Native run: `validation/beta/mac-native-final-plan/execution-v14/`.
- Typing comparison: `validation/beta/native-typing-inputs/NATIVE_RESULTS.md`.
- Five-set comparison: `validation/beta/positive-subsampling-plan/independent-five-comparison-v7/`.
- Supported rerun: `validation/beta/positive-subsampling-plan/time-closeout-selection001-v16/`.
- Report polish: `validation/beta/mac-polish-positive-v16/`.
- Four-state network checks: `tmp/state-network-fields/verification.json`.
- Final simulated reviews and [closeout receipt](evidence/closeout.json): linked from the worked example.

## Deferred or research work

- E. coli/HierCC reference sourcing and validation remain parked.
- Linux x86_64 and physical HPC/SLURM verification await a supplied environment.
  Earlier Linux ARM64 fixture checks do not establish either.
- Dynamic, time-aware cluster cutoffs require independent epidemiological
  benchmarks. No validated automatic cutoff or geographic origin is claimed.
- Historical typing-reference differences, sampling bias and conditional dating
  uncertainty remain limitations of scientific interpretation.
- Real human usability feedback is a follow-up, not a blocker for this beta.

The next scientific work should use frozen hypotheses, labelled cohorts and
held-out comparisons. Passing software checks does not resolve these questions.
