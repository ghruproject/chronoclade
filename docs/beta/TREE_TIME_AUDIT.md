# Tree/time beta audit

> **Current status (18 October 2026):** The findings and gap table below are the historical initial audit, dated 17 October 2026, and are not a current completion checklist. The current V11 software suite passed (1,230 passed, 2 skipped; Ruff and strict docs passed). The previously listed tree/time software items have scoped evidence: five V7 positive selections and the V7 negative control are terminal; V9 optional UFBoot support has an independent scoped review. See the [current acceptance plan](ACCEPTANCE_PLAN.md) and V11 requirement audit at `validation/beta/REQUIREMENT_COMPLETION_AUDIT_V11.md`; terminal evidence includes `validation/beta/positive-subsampling-plan/tree-five-terminal-review-v7.json`, `validation/beta/positive-subsampling-plan/independent-five-comparison-v7/acceptance-receipt.json`, `validation/beta/linux-aarch64-acceptance-v7-independent-review/terminal-review-v7.json`, and `validation/beta/optional-support-v9/retry-independent-review.json`. Historical Kp benchmark equivalence and broader uncertainty validation remain open.

Date: 17 October 2026. This is a current source/evidence audit, not a fresh
published-data replication. Scope: independent tree/time commands.

## Existing implementation

- `tree_stage/workflow.py` binds exact selections, retrieves assemblies without
  choosing replacements, and records SHA256 acquisition evidence. A common
  highest-N50 input is selected as reference for alternatives; inference seed
  policies are fingerprinted.
- SKA, IQ-TREE and ClonalFrameML run via the native core. Strict raw/filtered
  relative-distance coherence screening is the default. An explicit diagnostic
  policy retains flagged genomes with recorded warnings.
- `tree_stage/temporal.py` identifies maximal monophyletic same-date clades from
  the date-independent tree. It shuffles raw date strings at cluster level and
  refits least-squares rooting for observed/permuted screens. Original partial
  dates are retained. Failed null fits cannot count as successful evidence.
- Fewer than three independent date clusters/distinct dates returns not-assessed.
  Assessment success/unsupported status is checked by result contracts.
- `time_stage/workflow.py` dates only supported saved trees unless the user
  explicitly overrides. The input seed is retained and polytomies preserved;
  exact tip sets and rooted nontrivial internal leafsets are checked.
- Native dating requests `--confidence --time-marginal only-final`; reports label
  90% marginal node-date regions separately from between-selection sensitivity.
  Joint fitting also requests covariation; its model assumptions need assessment
  on independent positive/negative biological benchmark datasets.

## Rechecked historical evidence

`validation/beta/historical-native-inventory.json` inventories the retained
ST258 tree-final and time-final manifests: five completed tree runs and five
completed dated runs, each with302 selected genomes. All tree assessments say
supported. Each time manifest records topology preservation.

This does not constitute new biological validation. The tree ensemble used
`coherence_policy=diagnostic`,20permutations, and a previous input/context
selection. The canonical beta path defaults to strict screening, so the existence
of these outputs cannot be used as a pass for strict execution or precise dates.
Full artifact byte validation and an independent Kp benchmark remain outstanding.
The confidence/date intervals must be inspected, not replaced by point-estimate
agreement across highly overlapping sets.

Supplemental historical studies exist under validation/st239_baines2015 and
validation/komori2024_st131. They use older commands and species/inputs distinct
from the primary Kp beta workflow. Their recorded findings are candidates for
rerun/compatibility checks, not substitutes for current Kp acceptance.

## Historical initial gaps and acceptance work (17 October 2026)

| Priority | Gap | Acceptance |
| --- | --- | --- |
| Blocker | Current Kp positive/negative benchmark not independently frozen/reproduced | Freeze paper accessions, expected clades and temporal outcomes before fits; verify source hashes; publish discrepancies |
| Important | Default native tree uses20permutations, giving only1/21minimum p-value | Document quick/beta validation settings; use100or justified larger count for validation; retain complete-null requirement |
| Important | Rapid branch-support option absent from independent tree CLI | Add optional supported IQ-TREE procedure, recorded source tree, seed and resume key; do not conflate support with corrected-tree/date confidence |
| Important | Same-date clustering may be weak for datasets using many unique day dates | Review sampling-design sensitivity using source study labels/precision; do not fabricate common dates |
| Important | Runtime exceptions and not-assessed outcomes need uniform user-facing categories | Four distinct states in summary and JSON; demonstrate a computational failure vs legitimate unsupported dates |
| Important | Fresh native supported/unsupported small fixture evidence incomplete | Run native tests in stage environment; identify stubs separately; inspect node interval and false-positive results |
| Important | Module root README/export entry points not uniform | Shared output guide for prepare/tree/time; actual Microreact/GrapeTree imports |
| Deferred | Topology/model/root uncertainty not integrated into date intervals | Explicit limitation; compare selections and alternative justified settings without relabelling range asCI |

## First corrective implementation

PERF-01: root-to-tip diagnostics copied full profiles/provider payload into each
point and again into baseline, scaled and signed fits. Retained ST258 root
comparison JSON is roughly747MiB and profile_analysis roughly948MiB. A shared
`date_distance_point` helper now selects identity/date/plot fields while keeping
canonical bounds and precision. It is used by root and matched-clock comparisons.
Focused tests compare unchanged fitted statistics, original input preservation,
and bounded evidence size despite deliberately large unused source payloads.
A fresh ST258 diagnostic-only benchmark writes a separate output; original
matrices, NJ trees, selections and historical results are retained.
