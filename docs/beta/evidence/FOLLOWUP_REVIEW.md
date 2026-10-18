# Follow-up simulated reader review

Reviewed 2026-10-18 after the presentation changes. Independent agent content/link review only: not human validation, browser-interaction QA, or a rerun of scientific inference. No code edits or Git mutations were made.

## Confirmed improvements

- `tmp/beta-worked-example/index.html` clearly walks through the real two-input Klebsiella workflow: preparation, 234-profile cgMLST context, a 28-genome selection, SNP tree, and unsupported dating. Its separate supported example is explicitly labelled **a different 63-genome benchmark**. This prevents the positive example from appearing to rescue dating of the negative selection.
- The index explicitly says **recorded real runs**, existing references/provider/assembly caches reused, and **no new inference for this presentation**. Both refreshed time reports likewise state that their presentation was refreshed from saved results and scientific inference was not rerun. The README says original large data are reached by local symlinks, so this is a local worked example rather than a portable standalone package.
- `tmp/beta-worked-example/unsupported-time/report.html` now has working **Worked example**, **SNP-tree report**, and **Evidence files** links. It still prominently states **28 exact selected genomes; dated tree unavailable**, unsupported dating, and p=0.0891 against threshold 0.05. There are no forced calendar estimates.
- `tmp/beta-worked-example/supported-time/report.html` now presents a visible **Fitted time tree** section with final TreeTime rate **1.02e-06 substitutions/site/year**, whole-tree root **1995.6604**, and conditional **90% max-posterior region 1994.2684–1998.0784**. It explicitly distinguishes this final fit from the screening regression and explains that intervals exclude uncertainty from genome selection, missing sampling, and residual recombination. This resolves the earlier need to download files merely to obtain the final rate/root/interval.
- The shared source `chronoclade/report_components/profile_summary.py`, including `subsample_summary`, now emits concrete tree and time commands and explicitly says to run them from the analysis block output directory with the ChronoClade environment active. It also explains assembly acquisition and the temporal gate. Source inspection confirms the missing-next-command fix exists.

## Link checks

All six worked-example index targets exist:

1. `native/prepare/report.html`
2. `native/cgmlst/index.html`
3. `native/cgmlst/cglin_3981e3a96c24/profile_report.html#subsamples`
4. `tree/report.html`
5. `unsupported-time/report.html`
6. `supported-time/report.html`

No missing local file targets or same-page anchor IDs were found in either refreshed time report. The unsupported upstream links resolve now that the index exists. Download links in the positive report also resolve. This verifies target existence, not browser rendering or export behavior.

## Final narrow recheck: next commands and refreshed cgMLST report

The two snapshot gaps above are now resolved in the delivered worked example:

- `tmp/beta-worked-example/cgmlst-report.html` contains the refreshed **Continue with the selected genomes** section, both copyable commands, the active-environment requirement, and the instruction to run from the analysis block output directory. It also explains selected assembly acquisition and the temporal gate.
- `tmp/beta-worked-example/index.html` now links this refreshed report at `#summary` and `#subsamples`. The command introduction explicitly identifies the working directory as `native/cgmlst/cglin_3981e3a96c24`.
- Rechecked **all 31 local links in the refreshed cgMLST report** and **all six local links in the index**, including fragment IDs in same-page and linked HTML targets: **zero missing file targets or anchors**. Downloads resolve to the original native evidence.

The original saved native cgMLST HTML remains a historical artifact; the worked-example entry point now leads to the refreshed presentation. This final narrow check used report text and local target existence only, without a new browser-interaction review or scientific rerun.

## Reader conclusion

The affected reader tasks pass this follow-up content review. The next tree/time commands are delivered with working-directory guidance, unsupported dating has upstream navigation, and final within-run dating uncertainty is readable without downloading evidence. Cached real execution, the 28-genome negative outcome, and the different 63-genome positive benchmark remain clearly distinguished. No outstanding issues remain within this narrow follow-up scope. Human usability and interactive rendering remain outside this simulated review.
