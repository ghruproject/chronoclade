# Published export usability audit

Date: 18 October 2026. This scoped target-app check closes the import portion of EXP-01 and EXP-02 for the named published exports. It does not complete beta usability acceptance.

The actual 95-profile native GrapeTree export was imported into local GrapeTree 3.0.0. It displayed the calculated tree and metadata legend. Independent reconciliation verified 95 matching profile/tree/metadata IDs, all 59,755 cells across 629 loci, and all 168 missing calls represented as zero. The export calculation itself exited 0 in 3.27 seconds. This checks values and identity for this import; it does not establish equivalence between GrapeTree’s MST and ChronoClade’s NJ methods.

The actual 63-tip dated ChronoClade tree and metadata were imported through the hosted Microreact upload viewer. All 63 IDs mapped, and the imported topology and branch lengths matched the source dated tree exactly. A Save Project button was visible, but no Save Project action was performed.

The independent receipt, verifier, input/output artifacts and screenshots are in `validation/beta/published-exports-v7/identity-value-reconciliation.json`. The screenshot hashes are recorded in the receipt. The checks do not establish a Linux target-app run, accessibility compliance or completion of either human walkthrough; those acceptance gates remain open.
