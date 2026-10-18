# Pinned Klebsiella LIN scheme: correction review

The parent reviewed and applied the correction package changes, then derived the
separate real `tmp/typing-setup-smoke/query_config.published.json`. Original
installed scheme/config/database bytes and historical paths remain unchanged.
The isolated implementation and pre-application patch remain in
`validation/beta/plincer-scheme-sandbox/` as audit evidence. Its 28 sandbox tests
and 39 portable helper/setup controls are historical isolated checks; the parent
then ran 76 active focused tests successfully. The version-2 broad regression
was still pending when the native controls below completed.

## Scientific source and exact defects

The official [Pasteur LIN specification](https://bigsdb.pasteur.fr/klebsiella/cgmlst-lincodes/)
was read on 2026-10-17. It identifies the 629-locus scheme and gives these
adjacent identical-allele ranges: bin 6 [619,622), bin 7 [622,625), bin 8
[625,627). Its corresponding percentage boundaries are 98.4102, 98.8871,
99.3641 and 99.6820. These equal the rounded 629-locus fractions, including
622/629×100 = 98.8871 and 625/629×100 = 99.3641, to four decimals.

The pinned upstream file at commit
`c6d3d85e7484aee4d713429b78402b0590ac8fdc` has SHA256
`c52db6d411e8c5e677f4124c290464c3c2e46b9857861e15fc71daac8c35414c`.
It contains two deviations from that published table:

- Bin 6 upper bound 798.8871 exceeds the percentage domain; expected 98.8871.
- Bin 7 upper/bin 8 lower bound 99.3461 differs from published 99.3641.

The second defect is structurally contiguous, so percentage range/contiguity
checks alone cannot identify it. The foundational paper is Hennart et al., 2022,
DOI 10.1093/molbev/msac135; direct DOI retrieval failed during this bounded
review. The scientific evidence here is the official Pasteur table, with its
citation to that paper, rather than an independently retrieved paper supplement.

## Actual numerical effects

`boundary_audit.py` executes the unchanged upstream `calculate_identity` and
`assign_bin` functions, extracted from the pinned source. It retains upstream
four-decimal rounding and tests every possible identical count at every callable
denominator from 599 through 629: 19,065 combinations. This is a boundary audit,
not a claim about observed dataset frequencies or inferred cluster thresholds.

The level-6 correction changes 76 combinations. The second correction changes
17 combinations: precisely four mismatches with 612–628 jointly called loci.
At full 629-locus coverage, it changes none. Examples:

| Identical / called | Upstream identity | Upstream inherited positions | Both published corrections |
| --- | --- | --- | --- |
| 622 / 629 | 98.8871% | 5 | 6 |
| 625 / 629 | 99.3641% | 7 | 7 |
| 624 / 628 | 99.3631% | 7 | 6 |
| 623 / 627 | 99.3620% | 7 | 6 |

These counts describe the prefix inherited from the closest defined reference;
unresolved suffix positions remain stars. The correction does not mint new
source-database LIN identifiers or change upstream closest-profile tie handling.
Missing calls are therefore consequential: a complete-profile-only test cannot
validate the second boundary. All exact outcomes are in `boundary-results.json`.

## Proposed provenance and safety contract

The shared helper is used by setup and typing readiness. It requires exactly
10 ordered named bins, finite bounds within 0–100, no gaps/overlaps, full domain
coverage and the Pasteur scheme-18 identity. Unknown malformed files fail closed.

Only the exact original SHA and its expected one-token anomaly permit the
known correction. Setup preserves downloaded bytes as `scheme.upstream.toml`,
writes a separate working `scheme.published.toml`, and records original/working hashes,
pinned source URL, official threshold source, exact field transformations and
policy in `scheme.published.correction.json`. Readiness rejects malformed references and
known-corrected bytes lacking their exact original/receipt. Query provenance
includes the verified receipt, binding both source and transformation into reuse.

Two isolated review variants are retained:

- `level6-only`: repairs the out-of-range bound and explicitly records
  `published_boundary_discrepancy_remaining`; it must not be described as fully
  scientifically validated.
- `published-kpsc-629`: additionally repairs both occurrences of 99.3461 and
  records matching published numeric boundaries. This still does not imply
  validated caller concordance or reference-release equivalence.

The focused patch uses `published-kpsc-629` at setup, following parent review;
no extra user-facing switch is proposed. Readiness rejects the diagnostic-only
level6 variant and any unknown valid-looking scheme that differs from the
published rounded 629-locus boundary vector. The published variant is now available through the separate active descriptor;
the original scheme and original descriptor remain unchanged. Retain historical results under their
original reference hashes; no relabelling of old files.

## Acceptance before applying or declaring typing ready

1. Review official boundaries and choose the published variant explicitly. Run
   malformed/unknown-file, tamper, missing-receipt, exact-byte preservation and
   all-boundary tests. The isolated 28-test suite already exercises readiness
   against the actual installed bad reference and a separately corrected copy.
2. Update the existing setup test's `controlled-fixture` placeholder to a valid
   scheme-18 TOML fixture. Any provenance tests constructing fake scheme files
   likewise need a structurally valid fixture; do not weaken validation to keep
   invalid fixtures passing.
3. Replay the eight known native Klebsiella profiles through the real native
   assignment worker under the original and separate corrected schemes, retaining
   full best-match lists, identical/callable counts, prefixes and outputs. Known
   exact-reference matches should remain unchanged. Explicitly report changed
   incomplete-prefix assignments and their triggering rounded identities.
4. Add synthetic near-reference profiles at seven mismatches and four mismatches
   with 628 and 629 callable loci; freeze their expected inherited prefix lengths
   independently of output under test. Include missing-data ties and verify the
   shared-prefix rule remains unchanged. Neither reference ordering nor ties may
   be relaxed.
5. Re-run the ordinary native-public bridge and whole-ST CLI after any parent
   memoization equivalence acceptance. All input/public profiles must use the
   same corrected working scheme/reference provenance before LIN8 scope
   filtering. An old source/reference-bound checkpoint must not silently resume.

This correction does not address whole-ST throughput by itself. The unchanged
30-minute attempt remains a failed end-to-end acceptance item, with all acquired
source files preserved.

## Historical-path preservation and safe migration

The revised patch never writes `databases/plincer/scheme.toml`. Its bytes, path
and modification time remain unchanged. Both new-root setup and generated query
configuration use `scheme.published.toml`; the correction receipt has the distinct
name `scheme.published.correction.json`. The immutable source archive remains
`scheme.upstream.toml`. Unknown/pre-existing conflicting working files fail rather
than being replaced.

The parent reviewed/applied the patch and ran this exact local command to derive
a separate descriptor from the existing installation:

```sh
.pixi/envs/default/bin/python -m chronoclade.plincer_scheme \
  --config tmp/typing-setup-smoke/query_config.json \
  --out tmp/typing-setup-smoke/query_config.published.json
```

This does not call APIs, authenticate, install tools, re-index alleles or rebuild
profiles. It creates only the small published scheme, upstream archive, receipt
and new configuration. Existing profile, allele, caller-index and tool paths are
reused. Relative paths are resolved against the source descriptor so relocating
the new descriptor cannot silently change their meaning. Existing output
descriptors are refused; the original config is never overwritten. The new
config records its source-config SHA256. Historical run descriptors retain their
original bytes/path, and new corrected assignments have distinct reference and
execution provenance.

The command was exercised against isolated small copies of the actual scheme
and descriptor, retaining actual unchanged large-resource paths; receipt
`migration-command-control.json` confirms exit 0 and unchanged active scheme and
config SHA256. The parent subsequently derived the real separate descriptor and independently
verified all 637 immutable resources unchanged. Migration means a new descriptor
and new working file; it does not replace the historical scheme/config paths.

## Separate pinned missing-locus limitation

The official Pasteur page permits up to 30 missing loci. The pinned upstream
`compare_profiles` defaults to `max_missing_loci=29`: an otherwise matching
629-locus pair with 29 missing calls is accepted, while 30 missing calls is
rejected. Both controls execute the unchanged upstream function and are saved in
`boundary-results.json`. The function can be called with a different argument,
but ChronoClade has not modified or monkeypatched that policy.

Therefore the audit's 599-called cases are mathematical `assign_bin` checks,
not profiles reachable through the default `closest_profiles` comparison. The
17 second-boundary changes use 612–628 called loci (at most 17 missing) and
remain reachable under the pinned policy. Fixing TOML boundaries does not fix
this missing-locus difference, nor establish full published-scheme compatibility.
The correction receipt explicitly records this limitation. A separate policy
decision and real typing controls are required before changing upstream missing
call handling.

The final actual-schema migration control hashes every immutable caller-index
resource plus profiles, allele database, metadata, original scheme and original
configuration before/after derivation and full corrected-config readiness. Hashes,
lengths and modification times remain identical. Transient SQLite `-shm`/`-wal`
files are excluded from immutable-reference invariance: SQLite read-only schema
validation can update shared-memory lock state. Tools and `urlopen` were
monkeypatched to raise in this safety test, proving the derivation/readiness path
did not invoke a build or API. Evidence: `actual-schema-migration-hash-control.json`.

The receipt correction ID is `plincer-kpsc-published-629-boundaries-v1`, covering
both documented discrepancies. There are no active persisted correction receipts
requiring migration from the earlier sandbox-only ID.

## Active version-2 native assignment controls

`validation/beta/plincer-corrected-native-controls.json` records real subprocess
commands, source identity and timings. Complete raw outputs remain in
`tmp/beta-plincer-corrected-controls/`. The source matches freeze v2 SHA256
`2064afdafa8b1294450f303f86b2e06e08aa4aa4989f707cb73350977908e013`
(126 Python files). The old intentionally malformed scheme was passed directly
to the real native worker for scientific comparison, rather than bypassing the
new configuration-readiness validator. No caller or database rebuild occurred.

All **8/8 native profiles produced identical complete output objects** under the
original and published schemes, including complete best-match list order. Two
profiles had genuine two-way and eleven-way best-match ties; their lists and
results were preserved. Each real-reference batch took approximately 4.5 seconds.
These eight known profiles all matched at 100% identity, so this observation
alone does not test the repaired boundaries.

Separate controlled perturbations froze one complete actual native reference
profile in a one-reference background and executed the same worker. Seven
mismatches across 629 called loci changed the inherited prefix from five to six
positions as expected. Four mismatches across 629 called loci retained seven
positions; four mismatches across 628 called loci changed seven to six positions.
Best-match lists remained identical between schemes. The old/new 29-missing
controls were accepted at 100% callable-locus identity; 30-missing controls
produced no matches under both schemes. These are controlled reference-profile
perturbations, not a full-catalogue nearest-neighbour or new-assembly validation.

The guarded ordinary full-ST replay script is prepared at
`validation/beta/native-wholest-published-replay.py`. It has not been launched.
It validates the parent-sealed original 13,421-record raw capture at
`tmp/beta-native-wholest-raw-source/provider-source.json` using the shared
`ProviderSource`, then invokes normal `cgmlst --context-source` with the published
descriptor, LIN8 scope and scale backend. It has 90-minute / 12-GiB sampled
process-family RSS / 4-GiB new allocated-output guards, source-v2 identity and
preflight memory/disk checks. No normalized context cache, manual candidate
filtering or context assembly download is involved. This still requires parent
approval after other native/regression jobs complete.
