# E. coli HierCC report and profile-only partition audit

Date: 17 October 2026. This is an engineering audit using a deliberately
synthetic five-profile fixture, not biological validation. The fixture is in
`validation/beta/ecoli-report-audit/`. Report changes are engineering label and
display fixes with focused regressions; no native caller, public retrieval, or
fitting was run on biological data.

## Findings

The normal `prepare` and `cgmlst` commands correctly preserve and use explicit
E. coli HierCC evidence when the requested HC level is present. With
`--hiercc-level HC10`, the fixture produced two blocks, ST95/HC10:23 (three
inputs) and ST95/HC10:91 (one input). The fifth input has no HC10 assignment and
is recorded as `missing_or_unresolved_partition_typing`; it is not assigned a
neighbour's code. The run is `partial`, and the top-level index says
`HierCC: HC10`. The block reports now colour and summarize profiles by their
reported HC10 assignments. They do not use generic 2% distance-cut groups as
the reader-facing group basis. Distance-cut diagnostics remain available as
underlying exploratory evidence.

Two reader-label defects in the explicit-HC path were corrected with focused
regression tests:

1. The block index now says `Analysis block HierCC level: HC10.` when an HC
   partition is selected; it uses the cgLIN depth caption only for the default
   lineage partition and detects `partition_by="st"` from partition evidence.
2. Distance-group notes no longer assume LIN blocks. HierCC diagnostics say
   they do not assign or replace HC codes, and ST-only reports explicitly state
   that no lineage/context refinement is included.

Separately, HierCC block reports now use the reported HC assignment at the
partition's explicit level for NJ/ordination colours and the group summary.
Generic distance-cut groups remain in the optional sensitivity/evidence
section, where the cutoff is labelled exploratory and not biologically
validated. ST-only reports use one neutral unassigned-lineage display category.

Profiles without any usable LIN/HierCC assignment cannot be analysed by the
lineage-default cgMLST path, even when every input has the same species, MLST
ST and compatible cgMLST profile schema. With no `--hiercc-level`, the
command uses its default cgLIN partition; with `--hiercc-level HC25`, it asks
for an HC level absent from every fixture profile. Both runs finish by writing
an `unavailable` manifest with zero blocks and list all five inputs as
unresolved. This is not an error in the native profile caller; it is the lineage-partition
contract. The opt-in ST profile-only path now supports relationship analysis
when HierCC has not been supplied or assigned.

With the lineage-default path, no automatic group is inferred from profile
distances in either no-lineage case. The opt-in `--partition-by st` path now
analyses the supplied profiles in a species/ST/profile-scope block while
keeping the missing HC explicit and displaying one neutral unassigned-lineage
category.

## Reproduction

From the repository root, run the following commands. The fixture catalogue
has four synthetic loci and is not an EnteroBase catalogue. `--no-enrich-metadata`
keeps this run local and deterministic.

```bash
pixi run chronoclade prepare \
  validation/beta/ecoli-report-audit/synthetic_profiles.json \
  --catalogues validation/beta/ecoli-report-audit/synthetic_catalogues.json \
  --out tmp/ecoli-report-audit-repro/prepared --no-enrich-metadata

pixi run chronoclade cgmlst tmp/ecoli-report-audit-repro/prepared/dataset.json \
  --out tmp/ecoli-report-audit-repro/hc10 --hiercc-level HC10 \
  --bootstrap-replicates 0 --subsamples 1 --min-overlap 0.5

pixi run chronoclade cgmlst tmp/ecoli-report-audit-repro/prepared/dataset.json \
  --out tmp/ecoli-report-audit-repro/default --bootstrap-replicates 0 \
  --subsamples 1 --min-overlap 0.5

pixi run chronoclade cgmlst tmp/ecoli-report-audit-repro/prepared/dataset.json \
  --out tmp/ecoli-report-audit-repro/hc25 --hiercc-level HC25 \
  --bootstrap-replicates 0 --subsamples 1 --min-overlap 0.5

pixi run chronoclade cgmlst tmp/ecoli-report-audit-repro/prepared/dataset.json \
  --out tmp/ecoli-report-audit-repro/st --partition-by st \
  --bootstrap-replicates 0 --subsamples 1 --min-overlap 0.5
```

Observed output: `hc10/cgmlst.json` is `partial`, with two blocks and one
unresolved input (`eco5`); `default/cgmlst.json` and `hc25/cgmlst.json` are
`unavailable`, with no blocks and five unresolved inputs each. The opt-in ST
run is `complete` with one species/ST/profile-scope block and no context. Its
index and block report say “Supplied inputs only; no public context retrieval
or lineage refinement”; the block report shows one “Lineage unassigned” display
category. For HC10, the reports display the actual reported HierCC groups; the
distance-cut sensitivity remains under optional diagnostics. The prepared manifest contains a creation timestamp, so its file hash changes
on each run; source fixture hashes below are stable.

Fixture checksums:

| File | SHA-256 |
| --- | --- |
| `synthetic_profiles.json` | `89c17066e776b4674307a145a6ee0e59425aaca321b6986266a7858ce3627803` |
| `synthetic_catalogues.json` | `70e90c43bdca6141a4b59183766502aa6a556d60b5c0b5560cfdf91266bdac26` |

## Approved profile-only ST mode; implementation verified on fixture

Parent approved an explicit `--partition-by lineage|st` option, keeping
`lineage` as the default. `st` is a supplied-profiles-only analysis partition,
keyed by normalized species, declared MLST scheme/ST namespace, and compatible
cgMLST scheme/version/database fingerprint. Missing or conflicting HierCC
remains visible as input metadata and does not create labels or block
membership. The analysis uses observed cgMLST profile distances for
relationships. In ST-only mode, one neutral “Lineage unassigned” category is
used for display; distance-cutoff groups remain underlying exploratory
evidence, not the reader-facing group basis.

`st` mode rejects `--fetch-context` and all prepared context records. The
report states that it analyses supplied input profiles partitioned by
species/ST; no lineage-compatible context is implied. It retains explicit
missing-HC statuses. Any generic distance groups must not be renamed as
HierCC/cgLIN. It displays one neutral unassigned-lineage category, not
cutoff-derived groups. This is a profile relationship analysis, not evidence
of lineage assignment or a substitute for a later lineage-scoped comparison.

The core implementation binds `partition_by` in `cgmlst.json` and resume
identity, records membership and exclusions, rejects prepared context in ST
mode, and prevents public-context selection claims. Report captions identify
ST and HierCC partitions correctly. HierCC reports use the block's explicit HC
level and reported assignments for their reader-facing group display; ST mode
uses one neutral unassigned-lineage category. Focused unit tests and the normal
CLI fixture run verify these behaviors. This does not establish biological
validity of either grouping or native caller output.

## Scope and readiness

This synthetic check verifies ordinary command plumbing and report labelling
only. It does not validate native cgMLST allele calls, HierCC assignments,
biological cluster boundaries, or the published E. coli benchmark. Exact
caller concordance for the separate ST131 sample and the frozen ST95 published
panel remain separate acceptance evidence. The ST95 study panel cannot produce
a lineage-scoped cgMLST block when its profiles lack compatible HC evidence;
the opt-in ST mode can support a clearly labelled profile-only analysis, not
turn its published target into a positive temporal result.
