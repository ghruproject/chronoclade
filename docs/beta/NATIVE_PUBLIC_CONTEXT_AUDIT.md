# Native queries and public Pathogenwatch context

17 October 2026. **Implemented canonical comparison and bounded acceptance.**
The original eight-assembly native validation and native-only four-module run
were insufficient. The ordinary prepare/cgMLST CLI now also completes a mixed
three-native-query/five-real-public-profile benchmark using the downloaded
references. Live whole-ST throughput remains a separate acceptance gate.

## Implemented canonical comparison

`prepare --typing-config` retains a safe configuration path/SHA256 descriptor.
`cgmlst --fetch-context` reuses that verified descriptor; the same ordinary
`--typing-config` option can point to local references after moving a dataset.
No hand-written certificate, scheme-label rewrite or downloaded context assembly
is required.

The implementation validates the installed caller's Pasteur host, seqdef database,
scheme 18 and complete 629-locus metadata. It translates called numeric alleles
through the caller's full SQLite exact-hit DNA-hash index (`reverse=0`). Its
indexed FASTAs are only BLAST representatives and cannot be used as full allele
catalogues. Novel lower-case DNA SHA1 identifiers retain their sequence identity;
missing and ambiguous calls remain unavailable. Unknown/deleted numeric IDs are
excluded explicitly with per-locus evidence, not counted as new mutations.
Conflicting sequences for the same locus/ID cause failure.

The derivative comparison scope is an explicit canonical Pasteur locus/DNA-hash
protocol; it does not pretend the historical database snapshots were equal.
Original native/public calls, job IDs and reference evidence are saved in
normalization evidence under the prepared sources. The original input bundle
also remains intact. Actual current cgLIN assignments are regenerated with the
configured upstream plincer reference, independently of allele normalization.
The batch wrapper loads reference profiles once and uses upstream classifier
functions; it does not introduce another clustering method.

Crucially, this route retrieves the accessible ST candidate profiles **before**
current-reference LIN filtering. Unknown old LIN reference snapshots cannot
silently exclude a genome that would join the requested group after reassignment.
Every profile-unavailable or out-of-scope reassignment is counted with source IDs.

Reproduce the bounded real-data CLI acceptance:

```sh
.pixi/envs/default/bin/python validation/beta/native-typing-inputs/run_native_public_bridge.py
```

Only the external catalogue transport is replayed from five frozen real public
records. Production prepare, automatic reference resolution, normalization,
upstream plincer reassignment, cgMLST comparisons, selection and reports run
normally. Results: **3 native inputs + 5 public comparisons; one complete
629-locus matrix; all 28 pairwise mismatch counts, callable denominators and
distances exactly equal the independently frozen all-native control; all eight
cgLIN codes agree; no context assemblies requested.** Evidence is
`validation/beta/native-typing-inputs/native-public-bridge-results.json` and
`tmp/beta-native-public-bridge-v2/`. This is not claimed as a fresh live API run
or a test of all future Pasteur correction histories.

Focused normalization/provider/prepare checks: **51 passed** in
`validation/beta/native-typing-inputs/normalization-unit-tests.txt`. Controls cover
hash/number synonyms, missing calls, unknown numeric IDs, inconsistent ID/DNA,
unverified source schemes, mismatched locus panels, configuration changes and
old-LIN-prefilter exclusion. Broader release checks remain in the main plan.

The original diagnostic below remains valid for raw incompatible representations
when no normalization reference is supplied.

## Reproduction against the current program

Run from the worktree root:

```sh
.pixi/envs/default/bin/python validation/beta/native-typing-inputs/audit_native_public_context.py
```

The reproducer uses the same exact eight source identities and hashes as the
native validation. Three Greek query assemblies retain their native typing; the
five comparison records retain their historical Pathogenwatch typing. It invokes
the ordinary `prepare` CLI with both explicitly declared catalogues, then inspects
the actual shared context-matching and partition functions. It performs no API
calls, reference changes, retyping, tree inference or scheme-label rewrites. It
also invokes the actual production provider with the five frozen public rows
and no patched functions, checking the explicit refusal and saved audit.
Outputs are in `tmp/beta-native-public-context-audit`; the compact evidence is
`validation/beta/native-typing-inputs/native-public-context-audit.json`.

Observed:

| Check | Result |
| --- | ---: |
| Prepared records | 8 |
| Native Greek inputs | 3 |
| Historical public comparisons | 5 |
| Distinct prepared cgMLST matrices | 2 |
| Compatible native/public cgMLST pairs | 0 of 15 |
| Compatible native/public cgLIN pairs | 0 of 15 |
| Retrieval lineage matches, every LIN5–10 depth | 0 |
| Actual analysis block membership | 3 inputs; 0 context |
| Public records unmatched to an input block | 5 |

Preparation correctly retains separate typed matrices. It does not legitimise
cross-matrix distances. `context_partitions._group_key` includes lineage scope,
reference version and fingerprint; `profile_context_provider._matching` requires
that scope and matching database identities before comparing a prefix. The later
public-profile gate additionally checks cgMLST scope and references. These checks
must remain fail closed.

Before the fix, there was also a usability defect: lineage incompatibility can produce an apparent
zero-match result before the cgMLST compatibility gate is reached. `counts` and
`matched` use `_matching` directly, and the retrieval audit does not separately
count numeric-prefix candidates rejected for scope incompatibility. A report
with zero context can therefore look like absence of public relatives. The normal
live `cgmlst --fetch-context` route calls this exact function. The reproducer
exercises its frozen matching predicates; it does not claim an additional live
API retrieval was performed. The fixed provider now refuses this case with
`failed_incompatible_typing_scope`, saving five rejected candidate records in
`public-retrieval_retrieval/context_compatibility.json` before publishing a
successful context dataset. A genuinely empty compatible pool still succeeds.

## What is known, and what is not

| Item | Native | Historical provider |
| --- | --- | --- |
| cgMLST namespace | `klebsiella_1` | `pathogenwatch:cgmlst-573-2` |
| Declared scope/version | `2026-10-05` | `cgmlst-573-2` (job identity, not a database release) |
| Canonical loci | 629, explicitly complete | 628 observed, explicitly incomplete |
| Reference checksum | Pinned | Unknown |
| cgLIN scheme | `scgMLST629_S` | `scgMLST629_S` |
| cgLIN scope/version | `2026-10-05` | `lincodes-3390273-2` (job identity) |
| cgLIN reference checksum | Pinned | Unknown |

The native-only `KP1_0515_S` is missing in all eight assemblies. An export with no
observed call at that locus is not evidence that it is outside the original
scheme. Backfilling it as missing would be legitimate only after independently
establishing the provider's canonical locus universe.

All eight literal cgLIN codes agree numerically. Of 5,024 shared-locus
observations, 14 old novel SHA1 identifiers match the DNA of a currently numbered
allele exactly. This supports the particular paired observations; it does not
prove every historical numbered allele has the same sequence or every public
LIN assignment belongs to the current reference's code namespace.

## Primary upstream evidence

The inspected pinned sources are the installed native tools:

- MLST commit `4067e2aecbedc21ca2921fc0e6e149a0972119ff`:
  [hashing and output implementation](https://github.com/pathogenwatch-oss/mlst/blob/4067e2aecbedc21ca2921fc0e6e149a0972119ff/src/mlst.js).
  `hashHit` hashes lower-case, locus-oriented DNA; `formatOutput` also produces
  `raw_code`, a per-locus checksum vector, in addition to numeric/novel allele IDs.
- Typing databases commit `d36e3c97e18bc0dc4d4eaa1136d1f15031f4384d`:
  [scheme configuration](https://github.com/pathogenwatch-oss/typing-databases/blob/d36e3c97e18bc0dc4d4eaa1136d1f15031f4384d/config/schemes.json).
  `klebsiella_1` identifies Pasteur `pubmlst_klebsiella_seqdef`, scheme `18`.
  Its human-readable name says `scgMLST639_S`; the actual downloaded ordered
  metadata contains 629 loci. Do not infer size from that title.
- LIN assigner commit `c6d3d85e7484aee4d713429b78402b0590ac8fdc`:
  [README](https://github.com/pathogenwatch-oss/klebsiella-lincodes/blob/c6d3d85e7484aee4d713429b78402b0590ac8fdc/README.md)
  specifies input from the Pasteur/`klebsiella_1` caller.

These establish the native route's source scheme and sequence-hash policy. They
alone do not resolve the public job ID to an exact historical reference snapshot.
The cached public CSV contains Genome ID, Genome Name, Gene, Allele ID, Start,
End, Contig and Direction; it does not expose the caller's `raw_code` DNA vector
or a per-genome reference version. No credentials are copied into this evidence.

Additional primary sources materially narrow the problem:

- The [Pathogenwatch Klebsiella paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8634497/)
  states that cgMLST uses Institut Pasteur's allele and profile definitions.
- [Pasteur's cgLIN definition](https://bigsdb.web.pasteur.fr/klebsiella/cgmlst-lincodes/)
  identifies the KpSC 629-locus scheme and a single source database for cgLIN
  nomenclature. A common source namespace need not use identical snapshots forever.
- The [BIGSdb curator guide](https://bigsdb.readthedocs.io/en/latest/curator_guide/0030_updating_and_deleting_alleles.html)
  prohibits editing an allele's sequence in place, but documents deletion,
  re-addition and optionally retiring IDs. Generic software documentation alone
  does not prove that every historical Pasteur ID was never reused.

**This is an adapter/evidence gap, not demonstrated biological incompatibility.**
ChronoClade uses execution-job identifiers as public `scheme_version` and an
internal tool shortname as native `scheme_id`. These preserve different provenance
but need not imply different biological nomenclatures. Additive reference updates
do not inherently invalidate distances. A validated canonical source-and-sequence
representation should retain release/checksum data as provenance while using
proved semantic compatibility for comparison. It must not rewrite historical
metadata as though the reference releases were identical.

## Scientific contract used by the canonical bridge

A bridge must be an explicit, validated transformation with its own versioned
provenance, not an alias added to `compatible_typing`.

1. **Identity:** bind both inputs to exact source IDs, checksummed assemblies
   where available and frozen response/file hashes. Never substitute a nearest
   public genome for a newly sequenced query.
2. **Scheme:** establish the authoritative host/database/scheme identity on both
   sides (for example Pasteur seqdef scheme 18), including the complete ordered
   locus definition and any rename/removal history. Observed locus overlap is
   insufficient. Retain unknown provider reference versions as unknown.
3. **Alleles:** prefer matching per-locus DNA digests from both callers. Convert
   numbered IDs through a checksummed locus-specific allele sequence catalogue;
   resolve historical novel SHA1 values by the documented hash policy. Never
   compare numbers from unproved namespaces or count a proven hash/number synonym
   as a mutation. Reject duplicate-ID/different-sequence conflicts. Keep missing
   and ambiguous status separate from a called allele. If the provider only
   supplies historical numbers, establish their stable source nomenclature and
   check observed locus-specific IDs against the current catalogue; obtain
   historical snapshots or correction/retirement evidence for unresolved IDs.
   Byte-identical whole-database references are not inherently required when
   sequence identities and source nomenclature are proved.
4. **Lineages:** public numeric prefixes can be provisional discovery candidates,
   not proof of cross-reference equality. Prefer reassigning complete provider
   profiles with the same pinned native LIN reference after profile normalisation.
   Check locus order, assigner coverage, unresolved/provisional results and code
   changes. Preserve original assignments alongside reassignments.
5. **Distances:** create a new, explicit comparison representation whose complete
   locus universe and digest policy are declared. Preserve original native and
   public matrices/provenance. Report jointly callable loci; do not turn missing
   into an allele, or certify a full universe from an observed union.
6. **Controls:** validate eight exact paired assemblies as a starting panel plus
   unseen assemblies and context. Include novel-to-numbered identity, genuine
   one-locus differences, incompatible scheme/locus maps, stale/reused allele IDs,
   callability disagreements and conflicting lineage references. Evidence must
   demonstrate both acceptance of valid translation and refusal of invalid data.
7. **Normal end-to-end run:** new assembly → native typing → public profile/metadata
   discovery → translation → ranked nearest neighbours → ordinary cgMLST report
   and selections. Verify no context assemblies are downloaded during this stage.
   Compare distances and rankings with an independently generated sequence-based
   shared-locus baseline on a manageable checksummed panel.

Two particularly strong evidence routes are: obtain public per-locus DNA hashes plus
canonical scope metadata, or obtain its historical pinned allele/LIN references
and type queries against those. Independently retyping every public context
assembly would give a homogeneous benchmark, but violates the intended
profile-first context retrieval workflow and cannot replace this acceptance gate.

## Immediate code change and remaining bridge

The implemented immediate fix preserves all scientific checks and adds explicit
retrieval diagnostics for numeric-prefix candidates rejected by lineage or
profile scope. When public candidates exist but none can be compared because
scope/reference compatibility is unproved, stop with an actionable message and
save the incompatibility audit. Do not publish a successful zero-context run
that implies no close public genomes exist. This changes failure
reporting, not biological comparability. The parent approved the narrow shared-code
change. Fourteen provider tests pass, including lineage mismatch, profile mismatch
and genuinely zero compatible-prefix matches. The real frozen eight-record
reproducer confirms the refusal. Those refusal-only controls do not establish the bridge; the separate real CLI
acceptance above supplies that evidence.
