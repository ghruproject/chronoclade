# Typing new query assemblies

Prepare the input assemblies first, then use the cgMLST module to discover
public comparisons. Species and seven-locus MLST ST identify the candidate
Pathogenwatch pool; compatible cgMLST profiles define the comparisons. LIN and
HierCC labels are optional in default profile mode. The
cgMLST report analyses the retained profiles and selects representatives for
the tree module. Context assemblies are acquired for the selected SNP-tree
analysis. A LIN code alone does not establish how close two genomes are.

The main walkthrough below targets Klebsiella. S. aureus also supports native
cgMLST calling without a lineage assignment, as described below. E. coli HierCC
requires a compatible assignment reference database.

## Classical MLST for public discovery

`prepare` reuses supplied positive numeric STs and attempts missing classical
STs from local input assemblies independently of cgMLST or LIN/HierCC. Install
the pinned classical caller and bundled schemes in their separate environment:

```bash
pixi install --frozen -e mlst
```

The ordinary `prepare` command finds the installed sibling environment; no
manual intermediate typing file is needed. The pinned `mlst` 2.32.2 mappings
are `klebsiella` for K. pneumoniae, `ecoli_achtman_4` for E. coli and `saureus`
for S. aureus. All target Pathogenwatch's primary `mlst` field. The local scheme
named `ecoli` is Pasteur eight-locus MLST and must not be used as the primary
Achtman ST. Explicit alternative-scheme input cannot be silently relabelled for
primary discovery; alternative-scheme searches are not implemented here.

Only complete, exact, unambiguous calls assign an ST. Missing tools, references,
assemblies or unresolved calls stay explicit in the preparation audit. Compatible
local cgMLST profiles remain usable without an ST; public same-ST discovery for
that input waits for a resolved primary ST. See [preparation](prepare.md#live-query-preparation)
for the separate evidence and requirements.

## Install the cgMLST and lineage tools

Install [uv](https://docs.astral.sh/uv/getting-started/installation/) first, then run:

```bash
pixi run chronoclade setup-typing --species klebsiella
```

This installs pinned releases of Pathogenwatch's own tools in an isolated local
directory. `plincer` and `hclink` use separate uv environments. The cgMLST caller
is a Node program, so setup also installs a checksum-verified Node runtime and
its native dependencies. BLAST is provided by ChronoClade's Pixi environment.
The tools do not upload query assemblies.

Software installation does not install reference data. The setup manifest records
which databases are missing, and its query configuration cannot be used until
those resources are ready. The default directory is
`~/.local/share/chronoclade/typing`; use `--destination` for a shared cluster path.

To download and index publicly available cgMLST alleles as well:

```bash
pixi run chronoclade setup-typing --species klebsiella --prepare-db
```

Database downloads and indexing can be substantial. Without cached Pasteur
authorization, Klebsiella cgMLST downloads use the upstream public,
redistributable route. To use authenticated Pasteur cgMLST records, supply the
protected client-key file and first place a valid `klebsiella:access` or
`klebsiella:session` token in the setup's existing private
`credential-cache/plincer.json`. `setup-typing` copies that cache to a separate
private downloader cache and does not start browser authorization or pass account
passwords. The setup report records the route used and the actual downloaded
scheme/version; it does not infer a current cutoff from the downloader's public
mode. Current LIN and HierCC reference builds need their own upstream
credentials; a PathogenWatch API token does not grant this access. Supply
protected files when building the references:

```bash
pixi run chronoclade setup-typing --species klebsiella --prepare-db \
  --pasteur-secrets /protected/pasteur-client-key.txt
pixi run chronoclade setup-typing --species ecoli --prepare-db \
  --enterobase-key-file /protected/ecoli-api.key
```

The installed tools and completed databases can be reused without downloading
them for each query. Install the native tools separately on each operating system;
the completed reference files can be shared with a SLURM cluster. Credentials are not included in the
query configuration or analysis outputs.

## S. aureus cgMLST without a lineage reference

Install the classical MLST environment above, then prepare the native cgMLST
caller and its public reference through the normal setup command:

```bash
pixi run chronoclade setup-typing --species saureus --prepare-db \
  --destination /shared/chronoclade/saureus-typing
```

This uses Pathogenwatch's pinned native `saureus_1` caller and PubMLST S. aureus
scheme 2, with 2,208 declared loci. It does not require or install plincer or
HierCC. The generated query configuration is caller-only (`profiles-only`).
The setup report distinguishes a ready cgMLST reference from lineage assignment,
which is not required and remains unassigned.

Unauthenticated PubMLST access includes records submitted on or before
2024-12-31. Download receipts retain this restriction and the downloader-recorded reference
version. If the server omits `last_updated`, the upstream downloader uses the
retrieval date; that date does not imply current allele completeness. Content
hashes identify the actual reference. A reference can be ready to call those known alleles without including
every current allele; novel or missing calls remain explicit. Reference setup
requires network access and can take substantial time. Preserve the full declared
locus universe rather than restricting it to columns observed in provider exports.

Provide local assemblies and declared species in the ordinary assembly input
CSV; omit ST to exercise automatic classical MLST. Pass the generated
`/shared/chronoclade/saureus-typing/query_config.json` to `prepare --typing-config`,
then use the normal profile-mode cgMLST command. Compatible public comparison
requires checked scheme/locus and allele nomenclature correspondence; neither a
species name nor matching locus counts establishes it alone.

## Call and refine

Use the query configuration written by setup, once its readiness checks pass:

```bash
pixi run --frozen --environment default chronoclade prepare focal.csv \
  --input-kind assemblies \
  --typing-config ~/.local/share/chronoclade/typing/query_config.json \
  --out prepared

pixi run --frozen --environment default chronoclade cgmlst prepared/dataset.json \
  --fetch-context \
  --typing-config ~/.local/share/chronoclade/typing/query_config.json \
  --out cgmlst
```

The input CSV needs the species. A supplied positive numeric classical ST is
reused; a missing ST is called from an available local input assembly using the
separate classical caller above. Native cgMLST typing does not calculate that ST.
Preparation types new cgMLST profiles with the configured caller and reference;
LIN/HierCC assignment remains optional for default profile analysis.
During context discovery, the same reference configuration normalises compatible
public typing where required. Public profiles and metadata are retrieved at this
stage; representative SNP-tree selection is a separate step. See the
[cgMLST module](cgmlst.md) for explicit context scope and selection controls.

`type-queries` remains available when a separate native-typing output is useful:

```bash
pixi run --frozen --environment default chronoclade type-queries focal.csv \
  --typing-config ~/.local/share/chronoclade/typing/query_config.json \
  --output query_typing
```

Alternatively, use `--query-typing query_typing/query_typing.json` to import
existing results whose sample identifiers, species and assembly SHA256 hashes
match the focal files. `--public-typing` imports a frozen JSON list or an
`assignments` envelope keyed by exact Pathogenwatch `source_genome_id`. It can
contain cgMLST profiles, cgLIN assignments or HierCC assignments. Importing a
query file alone does not make public profiles reference-compatible; context
selection still requires compatible public typing.

## Profiles without a lineage reference

When the cgMLST caller index is ready but the LIN/HierCC background is unavailable,
normal `prepare --typing-config CONFIG` automatically uses profile-only calling
and records why lineage assignment was unavailable. No separate intermediate file
is required. To retain a standalone typing artifact, explicitly request profiles:

```bash
pixi run chronoclade type-queries focal.csv --profiles-only \
  --typing-config ~/.local/share/chronoclade/typing/query_config.json \
  --output focal-profiles
pixi run chronoclade prepare focal.csv --input-kind assemblies \
  --query-typing focal-profiles/query_typing.json --out prepared
pixi run chronoclade cgmlst prepared/dataset.json --out profile-cgmlst
```

This uses the same pinned caller, parser, assembly identity checks and database
provenance as full typing. Missing, novel and ambiguous calls stay distinct.
The saved mode is `profiles-only`; lineage assignment is `not_requested`, and
no LIN/HierCC group is inferred. Full typing still requires a verified assignment
reference. Results from the two modes cannot silently reuse each other's cache.

Default `--partition-by profiles` groups compatible profiles by species and
profile scope, independently of ST or lineage labels. Add `--fetch-context` when
known species/ST and provider access allow compatible public profiles to be
retrieved. The legacy `--partition-by st` route remains an explicit input-only
analysis. Choose `--partition-by lineage` only when verified HC/LIN assignments
are available and you intend those labels to restrict analysis.

## Comparing native calls with provider profiles

Pass the same ready native configuration to `cgmlst --fetch-context
--typing-config CONFIG` when comparing newly called inputs with provider exports.
For the explicitly supported Klebsiella, E. coli and S. aureus routes, ChronoClade
checks the biological scheme identity and translates locus-specific allele IDs
through the installed reference's DNA hashes. This produces a shared categorical
comparison namespace; it does not align alleles or turn allele mismatches into
SNP distances. The raw source profiles and normalization audit remain preserved.

Every declared reference locus stays in the panel. Provider allele numbers absent
from the verified reference remain unsupported/missing calls with explicit audit
reasons; they are never treated as matches. Coverage requirements still apply.
A matching species name, locus count or database date cannot substitute for this
compatibility check. The reference's actual contents are bound to resume evidence.
Lineage labels remain optional annotations in default profile mode.

## Klebsiella LIN

The [Pathogenwatch cgMLST caller](https://github.com/pathogenwatch-oss/mlst) uses
the `klebsiella_1` analysis scheme. Its output feeds
[plincer](https://github.com/pathogenwatch-oss/klebsiella-lincodes), which links to
the Pasteur `scgMLST629_S` LIN nomenclature. These are distinct names for the
caller and the lineage system. Ambiguous loci are excluded from assignment.
Partial or provisional LIN prefixes remain partial; tied reference profiles
cannot supply a more specific code than their common prefix.

Comparisons use the entire prefix at `--cglin-depth`, with a compatible scheme
and database version. Unversioned codes may compare within one frozen export,
identified as export provenance rather than a database release. Independent
unversioned exports are not silently combined.

## E. coli HierCC

`setup-typing --species ecoli --prepare-db` now installs the public
[EnteroBase Escherichia.cgMLSTv1 allele catalogue](https://enterobase.warwick.ac.uk/schemes/Escherichia.cgMLSTv1/)
through a frozen, resumable directory index. The catalogue currently has
2,513 loci; setup derives the count from the snapshot rather than assuming a
fixed count. It does not download the large profile table merely to read its
header. The original compressed allele files remain unchanged; a separate caller
view records numeric-header normalization, excluded non-ACGT alleles and hashes.
The caller database version is `snapshot:<reference SHA256>`, not an inferred
release date.

This installs **allele definitions and the cgMLST caller index**, not a HierCC
background. The public profile table has no HierCC columns, and the advertised
`hierBG` files returned access denied during testing. Upstream hclink database
builds still require their documented EnteroBase API access. Setup reports caller
and lineage readiness separately. Full native E. coli typing remains blocked
until a compatible HierCC reference and verified caller/reference compatibility
are available; a freshly downloaded caller snapshot is not silently declared
compatible with an older dated hclink database.


## Current module handoff

Use `type-queries` to create a saved typing result, then prepare the input dataset
and retrieve context through cgMLST. The independent workflow is documented in
[the tool suite guide](independent-suite.md); E. coli/HierCC remains outside the
accepted beta scope described above.

The `ecoli` discovery shortcut uses Pathogenwatch's `mlst` search field; `mlst2`
is an explicit alternative field. Always check that the chosen ST refers to the
same seven-locus scheme as your input. The `ecoli_1` cgMLST caller pairs with
[hclink](https://github.com/pathogenwatch-oss/hclink) and the EnteroBase
Escherichia/Shigella cgMLST reference. Its HierCC assignments retain every
returned level. `--hiercc-level` selects one level for context prioritisation;
missing levels remain unresolved. HClink references must be at least as recent
as the caller database. The canonical locus order is checked before assignment.

HierCC identifiers describe clusters at named thresholds, not LIN prefixes.
Do not compare HC1100 identifiers with HC10 identifiers or with Klebsiella codes.
The current country-group atlas uses Klebsiella cgLIN; E. coli HierCC currently
supports query assignment and comparison-pool refinement, rather than a HierCC
country-group atlas. Country counts and the location network are organism-neutral.

For cluster execution, use the generated stage job scripts documented in the
[independent workflow](independent-suite.md#local-and-slurm-jobs). Stage the
saved input and referenced files on the cluster; local validation does not
establish that a physical cluster run has completed.

## Evidence and limits

Typing outputs preserve the raw cgMLST and assignment JSON, assembly hash, tool
versions, database hashes, ordered loci, novel allele hashes, missing calls and
assignment status. cgMLST distances use shared callable loci and require at
least 90% overlap across the scheme; they are not SNP distances. Duplicate or
ambiguous calls cannot establish similarity. Conflicting fingerprints, schemes
or versions block comparison.

`context_selection.json` records typing-priority counts, background counts,
per-query comparisons and reasons for same-ST fallback. Unsupported, missing or
unresolved typing never becomes a made-up group. The candidate manifest retains
`pool_selection_reason` separately from the final SKA selection reason. This is
context within the analysed pool, not an exhaustive nearest-neighbour search.
