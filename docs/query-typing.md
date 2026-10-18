# Typing new query assemblies

Prepare the input assemblies first, then use the cgMLST module to discover
public comparisons. Species and seven-locus MLST ST identify the candidate
Pathogenwatch pool; compatible cgMLST profiles and LIN groups refine it. The
cgMLST report analyses the retained profiles and selects representatives for
the tree module. Context assemblies are acquired for the selected SNP-tree
analysis. A LIN code alone does not establish how close two genomes are.

The walkthrough below targets Klebsiella. E. coli HierCC requires a compatible assignment reference database.

## Install the tools

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

The input CSV needs the species and a verified seven-locus MLST ST for public
context discovery. Native cgMLST typing does not calculate that ST. Preparation
types new input assemblies using the configured caller and lineage reference.
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
request only the allele profiles explicitly:

```bash
pixi run chronoclade type-queries focal.csv --profiles-only \
  --typing-config ~/.local/share/chronoclade/typing/query_config.json \
  --output focal-profiles
pixi run chronoclade prepare focal.csv --input-kind assemblies \
  --query-typing focal-profiles/query_typing.json --out prepared
pixi run chronoclade cgmlst prepared/dataset.json --partition-by st \
  --out input-only-cgmlst
```

This uses the same pinned caller, parser, assembly identity checks and database
provenance as full typing. Missing, novel and ambiguous calls stay distinct.
The saved mode is `profiles-only`; lineage assignment is `not_requested`, and
no LIN/HierCC group is inferred. Full typing still requires a verified assignment
reference. Results from the two modes cannot silently reuse each other's cache.

`--partition-by st` is an explicit analysis of supplied inputs grouped by species,
ST and compatible profile scope. It does not discover public context or replace
an HC/LIN assignment. Standard lineage partitioning remains the default.

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
