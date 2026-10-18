# Independent preparation

`chronoclade prepare` resolves live query collections, accessions or assembly
metadata, imports frozen canonical profiles, or validates and republishes an
existing prepared dataset. It writes `dataset.json`, typed categorical allele
matrices, inspection CSVs, complete imported lineage and metadata evidence,
`readiness.json`, `report.html` and a checksummed `prepare.json` stage result.
The complete output directory is published atomically; an existing output is never
replaced. Choose a new output directory for a revised input.

```console
chronoclade prepare profiles.json --catalogues catalogues.json --out prepared/
chronoclade prepare prepared/dataset.json --input-kind dataset --out imported/
chronoclade cgmlst prepared/dataset.json --out cgmlst/
```

Frozen profile and dataset imports are offline. Live preparation uses the same
portable bundle contract, with exact query identity, provider response snapshots,
metadata precedence/conflicts and per-sample retrieval/readiness evidence. It does
not discover public context or download context assemblies. Linked local input
assemblies are copied under content hashes and remain portable.

## Live query preparation

```console
chronoclade prepare https://pathogen.watch/collections/jX5cwsoyJ1KDquMqssUzAD \
  --input-kind collection --out prepared-collection/
chronoclade prepare accessions.csv --input-kind accessions \
  --species "Klebsiella pneumoniae" --metadata overrides.csv --out prepared-accessions/
chronoclade prepare assemblies.csv --input-kind assemblies \
  --typing-config ~/.local/share/chronoclade/typing/query_config.json --out prepared-assemblies/
```

Collection identifiers accept a short UUID or a full HTTPS Pathogenwatch collection
URL, including its optional slug. Membership pages must reconcile exactly with the
advertised count. An accession input is a unique one-ID-per-line list or CSV with
`accession` or `source_genome_id`; CSV metadata is applied as an explicit override.
Public accession search requires species, a complete bounded response and exactly
one exact accession match. Ambiguous identities fail without choosing a neighbour.
An assembly input is a CSV with `sample_id`, `assembly` and `species` (or
`--species`). Assembly paths are relative to the CSV. Additional `--metadata`
overrides join by exact ID and take precedence over provider values; differences
remain in the conflicts table.

Existing server-advertised cgMLST and cgLIN analysis exports are retrieved for the
resolved input IDs. The configured Pathogenwatch key is read from the environment
or the protected local configuration; it is never accepted as a command-line
argument or written into snapshots. Public searches remain credential-free.
`--query-typing`, `--public-typing` and `--cglin-export` accept validated frozen
assignments. Query typing must match exact sample identity, species and any known
input assembly SHA256. Missing profiles and lineages remain separate checks.

`--typing-config` enables the pinned native upstream caller and assigner. A
native assembly call needs `blastn` and `makeblastdb` on the executable path:
run it with `pixi run --frozen --environment default chronoclade prepare ...`.
The lean `-e prepare` environment handles existing profiles/provider typing
and does not include those executables for native cgMLST calling.

For public species/ST discovery, preparation reuses a supplied positive numeric
classical MLST ST. When that ST is missing and an input assembly is available
locally, it attempts an independent classical MLST call for supported species.
Install the separate pinned environment once:

```bash
pixi install --frozen -e mlst
```

The normal default/prepare environment discovers the installed sibling `mlst`
environment automatically. The separate environment keeps its BioPerl dependencies
from changing the native inference environment. An unavailable caller is reported
with its setup command. The bundled reference is a dated database snapshot;
unknown recent alleles or profiles may remain unresolved.

The supported primary mappings in the pinned `mlst` 2.32.2 package are:

| Species | Local caller scheme | Pathogenwatch search namespace |
| --- | --- | --- |
| Klebsiella pneumoniae | `klebsiella` | Primary `mlst`, Pasteur seven-locus scheme |
| Escherichia coli | `ecoli_achtman_4` | Primary `mlst`, Achtman seven-locus scheme |
| Staphylococcus aureus | `saureus` | Primary `mlst`, PubMLST seven-locus scheme |

The caller must return the exact expected scheme, seven expected loci, a numeric
ST and complete unambiguous allele calls for the exact input assembly. The audit
preserves caller/database and assembly hashes and raw result evidence. Novel,
partial, missing or ambiguous calls do not receive an inferred ST. No local
assembly means the missing-ST reason remains explicit. Supplied numeric STs are
reused rather than overwritten; an explicitly incompatible supplied namespace
cannot seed a primary public search. In particular, the local `ecoli` scheme is
the alternative Pasteur eight-locus scheme, not the primary Achtman scheme.
Alternative-scheme public discovery is not implemented by this workflow.

Missing ST does not prevent compatible local profile analysis, but that input
cannot seed its own public same-ST search until resolved. Classical ST, cgMLST
profile IDs and LIN/HierCC assignments are different evidence.
With a complete assignment configuration, an existing compatible profile missing
its lineage can go directly to the assigner without acquiring an assembly.
If the cgMLST caller and its reference are ready but optional LIN/HierCC resources
are unavailable, preparation automatically falls back to profile-only calling.
Its audit records the fallback and the unavailable assignment reason. Missing
profiles are called from **input** assemblies; an existing profile is not recalled
solely to obtain an optional lineage label. A missing/unready cgMLST configuration
publishes an incomplete bundle with an explicit reference-database reason rather
than inventing calls. This fallback does not make lineage-restricted analysis ready.
Typing-tool installation and database readiness are separate audit fields. Existing
frozen typing does not require a Pasteur key. Fresh Klebsiella lineage assignment
requires its configured Pasteur reference database; E. coli assignment requires
its scheme-specific HierCC reference database.

For native Klebsiella reference setup, follow the [one-time Pasteur authorisation](pasteur-auth.md)
walkthrough. Client-key files, browser consent and reference installation are
separate steps; the command saves access privately without requiring a password
in your typing configuration.

Live metadata enrichment uses exact linked ENA run/BioSample accessions. It fills
missing country, collection date, host and isolation source only when returned
records agree. Existing values retain precedence, source raw values and conflicts
remain recorded, and unrelated identities are rejected. It never derives a region
or date precision from an upload date. `--no-enrich-metadata` disables enrichment.
Missing metadata and request failures preserve the sample.

Location metadata uses separate `country`, `region` and `city` columns. Supply
only the detail you have; ChronoClade does not infer a city from a region or a
legacy NUTS2 code. Reports omit entirely unreported region, city and host fields,
including their colour controls. Where a field is partly reported, missing
values remain visible as Unknown.

Provider cgMLST exports list observed loci only. Their automatically generated
catalogue is explicitly incomplete; no full-universe call fraction is claimed.
Supply `--catalogues` with a matching full scheme catalogue to assess completeness.
Native callers provide the ordered full locus universe. Source responses and exports
are copied into the checksummed bundle after credential checks.

## Prepare local assemblies in parallel on Slurm

Use `slurm-prepare` to distribute an assembly CSV across bounded, independent
workers. Each worker runs the ordinary `prepare` workflow, including missing
classical ST and cgMLST calling. Configure the typing reference and install the
separate classical MLST environment first; use shared paths visible to workers.

```bash
pixi run --frozen chronoclade slurm-prepare assemblies.csv \
  --typing-config /shared/typing/query_config.json \
  --out prepare-jobs --results prepared \
  --account YOUR_ACCOUNT --partition YOUR_PARTITION \
  --cpus 1 --memory-gb 4 --max-parallel 4 --walltime 01:00:00
bash prepare-jobs/submit.sh
```

Generation writes scripts and a checksummed plan without submitting. The submit
script schedules one array task per assembly and an `afterany` collector. Four
active workers in this example reserve four CPUs and 16 GiB in total, excluding
the collector. Metadata enrichment is disabled by default for this route; add
`--enrich-metadata` deliberately if required. Collection/accession inputs use
ordinary preparation; this parallel command accepts local assembly CSV input.

The collector validates completed worker artifacts and preserves the full original
cohort and input order. It does not run typing or assemble scientific intermediate
files by hand. If any worker is missing or failed, the collector exits nonzero and
does not publish a partial final dataset. Inspect task logs/status and
`collection-status.json`; retry the saved submission to reuse checked successes,
or submit only failed zero-based indices:

```bash
sbatch --array=0,2 prepare-jobs/array.slurm
# After those retries finish:
sbatch prepare-jobs/collect.slurm
```

For incomplete workers, the collector can make one read-only `sacct` query
(five-second timeout, at most 1,000 job IDs). It trusts only task-start IDs bound
to the saved plan or the submission ledger's plan checksum. Recognized scheduler
states and exit codes are retained in the collection status, helping distinguish
cancellation, timeout and out-of-memory termination. Missing or delayed accounting
is explicit; it never substitutes for checked scientific outputs or changes which
tasks require retry.

Keep the saved plan, source assemblies, reference contents, software and shared
paths unchanged when retrying. Changed inputs or resources require a fresh plan;
verified successful task results are reused without retyping. After collection,
continue with `prepared/dataset.json`, the normal context preflight and cgMLST
commands. No manual profile JSON or dataset construction is required.

## Resume interrupted native preparation

If live/local preparation using `--typing-config` is interrupted, repeat the
same command with `--resume`:

```console
pixi run --frozen --environment default chronoclade prepare assemblies.csv \
  --input-kind assemblies --typing-config /path/to/query_config.json \
  --out prepared-assemblies --resume
```

Preparation retains its work beside the requested output, for example
`prepared-assemblies.work/`. Completed sample typing is checkpointed; a sample
interrupted before completion is rerun in a new attempt directory. The error
message identifies the retained work and the resume command. Keep the input,
configuration, software and references unchanged when resuming. Changed or
corrupt recovery state is refused rather than silently mixed with new results.
Provider metadata may be retrieved again before exact cohort checks complete.

A completed output remains immutable, even with `--resume`. Frozen profile and
dataset imports do not use this native recovery mechanism. The standalone
`type-queries` command also accepts `--resume` for its own typing output.

## Frozen input contract

Supply a JSON list of canonical records, or an object with `records` plus optional
`lineages`, `crosswalk`, `provenance`, `conflicts`, `exclusions`, `retrieval` and
`parameters`. Evidence tables are lists; parameters is an object. Explicit evidence
tables override the frozen-record adapter's automatic evidence. Every evidence row
must refer to an existing `sample_id`. Unknown wrapper fields are rejected. The
legacy resolver's `queries`/`context` output requires explicit conversion to this
contract; the new command does not run the legacy resolver.

```json
[
  {
    "sample_id": "input-001",
    "role": "input",
    "species": "Klebsiella pneumoniae",
    "mlst_st": "39",
    "collection_date": "2024-02",
    "country": "Greece",
    "cgmlst_scheme": "example:cgmlst",
    "cgmlst_scheme_version": "v1",
    "cgmlst_status": "resolved",
    "cgmlst_profile": {"locus-a": "001", "locus-b": "7"},
    "cglin_scheme": "example:cglin",
    "cglin_scheme_version": "lin-v1",
    "cglin_status": "resolved",
    "cglin_raw": "1.2.3.4.5.6.7"
  }
]
```

Records require safe unique `sample_id` values. `role` is `input` or `context`;
the existing `origin` values `local`, `query`, `focal`, `input`, `context` are
also accepted. A record with neither role nor origin defaults to `input`.
Explicitly supplied context records are preserved without discovery or selection.
`--species` fills missing species values and rejects conflicting values. It cannot
override an existing prepared bundle's species. Date precision, null metadata,
ambiguous identities and complete lineage assignment strings/maps remain saved.
Profile allele values are categorical; missing calls never become allele matches.

An explicit catalogue JSON list is required for records, even when profiles are
unavailable. An object containing only `catalogues` is also accepted.

```json
[
  {
    "scheme_id": "example:cgmlst",
    "scheme_version": "v1",
    "loci": ["locus-a", "locus-b", "locus-never-called"],
    "source": "frozen-schema.json",
    "complete": true,
    "database_version": null,
    "database_sha256": null
  }
]
```

`complete` must be explicitly `true` or `false`. Do not claim completeness from
a union of observed calls. Every declared locus is retained, including loci absent
from every sample. Catalogue scheme/version and any database SHA256 must match
the records; incompatible called profiles fail rather than being discarded.
Missing profiles remain as samples, including all-null matrix rows when their
scope matches a supplied catalogue. Credentials are rejected before source
snapshots or outputs are published.

## Readiness and public profile context

Dataset publication has a separate completion status from typing readiness. The
audit marks available profiles, absent loci, incomplete catalogue universes and
lineage assignment independently. Klebsiella readiness expects a resolved scoped
cgLIN assignment, including the authentic provider status `complete`; partial and
provisional assignments remain incomplete. E. coli expects resolved HierCC. Other species have lineage readiness
explicitly unassessed. The readiness audit can therefore describe a profile with
missing lineage as incomplete for lineage assignment. That does not prevent
default profile-mode analysis when the profile itself is comparable; inspect
profile and lineage readiness separately. Missing metadata is counted and retained; it does not exclude valid samples or by itself
make existing typing incomplete.

Frozen-import audits say tools/databases were not checked because those imports
do not need them. Live audits distinguish configured, unavailable and ready
references, and whether native typing or assignment actually ran. Dataset
publication completion does not imply typing completeness.

The cgMLST provider hook `discover_profile_context` consumes a prepared bundle and
publishes a new prepared bundle for profile analysis. It searches the complete
accessible public same-ST pool and checks species and profile scheme/version/database
compatibility. The default profile mode does not require a cgLIN/HierCC label.
Explicit lineage mode additionally records available LIN5/6/7 counts and restricts
to the requested LIN depth (default 5) or explicit E. coli HierCC level; lineage
scheme/version/database scope must then be compatible. The hook
requests every matching profile, without a candidate limit or context assembly
acquisition. It deduplicates source identity, preserving different samples with
identical profiles; exact inputs are excluded. Unavailable or incompatible profiles
remain audited exclusions.

`readiness.json` and dataset parameters record discovered, requested, retrieved
and usable counts separately, pagination evidence, export status, exclusions,
source dataset SHA256 and the explicit policy. Input metadata, profile calls and
lineage/crosswalk/retrieval evidence are preserved. Accessible public matches are
reported; global completeness is never claimed. Raw provider snapshots are
checksummed under `sources/provider/`.
