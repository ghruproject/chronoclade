# Command line

## Independent module commands

Use these commands for the four-module workflow:

| Command | Input | Purpose |
| --- | --- | --- |
| `prepare SOURCE --out prepared` | Collection, accessions, assemblies or frozen profiles | Resolve input typing and metadata |
| `cgmlst prepared/dataset.json --out cgmlst` | Prepared bundle | Analyse compatible profiles and save selections |
| `tree SELECTION --out tree` | Selection JSON or `ensemble.json` | Build corrected SNP trees and assess dates |
| `time TREE --out time` | `tree.json` or `tree_ensemble.json` | Date supported trees and compare alternatives |

`cgmlst --fetch-context` retrieves matching public profiles and metadata;
assembly acquisition happens in `tree`. `--subsamples` defaults to 1 and requests
alternative sets from the same block. `--context-size` defaults to 50, with
mandatory neighbours retained even when they exceed that budget.
`--representative-reserve` defaults to 0; set it explicitly to retain room for
optional representatives beyond required neighbours. See [cgMLST](cgmlst.md),
[tree](tree-stage.md) and [time](time-stage.md) for stage options and examples.

Report and diagnostic utilities:

| Command | Purpose |
| --- | --- |
| `refresh-cgmlst-report CGMLST_JSON` | Refresh saved report views; retain numerical evidence and selections. Add `--recompute-diagnostics` to refit root/date diagnostics |
| `time-compare ENSEMBLE --run TIME_JSON ... --out comparison` | Compare checked completed time runs |
| `compare-clocks TREE_ENSEMBLE --out comparison` | Compare cgMLST and corrected SNP diagnostics on identical selected genomes/dates |
| `calibrate-lin ENSEMBLE --folds 3 --output calibration` | Experimental held-out-date checks for finer LIN groups; does not set biological cutoffs |
| `export SOURCE --viewer grapetree\|microreact --out DIRECTORY` | Write checked local viewer files from a saved module result; no upload |
| `list-schemes [--inspect]` | List registered public EnteroBase cgMLST references; optionally check indexes |
| `setup-scheme SCHEME --destination DIRECTORY` | Install/resume one verified public allele catalogue; caller and HierCC setup remain separate |
| `authorize-pasteur --pasteur-secrets KEYFILE` | Begin one-time browser authorisation; use `--complete` to enter the verifier privately |
| `job ...` | Save a local/SLURM stage script and resource manifest without submitting |

Run `chronoclade COMMAND --help` for the installed options. The older `run` and
`prepare-context` commands below have different context limits and output paths.
See [Pasteur authorisation](pasteur-auth.md) for the tool-install, consent and
authenticated reference-build sequence.

## Viewer exports

`export` reads a saved result and creates a **new** directory. It checks sample
identities and source evidence before publishing files, and never uploads data.
Keep the source bundle and its referenced files available while exporting.

```bash
# Export prepared categorical profiles for GrapeTree.
pixi run chronoclade export prepared/dataset.json \
  --viewer grapetree --out exports/grapetree

# Export the saved cgMLST NJ tree and matching metadata for Microreact.
pixi run chronoclade export cgmlst/cgmlst.json \
  --viewer microreact --out exports/cgmlst-microreact

# Export one recombination-adjusted tree or supported dated result.
pixi run chronoclade export tree/tree.json \
  --viewer microreact --out exports/tree-microreact
pixi run chronoclade export time/time.json \
  --viewer microreact --out exports/time-microreact
```

For a cgMLST result with several blocks or comparable cohorts, add `--block
BLOCK_ID` and `--cohort COHORT_ID` using the identifiers in the saved manifest.
The command asks for an explicit choice when selection is ambiguous. These
options do not apply to `dataset.json`, `tree.json` or `time.json`.

GrapeTree accepts a prepared dataset or cgMLST result with one compatible profile
catalogue. Microreact needs a saved tree: use `cgmlst.json`, an individual
`tree.json`, or an individual `time.json` containing a dated tree. For ensembles,
export each individual result linked by the ensemble manifest. If dating is
unsupported, export its tree-stage result instead.

Open the generated `README.txt` for import instructions. `export.json` records
source checksums, exported sample count and output checksums. Existing export
directories are never overwritten. See [viewer export files](outputs.md#viewer-export-files)
for identifiers, missing calls and file layouts.

## Native query typing

`setup-typing` installs pinned Pathogenwatch tools with uv-managed Python
environments and a native Node runtime. `--species klebsiella` and
`--species ecoli` can be combined. `--prepare-db` downloads/indexes cgMLST
data. Klebsiella uses the public route by default; with `--pasteur-secrets` and
a pre-existing Pasteur access/session token in the private Plincer cache, setup
uses an isolated authenticated downloader cache without starting browser
authorization. Protected `--pasteur-secrets` and `--enterobase-key-file` also
enable the corresponding lineage reference builds. Software-only installation writes
a configuration with `ready: false` until the databases are available.

`type-queries METADATA --typing-config CONFIG --output DIRECTORY` runs cgMLST
and the organism's LIN/HierCC assigner locally. `prepare-context` can use the same
`--typing-config` to type focal and bounded public assemblies, or import
`--query-typing` and `--public-typing` files. `--cglin-depth` and
`--hiercc-level` select compatible groups within the initial same-ST pool.
See [query typing](query-typing.md) for examples and database requirements.

## `chronoclade preflight`

Reports whether the native genome-analysis and temporal-analysis tools are
available in the active environment. Profile-first `fast` mode does not invoke
these assembly/tree tools.

```bash
pixi run chronoclade preflight
```

## `chronoclade validate`

Checks metadata, assembly paths, dates, lineage sizes and reference selection.

```bash
pixi run chronoclade validate metadata.csv --min-samples 10
```

## `chronoclade run`

Runs the cumulative profile-first workflow and writes an index for completed
stages. The default `fast` mode stops after profile analysis.

```bash
pixi run chronoclade run metadata.csv [OPTIONS]
```

| Option | Default | Meaning |
| --- | ---: | --- |
| `--output`, `-o` | `chronoclade_results` | Output directory |
| `--threads`, `-t` | `4` | Total CPU budget |
| `--lineage-jobs` | `2` | Maximum concurrent lineages |
| `--randomisation-jobs` | `4` | Maximum TreeTime permutations per lineage |
| `--mode` | `fast` | `fast` profile analysis, `full` corrected genomic analysis, or `finish` temporal analysis |
| `--collection` | none | Pathogenwatch collection UUID or full URL as query input |
| `--accessions` | none | Query accession list or CSV |
| `--species` | none | Declared species for new assemblies or accession lookup |
| `--catalogue` | none | Frozen public metadata catalogue |
| `--public-typing` | none | Frozen public cgMLST and group assignments |
| `--query-typing` | none | Verified query typing JSON |
| `--typing-config` | none | Native caller and reference configuration for query assemblies |
| `--cglin-export` | none | Frozen public cgLIN assignments |
| `--profile-limit` | `500` | Maximum public records in profile analysis; query records are additional |
| `--context-size` | `50` | Maximum selected context assemblies per lineage; query assemblies are additional |
| `--nearest-per-query` | `3` | Nearest profile relatives considered per query during context selection |
| `--profile-bootstraps` | `10` | Requested locus-bootstrap replicates for descriptive group support |
| `--group-distance` | `0.02` | Exploratory complete-linkage cgMLST mismatch fraction |
| `--date-randomisations` | `100` | Number of tip-date permutations in `finish` |
| `--date-randomisation-method` | `root-to-tip` | `root-to-tip` screen or `full-tree` TreeTime refits in `finish` |
| `--temporal-p-value` | `0.05` | Temporal gate threshold |
| `--min-samples` | `10` | Minimum genomes per lineage for whole-genome workflow |
| `--seed` | `20260818` | Selection and analysis seed |
| `--context-manifest` | none | Frozen manifest from `prepare-context` |
| `--force` | false | Rerun completed stages |
| `--dry-run` | false | Validate and print the plan only |

`--mode fast` compares the query set with the available public typing context.
It does not download context assemblies or build a corrected tree. Its report
shows query/context profile coverage, exclusions, country figures, PCoA and
neighbour-joining views where possible, descriptive genetic groups, separate
year-persistence and time/place-concentration summaries, an exploratory
allele-distance root-to-tip diagnostic, and nearest profile relatives.

`--mode full` selects up to `--context-size` public context assemblies per
lineage, retains all query genomes, and runs the corrected genomic analysis.
It stops before date tests and time scaling. `--mode finish` runs the same
corrected analysis and then performs the temporal tests; TreeTime creates a
dated tree only when the configured temporal gate passes. If there are fewer
than three usable distinct collection dates, `finish` retains the corrected
report and records that temporal analysis was not assessed.

Supply one query source: a metadata CSV, `--collection`, or `--accessions`.
`--dry-run` writes a plan without fetching typing exports, downloading
assemblies or running analyses. `--public-typing` and `--query-typing` import
local frozen files; they do not grant access to protected production databases.

## `chronoclade prepare-context`

Builds a reproducible public context set for one species and ST.

```bash
pixi run chronoclade prepare-context focal_metadata.csv \
  --scheme SCHEME \
  --st ST \
  --output DIRECTORY [OPTIONS]
```

Use `pixi run chronoclade prepare-context --help` for metadata filters and
screening controls. Start with `--dry-run` to inspect the candidate pool before
downloading assemblies.

Pathogenwatch is the built-in provider, currently supporting *K. pneumoniae*
with `--scheme klebsiella`. Other organisms and schemes return an explicit
unsupported-provider error. `--catalogue FILE` imports a verified frozen JSON
catalogue without live metadata access. Repeated runs reuse the
frozen metadata; `--refresh-catalogue` explicitly fetches a new snapshot.
`--cglin-export FILE` imports a validated CSV/TSV/JSON export and
`--focal-crosswalk FILE` imports sample IDs and strongly verified public IDs.
Without focal assignments, the report states that limitation.

Selection is recomputed when called; downloads resume from checksum-verified
source-ID caches.

Country figures can be regenerated offline from `context_catalogue.json`
without running phylogeny. See [context geography](context-geography.md).
