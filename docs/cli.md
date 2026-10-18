# Command line

## Independent module commands

Use these commands for the four-module workflow:

| Command | Input | Purpose |
| --- | --- | --- |
| `prepare SOURCE --out prepared` | Collection, accessions, assemblies or frozen profiles | Resolve input typing and metadata |
| `cgmlst prepared/dataset.json --out cgmlst` | Prepared bundle | Analyse compatible profiles and save selections |
| `tree SELECTION --out tree` | Selection JSON or `ensemble.json` | Build corrected SNP trees and assess dates |
| `time TREE --out time` | `tree.json` or `tree_ensemble.json` | Date supported trees and compare alternatives |

`cgmlst` defaults to `--partition-by profiles`: species and compatible profile
scope determine blocks, with no cgLIN/HierCC prerequisite. Explicit
`--partition-by lineage` requires verified labels; legacy `--partition-by st`
remains input-only. `cgmlst --fetch-context` retrieves compatible public profiles
and metadata from known species/ST candidate pools;
assembly acquisition happens in `tree`. `--subsamples` defaults to 1 and requests
alternative sets from the same block. `--context-size` defaults to 50, with
all inputs and mandatory neighbours, including ties, retained even when the
context count exceeds that budget. The budget does not limit public retrieval.
Use `cgmlst DATASET --context-preflight --out context-plan` to save a count/cost
report without details, profiles, assemblies or analysis. `--context-max-candidates N`
works with that preview or `--fetch-context`. All species/ST pool counts and
guards are checked before detailed retrieval begins; the guard never truncates
pools and defaults to unlimited. Repeat an interrupted fetch command with
`--resume` to reuse checked metadata and profile batches.
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
| `slurm-prepare ASSEMBLIES.csv --typing-config CONFIG ...` | Generate bounded per-assembly preparation workers and a strict full-cohort collector; no submission until its submit script runs |
| `job ...` | Save a local/SLURM stage script and resource manifest without submitting |
| `slurm-ensemble tree\|time ...` | Generate bounded selection job arrays and checked collectors; run the generated submit script to submit |

Run `chronoclade COMMAND --help` for the installed options. See
[Pasteur authorisation](pasteur-auth.md) for the tool-install, consent and
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
and the organism's LIN/HierCC assigner locally. Its saved output can be used when
preparing inputs. See [query typing](query-typing.md) for examples and database
requirements.

## `chronoclade preflight`

Reports whether native genome-analysis and temporal-analysis tools are
available in the active environment. Use `--stage` to check a stage's
requirements where supported.

```bash
pixi run chronoclade preflight
```
