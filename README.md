# ChronoClade

ChronoClade compares bacterial genomes with public context using cgMLST profiles,
then builds recombination-adjusted SNP trees and assesses whether sampling dates
support a time tree. Four commands save separate datasets, reports and manifests:

| Command | Purpose | Next-stage input |
| --- | --- | --- |
| `prepare` | Resolve input genomes, typing and metadata | `prepared/dataset.json` |
| `cgmlst` | Retrieve matching context profiles, analyse relationships and select representative sets | `cgmlst/BLOCK/selection_bundle/ensemble.json` |
| `tree` | Acquire the selected assemblies, build corrected SNP trees and assess temporal signal | `tree/tree_ensemble.json` |
| `time` | Date supported trees and compare alternative selections | `time/time_ensemble.json` and its linked report |

Categorical allele-profile distance is the recommended cgMLST default. An optional
direct-DNA mode compares nucleotide sequences after aligning each locus separately
within the same cgMLST genes.

## Install

[Pixi](https://pixi.sh/) installs the Python package and native dependencies on
macOS or Linux. Windows users can use WSL2.

```bash
git clone https://github.com/ghruproject/chronoclade.git
cd chronoclade
pixi install
pixi run chronoclade preflight
```

Large cgMLST blocks require RapidNJ. On supported macOS hosts, run
`pixi run chronoclade setup-scale`; Apple Silicon currently requires Rosetta.
See [installation](docs/installation.md) and [large-block analysis](docs/cgmlst-scale.md).
Start with the [Klebsiella walkthrough](docs/mac-klebsiella-walkthrough.md)
or the [small offline example](examples/beta/README.md).

Separate Pixi environments are available for each stage.

## Prepare your input genomes

Supply a Pathogenwatch collection UUID or full collection URL, an accession list,
assembly metadata, or a frozen profile bundle. For example:

```bash
pixi run chronoclade prepare accessions.csv --input-kind accessions \
  --species "Klebsiella pneumoniae" --metadata study-metadata.csv --out prepared
```

Or start from a collection:

This collection contains 23 *Staphylococcus aureus* inputs. Keep these genomes:
they are the study inputs, not unwanted context. Compatible supplied allele
profiles can be analysed without Klebsiella cgLIN or E. coli HierCC references.
For Klebsiella on
Linux x86_64 and a Slurm cluster, use the [Linux walkthrough](docs/linux-klebsiella-walkthrough.md).

```bash
pixi run chronoclade prepare \
  https://next.pathogen.watch/collections/jX5cwsoyJ1KDquMqssUzAD-ellington-et-al-2015 \
  --input-kind collection --out prepared-collection
```

Preparation retrieves available typing and metadata for the input genomes. New
assemblies require a declared species, configured typing tools and ready reference
databases. Preparation reuses a supplied positive numeric MLST ST. If it is
missing, preparation attempts classical seven-locus MLST from a local input
assembly for Klebsiella pneumoniae, E. coli or S. aureus, using the provider's
primary scheme. Install its separate pinned environment once:

```bash
pixi install --frozen -e mlst
```

The normal `prepare` command discovers that installed sibling environment.
Missing assemblies, unavailable tools or incomplete/ambiguous calls remain
explicit; compatible local cgMLST profiles can still be analysed. This classical
ST call is separate from cgMLST and optional LIN/HierCC assignment. Installation
and reference access are separate requirements. Use
`setup-typing` to install the upstream cgMLST/LIN/HierCC tools; see
[input preparation](docs/prepare.md) and [query typing](docs/query-typing.md).

On Slurm, `slurm-prepare` distributes a local assembly CSV across bounded
workers, supports checked retries and collects every original input into one
prepared dataset. See [parallel assembly preparation](docs/prepare.md#prepare-local-assemblies-in-parallel-on-slurm)
and the [Linux walkthrough](docs/linux-klebsiella-walkthrough.md).

Pathogenwatch credentials come from `PATHOGENWATCH_API_KEY` or a protected local
configuration file. Keep credentials outside input files and version control.
See [Pathogenwatch access](docs/pathogenwatch.md).

## Analyse profiles and choose representative sets

```bash
pixi run chronoclade cgmlst prepared/dataset.json --fetch-context \
  --subsamples 1 --context-size 10 \
  --out cgmlst
```

The default `--partition-by profiles` groups genomes by species and compatible
cgMLST profile scope. LIN, HierCC and MLST ST are not prerequisites for local
profile analysis. Public discovery uses known species/ST to bound the candidate
pool, then checks profile compatibility. Distance comparisons determine nearest
relatives; available compatible lineage labels remain annotations.
Use `--partition-by lineage --lin-level 5` or an explicit HierCC level only when
you intend to restrict analysis by verified lineage assignments. Missing labels
remain unknown and are never copied from neighbours. Neither NJ distances nor
lineage prefixes define a universally validated close-relative cutoff. See
[comparison scopes](docs/cgmlst.md).
The command retrieves compatible context profiles and metadata;
context assemblies are acquired by the later tree stage. Omit `--fetch-context`
when the prepared bundle already contains the desired public comparisons.
For a count/cost preview without downloads or analysis, run:

```bash
pixi run chronoclade cgmlst prepared/dataset.json --context-preflight \
  --context-max-candidates 300 --out context-plan
```

The normal fetch route also previews every pool before detailed retrieval.
The guard refuses any advertised same-ST pool above 300 without truncation;
its default is unlimited. Review `context-plan/index.html` before increasing
the guard or planning Slurm resources. Repeat interrupted retrieval with
`--resume` to reuse checked search, metadata and profile batches.

Open `cgmlst/index.html` for the partition summary and links to each block report.
The report presents the combined input/context dataset, nearest relatives, the full
NJ tree, NJ-distance MDS, exploratory date regression, country/host/region/city state-change
networks, and the selections for SNP-tree analysis. It retains the existing report
style and uses accessions where available.

Tables have search, sorting and a 10-row default. NJ opens as a complete overview,
with zoom to inputs, a group or a sample. NJ and MDS share group highlighting;
their selectors show the 20 largest groups plus all input-containing groups, with
search across the complete group list. Default profile mode displays
distance-derived genetic groups, with verified lineage annotations shown
separately where available. Explicit lineage mode can display compatible
reported cgLIN prefixes. These labels are not validated close-relative cutoffs.
Nearest relatives rank by actual distances,
with coverage and ties retained. Additional automatic cluster selection that is
distinct and time-aware remains research, pending independent epidemiological
benchmarks, and is not a beta blocker.

Every comparable profile remains in the guide analysis. Subsampling applies to the
next stage. All inputs and required neighbours are retained; `--context-size` is a
context budget that can expand for mandatory neighbours, including boundary
ties. The example requests one selection and ten context genomes per block.
This budget limits the downstream assembly selection, not public profile retrieval
or the full guide analysis. Use `--representative-reserve` to request optional
representatives beyond mandatory neighbours. Multiple sets vary their optional
representatives to assess selection sensitivity.

See [cgMLST options and distance modes](docs/cgmlst.md).

## Build SNP trees, then time trees

Replace `BLOCK` with a block directory listed in `cgmlst/index.html`:

```bash
pixi run chronoclade tree cgmlst/BLOCK/selection_bundle/ensemble.json \
  --out tree --threads 4 --randomisations 100 --ensemble-jobs 1

pixi run chronoclade time tree/tree_ensemble.json --out time --ensemble-jobs 1
```

Tree acquires assemblies for the exact selected IDs, runs SKA2, IQ-TREE and
ClonalFrameML, and assesses temporal signal. Alternative selections share a
verified assembly cache, mapping reference and recorded seed policy. Time reads
the saved results and uses TreeTime when the temporal assessment supports dating.
It also writes an ensemble comparison report. A single selection manifest is
accepted too; that route produces `tree.json` and `time.json`.

Differences between alternative selections describe subsampling sensitivity.
They are not confidence intervals or proof of convergence. cgMLST date regressions
are exploratory; an improved regression alone does not establish temporal signal.
State-change networks describe reconstructed metadata changes on a rooted tree
and do not prove transmission.

See the [complete workflow](docs/independent-suite.md), [tree stage](docs/tree-stage.md)
and [time stage](docs/time-stage.md). The same commands run locally or through
[generated SLURM job scripts](docs/independent-suite.md#local-and-slurm-jobs).
Use [`slurm-ensemble`](docs/independent-suite.md#parallel-selections-on-slurm)
to distribute selections as bounded job arrays and collect their checked results.

Try the [small offline profile example](examples/beta/README.md) to check
preparation, cgMLST reports and exports without credentials or assemblies.
It uses invented data and does not validate biological inference.
For the real local Apple Silicon assembly-to-report route, see the
[Klebsiella walkthrough](docs/mac-klebsiella-walkthrough.md). It records which
references, context data and assemblies the accepted run reused.

## Klebsiella beta scope

The Klebsiella pneumoniae beta has been exercised on Apple Silicon macOS and
Linux x86_64. Linux Slurm validation includes native classical ST and cgMLST
calling, parallel preparation with failure/retry recovery, public-context
planning and retrieval, and selected SNP-tree/time runs. Temporal outcomes
depend on the cohort: unsupported dating remains a valid completed result.
See the [Linux walkthrough](docs/linux-klebsiella-walkthrough.md) for the tested
inputs, reference limits and commands. These runs do not establish biological
accuracy or validate every species, platform or cluster configuration.

## Reports and output folders

Each module writes its own report, machine-readable stage result and output guide.
New prepare, tree and time outputs include `README.md` with the next command and
recovery notes. cgMLST runs
also include root/block `README.md` files and a `files/index.html` browser grouping
tables, trees, figures, networks, subsamples and methods. Original evidence stays
at its recorded paths so hashes and downstream manifests remain valid.

Report figures and interactive views are embedded in the HTML and can be opened
without a server. Keep the output folder for downloads and subsequent stages.
Large distance matrices can occupy several gigabytes; avoid retaining unnecessary
copies of complete runs. See [output files](docs/outputs.md) and
[reading the reports](docs/report.md).

Use `chronoclade export SOURCE --viewer grapetree --out DIRECTORY`
to create checked local GrapeTree files, or choose `--viewer microreact` for a
saved tree and metadata. The command does not upload data. See
[export commands](docs/cli.md#viewer-exports) and
[export file descriptions](docs/outputs.md#viewer-export-files).

## Documentation and development

The [documentation site](https://ghruproject.github.io/chronoclade/) includes
installation, methods, command reference and the [Klebsiella example](docs/mac-klebsiella-walkthrough.md).
Documentation changes on a pull request reach the published site after merging
into `main` and successful deployment.

```bash
pixi run lint
pixi run test
pixi run docs
```

ChronoClade is under active development. Interpret results against the sampling
frame, metadata quality and local epidemiology.
