# Klebsiella on an Apple Silicon Mac

This walkthrough takes a small *Klebsiella pneumoniae* dataset from local
assemblies through public cgMLST context, a recombination-adjusted SNP tree and
a decision about whether dating is supported. It uses the normal live-data
workflow. The final section records a separate two-isolate Mac validation; its
cached inputs are not needed for your own run.

## Install and prepare references

Install [Pixi](https://pixi.sh/latest/). Use the branch or commit named in the beta publication receipt until the
changes are merged into the default branch.
Install the locked environments and check the native programs:

```sh
git clone https://github.com/ghruproject/chronoclade.git
cd chronoclade
pixi install --frozen
pixi run --frozen --environment default chronoclade preflight --stage all
```

Preflight checks software, not typing references, credentials or your data. For
native assembly typing, follow the [Pasteur authorisation](pasteur-auth.md) and
[typing setup](installation.md#pathogenwatch-credentials) instructions. Prepare
the Klebsiella references and check their report and paths before use; setup can
take hours and is separate from installing ChronoClade.

Authenticated Pathogenwatch downloads use `PATHOGENWATCH_API_KEY` or the
protected file `~/.config/chronoclade/pathogenwatch.json` with permissions
`0600`. Follow [Pathogenwatch credentials](installation.md#pathogenwatch-credentials);
never put a key in a CSV, command history or report. Tree checks for a configured
credential before its assembly cache, so even an all-cache-hit run needs a valid
key.

Set paths for this analysis. Keep the run output together and choose a fresh
folder for each attempt:

```sh
RUN=results/klebsiella-study
TYPING_CONFIG="$HOME/.local/share/chronoclade/typing/query_config.json"
mkdir -p "$RUN"
```

## Make an input CSV

Create `$RUN/assemblies.csv`, with one row per isolate. Use an accession or
stable sample ID, an absolute path to the local assembly FASTA, the species, a
verified seven-gene MLST ST, and collection date/country where known:

```csv
sample_id,assembly,species,mlst_st,collection_date,country
ERR1234567,/absolute/path/to/assembly.fasta,Klebsiella pneumoniae,258,2019-05-14,Greece
```

Replace the example values with your own. Assembly paths may be relative to the
CSV. Include `source_genome_id` only with the exact provider UUID, so ChronoClade
can exclude that exact input from public context. Use collection dates, not
upload dates; unknown metadata can be left blank.

The seven-gene `mlst_st` is required for public species/ST context discovery.
The native cgMLST caller does not calculate it, and a cgMLST profile ID is not a
substitute. Resolve a missing ST with a verified source before requesting
context. Native assembly typing also needs BLAST, the cgMLST caller and
compatible references; use the full `default` Pixi environment below.

## Prepare and inspect the assemblies

```sh
pixi run --frozen --environment default chronoclade prepare \
  "$RUN/assemblies.csv" --input-kind assemblies \
  --typing-config "$TYPING_CONFIG" --out "$RUN/prepare"
```

Read `prepare/readiness.json` and `prepare/report.html`. Check that each intended
sample is present, its typing status is ready, and its metadata and locus
coverage are suitable for the next stage. Readiness does not mean every locus
was called or every metadata field was supplied. If preparation reports a
missing ST, caller, reference or BLAST executable, resolve that specific issue
before continuing.

## Retrieve public context and choose the cohort

Request the same-ST public profile pool, refined to the input-associated LIN8
scope, then analyse it in LIN5 blocks. The context budget is 50; all inputs and
all nearest-neighbour ties are mandatory, so the final selected cohort may be
larger.

```sh
pixi run --frozen --environment default chronoclade cgmlst \
  "$RUN/prepare/dataset.json" --fetch-context \
  --typing-config "$TYPING_CONFIG" \
  --lin-level 5 --context-lin-level 8 --display-lin-level 8 \
  --context-size 50 --nearest-per-query 3 --subsamples 1 \
  --seed 20261018 --out "$RUN/cgmlst"
```

LIN5 defines analysis blocks; LIN8 scopes context and reported group labels.
Neither is a genomic-distance cutoff. Neighbours are ranked from allele
comparisons within each block, and every tied nearest profile is retained, even
if that expands the budget. Provider results are accessible records, not
complete worldwide surveillance.

Open `cgmlst/index.html` and the relevant block report. Check retrieved and
included counts, shared-locus denominators, missing inputs, ties and selection.
Tables show 10 rows by default and support search and paging. Use the generated
`selection_bundle/ensemble.json` named by the output README for the tree. Keep
stage directories together so relative links and manifests remain valid.

## Build the SNP tree and assess time support

Pass the generated ensemble path from the previous stage. The tree stage
acquires assemblies for selected public genomes as needed; configure a valid
Pathogenwatch credential even if you expect cache hits.

```sh
SELECTION="$RUN/cgmlst/<block>/selection_bundle/ensemble.json"
pixi run --frozen --environment default chronoclade tree "$SELECTION" \
  --threads 4 --randomisations 100 --randomisation-jobs 4 \
  --support-replicates 0 --coherence-policy strict --out "$RUN/tree"

pixi run --frozen --environment default chronoclade time \
  "$RUN/tree/tree_ensemble.json" --out "$RUN/time"
```

Replace `<block>` with the directory named in the cgMLST output README. Strict
coherence is the default; inspect failures rather than switching to diagnostic
mode to force a result. Review the corrected tree, callable-site denominators,
recombination mask and temporal assessment in `tree/README.md`.

The time stage uses the saved tree; it does not infer a new SNP tree. An
`unsupported` assessment is a completed negative result, with no dated tree or
node-date intervals. Do not override it. If dating is supported, TreeTime's 90%
marginal max-posterior regions are conditional on the model and omit uncertainty
from missing genomes, residual recombination and the sampling frame.

## Outputs and troubleshooting

Open the preparation, cgMLST, tree and time reports in a browser. Follow each
stage README for the next output path and keep stage folders together. Completed
outputs are immutable; use a new folder for revised settings and inspect a
failed stage log before removing its work folder.

- **Native caller or BLAST missing:** use the full `default` Pixi environment;
  check typing setup status and reference paths in the configuration.
- **Context search cannot start:** verify the numeric seven-gene ST and
  configured Pathogenwatch access.
- **Context counts differ:** inspect the saved retrieval audit, scope and exact
  identity exclusions; provider results can change over time.
- **Tree coherence fails:** inspect raw and corrected distance audits. Strict
  mode stops for review.
- **Time is unsupported:** report the negative assessment; there is no dated
  tree or confidence region for that run.

## Recorded local check

A separate V14 validation followed two Greek ST258 isolates, ERR8771929 and
ERR8772239, collected on 2019-05-14. It reused a frozen provider snapshot,
configured native references and an assembly cache, so it is evidence for this
local replay rather than a fresh-download recipe. The frozen source had 13,421
accessible records; LIN8 refinement retained 386 source rows and identity
reconciliation produced 232 unique public comparisons. Nearest-neighbour ties
expanded the requested context budget to 26 contexts plus the two inputs. All
26 selected context assemblies were cached.

Both inputs were ready with 628/629 loci called. The tree completed 100/100 date
permutations over 24 date clusters. Dating was unsupported (p=0.0891 at the
0.05 threshold), so no dated tree or node intervals were produced. The recorded
run and independent evidence are under
`validation/beta/mac-native-final-plan/execution-v14/`; the validation manifest
and cached FASTAs are local files and are not required for the normal workflow
above.

A separate supported 63-genome benchmark was rerun with the final V16 build,
using its saved corrected SNP tree. It retained all 63 genomes and the same
rooted topology, and reproduced the MRCA estimate of 1995.6604 with a 90%
max-posterior region of 1994.268418–1998.078415. This checks the supported
dating route and its uncertainty output; it is not a fresh SNP-tree run.
The source tree, old result and new result remained separately checksummed.
Evidence: `validation/beta/positive-subsampling-plan/time-closeout-selection001-v16/`.
