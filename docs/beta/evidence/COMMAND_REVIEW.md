# Simulated first-time command walkthrough

Date: 2026-10-18. Host: Darwin arm64. Repository: `/Users/nfareed/.codex/worktrees/pathogenwatch-context/chronoclade`.

## Outcome

The documented small offline example worked without credentials or manual data repair: prepare, cgMLST, GrapeTree export and Microreact export succeeded. The generated stage artifacts retained all three invented genomes (one input, two comparison records). Both the default environment and the documented lean prepare/cgmlst environments worked. The normal tree handoff accepted the generated ensemble but correctly stopped at missing assemblies; time retained that failure without inventing a dated result. SLURM script generation worked.

This is agent-simulated command acceptance, not evidence that a real first-time bioinformatician can use the software independently. It does not establish biological correctness or completion of the real Mac assembly-to-report workflow.

## Method and installation limits

Read README, the Mac Klebsiella walkthrough, linked installation/independent-suite/tree-stage/time-stage/CLI documentation, and examples/beta README and supplied inputs. Did not inspect previous validation notes or request help from another agent. The README/Mac walkthrough themselves contain historical validation claims; those claims were not used as proof. One ancillary read displayed the opening package metadata from pyproject.toml; no implementation or tests were read, and no inference depended on it.

Started with a new output folder under `tmp/beta-first-user-review/run`. Used the existing reviewed checkout and existing Pixi package/runtime caches. `pixi install --frozen --environment prepare` and later the cgmlst equivalent both returned success, but this was not a fresh clone or clean-machine dependency download. The default runtime already existed and was reused. No large reference databases were downloaded, no native phylogenetic inference was attempted beyond acquisition validation, and no credential values were displayed or changed.

## Exact executable workflow and observations

All commands below ran from the repository root. Exit codes apply to the named ChronoClade/Pixi command. The mkdir step was successful. Paths replace the example's `tmp/beta-example` only to guarantee a fresh output location.

1. Exit `0`. Installed/reconciled lean prepare environment; downloaded caches were already available.

```sh
pixi install --frozen --environment prepare
```

2. Exit `0`. Python profile workflow available; explicitly warns references/credentials separate.

```sh
pixi run -e prepare chronoclade preflight --stage prepare
```

3. Exit `0`. Fresh output directory; existing repository/environment reused.

```sh
mkdir -p tmp/beta-first-user-review/run
```

4. Exit `0`. Prepared three records, readiness ready; each assembly unavailable.

```sh
pixi run chronoclade prepare examples/beta/profiles.json --catalogues examples/beta/catalogues.json --out tmp/beta-first-user-review/run/prepared
```

5. Exit `0`. One block; 1 input, 2 contexts, 3/3 profiles; complete report, selection bundle.

```sh
pixi run chronoclade cgmlst tmp/beta-first-user-review/run/prepared/dataset.json --out tmp/beta-first-user-review/run/cgmlst --bootstrap-replicates 0
```

6. Exit `0`. Three genomes; profiles.tsv and metadata.tsv with import instructions.

```sh
pixi run chronoclade export tmp/beta-first-user-review/run/cgmlst/cgmlst.json --viewer grapetree --out tmp/beta-first-user-review/run/grapetree
```

7. Exit `0`. Three genomes; saved NJ Newick and matching metadata with import instructions.

```sh
pixi run chronoclade export tmp/beta-first-user-review/run/cgmlst/cgmlst.json --viewer microreact --out tmp/beta-first-user-review/run/microreact
```

8. Exit `2`. Expected boundary: 0 completed, 1 failed selection; all three records lack local assemblies/exact provider source IDs.

```sh
pixi run chronoclade tree tmp/beta-first-user-review/run/cgmlst/cglin_7bbf5faa1480/selection_bundle/ensemble.json --out tmp/beta-first-user-review/run/tree --threads 4 --randomisations 100 --ensemble-jobs 1
```

9. Exit `2`. Preserved failed tree selection and reason; comparison report written, no dates fabricated.

```sh
pixi run chronoclade time tmp/beta-first-user-review/run/tree/tree_ensemble.json --out tmp/beta-first-user-review/run/time --ensemble-jobs 1
```

10. Exit `0`. ska, iqtree, ClonalFrameML, treetime and Profile shown available; native platform supported.

```sh
pixi run chronoclade preflight --stage all
```

11. Exit `0`. Created job.json and run.slurm; placeholders intentionally retained; not submitted.

```sh
pixi run chronoclade job --slurm --account YOUR_ACCOUNT --partition YOUR_PARTITION --cpus 8 --memory-gb 32 --out tmp/beta-first-user-review/run/tree-job tree -- tmp/beta-first-user-review/run/cgmlst/cglin_7bbf5faa1480/selection_bundle/selections/selection-001.json --out tmp/beta-first-user-review/run/tree-001 --threads 8
```

12. Exit `0`. Same documented fixture worked in lean prepare environment.

```sh
pixi run -e prepare chronoclade prepare examples/beta/profiles.json --catalogues examples/beta/catalogues.json --out tmp/beta-first-user-review/run/lean-prepared
```

13. Exit `0`. Installed/reconciled lean cgmlst environment.

```sh
pixi install --frozen --environment cgmlst
```

14. Exit `0`. Same three profiles, one input and two contexts; complete lean-environment report.

```sh
pixi run -e cgmlst chronoclade cgmlst tmp/beta-first-user-review/run/lean-prepared/dataset.json --out tmp/beta-first-user-review/run/lean-cgmlst --bootstrap-replicates 0
```

## What I inspected in the generated outputs

- Prepared README and readiness JSON: `ready`, IDs fixture_1/2/3, one input, two contexts, three called loci per record; assemblies explicitly unavailable. This makes the distinction between profile readiness and tree readiness visible.
- Root cgMLST README links the `cglin_7bbf5faa1480` block. The block README gives `selection_bundle/ensemble.json` as the tree input and describes categorized files.
- GrapeTree profiles.tsv has three safe sample IDs and three numeric loci; generated README explains ID metadata matching, missing value 0 and allele_key.csv.
- Microreact Newick has three matching safe IDs; README instructs selecting ID as shared key and preserves source identifiers in id_map.csv.
- Tree ensemble records 0 completed/1 failed and a reason pointing to assembly_audit.json. That audit states for every fixture: `No linked local assembly or exact Pathogenwatch source_genome_id`.
- Time ensemble records the same failure; no successful native tree/time result exists.
- Generated job.json names run.slurm and records the requested resources and repository working directory. run.slurm quotes the working directory and runs `exec chronoclade tree ...` with the chosen selection. It needs a suitable executable environment on the worker, as documented.

Report HTMLs were produced and output navigation/readme text inspected. Browser rendering, responsiveness and interactions were outside this CLI walkthrough; this review must not be counted as visible report acceptance. No external viewer import/upload was performed.

## Friction and manual intervention

1. **Failure diagnosis is buried.** Tree stdout reports `Completed: 0; failed: 1` and the ensemble path without the acquisition reason. I needed to open tree_ensemble.json, then its deeply nested assembly_audit.json, to discover missing assemblies. The reason is accurate and the fixture README explicitly warns there are no assemblies, so this is an expected boundary rather than a broken fixture. Printing the first failure reason plus audit path would help a new user.
2. **Failed tree output still gives a next time command.** Generated tree README says to date the saved ensemble and gives an executable command even when no selection completed. It also says to inspect failures, but a user following the prominent next command reaches another exit 2 and a comparison report. A state-specific recovery instruction would be clearer.
3. **Block replacement requires navigation.** The Mac guide says replace `<block>` with a directory named by the cgMLST output README. The root README links a block rather than printing a complete ready-to-copy tree command. I combined that directory with the relative path from the block README. This is manageable for a command-line bioinformatician but still a manual handoff.
4. **The lightweight introduction initially activates default.** The synthetic example's unqualified `pixi run` uses the full default environment. The linked installation guide explains lean environments, and I verified those separately. A user seeking only the small example must discover the lean route separately.
5. **Installation claims need restraint.** Environment reconciliation succeeded using existing caches/runtimes; fresh download, clone availability, and dependency solver behavior on a new host remain untested.

No manual source-data correction, hidden reference setup, identity remapping or fabricated assembly input was needed. I deliberately stopped at the fixture's documented assembly boundary. The generated SLURM script was not submitted. A mistaken inspection command looked for `tree-job/job.sh`; cat returned missing file, then job.json correctly identified `run.slurm`. This was my filename assumption, not a documentation error.

## Acceptance boundaries

Supported here: documented offline profile commands, output handoff paths, both local viewer export formats, lean runtime use, expected missing-assembly failure, preservation of that failure through time, and SLURM script preparation.

Not supported here: fresh-machine installation, fresh assembly typing/reference builds, Pasteur authorisation, live public retrieval, authenticated downloads/cache use, successful SNP inference or temporal inference, real scheduler execution, browser report usability, or real human usability. The provided synthetic fixture cannot test those, as its own README states.
