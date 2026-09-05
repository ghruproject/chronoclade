# Running an analysis

## Inspect the plan

```bash
pixi run chronoclade run metadata.csv --dry-run
```

The plan lists each species/lineage group, its sample count, number of distinct
dates, chosen reference and whether it meets the minimum sample requirement.

## Run

```bash
pixi run chronoclade run metadata.csv \
  --output chronoclade_results \
  --threads 8 \
  --lineage-jobs 2 \
  --randomisation-jobs 4 \
  --date-randomisations 100
```

`--threads` is the total CPU budget. ChronoClade divides it across concurrent
lineages and caps TreeTime permutation workers so the run does not request more
CPUs than the total.

One hundred permutations give a minimum attainable corrected p-value of
1 / 101 = 0.0099. Larger tests give finer resolution but run TreeTime more
times. Choose the number before examining the result.

Completed stages are reused when their command and direct-input SHA-256
fingerprints still match. A changed assembly, metadata file, upstream tree or
option reruns the affected stage. `--force` reruns every stage.

## Add public context

The context step can be inspected without downloading assemblies:

```bash
pixi run chronoclade prepare-context focal_metadata.csv \
  --scheme ecoli_achtman_4 \
  --st 131 \
  --output context/ST131 \
  --candidate-pool 500 \
  --max-context 150 \
  --nearest-per-focal 3 \
  --dry-run
```

Remove `--dry-run` to fetch and screen the candidates. The command writes
`combined_metadata.csv` and `context_manifest.tsv` for the final analysis:

```bash
pixi run chronoclade run context/ST131/combined_metadata.csv \
  --context-manifest context/ST131/context_manifest.tsv \
  --output chronoclade_results
```

The context manifest records the source snapshot, candidate attrition,
screening distances and reasons for retaining each public genome.

## Demonstration data

The synthetic demonstration exercises both outcomes of the temporal gate:

```bash
pixi run python examples/demo/generate_demo.py --output demo_run/input
pixi run chronoclade run demo_run/input/metadata.csv \
  --output demo_run/results \
  --threads 8 \
  --date-randomisations 20
```

Open `demo_run/results/index.html`. The demonstration checks software behaviour;
it is not evidence that the method performs well on a real outbreak dataset.
