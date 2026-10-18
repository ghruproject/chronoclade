# Large-profile performance acceptance

## Current scope

This fresh check measures the **actual scale analysis engine** used by
`chronoclade cgmlst --backend scale` on every record in the frozen biological
ST258 dataset. It does not time the whole native command: portable selection
copies, stage report enrichment, extra clock-root optimisation and LIN-colour
regeneration are excluded. The normal native command was separately verified
on the 799-genome LIN8 context union.

The full biological workload is 7,815 profiles (41 inputs, 7,774 context records),
from the frozen prepared manifest
`tmp/st258-native-e2e/prepared/dataset.json`. Its observed Pathogenwatch export
contains 628 loci and declares canonical completeness unknown. It is not a
synthetic workload, nor does this computational test validate an epidemiological
cutoff, geographic inference, reference-database equivalence or temporal clock.

The former 10,000-profile synthetic benchmark remains historical evidence. Its
structured 20-template data, metadata and older output scope differ from this
full biological run. Timings should not be compared as if only genome count
changed.

## Method and completion gates

Run:

```bash
pixi run python validation/beta/performance-cgmlst-engine.py --inspect
pixi run python validation/beta/performance-cgmlst-engine.py
```

The harness refuses existing output, freezes the prepared manifest identity and
sample-ID digest, and runs a fresh child process. It sets numerical-library thread
limits to one, observes available resources, and waits when its guard fails.
The default guard requires 8 GiB of available/reclaimable memory, approximately
two idle CPU cores, and the conservative disk estimate plus a 10 GiB free-space
reserve. On macOS the memory figure is free + inactive + speculative pages,
which is a reclaimable estimate rather than physically free RAM.

The parent approved a **4 GiB cap on new output**, including transient files.
The harness monitors logical and allocated disk bytes once per second; an excess
interrupts only its own worker/process group and preserves the failed evidence.
It never deletes old runs or changes the scientific backend. A one-second sampler
can miss shorter peaks; that limitation is retained in the evidence.

Completion requires:

- All 7,815 input identities remain in the complete comparable pool, NJ tree,
  direct-distance PCoA and NJ-distance MDS, without duplicate or missing IDs.
- Both coordinate outputs contain finite values for every expected genome.
- All 30,533,205 pairwise comparisons are represented and no profiles excluded.
- A reproducible sample of 1,024 pairs matches the separate categorical reference
  implementation, including raw differences and jointly called loci.
- The actual wall time, OS-reported Python and child peak RSS, sampled concurrent
  process-family RSS, final/peak logical disk and allocated disk are recorded.
- The driver and worker complete successfully without breaching the resource cap.

The OS peak RSS values and sampled concurrent-family peak answer different
questions. Separate process peaks cannot be added and called an observed
simultaneous peak. Resident memory also differs from virtual allocation and
filesystem page cache.

## Updated storage accounting

For one fully comparable cohort of `n` profiles:

| Storage | Approximate bytes |
| --- | --- |
| Global distance and three callable-count arrays | `20 × n²` |
| Cohort square distance plus condensed vector | `12 × n²` |
| Saved NJ patristic distance matrix | `8 × n²` |
| Retained matrix subtotal | `40 × n²` |
| Portable selection copy, excluded from this engine run | another `20 × n²` |

The previous scale documentation's `32 × n²` engine estimate omitted the newer
saved NJ patristic matrix. The engine also writes transient streamed PHYLIP for
RapidNJ before it creates the NJ patristic matrix. The harness allows `24 × n²`
for this text scratch and 256 MiB for reports: conservative peak estimate
`56 × n² + 256 MiB`. The file is removed by the normal backend after RapidNJ
succeeds; the benchmark does not remove numerical evidence to meet its budget.

For 7,815 profiles, retained arrays alone are approximately 2.44 GB decimal
(2.28 GiB), and the conservative peak new-disk estimate is 3.69 GB decimal
(3.44 GiB). Thus a strict 2 GiB cap could not hold this current full engine run.
Actual logical/allocated sizes, rather than these estimates, determine the
recorded result. A full portable selection bundle would increase the footprint
and is outside this engine-only timing.

## Evidence and status

**Passed on 17 October 2026.** The fresh biological run completed on an Apple M1
Mac with 32 GiB RAM, macOS 14.6.1 and the pinned RapidNJ executable. The environment
records ChronoClade 0.1.0 (local beta working tree), Python 3.11.15, NumPy 2.4.6,
SciPy 1.17.1, Biopython 1.88 and Matplotlib 3.11.1. Package versions and engine
source hashes were recorded immediately after the completed run; they are labelled
as post-run observations in the consolidated evidence. Numerical
library thread limits were one. Two IQ-TREE jobs and other application activity
were present; this is an observed run under concurrent load, not an isolated
machine performance claim.

| Measurement | Observed |
| --- | --- |
| Scale engine wall time | 198.64 s |
| Worker data loading, engine and verification | 211.08 s |
| Driver wall time after resource guard passed | 212.30 s |
| OS peak Python RSS | 2.79 GiB |
| OS peak single child RSS (RapidNJ) | 0.35 GiB |
| Sampled concurrent process-family peak RSS | 2.79 GiB |
| Sampled peak allocated new disk, including PHYLIP | 3.01 GiB |
| Final allocated new disk | 2.39 GiB |
| Final logical file sizes | 2.38 GiB |

The full pool, NJ tree and **both** coordinate outputs retain all 7,815 unique
source IDs with finite coordinates. All 30,533,205 comparisons are available;
there are no profile exclusions. All 1,024 independently checked pairs matched
raw differences, shared-locus counts and normalized distances. The network
contains 41 country nodes and evaluates five candidate roots with its stated
fixed-tree uncertainty scope. The 4 GiB new-disk cap was not breached.

The 36.72 s compiled distance phase and 68.94 s network phase are included in the
engine total. They are not extra times to add to that total. Full CLI staging,
portable selections and biological inference validation remain outside this
resource result.

- `validation/beta/performance-results.json`: consolidated measured result, source identity, resource guard and exact byte counts.
- `validation/beta/performance-inspection.json`: source identity and initial plan.
- `validation/beta/performance-resource-preflight.json`: latest live resource guard.
- `validation/beta/performance-harness-tests.txt`: negative controls for disk
  preflight, missing/duplicate/nonfinite coordinate identities, and allocated disk.
- `validation/beta/performance-driver.log`: complete driver outcome.
- `tmp/beta-performance-st258-7815-engine/README.md`: generated output navigation.
- `tmp/beta-performance-st258-7815-engine/worker-result.json`: engine timing,
  correctness checks, backend and OS memory peaks after successful completion.
- `tmp/beta-performance-st258-7815-engine/monitor.json`: measured wall time and
  sampled process-family/disk resource outcomes after worker completion.

## Desktop viewer interaction check

The generated complete-tree and ordination views were opened in Chrome through
the supported browser automation. The tree displayed **7,815 tips**, a compact
legend of 25 of 328 diagnostic groups, working input/group zoom, and accession
lookup. Locating `ERR8772835` produced five tips at readable spacing and displayed
its country, host and collection date. Group selection enabled group zoom.

The full ordination allowed switching between NJ-tree and direct-profile
coordinates and country colours. Selecting a point displayed accessions and
metadata for seven overlapping/nearby samples. Both views had no captured
JavaScript console errors. These are actual functional desktop observations;
automation call timings include transport overhead and are not reported as
browser rendering benchmarks or guarantees for every device/browser. Mobile
layout and a complete user-report walkthrough are separate acceptance checks.
The engine views use optional distance-cut groups; the normal stage separately
replaces those display groups with reported LIN prefixes.

Evidence: `validation/beta/performance-browser-check.json`,
`performance-tree-locate.png`, and `performance-ordination-country.png`.
