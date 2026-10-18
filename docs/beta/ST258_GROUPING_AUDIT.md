# ST258 relatedness scope audit

Date: 2026-10-17. Status: local diagnostic; no production grouping policy changed.

## Finding

The current report's broad search block and displayed groups are not satisfactory definitions of close relatives. The block is a useful diversity background, but conflating its membership with proximity makes the report misleading.

For this frozen example, level 8 is a promising **descriptive neighbourhood scope**: the union of the 12 level-8 prefixes represented among the inputs contains 758 public genomes, versus 7,774 in the broad block. It still includes every public genome tied for nearest to any input, and all 211 public genomes required by the current top-three-with-boundary-ties subsampling policy. This is a dataset-specific observation, **not** validation of level 8 as an outbreak definition, nor evidence of time-calibrated cutoffs.

## Implementation at original audit

- `chronoclade/context_partitions.py:partition_records` defaults to level 5 and allows only levels 5, 6 and 7. Its keys include species, MLST ST and typing namespace/version. It creates disjoint reports, not per-query nearest-distance groups.
- `chronoclade/report_components/lineage_groups.py:assign_lineage_groups` independently fixes the displayed prefix depth at 7. Therefore partition scope and plotting groups differ.
- `chronoclade/cglin.py` parses up to ten components, but its default exported grouping depths are only 5–7.
- Actual nearest-neighbour distances are independent of the plotting groups: categorical allele mismatch fraction over jointly called loci. The persisted nearest table retains ties.
- This audit reads prepared metadata/lineages and memory maps four existing binary matrices. It does not load the approximately 1 GB `profile_analysis.json`, rewrite scientific evidence, download assemblies or allocate another dense distance matrix.

## Current implementation status (2026-10-17)

The limits above describe the implementation when this audit was first run; they are
not current CLI limits. The current code supports explicit context and display
depths from 5 through 10, with separate controls for context selection, analysis
partition and displayed lineage depth. The defaults remain analysis depth 5 and
display depth 7. These are software defaults, not a newly validated biological
resolution. Automatic validated resolution selection has not been implemented.
The frozen numerical results below remain the original audit evidence and were not
recomputed for this documentation correction.

## Frozen dataset and depth comparison

Source: `tmp/st258-native-e2e/prepared/` and `cgmlst/cglin_37d98c4c5fe0/binary_distances/`.
There are 41 input genomes and 7,774 deduplicated public comparisons. These are the frozen accessible ST258 profiles, not all global CG258 genomes or all nearby STs. Claims of worldwide completeness are unsupported.

| Prefix depth | Groups in complete block | Input-containing groups | Public genomes in their union | Nearest pairs inside the particular input's own group | Extra genomes needed beyond union to keep every nearest |
|---|---:|---:|---:|---:|---:|
| 5 | 1 | 1 | 7,774 | 842/842 | 0 |
| 6 | 92 | 1 | 7,428 | 842/842 | 0 |
| 7 | 512 | 3 | 6,274 | 827/842 | 0 |
| 8 | 1,876 | 12 | 758 | 805/842 | 0 |
| 9 | 3,723 | 25 | 180 | 553/842 | 43 |
| 10 | 3,723 | 25 | 180 | 553/842 | 43 |

A nearest pair means one input–public pair; the same public genome can be nearest to several inputs. Thus 842 is not a count of unique public genomes. Level-9/10 input-group unions plus explicit nearest exceptions contain 223 public genomes. Cross-group nearest matches can already belong to a *different input's* group, explaining why own-group retention and union retention differ.

All 7,815 imported codes are marked complete, and every tenth component is zero. Levels 9 and 10 therefore provide identical partitioning here. Do not advertise additional level-10 resolution without investigating source exports and database provenance.

### Input-associated level-8 neighbourhoods

All prefixes start `0.0.105.6.0.0`; the abbreviated suffixes below are components 7–8.

| Suffix | Inputs | Public | Input–public mismatch range within prefix | Nearest pairs inside / total |
|---|---:|---:|---:|---:|
| `0.5` | 6 | 484 | 0–10 | 402/402 |
| `0.50` | 12 | 189 | 0–7 | 227/227 |
| `0.58` | 1 | 24 | 1–4 | 11/11 |
| `0.106` | 4 | 3 | 0–4 | 4/4 |
| `0.114` | 1 | 0 | — | 0/7 |
| `0.115` | 2 | 1 | 0 | 2/2 |
| `0.120` | 7 | 22 | 0–4 | 98/98 |
| `0.122` | 1 | 0 | — | 0/14 |
| `0.123` | 1 | 2 | 1–2 | 1/1 |
| `17.62` | 4 | 33 | 1–4 | 60/60 |
| `17.195` | 1 | 0 | — | 0/1 |
| `98.0` | 1 | 0 | — | 0/15 |

These ranges are actual input–public pair distances. Exact all-member mismatch diameters are 13, 10, 6, 4, 0, 0, 4, 0, 3, 5, 0 and 0 respectively in the table order (singletons have diameter zero). Even the largest level-8 group is not uniformly within four mismatches of every other member. Group membership does not guarantee every member is equally close to every input. The four zero-public prefixes must still display their external nearest matches, rather than “no relatives”.

## Actual distance evidence

Across all 318,734 input–public pairs, mismatch count quantiles are min 0, Q1 13, median 17, Q3 22, 95th percentile 28, max 54. Jointly called loci span 603–628; median 628. The export contains 628 observed loci; canonical 629-locus completeness is not verified.

Across 820 input–input pairs, mismatch count quantiles are min 0, Q1 6, median 9, Q3 15, 95th percentile 25, max 34. The inputs are not one uniformly tight group.

For the 41 individual inputs, the nearest public raw mismatch counts are:

| Mismatches to nearest public genome(s) | Inputs |
|---|---:|
| 0 | 18 |
| 1 | 15 |
| 2 | 6 |
| 4 | 1 |
| 7 | 1 |

The production nearest table and subsampling rank by raw allele mismatches, retaining every equal-count tie and reporting fractions and called-locus denominators. This audit reproduces all 842 saved nearest pairs exactly (190 distinct public genomes). A fraction-ranked comparison gives 811 pairs and 186 distinct genomes: coverage differences account for the discrepancy. The ranking policy must remain explicit and consistent; this is not a production mismatch.

The 37 pairs crossing a query's own level-8 prefix comprise 8 pairs at two mismatches, 14 at four, and 15 at seven. None is a zero-mismatch pair. Their jointly called loci are 628 for 19 pairs, 627 for 15, and 624, 625 and 626 for one pair each. These pairs therefore are not simply extremely low-coverage noise. Their underlying cause cannot be established from prefix labels alone; nomenclature history, reference matching and scheme/export compatibility need checking.

An instructive source discrepancy: `ERR8774297` and `SRR10028403` have the same full imported code `0.0.105.6.0.0.0.50.0.0`, but differ at three of the directly compared export loci. Their upstream records associate different cgSTs (191 and 22158), each at stated 100% identity to a reference, with 622 and 620 identical loci respectively. Thus matching imported full codes must not be presented as proof of identical observed profiles. Missing/reference-profile matching is a plausible contributor; reconstructing reference profiles and confirming a common reference snapshot is needed to determine the mechanism.

### Why prefix depth is not a guaranteed pairwise cutoff

Two concrete maximum-distance examples are retained in the machine-readable audit:

| Shared level-8 prefix suffix | Genome A | Genome B | Full code suffixes (components 8–10) | Raw mismatches / jointly called loci |
|---|---|---|---|---:|
| `0.5` | ERR1761546, Germany, 2012 | SRR5931295, Australia, 2016-10-11 | `5.106.0` / `5.93.0` | 13 / 627 |
| `0.50` | ERR11592230, Greece, 2017 | ERR11818272, Greece, 2022 | `50.59.0` / `50.69.0` | 10 / 628 |

These assignments share `scgMLST629_S` / `lincodes-3390273-2`, but their reference-database versions/fingerprints are unavailable. The actual profile catalogue declares `pathogenwatch:cgmlst-573-2`, 628 observed export loci, `complete=False`, and no database version/fingerprint. A common advertised analysis job is not proof of identical reference snapshots or a fully reconstructed canonical scheme.

The published LIN procedure copies a prefix from a nearest *previously encoded* profile. It does not describe complete-linkage checking against every existing member. A chain of genomes that each acquire one new difference can repeatedly inherit the same early prefix while its endpoints grow far apart. **This is an inference from the assignment procedure**, not proof of why these particular observed pairs received their codes. Therefore the nominal level-8 threshold must not be advertised as “every pair is within two alleles”. Confirming the observed mechanism requires source reference profiles and assignments; changing the report label alone does not validate a cutoff.

## Primary-source interpretation

Pasteur describes LIN as a stable, sequential, proximity-based nomenclature, rather than a fresh distance clustering of the present cohort. Its deeper bins are intended for epidemiological use; partial profiles and reference assignment affect interpretation. Its thresholds belong to the 629-locus scheme. Source: [Pasteur LIN specification](https://bigsdb.pasteur.fr/klebsiella/cgmlst-lincodes/).

The foundational study reports small observed differences among previously linked outbreak isolates but presents shallow classification as tentative surveillance evidence. This supports benchmarking against independent links rather than declaring a universal cutoff from group size or temporal regression. Source: [Hennart et al., 2022](https://doi.org/10.1093/molbev/msac135).

Pathogenwatch's own assigner consumes cgMLST and uses a reference catalogue. Its current README documents this dependency. Source: [plincer](https://github.com/pathogenwatch-oss/klebsiella-lincodes). The saved local upstream `classify.py` also explicitly skips missing calls and can return multiple reference matches. Its version is not assumed identical to the frozen online export.

## Recommended beta behaviour and acceptance gates

1. **Separate three concepts.** Preserve the broad search/background block, identify input-associated finer reported lineages, and list individual nearest relatives from actual distances. Use these exact distinctions in headings and denominators. Do not rename the entire block “close relatives”.
2. **Expose depth evidence.** Show a compact depth 5–10 summary, including inputs/public counts, own-group nearest retention and union nearest retention. Do not silently choose an arbitrary sample-count cap or universal biological threshold.
3. **Prototype, then validate a level-8 view.** For this ST258 example show the 12 input-associated neighbourhoods, their own nearest outsiders and the 758-public combined scope, with the broad background still available. An explicit user-chosen depth is safe sooner than an automatic “best” depth.
4. **Preserve nearest exceptions.** Any narrowed scope must retain every qualifying global nearest tie from the eligible background, annotate cross-prefix matches, and state that “nearest available” can still be distant. Narrowing must not redefine nearest by construction.
5. **Do not create extra trees as the first fix.** A shared full tree/ordination with linked subgroup highlighting and consistent filtered metadata summaries can support these views without twelve separate sprawling reports. Distinguish view filtering from recomputing an analysis.
6. **Validate before automatic defaults.** Evaluate independent linked/unlinked collections, an additional ST and an external study; report sensitivity to sampling, missing calls and scheme versions. Precision/recall require independent epidemiological labels. Levels 8/9 are candidate nomenclature resolutions, not validated epidemiological labels.
7. **Resolve full-code discordance.** Compare same-full-code pairs with source cgST profiles/alleles and upstream assignments, including repeated reference snapshots and incomplete profiles. Assert that report language never claims identity solely from a code.
8. **Time calibration remains research.** Neither root-to-tip slope nor a static prefix establishes a time-dependent inclusion rule. Compare any time-aware candidate against nearest-distance and reported-lineage baselines on temporally held-out linked/unlinked data; accept a negative result.

### Concrete regression checks

- Frozen ST258 check: 41 inputs; 7,774 contexts; level counts and 12 level-8 input groups exactly as above.
- Every input remains represented, including singleton or zero-public-prefix groups.
- Level-8 union retains all 842 nearest pairs and all 211 mandatory top-three public neighbours including boundary ties; level-9 union plus 43 nearest-only exceptions also retains all nearest-minimum pairs, but needs 63 exceptions (243 total public genomes) to retain all 211 mandatory top-three neighbours. An “own prefix only” implementation fails intentionally.
- Compare distances before and after scope selection by identity, not row position; no biological evidence is overwritten.
- Missing/provisional prefixes are unassigned at unresolved depths, not padded with zero or borrowed from a neighbour.
- Coverage and namespace incompatibility produce explicit warnings/exclusions; zero distance with missing calls is never labelled whole-genome identity.
- Country, time and host breakdowns use the selected scope denominator, and show missing data rather than inventing values.
- Identical input and cache produce identical selection; source-accessible totals, deduplicated totals and narrowed totals remain distinct.

## Reproduce

```sh
pixi run python validation/beta/st258_grouping_audit.py --run tmp/st258-native-e2e
```

Machine-readable evidence: `validation/beta/st258_grouping_results.json` (about 630 KB, no matrix copy). It includes per-input nearest labels, coverage, prefix crossovers and per-group counts. Inputs are local frozen artifacts; no credential or network access is needed. Numerical summaries are verified against the frozen binary evidence, not a new assembly-typing concordance test.

Verification: the frozen run diagnostic matched the persisted nearest table (842 pairs) and persisted mandatory top-three set (211 distinct public genomes) exactly. Three small fixture tests passed, including raw ties with differing coverage, cross-prefix nearest exceptions and top-three exceptions; Ruff passed for both diagnostic files.

## Local explicit-LIN8 pilot

`validation/beta/st258_relatedness_pilot.py` builds a new prepared bundle and report under `tmp/beta-st258-relatedness/`; it refuses to overwrite an existing pilot. The original ST258 outputs remain unchanged. It explicitly calls the scale/RapidNJ backend, then reuses the standard metadata, lineage, tree, ordination, root-diagnostic, network, report and selection components.

The run completed in 40.96 seconds on this workstation, retaining all 799 genomes and showing the 12 reported LIN8 neighbourhoods. It recomputed the scope-specific NJ tree, ordinations and country/host state networks. All 842 nearest pairs and all 211 mandatory public neighbours matched the original evidence. Five distinct new 302-genome representative sets were generated; SNP trees and dated trees were not rerun.

Entry point: `tmp/beta-st258-relatedness/cgmlst/profile_report.html` (about 4.85 MB, with all image sources embedded and two embedded interactive frames). `scope_acceptance.json` records the counts and scope. Source/profile namespace uncertainty remains explicit. This is a pilot for review, not a new production default or validated outbreak rule.

The shared lineage-display API gained an optional validated integer `level` (5–10), preserving default 7. Six lineage tests, plus group-tree and refresh checks (14 total), passed. The diagnostic and pilot scripts passed Ruff. Browser inspection is a separate presentation check; these checks establish generated content and data invariants.
