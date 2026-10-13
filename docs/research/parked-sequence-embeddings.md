# Parked research: sequence embeddings for bacterial clustering

**Status:** parked on 11 October 2026. The experiments below are exploratory
research records. They do not justify an embedding workflow or a clustering
default in ChronoClade. The ESM2 workflow and user options have been removed
from the active codebase. Reopen this work only as a separate research task with
a defined validation question and the archived source snapshot.

## Why this is parked

The work established that several sequence models can produce reproducible,
non-identical pairwise distances from the pilot sequences. It did not show that
those distances find more useful within-lineage clusters than allele profiles
or direct sequence comparisons. A model can preserve broad distance ranks while
changing nearest neighbours or tree splits. Protein embeddings also cannot
distinguish DNA alleles that translate to identical proteins. There is not yet a
validated rule for turning embedding distances into biological groups.

ChronoClade's current core analysis is based on compatible cgMLST allele
profiles, with a direct DNA comparison available as a distinct analysis mode.
The experiments here remain useful background if a future study asks whether
learned sequence representations reveal a defined property that those methods
miss.

## Research questions and what was tested

Three saved comparisons were used:

1. A natural paired comparison of cgMLST allele distances, direct DNA distances,
   two nucleotide model families and two ESM2 protein checkpoints on the same
   20 ST39 genomes and 614 loci.
2. A sensitivity study of model layer and pooling choices on those same natural
   profiles, including repeated locus subsampling.
3. A controlled mutation simulation using real gene templates, held-out genes,
   and simple synthetic mutation histories. Model configuration was selected
   on discovery data before held-out scoring.

These answer implementation and measurement questions. They do not establish
outbreak-cluster performance, broad phylogenetic accuracy, or epidemiological
utility.

## Findings from the natural 20-genome comparison

The main saved comparison is
`tmp/dna-embedding-benchmark/comparison/comparison.json`.
Its paired input was frozen from
`tmp/paired-cgmlst-esm2/prepared-snapshot/dataset.json`
and the sequence panel at `tmp/paired-cgmlst-esm2/catalogue/alleles.csv`.
There were 20 genomes, 614 matched loci, and 13 distinct DNA profiles. The
profile collection had 22 variable loci. The comparison therefore contains 190
pair rows, but only 13 unique sequence profiles and dependent pairwise
observations.

The direct DNA and cgMLST distance ranks were strongly associated (Spearman
rho 0.878 across 190 pairs). Their unrooted NJ split comparison had Jaccard 1.0
in this small cohort; nearest-neighbour mean Jaccard was 0.938, with exact sets
for 18 of 20 genomes. This is evidence that both views captured similar broad
relationships in this particular panel, not evidence that either is generally
superior.

DNA model representations tracked the direct DNA rankings more closely than
the protein embeddings did. DNABERT-2 and Nucleotide Transformer pairwise rank
correlations with direct DNA were 0.941 and 0.883. Their nearest-neighbour sets
matched the direct DNA sets for all 20 genomes under the saved tie rule, though
their unrooted NJ split Jaccards against cgMLST were 0.167 and 0.313. The
protein ESM2 8M and 35M distances had rank correlations with cgMLST of 0.494 and
0.555; each resolved 58 of the 67 zero-distance pairs produced by protein
identity. Their nearest-neighbour average Jaccard against cgMLST was 0.254,
with exact sets for 2 of 20 genomes; their NJ split Jaccard against cgMLST was
0.083. Distance correlations are descriptive: pairs share genomes and no
inferential p-values were calculated.

The paired natural analysis is also saved in
`tmp/paired-cgmlst-esm2/comparison-35M/comparison.json`.
All 20 genomes were included. The scheme contained 629 loci and the fixed
protein panel included 614 (97.6%); cgMLST panel calls were complete. Across 190
pairs, full-scheme cgMLST versus ESM2 rank rho was 0.376; fixed-panel cgMLST
versus ESM2 was 0.555. The full-vs-panel cgMLST rho was 0.832, showing that
restricting to the protein panel changes the comparison independently of
protein representation. The panel cgMLST / ESM2 NJ split Jaccard was 0.083.
There were 67 zero-distance ESM2 pairs and 58 pairs with differing panel DNA
alleles but identical protein representations. This is expected information
loss from translation; it is not a model bug.

ESM2 used checkpoint `esm2_t12_35M_UR50D`, fair-esm 2.0.0, representation
layer 12, 480 dimensions, residue mean excluding BOS/EOS/PAD, float32 on Apple
MPS. The local checkpoint SHA256 was
`7f21e80e61d16a71735163ef555d3009afb0c98da74c48e29df08606973cc55e`.
The 8M and 35M implementations produced nearly identical distance ranks in the
natural pilot (rho 0.975), which demonstrates consistency between those two
protein checkpoints on this input, not biological correctness.

## Findings from representation and pooling sensitivity

`tmp/dna-embedding-benchmark/pooling/comparison/results.json`
compares six DNABERT-2 layer/pooling variants and six Nucleotide Transformer
layer/pooling variants, plus DNA and cgMLST. The cohort and panel match the
natural comparison: 20 genomes, 614 loci, 13 unique profiles and 22 variable
loci. The protocol used 100 seeded subsets of 491 of 614 loci, sampled without
replacement with seed `20261011`; the same subsets were applied to every
representation. Cuts at k=2, 3 and 4 were declared resolution probes only, not
biological cluster definitions.

Removing invariant loci changed only distance scale in every method. In other
words, the pilot's ranking was driven by the 22 variable loci, while the many
invariant loci diluted normalized distance magnitudes. This is a property of
the cohort and chosen aggregation, not a general recommendation to remove
invariant loci.

At k=3, the direct DNA partition had adjusted Rand index (ARI) 1.0 against the
cgMLST partition, while the selected Nucleotide Transformer variant had ARI
0.383. Across locus subsets, that embedding's k=3 partition had median ARI
0.383 against its full-panel partition, ranging from 0.149 to 1.0. The
Nucleotide Transformer variants were not interchangeable: their k=3 ARIs
against one another and against DNA varied by layer and pooling. Chord and
cosine distances happened to give identical partitions for each model variant
at the tested cuts, but this does not settle the metric choice for other data.
These cuts are arbitrary and the small number of unique profiles limits what
the ARIs can say.

The output has full metric details and the sampled partition labels. Do not
interpret a high two-axis plot separation or one favorable cut as validation.

The compared representations were not one fixed embedding recipe. DNABERT-2
used hidden layers 6, 11 and 12 with either an equal mean over DNA tokens or a
mean weighted by each token's nucleotide span. Nucleotide Transformer 50M used
halfway, penultimate and final layers with token-mean or nucleotide-span-mean
pooling. The pooler output was not used. Each genome-pair distance averaged
per-locus cosine distances while preserving locus identity. This grid is why a
future evaluation must predeclare a small set of settings and choose them only
on discovery data.

## Findings from the synthetic mutation study

`tmp/dna-embedding-benchmark/simulation/comparison/results.json`
records the controlled study. It generated sequences from 12 real gene
templates, used six disjoint genes for discovery and six for held-out scoring,
and evaluated 86 sequence sets. Mutations were substitutions only; indels,
recombination and homoplasy were absent. The primary target was the equal mean
of per-locus Hamming fractions; raw SNP counts were retained separately.

The model, hidden layer and pooling were selected using discovery ladder rank
correlation, then frozen before held-out scoring. The chosen setting was
Nucleotide Transformer 50M, halfway layer, base-weighted token mean with cosine
distance. On the 43 held-out sets, it averaged Spearman rho 0.947 against the
per-locus Hamming target, raw-SNP rho 0.939, and nearest-neighbour Jaccard
0.872. Direct SNP distance was the reference comparison (rho 1.0 to the
per-locus Hamming target by construction, raw-SNP rho 0.985, neighbour Jaccard
1.0). A gene-mismatch proxy scored lower (mean rho 0.480 and neighbour Jaccard
0.291). These results show that some embeddings tracked the simulated
substitution target in this design; they do not show that embeddings improve on
direct SNPs.

The study deliberately has limitations: only substitutions were simulated,
mutation histories were simplified, gene templates were few, tree experiments
had three replicates per mixed condition, and a gene mismatch proxy is not an
actual cgMLST typing workflow. The selected setting was not chosen to optimize
the held-out data, but selection on the discovery data still defines a narrow
simulation target. Treat the held-out scores as a pilot result, not as a
performance guarantee.

## Reproduction and retained evidence

The benchmark workspace is under `tmp/`, which is temporary and may be ignored
by Git. Keep the artifacts intact when reopening this work. The JSON reports
contain source checksums; the central saved artifacts are:

| Analysis | Reproduction script | Main result | Frozen inputs / key provenance |
| --- | --- | --- | --- |
| Natural multi-method comparison | `tmp/dna-embedding-benchmark/compare.py` | `tmp/dna-embedding-benchmark/comparison/comparison.json` | Paired 8M comparison SHA256 `e4fad7dc836b9cc6a2edff678bc70804afb35f65a21c7734a7a5ad528b298d5c`; arrays SHA256 `74defe3922011c7e05442341404f2f2653c6afd636bad0f7f40ad25d1dc2bf75`; DNABERT-2 checkpoint SHA256 `7ff39ec77a484dd01070a41bfd6e95cdd7247bec80fe357ab43a4be33687aeba`, repo `zhihan1996/DNABERT-2-117M`, revision `7bce263b15377fc15361f52cfab88f8b586abda0`; Nucleotide Transformer revision `81b29e5786726d891dbf929404ef20adca5b36f1`, weights SHA256 `17e75af297556ea56828716d8aa539e8f12b8b625547204b74171ab91aa33569`; ESM2 natural comparison used `esm2_t12_35M_UR50D`, checkpoint SHA256 `7f21e80e61d16a71735163ef555d3009afb0c98da74c48e29df08606973cc55e`. |
| Layer/pooling/locus sensitivity | `tmp/dna-embedding-benchmark/pooling/analyse.py` | `tmp/dna-embedding-benchmark/pooling/comparison/results.json` | Script SHA256 `69f09317751557b6b60b0dc2ad34501c387d29ecc246f24e91a6edb44127409a`; exact vector inputs and hashes are listed in `source_sha256` inside the result JSON. |
| Synthetic mutation study | `tmp/dna-embedding-benchmark/simulation/analyse.py` | `tmp/dna-embedding-benchmark/simulation/comparison/results.json` | Script SHA256 `ac50b86dc9fd60b69c17533794d0376ceb974646df142d3c60d59d465b204e3a`; dataset SHA256 `91559d5b0653c145e105ff83b55c993469a79ef78a53383c3de7dd271597659e`; simulation seed `20261012`; exact model-vector hashes are in `source_sha256`. |
| Paired natural report | Existing paired comparison pipeline; result records all method settings | `tmp/paired-cgmlst-esm2/comparison-35M/comparison.json` | Prepared dataset ID `sha256:8ca3e097ae2ada8aa1f05cc1b25727f5829bba920463a859a4271adc3ebac509`; the result contains model, panel, topology and input provenance. |

The codebase archive is at
`tmp/parked-embeddings/source-snapshot/`
with a file-hash inventory in
`tmp/parked-embeddings/source-snapshot-manifest.json`.
It preserves the removed implementation, reports, fixtures and tests for
review. The snapshot is a recovery aid, not active ChronoClade functionality.
The model checkpoint files remain in the benchmark workspace; no installation
or runtime dependency should fetch them.

## Revisit only with a task-specific validation plan

Before writing more embedding code, specify what a successful embedding would
add beyond the two supported comparison modes. A useful next study would:

1. Define the intended target: ancestry, within-lineage genetic grouping,
   functional similarity, or another explicit biological property. Do not use
   the word “cluster” without defining what the groups are meant to represent.
2. Freeze a larger set of lineages and independent genomes with usable
   compatible allele profiles, sequence calls, and a suitable SNP reference.
   Keep train/discovery and evaluation genomes or whole lineages separate to
   avoid choosing models on the final evaluation.
3. Compare identical input scope using cgMLST categorical allele mismatch,
   direct DNA distances over the same cgMLST genes, and a predeclared small set
   of embedding configurations. Report missingness, callable bases, locus
   weighting, synonymous changes and all model/checkpoint provenance.
4. Measure the property the task needs: held-out group recovery, pairwise
   neighbour stability, agreement with independent SNP phylogeny, or external
   functional labels. Do not select an embedding by PCA appearance or by its
   agreement with an outcome that was used to tune it.
5. Test stability under resampling genomes and loci, lineage composition,
   sequence lengths, model layers/pooling, synonymous-only changes, and plausible
   mutation processes including recombination where relevant. Keep direct DNA
   as the baseline and distinguish a useful functional grouping from a lineage
   or transmission cluster.
6. If the question is specifically whether embeddings produce better genetic
   clusters, require reproducible gains over both allele and direct DNA methods
   on held-out lineages and an independently defined reference. Otherwise keep
   embeddings as research-only and do not restore a routine report option.

No acceptance threshold is set here; the scientific use case and evaluation
cohort need to be agreed before a new study starts.
