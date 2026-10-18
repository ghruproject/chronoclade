# ST258 comparison: completed local run

18 October 2026. Start with `tmp/science-st258/comparison/report.html`.
The report uses the existing ChronoClade style, embedded figures and searchable,
sortable tables opening at ten rows. The evidence-folder README explains navigation.
No GitHub push was performed.

## What was completed

All 7,815 cached ST258 genomes were compared: 41 Greek inputs and 7,774 public
comparisons. The 628 observed loci were aligned independently using MAFFT,
including equal-length variable alleles, then reconstructed as a 512,056-column
concatenated alignment. The calculation retains joint callable sites/loci and
indel evidence. It does not imply that these loci are recombination-free.

Historical profile fingerprints verified the catalogue bridge for every genome.
DNA remained unavailable for 8,716 individual called cells (0.178%), rather than
being guessed. All three representations used the same 30,532,764 eligible pairs
out of 30,533,205 possible pairs. Coverage requires at least 597 joint loci.

## Main results

- All 41 inputs share at least one closest public comparison between categorical
  allele and aligned locus-SNP distances. Exact tied-neighbour sets agree for
  19/41 inputs. Restricting categorical alleles to the DNA-comparable loci raises
  agreement to 31/41. Missing DNA explains part, but not all, of the difference.
- Small single-linkage cutoffs can produce very broad groups. At two allele
  differences, 603 contexts are directly within the cutoff of an input, while
  chains link 5,089 contexts. Available within-group distances reach 24 alleles.
  Joining cutoff and pairwise maximum distance are different quantities.
- Context removal changes some boundaries substantially. With 25% of context
  retained, at the illustrative two-allele cutoff, mean per-input retained-member
  recovery ranges from 59.5% to 81.7% across five alternatives. At the independently
  measured two-SNP cutoff, the range is 66.7% to 85.5%. These equal numeric cutoffs
  do not imply equivalent biological boundaries.
- Date/country summaries describe groups; they do not determine them or establish
  a clock, transmission, origin or incidence. Random removal is a sensitivity
  analysis, not a model of surveillance bias. No universal cutoff was selected.

## Verification

Independent code verified all 32 full-dataset threshold partitions by graph
traversal, all 123 focal nearest sets and coverage, direct/path-linked context
counts, all 15 seeded sampling draws, and 288 sampling-score rows across nine
retained panels. The other 192 sampling-score rows were not independently graph
recomputed. Native categorical recomputation matched all 318,734 input-context
pairs. A separate character-level oracle matched all five aligned evidence counts
for 68 selected full-matrix pairs, covering all 41 inputs, missing DNA, gaps,
low coverage, diagonal entries, and nearest ties; 129 nearest-evidence rows match.

Receipts:

- `validation/science/st258/native-focal-audit.json`
- `validation/science/st258/comparison-audit.json`
- `tmp/science-st258/full-matrix-independent-audit.json`
- `tmp/science-st258/comparison/receipt.json`

## Laptop cost and remaining scope

Per-locus alignment and concatenation took 1,199.8 seconds (~20 minutes).
The full SNP matrix took 325.8 seconds (~5.4 minutes), with peak main-process
RSS 1.64 GB (not a system-wide memory measurement). Comparison and all sampling
runs took 48.3 seconds. Large evidence stays local rather than being copied into
reports or Git. Reusable alignment and threshold components are in ChronoClade;
the full-panel comparison is a reproducible research runner, not yet a new
polished CLI workflow or a validated epidemiological cutoff-selection method.
