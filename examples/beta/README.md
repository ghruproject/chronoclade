# Synthetic profile walkthrough

Three invented Klebsiella records and three invented loci, for installation,
module handoffs and viewer exports. These are not real typing schemes or
biological evidence. One input and two comparison records are already included;
no online fetch or credentials are needed. There are no assemblies.

From the repository root:

```sh
pixi run chronoclade prepare examples/beta/profiles.json --catalogues examples/beta/catalogues.json --out tmp/beta-example/prepared
pixi run chronoclade cgmlst tmp/beta-example/prepared/dataset.json --out tmp/beta-example/cgmlst --bootstrap-replicates 0
pixi run chronoclade export tmp/beta-example/cgmlst/cgmlst.json --viewer grapetree --out tmp/beta-example/grapetree
```

Open each output README and the cgMLST index.html. All three genomes should
remain in the analysed profile pool. Use a new output location for each exercise.
This fixture cannot exercise native typing, SNP trees, dated trees or realistic
relatedness. Full-stack and scientific acceptance use separate datasets.

## Synthetic direct DNA comparison

The separate `snp-*` files contain four invented local profiles and three
invented DNA loci. They exercise the direct DNA CLI without assemblies,
credentials, downloads or a real typing reference. The species is a scope label;
`beta:synthetic-snp`, database `synthetic-v1` and ST `synthetic` are invented.
These records make no real typing or biological relatedness claims. The original
allele-only example above stays independent.

Use an installation with MAFFT available (the repository Pixi environment
includes it). From the repository root, choose a fresh output directory:

```sh
pixi run chronoclade prepare examples/beta/snp-profiles.json --catalogues examples/beta/snp-catalogues.json --out tmp/beta-snp/prepared
pixi run chronoclade cgmlst tmp/beta-snp/prepared/dataset.json --out tmp/beta-snp/alleles --partition-by st --min-overlap 0.66 --bootstrap-replicates 0
pixi run chronoclade cgmlst tmp/beta-snp/prepared/dataset.json --out tmp/beta-snp/direct-dna --partition-by st --min-overlap 0.66 --bootstrap-replicates 0 --distance-method direct-dna --allele-catalogue examples/beta/snp-allele-catalogue.json --distance-threshold 0.04
```

An installed wheel uses the same commands with `chronoclade` in place of
`pixi run chronoclade`. ST partitioning here uses the supplied profiles only;
there is no context retrieval or lineage assignment. The overlap of `0.66`
explicitly permits the sample with only two of three loci called. The DNA
threshold `0.04` is an exercise setting with no biological interpretation.

Open `tmp/beta-snp/direct-dna/index.html`, then the block's
`pairwise_distances.csv` and `locus_alignment/` evidence. Both modes retain all
four profiles. Real MAFFT runs for all three variable loci: `snp-a` and
`missing-c` have equal-length variants; `indel-b` has 30- and 29-base variants.
The reconstructed concatenation is 90 columns wide, with one gap for the
single-base deletion and 30 `N` columns for `snp_missing`'s absent `missing-c`
call. Missingness and gaps never count as substitutions.

Expected pairwise evidence is frozen in `snp-expected.json`:

| Pair | Allele differences / shared loci | SNP differences / jointly called ACGT sites | Gap/base columns | Indel tracts |
| --- | --- | --- | --- | --- |
| snp_query / snp_equal | 1 / 3 | 1 / 90 | 0 | 0 |
| snp_query / snp_indel | 2 / 3 | 1 / 89 | 1 | 1 |
| snp_query / snp_missing | 1 / 2 | 1 / 60 | 0 | 0 |
| snp_equal / snp_indel | 3 / 3 | 2 / 89 | 1 | 1 |
| snp_equal / snp_missing | 1 / 2 | 2 / 60 | 0 | 0 |
| snp_indel / snp_missing | 2 / 2 | 1 / 59 | 1 | 1 |

The direct DNA distance is each row's SNP numerator divided by its ACGT
site denominator. These are substitutions within an invented gene panel,
not genome-wide SNPs, recombination-corrected distances or transmission evidence.
The report includes a neighbour-joining summary of these gene-panel distances.
This example does not exercise the assembly-based `tree` or dated `time` stages.

`snp-allele-catalogue.json` references `snp-alleles.csv` by a relative path and
validates its SHA256. The same digest identifies this entire invented database:

```text
7967e0ab38ece556f731e6fbe5787621084706237392429c81fbf07f6677e3fd
```

`snp-expected.json` also records the hashes of the CSV, profiles and both
catalogues, so the fixture can be copied and checked independently. No test-suite
module is needed at runtime.
