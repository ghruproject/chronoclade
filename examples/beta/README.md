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
