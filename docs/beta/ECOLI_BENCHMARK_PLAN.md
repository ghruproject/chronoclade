# E. coli published benchmark freeze

Freeze date: 17 October 2026. Status: Australian metadata and relationship targets frozen; eight assembly inputs downloaded and identity-checked for a bounded native caller pilot. No ChronoClade biological analysis or temporal fit has been run. This benchmark is a method comparison, not a claim of exact reproduction.

## Main One Health benchmark

Use Watt et al., *Parameters for one health genomic surveillance of Escherichia coli from Australia*, Nature Communications, first published online 2 January 2025 (DOI `10.1038/s41467-024-55103-2`). The 2025 publisher correction (DOI `10.1038/s41467-025-58600-0`) states that the initially supplied Supplementary Data 1 was an earlier, incorrect version and replaces it. All frozen metadata uses the corrected workbook. The superseded file is retained only for audit and must not be used as an input.

The paper reports 5,471 Australian genomes from 1986–2022, with six source niches: human (2,996), wild animal (870), livestock (649), companion animal (375), environment (292), and food (289). The corrected Supplementary Data 1 contains exactly 5,471 identified sample rows and provides BioSample and BioProject identifiers for all rows, assembly accessions for 1,223, and SRA accessions for 4,248. The workbook yields 794 distinct entries in its ST field, including one shared novel-ST label, while the paper reports 827 STs. This discrepancy is unresolved and must be investigated before using the workbook to accept a full-study ST-count reproduction.

The source calls sample origin a `Revised_Source_Niche` and `Revised_Source_Type`; those fields are source metadata, not a uniform biological host taxonomy. “Food” and “Environmental” are source niches, not hosts. The supplement records Australian state, collection year, and source labels. It does not provide city or finer date precision, so use year-only dates and state-only geography. Missing state values remain missing. Do not derive a collection day or substitute a city.

## Declared comparison targets

The paper typed with chewBBACA 2.8.5 and the public EnteroBase `Escherichia.cgMLSTv1` scheme (accessed 19 July 2022), through Coreugate 2.0.5; its distance tree used RapidNJ 2.3.3 defaults. It treated a pairwise cgMLST distance of at most 40 alleles as a putative genomic cluster. For ten selected STs, it used SKA 1.0 and called pairs at most 100 SNPs apart cross-source genomic clusters. These were One Health relatedness screens, not assertions of direct transmission. Any ChronoClade comparison must state which metric and denominator it uses and must not treat those published cutoffs as universal thresholds.

Frozen metadata and pairwise targets are in the repository directory `validation/beta/ecoli-onehealth/`:

- `panel_st95_network_cluster_6.csv` contains 30 exact sample records in the published Supplementary Data 2 ST95 network cluster 6. `expected_published_pairs_st95_cluster_6.csv` contains its 224 reported pairwise edges at no more than 100 SNPs. The supplement resolves this panel as 15 food, 11 human, and 4 livestock sources over 2010–2016. The main-text cluster discussion has an inconsistent food count (16 rather than 15); preserve the workbook count and flag the discrepancy. The numeric published cluster ID is retained as supplied rather than relabelled as a named subclade.
- `native_panel_8_assemblies.csv` and `native_panel_expected_published_pairs.csv` declare a bounded caller pilot before any native calls: eight ST95 genomes with exact assembly accessions, five published cgMLST cluster labels (40-allele scheme), three source niches (human, food, environmental), and two states (ACT, NSW). The 28 reported pairwise ST95 SNP/allele distances are frozen. It tests published allele-profile and group relationships across source categories; it is not a locus-by-locus gold standard because the paper does not publish the underlying eight allele profiles.
- `assemblies/` contains only the eight compressed assemblies in that pilot, totaling 13,003,457 bytes. `native_assembly_manifest.json` records the exact ENA assembly and BioSample cross-check, assembly URL, checksum, year-only date, state, and source labels. No raw reads or other study genomes were downloaded.
- `freeze_manifest.json` and `frozen_counts.json` record source hashes, exact panel counts, and data limits. `freeze_metadata.py` recreates the tables from the downloaded official workbooks. `acquire_native_assemblies.py` acquires only the predeclared eight ENA FASTAs and enforces a 150 MB aggregate cap.

Full-study targets for later report comparison are the six source-niche counts and published state/year summaries; the ten cross-source STs listed in the article are ST131, ST963, ST1193, ST95, ST69, ST80, ST117, ST457, ST648, and ST57. The paper specifically reports that ST648 formed cross-source cgMLST groupings but no cross-source SKA clusters at the 100-SNP threshold. Compare only after resolving the 827-versus-794 ST count discrepancy and recording actual input coverage. The study is a convenience collection with uneven surveillance, not a prevalence sample.

## Temporal targets

The Australian study is a valid published negative temporal target. The authors report poor strict-clock signal for ST131 (`R²=0.0003`), ST1193 (`R²=0.02`), and ST95 (`R²=0.09`); their BEAST evidence for ST131 and ST95 did not converge to a stationary distribution. ChronoClade should return an explicit unsupported result for those evaluated subsets when the corresponding signal assessment agrees. This negative is a valid target, not a failed positive benchmark and not proof that all Australian E. coli lineages lack temporal signal.

Use Ludden et al., *Genomic surveillance of Escherichia coli ST131 identifies local expansion and serial replacement of subclones* (Microbial Genomics, 2020; DOI `10.1099/mgen.0.000352`) as the positive temporal candidate. The official supplement lists 794 unique run accessions and year-of-isolation values spanning 1967–2014, including 104 Irish records; the paper states 69 came from the LTCF outbreak collection. The paper reports root-to-tip `R²=0.4`, a time-scaled analysis of 794 genomes using 8,567 core SNPs after masking mobile elements, recombination, and sites missing in more than 1% of genomes, and dated targets including Clade C emergence around 1985, C2 expansion around 1994, and the C2_8 LTCF lineage dated to 2001–2005 with local expansion from 2006. The archived article method used BEAST 2.4.7; ChronoClade results should be compared as a method comparison, not expected to reproduce its dates exactly.

The positive study’s supplementary table provides exact run accessions, clade, country, year, and source type; the article’s described 69-member LTCF subset is not directly marked as a cohort column in that table. Recover a defensible accession-to-LTCF membership crosswalk and verify availability of all selected records before freezing a smaller LTCF-only input. Its supplementary PDF does not provide a machine-readable original time tree or alignment, so an exact historical reproduction remains unavailable from the frozen sources. The complete accession/year candidate is suitable for a future fresh analysis after that crosswalk and sequence/reference identity checks.

## Readiness and boundaries

Ready now: corrected Australian source metadata, the ST95 30-record target panel, and an exact-identity eight-assembly caller pilot. Pending: native caller results and any locus-level profile concordance; reconciliation of the reported ST total and one ST95 source-count difference; full-study sequence retrieval/selection; and a positive temporal accession freeze. No fitting, large genome download, or core application source modification is part of this freeze.

## Source-consistency audit of the native eight-sample panel

The exact eight-sample panel is retained unchanged, with an additional source
consistency caveat. The publisher's supplementary-file description maps the
pairwise-distance workbook (`41467_2024_55103_MOESM4_ESM.xlsx`) to official
Supplementary Data 2, which contains SNP and cgMLST distances for pairwise
comparisons in the ten detailed STs. In the predeclared ST95 panel, five of 28
published pairs have `Allelic_Distance <= 40`; four have matching corrected
Supplementary Data 1 `cgMLST cluster (40 allelic distance)` values. The fifth
pair, TEA_70_4_DC20 and TEA_70_5_R4, has published `Allelic_Distance = 0` and
`SNPs = 93`, but corrected metadata labels 260 and 895. These two cluster-label
values are also present in the superseded workbook. The pairwise supplement's
own `Cluster_*` columns both show 8 for this pair; those fields are not the
Supplementary Data 1 cgMLST-40 labels and do not resolve the discrepancy.

The audit extended this comparison across all 441,244 pair records for the ten
STs in Supplementary Data 2. Of 70,218 direct pairs at `Allelic_Distance <=
40`, 51,226 have different corrected Data 1 cgMLST-40 labels, 18,992 have the
same label, and none lack either label. Across the ST95 rows specifically,
1,697 of 30,381 are at or below 40; 440 have different labels and 1,257 have
matching labels. There are 369 zero-distance pairs across the ten STs; just
one has different labels, the ST95 pair above. This audit compares literal
pair records and sample labels; it does not reconstruct or substitute a new
clustering algorithm.

The article reports that genomes retained for analysis had at least 95% of
2,513 cgMLST loci present. The pairwise table does not provide the number of
comparable loci for each pair. Coreugate's author-maintained README says
pairwise distances ignore missing data, but the precise denominator for this
pair is not published. The correction notice only says the originally supplied
Supplementary Data 1 was an earlier incorrect version and replaces it; it does
not discuss this pairwise table or explain the cluster-label mismatch. The two
files therefore remain a conflicting source of published evidence. Preserve
the frozen distance and labels, do not relabel the samples or select the more
favorable result, and do not use these labels as validated acceptance truth.
This finding limits the benchmark's group-label comparison; it does not
invalidate the predeclared sample panel or establish which source field is
wrong.

Counts, exact pair identifiers, source checksums, per-ST results, and caveats
are frozen in
`validation/beta/ecoli-onehealth/source-consistency.json` (repository-local evidence).
The methods and correction were checked against the [publisher article](https://www.nature.com/articles/s41467-024-55103-2),
[publisher correction](https://www.nature.com/articles/s41467-025-58600-0),
[publisher supplement description](https://www.nature.com/articles/s41467-024-55103-2#Sec32),
and [Coreugate README](https://github.com/MDU-PHL/Coreugate). No native caller
was rerun and no ChronoClade fit was performed for this audit.

## Primary sources and frozen files

- [Watt et al., Nature Communications](https://www.nature.com/articles/s41467-024-55103-2) and [publisher correction](https://www.nature.com/articles/s41467-025-58600-0).
- [Ludden et al., Microbial Genomics](https://doi.org/10.1099/mgen.0.000352), [Europe PMC full text](https://europepmc.org/articles/PMC7276707), and the paper-declared ENA project [PRJEB2974](https://www.ebi.ac.uk/ena/browser/view/PRJEB2974).
- The paper-linked [APG-OHEC-Retro-M1 v1.0.1 source](https://github.com/maxlcummins/APG-OHEC-Retro-M1/tree/v1.0.1) resolves through GitHub to commit `72411c503d9187db41a60750cb995e3d5c98dcf2`; its tag response is retained and checksummed. The paper also cites the [pipelord analysis pipeline](https://github.com/maxlcummins/pipelord).
- Official source workbooks, correction page, full-text XML and Ludden supplementary files are preserved under `validation/beta/ecoli-onehealth/sources/`; file SHA256 values are in `freeze_manifest.json` (Australian source) and the source checksums below in `ECOLI_PAPER_CANDIDATES.md`.
