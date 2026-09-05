# Real-data validation: *Staphylococcus aureus* ST239

This validation uses the longitudinal ST239 collection from Baines et al.
(2015), *Convergent Adaptation in the Dominant Global Hospital Clone ST239 of
Methicillin-Resistant Staphylococcus aureus*.

The published study analysed 123 genomes sampled between 1980 and 2012. It
reported two distinct Australian ST239 clades, consistent with long-term local
expansion plus the introduction of another clade. The study estimated a mean
substitution rate of 1.6 × 10⁻⁶ substitutions per site per year and dated the
emergence of ST239 to 1946.

## Dataset used here

The ChronoClade run used the 72 study samples that had a high-quality ST239
assembly in the bundled AllTheBacteria 2025-05 snapshot:

- 69 Australian isolates;
- one New Zealand isolate;
- two external comparators from Turkey and Hungary;
- 17 distinct collection years spanning 1980–2012; and
- 2,462,750 complete sites after reference mapping.

`input_manifest.csv` records the public sample accession, original isolate
name, collection year, country, analysis role and assembly URL. Assemblies are
not committed to the repository.

## ChronoClade result

The complete native SKA2 → IQ-TREE → ClonalFrameML → TreeTime workflow was run
with 100 date randomisations.

| Result | ChronoClade | Published analysis |
|---|---:|---:|
| Root-to-tip R² | 0.94 | 0.8261 for the full collection |
| Clock rate | 1.693 × 10⁻⁶/site/year | 1.6 × 10⁻⁶/site/year |
| Estimated ST239 root | 1946 | 1946 |
| Date-randomisation p-value | 0.0099 | Not reported in this study |
| Working scenario | Persistence plus additional introductions | Two independently circulating Australian clades |

The workflow identified two longitudinal focal groups. Their median
recombination-filtered distance was 124 SNPs within groups and 454 SNPs between
groups. It therefore produced the cautious working interpretation “consistent
with local persistence plus additional introductions”. Confidence remained low
because this initial run included only two external contextual genomes and did
not include coded patient identifiers or corrected-tree branch support.

This is strong validation of the temporal-signal gate and time-tree estimates,
but it is not yet a complete validation of contextual-neighbour retrieval.

## Preserved outputs

The compact `results/` bundle contains the stage-by-stage HTML report, temporal
diagnostics, dated tree, recombination-corrected tree, distance matrices and
machine-readable evidence. Figures are preserved as SVG and high-resolution PNG,
the temporal statistics and node dates are also available as CSV, and
`supporting_results.zip` packages the reader-facing evidence in one file. Large
assemblies and intermediate alignments are deliberately excluded.

## Sources

- Baines SL et al. 2015. <https://doi.org/10.1128/mBio.00080-15>
- ENA study ERP009308 / PRJEB8247: <https://www.ebi.ac.uk/ena/browser/view/PRJEB8247>
- Harris SR et al. 2010 ST239 study: <https://doi.org/10.1126/science.1182395>
- BactDating ST239 vignette: <https://xavierdidelot.github.io/BactDating/articles/Staph.html>
