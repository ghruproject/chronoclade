# Bounded Klebsiella relatedness source review

Reviewed 2026-10-18. Source selection and label extraction only: no allele profiles, LIN assignments, measured distances, or newly inferred trees were inspected. The proposal below precedes those measurements.

## Recommendation and evidential ceiling

Use **Zhou2017 as the first tractable comparison**, and reserve **NIH2011 / Snitkin2012 for an external, limited positive-control follow-up**. These studies do not together constitute a clean independent binary transmission benchmark. Report agreement with author epidemiological categories, and state their ascertainment limits. Do not call a preferred LIN level a validated transmission boundary from this evidence.

### Development panel: Zhou et al. 2017

[Primary article, including Table 1 and methods](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.00371/full), DOI `10.3389/fmicb.2017.00371`.

Table 1 identifies 26 outbreak isolates and five sporadic controls, all author-reported ST11. Methods define sporadic cases as having no direct epidemiological relevance to other cases; outbreak assignment used epidemiological information **and PFGE/MLST**. Thus labels predate the new cgMLST measurements but are not purely independent clinical transmission truth. Patient-level links are unavailable. Controls are A1502, A1732, A1760, A1763, A1771. The outbreak occurred in Beijing SICU in October–December 2011. Its original scheme contained 1,143 targets; its numerical cutoff must not transfer to the native 629-locus scheme. [Source](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.00371/full).

**Verified archive crosswalk:** all 31 exact strain names resolved to versioned ENA assemblies. Eight belong to `PRJNA313004`; 23 have other individual BioProjects. The case-only alias `TR231-M` → `TR231-m` is explicit. The paper's accession-range text is insufficient for automated acquisition; use the individual verified manifest. Archive sequence identity and quality remain to be checked after downloading.

Frozen material is in the local `validation/science/zhou2017/` directory (start with its `README.md`): full article HTML, lossless extracted Table 1 rows, ENA query result, 31 assembly XML records, checksums, and an enriched manifest with exact author labels, raw dates, normalized day dates, hospital pseudonyms, author ST/PFGE, sample/project/accession joins, and FASTA URLs. Patient IDs remain `unknown`. These local evidence files are not bundled with the documentation site.

### Reserved follow-up: Snitkin et al. 2012 / NIH2011

[Primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC3521604/), DOI `10.1126/scitranslmed.3004129`; [NCBI umbrella project](https://www.ncbi.nlm.nih.gov/bioproject/188534).

The study describes 18 affected patients and seven isolates from four body sites of patient 1. Same-patient isolate identity would provide independently recorded within-host positives, once sample-to-patient mapping is verified. The transmission map and major clusters were inferred using genetic information with patient traces; they are not independent labels. Absence of ward overlap is not proof of absence of transmission. No explicit independent negative-control set was established in this review. [Source](https://pmc.ncbi.nlm.nih.gov/articles/PMC3521604/).

[Current candidate assemblies](snitkin2012-current-assembly-candidates.csv) list 20 GenBank assemblies returned by the [live NCBI report](https://api.ncbi.nlm.nih.gov/datasets/v2/genome/bioproject/PRJNA188534/dataset_report?page_size=1000), after excluding paired RefSeq duplicates. KPNIH1 and KPNIH10 now resolve to version 2 complete genomes; historical paper WGS versions must be pinned separately. Do not equate strain numbering with patient numbering. The paper lists 20 WGS accessions, while later reuse may contain additional assemblies; this crosswalk does not resolve all later reuse.

Supplement: `NIHMS424107-supplement-Sup_Materials.pdf`, linked by the primary article. Tables S1/S2 and patient traces are the intended sample-mapping evidence. Download attempts encountered PMC download-challenge HTML and publisher HTTP 403; supplement contents were not verified. **Release gate:** no NIH patient-pair labels until that mapping and exact assembly versions are frozen. Treat all candidate pair labels as unknown now. This is a reserved external sensitivity check, not a completed held-out validation set.

## Proposed label policy before measurements

These are proposed analysis rules, not additional facts asserted by either paper.

| Category | Frozen inclusion rule | Interpretation |
| --- | --- | --- |
| Zhou author outbreak co-membership | Both isolates carry Table 1 `Outbreak` | Study outbreak association, not direct transmission |
| Zhou outbreak-versus-sporadic comparison | One Table 1 `Outbreak`, one `Sporadic` | Author-reported lack of direct epidemiological association, with incomplete ascertainment |
| Zhou sporadic-versus-sporadic | Keep separate and report individually | Optional descriptive controls; do not inflate the main negative sample count |
| NIH within-host | Same verified patient ID from source mapping | Within-host relatedness only; cannot estimate transmission specificity |
| NIH interpatient | Unknown for independent binary evaluation | Genetically informed published reconstructions can be annotated separately |
| Other genomes, including Bowers context | Unknown unless explicit source evidence is separately curated | Shared location/date/ST/clade never supplies an outbreak label |

For 26 and five retained Zhou isolates, the first two strata contain 325 and 130 unordered pairs. Those are dependent pairs from one outbreak, not 455 independent observations. Any uncertainty analysis must respect isolate/study dependence; an apparent high pairwise score cannot support broad generalization.

Freeze inclusion/QC rules and the full LIN-level comparison grid before typing. Recompute all measured allele distances using one pinned native scheme, record jointly callable loci and missingness, and retain denominator information. Evaluate LIN level co-membership and measured pairwise distance separately: a shared hierarchical level need not impose a bound on every pair's distance. Report every level's results before naming any preferred level.

Construct SNP trees as a separate descriptive comparison after the manifest freeze, with fixed reference, masks, recombination handling, and missingness rules recorded. Annotate the frozen author categories. Tree clustering cannot retroactively create positive/negative labels, and matching the published tree cannot independently validate a transmission decision. Do not tune on NIH follow-up results and still describe it as held out.
