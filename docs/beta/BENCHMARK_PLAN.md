# Preregistered beta benchmark targets

Date: 17 October 2026. Status: targets declared; BIO-TIME-01 metadata frozen,
sequence/reference and exact historical alignment recovery pending.
Targets below are declared before any new biological benchmark fit. Existing
CCRE ST258 outputs are engineering/sensitivity evidence, not a reproduction of
these publications.

## Tests and allowed conclusions

| ID | Source | Target | Independent acceptance evidence |
| --- | --- | --- | --- |
| BIO-KP-01 | David et al.,2019; https://pmc.ncbi.nlm.nih.gov/articles/PMC7244338/ | Closest genomic neighbours and local/hospital versus wider context | Exact study accessions; verify chosen ST and source site/country; reproduce selected published nearest-neighbour/geographic findings with declared denominator; no importation counts inferred |
| BIO-KP-02 | International and regional spread,2024; https://www.nature.com/articles/s41467-024-49349-z | ST258/512 structure and context selection across surveys | Freeze actual analysed subset and source reference before comparison; compare selected clades and nearest relationships, retaining distinct MLST scope unless explicit combined-study contract |
| BIO-KP-03 | Pavia study,2022; https://www.nature.com/articles/s41564-022-01263-0 | Human/animal/environment context summaries | Freeze source host categories and locations; retain unknown/clinical ascertainment bias; do not interpret host proportions as population prevalence |
| BIO-TIME-01 | Bowers et al.,2015; https://doi.org/10.1371/journal.pone.0133727 | Longitudinal ST258 inference as primary positive candidate | Paper used101ST258 genomes for dated analysis and1425coreSNPs after an explicit215kb mask. Metadata freeze resolves 167 study records and 101 combined-lineage candidates (97 ST258 + 3 ST512 + 1 ST1199); recover reference sequence/mask and exact historical alignment before running; inspect temporal support, important clades, fitted rate and root uncertainty. Different inference/filtering is a method comparison, not exact reproduction |
| BIO-TIME-02 | Frozen date-shuffled and date-homogeneous derivative of a verified positive dataset | Negative controls for temporal gate | Dates generated with preregistered seed before fitting; preserve exact genotypes/tree and document derivative status. Negative controls are statistical software checks, not independent observed biological datasets |
| BIO-TIME-03 | Independently justified biological dataset with insufficient temporal information | Legitimate unsupported dating | Select using sampling design/literature before fitting, not because ChronoClade returned a lowR². Small time span does not guarantee no biological signal; record actual conclusion |

## Input freeze contract

For each chosen study, save source citation, metadata/supplement checksums,
exact accessions, source role/date/location/host labels, exclusions and reasons,
reference identity, and available original alignment/tree. Download only after
identity resolution. A strain-name match alone is not an exact sequence identity.
Do not substitute a different assembly of the same BioSample without recording
that difference.

Before a run, declare retained study subset, expected labelled clades/relationships,
metrics and tolerances appropriate to method differences. Use separate known
related/unrelated labels for any threshold validation. Do not select tolerances
post hoc to obtain a pass. Missing original inputs or inaccessible supplements
must be recorded; illustrative reruns are not promoted to reproduction claims.

## Temporal interpretation

A published BEAST result is not automatically a gold-standard proof of temporal
signal. ChronoClade uses a different topology/filtering/dating pipeline and must
run its own documented screen. Check numerical units and ascertainment/site counts.
Compare ancestral dates for the same descendant set, not whole trees with different
ancestors. Model/date intervals, topology uncertainty and subsampling sensitivity
are reported separately. The first target is a candidate, not a promised pass.

## Resource and walkthrough acceptance

The frozen small profile fixture under tmp/beta-usability provides a regression
and file-handoff check only; it cannot establish native typing correctness or
biological validation. Small native positive/negative tests and a checked
published-data subset will provide separate evidence. Existing7815-profileST258
and synthetic10k runs establish scale candidates; record fresh export/report
resource measurements after the first fixes without duplicating all matrices.

## BIO-TIME-01 metadata freeze: Bowers 2015

Completed locally before any fit. Reproduction script and evidence are in
`validation/beta/bowers2015/`; read its `README.md` and `freeze_manifest.json`.

| Item | Frozen result | Remaining gap |
| --- | --- | --- |
| Official S2 isolate table | 167 isolates, source XLSX row and historical MLST retained | No accession/date columns in the supplement itself |
| Official repository crosswalk | All 167 uniquely matched to 167 runs / 167 BioSamples using exact ENA `sample_alias` within article-declared PRJNA252957 | This is a metadata link; sequence checksums/assembly versions not yet frozen |
| Date provenance | ENA collection-date field, all 167 year precision; full set 1996–2012 | Never parse years from isolate labels or invent observation days |
| Published combined-lineage candidate | 101 genomes: 97 ST258 + 3 ST512 + 1 ST1199; 2003–2012 | Original BEAST tip manifest/reference inclusion still unresolved; strict ST258-only would omit 4 candidates |
| Metadata exclusions | None: 167/167 records resolved | Candidate analysis scope excludes 66 other/novel MLST records by declared combined-lineage criterion, not by a fitted result |
| Reference | Published NJST258_1 accession CP006923 recorded | Exact accession version/checksum and 215 kb mask coordinates unresolved; no reference sequence downloaded |
| Original analysis artifacts | Official S1 expanded phylogeny PDF frozen | Figure is not Newick or the 101-genome BEAST tree; machine-readable historical SNP alignment/tree not listed in PLOS supplements |

Files:

- `sources/S2_isolates.xlsx`, `sources/ena_PRJNA252957.tsv`,
  `sources/article.html`, `sources/S1_published_phylogeny.pdf`.
- `resolved_accessions_dates.csv`: 167 rows, repository links and year intervals.
- `temporal_candidate_101.csv`: 101-row combined-lineage candidate, retaining actual MLST labels.
- `freeze_metadata.py`: standard-library downloader/parser with pinned source checksums.
- `freeze_manifest.json`: exact source/output checksums and unresolved prerequisites.

A repeated cached freeze regenerated identical derived-table/manifest checksums.
No read/assembly/reference sequences were downloaded, no model fitted and no
biological benchmark pass claimed. The sequence/mask/alignment gaps must be
resolved or explicitly change the goal from historical reproduction to a method
comparison before running. No acceptance tolerance has been selected from a fit.
