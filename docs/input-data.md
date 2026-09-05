# Input data

ChronoClade requires assembled genomes and a CSV metadata table. Each row names
one assembly. Paths are resolved relative to the metadata file.

```csv
sample_id,assembly,collection_date,location,species,lineage,origin,is_reference,patient_id
KPN001,assemblies/KPN001.fasta,2023-03,KIMS,Klebsiella_pneumoniae,ST15,local,true,P001
KPN002,assemblies/KPN002.fasta,2024-01-17,KIMS,Klebsiella_pneumoniae,ST15,local,false,P002
KPN_REF,context/KPN_REF.fasta,2021,Philippines,Klebsiella_pneumoniae,ST15,context,false,
```

## Required columns

| Column | Content |
| --- | --- |
| `sample_id` | Unique tree-safe identifier using letters, numbers, `.`, `_` or `-` |
| `assembly` | FASTA or compressed FASTA path |
| `collection_date` | `YYYY`, `YYYY-MM`, `YYYY-MM-DD`, `YYYY-XX-XX`, or a TreeTime numeric/range date |
| `location` | Hospital, region or country used for exploratory state reconstruction |
| `species` | Confirmed species or species-complex label |
| `lineage` | Analysis unit, usually an ST or genomic cluster |
| `origin` | `local`, `retrospective` or `context` |

## Optional columns

`is_reference` may mark one assembly per lineage with `true`. If it is absent or
all values are false, ChronoClade selects the assembly with the highest N50.

`patient_id` supports summaries that retain one isolate per patient. Use a coded
identifier. The analysis can run with missing patient identifiers. In that case,
ChronoClade cannot assess repeated sampling at patient level.

## Choosing the analysis unit

Do not mix distant bacterial populations in one clock analysis. An MLST
sequence type is a useful starting point when it captures the clone under
investigation. A species-wide collection, a polyphyletic ST, or a lineage with
deep internal structure may need to be divided before temporal analysis.

Keep focal, retrospective and public context genomes in the same metadata file.
ChronoClade uses `origin` to separate the locally sampled collection from the
context used to interpret it.

## Dates

Use the most precise verified date available. Month-level sampling is sufficient
for many two-year surveillance collections. Year-only dates are accepted, but
the wider uncertainty should temper the interpretation of short branches and
recent ancestral dates.

Validate the table before starting expensive phylogenetic steps:

```bash
pixi run chronoclade validate metadata.csv --min-samples 10
```

Use de-identified metadata. Do not commit patient-level data or unpublished
assemblies to a public repository.
