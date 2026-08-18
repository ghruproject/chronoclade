# Bundled AllTheBacteria context metadata

`atb_context_202505.parquet` is a compact, analysis-ready subset of the
AllTheBacteria 2025-05 metadata release. It contains high-quality downloadable
assemblies assigned to these species/scheme pairs:

- *Escherichia coli* / `ecoli_achtman_4`
- *Klebsiella pneumoniae* / `klebsiella`

The fields are limited to those used for same-ST discovery, metadata filtering,
context reporting and assembly retrieval. `atb_context_202505.json` records the
source URLs, selection rule, row counts and SHA-256 checksum.

The source data are from the public AllTheBacteria/ENA dataset and are
distributed under the MIT licence. Cite:

> Hunt M, Lima L, Anderson D, et al. AllTheBacteria - all bacterial genomes
> assembled, available, and searchable. bioRxiv 2024.03.08.584059.

Regenerate the snapshot with:

```bash
pixi run python scripts/build_atb_context_snapshot.py \
  --sqlite ~/.atbfetcher/atb.metadata.202505.sqlite \
  --mlst ~/.atbfetcher/mlst.parquet \
  --output beyondmlst/data/atb_context_202505.parquet \
  --manifest beyondmlst/data/atb_context_202505.json
```
