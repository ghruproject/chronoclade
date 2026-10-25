# Compressed storage and reuse

ChronoClade can keep large reusable inputs and intermediate files compressed at
rest. Compression is lossless: scientific paths, manifest checksums, dataset IDs
and native command fingerprints continue to describe the original bytes.
Completed analyses remain reusable without downloading assemblies or rerunning
native tools just because their files were compressed.

## Enable compression for new work

Put the global option **before** the command:

```bash
pixi run chronoclade --storage compact prepare profiles.json --out prepared
pixi run chronoclade --storage compact cgmlst prepared/dataset.json --out cgmlst
pixi run chronoclade --storage compact tree selections/ensemble.json --out tree
pixi run chronoclade --storage compact time tree/tree_ensemble.json --out time
```

Alternatively, `export CHRONOCLADE_STORAGE=compact` enables this policy for CLI
commands in that shell. The default is `plain`; both policies can read managed
compressed inputs. Successful stage commands compact their own output directory;
cgMLST context acquisition also compacts its completed context bundle and
retrieval checkpoints. Failed commands keep their partial outputs available for
inspection and recovery.

Generated jobs preserve the choice:

```bash
pixi run chronoclade --storage compact slurm-ensemble tree selections/ensemble.json \
  --out tree-jobs --results tree-results
pixi run chronoclade --storage compact slurm-prepare assemblies.csv \
  --typing-config typing.json --out prepare-jobs --results prepared
pixi run chronoclade --storage compact job tree --out tree-job --slurm -- \
  selections/selection-001.json --out tree
```

Array plans record `storage`; workers use that saved setting even if the submit
shell changes. Older plans default to plain storage. Each successful worker
compacts only its own result. Collectors read checked compressed workers and
release temporary restorations when finished. No storage operation submits jobs.
The shared assembly download cache is retained; compact it separately when idle.

## Compact existing outputs

Wait until all writers using the directory have stopped. In particular, do not
compact inputs while older ChronoClade installations or external tools are reading
them. A running `job.json` prevents packing, but this check cannot identify every
external process or scheduler job.

```bash
# Preview eligible files without modifying them.
pixi run chronoclade storage pack chronoclade_results --dry-run

# Compress large reusable payloads with large-window Zstandard.
pixi run chronoclade storage pack chronoclade_results
pixi run chronoclade storage status chronoclade_results
pixi run chronoclade storage verify chronoclade_results

# Optional gzip format; use smaller files too if desired.
pixi run chronoclade storage pack another-idle-run --codec gzip --min-size-mb 0.1
```

The default threshold is 1 MiB, and files must shrink by at least 5%. Eligible
payloads are FASTA (`.fasta`, `.fa`, `.fna`), PHYLIP (`.phy`), SKA (`.skf`), NumPy
(`.npy`), JSONL (`.jsonl`) and the named provider checkpoints `public_N.json`,
`matched_N.json`, `input.json` and `collection_snapshot.json`. Reports, control
manifests, generic JSON indexes, CSV/TSV exports, logs and reference databases
remain directly accessible. Directory and file symlinks are not traversed.

Packing does not discard evidence or remove failed runs. A completed archive is
fully decompressed and checked against the source before the original file is
removed. Files that change during packing are retained and reported as errors.
Repeated packing reuses an unchanged, verified archive. Switching `--codec` does
not transcode files already packed.

## Reading, restoring and exporting

Managed compression uses adjacent files:

```text
core_alignment.fasta.cc.zst          # or .cc.gz
core_alignment.fasta.cc-storage.json # storage receipt
core_alignment.fasta.cc-lock         # small advisory lock
```

Keep the archive and receipt together. The receipt records the logical name,
codec, original bytes and SHA256, stored bytes and SHA256, and ordinary file
permissions. Scientific manifests still refer to `core_alignment.fasta` and its
original checksum; they are not rewritten by packing.

Text readers stream the compressed data. NumPy memory maps and native programs
need actual files, so ChronoClade restores exact bytes atomically at the original
path. CLI sessions hold read leases for those consumers and remove unchanged
**newly restored** copies at exit. A pre-existing restored file is retained unless
a successful compact policy explicitly packs that output. Changed copies are
never silently removed. Compressed downloads can be copied directly into an
independent job snapshot without expanding the shared cache. No mutable
hardlinks are introduced.

Storage verification checks both compressed and logical checksums; it also
checks a materialized copy if one exists. Corrupt archives, malformed receipts,
unsafe destinations and checksum mismatches fail explicitly. Interrupted packing
can leave an unreferenced archive or temporary file, but the original is kept
until the archive and receipt are committed. Investigate leftovers before removal.
Advisory locks coordinate cooperating readers and writers; they do not stop
external programs that ignore them.

For external tools or older ChronoClade versions, restore the files first:

```bash
pixi run chronoclade storage restore chronoclade_results
```

Restore retains archives, so it temporarily uses **both** sizes on disk. Run
`storage pack` again to evict eligible unchanged copies. Explicit user FASTA
`.gz` inputs can be streamed by the alignment reader; arbitrary renamed archives
are not a substitute for managed storage receipts.

## Resource and compatibility notes

Zstandard uses level 3 with a 128 MiB window and long-distance matching to exploit
similarity across aligned sequences. gzip uses level 1. Compression is sequential
and may need roughly a window of memory plus codec overhead. Native/mmap work
needs free disk space for restored inputs, and packing briefly needs space for
both the source and archive. Interrupted processes can leave restored copies;
verification followed by packing reclaims eligible copies safely.

Reuse avoids expensive inference and downloads, but validation and decompression
still cost I/O and CPU. Compression does not change analysis thresholds, selected
genomes, outlier rules, dates or results. Representation changes alone do not
invalidate hashes. Normal software-version/provenance checks still apply: a
pipeline upgrade can invalidate computation caches independently of compression.
Do not change the software beneath active jobs.

Python callers can use `chronoclade.storage.StorageSession()` around an analysis
and `storage.pack(output)` after successful completion. Without a session,
materialized files remain available until explicitly packed; keep the session
open for the lifetime of any native consumers or memory maps.
