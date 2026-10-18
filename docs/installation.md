# Installation

## Requirements

ChronoClade runs on macOS and Linux. Windows users can use WSL2. Pixi installs
the Python package and the compiled phylogenetic programs from the checked-in
lock file.

Install [Pixi](https://pixi.sh/latest/), then clone the repository:

```bash
git clone https://github.com/ghruproject/chronoclade.git
cd chronoclade
pixi install
```

Check the installed programs:

```bash
pixi run chronoclade preflight
```

The table should report paths for `ska`, `iqtree`, `ClonalFrameML` and
`treetime`. Run ChronoClade through `pixi run` so these programs remain on the
same executable path.

## Separate command environments

For existing profiles or provider typing, `pixi run -e prepare chronoclade
prepare ...` avoids the native SNP-tree stack. The lean environment does not
include BLAST. When preparation will call cgMLST from assemblies, use
`pixi run --frozen --environment default chronoclade prepare ...`; this activates
the installed `blastn` and `makeblastdb` executables required by the caller.
Use the same full environment when missing provider typing may trigger native
calling. Calling the environment's `chronoclade` executable by absolute path
does not activate its external tools.

`pixi run -e cgmlst chronoclade cgmlst ...` analyses existing profiles.
Use `-e tree` and `-e time` for those stages. The default environment includes
development and verification tools.

Check the intended stage with `pixi run -e prepare chronoclade preflight --stage
prepare` (or `cgmlst`, `tree`, `time`). This checks stage software only; it does
not claim that typing references, credentials or input assemblies are ready.
Without `--stage`, preflight checks the complete native tool suite.

On Linux, the cgMLST environment includes RapidNJ. On Apple Silicon, install the
pinned macOS backend with `pixi run chronoclade setup-scale`; the installer records
its revision, binary checksum and compiler. The upstream backend currently runs
through Rosetta on arm64.

## Updating an installation

```bash
git pull
pixi install
pixi run chronoclade version
```

Pixi will reuse downloaded packages where possible. Changes to `pixi.lock`
alter the resolved environment and should be reviewed like source-code changes.

## Build the documentation

```bash
pixi run docs
```

The static site is written to `site/`. Use `pixi run docs-serve` while editing
Markdown files.

## Pathogenwatch credentials

Public metadata discovery is anonymous and explicitly excludes private genomes.
Optional native cgMLST/LIN/HierCC typing is installed with
`pixi run chronoclade setup-typing`. It uses uv environments for Python tools and
a pinned Node runtime for the cgMLST caller. Reference databases have separate
readiness and authentication requirements; see [query typing](query-typing.md).
Assembly calling also requires BLAST on the executable path; use the full Pixi
environment described above. The pinned Node runtime alone does not supply it.
For current Pasteur references, follow the [one-time browser authorisation](pasteur-auth.md)
before building with `setup-typing --prepare-db --pasteur-secrets KEYFILE`.

Authenticated sequence/cgLIN downloads use `PATHOGENWATCH_API_KEY`, or
`~/.config/chronoclade/pathogenwatch.json` containing an `api_key` field with
permissions `0600`. Never commit this file. Offline fixture tests need no
credentials. The same Pixi CLI runs on macOS/Linux laptops and Linux SLURM
workers; SLURM is optional. See the [independent workflow](independent-suite.md#local-and-slurm-jobs).
