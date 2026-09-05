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
