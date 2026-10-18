# Publication scope review — 18 October 2026

Read-only Git audit of `codex/st258-native-e2e`. No staging, commit, push, fetch, checkout, deletion or source edit performed by this reviewer.

## Recommended publication scope

- Include all 91 currently modified tracked files: implementation, regression tests, README, user documentation and MkDocs configuration. The tracked diff is 6,283 insertions and 919 deletions. Parent review should remove the E. coli acceptance links from public MkDocs navigation to respect the parked scope; preserve the E. coli implementation and regression tests.
- Include all new `chronoclade/**/*.py` modules. Existing changed modules import them; committing only tracked changes would produce a broken checkout.
- Include all new `tests/test_*.py` and all eight files in `tests/fixtures/`. The fixtures total about 0.01 MiB. The only binary among the proposed small untracked publication set is the small gzipped provider CSV fixture.
- Include `examples/beta/{README.md,profiles.json,catalogues.json}`, `docs/mac-klebsiella-walkthrough.md`, `docs/pasteur-auth.md`, and the beta documentation needed for existing cross-links. Keep parked E. coli audit documents out of public navigation/marketing claims. Retaining their historical text in the repository preserves evidence; do not describe it as release acceptance.
- `uv.lock` is a small, portable PyPI lock (~0.04 MiB), with editable source `.` and public registry references. Include it if the documented uv installation workflow is intended to be supported; it contains no local absolute source dependency.
- These untracked source/tests/docs/examples/lock files together total 111 files and 847,329 bytes before any parent documentation adjustments.

## Exclude from this push

- All `validation/beta/`: 7,024 untracked files / 281.28 MiB at audit time. This directory contains generated HTML reports, raw provider snapshots, repeated staging copies of source trees, execution logs, local authorization manifests, and study supplements. Largest file is a 53.5 MiB XLSX; two repeated fixture tarballs are 13.06 MiB each. Do not use `git add .`.
- If public evidence receipts are required, curate compact final summaries with path/hash/status evidence into a separate documented selection after inspecting them. The existing documentation references many local validation paths. Excluding bulk validation means those are local audit references, not files obtainable from a fresh clone; label that limitation or supply selected receipts rather than silently promising downloadable evidence.
- `runtime-provenance.json`: generated single-runtime evidence containing absolute local checkout/module paths, not a portable project configuration.
- `.impeccable/critique/`: internal review artifact, unrelated to runtime or user documentation.
- Already ignored `tmp/`, `.pixi/`, caches, build/site outputs, assemblies and previous result directories remain excluded. Existing ignore rules do not cover validation/beta wholesale, credential filenames, or runtime-provenance.json; explicit staging exclusions are necessary. Adding narrow ignore rules for generated validation/runtime provenance would prevent accidental future inclusion, without deleting local evidence.

## Credential and dependency review

A filename-safe scan of tracked and untracked text found no PEM/OpenSSH private key blocks and no GitHub token patterns. Credential-value assignment patterns appeared only in `tests/test_typing_setup.py` and its four generated validation staging copies. Inspection of the source test confirms clearly synthetic fixture values. No credential values were printed by the scan. This is a targeted pattern scan, not a claim that arbitrary secrets can be detected reliably. Keep all runtime credential caches and provider authorization material local.

The installed ChronoClade executable responds to `--help`. Installed requirements of ChronoClade, Biopython, Matplotlib, NumPy, Rich and Typer satisfy their package metadata constraints. Static internal import resolution found no missing ChronoClade modules when all new source files are included. `plincer` and `toml` imports in `pasteur_assignment_worker.py` occur inside `main()` for the separately installed plincer environment; they are not missing core dependencies. `python -m chronoclade` is not a supported entry point (there is no `__main__.py`); use the installed `chronoclade` command as documented.

No push blocker found in the proposed source + docs + small-fixture scope. Scientific acceptance claims and complete tests remain the parent reviewers' responsibility. Before committing, inspect the staged file list and staged size, verify no validation/runtime/credential paths entered it, and run staged diff whitespace checks.
