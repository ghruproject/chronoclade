# Installation and usability acceptance evidence

Date: 17 October 2026. Scope: one isolated Apple Silicon install and the
synthetic profile walkthrough. No human usability sessions are represented here.

## Route under test

The recommended route remains Pixi from a checkout. On macOS and Linux, install
Pixi, clone ChronoClade, then run `pixi install`; use the named `prepare`,
`cgmlst`, `tree` and `time` environments when following a staged workflow. Run
the CLI through `pixi run` so its Python and native tools come from the same
environment. Check software for the stage you intend to run with
`chronoclade preflight --stage prepare|cgmlst|tree|time`; `--stage all` checks
the complete native stack. These are software checks. Keep typing
references/databases, credentials and input readiness as separate setup checks;
an offline frozen-profile analysis does not require them. The current
[installation guide](../installation.md) describes these boundaries but does
not yet give one tested clean-install route for both Apple Silicon macOS and
Linux.

This is the same basic installation pattern used by the local BactScout
checkout: Pixi is the documented route, and the README provides short task
commands. BactScout additionally combines preflight checks with database and
resource checks, writes per-sample logs and a consolidated summary, and makes
WARN/FAIL thresholds configurable. ChronoClade already has named environments,
stage manifests and output READMEs; a single preflight/readiness view, consistent
stage/count/elapsed-time progress, durable user-facing failure guidance and an
end-of-run next-command summary remain usability gaps. This is a workflow
comparison, not a claim that the scientific tools or outputs are equivalent.

## Fresh Apple Silicon run

The source and fixtures were copied to a new `/tmp/chronoclade-beta-macos-*`
directory with `.pixi`, prior `tmp` outputs, Git metadata and caches excluded.
Pixi then created new environments there. Host: macOS 14.6.1, Darwin arm64;
Pixi 0.72.2; environment Python 3.11.17; package version 0.1.0. The copied
workspace's checked-in `pixi.lock` was used with `--frozen`.

After that fresh environment install, the frozen temporary source snapshot was
updated by copying only `chronoclade/cli.py` and `chronoclade/prepare_cli.py`
from the current worktree. The environment was not rebuilt. Their source and
temporary-copy SHA256 values match in `source-snapshot-update.txt`; the
additional preflight and prepare observations below therefore verify the
current two-file change against the newly installed dependency environment.

```sh
pixi install --frozen -e prepare
pixi install --frozen -e cgmlst
pixi install --frozen
pixi run chronoclade version
pixi run chronoclade preflight --stage all
pixi run -e prepare chronoclade preflight --stage prepare
pixi run -e cgmlst chronoclade preflight --stage cgmlst
pixi run -e prepare chronoclade prepare examples/beta/profiles.json \
  --input-kind profiles --catalogues examples/beta/catalogues.json \
  --out tmp/install-acceptance/prepared
pixi run -e cgmlst chronoclade cgmlst \
  tmp/install-acceptance/prepared/dataset.json \
  --out tmp/install-acceptance/cgmlst --lin-level 5 --subsamples 1 \
  --context-size 0 --bootstrap-replicates 0
pixi run chronoclade export tmp/install-acceptance/cgmlst/cgmlst.json \
  --viewer grapetree --out tmp/install-acceptance/grapetree
```

All install, stage preflight and analysis commands exited successfully. The
`prepare` and `cgmlst` preflights passed inside their lean stage environments
without reporting absent native SNP-tree programs. `cgmlst` notes that
RapidNJ is additionally required for large profile analyses. The full default
environment's `--stage all` check reported `ska`, `iqtree`, `ClonalFrameML`,
`treetime`, the profile backend and native platform as available/supported.
Each preflight states it checks software only. The fixture contains three
invented Klebsiella profiles across three invented loci, with one input and two
context records. Preparation reported `readiness: ready`, and cgMLST completed
one block. The GrapeTree bundle contains three profile rows and three matching
metadata rows with shared `CC_...` IDs.

Readiness has a narrow scope: the fixture's frozen profiles and lineage values
are usable by cgMLST. Its readiness rows still say `assembly_status: unavailable`,
and the dataset lacks several geographic/source fields. The status does not
authorize assembly-tree or time-tree commands. The updated preparation command
prints: “Preparation readiness describes frozen typing and metadata; tree
assembly availability and temporal signal are assessed in later stages.” The
prepare README gives the cgMLST command; the cgMLST README links to the group
report. The export README explains IDs and per-locus allele recoding. No
tree/time run was attempted on a profile-only fixture.

The missing-input check used the same installed CLI with a nonexistent profile
file. It exited 2 and printed `Error: Cannot read preparation JSON` with the
missing path and OS reason. The updated stage preflights passed for `prepare`
and `cgmlst`, and the full `--stage all` preflight passed in the default
environment. Neither command inspects Pathogenwatch credentials or
reference-database readiness; those have separate setup commands and reports.

The small command outputs and error checks are preserved in
the repository directory `validation/beta/install-macos/`.
The successful report/export outputs remain in the temporary checkout named in
`checkout-path.txt`. Reproduce with the commands above from a clean temporary
copy; this does not modify the shared source checkout's Pixi environment.

## Linux status

No Linux fixture run was possible in this acceptance pass. Docker is installed,
but its configured API socket had no daemon. Podman 5.2.2 had two existing
AppleHV machines configured (4 CPUs/8 GiB each); `podman machine start
podman-machine-default` reported success after a slow boot, but the machine
returned to `stopped` and its socket refused connections before a container
could run. Startup created/refreshed that machine's configured sparse backing
file: 100 GiB logical size, 5.3 GiB allocated on disk. The machine remains
stopped, its backing file is preserved, and the other existing machine was not
started. No new VM was configured.

Concrete prerequisite for the remaining check: a running Linux x86_64 host or
an available AMD64-capable Docker/Podman runtime. In that environment, repeat
the Pixi installation and the exact fixture commands above, record `uname -m`,
Pixi/package versions and `preflight`, and inspect the generated readiness and
GrapeTree bundle. `pixi install --platform linux-64` on this Mac would only
resolve/install a target environment; it cannot substitute for executing the
commands on Linux.

## Acceptance decision

| Requirement | Status | Evidence / remaining work |
| --- | --- | --- |
| Fresh Apple Silicon Pixi install | Verified for current locked source | Fresh temporary checkout, separate stage environments and full default environment installed; logs in `install-macos/`. |
| Small profile workflow and export | Verified as an installation/navigation smoke test | `prepare → cgmlst → GrapeTree export`; three invented records only. No biological validity claim. |
| Stage software preflight | Verified on Apple Silicon | `prepare` and `cgmlst` pass in their own environments; `all` passes in default. Evidence logs record the current CLI source refresh. |
| CLI error and next action | Partially verified | Missing-input exit 2 is readable; preparation now states readiness scope and stage README gives the next command. Consistent progress and end-of-run summary remain open. |
| Separate tool and reference readiness | Partial | Software preflight explicitly excludes references, credentials and input readiness. Native typing databases and credentials were not installed or tested in this fresh environment. |
| Linux x86_64 install and fixture | Blocked by runtime availability | Requires a running x86_64 Linux host or healthy AMD64-capable container runtime. |
| Two independent human walkthroughs | Not done; required beta gate remains open | Use [HUMAN_WALKTHROUGH.md](HUMAN_WALKTHROUGH.md). Agent execution does not count. |

## Frozen-source wheel install on macOS arm64

On 17 October 2026, built `chronoclade-0.1.0-py3-none-any.whl` from the
current frozen checkout with uv 0.10.6 and installed it with the declared
runtime dependencies into a fresh Python 3.11.14 environment at
`/tmp/beta-final-install/venv`. The active module resolved to
`/tmp/beta-final-install/venv/lib/python3.11/site-packages/chronoclade/`, and
an assertion confirmed that this path is outside the working checkout. The
wheel SHA-256 is recorded in
`validation/beta/install-final/wheel-sha256.txt` (repository-local evidence);
the dependency freeze, command outputs, and installed workflow artifacts are
in the same evidence directory. The final fresh-output workflow logs and
exit-status record are under `validation/beta/install-final/external-fresh/`;
the top-level collision-refusal logs are explicitly separated as negative
controls. The wheel itself is in
`tmp/beta-final-install/wheelhouse/`.

Installed `version`, CLI help, `list-schemes`, and `preflight --stage
prepare|cgmlst` all ran successfully. The software preflights reported the
Python profile workflow available; they do not establish reference or native
tool readiness. `list-schemes` showed registered references without inspecting
public indexes or downloading catalogues. The small offline workflow completed
`prepare → cgmlst --partition-by st → GrapeTree export` using the synthetic
five-profile E. coli fixture. It produced one ST-only block with no unresolved
inputs and five profile and metadata rows in the export. This is an installed
package/CLI smoke test, not biological validation.

The fresh venvs, wheel, and generated evidence occupied about 234 MiB in total
at the time of measurement, below the 2 GiB cap. No native caller, allele
reference/index, large genome data, or SLURM job was installed or exercised.
This evidence does not replace the separate full native/index or real-SLURM
acceptance.

## Independent frozen-source v2 wheel rerun

On 17 October 2026, built a new wheel from the current source matching
`validation/beta/source-freeze-v2.json` (SHA-256
`2064afdafa8b1294450f303f86b2e06e08aa4aa4989f707cb73350977908e013`, 126
Python files). Fingerprints captured before and after the build/test match the
freeze. The wheel SHA-256, install logs, source fingerprints and exact CLI
evidence are under `validation/beta/install-v2/`; the wheel and new environment
are under `tmp/beta-install-v2/` and `/tmp/beta-install-v2/`.

The wheel installed with its declared runtime dependencies into a new Python
3.11.14 uv environment. The imported module resolved to
`/private/tmp/beta-install-v2/venv/lib/python3.11/site-packages/chronoclade/`,
outside the checkout. Its installed-source fingerprint matched the frozen
source hash both before and after the workflow. Version, help, local
`list-schemes`, and `prepare`/`cgmlst` software preflights exited successfully.

From `/tmp`, the installed CLI completed synthetic `prepare → cgmlst
--partition-by st → GrapeTree export`; all three commands exited 0. It analyzed
all five input profiles in one complete ST95 supplied-profile block and
exported five profile rows plus five metadata rows. Preparation readiness was
incomplete for the set because one sample lacks known-scope HierCC evidence;
the profile-only ST analysis still completed. The fixture is an engineering
smoke test, not biological validation. New venv, wheelhouse and evidence used
about 117 MiB, within the 2 GiB cap. No native tools, large references or
assemblies were installed or run.
