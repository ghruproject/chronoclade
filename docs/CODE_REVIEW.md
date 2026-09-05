# Code-quality review

This audit covers the ChronoClade branch prepared for the first public
documentation release. It concentrates on maintainability and failure modes.
Scientific validation of the ST239 result is recorded in the worked example.

## Review method

The review compared the branch with `main`, measured every Python module and
function, traced the external-command path, examined cross-module data shapes,
and ran the native integration test. Generated SVG, JSON and tabular validation
outputs were treated as results rather than source code.

The initial branch added 24,000 lines because it included a real validation
bundle. The useful source-code signal was concentrated in four modules:
`report.py` had 1,356 lines, `evidence.py` 864, `context.py` 820 and
`workflow.py` 707.

## Findings and disposition

| Priority | Finding | Disposition |
| --- | --- | --- |
| Blocker | Report generation mixed plotting, TreeTime parsing, HTML assembly and archive creation in a 1,356-line module. | Fixed |
| High | One 278-line lineage function owned all paths, seven external stages, temporal gating, evidence and serialization. | Fixed |
| High | A bespoke Newick parser sat in the public-health topology path despite Biopython already being required. | Fixed |
| High | Existing stage files were reused without a content fingerprint for their inputs and command options. | Fixed |
| Medium | Matplotlib SVGs and ZIP entry metadata made report bundles change between identical builds. | Fixed |
| Medium | Saved analysis records cross module boundaries as `dict[str, object]`; the JSON contract has no typed model. | Open |
| Medium | Corrected-tree branch support does not enter candidate-group assignment. | Open scientific limitation |

## Structural changes

Temporal assessment, TreeTime output parsing and figure generation now live in
`temporal_report.py`. HTML assembly and supporting-file packaging remain in
`report.py`. Neither module exceeds 750 lines. The main report function fell
from 394 to 185 lines; its remaining length is dominated by the standalone HTML
document rather than control flow.

`workflow.py` now plans lineages, allocates CPUs and coordinates independent
runs. `lineage.py` owns one lineage's files and execution stages. A frozen
`LineageFiles` model gives every stage the same path contract. The core
phylogeny, observed clock, date permutations, dated tree, location model and
report serialization have separate functions. The lineage runner fell from 278
to 59 lines, and its largest helper is 73 lines.

The custom Newick classes were deleted. Biopython now parses the corrected tree,
including standard Newick support values and comments, before ChronoClade
applies its focal/context grouping rule. This removes about 75 lines of parsing
code from a scientifically sensitive path.

Shared helper names imported across the report boundary are public rather than
underscore-prefixed internals. `WorkflowError` moved to a small shared errors
module, avoiding a circular dependency between planning and lineage execution.
The SVG backend now uses a fixed hash seed and omits creation timestamps. ZIP
entries use a fixed timestamp, permissions and order. Identical report inputs
therefore produce byte-identical dated-tree figures and support bundles.

## Stage reuse

Each external stage now writes a fingerprint beside its log. The fingerprint
contains the command and SHA-256 hashes of the direct inputs. ChronoClade reuses
stage outputs only when those values match. Changed assemblies, metadata,
upstream trees or options remove the expected outputs and rerun the affected
stage. The fingerprint is written with an atomic replacement after the command
completes and all expected files exist.

## Typed result boundary

Sixty-seven annotations still use `dict[str, object]`. These dictionaries mirror
the public JSON outputs, so replacing them piecemeal would add casts and adapter
code. A follow-up change should define versioned models for temporal assessment,
context evidence, topology groups and the lineage report. Serialization should
happen at the package boundary. Existing JSON keys need fixture tests before
that migration because users may already parse them.

## Verification

The revised code passes Ruff and all 41 tests. The integration test runs SKA2,
IQ-TREE, ClonalFrameML and TreeTime on a small lineage, reaches the confidence
aware dated-tree branch, performs location reconstruction, and checks the report
bundle. The Material documentation builds with strict link and navigation
validation. The ST239 example preparation script was also exercised against all
72 recorded accessions using a temporary assembly directory.

The review bar is met for this branch. Versioned typed result models remain the
first maintainability task before the public JSON interface is declared stable.
