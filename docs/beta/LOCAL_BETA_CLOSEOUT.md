# Local Klebsiella beta closeout

18 October 2026. This is the final local package check after the per-locus SNP
alignment correction. The checks were completed locally before publication;
the user subsequently authorised a GitHub pull request and merge. No versioned
release or tag is part of this closeout. The package version remains `0.1.0`.

**Decision: local beta closeout passed.** The installed-package run, independent
handoff audit, SNP example, focused tests, lint, strict documentation build and
browser review completed successfully within the scope below.

## Supported beta scope

The supported workflow is exploratory *Klebsiella pneumoniae* analysis on Apple
Silicon macOS: **prepare → cgmlst → tree → time**.

- Categorical allele profiles are the cgMLST default.
- The explicit `--distance-method direct-dna` alternative aligns allele DNA
  independently per locus with MAFFT, then reconstructs the gene-panel alignment.
  It counts joint A/C/G/T substitutions; gaps and unavailable calls remain
  separate evidence. These locus SNPs are not recombination-adjusted whole-genome
  SNPs. The later tree stage provides the recombination-adjusted SNP analysis.
- Reported compatible LIN groups describe lineages and context scope. Closest
  relatives are measured from distances, coverage and ties. Neither a LIN prefix
  nor a distance threshold is an automatically validated outbreak boundary.
- Representatives retain inputs and mandatory neighbours, with reproducible
  optional selections. Their variability is separate from node-date uncertainty.
- Time supports an explicit unsupported result as well as supported dating with
  conditional node-date intervals. A negative temporal result is not a software
  failure, and it is not overridden to force dates.

Automatic discovery of epidemiologically meaningful clusters and time-aware
cutoffs remains research. It does not block this beta. The ReporTree reproduction
and ST258 comparison are [separate scientific work](../research/klebsiella-grouping/REPORTREE_REFERENCE.md).
E. coli/HierCC, Linux x86_64 and physical HPC/SLURM verification remain deferred.

## Current installed-package checks

The wheel was built from this working tree and installed into a new Python 3.11
virtual environment with fresh Python dependencies, without system-site-packages.
Commands ran from outside the checkout, and imports resolved to the installed
package. All 144 packaged files match both the source and installed copy, including
all 138 Python modules and the new alignment components.

Native tools were reused from the locked default Pixi environment. This is not a
new native-tool/reference installation or a new provider download benchmark.
The lean cgMLST Pixi environment also installed from its frozen lock and reports
MAFFT 7.526.

| Check | Result and scope |
| --- | --- |
| Focused alignment/distance checks | 45 passed in the activated default Pixi environment |
| Original bundled allele example | Installed prepare/cgmlst retain all 3 synthetic records |
| Bundled aligned-SNP example | Installed prepare and both distance modes retain 4 records; all 6 expected pairwise counts match |
| Actual MAFFT alignment | Three variable loci, including equal/unequal-length variants; 90-column concatenation, 30 missing N columns, deletion recorded separately |
| Four-module eight-assembly run | Passed in about 3½ minutes; exact IDs/assemblies retained, recombination evidence verified, 100/100 permutations, unsupported time result (p=0.376), no override |
| Reports and package handoff | Browser checked entry page, tree figures, drawn interactive network and SNP report; Ruff, strict docs and diff checks passed |

The original direct invocation of the test Python without Pixi activation could
not find MAFFT. The documented activated environment passed all 45 checks. Installing
the Python wheel alone does not install external bioinformatics tools.

## Examples and navigation

For a fresh clone, begin with the bundled profile and SNP examples in
`examples/beta/README.md`, alongside the [profile walkthrough](USER_WALKTHROUGH.md).
They need no credentials, use explicitly invented data, and do not validate typing
or temporal biology. Their expected pair counts and catalogue hashes are supplied.

For your own biological data, follow the [Mac walkthrough](../mac-klebsiella-walkthrough.md).
For the earlier complete native preparation/context run and separate supported
dating example, see the [recorded four-stage example](WORKED_EXAMPLE.md).

The local entry page is `tmp/beta-local-closeout/index.html`. It links the current
installed-package outputs, the optional SNP example, separate earlier supported
dating evidence and each output README. The eight-genome full-stage fixture uses
three inputs, five frozen comparisons, existing profiles and exact local assemblies.
It does not perform new native typing or a live public-context search.

Package/run/audit receipts live under `validation/beta/local-beta-closeout/`;
large generated evidence remains local. The compact checked-in summary is
[evidence/local-beta-closeout.json](evidence/local-beta-closeout.json).

## Acceptance boundaries

Earlier native typing and public-context checks remain separately recorded in
[the acceptance matrix](ACCEPTANCE_PLAN.md). This closeout exercises the installed
current code without claiming to repeat every previous biological benchmark.
The previous 1,260-test baseline is historical; the 45 focused checks above apply
to this change. The current four-stage run records its own inputs and results.

The approved usability gate comprises independent agent walkthroughs, explicitly
labelled simulated evidence. Real-human usability remains a beta follow-up, not
an established finding or a hidden release blocker. No universal cutoff,
transmission route or geographic origin is claimed.
