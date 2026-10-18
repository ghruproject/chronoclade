Current starting points: [Mac walkthrough](../mac-klebsiella-walkthrough.md) and
[recorded Klebsiella example](WORKED_EXAMPLE.md). This page retains earlier fixture
review instructions; simulated acceptance is recorded [separately](SIMULATED_REVIEWS.md).

# Independent user walkthrough

Use this script with two basic bioinformatics users, independently and without
developer coaching. It records usability evidence for
[UX-02](ACCEPTANCE_PLAN.md). No participant findings are
pre-filled; UX-02 remains blocked until two people complete the exercise.

## Pick a route

**Short synthetic route:** profile-only data exercises `prepare`, readiness,
`cgmlst`, its summary/groups/nearest evidence and selection, then GrapeTree
export. It has no assemblies and cannot exercise assembly tree/time. The
three-record/three-locus fixture used for export validation is synthetic; its
target import checks are recorded in
`validation/beta/grapetree-import-validation.md`.

The ready-to-run synthetic fixture is in `examples/beta/profiles.json` and
`examples/beta/catalogues.json`; its README explains its limits. No editing of
Python tests or manual construction of JSON is required. There is no complete
biological full-stack fixture shipped with the repository or release. See
[prepare input requirements](../prepare.md) and [the independent workflow](../independent-suite.md).

From the repository root, run these commands against the bundled fixture. Use
fresh output directories. Its single block and cohort are selected automatically.

```bash
pixi run -e prepare chronoclade prepare examples/beta/profiles.json \
  --input-kind profiles \
  --catalogues examples/beta/catalogues.json \
  --out tmp/user-walkthrough/prepared
cat tmp/user-walkthrough/prepared/readiness.json

pixi run -e cgmlst chronoclade cgmlst \
  tmp/user-walkthrough/prepared/dataset.json \
  --out tmp/user-walkthrough/cgmlst --lin-level 5 \
  --subsamples 1 --context-size 0 --bootstrap-replicates 0

pixi run chronoclade export tmp/user-walkthrough/cgmlst/cgmlst.json \
  --viewer grapetree \
  --out tmp/user-walkthrough/grapetree
```

Open each output `README.md`, then the prepare report and cgMLST `index.html`
and block report. Find the readiness explanation, group summary, nearest
relative counts and shared-locus denominator, ties, and the saved selection
under `selection_bundle/selections/selection-001.json`. Check that the export
profile and metadata IDs match and that `allele_key.csv` and the crosswalk make
the recoding/source IDs understandable. For larger datasets with several blocks or cohorts, use the explicit
`--block`/`--cohort` values listed by the command. Do not run tree/time on
profile-only data.

## Local biological four-module prerelease

For a separate full workflow usability exercise, this checkout has a local-only
packet at `tmp/beta-human-walkthrough/`. Start with its `README.md` and hand
participants `participant-task.md`, the input folder, and a fresh session path;
keep its `assessor-notes.md` separate. This is a temporary working-copy packet,
not a downloadable or published release. The supplied frozen profiles and
lineages plus all eight exact assemblies allow ordinary prepare, cgMLST, tree,
and time steps without setting up Pathogenwatch credentials or Pasteur typing
references. It does not fetch public context. Tree and time still require the
default Pixi native software environment.

The local packet now includes current V7 default-Pixi software preflight evidence
and a separately scoped, independently reviewed V7 Linux ARM64 four-module
run on the same supplied pretyped fixture. These are automated checks, not
participant sessions. It also retains the historical V4 CLI/source identity, a
prior fresh V3 native eight-genome reference run, a separate scoped read-only link audit, and distinct
V4 profile/wheel/source-check evidence. The V3 strict runner remains recorded
as failed at its time-report link verifier. The packet does not claim that the
full native four-module workflow has been rerun under V4. No participant
results are prefilled; UX-02 remains open until two independent people finish
and their observations are recorded.

**Biological full-stack route:** use a separately agreed biological input with
matching complete cgMLST catalogue, compatible lineage assignments, exact
assemblies for every selected sample, dates, and required native tools and
reference databases. Confirm the prerequisites using
[installation](../installation.md), [prepare](../prepare.md) and
[tree-stage](../tree-stage.md). Local assemblies and frozen profiles need no
Pathogenwatch credential. If a source requires a key, configure it locally per
[the credential instructions](../installation.md#pathogenwatch-credentials);
never paste it into logs or notes. Do not treat the synthetic route as
biological readiness.

For that route, use new output folders and exact IDs. For E. coli, supply the
required compatible HierCC level instead of assuming the cgLIN example.
For new assemblies, supply a verified seven-gene MLST ST in the metadata as well
as species. The current native caller produces cgMLST profiles and lineage
assignments; it does not calculate seven-gene MLST ST. Without that ST, public
context discovery cannot define its search pool. Do not substitute a cgMLST
profile identifier or infer ST from a nearest neighbour.

For an accession CSV, the preparation command is:

```bash
pixi run --frozen --environment default chronoclade prepare accessions.csv --input-kind accessions \
  --species "Klebsiella pneumoniae" --catalogues CATALOGUES --out prepared
```

For a new-assembly CSV, use this alternative, with the installed typing config:

```bash
pixi run --frozen --environment default chronoclade prepare assemblies.csv --input-kind assemblies \
  --species "Klebsiella pneumoniae" --typing-config TYPING_CONFIG \
  --catalogues CATALOGUES --out prepared
```

Choose one input route. Check readiness before continuing. For live Klebsiella
assembly typing, the full environment supplies the caller's BLAST dependencies;
the lean `prepare` environment is for existing profiles/provider typing.
For live Klebsiella
context, use the same verified native reference when the saved public and input
typing scopes require normalization:

```bash
pixi run chronoclade preflight
pixi run -e cgmlst chronoclade cgmlst prepared/dataset.json --out cgmlst \
  --fetch-context --typing-config TYPING_CONFIG \
  --lin-level 5 --subsamples 1 --context-size 50
pixi run -e tree chronoclade tree \
  cgmlst/BLOCK/selection_bundle/selections/selection-001.json --out tree-001
pixi run -e time chronoclade time tree-001/tree.json --out time-001
pixi run chronoclade export tree-001/tree.json --viewer microreact \
  --out microreact-export
```

For E. coli with verified compatible HierCC assignments, replace the cgMLST
command with the following, and use `--species "Escherichia coli"` with the
matching typing configuration during preparation:

```bash
pixi run -e cgmlst chronoclade cgmlst prepared/dataset.json --out cgmlst \
  --fetch-context --hiercc-level HC_LEVEL --subsamples 1 --context-size 50
```

Choose `HC_LEVEL` from the supplied assignments. Unknown HierCC reference scope
is not proof of compatibility. A profiles-only E. coli analysis instead uses
`--partition-by st` without `--fetch-context`; it does not add public context.

Replace `CATALOGUES`, `TYPING_CONFIG`, `HC_LEVEL` and `BLOCK` with the actual
reference paths, declared HC level and generated block. Tree uses the saved selection and never substitutes missing
assemblies. Read `assessment.json`: `unsupported` is a completed assessment
without a dated tree; a failed computation is an error. Do not override an
unsupported gate to manufacture a dated result. For the Microreact bundle,
check every Newick tip has matching metadata and a source mapping. Only the
exact tiny synthetic app imports have been checked; larger-study imports remain
untested.

## Record the sessions

Give each participant the route matching the data available. Ask them to think
aloud, but record help requests and exact errors before intervening. Do not
record credentials, patient details or private accessions here.

| Field | Participant 1 | Participant 2 |
| --- | --- | --- |
| Role/background (with consent), date, OS/architecture |  |  |
| Route and starting files |  |  |
| Tasks completed without help |  |  |
| Help requested and response |  |  |
| Failures/errors and recovery |  |  |
| Interpretation questions/misunderstandings |  |  |
| Output paths and IDs checked |  |  |
| Suggested improvement |  |  |

Ask each participant: What does ready mean, and does it guarantee the full
analysis can run? What do the nearest count and shared-locus denominator mean?
Why are selection members included, and what happens if a selected assembly
is missing? How do unsupported and failed tree/time results differ? How did you
verify viewer IDs link back to profile/tree tips and source samples?

**Pass criteria:** both users finish their assigned route unaided and correctly
explain readiness, nearest-relative evidence, selection, relevant tree/time
handoffs and export ID mapping. Record all friction even if they ultimately
succeed. A synthetic pass tests navigation and format understanding only; it
does not validate native typing, assembly acquisition, biological results,
temporal inference or large-study viewer performance.
