# Optional human beta feedback

The user moved this activity to beta follow-up on 18 October 2026 because no
participants are currently available. The completed agent checks are documented
in [simulated reviews](SIMULATED_REVIEWS.md).

Use this package with two people who can run command-line bioinformatics tools
but are not phylogenetics specialists. Run the sessions separately. The
participant should get the task sheet and fixture, not coaching. An observer may
record the screen and exact help requests with consent; do not collect
credentials, patient information or private accessions.

This is a workflow and interpretation exercise using invented data. It cannot
validate typing databases, biology, assembly retrieval, SNP trees, time trees or
performance at study scale. Any future sessions should be recorded separately from simulated agent checks.

## Before each session

Prepare a clean checkout and complete the software installation described in
[installation](../installation.md). Check the software needed for this profile-only exercise with
`pixi run -e prepare chronoclade preflight --stage prepare` and
`pixi run -e cgmlst chronoclade preflight --stage cgmlst`. Assembly-tree and dating
tools are not prerequisites for this synthetic route. Do not prepare the participant's
outputs in advance. Give them `examples/beta/profiles.json`,
`examples/beta/catalogues.json`, this task list and a fresh output path. Start
the clock when they begin installation or, if installation is outside the
session, record which steps were already completed.

The fixture contains three invented Klebsiella profiles and three invented
loci: one input and two context records. It includes no assemblies and needs no
network access or credentials. Ask the participant to use fresh output paths.

## Participant task sheet

Please work through these tasks without asking the developer what a command or
result means. You may use the project documentation. Say aloud what you expect
before each command and what the output tells you afterward.

1. **Check the installation.** Find the version and check the required software
   for the prepare and cgMLST stages. Explain what a software preflight does
   and does not verify.
2. **Prepare the fixture.** Produce the prepared dataset. Find its readiness
   result and decide which later steps that result permits.
3. **Analyse the profiles.** Run cgMLST with the documented small-example
   settings. Find the summary, group and nearest-relative evidence in the saved
   report.
4. **Export for GrapeTree.** Produce an export, find how to load it, and explain
   how the profile rows, metadata rows and source IDs are linked.
5. **Recover from an error.** Re-run preparation with a deliberately incorrect
   input path, read the message, then return to the valid fixture and identify
   the next command.
6. **Explain the limits.** Say whether this synthetic result validates
   biological relatedness, native typing, assembly-tree inference or dating,
   and identify what evidence would be required for those steps.

Suggested command route (participants may discover these commands from the
installation and walkthrough docs):

```sh
pixi run chronoclade version
pixi run -e prepare chronoclade preflight --stage prepare
pixi run -e prepare chronoclade prepare examples/beta/profiles.json \
  --input-kind profiles --catalogues examples/beta/catalogues.json \
  --out tmp/user-1/prepared
cat tmp/user-1/prepared/readiness.json
pixi run -e cgmlst chronoclade preflight --stage cgmlst
pixi run -e cgmlst chronoclade cgmlst tmp/user-1/prepared/dataset.json \
  --out tmp/user-1/cgmlst --lin-level 5 --subsamples 1 \
  --context-size 0 --bootstrap-replicates 0
pixi run chronoclade export tmp/user-1/cgmlst/cgmlst.json \
  --viewer grapetree --out tmp/user-1/grapetree
```

For the error task, replace `examples/beta/profiles.json` with
`examples/beta/no-such-file.json`; record the exit code and message. Do not
delete or change the fixture.

## Expected outcomes for the observer

- `version` reports the installed package version. The prepare and cgMLST
  preflights pass in their matching environments; they check software only.
  References, credentials and specific input readiness are assessed separately.
- Preparation completes with three profile-ready rows and one input/two context
  roles. Its output says: “Preparation readiness describes frozen typing and
  metadata; tree assembly availability and temporal signal are assessed in
  later stages.” The fixture's assemblies are unavailable, so tree/time cannot
  run from it.
- cgMLST completes one block. The report and its README expose a group summary,
  nearest evidence and the number of jointly compared loci. These are synthetic
  values, not pathogen epidemiology.
- The GrapeTree bundle has three profile records and three metadata records.
  Each file uses the same export ID; the metadata crosswalk returns to the
  ChronoClade sample ID. Alleles are recoded categories and `0` represents a
  missing call.
- The incorrect path exits nonzero with an actionable file error. The prepared
  output README states the next cgMLST command.
- The participant distinguishes successful profile preparation from database
  readiness and from assembly/tree/time eligibility. A temporal result marked
  `unsupported` means assessment completed without supported dating; it differs
  from a failed computation.

## Session record

Use one copy per participant. Record observations in their words where possible.

| Field | Record |
| --- | --- |
| Participant code and consent | |
| Bioinformatics experience / phylogenetics experience | |
| Date, OS, architecture | |
| Installation steps completed by participant | |
| Start/end time and total duration | |
| Tasks completed independently | |
| Exact command or screen where help was requested | |
| Errors encountered and whether participant recovered | |
| Readiness interpretation (quote or close paraphrase) | |
| Nearest/group interpretation and denominator | |
| Export ID mapping explained correctly? | |
| Unsupported vs failed explained correctly? | |
| Places where participant hesitated or guessed | |
| Most useful documentation or output | |
| Suggested change | |

Observer notes:

```text


```

## Pass rule and follow-up

The usability gate passes only when both independent participants complete all
six tasks without developer help and correctly explain readiness scope, the
synthetic nature of the evidence, nearest-relative denominators, export IDs,
and the distinction between unsupported and failed temporal assessment. Keep
friction notes even when a participant eventually succeeds. If a participant
fails a task or misunderstands an output, retain that as a finding and address
it before repeating the affected task. Do not mark this gate complete from an
agent run, a developer demonstration or the existing synthetic install check.
