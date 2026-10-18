# ChronoClade beta specification

Approved locally by the user on 17 October 2026. This is the working specification
for the beta audit, implementation and acceptance evidence. It supersedes older
roadmap statements where they conflict. Work remains local until the user approves
publication: no pushes, pull requests or releases.

Amendment, 18 October 2026: the user authorised final documentation and two
independent agent walkthroughs, then a GitHub push if the checks pass. A versioned
release is not requested. Human usability feedback moves to beta follow-up; agent
walkthroughs are labelled simulated evidence, not human validation.

## Objective and working rules

Develop ChronoClade into a usable, scientifically validated beta for
bioinformaticians who can run command-line tools but are not specialists in
phylogenetics.

Preserve the existing report style and approved state-change network. Use
parallel subagents for independent tasks where useful, with shared components
rather than duplicated implementations. Essential workflow requirements must
pass without undocumented manual repairs. Experimental research must be labelled
and must not silently determine routine results.

Distinguish correctness, scientific validity, usability and performance evidence.
A valid unsupported result is a successful scientific outcome; a plausible plot
or a completed command does not establish scientific correctness.

## 1. Scope and acceptance plan

Audit the current implementation against this specification. Record:

- Implemented and verified.
- Implemented but insufficiently tested.
- Missing or defective.
- Experimental or dependent on external access.

Prioritise a complete, validated Klebsiella pneumoniae workflow. Identify the
validation status of E. coli/HierCC and other supported inputs.

Before changing scientific methods, specify the proposed method, reference,
benchmark and acceptance criteria. Do not choose thresholds because they make
an example look better.

## 2. Installation and execution

Use BactScout as a reference for reproducible installation, short commands,
preflight checks, configurable warnings and clear output summaries:
https://github.com/ghruproject/bactscout

Provide:

- One recommended installation route, tested on Apple Silicon macOS and Linux.
- Separate software/reference installation and readiness checks.
- Clear handling of credentials, downloads, unavailable references and interrupted setup.
- Progress showing stage, sample counts, elapsed time and failures.
- An end-of-run summary linking to reports and the next command.
- Logs and resumable execution without requiring Python traceback inspection.

Test the walkthrough in a fresh output directory and disclose environment/cache
reuse. Two independent agents should attempt the documented commands and report
reading without developer assistance; record friction. These simulated checks
satisfy the amended beta gate. Recruit representative human users for later
feedback when available; do not describe agent checks as human usability evidence.

## 3. Preparation and native typing

Accept assemblies, accession lists with metadata, and Pathogenwatch collection
identifiers or URLs. Resolve identities, retrieve existing typing/metadata where
available, and run configured native typing when necessary.

Report missing profiles, ambiguous identities, insufficient locus coverage,
conflicting metadata and unavailable databases.

Typing benchmark:

- Freeze Pathogenwatch assemblies and exported cgMLST/cgLIN results.
- Record exact identities, assembly checksums, schemes, tool and reference versions.
- Run those assemblies through the native typing route.
- Compare allele calls locus by locus, including missing calls, and compare lineages.
- Publish concordance tables; investigate every unexplained discrepancy.
- Include incomplete/difficult assemblies, rather than only easy matches.

Exact LIN comparisons require compatible reference databases. Database-version
changes and provisional assignments must remain distinguishable from software
errors. Do not copy a nearest neighbour's LIN code to an unassigned query.

## 4. cgMLST analysis

The report should answer:

- What genomes and metadata entered the analysis?
- Which public genomes are closest to the inputs?
- Where, when and in which hosts are related genomes represented?
- Are inputs concentrated in particular reported lineage groups or distributed across several?
- What representatives were selected for deeper analysis, and why?

Use species/ST for candidate discovery and compatible LIN/HierCC assignments for
analysis blocks. Distinguish discovered, retrieved, excluded and analysed records.
Retrieve context profiles/metadata before downloading context assemblies.

Categorical allele distance is the default. Direct nucleotide comparisons within
the cgMLST genes remain an explicit alternative with accurate units and coverage.
Retain the full comparable pool in NJ and ordination. Identify neighbours using
actual distances and jointly compared loci; retain ties.

Describe local versus wider representation with explicit denominators. Geographic
origins and state-change links are investigation hypotheses, not demonstrated
introductions or transmission routes.

## 5. Additional cluster detection and the ST258 amendment

Keep reported LIN/HierCC groups as the baseline, but do not equate a broad context
search block with a close-relative cluster.

The user explicitly finds the current ST258 grouping too broad. Audit the input
subgroups, available finer LIN levels, actual allele-distance distributions,
missing-locus coverage, and public-neighbour composition. Test whether finer
reported prefixes and/or a documented distance-based view provide a useful
close-relative overview while retaining a wider background for context.

Do not assume that LIN7 is the appropriate display level for every lineage.
Separate search scope, displayed grouping, nearest-neighbour definition and
subsampling scope. Label their denominators and purposes independently.

Review published approaches for distinct genetic groups and relatedness
thresholds, including time. Compare a bounded number of justified candidates
against existing reported groups.

Evaluate using independent epidemiology, known related/unrelated genomes and
held-out studies or time periods. Avoid circular validation against the same
threshold being tested. Report false inclusions, missed relationships,
missing-data/sampling sensitivity, and performance relative to simpler methods.

Do not infer a universal cutoff from a root-to-tip slope. If no candidate improves
reliably, retain the established grouping and document that result. This research
outcome must not indefinitely block the core beta.

## 6. Representative subsampling

Retain inputs, explicit pins and required close neighbours. Select additional
representatives using genetic diversity, time, location and host where available.
Alternative sets vary optional representatives reproducibly.

Report requested/actual counts, mandatory/optional reasons, strata coverage,
unselected-to-selected distances, alternative overlap and cumulative coverage.

Test effects on relevant clades, nearest relationships, clock estimates and
shared-target ancestral dates. Subsampling sensitivity is separate from formal
statistical uncertainty.

## 7. Tree and temporal assessment

Build recombination-adjusted SNP trees for exact selected genomes; do not promise
that every recombination event has been removed. Validate identity, alignment,
callable sites, filtering, distance units and reference choice. Compare important
relationships with published analyses.

Evaluate optional rapid branch support, identifying its source tree and effects
of recombination filtering. Branch support and temporal signal are separate.

Use explicit assessment states:

- Supported: assessment completed and supports proceeding.
- Unsupported: assessment completed but does not support dating.
- Not assessed: insufficient evidence or an unmet prerequisite.
- Failed: computational/validation error prevented assessment.

Review date randomisation, including clustered permutations where sampling
structure warrants them. Record the method and evidence behind each decision.

## 8. Time trees and uncertainty

Date supported trees, or record an explicit user override. Report node-date
intervals, rate uncertainty where available, date precision, clock outliers and
limitations. State the interval level/method accurately; TreeTime's documented
confidence output ordinarily reports 90% marginal node-date regions.

Distinguish conditional node-date uncertainty, between-selection variation, and
topology/root/model uncertainty. Include positive and negative temporal benchmarks.
Insufficient signal must produce a useful unsupported result, not false precision.

## 9. Reports and interoperability

Preserve the established style and approved interactive state network.
Report order: dataset, relationships, temporal diagnostics, state network,
subsamples, methods/downloads.

Provide self-contained HTML, searchable/sortable tables, usable tree navigation,
metadata tooltips and compact legends. Hide empty metadata sections. Root-sensitivity
comparison tables/figures and discussion are omitted from the reader-facing report.

Every output directory needs a README, clear entry points and explained filenames.
Provide exports that actually load into GrapeTree and Microreact. Stable IDs must
match between trees, profiles and metadata; display labels remain traceable.
Test imports in the target software. Matching file extensions are insufficient.

## 10. Benchmarks and beta delivery

Use three testing levels:

| Level | Purpose |
| --- | --- |
| Small frozen dataset | Fast regression, typing concordance, user walkthrough |
| Published biological datasets | Relationships, metadata summaries, temporal conclusions |
| Large profile dataset | Completeness, runtime, memory, disk and report responsiveness |

Choose accessible published inputs and explicit findings. Candidate sources:

- David et al., European hospital spread: https://pmc.ncbi.nlm.nih.gov/articles/PMC7244338/
- Regional/international European spread: https://www.nature.com/articles/s41467-024-49349-z
- Pavia clinical/non-clinical collection: https://www.nature.com/articles/s41564-022-01263-0

Different datasets may validate geography, host and temporal questions. Define
reproducible target findings before comparison; record differing methods/inputs
instead of promising identical results.

Deliver an acceptance matrix linking requirements to commands/artifacts/status;
a reproducible small walkthrough; scientific/resource benchmarks; and a visual
presentation using actual generated reports/results to demonstrate the four
modules, handoffs, validation and remaining limits. Supply release blockers and
explicit research deferrals.

## Method references

- Pasteur LIN definitions and reference dependence:
  https://bigsdb.pasteur.fr/klebsiella/cgmlst-lincodes/
- Murray et al., structured sampling and temporal tests:
  https://doi.org/10.1111/2041-210X.12466
- TreeTime date inference and interval semantics:
  https://treetime.readthedocs.io/en/latest/tutorials/timetree.html
- Microreact shared identifier requirements:
  https://docs.microreact.org/instructions/creating-a-microreact-project/metadata-column-types

## Completion rule

A beta requirement is verified only with a reproducible command, frozen input
identity, inspected result and evidence location. Unit tests do not establish
biological validity; documentation does not establish a working installation;
script generation does not establish a successful real SLURM execution.

## Scope addition: E. coli and public EnteroBase schemes

Requested 17 October 2026: implement and test E. coli support using the public
EnteroBase cgMLST scheme, select a published reproducible benchmark, and reuse
shared reference machinery for other publicly available cgMLST schemes. Keep
catalogue availability, native calling, HierCC assignment, public context and
biological validation as separately evidenced support levels. Preserve the
Klebsiella acceptance work and local-only publication constraint.

See [EnteroBase support and acceptance sequence](ENTEROBASE_SUPPORT.md).
