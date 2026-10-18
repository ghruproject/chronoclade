# Reproducible ST147 pilot on a laptop and SLURM

**Historical pilot prototype.** This page records a reproducibility design for a
Pathogenwatch context pilot and its preserved fixtures. It is not the current
ChronoClade usage guide. Use the [independent four-module workflow](independent-suite.md)
and [Apple Silicon Klebsiella walkthrough](mac-klebsiella-walkthrough.md) for
current commands. The pilot's context selection, scripts and acceptance checks
are prototype-specific.

Pathogenwatch was the public context provider for this pilot. It froze metadata
and annotations once, acquired a bounded pool and recorded lineage, recombination
and temporal-signal checks. A failed dating gate remains a valid result.

## Inputs and scope

The prototype's intended source was a public Klebsiella pneumoniae ST147
catalogue (`organismId=573`, MLST field `mlst`) with public QC-pass/fail records
and detailed metadata. It used cgLIN assignments for the frozen UUID set with the
`scgMLST629_S` scheme and recorded retrieval time, source bytes, hashes and
provider limitations.

Its focal input schema included assemblies, dates, country/location, species,
ST147, origin and a reference. Public-source stand-ins were to be labelled as
such and excluded from context selection by verified identity. A crosswalk tied
sample IDs to source UUIDs and BioSamples where available; ST alone was not used
to infer cgLIN. Date precision and missing metadata were retained.

## Prototype workflow record

The prototype used the repository's Pixi environment and one Python
interpreter for its stages. Credentials were kept in protected configuration or
`PATHOGENWATCH_API_KEY`, outside command arguments and inputs. Its public-source
stand-in option marked outputs as demonstrations.

The minimum sample count is four, permitting this bounded demonstration rather
than changing temporal significance. The pilot uses p=0.05 and the existing
coherence/divergence gates. Record the resulting gate outcomes explicitly.
Native tools must be installed and supported on the execution platform; metadata,
selection and geography replay can be exercised separately on other laptops.

The prototype driver recorded argument arrays, input hashes, Python/platform,
threads, seed, status and elapsed time. Context provenance included tool
versions, candidate and selected IDs, reasons, a per-record download ledger,
checksums and rule-specific losses. Its acceptance notes identified resumability,
interrupted downloads, exact source-ID-to-FASTA mapping and cache invalidation as
checks; one successful download did not establish bulk-endpoint reliability.

## Cluster design notes

The prototype staged the complete context, referenced assemblies and checkout
on the cluster. Manifest paths had to resolve there, and the same locked
environment was installed. No heavy biological computation ran inside web
requests. Current stage job generation is documented in the
[independent workflow](independent-suite.md#local-and-slurm-jobs).



The prototype batch script requested 4 CPUs, 16 GiB and 12 hours and recorded a
SLURM job log. Account and partition were site-specific. Its design called for
preserving the job ID, resources, software environment and input/output hashes;
script review alone did not establish a completed cluster run.

## Validation and offline reproduction



The historical end-to-end offline fixture froze raw public responses, imports complete and
partial cgLIN, includes duplicate BioSamples, Unknown country, an undated record,
QC failure, singleton and a deliberately failed candidate download. It exercises
the prototype context-preparation implementation with download/SKA transports
substituted, checks manifests and geography, and replays selected IDs/reasons plus
every country table row without live credentials.
These fixtures are plumbing/contract tests, not evidence for biological
accuracy.

The prototype's opt-in upstream contract smoke test checked public search and
one UUID/details join without credentials. It did not freeze all detailed metadata
or verify download/cgLIN capability. Protected export/download checks were
separate from ordinary CI, which did not depend on upstream uptime.

The offline replay checks the annotated snapshot hash, reconstructs selection from
stored SKA distances and settings, and compares exact selected source IDs/reasons
and source table counts. Changing a snapshot, cgLIN scheme/version, selection
settings or focal assembly requires affected preparation/selection stages to be
rerun, with a changed selection fingerprint. Never silently reuse another run's
manifest. Regeneration of country figures does not require phylogeny.

Before release, hand-check a subset of raw-export groups/countries/focal memberships;
reconcile eligible genome/sample-unit counts with duplicate audit, valid groups and
unresolved coverage. Inspect `context/context_geography/index.html`, the integrated
report, exported SVG/PNG and supporting archive at readable size. The full public
catalogue denominator includes undated QC-pass records and is never replaced by
the selected tree tips. Keep separate public, selected and focal denominators.

The figures measure composition of available sequenced records, not prevalence,
incidence, migration direction or transmission. Submission/surveillance bias,
incomplete or provisional cgLIN, organism-specific capability, metadata conflicts
and date uncertainty remain material limitations. Country is not country of
acquisition. A same-ST catalogue is not the full global lineage across all STs.

## Historical public freeze for offline figures

The historical compact public fixture retained 5,647 QC-passing deduplicated
units, source UUIDs, sample-unit/raw-record counts, normalised country/date/QC
fields, cgLIN annotations, three public focal stand-ins and source hashes. It
excluded raw detailed responses and assemblies, so it served as a geography and
coverage freeze rather than a complete metadata source.

The associated offline replay checked the payload hash and regenerated figures,
tables and HTML without network access, assemblies or phylogenetic analyses. The
real-freeze checks reconciled focal country denominators, assignment coverage,
raw/sample-unit counts, source identity and focal memberships.

Metadata and manifest assembly/catalogue paths were written relative to their
files. A regression fixture checked that the complete frozen directory could move
between laptop and cluster without rewriting absolute laptop paths.
