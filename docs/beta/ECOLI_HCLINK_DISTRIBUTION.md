# E. coli HierCC reference distribution audit

Checked 17 October 2026, read-only. No credentials were requested; the existing
configured Pathogenwatch API key was used only for one public-genome HClink
export, and its value was never printed or saved. No full profile/allele/HierCC
resource or OCI image was downloaded, and no package code was changed.
Machine-readable statuses and exact paths are in
`validation/beta/ecoli-hclink-distribution.json` (repository-local evidence).

## Finding

No anonymous static hclink E. coli database download was identified among the
documented routes checked. The upstream [`hclink` README](https://github.com/pathogenwatch-oss/hclink)
and [`schemes.json`](https://github.com/pathogenwatch-oss/hclink/blob/main/schemes.json)
show two distinct inputs: public EnteroBase cgMLST profiles/allele FASTAs, and a
separate HierCC profile endpoint. Upstream `hclink build` fetches HierCC batches
using HTTP Basic authentication with an EnteroBase API key, then downloads the
public profile and allele files. `write_db` combines those sources locally into
the search database. The public `profiles.list.gz` header checked here has an ST
column and 2,513 loci, with no HierCC fields; only the first compressed chunk
was read.

The only upstream prebuilt artifact route found is a versioned GitLab container:
the [tag workflow](https://github.com/pathogenwatch-oss/hclink/blob/main/.github/workflows/build-tagged-image.yml)
logs in with a GitLab token secret and pushes the versioned image;
[`Dockerfile.pw`](https://github.com/pathogenwatch-oss/hclink/blob/main/Dockerfile.pw)
consumes that image. The anonymous registry pull/token probes returned 401 and
403. GitHub releases and tags returned empty lists; PyPI returned 404. The
workflow's current `docker build` line provides neither `--build-arg SPECIES`
nor the EnteroBase `--secret id=api` mounted by the current `Dockerfile`, so
the workflow source does not by itself demonstrate that a fresh tagged build
succeeds. No credential or bypass was attempted.

Pathogenwatch documentation describes browsing public genomes and cgMLST
context search. Its clustering guide says its profiles are functionally similar
to EnteroBase HierCC profiles, with some divergence; it does not specify a
standalone local hclink database package. The recent local provider response
advertises an E. coli `hclink-562-4` (HierCC) job alongside `cgmlst-562-1`, but
only for one public genome. I requested that exact ID through the documented
`POST /api/downloads/{download}` endpoint using the existing configured API
credential. The returned single-row HClink CSV identifies `SAMD00110466`
(`iLGx7naVnHXTvwPF8DEgAo`), closest ST `246645`, distance `35.7`, empty values
at HC0, HC2, HC5, HC10 and HC20, then HC50 through HC2000 all `7`, and HC2350
`1`. These are the service's returned values; none were inferred locally.
The result has no hclink database release/date, profile scope version, or
reference fingerprint. Since its closest ST distance is nonzero, retain it as
that service job's nearest-profile assignment, not as a verified direct
EnteroBase label for this genome. The paired cgMLST export covers 2,500
observed loci and has no reference fingerprint. The raw compressed response,
decoded assignment and request provenance are stored under
`validation/beta/ecoli-hclink-distribution/`.

Ordinary `prepare` now uses the advertised Pathogenwatch route and has imported
the server-produced assignment for the exact public genome
`DRR120871`/`SAMD00110466`; the request, response checksum and imported fields
are recorded in `validation/beta/ecoli-hclink-final-live/verified-run/acceptance-evidence.json`.
This is one provider-labelled record. The route does not expose the indexed
HClink database, its release fingerprint, or the full EnteroBase profile-to-HC
map, so it cannot support local assignment of arbitrary new genomes. Native
E. coli profile-only results remain separate; do not infer HierCC labels from
ST or profile similarity.

## Primary sources

- [HClink README](https://github.com/pathogenwatch-oss/hclink), [database build source](https://github.com/pathogenwatch-oss/hclink/blob/main/src/hclink/build.py), [Dockerfile](https://github.com/pathogenwatch-oss/hclink/blob/main/Dockerfile), and [tagged-image workflow](https://github.com/pathogenwatch-oss/hclink/blob/main/.github/workflows/build-tagged-image.yml).
- [Pathogenwatch cgMLST clustering documentation](https://cgps.gitbook.io/pathogenwatch/technical-descriptions-of-analysis-tools/trees-clustering-and-context-search/cgmlst-clustering).
