# E. coli native reference acceptance

## Scope

The public [EnteroBase Escherichia.cgMLSTv1 catalogue](https://enterobase.warwick.ac.uk/schemes/Escherichia.cgMLSTv1/)
contains 2,513 locus allele files in the inspected snapshot, totalling
286,755,118 compressed bytes. These are allele definitions, not an independently
validated HierCC reference. No fixed 2,512-locus assumption is used.

`setup-typing --species ecoli --prepare-db` uses the shared frozen-index downloader
and builds the pinned Pathogenwatch cgMLST caller index. The previous upstream
public downloader assembled the directory URL without a slash and fetched the
758 MB profile archive just to derive its header. This route instead reads the
small directory listing and downloads only explicitly listed allele files.

Each original gzip file and its source evidence remain unchanged. The adapter
creates a separate `<locus>.fa.gz` view with numeric allele IDs, recording all
non-ACGT/non-numbered exclusions and rejecting duplicate or zero allele IDs.
The view retains the ordered snapshot locus list. The caller version is
`snapshot:<content reference SHA256>`; the download date is not presented as a
biological database release date.

The standard setup manifest separates `cgmlst_ready`, `lineage_ready` and
`ready`. Passing caller readiness must not imply that new assemblies can receive
verified HierCC assignments.

## Current HierCC reference gap

The two `hierBG` files below are not required by hclink. Its fresh build route
downloads public cgMLST profiles and HierCC labels separately, then constructs
a prepared assignment database. Using an existing compatible prepared database
does not require rebuilding it with an EnteroBase API key.

The advertised `hierBG.list.gz` and `hierBG.npy` returned HTTP 403 during bounded
reads. The accessible `profiles.list.gz` header has ST plus 2,513 loci and no
HierCC assignments. The official hclink repository has no public GitHub release
assets. Its documented GitLab image registry rejected a normal anonymous pull
authorization request. These observations are recorded in
`validation/beta/ecoli-reference-probe.json`; they are availability checks, not
proof that a compatible reference cannot be distributed by its maintainers.

The existing authenticated hclink builder remains available with a legitimate
EnteroBase API key. Its compatibility guard currently expects dated databases;
a public content-addressed caller snapshot is deliberately not declared compatible
with a dated hclink reference without a separately verified reference relation.
No query's own known HierCC label is used as its assignment reference.

## Verification and live evidence

`tests/test_enterobase_caller.py` verifies source preservation, numeric allele
normalization, explicit exclusion counts, duplicate/zero-ID rejection,
empty-catalogue rejection, frozen-reference change detection, normalized-file
reuse, the public setup route and independent lineage readiness. Combined caller/setup/query/reference testing is recorded in
`validation/beta/ecoli-native-focused-tests.txt`; the final frozen-source suite
is recorded separately by the parent acceptance run.

The original 4 GiB storage approval was revised after the complete source
manifest showed 9,663,183 alleles / 12,713,186,930 decompressed FASTA bytes.
No decompressed full-catalogue dump is written.

The live harness `validation/beta/ecoli-public-caller-build.py` runs the normal
setup helper against the public source in a separate directory. It reuses the
already installed pinned native tools through explicit absolute paths, limits
new allocated output to 8 GiB, checks at least 25 GiB free disk before starting,
and enforces a 12 GiB sampled process-family RSS cap and records actual wall time and output allocation.
The output is `tmp/typing-setup-ecoli-public/`. The live build completed; results and limitations are recorded below.
Fixture tests alone do not meet native biological acceptance.

## Proposed caller/reference compatibility bridge

The upstream hclink documentation permits pure SHA1 allele profiles even when the
caller catalogue is newer. Its assignment command reads `code` from the input
JSON. Once a genuine HierCC background is available, an explicit checksum mode
can construct that code from the caller's sequence checksums rather than numeric
allele IDs. It must validate the exact locus-ID sets, reorder the calls into the
reference's verified locus order, blank missing and ambiguous calls, and preserve
both reference identities. A changed numeric allele namespace must never be
assumed compatible merely because the organism or locus count matches.

This is a proposed compatibility path, not an implemented or validated HierCC
workflow. It needs held-out reference-profile tests and exact known public query
assignments after a compatible background has been obtained.


## Completed live reference build

The full 2,513-locus catalogue was installed and indexed without truncation.
Source reference identity:
`8365f5eea5b1fc6568aef5115dbcbbdf8406284ec8dc3c15cbe7df178c24ef4f`.
Normalization retained 9,660,329 numbered ACGT alleles and recorded 2,854 excluded
records. Caller readiness passed; HierCC/full-typing readiness remained false.

The actual download/normalization/index/readiness helper run took 1,619.56 seconds
(27.0 minutes) on the local M1 under concurrent application/test load. Peak sampled
family RSS after monitor attachment was 2,410,479,616 bytes (2.24 GiB); peak newly
allocated disk was 4,900,978,688 bytes (4.56 GiB), including SQLite transient files.
The observed output directory also held the one small ST131 query export/assembly
downloaded during the build; this is the actual staged directory footprint.
RSS monitoring began after the download/part of normalization, so this is not an
unqualified full-run peak RSS claim. Worker OS-reported peak RSS is recorded
separately. The original observer was replaced to implement the approved larger
budget without restarting the unchanged worker. No old reference or run was deleted.

Evidence: `tmp/typing-setup-ecoli-public/caller-build-result.json` and
`validation/beta/ecoli-public-caller-build-revised.json`. The guarded run reused
already installed pinned tools through verified explicit paths; it is a live
reference/index acceptance run, not a claim of a clean-machine whole-software
installation benchmark.

## Exact Pathogenwatch assembly comparison

The exact public ST131 assembly for SAMD00110466 / DRR120871 was acquired through
ChronoClade's identity-checked assembly downloader and its server-advertised
cgMLST export. Ordinary `type-queries --profiles-only` completed in 26.11 seconds
(native caller 15.9 seconds). Ordinary full mode refused the missing HierCC
reference with exit status 2; it did not classify a genome without a background.

The export has 2,500 called tokens. Native calls comprise 2,468 known numeric
alleles, 30 novel SHA1 alleles and 15 missing loci across the complete 2,513-locus
universe. Of these, 2,496 tokens agree exactly: 2,467 numeric and 29 SHA1 tokens.
Two called tokens differ because the matching genomic segment boundaries differ:
the frozen public b2306 SHA1 matches the native contig hit with two bases removed
(772 versus 774 bp), while public b1858 matches the native hit with one base added
(756 versus 755 bp). This was verified by exact lowercase, strand-oriented SHA1
matches in the unchanged query assembly. The diagnostic did not alter any calls.

Two public-called loci, b0967 and b3473, remain missing natively; 13 loci are missing
in both outputs. This is **partial concordance**, not a perfect profile match or
proof that the old export and current reference releases can be joined. The
initial known-only comparison is explicitly marked superseded; authoritative
counts, source hashes, boundary diagnostics and execution provenance are in
`validation/beta/ecoli-caller-concordance.json`.

## Published eight-assembly ST95 pilot

Eight independently frozen, ENA-checked assemblies completed the ordinary local
CLI chain: `type-queries --profiles-only` (120.26 seconds), `prepare --query-typing`
(1.36 seconds for the final metadata-corrected dataset), then `cgmlst --partition-by st --backend dense
--bootstrap-replicates 0` (3.16 seconds for the final report). No HC assignment was fabricated and no
public context was fetched. All eight inputs were analysed in one supplied-profile
ST95 block. Actual report/numerical provenance uses the final frozen package source.

Across the 28 independently supplied published allele-distance targets, current
native calls (including novel SHA1 alleles) have Pearson r=0.99923 and Spearman
rho=0.98725. Absolute distances are not identical: 2/28 exact pair counts, median
absolute difference 11 alleles and maximum 21, with 2,496–2,513 jointly called loci.
All 28 pair memberships at the paper's stated 40-allele comparison threshold agree.
These are descriptive comparisons of dependent genome pairs, not independent
validation of a universal threshold or exact replication of the study's method.
Published SNP counts retain separate units. The corrected metadata also gives
different cluster labels to one pair with published/native allele distance zero;
those labels are preserved independently and are not treated as gold-standard HC.

`validation/beta/ecoli-native-panel-pairs.json` records every pair. The generated
standard report has embedded images and inline interactive views; its offline
asset check is in `validation/beta/ecoli-st95-report-artifact-check.json`.
Direct browser inspection was blocked by the browser's file-protocol policy;
no alternate browser route was attempted and no visual acceptance is claimed.
The preliminary pilot's mapping of source niche into `host` was corrected before
final acceptance. The authoritative inputs preserve `source_niche`, `source_type`
and `isolation_source`; only explicitly human samples have host Human. Food and
water sources retain missing host (4 known / 4 missing), without inferring an
animal from a food product. Raw source fields remain in the immutable prepared
source JSON/CSV evidence even when they are not displayed as ordinary sample
columns. The old report remains preliminary; no typed profiles were recomputed.

Authoritative report:
`tmp/typing-setup-ecoli-public/published-panel/final-cgmlst/st_02baf3f94bff/profile_report.html`.
Every one of the 28 comparison pair allele counts and jointly called locus counts
was verified against the normal cgMLST module's `pairwise_distances.csv`; all
profile-matrix entries remained identical after the metadata correction.
Correction provenance is in `validation/beta/ecoli-st95-metadata-correction.json`.
