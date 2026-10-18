# EnteroBase species expansion

This expands the beta scope requested on 17 October 2026. It does not replace the
Klebsiella acceptance work. A downloadable scheme, a working native caller,
lineage assignment, public-context retrieval and a validated report are separate
acceptance gates. Do not advertise full species support when only a catalogue
has been downloaded.

## Public catalogue check

Live index observations on 17 October 2026, with exact index hashes and listed
compressed byte totals in `validation/beta/enterobase/scheme-inventory.json`:

| EnteroBase directory | Public index | Loci listed | Compressed allele files |
| --- | --- | ---: | ---: |
| Escherichia.cgMLSTv1 | Available | 2,513 | 286,755,118 bytes |
| Salmonella.cgMLSTv2 | Available | 3,002 | 306,270,843 bytes |
| Streptococcus.cgMLSTv1 | Available | 372 | 118,452,752 bytes |
| VIBwgMLST.cgMLSTv1 | Available | 1,128 | 265,136,733 bytes |
| Yersinia.cgMLSTv1 | Available | 1,553 | 48,235,693 bytes |
| clostridium.cgMLSTv1 | Available | 2,556 | 20,142,252 bytes |
| Photorhabdus.cgMLSTv1 | HTTP 403 | Not verified | Not downloaded |
| klebsiella.pasteurcgmlst | HTTP 403 | Not verified | Not downloaded |

These are directory identities, not assertions about every species within a
genus. The existing authenticated Pasteur Klebsiella reference route remains
separate. Derive loci from the downloaded scheme; do not hardcode an expected
published locus count or equate versions of different schemes.

Source: [EnteroBase scheme directory](https://enterobase.warwick.ac.uk/schemes/).

## User commands

```bash
chronoclade list-schemes
chronoclade list-schemes --inspect
chronoclade setup-scheme Escherichia.cgMLSTv1 --destination refs/ecoli
```

Listing is offline unless `--inspect` is supplied. Inspection returns exit 2 if
any registered index is inaccessible; available entries remain displayed.
Installation fetches only the chosen scheme's allele FASTA files, not the large
profile table. It freezes the index, validates gzip/FASTA/DNA, records file and
reference checksums, and resumes the same snapshot only. Changed cached data
fails validation. Choose a new directory to obtain a newer snapshot.

`reference.json` means the allele catalogue is ready. It is not a native caller
index or a HierCC assignment database. Native caller setup is integrated through
`setup-typing`; its readiness report must preserve this distinction.

## E. coli first

Use the public Escherichia.cgMLSTv1 catalogue with the pinned Pathogenwatch native
caller, then verify exact sample identities and profile concordance. Its public
profile table has 2,513 locus columns plus ST, without HierCC columns. Those public profiles support direct distance analysis. Official HierCC labels
are separate: upstream hclink downloads them independently when building its
reference. The optional `hierBG.list.gz` and `hierBG.npy` files returned HTTP 403
in our probe; hclink does not require those particular files. A compatible
prepared hclink database would avoid a fresh API build. Do not manufacture
HierCC codes or silently use cgLIN labels as substitutes. Retrieval of existing
Pathogenwatch assignments is a different route from assigning a new assembly.

Benchmark candidate: [Parameters for one health genomic surveillance of
Escherichia coli from Australia](https://www.nature.com/articles/s41467-024-55103-2).
The preregistration must account for its 2025 correction and freeze source
metadata before choosing a manageable panel. It provides host/location and
cgMLST-versus-SNP comparisons and a useful negative temporal-signal expectation.
See `ECOLI_BENCHMARK_PLAN.md` for the frozen protocol when completed.

## Acceptance sequence

1. Validate catalogue installation and interrupted-download recovery.
2. Build and exercise the native E. coli caller on independently known profiles.
3. Resolve and verify the HierCC reference dependency, retaining source identity.
4. Run the normal prepare and cgmlst workflow on a preregistered published panel.
5. Compare tree and temporal conclusions with that study, recording differences
   in samples, reference snapshots, loci and methods.
6. Add further listed species through the same shared infrastructure, each with
   its own caller, group assignment, context and biological acceptance evidence.

The shared catalogue machinery is reusable; tests for one organism do not prove
biological acceptance of the others. No GitHub publication is authorised yet.

### Profile-only analysis while lineage assignment is unavailable

The explicit `cgmlst --partition-by st` route analyses supplied inputs by
species/ST/compatible profile scope. It refuses public context and automatic
context retrieval and does not create a lineage code. Default lineage partitioning
remains unchanged. The new mode has focused membership, exclusion, context-refusal
and interrupted-resume tests. Native profile-only calling through `type-queries --profiles-only` has completed
the ordinary CLI pipeline on eight published ST95 assemblies. Exact counts and
limitations are in `ECOLI_NATIVE_ACCEPTANCE.md`; this is not HierCC acceptance.
