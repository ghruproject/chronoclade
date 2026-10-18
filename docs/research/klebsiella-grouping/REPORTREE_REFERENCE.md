# ReporTree reference reproduction — 18 October 2026

The user requested one clear existing-tool example before introducing more
clustering approaches. We ran ReporTree itself on its complete Neisseria worked
example, then used its supplied-partition stability mode on cached ST258 results.
No GitHub push was performed.

## Reference reproduction

- Source: https://github.com/insapathogenomics/ReporTree/wiki/5.-Examples
- Data: 3,791 genomes, 822 loci, 21 countries; all supplied inputs used.
- Exact logged example command, including MSTree, nomenclature and metadata.
- Fresh version2.7.0, commit e5a0e1cb8192b314a97d089f6a18b995eeff94d8.
- Reference outputs version2.0.0, dated2023-03-28. No source patches.
- Adjusted Wallace >=0.99 for >=9 consecutive adjacent comparisons; default
  zero-partition redundancy reduction. This is genetic threshold stability,
  not temporal persistence or context-removal robustness.
- Fresh full run:93.5seconds, after isolated environment setup.

Independent verification passed: all14 stable-region boundaries/lengths,
all472 adjacent point scores (absolute tolerance1e-8), all14 reported partitions
by actual sample memberships, and620 group/country/date summaries. The reference
has473 levels (0–472); the logged range0..473 is end-exclusive. Cluster names and
country display ordering were not confused with membership identity.

At example thresholds39 and78, respectively, the full data contain329 and181
groups including singletons (160 and101 groups containing >=2 genomes). These
are illustrative genetic resolutions, not transferable Klebsiella thresholds.

## Applying the same stability test to ST258

Input: previously calculated exact single-linkage partitions for all7,815 cached
ST258 genomes, retaining the common pair-coverage mask. Every integer threshold
through the final merge was exported. Exported memberships were checked against
every previously saved selected partition in that range. ReporTree performed the
stability scoring and metadata reporting; it did not rerun GrapeTree for ST258.
The existing distance policies are not claimed to be identical to GrapeTree's.

- Alleles: no qualifying range through final merge at20.
- Aligned locus SNPs: qualifying transitions29->28 through37->36, covering levels
  28–37; group sizes at level28 are7,814 and1.

This shows why a plateau alone cannot identify useful fine-scale clusters. Its
absence under the example's settings does not prove biological subgroups absent.
We have reproduced a concrete method, not established an outbreak cutoff.

## Outputs and reproducibility

Start: `tmp/reportree-reference/report/report.html`.
Navigation/commands: `tmp/reportree-reference/README.md`.
Independent audit: `tmp/reportree-reference/verification.json`.
Version/input/command receipt: `tmp/reportree-reference/reproduced/execution.json`.
ST258 commands: `tmp/reportree-reference/st258/execution.json`.
Code: `research/reportree_reference/`.

Browser review confirmed the embedded figures, ten-row table controls and search.
All local evidence links resolve; the report has no external rendering assets.
The current-tool successful run and supplied expected outputs are kept separate.
An initial missing-dependency failure was retained; expected outputs and completion
markers must be checked because this upstream tool can exit0 after a child error.
