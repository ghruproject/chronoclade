# Simulated first-reader report review

Reviewed 2026-10-18. This is an independent agent reading of delivered report contents, not human validation or usability testing. No source changes, Git mutations, or heavy analyses were performed. I read the HTML text, report-linked evidence where noted, and checked local hyperlink/anchor targets. I did not exercise rendered interactive controls, exports, or mobile/browser layout; those remain unverified by this review.

## Sources

- **Prepare:** `<local-native-run>/prepare/report.html`
- **Context preparation:** `<local-native-run>/cgmlst_context/report.html`
- **cgMLST index:** `<local-native-run>/cgmlst/index.html`
- **cgMLST block:** `<local-native-run>/cgmlst/cglin_3981e3a96c24/profile_report.html`
- **Native tree:** `<local-native-run>/tree/runs/selection-001/jobs/aaf534cc4c9308395256d3b10b6875ddb712abebc0033687209bbdc7ad411a37/report.html`
- **Polished tree:** `tmp/beta-mac-polish-v16/tree/report.html`
- **Unsupported time run:** `<local-native-run>/time/runs/selection-001/jobs/127fadf0f74e28883dc7c20779d1a3f50cfa881a663c4084eabbe32a4c34384a/report.html`
- **Time comparison:** `<local-native-run>/time/comparisons/b2007ab806deee1252185c4a51f80eeb28ac72266715721fbed1652070382021/report.html`
- **Positive time run:** `tmp/beta-positive-time-closeout-v16/selection-001/jobs/99f3fdae583c04ed42c0f64e626411715f3fdcd62c417bed37f6c4f5617434c8/report.html`

## Answers obtainable as a reader

### What was input?

The prepare report identifies two input assemblies, **ERR8771929 and ERR8772239**, both ready, with 628/629 called loci. Local assemblies were supplied; provider retrieval was not requested. The report says preparation discovers no public context. The selected-sample tree table makes their metadata explicit: both Greece, 2019-05-14, host unknown. The cgMLST block identifies the organism as *Klebsiella pneumoniae*, ST258, analysis at cgLIN5 prefix `0.0.105.6.0`.

### What public context was added?

The cgMLST index and block explain the scope and funnel: **13,421 same-ST candidate source records → 386 prefix-matching source records → 232 unique public comparisons**. Context is the union of input LIN8 prefixes; analysis is at LIN5 and colours are at LIN8. All 232 public profiles are usable and comparable in this block, producing 234 genomes with the two inputs. The report explicitly warns that prefix membership does not establish a tight pairwise cluster. Context assemblies are acquired only by the separate tree stage. Context preparation reports 230 ready and four incomplete frozen typing records, whereas the cgMLST block uses 234 usable profiles; these are evidently different readiness concepts, but their relationship is not explained directly to the reader.

### Closest relatives and distance denominator

The cgMLST block gives a useful at-a-glance table and a full tied-neighbour table:

- **ERR8771929:** 17 equally closest public profiles, **0 allele differences over 626–628 jointly called loci**. Countries: Germany, Greece, Netherlands. The complete accession list is ERR5432174, ERR5432269, ERR5432346, ERR5432342, ERR5432273, ERR5432275, ERR5432271, ERR3712846, ERR5432270, ERR10669789, ERR8771775, ERR5432268, ERR5432344, ERR11818296, SRR16388401, ERR8772882, ERR8771776.
- **ERR8772239:** nine equally closest public profiles, **1 allele difference over 627–628 jointly called loci**. Countries: Greece, Slovakia, United Kingdom, Unknown. Accessions: ERR1204918, ERR1334451, ERR1289745, ERR1289762, SRR10991758, ERR1334607, ERR2523242, ERR2523241, ERR11818327.

The search universe is the 232 compatible public profiles in this block, rather than the entire public catalogue. Each jointly called gene contributes one mismatch if allele identifiers differ; missing calls are excluded. The report states raw differences determine ranking and retains ties. The normalized values and shared-locus denominator are provided alongside the ranking.

The polished tree clearly presents a different comparison: nearest **selected assembly context** using corrected SNPs, restricted to the exact 28-genome set. ERR8771929 is closest to ERR8771775 at **6 SNPs / 5,171,197 callable sites**; ERR8772239 to SRR10991758 at **37 SNPs / 5,171,197 sites**. This makes the allele-versus-SNP distinction and denominator clear without implying a search over unselected assemblies.

### Geography

The full profile-pool country table separates inputs and context. Combined totals are Greece 120 (two inputs + 118 context), Bulgaria 28, United Kingdom 19, Germany 15, Unknown 13, Switzerland 9, United States 8, Australia 4, Montenegro 4, Netherlands 4, Serbia 4, Norway 2, and Canada, Slovakia, Slovenia, Sweden one each. Host is Human 36 and Unknown 198.

The selected assembly reports separately show **Greece 12, Germany 11, Netherlands 1, Slovakia 1, United Kingdom 1, Unknown 2**; Human 7 and Unknown 21. These sum to 28 and the cgMLST selection comparison makes the change from the full pool inspectable. Networks explicitly describe inferred country-state changes on a tree and warn that arrows do not demonstrate transmission. The absent host reconstruction is explained by there being only one distinct known host state.

### Chosen subsample and next action

There is **one saved selection, selection-001: all two inputs + 26 public comparisons = 28 genomes**. All 26 context genomes are marked mandatory nearest neighbours, consistent with retaining all 17 + nine tied nearest public profiles. The report explains boundary ties and provides a manifest download. It explicitly says the cgMLST step has not run the downstream tree/time stages.

I can identify the next stage as the **tree module**, and the linked block README identifies `selection_bundle/ensemble.json` as its input. However, **no copyable next tree/time command appears in the cgMLST report or its linked README**. A new CLI user must consult separate documentation to know the actual invocation. I did not invent a command from repository implementation.

### Is dating supported, and what uncertainty is given?

For the real 28-genome beta selection, **no**. Both the polished tree and unsupported time report state that observed root-to-tip fit did not exceed randomized-date expectations: **p = 0.0891, threshold 0.05**. The time report states dated tree unavailable. The comparison report shows **0/1 requested runs produced dates**, 25 dated metadata records and three undated, and unavailable rate/MRCA/interval estimates. Thus the reports do not claim calendar dates or transmission times for this selection. With only one requested selection there is also no empirical between-selection spread to interpret.

The separate positive example has **63 genomes (12 inputs + 51 context)**, years shown on the time tree, and a supported temporal-screen result. Its saved assessment reports observed **R² 0.39**, screening rate **1.265e-6**, **p = 0.00990099**, threshold 0.05, 100 successful permutations, and 53 clusters. It explicitly limits this to a fixed-topology screen: passing does not establish model adequacy or eliminate sampling bias.

The positive report labels a separate figure **TreeTime node-date uncertainty; 90% max-posterior regions**, and provides `node_dates.csv`, `timetree_confidence.csv`, and `molecular_clock.txt`. Reading these linked downloads gives final TreeTime rate **1.02e-6 substitutions/site/year ± 9.8e-8 (one standard deviation)**, whole-tree root point estimate **1995.6604 (1995-08-30)**, and root 90% region **1994.268418–1998.078415**. Median internal-node interval width is about 1.72 years; maximum about 3.81. These are within-run uncertainty summaries; they are not evidence of robustness across alternative selections. The final TreeTime rate differs from the screening rate, but their roles are not explained together in reader-facing prose.

## Findings

1. **Next-command gap, moderate:** the cgMLST subsample section is sufficient to identify and download the chosen manifest but insufficient to execute the next module without consulting separate documentation. Give the reader a concrete tree command using the saved ensemble and subsequently the time command.
2. **Unsupported-time navigation gap, minor:** the unsupported time report has only its two download links and no upstream SNP-tree/cgMLST link. The native tree also lacks upstream navigation; the polished tree adds a working cgMLST link. The positive time report adds a working SNP-tree link. Apply this navigation consistently, including failure reports.
3. **Dating interpretation needs a compact numerical summary, moderate:** the positive report's headline is clear and uncertainty figure is explicitly labelled, but the exact p-value and screen limitation live inside a lengthy saved JSON assessment; the final clock rate/root interval require downloads. A short summary could distinguish screen evidence from final clock fit and identify the root estimate as the ancestor of this selected set. The existing downloads preserve the evidence.
4. **Readiness terminology, minor:** context preparation says four incomplete typing records while cgMLST says all 234 profiles are usable. This need not be inconsistent, but a first reader may infer four exclusions. Explain that missing resolved lineage readiness differs from profile comparability if that is the intended distinction.
5. **Suspicious historical date visible, minor/data quality:** full-pool collection-year summaries include one record in **1800**, and one LIN8 group spans 1800–2025. This remains visible rather than silently corrected, which is good, but merits an explicit suspicious-date flag or review cue. The selected 28 have a 2012–2022 metadata interval, so the comparison audit does not suggest this record entered that tree.
6. **Precision/readability, minor:** call overlap displays long floating-point fractions, and machine keys such as `not_supported` appear in the temporal section even though the header has plain language. Rounded overlap percentages and plain labels would be easier to read.

## Link check and limits

No missing local file targets or missing same-page anchor IDs were found in the prepare report, cgMLST index/block/file-browser pages, native tree, polished tree, unsupported time run/comparison, or positive time report. Embedded data-image links were recognized but not browser-exercised. This checks target existence, not that every download renders correctly or an interactive control behaves correctly. The positive uncertainty CSV and molecular-clock download were read directly. No prior validation notes were used to establish the answers above.

Overall the report contents support the principal biological reading: two inputs gain a documented restricted public background; all tied nearest profile relatives are visible with callable-locus denominators; the assembly comparison explicitly narrows to the selected genomes; the negative run declines dating; and the positive example preserves labelled within-run uncertainty. The main first-user task gap is moving from a saved cgMLST selection to an executable next command.
