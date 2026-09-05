---
name: ChronoClade Lineage Reports
description: A public-health case-conference evidence wall for defensible genomic review.
colors:
  ink: "#17191f"
  muted-ink: "#5f6470"
  hairline: "#c9ccd4"
  paper: "#ffffff"
  outer-wash: "#eef0f4"
  evidence-cobalt: "#2855a6"
  pass-green: "#14734f"
  decision-vermillion: "#d63c2f"
  interval-vermillion: "#e2aaa3"
  panel-wash: "#f4f5f7"
  evidence-hover: "#edf3ff"
  caution-wash: "#fff8dc"
  focus-gold: "#ffb800"
  selection-yellow: "#ffdf78"
  scroll-track: "#e5e7ec"
  footer-rule: "#555b67"
typography:
  display:
    fontFamily: 'Archivo, "Helvetica Neue", sans-serif'
    fontSize: "clamp(32px, 5vw, 68px)"
    fontWeight: 700
    lineHeight: 1.12
    letterSpacing: "-0.035em"
  headline:
    fontFamily: 'Archivo, "Helvetica Neue", sans-serif'
    fontSize: "clamp(27px, 3.2vw, 44px)"
    fontWeight: 700
    lineHeight: 1.05
    letterSpacing: "-0.025em"
  title:
    fontFamily: 'Archivo, "Helvetica Neue", sans-serif'
    fontSize: "21px"
    fontWeight: 700
  body:
    fontFamily: 'Archivo, "Helvetica Neue", sans-serif'
    fontSize: "16px"
    fontWeight: 400
    lineHeight: 1.58
  label:
    fontFamily: 'Archivo, "Helvetica Neue", sans-serif'
    fontSize: "11px"
    fontWeight: 700
    letterSpacing: "0.08em"
  decision:
    fontFamily: 'Archivo, "Helvetica Neue", sans-serif'
    fontSize: "12px"
    fontWeight: 700
    letterSpacing: "0.08em"
spacing:
  micro: "4px"
  compact: "8px"
  control-y: "10px"
  control-x: "12px"
  cluster: "18px"
  section: "22px"
  gutter: "24px"
  panel: "30px"
  panel-wide: "42px"
  stage-wide: "44px"
components:
  report-identity:
    backgroundColor: "{colors.paper}"
    textColor: "{colors.ink}"
    width: "100%"
  stage-navigation:
    backgroundColor: "{colors.paper}"
    textColor: "{colors.ink}"
    padding: "13px 18px"
  stage-rail:
    backgroundColor: "{colors.panel-wash}"
    textColor: "{colors.evidence-cobalt}"
    padding: "38px 24px"
    width: "248px"
  decision-review:
    backgroundColor: "transparent"
    textColor: "{colors.evidence-cobalt}"
    typography: "{typography.decision}"
    padding: "10px 12px"
  decision-proceed:
    backgroundColor: "transparent"
    textColor: "{colors.pass-green}"
    typography: "{typography.decision}"
    padding: "10px 12px"
  evidence-manifest:
    backgroundColor: "{colors.paper}"
    textColor: "{colors.ink}"
    padding: "14px 0"
  measure-strip:
    backgroundColor: "{colors.paper}"
    textColor: "{colors.ink}"
    padding: "16px"
  confidence-grid:
    backgroundColor: "{colors.paper}"
    textColor: "{colors.ink}"
    padding: "16px"
  figure-sheet:
    backgroundColor: "{colors.paper}"
    textColor: "{colors.muted-ink}"
    width: "100%"
  evidence-table:
    backgroundColor: "{colors.paper}"
    textColor: "{colors.ink}"
    width: "100%"
  report-completion:
    backgroundColor: "{colors.ink}"
    textColor: "{colors.paper}"
    padding: "30px 42px"
    width: "100%"
---

# Design System: ChronoClade Lineage Reports

## Overview

**Creative North Star: "The Public-Health Case Conference Evidence Wall"**

ChronoClade reports feel like an expert case conference assembled on a white working wall: numbered sheets establish the order of review, large scientific figures carry the evidence, and hairline rules keep every claim attached to its source. The world is rigorous, direct, and deliberately more editorial than dashboard-like.

Near-black ink supplies authority, cobalt marks evidence and continued investigation, green marks a passed gate, and vermillion marks a failed gate or stop decision. Pale interval vermillion carries dated-tree uncertainty without competing with the stronger status marks. The report stays dense enough for expert scrutiny while preserving a clear reading sequence, explicit gates, and ready access to underlying files.

**Key Characteristics:**

- White paper on a pale outer wash, divided by near-black and grey rules.
- Embedded Archivo throughout, with italic spaced taxon names, oversized compact headings, and tabular numerals for measures.
- A static numbered four-stage rail that keeps the analytical sequence visible without tracking the reader.
- Large figures—including dated phylogenies with visible 90% internal-node intervals—paired with open evidence manifests and downloadable source files.
- Flat, square, print-capable surfaces without decorative chrome or shadow.

## Colors

The palette separates evidence, decisions, supporting structure, and caution without asking colour to carry meaning alone.

### Primary

- **Evidence Cobalt** (`{colors.evidence-cobalt}`): Marks exploratory stages, evidence links, charts, scrollbars, and review states.

### Secondary

- **Pass Green** (`{colors.pass-green}`): Marks passed gates and explicit proceed decisions.
- **Decision Vermillion** (`{colors.decision-vermillion}`): Marks failed gates, stop decisions, and observed values in statistical figures.
- **Interval Vermillion** (`{colors.interval-vermillion}`): Shows 90% inferred node-date intervals behind the stronger vermillion root and decision marks.

### Tertiary

- **Focus Gold** (`{colors.focus-gold}`): Supplies the visible keyboard-focus outline.
- **Selection Yellow** (`{colors.selection-yellow}`): Makes selected report text unmistakable.
- **Caution Wash** (`{colors.caution-wash}`): Backs interpretation boundaries and guardrails without competing with decision marks.

### Neutral

- **Case Ink** (`{colors.ink}`): Primary copy, strong rules, identity blocks, and the report-completion footer.
- **Muted Ink** (`{colors.muted-ink}`): Questions, captions, units, metadata, and table labels.
- **Hairline Grey** (`{colors.hairline}`): Internal dividers within measures, downloads, navigation, and tables.
- **Paper White** (`{colors.paper}`): The report surface, figures, cards, and sticky table cells.
- **Outer Wash** (`{colors.outer-wash}`): Separates the bounded report sheet from the browser viewport.
- **Panel Wash** (`{colors.panel-wash}`): Stage rail sheets and paired logic annotations.
- **Evidence Hover** (`{colors.evidence-hover}`): Gives evidence navigation and the mobile table notice a quiet cobalt tint.
- **Scroll Track** (`{colors.scroll-track}`): Keeps horizontal evidence-table scrolling visible but recessive.
- **Footer Rule** (`{colors.footer-rule}`): Separates the label column inside the dark report-completion footer.

### Named Rules

**The Evidence/Decision Rule.** Use green for pass/proceed, vermillion for fail/stop, and cobalt for review or unresolved evidence. Observed statistical marks may remain vermillion when they are not status indicators.

**The Redundancy Rule.** Every colour-coded state also carries a readable label such as “PROCEED”, “CONTINUE TO TEST”, or “DO NOT TIME-SCALE”.

**The Interval Is Uncertainty Rule.** Use pale interval vermillion only for dated-node ranges; keep roots, observed marks, and gate decisions in the stronger decision vermillion.

## Typography

**Display Font:** Archivo (with Helvetica Neue and sans-serif fallbacks)  
**Body Font:** Archivo (with Helvetica Neue and sans-serif fallbacks)  
**Label Font:** Archivo (with Helvetica Neue and sans-serif fallbacks)

**Character:** A single embedded variable sans serif makes the report self-contained and gives scientific content a direct, institutional voice. Hierarchy comes from scale, weight, tracking, and case rather than from ornamental type pairing.

### Hierarchy

- **Display** (700, fluid 32–68px, 1.12): Species and lineage identity; set the reader-facing taxon name in italics with source underscores converted to spaces, then place the lineage on its own line.
- **Headline** (700, fluid 27–44px, 1.05): The scientific question that opens each stage.
- **Title** (700, 21px): Section headings and prominent measure values.
- **Body** (400, 16px, 1.58): Explanations and findings, held to 75 characters per line; identity metadata is held to 70 characters.
- **Question** (400, 18px): A muted explanatory line immediately beneath each stage question.
- **Label** (700, 11px, 0.08em tracking, uppercase): Measure names, file types, and decision metadata.
- **Rail number** (800, 74px, 0.8 line-height, -0.06em tracking): The stage number; it compresses to 48px on narrow screens.

### Named Rules

**The Question-First Rule.** Every analytical stage begins with a plain-language scientific question before statistics or files.

**The Embedded-Type Rule.** Keep Archivo embedded in every generated report so offline, shared, and printed copies retain the same hierarchy.

**The Taxon Display Rule.** Present scientific taxon names as readable spaced names in italics; preserve machine identifiers only in source data and filenames.

## Layout

The report is a single white sheet capped at 1480px and centred on the outer wash. Its primary grid pairs a 248px stage rail with a flexible evidence field. Each stage is a full-width horizontal sheet divided from the next by a strong two-pixel ink rule; the body uses generous 40px 44px 50px padding.

The first stage gives the scientific figure the majority of the evidence field: a narrow question-and-measure column sits beside the figure, with a compact evidence manifest at its right. Later stages use a flexible figure-plus-270px manifest pattern. Standard measure strips form three equal ruled cells; time-tree confidence uses six cells in a 3×2 ruled grid. Comparison logic uses two equal sheets joined by a directional arrow.

At 1100px, the first-stage split and all figure-manifest pairs stack into the semantic sequence figure → meaning or verdict → downloads. At 800px, the rail narrows to 76px, stage labels turn vertically, stage padding reduces to 30px 18px 38px, and the four sticky stage links become equal compact columns without horizontal navigation scrolling. Standard measures and all six confidence cells stack, and download rows become two-column. Wide tables keep a 680px minimum width inside labelled, keyboard-focusable scroll regions; a visible instruction announces horizontal scrolling and the first column remains sticky. Compact two-column tables instead become bordered evidence rows, with the header hidden and the first cell acting as the row label.

Print removes the sticky contents rail, returns the canvas and footer to white, keeps stages together where possible, and exposes detail content so the inferential order and evidence remain legible on paper.

### Named Rules

**The Ordered Wall Rule.** Preserve the Explore → Test → Time scale → Interpret sequence and its numbered rail in every lineage report.

**The Evidence-Beside-Claim Rule.** Pair each major figure with its open evidence manifest at wide widths; on narrow screens preserve the semantic figure → meaning or verdict → downloads order.

## Elevation & Depth

The system is entirely flat. It uses no shadows, gradients, blur, or simulated lift. Depth and grouping come from the contrast between the outer wash and paper, one- and two-pixel rules, and sticky positioning for the stage navigation and mobile table key.

### Named Rules

**The Flat Evidence Rule.** Never introduce box shadows or floating cards; a report surface earns hierarchy through rules, contrast, and position.

## Shapes

Sheets, verdicts, figures, manifests, measures, and tables are square. Borders are one-pixel hairlines for internal structure and two-pixel ink rules for major stage or decision boundaries. The system does not use pills, rounded cards, decorative clipping, or soft container silhouettes.

### Named Rules

**The Square Sheet Rule.** Keep report containers and state marks at square corners; rounded UI would weaken the case-file character.

## Components

### Report Identity

- **Character:** A compact dark institutional mark anchors an oversized species and lineage heading.
- **Structure:** The 248px identity mark aligns with the stage rail; the flexible copy field holds an italic, reader-facing taxon name with spaces instead of underscores, the lineage on its own line, report metadata, and the overall temporal decision.
- **Responsive behavior:** The grid inherits the 76px narrow rail while identity copy reduces to 22px 20px padding.

### Stage Rail

- **Default:** A pale square sheet with a 74px cobalt number and an uppercase 12px label.
- **Reading state:** The rail stays visually stable while scrolling; it does not highlight or track the section in view.
- **Responsive behavior:** At 800px the rail becomes 76px wide, the number becomes 48px, and the label runs vertically.

### Stage Navigation

- **Style:** A sticky horizontal rail of equal-width text links, offset by the stage rail on wide screens.
- **States:** Evidence-hover tint on hover; strong gold outline with a three-pixel offset on keyboard focus.
- **Responsive behavior:** It spans the full viewport; at 800px the four links divide the available width equally, reduce to 11px type with compact horizontal padding, and remain fully visible without horizontal scrolling.

### Decision Marks

- **Shape:** Square transparent label with a two-pixel current-colour border and 10px 12px padding.
- **Review:** Cobalt, used for exploratory continuation.
- **Proceed / Stop:** Vermillion, used when the gate changes the analytical path.
- **Accessibility:** Wording states the decision independently of colour.

### Figure + Evidence Manifest

- **Figure:** Full-width scientific image on white with a one-pixel ink frame and a muted 13px caption.
- **Manifest:** Native open disclosure bounded by ink rules; each file row gives title, description, and uppercase format across hairline dividers.
- **Responsive behavior:** The side-by-side grid stacks at 1100px in semantic figure → meaning or verdict → downloads order.

### Dated Phylogeny

- **Structure:** Time-scaled branches are near-black, sampled genomes are cobalt points, and the inferred root is a strong vermillion point with a labelled interval.
- **Uncertainty:** Every inferred internal node receives a pale vermillion horizontal bar showing its 90% maximum-posterior date region.
- **Detail region:** On narrow screens the full-resolution tree keeps a wide intrinsic canvas inside a labelled, keyboard-focusable horizontal scroll region with a visible gold focus outline.
- **Evidence:** The interval-bearing SVG is the displayed figure; the original TreeTime SVG, dated tree, node-date table, confidence summary, and raster export remain in the adjacent manifest.

### Measure Strip

- **Style:** Three equal white cells inside one ink frame, separated by hairlines.
- **Typography:** Uppercase muted labels, 21px values with tabular numerals, and muted 12px units or interpretation.
- **Responsive behavior:** Cells stack at 800px with horizontal dividers.

### Time-Tree Confidence Grid

- **Style:** A square 3×2 ruled grid using the same label, value, unit, and tabular-number treatment as the standard measure strip.
- **Measures:** Clock rate, final dated-tree fit, root estimate and 90% region, internal-node count, median internal-node interval width, and widest interval. Earlier stages label root-to-tip values as exploratory so readers do not confuse the diagnostic fit with the refitted dated-tree result.
- **Responsive behavior:** All six cells become a single column at 800px, with one hairline between each measure.

### Evidence Tables

- **Style:** Collapsed square grid with 14px body copy, uppercase 11px muted headers, and hairline row dividers.
- **Responsive behavior:** Wide tables scroll within labelled, keyboard-focusable flat regions; at 800px they announce the gesture and hold the first column on paper white. Compact two-column tables do not scroll: they stack each comparison and finding as one bordered evidence row with the first cell acting as its label.

### Interpretation Boundary

- **Style:** A pale caution wash behind plain-language limits and review guardrails.
- **Role:** It is an annotation within the evidence wall, not a warning card or a substitute for the decision label.

### Report Completion

- **Style:** A near-black full-width footer aligned to the stage rail, with a compact “REPORT COMPLETE” label and a direct statement of what the supporting-results bundle preserves.
- **Role:** Close the evidence package cleanly without appending irrelevant method explainers to the visual report.

### Motion

- **Stage navigation:** Smooth anchored scrolling preserves orientation between the four stages.
- **Reduced motion:** When the reader requests reduced motion, anchored scrolling switches to immediate movement.

## Do's and Don'ts

### Do:

- **Do** begin each analytical stage with its scientific question and explicit decision state.
- **Do** keep figures large and place their exact source files in an open adjacent evidence manifest.
- **Do** use square sheets, hairline dividers, and strong stage boundaries to organise dense evidence.
- **Do** preserve labelled states, keyboard focus, horizontal-scroll guidance, sticky first columns, and print order.
- **Do** keep taxon names reader-facing, spaced, and italic while retaining machine identifiers in data outputs.
- **Do** distinguish exploratory root-to-tip fit from the final refitted dated-tree fit in every measure label.
- **Do** preserve the figure → meaning or verdict → downloads sequence when evidence layouts stack.
- **Do** embed Archivo and keep figures and report assets self-contained for offline review.
- **Do** expose dated-tree uncertainty in both the phylogeny and the six-measure confidence grid.
- **Do** end with the compact report-completion statement and supporting-results provenance.

### Don't:

- **Don't** convert stages into generic dashboard cards or detach statistics from the question they answer.
- **Don't** use cobalt and vermillion interchangeably: cobalt is evidence; vermillion is a decision.
- **Don't** add rounded corners, shadows, gradients, glass effects, or decorative chrome.
- **Don't** collapse evidence access into an unlabeled icon or hide the only copy of a supporting result.
- **Don't** rely on colour alone or let narrow tables silently clip evidence.
- **Don't** force compact two-column evidence tables into horizontal scrolling on narrow screens.
- **Don't** append irrelevant method explainers to the visual report.
