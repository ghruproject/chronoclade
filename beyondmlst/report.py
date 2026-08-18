"""Human-readable temporal-signal reports."""

from __future__ import annotations

import math
from html import escape
from pathlib import Path


def assess_temporal_signal(
    temporal: dict[str, object], *, p_value_threshold: float
) -> dict[str, object]:
    """Return a stable verdict and reader-facing explanation."""

    observed = temporal["observed"]
    if not isinstance(observed, dict):
        raise ValueError("temporal observed metrics must be a mapping")
    rate = float(observed["rate"])
    requested = int(temporal["requested_randomisations"])
    successful = int(temporal["successful_randomisations"])
    p_value = float(temporal["p_value_r_squared"])

    if requested == 0:
        return {
            "code": "not_assessed",
            "label": "Temporal signal not assessed",
            "reason": "No date randomisations were requested, so no dated tree was generated.",
            "supported": False,
        }
    if successful != requested:
        return {
            "code": "not_assessed",
            "label": "Temporal signal assessment incomplete",
            "reason": (
                f"Only {successful} of {requested} date randomisations completed successfully. "
                "The temporal-signal gate was not evaluated."
            ),
            "supported": False,
        }
    if rate <= 0:
        return {
            "code": "not_supported",
            "label": "Temporal signal not supported",
            "reason": (
                "The estimated evolutionary rate was not positive. Sampling dates do not "
                "support reliable time scaling for this lineage."
            ),
            "supported": False,
        }
    if p_value > p_value_threshold:
        return {
            "code": "not_supported",
            "label": "Temporal signal not supported",
            "reason": (
                f"The observed root-to-tip fit was not stronger than randomised dates at the "
                f"predefined threshold (p={p_value:.3g}; threshold={p_value_threshold:.3g})."
            ),
            "supported": False,
        }
    return {
        "code": "supported",
        "label": "Temporal signal supported",
        "reason": (
            "The lineage has a positive estimated rate and its observed root-to-tip fit was "
            "stronger than expected after randomly permuting sampling dates."
        ),
        "supported": True,
    }


def _metric(temporal: dict[str, object], name: str) -> float:
    observed = temporal["observed"]
    if not isinstance(observed, dict):
        raise ValueError("temporal observed metrics must be a mapping")
    return float(observed[name])


def _format_rate(value: float) -> str:
    return f"{value:.3g} substitutions/site/year"


def _histogram_panel(
    values: list[float],
    *,
    observed: float,
    x: float,
    title: str,
    x_label: str,
    formatter: str,
) -> str:
    panel_width = 430.0
    plot_x = x + 52
    plot_y = 72.0
    plot_width = 340.0
    plot_height = 220.0
    title_text = escape(title)
    label_text = escape(x_label)

    if not values:
        return f"""
        <g>
          <text x="{x + panel_width / 2}" y="36" text-anchor="middle" class="title">{title_text}</text>
          <rect x="{plot_x}" y="{plot_y}" width="{plot_width}" height="{plot_height}" class="plot-bg"/>
          <text x="{x + panel_width / 2}" y="185" text-anchor="middle" class="empty">No completed permutations</text>
          <text x="{x + panel_width / 2}" y="330" text-anchor="middle" class="axis-label">{label_text}</text>
        </g>
        """

    lower = min([*values, observed])
    upper = max([*values, observed])
    if math.isclose(lower, upper):
        padding = abs(lower) * 0.1 or 1.0
        lower -= padding
        upper += padding
    else:
        padding = (upper - lower) * 0.05
        lower -= padding
        upper += padding

    bin_count = min(20, max(5, round(math.sqrt(len(values)))))
    bin_width = (upper - lower) / bin_count
    counts = [0] * bin_count
    for value in values:
        index = min(bin_count - 1, max(0, int((value - lower) / bin_width)))
        counts[index] += 1
    maximum = max(counts) or 1
    bar_width = plot_width / bin_count
    bars: list[str] = []
    for index, count in enumerate(counts):
        height = plot_height * count / maximum
        bars.append(
            f'<rect x="{plot_x + index * bar_width + 1:.2f}" '
            f'y="{plot_y + plot_height - height:.2f}" width="{max(1, bar_width - 2):.2f}" '
            f'height="{height:.2f}" class="bar"/>'
        )

    observed_x = plot_x + ((observed - lower) / (upper - lower)) * plot_width
    lower_text = format(lower, formatter)
    upper_text = format(upper, formatter)
    observed_text = format(observed, formatter)
    return f"""
    <g>
      <text x="{x + panel_width / 2}" y="36" text-anchor="middle" class="title">{title_text}</text>
      <rect x="{plot_x}" y="{plot_y}" width="{plot_width}" height="{plot_height}" class="plot-bg"/>
      {"".join(bars)}
      <line x1="{observed_x:.2f}" y1="{plot_y - 8}" x2="{observed_x:.2f}" y2="{plot_y + plot_height}" class="observed"/>
      <text x="{observed_x:.2f}" y="{plot_y - 16}" text-anchor="middle" class="observed-label">Observed {escape(observed_text)}</text>
      <line x1="{plot_x}" y1="{plot_y + plot_height}" x2="{plot_x + plot_width}" y2="{plot_y + plot_height}" class="axis"/>
      <text x="{plot_x}" y="{plot_y + plot_height + 20}" text-anchor="start" class="tick">{escape(lower_text)}</text>
      <text x="{plot_x + plot_width}" y="{plot_y + plot_height + 20}" text-anchor="end" class="tick">{escape(upper_text)}</text>
      <text x="{x + panel_width / 2}" y="330" text-anchor="middle" class="axis-label">{label_text}</text>
    </g>
    """


def write_randomisation_plot(temporal: dict[str, object], output: Path) -> Path:
    """Write an accessible two-panel SVG comparing observed and permuted metrics."""

    raw_randomised = temporal.get("randomised", [])
    randomised = [value for value in raw_randomised if isinstance(value, dict)]
    r_squared = [float(value["r_squared"]) for value in randomised]
    rates = [float(value["rate"]) for value in randomised]
    observed_r_squared = _metric(temporal, "r_squared")
    observed_rate = _metric(temporal, "rate")
    p_value = float(temporal["p_value_r_squared"])
    requested = int(temporal["requested_randomisations"])
    successful = int(temporal["successful_randomisations"])

    left = _histogram_panel(
        r_squared,
        observed=observed_r_squared,
        x=20,
        title="Root-to-tip fit under randomised dates",
        x_label="R²",
        formatter=".3g",
    )
    right = _histogram_panel(
        rates,
        observed=observed_rate,
        x=470,
        title="Clock rate under randomised dates",
        x_label="Substitutions/site/year",
        formatter=".2e",
    )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="920" height="390" viewBox="0 0 920 390" role="img" aria-labelledby="title description">
  <title id="title">Date-randomisation comparison</title>
  <desc id="description">Histograms compare the observed root-to-tip R squared and clock rate with values from permuted sampling dates.</desc>
  <style>
    text {{ font-family: Inter, ui-sans-serif, system-ui, sans-serif; fill: #243340; }}
    .title {{ font-size: 16px; font-weight: 700; }}
    .axis-label {{ font-size: 13px; font-weight: 600; }}
    .tick, .empty {{ font-size: 12px; fill: #607080; }}
    .plot-bg {{ fill: #f5f8fa; stroke: #d6e0e6; }}
    .bar {{ fill: #6285a3; opacity: 0.82; }}
    .observed {{ stroke: #c6403d; stroke-width: 3; }}
    .observed-label {{ font-size: 12px; font-weight: 700; fill: #a52f2c; }}
    .axis {{ stroke: #607080; stroke-width: 1; }}
    .footer {{ font-size: 13px; fill: #40515e; }}
  </style>
  {left}
  {right}
  <text x="460" y="372" text-anchor="middle" class="footer">{successful} of {requested} permutations completed · root-to-tip randomisation p={p_value:.3g}</text>
</svg>
"""
    output.write_text(svg, encoding="utf-8")
    return output


def _file_link(relative_path: str, label: str) -> str:
    return f'<a href="{escape(relative_path, quote=True)}">{escape(label)}</a>'


def _display_summary(value: object) -> str:
    if not isinstance(value, dict) or not value.get("comparisons"):
        return "No comparisons"
    return (
        f"n={int(value['comparisons'])}; median {float(value['median']):.3g}; "
        f"range {int(value['minimum'])}–{int(value['maximum'])} SNPs"
    )


def _public_health_visual(report: dict[str, object], directory: Path) -> str:
    raw_public_health = report.get("public_health", {})
    public_health = raw_public_health if isinstance(raw_public_health, dict) else {}
    raw_scenario = public_health.get("scenario", {})
    scenario = raw_scenario if isinstance(raw_scenario, dict) else {}
    code = str(scenario.get("code", "indeterminate"))
    label = str(scenario.get("label", "Indeterminate"))
    confidence = str(scenario.get("confidence", "low"))
    reasons = scenario.get("reasons", [])
    actions = scenario.get("recommended_follow_up", [])
    evidence = scenario.get("evidence", [])
    if not isinstance(reasons, list):
        reasons = []
    if not isinstance(actions, list):
        actions = []
    if not isinstance(evidence, list):
        evidence = []
    evidence_rows = []
    for item in evidence:
        if not isinstance(item, dict):
            continue
        evidence_rows.append(
            "<tr>"
            f"<td>{escape(str(item.get('id', '')).replace('_', ' ').title())}</td>"
            f'<td><span class="evidence {escape(str(item.get("status", "unavailable")))}">'
            f"{escape(str(item.get('status', 'unavailable')).title())}</span></td>"
            f"<td>{escape(str(item.get('finding', '')))}</td>"
            "</tr>"
        )

    raw_distances = public_health.get("distance_summary", {})
    distances = raw_distances if isinstance(raw_distances, dict) else {}
    raw_categories = distances.get("categories", {})
    categories = raw_categories if isinstance(raw_categories, dict) else {}
    category_labels = (
        ("focal_focal", "Focal–focal"),
        ("focal_context", "Focal–context"),
        ("same_month", "Focal pairs in the same month"),
        ("same_year", "Focal pairs in the same year"),
        ("between_years", "Focal pairs between years"),
        ("within_candidate_groups", "Within candidate local groups"),
        ("between_candidate_groups", "Between candidate local groups"),
        ("same_patient", "Within-patient focal pairs"),
        ("different_patients", "Between-patient focal pairs"),
    )
    distance_rows = "".join(
        f"<tr><td>{escape(label_text)}</td><td>{escape(_display_summary(categories.get(key)))}</td></tr>"
        for key, label_text in category_labels
    )

    raw_topology = public_health.get("topology", {})
    topology = raw_topology if isinstance(raw_topology, dict) else {}
    raw_groups = topology.get("groups", [])
    groups = raw_groups if isinstance(raw_groups, list) else []
    group_rows = []
    for group in groups:
        if not isinstance(group, dict):
            continue
        sample_ids = group.get("sample_ids", [])
        if not isinstance(sample_ids, list):
            sample_ids = []
        locations = group.get("locations", [])
        if not isinstance(locations, list):
            locations = []
        nearest = group.get("nearest_context")
        nearest_text = "—"
        if isinstance(nearest, dict):
            nearest_text = (
                f"{nearest.get('sample_id', '')} ({nearest.get('clonal_snps', '—')} SNPs)"
            )
        group_rows.append(
            "<tr>"
            f"<td>{escape(str(group.get('group_id', '')))}</td>"
            f"<td>{int(group.get('sample_count', 0))}</td>"
            f"<td>{int(group.get('patient_count', 0)) or '—'}</td>"
            f"<td>{escape(str(group.get('first_collection_date', '') or '—'))} – "
            f"{escape(str(group.get('last_collection_date', '') or '—'))}</td>"
            f"<td>{escape(', '.join(str(value) for value in locations) or '—')}</td>"
            f"<td>{escape(_display_summary(group.get('within_group_clonal_snps')))}</td>"
            f"<td>{escape(nearest_text)}</td>"
            f"<td><details><summary>View IDs</summary>{escape(', '.join(str(value) for value in sample_ids))}</details></td>"
            "</tr>"
        )
    group_table = (
        "<table><thead><tr><th>Candidate group</th><th>Focal genomes</th><th>Patients</th>"
        "<th>Date span</th><th>Locations</th><th>Within-group clonal distance</th>"
        "<th>Nearest final context</th><th>Samples</th></tr></thead>"
        f"<tbody>{''.join(group_rows)}</tbody></table>"
        if group_rows
        else "<p>No topology-defined focal groups were available.</p>"
    )

    raw_sensitivity = public_health.get("patient_sensitivity", {})
    sensitivity = raw_sensitivity if isinstance(raw_sensitivity, dict) else {}
    excluded = sensitivity.get("excluded_repeated_patient_samples", [])
    if not isinstance(excluded, list):
        excluded = []
    sensitivity_text = str(
        sensitivity.get("interpretation", "Patient-level sensitivity was not available.")
    )
    if not bool(sensitivity.get("patient_metadata_complete")):
        sensitivity_text += " Patient identifiers are incomplete for the focal samples."
    if not excluded:
        exclusion_text = (
            "No repeated-patient samples would be excluded in the deterministic "
            "one-isolate-per-patient view."
        )
    elif len(excluded) == 1:
        exclusion_text = (
            "One repeated-patient sample would be excluded in the deterministic "
            "one-isolate-per-patient view."
        )
    else:
        exclusion_text = (
            f"{len(excluded)} repeated-patient samples would be excluded in the deterministic "
            "one-isolate-per-patient view."
        )

    heatmap = directory / "clonal_snp_heatmap.svg"
    heatmap_visual = (
        '<img src="clonal_snp_heatmap.svg" alt="Heatmap of recombination-filtered pairwise SNP counts">'
        if heatmap.is_file()
        else '<p class="missing">The clonal SNP heatmap was not available.</p>'
    )
    evidence_table = (
        "<table><thead><tr><th>Evidence</th><th>Status</th><th>Finding</th></tr></thead>"
        f"<tbody>{''.join(evidence_rows)}</tbody></table>"
        if evidence_rows
        else "<p>No structured evidence ledger was available.</p>"
    )
    return f"""
  <section class="scenario {escape(code)}">
    <span class="eyebrow">Working public-health interpretation</span>
    <strong>{escape(label)}</strong>
    <p><b>Confidence: {escape(confidence)}.</b> This automated summary is provisional and requires epidemiological review.</p>
    <p><b>Rule applied:</b> {escape(str(scenario.get("decision_rule", "No scenario rule was available.")))}</p>
    <ul>{"".join(f"<li>{escape(str(reason))}</li>" for reason in reasons)}</ul>
    <p class="guardrail">{escape(str(scenario.get("guardrail", "This analysis does not establish direct transmission.")))}</p>
  </section>

  <section class="card action">
    <h2>Recommended follow-up</h2>
    <ul>{"".join(f"<li>{escape(str(action))}</li>" for action in actions)}</ul>
  </section>

  <section class="card">
    <h2>Evidence ledger</h2>
    <p>Every part of the working interpretation is exposed here rather than combined into an opaque score.</p>
    {evidence_table}
  </section>

  <section class="card">
    <h2>Candidate local groups</h2>
    <p>{escape(str(topology.get("interpretation", "No topology summary was available.")))}</p>
    {group_table}
  </section>

  <section class="card">
    <h2>Recombination-filtered genomic distances</h2>
    <p>Counts use only pair-specific A/C/G/T sites after ClonalFrameML recombination filtering. No universal SNP threshold has been applied.</p>
    {heatmap_visual}
    <table><thead><tr><th>Comparison</th><th>Clonal SNP summary</th></tr></thead><tbody>{distance_rows}</tbody></table>
    <p>{_file_link("clonal_pairwise_distances.tsv", "Download pairwise SNPs and callable sites")} · {_file_link("clonal_snp_matrix.tsv", "Download SNP matrix")} · {_file_link("pairwise_callable_sites.tsv", "Download callable-site matrix")}</p>
  </section>

  <section class="card">
    <h2>Patient-level sensitivity</h2>
    <p>{escape(sensitivity_text)}</p>
    <p>{escape(exclusion_text)}</p>
  </section>
"""


def write_lineage_report(
    report: dict[str, object], *, directory: Path, p_value_threshold: float
) -> Path:
    """Write a standalone, reader-facing HTML report for one lineage."""

    temporal = report["temporal_signal"]
    if not isinstance(temporal, dict):
        raise ValueError("lineage report lacks temporal metrics")
    assessment = assess_temporal_signal(temporal, p_value_threshold=p_value_threshold)
    status = str(assessment["code"])
    raw_context = report.get("context", {})
    context = raw_context if isinstance(raw_context, dict) else {}
    local_count = int(context.get("local_samples", 0))
    context_count = int(context.get("context_samples", 0))
    context_locations = context.get("context_locations", [])
    if not isinstance(context_locations, list):
        context_locations = []
    public_health_visual = _public_health_visual(report, directory)
    observed_rate = _metric(temporal, "rate")
    observed_r_squared = _metric(temporal, "r_squared")
    p_value = float(temporal["p_value_r_squared"])
    requested = int(temporal["requested_randomisations"])
    successful = int(temporal["successful_randomisations"])

    root_plot = directory / "clock" / "root_to_tip_regression.svg"
    timetree_plot = directory / "timetree" / "timetree.svg"
    randomisation_plot = directory / "date_randomisation.svg"
    write_randomisation_plot(temporal, randomisation_plot)

    root_visual = (
        '<img src="clock/root_to_tip_regression.svg" alt="TreeTime root-to-tip regression plot">'
        if root_plot.is_file()
        else '<p class="missing">Root-to-tip plot was not available. See the clock output files.</p>'
    )
    timetree_visual = ""
    if bool(assessment["supported"]) and timetree_plot.is_file():
        timetree_visual = """
        <section class="card">
          <h2>Time-scaled phylogeny</h2>
          <p>This visual is shown because the predefined temporal-signal gate passed. Grey
          horizontal bars show TreeTime 90% confidence intervals for inferred node dates.</p>
          <img src="timetree/timetree.svg" alt="TreeTime time-scaled phylogeny">
        </section>
        """

    outliers = directory / "clock" / "outliers.tsv"
    outlier_note = (
        "Automatic clock-outlier filtering was disabled for the temporal-signal gate, so all "
        "supplied genomes contributed. Review clock deviations in the root-to-tip table."
    )
    if outliers.is_file():
        rows = [line for line in outliers.read_text(encoding="utf-8").splitlines() if line]
        count = max(0, len(rows) - 1)
        outlier_note = (
            f"TreeTime reported {count} potential clock outlier(s). Review the outlier table."
        )

    outputs = [
        _file_link("clock/rtt.csv", "Root-to-tip data"),
        _file_link("clock/molecular_clock.txt", "Clock statistics"),
        _file_link("temporal_signal.json", "Date-randomisation results"),
        _file_link("clonalframeml.labelled_tree.newick", "Recombination-corrected tree"),
        _file_link("clonalframeml.importation_status.txt", "Inferred recombination events"),
        _file_link("public_health_evidence.json", "Machine-readable public-health evidence"),
        _file_link("clonal_pairwise_distances.tsv", "Clonal SNPs and callable sites"),
        _file_link("report.json", "Machine-readable lineage report"),
    ]
    if outliers.is_file():
        outputs.append(_file_link("clock/outliers.tsv", "Clock outliers"))
    if bool(assessment["supported"]):
        outputs.append(_file_link("timetree/timetree.nexus", "Dated tree"))
    if bool(context.get("manifest_available")):
        outputs.append(_file_link("context_manifest.tsv", "Frozen context manifest"))

    raw_nearest = context.get("nearest_screening_contexts", [])
    nearest = raw_nearest if isinstance(raw_nearest, list) else []
    nearest_rows: list[str] = []
    for value in nearest:
        if not isinstance(value, dict):
            continue
        distance = value.get("min_ska_distance")
        distance_text = "—" if distance is None else f"{float(distance):.3g}"
        nearest_rows.append(
            "<tr>"
            f"<td>{escape(str(value.get('sample_id', '')))}</td>"
            f"<td>{escape(str(value.get('country', '') or 'Unknown'))}</td>"
            f"<td>{escape(str(value.get('collection_date', '') or 'Unknown'))}</td>"
            f"<td>{escape(str(value.get('nearest_focal', '')))}</td>"
            f"<td>{escape(distance_text)}</td>"
            "</tr>"
        )
    if context_count:
        nearest_table = (
            "<table><thead><tr><th>Context genome</th><th>Country</th><th>Date</th>"
            "<th>Nearest focal isolate in screen</th><th>SKA screening distance</th>"
            "</tr></thead><tbody>" + "".join(nearest_rows) + "</tbody></table>"
            if nearest_rows
            else "<p>The context set was supplied without a beyondMLST selection manifest.</p>"
        )
        context_visual = f"""
        <section class="card">
          <h2>Public contextual genomes</h2>
          <p>{escape(str(context.get("interpretation", "")))}</p>
          <p><strong>{context_count}</strong> context genomes cover <strong>{len(context_locations)}</strong>
          reported locations. These are the nearest genomes within the downloaded screening pool,
          not necessarily the nearest in all public data. SKA distances were used only to select
          candidates; final relatedness must be interpreted from the recombination-corrected tree.</p>
          {nearest_table}
        </section>
        """
    else:
        context_visual = f"""
        <section class="card caveat">
          <h2>Public contextual genomes</h2>
          <p>{escape(str(context.get("interpretation", "No contextual genomes were supplied.")))}</p>
          <p>The report can describe structure among local isolates, but it cannot distinguish
          persistence from repeated introductions without an appropriate same-lineage context set.</p>
        </section>
        """

    output = directory / "report.html"
    output.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>beyondMLST analysis report — {escape(str(report["species"]))} {escape(str(report["lineage"]))}</title>
  <style>
    :root {{ color-scheme: light; --ink:#1d2b34; --muted:#5b6b76; --line:#d8e1e7; --paper:#fff; --wash:#f2f6f8; --accent:#265f7d; --good:#157347; --bad:#a33a36; --unknown:#8a6420; }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; background:var(--wash); color:var(--ink); font:16px/1.55 Inter, ui-sans-serif, system-ui, sans-serif; }}
    main {{ width:min(1120px, calc(100% - 32px)); margin:32px auto 64px; }}
    header {{ background:#17394d; color:white; padding:30px 34px; border-radius:16px; box-shadow:0 8px 24px #17394d18; }}
    header p {{ margin:6px 0 0; color:#d5e5ed; }}
    h1 {{ margin:0; font-size:clamp(28px,4vw,42px); line-height:1.15; }}
    h2 {{ margin:0 0 10px; font-size:22px; }}
    p {{ margin:8px 0; }}
    .verdict {{ margin-top:18px; border-left:7px solid; padding:18px 22px; border-radius:10px; background:var(--paper); }}
    .verdict.supported {{ border-color:var(--good); }}
    .verdict.not_supported {{ border-color:var(--bad); }}
    .verdict.not_assessed {{ border-color:var(--unknown); }}
    .verdict strong {{ display:block; font-size:23px; }}
    .scenario {{ margin-top:18px; border:2px solid var(--line); border-left:9px solid var(--unknown); padding:22px 26px; border-radius:13px; background:var(--paper); }}
    .scenario.persistent_local_lineage {{ border-left-color:var(--good); }}
    .scenario.multiple_introductions {{ border-left-color:#7b4bb7; }}
    .scenario.mixed {{ border-left-color:#c46b19; }}
    .scenario strong {{ display:block; margin-top:3px; font-size:25px; }}
    .eyebrow {{ color:var(--muted); font-size:12px; font-weight:800; letter-spacing:.07em; text-transform:uppercase; }}
    .guardrail {{ border-top:1px solid var(--line); margin-top:16px; padding-top:13px; color:var(--muted); }}
    .action {{ border-left:7px solid var(--accent); }}
    .evidence {{ display:inline-block; border-radius:999px; padding:3px 8px; font-size:12px; font-weight:800; }}
    .evidence.supported {{ color:#0d623d; background:#dff4e9; }}
    .evidence.suggestive {{ color:#76510d; background:#fff0c8; }}
    .evidence.contradicted {{ color:#8a2e2b; background:#fde3e2; }}
    .evidence.unavailable {{ color:#53616b; background:#e8edf0; }}
    .metrics {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(175px,1fr)); gap:12px; margin:18px 0; }}
    .metric, .card {{ background:var(--paper); border:1px solid var(--line); border-radius:13px; box-shadow:0 4px 14px #17394d0a; }}
    .metric {{ padding:15px 17px; }}
    .metric span {{ display:block; color:var(--muted); font-size:13px; text-transform:uppercase; letter-spacing:.04em; }}
    .metric strong {{ display:block; margin-top:4px; font-size:19px; }}
    .card {{ margin-top:18px; padding:24px; overflow:auto; }}
    .grid {{ display:grid; grid-template-columns:1fr; gap:18px; }}
    .grid .card {{ margin-top:0; }}
    img {{ display:block; width:100%; height:auto; margin-top:14px; border:1px solid var(--line); border-radius:9px; background:white; }}
    a {{ color:var(--accent); font-weight:600; }}
    ul {{ margin-bottom:0; }}
    table {{ width:100%; border-collapse:collapse; margin-top:16px; }}
    th,td {{ padding:10px; text-align:left; border-bottom:1px solid var(--line); }}
    th {{ color:var(--muted); font-size:12px; text-transform:uppercase; }}
    .caveat {{ background:#fff8df; border-color:#e8d185; }}
    .missing {{ color:var(--muted); font-style:italic; padding:30px 0; }}
    footer {{ color:var(--muted); margin-top:22px; font-size:14px; }}
    @media (max-width:600px) {{ main {{ width:min(100% - 20px,1120px); margin-top:10px; }} header,.card {{ padding:20px; }} .grid {{ grid-template-columns:1fr; }} }}
  </style>
</head>
<body>
<main>
  <header>
    <h1>{escape(str(report["species"]))} · {escape(str(report["lineage"]))}</h1>
    <p>beyondMLST public-health evidence report</p>
  </header>

  {public_health_visual}

  <section class="verdict {escape(status)}">
    <span class="eyebrow">Temporal analysis</span>
    <strong>{escape(str(assessment["label"]))}</strong>
    <p>{escape(str(assessment["reason"]))}</p>
  </section>

  <section class="metrics" aria-label="Key statistics">
    <div class="metric"><span>Genomes</span><strong>{int(report["sample_count"])}</strong></div>
    <div class="metric"><span>Local genomes</span><strong>{local_count or "—"}</strong></div>
    <div class="metric"><span>Context genomes</span><strong>{context_count}</strong></div>
    <div class="metric"><span>Distinct dates</span><strong>{int(report["distinct_dates"])}</strong></div>
    <div class="metric"><span>Clock rate<br>subs/site/year</span><strong>{observed_rate:.3g}</strong></div>
    <div class="metric"><span>Root-to-tip R²</span><strong>{observed_r_squared:.3g}</strong></div>
    <div class="metric"><span>Randomisation p</span><strong>{p_value:.3g}</strong></div>
    <div class="metric"><span>Permutations</span><strong>{successful}/{requested}</strong></div>
  </section>

  {context_visual}

  <div class="grid">
    <section class="card">
      <h2>Root-to-tip regression</h2>
      <p>Each point is a sampled genome. The fitted line asks whether genetic divergence accumulates with sampling time. R² alone is descriptive and does not establish temporal signal.</p>
      {root_visual}
    </section>
    <section class="card">
      <h2>Date-randomisation test</h2>
      <p>The observed result is compared with analyses in which collection dates were permuted among genomes. A convincing temporal pattern should outperform this null distribution.</p>
      <img src="date_randomisation.svg" alt="Observed TreeTime statistics compared with randomised sampling dates">
    </section>
  </div>

  {timetree_visual}

  <section class="card">
    <h2>Clock diagnostics</h2>
    <p>{escape(outlier_note)}</p>
  </section>

  <section class="card caveat">
    <h2>Interpretation boundary</h2>
    <p>Temporal signal means that these sampling dates contain information about evolutionary rate and node timing under this model. It does not prove direct transmission, local circulation or a definitive number of introductions. Those interpretations also depend on recombination filtering, contextual sampling and epidemiological metadata.</p>
  </section>

  <section class="card">
    <h2>Download and audit</h2>
    <ul>{"".join(f"<li>{link}</li>" for link in outputs)}</ul>
  </section>

  <footer>Generated by beyondMLST. Review the complete diagnostics before using dated estimates.</footer>
</main>
</body>
</html>
""",
        encoding="utf-8",
    )
    return output


def write_summary_report(summary: dict[str, object], *, output: Path) -> Path:
    """Write a compact landing page linking all lineage reports."""

    rows: list[str] = []
    lineages = summary.get("lineages", [])
    if not isinstance(lineages, list):
        raise ValueError("summary lineages must be a list")
    for lineage in lineages:
        if not isinstance(lineage, dict):
            continue
        status = str(lineage.get("temporal_status", lineage.get("status", "unknown")))
        label = status.replace("_", " ").title()
        metrics = lineage.get("temporal_signal", {})
        rate = "—"
        r_squared = "—"
        p_value = "—"
        if isinstance(metrics, dict) and "observed" in metrics:
            rate = _format_rate(_metric(metrics, "rate"))
            r_squared = f"{_metric(metrics, 'r_squared'):.3g}"
            p_value = f"{float(metrics['p_value_r_squared']):.3g}"
        slug = str(lineage["slug"])
        lineage_context = lineage.get("context", {})
        context_count = (
            int(lineage_context.get("context_samples", 0))
            if isinstance(lineage_context, dict)
            else 0
        )
        lineage_public_health = lineage.get("public_health", {})
        public_health = lineage_public_health if isinstance(lineage_public_health, dict) else {}
        lineage_scenario = public_health.get("scenario", {})
        scenario = lineage_scenario if isinstance(lineage_scenario, dict) else {}
        scenario_label = str(scenario.get("label", "Not assessed"))
        scenario_code = str(scenario.get("code", "indeterminate"))
        report_link = (
            f'<a href="{escape(slug, quote=True)}/report.html">Open report</a>'
            if "temporal_signal" in lineage
            else "Not analysed"
        )
        rows.append(
            f"<tr><td>{escape(str(lineage['species']))}</td><td>{escape(str(lineage['lineage']))}</td>"
            f"<td>{int(lineage['sample_count'])}</td><td>{context_count}</td>"
            f'<td><span class="status {escape(scenario_code)}">{escape(scenario_label)}</span></td>'
            f'<td><span class="status {escape(status)}">{escape(label)}</span></td>'
            f"<td>{escape(rate)}</td><td>{escape(r_squared)}</td><td>{escape(p_value)}</td><td>{report_link}</td></tr>"
        )

    report_path = output / "index.html"
    report_path.write_text(
        f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>beyondMLST analysis report</title>
<style>
body{{margin:0;background:#f2f6f8;color:#1d2b34;font:16px/1.5 Inter,ui-sans-serif,system-ui,sans-serif}}main{{width:min(1200px,calc(100% - 32px));margin:32px auto}}header{{background:#17394d;color:white;padding:30px;border-radius:16px}}h1{{margin:0}}header p{{color:#d5e5ed}}section{{background:white;border:1px solid #d8e1e7;border-radius:13px;margin-top:18px;padding:22px;overflow:auto}}table{{width:100%;border-collapse:collapse;min-width:900px}}th,td{{padding:12px;text-align:left;border-bottom:1px solid #e4eaee}}th{{font-size:13px;text-transform:uppercase;color:#5b6b76}}a{{color:#265f7d;font-weight:700}}.status{{font-weight:700}}.supported{{color:#157347}}.not_supported{{color:#a33a36}}.not_assessed{{color:#8a6420}}footer{{margin-top:20px;color:#5b6b76}}
</style></head><body><main><header><h1>beyondMLST analysis report</h1><p>Recombination-aware temporal and contextual bacterial phylogenetics</p></header>
<section><h2>Lineage results</h2><table><thead><tr><th>Species</th><th>Lineage</th><th>Genomes</th><th>Context</th><th>Working interpretation</th><th>Temporal result</th><th>Rate</th><th>R²</th><th>p</th><th></th></tr></thead><tbody>{"".join(rows)}</tbody></table></section>
<section><h2>How to interpret this report</h2><p>The working interpretation combines recombination-filtered distances, corrected topology, longitudinal sampling and public context. It is provisional and must be reviewed with epidemiological information. A supported temporal signal permits time scaling under the fitted model but is not required for the genomic scenario assessment.</p></section>
<footer>Generated by beyondMLST.</footer></main></body></html>
""",
        encoding="utf-8",
    )
    return report_path
