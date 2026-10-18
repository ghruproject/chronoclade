"""Result-only offline reports preserving the established full-report style."""

from collections import Counter
from html import escape
from pathlib import Path
import json
from math import isfinite
from os.path import relpath
from urllib.parse import quote

from chronoclade.report_components.tables import write_report_html
from chronoclade.report_components.styles import report_styles
from chronoclade.location_network.widget import interactive_network_html, widget_assets
from chronoclade.location_network.state_selector import state_network_selector, state_network_panel
from chronoclade.report_components.profile_summary import _reported
from chronoclade.sample_labels import sample_label


def _root_to_tip_annotation(assessment):
    rate = ((assessment.get("temporal") or {}).get("observed") or {}).get("rate")
    nonpositive = isinstance(rate, (int, float)) and not isinstance(rate, bool) and rate <= 0
    return (assessment.get("supported") is not True or nonpositive,
            'The plotted "root date" is an extrapolated regression intercept, '
            'not an accepted MRCA estimate.')


def _figure_title(name, label, assessment):
    if Path(name).name != "root_to_tip_regression.svg":
        return label
    unsupported, explanation = _root_to_tip_annotation(assessment)
    return label + ". " + ("Dating unsupported. " if unsupported else "") + explanation


def _figure_html(name, label, assessment):
    """Distinguish a raw regression intercept from a supported dating result."""
    caption = escape(label)
    alternative = _figure_title(name, label, assessment)
    if Path(name).name == "root_to_tip_regression.svg":
        unsupported, explanation = _root_to_tip_annotation(assessment)
        if unsupported:
            caption += (
                '<p class="diagnostic-warning" role="alert"><strong>Diagnostic only: '
                'dating unsupported.</strong> ' + escape(explanation)
                + ' Do not use it as an ancestor date or compare it with published '
                'dating results.</p>'
            )
        else:
            caption += '<p>' + escape(explanation) + ' Final dating results are reported separately.</p>'
    image = f'<img src="{escape(name, quote=True)}" alt="{escape(alternative, quote=True)}">'
    if Path(name).name == "root_to_tip_regression.svg":
        image = (f'<a href="{escape(name, quote=True)}" '
                 f'data-figure-title="{escape(alternative, quote=True)}">{image}</a>')
    return f'<figure>{image}<figcaption>{caption}</figcaption></figure>'


def _temporal_screen_summary(assessment):
    """Display saved permutation evidence without changing the dating decision."""
    temporal = assessment.get("temporal") or {}
    successful = temporal.get("successful_randomisations")
    requested = temporal.get("requested_randomisations")
    p_value = temporal.get("p_value_r_squared")
    threshold = assessment.get("p_value_threshold")
    def count(value):
        return type(value) is int and value >= 0

    def number(value):
        return type(value) in (int, float) and isfinite(value)
    parts = []
    if number(p_value) and 0 <= p_value <= 1:
        value = format(p_value, ".6g")
        if count(successful):
            denominator = successful + 1
            numerator = p_value * denominator
            if float(numerator).is_integer():
                value += f" ({int(numerator)}/{denominator})"
        parts.append("Empirical root-to-tip p: " + value)
    if number(threshold):
        parts.append("threshold: " + format(threshold, ".6g"))
    if count(successful) and count(requested):
        parts.append(f"successful permutations: {successful}/{requested}")
    if not parts:
        return ""
    return '<p class="temporal-screen-summary">' + escape("; ".join(parts)) + '.</p>'


def write_report(
    path, *, title, result, rows=(), network=None, host_network=None, figures=(), extra="", overview=None, report_label="TREE ANALYSIS REPORT",
    figure_heading="Tree figures", assessment_first=False,
    assessment_title=None, navigation=(), state_networks=None
):
    path = Path(path)
    optional = [field for field in ("region", "city", "host")
                if any(_reported(row.get(field)) for row in rows)]
    metadata_fields = ["role", "collection_date", "country", *optional]
    counts = "".join(
        f"<h3>{escape(field.title())}</h3><p>"
        + "; ".join(
            f"{escape(str(key))}: {count}"
            for key, count in sorted(Counter(row.get(field) or "Unknown" for row in rows).items())
        )
        + "</p>"
        for field in ("country", *optional)
    )
    table = (
        "<table><thead><tr><th>Genome</th>"
        + "".join(f"<th>{escape(field.replace('_', ' ').title())}</th>" for field in metadata_fields)
        + "</tr></thead><tbody>"
        + "".join(
            "<tr>"
            + f'<td title="{escape(row["sample_id"], quote=True)}">{escape(sample_label(row))}</td>'
            + "".join(
                f"<td>{escape(str(row.get(field) or 'Unknown'))}</td>"
                for field in metadata_fields
            )
            + "</tr>"
            for row in rows
        )
        + "</tbody></table>"
    )
    assessment = result["temporal_assessment"]
    plots = "".join(
        _figure_html(name, label, assessment)
        for name, label in figures
        if (path.parent / name).is_file()
    )
    labels = {row["sample_id"]: sample_label(row) for row in rows}
    nearest = result.get("nearest_contexts_by_corrected_snps", [])
    if nearest:
        plots += "<h3>Nearest selected context by corrected SNP distance</h3><p>The search denominator is this exact selected assembly set. Distances are SNP counts and pair-specific callable sites.</p><div class=\"table-scroll\"><table><thead><tr><th>Input</th><th>Nearest context</th><th>Corrected SNPs</th><th>Callable sites</th></tr></thead><tbody>"
        for pair in nearest:
            plots += (
                "<tr>"
                + "".join(
                    "<td>" + escape(str(value)) + "</td>"
                    for value in (
                        labels.get(pair["focal_sample"], pair["focal_sample"]),
                        labels.get(pair["context_sample"], pair["context_sample"]),
                        pair["clonal_snps"],
                        pair["callable_sites"],
                    )
                )
                + "</tr>"
            )
        plots += "</tbody></table></div>"
    downloads = "".join(
        f'<li><a href="{escape(p.relative_to(path.parent).as_posix(), quote=True)}"'
        + (f' data-figure-title="{escape(_figure_title(p.name, "Root-to-tip temporal diagnostic", assessment), quote=True)}"'
           if p.name == "root_to_tip_regression.svg" else "")
        + f' download>{escape(p.name)}</a></li>'
        for p in sorted(path.parent.rglob("*"))
        if p.is_file() and p != path and p.name != "job.json"
    )
    coherence = result.get("coherence_screening") or {}
    coherence_warning = (
        '<p role="alert"><strong>Coherence warning:</strong> '
        + escape(coherence["warning"]) + "</p>"
        if coherence.get("warning") else ""
    )
    confidence = result.get("clock_confidence") or {}
    if confidence:
        def number(key, digits=".4f"):
            value = confidence.get(key)
            return format(value, digits) if isinstance(value, (int, float)) else "Unavailable"
        extra = (
            '<h3>Fitted time tree</h3><p>TreeTime rate: <strong>'
            + number("rate", ".4g") + '</strong> substitutions/site/year. '
            'This is the final time-tree fit; the screening regression is a separate diagnostic.</p>'
            '<p>Whole-tree root year: <strong>' + number("root_numeric_date")
            + '</strong>; conditional 90% max-posterior region: '
            + number("root_lower_bound") + '–' + number("root_upper_bound") + '.</p>'
            '<p>These within-run intervals do not include uncertainty from selection '
            'of genomes, missing sampling or residual recombination.</p>' + extra
        )
    network_html = ""
    assets = ""
    if network is not None:
        basis = "Dated-tree" if result.get("dating_status") == "dated" else "Corrected-tree"
        country_panel = f"<h3>{basis} location history</h3>" + interactive_network_html(
            network, "selected-tree-network"
        )
        optional_networks = dict(state_networks or {})
        if host_network:
            optional_networks.setdefault("host", host_network)
        panels = {
            field: state_network_panel(data, basis=basis, widget_id=f"selected-tree-{field}-network")
            for field, data in optional_networks.items() if field != "country"
        }
        if not network.get("nodes") and not network.get("directed_edges"):
            country_panel = None
        network_html = state_network_selector(country_panel, panels=panels)
        assets = widget_assets()
    sections = []
    if rows:
        sections.append(
            ("Selected samples", counts + '<div class="table-scroll">' + table + "</div>")
        )
    if plots:
        sections.append((figure_heading, plots))
    if network_html:
        sections.append(("State-change network", network_html))
    sections.extend(
        [
            (
                assessment_title or ("Ensemble assessment"
                if assessment["code"] == "sensitivity"
                else "Temporal assessment"),
                coherence_warning + ("" if assessment_title else f"<p><strong>{escape(assessment['code'])}</strong>: {escape(assessment['reason'])}</p>") + f"<details><summary>Saved assessment</summary><pre>{escape(json.dumps(assessment, indent=2))}</pre></details>"
                + extra,
            ),
            ("Downloads", "<ul>" + downloads + "</ul>"),
        ]
    )
    if assessment_first:
        sections.insert(0, sections.pop(-2))
    if overview is None:
        overview = f"{len(result['selected_sample_ids'])} exact selected genomes"
        if result.get("branch_units"):
            overview += " · " + result["branch_units"]
        elif result.get("dating_status") == "unsupported":
            overview += " · dated tree unavailable"
    supported = assessment.get("supported", False)
    verdict_title = (
        "Dating unsupported"
        if result.get("dating_status") == "unsupported"
        else "Temporal signal supported"
        if supported
        else "Temporal signal not supported"
    )
    verdict_class = "supported" if supported else "not_supported"
    if assessment["code"] == "sensitivity":
        verdict_title = "Selection sensitivity assessment"
        verdict_class = ""
    if assessment_title is not None:
        verdict_title = assessment_title
    verdict = (
        f'<div class="verdict {verdict_class}">'
        f"<strong>{escape(verdict_title)}</strong><p>{escape(assessment['reason'])}</p>"
        + _temporal_screen_summary(assessment) + "</div>"
    )
    report_links = "".join(
        f'<a href="{escape(quote(relpath(target, path.parent), safe="/"), quote=True)}">{escape(label)}</a>'
        for label, target in navigation if Path(target).is_file()
    )
    navigation_html = (
        '<nav class="contents" aria-label="Report navigation">' + report_links
        + ('<a href="#state-network">State network</a>' if network_html else '')
        + '<a href="#evidence">Evidence files</a></nav>'
    )
    body = []
    for index, (label, content) in enumerate(sections, 1):
        anchor = (' id="evidence"' if label == "Downloads" else
                  ' id="state-network"' if label == "State-change network" else "")
        body.append(
            f'<section class="stage"{anchor}><div class="stage-index"><span>{index}</span></div><div class="stage-body"><h2>{escape(label)}</h2>{content}</div></section>'
        )
    body = "".join(body)
    write_report_html(path,
        f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{escape(title)}</title><style>{report_styles()}pre{{white-space:pre-wrap;overflow-wrap:anywhere;max-height:500px;overflow:auto}}</style></head><body><main class="shell"><header class="identity"><div class="identity-mark">ChronoClade · {escape(report_label)}</div><div class="identity-copy"><h1>{escape(title)}</h1><p>{escape(overview)}</p>{verdict}</div></header>{navigation_html}{body}</main>{assets}</body></html>',
        encoding="utf-8",
    )
    return path
