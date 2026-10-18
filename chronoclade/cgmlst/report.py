"""The existing report shell presents saved profile partitions and selections."""

import html
from pathlib import Path
from chronoclade.report_components.tables import write_report_html
from chronoclade.report_components.styles import report_styles
from chronoclade.report_components.profile_summary import _read_json, _selection_path


def _tree_eligibility(row, output):
    """Read the exact saved selections rather than assuming the pool was selected."""
    path = Path(output) / row["ensemble"]
    ensemble = _read_json(path) or {}
    counts = []
    for run in ensemble.get("selections", []):
        selection_path = _selection_path(path.parent, run)
        selection = _read_json(selection_path) if selection_path else None
        if selection is None:
            return "Unknown: inspect selection manifests"
        selected = set(selection.get("selected_sample_ids") or [])
        if not selected:
            selected = set(selection.get("query_ids") or []) | set(selection.get("selected_context_ids") or [])
        counts.append(len(selected))
    if not counts:
        return "Unknown: no readable selections"
    eligible = sum(count >= 3 for count in counts)
    return (f"Eligible: {eligible}/{len(counts)} selections (minimum 3 genomes)" if eligible == len(counts)
            else f"Blocked: {len(counts) - eligible}/{len(counts)} selections have fewer than 3 genomes")


def write_partition_index(
    products, audit, *, records_count, lin_level, hiercc_level, status, output,
    partition_by="lineage", display_lin_level=7, profile_context=None,
):
    items = "".join(
        f"<tr><td>{html.escape(row['label'])}</td><td>{row['input_count']}</td>"
        f'<td>{row["context_count"]}</td><td>{html.escape(_tree_eligibility(row, output))}</td>'
        f'<td><a href="{html.escape(row["report"])}">Report</a> · '
        f'<a href="{html.escape(row["ensemble"])}">Selections</a></td></tr>'
        for row in products
    )
    unresolved_scope = audit.get("unresolved_context_scope_inputs") or []
    scope_level = audit.get("context_lin_level")
    reason_labels = {
        "missing_or_unresolved_context_scope_typing":
            f"Requested LIN{scope_level} public context is unresolved"
    }
    scope_notice = (
        f"<p>The requested LIN{scope_level} context cannot be defined for these inputs. "
        "Public relatives could not be assessed under that requested scope; "
        "zero retained contexts does not mean no relatives exist.</p>"
        if unresolved_scope else ""
    )
    unresolved = "".join(
        f"<li>{html.escape(row['sample_id'])}: {html.escape(reason_labels.get(row['reason'], row['reason']))}</li>"
        for row in audit["unresolved_inputs"] + audit.get("unresolved_context_scope_inputs", [])
    )
    counts = "".join(
        f"<tr><td>{row['level']}</td><td>ST{html.escape(row['mlst_st'])}</td>"
        f"<td>{html.escape('.'.join(map(str, row['prefix'])))}</td>"
        f"<td>{row['input_count']}</td><td>{row['context_count']}</td></tr>"
        for row in audit["available_level_counts"]
    )
    count_table = (
        (
            "<h2>Available LIN depths</h2><p>Counts use the full prepared dataset before "
            "the optional context-prefix scope and analysis partition filters.</p>"
            "<table><thead><tr><th>Level</th><th>ST</th>"
            "<th>Prefix</th><th>Input</th><th>Context</th></tr></thead><tbody>"
            + counts
            + "</tbody></table>"
        )
        if counts
        else ""
    )
    from chronoclade.report_components.relatedness import scope_summary
    rendered_scope = scope_summary(
        {"partition_by": partition_by, "lin_level": lin_level,
         "context_lin_level": audit.get("context_lin_level"),
         "hiercc_level": hiercc_level, "display_lin_level": display_lin_level},
        audit, provider=profile_context,
    )
    partition_by = audit.get("partition_by", partition_by)
    if audit.get("policy") == "supplied-inputs-species-st-profile-scope":
        partition_by = "st"
    if partition_by == "st":
        context_description = html.escape(
            audit.get("context_scope", "Supplied inputs only; no public context retrieval or lineage refinement")
        )
    elif partition_by == "profiles":
        context_description = (
            "Inputs are grouped by species and compatible cgMLST scheme. "
            "cgLIN/HierCC codes are optional annotations; nearest neighbours use measured profile distances. "
            "Context counts describe compatible prepared profiles before assembly subsampling."
        )
    else:
        context_description = (
            "Context counts describe compatible profiles in the prepared dataset, before assembly subsampling. "
            "Public retrieval provenance is saved with that dataset. Context assemblies are acquired only by the tree stage."
        )
    report_path = output / "index.html"
    write_report_html(report_path,
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>ChronoClade · cgMLST</title><style>{report_styles()}</style></head>"
        '<body><main class="shell"><header class="identity"><div class="identity-mark">ChronoClade · cgMLST ANALYSIS REPORT</div><div class="identity-copy">'
        '<h1>Profile analysis blocks</h1></div></header><section class="stage"><div class="stage-body"><h2>Prepared dataset</h2>'
        f"<p>{records_count} samples; {len(products)} disjoint blocks. Status: {status}.</p>"
        f"<p>{context_description}</p>"
        f"{rendered_scope}{scope_notice}"
        "<table><thead><tr><th>Block</th><th>Input</th><th>Context</th><th>Tree eligibility</th><th>Outputs</th></tr></thead>"
        f"<tbody>{items}</tbody></table>{count_table}<h2>Unresolved inputs</h2><ul>{unresolved or '<li>None</li>'}</ul>"
        '<p><a href="files/index.html">Browse output files by block</a> · '
        '<a href="README.md">Output folder guide</a> · '
        '<a href="partitions.json">Partition counts and exclusions</a></p></div></section></main></body></html>'
    )
    return report_path
