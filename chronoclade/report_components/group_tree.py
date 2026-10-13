"""Offline, full-topology NJ viewer for direct-distance genetic groups."""
from __future__ import annotations

import json
import re
from html import escape
from pathlib import Path

from chronoclade.report_components.styles import report_styles
from chronoclade.sample_labels import display_label

_PALETTE = (
    "#2166ac", "#b2182b", "#1b7837", "#762a83", "#e08214", "#008080",
    "#c51b7d", "#4d4d4d", "#a6d854", "#3288bd", "#f46d43", "#66c2a5",
)
_INPUT_ORIGINS = {"input", "local", "query", "focal"}


def _group_name(value):
    value = str(value or "Ungrouped")
    match = re.search(r"(?:^|_)group_(\d+)$", value, re.IGNORECASE)
    return f"Group {match.group(1)}" if match else value


def _safe_json(value):
    return (json.dumps(value, ensure_ascii=False, separators=(",", ":"))
            .replace("<", "\\u003c").replace(">", "\\u003e")
            .replace("&", "\\u0026"))


def write_group_tree(tree, records, labels, groups, output_path, *, distance_units=None) -> str:
    """Write an offline interactive rectangular phylogram containing every input tip.

    Layout traverses the supplied Bio.Phylo tree without pruning or editing it. The
    SVG retains each signed branch length; negative horizontal segments are dashed.
    """
    records = list(records)
    labels = labels or {}
    groups = groups or {}
    by_id = {str(row.get("sample_id")): row for row in records if row.get("sample_id") is not None}
    if len(by_id) != len(records):
        raise ValueError("Tree records require unique sample_id values")
    tips = tree.get_terminals()
    tip_ids = [tip.name for tip in tips]
    if any(name is None for name in tip_ids) or len(set(tip_ids)) != len(tip_ids):
        raise ValueError("Group tree requires uniquely named tips")
    if set(tip_ids) != set(by_id):
        raise ValueError("Group tree records must match every tree tip")

    group_names = sorted({_group_name(groups.get(name)) for name in tip_ids})
    colors = {name: _PALETTE[index % len(_PALETTE)] for index, name in enumerate(group_names)}
    # Explicit iterative traversal avoids recursion limits on large NJ trees.
    nodes, depths, parent = [], {}, {}
    stack = [(tree.root, 0.0, None)]
    while stack:
        node, depth, parent_node = stack.pop()
        nodes.append(node)
        depths[node] = depth
        if parent_node is not None:
            parent[node] = parent_node
        for child in reversed(node.clades):
            length = child.branch_length
            if length is None:
                length = 0.0
            length = float(length)
            if not (float("-inf") < length < float("inf")):
                raise ValueError("Group tree requires finite branch lengths")
            stack.append((child, depth + length, node))

    tip_order = [node for node in nodes if not node.clades]
    row_by_node = {node: index for index, node in enumerate(tip_order)}
    # Internal y positions are the midpoint of their first and last descendant tips.
    y_positions = {}
    for node in reversed(nodes):
        if not node.clades:
            y_positions[node] = row_by_node[node]
        else:
            child_y = [y_positions[child] for child in node.clades]
            y_positions[node] = (min(child_y) + max(child_y)) / 2

    min_depth = min(depths.values(), default=0.0)
    max_depth = max(depths.values(), default=0.0)
    span = max_depth - min_depth
    chart_left, chart_width = 48.0, 900.0
    label_width = max(160, min(420, max((len(str(labels.get(name, name))) for name in tip_ids), default=10) * 7))
    width = chart_left + chart_width + label_width + 20
    unit = escape(str(distance_units), quote=True) if distance_units else ""
    metadata = []
    for name in tip_ids:
        row = by_id[name]
        label = display_label(str(labels.get(name) or name))
        # Public accession aliases supplied by the caller take precedence over internal IDs.
        metadata.append({
            "id": name,
            "label": label,
            "group": _group_name(groups.get(name)),
            "country": str(row.get("country") or "Unknown country"),
            "host": str(row.get("host") or "Unknown host"),
            "date": str(row.get("collection_date") or row.get("date") or row.get("collection_year") or "Unknown date"),
            "origin": str(row.get("origin") or row.get("role") or "context"),
        })
    data = {"tips": metadata, "groups": group_names, "colors": colors,
            "units": str(distance_units or "")}
    json_data = _safe_json(data)

    def x_for(node):
        if span == 0:
            return chart_left + chart_width * .05
        return chart_left + (depths[node] - min_depth) / span * chart_width

    # A modest initial row spacing keeps even a few hundred tips tractable, while
    # the canvas itself scrolls and offers a label-reading spacing mode.
    row_space = 20
    top = 28
    height = top * 2 + max(1, len(tip_order) - 1) * row_space
    index_for_node = {node: index for index, node in enumerate(tip_order)}
    svg_parts = [f'<svg id="tree" role="img" aria-label="Complete NJ group tree" preserveAspectRatio="none" viewBox="0 0 {width:.1f} {height:.1f}" width="{width:.1f}" height="{height:.1f}"><g id="layout">']
    # Branches are separate line elements to retain exact signed lengths and make
    # negative edges visually explicit without changing topology or values.
    for node in nodes:
        if node.clades:
            low = top + min(y_positions[child] for child in node.clades) * row_space
            high = top + max(y_positions[child] for child in node.clades) * row_space
            svg_parts.append(f'<line class="branch" x1="{x_for(node):.2f}" y1="{low:.2f}" x2="{x_for(node):.2f}" y2="{high:.2f}"/>')
        if node in parent:
            ancestor = parent[node]
            x1, x2 = x_for(ancestor), x_for(node)
            y = top + y_positions[node] * row_space
            negative = (node.branch_length or 0) < 0
            css = "branch negative" if negative else "branch"
            svg_parts.append(f'<line class="{css}" x1="{x1:.2f}" y1="{y:.2f}" x2="{x2:.2f}" y2="{y:.2f}" data-length="{float(node.branch_length or 0):.12g}"/>')
        if not node.clades:
            index = index_for_node[node]
            item = metadata[index]
            shape = "input" if item["origin"] in _INPUT_ORIGINS else "context"
            color = colors[item["group"]]
            y2 = top + y_positions[node] * row_space
            x2 = x_for(node)
            marker = (f'<path class="tip-marker {shape}" d="M {x2-5:.2f} {y2-5:.2f} L {x2+5:.2f} {y2:.2f} L {x2-5:.2f} {y2+5:.2f} Z" fill="{color}"/> '
                      if shape == "input" else
                      f'<circle class="tip-marker {shape}" cx="{x2:.2f}" cy="{y2:.2f}" r="4.5" fill="{color}"/> ')
            band = (f'<rect class="tip-band" x="{x2:.2f}" y="{y2-9:.2f}" '
                    f'width="{width-x2-8:.2f}" height="18" fill="{color}"/>')
            svg_parts.append(f'<g class="tip" tabindex="0" role="button" data-tip="{index}" aria-label="{escape(item["label"], quote=True)}">{band}{marker}<text x="{x2+9:.2f}" y="{y2+4:.2f}">{escape(item["label"])}</text></g>')
    svg_parts.append("</g></svg>")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lineage = all(group.startswith('cgLIN ') for group in group_names)
    note = ("Groups share reported cgLIN prefixes; no additional distance cutoff defines these colours."
            if lineage else "Groups defined from direct profile distances, not tree clades.")
    document = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        '<title>Genetic groups on NJ tree</title><style>' + report_styles() + '''
body{background:white;padding:8px}.tree-controls{display:flex;gap:8px;align-items:center;flex-wrap:wrap}
button{font:inherit;padding:6px 9px;background:white;border:1px solid #aeb4bf;border-radius:4px;cursor:pointer}
.legend{display:flex;gap:6px 10px;flex-wrap:wrap;margin:8px 0}.legend button[aria-pressed="true"]{outline:2px solid #17191f}
.swatch{display:inline-block;width:12px;height:12px;margin-right:5px;vertical-align:-1px}
.tree-canvas{overflow:auto;max-height:68vh;border:1px solid #d0d3da;background:white}
.tree-canvas svg{display:block;width:100%;height:400px;max-width:none;font:12px 'Archivo',sans-serif}.tree-canvas.read-labels svg{width:max(900px,100%);height:auto}
.branch{stroke:#68707c;stroke-width:1;fill:none}.branch.negative{stroke:#888;stroke-dasharray:4 3}
.branch{vector-effect:non-scaling-stroke;stroke-width:.7}.tip-band{opacity:.12}.tip.selected .tip-band{opacity:.4}
.tree-canvas:not(.read-labels) .tip text{display:none}
.tip text{fill:#17191f}.tip{cursor:pointer}.tip:focus{outline:none}.tip:focus .tip-marker{stroke:#ffb800;stroke-width:3}.tip.selected .tip-marker{stroke-width:2}
.tip.dimmed{opacity:.18}.tip-marker{stroke:#222;stroke-width:.65}.internal-note{color:#5f6470;font-size:.9em}
#tip-details{min-height:2.8em;padding:8px;margin:8px 0;border:1px solid #ccd0d8;white-space:pre-wrap}
</style></head><body><h2>Genetic groups on NJ tree</h2>
<div class="tree-controls"><button id="fit" type="button">Fit all</button><button id="read" type="button">Read labels</button><button id="clear" type="button">Clear group highlight</button><label><input id="scroll-group" type="checkbox"> Scroll to first selected group tip</label></div>
<div id="legend" class="legend" aria-label="Genetic group counts"></div><div class="tree-canvas" id="canvas">''' + "".join(svg_parts) + '''</div>
<p id="tip-details" role="status" aria-live="polite">Hover, focus, or select a tip to inspect its metadata.</p>
<p class="internal-note">''' + note + (f' Branch lengths are in {unit}.' if unit else '') + ''' Triangles mark input genomes; circles mark comparisons. Colours show direct-distance group assignment.</p>
<script type="application/json" id="tree-data">''' + json_data + '''</script><script>
(()=>{
const data=JSON.parse(document.getElementById('tree-data').textContent),svg=document.getElementById('tree'),canvas=document.getElementById('canvas');
const details=document.getElementById('tip-details'),legend=document.getElementById('legend'),tips=[...svg.querySelectorAll('.tip')];
const counts=new Map(data.groups.map(g=>[g,0]));data.tips.forEach(t=>counts.set(t.group,(counts.get(t.group)||0)+1));
let selected=null;
function showTip(index){const t=data.tips[index];details.textContent=`${t.label} | Group: ${t.group} | Country: ${t.country} | Host: ${t.host} | Date: ${t.date}`}
function send(group){if(window.parent!==window)window.parent.postMessage({type:'chronoclade-group-highlight',group},'*')}
function highlight(group,notify){selected=group;tips.forEach(el=>{const active=group&&data.tips[Number(el.dataset.tip)].group===group;el.classList.toggle('dimmed',!!group&&!active);el.classList.toggle('selected',!!active)});legend.querySelectorAll('button').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.group===group)));if(notify)send(group)}
data.groups.forEach(group=>{const button=document.createElement('button');button.type='button';button.dataset.group=group;button.setAttribute('aria-pressed','false');const swatch=document.createElement('span');swatch.className='swatch';swatch.style.backgroundColor=data.colors[group];button.append(swatch,document.createTextNode(`${group} (${counts.get(group)})`));button.addEventListener('click',()=>{const next=selected===group?null:group;highlight(next,true);if(next&&document.getElementById('scroll-group').checked){const first=tips.find(el=>data.tips[Number(el.dataset.tip)].group===next);if(first)first.scrollIntoView({block:'center'})}});legend.append(button)});
tips.forEach(el=>{const index=Number(el.dataset.tip);el.addEventListener('pointerenter',()=>showTip(index));el.addEventListener('focus',()=>showTip(index));el.addEventListener('click',()=>{showTip(index);el.classList.toggle('selected')})});
document.getElementById('clear').addEventListener('click',()=>highlight(null,true));
function spacing(space,mode){const n=Math.max(1,data.tips.length),height=56+(n-1)*space,view=svg.getAttribute('viewBox').split(' ');view[3]=String(height);svg.setAttribute('viewBox',view.join(' '));svg.setAttribute('height',String(height));canvas.classList.toggle('read-labels',mode==='read');document.getElementById('layout').setAttribute('transform',`translate(0 28) scale(1 ${space/20}) translate(0 -28)`)}
function fitAll(){const n=Math.max(1,data.tips.length),space=Math.max(1,Math.floor((400-56)/Math.max(1,n-1)));spacing(space,'fit')}
document.getElementById('fit').addEventListener('click',fitAll);
document.getElementById('read').addEventListener('click',()=>spacing(24,'read'));
window.addEventListener('message',event=>{if(event.source!==window.parent||!event.data||event.data.type!=='chronoclade-group-highlight')return;const group=event.data.group===null?null:String(event.data.group);highlight(data.groups.includes(group)?group:null,false)});
fitAll();
})();
</script></body></html>'''
    )
    if lineage:
        document = document.replace('Genetic groups on NJ tree', 'cgLIN groups on NJ tree')
        document = document.replace('Colours show direct-distance group assignment.', 'Colours show reported cgLIN group assignment.')
    output_path.write_text(document, encoding="utf-8")
    return str(output_path)


def tree_viewer(cohorts, directory) -> str:
    """Return report markup embedding generated group-tree pages when available."""
    directory = Path(directory)
    root = directory.resolve()
    if isinstance(cohorts, dict):
        cohorts = cohorts.values()
    frames = []
    for cohort in cohorts or []:
        if not isinstance(cohort, dict):
            continue
        raw_path = cohort.get("tree_groups_html")
        if not raw_path:
            continue
        path = Path(raw_path)
        if not path.is_absolute():
            path = directory / path
        path = path.resolve()
        try:
            src = path.relative_to(root).as_posix()
        except ValueError:
            continue
        if not path.is_file():
            continue
        safe_src = escape(src, quote=True)
        cohort_id = escape(str(cohort.get("cohort_id") or ""), quote=True)
        frames.append(
            f'<section><iframe data-group-cohort="{cohort_id}" title="NJ tree genetic groups" src="{safe_src}" '
            'loading="lazy" style="width:100%;height:700px;border:1px solid #c9ccd4" '
            'sandbox="allow-scripts allow-same-origin"></iframe>'
            f'<p><a href="{safe_src}">Open the complete interactive tree</a></p></section>'
        )
    return "".join(frames)


def linked_group_highlights() -> str:
    """Return a parent-page router for group highlights between same-cohort frames."""
    return '''<script data-group-links>
(()=>{window.addEventListener('message',event=>{
 if(!event.data||event.data.type!=='chronoclade-group-highlight')return;
 const source=[...document.querySelectorAll('iframe[data-group-cohort]')].find(frame=>frame.contentWindow===event.source);
 if(!source)return;
 const group=event.data.group;
 if(group!==null&&(typeof group!=='string'||group.length>200))return;
 const cohort=source.dataset.groupCohort;
 for(const frame of document.querySelectorAll('iframe[data-group-cohort]'))
  if(frame!==source&&frame.dataset.groupCohort===cohort)frame.contentWindow.postMessage({type:'chronoclade-group-highlight',group},'*');
 });})();
</script>'''
