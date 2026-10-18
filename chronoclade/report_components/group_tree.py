"""Offline, full-topology NJ viewer for direct-distance genetic groups."""
from __future__ import annotations

import json
from html import escape
from pathlib import Path

from chronoclade.report_components.group_controls import (
    INPUT_ORIGINS, group_colors, group_label_map, group_name, group_selector_config,
    group_selector_script, group_selector_styles,
)
from chronoclade.report_components.styles import report_styles
from chronoclade.sample_labels import display_label
from chronoclade.report_components.profile_summary import _reported


def _safe_json(value):
    return (json.dumps(value, ensure_ascii=False, separators=(",", ":"))
            .replace("<", "\\u003c").replace(">", "\\u003e")
            .replace("&", "\\u0026"))


def grouping_caption(group_ids, basis=None, threshold=None):
    """Describe the supplied grouping provenance, never infer clustering from labels."""
    if basis == "unassigned_lineage" or set(group_ids) <= {"Lineage unassigned", "Ungrouped"}:
        return "NJ tree", "Lineage groups are unavailable or not used here. This display category does not infer genetic clusters."
    if basis == "reported_cglin" or all(group.startswith("cgLIN ") for group in group_ids):
        return "cgLIN groups on NJ tree", "Colours show reported cgLIN prefixes; unassigned tips have no prefix grouping. No additional distance cutoff defines these colours."
    if basis == "reported_hiercc" or all(group.startswith("HierCC ") for group in group_ids):
        return "HierCC groups on NJ tree", "Colours show reported HierCC assignments; unassigned tips have no assignment. No additional distance cutoff defines these colours."
    if basis == "distance_threshold":
        detail = f" (mismatch fraction {escape(str(threshold))})" if threshold is not None else ""
        return "Genetic groups on NJ tree", "Groups use the declared direct-profile distance threshold" + detail + ", not tree clades."
    return "Groups on NJ tree", "Colours identify supplied group labels; no grouping method is declared in this view."


def write_group_tree(tree, records, labels, groups, output_path, *, distance_units=None,
                     group_basis=None, group_threshold=None) -> str:
    """Write every branch and tip with an independently navigable viewport.

    The supplied topology and signed branch lengths are never modified. Only the
    viewport changes; labels and markers retain their screen size when zooming.
    """
    records = list(records)
    labels, groups = labels or {}, groups or {}
    by_id = {str(row.get("sample_id")): row for row in records if row.get("sample_id") is not None}
    if len(by_id) != len(records):
        raise ValueError("Tree records require unique sample_id values")
    tips = tree.get_terminals()
    tip_ids = [tip.name for tip in tips]
    if any(name is None for name in tip_ids) or len(set(tip_ids)) != len(tip_ids):
        raise ValueError("Group tree requires uniquely named tips")
    if set(tip_ids) != set(by_id):
        raise ValueError("Group tree records must match every tree tip")
    nodes, depths, parent = [], {}, {}
    stack = [(tree.root, 0.0, None)]
    while stack:
        node, depth, ancestor = stack.pop()
        nodes.append(node)
        depths[node] = depth
        if ancestor is not None:
            parent[node] = ancestor
        for child in reversed(node.clades):
            length = float(child.branch_length or 0)
            if not (float("-inf") < length < float("inf")):
                raise ValueError("Group tree requires finite branch lengths")
            stack.append((child, depth + length, node))
    tip_order = [node for node in nodes if not node.clades]
    positions = {node: index for index, node in enumerate(tip_order)}
    for node in reversed(nodes):
        if node.clades:
            child_y = [positions[child] for child in node.clades]
            positions[node] = (min(child_y) + max(child_y)) / 2
    minimum, maximum = min(depths.values()), max(depths.values())
    span = maximum - minimum
    width, row_space, top = 1200, 20, 28
    height = top * 2 + max(1, len(tip_order) - 1) * row_space

    def x_for(node):
        return 48 + ((depths[node] - minimum) / span * 900 if span else 45)

    optional_fields = [field for field in ("region", "city", "host")
                       if any(_reported(row.get(field)) for row in records)]
    metadata = []
    for node in tip_order:
        name, row = node.name, by_id[node.name]
        raw_group = str(groups.get(name) or "Ungrouped")
        metadata.append({
            "id": name, "label": display_label(str(labels.get(name) or name)),
            "group": group_name(raw_group), "group_id": raw_group,
            "country": str(row.get("country") or "Unknown country"),
            "date": str(row.get("collection_date") or row.get("date") or row.get("collection_year") or "Unknown date"),
            "origin": str(row.get("origin") or row.get("role") or "context"),
            "x": x_for(node), "y": top + positions[node] * row_space,
            **{field: str(row[field]).strip() if _reported(row.get(field))
               else f"Unknown {field}" for field in optional_fields},
        })
    group_ids = sorted({item["group_id"] for item in metadata})
    group_labels = group_label_map(group_ids)
    for item in metadata:
        item["group"] = group_labels[item["group_id"]]
    colors = group_colors(group_ids, labels=group_labels)
    config = group_selector_config(metadata, group_key="group_id", labels=group_labels)
    config["labels"] = group_labels
    data = {"metadata_fields": optional_fields, "tips": metadata, "groups": group_ids, "colors": colors,
            "displayColors": {group_labels[key]: value for key, value in colors.items()},
            "selector": config, "width": width, "height": height, "rowSpace": row_space,
            "units": str(distance_units or "")}
    svg = [f'<svg id="tree" role="img" aria-label="Complete NJ group tree" viewBox="0 0 {width} 400" width="{width}" height="400"><g id="branches">']
    for node in nodes:
        if node.clades:
            low = top + min(positions[child] for child in node.clades) * row_space
            high = top + max(positions[child] for child in node.clades) * row_space
            svg.append(f'<line class="branch" x1="{x_for(node):.6f}" y1="{low:.4f}" x2="{x_for(node):.6f}" y2="{high:.4f}"/>')
        if node in parent:
            ancestor = parent[node]
            y = top + positions[node] * row_space
            css = "branch negative" if (node.branch_length or 0) < 0 else "branch"
            svg.append(f'<line class="{css}" x1="{x_for(ancestor):.6f}" y1="{y:.4f}" x2="{x_for(node):.6f}" y2="{y:.4f}" data-length="{float(node.branch_length or 0):.12g}"/>')
    svg.append('</g><g id="tip-layer">')
    for index, item in enumerate(metadata):
        shape = "input" if item["origin"] in INPUT_ORIGINS else "context"
        color = colors[item["group_id"]]
        marker = (f'<path class="tip-marker input" d="M -5 -5 L 5 0 L -5 5 Z" fill="{color}"/>'
                  if shape == "input" else f'<circle class="tip-marker context" r="3" fill="{color}"/>')
        svg.append(f'<g class="tip" tabindex="-1" role="button" data-tip="{index}" aria-label="{escape(item["label"], quote=True)}"><rect class="tip-band" x="0" y="-9" width="240" height="18" fill="{color}"/>{marker}<text x="9" y="4">{escape(item["label"])}</text></g>')
    svg.append('</g></svg>')
    title, note = grouping_caption(group_ids, group_basis, group_threshold)
    unit = f' Branch lengths are in {escape(str(distance_units))}.' if distance_units else ''
    document = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
                '<meta name="viewport" content="width=device-width,initial-scale=1">'
                f'<title>{title}</title><style>' + report_styles() + group_selector_styles() + TREE_CSS +
                f'</style></head><body><h2>{title}</h2>' + TREE_CONTROLS +
                '<div id="legend" aria-label="Genetic group counts"></div>'
                '<div class="tree-canvas" id="canvas" tabindex="0" aria-label="Tree viewport. Drag to pan; scroll to move; control or command scroll to zoom.">' + ''.join(svg) + '</div>'
                '<p id="tip-details" role="status" aria-live="polite">Hover, focus, or select a tip to inspect its metadata.</p>'
                '<p class="internal-note">' + note + unit + ' Triangles mark input genomes; circles mark comparisons. '
                'Drag to pan; scroll to move vertically, Shift-scroll horizontally; Ctrl/Command-scroll to zoom. '
                'Labels appear at readable spacing. Every branch and tip remains available.</p>'
                '<script type="application/json" id="tree-data">' + _safe_json(data) + '</script><script>' +
                group_selector_script() + TREE_SCRIPT + '</script></body></html>')
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(document, encoding="utf-8")
    return str(output_path)


TREE_CSS = '''
body{background:#fff;padding:8px;margin:0}h2{margin:0 0 6px;font-size:20px}.tree-controls{display:flex;gap:5px 8px;align-items:center;flex-wrap:wrap}
button,input{font:inherit}button{padding:5px 8px;background:#fff;border:1px solid #858585;border-radius:3px;cursor:pointer}button:disabled{color:#666;cursor:default;border-color:#ccc}
.tree-controls label{font-size:12px}.tree-controls input[type=search]{width:180px;max-width:100%;padding:5px;border:1px solid #858585}
#view-status{font-size:12px;color:#454545}.tree-canvas{height:400px;min-height:400px;overflow:hidden;position:relative;border:1px solid #aaa;background:#fff;touch-action:none;cursor:grab}
.tree-canvas.dragging{cursor:grabbing}.tree-canvas svg{display:block;width:100%;height:100%;max-width:none;font:12px 'Archivo',sans-serif;user-select:none}
.branch{stroke:#68707c;stroke-width:.7;fill:none;vector-effect:non-scaling-stroke}.branch.negative{stroke:#888;stroke-dasharray:4 3}
.tip-band{opacity:0}.tip.selected .tip-band{opacity:.16}.tip text{fill:#17191f;paint-order:stroke;stroke:white;stroke-width:2px;stroke-linejoin:round}.tip{cursor:pointer}.tip.dimmed{opacity:.18}.tip-marker{stroke:#222;stroke-width:.65}.tip:focus{outline:none}.tip:focus .tip-marker,.tip.located .tip-marker{stroke:#000;stroke-width:2}
#tip-details{font-size:12px;line-height:1.4;min-height:2.8em;max-height:4.2em;overflow:auto;padding:6px;margin:6px 0;border:1px solid #aaa;white-space:pre-wrap}.internal-note{color:#454545;font-size:11px;line-height:1.35;margin:6px 0}
::selection{background:#17191f;color:white}
'''
TREE_CONTROLS = '''<div class="tree-controls">
<button id="fit" type="button">Fit all</button><button id="zoom-in" type="button" aria-label="Zoom in">Zoom +</button><button id="zoom-out" type="button" aria-label="Zoom out">Zoom −</button>
<button id="read" type="button">Read labels</button><button id="zoom-inputs" type="button">Zoom to inputs</button><button id="zoom-group" type="button" disabled>Zoom to selected group</button><button id="clear" type="button">Clear highlight</button>
<label>Find sample <input id="tip-search" type="search" placeholder="Accession or sample name" list="tip-matches"></label><datalist id="tip-matches"></datalist><button id="locate-tip" type="button">Locate</button><span id="view-status" role="status"></span></div>'''
TREE_SCRIPT = '''
(()=>{
const data=JSON.parse(document.getElementById('tree-data').textContent),svg=document.getElementById('tree'),canvas=document.getElementById('canvas'),branches=document.getElementById('branches');
const details=document.getElementById('tip-details'),tips=[...svg.querySelectorAll('.tip')],status=document.getElementById('view-status'),groupZoom=document.getElementById('zoom-group');
let selected=null,located=null,mode='overview',view={x:0,y:0,w:data.width,h:data.height},frame=null;
const selector=createGroupSelector({container:document.getElementById('legend'),groups:data.groups,colors:data.colors,config:data.selector,onSelect:group=>highlight(group,true)});
const inputIndexes=data.tips.map((t,i)=>['input','local','query','focal'].includes(t.origin)?i:-1).filter(i=>i>=0);
document.getElementById('zoom-inputs').disabled=!inputIndexes.length;
function showTip(index){const t=data.tips[index];details.textContent=[t.label,`Group: ${t.group}`,`Country: ${t.country}`,...data.metadata_fields.map(f=>`${f[0].toUpperCase()+f.slice(1)}: ${t[f]}`),`Date: ${t.date}`].join(' | ')}
function send(group){if(window.parent!==window)window.parent.postMessage({type:'chronoclade-group-highlight',group},'*')}
function highlight(group,notify){selected=group;selector.setSelected(group);groupZoom.disabled=!group;tips.forEach((el,i)=>{const active=!!group&&data.tips[i].group_id===group;el.classList.toggle('dimmed',!!group&&!active);el.classList.toggle('selected',active)});if(notify)send(group)}
function bounds(){view.w=Math.min(data.width,Math.max(80,view.w));view.h=Math.min(data.height,Math.max(100,view.h));view.x=Math.max(0,Math.min(data.width-view.w,view.x));view.y=Math.max(0,Math.min(data.height-view.h,view.y))}
function render(){frame=null;bounds();const w=canvas.clientWidth,h=canvas.clientHeight,sx=w/view.w,sy=h/view.h,readable=data.rowSpace*sy>=14;svg.setAttribute('viewBox',`0 0 ${w} ${h}`);branches.setAttribute('transform',`matrix(${sx} 0 0 ${sy} ${-view.x*sx} ${-view.y*sy})`);
tips.forEach((el,i)=>{const t=data.tips[i],x=(t.x-view.x)*sx,y=(t.y-view.y)*sy,visible=y>=-15&&y<=h+15&&x>=-260&&x<=w+15;el.style.display=visible?'':'none';el.setAttribute('tabindex',visible&&readable?'0':'-1');if(!visible)return;el.setAttribute('transform',`translate(${x} ${y})`);el.querySelector('text').style.display=readable?'':'none';el.querySelector('.tip-band').style.display=readable?'':'none';const marker=el.querySelector('.tip-marker');if(marker.tagName.toLowerCase()==='circle')marker.setAttribute('r',readable?'3':'1.1');el.classList.toggle('located',i===located)});
status.textContent=mode==='overview'?`Full tree · ${data.tips.length.toLocaleString()} tips`:`${readable?'Read labels':'Zoomed tree'} · ${Math.max(1,Math.round(data.tips.length*view.h/data.height)).toLocaleString()} tips in vertical range`;}
function schedule(){if(frame===null)frame=requestAnimationFrame(render)}
function fitAll(){mode='overview';view={x:0,y:0,w:data.width,h:data.height};schedule()}
function zoom(factor,ax=.5,ay=.5,horizontal=false){mode='focused';const oldW=view.w,oldH=view.h;if(horizontal)view.w=Math.min(data.width,Math.max(80,oldW/factor));view.h=Math.min(data.height,Math.max(100,oldH/factor));view.x+=(oldW-view.w)*ax;view.y+=(oldH-view.h)*ay;schedule()}
function focusTips(indexes){if(!indexes.length)return;mode='focused';const xs=indexes.map(i=>data.tips[i].x),ys=indexes.map(i=>data.tips[i].y),low=Math.min(...ys),high=Math.max(...ys);const height=Math.max(100,(high-low)*1.2+60);view={x:Math.max(0,Math.min(...xs)-120),y:(low+high-height)/2,w:Math.max(400,Math.max(...xs)-Math.min(...xs)+380),h:height};schedule()}
function readLabels(){mode='labels';const center=view.y+view.h/2;view.h=Math.min(data.height,canvas.clientHeight/24*data.rowSpace);view.y=center-view.h/2;schedule()}
function locateTip(){const query=document.getElementById('tip-search').value.trim().toLocaleLowerCase();if(!query)return;const exact=data.tips.findIndex(t=>t.label.toLocaleLowerCase()===query||t.id.toLocaleLowerCase()===query);const index=exact>=0?exact:data.tips.findIndex(t=>`${t.label} ${t.id}`.toLocaleLowerCase().includes(query));if(index<0){details.textContent='No matching sample. Search an accession or sample name.';return}located=index;focusTips([index]);showTip(index)}
document.getElementById('fit').addEventListener('click',fitAll);document.getElementById('zoom-in').addEventListener('click',()=>zoom(2));document.getElementById('zoom-out').addEventListener('click',()=>zoom(.5));document.getElementById('read').addEventListener('click',readLabels);document.getElementById('zoom-inputs').addEventListener('click',()=>focusTips(inputIndexes));groupZoom.addEventListener('click',()=>focusTips(data.tips.map((t,i)=>t.group_id===selected?i:-1).filter(i=>i>=0)));document.getElementById('clear').addEventListener('click',()=>highlight(null,true));
const search=document.getElementById('tip-search'),matches=document.getElementById('tip-matches');search.addEventListener('input',()=>{const query=search.value.trim().toLocaleLowerCase();matches.replaceChildren();if(query.length<2)return;data.tips.filter(t=>`${t.label} ${t.id}`.toLocaleLowerCase().includes(query)).slice(0,20).forEach(t=>{const option=document.createElement('option');option.value=t.label;matches.append(option)})});search.addEventListener('keydown',event=>{if(event.key==='Enter')locateTip()});document.getElementById('locate-tip').addEventListener('click',locateTip);
svg.addEventListener('pointerover',event=>{const tip=event.target.closest('.tip');if(tip)showTip(Number(tip.dataset.tip))});svg.addEventListener('focusin',event=>{const tip=event.target.closest('.tip');if(tip)showTip(Number(tip.dataset.tip))});svg.addEventListener('click',event=>{const tip=event.target.closest('.tip');if(tip){located=Number(tip.dataset.tip);showTip(located);schedule()}});
let drag=null;canvas.addEventListener('pointerdown',event=>{if(event.button!==0)return;drag={x:event.clientX,y:event.clientY,vx:view.x,vy:view.y};canvas.setPointerCapture(event.pointerId);canvas.classList.add('dragging')});canvas.addEventListener('pointermove',event=>{if(!drag)return;mode='focused';view.x=drag.vx-(event.clientX-drag.x)*view.w/canvas.clientWidth;view.y=drag.vy-(event.clientY-drag.y)*view.h/canvas.clientHeight;schedule()});function stopDrag(){drag=null;canvas.classList.remove('dragging')}canvas.addEventListener('pointerup',stopDrag);canvas.addEventListener('pointercancel',stopDrag);
canvas.addEventListener('wheel',event=>{event.preventDefault();if(event.ctrlKey||event.metaKey){const rect=canvas.getBoundingClientRect();zoom(Math.exp(-event.deltaY*.005),(event.clientX-rect.left)/rect.width,(event.clientY-rect.top)/rect.height,true)}else{mode='focused';view.x+=(event.shiftKey?event.deltaY:event.deltaX)*view.w/canvas.clientWidth;if(!event.shiftKey)view.y+=event.deltaY*view.h/canvas.clientHeight;schedule()}},{passive:false});
canvas.addEventListener('keydown',event=>{const directions={ArrowDown:[0,60],ArrowUp:[0,-60],ArrowLeft:[-60,0],ArrowRight:[60,0]};if(directions[event.key]){event.preventDefault();mode='focused';view.x+=directions[event.key][0]*view.w/canvas.clientWidth;view.y+=directions[event.key][1]*view.h/canvas.clientHeight;schedule()}});
window.addEventListener('message',event=>{if(event.source!==window.parent||!event.data||event.data.type!=='chronoclade-group-highlight')return;const group=event.data.group===null?null:String(event.data.group);const match=data.groups.includes(group)?group:data.groups.find(key=>data.selector.labels[key]===group);highlight(match||null,false)});
new ResizeObserver(schedule).observe(canvas);fitAll();
})();
'''

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
