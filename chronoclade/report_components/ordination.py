"""Shared offline ordination views; dataset markers stay distinct in every view."""
from __future__ import annotations

import json
import re
from collections import Counter
from html import escape
from pathlib import Path

import numpy as np

from chronoclade.metadata_dates import sample_date_interval
from chronoclade.sample_labels import sample_labels

_INPUT_ORIGINS = {"input", "local", "query", "focal"}


def write_ordination_views(output, prefix, coordinates, records, *,
                           title="Distance ordination", axis_labels=("Axis 1", "Axis 2"),
                           groups=None, coordinate_views=None, default_basis=None):
    """Write consistent dataset/country/date/group views using identical coordinates.

    Collection dates are midpoint decimal years of their recorded intervals.
    Missing dates remain visible as grey points, never zero-valued dates.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    from chronoclade.report_components.styles import report_styles

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if not re.fullmatch(r"[A-Za-z0-9_-]+", prefix):
        raise ValueError("Ordination prefix must be a safe artifact name")
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.shape != (len(records), 2) or not np.isfinite(coordinates).all():
        raise ValueError("Ordination requires two finite coordinates per record")
    ids = [row["sample_id"] for row in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Ordination sample IDs must be unique")
    bases = {} if coordinate_views else {
        "ordination": {"coordinates": coordinates, "label": "Ordination",
                       "axis_labels": axis_labels, "note": ""}
    }
    for key, view in (coordinate_views or {}).items():
        if not re.fullmatch(r"[A-Za-z0-9_-]+", str(key)):
            raise ValueError("Coordinate basis names must be safe identifiers")
        basis_coordinates = np.asarray(view["coordinates"], dtype=float)
        if basis_coordinates.shape != (len(records), 2) or not np.isfinite(basis_coordinates).all():
            raise ValueError("Each coordinate basis requires two finite coordinates per record")
        basis_axes = tuple(view.get("axis_labels", ("Axis 1", "Axis 2")))
        if len(basis_axes) != 2:
            raise ValueError("Coordinate basis axis_labels must contain two labels")
        bases[str(key)] = {"coordinates": basis_coordinates, "label": str(view.get("label", key)),
                           "axis_labels": basis_axes, "note": str(view.get("note", ""))}
    if default_basis is None:
        default_basis = next(iter(coordinate_views or {}), "ordination")
    if default_basis not in bases:
        raise ValueError(f"Unknown default coordinate basis: {default_basis}")
    coordinates = bases[default_basis]["coordinates"]
    axis_labels = bases[default_basis]["axis_labels"]
    roles = ["input" if row.get("role", row.get("origin")) in _INPUT_ORIGINS else "context"
             for row in records]
    markers = {"input": "^", "context": "o"}
    names = {"input": "Input genomes", "context": "Public comparisons"}
    dates = [sample_date_interval(row) for row in records]
    years = np.array([value["midpoint_year"]
                      if value["status"] == "valid" else np.nan for value in dates])
    categories = {
        "dataset": [names[role] for role in roles],
        "country": [str(row.get("country") or "Unknown country") for row in records],
        "region": [str(row.get("region") or "Unknown region") for row in records],
        "host": [str(row.get("host") or "Unknown host") for row in records],
    }
    if groups is not None:
        def group_label(ident):
            value = str(groups.get(ident, "Ungrouped"))
            match = re.search(r"(?:^|_)group_(\d+)$", value, re.IGNORECASE)
            return f"Group {match.group(1)}" if match else value
        categories["group"] = [group_label(ident) for ident in ids]
    paths, views = {}, []
    for view in [*categories, "date"]:
        fig, ax = plt.subplots(figsize=(9, 6))
        if view == "date":
            finite = np.isfinite(years)
            lower, upper = (float(years[finite].min()), float(years[finite].max())) if finite.any() else (0, 1)
            norm = plt.Normalize(lower, upper if upper > lower else lower + 1)
            for role, marker in markers.items():
                zorder = 3 if role == "input" else 2
                selected = np.array([r == role for r in roles])
                known, missing = selected & finite, selected & ~finite
                if known.any():
                    scatter = ax.scatter(*coordinates[known].T, c=years[known], cmap="viridis",
                                         norm=norm, marker=marker, alpha=.8, s=35,
                                         edgecolors="#222", linewidths=.3, zorder=zorder)
                if missing.any():
                    ax.scatter(*coordinates[missing].T, c="#999", marker=marker, s=35,
                               edgecolors="#222", linewidths=.3, zorder=zorder)
            if finite.any():
                fig.colorbar(scatter, ax=ax, label="Collection year (interval midpoint)")
            handles = [Line2D([], [], color="#999", marker="o", linestyle="None", label="Date unknown")]
        else:
            values = categories[view]
            unique = sorted(set(values))
            palette = plt.get_cmap("tab20", max(1, len(unique)))
            colours = {value: palette(i) for i, value in enumerate(unique)}
            if view == "dataset":
                colours = {"Input genomes": "#2166ac", "Public comparisons": "#666666"}
            counts = Counter(values)
            shown = sorted(unique, key=lambda value: (-counts[value], value))[:20]
            handles = [Line2D([], [], color=colours[value], marker="o", linestyle="None", label=value)
                       for value in shown]
            # One scatter per role avoids quadratic scans and thousands of artists.
            for role, marker in markers.items():
                zorder = 3 if role == "input" else 2
                selected = [i for i, observed in enumerate(roles) if observed == role]
                if selected:
                    ax.scatter(*coordinates[selected].T,
                               c=[colours[values[i]] for i in selected], marker=marker,
                               alpha=.8, s=35, edgecolors="#222", linewidths=.3, zorder=zorder)
            if len(unique) > len(shown):
                ax.text(0, -.15,
                        f"Legend shows {len(shown)} of {len(unique)} categories by frequency; all genomes are plotted. "
                        "Full labels are in the downloaded tables.",
                        transform=ax.transAxes, fontsize=8, wrap=True)
        handles += [Line2D([], [], color="#222", marker=marker, linestyle="None", label=names[role])
                    for role, marker in markers.items() if role in roles]
        if view == "date":
            ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.01),
                      ncol=3, fontsize=8, frameon=False)
        else:
            ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1),
                      fontsize=8, frameon=False)
        ax.set(xlabel=axis_labels[0], ylabel=axis_labels[1], title=title)
        if view == "date":
            ax.set_title(title, pad=30)
        fig.tight_layout()
        filename = f"{prefix}_pcoa" + ("" if view == "dataset" else f"_{view}") + ".svg"
        fig.savefig(output / filename, bbox_inches="tight")
        plt.close(fig)
        key = "pcoa_figure" if view == "dataset" else f"pcoa_{view}_figure"
        paths[key] = str(output / filename)
        views.append((view, filename))
    labels = sample_labels(records)
    rows = []
    for i, row in enumerate(records):
        item = {
            "id": ids[i], "label": labels.get(ids[i], ids[i]), "dataset": names[roles[i]],
            "role": roles[i], "country": categories["country"][i],
            "region": categories["region"][i], "host": categories["host"][i],
            "date": (f'{dates[i]["start"]} to {dates[i]["end"]}'
                     if dates[i].get("status") == "valid" else "Unknown date"),
            "year": (float(years[i]) if np.isfinite(years[i]) else None),
            "group": (categories["group"][i] if "group" in categories else "Ungrouped"),
            "group_id": str(groups.get(ids[i], "Ungrouped")) if groups is not None else "Ungrouped",
        }
        item["coordinates"] = {key: value["coordinates"][i].tolist() for key, value in bases.items()}
        rows.append(item)
    viewer_data = {
        "title": title, "bases": {key: {"label": value["label"],
                                           "axis_labels": value["axis_labels"],
                                           "note": value["note"]} for key, value in bases.items()},
        "default_basis": default_basis,
        "default_view": "group" if groups is not None else "dataset",
        "views": list(categories) + ["date"], "rows": rows,
    }
    basis_options = "".join(
        f'<option value="{escape(key, quote=True)}"{" selected" if key == default_basis else ""}>'
        f'{escape(value["label"])}</option>' for key, value in bases.items())
    group_name = 'cgLIN group' if groups and all(str(value).startswith('cgLIN ') for value in groups.values()) else 'Group'
    view_options = "".join(f'<option value="{escape(view, quote=True)}">{escape(group_name if view == "group" else view.title())}</option>'
                           for view in viewer_data["views"])
    view_options = view_options.replace(
        f'value="{viewer_data["default_view"]}"',
        f'value="{viewer_data["default_view"]}" selected',
    )
    payload = json.dumps(viewer_data, ensure_ascii=False, separators=(",", ":"))
    payload = payload.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    document = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
                '<meta name="viewport" content="width=device-width,initial-scale=1">'
                f'<title>{escape(title)}</title><style>{report_styles()} '
                'body{background:white;padding:8px} .ordination-controls{display:flex;gap:14px;flex-wrap:wrap;align-items:center}'
                'fieldset{border:0;border-bottom:1px solid #ccc;padding:8px 0}label{font-weight:700}'
                'select,button{font:inherit;padding:5px}.plot-wrap{position:relative;width:100%;overflow:hidden}'
                'svg{display:block;width:100%;height:min(52vh,400px);min-height:320px;background:#fff}'
                '.legend{display:flex;gap:6px 12px;flex-wrap:wrap;max-height:150px;overflow-y:auto;margin:8px 0}.legend button{border:1px solid #bbb;background:white;cursor:pointer}'
                '.legend button[aria-pressed="true"]{outline:2px solid #17191f}.swatch{display:inline-block;width:12px;height:12px;margin-right:5px;vertical-align:-1px}'
                '#details{min-height:3em;max-height:180px;overflow-y:auto;padding:8px;border:1px solid #ccd0d8;white-space:pre-wrap}.point:focus{outline:none;stroke:#ffb800;stroke-width:3px}'
                '</style></head><body>'
                '<fieldset><legend>Ordination view</legend><div class="ordination-controls">'
                f'<label>Coordinates <select id="basis">{basis_options}</select></label>'
                f'<label>Colour by <select id="colour">{view_options}</select></label>'
                '<button id="clear" type="button">Clear highlight</button></div></fieldset>'
                '<div id="legend" class="legend" aria-label="Clickable category legend"></div>'
                '<div class="plot-wrap"><svg id="plot" role="img" aria-label="Interactive ordination scatter plot" tabindex="0"></svg></div>'
                '<div id="details" role="status" aria-live="polite">Focus, hover over, or select a point to inspect its metadata.</div>'
                '<p id="note" class="muted"></p><p class="muted">Triangles mark input genomes; circles mark public comparisons. '
                'Point colour reflects the selected field. Dates use collection-interval midpoints.</p>'
                f'<script type="application/json" id="ordination-data">{payload}</script>'
                '<script>' + _ORDINATION_SCRIPT + '</script></body></html>')
    viewer = output / f"{prefix}_ordination.html"
    viewer.write_text(document)
    paths["pcoa_views_html"] = str(viewer)
    return paths


_ORDINATION_SCRIPT = r"""
(()=>{
const data=JSON.parse(document.getElementById('ordination-data').textContent);
const svg=document.getElementById('plot'), details=document.getElementById('details');
const basis=document.getElementById('basis'), colour=document.getElementById('colour');
const legend=document.getElementById('legend'), note=document.getElementById('note');
const ns='http://www.w3.org/2000/svg', palette=['#2166ac','#b2182b','#1b7837','#762a83','#e08214','#008080','#c51b7d','#4d4d4d','#a6d854','#3288bd','#f46d43','#66c2a5'];
let points=[], highlight=null, cursor=0;
const svgEl=(tag,attrs={})=>{const e=document.createElementNS(ns,tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,String(v));return e};
const fieldValue=(r,f)=>f==='dataset'?r.dataset:f==='date'?(r.year===null?'Unknown date':String(Math.floor(r.year))):r[f];
const dateColour=(year,lo,hi)=>{const t=hi===lo?0.5:Math.max(0,Math.min(1,(year-lo)/(hi-lo))),a=[68,1,84],b=[253,231,37];return '#'+a.map((x,i)=>Math.round(x+(b[i]-x)*t).toString(16).padStart(2,'0')).join('')};
function colourLookup(field){
 if(field==='date'){const years=data.rows.filter(r=>r.year!==null).map(r=>r.year),lo=years.length?Math.min(...years):0,hi=years.length?Math.max(...years):1;return {colours:data.rows.map(r=>r.year===null?'#999999':dateColour(r.year,lo,hi)),lo,hi}}
 const values=[...new Set(data.rows.map(r=>fieldValue(r,field)))].sort(),byValue=new Map();
 values.forEach((value,i)=>byValue.set(value,field==='dataset'?(value==='Input genomes'?'#2166ac':'#666666'):palette[i%palette.length]));
 return {colours:data.rows.map(r=>byValue.get(fieldValue(r,field))),byValue};
}
const sendGroup=()=>{if(window.parent!==window)window.parent.postMessage({type:'chronoclade-group-highlight',group:highlight?.field==='group'?highlight.value:null},'*')};
const resetDetails=()=>{details.textContent='Focus, hover over, or select a point to inspect its metadata.'};
const showRows=(rows)=>{details.textContent=`Samples at/near this point (${rows.length})\n\n`+rows.map(r=>`Accession: ${r.label}\nDataset: ${r.dataset}\nCountry: ${r.country}\nHost: ${r.host}\nCollection date: ${r.date}\nGroup: ${r.group}\nColour category: ${fieldValue(r,colour.value)}`).join('\n\n')};
function renderLegend(){legend.replaceChildren();const f=colour.value;
 if(f==='date'){const years=data.rows.filter(r=>r.year!==null).map(r=>r.year),gradient=document.createElement('span');gradient.className='swatch';gradient.style.cssText='width:110px;background:linear-gradient(90deg,#440154,#fde725)';legend.append(gradient,document.createTextNode(`${years.length?Math.min(...years).toFixed(2):'No dates'} – ${years.length?Math.max(...years).toFixed(2):''}`));const missing=data.rows.filter(r=>r.year===null).length;if(missing){const b=document.createElement('button');b.type='button';b.setAttribute('aria-pressed',String(highlight?.field==='date'&&highlight.value==='Unknown date'));const sw=document.createElement('span');sw.className='swatch';sw.style.backgroundColor='#999';b.append(sw,document.createTextNode(`Unknown date (${missing})`));b.addEventListener('click',()=>{highlight=highlight?.field==='date'&&highlight.value==='Unknown date'?null:{field:'date',value:'Unknown date'};resetDetails();renderLegend();draw();if(f==='group')sendGroup()});legend.append(b)}return}
 const counts=new Map();for(const r of data.rows){const v=fieldValue(r,f);counts.set(v,(counts.get(v)||0)+1)}const lookup=colourLookup(f);
 [...counts].sort((a,b)=>b[1]-a[1]||String(a[0]).localeCompare(String(b[0]))).forEach(([v,n])=>{const b=document.createElement('button');b.type='button';b.setAttribute('aria-pressed',String(highlight?.field===f&&highlight.value===v));const sw=document.createElement('span');sw.className='swatch';sw.style.backgroundColor=lookup.byValue.get(v);b.append(sw,document.createTextNode(`${v} (${n})`));b.addEventListener('click',()=>{highlight=highlight?.field===f&&highlight.value===v?null:{field:f,value:v};resetDetails();renderLegend();draw();if(f==='group')sendGroup()});legend.append(b)});
}
function draw(){const rect=svg.getBoundingClientRect(),w=Math.max(500,rect.width||900),h=Math.max(320,rect.height||400),pad={l:74,r:24,t:20,b:54},key=basis.value,axes=data.bases[key].axis_labels;
 svg.replaceChildren();svg.setAttribute('viewBox',`0 0 ${w} ${h}`);const coords=data.rows.map(r=>r.coordinates[key]),xs=coords.map(p=>p[0]),ys=coords.map(p=>p[1]);let xmin=xs.length?Math.min(...xs):0,xmax=xs.length?Math.max(...xs):1,ymin=ys.length?Math.min(...ys):0,ymax=ys.length?Math.max(...ys):1;if(xmin===xmax){xmin-=1;xmax+=1}if(ymin===ymax){ymin-=1;ymax+=1}const xpad=(xmax-xmin)*.06,ypad=(ymax-ymin)*.06;xmin-=xpad;xmax+=xpad;ymin-=ypad;ymax+=ypad;const X=x=>pad.l+(x-xmin)/(xmax-xmin)*(w-pad.l-pad.r),Y=y=>h-pad.b-(y-ymin)/(ymax-ymin)*(h-pad.t-pad.b);
 const xticks=Array.from({length:5},(_,i)=>xmin+(xmax-xmin)*i/4),yticks=Array.from({length:5},(_,i)=>ymin+(ymax-ymin)*i/4),fmt=n=>Number(n.toPrecision(4)).toString();
 for(const tick of xticks){const x=X(tick);svg.append(svgEl('path',{d:`M${x},${pad.t}V${h-pad.b}`,stroke:'#e4e6eb'}));const text=svgEl('text',{x,y:h-pad.b+18,'text-anchor':'middle',fill:'#555','font-size':11});text.textContent=fmt(tick);svg.append(text)}
 for(const tick of yticks){const y=Y(tick);svg.append(svgEl('path',{d:`M${pad.l},${y}H${w-pad.r}`,stroke:'#e4e6eb'}));const text=svgEl('text',{x:pad.l-8,y:y+4,'text-anchor':'end',fill:'#555','font-size':11});text.textContent=fmt(tick);svg.append(text)}
 const axis=svgEl('path',{d:`M${pad.l},${pad.t}V${h-pad.b}H${w-pad.r}`,fill:'none',stroke:'#777'});svg.append(axis);const xlabel=svgEl('text',{x:(pad.l+w-pad.r)/2,y:h-8,'text-anchor':'middle',fill:'#222'});xlabel.textContent=axes[0];svg.append(xlabel);const ylabel=svgEl('text',{x:17,y:h/2,transform:`rotate(-90 17 ${h/2})`,'text-anchor':'middle',fill:'#222'});ylabel.textContent=axes[1];svg.append(ylabel);
 const lookup=colourLookup(colour.value);points=[];const pointOrder=data.rows.map((r,i)=>i).sort((a,b)=>(data.rows[a].role==='input'?1:0)-(data.rows[b].role==='input'?1:0));pointOrder.forEach(i=>{const r=data.rows[i],x=X(coords[i][0]),y=Y(coords[i][1]),v=fieldValue(r,colour.value),c=lookup.colours[i],active=!highlight||(highlight.field===colour.value&&(colour.value==='date'?((highlight.value==='Unknown date')===(r.year===null)):highlight.value===v)),p=r.role==='input'?svgEl('path',{d:`M${x},${y-7} L${x+7},${y+6} L${x-7},${y+6} Z`}):svgEl('circle',{cx:x,cy:y,r:5});p.classList.add('point');p.dataset.index=i;p.setAttribute('fill',c);p.setAttribute('opacity',active?'1':'.14');p.setAttribute('stroke',active?'#222':'#ddd');p.setAttribute('stroke-width',highlight&&active?'2':'.7');p.setAttribute('tabindex','0');p.setAttribute('role','button');p.setAttribute('aria-label',`${r.label}, ${r.dataset}, ${r.country}, ${r.host}, ${r.date}, ${r.group}`);const nearby=()=>showRows(data.rows.filter(q=>Math.hypot(X(q.coordinates[key][0])-x,Y(q.coordinates[key][1])-y)<12));p.addEventListener('focus',nearby);p.addEventListener('click',nearby);p.addEventListener('mouseenter',nearby);svg.append(p);points.push(p)});
 note.textContent=data.bases[key].note||'';
}
colour.addEventListener('change',()=>{highlight=null;resetDetails();renderLegend();draw();sendGroup()});basis.addEventListener('change',()=>{resetDetails();draw()});document.getElementById('clear').addEventListener('click',()=>{highlight=null;resetDetails();renderLegend();draw();sendGroup()});
svg.addEventListener('keydown',e=>{if(e.key==='ArrowRight'||e.key==='ArrowDown'){e.preventDefault();cursor=(cursor+1)%points.length;points[cursor].focus()}else if(e.key==='ArrowLeft'||e.key==='ArrowUp'){e.preventDefault();cursor=(cursor+points.length-1)%points.length;points[cursor].focus()}});
window.addEventListener('message',e=>{if(e.source!==window.parent||e.data?.type!=='chronoclade-group-highlight')return;const group=e.data.group;if(group!==null&&(typeof group!=='string'||!data.rows.some(r=>r.group===group)))return;if(group!==null){colour.value='group';highlight={field:'group',value:group}}else{highlight=null}resetDetails();renderLegend();draw()});
renderLegend();draw();window.addEventListener('resize',draw);
})();
"""


def ordination_viewer(cohorts, directory):
    """Embed only existing contained report assets, with an offline fallback link."""
    directory = Path(directory).resolve()
    frames = []
    for cohort in cohorts or []:
        reference = cohort.get("pcoa_views_html")
        if not reference:
            continue
        path = Path(reference)
        path = (directory / path).resolve() if not path.is_absolute() else path.resolve()
        if not path.is_file() or not path.is_relative_to(directory):
            continue
        relative = escape(path.relative_to(directory).as_posix(), quote=True)
        frames.append(f'<iframe src="{relative}" title="Ordination metadata views" '
                      f'data-group-cohort="{escape(str(cohort.get("cohort_id", "")), quote=True)}" '
                      'style="width:100%;height:760px;border:0" loading="lazy"></iframe>'
                      f'<p><a href="{relative}">Open ordination views</a></p>')
    return "".join(frames)
