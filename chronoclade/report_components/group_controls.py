"""Shared group colours and bounded selectors for standalone report views."""
from __future__ import annotations

import re
from collections import Counter

PALETTE = (
    "#2166ac", "#b2182b", "#1b7837", "#762a83", "#e08214", "#008080",
    "#c51b7d", "#4d4d4d", "#a6d854", "#3288bd", "#f46d43", "#66c2a5",
)
INPUT_ORIGINS = {"input", "local", "query", "focal"}


def group_name(value):
    value = str(value or "Ungrouped")
    match = re.search(r"(?:^|_)group_(\d+)$", value, re.IGNORECASE)
    return f"Group {match.group(1)}" if match else value


def group_label_map(group_ids):
    """Keep familiar short names, distinguishing any aliases that collide."""
    labels = {key: group_name(key) for key in group_ids}
    counts = Counter(labels.values())
    return {key: f"{label} [{key}]" if counts[label] > 1 else label
            for key, label in labels.items()}


def group_colors(group_names, *, labels=None):
    return {name: PALETTE[index % len(PALETTE)]
            for index, name in enumerate(sorted(set(group_names), key=lambda name: ((labels or {}).get(name, name), name)))}


def group_selector_config(rows, *, group_key="group", role_key="origin", labels=None):
    rows = list(rows)
    counts = Counter(row[group_key] for row in rows)
    ordered = sorted(counts, key=lambda name: (-counts[name], (labels or {}).get(name, name), name))
    inputs = {row[group_key] for row in rows if row.get(role_key) in INPUT_ORIGINS}
    compact = set(ordered[:20]) | inputs
    return {"counts": dict(counts), "inputGroups": sorted(inputs),
            "compactGroups": [name for name in ordered if name in compact], "limit": 20}


def group_selector_styles():
    return '''.group-tools{display:flex;gap:6px 12px;flex-wrap:wrap;align-items:center;margin:6px 0}
.group-tools input[type=search]{font:inherit;padding:5px;min-width:160px;max-width:100%;border:1px solid #858585;border-radius:3px}
.group-tools label{font-size:12px}.group-summary{font-size:12px;color:#454545}
.group-options{display:flex;gap:5px 8px;flex-wrap:wrap;max-height:140px;overflow-y:auto;align-content:start;padding:3px;margin:4px 0;scrollbar-color:#666 #eee}
.group-options button{font:inherit;font-size:12px;padding:4px 6px;border:1px solid #aaa;background:#fff;border-radius:3px;cursor:pointer}
.group-options button[aria-pressed=true]{outline:2px solid #17191f;outline-offset:0}
.group-options .swatch{display:inline-block;width:10px;height:10px;margin-right:5px;vertical-align:-1px}
button:focus-visible,input:focus-visible{outline:2px solid #17191f;outline-offset:2px}'''


def group_selector_script():
    return '''function createGroupSelector({container,groups,colors,config,onSelect}){
let selected=null;const ordered=[...groups].sort((a,b)=>(config.counts[b]||0)-(config.counts[a]||0)||(config.labels?.[a]||a).localeCompare(config.labels?.[b]||b)||a.localeCompare(b));
container.replaceChildren();const tools=document.createElement('div');tools.className='group-tools';
const searchLabel=document.createElement('label');searchLabel.textContent='Search all groups ';const search=document.createElement('input');search.type='search';search.placeholder='Group name';search.setAttribute('aria-label','Search all groups');searchLabel.append(search);
const showLabel=document.createElement('label'),all=document.createElement('input');all.type='checkbox';showLabel.append(all,document.createTextNode(' Show all groups'));
const summary=document.createElement('span');summary.className='group-summary';summary.setAttribute('role','status');tools.append(searchLabel,showLabel,summary);
const options=document.createElement('div');options.className='group-options';options.setAttribute('aria-label','Genetic group counts');container.append(tools,options);
function render(){const query=search.value.trim().toLocaleLowerCase();const compact=new Set(config.compactGroups);const visible=ordered.filter(group=>query?`${config.labels?.[group]||group} ${group}`.toLocaleLowerCase().includes(query):all.checked||compact.has(group)||group===selected);options.replaceChildren();
for(const group of visible){const button=document.createElement('button');button.type='button';button.dataset.group=group;button.setAttribute('aria-pressed',String(group===selected));const swatch=document.createElement('span');swatch.className='swatch';swatch.style.backgroundColor=colors[group];button.append(swatch,document.createTextNode(`${config.labels?.[group]||group} (${config.counts[group]||0})${config.inputGroups.includes(group)?' · input':''}`));button.addEventListener('click',()=>{selected=selected===group?null:group;render();onSelect(selected)});options.append(button)}
summary.textContent=`${visible.length} of ${groups.length} groups${query?' match':all.checked?'':' · largest 20 + input groups'}`;if(!visible.length){options.textContent='No matching groups.'}}
search.addEventListener('input',render);all.addEventListener('change',render);render();
return{setSelected(group){selected=groups.includes(group)?group:null;render()},destroy(){container.replaceChildren()}};
}
'''
