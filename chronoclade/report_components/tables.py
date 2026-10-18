"""Shared, offline table controls; source rows and downloaded evidence stay intact."""

from pathlib import Path

_TABLE_CSS = """
<style data-chronoclade-tables>
.table-tools,.table-pages{display:flex;align-items:center;flex-wrap:wrap;gap:12px;margin:14px 0;font-size:13px;color:var(--muted,#5f6470)}
.table-tools label{display:flex;align-items:center;gap:8px}
.table-tools input,.table-tools select,.table-pages button{font:inherit;color:var(--ink,#17191f);background:#fff;border:1px solid var(--line,#c9ccd4);padding:8px;min-height:36px;border-radius:0}
.table-tools input{width:220px;max-width:100%;caret-color:var(--blue,#2855a6)}
.table-pages button:disabled{color:var(--muted,#5f6470);background:var(--wash,#f4f5f7);cursor:default}
.table-tools input:focus-visible,.table-tools select:focus-visible,.table-pages button:focus-visible,.table-sort:focus-visible{outline:3px solid #ffb800;outline-offset:2px}
.table-sort{font:inherit;letter-spacing:inherit;text-transform:inherit;text-align:left;color:inherit;border:0;background:transparent;padding:4px 0;cursor:pointer;display:inline-flex;align-items:center;gap:6px}
.table-sort:hover{color:var(--blue,#2855a6)}
.table-sort svg{width:12px;height:12px;flex:none;fill:none;stroke:currentColor;stroke-width:1.5}
.table-sort .sort-up,.table-sort .sort-down{opacity:.35}
th[aria-sort=ascending] .sort-up,th[aria-sort=descending] .sort-down{opacity:1}
.table-pages [role=status]{margin-right:auto;font-variant-numeric:tabular-nums}
.table-mobile-sort{display:none!important}
tr[data-table-row][hidden],tr[data-table-empty][hidden]{display:none!important}
@media(max-width:800px){.table-tools{align-items:stretch}.table-tools label{flex-wrap:wrap}.table-tools input{width:180px}.table-mobile-sort{display:flex!important}}
@media print{.table-tools,.table-pages{display:none!important}.table-sort svg{display:none}.table-sort{padding:0}tr[data-table-row][hidden]{display:table-row!important}tr[data-table-empty]{display:none!important}}
</style>
"""

_TABLE_SCRIPT = r"""
<script data-chronoclade-tables>
(() => {
  const collator = new Intl.Collator(undefined, {numeric:true,sensitivity:'base'});
  let serial = 0;
  const missing = /^(?:|unknown|unavailable|not reported|none|null|n\/a|na|—|–|-)$/i;
  const numeric = value => {
    const text = value.replace(/,/g,'').replace(/%$/,'').trim();
    return /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(text) ? Number(text) : null;
  };
  function enhance(table) {
    if (table.dataset.tableReady || table.dataset.noControls !== undefined) return;
    const heads = [...(table.tHead?.rows[0]?.cells || [])];
    const body = table.tBodies[0];
    if (!body || !heads.length || table.tHead.rows.length !== 1 || table.tBodies.length !== 1) return;
    const rows = [...body.rows];
    if (rows.length < 2 || heads.some(h => h.colSpan !== 1 || h.rowSpan !== 1) ||
        rows.some(r => r.cells.length !== heads.length || [...r.cells].some(c=>c.colSpan!==1 || c.rowSpan!==1))) return;
    table.dataset.tableReady = 'true';
    if(!table.id){let id;do{id='cc-table-'+(++serial);}while(document.getElementById(id));table.id=id;}
    const labels = heads.map(h=>h.textContent.trim());
    const records = rows.map((row,index)=>{const values=[...row.cells].map(c=>c.textContent.trim());return {row,index,values,text:values.join(' ').toLocaleLowerCase()};});
    rows.forEach(r=>r.dataset.tableRow='');
    let column = -1, direction = 1, page = 0, size = 10, query = '';
    const large = rows.length > 10;
    let status, previous, next, sortMenu;
    const element = (tag,text) => {const e=document.createElement(tag);if(text)e.textContent=text;return e;};
    const linkControl = e => e.setAttribute('aria-controls',table.id);
    const empty = element('tr');empty.dataset.tableEmpty='';empty.hidden=true;
    const emptyCell=element('td','No matching rows');emptyCell.colSpan=heads.length;empty.append(emptyCell);
    body.append(empty);
    const compare = (a,b) => {
      const left=a.values[column],right=b.values[column];
      const lm=missing.test(left),rm=missing.test(right);
      if(lm!==rm)return lm?1:-1;
      const ln=numeric(left),rn=numeric(right);
      const diff=ln!==null&&rn!==null ? ln-rn : collator.compare(left,right);
      return direction*diff || a.index-b.index;
    };
    function render() {
      const ordered=column<0?[...records]:[...records].sort(compare);
      const terms=query.split(/\s+/).filter(Boolean);
      const selected=ordered.filter(r=>terms.every(term=>r.text.includes(term)));
      const pages=Math.max(1,Math.ceil(selected.length/size));page=Math.min(page,pages-1);
      const start=page*size,end=Math.min(start+size,selected.length);
      const shown=new Set(selected.slice(start,end));
      const fragment=document.createDocumentFragment();
      ordered.forEach(r=>{r.row.hidden=large&&!shown.has(r);fragment.append(r.row);});
      fragment.append(empty);body.append(fragment);empty.hidden=selected.length!==0;
      heads.forEach((h,i)=>{h.setAttribute('aria-sort',column===i?(direction===1?'ascending':'descending'):'none');
        h.querySelector('button').setAttribute('aria-label','Sort by '+labels[i]+(column===i&&direction===1?', descending':', ascending'));});
      if(sortMenu)sortMenu.value=column<0?'':String(column);
      if(status){status.textContent=selected.length ? `${start+1}–${end} of ${selected.length} rows` : '0 matching rows';
        if(query)status.textContent+=` (filtered from ${rows.length})`;
        previous.disabled=page===0;next.disabled=page>=pages-1;}
    }
    function sort(index){direction=column===index?-direction:1;column=index;page=0;render();}
    heads.forEach((h,i)=>{
      const button=element('button');button.type='button';button.className='table-sort';
      while(h.firstChild)button.append(h.firstChild);
      button.insertAdjacentHTML('beforeend','<svg aria-hidden="true" viewBox="0 0 12 16"><path class="sort-up" d="M2 6l4-4 4 4"/><path class="sort-down" d="M2 10l4 4 4-4"/></svg>');
      linkControl(button);h.append(button);button.addEventListener('click',()=>sort(i));
    });
    if(large){
      const tools=element('div');tools.className='table-tools';
      const searchLabel=element('label','Search this table');const search=element('input');search.type='search';linkControl(search);
      searchLabel.append(search);tools.append(searchLabel);
      search.addEventListener('input',()=>{query=search.value.trim().toLocaleLowerCase();page=0;render();});
      const sizeLabel=element('label','Rows');const select=element('select');linkControl(select);
      [10,25,50,100,'All'].forEach(v=>{const o=element('option',String(v));o.value=String(v);select.append(o);});select.value='10';
      sizeLabel.append(select);tools.append(sizeLabel);select.addEventListener('change',()=>{size=select.value==='All'?rows.length:Number(select.value);page=0;render();});
      const sortLabel=element('label','Sort by');sortLabel.className='table-mobile-sort';sortMenu=element('select');linkControl(sortMenu);
      const original=element('option','Original order');original.value='';sortMenu.append(original);
      labels.forEach((label,i)=>{const o=element('option',label);o.value=String(i);sortMenu.append(o);});sortLabel.append(sortMenu);tools.append(sortLabel);
      sortMenu.addEventListener('change',()=>{column=sortMenu.value===''?-1:Number(sortMenu.value);direction=1;page=0;render();});
      const reverse=element('button','Reverse sort');reverse.type='button';reverse.className='table-mobile-sort table-sort';linkControl(reverse);tools.append(reverse);
      reverse.addEventListener('click',()=>{if(column>=0){direction*=-1;page=0;render();}});
      const wrapper=table.parentElement.classList.contains('table-scroll')?table.parentElement:table;
      wrapper.before(tools);
      const pages=element('div');pages.className='table-pages';status=element('span');status.setAttribute('role','status');status.setAttribute('aria-live','polite');
      previous=element('button','Previous');next=element('button','Next');
      [previous,next].forEach(b=>{b.type='button';linkControl(b);});
      previous.addEventListener('click',()=>{page--;render();});next.addEventListener('click',()=>{page++;render();});pages.append(status,previous,next);wrapper.after(pages);
    }
    render();
  }
  document.querySelectorAll('table').forEach(enhance);
  new MutationObserver(changes=>changes.forEach(change=>change.addedNodes.forEach(node=>{
    if(node.nodeType!==1)return;
    if(node.matches('table'))enhance(node);
    node.querySelectorAll('table').forEach(enhance);
  }))).observe(document.body,{childList:true,subtree:true});
})();
</script>
"""


def enhance_report_tables(document: str) -> str:
    """Progressively enhance plain tables; no JavaScript still shows every row."""
    if 'data-chronoclade-tables' in document:
        return document
    document = document.replace('</head>', _TABLE_CSS + '</head>', 1)
    before, closing, after = document.rpartition('</body>')
    return before + _TABLE_SCRIPT + closing + after if closing else document


def write_report_html(path, document: str, *, encoding="utf-8"):
    """Write an offline report with the common table behaviour."""
    from chronoclade.report_components.portable import embed_report_assets
    from chronoclade.report_components.downloads import enhance_report_downloads

    document = enhance_report_downloads(document, Path(path).parent)
    document = embed_report_assets(document, Path(path).parent)
    return Path(path).write_text(enhance_report_tables(document), encoding=encoding)
