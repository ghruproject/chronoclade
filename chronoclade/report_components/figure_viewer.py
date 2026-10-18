"""Offline full-size viewing for embedded report figures."""

_VIEWER = r'''<style id="cc-figure-viewer-style">
#cc-figure-viewer {box-sizing:border-box;width:96vw;max-width:96vw;height:94vh;
 max-height:94vh;padding:0;border:1px solid #c9ccd4;background:#fff;color:#17191f;}
#cc-figure-viewer::backdrop {background:rgb(23 25 31 / 55%);}
#cc-figure-viewer[open] {display:grid;grid-template-rows:auto minmax(0,1fr);}
#cc-figure-viewer .cc-figure-toolbar {display:flex;align-items:center;gap:1rem;
 flex-wrap:wrap;padding:1rem;border-bottom:1px solid #c9ccd4;font:inherit;}
#cc-figure-viewer .cc-figure-title {flex:1;min-width:0;overflow-wrap:anywhere;}
#cc-figure-viewer .cc-figure-scroll {overflow:auto;min-height:0;}
#cc-figure-viewer img {display:block;max-width:none;max-height:none;margin:0;}
#cc-figure-viewer button,#cc-figure-viewer a {font:inherit;color:#2855a6;
 background:#fff;border:1px solid #c9ccd4;padding:.5rem .75rem;cursor:pointer;}
#cc-figure-viewer :focus-visible {outline:3px solid #ffb800;outline-offset:2px;}
@media(max-width:600px) {#cc-figure-viewer .cc-figure-title {flex-basis:100%;}}
@media print {#cc-figure-viewer {display:none;}}
</style>
<script id="cc-figure-viewer-script">
(()=>{
 let viewer, image, title, save, opener, downloadUrl;
 function build(){
  viewer=document.createElement('dialog');
  viewer.id='cc-figure-viewer';
  viewer.setAttribute('aria-label','Full-size figure');
  const toolbar=document.createElement('div');toolbar.className='cc-figure-toolbar';
  title=document.createElement('strong');title.className='cc-figure-title';
  save=document.createElement('a');save.textContent='Save figure';
  const close=document.createElement('button');close.type='button';
  close.textContent='Close figure';close.autofocus=true;
  close.addEventListener('click',()=>viewer.close());
  toolbar.append(title,save,close);
  const scroll=document.createElement('div');scroll.className='cc-figure-scroll';
  scroll.tabIndex=0;scroll.setAttribute('aria-label','Full-size image; scroll to explore');
  image=document.createElement('img');scroll.append(image);
  viewer.append(toolbar,scroll);document.body.append(viewer);
  viewer.addEventListener('close',()=>{
   image.removeAttribute('src');save.removeAttribute('href');
   if(downloadUrl){URL.revokeObjectURL(downloadUrl);downloadUrl=null;}
   if(opener?.isConnected)opener.focus();
  });
 }
 document.addEventListener('click',event=>{
  const link=event.target.closest?.('a[data-source-asset]');
  if(!link || !/^data:image\//i.test(link.getAttribute('href')||''))return;
  if(event.button!==0 || event.ctrlKey || event.metaKey || event.shiftKey || event.altKey)return;
  if(typeof HTMLDialogElement==='undefined')return;
  event.preventDefault();if(!viewer)build();opener=link;
  const source=link.getAttribute('data-source-asset');
  let name=source.split('/').pop().split(/[?#]/)[0];
  try{name=decodeURIComponent(name);}catch{}
  const label=link.querySelector('img')?.alt || link.closest('figure')?.querySelector('img')?.alt || name;
  title.textContent=label;image.alt=label;image.src=link.href;
  const [header,encoded]=link.href.split(',',2);
  const bytes=Uint8Array.from(atob(encoded),letter=>letter.charCodeAt(0));
  downloadUrl=URL.createObjectURL(new Blob([bytes],{type:header.slice(5).split(';')[0]}));
  save.href=downloadUrl;save.download=name;viewer.showModal();
 });
})();
</script>'''


def add_figure_viewer(document: str) -> str:
    """Enhance a complete report once, leaving fragments and evidence links alone."""
    if 'id="cc-figure-viewer-script"' in document or '</body>' not in document:
        return document
    before, closing, after = document.rpartition('</body>')
    return before + _VIEWER + '\n' + closing + after
