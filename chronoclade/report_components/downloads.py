"""Make small explicit evidence downloads work when reports open via file://."""

import base64
from html import escape, unescape
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

_ANCHOR = re.compile(r'<a\b[^>]*>', re.IGNORECASE)
_HREF = re.compile(r'\bhref="([^"]*)"', re.IGNORECASE)
_DOWNLOAD = re.compile(r'\sdownload(?:\s|=|>)', re.IGNORECASE)
_TYPES = {'.json': 'application/json', '.nwk': 'text/plain', '.newick': 'text/plain',
          '.nex': 'text/plain', '.csv': 'text/csv', '.tsv': 'text/tab-separated-values'}
_LIMIT = 1024 * 1024

_SCRIPT = '''<script id="cc-offline-downloads">
(function(){
 document.addEventListener('click',function(event){
  const link=event.target.closest?.('a[data-cc-download]');
  if(!link||event.button!==0||event.ctrlKey||event.metaKey||event.shiftKey||event.altKey)return;
  event.preventDefault();
  const bytes=Uint8Array.from(atob(link.dataset.ccDownload),letter=>letter.charCodeAt(0));
  const url=URL.createObjectURL(new Blob([bytes],{type:link.dataset.ccDownloadType}));
  const save=document.createElement('a');save.href=url;save.download=link.dataset.ccDownloadName;
  document.body.append(save);save.click();save.remove();
  setTimeout(()=>URL.revokeObjectURL(url),1000);
 });
})();
</script>'''


def enhance_report_downloads(document, directory):
    """Retain source hrefs; embed bounded text evidence for browser Blob saves.

    Image/viewer navigation and large matrices are unchanged. The frozen files
    remain authoritative; embedded payloads are exact byte copies at generation.
    """
    directory = Path(directory).resolve()

    def anchor(match):
        tag = match.group()
        href = _HREF.search(tag)
        if not href or not _DOWNLOAD.search(tag) or 'data-cc-download=' in tag:
            return tag
        url = urlsplit(unescape(href[1]))
        if url.scheme or url.netloc or not url.path:
            return tag
        path = (directory / unquote(url.path)).resolve()
        kind = _TYPES.get(path.suffix.lower())
        if (not kind or not path.is_relative_to(directory) or not path.is_file()
                or path.stat().st_size > _LIMIT):
            return tag
        payload = base64.b64encode(path.read_bytes()).decode('ascii')
        return (tag[:-1] + f' data-cc-download="{payload}"'
                f' data-cc-download-type="{kind}"'
                f' data-cc-download-name="{escape(path.name, quote=True)}">')

    document = _ANCHOR.sub(anchor, document)
    if 'data-cc-download=' in document and 'id="cc-offline-downloads"' not in document:
        before, closing, after = document.rpartition('</body>')
        if closing:
            document = before + _SCRIPT + closing + after
    return document
