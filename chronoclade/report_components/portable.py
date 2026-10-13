"""Embed report display assets so a single HTML file can be viewed offline."""

import base64
import json
import re
from html import escape, unescape
from pathlib import Path
from urllib.parse import unquote, urlsplit

_MEDIA = {'.svg': 'image/svg+xml', '.png': 'image/png', '.jpg': 'image/jpeg',
          '.jpeg': 'image/jpeg', '.webp': 'image/webp', '.gif': 'image/gif'}
_ATTRIBUTE = re.compile(r'\b(src|href)="([^"]*)"')
_FRAME = re.compile(r'<iframe\b([^>]*?)\bsrc="([^"]*)"([^>]*)>', re.IGNORECASE)
_PAYLOAD = re.compile(r'(<script\b[^>]*type="application/json"[^>]*>)(.*?)(</script>)', re.DOTALL)


def embed_report_assets(document: str, directory: Path) -> str:
    """Embed only local display files within the report directory.

    Evidence downloads and navigation remain relative links. Missing display assets
    fail generation instead of producing a completed report with broken pictures.
    """
    directory = Path(directory).resolve()

    def local_file(raw):
        url = urlsplit(unescape(raw))
        if url.scheme or url.netloc or not url.path or url.path.startswith('/'):
            return None
        path = (directory / unquote(url.path)).resolve()
        if not path.is_relative_to(directory):
            return None
        return path

    def attribute(match):
        kind, raw = match.groups()
        path = local_file(raw)
        if path is None or path.suffix.lower() not in _MEDIA:
            return match.group()
        if not path.is_file():
            raise ValueError(f'Report display asset is missing: {raw}')
        encoded = base64.b64encode(path.read_bytes()).decode('ascii')
        return (f'{kind}="data:{_MEDIA[path.suffix.lower()]};base64,{encoded}"'
                f' data-source-asset="{escape(raw, quote=True)}"')

    def frame(match):
        before, raw, after = match.groups()
        path = local_file(raw)
        if path is None or path.suffix.lower() != '.html':
            return match.group()
        if not path.is_file():
            raise ValueError(f'Report interactive view is missing: {raw}')
        content = embed_report_assets(path.read_text(encoding='utf-8'), path.parent)
        return f'<iframe{before}srcdoc="{escape(content, quote=True)}"{after}>'

    def transform(value):
        if isinstance(value, str):
            return _ATTRIBUTE.sub(attribute, _FRAME.sub(frame, value))
        if isinstance(value, list):
            return [transform(item) for item in value]
        if isinstance(value, dict):
            return {key: transform(item) for key, item in value.items()}
        return value

    def payload(match):
        try:
            value = json.loads(match[2])
        except json.JSONDecodeError:
            return match.group()
        encoded = json.dumps(transform(value), ensure_ascii=True, separators=(",", ":")).replace('<', '\\u003c')
        return match[1] + encoded + match[3]

    # JSON contains HTML for dynamically selected network views.
    document = _PAYLOAD.sub(payload, document)
    return _ATTRIBUTE.sub(attribute, _FRAME.sub(frame, document))
