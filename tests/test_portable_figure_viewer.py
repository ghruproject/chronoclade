"""Portable reports retain inline images, evidence links and network export code."""

import base64
import re

import pytest

from chronoclade.report_components.figure_viewer import add_figure_viewer
from chronoclade.report_components.portable import embed_report_assets


def test_viewer_preserves_network_export_script_and_is_idempotent():
    document = '<html><body><script>const exported="</body>";</script><p>End</p></body></html>'
    result = add_figure_viewer(document)
    assert '<script>const exported="</body>";</script><p>End</p>' in result
    assert result.index('<p>End</p>') < result.index('id="cc-figure-viewer-script"')
    assert add_figure_viewer(result) == result
    fragment = '<a href="evidence.csv">Download</a>'
    assert add_figure_viewer(fragment) == fragment


def test_embedded_image_bytes_and_navigation_remain_portable(tmp_path):
    svg = b'<svg xmlns="http://www.w3.org/2000/svg"><text>Example</text></svg>'
    (tmp_path / 'plot.svg').write_bytes(svg)
    document = ('<html><body><figure><a href="plot.svg">Open full size</a>'
                '<img alt="Example" src="plot.svg"></figure>'
                '<a href="data.csv">Evidence</a><a href="../index.html">All groups</a>'
                '</body></html>')
    result = embed_report_assets(document, tmp_path)
    payloads = re.findall(r'data:image/svg\+xml;base64,([^\"]+)', result)
    assert len(payloads) == 2
    assert all(base64.b64decode(payload) == svg for payload in payloads)
    assert 'href="data.csv"' in result
    assert 'href="../index.html"' in result
    assert 'data-source-asset="plot.svg"' in result
    assert embed_report_assets(result, tmp_path) == result


def test_missing_image_fails_report_generation(tmp_path):
    with pytest.raises(ValueError, match='Report display asset is missing'):
        embed_report_assets('<html><body><img src="missing.svg"></body></html>', tmp_path)


def test_interactive_frame_is_embedded_with_its_image(tmp_path):
    (tmp_path / 'network.svg').write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
    (tmp_path / 'network.html').write_text(
        '<html><body><img src="network.svg"></body></html>'
    )
    result = embed_report_assets(
        '<html><body><iframe src="network.html"></iframe></body></html>', tmp_path
    )
    assert 'srcdoc=' in result
    assert 'src="network.html"' not in result
    assert 'data:image/svg+xml;base64,' in result
