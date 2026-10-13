"""A report carries its display assets even when copied away from its bundle."""

import base64
import json
import re

import pytest

from chronoclade.report_components.tables import write_report_html


def test_single_file_embeds_figures_frames_and_dynamic_network_views(tmp_path):
    svg = b'<svg xmlns="http://www.w3.org/2000/svg"><text>Country</text></svg>'
    (tmp_path / 'country.svg').write_bytes(svg)
    (tmp_path / 'ordination.html').write_text('<html><body><button>Date</button></body></html>')
    payload = json.dumps([{'all': '<img src="country.svg">'}])
    document = ('<html><head></head><body><img src="country.svg">'
                '<a href="country.svg">Full size</a><a href="distances.csv">Evidence</a>'
                '<iframe src="ordination.html"></iframe>'
                f'<script type="application/json" data-network-views>{payload}</script></body></html>')
    output = tmp_path / 'report.html'
    write_report_html(output, document)
    result = output.read_text()
    assert 'src="country.svg"' not in result
    assert 'src="ordination.html"' not in result
    assert 'srcdoc="' in result and 'Date' in result
    assert 'href="distances.csv"' in result
    assert base64.b64encode(svg).decode() in result
    dynamic = json.loads(re.search(r'data-network-views>(.*?)</script>', result).group(1))
    assert 'src="data:image/svg+xml;base64,' in dynamic[0]['all']
    isolated = tmp_path / 'copy' / 'report.html'
    isolated.parent.mkdir()
    isolated.write_text(result)
    assert isolated.read_text() == result


def test_missing_local_display_asset_fails_instead_of_broken_image(tmp_path):
    with pytest.raises(ValueError, match='display asset is missing'):
        write_report_html(tmp_path / 'report.html', '<html><body><img src="missing.svg"></body></html>')


def test_table_script_follows_export_code_containing_html_body_tag(tmp_path):
    output = tmp_path / 'report.html'
    write_report_html(output, '<html><head></head><body><script>const closing="</body>";</script></body></html>')
    result = output.read_text()
    assert 'const closing="</body>";</script>' in result
    assert result.index('const closing') < result.index('<script data-chronoclade-tables>')
