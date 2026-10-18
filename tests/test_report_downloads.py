"""Small offline evidence downloads keep exact bytes and navigation fallbacks."""

import base64
import json
import re

from chronoclade.report_components.downloads import enhance_report_downloads
from chronoclade.report_components.group_tree import grouping_caption
from chronoclade.report_components.tables import write_report_html


def test_explicit_text_download_retains_href_and_embeds_exact_bytes(tmp_path):
    content = b'{"sample": "q", "status": "complete"}\n'
    (tmp_path / 'selection.json').write_bytes(content)
    document = '<html><body><a href="selection.json" download>Download manifest</a></body></html>'
    result = enhance_report_downloads(document, tmp_path)
    assert 'href="selection.json" download' in result
    payload = re.search(r'data-cc-download="([^"]+)"', result)[1]
    assert base64.b64decode(payload) == content
    assert 'data-cc-download-name="selection.json"' in result
    assert 'URL.createObjectURL' in result and 'save.download=' in result
    assert enhance_report_downloads(result, tmp_path) == result


def test_large_missing_nontext_and_outside_targets_keep_native_fallbacks(tmp_path):
    (tmp_path / 'large.csv').write_bytes(b'x' * (1024 * 1024 + 1))
    (tmp_path / 'alignment.fasta').write_text('>q\nACT\n')
    (tmp_path / 'picture.svg').write_text('<svg/>')
    (tmp_path.parent / 'outside.json').write_text('{}')
    links = ''.join(f'<a href="{path}" download>Evidence</a>' for path in [
        'large.csv', 'alignment.fasta', 'picture.svg', '../outside.json', 'missing.json',
        'https://example.com/file.json',
    ])
    document = '<html><body>' + links + '</body></html>'
    assert enhance_report_downloads(document, tmp_path) == document


def test_download_script_does_not_corrupt_network_html_export_source(tmp_path):
    (tmp_path / 'tree.nwk').write_text('(a:1,b:1);')
    document = ('<html><head></head><body><script>const network="</body></html>";</script>'
                '<a href="tree.nwk" download>Tree</a></body></html>')
    path = tmp_path / 'report.html'
    write_report_html(path, document)
    result = path.read_text()
    assert '<script>const network="</body></html>";</script>' in result
    assert result.index('id="cc-offline-downloads"') > result.index('<a href="tree.nwk"')
    assert 'data-cc-download=' in result


def test_profile_distance_tree_caption_declares_its_method_and_cutoff():
    title, note = grouping_caption(['Distance group 1'], 'exploratory_profile_distance', .02)
    assert title == 'Genetic groups on NJ tree'
    assert 'direct-profile distance threshold' in note and '0.02' in note
    assert 'no grouping method' not in note


def test_normal_profile_refresh_preserves_scope_and_numerical_evidence(tmp_path):
    from chronoclade.cgmlst.refresh import refresh_clock_reports
    from chronoclade.cgmlst.workflow import run_cgmlst
    from chronoclade.datasets import LocusCatalogue, from_profile_records, write_dataset

    rows = [{'sample_id': f'q{i}', 'species': 'Staphylococcus aureus', 'role': 'input',
             'cgmlst_scheme': 'test', 'cgmlst_scheme_version': '1',
             'cgmlst_profile': {'a': str(i + 1), 'b': '1'}, 'cgmlst_status': 'resolved'}
            for i in range(3)]
    dataset = write_dataset(from_profile_records(rows, [
        LocusCatalogue('test', '1', ('a', 'b'), 'test'),
    ]), tmp_path / 'prepared')
    stage = run_cgmlst(dataset, tmp_path / 'analysis', bootstrap_replicates=0)
    manifest = json.loads(stage.manifest_path.read_text())
    directory = stage.manifest_path.parent / manifest['blocks'][0]['block_id']
    numerical = directory / 'numerical_analysis.json'
    before = numerical.read_bytes()
    refresh_clock_reports(stage.manifest_path)
    index = stage.report_path.read_text()
    assert 'Analysis: species and compatible cgMLST scheme' in index
    assert 'Analysis: LIN' not in index
    assert 'optional annotations' in index
    assert numerical.read_bytes() == before
