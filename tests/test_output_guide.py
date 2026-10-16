"""Output navigation keeps checksum-bound artifacts at their original paths."""
import re
from urllib.parse import unquote

from chronoclade.report_components.output_guide import write_output_guide


def test_grouped_navigation_preserves_files_and_resolves_links(tmp_path):
    originals = {
        'profile_report.html': '<html>Report</html>',
        'cohort_nj.nwk': '(a:1,b:1);',
        'nearest_neighbours.csv': 'a,b\n1,2\n',
        'selection_bundle/ensemble.json': '{}',
        'binary_distances/distance.npy': 'frozen',
        'location network.html': '<html>Network</html>',
    }
    for relative, content in originals.items():
        path = tmp_path / relative
        path.parent.mkdir(exist_ok=True)
        path.write_text(content)
    write_output_guide(tmp_path)
    write_output_guide(tmp_path)  # Repeat refresh must not recursively list guides.
    for relative, content in originals.items():
        assert (tmp_path / relative).read_text() == content
    for page in (tmp_path / 'files').rglob('*.html'):
        for href in re.findall(r'<a href="([^"]+)"', page.read_text()):
            assert (page.parent / unquote(href)).is_file(), href
    network_index = (tmp_path / 'files/networks/index.html').read_text()
    assert 'location%20network.html' in network_index
    assert 'Archivo' in network_index
    assert 'frozen numerical checkpoint' in (tmp_path / 'README.md').read_text()
