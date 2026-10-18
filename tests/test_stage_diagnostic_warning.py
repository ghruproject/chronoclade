from html.parser import HTMLParser

import pytest

from chronoclade.report_components.stage import write_report


class Figures(HTMLParser):
    def __init__(self):
        super().__init__()
        self.images = []
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == 'img':
            self.images.append(dict(attrs))
        if tag == 'a':
            self.links.append(dict(attrs))


@pytest.mark.parametrize('supported,rate,warns', [
    (False, -5.474e-8, True),
    (False, 3.4e-7, True),
    (True, 0, True),
    (True, -1e-8, True),
    (True, 3.4e-7, False),
])
def test_raw_clock_intercept_cannot_be_confused_with_accepted_mrca(
        tmp_path, supported, rate, warns):
    clock = tmp_path/'clock'
    clock.mkdir()
    raw = b'<svg xmlns="http://www.w3.org/2000/svg"><text>root date: 2193.7</text></svg>'
    (clock/'root_to_tip_regression.svg').write_bytes(raw)
    (tmp_path/'corrected_tree.svg').write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
    report = tmp_path/'report.html'
    write_report(report, title='Saved tree', result={
        'selected_sample_ids': [],
        'temporal_assessment': {'code': 'supported' if supported else 'not_supported',
                                'supported': supported, 'reason': 'screening fixture',
                                'temporal': {'observed': {'rate': rate}}}},
        figures=[('clock/root_to_tip_regression.svg', 'Root-to-tip diagnostic'),
                 ('corrected_tree.svg', 'Corrected biological tree')])
    html = report.read_text()
    assert ('class="diagnostic-warning" role="alert"' in html) == warns
    assert ('Do not use it as an ancestor date' in html) == warns
    parser = Figures()
    parser.feed(html)
    root = next(image for image in parser.images
                if image.get('data-source-asset') == 'clock/root_to_tip_regression.svg')
    assert 'regression intercept, not an accepted MRCA estimate' in root['alt']
    assert ('Dating unsupported.' in root['alt']) == warns
    tree = next(image for image in parser.images
                if image.get('data-source-asset') == 'corrected_tree.svg')
    assert tree['alt'] == 'Corrected biological tree'
    links = [link for link in parser.links
             if link.get('data-source-asset') == 'clock/root_to_tip_regression.svg']
    assert len(links) == 2  # Image and Downloads open the same warned viewer.
    assert all('regression intercept, not an accepted MRCA estimate'
               in link['data-figure-title'] for link in links)
    assert all(('Dating unsupported.' in link['data-figure-title']) == warns
               for link in links)
    assert (clock/'root_to_tip_regression.svg').read_bytes() == raw


@pytest.mark.parametrize('p_value,successful,expected', [
    (5/101, 100, '0.049505 (5/101)'),
    (77/101, 100, '0.762376 (77/101)'),
    (0.0495, 100, '0.0495'),  # Rounded p must not invent an exact numerator.
    (0.5, 9, '0.5 (5/10)'),
])
def test_saved_permutation_resolution_is_visible_at_decision(
        tmp_path, p_value, successful, expected):
    path = tmp_path/'report.html'
    assessment = {'code': 'not_supported', 'supported': False,
                  'reason': 'Negative rate stays unsupported despite displayed p.',
                  'p_value_threshold': 0.05,
                  'temporal': {'p_value_r_squared': p_value,
                               'successful_randomisations': successful,
                               'requested_randomisations': 100}}
    write_report(path, title='Permutation resolution', result={
        'selected_sample_ids': [], 'temporal_assessment': assessment})
    html = path.read_text()
    assert 'Empirical root-to-tip p: ' + expected in html
    assert 'threshold: 0.05' in html
    assert f'successful permutations: {successful}/100' in html
    if p_value == 0.0495:
        assert '(5/101)' not in html
    verdict = html.split('<div class="verdict ', 1)[1].split('</div>', 1)[0]
    assert 'Temporal signal not supported' in verdict
    assert 'temporal-screen-summary' in verdict
    assert assessment['supported'] is False
