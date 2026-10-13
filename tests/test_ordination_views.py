from pathlib import Path
import json
import re

import numpy as np
import pytest

from chronoclade.report_components.ordination import ordination_viewer, write_ordination_views


def test_shared_coordinates_metadata_views_keep_missing_dates(tmp_path):
    records = [
        dict(sample_id='input', role='input', country='Greece', collection_date='2020'),
        dict(sample_id='context', role='context', country='Italy', date_start='2021-01-01',
             date_end='2022-12-31', date_precision='interval', region='Lazio', host='human'),
        dict(sample_id='missing', role='context'),
    ]
    paths = write_ordination_views(tmp_path, 'cohort', np.array([[0, 0], [1, 2], [3, 2]]), records)
    assert len(paths) == 6
    assert all(Path(path).exists() for path in paths.values())
    html = Path(paths['pcoa_views_html']).read_text()
    assert 'Colour by' in html and 'Country' in html and 'Region' in html and 'Host' in html and 'Date' in html
    region = Path(paths['pcoa_region_figure']).read_text()
    assert 'Lazio' in region and 'Unknown region' in region
    date = Path(paths['pcoa_date_figure']).read_text()
    assert 'Date unknown' in date and 'Input genomes' in date and 'Public comparisons' in date
    assert 'iframe' in ordination_viewer([paths], tmp_path)
    assert ordination_viewer([dict(pcoa_views_html='../escape.html')], tmp_path) == ''


def test_invalid_coordinates_and_names_rejected(tmp_path):
    with pytest.raises(ValueError, match='two finite'):
        write_ordination_views(tmp_path, 'safe', [[np.nan, 0]], [dict(sample_id='a')])
    with pytest.raises(ValueError, match='safe artifact'):
        write_ordination_views(tmp_path, '../bad', [[0, 0]], [dict(sample_id='a')])


def test_many_categories_keep_all_points_and_bound_legend(tmp_path):
    rows = [dict(sample_id=f's{i}', role='context') for i in range(25)]
    paths = write_ordination_views(tmp_path, 'many', np.arange(50).reshape(25,2), rows,
                                  groups={row['sample_id']:f'g{i}' for i,row in enumerate(rows)})
    svg = Path(paths['pcoa_group_figure']).read_text()
    assert 'Legend shows 20 of 25 categories' in svg
    assert 'all genomes are plotted' in svg


def test_interactive_viewer_exposes_sample_metadata_and_coordinate_bases(tmp_path):
    rows = [
        dict(sample_id='PW_local', role='input', origin='query', country='UK', host='human',
             collection_date='2020', group='cohort_1_group_1'),
        dict(sample_id='public-1', role='context', accession='SRR123456', country='France',
             host='bird', collection_date='2021', group='cohort_1_group_2'),
    ]
    paths = write_ordination_views(
        tmp_path, 'interactive', [[0, 0], [1, 1]], rows,
        groups={'PW_local': 'cohort_1_group_1', 'public-1': 'cohort_1_group_2'},
        coordinate_views={
            'nj': {'coordinates': [[2, 3], [4, 5]], 'label': 'Neighbour joining',
                   'axis_labels': ('NJ1', 'NJ2'), 'note': 'Tree-derived view'},
            'direct': {'coordinates': [[8, 9], [6, 7]], 'label': 'Direct distances',
                       'axis_labels': ('D1', 'D2'), 'note': 'Distance view'},
        },
        default_basis='nj',
    )
    html = Path(paths['pcoa_views_html']).read_text()
    payload = re.search(r'<script type="application/json" id="ordination-data">(.*?)</script>',
                        html).group(1)
    data = json.loads(payload)
    assert data['default_basis'] == 'nj'
    assert list(data['bases']) == ['nj', 'direct']
    assert data['default_view'] == 'group'
    assert data['rows'][1]['label'] == 'SRR123456'
    assert data['rows'][1]['group'] == 'Group 2'
    assert data['rows'][1]['group_id'] == 'cohort_1_group_2'
    assert data['rows'][0]['coordinates']['direct'] == [8.0, 9.0]
    assert 2020.49 < data['rows'][0]['year'] < 2020.51
    assert data['bases']['nj']['axis_labels'] == ['NJ1', 'NJ2']
    static_svg = Path(paths['pcoa_figure']).read_text()
    assert 'NJ1' in static_svg and 'NJ2' in static_svg
    assert 'id="basis"' in html and 'id="colour"' in html
    assert 'mouseenter' in html and 'showRows' in html and 'textContent' in html
    assert 'Colour category' in html and 'Samples at/near this point' in html
    assert 'aria-pressed' in html and 'Clear highlight' in html
    assert 'dateColour(r.year' in html and 'Math.max(0,Math.min(1' in html
    assert 'linear-gradient(90deg' in html and 'max-height:150px' in html and 'max-height:180px' in html
    assert 'toPrecision(4)' in html and '#e4e6eb' in html


def test_interactive_payload_escapes_markup_and_validates_coordinate_bases(tmp_path):
    rows = [dict(sample_id='a', role='context', country='</script><img src=x>')]
    path = write_ordination_views(tmp_path, 'safe', [[0, 0]], rows)['pcoa_views_html']
    html = Path(path).read_text()
    assert '</script><img src=x>' not in html
    assert r'\u003c/script\u003e' in html
    with pytest.raises(ValueError, match='coordinate basis'):
        write_ordination_views(tmp_path, 'badbasis', [[0, 0]], rows,
                               coordinate_views={'nj': {'coordinates': [[np.inf, 0]]}})
