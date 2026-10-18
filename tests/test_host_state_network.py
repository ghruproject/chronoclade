"""Host reconstruction must use host observations, including explicit missingness."""
from io import StringIO
import json
import re

from Bio import Phylo

from chronoclade.cgmlst.workflow import _add_state_networks
from chronoclade.location_network.reconstruction import build_location_network, build_state_network
from chronoclade.location_network.widget import interactive_network_html


def fixture():
    tree = Phylo.read(StringIO('((a:1,b:1):1,(c:1,d:1):1);'), 'newick')
    rows = [dict(sample_id=ident, country='UK', host=host, role=role)
            for ident, host, role in [('a', 'Human', 'input'), ('b', 'Human', 'context'),
                                     ('c', 'Chicken', 'context'), ('d', 'Unknown', 'context')]]
    return tree, rows


def test_hosts_reconstructed_independently_with_unknown_observations_retained():
    tree, rows = fixture()
    country = build_location_network(tree, rows, tree_basis='cgmlst')
    host = build_state_network(tree, rows, state_field='host', tree_basis='cgmlst')
    assert country['optimum_changes'] == 0
    assert host['optimum_changes'] == 1
    assert {node['state'] for node in host['nodes']} == {'Human', 'Chicken'}
    assert host['unknown_state_count'] == 1
    unknown = next(node for node in host['reconstruction']['nodes'] if node['name'] == 'd')
    assert unknown['observed_state'] is None and unknown['state'] is None
    assert unknown['unknown_state']
    assert host['status'] == 'available'
    assert all(edge['source'] in {'Human', 'Chicken'} and edge['target'] in {'Human', 'Chicken'}
               for edge in host['directed_edges'])
    assert rows[0]['country'] == 'UK'


def test_missing_hosts_have_explicit_unavailable_status():
    tree, rows = fixture()
    for row in rows:
        row['host'] = None
    host = build_state_network(tree, rows, state_field='host', tree_basis='cgmlst')
    assert host['status'] == 'insufficient_metadata'
    assert host['unknown_state_count'] == 4
    assert host['nodes'] == [] and host['directed_edges'] == []


def test_workflow_host_views_and_widget_exports(tmp_path):
    tree, rows = fixture()
    location = build_location_network(tree, rows, tree_basis='cgmlst')
    location.update(cohort_id='CG_1', id='all', label='All inputs')
    analysis = dict(location_network=[location], paths={})
    _add_state_networks(analysis, rows, tmp_path)
    host = analysis['state_networks']['host'][0]
    assert json.loads((tmp_path / 'CG_1_host_network.json').read_text())['optimum_changes'] == 1
    markup = interactive_network_html(host, 'hosts')
    assert 'Highlight host' in markup and 'All hosts' in markup
    assert 'Host colour key' in markup and 'Arrows show host-state changes' in markup
    assert 'Save PNG' in markup and 'Save interactive HTML' in markup
    payload = json.loads(re.search(r'data-network-data>(.*?)</script>', markup, re.S).group(1))
    assert payload['stateLabel'] == 'host'
    assert all('country' not in node for node in payload['nodes'])
    assert {node['id'] for node in payload['nodes']} == {'Human', 'Chicken'}


def test_host_case_variants_and_metadata_summary_are_consistent():
    from chronoclade.cgmlst.workflow import add_metadata_summary
    tree, rows = fixture()
    rows[1]['host'] = 'human'
    rows[0]['collection_date'] = '2021-03-12'
    rows[1]['collection_date'] = '2021'
    rows[2]['collection_date'] = '2022-05'
    analysis = add_metadata_summary({}, rows)
    assert analysis['metadata_summary']['host'] == {'total': 4, 'reported': 3, 'missing': 1}
    assert analysis['metadata_summary']['collection_date'] == {'total': 4, 'reported': 3, 'missing': 1}
    human = next(row for row in analysis['metadata_composition']
                 if row['field'] == 'Host' and row['value'] == 'Human')
    assert human == {'field': 'Host', 'value': 'Human', 'input': 1, 'context': 1, 'total': 2}
    host = build_state_network(tree, rows, state_field='host', tree_basis='cgmlst')
    assert {node['state'] for node in host['nodes']} == {'Human', 'Chicken'}


def test_nearest_host_comparison_cannot_fall_back_to_context_country():
    tree, rows = fixture()
    comparison = dict(query_id='a', context_id='outside_tree', country='France', status='matched')
    host = build_state_network(tree, rows, [comparison], state_field='host', tree_basis='cgmlst')
    assert host['nearest_edges'][0]['source'] == 'Human'
    assert host['nearest_edges'][0]['target'] == 'Unknown'


def test_optional_metadata_networks_use_subgroup_tree_and_rows(tmp_path):
    from chronoclade.cgmlst.workflow import _add_state_networks

    tree = Phylo.read(StringIO('((a:1,b:1):1,(c:1,d:1):1);'), 'newick')
    rows = [
        dict(sample_id='a', country='UK', region='North', city='Oxford', host=None, role='input'),
        dict(sample_id='b', country='UK', region='North', city='Oxford', host=None, role='context'),
        dict(sample_id='c', country='UK', region='South', city='Leeds', host=None, role='input'),
        dict(sample_id='d', country='UK', region='South', city='Leeds', host=None, role='context'),
    ]
    location = build_location_network(tree, rows, tree_basis='cgmlst')
    location.update(cohort_id='CG_1', id='all', label='All inputs')
    view_tree = Phylo.read(StringIO('(a:1,b:1);'), 'newick')
    view = build_location_network(view_tree, rows[:2], tree_basis='cgmlst')
    view.update(cohort_id='CG_1', id='lin_1', label='LIN5: 1 (1 input)')
    view['input_ids'] = ['a']
    location['views'] = [view]

    analysis = dict(location_network=[location], paths={}, nearest_neighbours=[])
    _add_state_networks(analysis, rows, tmp_path)

    assert analysis['state_networks']['host'] == []
    for field, expected in [('region', 'UK / North'), ('city', 'UK / North / Oxford')]:
        network = analysis['state_networks'][field][0]
        assert network['sample_count'] == 4
        assert network['state_field'] == field
        subgroup = network['views'][0]
        assert subgroup['sample_count'] == 2
        assert subgroup['input_ids'] == ['a']
        assert {item['state'] for item in subgroup['state_definitions']} == {expected}
        assert all(item['state'] not in {'South', 'UK / South / Leeds'}
                   for item in subgroup['state_definitions'])
        audit = json.loads((tmp_path / f'CG_1_{field}_network.json').read_text())
        assert audit == network


def test_profile_selector_exposes_only_observed_optional_fields(tmp_path, monkeypatch):
    import chronoclade.location_network.state_selector as selector
    import chronoclade.profile_report as report
    from chronoclade.profile_report import _state_network_views

    captured = {}
    rendered_fields = []

    def render_networks(value, directory, cohorts=None, *, state_field='country', include_assets=True):
        rendered_fields.append((state_field, value))
        return f'<panel data-state="{state_field}"></panel>'

    def choose(country_panel, host_panel=None, *, panels=None):
        captured['fields'] = set(panels or {})
        captured['panels'] = panels or {}
        captured['country'] = country_panel
        return '<selector></selector>'

    monkeypatch.setattr(report, '_location_network', render_networks)
    monkeypatch.setattr(selector, 'state_network_selector', choose)
    data = {
        'state_networks': {
            'country': [{'nodes': [{'country': 'UK'}], 'edges': []}],
            'host': [],
            'region': [{
                'state_field': 'region', 'state_label': 'region',
                'known_state_count': 1, 'sample_count': 4,
                'unknown_state_count': 0, 'status': 'insufficient_metadata',
                'state_plural': 'regions', 'cohort_id': 'CG_1',
            }],
            'city': [],
        }
    }

    html = _state_network_views(data, tmp_path)
    assert html == '<selector></selector>'
    assert captured['fields'] == {'region'}
    assert captured['country'] == '<panel data-state="country"></panel>'
    assert captured['panels']['region'] == '<panel data-state="region"></panel>'
    assert [field for field, _ in rendered_fields] == ['country', 'region']
