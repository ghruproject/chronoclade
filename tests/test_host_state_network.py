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
