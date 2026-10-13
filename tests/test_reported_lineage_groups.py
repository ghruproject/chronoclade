from chronoclade.cglin import DEFAULT_SCHEME
from chronoclade.report_components.lineage_groups import assign_lineage_groups
from chronoclade.report_components.profile_overview import grouping_rule


def record(name, code):
    return {'sample_id': name, 'cglin_raw': code, 'cglin_scheme': DEFAULT_SCHEME,
            'cglin_scheme_version': 'test', 'cglin_status': 'complete'}


def test_display_membership_uses_codes_not_allele_distance_cutoff():
    rows = [record('a', '0,0,107,0,0,0,0,0'), record('b', '0,0,107,0,0,0,0,1'),
            record('c', '0,0,107,0,0,0,1,0')]
    rows[0]['cgmlst_profile'] = {'x': 1}
    rows[1]['cgmlst_profile'] = {'x': 99999}
    data = {'partition': {'kind': 'cglin'}, 'distance_threshold': .02}
    assignments = assign_lineage_groups(data, rows)
    assert assignments['a'] == assignments['b']
    assert assignments['a'] != assignments['c']
    assert len(data['lineage_groups']) == 2
    assert data['display_group_level'] == 7
    assert 'first 7 cgLIN positions match' in grouping_rule(data)


def test_unassigned_code_stays_unassigned():
    data = {'partition': {'kind': 'cglin'}}
    rows = [record('a', '0,0,107,0,0,0,0'), record('b', '0,0,107,0,0')]
    assignments = assign_lineage_groups(data, rows)
    assert assignments['b'] == 'cgLIN 7: unassigned'
    assert assignments['a'] != assignments['b']


def test_non_lin_input_is_not_given_a_fabricated_lin_group():
    data = {'partition': {'kind': 'hiercc'}}
    assert assign_lineage_groups(data, [record('a', '0,0,107,0,0,0,0')]) is None
    assert 'lineage_groups' not in data
