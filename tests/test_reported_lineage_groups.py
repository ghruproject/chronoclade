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


def test_hiercc_display_uses_the_partition_level_and_reported_codes():
    def hiercc(name, hc, status='resolved'):
        return {'sample_id': name, 'hiercc_scheme': 'enterobase:ecoli-hiercc',
                'hiercc_scheme_version': 'v1', 'hiercc_status': status,
                'hiercc_codes': {'HC10': hc} if hc is not None else {}}

    rows = [hiercc('a', '23'), hiercc('b', '23'), hiercc('c', '91'), hiercc('d', None)]
    data = {'partition': {'kind': 'hiercc', 'level': 'HC10'}}
    assignments = assign_lineage_groups(data, rows)
    assert assignments['a'] == assignments['b'] == 'HierCC HC10: 23'
    assert assignments['c'] == 'HierCC HC10: 91'
    assert assignments['d'] == 'HierCC HC10: unassigned'
    assert data['display_group_basis'] == 'reported_hiercc'
    assert data['display_group_level'] == 'HC10'
    assert len(data['lineage_groups']) == 3
    rule = grouping_rule(data)
    assert 'reported HierCC assignments at HC10' in rule
    assert 'no profile-distance cutoff' in rule


def test_st_only_display_has_one_neutral_unassigned_lineage_category():
    data = {'partition': {'kind': 'st'}}
    rows = [record('a', '0,0,107,0,0,0,0'), record('b', '')]
    assignments = assign_lineage_groups(data, rows)
    assert assignments == {'a': 'Lineage unassigned', 'b': 'Lineage unassigned'}
    assert len(data['lineage_groups']) == 1
    assert data['lineage_groups'][0]['method'] == 'lineage not used in ST-only mode'
    rule = grouping_rule(data)
    assert 'No lineage or genetic clusters are inferred' in rule


def test_unknown_partition_kind_does_not_create_a_lineage_group():
    data = {'partition': {'kind': 'other'}}
    assert assign_lineage_groups(data, [record('a', '0,0,107,0,0,0,0')]) is None
    assert 'lineage_groups' not in data


def test_explicit_level_eight_preserves_reported_fine_prefixes():
    rows = [record('a', '0,0,105,6,0,0,0,5'), record('b', '0,0,105,6,0,0,0,50')]
    data = {'partition': {'kind': 'cglin'}}
    assigned = assign_lineage_groups(data, rows, level=8)
    assert assigned['a'] != assigned['b']
    assert assigned['a'] == 'cgLIN 8: 0.0.105.6.0.0.0.5'
    assert data['display_group_level'] == 8
    assert 'first 8 cgLIN positions match' in grouping_rule(data)


def test_level_eight_does_not_pad_short_or_missing_codes():
    rows = [record('short', '0,0,105,6,0,0,0'), record('missing', ''),
            record('failed', '0,0,105,6,0,0,0,0')]
    rows[2]['cglin_status'] = 'failed'
    assigned = assign_lineage_groups({'partition': {'kind': 'cglin'}}, rows, level=8)
    assert set(assigned.values()) == {'cgLIN 8: unassigned'}


def test_display_level_rejects_invalid_depths():
    import pytest
    for level in [4, 11, 7.5, '8', True]:
        with pytest.raises(ValueError, match='integer between 5 and 10'):
            assign_lineage_groups({'partition': {'kind': 'cglin'}}, [], level=level)
