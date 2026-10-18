"""Reporting annotations never redefine distance ranking, identities or code scope."""
from copy import deepcopy

import pytest

from chronoclade.cgmlst.report import write_partition_index
from chronoclade.profile_report import write_profile_report
from chronoclade.query_typing import _assignment_result
from chronoclade.report_components.relatedness import annotate_neighbours, prefix_relationship, scope_summary


def row(ident, code='0,0,105,6,0,0,0,5,1,0', **extra):
    return dict(sample_id=ident, role='input' if ident == 'q' else 'context',
                cglin_raw=code, cglin_status='complete', cglin_scheme='scgMLST629_S',
                cglin_scheme_version='v1', **extra)


def test_same_different_and_unknown_prefix_depths_without_code_mutation():
    q = row('q')
    different = row('c', '0,0,105,6,0,0,1,5,1,0')
    partial = row('p', '0,0,105,6,0,0,0,5,*,*')
    before = deepcopy([q, different, partial])
    assert prefix_relationship(q, q, 7) == ('Same LIN7 prefix', '10')
    assert prefix_relationship(q, different, 7) == ('Different LIN7 prefixes', '6')
    assert prefix_relationship(q, partial, 7) == ('Same LIN7 prefix', '≥8')
    assert prefix_relationship(q, partial, 9) == ('Unresolved LIN9 prefix', '≥8')
    assert [q, different, partial] == before


def test_incompatible_versions_fingerprints_and_external_provisional_are_not_claimed():
    q = row('q')
    assert prefix_relationship(q, dict(row('c'), cglin_scheme_version='v2'), 7)[0] == 'Incompatible typing scopes'
    assert prefix_relationship(q, dict(row('c'), cglin_database_sha256='a' * 64), 7) == (
        'Reference compatibility unverified', 'Not assessed'
    )
    assert prefix_relationship(dict(q, cglin_database_sha256='b' * 64),
                               dict(row('c'), cglin_database_sha256='a' * 64), 7) == (
        'Incompatible reference versions', 'Not comparable'
    )
    assert prefix_relationship(q, dict(row('c'), cglin_provisional=True), 7)[0] == 'Unresolved LIN7 prefix'
    assert prefix_relationship(q, {}, 7)[0] == 'Unresolved typing scope'


def test_proven_conservative_prefix_honours_known_depth_only():
    raw = {'LINcode': list('0010000510'), 'matches': [{'LINcode': list('0010000520')}]}
    assignment = _assignment_result(raw, {'kind':'plincer', 'scheme':'scgMLST629_S', 'scheme_version':'v1'})
    q = row('q', '0,0,1,0,0,0,0,5,1,0')
    other = dict(assignment, sample_id='c', role='context')
    assert prefix_relationship(q, other, 8) == ('Same LIN8 prefix', '≥8')
    assert prefix_relationship(q, other, 9) == ('Unresolved LIN9 prefix', '≥8')


def test_annotation_preserves_order_distances_and_ties():
    neighbours = [dict(query_id='q', context_id=c, status='matched', distance=0, allele_differences=0) for c in ['c2','c1']]
    before = deepcopy(neighbours)
    result = annotate_neighbours(neighbours, [row('q'),row('c1'),row('c2')], level=7, kind='cglin')
    assert neighbours == before
    partial = row('p', '0,0,105,6,0,0,0,5,*,*')
    shown = annotate_neighbours([dict(neighbours[0], context_id='p')], [row('q'),partial], level=7, kind='cglin')
    assert shown[0]['lin_relationship'] == 'Same LIN7 prefix · at least 8 shared levels'
    assert [{k:v for k,v in r.items() if k in before[i]} for i,r in enumerate(result)] == before
    assert [r['context_id'] for r in result] == ['c2','c1']


@pytest.mark.parametrize('kind,parameters,words', [
    ('cglin', {'lin_level':5,'display_lin_level':7}, ['Search: not run', 'Analysis: LIN5', 'Colours: LIN7', 'background']),
    ('cglin', {'lin_level':5,'context_lin_level':8,'display_lin_level':9}, ['LIN8 input-prefix union','Analysis: LIN5','Colours: LIN9']),
    ('st', {'partition_by':'st'}, ['Search: not run','supplied profiles by species/ST','lineage unassigned']),
    ('hiercc', {'hiercc_level':'HC1100'}, ['Analysis: HC1100','Colours: HC1100']),
])
def test_scope_modes_are_explicit_and_do_not_invent_live_search(kind, parameters, words):
    audit = {'context_scope_counts':{'prepared_context':20,'scope_context':4}}
    text = scope_summary(parameters, audit)
    for word in words:
        assert word in text
    assert '20 prepared public comparisons' in text and '4 within context scope' in text


def test_saved_provider_search_and_cohort_denominators_are_distinct():
    text = scope_summary({'lin_level':5,'display_lin_level':7},
                         {'context_scope_counts':{'prepared_context':3,'scope_context':3}},
                         provider={'policy':{'lin_level':5},'pools':[{'discovered':100}], 'accessible_matches':8},
                         records=[row('q'),row('a'),row('b'),row('c')],
                         analysis={'cohorts':[{'sample_ids':['q','a','b']}]})
    for phrase in ['Saved search: species/ST → LIN5','100 ST candidate records','8 prefix-matching source records', '3 prepared public comparisons / within context scope / public comparisons in this block','2 in NJ comparable cohorts']:
        assert phrase in text


@pytest.mark.parametrize('direct', [False, True])
def test_rendered_nearest_annotations_keep_accessions_and_mode_units(tmp_path, direct):
    records = [row('q'),row('c', '0,0,105,6,0,0,1,5,1,0')]
    neighbour = dict(query_id='q', context_id='c', status='matched', distance=0,
                     allele_differences=0, shared_called_loci=628, snp_differences=0, shared_called_sites=1000)
    analysis = {'partition':{'kind':'cglin'},'display_group_level':7,
                'sample_labels':{'q':'SAMEA111','c':'SRR222'},'nearest_neighbours':[neighbour]}
    if direct:
        analysis['distance_method'] = 'direct-dna'
    before = deepcopy(analysis)
    text = write_profile_report(analysis, directory=tmp_path,
            provenance={'relatedness_parameters':{'lin_level':5,'display_lin_level':7},
                        'relatedness_audit':{'context_scope_counts':{'prepared_context':1,'scope_context':1}}},
            relatedness_records=records).read_text()
    assert 'Different LIN7 prefixes · 6 shared levels' in text and 'LIN relationship' in text
    assert 'Shared resolved LIN depth' not in text and 'Reported prefix relation' not in text
    assert 'SAMEA111' in text and 'SRR222' in text
    assert ('SNP differences / compared sites' if direct else 'Allele differences / shared loci') in text
    assert 'Comparison scope.' in text and 'Analysis: LIN5' in text and 'Colours: LIN7' in text
    assert analysis == before


def test_stage_default_scope_is_present_without_explicit_context_level(tmp_path):
    audit = {'unresolved_inputs':[], 'available_level_counts':[],
             'context_scope_counts':{'prepared_context':50,'scope_context':50}}
    path = write_partition_index([], audit, records_count=51, lin_level=5, hiercc_level=None,
                                 status='complete',output=tmp_path,display_lin_level=8)
    text = path.read_text()
    assert 'Analysis: LIN5' in text and 'Colours: LIN8' in text
    assert '50 prepared public comparisons' in text and 'tight pairwise cluster' in text


def test_refresh_scope_keeps_full_block_denominator_beside_smaller_cohort():
    text = scope_summary({'lin_level':5}, {}, records=[row('q'),row('a')],
                         analysis={'cohorts':[{'sample_ids':['q','a']}]}, block_context_count=3)
    assert '3 public comparisons in this block' in text and '1 in NJ comparable cohorts' in text


def test_equal_denominators_coalesce_without_hiding_different_counts():
    provider = {'policy':{'lin_level':8},'pools':[{'discovered':13421}],
                'accessible_matches':1298,'context_deduplication':{'counting_unit':'accession units'}}
    records = [row('q'),row('a'),row('b')]
    text = scope_summary({'lin_level':5,'context_lin_level':8,'display_lin_level':8},
                         {'context_scope_counts':{'prepared_context':2,'scope_context':2}},
                         provider=provider,records=records,
                         analysis={'cohorts':[{'sample_ids':['q','a','b']}]})
    assert '13,421 ST candidate records → 1,298 prefix-matching source records → 2 unique public comparisons.' in text
    assert 'All 2 are included and comparable in this block.' in text
    assert '2 within context scope' not in text and '2 in NJ comparable cohorts' not in text
    different = scope_summary({'lin_level':5}, {'context_scope_counts':{'prepared_context':3,'scope_context':2}},
                              records=records,analysis={'cohorts':[{'sample_ids':['q','a']}]})
    assert '3 prepared public comparisons' in different
    assert '2 within context scope / public comparisons in this block' in different
    assert '1 in NJ comparable cohorts' in different
