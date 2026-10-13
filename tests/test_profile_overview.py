import json

from chronoclade.report_components.profile_overview import (
    genetic_group_summary,
    grouping_rule,
    input_scope,
    nearest_summary,
    selection_handoff,
)


def test_scope_uses_explicit_original_and_block_counts():
    assert '7 of those inputs' in input_scope({}, {'input_partition':
        {'total_inputs': 10, 'block_inputs': 7, 'blocks': 4}})
    assert 'remaining 3' in input_scope({}, {'input_partition':
        {'total_inputs': 10, 'block_inputs': 7, 'blocks': 4}})
    assert input_scope({}, {}) == ''
    assert 'All 7' in input_scope({}, {'input_partition':
        {'total_inputs': 7, 'block_inputs': 7, 'blocks': 1}})


def test_nearest_summary_retains_minimum_raw_ties_and_denominator_range():
    rows = [{'query_id': 'q', 'context_id': ident, 'status': 'matched',
             'allele_differences': difference, 'shared_called_loci': loci,
             'country': country} for ident, difference, loci, country in
            [('PW_a', 1, 600, 'Greece'), ('PW_b', 1, 620, 'Italy'), ('PW_c', 2, 629, 'France')]]
    html = nearest_summary({'nearest_neighbours': rows,
        'sample_labels': {'PW_a': 'ERR123', 'PW_b': 'ERR124', 'PW_c': 'ERR125'}})
    assert '1 allele difference; 600–620 compared loci' in html
    assert '2 closest public records' in html
    assert 'ERR123' in html and 'ERR124' in html
    assert 'ERR125' not in html and 'France' not in html and 'PW_' not in html
    for row in rows:
        row['snp_differences'] = row['allele_differences'] * 3
        row['shared_called_sites'] = row['shared_called_loci'] * 100
    dna = nearest_summary({'distance_method': 'direct-dna', 'nearest_neighbours': rows})
    assert '3 SNP differences; 60000–62000 compared sites' in dna


def test_groups_join_analytical_ids_without_merging_accession_aliases(tmp_path):
    data = {'sample_metadata': [
        {'sample_id': 'q', 'role': 'input', 'country': 'Greece', 'host': 'Human', 'collection_date': '2019'},
        {'sample_id': 'PW_a', 'role': 'context', 'country': 'Italy', 'host': 'human', 'collection_date': '2021'},
        {'sample_id': 'PW_b', 'role': 'context', 'country': '<script>', 'host': None},
    ], 'sample_labels': {'PW_a': 'ERR123', 'PW_b': 'ERR123'},
        'genetic_groups': [{'group_id': 'cohort_1_group_1', 'sample_ids': ['q', 'PW_a', 'PW_b']}]}
    html = genetic_group_summary(tmp_path, data)
    assert 'Group 1' in html and 'Human 2' in html and 'Unknown 1' in html
    assert '2019–2021; 1 undated' in html
    assert '<td>1</td><td>2</td><td>3</td>' in html
    assert '&lt;script&gt; 1' in html and '<script>' not in html


def test_handoff_counts_union_without_claiming_downstream_analysis(tmp_path):
    bundle = tmp_path / 'selection_bundle'
    (bundle / 'selections').mkdir(parents=True)
    (bundle / 'ensemble.json').write_text(json.dumps({'selections': [{'path': 'selections/a.json'}]}))
    (bundle / 'selections/a.json').write_text(json.dumps({'query_ids': ['q'],
        'available_context_ids': ['a', 'b', 'c'], 'selected_context_ids': ['a', 'b']}))
    html = selection_handoff(tmp_path, {'coverage': {'total_inputs': 1}})
    assert '4 genomes in this block → 3 selected' in html
    assert 'all 1 inputs + 2 public comparisons' in html
    assert 'href="#subsamples"' in html
    assert selection_handoff(tmp_path / 'missing', {}) == ''


def test_group_rule_exposes_saved_threshold_and_correct_denominator():
    html = grouping_rule({'distance_threshold': .02, 'genetic_groups': [{}, {}, {}]})
    assert '3 coloured genetic groups' in html
    assert 'no more than 2% of jointly called cgMLST loci' in html
    assert 'not inferred from the MDS plot' in html and 'not Pathogenwatch LIN codes' in html
    dna = grouping_rule({'distance_method': 'direct-dna', 'distance_threshold': .000005,
                         'genetic_groups': [{}]})
    assert '0.0005% of jointly compared nucleotide sites' in dna
    assert grouping_rule({}) == ''
