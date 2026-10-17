"""Exact panel and unit contracts for matched clock comparisons."""
from io import StringIO
from types import SimpleNamespace

import numpy as np
import pytest
from Bio import Phylo

from chronoclade.cgmlst.clock_comparison import (
    _check_binding, _fit, compare_panel, exact_panel_nj, matched_dates,
)
from chronoclade.errors import WorkflowError


def test_fresh_nj_uses_exact_ordered_submatrix():
    evidence = SimpleNamespace(sample_ids=('unused', 'a', 'b', 'c'), arrays={
        'distance': np.array([[0, .9, .8, .7], [.9, 0, .1, .2], [.8, .1, 0, .3], [.7, .2, .3, 0]]),
        'scheme_loci': np.full((4, 4), 100),
    })
    tree, loci = exact_panel_nj(evidence, ['c', 'a', 'b'])
    assert loci == 100
    assert {t.name for t in tree.get_terminals()} == {'a', 'b', 'c'}
    assert tree.distance('a', 'b') == pytest.approx(.1)
    assert tree.distance('a', 'c') == pytest.approx(.2)
    with pytest.raises(WorkflowError, match='unique'):
        exact_panel_nj(evidence, ['a', 'a', 'b'])
    evidence.arrays['distance'][1, 2] = np.nan
    with pytest.raises(WorkflowError, match='unavailable'):
        exact_panel_nj(evidence, ['a', 'b', 'c'])


def test_date_panel_matches_before_any_fit_and_keeps_original_dates():
    rows = [{'sample_id': name, 'collection_date': year} for name, year in
        [('a', '1899'), ('b', '1900'), ('c', '2000'), ('d', '2001'), ('e', 'unknown')]]
    eligible, audit = matched_dates(rows, ['a', 'b', 'c', 'd', 'e'])
    assert [r['sample_id'] for r in eligible] == ['b', 'c', 'd']
    assert audit['date_filter']['minimum_collection_year'] == 1900
    assert rows[0]['collection_date'] == '1899'
    with pytest.raises(WorkflowError, match='ordered'):
        matched_dates(rows, ['b', 'a', 'c', 'd', 'e'])


def test_distinct_units_and_diagnostic_copy(tmp_path):
    tree = Phylo.read(StringIO('(a:-.1,b:.2,c:.4);'), 'newick')
    rows = [{'sample_id': name, 'collection_date': year} for name, year in
        [('a', '2000'), ('b', '2001'), ('c', '2002')]]
    before = tree.find_any(name='a').branch_length
    result = compare_panel(tree, rows, tmp_path, 'cgmlst', 'cgMLST mismatch fraction', 100)
    assert tree.find_any(name='a').branch_length == before
    assert result['branch_length_policy']['negative_edge_count'] == 1
    assert result['branch_length_policy']['negative_edge_total_magnitude'] == pytest.approx(.1)
    assert result['clock_compatible']['normalization_loci'] == 100
    assert 'scheme-equivalent' in result['clock_compatible']['slope_units']
    snp = _fit(tree, rows, 'substitutions/site')
    assert snp['slope_units'] == 'substitutions/site per year'
    assert 'normalization_loci' not in snp
    incomplete = _fit(tree, rows, 'cgMLST mismatch fraction', 100, catalogue_complete=False)
    assert 'observed export loci' in incomplete['slope_units']


def test_exact_source_and_metadata_binding():
    entry = {'selection_id': 'selection-001', 'sha256': 'selection-sha'}
    rows = [{'sample_id': 'a', 'collection_date': '2000'}]
    result = {'selection_id': 'selection-001', 'selected_sample_ids': ['a'],
        'source': {'selection': {'sha256': 'selection-sha'}, 'dataset': {'sha256': 'dataset-sha'}}}
    _check_binding(result, entry, ['a'], rows, {'a': rows[0]}, 'dataset-sha')
    with pytest.raises(WorkflowError, match='prepared dataset'):
        _check_binding(result, entry, ['a'], [{'sample_id': 'a', 'collection_date': '2001'}],
                       {'a': rows[0]}, 'dataset-sha')
    result['source']['selection']['sha256'] = 'different-sha'
    with pytest.raises(WorkflowError, match='exact saved selection'):
        _check_binding(result, entry, ['a'], rows, {'a': rows[0]}, 'dataset-sha')
