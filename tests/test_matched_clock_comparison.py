"""Exact panel and unit contracts for matched clock comparisons."""
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from Bio import Phylo
import chronoclade

from chronoclade.cgmlst.clock_comparison import (
    _check_binding, _comparison_computed_with, _fit, _report_slope_units, compare_panel, exact_panel_nj, matched_dates,
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


def _clock_profile(sample_id, calls, database=''):
    loci = [f'locus_{index}' for index in range(len(calls))]
    return {
        'sample_id': sample_id,
        'role': 'context',
        'origin': 'context',
        'species': 'Klebsiella pneumoniae',
        'cgmlst_scheme': 'scheme',
        'cgmlst_scheme_version': '1',
        'cgmlst_loci': loci,
        'cgmlst_profile': dict(zip(loci, calls)),
        'cgmlst_database_sha256': database,
    }


def _dense_from_binary(binary):
    from chronoclade.selections import DistanceEvidence

    rows = tuple({
        'sample_id_1': pair['sample_a'],
        'sample_id_2': pair['sample_b'],
        'distance': str(pair['distance']),
        'allele_differences': str(pair['allele_differences']),
        'shared_called_loci': str(pair['shared_called_loci']),
        'call_overlap': str(pair['call_overlap']),
    } for pair in binary.iter_pairs())
    return DistanceEvidence(
        rows,
        method='alleles',
        distance_definition='fraction of mismatching jointly called cgMLST loci',
    )


def test_dense_and_binary_evidence_build_the_same_selected_nj(tmp_path):
    from chronoclade.cgmlst.scale_distance import write_distances
    from chronoclade.matrix_distances import BinaryDistanceEvidence

    rows = [
        _clock_profile('a', [1, 1, 1, 1]),
        _clock_profile('b', [1, 1, 2, 1]),
        _clock_profile('c', [2, 1, 2, 1]),
    ]
    written, _, _ = write_distances(rows, tmp_path / 'binary', min_overlap=1)
    binary = BinaryDistanceEvidence(written.manifest_path)
    dense = _dense_from_binary(binary)
    selected = ['c', 'a', 'b']

    binary_tree, binary_loci = exact_panel_nj(binary, selected, expected_loci=4)
    dense_tree, dense_loci = exact_panel_nj(dense, selected, expected_loci=4)

    assert binary_loci == dense_loci == 4
    for left in selected:
        for right in selected:
            assert dense_tree.distance(left, right) == pytest.approx(binary_tree.distance(left, right))


def test_dense_and_binary_evidence_reject_the_same_unavailable_selected_pair(tmp_path):
    from chronoclade.cgmlst.scale_distance import write_distances
    from chronoclade.matrix_distances import BinaryDistanceEvidence

    rows = [
        _clock_profile('a', [1, 1, 1, 1], database='a' * 64),
        _clock_profile('b', [1, 1, 2, 1], database='a' * 64),
        _clock_profile('c', [2, 1, 2, 1], database='b' * 64),
        _clock_profile('d', [2, 2, 2, 1], database='b' * 64),
    ]
    written, _, _ = write_distances(rows, tmp_path / 'binary', min_overlap=1)
    binary = BinaryDistanceEvidence(written.manifest_path)
    dense = _dense_from_binary(binary)

    for evidence in (binary, dense):
        with pytest.raises(WorkflowError, match='unavailable distance pairs'):
            exact_panel_nj(evidence, ['a', 'b', 'c'], expected_loci=4)


def test_dense_evidence_validates_selected_identity_and_catalogue_loci(tmp_path):
    from chronoclade.cgmlst.scale_distance import write_distances
    from chronoclade.matrix_distances import BinaryDistanceEvidence

    rows = [_clock_profile(name, calls) for name, calls in [
        ('a', [1, 1, 1, 1]), ('b', [1, 1, 2, 1]), ('c', [2, 1, 2, 1]),
    ]]
    written, _, _ = write_distances(rows, tmp_path / 'binary', min_overlap=1)
    dense = _dense_from_binary(BinaryDistanceEvidence(written.manifest_path))

    binary = BinaryDistanceEvidence(written.manifest_path)
    for evidence in (binary, dense):
        with pytest.raises(WorkflowError, match='present in full-pool distances'):
            exact_panel_nj(evidence, ['a', 'b', 'absent'], expected_loci=4)
        with pytest.raises(WorkflowError):
            exact_panel_nj(evidence, ['a', 'b', 'c'], expected_loci=5)
    with pytest.raises(WorkflowError, match='call_overlap disagrees'):
        exact_panel_nj(dense, ['a', 'b', 'c'], expected_loci=5)


def test_dense_evidence_rejects_conflicting_duplicate_pair_and_inconsistent_distance(tmp_path):
    from dataclasses import replace
    from chronoclade.cgmlst.scale_distance import write_distances
    from chronoclade.matrix_distances import BinaryDistanceEvidence

    rows = [_clock_profile(name, calls) for name, calls in [
        ('a', [1, 1, 1, 1]), ('b', [1, 1, 2, 1]), ('c', [2, 1, 2, 1]),
    ]]
    written, _, _ = write_distances(rows, tmp_path / 'binary', min_overlap=1)
    dense = _dense_from_binary(BinaryDistanceEvidence(written.manifest_path))
    duplicate = replace(dense, records=(*dense.records, dense.records[0]))
    with pytest.raises(WorkflowError, match='Duplicate dense distance evidence'):
        exact_panel_nj(duplicate, ['a', 'b', 'c'], expected_loci=4)
    changed = dict(dense.records[0], distance='0.75')
    inconsistent = replace(dense, records=(changed, *dense.records[1:]))
    with pytest.raises(WorkflowError, match='disagree with mismatch/callable counts'):
        exact_panel_nj(inconsistent, ['a', 'b', 'c'], expected_loci=4)


def test_clock_evidence_loader_accepts_dense_and_binary_schemas(tmp_path):
    import json
    from chronoclade.cgmlst.clock_comparison import _read_clock_distance_evidence
    from chronoclade.cgmlst.scale_distance import write_distances
    from chronoclade.matrix_distances import BinaryDistanceEvidence

    rows = [_clock_profile(name, calls) for name, calls in [
        ('a', [1, 1, 1, 1]), ('b', [1, 1, 2, 1]), ('c', [2, 1, 2, 1]),
    ]]
    written, _, _ = write_distances(rows, tmp_path / 'binary', min_overlap=1)
    binary = BinaryDistanceEvidence(written.manifest_path)
    source = {
        'method': 'alleles',
        'distance_definition': 'fraction of mismatching jointly called cgMLST loci',
    }
    dense_path = tmp_path / 'dense-distance-evidence.json'
    dense_path.write_text(json.dumps({
        **source,
        'records': [dict(row) for row in _dense_from_binary(binary).records],
    }))
    dense_loaded = _read_clock_distance_evidence(dense_path, source)
    binary_loaded = _read_clock_distance_evidence(
        written.manifest_path,
        {'method': 'cgmlst', 'distance_definition': source['distance_definition']},
    )
    assert len(dense_loaded.records) == 3
    assert isinstance(binary_loaded, BinaryDistanceEvidence)


def test_clock_evidence_loader_rejects_unknown_declared_schema(tmp_path):
    import json
    from chronoclade.cgmlst.clock_comparison import _read_clock_distance_evidence

    source = {
        'method': 'alleles',
        'distance_definition': 'fraction of mismatching jointly called cgMLST loci',
    }
    path = tmp_path / 'unknown-schema.json'
    path.write_text(json.dumps({
        **source,
        'schema': 'vendor.custom-distance-v1',
        'records': [{'sample_id_1': 'a', 'sample_id_2': 'b'}],
    }))
    with pytest.raises(WorkflowError, match='Unsupported distance evidence schema'):
        _read_clock_distance_evidence(path, source)


@pytest.mark.parametrize('field', ['allele_differences', 'shared_called_loci'])
def test_dense_distance_evidence_rejects_boolean_count_values(tmp_path, field):
    from dataclasses import replace

    from chronoclade.matrix_distances import BinaryDistanceEvidence
    from chronoclade.cgmlst.scale_distance import write_distances

    rows = [_clock_profile(name, calls) for name, calls in [
        ('a', [1, 1, 1, 1]), ('b', [1, 1, 2, 1]), ('c', [2, 1, 2, 1]),
    ]]
    written, _, _ = write_distances(rows, tmp_path / field, min_overlap=1)
    dense = _dense_from_binary(BinaryDistanceEvidence(written.manifest_path))
    changed_row = dict(dense.records[0], **{field: True})
    changed = replace(dense, records=(changed_row, *dense.records[1:]))
    with pytest.raises(WorkflowError, match=f'invalid {field}'):
        exact_panel_nj(changed, ['a', 'b', 'c'], expected_loci=4)


def test_matched_clock_report_uses_display_units_and_records_execution_identity():
    fitted = {
        'slope_units': 'Estimated allele differences/year on 629 observed loci',
        'display_slope_units': 'Estimated allele differences/year',
    }
    assert _report_slope_units(fitted) == 'Estimated allele differences/year'
    assert _report_slope_units({'slope_units': 'legacy units'}) == 'legacy units'
    computed = _comparison_computed_with()
    assert computed['component'] == 'matched-clock-comparison'
    assert computed['method_revision'] == 'matched-root-diagnostics-v1'
    package_root = Path(chronoclade.__file__).resolve().parent
    assert computed['software']['source_file_count'] == sum(1 for _ in package_root.rglob('*.py'))
    from chronoclade.execution_provenance import execution_provenance

    expected = execution_provenance(computed['component'], computed['method_revision'])
    assert computed['software'] == expected['software']
