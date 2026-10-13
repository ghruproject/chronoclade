import io

import numpy as np
import pytest
from Bio import Phylo

from chronoclade.report_components.tree_ordination import nj_distances, ordination_bases, pcoa


def test_patristic_distances_match_tree_and_requested_tip_order(tmp_path):
    tree = Phylo.read(io.StringIO('((a:1,b:2):3,(c:4,d:5):6);'), 'newick')
    names = ['d', 'a', 'c', 'b']
    matrix, audit = nj_distances(tree, names, output=tmp_path / 'distances.npy', batch_size=1)
    expected = np.array([[tree.distance(a, b) for b in names] for a in names])
    np.testing.assert_allclose(matrix, expected)
    np.testing.assert_allclose(np.load(tmp_path / 'distances.npy'), expected)
    assert audit['negative_branch_count'] == 0


def test_negative_branches_are_clipped_only_for_ordination():
    tree = Phylo.read(io.StringIO('((a:-2,b:1):3,c:4);'), 'newick')
    original = [node.branch_length for node in tree.find_clades()]
    matrix, audit = nj_distances(tree, ['a', 'b', 'c'])
    np.testing.assert_allclose(matrix, [[0, 1, 7], [1, 0, 8], [7, 8, 0]])
    assert audit['negative_branch_count'] == 1
    assert audit['negative_branch_total_magnitude'] == 2
    assert original == [node.branch_length for node in tree.find_clades()]


def test_mds_preserves_line_distances_and_exposes_both_bases():
    tree = Phylo.read(io.StringIO('(a:0,b:1,c:2);'), 'newick')
    direct = np.array([[0., 0.], [1., 0.], [2., 0.]])
    coordinates, audit, bases = ordination_bases(tree, [{'sample_id': x} for x in 'abc'], direct)
    assert coordinates.shape == (3, 2)
    assert np.isfinite(coordinates).all()
    assert audit['basis'] == 'nj_tree'
    np.testing.assert_array_equal(bases['direct']['coordinates'], direct)
    assert 'negative branches' in bases['nj_tree']['note']
    line = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    points, _ = pcoa(line)
    recovered = np.linalg.norm(points[:, None] - points[None, :], axis=2)
    np.testing.assert_allclose(recovered, line, atol=1e-7)


def test_tree_tip_mismatch_is_rejected():
    tree = Phylo.read(io.StringIO('(a:1,b:1);'), 'newick')
    with pytest.raises(ValueError, match='matching every tree tip'):
        nj_distances(tree, ['a', 'c'])
