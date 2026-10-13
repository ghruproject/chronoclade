import io

import pytest
from Bio import Phylo

from chronoclade.profile_rooting import project_outgroup_root


def tree(text):
    return Phylo.read(io.StringIO(text), 'newick')


def test_external_attachment_preserves_original_ingroup_distances_and_does_not_mutate():
    original = tree('(a:1,(b:2,(c:3,d:4):5):6);')
    augmented = tree('((a:1,b:2):2,(c:3,d:4):3,out:7);')
    before = {(a, b): original.distance(a, b) for a in 'abcd' for b in 'abcd'}
    rooted, audit = project_outgroup_root(original, augmented, 'out')
    assert sorted(len(child.get_terminals()) for child in rooted.root.clades) == [2, 2]
    for pair, distance in before.items():
        assert rooted.distance(*pair) == pytest.approx(distance)
        assert original.distance(*pair) == pytest.approx(distance)
    assert audit['dates_used_to_choose_root'] is False
    assert audit['ingroup_split_sizes'] == [2, 2]


def test_root_projection_rejects_incompatible_basal_splits():
    with pytest.raises(ValueError, match='absent'):
        project_outgroup_root(tree('((a:1,b:1):2,(c:1,d:1):2);'),
                              tree('((a:1,c:1):2,(b:1,d:1):2,out:3);'), 'out')


def test_root_projection_rejects_negative_basal_length_and_internal_outgroup():
    with pytest.raises(ValueError, match='negative'):
        project_outgroup_root(tree('((a:1,b:1):1,c:1);'),
                              tree('((a:1,b:1):-1,c:2,out:3);'), 'out')
    with pytest.raises(ValueError, match='external outgroup'):
        project_outgroup_root(tree('(a:1,b:1);'), tree('(a:1,b:1);'), 'a')
