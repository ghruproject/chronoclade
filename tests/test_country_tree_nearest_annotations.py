from io import StringIO
from xml.etree import ElementTree

from Bio import Phylo
import pytest

from chronoclade.location_network.tree import draw_country_tree
from chronoclade.profile_network import build_profile_network


@pytest.mark.parametrize('nearest,shown,marked', [
    ([], None, 0),
    (['absent'], None, 0),
    (['c'], ['q'], 0),
    (['c'], None, 1),
])
def test_nearest_legend_requires_visible_marked_tips(tmp_path, nearest, shown, marked):
    rows = [dict(sample_id='q', origin='query', country='Greece'),
            dict(sample_id='c', origin='context', country='Italy')]
    tree = Phylo.read(StringIO('(q:0.1,c:0.1);'), 'newick')
    network = build_profile_network(tree, rows)
    path = tmp_path / 'country.svg'
    label = 'Mandatory cgMLST-guide neighbours (ties retained)'
    draw_country_tree(network, rows, path, nearest_ids=nearest, display_ids=shown,
                      nearest_label=label)
    text = path.read_text()
    assert (label in text) == bool(marked)
    assert 'Outlined tips: input genomes' in text
    assert 'Hollow ancestors: multiple optimal countries' in text
    # The orange star appears once per marked tip, plus its conditional legend.
    svg = ElementTree.fromstring(text)
    stars = [element for element in svg.iter()
             if element.tag.endswith('}use')
             and 'stroke: #e66101' in element.attrib.get('style', '')]
    assert len(stars) == marked + bool(marked)
