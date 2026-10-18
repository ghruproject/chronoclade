"""Legacy presentation must not change saved numerical or grouping evidence."""

from copy import deepcopy
from io import StringIO
import json

from Bio import Phylo
import pytest

from chronoclade.profile_report import _root_to_tip
from chronoclade.temporal_diagnostics import diagnostic_display_units
from chronoclade.report_components.group_tree import write_group_tree


@pytest.mark.parametrize('panel', ['2513 declared loci','629 observed export loci','600 fixed callable loci'])
def test_legacy_allele_units_are_concise_without_mutating_saved_diagnostic(panel):
    detailed=('allele-distance units on '+panel if 'fixed' in panel else
              'scheme-equivalent allele-distance units on '+panel)
    diagnostic={'cohort_id':'cohort_1','distance_units':detailed,'slope_units':detailed+' per year',
                'normalization_loci':int(panel.split()[0]),'fraction_slope':0.001,
                'midpoint_date_slope':2.513,'slope':2.513,'points':[{'root_to_tip':5.0}]}
    before=deepcopy(diagnostic)
    assert diagnostic_display_units(diagnostic)=='Estimated allele differences'
    assert diagnostic_display_units(diagnostic,slope=True)=='Estimated allele differences/year'
    html=_root_to_tip([diagnostic])
    assert 'Estimated allele differences/year' in html and detailed not in html
    assert diagnostic==before


def test_direct_dna_and_unrecognized_units_remain_distinct_and_escaped():
    for units in ('Nucleotide differences/year','Substitutions/site/year','custom <units>/year'):
        diagnostic={'slope_units':units,'distance_units':units,'midpoint_date_slope':0.1}
        assert diagnostic_display_units(diagnostic,slope=True)==units
        assert 'Estimated allele' not in _root_to_tip([diagnostic])
    assert 'custom &lt;units&gt;/year' in _root_to_tip([{'slope_units':'custom <units>/year'}])
    assert diagnostic_display_units({'slope_units':'legacy','display_slope_units':'Explicit unit'},slope=True)=='Explicit unit'


@pytest.mark.parametrize('groups,basis,threshold,expected,forbidden',[
    ({'a':'Lineage unassigned','b':'Lineage unassigned'},'unassigned_lineage',None,
     'does not infer genetic clusters','Groups use the declared'),
    ({'a':'cgLIN 8: 0.0.1','b':'cgLIN 8: unassigned'},'reported_cglin',None,
     'reported cgLIN prefixes','Groups use the declared'),
    ({'a':'HierCC HC200: 1','b':'HierCC HC200: unassigned'},'reported_hiercc',None,
     'reported HierCC assignments','Groups use the declared'),
    ({'a':'cohort_1_group_1','b':'cohort_1_group_2'},'distance_threshold',0.02,
     'mismatch fraction 0.02','reported cgLIN'),
    ({'a':'label-a','b':'label-b'},None,None,
     'no grouping method is declared','direct profile distances'),
])
def test_tree_caption_uses_actual_grouping_basis_and_preserves_data(tmp_path,groups,basis,threshold,expected,forbidden):
    tree=Phylo.read(StringIO('(a:-0.1,b:0.2);'),'newick')
    original=[(node.name,node.branch_length) for node in tree.find_clades()]
    records=[{'sample_id':'a','role':'input'},{'sample_id':'b','role':'context'}]
    before=deepcopy(groups)
    path=write_group_tree(tree,records,{},groups,tmp_path/'tree.html',
                          group_basis=basis,group_threshold=threshold)
    html=open(path).read()
    assert expected in html and forbidden not in html
    assert [(node.name,node.branch_length) for node in tree.find_clades()]==original
    assert groups==before
    raw=html.split('<script type="application/json" id="tree-data">')[1].split('</script>')[0]
    data=json.loads(raw)
    assert {row['id']:row['group_id'] for row in data['tips']}==groups
