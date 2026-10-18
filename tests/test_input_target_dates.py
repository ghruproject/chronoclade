"""Explicit input ancestors reuse saved dated-tree evidence without inference."""

from copy import deepcopy
from io import StringIO

from Bio import Phylo

from chronoclade.stage_artifacts import reference
from chronoclade.time_stage import comparison_evidence as evidence


def fixture(tmp_path):
    tree = Phylo.read(StringIO('((a:1,b:1,c:1)INPUT_ANCESTOR:1,(d:1,e:1)OTHER:1)ROOT;'), 'newick')
    path = tmp_path/'dated.nexus'
    Phylo.write(tree, path, 'nexus')
    dates = tmp_path/'node_dates.csv'
    dates.write_text('node,numeric_date,lower_bound,upper_bound\nROOT,2000,1999,2001\nINPUT_ANCESTOR,2012,2010,2014\n')
    result = {'dated_tree':reference(tmp_path,path), 'clock_confidence':{'interval_level':'90% max-posterior region'},
              'artifacts':{'node_dates.csv':reference(tmp_path,dates)}}
    return tmp_path/'time.json',result


def test_distinct_input_target_has_own_saved_ancestor_and_preserves_common_numbers(tmp_path):
    path,result=fixture(tmp_path)
    common=evidence.anchor_date_evidence(path,result,['a','b','d'])
    before=deepcopy(common)
    target=evidence.input_target_date_evidence(path,result,['a','b'],['a','b','d'],common)
    assert common==before and common['estimate']==2000
    assert target['estimate']==2012 and (target['lower_bound'],target['upper_bound'])==(2010,2014)
    assert target['target_sample_ids']==['a','b'] and target['target_count']==2
    assert target['descendant_sample_ids']==['a','b','c']
    assert target['extra_descendant_sample_ids']==['c']
    assert target['same_target_as_common_anchors'] is False


def test_same_targets_reuse_common_without_a_second_read(tmp_path,monkeypatch):
    path,result=fixture(tmp_path)
    common=evidence.anchor_date_evidence(path,result,['a','b'])
    before=deepcopy(common)
    monkeypatch.setattr(evidence,'anchor_date_evidence',lambda *a,**k: (_ for _ in ()).throw(AssertionError('duplicate date read')))
    target=evidence.input_target_date_evidence(path,result,['b','a'],['a','b'],common)
    assert target['same_target_as_common_anchors'] is True
    assert target['estimate']==common['estimate'] and common==before


def test_unsupported_and_missing_dates_remain_unavailable(tmp_path):
    path,result=fixture(tmp_path)
    result['dated_tree']=None
    common=evidence.anchor_date_evidence(path,result,['a','b','d'])
    unavailable=evidence.input_target_date_evidence(path,result,['a','b'],['a','b','d'],common)
    assert unavailable['estimate'] is None and unavailable['status']=='unavailable'
    path,result=fixture(tmp_path)
    result['artifacts']={}
    common=evidence.anchor_date_evidence(path,result,['a','b','d'])
    missing=evidence.input_target_date_evidence(path,result,['a','b'],['a','b','d'],common)
    assert missing['estimate'] is None and missing['reason']=='No saved node-date table'
    assert missing['extra_descendant_sample_ids']==['c']


def test_missing_bounds_are_not_invented(tmp_path):
    path,result=fixture(tmp_path)
    dates=tmp_path/'node_dates.csv'
    dates.write_text('node,numeric_date,lower_bound,upper_bound\nROOT,2000,,\nINPUT_ANCESTOR,2012,,\n')
    result['artifacts']['node_dates.csv']=reference(tmp_path,dates)
    common=evidence.anchor_date_evidence(path,result,['a','b','d'])
    target=evidence.input_target_date_evidence(path,result,['a','b'],['a','b','d'],common)
    assert target['estimate']==2012 and target['lower_bound'] is target['upper_bound'] is None
    assert target['interval_available'] is False


def render_result(tmp_path, *, same=False, dated=True):
    from chronoclade.time_stage.comparison import _spread
    from chronoclade.time_stage.comparison_report import write_comparison_report
    path,result=fixture(tmp_path)
    anchors=['a','b'] if same else ['a','b','d']
    common=evidence.anchor_date_evidence(path,result,anchors)
    target=evidence.input_target_date_evidence(path,result,['a','b'],anchors,common)
    row={'selection_id':'one','status':'dated' if dated else 'unsupported',
         'common_anchor_mrca':common,'common_anchor_mrca_date':common['estimate'],
         'input_target_mrca':target,'input_target_mrca_date':target['estimate']}
    comparison={'runs':[row],'selection_coverage':{'interpretation':'fixture','cumulative_unique_genomes':[],
                                                 'pairwise_overlap':[],'runs':[]},
                'common_anchor_ids':anchors,'target_ids':['a','b'],'input_target_ids':['a','b'],
                'dated_runs':int(dated),'sensitivity':{'common_anchor_mrca_date':_spread([common['estimate']]),
                                                      'input_target_mrca_date':_spread([target['estimate']])},
                'target_pair_agreement':[],'shared_anchor_splits':{'interpretation':'fixture','splits':[]},
                'interpretation':'Between-run range is not a confidence interval.'}
    write_comparison_report(tmp_path,comparison)
    return (tmp_path/'report.html').read_text()


def test_report_distinguishes_input_target_without_adding_a_figure(tmp_path):
    html=render_result(tmp_path)
    assert 'shared-anchor ancestor uses 3 genomes' in html
    assert 'including inputs and any shared required public comparisons' in html
    assert 'Input-genome ancestor (2 genomes)' in html
    assert '2012' in html and '90% region lower' in html
    assert 'Additional descendants' in html
    assert '<img' not in html


def test_report_deduplicates_identical_targets_and_hides_unsupported_panel(tmp_path):
    html=render_result(tmp_path,same=True)
    assert 'These are the same target' in html
    assert 'Input-genome ancestor (' not in html
    html=render_result(tmp_path,dated=False)
    assert 'Input-genome ancestor (' not in html


def test_empty_input_target_has_consistent_unavailable_interval_fields(tmp_path):
    path,result=fixture(tmp_path)
    common=evidence.anchor_date_evidence(path,result,['a','b','d'])
    empty=evidence.input_target_date_evidence(path,result,[],['a','b','d'],common)
    assert empty['status']=='unavailable' and empty['estimate'] is None
    assert empty['lower_bound'] is empty['upper_bound'] is None
    assert empty['interval_level']=='90% max-posterior region'
    assert empty['target_sample_ids']==[] and empty['target_count']==0
