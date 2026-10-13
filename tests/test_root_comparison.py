import json
from chronoclade.report_components.root_comparison import add_root_comparisons
from chronoclade.profile_report import write_profile_report


def test_clock_comparison_keeps_primary_tree_network_and_embeds_evidence(tmp_path):
    path=tmp_path/'tree.nwk'
    path.write_text('((a:1,b:2):1,(c:3,d:4):1);')
    before=path.read_bytes()
    rows=[{'sample_id':s,'collection_date':str(2010+i)} for i,s in enumerate('abcd')]
    analysis={'distance_units':'mismatch fraction','cohorts':[{'cohort_id':'one',
      'sample_ids':list('abcd'),'tree_path':'tree.nwk'}],
      'root_to_tip':[{'cohort_id':'one','root':'outgroup','normalization_loci':629}],
      'location_network':[{'cohort_id':'one','method':'unchanged'}]}
    add_root_comparisons(analysis,rows,tmp_path)
    result=analysis['root_comparisons'][0]
    assert path.read_bytes()==before
    assert analysis['location_network']==[{'cohort_id':'one','method':'unchanged'}]
    assert len(result['clock_compatible']['points'])==len(result['baseline']['points'])==4
    assert result['clock_compatible']['residual_mse']<=result['baseline']['residual_mse']+1e-10
    assert result['slope_change']==result['clock_compatible']['slope']-result['baseline']['slope']
    assert result['audit']['pairwise_distances_preserved']
    assert json.loads((tmp_path/'one_root_comparison.json').read_text())['audit']['status']=='rooted'
    html=write_profile_report(analysis,directory=tmp_path).read_text()
    assert 'Clock-compatible (date-optimised)' in html
    assert 'does not establish temporal signal' in html
    assert 'data:image/svg+xml' in html


def test_comparison_reports_insufficient_dates(tmp_path):
    (tmp_path/'tree.nwk').write_text('(a:1,b:1,c:1);')
    analysis={'cohorts':[{'cohort_id':'one','sample_ids':list('abc'),'tree_path':'tree.nwk'}]}
    add_root_comparisons(analysis,[{'sample_id':s} for s in 'abc'],tmp_path)
    assert analysis['root_comparisons'][0]['audit']['status']=='skipped'
    html=write_profile_report(analysis,directory=tmp_path).read_text()
    assert 'Comparison unavailable' in html


def test_clock_comparison_has_real_cgmlst_workflow_outputs(tmp_path):
    from chronoclade.datasets import LocusCatalogue, from_profile_records, write_dataset
    from chronoclade.cgmlst.workflow import run_cgmlst
    rows=[]
    for index,name in enumerate(('q','c1','c2','c3')):
        rows.append(dict(sample_id=name,origin='local' if name=='q' else 'context',
            species='Klebsiella pneumoniae',mlst_st='39',collection_date=str(2010+index),
            country='Greece',cgmlst_scheme='test',cgmlst_scheme_version='1',
            cgmlst_profile={'a':str(index+1),'b':'1','c':str(1+index//2)},
            cgmlst_status='resolved',cglin_scheme='scgMLST629_S',
            cglin_scheme_version='1',cglin_raw='0,0,1,0,1,0,0',cglin_status='complete'))
    manifest=write_dataset(from_profile_records(rows,[LocusCatalogue('test','1',('a','b','c'),'test')]),tmp_path/'prepared')
    result=run_cgmlst(manifest,tmp_path/'run',bootstrap_replicates=0,context_size=2)
    data=json.loads(result.manifest_path.read_text())
    block=result.manifest_path.parent/data['blocks'][0]['block_id']
    analysis=json.loads((block/'profile_analysis.json').read_text())
    comparison=analysis['root_comparisons'][0]
    assert comparison['audit']['status']=='rooted'
    assert (block/comparison['tree_path']).is_file()
    assert (block/comparison['audit_path']).is_file()
    assert 'Clock-compatible (date-optimised)' in (block/'profile_report.html').read_text()
    assert comparison['baseline']['distance_units']==comparison['clock_compatible']['distance_units']
