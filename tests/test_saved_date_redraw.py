"""Saved-point rendering must not become another diagnostic calculation."""

from copy import deepcopy

import matplotlib.pyplot as plt
import pytest

from chronoclade.report_components.date_diagnostic import write_saved_date_figure, redraw_saved_date_figures


def diagnostic(**changes):
    return {'cohort_id':'cohort_1','distance_field':'root_to_tip',
            'distance_units':'scheme-equivalent allele-distance units on 2513 declared loci',
            'slope':2.0,'intercept':999999999999.0,'reference_year':2000.0,'centered_intercept':4.0,
            'interpretation':'Exploratory, arbitrary-root NJ diagnostic.',
            'points':[{'year_min':1999.0,'year_max':1999.0,'root_to_tip':2.0},
                      {'year_min':2001.0,'year_max':2001.0,'root_to_tip':6.0}],**changes}


def test_saved_centered_fit_and_units_are_drawn_without_mutation(tmp_path,monkeypatch):
    row=diagnostic()
    before=deepcopy(row)
    figure,axis=plt.subplots()
    original_close=plt.close
    monkeypatch.setattr(plt,'subplots',lambda **k:(figure,axis))
    monkeypatch.setattr(plt,'close',lambda fig:None)
    path=tmp_path/'date.svg'
    assert write_saved_date_figure(row,path)==str(path)
    line=next(line for line in axis.lines if line.get_label()=='Saved regression')
    assert list(line.get_ydata())==[2.0,6.0]  # Uncentered intercept is deliberately inconsistent.
    assert axis.get_ylabel()=='Estimated allele differences'
    assert 'arbitrary tip root' in axis.get_title()
    assert row==before
    original_close(figure)


@pytest.mark.parametrize('changes,title',[
    ({'interpretation':'saved profile diagnostic','cohort_root':'input-id'},'declared input tip root'),
    ({'interpretation':'root identity saved'},'saved root'),
])
def test_root_description_is_not_invented(tmp_path,changes,title):
    path=tmp_path/'figure.svg'
    write_saved_date_figure(diagnostic(**changes),path)
    assert title in path.read_text()


def test_legacy_fallback_and_dna_units_remain_valid(tmp_path,monkeypatch):
    row=diagnostic(distance_units='Nucleotide differences',centered_intercept=None,reference_year=None,intercept=-3996.0)
    figure,axis=plt.subplots()
    close=plt.close
    monkeypatch.setattr(plt,'subplots',lambda **k:(figure,axis))
    monkeypatch.setattr(plt,'close',lambda fig:None)
    write_saved_date_figure(row,tmp_path/'dna.svg')
    assert axis.get_ylabel()=='Nucleotide differences'
    assert list(next(line for line in axis.lines if line.get_label()=='Saved regression').get_ydata())==[2.0,6.0]
    close(figure)


def test_missing_stored_points_preserve_existing_asset_and_analysis(tmp_path):
    path=tmp_path/'old.svg'
    path.write_text('original saved asset')
    row=diagnostic()
    row.pop('points')
    analysis={'cohorts':[{'cohort_id':'cohort_1','root_to_tip_figure':'old.svg'}],'root_to_tip':[row]}
    before=deepcopy(analysis)
    audit=redraw_saved_date_figures(analysis,tmp_path)
    assert audit[0]['status']=='skipped_missing_or_invalid_saved_points'
    assert path.read_text()=='original saved asset' and analysis==before
    row['points']=[]
    write_saved_date_figure(row,path)
    assert 'No valid collection dates available' in path.read_text()


def test_redraw_cannot_escape_report_block(tmp_path):
    analysis={'cohorts':[{'cohort_id':'cohort_1','root_to_tip_figure':'../outside.svg'}],
              'root_to_tip':[diagnostic()]}
    with pytest.raises(ValueError,match='within its report block'):
        redraw_saved_date_figures(analysis,tmp_path)
    assert not (tmp_path.parent/'outside.svg').exists()


def test_normal_refresh_redraws_legacy_svg_without_changing_numeric_evidence(tmp_path,monkeypatch):
    import json
    from pathlib import Path
    from chronoclade.cgmlst import refresh
    from test_cgmlst_refresh import completed_run

    manifest=completed_run(tmp_path)
    main=json.loads(manifest.read_text())
    block=manifest.parent/main['blocks'][0]['block_id']
    analysis_path=block/'profile_analysis.json'
    analysis=json.loads(analysis_path.read_text())
    for row in analysis['root_to_tip']:
        row.pop('display_distance_units',None)
        row.pop('display_slope_units',None)
    analysis_path.write_text(json.dumps(analysis))
    checkpoint=block/'numerical_analysis.json'
    protected=[checkpoint,*[block/name for name in json.loads(checkpoint.read_text())['artifacts']]]
    saved={p:p.read_bytes() for p in protected}
    figure=block/analysis['cohorts'][0]['root_to_tip_figure']
    figure.write_text('legacy figure label')
    def no_fit(*a,**k):
        raise AssertionError('presentation refresh must not fit diagnostics')
    monkeypatch.setattr(refresh,'add_root_comparisons',no_fit)
    result=refresh.refresh_clock_reports(manifest)
    after=json.loads(analysis_path.read_text())
    assert after['root_to_tip']==analysis['root_to_tip']
    assert all(p.read_bytes()==value for p,value in saved.items())
    assert result['computed_with']==main['computed_with']
    assert result['clock_diagnostic_refresh']['recomputed'] is False
    audit=result['clock_diagnostic_refresh']['blocks'][0]['date_figure_presentation']
    assert audit[0]['status']=='redrawn_saved_evidence'
    assert 'Estimated allele differences' in figure.read_text()
    assert 'arbitrary tip root' in figure.read_text()
    assert Path(figure).is_file()
