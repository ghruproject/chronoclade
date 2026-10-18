"""Clock report refresh must preserve its frozen tree and selection evidence."""

import json
from pathlib import Path

import pytest

from chronoclade.cgmlst.refresh import refresh_clock_reports
from chronoclade.errors import WorkflowError
from chronoclade.datasets import LocusCatalogue, from_profile_records, write_dataset
from chronoclade.cgmlst.workflow import run_cgmlst


def completed_run(tmp_path, *, display_lin_level=7):
    rows = [dict(sample_id=name, origin='local' if i==0 else 'context',
        species='Klebsiella pneumoniae', mlst_st='39', collection_date=str(2010+i),
        country='Greece', cgmlst_scheme='test', cgmlst_scheme_version='1',
        cgmlst_profile={'a':str(i+1),'b':'1','c':str(1+i//2)}, cgmlst_status='resolved',
        cglin_scheme='scgMLST629_S', cglin_scheme_version='1',
        cglin_raw='0,0,1,0,1,0,0', cglin_status='complete')
        for i,name in enumerate(('q','c1','c2','c3'))]
    dataset = write_dataset(from_profile_records(rows,[LocusCatalogue('test','1',('a','b','c'),'test')]),tmp_path/'prepared')
    return run_cgmlst(dataset,tmp_path/'run',bootstrap_replicates=0,context_size=2,
                      display_lin_level=display_lin_level,partition_by='lineage').manifest_path


def test_refresh_preserves_numeric_tree_and_selection(tmp_path):
    path = completed_run(tmp_path, display_lin_level=8)
    manifest = json.loads(path.read_text())
    directory = path.parent / manifest['blocks'][0]['block_id']
    analysis = json.loads((directory/'profile_analysis.json').read_text())
    tree = directory / analysis['cohorts'][0]['tree_path']
    selection = directory/'selection_bundle/ensemble.json'
    coordinates = [directory / analysis['cohorts'][0][key]
                   for key in ('pcoa_csv', 'direct_pcoa_csv')]
    frozen = [tree, selection, *coordinates]
    snapshots = [file.read_bytes() for file in frozen]
    result = refresh_clock_reports(path)
    assert snapshots == [file.read_bytes() for file in frozen]
    assert result['clock_diagnostic_refresh']['blocks'][0]['sample_count']==4
    refreshed = json.loads((directory/'profile_analysis.json').read_text())
    assert refreshed['display_group_level'] == 8
    report = Path(path.parent/result['blocks'][0]['report']).read_text()
    assert 'How much does the root affect the slope?' not in report
    assert 'Rooting method and sensitivity checks' not in report
    assert 'Archivo' in report
    assert '1 input genomes · 3 public comparisons' in report
    assert 'Browse output files by type' in report
    assert (directory / 'README.md').is_file()
    assert (path.parent / 'files/index.html').is_file()


def test_refresh_rejects_changed_coordinates(tmp_path):
    path = completed_run(tmp_path)
    manifest = json.loads(path.read_text())
    directory = path.parent / manifest['blocks'][0]['block_id']
    analysis = json.loads((directory/'profile_analysis.json').read_text())
    coordinates = directory / analysis['cohorts'][0]['pcoa_csv']
    coordinates.write_text(coordinates.read_text().replace('0.', '1.'))
    with pytest.raises(WorkflowError, match='numerical checkpoint'):
        refresh_clock_reports(path)


def test_refresh_rejects_changed_source_tree(tmp_path):
    path = completed_run(tmp_path)
    manifest=json.loads(path.read_text())
    directory=path.parent/manifest['blocks'][0]['block_id']
    analysis=json.loads((directory/'profile_analysis.json').read_text())
    (directory/analysis['cohorts'][0]['tree_path']).write_text('(q:10,c1:10,c2:10,c3:10);')
    with pytest.raises(WorkflowError,match='numerical checkpoint'):
        refresh_clock_reports(path)


def test_presentation_refresh_preserves_computation_identity_and_skips_fits(monkeypatch, tmp_path):
    from chronoclade.cgmlst import refresh
    path = completed_run(tmp_path)
    before = json.loads(path.read_text())['computed_with']
    def forbidden(*args, **kwargs):
        raise AssertionError('presentation refresh must not recompute clock fits')
    monkeypatch.setattr(refresh, 'add_root_comparisons', forbidden)
    result = refresh.refresh_clock_reports(path)
    assert result['computed_with'] == before
    assert result['rendered_with']['component'] == 'cgmlst-report'
    assert result['clock_diagnostic_refresh']['recomputed'] is False
    assert result['clock_diagnostic_refresh']['computed_with'] is None


def test_explicit_diagnostic_refresh_records_separate_identity(monkeypatch, tmp_path):
    from chronoclade.cgmlst import refresh
    path = completed_run(tmp_path)
    before = json.loads(path.read_text())['computed_with']
    original = refresh.add_root_comparisons
    calls = []
    def observed(*args, **kwargs):
        calls.append(True)
        return original(*args, **kwargs)
    monkeypatch.setattr(refresh, 'add_root_comparisons', observed)
    result = refresh.refresh_clock_reports(path, recompute_diagnostics=True)
    assert calls == [True]
    assert result['computed_with'] == before
    assert result['clock_diagnostic_refresh']['recomputed'] is True
    assert result['clock_diagnostic_refresh']['computed_with']['component'] == 'cgmlst-clock-diagnostics'
