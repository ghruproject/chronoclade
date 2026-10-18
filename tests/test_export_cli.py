"""Supported viewer commands retain checked source identity and reject stale data."""

import csv
import json

from Bio import Phylo
import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.errors import WorkflowError
from chronoclade.export_cli import export_saved
from test_cgmlst_refresh import completed_run
from test_tree_stage import selection_fixture, native_stubs as _native_stubs


@pytest.fixture
def native_stubs(monkeypatch):
    return _native_stubs.__wrapped__(monkeypatch)


def test_corrected_tree_export_and_unsupported_time_handoff(tmp_path, native_stubs):
    from chronoclade.tree_stage.workflow import run_tree
    from chronoclade.time_stage.workflow import run_time

    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    result = run_tree(selection, output=tmp_path / 'tree')
    evidence = export_saved(tmp_path / 'tree/tree.json', tmp_path / 'tree-viewer', 'microreact')
    assert evidence['scope']['branch_units'] == 'substitutions/site'
    metadata = list(csv.DictReader((tmp_path / 'tree-viewer/metadata.csv').open()))
    assert {row['chronoclade_sample_id'] for row in metadata} == {'q', 'c', 'd'}
    guide = (tmp_path / 'tree/README.md').read_text()
    assert result['report']['path'] in guide
    assert str(tmp_path / 'tree/tree.json') in guide
    dated = run_time(tmp_path / 'tree/tree.json', output=tmp_path / 'time')
    assert dated['report']['path'] in (tmp_path / 'time/README.md').read_text()
    with pytest.raises(WorkflowError, match='export the saved tree-stage'):
        export_saved(tmp_path / 'time/time.json', tmp_path / 'time-viewer', 'microreact')
    assert not (tmp_path / 'time-viewer').exists()


def test_cli_exports_prepared_profiles_with_exact_scope_and_matching_metadata(tmp_path):
    run = completed_run(tmp_path)
    dataset = tmp_path / 'prepared/dataset.json'
    target = tmp_path / 'grapetree'
    result = CliRunner().invoke(app, ['export', str(dataset), '--viewer', 'grapetree',
                                     '--out', str(target)])
    assert result.exit_code == 0, result.output
    metadata = list(csv.DictReader((target / 'metadata.tsv').open(), delimiter='\t'))
    with (target / 'profiles.tsv').open() as handle:
        profiles = list(csv.reader(handle, delimiter='\t'))
    assert [row[0] for row in profiles[1:]] == [row['ID'] for row in metadata]
    assert {row['chronoclade_sample_id'] for row in metadata} == {'q', 'c1', 'c2', 'c3'}
    evidence = json.loads((target / 'export.json').read_text())
    assert evidence['sample_count'] == 4
    assert evidence['scope']['locus_count'] == 3
    assert 'cgMLST' in result.output or 'grapetree' in result.output
    assert run.is_file()


def test_export_cgmlst_uses_complete_saved_cohort_not_assembly_subsample(tmp_path):
    run = completed_run(tmp_path)
    target = tmp_path / 'microreact'
    evidence = export_saved(run, target, 'microreact')
    metadata = list(csv.DictReader((target / 'metadata.csv').open()))
    tree = Phylo.read(target / 'tree.newick', 'newick')
    assert len(tree.get_terminals()) == len(metadata) == evidence['sample_count'] == 4
    assert {tip.name for tip in tree.get_terminals()} == {row['ID'] for row in metadata}
    assert evidence['scope']['cohort_id'] == 'cohort_1'
    assert evidence['scope']['checkpoint_sha256']


def test_export_rejects_changed_source_tree_without_partial_output(tmp_path):
    run = completed_run(tmp_path)
    manifest = json.loads(run.read_text())
    directory = run.parent / manifest['blocks'][0]['block_id']
    saved = json.loads((directory / 'numerical_analysis.json').read_text())
    tree = directory / saved['analysis']['cohorts'][0]['tree_path']
    tree.write_text('(q:10,c1:10,c2:10,c3:10);')
    with pytest.raises(WorkflowError, match='sha256'):
        export_saved(run, tmp_path / 'bad', 'microreact')
    assert not (tmp_path / 'bad').exists()


def test_export_rejects_overwrite_and_does_not_invent_prepared_tree(tmp_path):
    completed_run(tmp_path)
    dataset = tmp_path / 'prepared/dataset.json'
    target = tmp_path / 'existing'
    target.mkdir()
    (target / 'keep.txt').write_text('preserved')
    with pytest.raises(WorkflowError, match='already exists'):
        export_saved(dataset, target, 'grapetree')
    assert (target / 'keep.txt').read_text() == 'preserved'
    result = CliRunner().invoke(app, ['export', str(dataset), '--viewer', 'microreact',
                                     '--out', str(tmp_path / 'no-tree')])
    assert result.exit_code == 2 and 'needs a tree' in result.output
    assert not (tmp_path / 'no-tree').exists()


def test_export_requires_explicit_cohort_choice_when_ambiguous(tmp_path):
    run = completed_run(tmp_path)
    manifest = json.loads(run.read_text())
    directory = run.parent / manifest['blocks'][0]['block_id']
    checkpoint = directory / 'numerical_analysis.json'
    saved = json.loads(checkpoint.read_text())
    saved['analysis']['cohorts'].append(dict(saved['analysis']['cohorts'][0], cohort_id='other'))
    checkpoint.write_text(json.dumps(saved))
    with pytest.raises(WorkflowError, match='Choose --cohort'):
        export_saved(run, tmp_path / 'ambiguous', 'microreact')
    assert not (tmp_path / 'ambiguous').exists()


def test_cli_rejects_malformed_block_entries_with_useful_error(tmp_path):
    source = tmp_path / 'cgmlst.json'
    source.write_text(json.dumps({'schema': 'chronoclade.cgmlst.analysis',
        'schema_version': 1, 'status': 'complete', 'blocks': [None]}))
    target = tmp_path / 'viewer'
    result = CliRunner().invoke(app, ['export', str(source), '--viewer', 'microreact',
                                     '--out', str(target)])
    assert result.exit_code == 2 and 'list of objects' in result.output
    assert not target.exists()
