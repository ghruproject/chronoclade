"""Computational failures must reach shell callers; unsupported dating is success."""
import json

import pytest
from typer.testing import CliRunner

from chronoclade.cli import app


@pytest.mark.parametrize('stage,schema,runner', [
    ('tree', 'chronoclade.selection-ensemble', 'run_tree_ensemble'),
    ('time', 'chronoclade.tree-ensemble', 'run_time_ensemble'),
])
@pytest.mark.parametrize('failed,expected', [(0, 0), (1, 2)])
def test_failed_ensemble_returns_nonzero(tmp_path, monkeypatch, stage, schema, runner, failed, expected):
    source = tmp_path/'source.json'
    source.write_text(json.dumps({'schema': schema}))
    result = {'schema': f'chronoclade.{stage}-ensemble', 'completed_runs': 2-failed,
              'failed_runs': failed, 'report': {'path': 'report.html'}}
    monkeypatch.setattr(f'chronoclade.execution.ensemble.{runner}', lambda *a, **kw: result)
    monkeypatch.setattr('chronoclade.tree_stage.cli._api_key', lambda environment: '')
    actual = CliRunner().invoke(app, [stage, str(source), '--out', str(tmp_path/'out')])
    assert actual.exit_code == expected
    assert f'failed: {failed}' in actual.output


def test_legitimate_unsupported_time_is_success(tmp_path, monkeypatch):
    source = tmp_path/'tree.json'
    source.write_text(json.dumps({'schema': 'chronoclade.tree'}))
    monkeypatch.setattr('chronoclade.time_stage.workflow.run_time', lambda *a, **kw: {'schema': 'chronoclade.time', 'dating_status': 'unsupported'})
    result = CliRunner().invoke(app, ['time', str(source), '--out', str(tmp_path/'out')])
    assert result.exit_code == 0
    assert 'unsupported' in result.output


def test_failed_tree_prints_actionable_reason_and_audit_path(tmp_path, monkeypatch):
    source = tmp_path / 'ensemble.json'
    source.write_text(json.dumps({'schema': 'chronoclade.selection-ensemble'}))
    reason = 'Selected assembly acquisition failed. See /analysis/[study]/assembly_audit.json'
    monkeypatch.setattr('chronoclade.execution.ensemble.run_tree_ensemble', lambda *a, **kw: {
        'schema': 'chronoclade.tree-ensemble', 'completed_runs': 0, 'failed_runs': 1,
        'runs': [{'selection_id': 'selection-001', 'status': 'failed', 'reason': reason}],
    })
    monkeypatch.setattr('chronoclade.tree_stage.cli._api_key', lambda environment: '')
    result = CliRunner().invoke(app, ['tree', str(source), '--out', str(tmp_path / 'out')])
    assert result.exit_code == 2
    # Rich may wrap long paths; markup-like filenames must still be printed literally.
    assert 'Failed selection-001:' in result.output
    assert '/analysis/[study]/assembly_audit.json' in result.output.replace('\n', '')
