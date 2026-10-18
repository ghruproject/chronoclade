"""Current numerical cache identity is separate from reading saved history."""
import json

import pytest

from chronoclade.stage_artifacts import load_tree_result
from chronoclade.tree_stage import workflow as tree
from chronoclade.time_stage import workflow as time
from test_tree_stage import selection_fixture, native_stubs as tree_stubs


@pytest.fixture
def native_stubs(monkeypatch):
    return tree_stubs.__wrapped__(monkeypatch)


def test_tree_same_version_method_change_recomputes_and_history_loads(monkeypatch, tmp_path, native_stubs):
    selection, _ = selection_fixture(tmp_path, dates=('', '', ''))
    first = tree.run_tree(selection, output=tmp_path/'tree')
    historic = tmp_path/'tree/historic-tree.json'
    historic.write_text(json.dumps(first))
    original = tree.execution_provenance
    monkeypatch.setattr(tree, 'execution_provenance',
                        lambda component, revision, **kw: original(component, revision+'-changed', **kw))
    second = tree.run_tree(selection, output=tmp_path/'tree')
    assert native_stubs == ['corrected', 'corrected']
    assert first['fingerprint'] != second['fingerprint']
    assert first['software'] == second['software']
    assert load_tree_result(historic)['computed_with'] == first['computed_with']


def test_time_version_change_invalidates_completed_reuse(monkeypatch, tmp_path, native_stubs):
    import chronoclade.stage_artifacts as contracts
    selection, _ = selection_fixture(tmp_path, dates=('', '', ''))
    tree.run_tree(selection, output=tmp_path/'tree')
    first = time.run_time(tmp_path/'tree/tree.json', output=tmp_path/'time')
    assert time.run_time(tmp_path/'tree/tree.json', output=tmp_path/'time') == first
    monkeypatch.setattr(contracts, '__version__', 'test.changed')
    second = time.run_time(tmp_path/'tree/tree.json', output=tmp_path/'time')
    assert first['fingerprint'] != second['fingerprint']
    assert second['computed_with']['software']['version'] == 'test.changed'
    assert native_stubs == ['corrected']
