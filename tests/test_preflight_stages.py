"""Lean environments must not fail checks for unrelated native stages."""

from typer.testing import CliRunner

from chronoclade.cli import app


def test_lean_prepare_preflight_and_missing_tree_tools(monkeypatch):
    monkeypatch.setattr('chronoclade.native_tools.tool_status', lambda: {'ska': None, 'iqtree': None, 'ClonalFrameML': None, 'treetime': None})
    monkeypatch.setattr('chronoclade.native_tools.native_platform_supported', lambda: False)
    runner = CliRunner()
    lean = runner.invoke(app, ['preflight', '--stage', 'prepare'])
    assert lean.exit_code == 0
    assert 'Software check only' in lean.output
    assert 'ska' not in lean.output
    native = runner.invoke(app, ['preflight', '--stage', 'tree'])
    assert native.exit_code == 2
    assert 'missing' in native.output


def test_time_preflight_only_checks_time_tool(monkeypatch):
    monkeypatch.setattr('chronoclade.native_tools.tool_status', lambda: {'ska': None, 'iqtree': None, 'treetime': '/bin/treetime'})
    monkeypatch.setattr('chronoclade.native_tools.native_platform_supported', lambda: True)
    result = CliRunner().invoke(app, ['preflight', '--stage', 'time'])
    assert result.exit_code == 0
    assert '/bin/treetime' in result.output
    assert 'iqtree' not in result.output


def test_invalid_stage_fails_with_next_choices():
    result = CliRunner().invoke(app, ['preflight', '--stage', 'invalid'])
    assert result.exit_code == 2
    assert 'prepare, cgmlst, tree or time' in result.output
