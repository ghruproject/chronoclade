"""Keep one public analysis workflow, without aliases for retired entry points."""

import pytest
from typer.main import get_command
from typer.testing import CliRunner

from chronoclade.cli import app


def test_four_analysis_commands_are_registered_without_retired_interfaces():
    commands = get_command(app).commands
    assert {"prepare", "cgmlst", "tree", "time"} <= commands.keys()
    assert not {"run", "prepare-context", "validate"} & commands.keys()


@pytest.mark.parametrize("command", ["prepare", "cgmlst", "tree", "time"])
def test_analysis_help_is_available(command):
    result = CliRunner().invoke(app, [command, "--help"])
    assert result.exit_code == 0, result.output


@pytest.mark.parametrize("command", ["run", "prepare-context", "validate"])
def test_retired_commands_are_rejected_without_creating_outputs(tmp_path, monkeypatch, command):
    monkeypatch.chdir(tmp_path)
    result = CliRunner().invoke(app, [command, "--help"])
    assert result.exit_code == 2
    assert "No such command" in result.output
    assert not list(tmp_path.iterdir())
