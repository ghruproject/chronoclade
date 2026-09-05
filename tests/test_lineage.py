import sys
from pathlib import Path

from chronoclade.lineage import LineageFiles, _run_command, _run_dated_tree


COPY_INPUT = (
    "from pathlib import Path; import sys; "
    "Path(sys.argv[2]).write_text(Path(sys.argv[1]).read_text())"
)


def test_stage_reruns_when_direct_input_changes(tmp_path: Path) -> None:
    source = tmp_path / "source.txt"
    output = tmp_path / "output.txt"
    log = tmp_path / "logs" / "copy.log"
    source.write_text("first", encoding="utf-8")
    command = [sys.executable, "-c", COPY_INPUT, str(source), str(output)]

    _run_command(command, log=log, expected=output, force=False, inputs=(source,))
    assert output.read_text(encoding="utf-8") == "first"

    source.write_text("second", encoding="utf-8")
    _run_command(command, log=log, expected=output, force=False, inputs=(source,))

    assert output.read_text(encoding="utf-8") == "second"
    assert (tmp_path / "logs" / "copy.fingerprint.json").is_file()


def test_dated_tree_keeps_the_root_used_for_temporal_testing(
    tmp_path: Path, monkeypatch
) -> None:
    files = LineageFiles.in_directory(tmp_path)
    observed: dict[str, object] = {}

    def capture(command, **kwargs):
        observed["command"] = command
        observed["inputs"] = kwargs["inputs"]

    monkeypatch.setattr("chronoclade.lineage._run_command", capture)

    _run_dated_tree(files, sequence_length=1234, force=False)

    command = observed["command"]
    assert isinstance(command, list)
    assert str(files.clock_dir / "rerooted.newick") in command
    assert "--keep-root" in command
    assert "--reroot" not in command
    assert observed["inputs"] == (files.clock_dir / "rerooted.newick", files.metadata)
