import sys
from pathlib import Path

from chronoclade.native_pipeline import LineageFiles, _run_command, _run_dated_tree


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


def test_dated_tree_keeps_the_root_used_for_temporal_testing(tmp_path: Path, monkeypatch) -> None:
    files = LineageFiles.in_directory(tmp_path)
    observed: dict[str, object] = {}

    def capture(command, **kwargs):
        observed["command"] = command
        observed["inputs"] = kwargs["inputs"]

    monkeypatch.setattr("chronoclade.native_pipeline._run_command", capture)

    _run_dated_tree(files, sequence_length=1234, force=False)

    command = observed["command"]
    assert isinstance(command, list)
    assert str(files.clock_dir / "rerooted.newick") in command
    assert "--keep-root" in command
    assert "--reroot" not in command
    assert observed["inputs"] == (files.clock_dir / "rerooted.newick", files.metadata)


def test_unknown_locations_are_missing_traits_not_countries(tmp_path: Path) -> None:
    import csv
    from chronoclade.native_pipeline import _write_lineage_inputs
    from chronoclade.metadata import Sample

    samples = [
        Sample(
            str(i),
            tmp_path / f"{i}.fasta",
            "2020",
            location,
            "Klebsiella pneumoniae",
            "ST147",
            "local",
        )
        for i, location in enumerate(("Italy", "Unknown", "Public_context", "?", "N/A", "Missing"))
    ]
    _write_lineage_inputs(tmp_path, samples)
    states = list(csv.DictReader((tmp_path / "states.csv").open()))
    metadata = list(csv.DictReader((tmp_path / "metadata.csv").open()))
    assert [row["location"] for row in states] == ["Italy", "?", "?", "?", "?", "?"]
    assert [row["location"] for row in metadata] == [sample.location for sample in samples]


def test_lineage_metadata_csv_retains_optional_city_and_region(tmp_path: Path) -> None:
    import csv
    from chronoclade.native_pipeline import _write_lineage_inputs
    from chronoclade.metadata import Sample

    sample = Sample(
        "F1", tmp_path / "sample.fasta", "2020", "Greece", "Klebsiella pneumoniae",
        "ST147", "local", country="Greece", region="Attica", city="Athens",
    )
    _, metadata_path, _ = _write_lineage_inputs(tmp_path, [sample])
    row = next(csv.DictReader(metadata_path.open(newline="")))
    assert (row["country"], row["region"], row["city"]) == ("Greece", "Attica", "Athens")


def test_future_collection_dates_excluded_and_clock_uses_validated_dates(tmp_path, monkeypatch):
    from datetime import date
    from chronoclade.native_pipeline import _dated_members, _run_observed_clock
    from chronoclade.metadata import Sample

    valid = Sample("valid", tmp_path / "valid.fa", "2020", "UK", "Klebsiella pneumoniae", "ST147", "local")
    future = Sample("future", tmp_path / "future.fa", str(date.today().year + 10),
                    "UK", "Klebsiella pneumoniae", "ST147", "context")
    assert _dated_members([valid, future]) == [valid]
    assert future.collection_date == str(date.today().year + 10)
    files = LineageFiles.in_directory(tmp_path)
    files.metadata.write_text("raw original metadata")
    clock_dates = tmp_path / "clock_dates.csv"
    clock_dates.write_text("sample_id,collection_date\nvalid,2020\nfuture,\n")
    calls = []
    monkeypatch.setattr("chronoclade.native_pipeline._run_command",
                        lambda command, **kwargs: calls.append((command, kwargs["inputs"])))
    _run_observed_clock(files, tmp_path / "tree.nwk", 100, False)
    _run_dated_tree(files, 100, False)
    assert all(str(clock_dates) in command and clock_dates in inputs
               for command, inputs in calls)
    assert files.metadata.read_text() == "raw original metadata"
