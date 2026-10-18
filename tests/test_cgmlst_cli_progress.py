"""Check progress UX with a small real offline cgMLST CLI run."""

from __future__ import annotations

import json
from pathlib import Path
from unittest import mock

from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.cgmlst import workflow


def fixture_dataset(directory: Path) -> Path:
    from chronoclade.datasets import LocusCatalogue, from_profile_records, write_dataset

    rows = []
    for sample_id, origin, code, first_allele in (
        ("query", "local", "0,0,1,0,1,0,0,0,0,0", "1"),
        ("context-a", "context", "0,0,1,0,1,0,0,0,0,0", "1"),
        ("context-b", "context", "0,0,1,0,1,0,0,0,0,0", "2"),
    ):
        rows.append(
            {
                "sample_id": sample_id,
                "origin": origin,
                "species": "Klebsiella pneumoniae",
                "mlst_st": "39",
                "collection_date": "2020",
                "country": "Greece",
                "cgmlst_scheme": "offline-test",
                "cgmlst_scheme_version": "1",
                "cgmlst_profile": {"locus-a": first_allele, "locus-b": "1"},
                "cgmlst_status": "resolved",
                "cglin_scheme": "scgMLST629_S",
                "cglin_scheme_version": "1",
                "cglin_raw": code,
                "cglin_status": "complete",
            }
        )
    dataset = from_profile_records(
        rows, [LocusCatalogue("offline-test", "1", ("locus-a", "locus-b"), "offline fixture")]
    )
    return write_dataset(dataset, directory / "prepared")


def scientific_result(output: Path) -> tuple[list[list[str]], list[dict]]:
    manifest = json.loads((output / "cgmlst.json").read_text())
    block = manifest["blocks"][0]
    bundle = output / block["block_id"] / "selection_bundle"
    ensemble = json.loads((bundle / "ensemble.json").read_text())
    selected = []
    for run in ensemble["selections"]:
        selection = json.loads((bundle / run["path"]).read_text())
        selected.append(selection["selected_sample_ids"])
    distances = json.loads((bundle / "distance_evidence.json").read_text())["records"]
    return selected, distances


def test_normal_cli_progress_and_quiet_leave_science_identical(tmp_path):
    manifest = fixture_dataset(tmp_path)
    runner = CliRunner()

    def invoke(output: Path, quiet: bool):
        args = [
            "cgmlst",
            str(manifest),
            "--out",
            str(output),
            "--subsamples",
            "1",
            "--context-size",
            "1",
            "--nearest-per-query",
            "1",
            "--bootstrap-replicates",
            "0",
        ]
        if quiet:
            args.insert(0, "--quiet")
        emit_progress = workflow.progress

        def check_final_summary(stage, message):
            if message.startswith("Stage summary:"):
                assert (output / "cgmlst.json").is_file()
                assert (output / "index.html").is_file()
                assert (output / "README.md").is_file()
            emit_progress(stage, message)

        with mock.patch.object(workflow, "progress", side_effect=check_final_summary):
            result = runner.invoke(app, args)
        assert result.exit_code == 0, result.output
        return result

    normal = invoke(tmp_path / "normal", quiet=False)
    quiet = invoke(tmp_path / "quiet", quiet=True)
    for message in (
        "Block 1/1: 1 inputs, 2 contexts, 3/3 profiles available",
        "Starting block 1/1 numerical analysis",
        "Completed block 1/1 numerical analysis",
        "root/date diagnostics",
        "Starting block 1/1 selection bundle",
        "requested context budget 1; mandatory neighbours may expand it",
        "Completed block 1/1 selection bundle",
        "report complete",
        "Stage summary:",
    ):
        assert message in normal.output
    assert "[cgmlst]" not in quiet.output
    assert "cgMLST manifest:" in quiet.output
    assert "cgMLST report:" in quiet.output
    assert scientific_result(tmp_path / "normal") == scientific_result(tmp_path / "quiet")
