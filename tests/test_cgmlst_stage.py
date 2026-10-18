"""Frozen independent profile stages retain partitions, full trees and exact selections."""

import json

import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.context_partitions import partition_records, profile_records
from chronoclade.cgmlst.workflow import run_cgmlst
from chronoclade.datasets import LocusCatalogue, from_profile_records, write_dataset
from chronoclade.selection_manifest import load_selection_manifest


def fixture_dataset(tmp_path):
    records = []
    for ident, role, code, allele in [
        ("q", "local", "0,0,1,0,1,0,0,0,0,0", "1"),
        ("c1", "context", "0,0,1,0,1,0,0,0,0,0", "2"),
        ("c2", "context", "0,0,1,0,1,1,0,0,0,0", "3"),
        ("outside", "context", "0,0,1,0,2,0,0,0,0,0", "4"),
        ("unresolved", "local", None, "5"),
    ]:
        records.append(
            {
                "sample_id": ident,
                "origin": role,
                "species": "Klebsiella pneumoniae",
                "mlst_st": "39",
                "collection_date": "2020",
                "country": "Greece",
                "cgmlst_scheme": "test",
                "cgmlst_scheme_version": "1",
                "cgmlst_profile": {"a": allele, "b": "1"},
                "cgmlst_status": "resolved",
                "cglin_scheme": "scgMLST629_S",
                "cglin_scheme_version": "1",
                "cglin_raw": code,
                "cglin_status": "complete" if code else "missing",
            }
        )
    dataset = from_profile_records(
        records, [LocusCatalogue("test", "1", ("a", "b"), "frozen fixture")]
    )
    return dataset, write_dataset(dataset, tmp_path / "prepared")


def test_partitions_explicit_depth_full_membership_and_unresolved_audit(tmp_path):
    dataset, _ = fixture_dataset(tmp_path)
    records = profile_records(dataset)
    blocks, audit = partition_records(records)
    assert [row["sample_id"] for row in blocks[0]["records"]] == ["c1", "c2", "q"]
    assert audit["unmatched_context_ids"] == ["outside"]
    assert audit["unresolved_inputs"][0]["sample_id"] == "unresolved"
    assert {row["level"] for row in audit["available_level_counts"]} == {5, 6, 7, 8, 9, 10}
    narrowed, _ = partition_records(records, lin_level=6)
    assert [row["sample_id"] for row in narrowed[0]["records"]] == ["c1", "q"]
    records[1]["cglin_database_sha256"] = "a" * 64
    scoped, _ = partition_records(records)
    assert "c1" not in {row["sample_id"] for row in scoped[0]["records"]}


def test_real_stage_has_full_tree_and_validated_selections_without_assembly_retrieval(
    tmp_path, monkeypatch
):
    _, manifest = fixture_dataset(tmp_path)

    def forbidden(*args, **kwargs):
        raise AssertionError("Independent frozen profile analysis must not acquire assemblies")

    monkeypatch.setattr("chronoclade.profile_inputs.materialise_assemblies", forbidden)
    result = run_cgmlst(
        manifest,
        tmp_path / "cgmlst",
        replicates=2,
        context_size=1,
        nearest_per_query=1,
        bootstrap_replicates=0,
        partition_by="lineage",
    )
    record = json.loads(result.manifest_path.read_text())
    assert record["status"] == "partial"
    block = record["blocks"][0]
    assert block["input_count"] == 1 and block["context_count"] == 2
    directory = result.manifest_path.parent
    ensemble = json.loads((directory / block["ensemble"]).read_text())
    # Inspect the emitted exact saved selection via the public loader.
    for path in (directory / block["ensemble"]).parent.glob("selections/selection-*.json"):
        selection = load_selection_manifest(path)
        assert "q" in selection["selected_sample_ids"]
    analysis = json.loads((directory / block["block_id"] / "profile_analysis.json").read_text())
    assert analysis["records_count"] == 3
    comparison = analysis["root_comparisons"][0]
    assert comparison["audit"]["status"] == "skipped"
    assert "date" in comparison["audit"]["reason"]
    assert (directory / block["block_id"] / comparison["audit_path"]).is_file()
    assert ensemble
    newick_files = list((directory / block["block_id"]).glob("*.nwk"))
    if not newick_files:
        newick_files = list((directory / block["block_id"]).glob("*.newick"))
    assert newick_files
    tree = newick_files[0].read_text()
    assert all(ident in tree for ident in ("q", "c1", "c2"))
    assert result.report_path.is_file()
    index_html = result.report_path.read_text()
    assert "Available LIN depths" in index_html and "ST39" in index_html
    block_html = (directory / block["report"]).read_text()
    assert "cgMLST profile report" in block_html
    with pytest.raises(ValueError, match="empty"):
        run_cgmlst(manifest, tmp_path / "cgmlst")


def test_cgmlst_cli_errors_are_clear_and_do_not_publish_completion(tmp_path):
    result = CliRunner().invoke(
        app, ["cgmlst", str(tmp_path / "missing.json"), "--out", str(tmp_path / "out")]
    )
    assert result.exit_code == 2
    assert not (tmp_path / "out/cgmlst.json").exists()


def test_lineage_evidence_cannot_change_validated_identity_or_paths(tmp_path):
    from dataclasses import replace

    dataset, _ = fixture_dataset(tmp_path)
    lineage = dict(dataset.lineages[0])
    lineage["evidence"] = dict(
        lineage["evidence"], sample_id="escape", role="context", cgmlst_profile={"a": "999"}
    )
    changed = replace(dataset, lineages=(lineage, *dataset.lineages[1:]))
    rows = profile_records(changed)
    assert rows[0]["sample_id"] == "q" and rows[0]["role"] == "input"
    assert rows[0]["cgmlst_profile"]["a"] == "1"
    rows[0]["mlst_st"] = "39/../../escape"
    blocks, _ = partition_records(rows)
    assert all("/" not in block["block_id"] and ".." not in block["block_id"] for block in blocks)
    duplicated = replace(dataset, lineages=(*dataset.lineages, dataset.lineages[0]))
    with pytest.raises(ValueError, match="Repeated lineage"):
        profile_records(duplicated)


def test_cli_stage_help_is_lazy():
    import subprocess
    import sys

    script = (
        "import sys; from chronoclade.cli import app; from typer.testing import CliRunner; "
        "assert CliRunner().invoke(app,['--help']).exit_code==0; "
        "assert 'chronoclade.native_tools' not in sys.modules; "
        "assert 'chronoclade.native_pipeline' not in sys.modules; "
        "assert 'torch' not in sys.modules; assert 'esm' not in sys.modules"
    )
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize(
    "update,message", [({"kind": None}, "kind"), ({"evidence": ["bad"]}, "mapping")]
)
def test_malformed_generic_lineage_evidence_has_clear_error(tmp_path, update, message):
    from dataclasses import replace

    dataset, _ = fixture_dataset(tmp_path)
    changed = dict(dataset.lineages[0], **update)
    with pytest.raises(ValueError, match=message):
        profile_records(replace(dataset, lineages=(changed, *dataset.lineages[1:])))


def test_unknown_lineage_namespace_is_not_a_shared_partition(tmp_path):
    dataset, _ = fixture_dataset(tmp_path)
    rows = profile_records(dataset)
    for row in rows:
        row["cglin_scheme_version"] = "unknown"
    blocks, audit = partition_records(rows)
    assert blocks == []
    assert len(audit["unresolved_inputs"]) == 2


def test_distinct_mlst_namespaces_do_not_share_a_st_block(tmp_path):
    dataset, _ = fixture_dataset(tmp_path)
    rows = profile_records(dataset)
    rows[0]["mlst_scheme"] = "primary-scheme"
    rows[1]["mlst_scheme"] = "another-scheme"
    blocks, audit = partition_records(rows)
    assert "c1" in audit["unmatched_context_ids"]
    assert all(row["sample_id"] != "c1" for row in blocks[0]["records"])


def test_missing_mlst_metadata_can_compare_with_one_explicit_namespace(tmp_path):
    dataset, _ = fixture_dataset(tmp_path)
    rows = profile_records(dataset)
    rows[1]["mlst_scheme"] = "primary-scheme"
    blocks, audit = partition_records(rows)
    assert {r["sample_id"] for r in blocks[0]["records"]} == {"q", "c1", "c2"}
    assert next(r for r in blocks[0]["records"] if r["sample_id"] == "q")["mlst_scheme"] is None
    assert audit["ambiguous_mlst_namespace_ids"] == []


@pytest.mark.parametrize("packager", ["write_selection_ensemble", "write_profile_report"])
@pytest.mark.parametrize("backend", ["dense", "scale"])
def test_native_resume_reuses_completed_numerical_analysis(tmp_path, monkeypatch, packager, backend):
    import chronoclade.cgmlst.workflow as workflow

    if backend == "scale":
        monkeypatch.setattr("chronoclade.profile_analysis.MAX_RECORDS", 2)
    _, source = fixture_dataset(tmp_path)
    output = tmp_path / "interrupted"
    arguments = dict(bootstrap_replicates=0, context_size=1, nearest_per_query=1, partition_by="lineage")
    original = getattr(workflow, packager)

    def interrupted(*args, **kwargs):
        raise OSError("fixture packaging interruption")

    monkeypatch.setattr(workflow, packager, interrupted)
    with pytest.raises(OSError, match="packaging interruption"):
        workflow.run_cgmlst(source, output, **arguments)
    assert list(output.glob("*/numerical_analysis.json"))
    assert not (output / "cgmlst.json").exists()

    def forbidden(*args, **kwargs):
        raise AssertionError("Resume must not recalculate numerical analysis")

    monkeypatch.setattr(workflow, "analyse_profiles", forbidden)
    monkeypatch.setattr("chronoclade.cgmlst.scale.analyse_profiles_scale", forbidden)
    monkeypatch.setattr(workflow, packager, original)
    if packager == "write_profile_report":
        monkeypatch.setattr(workflow, "write_selection_ensemble", forbidden)
    with pytest.raises(ValueError, match="source or parameters differ"):
        workflow.run_cgmlst(source, output, resume=True, **{**arguments, "seed": 99})
    original_source = source.read_text()
    source.write_text(original_source + "\n")
    with pytest.raises(ValueError, match="source or parameters differ"):
        workflow.run_cgmlst(source, output, resume=True, **arguments)
    source.write_text(original_source)
    result = workflow.run_cgmlst(source, output, resume=True, **arguments)
    manifest = json.loads(result.manifest_path.read_text())
    assert result.report_path.is_file()
    assert manifest["status"] == "partial"
    block = manifest["blocks"][0]
    assert (output / block["ensemble"]).is_file()
    for path in (output / block["ensemble"]).parent.glob("selections/selection-*.json"):
        assert "q" in load_selection_manifest(path)["selected_sample_ids"]


def test_resume_rejects_changed_numerical_evidence(tmp_path, monkeypatch):
    import chronoclade.cgmlst.workflow as workflow

    _, source = fixture_dataset(tmp_path)
    output = tmp_path / "interrupted"
    arguments = dict(bootstrap_replicates=0)

    def interrupted(*args, **kwargs):
        raise OSError("fixture packaging interruption")

    monkeypatch.setattr(workflow, "write_selection_ensemble", interrupted)
    with pytest.raises(OSError):
        workflow.run_cgmlst(source, output, **arguments)
    checkpoint = next(output.glob("*/numerical_analysis.json"))
    saved = json.loads(checkpoint.read_text())
    evidence = checkpoint.parent / next(iter(saved["artifacts"]))
    evidence.write_text("changed")
    with pytest.raises(ValueError, match="numerical artifact changed"):
        workflow.run_cgmlst(source, output, resume=True, **arguments)
    assert not (output / "cgmlst.json").exists()


def test_cli_resume_repackages_saved_analysis(tmp_path, monkeypatch):
    import chronoclade.cgmlst.workflow as workflow

    _, source = fixture_dataset(tmp_path)
    output = tmp_path / "interrupted"
    original = workflow.write_profile_report

    def interrupted(*args, **kwargs):
        raise OSError("fixture packaging interruption")

    monkeypatch.setattr(workflow, "write_profile_report", interrupted)
    with pytest.raises(OSError):
        workflow.run_cgmlst(source, output, bootstrap_replicates=0, representative_reserve=1)
    monkeypatch.setattr(workflow, "write_profile_report", original)
    result = CliRunner().invoke(app, [
        "cgmlst", str(source), "--out", str(output), "--bootstrap-replicates", "0",
        "--representative-reserve", "1", "--resume",
    ])
    assert result.exit_code == 0, result.output
    assert json.loads((output / "cgmlst.json").read_text())["parameters"]["representative_reserve"] == 1


def test_resume_can_start_analysis_after_retrieval_interruption(tmp_path):
    _, source = fixture_dataset(tmp_path)
    result = run_cgmlst(source, tmp_path / "not_started", resume=True, bootstrap_replicates=0)
    assert result.manifest_path.is_file()
    unrecognised = tmp_path / "unrecognised"
    unrecognised.mkdir()
    (unrecognised / "foreign.txt").write_text("unrecognised output")
    with pytest.raises(ValueError, match="No native cgmlst resume checkpoint"):
        run_cgmlst(source, unrecognised, resume=True, bootstrap_replicates=0)


def test_resume_rejects_injected_valid_selection_bundle_from_other_parameters(tmp_path, monkeypatch):
    import shutil
    import chronoclade.cgmlst.workflow as workflow
    from chronoclade.selection_manifest import load_selection_ensemble, write_selection_ensemble

    dataset, source = fixture_dataset(tmp_path)
    output = tmp_path / "interrupted"

    def interrupted(*args, **kwargs):
        raise OSError("fixture packaging interruption")

    monkeypatch.setattr(workflow, "write_profile_report", interrupted)
    with pytest.raises(OSError):
        workflow.run_cgmlst(source, output, bootstrap_replicates=0, context_size=1, partition_by="lineage")
    checkpoint = next(output.glob("*/numerical_analysis.json"))
    directory = checkpoint.parent
    analysis = workflow._restore_artifact_paths(json.loads(checkpoint.read_text())["analysis"], directory)
    block = partition_records(profile_records(dataset))[0][0]
    analysis["partition"] = {key: value for key, value in block.items() if key != "records"}
    replacement = write_selection_ensemble(
        [row for row in block["records"] if row["role"] == "input"],
        [row for row in block["records"] if row["role"] == "context"],
        analysis, dataset_manifest=source, output=directory / "replacement_bundle",
        size=1, seed=99, partition_manifest=output / "partitions.json",
    )
    shutil.rmtree(directory / "selection_bundle")
    replacement.parent.rename(directory / "selection_bundle")
    ensemble = directory / "selection_bundle" / "ensemble.json"
    assert load_selection_ensemble(ensemble)["parameters"]["seed"] == 99
    with pytest.raises(ValueError, match="selection bundle differs"):
        workflow.run_cgmlst(source, output, bootstrap_replicates=0, context_size=1, partition_by="lineage", resume=True)
    assert not (output / "cgmlst.json").exists()


@pytest.mark.parametrize("level", [5, 6, 7, 8, 9, 10])
def test_partition_supports_finer_reported_prefixes_without_distance_inference(tmp_path, level):
    dataset, _ = fixture_dataset(tmp_path)
    records = profile_records(dataset)
    records = [r for r in records if r["sample_id"] != "unresolved"]
    query = next(r for r in records if r["sample_id"] == "q")
    prefix = query["cglin_raw"].split(",")
    same = dict(query, sample_id="same", role="context", origin="context")
    different = dict(same, sample_id="different")
    differing_code = list(prefix)
    differing_code[level - 1] = "99"
    different["cglin_raw"] = ",".join(differing_code)
    short = dict(same, sample_id="short", cglin_raw=",".join(prefix[:level - 1]))
    blocks, audit = partition_records([query, same, different, short], lin_level=level)
    assert len(blocks) == 1
    assert blocks[0]["level"] == level
    assert {r["sample_id"] for r in blocks[0]["records"]} == {"q", "same"}
    assert audit["unmatched_context_ids"] == ["different", "short"]


@pytest.mark.parametrize("level", [True, 4, 11, 8.0, "8"])
def test_partition_rejects_invalid_depths(tmp_path, level):
    dataset, _ = fixture_dataset(tmp_path)
    with pytest.raises(ValueError, match="integer between 5 and 10"):
        partition_records(profile_records(dataset), lin_level=level)


def test_cli_finer_scope_and_display_are_explicit_and_selection_independent(tmp_path, monkeypatch):
    import chronoclade.cgmlst.workflow as workflow

    original_selection = workflow.write_selection_ensemble
    selection_groupings = []

    def record_selection(queries, contexts, analysis, **kwargs):
        selection_groupings.append(analysis["genetic_groups"])
        return original_selection(queries, contexts, analysis, **kwargs)

    monkeypatch.setattr(workflow, "write_selection_ensemble", record_selection)
    _, source = fixture_dataset(tmp_path)
    runner = CliRunner()
    selections = []
    for display in (7, 10):
        output = tmp_path / f"display-{display}"
        result = runner.invoke(app, ["cgmlst", str(source), "--out", str(output),
            "--partition-by", "lineage", "--lin-level", "8", "--display-lin-level", str(display),
            "--bootstrap-replicates", "0", "--nearest-per-query", "1"])
        assert result.exit_code == 2, result.output  # Unresolved input remains audited.
        manifest = json.loads((output / "cgmlst.json").read_text())
        assert manifest["parameters"]["lin_level"] == 8
        assert manifest["parameters"]["display_lin_level"] == display
        assert manifest["parameters"]["selection_group_lin_level"] == 7
        assert manifest["blocks"][0]["context_count"] == 1
        block = output / manifest["blocks"][0]["block_id"]
        analysis = json.loads((block / "profile_analysis.json").read_text())
        assert analysis["display_group_level"] == display
        assert analysis["partition"]["level"] == 8
        selection_path = next((block / "selection_bundle/selections").glob("selection-*.json"))
        selections.append(load_selection_manifest(selection_path)["selected_sample_ids"])
    assert selections[0] == selections[1]
    assert selection_groupings[0] == selection_groupings[1]
    assert all(group["level"] == 7 for group in selection_groupings[0])


def test_context_prefix_union_preserves_cross_prefix_neighbours_in_wider_analysis_block(tmp_path):
    dataset, _ = fixture_dataset(tmp_path)
    rows = [r for r in profile_records(dataset) if r["sample_id"] in {"q", "c1"}]
    q = rows[0]
    q.update(cglin_raw="0,0,1,0,1,0,0,0,0,0", cgmlst_profile={"a": "1", "b": "1"})
    c = rows[1]
    c.update(cglin_raw="0,0,1,0,1,0,0,1,0,0", cgmlst_profile={"a": "1", "b": "1"})
    rows.append(dict(q, sample_id="q2", cglin_raw=c["cglin_raw"], role="input",
                     origin="local", cgmlst_profile={"a": "3", "b": "1"}))
    rows.append(dict(c, sample_id="outside", cglin_raw="0,0,1,0,1,0,0,2,0,0"))
    source = write_dataset(from_profile_records(rows,
        [LocusCatalogue("test", "1", ("a", "b"), "fixture")]), tmp_path / "union-prepared")
    output = tmp_path / "union"
    result = CliRunner().invoke(app, ["cgmlst", str(source), "--out", str(output),
        "--partition-by", "lineage", "--lin-level", "5", "--context-lin-level", "8", "--display-lin-level", "8",
        "--bootstrap-replicates", "0", "--nearest-per-query", "1"])
    assert result.exit_code == 0, result.output
    manifest = json.loads((output / "cgmlst.json").read_text())
    assert len(manifest["blocks"]) == 1
    block = manifest["blocks"][0]
    assert block["input_count"] == 2 and block["context_count"] == 1
    analysis = json.loads((output / block["block_id"] / "profile_analysis.json").read_text())
    matches = {(r["query_id"], r["context_id"]) for r in analysis["nearest_neighbours"]
               if r["status"] == "matched"}
    assert (q["sample_id"], c["sample_id"]) in matches  # Different own LIN8 prefixes.
    audit = json.loads((output / "partitions.json").read_text())
    assert audit["outside_context_scope_ids"] == ["outside"]
    assert audit["context_scope_counts"] == {
        "prepared_context": 2, "scope_context": 1, "outside_scope_context": 1}


@pytest.mark.parametrize("parameters, message", [
    ({"backend": "invalid"}, "backend must be"),
    ({"backend": "scale", "distance_method": "direct-dna"}, "Scale backend supports allele"),
])
def test_backend_rejects_invalid_combinations_before_reading_input(tmp_path, parameters, message):
    with pytest.raises(ValueError, match=message):
        run_cgmlst(tmp_path / "missing.json", tmp_path / "out", **parameters)
    assert not (tmp_path / "out").exists()


def test_missing_finer_context_assignment_is_reported_without_dropping_input(tmp_path):
    dataset, _ = fixture_dataset(tmp_path)
    q = next(r for r in profile_records(dataset) if r["sample_id"] == "q")
    q["cglin_raw"] = "0,0,1,0,1"
    blocks, audit = partition_records([q], lin_level=5, context_lin_level=8)
    assert blocks[0]["records"][0]["sample_id"] == "q"
    assert audit["unresolved_context_scope_inputs"] == [
        {"sample_id": "q", "reason": "missing_or_unresolved_context_scope_typing"}]
