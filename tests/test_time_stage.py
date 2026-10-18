"""Dating reads saved evidence, preserves products and reports unavailable runs."""

import json

import pytest

from chronoclade.errors import WorkflowError
from chronoclade.time_stage.workflow import run_time
from chronoclade.time_stage.comparison import compare_time_runs
from chronoclade.tree_stage.workflow import run_tree
from test_tree_stage import selection_fixture, native_stubs as _tree_native_stubs


@pytest.fixture
def native_stubs(monkeypatch):
    return _tree_native_stubs.__wrapped__(monkeypatch)


def test_unsupported_time_needs_no_native_tools_or_tree_rebuild(
    tmp_path, native_stubs, monkeypatch
):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    tree = run_tree(selection, output=tmp_path / "tree")
    monkeypatch.setattr(
        "chronoclade.time_stage.workflow.require_tools",
        lambda names: (_ for _ in ()).throw(AssertionError("no tools for unsupported evidence")),
    )
    result = run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    assert result["dating_status"] == "unsupported" and result["dated_tree"] is None
    report = (tmp_path / "time" / result["report"]["path"]).read_text()
    assert "<h1>Dating assessment</h1>" in report
    assert "dated tree unavailable" in report and "· years" not in report
    assert report.index("Dating unsupported") < report.index("Selected samples")
    assert "Biological tree and location history" not in report
    from urllib.parse import unquote
    from re import search

    upstream = search(r'<a href="([^"]+)">SNP-tree report</a>', report)
    assert upstream
    report_path = tmp_path / "time" / result["report"]["path"]
    assert (report_path.parent / unquote(upstream[1])).resolve() == (
        tmp_path / "tree" / tree["report"]["path"]
    ).resolve()
    assert run_time(tmp_path / "tree/tree.json", output=tmp_path / "time") == result
    assert native_stubs == ["corrected"]
    assert json.loads((tmp_path / "tree/tree.json").read_text()) == tree
    with pytest.raises(WorkflowError, match="three saved eligible"):
        run_time(tmp_path / "tree/tree.json", output=tmp_path / "override", allow_unsupported=True)
    assert not (tmp_path / "override/time.json").exists()


def test_ensemble_comparison_keeps_unsupported_and_absent_runs(tmp_path, native_stubs):
    selection, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "tree")
    run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    comparison = compare_time_runs(
        ensemble, [tmp_path / "time/time.json"], output=tmp_path / "comparison"
    )
    assert comparison["requested_runs"] == 2 and comparison["dated_runs"] == 0
    assert [run["status"] for run in comparison["runs"]] == ["unsupported", "unavailable"]
    assert comparison["sensitivity"]["rate"]["assessable_runs"] == 0
    assert "not confidence intervals" in comparison["interpretation"]
    assert comparison["target_pair_agreement"] == []
    report = (tmp_path / "comparison/report.html").read_text()
    assert "shared anchor genomes" in report and "· years" not in report
    assert "Selected samples" not in report
    assert "<h2>Ensemble assessment</h2>" in report
    assert '<div class="verdict not_supported">' not in report


def test_time_rejects_fabricated_gate_in_rehashed_manifest(tmp_path, native_stubs):
    from chronoclade.stage_artifacts import load_result

    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "tree")
    run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    manifest = tmp_path / "time/time.json"
    edited = json.loads(manifest.read_text())
    edited["temporal_assessment"]["supported"] = True
    manifest.write_text(json.dumps(edited))
    with pytest.raises(WorkflowError, match="differs from the saved tree"):
        load_result(manifest, "time")


@pytest.mark.integration
def test_native_treetime_dating_saved_tree_without_phylogeny_rebuild(
    tmp_path, native_stubs, monkeypatch
):
    import shutil
    from chronoclade.stage_artifacts import load_result

    if not shutil.which("treetime"):
        pytest.skip("Native TreeTime unavailable")
    import chronoclade.native_pipeline as lineage

    original = lineage._run_core_phylogeny

    def synthetic_tree(files, members, *args, **kwargs):
        result = original(files, members, *args, **kwargs)
        newick = (
            "("
            + ",".join(
                f"{sample.sample_id}:{(index + 1) * 0.001}" for index, sample in enumerate(members)
            )
            + ");"
        )
        files.tree.write_text(newick)
        files.starting_tree.write_text(newick)
        return result

    monkeypatch.setattr(lineage, "_run_core_phylogeny", synthetic_tree)
    selection, _ = selection_fixture(tmp_path)
    run_tree(selection, output=tmp_path / "tree", assess_temporal=False)
    original_dating = lineage._run_dated_tree
    observed_dating = []
    def dating(*args, **kwargs):
        observed_dating.append(kwargs)
        return original_dating(*args, **kwargs)
    monkeypatch.setattr(lineage, "_run_dated_tree", dating)
    result = run_time(tmp_path / "tree/tree.json", output=tmp_path / "time", allow_unsupported=True)
    assert observed_dating == [{"seed": 42, "keep_polytomies": True}]
    from Bio import Phylo
    dated = Phylo.read(tmp_path / "time" / result["dated_tree"]["path"], "nexus")
    assert len(dated.root.clades) == 3
    assert result["topology_preservation"] == {
        "method": "rooted_nontrivial_internal_leafsets_root_excluded",
        "source_clade_count": 0, "dated_clade_count": 0, "preserved": True,
    }
    assert result["dating_status"] == "dated" and result["dating_override"]
    assert result["clock_confidence"]["rate"] > 0
    assert load_result(tmp_path / "time/time.json", "time") == result
    assert native_stubs == ["corrected"]
    report = tmp_path / "time" / result["report"]["path"]
    assert "data-cc-network-library" in report.read_text()
    assert "<h3>Dated-tree location history</h3>" in report.read_text()
    assert 'src="timetree/timetree.svg"' not in report.read_text()


def test_comparison_retains_failed_job_reason(tmp_path, native_stubs):
    selection, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "tree")
    with pytest.raises(WorkflowError):
        run_time(tmp_path / "tree/tree.json", output=tmp_path / "time", allow_unsupported=True)
    failed = next((tmp_path / "time/jobs").glob("*/job.json"))
    result = compare_time_runs(ensemble, [failed], output=tmp_path / "compare")
    assert result["runs"][0]["status"] == "failed"
    assert "three saved eligible" in result["runs"][0]["reason"]


def test_new_tree_invalidates_time_but_preserves_previous_time_evidence(tmp_path, native_stubs):
    from chronoclade.stage_artifacts import load_result

    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "tree")
    first = run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    # The source snapshot remains valid while the root tree.json pointer advances.
    run_tree(selection, output=tmp_path / "tree", force=True)
    assert load_result(tmp_path / "time/time.json", "time") == first
    second = run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    assert second["fingerprint"] != first["fingerprint"]
    assert second["source"]["tree"]["sha256"] != first["source"]["tree"]["sha256"]


def test_time_preserves_tree_coherence_review(tmp_path, native_stubs):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    tree = run_tree(selection, output=tmp_path / "tree", coherence_policy="diagnostic")
    tree["coherence_screening"]["warning"] = "Review retained <genome> in the broad cohort."
    source = tmp_path / "tree/tree.json"
    source.write_text(json.dumps(tree))
    result = run_time(source, output=tmp_path / "time")
    assert result["coherence_screening"] == tree["coherence_screening"]
    report = (tmp_path / "time" / result["report"]["path"]).read_text()
    assert "Coherence warning:" in report
    assert "Review retained &lt;genome&gt;" in report
    assert "Review retained <genome>" not in report
    assert result["dating_status"] == "unsupported"


@pytest.mark.parametrize("seed,keep_polytomies", [(None, False), (20261013, True)])
def test_dated_tree_seed_and_polytomy_flags(tmp_path, monkeypatch, seed, keep_polytomies):
    from chronoclade.native_pipeline import LineageFiles, _run_dated_tree
    commands = []
    monkeypatch.setattr("chronoclade.native_pipeline._run_command", lambda command, **kwargs: commands.append(command))
    files = LineageFiles.in_directory(tmp_path)
    _run_dated_tree(files, 2000, False, seed=seed, keep_polytomies=keep_polytomies)
    command = commands[0]
    assert ("--keep-polytomies" in command) is keep_polytomies
    assert "--keep-root" in command
    if seed is None:
        assert "--rng-seed" not in command
    else:
        assert command[command.index("--rng-seed") + 1] == str(seed)


@pytest.mark.parametrize("prior_policy", ["unseeded", "unchecked_topology"])
def test_native_time_policy_invalidates_prior_unseeded_fingerprint(tmp_path, native_stubs, prior_policy):
    from chronoclade.stage_artifacts import fingerprint, reference
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    source = tmp_path / "tree/tree.json"
    run_tree(selection, output=source.parent, seed=20261013)
    output = tmp_path / "time"
    first = run_time(source, output=output)
    assert first["parameters"] == {"allow_unsupported": False, "dating_seed": 20261013, "keep_polytomies": True, "dating_inference_policy": "seeded_saved_topology_keep_polytomies", "topology_validation_policy": "rooted_nontrivial_internal_leafsets_root_excluded"}
    prior = {**first, "parameters": {"allow_unsupported": False}}
    if prior_policy == "unseeded":
        prior["fingerprint"] = fingerprint({"tree": reference(output, source)["sha256"], "tools": first["tools"], "allow_unsupported": False})
    else:
        prior["parameters"] = {key: value for key, value in first["parameters"].items() if key != "topology_validation_policy"}
        prior["fingerprint"] = fingerprint({"tree": reference(output, source)["sha256"], "tools": first["tools"], "parameters": prior["parameters"]})
    assert first["fingerprint"] != prior["fingerprint"]
    (output / "time.json").write_text(json.dumps(prior))
    second = run_time(source, output=output)
    assert second["fingerprint"] == first["fingerprint"]
    assert second["parameters"] == first["parameters"]
    assert json.loads((output / "time.json").read_text())["parameters"] == first["parameters"]


def test_topology_validation_ignores_node_labels_and_lengths():
    from io import StringIO
    from Bio import Phylo
    from chronoclade.time_stage.workflow import validate_rooted_topology
    source = Phylo.read(StringIO("((a:1,b:2)original:3,c:4)root;"), "newick")
    dated = Phylo.read(StringIO("(c:8,(b:7,a:6)renamed:5)dated;"), "newick")
    result = validate_rooted_topology(source, dated)
    assert result["source_clade_count"] == result["dated_clade_count"] == 1
    assert result["preserved"]


def test_same_tips_changed_rooted_topology_fails_before_publish(tmp_path, native_stubs, monkeypatch):
    from io import StringIO
    from Bio import Phylo
    selection, _ = selection_fixture(tmp_path)
    run_tree(selection, output=tmp_path / "tree", assess_temporal=False)
    monkeypatch.setattr("chronoclade.time_stage.workflow.require_tools", lambda names: None)
    def altered_dating(files, *args, **kwargs):
        # Source is a trifurcation; retaining exact tips while adding a clade is forbidden.
        files.time_tree.parent.mkdir(parents=True, exist_ok=True)
        Phylo.write(Phylo.read(StringIO("((q:1,c:1):1,d:2);"), "newick"), files.time_tree, "nexus")
    monkeypatch.setattr("chronoclade.native_pipeline._run_dated_tree", altered_dating)
    with pytest.raises(WorkflowError, match="changed the saved rooted topology"):
        run_time(tmp_path / "tree/tree.json", output=tmp_path / "time", allow_unsupported=True)
    assert not (tmp_path / "time/time.json").exists()
    failed = next((tmp_path / "time/jobs").glob("*/job.json"))
    assert json.loads(failed.read_text())["status"] == "failed"


def test_topology_validation_rejects_equal_counts_with_different_clades():
    from io import StringIO
    from Bio import Phylo
    from chronoclade.time_stage.workflow import validate_rooted_topology
    source = Phylo.read(StringIO("((a,b),c);"), "newick")
    altered = Phylo.read(StringIO("((a,c),b);"), "newick")
    with pytest.raises(WorkflowError, match="1 saved internal clade.*1 new internal clade"):
        validate_rooted_topology(source, altered)
