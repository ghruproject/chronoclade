"""Exact selection, immutable source evidence and independent tree jobs."""

import json
from io import StringIO
from pathlib import Path

import pytest
from Bio import Phylo

from chronoclade.datasets import PreparedDataset, samples_from_records, write_dataset
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import write_selection_ensemble
from chronoclade.stage_artifacts import load_tree_result
from chronoclade.tree_stage.workflow import run_tree


def selection_fixture(
    tmp_path, *, dates=("2018", "2019", "2020"), hosts=None,
    countries=None, regions=None, cities=None,
):
    rows = []
    for index, ident in enumerate(("q", "c", "d")):
        assembly = tmp_path / (ident + ".fasta")
        assembly.write_text(">contig\n" + "A" * 2000 + "\n")
        rows.append(
            {
                "sample_id": ident,
                "origin": "local" if ident == "q" else "context",
                "species": "Klebsiella pneumoniae",
                "collection_date": dates[index],
                "country": countries[index] if countries is not None else (
                    "Greece" if ident != "d" else "UK"
                ),
                "assembly_reference": str(assembly),
                **({"host": hosts[index]} if hosts is not None else {}),
                **({"region": regions[index]} if regions is not None else {}),
                **({"city": cities[index]} if cities is not None else {}),
            }
        )
    dataset = PreparedDataset(samples=samples_from_records(rows))
    manifest = write_dataset(dataset, tmp_path / "prepared")
    distance = tmp_path / "distances.csv"
    distance.write_text(
        "sample_id_1,sample_id_2,distance,allele_differences,shared_called_loci\nq,c,0.1,1,10\nq,d,0.2,2,10\nc,d,0.1,1,10\n"
    )
    ensemble = write_selection_ensemble(
        rows[:1],
        rows[1:],
        {"paths": {"pairwise_distances": str(distance)}},
        dataset_manifest=manifest,
        output=tmp_path / "selections",
        size=2,
        replicates=2,
    )
    return ensemble.parent / "selections/selection-001.json", ensemble


@pytest.fixture
def native_stubs(monkeypatch):
    import chronoclade.lineage as lineage

    monkeypatch.setattr("chronoclade.tree_stage.workflow.require_tools", lambda names: None)
    monkeypatch.setattr(
        "chronoclade.tree_stage.workflow.tools",
        lambda names: {name: "frozen-test-tool" for name in names},
    )
    calls = []

    def core(files, members, reference, threads, force, mode, coherence_policy="strict", seed=None):
        calls.append(mode)
        ids = [s.sample_id for s in members]
        newick = "(" + ",".join(f"{ident}:0.001" for ident in ids) + ");\n"
        files.starting_tree.write_text(newick)
        files.tree.write_text(newick)
        sequences = "".join(f">{ident}\n" + "A" * 2000 + "\n" for ident in ids)
        files.alignment.write_text(sequences)
        files.filtered_alignment.write_text(sequences)
        files.importations.write_text("Node\tBeg\tEnd\n")
        from chronoclade.coherence import screen_alignment
        screen_alignment(files.alignment, members, output=files.coherence_screen)
        screen_alignment(files.filtered_alignment, members, output=files.clonal_coherence_screen)
        return files.tree, files.filtered_alignment, None

    monkeypatch.setattr(lineage, "_run_core_phylogeny", core)
    return calls


def test_tree_exact_ids_hashes_and_independent_resume(tmp_path, native_stubs):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    result = run_tree(selection, output=tmp_path / "tree")
    assert result["selected_sample_ids"] == ["q", "c", "d"]
    assert result["undated_sample_ids"] == ["q", "c", "d"]
    assert not result["temporal_assessment"]["supported"]
    assert load_tree_result(tmp_path / "tree/tree.json") == result
    assert run_tree(selection, output=tmp_path / "tree") == result
    assert native_stubs == ["corrected"]
    (tmp_path / "q.fasta").write_text(">contig\n" + "C" * 2000 + "\n")
    changed = run_tree(selection, output=tmp_path / "tree")
    assert changed["fingerprint"] != result["fingerprint"]
    assert len(native_stubs) == 2
    original_job = tmp_path / "tree" / result["corrected_tree"]["path"]
    assert original_job.is_file()
    report = tmp_path / "tree" / changed["report"]["path"]
    assert "data-cc-network-library" in report.read_text()
    assert "Date randomisations" not in report.read_text()


def test_failed_replacement_preserves_last_manifest(tmp_path, native_stubs, monkeypatch):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    saved = run_tree(selection, output=tmp_path / "tree")

    def broken(*args, **kwargs):
        raise WorkflowError("native failure")

    monkeypatch.setattr("chronoclade.lineage._run_core_phylogeny", broken)
    with pytest.raises(WorkflowError, match="native failure"):
        run_tree(selection, output=tmp_path / "tree", force=True)
    assert load_tree_result(tmp_path / "tree/tree.json") == saved
    assert any(
        json.loads(path.read_text())["status"] == "failed"
        for path in (tmp_path / "tree/jobs").glob("*/job.json")
    )


def test_acquisition_failure_has_audit_without_replacement(tmp_path):
    from chronoclade.tree_stage.workflow import acquire_assemblies

    row = {"sample_id": "missing", "source_genome_id": "exact-id"}
    with pytest.raises(WorkflowError, match="no replacements"):
        acquire_assemblies(
            [row],
            tmp_path / "dataset.json",
            tmp_path,
            api_key="",
            base_url="https://pathogen.watch",
            threads=1,
        )
    assert json.loads((tmp_path / "assembly_audit.json").read_text())[0]["sample_id"] == "missing"


def test_tree_rejects_changed_artifact_and_source(tmp_path, native_stubs):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    result = run_tree(selection, output=tmp_path / "tree")
    artifact = tmp_path / "tree" / result["corrected_tree"]["path"]
    artifact.write_text("(q:1,c:1,outsider:1);")
    with pytest.raises(WorkflowError, match="sha256"):
        load_tree_result(tmp_path / "tree/tree.json")


def test_clustered_dates_preserve_exact_partial_intervals():
    import random
    from chronoclade.metadata import Sample
    from chronoclade.tree_stage.temporal import date_clusters, permute_cluster_dates

    tree = Phylo.read(StringIO("((a:1,b:1):1,c:1,d:1);"), "newick")
    samples = [
        Sample(i, Path("."), date, "", "sp", "lin", "local")
        for i, date in zip("abcd", ("2018", "2018", "2019-02", "2020-03-04"))
    ]
    groups = date_clusters(tree, samples)
    assert [group["sample_ids"] for group in groups] == [["a", "b"], ["c"], ["d"]]
    shuffled = permute_cluster_dates(groups, random.Random(3))
    assert shuffled["a"] == shuffled["b"]
    assert set(shuffled.values()) == {"2018", "2019-02", "2020-03-04"}


def test_clustered_assessment_requires_complete_null_refits(tmp_path, monkeypatch):
    from chronoclade.metadata import Sample
    from chronoclade.tree_stage.temporal import run_clustered_assessment

    tree = tmp_path / "tree.nwk"
    tree.write_text("(a:0.001,b:0.002,c:0.003);")
    clock = tmp_path / "clock.txt"
    clock.write_text("--rate: 0.001\n--r^2: 0.95\n")
    samples = [
        Sample(i, Path("."), value, "", "species", "lineage", "local")
        for i, value in zip("abc", ("2018", "2019-02", "2020-03-04"))
    ]
    monkeypatch.setattr(
        "chronoclade.tree_stage.temporal._run_randomised_clock",
        lambda **kwargs: {"rate": 0.0001, "r_squared": 0.1},
    )
    args = dict(
        tree=tree,
        sequence_length=2000,
        samples=samples,
        observed_clock=clock,
        randomisations=20,
        randomisation_jobs=2,
        seed=1,
        output=tmp_path / "assessment.json",
        p_value_threshold=0.05,
    )
    result = run_clustered_assessment(**args)
    assert result["supported"] and result["temporal"]["cluster_count"] == 3
    assert result["temporal"]["p_value_r_squared"] == 1 / 21
    monkeypatch.setattr(
        "chronoclade.tree_stage.temporal._run_randomised_clock", lambda **kwargs: None
    )
    incomplete = run_clustered_assessment(**args)
    assert not incomplete["supported"] and incomplete["code"] == "not_assessed"


@pytest.mark.parametrize("flag_stage", ["raw", "corrected", "both"])
@pytest.mark.parametrize("policy", ["strict", "diagnostic"])
def test_real_coherence_checks_preserve_flagged_cohort(tmp_path, monkeypatch, flag_stage, policy):
    """Raw recombinant-like divergence must reach correction; persistent divergence fails.

    Mirrors the observed Bowers101 raw-only flag without pretending that this
    synthetic alignment reproduces the historical biology. Both screens use
    actual pairwise SNP counts; only external binaries are stubbed.
    """
    from chronoclade.lineage import LineageFiles, _run_core_phylogeny
    from chronoclade.metadata import Sample
    from chronoclade.evidence import read_alignment

    files = LineageFiles.in_directory(tmp_path)
    members = [Sample(i, tmp_path / f"{i}.fa", "", "", "sp", "lin", "local") for i in "abcde"]
    ids = [s.sample_id for s in members]
    calls = []
    def command(cmd, **kwargs):
        calls.append(cmd[0])
        def alignment(path, flagged):
            path.write_text("".join(f">{i}\n" + ("C" if flagged and i == "e" else "A") * 200 + "\n" for i in ids))
        if cmd[:2] == ["ska", "map"]:
            alignment(files.alignment, flag_stage in {"raw", "both"})
        elif cmd[0] == "iqtree":
            files.starting_tree.write_text("(" + ",".join(f"{i}:0.1" for i in ids) + ");")
        elif cmd[0] == "ClonalFrameML":
            files.tree.write_text(files.starting_tree.read_text())
            alignment(files.filtered_alignment, flag_stage in {"corrected", "both"})
    monkeypatch.setattr("chronoclade.lineage._run_command", command)
    if policy == "strict" and flag_stage in {"corrected", "both"}:
        with pytest.raises(WorkflowError, match="outlier"):
            _run_core_phylogeny(files, members, members[0], 1, False, "corrected", coherence_policy=policy)
        assert "ClonalFrameML" in calls
    else:
        tree, alignment, _ = _run_core_phylogeny(files, members, members[0], 1, False, "corrected", coherence_policy=policy)
        assert "ClonalFrameML" in calls
        assert set(read_alignment(alignment)) == set(ids)
        assert {tip.name for tip in Phylo.read(tree, "newick").get_terminals()} == set(ids)
        for label, path in (("raw", files.coherence_screen), ("corrected", files.clonal_coherence_screen)):
            audit = json.loads(path.with_suffix(".json").read_text())
            assert audit["sample_count"] == len(ids)
            assert len(audit["sample_statistics"]) == len(ids)
            assert audit["flagged_samples"] == (["e"] if flag_stage in {label, "both"} else [])
            if audit["flagged_samples"]:
                assert audit["flagged_statistics"][0]["median_raw_snps"] == 200


def test_coherence_policy_binds_resume_and_report(tmp_path, native_stubs):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    strict = run_tree(selection, output=tmp_path / "tree")
    diagnostic = run_tree(selection, output=tmp_path / "tree", coherence_policy="diagnostic")
    assert strict["fingerprint"] != diagnostic["fingerprint"]
    assert diagnostic["coherence_screening"]["policy"] == "diagnostic"
    assert diagnostic["parameters"]["coherence_policy"] == "diagnostic"
    assert len(native_stubs) == 2
    assert "Raw and corrected coherence audits" in (tmp_path / "tree" / diagnostic["report"]["path"]).read_text()


def test_diagnostic_warning_and_flagged_statistics_in_saved_report(tmp_path, native_stubs, monkeypatch):
    import chronoclade.lineage as lineage
    core = lineage._run_core_phylogeny
    def flagged(*args, **kwargs):
        result = core(*args, **kwargs)
        files = args[0]
        for path in (files.coherence_screen, files.clonal_coherence_screen):
            audit_path = path.with_suffix(".json")
            audit = json.loads(audit_path.read_text())
            audit["flagged_samples"] = ["d"]
            audit["flagged_statistics"] = [{"sample_id": "d", "median_raw_snps": 840, "robust_z": 12.3, "ratio_to_cohort_median": 6.7}]
            audit_path.write_text(json.dumps(audit))
        return result
    monkeypatch.setattr(lineage, "_run_core_phylogeny", flagged)
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    result = run_tree(selection, output=tmp_path / "tree", coherence_policy="diagnostic")
    assert result["selected_sample_ids"] == ["q", "c", "d"]
    assert result["qc_exclusions"] == []
    assert result["coherence_screening"]["corrected"]["flagged_statistics"][0]["median_raw_snps"] == 840
    report = (tmp_path / "tree" / result["report"]["path"]).read_text()
    assert "Coherence warning:" in report
    assert "legitimate cohort diversity" in report
    assert "Flagged genomes: d" in report


@pytest.mark.parametrize("mode", ["fast", "corrected"])
@pytest.mark.parametrize("seed", [None, 20261013])
def test_iqtree_seed_is_explicit_only_when_configured(tmp_path, monkeypatch, mode, seed):
    from chronoclade.lineage import LineageFiles, _run_core_phylogeny
    from chronoclade.metadata import Sample
    files = LineageFiles.in_directory(tmp_path)
    members = [Sample(i, tmp_path / f"{i}.fa", "", "", "sp", "lin", "local") for i in "abc"]
    commands = []
    def command(cmd, **kwargs):
        commands.append(cmd)
    monkeypatch.setattr("chronoclade.lineage._run_command", command)
    monkeypatch.setattr("chronoclade.lineage.screen_alignment", lambda *args, **kwargs: {"flagged_samples": []})
    monkeypatch.setattr("chronoclade.lineage.run_phipack_screen", lambda **kwargs: {})
    _run_core_phylogeny(files, members, members[0], 1, False, mode, seed=seed)
    iq = next(cmd for cmd in commands if cmd[0] == "iqtree")
    if seed is None:
        assert "-seed" not in iq
    else:
        assert iq[iq.index("-seed") + 1] == str(seed)


def test_native_tree_passes_seed_and_binds_resume(tmp_path, native_stubs, monkeypatch):
    import chronoclade.lineage as lineage
    original = lineage._run_core_phylogeny
    observed = []
    def core(*args, **kwargs):
        observed.append(kwargs["seed"])
        return original(*args, **kwargs)
    monkeypatch.setattr(lineage, "_run_core_phylogeny", core)
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    first = run_tree(selection, output=tmp_path / "tree", seed=20261013)
    second = run_tree(selection, output=tmp_path / "tree", seed=20261014)
    assert observed == [20261013, 20261014]
    assert first["fingerprint"] != second["fingerprint"]
    assert first["parameters"]["seed"] == 20261013
    assert first["parameters"]["phylogeny_seed"] == 20261013
    assert first["parameters"]["phylogeny_seed_policy"] == "explicit_IQ_TREE"


def test_shared_input_reference_survives_context_variation(tmp_path, native_stubs):
    """Higher-N50 alternative contexts cannot displace the shared input reference."""
    rows = []
    for index, ident in enumerate(("q", "c", "d", "e")):
        assembly = tmp_path / f"{ident}.fasta"
        assembly.write_text(">contig\n" + "A" * (2000 + index * 2000) + "\n")
        rows.append({"sample_id": ident, "origin": "local" if ident == "q" else "context", "species": "Klebsiella pneumoniae", "collection_date": "", "country": "UK", "assembly_reference": str(assembly)})
    dataset = write_dataset(PreparedDataset(samples=samples_from_records(rows)), tmp_path / "prepared")
    distance = tmp_path / "distances.csv"
    distance.write_text("sample_id_1,sample_id_2,distance,allele_differences,shared_called_loci\n" + "".join(f"{a},{b},{(j-i)/10},{j-i},10\n" for i,a in enumerate("qcde") for j,b in enumerate("qcde") if j > i))
    ensemble_path = write_selection_ensemble(rows[:1], rows[1:], {"paths": {"pairwise_distances": str(distance)}}, dataset_manifest=dataset, output=tmp_path / "selections", size=2, replicates=5, nearest_per_query=1)
    ensemble = json.loads(ensemble_path.read_text())
    results = []
    for index, entry in enumerate(ensemble["selections"]):
        result = run_tree(ensemble_path.parent / entry["path"], output=tmp_path / f"tree-{index}")
        results.append(result)
    assert len({tuple(result["selected_sample_ids"]) for result in results}) >= 2
    assert {result["reference_sample_id"] for result in results} == {"q"}
    assert len({result["reference"]["sha256"] for result in results}) == 1
    assert {result["parameters"]["reference_selection_policy"] for result in results} == {"highest_N50_common_input"}


def test_native_tree_requires_input_for_reference(tmp_path, native_stubs, monkeypatch):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    from chronoclade.tree_stage.workflow import selection_rows
    dataset, rows = selection_rows(selection, json.loads(selection.read_text()))
    for row in rows:
        row["role"] = "context"
    monkeypatch.setattr("chronoclade.tree_stage.workflow.selection_rows", lambda *args: (dataset, rows))
    with pytest.raises(WorkflowError, match="requires input genomes"):
        run_tree(selection, output=tmp_path / "tree")
    assert native_stubs == []


def test_strict_raw_review_warning_and_policy_version(tmp_path, native_stubs, monkeypatch):
    """A recombinant raw outlier can pass the unchanged corrected gate, as in Bowers101."""
    import chronoclade.lineage as lineage
    import chronoclade.tree_stage.workflow as workflow

    core = lineage._run_core_phylogeny

    def raw_flagged(*args, **kwargs):
        result = core(*args, **kwargs)
        audit_path = args[0].coherence_screen.with_suffix(".json")
        audit = json.loads(audit_path.read_text())
        audit["flagged_samples"] = ["d"]
        audit_path.write_text(json.dumps(audit))
        return result

    monkeypatch.setattr(lineage, "_run_core_phylogeny", raw_flagged)
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    first = run_tree(selection, output=tmp_path / "tree")
    assert first["selected_sample_ids"] == ["q", "c", "d"]
    assert first["coherence_screening"]["corrected"]["flagged_samples"] == []
    warning = first["coherence_screening"]["warning"]
    assert "recombination-filtered screen passed" in warning
    assert "Diagnostic policy" not in warning
    assert first["parameters"]["coherence_screen_policy_version"] == "raw_review_corrected_gate_v2"
    monkeypatch.setattr(workflow, "COHERENCE_SCREEN_POLICY_VERSION", "future_policy")
    changed = run_tree(selection, output=tmp_path / "tree")
    assert first["fingerprint"] != changed["fingerprint"]
    assert len(native_stubs) == 2
