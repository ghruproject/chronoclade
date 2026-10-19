"""Explicit reference selection remains exact, input-bound and fingerprinted."""

import json

import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.errors import WorkflowError
from chronoclade.execution.ensemble import run_tree_ensemble
from chronoclade.execution.jobs import SlurmResources
from chronoclade.execution.slurm import write_slurm_ensemble, run_task, collect
from chronoclade.tree_stage.workflow import run_tree
from test_tree_stage import selection_fixture, native_stubs as _native_stubs


@pytest.fixture
def native_stubs(monkeypatch):
    return _native_stubs.__wrapped__(monkeypatch)


def test_explicit_reference_reaches_core_and_manifest(tmp_path, native_stubs, monkeypatch):
    import chronoclade.native_pipeline as native

    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    original = native._run_core_phylogeny
    chosen = []

    def core(files, members, reference, *args, **kwargs):
        chosen.append(reference.sample_id)
        return original(files, members, reference, *args, **kwargs)

    monkeypatch.setattr(native, "_run_core_phylogeny", core)
    default = run_tree(selection, output=tmp_path / "default")
    explicit = run_tree(selection, output=tmp_path / "explicit", reference_sample_id="q")
    assert chosen == ["q", "q"]
    assert explicit["reference_sample_id"] == "q"
    assert explicit["parameters"]["reference_selection_policy"] == "explicit_input"
    assert explicit["parameters"]["reference_sample_id"] == "q"
    assert default["parameters"]["reference_selection_policy"] == "highest_N50_common_input"
    assert default["fingerprint"] != explicit["fingerprint"]


@pytest.mark.parametrize("ident", ["missing", "c", ""])
def test_invalid_or_context_reference_rejected_before_native(tmp_path, monkeypatch, ident):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    monkeypatch.setattr(
        "chronoclade.tree_stage.workflow.tools", lambda *a: pytest.fail("No native tool work")
    )
    with pytest.raises(WorkflowError, match="exact selected input"):
        run_tree(selection, output=tmp_path / "tree", reference_sample_id=ident)
    assert not (tmp_path / "tree").exists()


def test_ensemble_and_slurm_reference_parity(tmp_path, native_stubs, monkeypatch):
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    _, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    normal = run_tree_ensemble(ensemble, output=tmp_path / "local", reference_sample_id="q")
    assert normal["parameters"]["reference_sample_id"] == "q"
    assert normal["parameters"]["reference_selection_policy"] == "explicit_input"
    plan = write_slurm_ensemble(
        "tree",
        ensemble,
        output=tmp_path / "jobs",
        results=tmp_path / "array",
        resources=SlurmResources(cpus=1),
        options={"reference_sample_id": "q"},
    )
    assert json.loads(plan.read_text())["options"]["reference_sample_id"] == "q"
    run_task(plan, 0)
    run_task(plan, 1)
    assert collect(plan) == 0
    saved = json.loads((tmp_path / "array/runs/selection-001/tree.json").read_text())
    assert saved["reference_sample_id"] == "q" and saved["parameters"]["reference_sample_id"] == "q"
    with pytest.raises(ValueError, match="exact selected input"):
        write_slurm_ensemble(
            "tree",
            ensemble,
            output=tmp_path / "badjobs",
            results=tmp_path / "badresults",
            resources=SlurmResources(cpus=1),
            options={"reference_sample_id": "c"},
        )


def test_normal_cli_reference_option(tmp_path, native_stubs):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    result = CliRunner().invoke(
        app,
        [
            "tree",
            str(selection),
            "--out",
            str(tmp_path / "cli"),
            "--reference-sample-id",
            "q",
            "--no-assess-temporal",
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads((tmp_path / "cli/tree.json").read_text())["reference_sample_id"] == "q"


def test_slurm_cli_reference_option(tmp_path):
    _, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    result = CliRunner().invoke(
        app,
        [
            "slurm-ensemble",
            "tree",
            str(ensemble),
            "--out",
            str(tmp_path / "jobs"),
            "--results",
            str(tmp_path / "results"),
            "--reference-sample-id",
            "q",
        ],
    )
    assert result.exit_code == 0, result.output
    assert (
        json.loads((tmp_path / "jobs/array.json").read_text())["options"]["reference_sample_id"]
        == "q"
    )


def test_explicit_reference_overrides_a_higher_n50_input(tmp_path, native_stubs):
    from chronoclade.datasets import PreparedDataset, samples_from_records, write_dataset
    from chronoclade.selection_manifest import write_selection_ensemble

    rows = []
    for ident, length in [("low", 2000), ("high", 4000), ("context", 2000)]:
        assembly = tmp_path / (ident + ".fa")
        assembly.write_text(">contig\n" + "A" * length + "\n")
        rows.append(
            {
                "sample_id": ident,
                "origin": "context" if ident == "context" else "local",
                "species": "Klebsiella pneumoniae",
                "assembly_reference": str(assembly),
            }
        )
    dataset = write_dataset(
        PreparedDataset(samples=samples_from_records(rows)), tmp_path / "prepared"
    )
    distances = tmp_path / "distances.csv"
    distances.write_text(
        "sample_id_1,sample_id_2,distance,allele_differences,shared_called_loci\nlow,high,.1,1,10\nlow,context,.1,1,10\nhigh,context,.1,1,10\n"
    )
    ensemble = write_selection_ensemble(
        rows[:2],
        rows[2:],
        {"paths": {"pairwise_distances": str(distances)}},
        dataset_manifest=dataset,
        output=tmp_path / "selection",
        size=1,
        replicates=1,
    )
    selection = ensemble.parent / "selections/selection-001.json"
    default = run_tree(selection, output=tmp_path / "default")
    explicit = run_tree(selection, output=tmp_path / "explicit", reference_sample_id="low")
    assert default["reference_sample_id"] == "high"
    assert explicit["reference_sample_id"] == "low"
