"""Storage must preserve scientific bytes, native reuse and independent evidence."""
from concurrent.futures import ThreadPoolExecutor, TimeoutError
import hashlib
import json
from pathlib import Path
import sys
import threading

import numpy as np
import pytest
from typer.testing import CliRunner

from chronoclade import storage
from chronoclade.artifacts import file_sha256
from chronoclade.cli import app
from chronoclade.native_pipeline import _run_command


@pytest.mark.parametrize("codec", ["gzip", "zstd"])
def test_roundtrip_streaming_and_permanent_restore(tmp_path, codec):
    path = tmp_path / "evidence.fasta"
    content = (">sample\n" + "ACGT" * 10000 + "\n").encode()
    path.write_bytes(content)
    path.chmod(0o640)
    row = storage.pack_file(path, codec=codec)
    assert row["stored_bytes"] < len(content) / 10
    assert not path.exists()
    assert file_sha256(path) == hashlib.sha256(content).hexdigest()
    with storage.open_binary(path) as stream:
        assert stream.read(2) == content[:2]  # Exit must verify unread bytes too.
    with storage.open_text(path) as stream:
        assert next(stream) == ">sample\n"
    assert not path.exists()  # Streaming never materializes.
    assert storage.inspect(tmp_path, verify=True)[0]["bytes"] == len(content)
    assert list(storage.logical_files(tmp_path)) == [path]
    storage.restore(tmp_path)
    assert path.read_bytes() == content
    assert path.stat().st_mode & 0o777 == 0o640
    storage.pack_file(path)  # Reuse the checked archive, regardless of default codec.
    assert not path.exists()
    assert storage.inspect(tmp_path)[0]["codec"] == codec


@pytest.mark.parametrize("tamper", ["archive", "logical_checksum", "logical_size", "materialized"])
def test_corrupt_evidence_fails_closed(tmp_path, tamper):
    path = tmp_path / "data.npy"
    path.write_bytes(b"data" * 10000)
    storage.pack_file(path)
    receipt = storage.receipt_path(path)
    data = json.loads(receipt.read_text())
    if tamper == "archive":
        archive = storage._archive(path, data["codec"])
        content = bytearray(archive.read_bytes())
        content[-1] ^= 1
        archive.write_bytes(content)
    elif tamper == "logical_checksum":
        data["sha256"] = "0" * 64
        receipt.write_text(json.dumps(data))
    elif tamper == "logical_size":
        data["bytes"] = 1
        receipt.write_text(json.dumps(data))
    else:
        storage.restore(tmp_path)
        path.write_bytes(b"changed")
    with pytest.raises(storage.StorageError, match="mismatch"):
        storage.inspect(tmp_path, verify=True)
    if tamper != "materialized":
        with pytest.raises(storage.StorageError):
            storage.materialize(path)
        assert not path.exists()
        assert not list(tmp_path.glob(".cc-storage-*"))


def test_interrupted_commit_keeps_original(tmp_path, monkeypatch):
    path = tmp_path / "a.fasta"
    content = b"ACGT" * 10000
    path.write_bytes(content)
    def fail(*args):
        raise OSError("simulated disk full")
    monkeypatch.setattr(storage, "_write_record", fail)
    with pytest.raises(OSError, match="disk full"):
        storage.pack_file(path)
    assert path.read_bytes() == content
    assert not list(tmp_path.glob(".cc-storage-*"))


def test_session_evicts_only_unchanged_restorations(tmp_path):
    paths = [tmp_path / f"{i}.npy" for i in range(2)]
    for path in paths:
        np.save(path, np.zeros(1000))
        storage.pack_file(path)
    with storage.StorageSession():
        arrays = [np.load(storage.materialize(p), mmap_mode="r") for p in paths]
        assert all(isinstance(a, np.memmap) for a in arrays)
        paths[1].write_bytes(b"replacement output")
    assert not paths[0].exists()
    assert arrays[0].sum() == 0  # Existing mmap survives unlink on supported POSIX platforms.
    assert paths[1].read_bytes() == b"replacement output"


def test_native_reader_lease_prevents_concurrent_pack(tmp_path):
    path = tmp_path / "a.fasta"
    path.write_bytes(b"ACGT" * 10000)
    storage.pack_file(path)
    started = threading.Event()
    def pack():
        started.set()
        return storage.pack_file(path)
    with ThreadPoolExecutor() as executor:
        with storage.StorageSession():
            storage.materialize(path)
            future = executor.submit(pack)
            assert started.wait(3)
            with pytest.raises(TimeoutError):
                future.result(timeout=0.1)
            assert path.is_file()
        assert future.result(timeout=5)["status"] in {"packed", "already_packed"}
    assert not path.exists()


def test_idle_candidates_and_dry_run(tmp_path):
    for name in ("a.fasta", "matrix.npy", "public_0.json", "dataset.json", "report.html", "index.json"):
        (tmp_path / name).write_bytes(b"a" * 10000)
    (tmp_path / "linked.fasta").symlink_to(tmp_path / "a.fasta")
    rows = storage.pack(tmp_path, min_bytes=0, dry_run=True)
    assert {Path(r["path"]).name for r in rows} == {"a.fasta", "matrix.npy", "public_0.json"}
    assert not list(tmp_path.glob("*.cc-storage.json"))
    (tmp_path / "job.json").write_text('{"status":"running"}')
    with pytest.raises(storage.StorageError, match="running"):
        storage.pack(tmp_path, min_bytes=0)
    with pytest.raises(storage.StorageError):
        storage.pack(tmp_path, codec="invalid", dry_run=True)


def test_symlink_destinations_rejected(tmp_path):
    path = tmp_path / "a.fasta"
    path.write_bytes(b"A" * 10000)
    external = tmp_path / "external"
    external.write_text("keep")
    storage.receipt_path(path).symlink_to(external)
    with pytest.raises(storage.StorageError, match="symlink"):
        storage.pack_file(path)
    assert external.read_text() == "keep"
    assert path.is_file()


def test_native_cache_reuses_packed_inputs_and_outputs(tmp_path, monkeypatch):
    source, output, log = [tmp_path / name for name in ("source.fasta", "out.fasta", "copy.log")]
    source.write_bytes(b"ACGT" * 10000)
    command = [sys.executable, "-c", "import shutil,sys; shutil.copyfile(*sys.argv[1:])", str(source), str(output)]
    _run_command(command, log=log, expected=output, force=False, inputs=(source,))
    fingerprint = log.with_suffix(".fingerprint.json").read_bytes()
    storage.pack_file(source)
    storage.pack_file(output)
    def forbid(*args, **kwargs):
        pytest.fail("Native command executed instead of reusing checked output")
    monkeypatch.setattr("chronoclade.native_pipeline.subprocess.run", forbid)
    with storage.StorageSession():
        _run_command(command, log=log, expected=output, force=False, inputs=(source,))
        assert output.read_bytes() == b"ACGT" * 10000
    assert not source.exists() and not output.exists()
    assert log.with_suffix(".fingerprint.json").read_bytes() == fingerprint


def test_cli_storage_lifecycle(tmp_path):
    path = tmp_path / "large.fasta"
    content = b"A" * 2**20
    path.write_bytes(content)
    runner = CliRunner()
    for flags in (["pack", str(tmp_path), "--dry-run"], ["pack", str(tmp_path)],
                  ["status", str(tmp_path), "--json"], ["verify", str(tmp_path)],
                  ["restore", str(tmp_path)]):
        result = runner.invoke(app, ["storage", *flags])
        assert result.exit_code == 0, (result.output, result.exception)
    assert path.read_bytes() == content
    result = runner.invoke(app, ["--storage", "invalid", "version"])
    assert result.exit_code == 2
    assert storage._SESSION.get() is None


def test_compact_policy_only_successful_registered_outputs(tmp_path):
    path = tmp_path / "large.fasta"
    path.write_bytes(b"A" * 2**20)
    with pytest.raises(RuntimeError):
        with storage.StorageSession("compact"):
            storage.completed_output(tmp_path)
            raise RuntimeError("failed command")
    assert path.exists()
    with storage.StorageSession("compact"):
        storage.completed_output(tmp_path)
    assert not path.exists()


def test_packed_download_cache_never_redownloads_or_restores_shared_blob(tmp_path):
    from chronoclade.pathogenwatch_download import download_assemblies
    content = b">contig\n" + b"ACGT" * 10000 + b"\n"
    calls = []
    def fetch(source):
        calls.append(source)
        return content
    arguments = dict(rows=[{"source_genome_id": "source"}], cache_dir=tmp_path / "cache",
                     api_key="not-used", fetch=fetch)
    paths, _ = download_assemblies(output=tmp_path / "first", **arguments)
    storage.pack(tmp_path / "cache", min_bytes=0)
    paths, ledger = download_assemblies(output=tmp_path / "second", **arguments)
    assert calls == ["source"]
    assert ledger[0]["status"] == "cached"
    assert paths["source"].read_bytes() == content
    assert not list((tmp_path / "cache").glob("*.fasta"))


def test_packed_prepared_bundle_keeps_identity_and_mmap(tmp_path):
    from chronoclade.datasets import LocusCatalogue, from_profile_records, write_dataset, load_dataset
    rows = [{"sample_id": f"sample{i}", "origin": "local", "species": "example",
             "cgmlst_scheme": "test", "cgmlst_scheme_version": "v1",
             "cgmlst_profile": {"a": "1", "b": "2"}, "cgmlst_status": "resolved"}
            for i in range(100)]
    catalogue = LocusCatalogue("test", "v1", ("a", "b"), "test.tsv")
    manifest = write_dataset(from_profile_records(rows, [catalogue]), tmp_path / "prepared")
    before = load_dataset(manifest)
    identity = manifest.read_bytes()
    storage.pack(manifest.parent, min_bytes=0)
    packed = [Path(row["path"]) for row in storage.inspect(manifest.parent)]
    assert any(p.suffix == ".npy" for p in packed)
    with storage.StorageSession():
        after = load_dataset(manifest)
        assert after.samples == before.samples
        assert after.sample_ids == before.sample_ids
        assert isinstance(after.profiles[0].codes, np.memmap)
        assert np.array_equal(after.profiles[0].codes, before.profiles[0].codes)
    assert all(not p.exists() for p in packed)
    assert manifest.read_bytes() == identity


def test_tree_resume_after_compaction_keeps_manifest_and_native_result(tmp_path, monkeypatch):
    # Use the existing synthetic stage fixture; this checks storage, not inference.
    from test_tree_stage import selection_fixture, native_stubs
    from chronoclade.tree_stage.workflow import run_tree
    from chronoclade.stage_artifacts import load_tree_result
    fixture = native_stubs.__wrapped__(monkeypatch)
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    with storage.StorageSession():
        result = run_tree(selection, output=tmp_path / "tree")
    saved = (tmp_path / "tree/tree.json").read_bytes()
    storage.pack(tmp_path, min_bytes=0)
    assert storage.inspect(tmp_path)
    with storage.StorageSession():
        assert load_tree_result(tmp_path / "tree/tree.json") == result
        assert run_tree(selection, output=tmp_path / "tree") == result
    assert fixture == ["corrected"]
    assert (tmp_path / "tree/tree.json").read_bytes() == saved


def test_packed_binary_distances_keep_values_and_mmap(tmp_path):
    from test_cgmlst_scale import record
    from chronoclade.cgmlst.scale_distance import write_distances
    from chronoclade.matrix_distances import BinaryDistanceEvidence
    evidence, _, _ = write_distances(
        [record(str(i), [1, 2, 3, 4]) for i in range(30)], tmp_path, min_overlap=1)
    expected = evidence.get("0", "1")
    manifest = evidence.manifest_path.read_bytes()
    storage.pack(tmp_path, min_bytes=0)
    with storage.StorageSession():
        restored = BinaryDistanceEvidence(evidence.manifest_path)
        assert restored.get("0", "1") == expected
        assert all(isinstance(a, np.memmap) for a in restored.arrays.values())
    assert evidence.manifest_path.read_bytes() == manifest
    assert all(not Path(row["path"]).exists() for row in storage.inspect(tmp_path))


def test_compact_slurm_plan_workers_and_collectors(tmp_path, monkeypatch):
    from chronoclade.execution.jobs import SlurmResources
    from chronoclade.execution.slurm import write_slurm_ensemble, run_task, collect
    from test_tree_stage import selection_fixture, native_stubs
    native = native_stubs.__wrapped__(monkeypatch)
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    _, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    real_pack = storage.pack
    monkeypatch.setattr(storage, "pack", lambda root: real_pack(root, min_bytes=0))
    with storage.StorageSession("compact"):
        plan = write_slurm_ensemble("tree", ensemble, output=tmp_path / "jobs",
                                   results=tmp_path / "results", resources=SlurmResources(cpus=1))
    assert json.loads(plan.read_text())["storage"] == "compact"
    run_task(plan, 0)
    run_task(plan, 1)
    assert storage.inspect(tmp_path / "results")
    assert collect(plan) == 0
    run_task(plan, 0)  # Resume a compressed task without native work.
    assert len(native) == 2
    assert storage._SESSION.get() is None


def test_cli_prepare_compacts_and_validates(tmp_path, monkeypatch):
    from test_prepare_stage import frozen
    source, catalogue, _ = frozen.__wrapped__(tmp_path)
    real_pack = storage.pack
    monkeypatch.setattr(storage, "pack", lambda root: real_pack(root, min_bytes=0))
    args = ["--storage", "compact", "prepare", str(source), "--catalogues", str(catalogue),
            "--out", str(tmp_path / "prepared")]
    runner = CliRunner()
    first = runner.invoke(app, args)
    assert first.exit_code == 0, (first.output, first.exception)
    saved = (tmp_path / "prepared/prepare.json").read_bytes()
    assert storage.inspect(tmp_path / "prepared")
    from chronoclade.prepare_stage import _validate_previous_stage
    from chronoclade.datasets import load_dataset
    with storage.StorageSession():
        _validate_previous_stage(tmp_path / "prepared")
        assert load_dataset(tmp_path / "prepared/dataset.json").sample_ids == ("query", "missing")
    assert (tmp_path / "prepared/prepare.json").read_bytes() == saved
    assert storage._SESSION.get() is None


@pytest.mark.parametrize("backend", ["dense", "scale"])
def test_cgmlst_resumes_packed_checkpoint_without_recomputation(tmp_path, monkeypatch, backend):
    from test_cgmlst_stage import fixture_dataset
    import chronoclade.cgmlst.workflow as workflow
    _, manifest = fixture_dataset(tmp_path)
    output = tmp_path / "cgmlst"
    packager = workflow.write_profile_report
    def interrupt(*args, **kwargs):
        raise OSError("report interrupted")
    monkeypatch.setattr(workflow, "write_profile_report", interrupt)
    with pytest.raises(OSError, match="report interrupted"):
        with storage.StorageSession():
            workflow.run_cgmlst(manifest, output, bootstrap_replicates=0, backend=backend)
    numerical = next(output.glob("*/numerical_analysis.json"))
    before = numerical.read_bytes()
    storage.pack(tmp_path, min_bytes=0)
    def forbid(*args, **kwargs):
        pytest.fail("Packed cgMLST checkpoint was recomputed")
    monkeypatch.setattr(workflow, "analyse_profiles", forbid)
    monkeypatch.setattr("chronoclade.cgmlst.scale.analyse_profiles_scale", forbid)
    monkeypatch.setattr(workflow, "write_selection_ensemble", forbid)
    monkeypatch.setattr(workflow, "write_profile_report", packager)
    with storage.StorageSession():
        reused = workflow.run_cgmlst(manifest, output, bootstrap_replicates=0, backend=backend, resume=True)
    assert reused.manifest_path.is_file()
    assert numerical.read_bytes() == before


def test_large_window_zstd_roundtrip(tmp_path):
    # Small fixtures cannot exercise the actual 128 MiB frame/window path.
    path = tmp_path / "large.fasta"
    block = b"ACGT" * (1024**2 // 4)
    with path.open("wb") as stream:
        for _ in range(129):
            stream.write(block)
    original = file_sha256(path)
    storage.pack_file(path)
    assert file_sha256(path) == original
    assert storage.inspect(tmp_path, verify=True)[0]["bytes"] == 129 * 1024**2


@pytest.mark.integration
def test_real_native_phylogeny_reuses_packed_cache(tmp_path, monkeypatch):
    import random
    import shutil
    from chronoclade.native_pipeline import LineageFiles, _write_lineage_inputs, _run_core_phylogeny
    from chronoclade.metadata import Sample
    from chronoclade.evidence import read_alignment
    if not all(shutil.which(name) for name in ("ska", "iqtree", "Profile")):
        pytest.skip("SKA, IQ-TREE and PHIPack required")
    randomizer = random.Random(83)
    reference = "".join(randomizer.choices("ACGT", k=10000))
    members = []
    for i in range(4):
        path = tmp_path / f"sample{i}.fasta"
        sequence = list(reference)
        for position in range(i * 100, i * 100 + 10):
            sequence[position] = {"A": "C", "C": "G", "G": "T", "T": "A"}[sequence[position]]
        path.write_text(f">contig\n{''.join(sequence)}\n")
        members.append(Sample(f"sample{i}", path, "", "Unknown", "example", "example", "local"))
    files = LineageFiles.in_directory(tmp_path)
    _write_lineage_inputs(tmp_path, members)
    options = dict(recombination_method="phipack", phipack_window=1000, phipack_step=500)
    with storage.StorageSession():
        tree, alignment, summary = _run_core_phylogeny(files, members, members[0], 1, False, **options)
        expected_tree = tree.read_bytes()
        expected_alignment = read_alignment(alignment)
    storage.pack(tmp_path, min_bytes=0)
    def forbidden(*args, **kwargs):
        pytest.fail("Native tool reran after compaction")
    monkeypatch.setattr("chronoclade.native_pipeline.subprocess.run", forbidden)
    with storage.StorageSession():
        reused_tree, reused_alignment, reused_summary = _run_core_phylogeny(
            files, members, members[0], 1, False, **options)
        assert reused_tree.read_bytes() == expected_tree
        assert read_alignment(reused_alignment) == expected_alignment
        assert reused_summary == summary


def test_materialize_normalizes_aliases_and_rejects_packing_active_input(tmp_path, monkeypatch):
    path = tmp_path / "a.npy"
    np.save(path, np.zeros(1000))
    storage.pack_file(path)
    monkeypatch.chdir(tmp_path)
    with storage.StorageSession():
        assert storage.materialize("a.npy") == path
        assert storage.materialize(path) == path
        with pytest.raises(storage.StorageError, match="active native/mmap input"):
            storage.pack_file(path)
    assert not path.exists()


def test_prepare_input_reader_materializes_packed_assembly(tmp_path):
    from chronoclade.profile_inputs import _canonical
    path = tmp_path / "a.fasta"
    content = b">contig\n" + b"ACGT" * 10000 + b"\n"
    path.write_bytes(content)
    storage.pack_file(path)
    with storage.StorageSession():
        row = _canonical({"sample_id": "q", "assembly": str(path), "species": "example"},
                         origin="local", base=tmp_path)
        assert Path(row["assembly"]).read_bytes() == content
        assert row["assembly_sha256"] == hashlib.sha256(content).hexdigest()
    assert not path.exists()
