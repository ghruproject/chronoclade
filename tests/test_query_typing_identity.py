"""Computed lineage evidence must not erase independently supplied identities."""

import hashlib
import json

import pytest

from chronoclade import pasteur_profile_normalization, prepare_provider, query_typing
from chronoclade.datasets import DatasetError, load_dataset
from chronoclade.pathogenwatch import PathogenwatchClient
from chronoclade.prepare_stage import run_prepare

UUID = "fMEWxZ3K4bX9kXoKo9iRQW"


@pytest.fixture
def native(tmp_path, monkeypatch):
    assembly = tmp_path / "query.fa"
    assembly.write_text(">query\nACGT\n")
    digest = hashlib.sha256(assembly.read_bytes()).hexdigest()
    inputs = tmp_path / "assemblies.csv"
    inputs.write_text(
        "sample_id,assembly,species,mlst_st,source_genome_id,assembly_sha256\n"
        f"SRR1561305,{assembly},Klebsiella pneumoniae,258,{UUID},{digest}\n"
    )
    config = tmp_path / "config.json"
    config.write_text("{}\n")
    caller = {
        "command": ["caller"],
        "scheme": "klebsiella_1",
        "scheme_version": "fixture-v1",
        "index_dir": tmp_path,
    }
    assigner = {
        "kind": "plincer",
        "command": ["plincer"],
        "canonical_loci": ["a", "b", "c", "d"],
        "scheme_file": "scheme",
        "profiles_file": "profiles",
        "alleles_db": "alleles",
        "scheme": "scgMLST629_S",
        "scheme_version": "fixture-lin-v1",
    }
    monkeypatch.setattr(
        query_typing,
        "load_typing_config",
        lambda *a, **k: {
            "organisms": {"Klebsiella pneumoniae": {"cgmlst": caller, "assignment": assigner}}
        },
    )
    monkeypatch.setattr(
        pasteur_profile_normalization, "reference_descriptor", lambda *a: {"fixture": True}
    )
    monkeypatch.setattr(prepare_provider, "assign_existing_profiles", lambda *a: set())
    monkeypatch.setattr(
        query_typing, "_provenance", lambda *a: {"databases": {"index_dir": {"sha256": "a" * 64}}}
    )
    raw_lineage = {
        "LINcode": ["0", "0", "105", "6", "0", "0", "0", "0", "0", "0"],
        "matches": [{"LINcode": ["0", "0", "105", "6", "0", "0", "0", "0", "0", "0"]}],
        "cgST": "123",
    }

    def native_run(command, payload, destination):
        if command[0] == "caller":
            return {
                "scheme": "klebsiella_1",
                "schemeSize": 4,
                "genes": ["a", "b", "c", "d"],
                "code": "1__" + "b" * 40 + "_3",
                "alleles": {"d": [{"id": 3}, {"id": 4}]},
            }
        return raw_lineage

    monkeypatch.setattr(query_typing, "_run", native_run)
    result = run_prepare(
        inputs,
        tmp_path / "native",
        input_kind="assemblies",
        typing_config=config,
        enrich_metadata=False,
        client=PathogenwatchClient(api_key=""),
    )
    return inputs, result, raw_lineage, assigner


def test_plincer_assignment_and_ordinary_native_prepare_preserve_identity(native):
    inputs, result, raw, assigner = native
    assignment = query_typing._assignment_result(raw, assigner)
    assert "source_genome_id" not in assignment
    assert assignment["cglin_raw"] == "0_0_105_6_0_0_0_0_0_0"
    dataset = load_dataset(result.dataset_manifest)
    assert dataset.crosswalk[0]["source_genome_id"] == UUID
    retained = json.loads(
        (result.dataset_manifest.parent / "sources/provider/live_profiles.json").read_text()
    )
    row = retained["records"][0]
    assert row["source_genome_id"] == UUID
    assert row["cgmlst_profile"] == {"a": "1"}
    assert row["cgmlst_novel_alleles"] == {"c": "b" * 40}
    assert row["cgmlst_missing_loci"] == ["b"]
    assert row["cgmlst_ambiguous_loci"] == ["d"]
    assert row["cglin_raw"] == assignment["cglin_raw"]
    assert UUID in inputs.read_text()


def legacy_file(native, tmp_path, **changes):
    _, result, _, _ = native
    source = result.dataset_manifest.parent / "sources/provider/native_typing/query_typing.json"
    envelope = json.loads(source.read_text())
    envelope["assignments"][0].update(source_genome_id="", **changes)
    path = tmp_path / "legacy-query-typing.json"
    path.write_text(json.dumps(envelope))
    return path, envelope


def test_frozen_native_legacy_empty_identity_does_not_clear_input(native, tmp_path):
    path, envelope = legacy_file(native, tmp_path)
    original = path.read_bytes()
    loaded = query_typing.load_query_typing(path)
    assert "source_genome_id" not in loaded[0]
    assert loaded[0] == {
        k: v for k, v in envelope["assignments"][0].items() if k != "source_genome_id"
    }
    result = run_prepare(
        native[0],
        tmp_path / "imported",
        input_kind="assemblies",
        query_typing=path,
        enrich_metadata=False,
        client=PathogenwatchClient(api_key=""),
    )
    dataset = load_dataset(result.dataset_manifest)
    original_dataset = load_dataset(native[1].dataset_manifest)
    assert dataset.crosswalk[0]["source_genome_id"] == UUID
    assert dataset.profiles[0].profile("SRR1561305") == original_dataset.profiles[0].profile(
        "SRR1561305"
    )
    assert dataset.lineages == original_dataset.lineages
    assert path.read_bytes() == original


@pytest.mark.parametrize("mismatch", ["hash", "species"])
def test_frozen_native_import_still_rejects_assembly_or_species_mismatch(
    native, tmp_path, mismatch
):
    path, _ = legacy_file(
        native, tmp_path, **({"assembly_sha256": "c" * 64} if mismatch == "hash" else {})
    )
    inputs = native[0]
    if mismatch == "species":
        inputs = tmp_path / "wrong-species.csv"
        inputs.write_text(
            native[0].read_text().replace("Klebsiella pneumoniae", "Escherichia coli")
        )
    with pytest.raises(DatasetError, match="species/assembly hash differs"):
        run_prepare(
            inputs,
            tmp_path / "rejected",
            input_kind="assemblies",
            query_typing=path,
            enrich_metadata=False,
            client=PathogenwatchClient(api_key=""),
        )
    assert not (tmp_path / "rejected").exists()


def test_frozen_native_nonempty_provider_identity_is_rejected_not_replaced(native, tmp_path):
    path, envelope = legacy_file(native, tmp_path)
    envelope["assignments"][0]["source_genome_id"] = "different-provider-identity"
    path.write_text(json.dumps(envelope))
    before = native[0].read_bytes(), path.read_bytes()
    with pytest.raises(DatasetError, match="Native query typing cannot supply a provider identity"):
        run_prepare(
            native[0],
            tmp_path / "rejected",
            input_kind="assemblies",
            query_typing=path,
            enrich_metadata=False,
            client=PathogenwatchClient(api_key=""),
        )
    assert (native[0].read_bytes(), path.read_bytes()) == before
    assert not (tmp_path / "rejected").exists()
