import json
from pathlib import Path
import sqlite3
from types import SimpleNamespace

import pytest

from chronoclade import pasteur_profile_normalization as m
from chronoclade.artifacts import file_sha256
from chronoclade.datasets import DatasetError
from chronoclade.pathogenwatch import content_hash


@pytest.fixture
def reference(tmp_path, monkeypatch):
    root = tmp_path / "reference"
    root.mkdir()
    config = root / "query.json"
    config.write_text('{"fixture":"immutable config"}')
    index = root / "index"
    index.mkdir()
    metadata = index / "klebsiella_1_metadata.json"
    metadata.write_text(json.dumps({"path": "details.json.gz"}))
    (index / "details.json.gz").write_bytes(b"controlled details content")
    allele = index / "allele.db"
    with sqlite3.connect(allele) as c:
        c.execute("CREATE TABLE alleles (value TEXT)")
        c.execute("INSERT INTO alleles VALUES ('reference')")
    assignment = {
        "kind": "plincer",
        "command": ["fixture-tool"],
        "scheme": "scgMLST629_S",
        "scheme_version": "2026-01-01",
        "tool_version": "pinned",
        "canonical_loci": ["gene_1", "gene_2"],
    }
    for field in ("scheme_file", "profiles_file", "alleles_db", "loci_file"):
        path = root / field
        path.write_text(field)
        assignment[field] = path
    caller = {
        "index_dir": index,
        "scheme": "klebsiella_1",
        "scheme_version": "2026-01-01",
        "tool_version": "pinned",
        "command": ["fixture-caller"],
        "index_loci": ["gene 1", "gene 2"],
    }
    tools = {"organisms": {"Klebsiella pneumoniae": {"cgmlst": caller, "assignment": assignment}}}
    monkeypatch.setattr(m, "load_typing_config", lambda path: tools)

    def provenance(tool, fields):
        return {
            "scheme": tool["scheme"],
            "scheme_version": tool["scheme_version"],
            "version": tool["tool_version"],
            "command": tool["command"],
            "databases": {
                field: {"path": str(tool[field]), "sha256": file_sha256(tool[field])}
                for field in fields
            },
        }

    monkeypatch.setattr(m, "_provenance", provenance)
    return config, tools


def test_descriptor_contains_actual_content_and_ordered_panels(reference):
    config, tools = reference
    descriptor = m.reference_descriptor(config)
    assert descriptor["schema"] == m.REFERENCE_SCHEMA and descriptor["schema_version"] == 2
    caller = descriptor["resources"]["Klebsiella pneumoniae"]["caller"]
    assert caller["canonical_loci_sha256"] == content_hash(["gene_1", "gene_2"])
    assert caller["indexed_loci_sha256"] == content_hash(["gene 1", "gene 2"])
    assignment = descriptor["resources"]["Klebsiella pneumoniae"]["assignment"]
    assert set(assignment["databases"]) == {
        "scheme_file",
        "profiles_file",
        "alleles_db",
        "loci_file",
    }
    assert descriptor["sha256"] == content_hash(
        {k: v for k, v in descriptor.items() if k != "sha256"}
    )
    assert (
        m.resolve_reference(SimpleNamespace(parameters={"native_typing_reference": descriptor}))
        == config
    )


@pytest.mark.parametrize(
    "field",
    [
        "scheme_file",
        "profiles_file",
        "alleles_db",
        "loci_file",
        "allele.db",
        "scheme_metadata",
        "scheme_details",
    ],
)
def test_inplace_reference_mutation_under_unchanged_config_is_rejected(reference, field):
    config, tools = reference
    before = m.reference_descriptor(config)
    config_bytes = config.read_bytes()
    organism = tools["organisms"]["Klebsiella pneumoniae"]
    if field in organism["assignment"]:
        path = organism["assignment"][field]
    else:
        path = Path(
            before["resources"]["Klebsiella pneumoniae"]["caller"]["resources"][field]["path"]
        )
    if field == "scheme_metadata":
        data = json.loads(path.read_text())
        data["modified"] = "same locus panel"
        path.write_text(json.dumps(data))
    else:
        path.write_bytes(path.read_bytes() + b"changed")
    assert config.read_bytes() == config_bytes
    assert m.reference_descriptor(config) != before
    with pytest.raises(DatasetError, match="content changed"):
        m.resolve_reference(SimpleNamespace(parameters={"native_typing_reference": before}))


def test_legacy_descriptor_requires_explicit_rebinding_and_records_it(reference):
    config, tools = reference
    legacy = {
        "config_path": str(config),
        "config_sha256": file_sha256(config),
        "species": ["Klebsiella pneumoniae"],
    }
    ds = SimpleNamespace(parameters={"native_typing_reference": legacy})
    with pytest.raises(DatasetError, match="legacy/unverified"):
        m.resolve_reference(ds)
    assert m.resolve_reference(ds, config) == config
    current = m.reference_descriptor(config)
    assert (
        m.reference_rebinding(legacy, current, explicit=True)["kind"]
        == "legacy_config_only_to_content_bound"
    )
    assert ds.parameters["native_typing_reference"] == legacy
    assert m.reference_rebinding(legacy, current, explicit=False) is None


def test_unsupported_descriptor_version_and_hash_rejected(reference):
    config, _ = reference
    descriptor = m.reference_descriptor(config)
    for changed in [dict(descriptor, schema_version=999), dict(descriptor, sha256="0" * 64)]:
        with pytest.raises(DatasetError, match="legacy/unverified"):
            m.resolve_reference(SimpleNamespace(parameters={"native_typing_reference": changed}))


def test_real_readonly_sqlite_reads_shm_and_empty_wal_do_not_change_descriptor(reference):
    config, tools = reference
    index = tools["organisms"]["Klebsiella pneumoniae"]["cgmlst"]["index_dir"]
    path = index / "allele.db"
    before = m.reference_descriptor(config)
    with sqlite3.connect(path.as_uri() + "?mode=ro", uri=True) as connection:
        assert connection.execute("SELECT COUNT(*) FROM alleles").fetchone()[0] == 1
    path.with_name("allele.db-shm").write_bytes(b"transient lock bytes")
    path.with_name("allele.db-wal").touch()
    assert m.reference_descriptor(config) == before


@pytest.mark.parametrize("where", ["caller", "assignment"])
def test_uncheckpointed_wal_rejected(reference, where):
    config, tools = reference
    org = tools["organisms"]["Klebsiella pneumoniae"]
    path = (
        org["cgmlst"]["index_dir"] / "allele.db"
        if where == "caller"
        else org["assignment"]["alleles_db"]
    )
    path.with_name(path.name + "-wal").write_bytes(b"uncheckpointed")
    with pytest.raises(DatasetError, match="uncheckpointed SQLite writes"):
        m.reference_descriptor(config)


@pytest.mark.parametrize("stage", ["before_matched", "completed"])
def test_context_resume_rejects_resource_mutation_without_config_change(
    reference, tmp_path, monkeypatch, stage
):
    from dataclasses import asdict
    import chronoclade.profile_context_provider as provider
    from chronoclade.datasets import LocusCatalogue
    from chronoclade.prepare_stage import run_prepare

    config, tools = reference
    species = "Klebsiella pneumoniae"
    row = {
        "sample_id": "q",
        "species": species,
        "role": "input",
        "mlst_st": "39",
        "cgmlst_scheme": "test",
        "cgmlst_scheme_version": "v1",
        "cgmlst_profile": {"a": "1"},
        "cglin_scheme": "scgMLST629_S",
        "cglin_scheme_version": "lin-v1",
        "cglin_status": "complete",
        "cglin_raw": "0,0,1,2,3,4,5,6,7,8",
    }
    source, cats = tmp_path / "input.json", tmp_path / "catalogues.json"
    source.write_text(json.dumps([row]))
    cats.write_text(json.dumps([asdict(LocusCatalogue("test", "v1", ("a",), "fixture"))]))
    prepared = run_prepare(source, tmp_path / "prepared", catalogues=cats)
    events = []

    class Normalizer:
        def __init__(self, *args):
            self.audits = []

        def normalize(self, values):
            if stage == "before_matched" and any(r["role"] == "context" for r in values):
                raise RuntimeError("classification interruption")
            return [dict(r, cglin_scheme_version="current") for r in values]

    monkeypatch.setattr(m, "PasteurProfileNormalizer", Normalizer)
    candidate = dict(
        row,
        sample_id="PW_public",
        source_genome_id="public",
        role="context",
        origin="context",
        cgmlst_profile={},
    )

    class Client:
        base_url = "https://pathogen.watch"

        def supported_organisms(self):
            return [{"fullName": species, "organismId": "573"}]

        def freeze_catalogue(self, path, **kwargs):
            events.append("catalogue")
            value = {
                "rows": [candidate],
                "provenance": {
                    "complete": True,
                    "expected_count": 1,
                    "query": {"organismId": "573", "mlst": ["39"]},
                    "organism_capabilities": {"fullName": species, "organismId": "573"},
                    "base_url": self.base_url,
                },
            }
            path.write_text(json.dumps(value))
            return value

    def exports(values, *args, download_names=None, **kwargs):
        events.append("lineages" if download_names else "profiles")
        return values if download_names else [dict(r, cgmlst_profile={"a": "1"}) for r in values]

    monkeypatch.setattr(provider, "_grouped_analysis_exports", exports)
    monkeypatch.setattr(provider, "load_catalogue", lambda p: json.loads(p.read_text()))
    output = tmp_path / "context"
    if stage == "before_matched":
        with pytest.raises(RuntimeError, match="classification interruption"):
            provider.discover_profile_context(
                prepared.dataset_manifest, output, client=Client(), typing_config=config
            )
        assert not (tmp_path / "context_retrieval/matched_0.json").exists()
    else:
        result = provider.discover_profile_context(
            prepared.dataset_manifest, output, client=Client(), typing_config=config
        )
        assert result.dataset_manifest.is_file()
        assert (tmp_path / "context_retrieval/matched_0.json").is_file()
    prior_events = list(events)
    config_sha = file_sha256(config)
    tools["organisms"][species]["assignment"]["profiles_file"].write_text(
        "changed actual profile bytes"
    )
    assert file_sha256(config) == config_sha
    with pytest.raises(DatasetError, match="resume inputs|resume inputs changed"):
        provider.discover_profile_context(
            prepared.dataset_manifest, output, client=Client(), typing_config=config, resume=True
        )
    assert events == prior_events


def test_stable_directory_omits_only_known_sqlite_transients(tmp_path):
    directory = tmp_path / "db"
    directory.mkdir()
    db = directory / "alleles.db"
    with sqlite3.connect(db) as connection:
        connection.execute("CREATE TABLE refs (id INTEGER)")
    before = m._directory_reference(directory)
    db.with_name(db.name + "-shm").write_bytes(b"lock state")
    db.with_name(db.name + "-wal").touch()
    assert m._directory_reference(directory) == before
    (directory / "unknown-shm").write_bytes(b"not an ignorable SQLite sidecar")
    assert m._directory_reference(directory) != before
    db.with_name(db.name + "-wal").write_bytes(b"uncheckpointed")
    with pytest.raises(DatasetError, match="uncheckpointed"):
        m._directory_reference(directory)


def test_caller_only_config_fails_with_explicit_full_readiness_error(tmp_path, monkeypatch):
    from chronoclade.query_typing import QueryTypingError

    config = tmp_path / "caller-only.json"
    config.write_text('{"ready":false,"organisms":{"Escherichia coli":{"cgmlst":{}}}}')

    def full_only(*args, **kwargs):
        assert not kwargs.get("profiles_only", False)
        raise QueryTypingError("Typing databases are not ready")

    monkeypatch.setattr(m, "load_typing_config", full_only)
    with pytest.raises(DatasetError, match="requires full caller and lineage readiness"):
        m.reference_descriptor(config)


@pytest.mark.parametrize("legacy_binding", [True, False])
def test_reference_migration_is_portable_through_prepared_publication(
    reference, tmp_path, legacy_binding
):
    from copy import deepcopy
    from chronoclade.datasets import (
        LocusCatalogue, from_profile_records, load_dataset, write_dataset,
    )
    from chronoclade.datasets.model import check_json

    config, _ = reference
    current = m.reference_descriptor(config)
    previous = (
        {"config_path": str(config), "config_sha256": file_sha256(config),
         "species": ["Klebsiella pneumoniae"]}
        if legacy_binding else dict(current, sha256="0" * 64)
    )
    original_previous, original_current = deepcopy(previous), deepcopy(current)
    files_before = {
        str(path): (file_sha256(path), path.stat().st_size, path.stat().st_mtime_ns)
        for path in config.parent.rglob("*") if path.is_file()
    }
    migration = m.reference_rebinding(previous, current, explicit=True)
    assert migration["action_basis"] == "explicit --typing-config for fresh analysis"
    assert "authorization" not in migration
    assert migration["previous_descriptor_sha256"] == content_hash(previous)
    assert migration["current_descriptor_sha256"] == current["sha256"]
    parameters = {"profile_context": {"native_reference_migration": migration}}
    check_json({"parameters": parameters})
    dataset = from_profile_records(
        [{"sample_id": "query", "role": "input", "species": "Klebsiella pneumoniae",
          "cgmlst_scheme": "fixture", "cgmlst_scheme_version": "v1",
          "cgmlst_profile": {"gene_1": "1"}}],
        [LocusCatalogue("fixture", "v1", ("gene_1",), "fixture")],
        parameters=parameters,
    )
    manifest = write_dataset(dataset, tmp_path / "published")
    loaded = load_dataset(manifest)
    assert loaded.parameters == parameters
    assert previous == original_previous and current == original_current
    assert files_before == {
        str(path): (file_sha256(path), path.stat().st_size, path.stat().st_mtime_ns)
        for path in config.parent.rglob("*") if path.is_file()
    }
    # Keep the credential guard strict: the original reserved key is rejected.
    rejected = dict(migration, authorization="explicit --typing-config for fresh analysis")
    with pytest.raises(DatasetError, match="credentials must be referenced"):
        check_json({"parameters": {"profile_context": {"native_reference_migration": rejected}}})
