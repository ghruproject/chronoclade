"""Canonical biological identity is separate from execution-job/reference metadata."""
import hashlib
import json
from pathlib import Path
import sqlite3

import pytest

from chronoclade.datasets import DatasetError
from chronoclade import pasteur_profile_normalization as module


@pytest.fixture
def reference(tmp_path, monkeypatch):
    loci = [f"gene_{i}" for i in range(629)]
    index = tmp_path / "index"
    index.mkdir()
    (index / "klebsiella_1_metadata.json").write_text(json.dumps({
        "host": "pasteur", "hostPath": "pubmlst_klebsiella_seqdef", "schemeId": "18", "path": "details.json.gz",
    }))
    (index / "details.json.gz").write_bytes(b"fixture scheme details")
    digest = hashlib.sha1(b"acgt").hexdigest()
    with sqlite3.connect(index / "allele.db") as connection:
        connection.execute("CREATE TABLE alleles (hash TEXT, gene TEXT, st INTEGER, reverse INTEGER)")
        connection.executemany("INSERT INTO alleles VALUES (?,?,?,0)",
                               [(digest, locus, 1) for locus in loci])
    command = tmp_path / "plincer"
    (tmp_path / "python").touch()
    config_path = tmp_path / "config.json"
    config_path.write_text("{}")
    assignment = {"command": [str(command)], "canonical_loci": loci, "kind": "plincer",
                  "scheme": "scgMLST629_S", "scheme_version": "2026-10-05"}
    for field in ("scheme_file", "profiles_file", "alleles_db", "loci_file"):
        assignment[field] = tmp_path / field
        assignment[field].write_text(field)
    config = {"organisms": {"Klebsiella pneumoniae": {
        "cgmlst": {"index_dir": index, "scheme": "klebsiella_1", "scheme_version": "2026-10-05",
                   "tool_version": "fixture", "command": [str(command)], "index_loci": loci}, "assignment": assignment,
    }}}
    monkeypatch.setattr(module, "load_typing_config", lambda _: config)
    monkeypatch.setattr(module, "_provenance", lambda *_: {
        "databases": {"reference": {"sha256": "a" * 64}},
    })

    def assign(command, data, output):
        requests = json.loads(Path(command[2]).read_text())
        return {"assignments": [{"LINcode": ["0"] * 10, "matches": [
            {"st": "1", "LINcode": ["0"] * 10}
        ]} for _ in requests]}

    monkeypatch.setattr(module, "_run", assign)
    return config_path, index, digest, loci


def row(ident, **changes):
    result = {"sample_id": ident, "species": "Klebsiella pneumoniae",
              "cgmlst_scheme": "pathogenwatch:cgmlst-573-2", "cgmlst_scheme_version": "cgmlst-573-2",
              "cgmlst_profile": {"gene_0": "1"}, "cgmlst_loci": ["gene_0"],
              "cgmlst_locus_universe_complete": False,
              "cglin_raw": "0,0,0,0,0,0,0,0,0,0", "cglin_scheme_version": "old-job"}
    return dict(result, **changes)


def test_numbered_and_former_novel_alleles_share_verified_dna_identity(reference, tmp_path):
    config, _, digest, loci = reference
    native = row("new", cgmlst_scheme="klebsiella_1", cgmlst_scheme_version="2026-10-05")
    public = row("old", cgmlst_profile={"gene_0": digest}, source_genome_id="public-identity")
    normalizer = module.PasteurProfileNormalizer(config, tmp_path / "out")
    output = normalizer.normalize([native, public])
    assert output[0]["cgmlst_profile"] == output[1]["cgmlst_profile"] == {"gene_0": digest}
    assert output[0]["cgmlst_scheme"] == module.CANONICAL_SCHEME
    assert output[0]["cgmlst_loci"] == loci
    assert len(output[0]["cgmlst_missing_loci"]) == 628
    assert output[0]["cglin_database_sha256"] == output[1]["cglin_database_sha256"]
    assert output[1]["source_genome_id"] == "public-identity"
    assert public["cgmlst_scheme_version"] == "cgmlst-573-2"  # caller data unmodified
    audit = json.loads(Path(normalizer.audits[0]["path"]).read_text())
    assert audit["records"][1]["original_typing"]["cgmlst_profile"]["gene_0"] == digest
    assert audit["records"][0]["original_typing"]["cgmlst_profile"]["gene_0"] == "1"


def test_unknown_number_is_excluded_and_novel_remains_called(reference, tmp_path):
    config, _, _, _ = reference
    normalizer = module.PasteurProfileNormalizer(config, tmp_path / "out")
    output = normalizer.normalize([row("q", cgmlst_profile={"gene_0": "777"},
                                            cgmlst_novel_alleles={"gene_1": "b" * 40})])[0]
    assert output["cgmlst_profile"] == {"gene_1": "b" * 40}
    assert "gene_0" in output["cgmlst_missing_loci"]
    audit = json.loads(Path(normalizer.audits[0]["path"]).read_text())
    assert audit["records"][0]["unsupported_calls"] == [{
        "locus": "gene_0", "allele": "777",
        "reason": "numbered_allele_unavailable_in_verified_reference",
    }]


@pytest.mark.parametrize("changes", [
    {"cgmlst_scheme": "another-laboratory"},
    {"cgmlst_profile": {"unknown_locus": "1"}},
    {"cgmlst_locus_universe_complete": True},
])
def test_normalization_refuses_unverified_namespaces_and_locus_panels(reference, tmp_path, changes):
    normalizer = module.PasteurProfileNormalizer(reference[0], tmp_path / "out")
    with pytest.raises(DatasetError):
        normalizer.normalize([row("q", **changes)])


def test_conflicting_same_locus_allele_dna_is_rejected(reference, tmp_path):
    config, index, _, _ = reference
    with sqlite3.connect(index / "allele.db") as connection:
        connection.execute("INSERT INTO alleles VALUES (?, 'gene_0', 1, 0)", ("b" * 40,))
    normalizer = module.PasteurProfileNormalizer(config, tmp_path / "out")
    with pytest.raises(DatasetError, match="Conflicting allele sequence"):
        normalizer.normalize([row("q")])


def test_reference_locator_checks_exact_config_bytes(reference):
    config = reference[0]
    descriptor = module.reference_descriptor(config)
    from types import SimpleNamespace
    dataset = SimpleNamespace(parameters={"native_typing_reference": descriptor})
    assert module.resolve_reference(dataset) == config
    config.write_text('{"changed": true}')
    with pytest.raises(DatasetError, match="content changed"):
        module.resolve_reference(dataset)
    assert module.resolve_reference(dataset, config) == config


def test_progress_changes_only_messages_not_normalized_outputs(reference,tmp_path,monkeypatch,capsys):
    from chronoclade.progress import configure_progress
    config=reference[0]
    original=module._assignment_progress
    monkeypatch.setattr(module,'_assignment_progress',lambda *args,**kwargs:None)
    expected=module.PasteurProfileNormalizer(config,tmp_path/'no_messages').normalize([row('PRIVATE_SAMPLE')])
    monkeypatch.setattr(module,'_assignment_progress',original)
    configure_progress()
    actual=module.PasteurProfileNormalizer(config,tmp_path/'visible').normalize([row('PRIVATE_SAMPLE')])
    messages=capsys.readouterr().err.splitlines()
    assert messages==[
        '[typing] Assigning cgLIN to 1 profile using the configured Pasteur reference',
        '[typing] cgLIN assignment complete: 1 resolved, 0 provisional, 0 unassigned']
    assert actual==expected
    assert 'PRIVATE_SAMPLE' not in '\n'.join(messages)
    assert str(tmp_path) not in '\n'.join(messages)
    configure_progress(quiet=True)
    quiet=module.PasteurProfileNormalizer(config,tmp_path/'quiet').normalize([row('PRIVATE_SAMPLE')])
    assert quiet==expected and capsys.readouterr().err==''
