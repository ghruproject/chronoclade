"""Native Pathogenwatch adapters preserve unresolved alleles and lineage levels."""

import json

import pytest

from chronoclade.query_typing import (
    QueryTypingError,
    _assignment_result,
    _database_hash,
    load_query_typing,
    parse_cgmlst_result,
)


def profile(**updates):
    return {
        "scheme": "ecoli_1",
        "genes": ["a", "b", "c", "d"],
        "code": "1__" + "a" * 40 + "_3",
        "schemeSize": 4,
        "alleles": {"d": [{"id": 3}, {"id": 4}]},
        **updates,
    }


def indexed_scheme(index, scheme="ecoli_1", version="2026-01-01", genes=("a", "b")):
    import gzip
    import sqlite3

    index.mkdir(parents=True, exist_ok=True)
    for existing in index.iterdir():
        if existing.is_file():
            existing.unlink()
    for gene in genes:
        (index / f"{gene}.fa.gz").write_bytes(gzip.compress(f">{gene}_1\nACGT\n".encode()))
    details = {
        "shortname": scheme,
        "genes": list(genes),
        "schemeSize": len(genes),
        "allelePaths": {gene: f"{gene}.fa.gz" for gene in genes},
        "lengths": {gene: {"1": 4} for gene in genes},
        "alleleDictionary": {"ACGT": {"4": 1}},
        "profiles": {"1_1": "1"},
    }
    (index / "metadata.json.gz").write_bytes(gzip.compress(json.dumps(details).encode()))
    (index / f"{scheme}_metadata.json").write_text(
        json.dumps({"shortname": scheme, "lastUpdated": version, "path": "metadata.json.gz"})
    )
    with sqlite3.connect(index / "allele.db") as database:
        database.execute("CREATE TABLE alleles (hash TEXT, gene TEXT, st INTEGER, reverse INTEGER)")
        database.execute("CREATE INDEX ix_hash ON alleles (hash)")
        database.execute("INSERT INTO alleles VALUES ('hash', ?, 1, 0)", (genes[0],))


def hclink_database(path, datestamp="2026-01-01T12:00:00", genes=("a", "b")):
    path.mkdir(parents=True, exist_ok=True)
    (path / "metadata.json").write_text(json.dumps({"datestamp": datestamp}))
    (path / "loci.json").write_text(json.dumps({"genes": list(genes)}))
    for filename in ("index.usearch", "alleles.db", "ST.txt.xz"):
        (path / filename).write_bytes(b"controlled-fixture")


def test_known_novel_missing_ambiguous_are_separate():
    row = parse_cgmlst_result(profile(), "ecoli_1", "2026-01-01")
    assert row["cgmlst_profile"] == {"a": "1"}
    assert row["cgmlst_novel_alleles"] == {"c": "a" * 40}
    assert row["cgmlst_missing_loci"] == ["b"]
    assert row["cgmlst_ambiguous_loci"] == ["d"]
    assert row["cgmlst_loci"] == ["a", "b", "c", "d"]
    assert row["cgmlst_called_fraction"] == 0.5
    assert row["cgmlst_known_fraction"] == 0.25


def test_modern_hit_object_is_not_mistaken_for_duplicate():
    row = parse_cgmlst_result(
        profile(alleles={"a": {"id": 1, "contig": "c", "start": 1}}), "ecoli_1", "2026-01-01"
    )
    assert row["cgmlst_profile"] == {"a": "1", "d": "3"}


@pytest.mark.parametrize(
    "updates",
    [
        {"scheme": "other"},
        {"genes": ["a", "a", "c", "d"]},
        {"code": "1_2"},
        {"code": "1__nonsense_3"},
        {"schemeSize": 5},
    ],
)
def test_invalid_profile_rejected(updates):
    with pytest.raises(QueryTypingError):
        parse_cgmlst_result(profile(**updates), "ecoli_1", "2026-01-01")


def test_plincer_disagreeing_best_matches_keep_only_shared_prefix():
    config = {"kind": "plincer", "scheme": "scgMLST629_S", "scheme_version": "2026-01-01"}
    row = _assignment_result(
        {
            "LINcode": [1, 2, 3, 4, 5, "*", "*", "*", "*", "*"],
            "cgST": "*new",
            "matches": [{"LINcode": [1, 2, 3, 9, 5, 6, 7, 8, 9, 10]}],
        },
        config,
    )
    assert row["cglin_raw"] == "1_2_3_*_*_*_*_*_*_*"
    assert row["cglin_resolved_depth"] == 3
    assert row["cglin_provisional"]
    assert row["cglin_group_5"] == ""


def test_hclink_only_returns_native_supported_levels():
    row = _assignment_result(
        {"hierCC": [["d200", "42"], ["d1100", "9"]], "hierccDistance": 122},
        {"kind": "hclink", "scheme": "Enterobase_ecoli", "scheme_version": "2026-01-01"},
    )
    assert row["hiercc_codes"] == {"HC200": "42", "HC1100": "9"}
    assert "HC0" not in row["hiercc_codes"]


def test_no_hiercc_match_stays_unassigned():
    row = _assignment_result(
        {"hierCC": []},
        {"kind": "hclink", "scheme": "Enterobase_ecoli", "scheme_version": "2026-01-01"},
    )
    assert row["hiercc_status"] == "unassigned"


def test_database_hash_depends_on_file_names_and_bytes(tmp_path):
    (tmp_path / "one").write_text("one")
    original = _database_hash(tmp_path)
    (tmp_path / "one").rename(tmp_path / "two")
    assert _database_hash(tmp_path) != original


@pytest.mark.parametrize(
    "assignments",
    [
        [{"sample_id": "../unsafe", "assembly_sha256": "a" * 64}],
        [{"sample_id": "a", "assembly_sha256": "unknown"}],
        [{"sample_id": "a", "assembly_sha256": "a" * 64}] * 2,
    ],
)
def test_import_requires_safe_unique_ids_and_assembly_identity(tmp_path, assignments):
    path = tmp_path / "typing.json"
    path.write_text(json.dumps({"schema_version": 1, "assignments": assignments}))
    with pytest.raises(QueryTypingError):
        load_query_typing(path)


def test_hclink_blank_levels_are_not_copied_from_neighbour():
    row = _assignment_result(
        {"hierCC": [["HC0", ""], ["HC200", "42"]]},
        {"kind": "hclink", "scheme": "Enterobase_ecoli", "scheme_version": "2026-01-01"},
    )
    assert row["hiercc_codes"] == {"HC200": "42"}


@pytest.mark.parametrize(
    "update",
    [
        {"cgmlst_profile": {"a": "novel"}},
        {"cgmlst_loci": ["a", "a"]},
        {"cgmlst_profile": {"other": "1"}},
        {"hiercc_codes": {"HC200": ""}},
        {"cglin_raw": "1_?_3"},
    ],
)
def test_import_rejects_unsafe_assignments(tmp_path, update):
    path = tmp_path / "typing.json"
    row = {
        "sample_id": "query",
        "species": "Escherichia coli",
        "assembly_sha256": "a" * 64,
        "cgmlst_scheme": "ecoli_1",
        "cgmlst_scheme_version": "2026-01-01",
        "cgmlst_profile": {"a": "1"},
        "cgmlst_loci": ["a"],
        **update,
    }
    path.write_text(json.dumps({"schema_version": 1, "assignments": [row]}))
    with pytest.raises(QueryTypingError):
        load_query_typing(path)


def test_hclink_database_date_is_not_allowed_to_precede_caller(tmp_path):
    import sys
    from chronoclade.query_typing import load_typing_config

    index = tmp_path / "index"
    indexed_scheme(index, version="2026-01-02", genes=("a", "b"))
    db = tmp_path / "db"
    hclink_database(db, datestamp="2026-01-01 12:00:00")
    config = {
        "organisms": {
            "Escherichia coli": {
                "cgmlst": {
                    "command": [sys.executable],
                    "scheme": "ecoli_1",
                    "scheme_version": "2026-01-02",
                    "tool_version": "8.0.0",
                    "index_dir": "index",
                },
                "assignment": {
                    "command": [sys.executable],
                    "kind": "hclink",
                    "scheme": "Escherichia.cgMLSTv1",
                    "scheme_version": "2026-01-01",
                    "tool_version": "4.0.1",
                    "reference_db": "db",
                    "loci_file": "db/loci.json",
                },
            }
        }
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    with pytest.raises(QueryTypingError, match="at least as recent"):
        load_typing_config(path)
    config["organisms"]["Escherichia coli"]["cgmlst"]["scheme_version"] = "2026-01-01"
    indexed_scheme(index, version="2026-01-01", genes=("a", "b"))
    path.write_text(json.dumps(config))
    assert (
        load_typing_config(path)["organisms"]["Escherichia coli"]["assignment"]["reference_db"]
        == db
    )


def test_ready_config_rejects_empty_cgmlst_index(tmp_path):
    import sys
    from chronoclade.query_typing import load_typing_config

    index = tmp_path / "empty-index"
    index.mkdir()
    db = tmp_path / "hclink"
    hclink_database(db)
    config = {
        "ready": True,
        "organisms": {
            "Escherichia coli": {
                "cgmlst": {
                    "command": [sys.executable],
                    "scheme": "ecoli_1",
                    "scheme_version": "2026-01-01",
                    "tool_version": "8.0.0",
                    "index_dir": str(index),
                },
                "assignment": {
                    "command": [sys.executable],
                    "kind": "hclink",
                    "scheme": "Escherichia.cgMLSTv1",
                    "scheme_version": "2026-01-01",
                    "tool_version": "4.0.1",
                    "reference_db": str(db),
                    "loci_file": "hclink/loci.json",
                },
            }
        },
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))

    with pytest.raises(QueryTypingError, match="missing valid ecoli_1 caller metadata"):
        load_typing_config(path)


def test_ready_config_checks_caller_database_signature_and_locus_order(tmp_path):
    import sys
    from chronoclade.query_typing import load_typing_config

    index = tmp_path / "index"
    indexed_scheme(index, genes=("a", "b"))
    db = tmp_path / "hclink"
    hclink_database(db, genes=("b", "a"))
    config = {
        "ready": True,
        "organisms": {
            "Escherichia coli": {
                "cgmlst": {
                    "command": [sys.executable],
                    "scheme": "ecoli_1",
                    "scheme_version": "2026-01-01",
                    "tool_version": "8.0.0",
                    "index_dir": "index",
                },
                "assignment": {
                    "command": [sys.executable],
                    "kind": "hclink",
                    "scheme": "Escherichia.cgMLSTv1",
                    "scheme_version": "2026-01-01",
                    "tool_version": "4.0.1",
                    "reference_db": "hclink",
                    "loci_file": "hclink/loci.json",
                },
            }
        },
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))

    with pytest.raises(QueryTypingError, match="caller index loci do not match"):
        load_typing_config(path)


def test_not_ready_config_reports_setup_instructions(tmp_path):
    from chronoclade.query_typing import load_typing_config

    path = tmp_path / "config.json"
    path.write_text(json.dumps({"ready": False, "organisms": {}}))
    with pytest.raises(QueryTypingError, match="setup-typing"):
        load_typing_config(path)


def test_locus_order_must_match_assignment_database():
    from chronoclade.query_typing import _check_locus_order

    _check_locus_order(["a_S", "b_S"], ["a_S", "b_S"])
    with pytest.raises(QueryTypingError, match="locus order"):
        _check_locus_order(["b_S", "a_S"], ["a_S", "b_S"])


def test_species_case_and_underscore_aliases_use_canonical_config_and_cache(tmp_path, monkeypatch):
    import hashlib
    from chronoclade import query_typing
    from chronoclade.metadata import Sample

    assembly = tmp_path / "query.fa"
    assembly.write_text(">query\nACGT\n")
    (tmp_path / "config.json").write_text("{}")
    provenance = {"cgmlst": {}, "assignment": {}, "execution": {"identity_sha256":"fixture"}}
    monkeypatch.setattr(query_typing, "execution_provenance", lambda *a: {"identity_sha256":"fixture"})
    monkeypatch.setattr(
        query_typing,
        "load_typing_config",
        lambda path: {
            "organisms": {
                "Klebsiella pneumoniae": {"cgmlst": {}, "assignment": {"kind": "plincer"}}
            }
        },
    )
    monkeypatch.setattr(query_typing, "_provenance", lambda tool, fields: {})
    reused = {
        "sample_id": "query",
        "species": "KLEBSIELLA PNEUMONIAE",
        "assembly_sha256": hashlib.sha256(assembly.read_bytes()).hexdigest(),
        "typing_provenance": provenance,
    }
    sample = Sample("query", assembly, "2026", "UK", "Klebsiella_pneumoniae", "ST1", "query")
    rows = query_typing.type_query_assemblies(
        [sample],
        tmp_path / "config.json",
        tmp_path / "output",
        existing_assignments={"query": reused},
    )
    assert rows[0]["typing_reused"] is True
    assert rows[0]["species"] == "Klebsiella_pneumoniae"
    assert query_typing.load_query_typing(tmp_path / "output/query_typing.json") == rows


def test_import_rejects_independent_unsupported_taxon(tmp_path):
    from chronoclade.query_typing import normalise_typing_species

    with pytest.raises(QueryTypingError, match="Unsupported"):
        normalise_typing_species("Klebsiella oxytoca")
    path = tmp_path / "typing.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "assignments": [
                    {
                        "sample_id": "query",
                        "species": "Salmonella_enterica",
                        "assembly_sha256": "a" * 64,
                    }
                ],
            }
        )
    )
    with pytest.raises(QueryTypingError, match="Unsupported"):
        load_query_typing(path)


def test_profiles_only_uses_shared_caller_and_never_assigns_or_reuses_full_mode(tmp_path, monkeypatch):
    import sys
    from chronoclade import query_typing
    from chronoclade.metadata import Sample

    index = tmp_path / 'index'
    version = 'snapshot:' + 'a'*64
    indexed_scheme(index, version=version, genes=('a','b','c','d'))
    config = tmp_path / 'config.json'
    config.write_text(json.dumps({'ready':False, 'organisms':{'Escherichia coli':{
        'cgmlst':{'command':[sys.executable], 'scheme':'ecoli_1', 'scheme_version':version,
                  'tool_version':'8.0.0', 'index_dir':'index'},
        'assignment':{'kind':'hclink', 'reference_db':'absent'},
    }}}))
    with pytest.raises(QueryTypingError, match='not ready'):
        query_typing.load_typing_config(config)
    calls = []
    def run(command, data, destination):
        calls.append(command)
        assert '--scheme=ecoli_1' in command
        assert b'ACGT' in data
        return profile()
    monkeypatch.setattr(query_typing, '_run', run)
    assembly = tmp_path / 'query.fa'
    assembly.write_text('>query\nACGT\n')
    sample = Sample('query', assembly, '2020', 'Australia', 'Escherichia coli', 'ST95', 'query')
    rows = query_typing.type_query_assemblies([sample], config, tmp_path/'first', profiles_only=True)
    row = rows[0]
    assert row['typing_mode'] == 'profiles-only'
    assert row['typing_status'] == 'profiles_only'
    assert row['lineage_status'] == row['hiercc_status'] == 'not_requested'
    assert row['cgmlst_profile'] == {'a':'1'}
    assert row['cgmlst_missing_loci'] == ['b']
    assert row['cgmlst_ambiguous_loci'] == ['d']
    assert 'assignment' not in row['typing_provenance']
    assert 'hiercc_codes' not in row
    assert len(calls) == 1
    assert load_query_typing(tmp_path/'first/query_typing.json') == rows
    reused = query_typing.type_query_assemblies([sample], config, tmp_path/'reuse', profiles_only=True,
                                               existing_assignments={'query':row})
    assert reused[0]['typing_reused'] is True
    assert len(calls) == 1
    # Even with matching tool/data provenance, another execution mode is not reusable.
    old_mode = dict(row, typing_mode='full')
    query_typing.type_query_assemblies([sample], config, tmp_path/'new-mode', profiles_only=True,
                                      existing_assignments={'query':old_mode})
    assert len(calls) == 2
    indexed_scheme(index, version=version, genes=('a','b','c','d'))
    (index/'a.fa.gz').write_bytes(__import__('gzip').compress(b'>a_1\nTGCA\n'))
    latest = query_typing.type_query_assemblies([sample], config, tmp_path/'changed-index', profiles_only=True,
                                      existing_assignments={'query':row})
    assert len(calls) == 3
    current = latest[0]
    legacy = dict(current, typing_provenance=dict(current['typing_provenance']))
    legacy['typing_provenance'].pop('execution')
    query_typing.type_query_assemblies([sample], config, tmp_path/'legacy-unverified', profiles_only=True,
                                      existing_assignments={'query':legacy})
    assert len(calls) == 4
    changed = dict(current['typing_provenance']['execution'], identity_sha256='changed-source')
    monkeypatch.setattr(query_typing, 'execution_provenance', lambda *a: changed)
    query_typing.type_query_assemblies([sample], config, tmp_path/'changed-execution', profiles_only=True,
                                      existing_assignments={'query':current})
    assert len(calls) == 5
    from chronoclade import pasteur_profile_normalization
    monkeypatch.setattr(pasteur_profile_normalization, "reference_descriptor",
                        lambda *args, **kwargs: pytest.fail("profiles-only import must not require a full native reference"))
    # The ordinary preparation importer accepts exact profiles-only output,
    # retains its complete locus universe, and leaves lineage evidence empty.
    from chronoclade.prepare_stage import run_prepare
    from chronoclade.datasets import load_dataset
    inputs = tmp_path/'queries.csv'
    inputs.write_text('sample_id,assembly,species,lineage,mlst_st,country,collection_date\n'
                      f'query,{assembly},Escherichia coli,ST95,95,Australia,2020\n')
    prepared = run_prepare(inputs, tmp_path/'prepared', input_kind='assemblies',
                           query_typing=tmp_path/'first/query_typing.json', enrich_metadata=False)
    dataset = load_dataset(prepared.dataset_manifest)
    assert dataset.profiles[0].catalogue.complete is True
    assert dataset.profiles[0].catalogue.loci == ('a','b','c','d')
    assert all(value["assignment"] is None for value in dataset.lineages)
    assert dataset.lineages[0]["evidence"]["hiercc_status"] == "not_requested"


def test_profiles_only_import_cannot_claim_hiercc(tmp_path):
    path = tmp_path/'typing.json'
    path.write_text(json.dumps({'schema_version':1, 'assignments':[{
        'sample_id':'q', 'typing_mode':'profiles-only', 'typing_status':'profiles_only',
        'lineage_status':'not_requested', 'hiercc_codes':{'HC10':'10'},
    }]}))
    with pytest.raises(QueryTypingError, match='cannot claim a lineage'):
        load_query_typing(path)


def test_native_failure_is_inside_progress_without_sensitive_arguments(tmp_path, monkeypatch):
    from contextlib import contextmanager
    import subprocess
    from chronoclade import query_typing
    events = []
    @contextmanager
    def task(stage, label):
        events.append((stage, label, 'start'))
        try:
            yield
        except QueryTypingError:
            events.append((stage, label, 'failed'))
            raise
    monkeypatch.setattr(query_typing, 'progress_task', task)
    monkeypatch.setattr(query_typing.subprocess, 'run', lambda *a, **kw:
        subprocess.CompletedProcess(a[0], 1, stdout=b'', stderr=b'sensitive detail'))
    with pytest.raises(QueryTypingError, match='exit code 1') as error:
        query_typing._run(['/native/tool', '--private-placeholder'], b'', tmp_path/'result.json')
    assert events == [('typing','tool','start'), ('typing','tool','failed')]
    assert 'private-placeholder' not in str(error.value)
    assert 'sensitive detail' not in str(error.value)
