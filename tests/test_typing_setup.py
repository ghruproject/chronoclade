import json
from pathlib import Path

import pytest

from chronoclade import typing_setup


def test_species_aliases_and_validation():
    assert typing_setup._normalise_species("K. pneumoniae") == ["klebsiella"]
    assert typing_setup._normalise_species(["ecoli", "Escherichia coli"]) == ["ecoli"]
    assert typing_setup._normalise_species("both") == ["klebsiella", "ecoli"]
    with pytest.raises(typing_setup.TypingSetupError, match="Unsupported"):
        typing_setup._normalise_species("salmonella")


def test_setup_writes_truthful_not_ready_manifest(tmp_path, monkeypatch):
    monkeypatch.setattr(typing_setup.shutil, "which", lambda name: "/usr/bin/uv")
    monkeypatch.setattr(
        typing_setup,
        "_install_mlst",
        lambda root: {
            "status": "installed",
            "version": typing_setup.UPSTREAM["mlst"]["version"],
            "commit": typing_setup.UPSTREAM["mlst"]["commit"],
            "command": ["node", "index.js"],
        },
    )
    monkeypatch.setattr(
        typing_setup,
        "_install_python_tool",
        lambda uv, root, name, python: {
            "status": "installed",
            "version": typing_setup.UPSTREAM[name]["version"],
            "commit": typing_setup.UPSTREAM[name]["commit"],
            "command": [name],
        },
    )

    manifest_path, report = typing_setup.setup_typing(tmp_path, species="ecoli")

    assert manifest_path == tmp_path / "typing_setup.json"
    assert json.loads(manifest_path.read_text()) == report
    assert report["species"] == ["ecoli"]
    assert report["tools"]["hclink"]["status"] == "installed"
    assert report["run_ready"] is False
    assert report["query_config_path"] == str(tmp_path / "query_config.json")
    query_config = json.loads(Path(report["query_config_path"]).read_text())
    assert query_config["ready"] is False
    assert query_config["organisms"]["Escherichia coli"]["cgmlst"]["scheme_version"] is None
    assert query_config["organisms"]["Escherichia coli"]["assignment"]["loci_file"].endswith(
        "/databases/hclink/loci.json"
    )
    db = report["databases"]["ecoli"]
    assert db["cgmlst_ready"] is False
    assert db["lineage_ready"] is False
    assert db["cgmlst_access_scope"] == "unverified"
    assert "api key" in db["lineage_database_authentication"].lower()


def test_prepare_db_is_explicit_and_records_builds(tmp_path, monkeypatch):
    monkeypatch.setattr(typing_setup.shutil, "which", lambda name: "/usr/bin/uv")
    monkeypatch.setattr(
        typing_setup,
        "_install_mlst",
        lambda root: {
            "status": "installed",
            "version": typing_setup.UPSTREAM["mlst"]["version"],
            "commit": typing_setup.UPSTREAM["mlst"]["commit"],
            "command": ["node", "index.js"],
        },
    )
    monkeypatch.setattr(
        typing_setup,
        "_install_python_tool",
        lambda uv, root, name, python: {
            "status": "installed",
            "version": typing_setup.UPSTREAM[name]["version"],
            "commit": typing_setup.UPSTREAM[name]["commit"],
            "command": [name],
        },
    )
    monkeypatch.setattr(
        typing_setup,
        "_prepare_public_cgmlst",
        lambda uv, root, species, tools, **kwargs: {
            "status": "built",
            "scheme": typing_setup.SPECIES[species]["mlst_scheme"],
        },
    )

    _, report = typing_setup.setup_typing(tmp_path, species="klebsiella", prepare_db=True)

    assert report["cgmlst_builds"]["klebsiella"] == {
        "status": "built",
        "scheme": "klebsiella_1",
    }
    assert "only when their credential files were supplied" in report["database_preparation"]


def test_setup_routes_pasteur_secrets_to_authenticated_cgmlst_and_lin_builds(tmp_path, monkeypatch):
    monkeypatch.setattr(typing_setup.shutil, "which", lambda _name: "/usr/bin/uv")
    monkeypatch.setattr(
        typing_setup,
        "_install_mlst",
        lambda _root: {
            "status": "installed", "version": "8.0.0", "commit": "fixture",
            "command": ["node", "index.js"],
        },
    )
    monkeypatch.setattr(
        typing_setup,
        "_install_python_tool",
        lambda _uv, _root, name, _python: {
            "status": "installed", "version": "1.0.0", "commit": "fixture",
            "command": [name],
        },
    )
    observed = {}

    def prepare(_uv, _root, species, _tools, *, pasteur_secrets=None):
        observed["species"] = species
        observed["secrets"] = pasteur_secrets
        return {"status": "built", "access_mode": "authenticated"}

    monkeypatch.setattr(typing_setup, "_prepare_public_cgmlst", prepare)
    monkeypatch.setattr(
        typing_setup,
        "_build_plincer_db",
        lambda _root, _tool, secrets: {"status": "built", "secret_path": str(secrets)},
    )
    secret_path = tmp_path / "pasteur-client-key.txt"

    _, report = typing_setup.setup_typing(
        tmp_path / "setup", species="klebsiella", prepare_db=True,
        pasteur_secrets=secret_path,
    )

    assert observed == {"species": "klebsiella", "secrets": secret_path}
    assert report["cgmlst_builds"]["klebsiella"]["access_mode"] == "authenticated"
    assert report["cgmlst_builds"]["klebsiella_lineage"]["status"] == "built"


def test_missing_uv_is_reported_without_writing(tmp_path, monkeypatch):
    monkeypatch.setattr(typing_setup.shutil, "which", lambda _name: None)
    with pytest.raises(typing_setup.TypingSetupError, match="uv executable is unavailable"):
        typing_setup.setup_typing(tmp_path)
    assert not (tmp_path / "typing_setup.json").exists()


def test_hclink_builder_reads_key_from_file_without_putting_it_in_command(tmp_path, monkeypatch):
    key_file = tmp_path / "key.txt"
    key_file.write_text("private-fixture-key")
    tool = {"source": str(tmp_path / "source"), "environment": str(tmp_path / "venv")}

    def run(args, **kwargs):
        assert all("private-fixture-key" not in str(arg) for arg in args)
        assert kwargs["env"]["CHRONOCLADE_ENTEROBASE_API_KEY_FILE"] == str(key_file)
        script = Path(args[1]).read_text()
        assert "private-fixture-key" not in script
        db = Path(kwargs["env"]["CHRONOCLADE_HCLINK_DB"])
        (db / "metadata.json").write_text(json.dumps({"datestamp": "2026-01-01T12:00:00"}))
        (db / "loci.json").write_text(json.dumps({"genes": ["a", "b"]}))
        (db / "index.usearch").write_bytes(b"controlled-fixture-index")
        (db / "alleles.db").write_bytes(b"fixture-alleles")
        (db / "ST.txt.xz").write_bytes(b"fixture-clusters")
        return ""

    monkeypatch.setattr(typing_setup, "_run", run)
    report = typing_setup._build_hclink_db(tmp_path, tool, key_file)
    assert report["scheme_version"] == "2026-01-01"
    assert "private-fixture-key" not in json.dumps(report)


def test_plincer_credential_cache_is_private_and_outside_reference_database(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    (source / "scheme.toml").write_bytes(
        (Path(__file__).parent / "fixtures" / "plincer_scheme_upstream.toml").read_bytes()
    )
    secret = tmp_path / "secrets.json"
    secret.write_text(
        '{"pasteur":{"consumer":{"TOKEN":"client","TOKEN SECRET":"secret"},'
        '"user":{"TOKEN":"account","TOKEN SECRET":"password"}}}'
    )

    def run(args, **kwargs):
        cache = Path(args[args.index("--secrets-cache-file") + 1])
        assert cache.parent == tmp_path / "credential-cache"
        assert cache.parent.stat().st_mode & 0o777 == 0o700
        normalized = Path(args[args.index("--secrets-file") + 1])
        assert normalized.stat().st_mode & 0o777 == 0o600
        assert json.loads(normalized.read_text())["pasteur"]["consumer"]["TOKEN"] == "client"
        db = tmp_path / "databases" / "plincer"
        (db / "profiles.json.xz").write_bytes(b"fixture")
        (db / "alleles.sqlite").write_bytes(b"fixture")
        (db / "metadata.json").write_text(
            json.dumps({"last_updated": "2026-01-01", "genes": ["a"]})
        )
        return ""

    monkeypatch.setattr(typing_setup, "_run", run)
    report = typing_setup._build_plincer_db(
        tmp_path, {"source": str(source), "command": ["plincer"]}, secret
    )
    assert report["status"] == "built"
    assert report["scheme_version"] == "2026-01-01"


def test_pasteur_client_letter_is_normalised_without_exposing_values(tmp_path):
    letter = tmp_path / "pasteur-letter.txt"
    letter.write_bytes(b"Dear colleague,\nClient id: example-client\nClient secret: example-secret\n")
    destination = tmp_path / "private" / "pasteur.json"

    typing_setup._normalise_pasteur_secrets(letter, destination)

    assert destination.parent.stat().st_mode & 0o777 == 0o700
    assert destination.stat().st_mode & 0o777 == 0o600
    assert json.loads(destination.read_text()) == {
        "pasteur": {
            "consumer": {"TOKEN": "example-client", "TOKEN SECRET": "example-secret"}
        }
    }


def test_plincer_token_cache_is_valid_json_and_private(tmp_path):
    cache = tmp_path / "plincer.json"

    typing_setup._ensure_plincer_token_cache(cache)

    assert json.loads(cache.read_text()) == {}
    assert cache.stat().st_mode & 0o777 == 0o600


def test_plincer_fails_early_without_web_authorization(tmp_path):
    secrets = tmp_path / "pasteur.json"
    cache = tmp_path / "plincer.json"
    secrets.write_text('{"pasteur":{"consumer":{"TOKEN":"client","TOKEN SECRET":"secret"}}}')
    cache.write_text("{}")

    with pytest.raises(typing_setup.TypingSetupError, match="prior Pasteur access/session token"):
        typing_setup._check_plincer_authorization_inputs(secrets, cache)


def test_public_scheme_indexer_receives_its_own_local_scheme_manifest(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    interpreter = tmp_path / "venvs" / "typing_databases" / "bin" / "python"
    interpreter.parent.mkdir(parents=True)
    interpreter.write_text("fixture")
    (interpreter.parent / "download_schemes").write_text("fixture")
    monkeypatch.setattr(typing_setup, "_clone_at", lambda *args: source)
    monkeypatch.setattr(typing_setup, "_python_env", lambda *args: interpreter)
    other = tmp_path / "databases" / "mlst-schemes" / "ecoli_1-selection.json"
    other.parent.mkdir(parents=True)
    other.write_text("other-scheme-preserved")
    calls = []

    def run(args, **kwargs):
        calls.append(args)
        if "--output-schemes-file" in args:
            folder = Path(args[args.index("--output-dir") + 1]) / "pasteur" / "klebsiella_1"
            folder.mkdir(parents=True)
            (folder / "metadata.json").write_text(json.dumps({"genes": ["a"]}))
            selected = Path(args[args.index("--output-schemes-file") + 1])
            selected.write_text(
                json.dumps(
                    {
                        "schemes": [
                            {
                                "shortname": "klebsiella_1",
                                "scheme_id": "18",
                                "db_path": str(folder),
                                "last_updated": "2026-01-01",
                            }
                        ]
                    }
                )
            )
        if "index" in args:
            database = Path(
                next(arg.split("=", 1)[1] for arg in args if arg.startswith("--database="))
            )
            manifest = json.loads((database / "schemes.json").read_text())
            record = manifest["schemes"][0]
            assert record["shortname"] == "klebsiella_1"
            assert record["db_path"] == "."
            assert (database / record["db_path"] / "metadata.json").is_file()
            index = Path(next(arg.split("=", 1)[1] for arg in args if arg.startswith("--index=")))
            index.mkdir()
            (index / "klebsiella_1_metadata.json").write_text("fixture-index")
        return ""

    monkeypatch.setattr(typing_setup, "_run", run)
    report = typing_setup._prepare_public_cgmlst(
        "uv",
        tmp_path,
        "klebsiella",
        {"mlst": {"npm": "npm", "source": str(source), "command": ["node", "index.js"]}},
    )
    assert report["status"] == "built"
    assert report["access_mode"] == "public_unauthenticated"
    assert "cutoff" not in report["scope"]
    assert other.read_text() == "other-scheme-preserved"


def test_authenticated_pasteur_download_uses_isolated_cached_credentials(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    interpreter = tmp_path / "venvs" / "typing_databases" / "bin" / "python"
    interpreter.parent.mkdir(parents=True)
    interpreter.write_text("fixture")
    (interpreter.parent / "download_schemes").write_text("fixture")
    monkeypatch.setattr(typing_setup, "_clone_at", lambda *args: source)
    monkeypatch.setattr(typing_setup, "_python_env", lambda *args: interpreter)
    credential_root = tmp_path / "credential-cache"
    credential_root.mkdir()
    source_cache = credential_root / "plincer.json"
    source_cache.write_text(json.dumps({"pasteur": {"klebsiella:access": {
        "TOKEN": "fixture-access-token", "TOKEN SECRET": "fixture-access-secret"
    }}}))
    original_cache = source_cache.read_bytes()
    secrets = tmp_path / "pasteur-letter.txt"
    secrets.write_text(
        "Client id: fixture-client\nClient secret: fixture-secret\n"
        "Account password: must-not-be-forwarded\n"
    )
    calls = []

    def run(args, **kwargs):
        calls.append(args)
        if "--output-schemes-file" in args:
            folder = Path(args[args.index("--output-dir") + 1]) / "pasteur" / "klebsiella_1"
            folder.mkdir(parents=True)
            (folder / "metadata.json").write_text('{"genes":["locus-a"]}')
            Path(args[args.index("--output-schemes-file") + 1]).write_text(json.dumps({
                "schemes": [{
                    "shortname": "klebsiella_1", "scheme_id": "18",
                    "db_path": str(folder), "last_updated": "2026-10-05",
                }]
            }))
        if "index" in args:
            database = Path(next(a.split("=", 1)[1] for a in args if a.startswith("--database=")))
            manifest = json.loads((database / "schemes.json").read_text())
            assert manifest["schemes"][0]["scheme_id"] == "18"
            index = Path(next(a.split("=", 1)[1] for a in args if a.startswith("--index=")))
            index.mkdir()
            (index / "fixture-index").write_text("complete")
        return ""

    monkeypatch.setattr(typing_setup, "_run", run)
    report = typing_setup._prepare_public_cgmlst(
        "uv", tmp_path, "klebsiella",
        {"mlst": {"npm": "npm", "source": str(source), "command": ["node", "index.js"]}},
        pasteur_secrets=secrets,
    )

    command = calls[1]
    assert command[0] == str(interpreter)
    assert command[1:3] == ["-c", typing_setup._CACHED_PASTEUR_DOWNLOADER_ADAPTER]
    assert command[3] == "--scheme"
    assert "--secrets-file" in command and "--secrets-cache-file" in command
    assert all("must-not-be-forwarded" not in arg for arg in command)
    private = credential_root / "typing-databases"
    normalized = json.loads((private / "pasteur-secrets.json").read_text())
    assert set(normalized["pasteur"]) == {"consumer"}
    assert (private / "pasteur-secrets.json").stat().st_mode & 0o777 == 0o600
    assert (private / "pasteur-cache.json").stat().st_mode & 0o777 == 0o600
    assert (private).stat().st_mode & 0o777 == 0o700
    assert source_cache.read_bytes() == original_cache
    assert report["scheme_version"] == "2026-10-05"
    assert report["access_mode"] == "authenticated"
    assert typing_setup._cgmlst_access_scope(tmp_path, "klebsiella") == "authenticated"
    assert "fixture-access-token" not in json.dumps(report)
    selection = tmp_path / "databases" / "mlst-schemes" / "klebsiella_1-selection.json"
    selection.write_text(selection.read_text() + " ")
    assert typing_setup._cgmlst_access_scope(tmp_path, "klebsiella") == "unverified"


def test_authenticated_pasteur_download_requires_existing_access_cache(tmp_path):
    credential_root = tmp_path / "credential-cache"
    credential_root.mkdir()
    (credential_root / "plincer.json").write_text('{"pasteur":{}}')
    secrets = tmp_path / "secrets.txt"
    secrets.write_text("Client id: fixture-client\nClient secret: fixture-secret\n")

    with pytest.raises(typing_setup.TypingSetupError, match="authorize-pasteur.*--destination"):
        typing_setup._prepare_cached_pasteur_download(tmp_path, secrets)


def test_cached_pasteur_auth_adapter_is_scoped_and_never_starts_login(monkeypatch):
    import sys
    from types import ModuleType

    class FakeKeyCache:
        def __init__(self, keys):
            self.keys = keys

        def has_key(self, key_type, host, database=None):
            return self.keys.get((host, database, key_type), False)

        def can_authenticate(self, host):
            return self.has_key("consumer", host) and self.has_key("user", host)

    app_calls = []
    package = ModuleType("download_schemes")
    package.__path__ = []
    downloader = ModuleType("download_schemes.download_schemes")
    downloader.app = lambda: app_calls.append(True)
    keycache = ModuleType("download_schemes.keycache")
    keycache.KeyCache = FakeKeyCache
    monkeypatch.setitem(sys.modules, "download_schemes", package)
    monkeypatch.setitem(sys.modules, "download_schemes.download_schemes", downloader)
    monkeypatch.setitem(sys.modules, "download_schemes.keycache", keycache)
    exec(typing_setup._CACHED_PASTEUR_DOWNLOADER_ADAPTER, {"__name__": "__main__"})

    assert app_calls == [True]
    cache_only = FakeKeyCache({("pasteur", None, "consumer"): True,
                               ("pasteur", "klebsiella", "access"): True})
    no_authorization = FakeKeyCache({("pasteur", None, "consumer"): True,
                                     ("pasteur", None, "user"): True})
    other_host = FakeKeyCache({("pubmlst", None, "consumer"): True,
                               ("pubmlst", None, "user"): True})
    assert cache_only.can_authenticate("pasteur") is True
    assert no_authorization.can_authenticate("pasteur") is False
    assert other_host.can_authenticate("pubmlst") is True


def test_hclink_readiness_requires_allele_lookup_and_cluster_info(tmp_path):
    database = tmp_path / "databases" / "hclink"
    database.mkdir(parents=True)
    (database / "metadata.json").write_text('{"datestamp":"2026-01-01T12:00:00"}')
    (database / "loci.json").write_text('{"genes":["locus-a","locus-b"]}')
    (database / "index.usearch").write_bytes(b"fixture")
    state = typing_setup._database_state(tmp_path, "ecoli")
    assert not state["lineage_ready"]
    assert str(database / "alleles.db") in state["missing_resources"]
    assert str(database / "ST.txt.xz") in state["missing_resources"]
    for name in ("alleles.db", "ST.txt.xz"):
        (database / name).write_bytes(b"fixture")
    assert typing_setup._database_state(tmp_path, "ecoli")["lineage_ready"]


def test_nonempty_but_unindexed_cgmlst_directory_is_not_ready(tmp_path):
    scheme = typing_setup.SPECIES["ecoli"]["mlst_scheme"]
    selection = tmp_path / "databases" / "mlst-schemes"
    selection.mkdir(parents=True)
    (selection / f"{scheme}-selection.json").write_text(
        json.dumps({"schemes": [{"shortname": scheme, "last_updated": "2026-01-01"}]})
    )
    index = tmp_path / "databases" / "mlst-index" / scheme
    index.mkdir(parents=True)
    (index / "unrelated.txt").write_text("not a caller index")

    state = typing_setup._database_state(tmp_path, "ecoli")

    assert state["cgmlst_ready"] is False
    assert state["ready"] is False
    assert str(index) in state["missing_resources"]
