"""Install pinned Pathogenwatch typing tools and record database readiness.

This bootstrap installs software only by default. Typing databases can require
separate credentials and must never be reported ready until their files exist.
"""

from __future__ import annotations

import json
import hashlib
import os
import platform
import shutil
import subprocess
import tarfile
import urllib.error
import urllib.request
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.progress import progress, progress_task
from chronoclade.plincer_scheme import prepare_scheme, validate_scheme, SchemeValidationError


UPSTREAM = {
    "mlst": {
        "url": "https://github.com/pathogenwatch-oss/mlst.git",
        "commit": "4067e2aecbedc21ca2921fc0e6e149a0972119ff",
        "version": "8.0.0",
    },
    "plincer": {
        "url": "https://github.com/pathogenwatch-oss/klebsiella-lincodes.git",
        "commit": "c6d3d85e7484aee4d713429b78402b0590ac8fdc",
        "version": "7.0.0",
    },
    "hclink": {
        "url": "https://github.com/pathogenwatch-oss/hclink.git",
        "commit": "e8282cd6fc813b90d0b1d65f6598b11e7d637e8a",
        "version": "4.0.1",
    },
    "typing_databases": {
        "url": "https://github.com/pathogenwatch-oss/typing-databases.git",
        "commit": "d36e3c97e18bc0dc4d4eaa1136d1f15031f4384d",
        "version": "3.1.0",
    },
}

NODE_VERSION = "22.19.0"
NODE_SHA256 = {
    ("Darwin", "arm64"): "c59006db713c770d6ec63ae16cb3edc11f49ee093b5c415d667bb4f436c6526d",
    ("Darwin", "x86_64"): "3cfed4795cd97277559763c5f56e711852d2cc2420bda1cea30c8aa9ac77ce0c",
    ("Linux", "aarch64"): "d32817b937219b8f131a28546035183d79e7fd17a86e38ccb8772901a7cd9009",
    ("Linux", "x86_64"): "d36e56998220085782c0ca965f9d51b7726335aed2f5fc7321c6c0ad233aa96d",
}

SPECIES = {
    "klebsiella": {
        "name": "Klebsiella pneumoniae",
        "mlst_scheme": "klebsiella_1",
        "mlst_scheme_id": "18",
        "lineage_tool": "plincer",
        "lineage_scheme": "scgMLST629_S",
        "lineage_database": "Pasteur BIGSdb",
        "credential": "Pasteur account and consumer key/secret for current LIN database builds",
    },
    "ecoli": {
        "name": "Escherichia coli",
        "mlst_scheme": "ecoli_1",
        "mlst_scheme_id": "Escherichia.cgMLSTv1",
        "lineage_tool": "hclink",
        "lineage_scheme": "Escherichia.cgMLSTv1",
        "lineage_database": "EnteroBase HierCC",
        "credential": "scheme-specific EnteroBase API key for HierCC reference database builds",
    },
}


class TypingSetupError(ValueError):
    """A setup request or upstream installation failed safely."""


def _run(
    args: Sequence[str],
    *,
    cwd: Path | None = None,
    timeout: int = 1800,
    env: Mapping[str, str] | None = None,
    log_path: Path | None = None,
) -> str:
    try:
        result = subprocess.run(
            list(args),
            cwd=cwd,
            env=dict(env) if env is not None else None,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise TypingSetupError(f"Setup command failed: {Path(args[0]).name}") from exc
    if log_path is not None:
        # Opt-in only for credential-free native indexing, never auth/download tools.
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(result.stdout + "\n" + result.stderr)
    if result.returncode:
        # Tool output can contain remote-service details; do not echo it into reports.
        raise TypingSetupError(
            f"Setup command {Path(args[0]).name} exited with status {result.returncode}"
        )
    return result.stdout.strip()


def _clone_at(name: str, destination: Path) -> Path:
    source = UPSTREAM[name]
    if destination.exists():
        current = _run(["git", "-C", str(destination), "rev-parse", "HEAD"])
        if current != source["commit"]:
            raise TypingSetupError(
                f"Existing {name} checkout is not the pinned upstream revision: {destination}"
            )
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    _run(
        [
            "git",
            "clone",
            "--quiet",
            "--filter=blob:none",
            "--no-checkout",
            "--depth=1",
            source["url"],
            str(destination),
        ],
        timeout=900,
    )
    _run(
        [
            "git",
            "-C",
            str(destination),
            "fetch",
            "--quiet",
            "--depth=1",
            "origin",
            source["commit"],
        ],
        timeout=900,
    )
    _run(["git", "-C", str(destination), "checkout", "--quiet", "--detach", source["commit"]])
    return destination


def _python_env(uv: str, directory: Path, python: str) -> Path:
    interpreter = directory / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    if not interpreter.is_file():
        directory.parent.mkdir(parents=True, exist_ok=True)
        _run([uv, "venv", "--python", python, str(directory)], timeout=900)
    return interpreter


def _install_python_tool(uv: str, root: Path, name: str, python: str) -> dict[str, Any]:
    source = UPSTREAM[name]
    checkout = _clone_at(name, root / "sources" / name)
    env_dir = root / "venvs" / name
    interpreter = _python_env(uv, env_dir, python)
    package = f"git+{source['url']}@{source['commit']}"
    requirements = [uv, "pip", "install", "--python", str(interpreter), package]
    if name == "hclink":
        # usearch 2.21 is the last tested build available for the current ARM
        # runtime; newer wheels currently fail loading their native extension.
        requirements.extend(["usearch==2.21.0", "typing-extensions>=4.15,<5"])
    _run(requirements, timeout=1800)
    executable = env_dir / "Scripts" / f"{name}.exe" if os.name == "nt" else env_dir / "bin" / name
    if not executable.is_file():
        raise TypingSetupError(f"{name} installed without its expected command")
    _run([str(executable), "--help"], timeout=120)
    return {
        "status": "installed",
        "version": source["version"],
        "commit": source["commit"],
        "source": str(checkout),
        "environment": str(env_dir),
        "command": [str(executable)],
    }


def _install_node22(root: Path) -> tuple[Path, Path]:
    system = platform.system()
    machine = platform.machine()
    expected = NODE_SHA256.get((system, machine))
    if expected is None:
        raise TypingSetupError(f"No pinned Node.js {NODE_VERSION} runtime for {system}/{machine}")
    platform_name = "darwin" if system == "Darwin" else "linux"
    arch_name = "arm64" if machine in {"arm64", "aarch64"} else "x64"
    filename = f"node-v{NODE_VERSION}-{platform_name}-{arch_name}.tar.gz"
    runtime = root / "runtimes" / filename.removesuffix(".tar.gz")
    node = runtime / "bin" / "node"
    npm = runtime / "bin" / "npm"
    if node.is_file() and npm.exists() and _run([str(node), "--version"]) == f"v{NODE_VERSION}":
        return node, npm
    root.mkdir(parents=True, exist_ok=True)
    archive = root / filename
    url = f"https://nodejs.org/dist/v{NODE_VERSION}/{filename}"
    try:
        with urllib.request.urlopen(url, timeout=60) as response, archive.open("wb") as output:
            shutil.copyfileobj(response, output)
        digest = hashlib.sha256(archive.read_bytes()).hexdigest()
        if digest != expected:
            raise TypingSetupError("Downloaded Node.js archive failed its official SHA256 check")
        runtime.parent.mkdir(parents=True, exist_ok=True)
        shutil.rmtree(runtime, ignore_errors=True)
        with tarfile.open(archive, "r:gz") as package:
            package.extractall(runtime.parent, filter="data")
        archive.unlink(missing_ok=True)
    except TypingSetupError:
        archive.unlink(missing_ok=True)
        raise
    except (OSError, tarfile.TarError, urllib.error.URLError) as exc:
        archive.unlink(missing_ok=True)
        raise TypingSetupError(f"Could not install pinned Node.js {NODE_VERSION}") from exc
    if not node.is_file() or not npm.exists():
        raise TypingSetupError(f"Pinned Node.js {NODE_VERSION} archive has no node/npm executable")
    return node, npm


def _install_mlst(root: Path) -> dict[str, Any]:
    source = UPSTREAM["mlst"]
    checkout = _clone_at("mlst", root / "sources" / "mlst")
    node_path, npm_path = _install_node22(root)
    node = str(node_path)
    npm = str(npm_path)
    version = _run([node, "--version"])
    marker = checkout / "node_modules" / "better-sqlite3" / "package.json"
    if not marker.is_file():
        npm_env = dict(os.environ)
        npm_env["PATH"] = str(node_path.parent) + os.pathsep + npm_env.get("PATH", "")
        _run([npm, "ci", "--prefix", str(checkout)], timeout=1800, env=npm_env)
    return {
        "status": "installed",
        "version": source["version"],
        "commit": source["commit"],
        "source": str(checkout),
        "node_version": version,
        "npm": npm,
        "node_sha256": NODE_SHA256[(platform.system(), platform.machine())],
        "command": [node, str(checkout / "index.js")],
    }


def _normalise_species(species: str | Sequence[str]) -> list[str]:
    values = [species] if isinstance(species, str) else list(species)
    if not values or any(not isinstance(item, str) for item in values):
        raise TypingSetupError("Select klebsiella, ecoli, or both")
    normalised = []
    for value in values:
        key = " ".join(value.casefold().replace("_", " ").replace(".", " ").split())
        if key in {"both", "all"}:
            return ["klebsiella", "ecoli"]
        if key in {"klebsiella", "klebsiella pneumoniae", "k pneumoniae", "kp"}:
            canonical = "klebsiella"
        elif key in {"ecoli", "e coli", "escherichia coli"}:
            canonical = "ecoli"
        else:
            raise TypingSetupError(f"Unsupported native typing species: {value}")
        if canonical not in normalised:
            normalised.append(canonical)
    return normalised


def _database_state(root: Path, species: str) -> dict[str, Any]:
    spec = SPECIES[species]
    # The native MLST index is generated separately from downloaded scheme data.
    index = root / "databases" / "mlst-index" / spec["mlst_scheme"]
    scheme_data = root / "databases" / "mlst-schemes"
    lineage = root / "databases" / spec["lineage_tool"]
    if species == "klebsiella":
        required = [
            lineage / "scheme.published.toml",
            lineage / "profiles.json.xz",
            lineage / "alleles.sqlite",
            lineage / "metadata.json",
        ]
    else:
        required = [
            lineage / "metadata.json",
            lineage / "index.usearch",
            lineage / "loci.json",
            lineage / "alleles.db",
            lineage / "ST.txt.xz",
        ]
    try:
        from chronoclade.query_typing import validate_cgmlst_index

        caller_version = _database_version(root, species, "cgmlst")
        if caller_version is None:
            raise TypingSetupError("cgMLST scheme version is missing")
        validate_cgmlst_index(index, spec["mlst_scheme"], caller_version)
        cgmlst_ready = True
    except (OSError, ValueError):
        # Incomplete or malformed caller indexes are unready, not fatal to the
        # setup report.
        cgmlst_ready = False
    lineage_missing = [
        path
        for path in required
        if not path.is_file() or path.stat().st_size == 0
    ]
    metadata_path = lineage / "metadata.json"
    loci_path = metadata_path if species == "klebsiella" else lineage / "loci.json"
    try:
        metadata = json.loads(metadata_path.read_text())
        loci_metadata = json.loads(loci_path.read_text())
        genes = loci_metadata.get("genes") if isinstance(loci_metadata, dict) else None
        if (
            not isinstance(genes, list)
            or not genes
            or any(not isinstance(gene, str) or not gene for gene in genes)
            or len(set(genes)) != len(genes)
        ):
            raise ValueError
        if species == "klebsiella":
            if not isinstance(metadata.get("last_updated"), str) or not metadata["last_updated"].strip():
                raise ValueError
        else:
            date.fromisoformat(str(metadata["datestamp"]).split("T", 1)[0].split(" ", 1)[0])
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        if metadata_path not in lineage_missing:
            lineage_missing.append(metadata_path)
    if species == "klebsiella" and (lineage / "scheme.published.toml").is_file():
        try:
            validate_scheme(lineage / "scheme.published.toml")
        except (OSError, SchemeValidationError):
            lineage_missing.append(lineage / "scheme.published.toml")
    lineage_ready = not lineage_missing
    missing = ([] if cgmlst_ready else [str(index)]) + [str(path) for path in lineage_missing]
    return {
        "ready": cgmlst_ready and lineage_ready,
        "cgmlst_ready": cgmlst_ready,
        "lineage_ready": lineage_ready,
        "cgmlst_index": str(index),
        "cgmlst_scheme_data": str(scheme_data),
        "lineage_database": str(lineage),
        "missing_resources": missing,
        "lineage_database_authentication": spec["credential"],
        "cgmlst_access_scope": _cgmlst_access_scope(root, species),
    }


def _database_version(root: Path, species: str, tool: str) -> str | None:
    try:
        if tool == "cgmlst":
            file = (
                root
                / "databases"
                / "mlst-schemes"
                / f"{SPECIES[species]['mlst_scheme']}-selection.json"
            )
            records = json.loads(file.read_text()).get("schemes", [])
            record = next(
                item for item in records if item.get("shortname") == SPECIES[species]["mlst_scheme"]
            )
            value = record.get("last_updated")
        elif species == "klebsiella":
            value = json.loads((root / "databases" / "plincer" / "metadata.json").read_text()).get(
                "last_updated"
            )
        else:
            value = json.loads((root / "databases" / "hclink" / "metadata.json").read_text()).get(
                "datestamp"
            )
            value = str(value).split("T", 1)[0].split(" ", 1)[0]
        return str(value).strip() or None
    except (OSError, ValueError, KeyError, StopIteration, TypeError):
        return None


def _query_config(root: Path, species: list[str], tools: Mapping[str, Any]) -> dict[str, Any]:
    organisms: dict[str, Any] = {}
    for name in species:
        spec = SPECIES[name]
        caller_version = _database_version(root, name, "cgmlst")
        lineage_version = _database_version(root, name, "lineage")
        caller_tool = tools["mlst"]
        lineage_tool = tools[spec["lineage_tool"]]
        caller = {
            "command": caller_tool["command"],
            "scheme": spec["mlst_scheme"],
            "scheme_version": caller_version,
            "tool_version": f"cgps-mlst {caller_tool['version']}@{caller_tool['commit']}",
            "index_dir": str(root / "databases" / "mlst-index" / spec["mlst_scheme"]),
        }
        if name == "klebsiella":
            assignment = {
                "kind": "plincer",
                "command": lineage_tool["command"],
                "scheme": spec["lineage_scheme"],
                "scheme_version": lineage_version,
                "tool_version": f"plincer {lineage_tool['version']}@{lineage_tool['commit']}",
                "scheme_file": str(root / "databases" / "plincer" / "scheme.published.toml"),
                "profiles_file": str(root / "databases" / "plincer" / "profiles.json.xz"),
                "alleles_db": str(root / "databases" / "plincer" / "alleles.sqlite"),
                "loci_file": str(root / "databases" / "plincer" / "metadata.json"),
            }
        else:
            assignment = {
                "kind": "hclink",
                "command": lineage_tool["command"],
                "scheme": spec["lineage_scheme"],
                "scheme_version": lineage_version,
                "tool_version": f"hclink {lineage_tool['version']}@{lineage_tool['commit']}",
                "reference_db": str(root / "databases" / "hclink"),
                "loci_file": str(root / "databases" / "hclink" / "loci.json"),
            }
        organisms[spec["name"]] = {"cgmlst": caller, "assignment": assignment}
    return {
        "schema_version": 1,
        "ready": all(
            _database_state(root, name)["ready"]
            and _database_version(root, name, "cgmlst") is not None
            and _database_version(root, name, "lineage") is not None
            for name in species
        ),
        "organisms": organisms,
    }


def _cgmlst_access_scope(root: Path, species: str) -> str:
    """Return a verified download route only while it matches the saved selection."""
    spec = SPECIES[species]
    data_dir = root / "databases" / "mlst-schemes"
    selection = data_dir / f"{spec['mlst_scheme']}-selection.json"
    evidence_path = data_dir / f"{spec['mlst_scheme']}-download.json"
    try:
        evidence = json.loads(evidence_path.read_text())
        selected = json.loads(selection.read_text())
        records = selected.get("schemes", [])
        if (
            evidence.get("schema") != "chronoclade.typing-database-download"
            or evidence.get("schema_version") != 1
            or evidence.get("scheme") != spec["mlst_scheme"]
            or str(evidence.get("scheme_id", "")) != str(spec.get("mlst_scheme_id", ""))
            or len(records) != 1
            or records[0].get("shortname") != spec["mlst_scheme"]
            or str(records[0].get("scheme_id", ""))
            != str(spec.get("mlst_scheme_id", ""))
            or evidence.get("scheme_version") != records[0].get("last_updated")
            or evidence.get("selection_sha256") != file_sha256(selection)
            or evidence.get("access_mode") not in {"authenticated", "public_unauthenticated"}
        ):
            return "unverified"
        if evidence.get("access_method") == "frozen_enterobase_index":
            record = records[0]
            reference_path = root / "databases" / "enterobase" / spec["mlst_scheme_id"] / "reference.json"
            if (species != "ecoli" or evidence.get("access_mode") != "public_unauthenticated"
                    or evidence.get("reference_manifest_sha256") != file_sha256(reference_path)
                    or evidence.get("reference_sha256") != record.get("reference_sha256")
                    or json.loads(reference_path.read_text()).get("reference_sha256") != record.get("reference_sha256")
                    or evidence.get("scheme_version") != "snapshot:" + record["reference_sha256"]):
                return "unverified"
        elif evidence.get("source_commit") != UPSTREAM["typing_databases"]["commit"]:
            return "unverified"
        return evidence["access_mode"]
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return "unverified"


_CACHED_PASTEUR_DOWNLOADER_ADAPTER = """\
from download_schemes.download_schemes import app
from download_schemes.keycache import KeyCache

_upstream_can_authenticate = KeyCache.can_authenticate

def _can_authenticate_with_cached_access(self, host):
    if host == "pasteur":
        return self.has_key("consumer", host) and (
            self.has_key("access", host, "klebsiella")
            or self.has_key("session", host, "klebsiella")
        )
    return _upstream_can_authenticate(self, host)

KeyCache.can_authenticate = _can_authenticate_with_cached_access
app()
"""


def _prepare_cached_pasteur_download(root: Path, secrets_file: Path) -> dict[str, Path]:
    """Prepare an isolated downloader cache using a pre-authorized Pasteur token."""
    credential_root = root / "credential-cache"
    credential_root.mkdir(parents=True, exist_ok=True, mode=0o700)
    credential_root.chmod(0o700)
    source_cache = credential_root / "plincer.json"
    try:
        cache = json.loads(source_cache.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise TypingSetupError(
            "Authenticated Pasteur cgMLST download requires the existing Plincer token cache"
        ) from exc
    if not isinstance(cache, dict):
        raise TypingSetupError("Existing Plincer token cache must contain a JSON object")
    pasteur = cache.get("pasteur")
    if not isinstance(pasteur, dict):
        pasteur = {}
    def valid_cached_key(name: str) -> bool:
        key = pasteur.get(f"klebsiella:{name}")
        return (
            isinstance(key, dict)
            and isinstance(key.get("TOKEN"), str)
            and bool(key["TOKEN"].strip())
            and isinstance(key.get("TOKEN SECRET"), str)
            and bool(key["TOKEN SECRET"].strip())
        )

    if not (valid_cached_key("access") or valid_cached_key("session")):
        raise TypingSetupError(
            "Authenticated Pasteur cgMLST download requires a cached klebsiella access or "
            "session token; run authorize-pasteur with the same --destination first"
        )
    cached_tokens = {
        f"klebsiella:{name}": pasteur[f"klebsiella:{name}"]
        for name in ("access", "session")
        if valid_cached_key(name)
    }

    downloader_root = credential_root / "typing-databases"
    downloader_root.mkdir(parents=True, exist_ok=True, mode=0o700)
    downloader_root.chmod(0o700)
    secrets_path = downloader_root / "pasteur-secrets.json"
    cache_path = downloader_root / "pasteur-cache.json"
    _normalise_pasteur_secrets(secrets_file, secrets_path, consumer_only=True)
    _write_private_json(cache_path, {"pasteur": cached_tokens})
    return {"secrets_file": secrets_path, "cache_file": cache_path}


def _prepare_public_cgmlst(
    uv: str,
    root: Path,
    species: str,
    tools: dict[str, Any],
    *,
    pasteur_secrets: Path | None = None,
) -> dict[str, Any]:
    """Fetch the public snapshot or, for Klebsiella, a cached-auth Pasteur snapshot."""
    if species == "ecoli":
        return _prepare_enterobase_cgmlst(root, species, tools)
    source = UPSTREAM["typing_databases"]
    authenticated = species == "klebsiella" and pasteur_secrets is not None
    auth_files = (
        _prepare_cached_pasteur_download(root, pasteur_secrets) if authenticated else None
    )
    checkout = _clone_at("typing_databases", root / "sources" / "typing_databases")
    env_dir = root / "venvs" / "typing_databases"
    interpreter = _python_env(uv, env_dir, "3.11")
    package = f"git+{source['url']}@{source['commit']}"
    _run([uv, "pip", "install", "--python", str(interpreter), package], timeout=1800)
    downloader = (
        env_dir / "Scripts" / "download_schemes.exe"
        if os.name == "nt"
        else env_dir / "bin" / "download_schemes"
    )
    if not downloader.is_file():
        raise TypingSetupError("Typing database downloader installed without its command")

    spec = SPECIES[species]
    data_dir = root / "databases" / "mlst-schemes"
    data_dir.mkdir(parents=True, exist_ok=True)
    selected_file = data_dir / f"{spec['mlst_scheme']}-selection.json"
    download_command = [
            str(downloader),
            "--scheme",
            spec["mlst_scheme"],
            "--output-dir",
            str(data_dir),
            "--output-schemes-file",
            str(selected_file),
            "--fail-on-error",
        ]
    if authenticated:
        # The pinned upstream predicate ignores cached access/session tokens.
        # Invoke its unchanged entry point through a narrow adapter that reports
        # Pasteur as authenticated only when the exact cached Klebsiella token exists.
        download_command = [
            str(interpreter),
            "-c",
            _CACHED_PASTEUR_DOWNLOADER_ADAPTER,
            *download_command[1:],
            "--secrets-file",
            str(auth_files["secrets_file"]),
            "--secrets-cache-file",
            str(auth_files["cache_file"]),
        ]
    _run(
        download_command,
        cwd=checkout,
        timeout=7200,
    )
    try:
        selected = json.loads(selected_file.read_text())
        records = selected.get("schemes", [])
        if len(records) != 1 or records[0].get("shortname") != spec["mlst_scheme"]:
            raise ValueError
        expected_scheme_id = spec.get("mlst_scheme_id")
        if expected_scheme_id is not None and str(records[0].get("scheme_id")) != expected_scheme_id:
            raise ValueError
        scheme_path = Path(records[0]["db_path"])
        if not scheme_path.is_absolute():
            scheme_path = (data_dir / scheme_path).resolve()
        scheme_version = str(records[0].get("last_updated", "")).strip()
        if not scheme_version or not scheme_path.is_dir():
            raise ValueError
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise TypingSetupError("cgMLST download did not produce the expected validated scheme") from exc

    access_mode = "authenticated" if authenticated else "public_unauthenticated"
    write_json(
        data_dir / f"{spec['mlst_scheme']}-download.json",
        {
            "schema": "chronoclade.typing-database-download",
            "schema_version": 1,
            "scheme": spec["mlst_scheme"],
            "scheme_id": str(records[0].get("scheme_id", "")),
            "scheme_version": scheme_version,
            "access_mode": access_mode,
            "access_method": "cached_access_or_session_token" if authenticated else "public_route",
            "source_commit": source["commit"],
            "selection_sha256": file_sha256(selected_file),
        },
    )

    index_dir = _index_cgmlst(root, spec, records[0], scheme_path, tools)
    return {
        "status": "built",
        "scheme": spec["mlst_scheme"],
        "data_path": str(scheme_path),
        "index_path": str(index_dir),
        "source": source["url"],
        "source_commit": source["commit"],
        "access_mode": access_mode,
        "scope": (
            "Authenticated Pasteur access via a pre-existing cached token; exact downloaded "
            "scheme identity and version are recorded."
            if authenticated
            else "Unauthenticated public/redistributable records; exact downloaded scheme "
            "identity and version are recorded."
        ),
        "scheme_version": scheme_version,
    }


def _index_cgmlst(root, spec, record, scheme_path, tools):
    # The indexer reads schemes.json beneath --database and then resolves each
    # db_path against that same directory. The downloader's selected file uses
    # the original absolute download path, so give this individual scheme its
    # own manifest and a local path rather than altering other scheme records.
    index_manifest = dict(record, db_path=".")
    (scheme_path / "schemes.json").write_text(
        json.dumps({"schemes": [index_manifest]}, indent=2) + "\n"
    )
    index_dir = root / "databases" / "mlst-index" / spec["mlst_scheme"]
    index_dir.parent.mkdir(parents=True, exist_ok=True)
    npm = tools["mlst"]["npm"]
    index_log = root / "logs" / f"{spec['mlst_scheme']}-index.log"
    progress("setup", f"{spec['mlst_scheme']}: building pinned native caller index; "
             "this can take several minutes for large allele catalogues")
    with progress_task("setup", f"{spec['mlst_scheme']} native caller indexing",
                       details=f"index directory {index_dir}; native log {index_log}"):
        _run(
            [
                npm,
                "run",
                "index",
                "--",
                f"--scheme={spec['mlst_scheme']}",
                f"--index={index_dir}",
                f"--database={scheme_path}",
            ],
            cwd=Path(tools["mlst"]["source"]),
            timeout=7200,
            log_path=index_log,
            env={
                **os.environ,
                "PATH": str(Path(tools["mlst"]["command"][0]).parent)
                + os.pathsep
                + os.environ.get("PATH", ""),
            },
        )
    if not index_dir.is_dir() or not any(index_dir.iterdir()):
        raise TypingSetupError("cgps-mlst indexing completed without an index")
    progress("setup", f"{spec['mlst_scheme']}: native caller index written to {index_dir}")
    return index_dir


def _prepare_enterobase_cgmlst(root, species, tools):
    from chronoclade.enterobase_caller import prepare_caller_scheme

    spec = SPECIES[species]
    try:
        record = prepare_caller_scheme(root, spec["mlst_scheme_id"], spec["mlst_scheme"])
    except (OSError, ValueError) as error:
        raise TypingSetupError("Public EnteroBase allele reference could not be prepared") from error
    data_dir = root / "databases" / "mlst-schemes"
    selection = data_dir / f"{spec['mlst_scheme']}-selection.json"
    write_json(selection, {"schemes": [record]})
    write_json(data_dir / f"{spec['mlst_scheme']}-download.json", {
        "schema": "chronoclade.typing-database-download", "schema_version": 1,
        "scheme": spec["mlst_scheme"], "scheme_id": spec["mlst_scheme_id"],
        "scheme_version": record["last_updated"], "access_mode": "public_unauthenticated",
        "access_method": "frozen_enterobase_index", "selection_sha256": file_sha256(selection),
        "reference_sha256": record["reference_sha256"],
        "reference_manifest_sha256": file_sha256(record["reference_manifest"]),
    })
    index_dir = _index_cgmlst(root, spec, record, Path(record["db_path"]), tools)
    return {"status": "built", "scheme": spec["mlst_scheme"],
            "data_path": record["db_path"], "index_path": str(index_dir),
            "source": "https://enterobase.warwick.ac.uk/schemes/" + spec["mlst_scheme_id"] + "/",
            "reference_sha256": record["reference_sha256"],
            "scheme_version": record["last_updated"], "access_mode": "public_unauthenticated",
            "scope": "Frozen public allele definitions; HierCC reference readiness is checked separately."}


def _build_plincer_db(root: Path, tool: Mapping[str, Any], secrets_file: Path) -> dict[str, Any]:
    if not secrets_file.is_file():
        raise TypingSetupError("Pasteur credentials file does not exist or is not a file")
    checkout = Path(tool["source"])
    output = root / "databases" / "plincer"
    output.mkdir(parents=True, exist_ok=True)
    scheme_file = output / "scheme.published.toml"
    try:
        scheme_correction = prepare_scheme(checkout / "scheme.toml", scheme_file)
    except (OSError, SchemeValidationError) as exc:
        raise TypingSetupError("plincer assignment scheme validation failed") from exc
    profiles = output / "profiles.json.xz"
    metadata = output / "metadata.json"
    alleles = output / "alleles.sqlite"
    scratch = output / "scratch"
    credential_cache = root / "credential-cache"
    credential_cache.mkdir(parents=True, exist_ok=True, mode=0o700)
    credential_cache.chmod(0o700)
    normalized_secrets = credential_cache / "pasteur-secrets.json"
    _normalise_pasteur_secrets(secrets_file, normalized_secrets)
    token_cache = credential_cache / "plincer.json"
    _ensure_plincer_token_cache(token_cache)
    _check_plincer_authorization_inputs(normalized_secrets, token_cache)
    command = [
        *tool["command"],
        "build",
        "--scheme-file",
        str(scheme_file),
        "--profiles-file",
        str(profiles),
        "--metadata-file",
        str(metadata),
        "--alleles-db",
        str(alleles),
        "--secrets-file",
        str(normalized_secrets),
        "--secrets-cache-file",
        str(token_cache),
        "--host-config-file",
        str(checkout / "host_config.json"),
        "--scratch-dir",
        str(scratch),
        "--clean",
    ]
    _run(command, cwd=output, timeout=12 * 3600)
    try:
        source_metadata = json.loads(metadata.read_text())
        scheme_version = str(source_metadata["last_updated"]).strip()
        if not scheme_version:
            raise ValueError
        ready = all(path.is_file() for path in (profiles, metadata, alleles, scheme_file))
        if not ready:
            raise ValueError
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise TypingSetupError("plincer completed without its required database files") from exc
    return {
        "status": "built",
        "scheme": SPECIES["klebsiella"]["lineage_scheme"],
        "scheme_version": scheme_version,
        "database_path": str(output),
        "scheme_correction": scheme_correction,
    }


def _normalise_pasteur_secrets(
    source: Path, destination: Path, *, consumer_only: bool = False
) -> None:
    """Convert a Pasteur-issued client-key letter or upstream JSON to plincer JSON."""
    try:
        raw = source.read_bytes()
        try:
            text = raw.decode("utf-8-sig")
        except UnicodeDecodeError:
            text = raw.decode("cp1252")
    except (OSError, UnicodeError) as exc:
        raise TypingSetupError("Cannot read Pasteur credentials file") from exc

    try:
        value = json.loads(text)
    except json.JSONDecodeError:
        import re

        fields = {}
        for line in text.splitlines():
            match = re.match(r"\s*(Client\s+(?:id|secret))\s*:?\s+(.+?)\s*$", line, re.I)
            if match:
                fields[match.group(1).casefold()] = match.group(2)
        client_id = fields.get("client id")
        client_secret = fields.get("client secret")
        if not client_id or not client_secret:
            raise TypingSetupError(
                "Pasteur credentials must be upstream JSON or a letter with Client id and Client secret"
            )
        value = {"pasteur": {"consumer": {"TOKEN": client_id, "TOKEN SECRET": client_secret}}}

    try:
        consumer = value["pasteur"]["consumer"]
        if not isinstance(consumer.get("TOKEN"), str) or not consumer["TOKEN"].strip():
            raise ValueError
        if not isinstance(consumer.get("TOKEN SECRET"), str) or not consumer["TOKEN SECRET"].strip():
            raise ValueError
    except (TypeError, KeyError, ValueError, AttributeError) as exc:
        raise TypingSetupError(
            "Pasteur JSON credentials require pasteur.consumer TOKEN and TOKEN SECRET"
        ) from exc
    if consumer_only:
        # The authenticated downloader uses its separate token cache. Never copy
        # account passwords or unrelated OAuth cache entries into that route.
        value = {"pasteur": {"consumer": consumer}}

    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    destination.parent.chmod(0o700)
    _write_private_json(destination, value)


def _write_private_json(path: Path, value: Mapping[str, Any]) -> None:
    """Write credential JSON without a permissive creation window or symlink follow."""
    flags = os.O_CREAT | os.O_WRONLY | os.O_TRUNC | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags, 0o600)
        os.fchmod(descriptor, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2)
            stream.write("\n")
    except OSError as exc:
        raise TypingSetupError("Cannot safely write protected Pasteur credential cache") from exc


def _ensure_plincer_token_cache(path: Path) -> None:
    if path.exists():
        if path.is_symlink() or not path.is_file():
            raise TypingSetupError("Pasteur token cache must be a regular file")
        path.chmod(0o600)
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise TypingSetupError("Pasteur token cache is unreadable or invalid") from exc
        if not isinstance(value, dict):
            raise TypingSetupError("Pasteur token cache is invalid")
    else:
        _write_private_json(path, {})


def _check_plincer_authorization_inputs(secrets_file: Path, token_cache: Path) -> None:
    """Avoid starting a database build without a web-login credential or prior access token."""
    try:
        secrets = json.loads(secrets_file.read_text(encoding="utf-8"))
        cache = json.loads(token_cache.read_text(encoding="utf-8"))
        user = secrets.get("pasteur", {}).get("user", {})
        cached = cache.get("pasteur", {})
        access = cached.get("klebsiella:access", {})
        session = cached.get("klebsiella:session", {})
        has_user = bool(user.get("TOKEN") and user.get("TOKEN SECRET"))
        has_access = any(
            isinstance(token, dict) and token.get("TOKEN") and token.get("TOKEN SECRET")
            for token in (access, session)
        )
    except (OSError, ValueError, TypeError, AttributeError) as exc:
        raise TypingSetupError("Pasteur authorization files are invalid") from exc
    if not has_user and not has_access:
        raise TypingSetupError(
            "Pasteur client credentials alone are insufficient: plincer also needs "
            "a prior Pasteur access/session token or account credentials to build the database. "
            "Run authorize-pasteur with the same --destination first."
        )


def _build_hclink_db(root: Path, tool: Mapping[str, Any], key_file: Path) -> dict[str, Any]:
    if not key_file.is_file():
        raise TypingSetupError("EnteroBase API key file does not exist or is not a file")
    try:
        key = key_file.read_text().strip()
    except OSError as exc:
        raise TypingSetupError("Cannot read EnteroBase API key file") from exc
    if not key:
        raise TypingSetupError("EnteroBase API key file is empty")
    output = root / "databases" / "hclink"
    output.mkdir(parents=True, exist_ok=True)
    script = root / "scripts" / "build_hclink_db.py"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text(
        "from pathlib import Path\n"
        "import os\n"
        "import hclink.hclink as hclink\n"
        "download_profiles = hclink.download_profiles\n"
        "def capture_loci(url, directory):\n"
        "    profiles = download_profiles(url, directory)\n"
        "    genes = hclink.extract_genes(profiles)\n"
        "    import json\n"
        "    Path(directory, 'loci.json').write_text(json.dumps({'genes': genes}) + '\\n')\n"
        "    return profiles\n"
        "hclink.download_profiles = capture_loci\n"
        "key = Path(os.environ['CHRONOCLADE_ENTEROBASE_API_KEY_FILE']).read_text().strip()\n"
        "hclink.build(key, hclink.Database.ECOLI, db_dir=Path(os.environ['CHRONOCLADE_HCLINK_DB']), exact=True)\n"
    )
    python = Path(tool["environment"]) / "bin" / "python"
    environment = dict(os.environ)
    # The key is read inside the isolated interpreter; it is never placed in argv,
    # output, or the setup manifest.
    environment["CHRONOCLADE_ENTEROBASE_API_KEY_FILE"] = str(key_file.expanduser().resolve())
    environment["CHRONOCLADE_HCLINK_DB"] = str(output)
    _run([str(python), str(script)], cwd=Path(tool["source"]), env=environment, timeout=24 * 3600)
    try:
        metadata = json.loads((output / "metadata.json").read_text())
        index = output / "index.usearch"
        db_date = str(metadata["datestamp"]).split("T", 1)[0].split(" ", 1)[0]
        loci = output / "loci.json"
        if (
            not index.is_file()
            or not db_date
            or not loci.is_file()
            or not loci.read_text().strip()
            or not (output / "alleles.db").is_file()
            or not (output / "ST.txt.xz").is_file()
        ):
            raise ValueError
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise TypingSetupError("hclink completed without its indexed database metadata") from exc
    return {
        "status": "built",
        "scheme": SPECIES["ecoli"]["lineage_scheme"],
        "scheme_version": db_date,
        "database_path": str(output),
        "library_datestamp": metadata["datestamp"],
    }


def setup_typing(
    destination: str | Path,
    species: str | Sequence[str] = "klebsiella",
    *,
    prepare_db: bool = False,
    pasteur_secrets: str | Path | None = None,
    enterobase_key_file: str | Path | None = None,
    uv_executable: str = "uv",
) -> tuple[Path, dict[str, Any]]:
    """Install pinned callers/assigners and write a setup manifest.

    Database downloads/builds are intentionally opt-in. `prepare_db` downloads
    cgMLST data publicly by default; Klebsiella can use an existing Pasteur
    access/session token when protected client credentials are supplied. Full
    lineage databases require credentials and are never reported ready when absent.
    """
    selected = _normalise_species(species)
    uv = shutil.which(uv_executable)
    if uv is None:
        raise TypingSetupError(f"uv executable is unavailable: {uv_executable}")
    root = Path(destination).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    tools: dict[str, Any] = {"mlst": _install_mlst(root)}
    if "klebsiella" in selected:
        tools["plincer"] = _install_python_tool(uv, root, "plincer", "3.11")
    if "ecoli" in selected:
        tools["hclink"] = _install_python_tool(uv, root, "hclink", "3.12")
    db_builds: dict[str, Any] = {}
    if prepare_db:
        # This is opt-in because the public scheme download and index can be large.
        # Assignment reference DBs remain separate because both require upstream auth.
        for name in selected:
            if tools["mlst"].get("status") != "installed":
                db_builds[name] = {
                    "status": "blocked",
                    "reason": tools["mlst"].get("blocker", "Pathogenwatch mlst is not installed."),
                }
                continue
            db_builds[name] = _prepare_public_cgmlst(
                uv,
                root,
                name,
                tools,
                pasteur_secrets=(
                    Path(pasteur_secrets).expanduser()
                    if name == "klebsiella" and pasteur_secrets is not None
                    else None
                ),
            )
        if "klebsiella" in selected and pasteur_secrets is not None:
            db_builds["klebsiella_lineage"] = _build_plincer_db(
                root, tools["plincer"], Path(pasteur_secrets).expanduser()
            )
        if "ecoli" in selected and enterobase_key_file is not None:
            db_builds["ecoli_lineage"] = _build_hclink_db(
                root, tools["hclink"], Path(enterobase_key_file).expanduser()
            )
        access_modes = {
            name: db_builds[name].get("access_mode")
            for name in selected
            if isinstance(db_builds.get(name), dict)
        }
        mode_summary = ", ".join(
            f"{name}={mode}" for name, mode in access_modes.items() if mode
        ) or "no cgMLST download completed"
        db_note = (
            f"cgMLST download access: {mode_summary}. LIN/HierCC reference databases "
            "were built only when their credential files were supplied."
        )
    else:
        db_note = "No typing databases were downloaded or built."
    databases = {name: _database_state(root, name) for name in selected}
    query_config_path = root / "query_config.json"
    query_config = _query_config(root, selected, tools)
    query_config_path.write_text(json.dumps(query_config, indent=2) + "\n")
    report = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "destination": str(root),
        "species": selected,
        "tools": tools,
        "databases": databases,
        "cgmlst_builds": db_builds,
        "run_ready": query_config["ready"]
        and all(item.get("status") == "installed" for item in tools.values()),
        "database_preparation": db_note,
        "query_config_path": str(query_config_path),
        "query_config_ready": query_config["ready"],
    }
    path = root / "typing_setup.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    return path, report
