"""Run Pathogenwatch's native cgMLST and lineage assigners on query assemblies.

Prepared databases are supplied explicitly; this module never builds databases,
handles API credentials, or substitutes a neighbour's unresolved lineage code.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import sqlite3
import subprocess
import uuid
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from chronoclade.cglin import DEFAULT_SCHEME, PREFIX_CONSENSUS_METHOD, normalise_assignment, parse_code
from chronoclade.metadata import SAFE_IDENTIFIER, Sample
from chronoclade.progress import progress, progress_task
from chronoclade.execution_provenance import execution_provenance
from chronoclade.artifacts import file_sha256, write_json
from chronoclade.prepare_recovery import (
    DirectoryLock, RecoveryError, canonical_digest, load_object, owned_directory,
)
from chronoclade.plincer_scheme import validate_scheme, SchemeValidationError


class QueryTypingError(ValueError):
    """A typing configuration or result cannot safely be used."""


def normalise_typing_species(value: Any) -> str:
    """Return a supported taxon key while retaining user labels in output rows."""
    if not isinstance(value, str):
        raise QueryTypingError("Query typing requires a supported species identity")
    key = " ".join(value.replace("_", " ").casefold().split())
    if key not in {"klebsiella pneumoniae", "escherichia coli", "staphylococcus aureus"}:
        raise QueryTypingError(f"Unsupported native typing species: {value}")
    return key


def _same_species(left: Any, right: Any) -> bool:
    try:
        return normalise_typing_species(left) == normalise_typing_species(right)
    except QueryTypingError:
        return False


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _database_hash(path: Path) -> str:
    """Bind stable database content, excluding only SQLite locking sidecars.

    SQLite readers may create, alter and remove SHM files concurrently. They are
    not allele definitions; nonempty WAL/journal files do contain pending data
    and must stop reference reuse instead of hiding uncheckpointed writes.
    """
    path = Path(path)

    def sqlite_file(resource):
        try:
            with resource.open('rb') as stream:
                return stream.read(16) == b'SQLite format 3\x00'
        except FileNotFoundError:
            return False

    def checkpointed(resource):
        if sqlite_file(resource):
            for suffix in ('-wal', '-journal'):
                sidecar = resource.with_name(resource.name + suffix)
                try:
                    size = sidecar.stat().st_size
                except FileNotFoundError:
                    continue
                if size:
                    raise QueryTypingError(
                        f'Typing database has uncheckpointed SQLite writes: {resource.name}')

    def stable_hash(resource):
        checkpointed(resource)
        try:
            before = resource.stat()
            value = _sha256(resource)
            after = resource.stat()
        except FileNotFoundError as error:
            raise QueryTypingError('Stable typing database resource disappeared while hashing') from error
        if (before.st_ino, before.st_size, before.st_mtime_ns) != (
                after.st_ino, after.st_size, after.st_mtime_ns):
            raise QueryTypingError('Typing database resource changed while hashing')
        checkpointed(resource)
        return value

    if path.is_file():
        return stable_hash(path)
    files = []
    for item in sorted(path.rglob('*')):
        suffix = next((s for s in ('-shm', '-wal', '-journal') if item.name.endswith(s)), None)
        if suffix:
            base = item.with_name(item.name[:-len(suffix)])
            if sqlite_file(base):
                if suffix != '-shm':
                    try:
                        size = item.stat().st_size
                    except FileNotFoundError:
                        continue  # Reader-owned SQLite sidecar vanished after rglob.
                    if size:
                        raise QueryTypingError(
                            f'Typing database has uncheckpointed SQLite writes: {base.name}')
                continue
        if item.is_file():
            files.append(item)
        elif not item.exists():
            raise QueryTypingError('Stable typing database resource disappeared while hashing')
    if not files:
        raise QueryTypingError(f"Empty typing database: {path}")
    digest = hashlib.sha256()
    for item in files:
        digest.update(item.relative_to(path).as_posix().encode())
        digest.update(b"\0")
        digest.update(stable_hash(item).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def validate_cgmlst_index(index_dir: Path, scheme: str, scheme_version: str) -> list[str]:
    """Check the files and metadata consumed by the pinned cgps-mlst caller."""
    root = index_dir.expanduser().resolve()
    if not root.is_dir():
        raise QueryTypingError("cgMLST index_dir must be a directory")
    scheme_metadata_path = root / f"{scheme}_metadata.json"
    try:
        scheme_metadata = json.loads(scheme_metadata_path.read_text())
    except (OSError, ValueError) as exc:
        raise QueryTypingError(
            f"cgMLST index is missing valid {scheme} caller metadata"
        ) from exc
    if not isinstance(scheme_metadata, dict) or scheme_metadata.get("shortname") != scheme:
        raise QueryTypingError("cgMLST index metadata does not match its configured scheme")
    indexed_version = scheme_metadata.get("lastUpdated", scheme_metadata.get("last_updated"))
    if not isinstance(indexed_version, str) or indexed_version != scheme_version:
        raise QueryTypingError("cgMLST index version does not match its configured scheme version")
    relative_details = scheme_metadata.get("path")
    if not isinstance(relative_details, str) or not relative_details.strip():
        raise QueryTypingError("cgMLST index metadata has no indexed scheme details path")
    details_path = Path(relative_details)
    if details_path.is_absolute():
        raise QueryTypingError("cgMLST scheme details path must be relative to its index")
    details_path = (root / details_path).resolve()
    try:
        details_path.relative_to(root)
    except ValueError as exc:
        raise QueryTypingError("cgMLST scheme details path escapes its index directory") from exc
    try:
        with gzip.open(details_path, "rt", encoding="utf-8") as handle:
            details = json.load(handle)
    except (OSError, EOFError, UnicodeDecodeError, ValueError) as exc:
        raise QueryTypingError("cgMLST index has no readable compressed scheme details") from exc
    if not isinstance(details, dict) or details.get("shortname") != scheme:
        raise QueryTypingError("cgMLST indexed scheme details do not match the configured scheme")
    genes = details.get("genes")
    if (
        not isinstance(genes, list)
        or not genes
        or any(not isinstance(gene, str) or not gene for gene in genes)
        or len(set(genes)) != len(genes)
        or details.get("schemeSize") != len(genes)
    ):
        raise QueryTypingError("cgMLST index has an invalid ordered locus universe")
    allele_paths = details.get("allelePaths")
    if not isinstance(allele_paths, Mapping) or set(allele_paths) != set(genes):
        raise QueryTypingError("cgMLST index allele FASTA paths do not cover its locus universe")
    for locus, value in allele_paths.items():
        if not isinstance(value, str) or not value.strip() or Path(value).is_absolute():
            raise QueryTypingError(f"cgMLST index has an invalid allele FASTA path for {locus}")
        allele_path = (root / value).resolve()
        try:
            allele_path.relative_to(root)
        except ValueError as exc:
            raise QueryTypingError("cgMLST allele FASTA path escapes its index directory") from exc
        if not allele_path.is_file() or allele_path.stat().st_size == 0:
            raise QueryTypingError(f"cgMLST index is missing allele FASTA data for {locus}")
    for field in ("lengths", "alleleDictionary", "profiles"):
        if not isinstance(details.get(field), Mapping):
            raise QueryTypingError(f"cgMLST indexed scheme details lack {field}")
    if not set(genes).issubset(details["lengths"]):
        raise QueryTypingError("cgMLST indexed scheme lengths do not cover every locus")

    database_path = details_path.parent / "allele.db"
    if not database_path.is_file() or database_path.stat().st_size == 0:
        raise QueryTypingError("cgMLST index is missing its allele lookup database")
    try:
        connection = sqlite3.connect(database_path.as_uri() + "?mode=ro", uri=True)
        try:
            connection.execute("SELECT hash, gene, st, reverse FROM alleles LIMIT 0")
            if connection.execute("SELECT COUNT(*) FROM alleles").fetchone()[0] < 1:
                raise QueryTypingError("cgMLST allele lookup database contains no indexed alleles")
        finally:
            connection.close()
    except (OSError, sqlite3.Error) as exc:
        raise QueryTypingError("cgMLST allele lookup database is unreadable or malformed") from exc
    return genes


def _path(value: Any, base: Path) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise QueryTypingError("Typing database paths must be nonempty strings")
    path = Path(value).expanduser()
    path = (base / path).resolve() if not path.is_absolute() else path.resolve()
    if not path.exists():
        raise QueryTypingError(f"Typing resource does not exist: {path}")
    return path


def _command(value: Any, base: Path) -> list[str]:
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(x, str) or not x for x in value)
    ):
        raise QueryTypingError("Tool command must be a nonempty JSON array of arguments")
    command = list(value)
    executable = shutil.which(command[0])
    if executable is None:
        candidate = base / command[0]
        if candidate.is_file():
            executable = str(candidate.resolve())
    if executable is None:
        raise QueryTypingError(f"Typing executable is unavailable: {command[0]}")
    command[0] = executable
    # Script paths can be relative to the configuration, not the process cwd.
    for index, arg in enumerate(command[1:], 1):
        if (base / arg).is_file():
            command[index] = str((base / arg).resolve())
    return command


def load_typing_config(path: Path, *, profiles_only: bool = False) -> dict[str, Any]:
    """Validate native tool and frozen database paths, resolving relative paths."""
    path = path.expanduser().resolve()
    try:
        config = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise QueryTypingError("Cannot read typing configuration JSON") from exc
    if not profiles_only and isinstance(config, dict) and config.get("ready") is False:
        raise QueryTypingError(
            "Typing databases are not ready. Run setup-typing and supply the prepared database files described in its setup instructions."
        )
    organisms = config.get("organisms") if isinstance(config, dict) else None
    if not isinstance(organisms, dict) or not organisms:
        raise QueryTypingError("Typing config requires a nonempty organisms mapping")
    species_keys = set()
    for species, organism in organisms.items():
        if not isinstance(species, str) or not isinstance(organism, dict):
            raise QueryTypingError("Invalid organism typing configuration")
        species_key = normalise_typing_species(species)
        if species_key in species_keys:
            raise QueryTypingError("Typing config contains duplicate aliases of the same species")
        species_keys.add(species_key)
        for name in (("cgmlst",) if profiles_only else ("cgmlst", "assignment")):
            tool = organism.get(name)
            if not isinstance(tool, dict):
                raise QueryTypingError(f"{species}: missing {name} configuration")
            for field in ("scheme", "scheme_version", "tool_version"):
                if (
                    not isinstance(tool.get(field), str)
                    or not tool[field].strip()
                    or tool[field] == "unknown"
                ):
                    raise QueryTypingError(f"{species}: {name} requires pinned {field}")
            tool["command"] = _command(tool.get("command"), path.parent)
        caller = organism["cgmlst"]
        caller["index_dir"] = _path(caller.get("index_dir"), path.parent)
        caller["index_loci"] = validate_cgmlst_index(
            caller["index_dir"], caller["scheme"], caller["scheme_version"]
        )
        if profiles_only:
            continue
        if species_key == 'staphylococcus aureus':
            raise QueryTypingError('Native Staphylococcus aureus typing supports profiles-only; lineage assignment is not configured')
        assignment = organism["assignment"]
        assignment["loci_file"] = _path(assignment.get("loci_file"), path.parent)
        assignment["canonical_loci"] = _canonical_loci(assignment["loci_file"])
        if [gene.replace(" ", "_") for gene in caller["index_loci"]] != assignment[
            "canonical_loci"
        ]:
            raise QueryTypingError(
                "cgMLST caller index loci do not match the lineage reference locus order"
            )
        kind = assignment.get("kind")
        expected_kind = "plincer" if species_key == "klebsiella pneumoniae" else "hclink"
        if kind != expected_kind:
            raise QueryTypingError("Native typing assigner does not match its configured species")

        if kind == "plincer":
            if caller["scheme"] != "klebsiella_1" or assignment["scheme"] != DEFAULT_SCHEME:
                raise QueryTypingError("plincer requires klebsiella_1 and scgMLST629_S schemes")
            for field in ("scheme_file", "profiles_file", "alleles_db"):
                assignment[field] = _path(assignment.get(field), path.parent)
                if not assignment[field].is_file() or assignment[field].stat().st_size == 0:
                    raise QueryTypingError(f"plincer {field} must be a nonempty file")
            try:
                validate_scheme(assignment["scheme_file"])
            except (OSError, SchemeValidationError) as exc:
                raise QueryTypingError("plincer assignment scheme validation failed") from exc
            try:
                metadata = json.loads(assignment["loci_file"].read_text())
            except (OSError, ValueError) as exc:
                raise QueryTypingError("plincer loci metadata is unreadable") from exc
            if not isinstance(metadata, dict) or metadata.get("last_updated") != assignment[
                "scheme_version"
            ]:
                raise QueryTypingError("plincer database version does not match its metadata")
        elif kind == "hclink":
            if (
                caller["scheme"] != "ecoli_1"
                or assignment["scheme"] != "Escherichia.cgMLSTv1"
            ):
                raise QueryTypingError("E. coli hclink requires the ecoli_1 cgMLST scheme")
            assignment["reference_db"] = _path(assignment.get("reference_db"), path.parent)
            for filename in ("metadata.json", "index.usearch", "loci.json", "alleles.db", "ST.txt.xz"):
                resource = assignment["reference_db"] / filename
                if not resource.is_file() or resource.stat().st_size == 0:
                    raise QueryTypingError(f"hclink reference database is missing nonempty {filename}")
            if assignment["loci_file"].resolve() != (
                assignment["reference_db"] / "loci.json"
            ).resolve():
                raise QueryTypingError("hclink loci file must belong to its configured reference database")
            try:
                caller_date = date.fromisoformat(caller["scheme_version"])
                assignment_date = date.fromisoformat(assignment["scheme_version"])
                metadata = json.loads((assignment["reference_db"] / "metadata.json").read_text())
                actual_date = datetime.fromisoformat(metadata["datestamp"]).date()
            except (ValueError, KeyError, OSError, TypeError) as exc:
                raise QueryTypingError(
                    "hclink requires ISO database dates and metadata.json datestamp"
                ) from exc
            if caller_date > assignment_date or actual_date != assignment_date:
                raise QueryTypingError(
                    "hclink database must match its declared date and be at least as recent as cgMLST"
                )
        else:
            raise QueryTypingError("Assignment kind must be plincer or hclink")
    return config


def _canonical_loci(path: Path) -> list[str]:
    try:
        raw = json.loads(path.read_text())
        genes = raw.get("genes") if isinstance(raw, dict) else raw
    except (OSError, ValueError) as exc:
        raise QueryTypingError("Cannot read assignment database ordered loci JSON") from exc
    if (
        not isinstance(genes, list)
        or not genes
        or any(not isinstance(gene, str) or not gene for gene in genes)
    ):
        raise QueryTypingError("Assignment database requires a nonempty ordered loci list")
    genes = [gene.replace(" ", "_") for gene in genes]
    if len(set(genes)) != len(genes):
        raise QueryTypingError("Assignment database locus IDs must be unique")
    return genes


def _check_locus_order(genes: list[str], canonical: list[str]) -> None:
    if [gene.replace(" ", "_") for gene in genes] != canonical:
        raise QueryTypingError(
            "cgMLST locus order does not match the assignment database; rebuild compatible databases before typing"
        )


def _provenance(tool: Mapping[str, Any], database_fields: Sequence[str]) -> dict[str, Any]:
    command = tool["command"]
    correction = None
    if "scheme_file" in database_fields:
        try:
            correction = validate_scheme(tool["scheme_file"])
        except (OSError, SchemeValidationError) as exc:
            raise QueryTypingError("plincer assignment scheme validation failed") from exc
    return {
        **({"scheme_transformation": correction} if correction is not None else {}),
        "version": tool["tool_version"],
        "command": command,
        "tool_files_sha256": {arg: _sha256(Path(arg)) for arg in command if Path(arg).is_file()},
        "databases": {
            field: {"path": str(tool[field]), "sha256": _database_hash(tool[field])}
            for field in database_fields
        },
        "scheme": tool["scheme"],
        "scheme_version": tool["scheme_version"],
    }


def _run(command: list[str], data: bytes, destination: Path) -> dict[str, Any]:
    try:
        with progress_task("typing", Path(command[0]).name):
            result = subprocess.run(command, input=data, capture_output=True, timeout=3600, check=False)
            if result.returncode:
                # Never expose subprocess stderr or arguments, which may include credentials.
                raise QueryTypingError(
                    f"Native typing tool {Path(command[0]).name} failed with exit code {result.returncode}"
                )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise QueryTypingError(
            f"Native typing tool could not complete: {Path(command[0]).name}"
        ) from exc
    try:
        parsed = json.loads(result.stdout)
    except (UnicodeDecodeError, ValueError) as exc:
        raise QueryTypingError("Native typing tool did not return a JSON object") from exc
    if not isinstance(parsed, dict):
        raise QueryTypingError("Native typing tool did not return a JSON object")
    destination.write_text(json.dumps(parsed, indent=2) + "\n")
    return parsed


def parse_cgmlst_result(raw: Mapping[str, Any], scheme: str, version: str) -> dict[str, Any]:
    """Keep known alleles, missing calls, novel hashes and duplicate hits distinct."""
    genes = raw.get("genes")
    code = raw.get("code")
    if raw.get("scheme") != scheme or not isinstance(genes, list) or not genes:
        raise QueryTypingError("cgMLST output has an incompatible scheme or no ordered loci")
    if any(not isinstance(gene, str) or not gene for gene in genes) or len(set(genes)) != len(
        genes
    ):
        raise QueryTypingError("cgMLST output locus IDs must be unique nonempty strings")
    if (
        not isinstance(code, str)
        or len(code.split("_")) != len(genes)
        or raw.get("schemeSize", len(genes)) != len(genes)
    ):
        raise QueryTypingError("cgMLST code length does not match its ordered loci")
    profile, novel, missing, ambiguous = {}, {}, [], []
    hits = raw.get("alleles", {})
    for gene, allele in zip(genes, code.split("_"), strict=True):
        gene_hits = hits.get(gene, [])
        if isinstance(gene_hits, list) and len(gene_hits) > 1:
            ambiguous.append(gene)
            continue
        if not allele:
            missing.append(gene)
        elif allele.isdigit() and int(allele) > 0:
            profile[gene] = str(int(allele))
        elif re.fullmatch(r"[0-9a-fA-F]{40}", allele):
            novel[gene] = allele.lower()
        else:
            raise QueryTypingError("cgMLST output contains an invalid allele token")
    return {
        "cgmlst_scheme": scheme,
        "cgmlst_scheme_version": version,
        "cgmlst_loci": genes,
        "cgmlst_profile": profile,
        "cgmlst_novel_alleles": novel,
        "cgmlst_missing_loci": missing,
        "cgmlst_ambiguous_loci": ambiguous,
        "cgmlst_called_fraction": (len(profile) + len(novel)) / len(genes),
        "cgmlst_known_fraction": len(profile) / len(genes),
        "cgmlst_code": code,
    }


def _assignment_result(raw: dict[str, Any], config: Mapping[str, Any]) -> dict[str, Any]:
    if config["kind"] == "plincer":
        code = raw.get("LINcode", [])
        if not isinstance(code, list):
            raise QueryTypingError("plincer LINcode must be an array")
        # Multiple best profiles can disagree: retain only their common prefix.
        codes = [code] + [match.get("LINcode", []) for match in raw.get("matches", [])]
        conservative = []
        for index, part in enumerate(code):
            if str(part) == "*" or any(
                len(other) <= index or str(other[index]) != str(part) for other in codes
            ):
                break
            conservative.append(str(part))
        consensus = None
        known_prefix, prefix_status = parse_code(conservative)
        if known_prefix and prefix_status in {"complete", "partial"}:
            consensus = {
                "method": PREFIX_CONSENSUS_METHOD,
                "raw_assignment_sha256": hashlib.sha256(
                    json.dumps(raw, sort_keys=True, separators=(",", ":")).encode()
                ).hexdigest(),
                "known_prefix": list(known_prefix),
                "best_match_codes": [list(parts) for parts in sorted({tuple(map(str, parts)) for parts in codes})],
            }
        if len(conservative) < 10 and conservative:
            conservative += ["*"] * (10 - len(conservative))
        result = normalise_assignment(
            {
                "LINcode": "_".join(conservative),
                "cgST": raw.get("cgST", ""),
                "cglin_source": "pathogenwatch_plincer",
                **({"cglin_prefix_consensus": consensus} if consensus else {}),
            },
            scheme=config["scheme"],
            scheme_version=config["scheme_version"],
        )
        # A computed lineage cannot assert or clear the input's provider identity.
        result.pop("source_genome_id", None)
        return result
    codes = raw.get("hierCC")
    if not isinstance(codes, list):
        raise QueryTypingError("hclink output requires a hierCC array")
    result = {}
    for pair in codes:
        if (
            not isinstance(pair, list)
            or len(pair) != 2
            or not re.fullmatch(r"(?:d|HC)\d+", str(pair[0]))
            or (str(pair[1]) != "" and not str(pair[1]).isdigit())
        ):
            raise QueryTypingError("Invalid hclink HierCC assignment")
        if str(pair[1]) == "":
            continue
        level = "HC" + re.sub(r"^(?:d|HC)", "", str(pair[0]))
        if level in result:
            raise QueryTypingError("Duplicate hclink HierCC level")
        result[level] = str(pair[1])
    return {
        "hiercc_codes": result,
        "hiercc_scheme": config["scheme"],
        "hiercc_scheme_version": config["scheme_version"],
        "hiercc_source": "pathogenwatch_hclink",
        "hiercc_status": "predicted" if result else "unassigned",
        "hiercc_distance": raw.get("hierccDistance", raw.get("hierCC_distance")),
    }


def _typing_completion(index, total, sample_id, result):
    called = len(result.get("cgmlst_profile", {})) + len(result.get("cgmlst_novel_alleles", {}))
    progress(
        "typing",
        f"Completed typing sample {index}/{total}: {sample_id}; "
        f"{called}/{len(result.get('cgmlst_loci', []))} loci called; "
        f"{len(result.get('cgmlst_missing_loci', []))} missing; "
        f"{len(result.get('cgmlst_ambiguous_loci', []))} ambiguous",
    )


def _type_query_assemblies_once(
    samples: Sequence[Sample],
    config_path: Path,
    output: Path,
    *,
    existing_assignments: Mapping[str, dict[str, Any]] | None = None,
    profiles_only: bool = False,
    progress_index: int = 1,
    progress_total: int | None = None,
    _config_cache: dict[str, Any] | None = None,
    _execution_cache: dict[str, Any] | None = None,
    _provenance_cache: dict[str, dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Type each assembly using pinned native tools; write query_typing.json."""
    config = _config_cache or (load_typing_config(config_path, profiles_only=True) if profiles_only
                               else load_typing_config(config_path))
    mode = "profiles-only" if profiles_only else "full"
    execution = _execution_cache or execution_provenance(
        "native-query-typing", "pinned-caller-parser-lineage-v1"
    )
    output.mkdir(parents=True, exist_ok=True)
    organisms = {
        normalise_typing_species(species): settings
        for species, settings in config["organisms"].items()
    }
    results, provenance = [], {}
    seen = set()
    total = progress_total or len(samples)
    for index, sample in enumerate(samples, start=progress_index):
        if not SAFE_IDENTIFIER.fullmatch(sample.sample_id) or sample.sample_id in seen:
            raise QueryTypingError("Query sample IDs must be unique safe identifiers")
        seen.add(sample.sample_id)
        assembly_hash = _sha256(sample.assembly)
        existing = (existing_assignments or {}).get(sample.sample_id)
        species_key = normalise_typing_species(sample.species)
        organism = organisms.get(species_key)
        if organism is None:
            raise QueryTypingError(f"No native typing configuration for {sample.species}")
        caller = organism["cgmlst"]
        assigner = organism.get("assignment") if not profiles_only else None
        if species_key in (_provenance_cache or {}):
            provenance[species_key] = _provenance_cache[species_key]
        elif species_key not in provenance:
            provenance[species_key] = {
                "cgmlst": _provenance(caller, ("index_dir",)), "execution": execution
            }
            if not profiles_only:
                fields = (
                    ("scheme_file", "profiles_file", "alleles_db", "loci_file")
                    if assigner["kind"] == "plincer"
                    else ("reference_db", "loci_file")
                )
                provenance[species_key]["assignment"] = _provenance(assigner, fields)
        if (
            existing
            and existing.get("assembly_sha256") == assembly_hash
            and _same_species(existing.get("species"), sample.species)
            and existing.get("typing_provenance") == provenance[species_key]
            and existing.get("typing_mode", "full") == mode
        ):
            results.append(
                dict(
                    existing, sample_id=sample.sample_id, species=sample.species, typing_reused=True
                )
            )
            continue
        progress("typing", f"Typing sample {index}/{total}: {sample.sample_id}")
        directory = output / hashlib.sha256(sample.sample_id.encode()).hexdigest()[:24]
        directory.mkdir(exist_ok=True)
        opener = gzip.open if sample.assembly.suffix == ".gz" else open
        with opener(sample.assembly, "rb") as handle:
            assembly = handle.read()
        raw = _run(
            caller["command"]
            + [f"--scheme={caller['scheme']}", f"--indexDir={caller['index_dir']}"],
            assembly,
            directory / "cgmlst.json",
        )
        parsed = parse_cgmlst_result(raw, caller["scheme"], caller["scheme_version"])
        if profiles_only:
            _check_locus_order(raw["genes"], caller["index_loci"])
            results.append({
                "sample_id": sample.sample_id, "species": sample.species,
                "assembly_sha256": assembly_hash, "typing_mode": mode,
                "typing_status": "profiles_only", "lineage_status": "not_requested",
                **({'cglin_status': 'not_requested'} if species_key == 'klebsiella pneumoniae'
                   else {'hiercc_status': 'not_requested'} if species_key == 'escherichia coli' else {}),
                **parsed, "cgmlst_locus_universe_complete": True,
                "typing_provenance": provenance[species_key],
                "cgmlst_database_sha256": provenance[species_key]["cgmlst"]["databases"]["index_dir"]["sha256"],
            })
            _typing_completion(index, total, sample.sample_id, results[-1])
            continue
        _check_locus_order(raw["genes"], assigner["canonical_loci"])
        # Upstream assigns numeric code even for duplicate loci. Blank ambiguous
        # loci before classification so they cannot drive lineage assignments.
        safe_raw = dict(raw)
        safe_raw["code"] = "_".join(
            "" if gene in parsed["cgmlst_ambiguous_loci"] else allele
            for gene, allele in zip(raw["genes"], raw["code"].split("_"), strict=True)
        )
        if parsed["cgmlst_ambiguous_loci"]:
            safe_raw["st"] = "unassigned"
        if assigner["kind"] == "plincer":
            argv = [
                "classify",
                "-",
                "--scheme-file",
                str(assigner["scheme_file"]),
                "--profiles-file",
                str(assigner["profiles_file"]),
                "--alleles-db",
                str(assigner["alleles_db"]),
            ]
        else:
            argv = ["assign", "-", "--exact", "--reference-db", str(assigner["reference_db"])]
        if not parsed["cgmlst_profile"] and not parsed["cgmlst_novel_alleles"]:
            assignment = {"typing_status": "unassigned"}
        else:
            assignment_raw = _run(
                assigner["command"] + argv,
                json.dumps(safe_raw).encode(),
                directory / "assignment.json",
            )
            assignment = _assignment_result(assignment_raw, assigner)
            assignment["typing_status"] = (
                "typed"
                if assignment.get("cglin_resolved_depth", 0) or assignment.get("hiercc_codes")
                else "unassigned"
            )
            assignment["typing_assignment_evidence"] = assignment_raw

        results.append(
            {
                "sample_id": sample.sample_id,
                "species": sample.species,
                "assembly_sha256": assembly_hash,
                "typing_status": "typed",
                "typing_mode": mode,
                **parsed,
                **assignment,
                "typing_provenance": provenance[species_key],
                "cgmlst_database_sha256": provenance[species_key]["cgmlst"]["databases"][
                    "index_dir"
                ]["sha256"],
                (
                    "cglin_database_sha256"
                    if assigner["kind"] == "plincer"
                    else "hiercc_database_sha256"
                ): hashlib.sha256(
                    json.dumps(
                        {
                            field: info["sha256"]
                            for field, info in provenance[species_key]["assignment"][
                                "databases"
                            ].items()
                        },
                        sort_keys=True,
                    ).encode()
                ).hexdigest(),
            }
        )
        _typing_completion(index, total, sample.sample_id, results[-1])
    envelope = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "typing_mode": mode,
        "assignments": results,
    }
    (output / "query_typing.json").write_text(json.dumps(envelope, indent=2) + "\n")
    return results


def load_query_typing(path: Path) -> list[dict[str, Any]]:
    """Load native typing output; assembly identity is checked by its consumer."""
    try:
        envelope = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise QueryTypingError("Cannot read query typing JSON") from exc
    return _validate_query_typing_envelope(envelope)


def _validate_query_typing_envelope(envelope: dict[str, Any]) -> list[dict[str, Any]]:
    """Validate decoded query typing data for both files and recovery envelopes."""
    if (
        not isinstance(envelope, dict)
        or envelope.get("schema_version") != 1
        or not isinstance(envelope.get("assignments"), list)
    ):
        raise QueryTypingError("Unsupported query typing schema")
    seen = set()
    for row in envelope["assignments"]:
        if (
            not isinstance(row, dict)
            or not isinstance(row.get("sample_id"), str)
            or not SAFE_IDENTIFIER.fullmatch(row["sample_id"])
            or row["sample_id"] in seen
        ):
            raise QueryTypingError("Query typing IDs must be unique safe identifiers")
        seen.add(row["sample_id"])
        if row.get("source_genome_id") not in (None, ""):
            raise QueryTypingError("Native query typing cannot supply a provider identity")
        if row.get("typing_mode") == "profiles-only" and (
            row.get("typing_status") != "profiles_only"
            or row.get("lineage_status") != "not_requested"
            or "hiercc_codes" in row or "cglin_raw" in row
        ):
            raise QueryTypingError("Profiles-only typing cannot claim a lineage assignment")
        if not re.fullmatch(r"[0-9a-f]{64}", str(row.get("assembly_sha256", ""))):
            raise QueryTypingError("Query typing requires an assembly SHA256")
        if "cgmlst_profile" in row and (
            not isinstance(row["cgmlst_profile"], dict)
            or not isinstance(row.get("cgmlst_loci"), list)
        ):
            raise QueryTypingError("Query cgMLST profile requires its declared locus universe")
        if not isinstance(row.get("species"), str) or not row["species"].strip():
            raise QueryTypingError("Query typing requires species identity")
        species_key = normalise_typing_species(row["species"])
        if "cgmlst_profile" in row:
            loci, profile = row["cgmlst_loci"], row["cgmlst_profile"]
            if (
                not loci
                or any(not isinstance(x, str) or not x for x in loci)
                or len(set(loci)) != len(loci)
            ):
                raise QueryTypingError("Invalid query typing locus universe")
            if not set(profile).issubset(loci) or any(
                not str(v).isdigit() or int(v) < 1 for v in profile.values()
            ):
                raise QueryTypingError("Invalid known query cgMLST alleles")
            for field in ("cgmlst_scheme", "cgmlst_scheme_version"):
                if not isinstance(row.get(field), str) or row[field] in {"", "unknown"}:
                    raise QueryTypingError("Query profile requires a pinned scheme")
        expected_scheme = {'klebsiella pneumoniae': 'klebsiella_1',
                           'escherichia coli': 'ecoli_1', 'staphylococcus aureus': 'saureus_1'}[species_key]
        if "cgmlst_scheme" in row and row["cgmlst_scheme"] != expected_scheme:
            raise QueryTypingError("Query cgMLST scheme does not match its species")
        if "cgmlst_novel_alleles" in row:
            novel = row["cgmlst_novel_alleles"]
            if (
                not isinstance(novel, dict)
                or not set(novel).issubset(row.get("cgmlst_loci", []))
                or any(not re.fullmatch(r"[0-9a-f]{40}", str(v)) for v in novel.values())
            ):
                raise QueryTypingError("Invalid novel query cgMLST allele checksums")
            if set(novel).intersection(row.get("cgmlst_profile", {})):
                raise QueryTypingError("Query locus cannot be both known and novel")
        if "cglin_raw" in row:
            if species_key != "klebsiella pneumoniae":
                raise QueryTypingError("cgLIN query assignment requires Klebsiella pneumoniae")
            if (
                row.get("cglin_scheme") != DEFAULT_SCHEME
                or not isinstance(row.get("cglin_scheme_version"), str)
                or row["cglin_scheme_version"] in {"", "unknown"}
            ):
                raise QueryTypingError("Query cgLIN requires a pinned compatible scheme")
            assignment = normalise_assignment(row)
            if assignment["cglin_code_status"] in {"malformed", "unsupported"}:
                raise QueryTypingError("Invalid query cgLIN assignment")
        if "hiercc_codes" in row:
            if species_key != "escherichia coli":
                raise QueryTypingError("HierCC query assignment requires Escherichia coli")
            codes = row["hiercc_codes"]
            if not isinstance(codes, dict) or any(
                not re.fullmatch(r"HC\d+", str(k)) or not str(v).isdigit() for k, v in codes.items()
            ):
                raise QueryTypingError("Invalid query HierCC assignment")
            if any(
                not isinstance(row.get(k), str) or row[k] in {"", "unknown"}
                for k in ("hiercc_scheme", "hiercc_scheme_version")
            ):
                raise QueryTypingError("Query HierCC requires a pinned scheme")
    # Legacy Plincer normalization emitted an empty identity. Do not let a
    # frozen typing import overwrite the independently resolved input identity.
    return [
        {key: value for key, value in row.items() if key != "source_genome_id"}
        for row in envelope["assignments"]
    ]


def _typing_recovery_identity(samples, config_path, profiles_only):
    """Bind a recovery cache to ordered inputs, tools, references, and runtime."""
    config = (load_typing_config(config_path, profiles_only=True) if profiles_only
              else load_typing_config(config_path))
    execution = execution_provenance("native-query-typing", "pinned-caller-parser-lineage-v1")
    organisms = {normalise_typing_species(key): value for key, value in config["organisms"].items()}
    identities = []
    provenance_by_species = {}
    for sample in samples:
        species_key = normalise_typing_species(sample.species)
        organism = organisms.get(species_key)
        if organism is None:
            raise QueryTypingError(f"No native typing configuration for {sample.species}")
        caller = organism["cgmlst"]
        if species_key not in provenance_by_species:
            provenance = {"cgmlst": _provenance(caller, ("index_dir",)), "execution": execution}
            if not profiles_only:
                assigner = organism["assignment"]
                fields = (("scheme_file", "profiles_file", "alleles_db", "loci_file")
                          if assigner["kind"] == "plincer" else ("reference_db", "loci_file"))
                provenance["assignment"] = _provenance(assigner, fields)
            provenance_by_species[species_key] = provenance
        provenance = provenance_by_species[species_key]
        identities.append({
            "sample_id": sample.sample_id,
            "species": species_key,
            "assembly_sha256": _sha256(sample.assembly),
            "typing_provenance": provenance,
        })
    return {
        "schema": "chronoclade.query-typing-recovery-identity",
        "schema_version": 1,
        "mode": "profiles-only" if profiles_only else "full",
        "config_sha256": file_sha256(Path(config_path).expanduser().resolve()),
        "execution": execution,
        "samples": identities,
    }


def type_query_assemblies(
    samples: Sequence[Sample],
    config_path: Path,
    output: Path,
    *,
    existing_assignments: Mapping[str, dict[str, Any]] | None = None,
    profiles_only: bool = False,
    resume: bool = False,
    attempt_output: Path | None = None,
) -> list[dict[str, Any]]:
    """Type samples with atomic per-sample checkpoints and explicit resume."""
    from chronoclade.metadata import SAFE_IDENTIFIER

    output = Path(output).expanduser()
    if output.is_symlink():
        raise QueryTypingError(f"Typing output path is a symlink and will not be adopted: {output}")
    output = output.resolve()
    ids = [sample.sample_id for sample in samples]
    if any(not SAFE_IDENTIFIER.fullmatch(ident) for ident in ids) or len(set(ids)) != len(ids):
        raise QueryTypingError("Query sample IDs must be unique safe identifiers")
    manifest_path = output / "typing-recovery.json"
    partial_path = output / "query_typing.partial.json"
    if not resume and any((output / name).exists() or (output / name).is_symlink() for name in (
        "typing-recovery.json", "query_typing.partial.json", "query_typing.json"
    )):
        raise QueryTypingError(
            f"Typing checkpoint already exists at {output}; use --resume or choose a fresh output"
        )
    identity = _typing_recovery_identity(samples, config_path, profiles_only)
    runtime_config = (load_typing_config(config_path, profiles_only=True) if profiles_only
                      else load_typing_config(config_path))
    provenance_by_species = {}
    for item in identity["samples"]:
        provenance_by_species.setdefault(item["species"], item["typing_provenance"])
    identity_sha = canonical_digest(identity)
    output.mkdir(parents=True, exist_ok=True)
    with DirectoryLock(output / ".typing-recovery.lock"):
        cached = {}
        if resume:
            if manifest_path.is_symlink() or partial_path.is_symlink():
                raise QueryTypingError(f"Typing checkpoint contains a symlink in {output}")
            manifest = load_object(manifest_path, "typing recovery manifest")
            if (manifest.get("schema") != "chronoclade.query-typing-recovery"
                    or manifest.get("schema_version") != 1
                    or manifest.get("identity") != identity
                    or manifest.get("identity_sha256") != identity_sha):
                raise QueryTypingError(
                    f"Typing checkpoint identity differs or is corrupt at {output}; use a fresh output directory"
                )
            if not partial_path.exists() and not manifest.get("checkpoint_initialized"):
                _write_typing_checkpoint(partial_path, identity_sha, [])
                manifest["checkpoint_initialized"] = True
                write_json(manifest_path, manifest)
            try:
                envelope = load_object(partial_path, "typing checkpoint")
                rows = _validate_typing_checkpoint_envelope(envelope, identity_sha)
            except Exception as error:
                raise QueryTypingError(f"Invalid typing checkpoint at {partial_path}: {error}") from error
            row_ids = [row["sample_id"] for row in rows]
            if row_ids != ids[:len(row_ids)]:
                raise QueryTypingError("Typing checkpoint rows are not the recorded completed cohort prefix")
            expected_by_id = {row["sample_id"]: row for row in identity["samples"]}
            for row in rows:
                expected = expected_by_id.get(row["sample_id"])
                if expected is None:
                    raise QueryTypingError("Typing checkpoint contains an unknown sample ID")
                if (row.get("assembly_sha256") != expected["assembly_sha256"]
                        or normalise_typing_species(row.get("species")) != expected["species"]
                        or row.get("typing_provenance") != expected["typing_provenance"]
                        or row.get("typing_mode", "full") != identity["mode"]):
                    raise QueryTypingError(
                        f"Typing checkpoint row identity/provenance differs for {row['sample_id']}"
                    )
                cached[row["sample_id"]] = row
        else:
            if manifest_path.exists() or partial_path.exists() or (output / "query_typing.json").exists():
                raise QueryTypingError(
                    f"Typing checkpoint already exists at {output}; use --resume or choose a fresh output"
                )
            manifest = {
                "schema": "chronoclade.query-typing-recovery",
                "schema_version": 1,
                "status": "incomplete",
                "identity": identity,
                "identity_sha256": identity_sha,
                "expected_sample_ids": ids,
                "completed_sample_ids": [],
                "checkpoint_initialized": False,
            }
            write_json(manifest_path, manifest)
            _write_typing_checkpoint(partial_path, identity_sha, [])
            manifest["checkpoint_initialized"] = True
            write_json(manifest_path, manifest)

        results = []
        for sample in samples:
            row = cached.get(sample.sample_id)
            newly_computed = row is None
            if newly_computed:
                if attempt_output is not None:
                    attempt_root = Path(attempt_output).expanduser()
                    if attempt_root.is_symlink():
                        raise QueryTypingError(f"Typing attempt path is a symlink: {attempt_root}")
                    attempt_root.mkdir(parents=True, exist_ok=True)
                    attempt = attempt_root / uuid.uuid4().hex
                    if attempt.exists() or attempt.is_symlink():
                        raise QueryTypingError(f"Typing attempt destination already exists: {attempt}")
                    attempt.mkdir()
                else:
                    try:
                        attempts_root = owned_directory(output, "attempts", create=True)
                        attempt = attempts_root / uuid.uuid4().hex
                        if attempt.exists() or attempt.is_symlink():
                            raise RecoveryError(f"Typing attempt destination already exists: {attempt}")
                        attempt.mkdir()
                    except RecoveryError as error:
                        raise QueryTypingError(str(error)) from error
                external = (existing_assignments or {}).get(sample.sample_id)
                computed = _type_query_assemblies_once(
                    [sample], config_path, attempt,
                    existing_assignments={sample.sample_id: external} if external else None,
                    profiles_only=profiles_only,
                    progress_index=len(results) + 1,
                    progress_total=len(samples),
                    _config_cache=runtime_config,
                    _execution_cache=identity["execution"],
                    _provenance_cache=provenance_by_species,
                )
                if len(computed) != 1:
                    raise QueryTypingError("Native typer returned an incomplete sample result")
                row = computed[0]
            else:
                row = dict(row, typing_reused=True)
            results.append(row)
            if newly_computed:
                _write_typing_checkpoint(partial_path, identity_sha, results)
                manifest["completed_sample_ids"] = [item["sample_id"] for item in results]
                write_json(manifest_path, manifest)
        write_json(output / "query_typing.json", {
            "schema_version": 1,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "typing_mode": identity["mode"],
            "assignments": results,
        })
        manifest["status"] = "complete"
        write_json(manifest_path, manifest)
        return results


def _typing_checkpoint_envelope(identity_sha, assignments):
    body = {"identity_sha256": identity_sha, "assignments": assignments}
    return {
        "schema": "chronoclade.query-typing-checkpoint",
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "identity_sha256": identity_sha,
        "assignments": assignments,
        "integrity_sha256": canonical_digest(body),
    }


def _write_typing_checkpoint(path, identity_sha, assignments):
    write_json(path, _typing_checkpoint_envelope(identity_sha, assignments))


def _validate_typing_checkpoint_envelope(envelope, identity_sha):
    if (envelope.get("schema") != "chronoclade.query-typing-checkpoint"
            or envelope.get("schema_version") != 1
            or envelope.get("identity_sha256") != identity_sha
            or not isinstance(envelope.get("assignments"), list)
            or envelope.get("integrity_sha256") != canonical_digest({
                "identity_sha256": identity_sha, "assignments": envelope.get("assignments")
            })):
        raise QueryTypingError("Checkpoint envelope identity or integrity check failed")
    return _validate_query_typing_envelope({
        "schema_version": 1, "assignments": envelope["assignments"]
    })


def typing_checkpoint_sample_ids(output: Path) -> list[str]:
    """Return checkpoint identities after validating its self-contained envelope."""
    output = Path(output)
    if output.is_symlink():
        raise QueryTypingError(f"Typing checkpoint directory is a symlink: {output}")
    partial_path = output / "query_typing.partial.json"
    if partial_path.is_symlink():
        raise QueryTypingError(f"Typing checkpoint file is a symlink: {partial_path}")
    if not partial_path.exists():
        manifest_path = output / "typing-recovery.json"
        if manifest_path.is_symlink():
            raise QueryTypingError(f"Typing recovery manifest is a symlink: {manifest_path}")
        manifest = load_object(manifest_path, "typing recovery manifest")
        identity = manifest.get("identity")
        identity_sha = manifest.get("identity_sha256")
        if (manifest.get("schema") != "chronoclade.query-typing-recovery"
                or manifest.get("schema_version") != 1
                or manifest.get("status") != "incomplete"
                or manifest.get("checkpoint_initialized") is not False
                or not isinstance(identity, dict)
                or identity_sha != canonical_digest(identity)):
            raise QueryTypingError("Missing typing checkpoint is not a valid pre-initialization state")
        return []
    envelope = load_object(partial_path, "typing checkpoint")
    identity_sha = envelope.get("identity_sha256")
    if not isinstance(identity_sha, str) or not re.fullmatch(r"[0-9a-f]{64}", identity_sha):
        raise QueryTypingError("Typing checkpoint has no valid identity digest")
    rows = _validate_typing_checkpoint_envelope(envelope, identity_sha)
    return [row.get("sample_id") for row in rows]


def validate_typing_checkpoint(samples, config_path, output, *, profiles_only=False):
    """Fail closed on cached rows before prepare adopts any profile assignments."""
    identity = _typing_recovery_identity(samples, config_path, profiles_only)
    output = Path(output)
    manifest_path = output / "typing-recovery.json"
    partial_path = output / "query_typing.partial.json"
    if output.is_symlink() or manifest_path.is_symlink() or partial_path.is_symlink():
        raise QueryTypingError(f"Typing checkpoint contains a symlink in {output}")
    manifest = load_object(manifest_path, "typing recovery manifest")
    saved_identity = manifest.get("identity")
    compatible = (
        isinstance(saved_identity, dict)
        and {key: value for key, value in saved_identity.items() if key != "samples"}
        == {key: value for key, value in identity.items() if key != "samples"}
        and saved_identity.get("samples", [])[:len(identity["samples"])] == identity["samples"]
    )
    if not compatible:
        raise QueryTypingError("Typing checkpoint identity differs from the resolved cohort")
    envelope = load_object(partial_path, "typing checkpoint")
    saved_identity_sha = manifest.get("identity_sha256")
    if envelope.get("identity_sha256") != saved_identity_sha:
        raise QueryTypingError("Typing checkpoint envelope differs from its manifest identity")
    rows = _validate_typing_checkpoint_envelope(envelope, saved_identity_sha)
    ids = [sample.sample_id for sample in samples]
    row_ids = [row.get("sample_id") for row in rows]
    if row_ids != ids[:len(row_ids)]:
        raise QueryTypingError("Typing checkpoint rows are not the recorded completed cohort prefix")
    for row, expected in zip(rows, identity["samples"], strict=False):
        if (row.get("sample_id") != expected["sample_id"]
                or row.get("assembly_sha256") != expected["assembly_sha256"]
                or normalise_typing_species(row.get("species")) != expected["species"]
                or row.get("typing_provenance") != expected["typing_provenance"]
                or row.get("typing_mode", "full") != identity["mode"]):
            raise QueryTypingError(f"Typing checkpoint row identity differs for {expected['sample_id']}")
