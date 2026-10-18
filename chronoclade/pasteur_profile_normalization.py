"""Canonical Pasteur KpSC profiles with original job/reference evidence retained.

Execution-job labels are provenance, not biological scheme identities. This
adapter is deliberately restricted to the documented Pasteur scheme 18 route.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re
import sqlite3

from chronoclade.datasets import DatasetError
from chronoclade.artifacts import file_sha256
from chronoclade.pathogenwatch import content_hash
from chronoclade.query_typing import _assignment_result, _provenance, _run, load_typing_config

CANONICAL_SCHEME = "pasteur:pubmlst_klebsiella_seqdef:18"
COMPARISON_VERSION = "locus-oriented-dna-sha1-v1"
SUPPORTED_SCOPES = {"klebsiella_1", "pathogenwatch:cgmlst-573-2", CANONICAL_SCHEME}


REFERENCE_SCHEMA = "chronoclade.native-typing-reference"
REFERENCE_VERSION = 2


def _assert_checkpointed_file(path):
    path = Path(path).resolve()
    for suffix in ("-wal", "-journal"):
        sidecar = path.with_name(path.name + suffix)
        if sidecar.is_file() and sidecar.stat().st_size:
            raise DatasetError(f"Native reference has uncheckpointed SQLite writes: {path.name}")
    if not path.is_file():
        raise DatasetError(f"Native reference resource is missing: {path.name}")
    return path


def _checkpointed_file(path):
    path = _assert_checkpointed_file(path)
    return {"path": str(path), "sha256": file_sha256(path), "bytes": path.stat().st_size}


def _directory_reference(path):
    """Use shared file/content hashes; omit only recognized SQLite lock sidecars."""
    path = Path(path).resolve()
    files = {}
    for item in sorted(path.rglob("*")):
        if not item.is_file():
            continue
        if item.name.endswith(("-shm", "-wal", "-journal")):
            base = item.with_name(item.name.rsplit("-", 1)[0])
            if base.is_file():
                with base.open("rb") as stream:
                    sqlite = stream.read(16) == b"SQLite format 3\x00"
                if sqlite:
                    if not item.name.endswith("-shm") and item.stat().st_size:
                        raise DatasetError(
                            f"Native reference has uncheckpointed SQLite writes: {base.name}"
                        )
                    continue
        files[item.relative_to(path).as_posix()] = _checkpointed_file(item)
    if not files:
        raise DatasetError("Native reference directory contains no stable resources")
    return {"path": str(path), "files": files, "sha256": content_hash(files)}


def reference_descriptor(path):
    """Bind canonicalization/assignment inputs by content, without credentials."""
    path = Path(path).expanduser().resolve()
    before = file_sha256(path)
    try:
        config = load_typing_config(path)
    except ValueError as error:
        raise DatasetError(
            "Native reference binding requires full caller and lineage readiness; import profiles-only query typing without --typing-config"
        ) from error
    organisms = {}
    for species, tools in sorted(config["organisms"].items()):
        caller, assigner = tools["cgmlst"], tools["assignment"]
        index = caller["index_dir"]
        metadata_path = index / f"{caller['scheme']}_metadata.json"
        metadata = json.loads(metadata_path.read_text())
        details_path = (index / metadata["path"]).resolve()
        if not details_path.is_relative_to(index.resolve()):
            raise DatasetError("Native caller metadata escapes its reference index")
        # These are the caller resources consumed by canonical profile translation;
        # native assembly-calling caches retain their separate full-index provenance.
        caller_files = {
            name: _checkpointed_file(value)
            for name, value in {
                "allele.db": index / "allele.db",
                "scheme_metadata": metadata_path,
                "scheme_details": details_path,
            }.items()
        }
        fields = (
            ("scheme_file", "profiles_file", "alleles_db", "loci_file")
            if assigner["kind"] == "plincer"
            else ("reference_db", "loci_file")
        )
        file_fields = [field for field in fields if Path(assigner[field]).is_file()]
        for field in file_fields:
            _assert_checkpointed_file(assigner[field])
        assignment = _provenance(assigner, file_fields)
        for field in fields:
            if field not in file_fields:
                assignment["databases"][field] = _directory_reference(assigner[field])
        organisms[species] = {
            "caller": {
                **_provenance(caller, ()),
                "resources": caller_files,
                "indexed_loci_sha256": content_hash(caller["index_loci"]),
                "canonical_loci_sha256": content_hash(assigner["canonical_loci"]),
            },
            "assignment": assignment,
        }
    if file_sha256(path) != before:
        raise DatasetError("Native typing configuration changed while binding its resources")
    descriptor = {
        "schema": REFERENCE_SCHEMA,
        "schema_version": REFERENCE_VERSION,
        "config_path": str(path),
        "config_sha256": before,
        "species": sorted(config["organisms"]),
        "resources": organisms,
    }
    descriptor["sha256"] = content_hash(descriptor)
    return descriptor


def reference_rebinding(previous, current, *, explicit):
    """Expose an authorized fresh reference binding; never relabel old evidence."""
    if not explicit or previous is None or previous == current:
        return None
    return {
        "kind": (
            "legacy_config_only_to_content_bound"
            if not isinstance(previous, dict) or previous.get("schema") != REFERENCE_SCHEMA
            else "explicit_content_reference_rebind"
        ),
        "action_basis": "explicit --typing-config for fresh analysis",
        "previous_descriptor_sha256": content_hash(previous),
        "current_descriptor_sha256": current["sha256"],
        "historical_descriptor_unchanged": True,
    }


def resolve_reference(dataset, explicit=None):
    if explicit is not None:
        return Path(explicit).expanduser().resolve()
    descriptor = dataset.parameters.get("native_typing_reference")
    if not descriptor:
        return None
    if (
        not isinstance(descriptor, dict)
        or descriptor.get("schema") != REFERENCE_SCHEMA
        or type(descriptor.get("schema_version")) is not int
        or descriptor.get("schema_version") != REFERENCE_VERSION
        or descriptor.get("sha256")
        != content_hash({k: v for k, v in descriptor.items() if k != "sha256"})
    ):
        raise DatasetError(
            "Prepared native reference has a legacy/unverified content binding; supply --typing-config for a fresh explicit binding"
        )
    locator = descriptor.get("config_path")
    if not isinstance(locator, str) or not locator:
        raise DatasetError("Prepared native reference has no valid configuration locator")
    path = Path(locator)
    try:
        changed = not path.is_file() or reference_descriptor(path) != descriptor
    except (OSError, ValueError) as error:
        raise DatasetError(
            "Prepared native reference is unavailable or changed; supply --typing-config for a fresh analysis"
        ) from error
    if changed:
        raise DatasetError(
            "Prepared native reference content changed; supply --typing-config for a fresh analysis, not an old numerical resume"
        )
    return path


class PasteurProfileNormalizer:
    """Verify a canonical locus panel, normalise alleles, then genuinely reassign LINs."""

    def __init__(self, config_path, output):
        self.config_path = Path(config_path).expanduser().resolve()
        config = load_typing_config(self.config_path)
        organism = next((v for k, v in config["organisms"].items()
                         if k.replace("_", " ").casefold() == "klebsiella pneumoniae"), None)
        if organism is None:
            raise DatasetError("Public Klebsiella normalisation requires its installed native typing references")
        self.caller, self.assigner = organism["cgmlst"], organism["assignment"]
        self.output = Path(output)
        self.output.mkdir(parents=True, exist_ok=True)
        metadata = json.loads((self.caller["index_dir"] / "klebsiella_1_metadata.json").read_text())
        if (metadata.get("host"), metadata.get("hostPath"), str(metadata.get("schemeId"))) != (
            "pasteur", "pubmlst_klebsiella_seqdef", "18"
        ):
            raise DatasetError("Native reference is not verified Pasteur Klebsiella scheme 18")
        self.loci = self.assigner["canonical_loci"]
        self.locus_set = set(self.loci)
        if len(self.loci) != 629 or len(set(self.loci)) != 629:
            raise DatasetError("Canonical Pasteur KpSC comparison requires the complete 629-locus panel")
        self.panel_digest = content_hash({"scheme": CANONICAL_SCHEME,
                                         "policy": COMPARISON_VERSION, "loci": self.loci})
        self.assignment_provenance = _provenance(self.assigner, (
            "scheme_file", "profiles_file", "alleles_db", "loci_file"
        ))
        self.assignment_digest = content_hash({key: value["sha256"] for key, value in
                                              self.assignment_provenance["databases"].items()})
        self.sequences = {}
        self.source_files = {}
        self.audits = []

    def _load_locus(self, locus):
        # The caller's FASTAs contain a BLAST representative subset, not every
        # allele. Its SQLite exact-hit index retains full oriented DNA SHA1s.
        if not self.sequences:
            path = self.caller["index_dir"] / "allele.db"
            if not path.is_file() or (path.with_name(path.name + "-wal").exists()
                    and path.with_name(path.name + "-wal").stat().st_size):
                raise DatasetError("Canonical allele index is missing or has uncheckpointed writes")
            self.source_files["allele.db"] = hashlib.sha256(path.read_bytes()).hexdigest()
            try:
                with sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True) as connection:
                    for digest, gene, number in connection.execute(
                        "SELECT hash, gene, st FROM alleles WHERE reverse=0"
                    ):
                        if gene not in self.locus_set or not re.fullmatch(r"[0-9a-f]{40}", digest):
                            raise DatasetError("Canonical allele index contains an invalid locus or DNA digest")
                        values = self.sequences.setdefault(gene, {})
                        ident = str(number)
                        if not ident.isdigit() or int(ident) < 1:
                            raise DatasetError("Canonical allele index contains an invalid allele ID")
                        if ident in values and values[ident] != digest:
                            raise DatasetError(f"Conflicting allele sequence definitions at {gene}, allele {ident}")
                        values[ident] = digest
            except sqlite3.Error as error:
                raise DatasetError("Canonical allele hash index is unreadable") from error
            if set(self.sequences) != set(self.loci):
                raise DatasetError("Canonical allele hash index does not cover all 629 loci")
        return self.sequences[locus]

    def normalize(self, rows):
        """Return a new comparison representation; never mutate original records."""
        adapted, requests, evidence = [], [], []
        for original in rows:
            row = deepcopy(original)
            if row.get("species", "").replace("_", " ").casefold() != "klebsiella pneumoniae":
                raise DatasetError("Pasteur profile normalisation currently supports Klebsiella pneumoniae")
            if row.get("cgmlst_scheme") not in SUPPORTED_SCOPES:
                raise DatasetError(f"Unverified cgMLST source namespace for {row['sample_id']}")
            regular, novel = row.get("cgmlst_profile") or {}, row.get("cgmlst_novel_alleles") or {}
            if set(regular) & set(novel):
                raise DatasetError("Called and novel alleles overlap during normalisation")
            calls = dict(regular, **novel)
            ambiguous = set(row.get("cgmlst_ambiguous_loci") or [])
            if (set(calls) | ambiguous | set(row.get("cgmlst_loci") or [])) - set(self.loci):
                raise DatasetError("Source calls lie outside canonical Pasteur scheme 18")
            if row.get("cgmlst_locus_universe_complete") is True and set(row.get("cgmlst_loci") or []) != set(self.loci):
                raise DatasetError("Declared complete source locus panel differs from canonical Pasteur scheme 18")
            normalized, unsupported = {}, []
            safe_calls = {}
            for locus, token in sorted(calls.items()):
                if locus in ambiguous:
                    continue
                if token is None or token == 0:
                    continue
                token = str(token)
                if token.isdigit() and int(token) > 0:
                    token = str(int(token))
                    digest = self._load_locus(locus).get(token)
                    if digest is None:
                        unsupported.append({"locus": locus, "allele": token,
                                            "reason": "numbered_allele_unavailable_in_verified_reference"})
                        continue
                elif re.fullmatch(r"[0-9a-fA-F]{40}", token):
                    digest = token.lower()
                elif token in {"", "0", "-", "?"}:
                    continue
                else:
                    raise DatasetError(f"Invalid allele token for {row['sample_id']} at {locus}")
                normalized[locus], safe_calls[locus] = digest, token
            raw = {"scheme": "klebsiella_1", "schemeSize": 629, "genes": self.loci,
                   "code": "_".join(safe_calls.get(locus, "") for locus in self.loci),
                   "st": "unassigned"}
            requests.append(raw)
            evidence.append({"sample_id": row["sample_id"],
                             "source_genome_id": row.get("source_genome_id"),
                             "original_typing": {k: v for k, v in original.items()
                                                 if k.startswith(("cgmlst_", "cglin_", "typing_"))},
                             "unsupported_calls": unsupported,
                             "normalized_called_loci": len(normalized)})
            row.update(cgmlst_scheme=CANONICAL_SCHEME, cgmlst_scheme_version=COMPARISON_VERSION,
                       cgmlst_database_version=None, cgmlst_database_sha256=self.panel_digest,
                       cgmlst_profile=normalized, cgmlst_novel_alleles={},
                       cgmlst_loci=self.loci, cgmlst_locus_universe_complete=True,
                       cgmlst_missing_loci=sorted(set(self.loci) - set(normalized) - ambiguous),
                       cgmlst_status="resolved" if normalized else "unassigned")
            adapted.append(row)
        if not adapted:
            return []
        batch = content_hash(requests)
        request_path, result_path = self.output / (batch + "-requests.json"), self.output / (batch + "-assignments.json")
        request_path.write_text(json.dumps(requests))
        command = self.assigner["command"]
        interpreter = Path(command[0]).parent / "python"
        if len(command) != 1 or not interpreter.is_file():
            raise DatasetError("Batched Pasteur reassignment requires the installed plincer virtual environment")
        _assignment_progress(adapted, starting=True)
        raw_assignments = _run([str(interpreter), str(Path(__file__).with_name("pasteur_assignment_worker.py")),
                                str(request_path), str(self.assigner["scheme_file"]),
                                str(self.assigner["profiles_file"]), str(self.assigner["alleles_db"])],
                               b"", result_path)
        if not isinstance(raw_assignments, dict) or len(raw_assignments.get("assignments", [])) != len(adapted):
            raise DatasetError("Batched Pasteur assignment did not reconcile with requested profiles")
        for row, raw in zip(adapted, raw_assignments["assignments"], strict=True):
            assignment = _assignment_result(raw, self.assigner)
            # normalise_assignment is also an export-row parser and emits an
            # empty source ID for local predictions. It cannot change identity.
            assignment.pop("source_genome_id", None)
            row.update(assignment)
            row.update(cglin_database_sha256=self.assignment_digest,
                       cglin_database_version=self.assigner["scheme_version"])
        _assignment_progress(adapted)
        audit = {"scheme": CANONICAL_SCHEME, "comparison_version": COMPARISON_VERSION,
                 "canonical_panel_sha256": self.panel_digest, "canonical_loci": self.loci,
                 "reference": reference_descriptor(self.config_path),
                 "assignment_reference": self.assignment_provenance,
                 "allele_fasta_sha256": self.source_files,
                 "records": evidence,
                 "policy": "Canonical Pasteur locus-specific DNA identities; genuine configured plincer reassignment; original job/release evidence retained"}
        (self.output / (batch + "-normalization.json")).write_text(json.dumps(audit, indent=2) + "\n")
        self.audits.append({"path": str(self.output / (batch + "-normalization.json")),
                            "records": len(adapted), "unsupported_calls": sum(len(r["unsupported_calls"]) for r in evidence)})
        return adapted


def _assignment_progress(rows, *, starting=False):
    """Log known batch size and actual statuses without modifying results."""
    from chronoclade.progress import progress

    if starting:
        noun = "profile" if len(rows) == 1 else "profiles"
        progress("typing", f"Assigning cgLIN to {len(rows)} {noun} using the configured Pasteur reference")
        return
    resolved = sum(row.get("cglin_status") == "complete" for row in rows)
    provisional = sum(row.get("cglin_status") == "provisional" for row in rows)
    progress("typing", f"cgLIN assignment complete: {resolved} resolved, {provisional} provisional, {len(rows) - resolved - provisional} unassigned")
