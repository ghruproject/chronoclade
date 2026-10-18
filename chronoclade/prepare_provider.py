"""Query-only provider integration and complete public profile context retrieval.

Public searches never authenticate. Only exact resolved query membership and
analysis exports use configured credentials. No context assembly is acquired.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

from chronoclade.datasets import DatasetError, LocusCatalogue
from chronoclade.artifacts import write_json
from chronoclade.datasets.model import check_json
from chronoclade.pathogenwatch import PathogenwatchClient, normalize_country
from chronoclade.profile_inputs import (
    _accession_input,
    _analysis_export_eligible,
    _analysis_exports,
    _apply_overrides,
    _canonical,
    _collection,
    _csv,
    _frozen_cglin,
    _grouped_analysis_exports,
    _join_accession,
    collection_id,
    materialise_assemblies,
)
from chronoclade.typing_scopes import typing_scope


def _save(path, value):
    check_json(value, location=path.name)
    path.write_text(json.dumps(value, indent=2) + "\n")


def _canonical_json_sha256(value):
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _save_accession_evidence(directory, filename, value):
    path = directory / filename
    _save(path, value)
    return {
        "path": path.name,
        "canonical_json_sha256": _canonical_json_sha256(value),
        "bytes": path.stat().st_size,
    }


def _has_profile(row):
    return bool(row.get("cgmlst_profile") or row.get("cgmlst_novel_alleles"))


def _resume_checkpoint_ids(cache_path):
    """An absent cache means typing never started; an existing bad cache is rejected."""
    from chronoclade.query_typing import typing_checkpoint_sample_ids

    cache_path = Path(cache_path)
    if cache_path.is_symlink():
        raise DatasetError("Typing checkpoint cache path is a symlink; refusing resume")
    if not cache_path.exists():
        return []
    if not cache_path.is_dir():
        raise DatasetError("Typing checkpoint cache is not a directory; refusing resume")
    return typing_checkpoint_sample_ids(cache_path)


def _has_lineage(row):
    kind = "cglin" if "klebsiella" in row.get("species", "").casefold() else "hiercc"
    return bool(
        typing_scope(row, kind)
        and row.get("cglin_raw" if kind == "cglin" else "hiercc_codes")
        and row.get(kind + "_status") in {"complete", "resolved", "predicted"}
        and not row.get(kind + "_provisional")
    )


def apply_metadata_overrides(rows, overrides, base):
    """Keep field-specific precedence rather than labelling every field as user supplied."""
    sources = {}
    for row in rows:
        fields = set()
        aliases = set(row.get("aliases") or []) | {row["sample_id"], row.get("source_genome_id")}
        for override in overrides:
            ident = (
                override.get("source_genome_id")
                or override.get("accession")
                or override.get("sample_id")
            )
            if ident in aliases:
                fields.update(
                    k
                    for k, value in override.items()
                    if k not in {"source_genome_id", "accession", "sample_id"}
                    and value is not None
                    and (not isinstance(value, str) or value.strip())
                )
        sources[row.get("source_genome_id") or row["sample_id"]] = fields
    updated = _apply_overrides(rows, overrides, base)
    for row in updated:
        row["prepare_user_metadata_fields"] = sorted(
            sources.get(row.get("source_genome_id") or row["sample_id"], set())
        )
    return updated


def _resolve_accessions(path, species, client, work=None):
    identifiers, overrides = _accession_input(path)
    candidates = []
    organisms = None
    evidence_dir = Path(work) / "accession_resolution" if work is not None else None
    evidence_records = []

    def persist_evidence():
        if evidence_dir is not None:
            _save(
                evidence_dir / "receipt.json",
                {
                    "schema": "chronoclade.accession-resolution-evidence",
                    "schema_version": 1,
                    "hash_definition": "sha256 of UTF-8 canonical decoded JSON; not raw transport bytes",
                    "records": evidence_records,
                },
            )

    if evidence_dir is not None:
        evidence_dir.mkdir(parents=True, exist_ok=True)

    for input_index, ident in enumerate(identifiers):
        prefix = f"{input_index:04d}"
        is_uuid = bool(re.fullmatch(r"[1-9A-HJ-NP-Za-km-z]{20,22}", ident))
        record = {
            "input_index": input_index,
            "identifier": ident,
            "search_mode": "exact_uuid" if is_uuid else "public_accession_search",
            "details": [],
        }
        evidence_records.append(record)
        if is_uuid:
            searches = [{"uuid": ident}]
            record.update(reported_count=1, returned_count=1, result_ids=[ident])
        else:
            organisms = organisms or client.supported_organisms()
            organism = next(
                (
                    r
                    for r in organisms
                    if r.get("fullName", "").casefold() == (species or "").casefold()
                ),
                None,
            )
            if organism is None:
                raise DatasetError("Accession queries require an exact supported --species")
            request = {
                "method": "POST",
                "endpoint": "/api/search/genomes",
                "body": {
                    "organismId": str(organism["organismId"]),
                    "searchText": ident,
                    "qc": [True, False],
                },
                "params": {"limit": 100, "sort": "id"},
            }
            if evidence_dir is not None:
                record["request"] = _save_accession_evidence(
                    evidence_dir, f"{prefix}-search-request.json", request
                )
            response = client.request_json(
                "POST",
                "/api/search/genomes",
                body=request["body"],
                params=request["params"],
            )
            count = response.get("meta", {}).get("count")
            searches = response.get("genomes")
            record.update(
                reported_count=count if type(count) is int else None,
                returned_count=len(searches) if isinstance(searches, list) else None,
                result_ids=[item.get("uuid") for item in searches if isinstance(item, dict)]
                if isinstance(searches, list)
                else [],
            )
            if evidence_dir is not None:
                record["response"] = _save_accession_evidence(
                    evidence_dir, f"{prefix}-search-response.json", response
                )
            record["status"] = "search_response_captured"
            persist_evidence()
            if type(count) is not int or not isinstance(searches, list) or count != len(searches):
                record["status"] = "rejected_incomplete_search"
                persist_evidence()
                raise DatasetError("Accession search is incomplete; supply an exact source UUID")
        record.setdefault("status", "exact_uuid")
        if evidence_dir is not None and is_uuid:
            persist_evidence()
        for search_index, search in enumerate(searches):
            detail = client.request_json(
                "GET", "/api/genomes/details", params={"id": search["uuid"]}
            )
            if evidence_dir is not None:
                detail_evidence = _save_accession_evidence(
                    evidence_dir, f"{prefix}-detail-{search_index:04d}.json", detail
                )
                detail_evidence["source_genome_id"] = search["uuid"]
                detail_evidence["result_index"] = search_index
                record["details"].append(detail_evidence)
                record["status"] = "detail_responses_captured"
                persist_evidence()
            if (
                detail.get("uuid") != search["uuid"]
                or detail.get("projectAccess", "PUBLIC") != "PUBLIC"
            ):
                raise DatasetError("Exact accession detail identity/public-access mismatch")
            row = _canonical({**search, **detail}, origin="local", species=species or "")
            row["organism_id"] = str(detail.get("organismId") or "")
            if row["source_genome_id"] not in {r["source_genome_id"] for r in candidates}:
                candidates.append(row)
    rows = [_join_accession(ident, candidates) for ident in identifiers]
    if overrides:
        rows = apply_metadata_overrides(rows, overrides, path.parent)
    return rows


def enrich_exact_ena(rows, *, fetch=None):
    """Fill missing metadata only from an exact run or BioSample-linked ENA row."""

    def fetch_remote(accession):
        fields = "run_accession,sample_accession,collection_date,country,host,isolation_source"
        url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urlencode(
            {"accession": accession, "result": "read_run", "fields": fields, "format": "json"}
        )
        with urlopen(url, timeout=30) as response:
            return json.load(response)

    fetch = fetch or fetch_remote
    provenance, conflicts, ledger = [], [], []
    for row in rows:
        links = sorted(
            set(row.get("run_accessions") or []) | set(row.get("biosample_accessions") or [])
        )
        if not links:
            ledger.append(
                {
                    "sample_id": row["sample_id"],
                    "operation": "ena_metadata",
                    "status": "unavailable",
                    "reason": "no_exact_linked_accession",
                }
            )
            continue
        found = []
        for accession in links:
            try:
                response = fetch(accession)
                if not isinstance(response, list):
                    raise TypeError("Unsupported ENA response")
                exact = [
                    item
                    for item in response
                    if accession in {item.get("run_accession"), item.get("sample_accession")}
                ]
                if len(exact) != len(response):
                    raise ValueError("ENA response contains an unrelated identity")
                found.extend(exact)
                ledger.append(
                    {
                        "sample_id": row["sample_id"],
                        "operation": "ena_metadata",
                        "accession": accession,
                        "status": "retrieved" if exact else "unavailable",
                        "records": len(exact),
                    }
                )
            except (OSError, ValueError, TypeError):
                ledger.append(
                    {
                        "sample_id": row["sample_id"],
                        "operation": "ena_metadata",
                        "accession": accession,
                        "status": "failed",
                        "reason": "exact_link_metadata_request_failed",
                    }
                )
        for field in ("collection_date", "country", "host", "isolation_source"):
            values = sorted({str(item[field]) for item in found if item.get(field)})
            for value in values:
                selected = len(values) == 1 and row.get(field) in {None, "", "Unknown"}
                provenance.append(
                    {
                        "sample_id": row["sample_id"],
                        "field": field,
                        "source": "ENA_exact_link",
                        "value": value,
                        "selected": selected,
                    }
                )
                if (
                    row.get(field) not in {None, "", "Unknown"}
                    and row[field] != value
                    or len(values) > 1
                ):
                    conflicts.append(
                        {
                            "sample_id": row["sample_id"],
                            "field": field,
                            "source": "ENA_exact_link",
                            "value": value,
                        }
                    )
            if len(values) == 1 and row.get(field) in {None, "", "Unknown"}:
                row[field] = normalize_country(values[0]) if field == "country" else values[0]
                if field == "collection_date":
                    for key in ("date_precision", "date_start", "date_end", "date_raw"):
                        row.pop(key, None)
    return provenance, conflicts, ledger


def catalogues_for_records(rows, explicit=()):
    """Keep complete universes strict; expand incomplete catalogues only by observed loci."""
    catalogues = list(explicit)
    scopes = {}
    for row in rows:
        scope = typing_scope(row, "cgmlst")
        if not scope:
            if _has_profile(row):
                raise DatasetError("Provider returned an unscoped called profile")
            continue
        key = (
            *scope,
            row.get("cgmlst_database_version") or None,
            row.get("cgmlst_database_sha256") or None,
        )
        scopes.setdefault(key, []).append(row)
    for (scheme, version, database_version, digest), members in scopes.items():
        matching = [index for index, catalogue in enumerate(catalogues)
                    if (catalogue.scheme_id, catalogue.scheme_version,
                        catalogue.database_version, catalogue.database_sha256)
                    == (scheme, version, database_version, digest)]
        if matching and all(catalogues[index].complete for index in matching):
            continue
        observed = set()
        for row in members:
            observed.update(row.get("cgmlst_profile") or {})
            observed.update(row.get("cgmlst_novel_alleles") or {})
            observed.update(row.get("cgmlst_loci") or [])
        loci = sorted(observed)
        if matching:
            for index in matching:
                catalogue = catalogues[index]
                extra = sorted(set(loci) - set(catalogue.loci))
                if extra and not catalogue.complete:
                    source = catalogue.source
                    if not source.startswith("Observed same-scope locus union; prior source: "):
                        source = "Observed same-scope locus union; prior source: " + source
                    catalogues[index] = replace(
                        catalogue, loci=(*catalogue.loci, *extra), source=source, complete=False,
                    )
            continue
        declared = [tuple(r["cgmlst_loci"]) for r in members if r.get("cgmlst_loci")]
        if loci:
            complete = bool(
                declared
                and all(set(v) == set(loci) for v in declared)
                and all(r.get("cgmlst_locus_universe_complete") is True for r in members)
            )
            catalogues.append(
                LocusCatalogue(
                    scheme,
                    version,
                    tuple(loci),
                    "provider_frozen_evidence",
                    database_version=database_version,
                    complete=complete,
                    database_sha256=digest,
                )
            )
    return tuple(catalogues)


def assign_existing_profiles(rows, config, output):
    """Run the configured assigner on compatible existing profiles, without an assembly."""
    from chronoclade.query_typing import (
        _assignment_result,
        _provenance,
        _run,
        normalise_typing_species,
    )

    organisms = {normalise_typing_species(k): v for k, v in config["organisms"].items()}
    assigned = set()
    for row in rows:
        organism = organisms.get(normalise_typing_species(row["species"]))
        if not organism or not _has_profile(row) or _has_lineage(row):
            continue
        caller, assigner = organism["cgmlst"], organism["assignment"]
        if typing_scope(row, "cgmlst") != (caller["scheme"], caller["scheme_version"]):
            continue
        caller_provenance = _provenance(caller, ("index_dir",))
        digest = caller_provenance["databases"]["index_dir"]["sha256"]
        if row.get("cgmlst_database_sha256") and row["cgmlst_database_sha256"] != digest:
            continue
        loci = assigner["canonical_loci"]
        calls = dict(row.get("cgmlst_profile") or {}, **(row.get("cgmlst_novel_alleles") or {}))
        if set(calls) - set(loci):
            continue
        ambiguous = set(row.get("cgmlst_ambiguous_loci") or [])
        raw = {
            "scheme": caller["scheme"],
            "schemeSize": len(loci),
            "genes": loci,
            "code": "_".join(
                str(calls.get(locus) or "") if locus not in ambiguous else "" for locus in loci
            ),
            "st": "unassigned",
        }
        directory = output / hashlib.sha256(row["sample_id"].encode()).hexdigest()[:24]
        directory.mkdir(parents=True, exist_ok=True)
        if assigner["kind"] == "plincer":
            args = [
                "classify",
                "-",
                "--scheme-file",
                str(assigner["scheme_file"]),
                "--profiles-file",
                str(assigner["profiles_file"]),
                "--alleles-db",
                str(assigner["alleles_db"]),
            ]
            fields = ("scheme_file", "profiles_file", "alleles_db", "loci_file")
        else:
            args = ["assign", "-", "--exact", "--reference-db", str(assigner["reference_db"])]
            fields = ("reference_db", "loci_file")
        assignment = _run(
            assigner["command"] + args, json.dumps(raw).encode(), directory / "assignment.json"
        )
        row.update(_assignment_result(assignment, assigner))
        row["typing_assignment_evidence"] = assignment
        row["typing_assignment_provenance"] = _provenance(assigner, fields)
        kind = "cglin" if assigner["kind"] == "plincer" else "hiercc"
        row[kind + "_database_sha256"] = hashlib.sha256(
            json.dumps(
                {
                    key: info["sha256"]
                    for key, info in row["typing_assignment_provenance"]["databases"].items()
                },
                sort_keys=True,
            ).encode()
        ).hexdigest()
        assigned.add(row["sample_id"])
    return assigned


def resolve_queries(
    input_value,
    work,
    *,
    input_kind,
    species=None,
    metadata=None,
    client=None,
    typing_config=None,
    query_typing=None,
    public_typing=None,
    cglin_export=None,
    enrich_metadata=True,
    native_checkpoint=None,
    native_attempt=None,
    resume_typing=False,
    cohort_validator=None,
):
    """Resolve and type inputs, without public context discovery."""
    client = client or PathogenwatchClient()
    work.mkdir(parents=True, exist_ok=True)
    provider = {"input_kind": input_kind, "context_discovery": "not_run"}
    if input_kind == "collection":
        rows, details = _collection(client, collection_id(str(input_value)), work)
        provider.update(details)
        rows = _analysis_exports(
            rows, details["downloads"], client, work / "query_exports", provider
        )
    elif input_kind == "accessions":
        rows = _resolve_accessions(Path(input_value).resolve(), species, client, work=work)
        evidence = work / "accession_resolution" / "receipt.json"
        provider["accession_resolution_evidence"] = {
            "path": "accession_resolution/receipt.json",
            "sha256": hashlib.sha256(evidence.read_bytes()).hexdigest(),
            "bytes": evidence.stat().st_size,
        }
    elif input_kind == "assemblies":
        path = Path(input_value).resolve()
        rows = [
            _canonical(r, origin="local", species=species or "", base=path.parent)
            for r in _csv(path)
        ]
        if any(not r.get("assembly") for r in rows):
            raise DatasetError("Assembly metadata requires an existing assembly on every row")
    else:
        raise DatasetError("Unsupported live preparation input kind")
    if metadata:
        rows = apply_metadata_overrides(rows, _csv(Path(metadata)), Path(metadata).resolve().parent)
    if not rows or len({r["sample_id"] for r in rows}) != len(rows):
        raise DatasetError("Inputs require nonempty unique exact sample identities")
    if species:
        for row in rows:
            if row.get("species") and row["species"].casefold() != species.casefold():
                raise DatasetError("Input species differs from --species")
            row["species"] = species
    if any(not r.get("species") for r in rows):
        raise DatasetError("Each input needs --species or a provider species identity")
    if cohort_validator is not None:
        cohort_validator(rows)
    provider_queries_requested = input_kind != "assemblies" or (
        bool(client.api_key)
        and any(
            _analysis_export_eligible(row, name)
            for row in rows
            for name in ("cgmlst", "klebsiella-lincodes", "hclink")
        )
    )
    rows = _grouped_analysis_exports(rows, client, work / "query_exports_grouped", provider)
    rows = _grouped_analysis_exports(
        rows, client, work / "query_lineages", provider, download_names={"klebsiella-lincodes"}
    )
    rows = _grouped_analysis_exports(
        rows, client, work / "query_hiercc", provider, download_names={"hclink"}
    )
    from chronoclade.provider_hiercc import without_internal_proof
    rows = without_internal_proof(rows)
    if public_typing:
        from chronoclade.public_typing import annotate_public_typing, load_public_typing

        rows = annotate_public_typing(rows, load_public_typing(Path(public_typing)))
    if cglin_export:
        from chronoclade.cglin import load_cglin_export

        rows = _frozen_cglin(rows, load_cglin_export(Path(cglin_export)))
    if query_typing:
        from chronoclade.query_typing import load_query_typing

        assignments = {r["sample_id"]: r for r in load_query_typing(Path(query_typing))}
        if set(assignments) - {r["sample_id"] for r in rows}:
            raise DatasetError("Query typing contains an identity outside the input cohort")
        for row in rows:
            assignment = assignments.get(row["sample_id"])
            if assignment:
                if assignment["species"].casefold() != row["species"].casefold() or (
                    row.get("assembly_sha256")
                    and row["assembly_sha256"] != assignment["assembly_sha256"]
                ):
                    raise DatasetError(
                        "Query typing species/assembly hash differs from the exact input"
                    )
                row.update(assignment)
    operations = {
        "provider_retrieval": "run_exact_queries"
        if provider_queries_requested
        else "not_run_local_assemblies",
        "context_discovery": "not_run",
        "typing": "not_run",
        "lineage_assignment": "not_run",
        "typing_tool_installation": "not_configured",
        "reference_database": "unavailable_not_configured",
    }
    ledger = []
    if typing_config:
        from chronoclade.query_typing import (
            QueryTypingError,
            load_typing_config,
            type_query_assemblies,
        )

        try:
            config = load_typing_config(Path(typing_config))
        except QueryTypingError:
            operations.update(
                reference_database="unavailable_configured_database_or_tool_not_ready",
                typing_tool_installation="configured_not_ready",
            )
            setup_manifest = Path(typing_config).resolve().parent / "typing_setup.json"
            if setup_manifest.is_file():
                setup = json.loads(setup_manifest.read_text())
                if all(
                    t.get("status") == "installed" for t in setup.get("tools", {}).values()
                ) and setup.get("tools"):
                    operations["typing_tool_installation"] = (
                        "installed_reference_database_unavailable"
                    )
        else:
            operations.update(
                reference_database="ready_validated_configuration",
                typing_tool_installation="configured_available",
            )
            if native_checkpoint is not None and resume_typing:
                from chronoclade.query_typing import validate_typing_checkpoint

                checkpoint_path = Path(native_checkpoint)
                cached_ids = _resume_checkpoint_ids(checkpoint_path)
                if cached_ids:
                    checkpoint_rows = [row for row in rows if row["sample_id"] in cached_ids]
                    if [row["sample_id"] for row in checkpoint_rows] != cached_ids:
                        raise DatasetError("Typing checkpoint contains samples outside the resolved cohort")
                    checkpoint_samples = materialise_assemblies(
                        checkpoint_rows, output=work / "query_assemblies_preflight", client=client
                    )
                    validate_typing_checkpoint(
                        checkpoint_samples, Path(typing_config), Path(native_checkpoint)
                    )
            assigned = assign_existing_profiles(rows, config, work / "native_lineage_assignment")
            if assigned:
                operations["lineage_assignment"] = "run_existing_compatible_profiles"
            need = [
                r
                for r in rows
                if r["sample_id"] not in assigned and (not _has_profile(r) or not _has_lineage(r))
            ]
            if need:
                samples = materialise_assemblies(
                    need, output=work / "query_assemblies", client=client
                )
                if native_checkpoint is None:
                    assignments = type_query_assemblies(
                        samples, Path(typing_config), work / "native_typing"
                    )
                else:
                    assignments = type_query_assemblies(
                        samples,
                        Path(typing_config),
                        Path(native_checkpoint),
                        # A prepare resume can fail before native typing ever
                        # creates its per-sample checkpoint directory. In that
                        # case this is the first native typing attempt, while
                        # an existing cache must still be resumed and validated.
                        resume=resume_typing and Path(native_checkpoint).exists(),
                        attempt_output=Path(native_attempt),
                    )
                    # Retain the established portable export as well as the
                    # prepare-owned recovery cache, which is removed on success.
                    write_json(work / "native_typing" / "query_typing.json", {
                        "schema_version": 1,
                        "created_at": datetime.now(timezone.utc).isoformat(),
                        "typing_mode": "full",
                        "assignments": assignments,
                    })
                by_id = {r["sample_id"]: r for r in assignments}
                for row, sample in zip(need, samples, strict=True):
                    row.update(by_id[row["sample_id"]])
                    row["assembly"] = str(sample.assembly)
                    # Caller supplies its ordered full universe, distinct from an observed export.
                    row["cgmlst_locus_universe_complete"] = True
                operations.update(
                    typing="run_missing_profile_or_lineage_inputs",
                    lineage_assignment="run_missing_profile_or_lineage_inputs",
                )
    for row in rows:
        row.update(role="input", identity_resolved=bool(row.get("source_genome_id")))
        ledger.append(
            {
                "sample_id": row["sample_id"],
                "operation": "prepare_typing",
                "status": "ready" if _has_profile(row) and _has_lineage(row) else "incomplete",
                "reason": None
                if _has_profile(row) and _has_lineage(row)
                else operations["reference_database"],
            }
        )
    provenance = [
        {
            "sample_id": r["sample_id"],
            "field": key,
            "source": "user"
            if key in r.get("prepare_user_metadata_fields", []) or input_kind == "assemblies"
            else "input_provider",
            "value": r[key],
            "selected": True,
        }
        for r in rows
        for key in (
            "country", "region", "city", "nuts2", "host", "isolation_source", "collection_date",
            "hospital", "hospital_code", "publication_isolate_id", "specimen_type",
            "clinical_specimen", "study_accession",
        )
        if r.get(key)
    ]
    conflicts = [
        {
            "sample_id": r["sample_id"],
            "field": c["field"],
            "source": "provider",
            "value": c["source_value"],
        }
        for r in rows
        for c in r.get("metadata_override_conflicts", [])
    ]
    if enrich_metadata:
        extra, clashes, ena_ledger = enrich_exact_ena(rows)
        provenance.extend(extra)
        conflicts.extend(clashes)
        ledger.extend(ena_ledger)
    _save(
        work / "resolved_queries.json",
        {
            "records": rows,
            "provenance": provenance,
            "conflicts": conflicts,
            "retrieval": ledger,
            "parameters": {"provider_prepare": provider, "prepare_operations": operations},
        },
    )
    return rows, provenance, conflicts, ledger, provider, operations
