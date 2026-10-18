"""Complete accessible public profile context, with explicit lineage scope and policy."""

from __future__ import annotations

import json
import shutil
from contextlib import ExitStack, nullcontext
from dataclasses import asdict
from pathlib import Path

from chronoclade.artifacts import file_sha256
from chronoclade.execution_provenance import execution_provenance
from chronoclade.cglin import parse_code, resolved_prefix
from chronoclade.context_partitions import LIN_LEVELS, validate_lin_level
from chronoclade.datasets import DatasetError, from_profile_records, load_dataset
from chronoclade.pathogenwatch import (
    PathogenwatchClient,
    content_hash,
    deduplicate_catalogue,
    load_catalogue,
)
from chronoclade.prepare_provider import _has_profile, _save, catalogues_for_records
from chronoclade.profile_inputs import _canonical, _frozen_cglin, _grouped_analysis_exports
from chronoclade.typing_scopes import compatible_typing
from chronoclade.progress import progress_task


def dataset_records(dataset, *, source_root=None):
    """Adapt authoritative typing without allowing evidence to overwrite identity/metadata."""
    from chronoclade.context_partitions import profile_records

    rows = profile_records(dataset, source_root=source_root)
    by_id = {row["sample_id"]: row for row in rows}
    allowed = {
        "source_genome_id",
        "numeric_source_id",
        "source_checksum",
        "source_length",
        "source_n50",
        "source_contigs",
        "qc_pass",
        "source_qc",
        "aliases",
        "accession",
        "run_accession",
        "runAccession",
        "sample_accession",
        "sampleAccession",
        "assembly_accession",
        "assemblyAccession",
        "biosample",
        "run_accessions",
        "biosample_accessions",
        "assembly_accessions",
        "identity_candidates",
        "identity_resolved",
        "identity_conflict",
    }
    for evidence in dataset.crosswalk:
        by_id[evidence["sample_id"]].update({k: v for k, v in evidence.items() if k in allowed})
    return rows


def _matching(query, candidate, kind, level):
    from chronoclade.provider_hiercc import grouping_scope
    left_reported, right_reported = grouping_scope(query), grouping_scope(candidate)
    from chronoclade.typing_scopes import typing_scope
    if kind == "hiercc" and not left_reported and not right_reported and (typing_scope(query, kind) is None or typing_scope(candidate, kind) is None):
        return False
    reported = kind == "hiercc" and left_reported is not None and left_reported == right_reported
    if ((not reported and not compatible_typing(query, candidate, kind))
        or query.get(kind + "_database_version") != candidate.get(kind + "_database_version")
        or query.get(kind + "_database_sha256") != candidate.get(kind + "_database_sha256")
        or (kind == "hiercc" and not reported and (left_reported or right_reported))
        or any(row.get(kind + "_status") in
               {"unassigned", "unsupported", "conflict", "failed", "unverified", "unavailable"}
               for row in (query, candidate))):
        return False
    if kind == "cglin":
        left, right = resolved_prefix(query, level), resolved_prefix(candidate, level)
        return left is not None and right is not None and left == right
    left, right = query.get("hiercc_codes") or {}, candidate.get("hiercc_codes") or {}
    return bool(left.get(level) and left.get(level) == right.get(level))


def _numeric_lineage_match(query, candidate, kind, level):
    """Discovery hint only: identical numbers never establish reference compatibility."""
    if kind == "cglin":
        left = parse_code(query.get("cglin_raw"))[0]
        right = parse_code(candidate.get("cglin_raw"))[0]
        return len(left) >= level and len(right) >= level and left[:level] == right[:level]
    left = query.get("hiercc_codes") or {}
    right = candidate.get("hiercc_codes") or {}
    return bool(left.get(level) and left.get(level) == right.get(level))


def _scope_evidence(row, kind):
    evidence = {key: row.get(kind + "_" + key) for key in (
        "scheme", "scheme_version", "database_version", "database_sha256"
    )}
    if kind == "hiercc":
        from chronoclade.provider_hiercc import grouping_scope
        evidence["provider_reported_scope"] = grouping_scope(row)
        evidence["provenance_limitation"] = row.get("hiercc_reported_grouping_limitation")
    return evidence


def _fail_incompatible_context(work, *, kind, candidates, scope=None):
    """Retain refusal evidence instead of presenting incompatibility as no relatives."""
    path = work / "context_compatibility.json"
    _save(path, {
        "status": "failed_incompatible_typing_scope",
        "kind": kind,
        "scope": scope,
        "rejected_candidates": candidates,
        "context_assemblies": "not_requested",
        "interpretation": "Compatibility is unproved; absence of public relatives is not established.",
    })
    raise DatasetError(
        f"Public {kind} candidates were found, but their typing scheme/reference scopes "
        "cannot be compared with the inputs. This does not mean there are no public relatives. "
        "Use profiles and lineage assignments from the same verified reference scope, or "
        "supply a validated reference translation; changing scheme labels is insufficient. "
        f"Compatibility audit: {path}"
    )


def _exports_complete(provenance):
    """Whether attempted public exports completed without omissions or failures."""
    if provenance.get("analysis_export_status", "complete") != "complete":
        return False
    if provenance.get("profile_export_unavailable"):
        return False
    return all(
        item.get("status") == "downloaded"
        for item in provenance.get("analysis_exports", [])
    )


def discover_profile_context(
    dataset_manifest,
    output,
    *,
    lin_level=5,
    hiercc_level=None,
    client=None,
    public_typing=None,
    cglin_export=None,
    catalogues=None,
    typing_config=None,
    context_source=None,
    resume=False,
):
    """Publish all accessible compatible block profiles; never fetch context assemblies."""
    from chronoclade.prepare_stage import (
        PrepareResult,
        _catalogues,
        _validate_previous_stage,
        run_prepare,
    )

    dataset = load_dataset(dataset_manifest, mmap_mode=None)
    from chronoclade.pasteur_profile_normalization import (
        PasteurProfileNormalizer, reference_descriptor, reference_rebinding, resolve_reference,
    )
    native_reference = resolve_reference(dataset, typing_config)
    native_descriptor = reference_descriptor(native_reference) if native_reference else None
    rebinding = reference_rebinding(dataset.parameters.get("native_typing_reference"), native_descriptor,
                                   explicit=typing_config is not None) if native_descriptor else None
    try:
        validate_lin_level(lin_level, label="Context LIN level")
    except ValueError as error:
        raise DatasetError(str(error)) from error
    rows = dataset_records(dataset, source_root=dataset_manifest)
    for row in rows:
        reference = row.get("assembly_reference")
        if reference and "://" not in reference:
            source_root = Path(dataset_manifest).resolve()
            if source_root.is_file():
                source_root = source_root.parent
            row["assembly_reference"] = str((source_root / reference).resolve())
    queries = [r for r in rows if r["role"] == "input"]
    if any(
        r.get("species", "").casefold() not in {"klebsiella pneumoniae", "escherichia coli"}
        for r in queries
    ):
        raise DatasetError(
            "Public lineage context supports Klebsiella pneumoniae and Escherichia coli"
        )
    if any("coli" in r.get("species", "").casefold() for r in queries) and not hiercc_level:
        raise DatasetError("E. coli public context requires an explicit --hiercc-level")
    raw_source = None
    if context_source is not None:
        from chronoclade.provider_sources import ProviderSource
        if public_typing or cglin_export:
            raise DatasetError("Frozen raw context source cannot be combined with normalized typing/lineage imports")
        if not native_reference or any(r["species"].casefold() != "klebsiella pneumoniae" for r in queries):
            raise DatasetError("Frozen raw source replay currently requires a verified Klebsiella native normalization reference")
        raw_source = ProviderSource(context_source)
        if client is not None and client.base_url.rstrip("/") != raw_source.base_url:
            raise DatasetError("Frozen provider source base URL differs from requested provider")
        for species, st in {(r["species"], str(r.get("mlst_st") or "")) for r in queries}:
            raw_source.for_scope(species, st)
    else:
        client = client or PathogenwatchClient()
    provider_base_url = raw_source.base_url if raw_source else getattr(client, "base_url", "https://pathogen.watch").rstrip("/")
    target = Path(output).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)

    def input_sha256(path):
        if path is None:
            return None
        candidate = Path(path).resolve()
        if not candidate.is_file():
            raise DatasetError(f"Resume input must be a file: {candidate}")
        return file_sha256(candidate)

    source_manifest = (
        Path(dataset_manifest) / "dataset.json"
        if Path(dataset_manifest).is_dir()
        else Path(dataset_manifest)
    )
    fingerprint = {
        "execution": execution_provenance("profile-context", "frozen-scope-export-v1"),
        "source_dataset_sha256": file_sha256(source_manifest),
        "lin_level": lin_level,
        "hiercc_level": hiercc_level,
        "public_typing_sha256": input_sha256(public_typing),
        "cglin_export_sha256": input_sha256(cglin_export),
        "catalogues_sha256": input_sha256(catalogues),
        "native_typing_reference": native_descriptor,
        "provider_base_url": provider_base_url,
        "raw_provider_source": raw_source.identity if raw_source else None,
    }
    retrieval = target.with_name(target.name + "_retrieval")
    if target.exists():
        if not resume:
            raise DatasetError(f"Preparation output already exists: {target}")
        try:
            _validate_previous_stage(target)
            manifest = target / "dataset.json"
            saved_dataset = load_dataset(manifest, mmap_mode=None)
            stage = json.loads((target / "prepare.json").read_text())
            readiness = json.loads((target / "readiness.json").read_text())
            saved_fingerprint = (saved_dataset.parameters.get("profile_context") or {}).get(
                "resume_fingerprint"
            )
            if (
                stage.get("status") != "complete"
                or stage.get("dataset_manifest") != "dataset.json"
                or saved_fingerprint != fingerprint
                or (readiness.get("provider_context") or {}).get("resume_fingerprint")
                != fingerprint
            ):
                raise DatasetError("Completed context preparation does not match resume inputs")
        except DatasetError:
            raise
        except (OSError, ValueError, KeyError, TypeError) as error:
            raise DatasetError("Completed context preparation failed validation") from error
        return PrepareResult(
            manifest,
            target / "report.html",
            target / "readiness.json",
            target / "prepare.json",
        )
    if resume:
        retrieval.mkdir(parents=True, exist_ok=True)
        state_path = retrieval / "resume.json"
        if state_path.exists():
            try:
                state = json.loads(state_path.read_text())
            except (OSError, ValueError) as error:
                raise DatasetError("Cannot read context retrieval resume state") from error
            if state.get("fingerprint") != fingerprint:
                raise DatasetError("Context retrieval resume inputs changed")
        elif any(retrieval.iterdir()):
            raise DatasetError("Context retrieval directory has unverified existing files")
        else:
            state_tmp = state_path.with_suffix(".json.tmp")
            _save(state_tmp, {"schema_version": 1, "fingerprint": fingerprint})
            state_tmp.replace(state_path)
    else:
        if retrieval.exists():
            raise DatasetError(
                f"Context retrieval cache already exists: {retrieval}; use resume to continue"
            )
        retrieval.mkdir(parents=True)
        state_path = retrieval / "resume.json"
        state_tmp = state_path.with_suffix(".json.tmp")
        _save(state_tmp, {"schema_version": 1, "fingerprint": fingerprint})
        state_tmp.replace(state_path)

    work_context = nullcontext(str(retrieval))

    with work_context as temporary, ExitStack() as phases:
        work = Path(temporary)
        if raw_source:
            raw_source.snapshot(work / "raw_source")
        normalizer = PasteurProfileNormalizer(native_reference, work / "normalization") if native_reference else None
        if normalizer:
            rows = normalizer.normalize(rows)
            queries = [r for r in rows if r["role"] == "input"]
        found, pools = {}, []
        scopes = {(r["species"], str(r.get("mlst_st") or "")) for r in queries}
        organisms = None
        for index, (species, st) in enumerate(sorted(scopes)):
            checkpoint_path = work / f"matched_{index}.json"
            scope = {"index": index, "species": species, "st": st}
            if checkpoint_path.exists():
                try:
                    checkpoint = json.loads(checkpoint_path.read_text())
                except (OSError, ValueError) as error:
                    raise DatasetError(f"Cannot read context scope checkpoint {index}") from error
                body = {k: v for k, v in checkpoint.items() if k != "checkpoint_sha256"}
                if (
                    checkpoint.get("fingerprint") != fingerprint
                    or checkpoint.get("scope") != scope
                    or checkpoint.get("complete") is not True
                    or checkpoint.get("checkpoint_sha256") != content_hash(body)
                    or not isinstance(checkpoint.get("matched"), list)
                ):
                    raise DatasetError(f"Context scope checkpoint {index} failed validation")
                from chronoclade.provider_hiercc import rebind_scope_exports
                try:
                    matched = rebind_scope_exports(
                        checkpoint["matched"], checkpoint["pool_audit"].get("exports", {}),
                        work / f"scope_{index}" / "lineages" / "export_capabilities.json",
                    )
                except (OSError, ValueError, KeyError, TypeError) as error:
                    raise DatasetError(f"Context scope checkpoint {index} raw lineage proof failed validation") from error
                for row in matched:
                    found[row["source_genome_id"]] = row
                pools.append(checkpoint["pool_audit"])
                continue

            if raw_source:
                envelope, raw_rows, raw_assignments = raw_source.for_scope(species, st)
                from chronoclade.progress import progress
                progress("context-source", f"Frozen {species} ST{st}: {len(raw_rows)} public records as of {envelope['provenance']['retrieved_at']}; no provider requests")
                organism = envelope["provenance"]["organism_capabilities"]
            else:
                public_path = work / f"public_{index}.json"
                if public_path.exists():
                    try:
                        envelope = load_catalogue(public_path)
                    except Exception as error:
                        raise DatasetError(f"Frozen context catalogue {index} failed validation") from error
                    provenance = envelope.get("provenance", {})
                    organism = provenance.get("organism_capabilities")
                    query = provenance.get("query", {})
                    if (
                        not organism
                        or organism.get("fullName", "").casefold() != species.casefold()
                        or str(query.get("organismId")) != str(organism.get("organismId"))
                        or str((query.get("mlst") or [None])[0]) != st
                        or provenance.get("base_url", "").rstrip("/") != client.base_url.rstrip("/")
                    ):
                        raise DatasetError(f"Frozen context catalogue {index} does not match its scope")
                else:
                    if organisms is None:
                        organisms = client.supported_organisms()
                    organism = next(
                        (r for r in organisms if r.get("fullName", "").casefold() == species.casefold()),
                        None,
                    )
                if not st.isdigit() or organism is None:
                    pools.append(
                        {
                            "species": species,
                            "st": st,
                            "status": "unavailable",
                            "reason": "missing_supported_species_or_exact_mlst_st",
                        }
                    )
                    continue
                if not public_path.exists():
                    envelope = client.freeze_catalogue(
                        public_path, organism_id=str(organism["organismId"]), st=st
                    )
            pool = [
                _canonical(
                    dict(r, sample_id="PW_" + r["source_genome_id"]),
                    origin="context",
                    species=species,
                )
                for r in (raw_rows if raw_source else envelope["rows"])
            ]
            export = {}
            scope_work = work / f"scope_{index}"
            scope_work.mkdir(exist_ok=resume)
            if raw_source:
                pool = [dict(row, **{key: value for key, value in raw_assignments.get(row["source_genome_id"], {}).items() if key.startswith("cgmlst_")}) for row in pool]
                export["raw_provider_source"] = raw_source.identity
                export["as_of"] = envelope["provenance"]["retrieved_at"]
                export["analysis_exports"] = [{"download": "cgmlst", "status": "downloaded", "origin": "validated_raw_source_replay", "requested_count": len(pool)}]
            else:
                pool = _grouped_analysis_exports(
                    pool, client, scope_work / "lineages", export,
                    download_names={"hclink"} if species.casefold() == "escherichia coli"
                    else {"klebsiella-lincodes"}, resume=resume,
                )
            if public_typing:
                from chronoclade.public_typing import annotate_public_typing, load_public_typing

                pool = annotate_public_typing(pool, load_public_typing(Path(public_typing)))
            if cglin_export:
                from chronoclade.cglin import load_cglin_export

                pool = _frozen_cglin(pool, load_cglin_export(Path(cglin_export)))
            relevant = [
                r for r in queries if r["species"] == species and str(r.get("mlst_st") or "") == st
            ]
            kind = "cglin" if "klebsiella" in species.casefold() else "hiercc"
            counts = {
                str(level): sum(any(_matching(q, r, kind, level) for q in relevant) for r in pool)
                for level in (LIN_LEVELS if kind == "cglin" else [hiercc_level])
            }
            matched = [
                r
                for r in pool
                if any(
                    _matching(q, r, kind, lin_level if kind == "cglin" else hiercc_level)
                    for q in relevant
                )
            ]
            selected_level = lin_level if kind == "cglin" else hiercc_level
            rejected_lineages = []
            for candidate in pool:
                comparable = any(_matching(q, candidate, kind, selected_level) for q in relevant)
                hinted_queries = [q for q in relevant if _numeric_lineage_match(
                    q, candidate, kind, selected_level
                )]
                incompatible_scope = hinted_queries and all(
                    _scope_evidence(q, kind) != _scope_evidence(candidate, kind)
                    for q in hinted_queries
                )
                if incompatible_scope and not comparable:
                    rejected_lineages.append({
                        "source_genome_id": candidate.get("source_genome_id"),
                        "candidate_scope": _scope_evidence(candidate, kind),
                        "input_scopes": [_scope_evidence(q, kind) for q in hinted_queries],
                        "reason": "numeric_lineage_match_not_verified_comparable",
                    })
            if rejected_lineages and not matched and not normalizer:
                _fail_incompatible_context(
                    work, kind=kind, candidates=rejected_lineages, scope=scope
                )
            if normalizer:
                # Old public assignments cannot safely exclude candidates under
                # the new reference. Acquire ST profiles, then reassign BEFORE
                # applying the requested current-reference LIN prefix.
                matched = pool
            requested = len(matched)
            if not raw_source:
                matched = _grouped_analysis_exports(matched, client, scope_work / "profiles", export, resume=resume)
            reassignment_exclusions = []
            if normalizer:
                available = [r for r in matched if _has_profile(r)]
                normalized = normalizer.normalize(available)
                counts = {str(level): sum(any(_matching(q, r, kind, level) for q in relevant)
                                         for r in normalized)
                          for level in LIN_LEVELS}
                matches_after_assignment = []
                for candidate in normalized:
                    if any(_matching(q, candidate, kind, selected_level) for q in relevant):
                        matches_after_assignment.append(candidate)
                    else:
                        reassignment_exclusions.append({
                            "source_genome_id": candidate["source_genome_id"],
                            "reason": "current_reference_lineage_outside_requested_scope_or_unassigned",
                        })
                reassignment_exclusions.extend({
                    "source_genome_id": candidate["source_genome_id"], "reason": "profile_unavailable"
                } for candidate in matched if not _has_profile(candidate))
                matched = matches_after_assignment
            for row in matched:
                row.pop("assembly", None)
                row.pop("assembly_reference", None)
                found[row["source_genome_id"]] = row
            pool_audit = {
                "species": species,
                "st": st,
                "public_search": envelope["provenance"],
                "level_counts": counts,
                "discovered": len(pool),
                "requested": requested,
                "retrieved": sum(_has_profile(r) for r in matched),
                "exports": export,
                "incompatible_numeric_lineage_candidates": rejected_lineages,
                "current_reference_reassignment": {
                    "enabled": normalizer is not None,
                    "prefilter": "all_ST_profiles_before_current_LIN_scope" if normalizer else "reported_compatible_lineage",
                    "retained": len(matched), "exclusions": reassignment_exclusions,
                    "profiles_before_lineage_scope": len(available) if normalizer else None,
                },
            }
            pools.append(pool_audit)
            if _exports_complete(export) and all(_has_profile(row) for row in matched):
                checkpoint = {
                    "complete": True,
                    "fingerprint": fingerprint,
                    "scope": scope,
                    "matched": matched,
                    "pool_audit": pool_audit,
                }
                checkpoint["checkpoint_sha256"] = content_hash(checkpoint)
                temporary_checkpoint = checkpoint_path.with_suffix(".json.tmp")
                temporary_checkpoint.write_text(json.dumps(checkpoint, sort_keys=True, indent=2))
                temporary_checkpoint.replace(checkpoint_path)
        phases.enter_context(progress_task(
            "context", "Preparing combined input and public-context profiles",
            details=f"{len(rows)} existing records; {len(found)} retrieved lineage matches",
        ))
        aliases = set()
        for row in rows:
            aliases.update(row.get("aliases") or [])
            for key in (
                "source_genome_id",
                "sample_id",
                "biosample",
                "sample_accession",
                "assembly_accession",
            ):
                if row.get(key):
                    aliases.add(row[key])
        context_candidates, exclusions = [], []
        for source_id, row in sorted(found.items()):
            if source_id in aliases or set(row.get("aliases") or []) & aliases:
                exclusions.append(
                    {"source_genome_id": source_id, "reason": "already_present_exact_identity"}
                )
            elif not _has_profile(row):
                exclusions.append({"source_genome_id": source_id, "reason": "profile_unavailable"})
            elif not any(
                compatible_typing(q, row, "cgmlst")
                and q.get("cgmlst_database_version") == row.get("cgmlst_database_version")
                and (q.get("cgmlst_database_sha256") or None)
                == (row.get("cgmlst_database_sha256") or None)
                for q in queries
            ):
                exclusions.append(
                    {"source_genome_id": source_id, "reason": "incompatible_profile_scope"}
                )
            else:
                row.update(role="context", identity_resolved=True)
                context_candidates.append(row)
        # Keep input records intact. Collapse only retrieved context after the
        # profile and lineage compatibility gates, with every source record
        # represented in the deterministic deduplication audit.
        context, context_deduplication = deduplicate_catalogue(context_candidates)
        incompatible_profiles = [e for e in exclusions if e["reason"] == "incompatible_profile_scope"]
        if incompatible_profiles and not context:
            _fail_incompatible_context(work, kind="cgmlst", candidates=incompatible_profiles)
        from chronoclade.provider_hiercc import source_snapshot, without_internal_proof
        reported_evidence = source_snapshot(rows + context)
        rows = without_internal_proof(rows)
        context = without_internal_proof(context)
        all_rows = rows + context
        explicit = [m.catalogue for m in dataset.profiles] + list(
            _catalogues(Path(catalogues)) if catalogues else []
        )
        catalogue_values = catalogues_for_records(all_rows, explicit)
        audit = {
            "policy": {"lin_level": lin_level, "hiercc_level": hiercc_level},
            "provider_reported_hiercc": reported_evidence,
            "resume_fingerprint": fingerprint,
            "raw_provider_source": raw_source.identity if raw_source else None,
            "source_mode": "frozen_raw_replay" if raw_source else "live_provider",
            "source_as_of": [dict(species=s["scope"]["scope"]["species"], st=s["scope"]["scope"]["st"], retrieved_at=s["scope"]["scope"]["retrieved_at"], base_url=s["scope"]["scope"]["base_url"], expected_count=s["scope"]["scope"]["expected_count"]) for s in raw_source.scopes] if raw_source else [],
            "pools": pools,
            "accessible_matches": len(found),
            "context_candidates": len(context_candidates),
            "context_deduplication": context_deduplication,
            "usable": len(context),
            "exclusions": exclusions,
            "source_dataset_sha256": file_sha256(
                Path(dataset_manifest) / "dataset.json"
                if Path(dataset_manifest).is_dir()
                else Path(dataset_manifest)
            ),
            "context_assemblies": "not_requested",
            "global_completeness": "not_claimed",
            "native_typing_reference": native_descriptor,
            "native_reference_migration": rebinding,
            "normalization": [dict(item, path=str(Path(item["path"]).relative_to(work)))
                              for item in normalizer.audits] if normalizer else [],
        }
        _save(work / "catalogues.json", [asdict(c) for c in catalogue_values])
        lineages = (
            list(from_profile_records(all_rows, catalogue_values).lineages) if normalizer
            else list(dataset.lineages) + (
                list(from_profile_records(context, catalogue_values).lineages) if context else []
            )
        )
        parameters = dict(dataset.parameters, profile_context=audit)
        if native_reference:
            parameters["native_typing_reference"] = native_descriptor
            if rebinding:
                parameters["native_typing_reference_migration"] = rebinding
        _save(
            work / "profiles.json",
            {
                "records": all_rows,
                "lineages": lineages,
                "provenance": list(dataset.provenance)
                + [
                    {
                        "sample_id": r["sample_id"],
                        "field": "source_genome_id",
                        "source": "public_pathogenwatch",
                        "value": r["source_genome_id"],
                    }
                    for r in context
                ],
                "crosswalk": list(dataset.crosswalk)
                + list(from_profile_records(context, catalogue_values).crosswalk)
                if context
                else list(dataset.crosswalk),
                "retrieval": list(dataset.retrieval)
                + list(from_profile_records(context, catalogue_values).retrieval)
                if context
                else list(dataset.retrieval),
                "exclusions": list(dataset.exclusions),
                "conflicts": list(dataset.conflicts),
                "parameters": parameters,
            },
        )
        source_root = Path(dataset_manifest).resolve()
        if source_root.is_file():
            source_root = source_root.parent
        source_copy = work / "source_dataset"
        source_copy.mkdir(exist_ok=True)
        raw_manifest = json.loads((source_root / "dataset.json").read_text())
        source_paths = {"dataset.json"}
        declared_hashes = {"dataset.json": file_sha256(source_root / "dataset.json")}

        def collect_paths(value):
            if isinstance(value, dict):
                if "path" in value and "sha256" in value:
                    source_paths.add(value["path"])
                    previous = declared_hashes.setdefault(value["path"], value["sha256"])
                    if previous != value["sha256"]:
                        raise DatasetError("Source dataset contains conflicting file hashes")
                for child in value.values():
                    collect_paths(child)
            elif isinstance(value, list):
                for child in value:
                    collect_paths(child)

        # Copy only validated authoritative dataset descriptors, excluding legacy
        # parameter paths and deferred assembly references.
        collect_paths(raw_manifest["tables"])
        collect_paths(raw_manifest["profiles"])
        if reported_evidence:
            source_audit = json.loads((source_root / "prepare.json").read_text())
            for artifact in source_audit["artifacts"]:
                if artifact["path"].startswith("sources/provider/"):
                    collect_paths(artifact)

        for name in sorted(source_paths):
            source_file = source_root / name
            destination = source_copy / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            expected_sha = declared_hashes[name]
            if file_sha256(source_file) != expected_sha:
                raise DatasetError(f"Source dataset evidence changed while copying: {name}")
            if destination.exists():
                if file_sha256(destination) != expected_sha:
                    raise DatasetError(f"Saved source dataset copy failed validation: {name}")
                continue
            temporary_copy = destination.with_name(destination.name + ".tmp")
            shutil.copyfile(source_file, temporary_copy)
            if file_sha256(temporary_copy) != expected_sha:
                temporary_copy.unlink(missing_ok=True)
                raise DatasetError(f"Source dataset copy failed validation: {name}")
            temporary_copy.replace(destination)
        if not catalogue_values:
            from chronoclade.datasets import write_dataset
            from chronoclade.prepare_stage import _frozen_dataset

            sample_only = _frozen_dataset(
                json.loads((work / "profiles.json").read_text()), (), None
            )
            manifest = write_dataset(sample_only, work / "sample_only")
            return run_prepare(manifest, target, _provider_sources=work)
        result = run_prepare(
            work / "profiles.json",
            target,
            catalogues=work / "catalogues.json",
            _provider_sources=work,
        )
    return result
