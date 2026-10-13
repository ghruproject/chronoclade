"""Complete accessible public profile context, with explicit lineage scope and policy."""

from __future__ import annotations

import json
import shutil
import uuid
from contextlib import nullcontext
from dataclasses import asdict
from pathlib import Path

from chronoclade.artifacts import file_sha256
from chronoclade.cglin import parse_code
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


def dataset_records(dataset):
    """Adapt authoritative typing without allowing evidence to overwrite identity/metadata."""
    from chronoclade.context_partitions import profile_records

    rows = profile_records(dataset)
    by_id = {row["sample_id"]: row for row in rows}
    allowed = {
        "source_genome_id",
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
    if (
        not compatible_typing(query, candidate, kind)
        or query.get(kind + "_database_version") != candidate.get(kind + "_database_version")
        or query.get(kind + "_database_sha256") != candidate.get(kind + "_database_sha256")
    ) or any(
        row.get(kind + "_status")
        in {"unassigned", "unsupported", "conflict", "failed", "unverified", "unavailable"}
        for row in (query, candidate)
    ):
        return False
    if kind == "cglin":
        left, right = (
            parse_code(query.get("cglin_raw"))[0],
            parse_code(candidate.get("cglin_raw"))[0],
        )
        return (
            len(left) >= level
            and len(right) >= level
            and left[:level] == right[:level]
            and not query.get("cglin_provisional")
            and not candidate.get("cglin_provisional")
        )
    left, right = query.get("hiercc_codes") or {}, candidate.get("hiercc_codes") or {}
    return bool(left.get(level) and left.get(level) == right.get(level))


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
    if lin_level not in {5, 6, 7}:
        raise DatasetError("Context LIN level must be explicitly 5, 6 or 7")
    rows = dataset_records(dataset)
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
    client = client or PathogenwatchClient()
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
        "source_dataset_sha256": file_sha256(source_manifest),
        "lin_level": lin_level,
        "hiercc_level": hiercc_level,
        "public_typing_sha256": input_sha256(public_typing),
        "cglin_export_sha256": input_sha256(cglin_export),
        "catalogues_sha256": input_sha256(catalogues),
        "provider_base_url": getattr(client, "base_url", "https://pathogen.watch").rstrip("/"),
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

    with work_context as temporary:
        work = Path(temporary)
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
                for row in checkpoint["matched"]:
                    found[row["source_genome_id"]] = row
                pools.append(checkpoint["pool_audit"])
                continue

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
                for r in envelope["rows"]
            ]
            export = {}
            scope_work = work / f"scope_{index}_{uuid.uuid4().hex}"
            scope_work.mkdir()
            pool = _grouped_analysis_exports(
                pool,
                client,
                scope_work / "lineages",
                export,
                download_names={"klebsiella-lincodes"},
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
                for level in ([5, 6, 7] if kind == "cglin" else [hiercc_level])
            }
            matched = [
                r
                for r in pool
                if any(
                    _matching(q, r, kind, lin_level if kind == "cglin" else hiercc_level)
                    for q in relevant
                )
            ]
            requested = len(matched)
            matched = _grouped_analysis_exports(matched, client, scope_work / "profiles", export)
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
        all_rows = rows + context
        explicit = [m.catalogue for m in dataset.profiles] + list(
            _catalogues(Path(catalogues)) if catalogues else []
        )
        catalogue_values = catalogues_for_records(all_rows, explicit)
        audit = {
            "policy": {"lin_level": lin_level, "hiercc_level": hiercc_level},
            "resume_fingerprint": fingerprint,
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
        }
        _save(work / "catalogues.json", [asdict(c) for c in catalogue_values])
        _save(
            work / "profiles.json",
            {
                "records": all_rows,
                "lineages": list(dataset.lineages)
                + list(from_profile_records(context, catalogue_values).lineages)
                if context
                else list(dataset.lineages),
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
                "parameters": dict(dataset.parameters, profile_context=audit),
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
