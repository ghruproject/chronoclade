"""Own a complete profile-stage output without acquiring assemblies or dating."""

from collections import Counter
from dataclasses import dataclass
import os
import hashlib
import json
import sys
from pathlib import Path

from chronoclade import __version__
from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import load_dataset
from chronoclade.profile_analysis import analyse_profiles
from chronoclade.profile_report import write_profile_report
from .report import write_partition_index
from chronoclade.selection_manifest import load_selection_ensemble, write_selection_ensemble
from chronoclade.context_partitions import partition_records, profile_records


@dataclass(frozen=True)
class CGMLSTResult:
    manifest_path: Path
    report_path: Path


def _relative_artifacts(value, root):
    """Keep legacy stage evidence paths portable in the saved analysis adapter."""
    if isinstance(value, dict):
        return {key: _relative_artifacts(item, root) for key, item in value.items()}
    if isinstance(value, list):
        return [_relative_artifacts(item, root) for item in value]
    if isinstance(value, str) and value.startswith(str(root) + os.sep):
        return Path(value).relative_to(root).as_posix()
    return value


def _restore_artifact_paths(value, root):
    """Resolve portable checkpoint references before analysis adapters read them."""
    if isinstance(value, dict):
        return {key: _restore_artifact_paths(item, root) for key, item in value.items()}
    if isinstance(value, list):
        return [_restore_artifact_paths(item, root) for item in value]
    if isinstance(value, str) and len(value) < 4096 and "\n" not in value:
        path = Path(value)
        try:
            if not path.is_absolute() and (root / path).is_file():
                return str(root / path)
        except OSError:
            pass
    return value


def _validate_resumed_ensemble(path, analysis, dataset_manifest, partition_path,
                               block_id, parameters):
    """Bind a portable valid bundle to this exact run and numerical checkpoint."""
    from chronoclade.matrix_distances import BinaryDistanceEvidence
    from chronoclade.selections import read_distance_evidence

    ensemble = load_selection_ensemble(path)
    source = ensemble["source"]
    expected_sources = {
        "dataset": file_sha256(dataset_manifest),
        "partition": file_sha256(partition_path),
    }
    if (ensemble.get("parameters") != parameters
            or source.get("partition_block_id") != block_id
            or any((source.get(key) or {}).get("sha256") != digest
                   for key, digest in expected_sources.items())):
        raise ValueError("cgmlst resume selection bundle differs from saved run parameters or sources")
    for run in ensemble["selections"]:
        selection = json.loads((path.parent / run["path"]).read_text())
        if selection.get("parameters") != parameters:
            raise ValueError("cgmlst resume selection parameters differ from saved run")
    evidence = read_distance_evidence(analysis)
    saved = json.loads((path.parent / source["distance_evidence"]["path"]).read_text())
    if isinstance(evidence, BinaryDistanceEvidence):
        expected = dict(evidence.descriptor)
        # Snapshot filenames may differ, but numerical array content must agree.
        for descriptor in (expected, saved):
            descriptor["paths"] = {
                key: {"sha256": spec["sha256"]}
                for key, spec in descriptor["paths"].items()
            }
    else:
        expected = {"method": evidence.method,
                    "distance_definition": evidence.distance_definition,
                    "records": list(evidence.records)}
    if saved != expected:
        raise ValueError("cgmlst resume selection distance evidence differs from numerical checkpoint")


def add_metadata_summary(analysis, rows):
    """Refresh report metadata coverage from the full, supplied profile pool."""
    import re
    from chronoclade.location_network.reconstruction import is_input_sample, observed_host

    counts = {}
    covered = Counter()
    for row in rows:
        role = "input" if is_input_sample(row) else "context"
        host = observed_host(row.get("host"))
        date = str(row.get("collection_date") or "").strip()
        year_match = re.match(r"^(\d{4})(?:$|[-/])", date)
        year = year_match.group(1) if year_match else None
        for field, value in (("Host", host), ("Collection year", year)):
            key = (field, value or "Unknown")
            counts.setdefault(key, Counter())[role] += 1
        covered["total"] += 1
        covered["host_reported"] += host is not None
        covered["date_reported"] += year is not None
    analysis["metadata_summary"] = {
        "host": {"total": covered["total"], "reported": covered["host_reported"],
                 "missing": covered["total"] - covered["host_reported"]},
        "collection_date": {"total": covered["total"], "reported": covered["date_reported"],
                            "missing": covered["total"] - covered["date_reported"]},
    }
    analysis["metadata_composition"] = [
        {"field": field, "value": value, "input": roles["input"], "context": roles["context"],
         "total": sum(roles.values())}
        for (field, value), roles in sorted(counts.items())
    ]
    return analysis


def _add_state_networks(analysis, rows, directory):
    """Reconstruct host observations on each exact saved country-view topology."""
    from io import StringIO
    from Bio import Phylo
    from chronoclade.location_network.reconstruction import build_state_network

    hosts = []
    for location in analysis.get("location_network", []):
        def reconstruct(model):
            newick = model.get("reconstruction", {}).get("newick")
            if not newick:
                return None
            tree = Phylo.read(StringIO(newick), "newick")
            network = build_state_network(tree, rows, analysis.get("nearest_neighbours", []),
                                          state_field="host", tree_basis="cgmlst")
            for field in ("cohort_id", "id", "label", "scope", "distance_method", "distance_units"):
                if field in model:
                    network[field] = model[field]
            return network
        host = reconstruct(location)
        if host is None:
            continue
        host["views"] = [view for model in location.get("views", [])
                         if (view := reconstruct(model)) is not None]
        path = directory / f"{host['cohort_id']}_host_network.json"
        write_json(path, host)
        host["audit_path"] = str(path)
        analysis["paths"][f"{host['cohort_id']}_host_network_audit"] = str(path)
        hosts.append(host)
    analysis["state_networks"] = {"country": analysis.get("location_network", []), "host": hosts}


def run_cgmlst(
    dataset_manifest,
    output,
    *,
    lin_level=5,
    hiercc_level=None,
    replicates=1,
    context_size=50,
    nearest_per_query=3,
    include=None,
    seed=42,
    bootstrap_replicates=None,
    min_overlap=0.9,
    distance_threshold=None,
    distance_method="alleles",
    allele_catalogue=None,
    representative_reserve=0,
    resume=False,
):
    dataset_manifest, output = Path(dataset_manifest).resolve(), Path(output).resolve()
    if distance_method not in {"alleles", "direct-dna"}:
        raise ValueError("distance_method must be alleles or direct-dna")
    if (distance_method == "direct-dna") != (allele_catalogue is not None):
        raise ValueError("Direct DNA requires --allele-catalogue; allele profiles do not use it")
    if distance_threshold is None:
        if distance_method == "direct-dna":
            raise ValueError("Direct DNA requires an explicit --distance-threshold in nucleotide mismatch fraction")
        distance_threshold = 0.02
    if not 0 <= distance_threshold <= 1:
        raise ValueError("distance_threshold must be a mismatch fraction in [0,1]")
    dataset = load_dataset(dataset_manifest)
    records = profile_records(dataset)
    blocks, audit = partition_records(records, lin_level=lin_level, hiercc_level=hiercc_level)
    if not any(row["role"] == "input" for row in records):
        raise ValueError("cgmlst requires at least one input sample")
    # Immutable namespaces avoid overwriting another stage or earlier completed result.
    if not resume and output.exists() and any(output.iterdir()):
        raise ValueError("cgmlst output must be empty; choose a separate output directory")
    if replicates < 1 or context_size < 0 or nearest_per_query < 1 or representative_reserve < 0:
        raise ValueError(
            "Selections require positive replicates/neighbours and nonnegative context size"
        )
    requested = {}
    selectable = [row for block in blocks for row in block["records"] if row["role"] == "context"]
    for identifier in include or []:
        matches = [
            row
            for row in selectable
            if identifier
            in {
                row["sample_id"],
                str(row.get("source_genome_id") or ""),
                str(row.get("accession") or ""),
            }
        ]
        if len(matches) != 1:
            raise ValueError(f"Pinned context {identifier!r} is not uniquely in a selected block")
        requested[identifier] = matches[0]["sample_id"]
    models = {}
    if distance_method == "direct-dna":
        from chronoclade.cgmlst.dna import DNAComparator
        from chronoclade.profile_analysis import MAX_RECORDS
        for block in blocks:
            if len(block["records"]) > MAX_RECORDS:
                raise ValueError(
                    f"Direct DNA currently supports up to {MAX_RECORDS} records per block; "
                    "no samples have been discarded. Use allele profiles for larger blocks."
                )
            models[block["block_id"]] = DNAComparator(allele_catalogue, block["records"])
    parameters = {
        "distance_method": distance_method,
        "allele_catalogue": str(Path(allele_catalogue).resolve()) if allele_catalogue else None,
        "lin_level": lin_level, "hiercc_level": hiercc_level,
        "replicates": replicates, "context_size": context_size,
        "representative_reserve": representative_reserve,
        "nearest_per_query": nearest_per_query, "include": include or [],
        "seed": seed, "bootstrap_replicates": bootstrap_replicates,
        "min_overlap": min_overlap, "distance_threshold": distance_threshold,
    }
    identity = {
        "schema": "chronoclade.cgmlst.resume", "schema_version": 1,
        "software_version": __version__, "parameters": parameters,
        "dataset_sha256": file_sha256(dataset_manifest),
        "records_sha256": hashlib.sha256(
            json.dumps(records, sort_keys=True, allow_nan=False).encode()
        ).hexdigest(),
        "catalogues": {key: model.provenance for key, model in models.items()},
    }
    checkpoint_path = output / "resume.json"
    # Retrieval may have failed before this stage created any output.
    if resume and (not output.exists() or not any(output.iterdir())):
        resume = False
    if resume:
        if not checkpoint_path.is_file():
            raise ValueError("No native cgmlst resume checkpoint exists in the output directory")
        if json.loads(checkpoint_path.read_text()) != identity:
            raise ValueError("cgmlst resume source or parameters differ from the saved run")
        if (output / "cgmlst.json").exists():
            raise ValueError("cgmlst run is already complete; choose a separate output directory")
    else:
        output.mkdir(parents=True, exist_ok=True)
        write_json(checkpoint_path, identity)
    partition_path = output / "partitions.json"
    partition_evidence = {
        "schema": "chronoclade.cgmlst.partitions",
        "schema_version": 1,
        "dataset_id": dataset.dataset_id,
        "dataset_sha256": file_sha256(dataset_manifest),
        **audit,
        "blocks": [
            {key: value for key, value in block.items() if key != "records"}
            | {"sample_ids": [row["sample_id"] for row in block["records"]]}
            for block in blocks
        ],
    }
    write_json(partition_path, partition_evidence)
    products = []
    for block in blocks:
        directory = output / block["block_id"]
        rows = block["records"]
        queries = [row for row in rows if row["role"] == "input"]
        contexts = [row for row in rows if row["role"] == "context"]
        from chronoclade.profile_analysis import MAX_RECORDS
        large = len(rows) > MAX_RECORDS
        actual_bootstraps = (0 if large else 30) if bootstrap_replicates is None else bootstrap_replicates
        if large:
            from chronoclade.cgmlst.scale import analyse_profiles_scale
            engine = analyse_profiles_scale
        else:
            engine = analyse_profiles
        distance_model = models.get(block["block_id"])
        numerical_path = directory / "numerical_analysis.json"
        if resume and numerical_path.is_file():
            saved = json.loads(numerical_path.read_text())
            for relative, digest in saved["artifacts"].items():
                path = directory / relative
                if not path.is_file() or file_sha256(path) != digest:
                    raise ValueError(f"cgmlst resume numerical artifact changed: {relative}")
            analysis = _restore_artifact_paths(saved["analysis"], directory)
        else:
            analysis = engine(
                rows,
                output=directory,
                seed=seed,
                min_overlap=min_overlap,
                bootstrap_replicates=actual_bootstraps,
                distance_threshold=distance_threshold,
                tree_limit=max(2, len(rows)),
                **({"distance_model": distance_model} if distance_model else {}),
            )
            if distance_model:
                analysis["paths"]["sequence_comparison_audit"] = str(distance_model.save_audit(directory))
            # Publish before report enrichment or selection packaging can fail.
            # Report widgets are regenerated; protect the underlying numerical evidence.
            artifacts = {
                path.relative_to(directory).as_posix(): file_sha256(path)
                for path in directory.rglob("*")
                if path.is_file() and path.suffix in {".csv", ".npy", ".nwk", ".newick", ".json"}
                and path.name != "profile_analysis.json"
            }
            write_json(numerical_path, {
                "analysis": _relative_artifacts(analysis, directory), "artifacts": artifacts,
            })
        if large:
            sys.setrecursionlimit(max(sys.getrecursionlimit(), 4 * len(rows) + 100))
        analysis.setdefault("distance_method", "alleles")
        analysis.setdefault("distance_definition", "fraction of mismatching jointly called cgMLST loci")
        analysis.setdefault("distance_units", "cgMLST mismatch fraction")
        analysis.setdefault("raw_difference_field", "allele_differences")
        analysis.update(
            species=block["species"],
            lineage="ST"
            + block["mlst_st"]
            + " · "
            + block["kind"]
            + " "
            + str(block["level"])
            + " · "
            + ".".join(map(str, block["prefix"])),
            partition={key: value for key, value in block.items() if key != "records"},
        )
        excluded_ids = {row["sample_id"] for row in analysis["exclusions"]}
        analysed_inputs = sum(row["sample_id"] not in excluded_ids for row in queries)
        analysis["status"] = "complete" if analysed_inputs == len(queries) else "partial"
        analysis["coverage"]["analysed_inputs"] = analysed_inputs
        analysis["coverage"]["total_inputs"] = len(queries)
        composition = Counter(
            (
                row["origin"],
                row.get("country") or "Unknown",
                row.get("region") or "Unknown",
                row.get("nuts2") or "Unknown",
            )
            for row in rows
        )
        analysis["metadata_geography"] = [
            dict(origin=key[0], country=key[1], region=key[2], nuts2=key[3], count=count)
            for key, count in sorted(composition.items())
        ]
        add_metadata_summary(analysis, rows)
        from chronoclade.report_components.root_comparison import add_root_comparisons
        add_root_comparisons(analysis, rows, directory)
        from chronoclade.report_components.lineage_groups import render_lineage_groups
        render_lineage_groups(analysis, rows, directory)
        _add_state_networks(analysis, rows, directory)
        block_include = [alias for alias, ident in requested.items()
                         if ident in {row["sample_id"] for row in contexts}]
        selection_parameters = {
            "size": context_size, "nearest_per_query": nearest_per_query,
            "representative_reserve": representative_reserve,
            "include": block_include, "seed": seed, "replicates": replicates,
        }
        ensemble_path = directory / "selection_bundle" / "ensemble.json"
        if resume and ensemble_path.is_file():
            _validate_resumed_ensemble(
                ensemble_path, analysis, dataset_manifest, partition_path,
                block["block_id"], selection_parameters,
            )
        else:
            ensemble_path = write_selection_ensemble(
                queries,
                contexts,
                {**analysis, 'genetic_groups': analysis.get('lineage_groups', analysis['genetic_groups'])},
                dataset_manifest=dataset_manifest,
                output=directory / "selection_bundle",
                replicates=replicates,
                size=context_size,
                representative_reserve=representative_reserve,
                nearest_per_query=nearest_per_query,
                include=block_include,
                seed=seed,
                partition_manifest=partition_path,
            )
        provenance = {
            "input_partition": {
                "total_inputs": sum(row['role'] == 'input' for row in records),
                "block_inputs": len(queries),
                "blocks": len(blocks),
                "assigned_inputs": sum(row['role'] == 'input' for item in blocks for row in item['records']),
            },
            "source_dataset": dataset.dataset_id,
            "distance_method": distance_method,
            "sequence_catalogue": distance_model.provenance if distance_model else None,
            "profile_context": dataset.parameters.get("profile_context"),
            "bounded_public_context_count": len(contexts),
            "context_scope": audit["context_scope"],
            "coverage": {
                label: {
                    "total": len(group),
                    "profiles_available": sum(
                        any(value is not None for value in row.get("cgmlst_profile", {}).values())
                        for row in group
                    ),
                }
                for label, group in (("queries", queries), ("context", contexts))
            },
        }
        analysis['input_partition'] = provenance['input_partition']
        report = write_profile_report(
            analysis,
            directory=directory,
            provenance=provenance,
            report_label="cgMLST profile report",
            group_index_href="../index.html",
        )
        write_json(directory / "profile_analysis.json", _relative_artifacts(analysis, directory))
        products.append(
            {
                "block_id": block["block_id"],
                "label": "ST"
                + block["mlst_st"]
                + " · "
                + block["kind"]
                + " "
                + str(block["level"])
                + " · "
                + ".".join(map(str, block["prefix"])),
                "input_count": len(queries),
                "analysed_input_count": analysed_inputs,
                "analysis_status": analysis["status"],
                "context_count": len(contexts),
                "backend": "scale" if large else "dense",
                "distance_method": distance_method,
                "distance_definition": analysis["distance_definition"],
                "distance_units": analysis["distance_units"],
                "bootstrap_replicates": actual_bootstraps,
                "report": report.relative_to(output).as_posix(),
                "ensemble": ensemble_path.relative_to(output).as_posix(),
                "report_sha256": file_sha256(report),
                "ensemble_sha256": file_sha256(ensemble_path),
            }
        )
    status = (
        "complete"
        if blocks and not audit["unresolved_inputs"] and all(
            product["analysis_status"] == "complete" for product in products)
        else "partial"
        if blocks
        else "unavailable"
    )
    report_path = write_partition_index(
        products,
        audit,
        records_count=len(records),
        lin_level=lin_level,
        hiercc_level=hiercc_level,
        status=status,
        output=output,
    )
    manifest = {
        "schema": "chronoclade.cgmlst.analysis",
        "schema_version": 1,
        "software_version": __version__,
        "status": status,
        "dataset_id": dataset.dataset_id,
        "distance_method": distance_method,
        "dataset": {
            "path": os.path.relpath(dataset_manifest, output),
            "sha256": file_sha256(dataset_manifest),
        },
        "parameters": parameters,
        "partitions": {"path": partition_path.name, "sha256": file_sha256(partition_path)},
        "blocks": products,
        "report": {"path": "index.html", "sha256": file_sha256(report_path)},
    }
    manifest_path = output / "cgmlst.json"
    write_json(manifest_path, manifest)
    return CGMLSTResult(manifest_path, report_path)
