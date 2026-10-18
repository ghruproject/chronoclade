"""Execute assembly phylogeny from a validated immutable selection."""

from pathlib import Path
import json
from html import escape

from Bio import Phylo

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.coherence import COHERENCE_SCREEN_POLICY_VERSION
from chronoclade.datasets import load_dataset
from chronoclade.errors import WorkflowError
from chronoclade.execution_provenance import execution_provenance
from chronoclade.evidence import read_alignment, pairwise_distances, nearest_contexts
from chronoclade.metadata import Sample, select_reference
from chronoclade.selection_manifest import load_selection_manifest, _resolve
from chronoclade.stage_artifacts import (
    reference,
    header,
    fingerprint,
    tools,
    require_tools,
    prepare_job,
    publish,
    fail_job,
    load_tree_result,
)
from chronoclade.report_components.stage import write_report


def selection_rows(selection_path, selection):
    dataset_path = _resolve(
        selection_path.parent, selection["source"]["dataset"], "source.dataset", external=True
    )
    dataset = load_dataset(dataset_path)
    by_id = {row["sample_id"]: dict(row) for row in dataset.samples}
    rows = [by_id[ident] for ident in selection["selected_sample_ids"]]
    return dataset_path, rows


def local_assembly(row, dataset_path):
    raw = row.get("assembly_reference")
    if not raw or "://" in raw:
        return None
    path = Path(raw).expanduser()
    path = path if path.is_absolute() else dataset_path.parent / path
    return path.resolve() if path.is_file() else None


def acquire_assemblies(rows, dataset_path, job, *, api_key, base_url, threads, cache_dir=None):
    from chronoclade.pathogenwatch_download import download_assemblies, validate_fasta
    from chronoclade.file_snapshots import copy_evidence_snapshot

    directory = job / "assemblies"
    directory.mkdir(exist_ok=True)
    paths, audit, pending = {}, [], []
    for row in rows:
        local = local_assembly(row, dataset_path)
        if local:
            original = local.read_bytes()
            content = validate_fasta(original)
            destination = directory / (fingerprint(row["sample_id"]) + ".fasta")
            if destination.exists() and destination.read_bytes() == content:
                copy_method = "existing_verified"
            else:
                destination.unlink(missing_ok=True)
                if content == original:
                    copy_method = copy_evidence_snapshot(local, destination)
                else:
                    destination.write_bytes(content)
                    copy_method = "decompressed_copy"
            paths[row["sample_id"]] = destination
            audit.append(
                {
                    "sample_id": row["sample_id"],
                    "status": "ready",
                    "method": "linked_local_assembly",
                    "copy_method": copy_method,
                    "sha256": file_sha256(destination),
                }
            )
        elif row.get("source_genome_id"):
            pending.append(row)
        else:
            audit.append(
                {
                    "sample_id": row["sample_id"],
                    "status": "failed",
                    "reason": "No linked local assembly or exact Pathogenwatch source_genome_id",
                }
            )
    if pending:
        if not api_key:
            audit.extend(
                {
                    "sample_id": row["sample_id"],
                    "status": "failed",
                    "reason": "Missing Pathogenwatch API key for selected missing assembly",
                }
                for row in pending
            )
        else:
            downloaded, download_audit = download_assemblies(
                pending,
                output=directory,
                cache_dir=Path(cache_dir) if cache_dir is not None else job / "assembly_cache",
                api_key=api_key,
                base_url=base_url,
                workers=threads,
                require_public_identity=True,
            )
            by_source = {entry["source_genome_id"]: entry for entry in download_audit}
            for row in pending:
                source = row["source_genome_id"]
                audit.append({**by_source[source], "sample_id": row["sample_id"]})
                if source in downloaded:
                    paths[row["sample_id"]] = downloaded[source]
    write_json(job / "assembly_audit.json", audit)
    if set(paths) != {row["sample_id"] for row in rows}:
        raise WorkflowError(
            f"Selected assembly acquisition failed; no replacements chosen. See {job / 'assembly_audit.json'}"
        )
    return paths, audit


def run_tree(
    selection_manifest,
    *,
    output,
    threads=1,
    randomisations=100,
    randomisation_jobs=1,
    temporal_p_value=0.05,
    seed=42,
    assess_temporal=True,
    api_key="",
    base_url="https://pathogen.watch",
    force=False,
    assembly_cache=None,
    coherence_policy="strict",
    support_replicates=0,
):
    from chronoclade.tree_stage.support import validate_replicates
    validate_replicates(support_replicates)
    if threads < 1 or randomisation_jobs < 1 or randomisations < 0 or not 0 < temporal_p_value <= 1:
        raise WorkflowError("Invalid tree resource or temporal parameters")
    if coherence_policy not in {"strict", "diagnostic"}:
        raise WorkflowError("Coherence policy must be strict or diagnostic")
    selection_path = Path(selection_manifest).expanduser().resolve()
    selection = load_selection_manifest(selection_path)
    dataset_path, rows = selection_rows(selection_path, selection)
    from chronoclade.progress import progress
    progress("tree", f"{len(rows)} exact selected genomes; output: {output}")
    if len(rows) < 3:
        raise WorkflowError("Tree needs at least three exact selected genomes")
    if not any(row["role"] == "input" for row in rows):
        raise WorkflowError("Native tree requires input genomes for its common reference")
    if len({row["species"] for row in rows}) != 1:
        raise WorkflowError("Tree selection spans incompatible species")
    tool_evidence = tools(
        ["ska", "iqtree", "ClonalFrameML"] + (["treetime"] if assess_temporal else [])
    )
    parameters = dict(
        threads=threads,
        randomisations=randomisations,
        randomisation_jobs=randomisation_jobs,
        temporal_p_value=temporal_p_value,
        seed=seed,
        phylogeny_seed=seed,
        phylogeny_seed_policy="explicit_IQ_TREE",
        reference_selection_policy="highest_N50_common_input",
        assess_temporal=assess_temporal,
        base_url=base_url,
        coherence_policy=coherence_policy,
        coherence_screen_policy_version=COHERENCE_SCREEN_POLICY_VERSION,
        support_replicates=support_replicates,
    )
    local_hashes = {
        row["sample_id"]: file_sha256(local_assembly(row, dataset_path))
        for row in rows
        if local_assembly(row, dataset_path)
    }
    computed_with = execution_provenance(
        "tree", "recombination-corrected-phylogeny-v2",
        software_version=header("tree")["software"]["version"],
    )
    digest = fingerprint(
        {
            "selection": file_sha256(selection_path),
            "computed_with": computed_with,
            "parameters": parameters,
            "tools": tool_evidence,
            "assemblies": local_hashes,
        }
    )
    output = Path(output).expanduser().resolve()
    saved = output / "tree.json"
    if saved.is_file() and not force:
        previous = json.loads(saved.read_text())
        if previous.get("fingerprint") == digest:
            return load_tree_result(saved)
    require_tools(["ska", "iqtree", "ClonalFrameML"])
    output, job = prepare_job(output, digest, "tree", force)
    try:
        from chronoclade.native_pipeline import (
            LineageFiles,
            _write_lineage_inputs,
            _run_core_phylogeny,
            _run_observed_clock,
            _dated_members,
            complete_alignment_sites,
            alignment_length,
            _write_csv,
        )
        from chronoclade.location_network.reconstruction import (
            build_location_network,
            build_optional_state_networks,
        )
        from chronoclade.location_network.tree import draw_country_tree
        from chronoclade.tree_stage.temporal import run_clustered_assessment

        paths, audit = acquire_assemblies(
            rows, dataset_path, job, api_key=api_key, base_url=base_url, threads=threads,
            cache_dir=assembly_cache,
        )
        members = [
            Sample(
                row["sample_id"],
                paths[row["sample_id"]],
                row.get("collection_date") or "",
                row.get("country") or "Unknown",
                row["species"],
                selection["selection_id"],
                "local" if row["role"] == "input" else "context",
                country=row.get("country") or "",
                region=row.get("region") or "",
                city=row.get("city") or "",
            )
            for row in rows
        ]
        files = LineageFiles.in_directory(job)
        _write_lineage_inputs(job, members)
        input_members = [sample for sample in members if sample.origin == "local"]
        selected_reference = select_reference(input_members)
        tree_path, alignment, _ = _run_core_phylogeny(
            files, members, selected_reference, threads, force, coherence_policy=coherence_policy, seed=seed
        )
        coherence = {"policy": coherence_policy, "screen_policy_version": COHERENCE_SCREEN_POLICY_VERSION, "raw_action": "review_and_retain", "corrected_action": "abort_on_flags" if coherence_policy == "strict" else "review_and_retain"}
        for label, audit_path in (("raw", files.coherence_screen), ("corrected", files.clonal_coherence_screen)):
            coherence[label] = json.loads(audit_path.with_suffix(".json").read_text())
        flagged = sorted(set(coherence["raw"]["flagged_samples"]) | set(coherence["corrected"]["flagged_samples"]))
        labels = {row["sample_id"]: row.get("label") or row.get("accession") or row["sample_id"] for row in rows}
        if flagged:
            action = (
                "Raw-distance review retained all selected genomes; the recombination-filtered screen passed. "
                if coherence_policy == "strict"
                else "Diagnostic policy retained all selected genomes despite relative-distance flags. "
            )
            coherence["warning"] = action + "These can reflect legitimate cohort diversity and do not establish incorrect genome identity. Review raw and corrected statistics. Flagged genomes: " + ", ".join(labels[ident] for ident in flagged)
        else:
            coherence["warning"] = None
        ids = selection["selected_sample_ids"]
        for actual, name in (
            (
                [tip.name for tip in Phylo.read(tree_path, "newick").get_terminals()],
                "corrected tree",
            ),
            (
                [tip.name for tip in Phylo.read(files.starting_tree, "newick").get_terminals()],
                "raw tree",
            ),
            (list(read_alignment(alignment)), "alignment"),
        ):
            if len(actual) != len(ids) or set(actual) != set(ids):
                raise WorkflowError(f"{name} does not preserve exact selected IDs")
        from chronoclade.recombination_report import write_recombination_evidence

        recombination = write_recombination_evidence(
            alignment=files.alignment,
            filtered_alignment=alignment,
            importations=files.importations,
            output_directory=job,
            reference=selected_reference.assembly,
        )
        from chronoclade.tree_stage.support import run_branch_support, support_html
        branch_support = run_branch_support(
            alignment=alignment, corrected_tree=tree_path, directory=job,
            output_root=output, replicates=support_replicates, threads=threads,
            seed=seed, force=force,
        )
        sites = complete_alignment_sites(alignment)
        dated = _dated_members(members)
        assessment = {
            "code": "not_assessed",
            "supported": False,
            "reason": "Temporal assessment disabled"
            if not assess_temporal
            else "Fewer than three usable distinct collection dates",
        }
        rooted = tree_path
        if (
            assess_temporal
            and len(dated) >= 3
            and len({sample.collection_date for sample in dated}) >= 3
        ):
            require_tools(["treetime"])
            _write_csv(
                job / "clock_dates.csv",
                ["sample_id", "collection_date"],
                ({"sample_id": s.sample_id, "collection_date": s.collection_date} for s in dated),
            )
            _run_observed_clock(files, tree_path, sites, force)
            assessment = run_clustered_assessment(
                tree=tree_path,
                sequence_length=sites,
                samples=dated,
                observed_clock=files.clock_file,
                randomisations=randomisations,
                randomisation_jobs=min(randomisation_jobs, threads),
                seed=seed,
                output=job / "assessment.json",
                p_value_threshold=temporal_p_value,
            )
            if (files.clock_dir / "rerooted.newick").is_file():
                rooted = files.clock_dir / "rerooted.newick"
            from chronoclade.temporal_report import (
                write_randomisation_plot,
                write_randomisation_csv,
            )

            write_randomisation_plot(assessment["temporal"], job / "date_randomisation.svg")
            write_randomisation_csv(assessment["temporal"], job / "date_randomisation.csv")
        write_json(job / "assessment.json", assessment)
        write_json(job / "selected_metadata.json", rows)
        biological_tree = Phylo.read(rooted, "newick")
        network = build_location_network(biological_tree, rows, tree_basis="sequence")
        network.update(
            branch_units="substitutions/site",
            tree_sha256=file_sha256(rooted),
            root_policy="least-squares clock fit"
            if rooted != tree_path
            else "native corrected-tree root; no temporal rooting",
        )
        write_json(job / "location_network.json", network)
        draw_country_tree(
            network,
            rows,
            job / "corrected_tree.svg",
            title="Recombination-adjusted selected-sample tree",
            distance_label="Substitutions/site",
            dating_note="Collection dates are metadata; this tree is not time scaled.",
        )
        distances = pairwise_distances(read_alignment(alignment), members)
        write_json(job / "snp_distances.json", distances)
        nearest = nearest_contexts(members, distances)
        write_json(job / "snp_nearest_contexts.json", nearest)
        optional_networks = build_optional_state_networks(
            biological_tree,
            rows,
            tree_basis="sequence",
            nearest_neighbours=nearest,
            provenance={
                "branch_units": network["branch_units"],
                "tree_sha256": network["tree_sha256"],
                "root_policy": network["root_policy"],
            },
        )
        state_networks = {"country": network}
        state_network_descriptors = {
            "country": reference(output, job / "location_network.json")
        }
        for field, state_network in optional_networks.items():
            path = job / f"{field}_state_network.json"
            write_json(path, state_network)
            state_networks[field] = state_network
            state_network_descriptors[field] = reference(output, path)
        host_network = optional_networks.get("host")
        result = {
            "fingerprint": digest,
            "computed_with": computed_with, "rendered_with": computed_with,
            "source": {
                "selection": reference(output, selection_path),
                "dataset": reference(output, dataset_path),
            },
            "selection_id": selection["selection_id"],
            "selected_sample_ids": ids,
            "parameters": parameters,
            "coherence_screening": coherence,
            "tools": tool_evidence,
            "reference_sample_id": selected_reference.sample_id,
            "recombination": recombination,
            "branch_support": branch_support,
            "nearest_contexts_by_corrected_snps": nearest,
            "reference": reference(output, selected_reference.assembly),
            "assembly_audit": audit,
            "raw_tree": reference(output, files.starting_tree),
            "corrected_tree": reference(output, tree_path),
            "rooted_tree": reference(output, rooted),
            "alignment": reference(output, alignment),
            "metadata": reference(output, job / "selected_metadata.json"),
            "assessment": reference(output, job / "assessment.json"),
            "alignment_columns": alignment_length(alignment),
            "complete_alignment_sites": sites,
            "branch_units": "substitutions/site",
            "root_policy": network["root_policy"],
            "state_networks": state_network_descriptors,
            "temporal_assessment": assessment,
            "dated_sample_ids": [s.sample_id for s in dated],
            "undated_sample_ids": [s.sample_id for s in members if s not in dated],
            "qc_exclusions": [],
            "report": reference(
                output,
                write_report(
                    job / "report.html",
                    title="Assembly phylogeny",
                    navigation=(("cgMLST report", selection_path.parent.parent.parent / "profile_report.html"),),
                    result={
                        "selected_sample_ids": ids,
                        "temporal_assessment": assessment,
                        "branch_units": "substitutions/site",
                        "nearest_contexts_by_corrected_snps": nearest,
                        "coherence_screening": coherence,
                    },
                    rows=rows,
                    extra=support_html(branch_support) + "<details><summary>Raw and corrected coherence audits</summary><pre>" + escape(json.dumps(coherence, indent=2)) + "</pre></details>",
                    network=network,
                    host_network=host_network,
                    state_networks=state_networks,
                    figures=[
                        (
                            "recombination_map.svg",
                            "ClonalFrameML inferred importation and filtering profile",
                        ),
                        ("corrected_tree.svg", "Corrected biological tree; substitutions/site"),
                        ("clock/root_to_tip_regression.svg", "Root-to-tip temporal diagnostic"),
                        ("date_randomisation.svg", "Clustered date-permutation screen"),
                    ],
                ),
            ),
        }
        if host_network is not None:
            result["host_state_network"] = state_network_descriptors["host"]
        return publish(output, job, "tree", result)
    except Exception as error:
        fail_job(job, "tree", error)
        raise
