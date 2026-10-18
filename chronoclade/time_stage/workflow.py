"""Date a saved corrected tree without constructing or changing phylogenies."""

from pathlib import Path
import shutil

from Bio import Phylo

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.errors import WorkflowError
from chronoclade.execution_provenance import execution_provenance
from chronoclade.selection_manifest import _resolve
from chronoclade.stage_artifacts import (
    header,
    load_tree_result,
    load_result,
    reference,
    fingerprint,
    tools,
    require_tools,
    prepare_job,
    publish,
    fail_job,
    read_dated_tree,
)
from chronoclade.report_components.stage import write_report


TOPOLOGY_VALIDATION_POLICY = "rooted_nontrivial_internal_leafsets_root_excluded"


def native_time_parameters(tree, allow_unsupported=False):
    """Bind dating to the saved tree seed while retaining its polytomies."""
    seed = tree.get("parameters", {}).get("seed")
    if type(seed) is not int:
        raise WorkflowError("Native dating requires an integer seed from the saved tree parameters")
    return {
        "allow_unsupported": allow_unsupported,
        "dating_seed": seed,
        "keep_polytomies": True,
        "dating_inference_policy": "seeded_saved_topology_keep_polytomies",
        "topology_validation_policy": TOPOLOGY_VALIDATION_POLICY,
    }


def validate_rooted_topology(source, dated):
    """Compare rooted internal leaf sets, ignoring labels, lengths and the root."""
    def clades(tree):
        tip_count = len(tree.get_terminals())
        return {
            leaves
            for node in tree.get_nonterminals()
            if 1 < len(leaves := frozenset(tip.name for tip in node.get_terminals())) < tip_count
        }

    source_clades, dated_clades = clades(source), clades(dated)
    if source_clades != dated_clades:
        raise WorkflowError(
            "Dated tree changed the saved rooted topology: "
            f"{len(source_clades - dated_clades)} saved internal clade(s) lost and "
            f"{len(dated_clades - source_clades)} new internal clade(s) introduced"
        )
    return {
        "method": TOPOLOGY_VALIDATION_POLICY,
        "source_clade_count": len(source_clades),
        "dated_clade_count": len(dated_clades),
        "preserved": True,
    }


def run_time(tree_manifest, *, output, force=False, allow_unsupported=False):
    tree_path = Path(tree_manifest).expanduser().resolve()
    tree = load_tree_result(tree_path)
    parameters = native_time_parameters(tree, allow_unsupported)
    seed = parameters["dating_seed"]
    output = Path(output).expanduser().resolve()
    tool_evidence = tools(["treetime"])
    computed_with = execution_provenance(
        "time", "saved-rooted-tree-dating-v1",
        software_version=header("time")["software"]["version"],
    )
    digest = fingerprint(
        {
            "tree": reference(output, tree_path)["sha256"],
            "computed_with": computed_with,
            "tools": tool_evidence,
            "parameters": parameters,
        }
    )
    if (output / "time.json").is_file() and not force:
        import json

        previous = json.loads((output / "time.json").read_text())
        if previous.get("fingerprint") == digest:
            return load_result(output / "time.json", "time")
    output, job = prepare_job(output, digest, "time", force)
    write_json(
        job / "job.json",
        {
            "stage": "time",
            "status": "running",
            "selection_id": tree["selection_id"],
            "source": {"tree": reference(job, tree_path)},
            "fingerprint": digest,
        },
    )
    try:
        frozen_tree = tree_path.with_name(
            "tree-" + reference(output, tree_path)["sha256"] + ".json"
        )
        if not frozen_tree.exists():
            write_json(frozen_tree, tree)
        load_tree_result(frozen_tree)
        rows = _read_rows(_resolve(tree_path.parent, tree["metadata"], "metadata"))
        assessment = tree["temporal_assessment"]
        result = {
            "fingerprint": digest,
            "source": {"tree": reference(output, frozen_tree)},
            "input_tree_sha256": reference(output, tree_path)["sha256"],
            "computed_with": computed_with,
            "rendered_with": computed_with,
            "selection_id": tree["selection_id"],
            "selected_sample_ids": tree["selected_sample_ids"],
            "temporal_assessment": assessment,
            "coherence_screening": tree.get("coherence_screening"),
            "parameters": parameters,
            "tools": tool_evidence,
            "dating_status": "unsupported",
            "dated_tree": None,
            "topology_preservation": None,
            "dating_override": bool(allow_unsupported and not assessment["supported"]),
            "branch_units": None,
            "clock_confidence": None,
        }
        figures = []
        network = None
        host_network = None
        state_networks = None
        if assessment["supported"] or allow_unsupported:
            if len(tree["dated_sample_ids"]) < 3:
                raise WorkflowError(
                    "Dating requires at least three saved eligible dated samples, even with an override"
                )
            require_tools(["treetime"])
            from chronoclade.native_pipeline import LineageFiles, _write_csv, _run_dated_tree
            from chronoclade.temporal_report import (
                write_timetree_confidence,
                write_timetree_figures,
            )

            files = LineageFiles.in_directory(job)
            files.clock_dir.mkdir(exist_ok=True)
            rooted = _resolve(tree_path.parent, tree["rooted_tree"], "rooted_tree")
            shutil.copyfile(rooted, files.clock_dir / "rerooted.newick")
            dated = set(tree["dated_sample_ids"])
            _write_csv(
                job / "clock_dates.csv",
                ["sample_id", "collection_date"],
                (
                    {
                        "sample_id": row["sample_id"],
                        "collection_date": row.get("collection_date") or "",
                    }
                    for row in rows
                    if row["sample_id"] in dated
                ),
            )
            _run_dated_tree(
                files, tree["complete_alignment_sites"], force,
                seed=seed, keep_polytomies=True,
            )
            tips = [tip.name for tip in Phylo.read(files.time_tree, "nexus").get_terminals()]
            if len(tips) != len(tree["selected_sample_ids"]) or set(tips) != set(
                tree["selected_sample_ids"]
            ):
                raise WorkflowError("Dated tree does not retain exact selected sample identifiers")
            from chronoclade.location_network.reconstruction import (
                build_location_network,
                build_optional_state_networks,
            )
            from chronoclade.location_network.tree import draw_country_tree

            biological = read_dated_tree(files.time_tree)
            result["topology_preservation"] = validate_rooted_topology(
                Phylo.read(rooted, "newick"), biological,
            )
            network = build_location_network(biological, rows, tree_basis="sequence")
            network.update(
                branch_units="years",
                tree_sha256=file_sha256(files.time_tree),
                root_policy="saved corrected tree root retained",
            )
            write_json(job / "location_network.json", network)
            optional_networks = build_optional_state_networks(
                biological,
                rows,
                tree_basis="sequence",
                nearest_neighbours=tree.get("nearest_contexts_by_corrected_snps"),
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
            result["state_networks"] = state_network_descriptors
            host_network = optional_networks.get("host")
            if host_network is not None:
                result["host_state_network"] = state_network_descriptors["host"]
            draw_country_tree(
                network,
                rows,
                job / "dated_country_tree.svg",
                title="Dated biological tree location history",
                distance_label="Years",
                dating_note="TreeTime time-scaled biological tree",
            )
            confidence = write_timetree_confidence(job)
            write_timetree_figures(job, confidence)
            result.update(
                dating_status="dated",
                branch_units="years",
                dated_tree=reference(output, files.time_tree),
                clock_confidence=confidence,
                root_policy="saved tree root retained",
            )
            figures = [
                ("dated_country_tree.svg", "Dated-tree location history; years"),
                (
                    "timetree_with_confidence.svg",
                    "TreeTime node-date uncertainty; 90% max-posterior regions",
                ),
            ]
        write_json(job / "assessment.json", assessment)
        write_json(job / "selected_metadata.json", rows)
        note = "<p>Dating override explicitly requested.</p>" if result["dating_override"] else ""
        result["report"] = reference(
            output,
            write_report(
                job / "report.html",
                title="Time-scaled phylogeny" if result["dated_tree"] else "Dating assessment",
                report_label="TIME TREE ANALYSIS REPORT",
                navigation=(("SNP-tree report", tree_path.parent / tree["report"]["path"]),),
                result=result,
                rows=rows,
                network=network,
                host_network=host_network,
                state_networks=state_networks,
                figures=figures,
                extra=note,
            ),
        )
        return publish(output, job, "time", result)
    except Exception as error:
        fail_job(job, "time", error)
        raise


def _read_rows(path):
    import json

    rows = json.loads(Path(path).read_text())
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise WorkflowError("Selected metadata must be a table of objects")
    return rows
