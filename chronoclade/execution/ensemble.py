"""Run checked selections independently and retain every requested outcome."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock
import os
import tempfile

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import _read, _resolve, load_selection_ensemble
from chronoclade.stage_artifacts import (
    fingerprint, header, load_result, load_tree_result, reference,
)


def _freeze(path):
    """Snapshot bytes beside their source so relative artifact references stay valid."""
    path = Path(path).resolve()
    digest = file_sha256(path)
    frozen = path.with_name(path.stem + "-" + digest + path.suffix)
    if frozen.exists():
        if file_sha256(frozen) != digest:
            raise WorkflowError("Existing immutable ensemble source was modified")
        return frozen
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(path.read_bytes())
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, frozen)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return frozen


def load_tree_ensemble(path):
    """Validate the index, exact selection bindings and completed tree results."""
    path = Path(path).expanduser().resolve()
    index = _read(path)
    if index.get("schema") != "chronoclade.tree-ensemble" or index.get("schema_version") != 1:
        raise WorkflowError("Expected a versioned tree ensemble")
    ensemble_path = _resolve(path.parent, index["source"]["ensemble"], "ensemble", external=True)
    ensemble = load_selection_ensemble(ensemble_path)
    expected = {entry["selection_id"]: entry for entry in ensemble["selections"]}
    runs = index.get("runs", [])
    if [run.get("selection_id") for run in runs] != list(expected):
        raise WorkflowError("Tree ensemble must retain every selection exactly once in order")
    for run in runs:
        if run.get("status") == "complete":
            result_path = _resolve(path.parent, run["result"], "tree result", external=True)
            result = load_tree_result(result_path)
            ident = run["selection_id"]
            if result["selection_id"] != ident or result["source"]["selection"]["sha256"] != expected[ident]["sha256"]:
                raise WorkflowError("Ensemble tree does not match the exact saved selection")
        elif run.get("status") not in {"pending", "running", "failed", "unavailable"}:
            raise WorkflowError("Unknown tree ensemble run status")
    return index, ensemble_path


def _initial(output, source, stage, ids, parameters, force):
    output = Path(output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = output / f"{stage}_ensemble.json"
    source_ref = reference(output, source)
    digest = fingerprint({"source": source_ref["sha256"], "parameters": parameters})
    if manifest.exists() and not force and _read(manifest).get("fingerprint") != digest:
        raise WorkflowError("Ensemble inputs or parameters changed; use a fresh output or --force")
    result = {
        **header(f"{stage}-ensemble"), "status": "running", "fingerprint": digest,
        "source": {"ensemble" if stage == "tree" else "tree_ensemble": source_ref},
        "parameters": parameters,
        "runs": [{"selection_id": ident, "status": "pending"} for ident in ids],
    }
    return output, manifest, result


def _execute(manifest, result, jobs, worker):
    """Publish progress atomically; one failing selection never removes others."""
    write_json(manifest, result)
    lock = Lock()

    def started(i, row):
        with lock:
            result["runs"][i] = {**row, "status": "running"}
            write_json(manifest, result)
        try:
            return worker(i, row)
        except Exception as error:
            return {**row, "status": "failed", "reason": str(error)}

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        tasks = {pool.submit(started, i, dict(row)): i for i, row in enumerate(result["runs"])}
        for future in as_completed(tasks):
            with lock:
                result["runs"][tasks[future]] = future.result()
                write_json(manifest, result)
    result["status"] = "complete"
    result["completed_runs"] = sum(row["status"] == "complete" for row in result["runs"])
    result["failed_runs"] = sum(row["status"] == "failed" for row in result["runs"])
    write_json(manifest, result)
    from chronoclade.stage_readme import write_stage_readme
    stage = 'tree' if manifest.name == 'tree_ensemble.json' else 'time'
    write_stage_readme(manifest.parent, stage, ensemble=True,
                       completed_runs=result['completed_runs'])
    return result


def run_tree_ensemble(ensemble_manifest, *, output, ensemble_jobs=1, force=False, **options):
    from chronoclade.tree_stage.workflow import run_tree

    if ensemble_jobs < 1:
        raise WorkflowError("Ensemble jobs must be positive")
    source = Path(ensemble_manifest).expanduser().resolve()
    ensemble = load_selection_ensemble(source)
    options.setdefault("coherence_policy", "strict")
    safe_parameters = {key: str(value) if isinstance(value, Path) else value
                       for key, value in options.items() if key != "api_key"}
    safe_parameters["ensemble_jobs"] = ensemble_jobs
    safe_parameters["phylogeny_seed"] = options.get("seed", 42)
    safe_parameters["phylogeny_seed_policy"] = "explicit_IQ_TREE"
    safe_parameters["reference_selection_policy"] = "highest_N50_common_input"
    output, manifest, result = _initial(
        output, source, "tree", [s["selection_id"] for s in ensemble["selections"]], safe_parameters, force,
    )
    # Namespace transport caches by provider; child scientific outputs remain independent.
    options.setdefault("assembly_cache", output / "assembly_cache" / fingerprint(options.get("base_url", "https://pathogen.watch")))

    def worker(i, row):
        selection = ensemble["selections"][i]
        path = _resolve(source.parent, selection, "selection")
        destination = output / "runs" / f"selection-{i + 1:03d}"
        try:
            run_tree(path, output=destination, force=force, **options)
            return {**row, "status": "complete", "result": reference(output, _freeze(destination / "tree.json"))}
        except Exception as error:
            reason = str(error)
            if options.get("api_key"):
                reason = reason.replace(options["api_key"], "[redacted]")
            return {**row, "status": "failed", "reason": reason}

    return _execute(manifest, result, ensemble_jobs, worker)


def register_tree_ensemble(ensemble_manifest, tree_manifests, *, output, force=False):
    """Bind already-computed checked trees without copying or rerunning them."""
    source = Path(ensemble_manifest).expanduser().resolve()
    ensemble = load_selection_ensemble(source)
    results = {}
    for path in map(Path, tree_manifests):
        tree = load_tree_result(path)
        ident = tree["selection_id"]
        if ident in results:
            raise WorkflowError("Duplicate tree for one selection")
        results[ident] = (path.resolve(), tree)
    expected = {entry["selection_id"]: entry for entry in ensemble["selections"]}
    if set(results) - set(expected):
        raise WorkflowError("Tree does not belong to this ensemble")
    registered_hashes = {ident: file_sha256(results[ident][0]) for ident in expected if ident in results}
    output, manifest, index = _initial(output, source, "tree", list(expected),
                                       {"operation": "register_checked_trees", "tree_hashes": registered_hashes}, force)
    for row in index["runs"]:
        ident = row["selection_id"]
        if ident not in results:
            row.update(status="unavailable", reason="No checked tree supplied")
            continue
        path, tree = results[ident]
        if tree["source"]["selection"]["sha256"] != expected[ident]["sha256"]:
            raise WorkflowError("Tree selection checksum differs from declared ensemble")
        row.update(status="complete", result=reference(output, _freeze(path)))
    index["status"] = "complete"
    index["completed_runs"] = len(results)
    index["failed_runs"] = 0
    write_json(manifest, index)
    return index


def run_time_ensemble(tree_ensemble, *, output, ensemble_jobs=1, force=False,
                      allow_unsupported=False, cluster_snp_cutoff=10, reuse_results=(),
                      collect_only=False, missing_outcomes=None):
    from chronoclade.time_stage.workflow import run_time, native_time_parameters
    from chronoclade.time_stage.comparison import compare_time_runs

    if missing_outcomes and not collect_only:
        raise WorkflowError("Missing task outcomes are only valid for collect-only execution")
    if ensemble_jobs < 1 or cluster_snp_cutoff < 0:
        raise WorkflowError("Invalid ensemble resource or clustering parameter")
    source = Path(tree_ensemble).expanduser().resolve()
    trees, selection_source = load_tree_ensemble(source)
    source = _freeze(source)
    reused = {}
    expected = {row["selection_id"]: row for row in trees["runs"]}
    if missing_outcomes:
        if (set(missing_outcomes) - set(expected)
                or any(not isinstance(value, dict) or set(value) - {"reason", "scheduler"}
                       for value in missing_outcomes.values())):
            raise WorkflowError("Missing task outcomes may only attach reason and scheduler metadata to known selections")
    dating_parameters = {
        ident: native_time_parameters(
            load_tree_result(_resolve(source.parent, row["result"], "tree result", external=True)),
            allow_unsupported,
        )
        for ident, row in expected.items() if row["status"] == "complete"
    }
    if force and reuse_results and not collect_only:
        raise WorkflowError("Choose checked result reuse or forced execution, not both")
    for path in map(Path, reuse_results):
        checked = load_result(path, "time")
        ident = checked["selection_id"]
        if ident in reused or ident not in expected:
            raise WorkflowError("Reused time results must belong uniquely to this ensemble")
        tree = expected[ident]
        if (tree["status"] != "complete"
                or checked["source"]["tree"]["sha256"] != tree["result"]["sha256"]
                or checked["parameters"]["allow_unsupported"] != allow_unsupported):
            raise WorkflowError("Reused dating result differs from the checked tree or gate parameters")
        if checked.get("parameters") != dating_parameters[ident]:
            raise WorkflowError(
                "Reused dating result differs from the native dating inference policy or seed; "
                f"selection {ident}: expected {dating_parameters[ident]!r}, "
                f"found {checked.get('parameters')!r}. "
                "Use results from the same inference method and seed, or rerun only this time task; "
                "the collector will not rerun inference."
            )
        reused[ident] = (_freeze(path), checked)
    output, manifest, result = _initial(
        output, source, "time", [r["selection_id"] for r in trees["runs"]],
        {"ensemble_jobs": ensemble_jobs, "allow_unsupported": allow_unsupported,
         "cluster_snp_cutoff": cluster_snp_cutoff,
         "dating_parameters_by_selection": dating_parameters,
         "collect_only": collect_only,
         "reused_time_hashes": {ident: file_sha256(path) for ident, (path, _) in reused.items()}}, force,
    )

    def worker(i, row):
        if row["selection_id"] in reused:
            path, checked = reused[row["selection_id"]]
            return {**row, "status": "complete", "dating_status": checked["dating_status"],
                    "execution": "reused_checksum_verified_result", "result": reference(output, path),
                    "reason": checked["temporal_assessment"]["reason"]}
        tree = trees["runs"][i]
        if tree["status"] != "complete":
            return {**row, "status": tree["status"], "reason": tree.get("reason", "Tree unavailable")}
        if collect_only:
            return {**row, "status": "failed", "reason": "No checked time result from array task; inspect task logs and retry",
                    **(missing_outcomes or {}).get(row["selection_id"], {})}
        path = _resolve(source.parent, tree["result"], "tree result", external=True)
        destination = output / "runs" / f"selection-{i + 1:03d}"
        try:
            dated = run_time(path, output=destination, force=force, allow_unsupported=allow_unsupported)
            checked = load_result(destination / "time.json", "time")
            return {**row, "status": "complete", "dating_status": checked["dating_status"],
                    "result": reference(output, destination / "time.json"),
                    "reason": dated["temporal_assessment"]["reason"]}
        except Exception as error:
            return {**row, "status": "failed", "reason": str(error)}

    result = _execute(manifest, result, ensemble_jobs, worker)
    paths = [_resolve(output, row["result"], "time result", external=True) for row in result["runs"] if row["status"] == "complete"]
    # A revised comparison gets its own namespace; published prior comparisons stay immutable.
    comparison_output = output / "comparisons" / fingerprint(result)
    if not (comparison_output / "comparison.json").exists():
        compare_time_runs(selection_source, paths, output=comparison_output,
                          cluster_snp_cutoff=cluster_snp_cutoff,
                          missing_runs={r["selection_id"]: r for r in result["runs"] if r["status"] != "complete"})
    result["comparison"] = reference(output, comparison_output / "comparison.json")
    result["report"] = reference(output, comparison_output / "report.html")
    write_json(manifest, result)
    from chronoclade.stage_readme import write_stage_readme
    write_stage_readme(output, 'time', ensemble=True, report=result['report']['path'])
    return result
