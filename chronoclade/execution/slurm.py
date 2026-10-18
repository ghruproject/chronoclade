"""Bounded SLURM arrays with checked, independent stage outputs and a reducer.

Generate on the shared filesystem in an activated environment. No scheduler calls
are made until the user runs submit.sh. Workers never edit the shared index.
"""

from __future__ import annotations

import os
from pathlib import Path
import shlex
import sys
from datetime import datetime, timezone
from time import perf_counter

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.errors import WorkflowError
from chronoclade.execution.ensemble import (
    _freeze, load_tree_ensemble, register_tree_ensemble, run_time_ensemble,
)
from chronoclade.selection_manifest import _read, _resolve, load_selection_ensemble
from chronoclade.stage_artifacts import reference


def _members(stage, source):
    if stage == "tree":
        ensemble = load_selection_ensemble(source)
        return [(s["selection_id"], _resolve(source.parent, s, "selection"))
                for s in ensemble["selections"]]
    if stage == "time":
        trees, _ = load_tree_ensemble(source)
        return [(s["selection_id"], _resolve(source.parent, s["result"], "tree", external=True))
                for s in trees["runs"] if s["status"] == "complete"]
    raise ValueError("SLURM ensemble stage must be tree or time")


def _script(path, directives, command, working_directory):
    # Environment binaries include native tools as well as Python. Avoid inherited
    # BLAS/OpenMP pools multiplying each native/temporal worker's CPU allocation.
    path.write_text(
        "#!/usr/bin/env bash\n" + "\n".join(directives) + "\nset -euo pipefail\n"
        + "export PATH=" + shlex.quote(str(Path(sys.executable).parent)) + ':"$PATH"\n'
        + "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\n"
        + "cd -- " + shlex.quote(str(working_directory)) + "\n" + command + "\n"
    )
    path.chmod(0o700)


def write_slurm_ensemble(stage, source, *, output, results, resources,
                         max_parallel=4, options=None):
    """Prepare an array and an afterany reducer; return the saved plan path."""
    if type(max_parallel) is not int or max_parallel < 1:
        raise ValueError("Maximum parallel tasks must be a positive integer")
    source = Path(source).expanduser().resolve()
    members = _members(stage, source)
    if not members:
        raise ValueError("No runnable selections in ensemble")
    if stage == "tree":
        cohorts = {tuple(_read(path)["selected_sample_ids"]) for _, path in members}
        if len(cohorts) < len(members):
            from chronoclade.progress import progress
            progress("slurm", f"{len(members)} requested selections contain {len(cohorts)} distinct ordered cohorts. "
                     "Repeated cohorts still run independently; agreement is not independent selection evidence.")
    if stage == "tree" and any(len(_read(path)["selected_sample_ids"]) < 3 for _, path in members):
        raise ValueError("Each tree selection needs at least three genomes")
    options = dict(options or {})
    allowed = ({"threads", "randomisations", "randomisation_jobs", "support_replicates", "seed",
                "assess_temporal", "coherence_policy", "temporal_p_value"} if stage == "tree"
               else {"allow_unsupported", "cluster_snp_cutoff"})
    if set(options) - allowed:
        raise ValueError("Unsupported array options: " + ", ".join(sorted(set(options) - allowed)))
    if stage == "tree":
        from chronoclade.tree_stage.support import validate_replicates
        validate_replicates(options.get("support_replicates", 0))
        if options.get("randomisations", 100) < 0:
            raise ValueError("Randomisations must be nonnegative")
        if options.get("coherence_policy", "strict") not in {"strict", "diagnostic"}:
            raise ValueError("Coherence policy must be strict or diagnostic")
        if not 0 < options.get("temporal_p_value", 0.05) <= 1:
            raise ValueError("Temporal p-value must be greater than zero and at most one")
        options.setdefault("threads", resources.cpus)
        options.setdefault("randomisation_jobs", resources.cpus)
        if any(type(options[k]) is not int or not 1 <= options[k] <= resources.cpus
               for k in ("threads", "randomisation_jobs")):
            raise ValueError("Threads and randomisation jobs must fit cpus per array task")
    output = Path(output).expanduser().resolve()
    results = Path(results).expanduser().resolve()
    if output == results or output.is_relative_to(results) or results.is_relative_to(output):
        raise ValueError("Job and result directories must be separate, non-nested paths")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Job output must be empty")
    if results.exists() and any(results.iterdir()):
        raise ValueError("Use a fresh results directory; retry the saved array to resume")
    output.mkdir(parents=True, exist_ok=True)
    (output / "logs").mkdir()
    tasks = [{"selection_id": ident, "input": reference(output, path),
              "output": str(results / "runs" / f"selection-{i + 1:03d}")}
             for i, (ident, path) in enumerate(members)]
    from dataclasses import asdict
    plan = output / "array.json"
    write_json(plan, {
        "schema": "chronoclade.slurm-ensemble", "schema_version": 1, "stage": stage,
        "source": reference(output, source), "results": str(results), "tasks": tasks,
        "options": options, "resources": asdict(resources), "max_parallel": max_parallel,
        "working_directory": str(Path.cwd()),
    })
    command = shlex.join([sys.executable, "-m", __name__, str(plan)])
    directives = [f"#SBATCH --cpus-per-task={resources.cpus}",
                  f"#SBATCH --mem={resources.memory_gb}G", f"#SBATCH --time={resources.walltime}"]
    for key in ("account", "partition"):
        if getattr(resources, key):
            directives.append(f"#SBATCH --{key}={getattr(resources, key)}")
    if resources.gpus:
        directives.append(f"#SBATCH --gpus={resources.gpus}")
    _script(output / "array.slurm", directives + [
        f"#SBATCH --array=0-{len(tasks) - 1}%{max_parallel}",
        "#SBATCH --output=" + shlex.quote(str(output / "logs/task-%A_%a.log")),
    ], "exec " + command + ' task "${SLURM_ARRAY_TASK_ID:?}"', Path.cwd())
    _script(output / "collect.slurm", directives + [
        "#SBATCH --output=" + shlex.quote(str(output / "logs/collect-%j.log")),
    ], "exec " + command + " collect", Path.cwd())
    submission_log = shlex.quote(str(output / "submissions.tsv"))
    _script(output / "submit.sh", [],
            "array=$(sbatch --parsable " + shlex.quote(str(output / "array.slurm")) + ")\n"
            + 'array=${array%%;*}\n[[ "$array" =~ ^[0-9]+$ ]]\n'
            + f'printf "array\\t%s\\t{file_sha256(plan)}\\n" "$array" >> ' + submission_log + "\n"
            + 'printf "Array job: %s\\n" "$array"\n'
            + 'collector=$(sbatch --parsable --dependency="afterany:$array" '
            + shlex.quote(str(output / "collect.slurm")) + ")\n"
            + 'collector=${collector%%;*}\n[[ "$collector" =~ ^[0-9]+$ ]]\n'
            + 'printf "collector\\t%s\\n" "$collector" >> ' + submission_log + "\n"
            + 'printf "Collector job: %s\\n" "$collector"', Path.cwd())
    return plan


def _load_plan(path):
    path = Path(path).resolve()
    plan = _read(path)
    if plan.get("schema") != "chronoclade.slurm-ensemble" or plan.get("schema_version") != 1:
        raise WorkflowError("Expected a versioned SLURM ensemble plan")
    source = _resolve(path.parent, plan["source"], "ensemble", external=True)
    members = _members(plan["stage"], source)
    if len(members) != len(plan["tasks"]):
        raise WorkflowError("Array task count differs from checked ensemble")
    for task, (ident, source_path) in zip(plan["tasks"], members, strict=True):
        checked = _resolve(path.parent, task["input"], "task input", external=True)
        if task["selection_id"] != ident or checked != source_path:
            raise WorkflowError("Array task differs from checked ensemble selection")
    return path, plan, source


def run_task(path, index):
    """One task owns one output; credentials are loaded only at execution time."""
    path, plan, _ = _load_plan(path)
    if not 0 <= index < len(plan["tasks"]):
        raise WorkflowError("Array task index out of range")
    (path.parent / f"task-{index}.json").unlink(missing_ok=True)
    task = plan["tasks"][index]
    started = perf_counter()
    status_path = path.parent / f"task-{index}.status.json"
    status = {"plan_sha256": file_sha256(path), "selection_id": task["selection_id"],
              "status": "running", "started_at": datetime.now(timezone.utc).isoformat(),
              "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
              "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
              "task_index": index, "output": task["output"]}
    write_json(status_path, status)
    from chronoclade.progress import progress_task, progress
    options = dict(plan["options"])
    try:
        allocation = os.environ.get("SLURM_CPUS_PER_TASK")
        if allocation and int(allocation) < plan["resources"]["cpus"]:
            raise WorkflowError("SLURM CPU allocation is smaller than the prepared task")
        if plan["stage"] == "tree":
            from chronoclade.tree_stage.workflow import run_tree as runner
            from chronoclade.pathogenwatch import load_api_key
            options.update(api_key=load_api_key(),
                           assembly_cache=Path(plan["results"]) / "assembly_cache")
        else:
            from chronoclade.time_stage.workflow import run_time as runner
            options.pop("cluster_snp_cutoff", None)
        with progress_task("slurm", f"{plan['stage']} task {index}: {task['selection_id']}",
                           details=f"output: {task['output']}"):
            runner(_resolve(path.parent, task["input"], "task input", external=True),
                   output=task["output"], **options)
        # Bind the output to this plan and task; failed retries remove old success.
        result = Path(task["output"]) / (plan["stage"] + ".json")
        write_json(path.parent / f"task-{index}.json", {
            "plan_sha256": file_sha256(path), "selection_id": task["selection_id"],
            "result": reference(path.parent, result),
        })
    except Exception as error:
        reason = str(error)
        if options.get("api_key"):
            reason = reason.replace(options["api_key"], "[redacted]")
        status.update(status="failed", reason=reason)
        progress("slurm", f"Task {index} failed: {reason}; inspect logs and retry array index {index}")
        raise WorkflowError(reason) from None
    else:
        status.update(status="complete")
    finally:
        status.update(finished_at=datetime.now(timezone.utc).isoformat(),
                      elapsed_seconds=round(perf_counter() - started, 3))
        write_json(status_path, status)


def collect(path):
    """Validate successful receipts, retain absent tasks, never rerun inference."""
    path, plan, source = _load_plan(path)
    manifests = []
    missing, statuses = {}, {}
    for index, task in enumerate(plan["tasks"]):
        receipt = path.parent / f"task-{index}.json"
        if not receipt.is_file():
            status_path = path.parent / f"task-{index}.status.json"
            saved_status = _read(status_path) if status_path.is_file() else {}
            statuses[index] = saved_status
            if (saved_status.get("plan_sha256") == file_sha256(path)
                    and saved_status.get("selection_id") == task["selection_id"]):
                detail = saved_status.get("reason", f"Last task state: {saved_status.get('status')}")
            else:
                detail = "No task start or success receipt; inspect scheduler state and logs"
            missing[task["selection_id"]] = f"{detail}; retry array index {index} and rerun collect.slurm"
            continue
        saved = _read(receipt)
        if saved.get("plan_sha256") != file_sha256(path) or saved.get("selection_id") != task["selection_id"]:
            raise WorkflowError("Task receipt does not match this array plan")
        result = _resolve(path.parent, saved["result"], "task result", external=True)
        if result != Path(task["output"]) / (plan["stage"] + ".json"):
            raise WorkflowError("Task receipt points outside its declared output")
        manifests.append(result)
    from chronoclade.execution.accounting import task_accounting, accounting_reason
    accounting = task_accounting(path, plan, statuses, identity_key="selection_id")
    outcomes = {}
    for index, scheduler in accounting.items():
        ident = plan["tasks"][index]["selection_id"]
        missing[ident] += "; " + accounting_reason(scheduler)
        outcomes[ident] = {"reason": missing[ident], "scheduler": scheduler}
    previous = Path(plan["results"]) / (plan["stage"] + "_ensemble.json")
    if previous.is_file():
        _freeze(previous)
    if plan["stage"] == "tree":
        result = register_tree_ensemble(source, manifests, output=plan["results"], force=True)
        for row in result["runs"]:
            if row["selection_id"] in missing:
                row.update(outcomes[row["selection_id"]])
        write_json(Path(plan["results"]) / "tree_ensemble.json", result)
        from chronoclade.stage_readme import write_stage_readme
        write_stage_readme(Path(plan["results"]), "tree", ensemble=True,
                           completed_runs=result["completed_runs"], executor_plan=path)
    else:
        result = run_time_ensemble(source, output=plan["results"], reuse_results=manifests,
                                   collect_only=True, force=True, missing_outcomes=outcomes, **plan["options"])
        from chronoclade.stage_readme import write_stage_readme
        write_stage_readme(Path(plan["results"]), "time", ensemble=True,
                           report=result["report"]["path"], executor_plan=path)
    from chronoclade.progress import progress
    completed = sum(row["status"] == "complete" for row in result["runs"])
    progress("slurm", f"Collected {plan['stage']}: {completed}/{len(result['runs'])} complete; results: {plan['results']}")
    for ident, reason in missing.items():
        progress("slurm", f"Incomplete {ident}: {reason}")
    return 0 if all(row["status"] == "complete" for row in result["runs"]) else 2


def main():
    from chronoclade.progress import configure_progress
    configure_progress()
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    parser.add_argument("operation", choices=["task", "collect"])
    parser.add_argument("index", type=int, nargs="?")
    args = parser.parse_args()
    try:
        if args.operation == "task":
            if args.index is None:
                parser.error("task requires an index")
            run_task(args.plan, args.index)
            return 0
        return collect(args.plan)
    except (ValueError, OSError, WorkflowError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
