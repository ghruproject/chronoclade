"""Apply upstream ReporTree stability analysis to cached ST258 partitions.

This exports ChronoClade's exact single-linkage history at every integer cut;
it does not claim GrapeTree MSTree and our distance policies are identical.
"""

import csv
import json
import os
from pathlib import Path
import subprocess
import time

from chronoclade.datasets import load_dataset
from research.klebsiella_grouping.st258_compare import dated, display_names


def export_history(sweep, labels, path):
    ids = sweep["sample_ids"]
    parent = list(range(len(ids)))
    representatives = {i: i for i in parent}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    history = iter(sweep["merge_history"])
    pending = next(history, None)
    last = int(sweep["merge_history"][-1]["distance"])
    columns = []
    for cut in range(last + 1):
        while pending is not None and pending["distance"] <= cut:
            children = [find(representatives[node]) for node in pending["children"]]
            target = min(children)
            for child in children:
                parent[child] = target
            representatives[pending["node_id"]] = target
            pending = next(history, None)
        current = [find(i) for i in range(len(ids))]
        columns.append(current)
        # Validate the export against every already saved selected partition.
        for partition in sweep["partitions"]:
            if partition["threshold"] != cut:
                continue
            expected = {frozenset(c) for c in partition["clusters"]}
            actual = {}
            for ident, group in zip(ids, current, strict=True):
                actual.setdefault(group, set()).add(ident)
            assert expected == {frozenset(c) for c in actual.values()}
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["sequence", *[f"single-{cut}x1.0" for cut in range(last + 1)]])
        for i, ident in enumerate(ids):
            writer.writerow([labels[ident], *[f"group_{column[i]}" for column in columns]])
    return last


def main():
    root = Path(__file__).resolve().parents[2]
    base = root / "tmp/reportree-reference"
    output = base / "st258"
    output.mkdir(exist_ok=True)
    if (output / "execution.json").exists():
        raise SystemExit("Preserve the completed analysis; refusing overwrite.")
    dataset = load_dataset(root / "tmp/st258-native-e2e/prepared/dataset.json")
    labels = display_names(dataset.samples)
    assert len(set(labels.values())) == len(labels)
    with (output / "metadata.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["sequence", "country", "date", "role"])
        for sample in dataset.samples:
            writer.writerow(
                [
                    labels[sample["sample_id"]],
                    sample.get("country") or "Unknown",
                    sample.get("collection_date", "") if dated(sample) else "",
                    sample["role"],
                ]
            )
    python = base / "env/bin/python"
    upstream = base / "upstream"
    env = os.environ.copy()
    env["PATH"] = f"{python.parent}:{base / 'cgmlst-dists'}:{env['PATH']}"
    receipt = {
        "upstream_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=upstream, text=True
        ).strip(),
        "dataset_id": dataset.dataset_id,
        "methods": {},
        "AdjustedWallace": 0.99,
        "minimum_consecutive_transitions": 9,
        "scope": "Upstream ReporTree on supplied ChronoClade single-linkage partitions; common pair-coverage mask retained. Not a GrapeTree rerun.",
    }
    for name in ("allele-profiles", "aligned-locus-snps"):
        folder = output / name
        folder.mkdir(exist_ok=True)
        sweep_path = root / f"tmp/science-st258/comparison/{name}-sweep.json"
        last = export_history(json.loads(sweep_path.read_text()), labels, folder / "partitions.tsv")
        command = [
            str(python),
            str(upstream / "reportree.py"),
            "-m",
            str(output / "metadata.tsv"),
            "-p",
            str(folder / "partitions.tsv"),
            "-out",
            str(folder / "ST258"),
            "--partitions2report",
            "stability_regions",
            "-AdjW",
            "0.99",
            "-n",
            "9",
            "--columns_summary_report",
            "n_sequence,first_seq_date,last_seq_date,timespan_days,country,n_country,role",
        ]
        start = time.monotonic()
        with (folder / "console.log").open("w") as log:
            result = subprocess.run(
                command, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT
            )
        expected = [
            folder / f"ST258_{suffix}.tsv"
            for suffix in (
                "metrics",
                "stableRegions",
                "metadata_w_partitions",
                "partitions_summary",
            )
        ]
        success = (
            result.returncode == 0
            and all(p.exists() for p in expected)
            and "ReporTree is done!" in (folder / "console.log").read_text()
        )
        receipt["methods"][name] = {
            "command": command,
            "last_joining_cutoff": last,
            "elapsed_seconds": time.monotonic() - start,
            "returncode": result.returncode,
            "expected_outputs_present": success,
        }
        if not success:
            (output / "failed-execution.json").write_text(json.dumps(receipt, indent=2) + "\n")
            raise SystemExit(1)
    (output / "execution.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
