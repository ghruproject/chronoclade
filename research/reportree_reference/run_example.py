"""Run the upstream Neisseria worked example without modifying its source."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import time


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    base = root / "tmp/reportree-reference"
    upstream = base / "upstream"
    example = upstream / "examples/Neisseria"
    output = base / "reproduced"
    output.mkdir(exist_ok=True)
    expected = [
        output / f"NG_822_{name}.tsv"
        for name in ("metadata_w_partitions", "metrics", "stableRegions", "partitions_summary")
    ]
    if (output / "execution.json").exists():
        if all(p.exists() for p in expected):
            raise SystemExit("Completed run already exists; preserve it rather than overwrite.")
        for name in ("execution.json", "console.log"):
            (output / name).rename(output / f"failed-attempt1-{name}")
    original = next(
        line
        for line in (example / "output/NG_822.log").read_text().splitlines()
        if line.startswith("../../reportree.py ")
    )
    args = shlex.split(original)
    args[0] = str(upstream / "reportree.py")
    args[args.index("--output") + 1] = str(output / "NG_822")
    python = base / "env/bin/python"
    command = [str(python), *args]
    env = os.environ.copy()
    env["PATH"] = f"{python.parent}:{base / 'cgmlst-dists'}:{env['PATH']}"
    start = time.monotonic()
    receipt = {
        "command": command,
        "cwd": str(example),
        "upstream_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=upstream, text=True
        ).strip(),
        "submodules": subprocess.check_output(
            ["git", "submodule", "status"], cwd=upstream, text=True
        ).splitlines(),
        "input_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((example / "input").iterdir())
            if p.is_file()
        },
        "reference_version": "2.0.0 (saved example log, 2023-03-28)",
    }
    (output / "command.json").write_text(json.dumps(receipt, indent=2) + "\n")
    with (output / "console.log").open("w") as log:
        result = subprocess.run(command, cwd=example, env=env, stdout=log, stderr=subprocess.STDOUT)
    success = result.returncode == 0 and all(p.exists() for p in expected)
    receipt.update(
        returncode=result.returncode,
        expected_outputs_present=success,
        elapsed_seconds=time.monotonic() - start,
    )
    receipt["packages"] = subprocess.check_output(
        ["uv", "pip", "freeze", "--python", str(python)], text=True
    ).splitlines()
    (output / "execution.json").write_text(json.dumps(receipt, indent=2) + "\n")
    raise SystemExit(0 if success else 1)


if __name__ == "__main__":
    main()
