#!/usr/bin/env python3
"""Build the documentation and add the preserved ST239 HTML report."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site"
EXAMPLE_RESULTS = ROOT / "validation" / "st239_baines2015" / "results"


def main() -> None:
    subprocess.run(["mkdocs", "build", "--strict"], cwd=ROOT, check=True)
    destination = SITE / "example-report"
    shutil.copytree(EXAMPLE_RESULTS, destination, dirs_exist_ok=True)
    if not (destination / "report.html").is_file():
        raise RuntimeError("ST239 example report was not copied into the documentation site")


if __name__ == "__main__":
    main()
