#!/usr/bin/env python3
"""Build the ChronoClade documentation site."""

from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    subprocess.run(["mkdocs", "build", "--strict"], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
