"""Independent file snapshots with copy-on-write support and portable fallback."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import sys

from chronoclade.errors import WorkflowError


def _clone_file(source: Path, destination: Path) -> bool:
    """Try an independent copy-on-write snapshot on supported filesystems."""
    try:
        if sys.platform == "darwin":
            import ctypes

            clonefile = ctypes.CDLL(None, use_errno=True).clonefile
            clonefile.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
            clonefile.restype = ctypes.c_int
            return clonefile(os.fsencode(source), os.fsencode(destination), 0) == 0
        if sys.platform.startswith("linux"):
            import fcntl

            with source.open("rb") as original, destination.open("xb") as cloned:
                fcntl.ioctl(cloned.fileno(), 0x40049409, original.fileno())  # FICLONE
            return True
    except (OSError, AttributeError):
        pass
    return False


def copy_evidence_snapshot(source: Path, destination: Path) -> str:
    """Copy evidence cheaply where possible, never sharing mutable inodes."""
    source, destination = Path(source), Path(destination)
    if os.path.lexists(destination):
        raise WorkflowError(f"Evidence destination already exists: {destination}")
    if _clone_file(source, destination):
        return "copy_on_write"
    destination.unlink(missing_ok=True)
    shutil.copyfile(source, destination)
    return "copy"


