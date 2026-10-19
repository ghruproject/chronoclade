"""Native executable discovery for stage preflight checks."""

import shutil
import sys

REQUIRED_TOOLS = ("ska", "iqtree", "Profile", "ClonalFrameML", "verticall", "minimap2", "treetime")


def tool_status() -> dict[str, str | None]:
    return {tool: shutil.which(tool) for tool in REQUIRED_TOOLS}


def native_platform_supported() -> bool:
    """Return whether native tools support this operating system."""
    return sys.platform.startswith(("linux", "darwin"))
