"""Lazy headless plotting setup shared by offline report writers."""

from threading import Lock

_setup_lock = Lock()
_pyplot_module = None


def headless_pyplot():
    """Select Agg before importing pyplot, including when first called by a worker."""
    global _pyplot_module
    with _setup_lock:
        if _pyplot_module is None:
            import matplotlib

            matplotlib.use("Agg", force=True)
            matplotlib.rcParams["svg.hashsalt"] = "chronoclade"
            from matplotlib import pyplot

            _pyplot_module = pyplot
        return _pyplot_module
