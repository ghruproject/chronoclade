"""Shared CLI progress without changing scientific artifacts or fingerprints."""

from contextlib import contextmanager
import logging
import sys
from threading import Event, Thread
from time import perf_counter

_LOG = logging.getLogger('chronoclade.progress')


class _CurrentStderr(logging.Handler):
    def emit(self, record):
        print(self.format(record), file=sys.stderr, flush=True)


def configure_progress(*, quiet=False):
    """Configure once per CLI invocation; library imports remain silent."""
    _LOG.handlers = [_CurrentStderr()]
    _LOG.propagate = False
    _LOG.setLevel(logging.WARNING if quiet else logging.INFO)


def progress(stage, message):
    _LOG.info('[%s] %s', stage, message)


@contextmanager
def progress_task(stage, label, *, details='', heartbeat_seconds=30):
    """Report a long operation and elapsed time, including a failed attempt."""
    if not _LOG.isEnabledFor(logging.INFO):
        yield
        return
    started = perf_counter()
    suffix = f'; {details}' if details else ''
    progress(stage, f'Starting {label}{suffix}')
    stopped = Event()

    def heartbeat():
        while not stopped.wait(heartbeat_seconds):
            progress(stage, f'{label} still running; elapsed {perf_counter()-started:.0f}s{suffix}')

    reporter = Thread(target=heartbeat, daemon=True)
    reporter.start()
    try:
        yield
    except BaseException:
        progress(stage, f'{label} stopped before completion; elapsed {perf_counter()-started:.1f}s{suffix}')
        raise
    else:
        progress(stage, f'Completed {label}; elapsed {perf_counter()-started:.1f}s')
    finally:
        stopped.set()
        reporter.join(timeout=1)
