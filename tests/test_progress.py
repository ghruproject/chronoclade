"""Progress is visible during CLI work and records failure without false success."""
import pytest

from chronoclade.progress import configure_progress, progress_task


@pytest.fixture(autouse=True)
def restore_quiet():
    yield
    configure_progress(quiet=True)


def test_progress_has_start_elapsed_and_completion(capsys):
    configure_progress()
    with progress_task('tree', 'IQ-TREE', details='log: job/iqtree.log'):
        pass
    output = capsys.readouterr().err
    assert '[tree] Starting IQ-TREE; log: job/iqtree.log' in output
    assert 'Completed IQ-TREE; elapsed' in output


def test_failure_does_not_log_success(capsys):
    configure_progress()
    with pytest.raises(ValueError):
        with progress_task('native', 'inference'):
            raise ValueError('deliberate failure')
    output = capsys.readouterr().err
    assert 'stopped before completion' in output
    assert 'Completed' not in output


def test_quiet_still_runs_operation(capsys):
    configure_progress(quiet=True)
    executed = []
    with progress_task('native', 'inference'):
        executed.append(True)
    assert executed == [True]
    assert capsys.readouterr().err == ''
