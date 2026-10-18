"""Caller-only S. aureus setup must not invent a lineage dependency."""
import json
from pathlib import Path

from chronoclade import typing_setup


def installed_caller(root):
    return {'status': 'installed', 'version': '8.0.0', 'commit': 'fixture',
            'command': ['node', 'index.js']}


def test_saureus_software_setup_reports_unready_profile_reference(tmp_path, monkeypatch):
    monkeypatch.setattr(typing_setup.shutil, 'which', lambda name: '/usr/bin/uv')
    monkeypatch.setattr(typing_setup, '_install_mlst', installed_caller)
    def forbidden(*args, **kwargs):
        raise AssertionError('S. aureus profile setup must not install a lineage assigner')
    monkeypatch.setattr(typing_setup, '_install_python_tool', forbidden)
    _, report = typing_setup.setup_typing(tmp_path, species='S. aureus')
    assert report['species'] == ['saureus']
    assert set(report['tools']) == {'mlst'}
    assert not report['run_ready'] and not report['cgmlst_ready']
    state = report['databases']['saureus']
    assert not state['cgmlst_ready'] and not state['lineage_ready']
    assert state['lineage_assignment'] == 'not_required_caller_only'
    assert state['lineage_database'] is None
    assert all('plincer' not in path and 'hclink' not in path for path in state['missing_resources'])
    config = json.loads(Path(report['query_config_path']).read_text())
    organism = config['organisms']['Staphylococcus aureus']
    assert organism['cgmlst']['scheme'] == 'saureus_1'
    assert organism['mode'] == 'profiles-only'
    assert 'assignment' not in organism


def test_saureus_ready_requires_validated_versioned_caller_index(tmp_path, monkeypatch):
    monkeypatch.setattr(typing_setup.shutil, 'which', lambda name: '/usr/bin/uv')
    monkeypatch.setattr(typing_setup, '_install_mlst', installed_caller)
    def prepare(uv, root, name, tools, **kwargs):
        assert name == 'saureus'
        assert typing_setup.SPECIES[name]['mlst_scheme_id'] == '2'
        directory = root/'databases/mlst-schemes'
        directory.mkdir(parents=True)
        (directory/'saureus_1-selection.json').write_text(json.dumps({'schemes': [
            {'shortname': 'saureus_1', 'scheme_id': '2', 'last_updated': '2024-12-31'}]}))
        return {'status': 'built', 'access_mode': 'public_unauthenticated'}
    monkeypatch.setattr(typing_setup, '_prepare_public_cgmlst', prepare)
    calls = []
    def validate(path, scheme, version):
        calls.append((scheme, version))
        assert scheme == 'saureus_1' and version == '2024-12-31'
        return ['SAUR0001']
    monkeypatch.setattr('chronoclade.query_typing.validate_cgmlst_index', validate)
    _, report = typing_setup.setup_typing(tmp_path, species='saureus', prepare_db=True)
    assert calls
    assert report['run_ready'] and report['cgmlst_ready']
    assert not report['databases']['saureus']['lineage_ready']
    config = json.loads(Path(report['query_config_path']).read_text())
    assert config['ready'] and config['cgmlst_ready']
    assert 'assignment' not in config['organisms']['Staphylococcus aureus']


def test_saureus_unpinned_or_invalid_index_is_not_ready(tmp_path, monkeypatch):
    def invalid(*args, **kwargs):
        raise ValueError('invalid index')
    monkeypatch.setattr('chronoclade.query_typing.validate_cgmlst_index', invalid)
    monkeypatch.setattr(typing_setup, '_database_version', lambda *args: '2024-12-31')
    state = typing_setup._database_state(tmp_path, 'saureus')
    assert not state['ready'] and not state['cgmlst_ready']
    assert len(state['missing_resources']) == 1
    assert state['missing_resources'][0].endswith('mlst-index/saureus_1')


def test_native_index_log_is_visible_before_command_completes(tmp_path):
    import sys
    import threading
    import time

    log = tmp_path/'index.log'
    release = tmp_path/'release'
    results = []
    errors = []
    script = (
        "from pathlib import Path; import sys,time; "
        "print('index stdout', flush=True); print('index stderr', file=sys.stderr, flush=True); "
        "release=Path(sys.argv[1]); "
        "exec('while not release.exists():\\n time.sleep(0.01)'); "
        "print('finished', flush=True)"
    )
    def run():
        try:
            results.append(typing_setup._run([sys.executable, '-c', script, str(release)],
                                            log_path=log, timeout=10))
        except BaseException as error:
            errors.append(error)
    worker = threading.Thread(target=run)
    worker.start()
    try:
        deadline = time.monotonic()+5
        while time.monotonic() < deadline:
            contents = log.read_text() if log.exists() else ''
            if 'index stdout' in contents and 'index stderr' in contents:
                break
            time.sleep(0.01)
        assert worker.is_alive()
        assert 'index stdout' in contents and 'index stderr' in contents
        assert 'finished' not in contents
    finally:
        release.touch()
        worker.join(timeout=10)
    assert not errors and not worker.is_alive()
    assert results == ['index stdout\nfinished']  # stderr must not contaminate stdout.
    assert 'finished' in log.read_text()
