import json
from types import SimpleNamespace
from urllib.parse import parse_qs, urlsplit

import pytest
from typer.testing import CliRunner
import typer

from chronoclade import pasteur_authorization as auth


@pytest.fixture
def installation(tmp_path):
    python = tmp_path / 'venvs/plincer/bin/python'
    python.parent.mkdir(parents=True)
    python.write_text('#!/bin/sh\n')
    python.chmod(0o700)
    letter = tmp_path / 'client.txt'
    letter.write_text('Client id: fixture-client\nClient secret: fixture-consumer-secret\n')
    return tmp_path, letter


def fake_oauth(_python, _private, operation, _credentials, pending, verifier=''):
    value = json.loads(pending.read_text())
    if operation == 'complete':
        assert verifier == 'fixture-verifier'
    return {'token':'fixture-' + operation, 'secret':'fixture-token-secret', 'nonce':value['nonce']}


def begin(installation, monkeypatch):
    root, letter = installation
    monkeypatch.setattr(auth, '_oauth', fake_oauth)
    url = auth.begin_authorization(root, letter)
    return root, root / 'credential-cache', url


def test_begin_uses_official_consent_and_protected_state(installation, monkeypatch):
    root, private, url = begin(installation, monkeypatch)
    parsed = urlsplit(url)
    assert parsed.scheme == 'https' and parsed.hostname == 'bigsdb.pasteur.fr'
    assert parse_qs(parsed.query)['oauth_token'] == ['fixture-begin']
    assert 'fixture-token-secret' not in url
    assert not (private / 'plincer.json').exists()
    assert private.stat().st_mode & 0o777 == 0o700
    for name in ('pasteur-authorization-consumer.json','pasteur-authorization-pending.json'):
        assert (private/name).stat().st_mode & 0o777 == 0o600
    assert 'fixture-consumer-secret' not in (private/'pasteur-authorization-pending.json').read_text()


def test_complete_preserves_other_keys_and_only_clears_old_grant_session(installation, monkeypatch):
    root, private, _ = begin(installation, monkeypatch)
    cache = private/'plincer.json'
    cache.write_text(json.dumps({'pubmlst':{'keep':1},'pasteur':{'other:access':{'keep':2},
                      'klebsiella:session':{'old':1},'klebsiella:request':{'old':2}}}))
    result = auth.complete_authorization(root, 'fixture-verifier')
    assert result == cache
    value = json.loads(cache.read_text())
    assert value['pubmlst'] == {'keep':1}
    assert value['pasteur']['other:access'] == {'keep':2}
    assert value['pasteur']['klebsiella:access']['TOKEN'] == 'fixture-complete'
    assert 'klebsiella:session' not in value['pasteur']
    assert 'klebsiella:request' not in value['pasteur']
    assert cache.stat().st_mode & 0o777 == 0o600
    assert not (private/'pasteur-authorization-pending.json').exists()


def test_missing_tools_gives_actionable_error_without_state(tmp_path):
    with pytest.raises(auth.PasteurAuthorizationError, match='setup-typing'):
        auth.begin_authorization(tmp_path, tmp_path/'missing.txt')
    assert not (tmp_path/'credential-cache').exists()


def test_invalid_key_does_not_replace_cache(installation):
    root, letter = installation
    private = root/'credential-cache'
    private.mkdir()
    cache = private/'plincer.json'
    cache.write_text('{"existing":"preserve"}')
    before = cache.read_bytes()
    letter.write_text('invalid fixture')
    with pytest.raises(auth.PasteurAuthorizationError, match='key file'):
        auth.begin_authorization(root, letter)
    assert cache.read_bytes() == before


def test_failed_exchange_keeps_cache_and_pending(installation, monkeypatch):
    root, private, _ = begin(installation, monkeypatch)
    cache = private/'plincer.json'
    cache.write_text('{"existing":"preserve"}')
    before = cache.read_bytes()
    pending = (private/'pasteur-authorization-pending.json').read_bytes()
    def fail(*args, **kwargs):
        raise auth.PasteurAuthorizationError('Safe failure')
    monkeypatch.setattr(auth, '_oauth', fail)
    with pytest.raises(auth.PasteurAuthorizationError, match='Safe failure'):
        auth.complete_authorization(root, 'fixture-verifier')
    assert cache.read_bytes() == before
    assert (private/'pasteur-authorization-pending.json').read_bytes() == pending


@pytest.mark.parametrize('change', ['expired','consumer','nonce'])
def test_pending_identity_lifetime_and_response_nonce(installation, monkeypatch, change):
    root, private, _ = begin(installation, monkeypatch)
    path = private/'pasteur-authorization-pending.json'
    value = json.loads(path.read_text())
    if change == 'expired':
        value['created_at'] -= auth.PENDING_SECONDS+1
    elif change == 'consumer':
        value['consumer_fingerprint'] = 'wrong-consumer'
    else:
        monkeypatch.setattr(auth,'_oauth',lambda *args: {'token':'new','secret':'private','nonce':'wrong'})
    path.write_text(json.dumps(value))
    with pytest.raises(auth.PasteurAuthorizationError):
        auth.complete_authorization(root, 'fixture-verifier')
    assert not (private/'plincer.json').exists()


def test_valid_begin_reuses_pending_request(installation, monkeypatch):
    root, private, url = begin(installation, monkeypatch)
    def forbidden(*args):
        raise AssertionError('An existing valid pending request must not be rotated')
    monkeypatch.setattr(auth, '_oauth', forbidden)
    assert auth.begin_authorization(root, installation[1]) == url
    assert (private/'pasteur-authorization-pending.json').exists()


def test_symbolic_link_cache_is_rejected_before_exchange(installation, monkeypatch):
    root, private, _ = begin(installation, monkeypatch)
    outside = root/'outside.json'
    outside.write_text('{"untouched":true}')
    (private/'plincer.json').symlink_to(outside)
    with pytest.raises(auth.PasteurAuthorizationError, match='symbolic'):
        auth.complete_authorization(root, 'fixture-verifier')
    assert outside.read_text() == '{"untouched":true}'


def test_worker_passes_verifier_only_on_stdin_and_sanitizes_failure(installation, monkeypatch):
    root, private, _ = begin(installation, monkeypatch)
    monkeypatch.undo()
    # Directly exercise the subprocess boundary without running a real OAuth request.
    def process(command, **kwargs):
        assert 'fixture-verifier' not in ' '.join(command)
        assert 'fixture-consumer-secret' not in ' '.join(command)
        assert kwargs['input'] == 'fixture-verifier'
        assert kwargs['capture_output'] and kwargs['timeout'] == 65
        return SimpleNamespace(returncode=2, stderr='REMOTE SECRET RESPONSE', stdout='PRIVATE TOKEN')
    monkeypatch.setattr(auth.subprocess, 'run', process)
    with pytest.raises(auth.PasteurAuthorizationError) as error:
        auth.complete_authorization(root, 'fixture-verifier')
    assert 'REMOTE' not in str(error.value) and 'PRIVATE' not in str(error.value)
    assert not (private/'plincer.json').exists()


def test_cli_complete_uses_hidden_prompt_and_never_prints_token(installation, monkeypatch):
    root, private, _ = begin(installation, monkeypatch)
    calls = []
    def prompt(_message, **kwargs):
        calls.append(kwargs)
        return 'fixture-verifier'
    monkeypatch.setattr(auth.typer,'prompt',prompt)
    app = typer.Typer()
    app.command('authorize-pasteur')(auth.authorize_pasteur_command)
    result = CliRunner().invoke(app, ['--destination', str(root), '--complete'])
    assert result.exit_code == 0, result.output
    assert calls == [{'hide_input':True}]
    for value in ('fixture-verifier','fixture-complete','fixture-token-secret','fixture-consumer-secret'):
        assert value not in result.output
    assert (private/'plincer.json').stat().st_mode & 0o777 == 0o600


def test_pending_remains_valid_for_one_hour(installation, monkeypatch):
    root, private, _ = begin(installation, monkeypatch)
    path = private/'pasteur-authorization-pending.json'
    value = json.loads(path.read_text())
    value['created_at'] -= 50*60
    path.write_text(json.dumps(value))
    auth.complete_authorization(root, 'fixture-verifier')
    assert (private/'plincer.json').exists()


def test_explicit_restart_creates_fresh_nonce_without_touching_access_cache(installation, monkeypatch):
    root, private, _ = begin(installation, monkeypatch)
    path = private/'pasteur-authorization-pending.json'
    old = json.loads(path.read_text())['nonce']
    cache = private/'plincer.json'
    cache.write_text('{"untouched":true}')
    auth.begin_authorization(root, installation[1], restart=True)
    assert json.loads(path.read_text())['nonce'] != old
    assert cache.read_text() == '{"untouched":true}'


def test_atomic_save_failure_preserves_existing_access_file(installation, monkeypatch):
    root, private, _ = begin(installation, monkeypatch)
    cache = private/'plincer.json'
    cache.write_text('{"existing":"preserve"}')
    def fail(*args):
        raise auth.TypingSetupError('private writer failed')
    monkeypatch.setattr(auth, '_write_private_json', fail)
    with pytest.raises(auth.PasteurAuthorizationError, match='safely save'):
        auth.complete_authorization(root, 'fixture-verifier')
    assert cache.read_text() == '{"existing":"preserve"}'
    assert (private/'pasteur-authorization-pending.json').exists()


def test_restart_recovers_corrupt_pending_file_without_manual_deletion(installation, monkeypatch):
    root, private, _ = begin(installation, monkeypatch)
    path = private/'pasteur-authorization-pending.json'
    path.write_text('interrupted invalid pending JSON')
    with pytest.raises(auth.PasteurAuthorizationError):
        auth.begin_authorization(root, installation[1])
    auth.begin_authorization(root, installation[1], restart=True)
    assert json.loads(path.read_text())['request_token'] == 'fixture-begin'
