"""One-time user-approved Pasteur OAuth, isolated in the installed plincer runtime."""

import hashlib
import json
import os
from pathlib import Path
import secrets
import subprocess
import tempfile
import time
from typing import Annotated
from urllib.parse import urlencode

import typer
from rich.console import Console

from chronoclade.typing_setup import (
    TypingSetupError,
    _normalise_pasteur_secrets,
    _write_private_json,
)

DEFAULT_DESTINATION = Path('~/.local/share/chronoclade/typing')
PENDING_SECONDS = 3600
WEB_URL = 'https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl'

# Runs only in the plincer venv; rauth/requests are not ChronoClade core dependencies.
# No browser interaction, account password or upstream automatic-login method.
_OAUTH_WORKER = r'''
import json, logging, os, sys
logging.disable(logging.CRITICAL)
from plincer.keycache import create_oauth_service
operation, credentials, pending_path, output_path = sys.argv[1:]
try:
    consumer = json.load(open(credentials))["pasteur"]["consumer"]
    pending = json.load(open(pending_path))
    service = create_oauth_service((consumer["TOKEN"], consumer["TOKEN SECRET"]),
                                  "https://bigsdb.pasteur.fr/api/db", "klebsiella")
    if operation == "begin":
        response = service.get_raw_request_token(params={"oauth_callback":"oob"},
                                                timeout=30, allow_redirects=False)
    elif operation == "complete":
        verifier = sys.stdin.read().strip()
        response = service.get_raw_access_token(pending["request_token"], pending["request_secret"],
                      params={"oauth_verifier":verifier}, timeout=30, allow_redirects=False)
    else:
        raise ValueError("Unsupported operation")
    if response.status_code != 200:
        raise ValueError("OAuth exchange rejected")
    value = response.json()
    result = {"token":value["oauth_token"], "secret":value["oauth_token_secret"],
              "nonce":pending["nonce"]}
    descriptor = os.open(output_path, os.O_WRONLY|os.O_CREAT|os.O_EXCL, 0o600)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(result, stream)
except Exception:
    sys.exit(2)
'''


class PasteurAuthorizationError(ValueError):
    """An actionable failure with no remote response or credential values."""


def _read_private(path, *, optional=False):
    if not path.exists() and not path.is_symlink() and optional:
        return {}
    if path.is_symlink() or not path.is_file():
        raise PasteurAuthorizationError('Authorization files must be regular files, not symbolic links')
    try:
        value = json.loads(path.read_text(encoding='utf-8'))
        if not isinstance(value, dict):
            raise ValueError
        path.chmod(0o600)
        return value
    except (OSError, ValueError):
        raise PasteurAuthorizationError('Authorization files are unreadable or invalid') from None


def _paths(destination):
    root = Path(destination).expanduser().resolve()
    python = root / 'venvs/plincer/bin/python'
    if not python.is_file() or not os.access(python, os.X_OK):
        raise PasteurAuthorizationError('plincer is not installed here; run setup-typing first with the same --destination')
    private = root / 'credential-cache'
    if private.is_symlink() or (private.exists() and not private.is_dir()):
        raise PasteurAuthorizationError('Credential cache must be a regular directory, not a symbolic link')
    private.mkdir(parents=True, exist_ok=True, mode=0o700)
    private.chmod(0o700)
    return python, private


def _atomic_private(path, value):
    if path.is_symlink():
        raise PasteurAuthorizationError('Refusing to replace a symbolic-link credential file')
    temporary = None
    try:
        fd, temporary = tempfile.mkstemp(prefix='.pasteur-', dir=path.parent)
        os.close(fd)
        _write_private_json(Path(temporary), value)
        os.replace(temporary, path)
    except (OSError, TypingSetupError):
        raise PasteurAuthorizationError('Cannot safely save private authorization state') from None
    finally:
        if temporary and Path(temporary).exists():
            Path(temporary).unlink()


def _fingerprint(value):
    try:
        consumer = value['pasteur']['consumer']
        if not all(isinstance(consumer[x], str) and consumer[x].strip() for x in ('TOKEN','TOKEN SECRET')):
            raise ValueError
    except (KeyError, TypeError, ValueError):
        raise PasteurAuthorizationError('Saved Pasteur consumer credentials are invalid') from None
    return hashlib.sha256(json.dumps(consumer, sort_keys=True).encode()).hexdigest()


def _validate_pending(pending, fingerprint):
    try:
        age = time.time() - pending['created_at']
        valid = (pending['consumer_fingerprint'] == fingerprint and 0 <= age <= PENDING_SECONDS
                 and all(isinstance(pending[x], str) and pending[x] for x in
                         ('nonce','request_token','request_secret')))
    except (KeyError, TypeError):
        valid = False
    if not valid:
        raise PasteurAuthorizationError('Pending authorization expired or differs from this consumer; start authorization again')


def _oauth(python, private, operation, credentials, pending, verifier=''):
    with tempfile.TemporaryDirectory(prefix='.pasteur-oauth-', dir=private) as temporary:
        output = Path(temporary) / 'response.json'
        try:
            process = subprocess.run([str(python), '-c', _OAUTH_WORKER, operation,
                                      str(credentials), str(pending), str(output)],
                                     input=verifier, text=True, capture_output=True, timeout=65, check=False)
            if process.returncode != 0:
                raise ValueError
            result = _read_private(output)
            if not all(isinstance(result.get(x), str) and result[x] for x in ('token','secret','nonce')):
                raise ValueError
            return result
        except (OSError, ValueError, subprocess.TimeoutExpired):
            raise PasteurAuthorizationError('Pasteur authorization exchange failed; check connectivity/account access and retry without sharing credentials') from None


def _url(token):
    return WEB_URL + '?' + urlencode({'db':'pubmlst_klebsiella_seqdef', 'page':'authorizeClient', 'oauth_token':token})


def begin_authorization(destination, pasteur_secrets, *, restart=False):
    """Request an official browser consent URL; existing access tokens stay intact."""
    python, private = _paths(destination)
    credentials = private / 'pasteur-authorization-consumer.json'
    pending_path = private / 'pasteur-authorization-pending.json'
    # Normalize privately before committing, so malformed input cannot replace saved state.
    with tempfile.TemporaryDirectory(prefix='.pasteur-consumer-', dir=private) as temporary:
        normalized = Path(temporary) / 'consumer.json'
        try:
            _normalise_pasteur_secrets(Path(pasteur_secrets).expanduser(), normalized, consumer_only=True)
        except TypingSetupError:
            raise PasteurAuthorizationError('Cannot read Pasteur key file; supply the issued client-key letter or supported JSON') from None
        value = _read_private(normalized)
    fingerprint = _fingerprint(value)
    try:
        existing = _read_private(pending_path, optional=True)
    except PasteurAuthorizationError:
        if not restart or pending_path.is_symlink() or not pending_path.is_file():
            raise
        existing = {}  # Explicit restart can replace corrupt pending state, never access tokens.
    if existing and not restart:
        try:
            _validate_pending(existing, fingerprint)
        except PasteurAuthorizationError:
            pass  # Explicit begin starts a fresh request; it never exchanges a stale one.
        else:
            _atomic_private(credentials, value)
            return _url(existing['request_token'])
    pending = {'created_at':time.time(), 'consumer_fingerprint':fingerprint, 'nonce':secrets.token_urlsafe(24)}
    _atomic_private(credentials, value)
    with tempfile.TemporaryDirectory(prefix='.pasteur-begin-', dir=private) as temporary:
        request_path = Path(temporary) / 'request.json'
        _write_private_json(request_path, pending)
        result = _oauth(python, private, 'begin', credentials, request_path)
    if result['nonce'] != pending['nonce']:
        raise PasteurAuthorizationError('Authorization request identity changed; start again')
    pending.update(request_token=result['token'], request_secret=result['secret'])
    _atomic_private(pending_path, pending)
    return _url(result['token'])


def complete_authorization(destination, verifier):
    """Exchange only the saved pending request and atomically preserve other cached keys."""
    if not isinstance(verifier, str) or not verifier.strip():
        raise PasteurAuthorizationError('A verification code is required; it is entered at a hidden prompt')
    python, private = _paths(destination)
    credentials = private / 'pasteur-authorization-consumer.json'
    pending_path = private / 'pasteur-authorization-pending.json'
    pending = _read_private(pending_path)
    _validate_pending(pending, _fingerprint(_read_private(credentials)))
    cache_path = private / 'plincer.json'
    cache = _read_private(cache_path, optional=True)
    if not isinstance(cache.get('pasteur', {}), dict):
        raise PasteurAuthorizationError('Existing Pasteur token cache is invalid; preserve it and resolve before authorizing')
    result = _oauth(python, private, 'complete', credentials, pending_path, verifier.strip())
    if result['nonce'] != pending['nonce'] or _read_private(pending_path) != pending:
        raise PasteurAuthorizationError('Pending authorization changed during exchange; start again')
    cache = _read_private(cache_path, optional=True)  # Preserve unrelated updates made during the exchange.
    if not isinstance(cache.get('pasteur', {}), dict):
        raise PasteurAuthorizationError('Existing Pasteur token cache became invalid; access was not saved')
    cache.setdefault('pasteur', {})['klebsiella:access'] = {'TOKEN':result['token'], 'TOKEN SECRET':result['secret']}
    # These refer to the former grant; other databases/hosts remain untouched.
    for key in ('klebsiella:request','klebsiella:session'):
        cache['pasteur'].pop(key, None)
    _atomic_private(cache_path, cache)
    pending_path.unlink()
    return cache_path


def authorize_pasteur_command(
    destination: Annotated[Path, typer.Option(help='The setup-typing installation destination')] = DEFAULT_DESTINATION,
    pasteur_secrets: Annotated[Path | None, typer.Option(help='Pasteur client-key letter or supported JSON; begin only')] = None,
    complete: Annotated[bool, typer.Option('--complete', help='Complete saved authorization using a hidden verifier prompt')] = False,
    restart: Annotated[bool, typer.Option('--restart', help='Begin a fresh request instead of reusing pending consent')] = False,
):
    """Authorize installed Pasteur typing tools once using your own browser login."""
    console = Console()
    try:
        if complete:
            if restart:
                raise PasteurAuthorizationError('--restart begins a fresh request; do not combine it with --complete')
            if pasteur_secrets is not None:
                raise PasteurAuthorizationError('Use --pasteur-secrets to begin; omit it with --complete')
            _, private = _paths(destination)
            _validate_pending(_read_private(private / 'pasteur-authorization-pending.json'),
                              _fingerprint(_read_private(private / 'pasteur-authorization-consumer.json')))
            verifier = typer.prompt('Pasteur verification code', hide_input=True)
            complete_authorization(destination, verifier)
            console.print('Pasteur access saved privately. Rerun setup-typing with --prepare-db and the same destination/key file.')
        else:
            if pasteur_secrets is None:
                raise PasteurAuthorizationError('Supply --pasteur-secrets with the issued key file to begin')
            url = begin_authorization(destination, pasteur_secrets, restart=restart)
            console.print('Open this official Pasteur URL, sign in and approve the application. Then run this command again with --complete:')
            console.print(url, markup=False, soft_wrap=True)
    except (PasteurAuthorizationError, TypingSetupError) as error:
        # Errors above are deliberately safe; never print subprocess stderr or OAuth responses.
        console.print(f'Error: {error}', style='bold red', markup=False)
        raise typer.Exit(2) from None
    except OSError:
        console.print('Error: Cannot read or safely save private Pasteur authorization files', style='bold red')
        raise typer.Exit(2) from None
