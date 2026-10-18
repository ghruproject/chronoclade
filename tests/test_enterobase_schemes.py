import gzip
import io
import json

import pytest

from chronoclade import enterobase_schemes as schemes

SCHEME = 'Escherichia.cgMLSTv1'


def index(data):
    return f'<a href="abc.fasta.gz">abc.fasta.gz</a> 2026-10-17 00:00 {len(data)}\n'.encode()


class Response(io.BytesIO):
    headers = {'Last-Modified': 'test'}


def source(monkeypatch):
    data = gzip.compress(b'>1\nACGT\n>2\nACGA\n', mtime=0)
    calls = []

    def fetch(url):
        calls.append(url)
        return Response(data if url.endswith('.fasta.gz') else index(data))

    monkeypatch.setattr(schemes, '_open_url', fetch)
    return data, calls


def test_install_resume_and_tamper(monkeypatch, tmp_path):
    data, calls = source(monkeypatch)
    first = schemes.install_scheme(SCHEME, tmp_path)
    assert first['files'][0]['allele_count'] == 2
    assert first['hiercc_status'] == 'not_installed'
    assert schemes.install_scheme(SCHEME, tmp_path)['reference_sha256'] == first['reference_sha256']
    assert len(calls) == 2
    (tmp_path/'alleles/abc.fasta.gz').write_bytes(data+b'changed')
    with pytest.raises(schemes.EnterobaseSchemeError, match='changed'):
        schemes.install_scheme(SCHEME, tmp_path)


@pytest.mark.parametrize('body', [b'', b'<a href="../abc.fasta.gz">bad</a> 20',
                                  b'<a href="abc.fasta.gz">a</a> 20\n<a href="abc.fasta.gz">a</a> 20'])
def test_bad_index(body):
    with pytest.raises(schemes.EnterobaseSchemeError):
        schemes.parse_index(body, SCHEME)


@pytest.mark.parametrize('body', [b'>1\nACGT\n>1\nAAAA\n', b'>1\n>2\nACGT\n', b'>1\nPROTEIN\n', b'ACGT\n'])
def test_bad_fasta(tmp_path, body):
    path = tmp_path/'a.gz'
    path.write_bytes(gzip.compress(body))
    with pytest.raises(schemes.EnterobaseSchemeError):
        schemes.validate_alleles(path)


def test_budget_and_snapshot_paths(monkeypatch, tmp_path):
    _, calls = source(monkeypatch)
    with pytest.raises(schemes.EnterobaseSchemeError, match='budget'):
        schemes.install_scheme(SCHEME, tmp_path, max_compressed_bytes=1)
    assert len(calls) == 1
    snapshot = json.loads((tmp_path/'index.json').read_text())
    snapshot['allele_files'][0]['filename'] = '../outside.fasta.gz'
    (tmp_path/'index.json').write_text(json.dumps(snapshot))
    with pytest.raises(schemes.EnterobaseSchemeError, match='entry'):
        schemes.install_scheme(SCHEME, tmp_path)
    assert len(calls) == 1


def test_interrupted_download_can_resume(monkeypatch, tmp_path):
    data, _ = source(monkeypatch)
    original = schemes._open_url
    monkeypatch.setattr(schemes, '_open_url', lambda url: Response(data[:-3]) if url.endswith('.gz') else original(url))
    with pytest.raises(schemes.EnterobaseSchemeError, match='Incomplete'):
        schemes.install_scheme(SCHEME, tmp_path)
    assert not list((tmp_path/'alleles').iterdir())
    monkeypatch.setattr(schemes, '_open_url', original)
    assert schemes.install_scheme(SCHEME, tmp_path)['status'] == 'alleles_ready'


def test_scheme_and_symlink_guard(monkeypatch, tmp_path):
    source(monkeypatch)
    schemes.install_scheme(SCHEME, tmp_path)
    with pytest.raises(schemes.EnterobaseSchemeError):
        schemes.install_scheme('Salmonella.cgMLSTv2', tmp_path)
    (tmp_path/'reference.json').unlink()
    (tmp_path/'reference.json').symlink_to(tmp_path/'outside')
    with pytest.raises(schemes.EnterobaseSchemeError, match='symlink'):
        schemes.install_scheme(SCHEME, tmp_path)


def test_cli_reference_commands(monkeypatch, tmp_path):
    from typer.testing import CliRunner
    from chronoclade.cli import app
    source(monkeypatch)
    runner = CliRunner()
    listing = runner.invoke(app, ['list-schemes'])
    assert listing.exit_code == 0
    assert SCHEME in listing.stdout
    installed = runner.invoke(app, ['setup-scheme', SCHEME, '--destination', str(tmp_path)])
    assert installed.exit_code == 0, installed.stdout
    assert 'Verified 1 loci' in installed.stdout
    assert 'HierCC assignment' in installed.stdout
    bad = runner.invoke(app, ['setup-scheme', '../bad', '--destination', str(tmp_path)])
    assert bad.exit_code == 2
