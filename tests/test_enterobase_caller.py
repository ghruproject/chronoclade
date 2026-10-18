"""Public caller adapter preserves source and separates HierCC readiness."""
import gzip
import json
from pathlib import Path

import pytest

from chronoclade import enterobase_caller, typing_setup
from chronoclade.artifacts import file_sha256, write_json


def frozen_reference(root, monkeypatch):
    source_dir = root / 'databases/enterobase/Escherichia.cgMLSTv1'
    source_dir.mkdir(parents=True)
    source = source_dir / 'b0001.fasta.gz'
    with gzip.open(source, 'wt') as handle:
        handle.write('>b0001_1\nACGT\n>b0001_2\nTGCA\n>b0001_3\nACNT\n')
    record = {'locus': 'b0001', 'path': source.name, 'sha256': file_sha256(source)}
    manifest = {'reference_sha256': 'a'*64, 'loci': ['b0001'], 'files': [record]}
    write_json(source_dir / 'reference.json', manifest)
    monkeypatch.setattr(enterobase_caller, 'install_scheme', lambda *args: manifest)
    return source


def test_numeric_adapter_preserves_source_and_records_exclusions(tmp_path, monkeypatch):
    source = frozen_reference(tmp_path, monkeypatch)
    digest = file_sha256(source)
    record = enterobase_caller.prepare_caller_scheme(tmp_path, 'Escherichia.cgMLSTv1', 'ecoli_1')
    folder = Path(record['db_path'])
    with gzip.open(folder / 'b0001.fa.gz', 'rt') as handle:
        assert handle.read() == '>1\nACGT\n>2\nTGCA\n'
    assert file_sha256(source) == digest
    receipt = json.loads((folder / 'normalisation.json').read_text())
    assert receipt['files'][0]['excluded_alleles'] == 1
    assert record['last_updated'] == 'snapshot:' + 'a'*64
    # Reusing a verified normalized file preserves its bytes.
    before = file_sha256(folder / 'b0001.fa.gz')
    enterobase_caller.prepare_caller_scheme(tmp_path, 'Escherichia.cgMLSTv1', 'ecoli_1')
    assert file_sha256(folder / 'b0001.fa.gz') == before
    source.write_bytes(b'changed')
    with pytest.raises(ValueError, match='Frozen allele file changed'):
        enterobase_caller.prepare_caller_scheme(tmp_path, 'Escherichia.cgMLSTv1', 'ecoli_1')


@pytest.mark.parametrize('text', ['>gene_1\nACGT\n>1\nTGCA\n', '>gene_0\nACGT\n', '>gene_1\nNNNN\n'])
def test_normalisation_refuses_ambiguous_numeric_ids_or_empty_catalogue(tmp_path, text):
    source = tmp_path / 'source.gz'
    with gzip.open(source, 'wt') as handle:
        handle.write(text)
    with pytest.raises(ValueError):
        enterobase_caller._normalise(source, tmp_path / 'view.gz')
    assert not (tmp_path / 'view.gz').exists()
    assert not (tmp_path / 'view.gz.part').exists()


def test_setup_uses_public_registry_without_profile_or_auth_download(tmp_path, monkeypatch):
    frozen_reference(tmp_path, monkeypatch)
    calls = []
    def index(root, spec, record, scheme_path, tools):
        calls.append((spec, record))
        path = root / 'databases/mlst-index/ecoli_1'
        path.mkdir(parents=True)
        return path
    monkeypatch.setattr(typing_setup, '_index_cgmlst', index)
    monkeypatch.setattr(typing_setup, '_clone_at', lambda *a: pytest.fail('old downloader used'))
    result = typing_setup._prepare_public_cgmlst('unused', tmp_path, 'ecoli', {})
    assert result['access_mode'] == 'public_unauthenticated'
    assert len(calls) == 1
    assert typing_setup._cgmlst_access_scope(tmp_path, 'ecoli') == 'public_unauthenticated'
    assert not typing_setup._database_state(tmp_path, 'ecoli')['lineage_ready']
    reference = tmp_path / 'databases/enterobase/Escherichia.cgMLSTv1/reference.json'
    reference.write_text('{}')
    assert typing_setup._cgmlst_access_scope(tmp_path, 'ecoli') == 'unverified'


def test_native_index_log_is_opt_in_and_exception_remains_concise(tmp_path, monkeypatch):
    import subprocess
    monkeypatch.setattr(typing_setup.subprocess, 'run', lambda *a, **kw:
        subprocess.CompletedProcess(a[0], 1, stdout='index progress', stderr='native detail'))
    log = tmp_path / 'logs/index.log'
    with pytest.raises(typing_setup.TypingSetupError, match='exited with status 1') as error:
        typing_setup._run(['index-command'], log_path=log)
    assert 'native detail' not in str(error.value)
    assert log.read_text() == 'index progress\nnative detail'
