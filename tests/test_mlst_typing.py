"""Classical ST recovery requires exact primary-scheme calls and saved evidence."""

from dataclasses import asdict
import json
from pathlib import Path
import subprocess

import pytest

from chronoclade.artifacts import file_sha256
from chronoclade.datasets import DatasetError, LocusCatalogue, load_dataset
from chronoclade.mlst_typing import PRIMARY_SCHEMES, fill_missing_mlst
from chronoclade.prepare_stage import run_prepare


@pytest.fixture
def caller(tmp_path, monkeypatch):
    from chronoclade import mlst_typing

    executable = tmp_path / 'bin' / 'mlst'
    executable.parent.mkdir()
    executable.write_text('fixture executable identity\n')
    database = tmp_path / 'db'
    for scheme, loci in PRIMARY_SCHEMES.values():
        directory = database / 'pubmlst' / scheme
        directory.mkdir(parents=True)
        (directory / (scheme + '.txt')).write_text('ST\t' + '\t'.join(sorted(loci)) + '\n')
        for locus in loci:
            (directory / (locus + '.tfa')).write_text(f'>{locus}_1\nACTG\n')
    (database / 'blast').mkdir()
    (database / 'blast' / 'mlst.fa.nhr').write_text('fixture index')
    monkeypatch.setattr(mlst_typing, '_resources', lambda: (
        executable, database, {'PATH': str(executable.parent)}, 'mlst fixture-1'))
    state = {'calls': [], 'change': lambda output: output, 'returncode': 0}

    def run(command, **kwargs):
        assert command[:1] == [str(executable)]
        assert '--threads' in command and command[command.index('--threads') + 1] == '1'
        assert kwargs['timeout'] == 600 and kwargs['check'] is False
        state['calls'].append(command)
        scheme = command[command.index('--scheme') + 1]
        loci = next(loci for expected, loci in PRIMARY_SCHEMES.values() if scheme == expected)
        output = [{'filename': command[-1], 'scheme': scheme, 'sequence_type': '149',
                   'alleles': {locus: str(index + 1) for index, locus in enumerate(sorted(loci))}}]
        raw = Path(command[command.index('--json') + 1])
        raw.write_text(json.dumps(state['change'](output)))
        return subprocess.CompletedProcess(command, state['returncode'], 'caller stdout\n', 'caller stderr\n')

    monkeypatch.setattr(mlst_typing.subprocess, 'run', run)
    return state, executable, database


def sample(tmp_path, *, species='Staphylococcus aureus', **changes):
    assembly = tmp_path / 'query.fasta'
    assembly.write_text('>contig\nACTGACTGACTG\n')
    return {'sample_id': 'query', 'role': 'input', 'species': species,
            'assembly_reference': assembly.name, **changes}


def test_supplied_numeric_st_is_reused_without_tool_or_assembly(tmp_path, monkeypatch):
    from chronoclade import mlst_typing

    def forbidden():
        raise AssertionError('Supplied ST must not probe a caller or database')
    monkeypatch.setattr(mlst_typing, '_resources', forbidden)
    supplied = {'sample_id': 'q', 'role': 'input', 'mlst_st': '149',
                'mlst_scheme': 'saureus', 'mlst_source': 'verified_metadata'}
    rows, audit = fill_missing_mlst([supplied], tmp_path, tmp_path / 'evidence')
    assert rows[0] == supplied
    assert audit == [{'sample_id': 'q', 'status': 'reused', 'st': '149',
                      'scheme': 'saureus', 'source': 'verified_metadata'}]


@pytest.mark.parametrize('species', list(PRIMARY_SCHEMES))
def test_exact_complete_primary_scheme_records_st_alleles_and_evidence(tmp_path, caller, species):
    row = sample(tmp_path, species=species)
    rows, audit = fill_missing_mlst([row], tmp_path, tmp_path / 'evidence')
    typed, result = rows[0], audit[0]
    scheme, loci = PRIMARY_SCHEMES[species]
    assert typed['mlst_st'] == '149' and typed['mlst_scheme'] == scheme
    assert typed['mlst_status'] == 'resolved' and typed['mlst_source'] == 'classical_mlst'
    assert set(typed['mlst_alleles']) == loci
    assert typed['mlst_assembly_sha256'] == file_sha256(tmp_path / row['assembly_reference'])
    assert typed['mlst_database_sha256'] == result['database']['sha256']
    assert result['tool_sha256'] == file_sha256(caller[1])
    assert result['search_field'] == 'mlst' and result['status'] == 'resolved'
    assert file_sha256(tmp_path / 'evidence' / result['raw_result']) == result['raw_result_sha256']
    assert len(result['database']['files']) == 8
    assert json.loads((tmp_path / 'evidence/audit.json').read_text())['samples'][0]['st'] == '149'
    assert row.get('mlst_st') is None  # Caller does not mutate supplied sample rows.


@pytest.mark.parametrize('st', [None, '-', '0', '~149', '149,150', '149?'])
def test_missing_novel_or_ambiguous_st_does_not_become_exact_assignment(tmp_path, caller, st):
    def change(output):
        output[0]['sequence_type'] = st
        return output
    caller[0]['change'] = change
    rows, audit = fill_missing_mlst([sample(tmp_path)], tmp_path, tmp_path / 'evidence')
    assert rows[0].get('mlst_st') is None
    assert audit[0]['status'] == 'unresolved'
    assert audit[0]['reason'] == 'novel_partial_missing_or_ambiguous_mlst_profile'


@pytest.mark.parametrize('defect', ['missing', 'novel', 'ambiguous', 'extra'])
def test_numeric_st_does_not_override_incomplete_or_inexact_alleles(tmp_path, caller, defect):
    def change(output):
        alleles = output[0]['alleles']
        if defect == 'missing':
            alleles.pop('arcC')
        elif defect == 'extra':
            alleles['unexpected'] = '1'
        else:
            alleles['arcC'] = '~1' if defect == 'novel' else '1,2'
        return output
    caller[0]['change'] = change
    rows, audit = fill_missing_mlst([sample(tmp_path)], tmp_path, tmp_path / 'evidence')
    assert rows[0].get('mlst_st') is None and audit[0]['status'] == 'unresolved'


@pytest.mark.parametrize('defect', ['scheme', 'identity', 'multiple', 'empty'])
def test_caller_scheme_and_assembly_identity_must_match_exactly(tmp_path, caller, defect):
    def change(output):
        if defect == 'scheme':
            output[0]['scheme'] = 'ecoli'  # Same species names cannot establish scheme identity.
        elif defect == 'identity':
            output[0]['filename'] = str(tmp_path / 'other.fasta')
        elif defect == 'multiple':
            return output + output
        else:
            return []
        return output
    caller[0]['change'] = change
    rows, audit = fill_missing_mlst([sample(tmp_path)], tmp_path, tmp_path / 'evidence')
    assert rows[0].get('mlst_st') is None and audit[0]['status'] == 'failed'
    assert audit[0]['reason'] == 'mlst_call_or_output_validation_failed'


def test_missing_tool_and_missing_assembly_are_distinct_nonfatal_audits(tmp_path, monkeypatch):
    from chronoclade import mlst_typing
    monkeypatch.setattr(mlst_typing, '_resources', lambda: None)
    assembly = sample(tmp_path)
    absent = {**assembly, 'sample_id': 'absent', 'assembly_reference': None}
    rows, audit = fill_missing_mlst([assembly, absent], tmp_path, tmp_path / 'evidence')
    assert all(row.get('mlst_st') is None for row in rows)
    assert audit[0]['reason'] == 'mlst_tool_unavailable_run_pixi_install_--frozen_-e_mlst'
    assert audit[1]['reason'] == 'local_assembly_required_for_missing_st'


def test_missing_database_and_caller_exit_do_not_claim_a_resolved_st(tmp_path, caller):
    row = sample(tmp_path)
    caller[0]['returncode'] = 1
    rows, audit = fill_missing_mlst([row], tmp_path, tmp_path / 'failed')
    assert rows[0].get('mlst_st') is None
    assert audit[0]['reason'] == 'mlst_caller_nonzero_exit' and audit[0]['exit_code'] == 1
    (caller[2] / 'pubmlst/saureus/arcC.tfa').unlink()
    rows, audit = fill_missing_mlst([row], tmp_path, tmp_path / 'missing-db')
    assert rows[0].get('mlst_st') is None
    assert audit[0]['reason'] == 'mlst_scheme_or_blast_database_unavailable'
    assert len(caller[0]['calls']) == 1


def frozen_profiles(tmp_path):
    row = sample(tmp_path)
    row.update(cgmlst_scheme='test', cgmlst_scheme_version='1', cgmlst_status='resolved',
               cgmlst_profile={'a': '1', 'b': '2'})
    source = tmp_path / 'profiles.json'
    source.write_text(json.dumps([row]))
    catalogues = tmp_path / 'catalogues.json'
    catalogues.write_text(json.dumps([asdict(LocusCatalogue('test', '1', ('a', 'b'), 'fixture'))]))
    return source, catalogues


def test_prepare_calls_missing_st_after_metadata_merge_and_preserves_audit(tmp_path, caller):
    source, catalogues = frozen_profiles(tmp_path)
    metadata = tmp_path / 'metadata.csv'
    metadata.write_text('sample_id,country\nquery,Greece\n')
    result = run_prepare(source, tmp_path / 'prepared', catalogues=catalogues, metadata=metadata)
    dataset = load_dataset(result.dataset_manifest)
    row = dataset.samples[0]
    assert row['country'] == 'Greece' and row['mlst_st'] == '149'
    assert row['mlst_scheme'] == 'saureus' and row['mlst_status'] == 'resolved'
    assert (result.dataset_manifest.parent / row['assembly_reference']).is_file()
    audit = json.loads(result.audit_path.read_text())
    assert audit['samples'][0]['mlst_st_ready'] is True
    assert audit['samples'][0]['mlst_st'] == '149'
    assert 'Classical MLST' in result.report_path.read_text()
    assert 'ST149 (saureus; resolved)' in result.report_path.read_text()
    classical = result.dataset_manifest.parent / 'sources/classical_mlst/audit.json'
    assert json.loads(classical.read_text())['samples'][0]['status'] == 'resolved'
    assert any(p['field'] == 'mlst_st' and p['source'] == 'classical_mlst' and p['selected']
               for p in dataset.provenance)


def test_frozen_metadata_supplies_verified_st_before_any_caller_probe(tmp_path, monkeypatch):
    from chronoclade import mlst_typing
    source, catalogues = frozen_profiles(tmp_path)
    metadata = tmp_path / 'metadata.csv'
    metadata.write_text('sample_id,mlst_st,mlst_scheme\nquery,22,saureus\n')

    def forbidden():
        raise AssertionError('Metadata-supplied ST must prevent caller discovery')
    monkeypatch.setattr(mlst_typing, '_resources', forbidden)
    result = run_prepare(source, tmp_path / 'prepared', catalogues=catalogues, metadata=metadata)
    dataset = load_dataset(result.dataset_manifest)
    assert dataset.samples[0]['mlst_st'] == '22' and dataset.samples[0]['mlst_scheme'] == 'saureus'
    assert dataset.parameters['classical_mlst']['samples'][0]['status'] == 'reused'
    assert any(p['field'] == 'mlst_st' and p['source'] == 'user' and p['selected']
               for p in dataset.provenance)


def test_frozen_metadata_rejects_nonnumeric_st(tmp_path):
    source, catalogues = frozen_profiles(tmp_path)
    metadata = tmp_path / 'metadata.csv'
    metadata.write_text('sample_id,mlst_st\nquery,ST22\n')
    with pytest.raises(DatasetError, match='positive numeric ST'):
        run_prepare(source, tmp_path / 'prepared', catalogues=catalogues, metadata=metadata)


def test_user_st_override_clears_native_alleles_hashes_and_tool_evidence(tmp_path, caller):
    source, catalogues = frozen_profiles(tmp_path)
    native = run_prepare(source, tmp_path / 'native', catalogues=catalogues)
    original = load_dataset(native.dataset_manifest).samples[0]
    assert original['mlst_st'] == '149' and original['mlst_alleles']
    assert original['mlst_database_sha256'] and original['mlst_tool_version']
    metadata = tmp_path / 'verified-metadata.csv'
    metadata.write_text('sample_id,mlst_st,mlst_scheme\nquery,22,saureus\n')
    revised = run_prepare(native.dataset_manifest, tmp_path / 'revised', metadata=metadata)
    row = load_dataset(revised.dataset_manifest).samples[0]
    assert row['mlst_st'] == '22' and row['mlst_scheme'] == 'saureus'
    assert row['mlst_status'] == 'supplied' and row['mlst_source'] == 'user_metadata'
    for field in ('mlst_alleles', 'mlst_database_sha256', 'mlst_tool_version',
                  'mlst_assembly_sha256', 'mlst_reason'):
        assert not row.get(field), f'{field} must not be attributed to the overridden ST'
    assert len(caller[0]['calls']) == 1  # The verified replacement is reused without another call.
    classical = json.loads((revised.dataset_manifest.parent / 'sources/classical_mlst/audit.json').read_text())
    assert classical['samples'][0]['status'] == 'reused'
    assert classical['samples'][0]['source'] == 'user_metadata'
