"""Parallel preparation preserves the whole cohort and never types in collection."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.datasets import load_dataset
from chronoclade.execution.jobs import SlurmResources
from chronoclade.execution import prepare_slurm as m
from test_query_typing import indexed_scheme, profile


@pytest.fixture
def prepared_plan(tmp_path, monkeypatch, request):
    from chronoclade import query_typing, mlst_typing
    monkeypatch.delenv('SLURM_CPUS_PER_TASK', raising=False)
    monkeypatch.setattr(m, 'execution_provenance', lambda *a: {'source': 'stable-test-fixture'})
    monkeypatch.setattr(mlst_typing, '_resources', lambda: None)
    if getattr(request, 'param', False):
        db = tmp_path / 'classical-db'
        scheme, loci = mlst_typing.PRIMARY_SCHEMES['escherichia coli']
        directory = db / 'pubmlst' / scheme
        directory.mkdir(parents=True)
        for name in [scheme + '.txt', *(locus + '.tfa' for locus in loci)]:
            (directory / name).write_text('controlled reference fixture')
        (db / 'blast').mkdir()
        (db / 'blast/mlst.fa.nhr').write_text('controlled BLAST fixture')
        executable = tmp_path / 'classical-caller'
        executable.write_text('#!' + sys.executable + '\nimport json,sys\n'
            + 'from pathlib import Path\n'
            + 'output=sys.argv[sys.argv.index("--json")+1]\n'
            + 'call={"filename":str(Path(sys.argv[-1]).resolve()),"scheme":' + repr(scheme)
            + ',"sequence_type":"95","alleles":' + repr({locus: '1' for locus in loci}) + '}\n'
            + 'Path(output).write_text(json.dumps([call]))\n')
        executable.chmod(0o700)
        monkeypatch.setattr(mlst_typing, '_resources',
                            lambda: (executable, db, dict(os.environ), 'fixture-mlst-1'))
    version = 'snapshot:' + 'a' * 64
    indexed_scheme(tmp_path / 'index', version=version, genes=('a', 'b', 'c', 'd'))
    config = tmp_path / 'typing.json'
    config.write_text(json.dumps({'ready': False, 'organisms': {'Escherichia coli': {'cgmlst': {
        'command': [sys.executable], 'scheme': 'ecoli_1', 'scheme_version': version,
        'tool_version': 'fixture-1', 'index_dir': str(tmp_path / 'index'),
    }}}}))
    source = tmp_path / 'queries.csv'
    with source.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['sample_id', 'assembly', 'species', 'country',
                                                   'collection_date', 'hospital'])
        writer.writeheader()
        for ident in ['second', 'first']:
            assembly = tmp_path / (ident + '.fa')
            assembly.write_text('>' + ident + '\nACGT\n')
            writer.writerow(dict(sample_id=ident, assembly=assembly.name, species='Escherichia coli',
                                 country='UK', collection_date='2020', hospital=ident + '-hospital'))
    calls = []
    def native(command, data, destination):
        calls.append(data)
        Path(destination).write_text(json.dumps(profile()))
        return profile()
    monkeypatch.setattr(query_typing, '_run', native)
    # This test concerns assembly tasks and reuse, not caller-reference authorization.
    monkeypatch.setattr('chronoclade.pasteur_profile_normalization.reference_descriptor',
                        lambda *a, **k: {'fixture': 'content-bound separately in reference tests'})
    plan = m.write_slurm_prepare(source, typing_config=config, output=tmp_path / 'jobs',
                                results=tmp_path / 'prepared', resources=SlurmResources(cpus=1),
                                max_parallel=2)
    return plan, source, calls


def test_array_script_scope_and_cli_generation(prepared_plan, tmp_path):
    plan, source, _ = prepared_plan
    assert '#SBATCH --array=0-1%2' in (plan.parent / 'array.slurm').read_text()
    for script in ['array.slurm', 'collect.slurm', 'submit.sh']:
        subprocess.run(['bash', '-n', str(plan.parent / script)], check=True)
    config = tmp_path / 'typing.json'
    outcome = CliRunner().invoke(app, ['slurm-prepare', str(source), '--typing-config', str(config),
        '--out', str(tmp_path / 'more-jobs'), '--results', str(tmp_path / 'more-prepared')])
    assert outcome.exit_code == 0, outcome.output
    assert 'Submit:' in outcome.output


def test_workers_and_collector_preserve_order_profiles_metadata_and_evidence(prepared_plan, monkeypatch):
    plan, _, calls = prepared_plan
    m.run_task(plan, 1)
    m.run_task(plan, 0)
    assert len(calls) == 2
    monkeypatch.setattr('chronoclade.mlst_typing.fill_missing_mlst',
                        lambda *a, **k: pytest.fail('Collector performed MLST inference'))
    monkeypatch.setattr('chronoclade.query_typing.type_query_assemblies',
                        lambda *a, **k: pytest.fail('Collector performed cgMLST inference'))
    assert m.collect(plan) == 0
    dataset = load_dataset(plan.parent.parent / 'prepared/dataset.json')
    assert dataset.sample_ids == ('second', 'first')
    assert [r['hospital'] for r in dataset.samples] == ['second-hospital', 'first-hospital']
    assert all(r['role'] == 'input' and r['mlst_status'] == 'unavailable' for r in dataset.samples)
    assert len(dataset.profiles) == 1
    assert dataset.profiles[0].profile('second')['a'] == '1'
    assert dataset.profiles[0].profile('first')['c'] == 'a' * 40
    assert len(dataset.profiles[0].catalogue.loci) == 4
    assert dataset.profiles[0].catalogue.complete
    assert len(dataset.parameters['classical_mlst']['samples']) == 2
    assert len(list((plan.parent.parent / 'prepared/sources/provider').rglob('cgmlst.json'))) == 2
    assert m.collect(plan) == 0


def test_retry_reuses_only_verified_successes_and_finishes_missing_tasks(prepared_plan):
    plan, _, calls = prepared_plan
    m.run_task(plan, 0)
    assert m.collect(plan) == 2
    assert not (plan.parent.parent / 'prepared').exists()
    m.run_task(plan, 0)
    assert len(calls) == 1
    m.run_task(plan, 1)
    assert len(calls) == 2
    assert m.collect(plan) == 0


def test_failed_native_worker_resumes_checkpoint_without_discarding_success(prepared_plan, monkeypatch):
    from chronoclade import query_typing
    plan, _, calls = prepared_plan
    m.run_task(plan, 0)
    native = query_typing._run
    monkeypatch.setattr(query_typing, '_run', lambda *a, **k: (_ for _ in ()).throw(ValueError('native failure')))
    with pytest.raises(ValueError, match='native failure'):
        m.run_task(plan, 1)
    assert json.loads((plan.parent / 'task-1.status.json').read_text())['status'] == 'failed'
    assert m.collect(plan) == 2
    monkeypatch.setattr(query_typing, '_run', native)
    m.run_task(plan, 0)
    m.run_task(plan, 1)
    assert len(calls) == 2
    assert m.collect(plan) == 0


@pytest.mark.parametrize('which', ['assembly', 'csv', 'task-input', 'task-result', 'execution'])
def test_changed_inputs_and_cached_artifacts_reject_unsafe_reuse(prepared_plan, which, monkeypatch):
    plan, source, _ = prepared_plan
    m.run_task(plan, 0)
    if which == 'assembly':
        (source.parent / 'second.fa').write_text('>second\nTGCA\n')
    elif which == 'csv':
        source.write_text(source.read_text() + '\n')
    elif which == 'task-input':
        p = plan.parent / 'inputs/000000.csv'
        p.write_text(p.read_text() + '\n')
    elif which == 'execution':
        monkeypatch.setattr(m, 'execution_provenance', lambda *a: {'source': 'changed'})
    else:
        p = plan.parent / 'runs/000000/report.html'
        p.write_text(p.read_text() + 'changed')
    with pytest.raises((ValueError, RuntimeError)):
        m.run_task(plan, 0)


@pytest.mark.parametrize('which', ['metadata', 'task-id', 'order', 'input-reference'])
def test_plan_edits_are_rejected_before_work(prepared_plan, which):
    plan, _, calls = prepared_plan
    value = json.loads(plan.read_text())
    if which == 'metadata':
        value['records'][0]['hospital'] = 'different metadata'
    elif which == 'task-id':
        value['tasks'][0]['sample_id'] = 'outsider'
    elif which == 'order':
        value['tasks'].reverse()
    else:
        value['tasks'][0]['input'] = value['tasks'][1]['input']
    plan.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='plan checksum changed'):
        m.run_task(plan, 0)
    assert not calls


def test_collector_recovers_after_publication_before_receipt(prepared_plan, monkeypatch):
    plan, _, _ = prepared_plan
    m.run_task(plan, 0)
    m.run_task(plan, 1)
    assert m.collect(plan) == 0
    (plan.parent / 'collection-result.json').unlink()
    monkeypatch.setattr('chronoclade.prepare_stage.run_prepare',
                        lambda *a, **k: pytest.fail('Published result was needlessly recomputed'))
    assert m.collect(plan) == 0
    assert (plan.parent / 'collection-result.json').is_file()


def test_interrupted_evidence_copy_is_not_adopted(prepared_plan, monkeypatch):
    plan, _, _ = prepared_plan
    m.run_task(plan, 0)
    m.run_task(plan, 1)
    copytree = m.shutil.copytree
    def interrupted(source, destination, **kwargs):
        Path(destination).mkdir()
        (Path(destination) / 'partial').write_text('incomplete')
        raise OSError('interrupted evidence copy')
    monkeypatch.setattr(m.shutil, 'copytree', interrupted)
    with pytest.raises(OSError, match='interrupted evidence copy'):
        m.collect(plan)
    assert not (plan.parent / 'collected-sources/task-000000').exists()
    monkeypatch.setattr(m.shutil, 'copytree', copytree)
    assert m.collect(plan) == 0


def test_changed_retained_evidence_is_rejected_before_publication(prepared_plan, monkeypatch):
    from chronoclade import prepare_stage
    plan, _, _ = prepared_plan
    m.run_task(plan, 0)
    m.run_task(plan, 1)
    prepare = prepare_stage.run_prepare
    monkeypatch.setattr(prepare_stage, 'run_prepare',
                        lambda *a, **k: (_ for _ in ()).throw(ValueError('stop before publication')))
    with pytest.raises(ValueError, match='stop before publication'):
        m.collect(plan)
    retained = next((plan.parent / 'collected-sources/task-000000').rglob('*.json'))
    retained.write_text(retained.read_text() + '\n')
    monkeypatch.setattr(prepare_stage, 'run_prepare', prepare)
    with pytest.raises(ValueError, match='Retained worker evidence differs'):
        m.collect(plan)
    assert not (plan.parent.parent / 'prepared').exists()


def test_parallel_collection_matches_normal_serial_prepare(prepared_plan):
    from chronoclade.prepare_stage import run_prepare
    plan, source, _ = prepared_plan
    m.run_task(plan, 0)
    m.run_task(plan, 1)
    assert m.collect(plan) == 0
    serial = run_prepare(source, source.parent / 'serial', input_kind='assemblies',
                         typing_config=source.parent / 'typing.json', enrich_metadata=False)
    left = load_dataset(plan.parent.parent / 'prepared/dataset.json')
    right = load_dataset(serial.dataset_manifest)
    assert left.sample_ids == right.sample_ids
    assert [{k: v for k, v in r.items() if k != 'assembly_reference'} for r in left.samples] == [
        {k: v for k, v in r.items() if k != 'assembly_reference'} for r in right.samples]
    assert [{k: v for k, v in row.items() if k != 'retained_evidence_directory'}
            for row in left.parameters['classical_mlst']['samples']] == right.parameters['classical_mlst']['samples']
    assert len(left.profiles) == len(right.profiles) == 1
    assert left.profiles[0].catalogue == right.profiles[0].catalogue
    for ident in left.sample_ids:
        assert left.profiles[0].profile(ident) == right.profiles[0].profile(ident)


@pytest.mark.parametrize('prepared_plan', [True], indirect=True)
def test_aggregated_mlst_audit_resolves_retained_raw_evidence(prepared_plan):
    from chronoclade.artifacts import file_sha256
    plan, _, _ = prepared_plan
    m.run_task(plan, 0)
    m.run_task(plan, 1)
    assert m.collect(plan) == 0
    root = plan.parent.parent / 'prepared'
    dataset = load_dataset(root / 'dataset.json')
    for audit in dataset.parameters['classical_mlst']['samples']:
        assert audit['status'] == 'resolved'
        retained = root / audit['retained_evidence_directory']
        original = json.loads((retained / 'audit.json').read_text())['samples'][0]
        assert audit['raw_result'] == original['raw_result']
        raw = retained / audit['raw_result']
        assert file_sha256(raw) == audit['raw_result_sha256']
        assert json.loads(raw.read_text())[0]['sequence_type'] == '95'



def test_failure_receipt_survives_unreadable_credential_configuration(prepared_plan, monkeypatch):
    plan, _, _ = prepared_plan
    monkeypatch.setattr('chronoclade.prepare_stage.run_prepare',
                        lambda *a, **k: (_ for _ in ()).throw(ValueError('original task failure')))
    monkeypatch.setattr('chronoclade.pathogenwatch.load_api_key',
                        lambda *a, **k: (_ for _ in ()).throw(RuntimeError('credential file unreadable')))
    with pytest.raises(ValueError, match='original task failure'):
        m.run_task(plan, 0)
    status = json.loads((plan.parent / 'task-0.status.json').read_text())
    assert status['status'] == 'failed'
    assert status['reason'] == 'original task failure'


def test_incomplete_collection_retains_scheduler_category(prepared_plan, monkeypatch):
    plan, _, _ = prepared_plan
    m.run_task(plan, 0)
    monkeypatch.setattr('chronoclade.execution.accounting.task_accounting', lambda *a, **k:
        {1: {'availability': 'available', 'job_id': '123_1', 'state': 'CANCELLED', 'exit_code': '0:15'}})
    assert m.collect(plan) == 2
    status = json.loads((plan.parent / 'collection-status.json').read_text())
    assert status['failures'][0]['scheduler']['state'] == 'CANCELLED'
    assert 'Scheduler CANCELLED' in status['failures'][0]['reason']
    assert not (plan.parent.parent / 'prepared').exists()
