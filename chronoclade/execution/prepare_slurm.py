"""Bounded per-assembly Slurm preparation with checked reuse and a strict collector."""
from __future__ import annotations

import csv
from dataclasses import asdict
from datetime import datetime, timezone
import os
from pathlib import Path
import shlex
import shutil
import sys
import tempfile

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import DatasetError, load_dataset
from chronoclade.errors import WorkflowError
from chronoclade.execution.slurm import _script
from chronoclade.execution_provenance import execution_provenance
from chronoclade.pathogenwatch import content_hash
from chronoclade.prepare_recovery import DirectoryLock
from chronoclade.profile_inputs import _canonical, _csv
from chronoclade.selection_manifest import _read, _resolve
from chronoclade.stage_artifacts import reference


def _mlst_binding(records):
    from chronoclade.mlst_typing import _resources, numeric_st, PRIMARY_SCHEMES
    needed = {' '.join(r['species'].replace('_', ' ').casefold().split()) for r in records
              if not numeric_st(r.get('mlst_st'))}
    if not needed:
        return {'status': 'not_needed_supplied_st'}
    resources = _resources()
    if resources is None:
        return {'status': 'unavailable'}
    executable, root, environment, version = resources
    files = {}
    for species in sorted(needed):
        spec = PRIMARY_SCHEMES.get(species)
        if spec:
            scheme, _ = spec
            for path in sorted((root / 'pubmlst' / scheme).glob('*')):
                if path.is_file():
                    files[str(path)] = file_sha256(path)
    for path in sorted((root / 'blast').glob('mlst.fa.*')):
        if path.is_file():
            files[str(path)] = file_sha256(path)
    dependencies = {}
    for name in ('blastn', 'any2fasta', 'perl'):
        located = shutil.which(name, path=environment.get('PATH'))
        dependencies[name] = ({'path': located, 'sha256': file_sha256(Path(located))}
                              if located else None)
    return {'status': 'available', 'dependencies': dependencies, 'executable': str(executable),
            'executable_sha256': file_sha256(executable), 'version': version, 'files': files}


def _identity(source, config, records, enrich_metadata, *, hash_assemblies=True):
    from chronoclade.pasteur_profile_normalization import reference_descriptor
    return {'input_sha256': file_sha256(source),
            'typing_reference': reference_descriptor(config),
            'assemblies': ([{'sample_id': r['sample_id'], 'sha256': file_sha256(Path(r['assembly']))}
                            for r in records] if hash_assemblies else None),
            'classical_mlst': _mlst_binding(records), 'enrich_metadata': enrich_metadata,
            'execution': execution_provenance('slurm-prepare', 'per-assembly-prepare-v1')}


def write_slurm_prepare(source, *, typing_config, output, results, resources,
                        max_parallel=4, species=None, enrich_metadata=False):
    """Freeze an assembly CSV, generate workers and collect all inputs in order."""
    if type(max_parallel) is not int or max_parallel < 1:
        raise ValueError('Maximum parallel tasks must be a positive integer')
    source = Path(source).expanduser().resolve()
    config = Path(typing_config).expanduser().resolve()
    raw = _csv(source)
    records = [_canonical(row, origin='local', species=species or '', base=source.parent) for row in raw]
    ids = [r['sample_id'] for r in records]
    if len(set(ids)) != len(ids):
        raise DatasetError('Assembly CSV requires unique sample identities')
    for row in records:
        if not row.get('species') or not row.get('assembly'):
            raise DatasetError('slurm-prepare requires an assembly CSV with local assemblies and species')
        if species and row['species'].casefold() != species.casefold():
            raise DatasetError('Assembly species differs from --species')
    identity = _identity(source, config, records, enrich_metadata)
    output, results = Path(output).expanduser().resolve(), Path(results).expanduser().resolve()
    if output == results or output.is_relative_to(results) or results.is_relative_to(output):
        raise ValueError('Job and result directories must be separate, non-nested paths')
    if output.exists() and any(output.iterdir()):
        raise ValueError('Use a fresh job directory, or retry the saved array')
    if results.exists():
        raise ValueError('Use a fresh results directory, or collect the saved array')
    output.mkdir(parents=True, exist_ok=True)
    (output / 'logs').mkdir()
    (output / 'inputs').mkdir()
    tasks = []
    for index, (original, row) in enumerate(zip(raw, records, strict=True)):
        task_input = output / 'inputs' / f'{index:06d}.csv'
        copied = dict(original, sample_id=row['sample_id'], species=row['species'], assembly=row['assembly'])
        with task_input.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(copied))
            writer.writeheader()
            writer.writerow(copied)
        tasks.append({'sample_id': row['sample_id'], 'input': reference(output, task_input),
                      'output': str(output / 'runs' / f'{index:06d}')})
    plan = output / 'prepare-array.json'
    write_json(plan, {'schema': 'chronoclade.slurm-prepare', 'schema_version': 1,
                      'source': reference(output, source), 'typing_config': str(config),
                      'records': records, 'species': species, 'identity': identity, 'results': str(results),
                      'tasks': tasks, 'resources': asdict(resources), 'max_parallel': max_parallel})
    plan.with_suffix('.sha256').write_text(file_sha256(plan) + '\n')
    command = shlex.join([sys.executable, '-m', __name__, str(plan)])
    directives = [f'#SBATCH --cpus-per-task={resources.cpus}',
                  f'#SBATCH --mem={resources.memory_gb}G', f'#SBATCH --time={resources.walltime}']
    for key in ('account', 'partition'):
        if getattr(resources, key):
            directives.append(f'#SBATCH --{key}={getattr(resources, key)}')
    _script(output / 'array.slurm', directives + [f'#SBATCH --array=0-{len(tasks)-1}%{max_parallel}',
            '#SBATCH --output=' + shlex.quote(str(output / 'logs/task-%A_%a.log'))],
            'exec ' + command + ' task "${SLURM_ARRAY_TASK_ID:?}"', Path.cwd())
    _script(output / 'collect.slurm', directives + [
            '#SBATCH --output=' + shlex.quote(str(output / 'logs/collect-%j.log'))],
            'exec ' + command + ' collect', Path.cwd())
    submission_log = shlex.quote(str(output / 'submissions.tsv'))
    _script(output / 'submit.sh', [],
            'array=$(sbatch --parsable ' + shlex.quote(str(output / 'array.slurm')) + ')\n'
            'array=${array%%;*}\n[[ "$array" =~ ^[0-9]+$ ]]\n'
            f'printf "array\\t%s\\t{file_sha256(plan)}\\n" "$array" >> ' + submission_log + '\n'
            'printf "Array job: %s\\n" "$array"\n'
            'collector=$(sbatch --parsable --dependency="afterany:$array" '
            + shlex.quote(str(output / 'collect.slurm')) + ')\n'
            'collector=${collector%%;*}\n[[ "$collector" =~ ^[0-9]+$ ]]\n'
            'printf "collector\\t%s\\n" "$collector" >> ' + submission_log + '\n'
            'printf "Collector job: %s\\n" "$collector"', Path.cwd())
    (output / 'README.md').write_text(
        '# Parallel assembly preparation\n\n'
        f'{len(tasks)} original inputs; at most {max_parallel} workers. Each worker runs ordinary prepare.\n\n'
        f'Submit: `bash {output / "submit.sh"}`. Inspect `logs/` and `task-*.status.json`.\n\n'
        'Retry the saved submission to reuse checksum-verified successes. For selected zero-based indices, '
        'use `sbatch --array=0,2 array.slurm`, then submit `collect.slurm` after those jobs finish. '
        'Collection never performs typing and requires every original input. Changed inputs, references '
        'or software require fresh job and result directories.\n')
    return plan


def _load(path, index=None):
    path = Path(path).resolve()
    if path.with_suffix('.sha256').read_text().strip() != file_sha256(path):
        raise DatasetError('Prepare plan checksum changed; use the original saved plan or create a fresh plan')
    plan = _read(path)
    if plan.get('schema') != 'chronoclade.slurm-prepare' or plan.get('schema_version') != 1:
        raise DatasetError('Expected a versioned Slurm prepare plan')
    source = _resolve(path.parent, plan['source'], 'assembly CSV', external=True)
    if index is not None and not 0 <= index < len(plan['tasks']):
        raise DatasetError('Prepare task index out of range')
    current = _identity(source, Path(plan['typing_config']), plan['records'],
                        plan['identity']['enrich_metadata'], hash_assemblies=False)
    expected = dict(plan['identity'], assemblies=None)
    if current != expected:
        raise DatasetError('Preparation inputs, references or execution identity changed; create a fresh plan')
    # The sealed plan and original CSV hash protect metadata. Workers verify
    # their own assembly; collection verifies the full cohort once, avoiding
    # quadratic assembly reads when thousands of independent tasks start.
    indices = range(len(plan['records'])) if index is None else [index]
    for position in indices:
        record, saved = plan['records'][position], plan['identity']['assemblies'][position]
        if (record['sample_id'] != saved['sample_id']
                or file_sha256(Path(record['assembly'])) != saved['sha256']):
            raise DatasetError('Original assembly content changed; create a fresh plan')
    if [t['sample_id'] for t in plan['tasks']] != [r['sample_id'] for r in plan['records']]:
        raise DatasetError('Prepare task identities differ from the original cohort')
    for task in plan['tasks']:
        _resolve(path.parent, task['input'], 'task CSV', external=True)
    return path, plan


def _checked_result(path, plan, index):
    from chronoclade.prepare_stage import _validate_previous_stage
    task = plan['tasks'][index]
    receipt = _read(path.parent / f'task-{index}.json')
    if receipt.get('plan_sha256') != file_sha256(path) or receipt.get('sample_id') != task['sample_id']:
        raise DatasetError('Prepare success receipt differs from its exact task and plan')
    result = _resolve(path.parent, receipt['dataset'], 'prepared task dataset', external=True)
    if result != Path(task['output']) / 'dataset.json':
        raise DatasetError('Prepare receipt points outside its declared task output')
    _validate_previous_stage(result.parent)
    dataset = load_dataset(result, mmap_mode=None)
    if dataset.sample_ids != (task['sample_id'],) or dataset.samples[0]['role'] != 'input':
        raise DatasetError('Prepared task did not retain its exact original input identity')
    return dataset, result


def run_task(path, index):
    from chronoclade.prepare_stage import run_prepare
    from chronoclade.progress import progress_task, progress
    path, plan = _load(path, index)
    if not 0 <= index < len(plan['tasks']):
        raise DatasetError('Prepare task index out of range')
    task = plan['tasks'][index]
    with DirectoryLock(path.parent / f'.task-{index}.lock'):
        receipt = path.parent / f'task-{index}.json'
        if receipt.is_file():
            _checked_result(path, plan, index)
            progress('slurm-prepare', f'Reused checked task {index}: {task["sample_id"]}')
            return
        status = {'plan_sha256': file_sha256(path), 'sample_id': task['sample_id'],
                  'status': 'running', 'started_at': datetime.now(timezone.utc).isoformat(),
                  'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
                  'slurm_array_job_id': os.environ.get('SLURM_ARRAY_JOB_ID'), 'task_index': index}
        status_path = path.parent / f'task-{index}.status.json'
        write_json(status_path, status)
        try:
            allocation = os.environ.get('SLURM_CPUS_PER_TASK')
            if allocation and int(allocation) < plan['resources']['cpus']:
                raise DatasetError('Slurm CPU allocation is smaller than the prepared task')
            destination = Path(task['output'])
            if destination.exists():
                raise DatasetError('Unreceipted task output exists; inspect it and use a fresh plan')
            with progress_task('slurm-prepare', f'Prepare {task["sample_id"]} (task {index})'):
                result = run_prepare(
                    _resolve(path.parent, task['input'], 'task CSV', external=True), destination,
                    input_kind='assemblies', typing_config=plan['typing_config'],
                    enrich_metadata=plan['identity']['enrich_metadata'],
                    resume=destination.with_name(destination.name + '.work').exists())
            write_json(receipt, {'plan_sha256': file_sha256(path), 'sample_id': task['sample_id'],
                                 'dataset': reference(path.parent, result.dataset_manifest)})
        except Exception as error:
            from chronoclade.pathogenwatch import load_api_key
            reason = str(error)
            try:
                key = load_api_key()
            except (OSError, ValueError, RuntimeError):
                key = None
            if key:
                reason = reason.replace(key, '[redacted]')
            status.update(status='failed', reason=reason)
            raise DatasetError(reason) from None
        else:
            status.update(status='complete')
        finally:
            status['finished_at'] = datetime.now(timezone.utc).isoformat()
            write_json(status_path, status)



def _tree_hashes(root):
    result = {}
    for entry in sorted(root.rglob('*')):
        if entry.is_symlink():
            raise DatasetError('Worker evidence contains a symlink; refusing collection')
        if entry.is_file():
            result[entry.relative_to(root).as_posix()] = file_sha256(entry)
    return result


def _retain_evidence(source, target):
    expected = _tree_hashes(source)
    if target.exists():
        if target.is_symlink() or _tree_hashes(target) != expected:
            raise DatasetError('Retained worker evidence differs from its checked source')
        return
    # Incomplete copies stay outside the directory later exported as evidence.
    temporary = Path(tempfile.mkdtemp(prefix='.prepare-evidence-', dir=target.parent.parent))
    try:
        copied = temporary / 'copy'
        shutil.copytree(source, copied)
        if _tree_hashes(copied) != expected:
            raise DatasetError('Worker evidence changed during collection copy')
        os.replace(copied, target)
    finally:
        shutil.rmtree(temporary, ignore_errors=True)


def _validate_final(destination, datasets, path, plan):
    from chronoclade.prepare_stage import _validate_previous_stage
    _validate_previous_stage(destination)
    final = load_dataset(destination / 'dataset.json', mmap_mode=None)
    if (final.sample_ids != tuple(t['sample_id'] for t in plan['tasks'])
            or final.parameters.get('parallel_prepare', {}).get('plan_sha256') != file_sha256(path)
            or final.parameters.get('native_typing_reference') != plan['identity']['typing_reference']):
        raise DatasetError('Collected preparation differs from its original cohort or reference binding')
    expected_tasks = final.parameters['parallel_prepare']['tasks']
    if [t['dataset_sha256'] for t in expected_tasks] != [file_sha256(p) for _, p in datasets]:
        raise DatasetError('Collected preparation differs from its checked per-assembly sources')
    for actual, (worker, manifest) in zip(final.samples, datasets, strict=True):
        expected = worker.samples[0]
        if ({k: v for k, v in actual.items() if k != 'assembly_reference'}
                != {k: v for k, v in expected.items() if k != 'assembly_reference'}):
            raise DatasetError('Collected sample metadata differs from its checked assembly task')
        if expected.get('assembly_reference'):
            if file_sha256(destination / actual['assembly_reference']) != file_sha256(
                    manifest.parent / expected['assembly_reference']):
                raise DatasetError('Collected assembly differs from its checked assembly task')
        for matrix in worker.profiles:
            matches = [m for m in final.profiles if asdict(m.catalogue) == asdict(matrix.catalogue)]
            if len(matches) != 1 or matches[0].profile(expected['sample_id']) != matrix.profile(expected['sample_id']):
                raise DatasetError('Collected allele calls differ from their checked assembly task')
    return final

def collect(path):
    """Publish the whole original cohort, with no ST or cgMLST inference."""
    from chronoclade.prepare_stage import run_prepare
    from chronoclade.progress import progress
    path, plan = _load(path)
    with DirectoryLock(path.parent / '.collect.lock'):
        datasets, failures, statuses = [], [], {}
        for index, task in enumerate(plan['tasks']):
            try:
                datasets.append(_checked_result(path, plan, index))
            except (OSError, ValueError, KeyError, WorkflowError) as error:
                status_path = path.parent / f'task-{index}.status.json'
                status = _read(status_path) if status_path.is_file() else {}
                reason = str(error)
                statuses[index] = status
                if (status.get('plan_sha256') == file_sha256(path)
                        and status.get('sample_id') == task['sample_id'] and status.get('reason')):
                    reason = status['reason'] + '; ' + reason
                failures.append({'sample_id': task['sample_id'], 'index': index, 'reason': reason})
        if failures:
            from chronoclade.execution.accounting import task_accounting, accounting_reason
            accounting = task_accounting(path, plan, statuses, identity_key='sample_id')
            for failure in failures:
                failure['scheduler'] = accounting[failure['index']]
                failure['reason'] += '; ' + accounting_reason(failure['scheduler'])
            write_json(path.parent / 'collection-status.json', {'status': 'incomplete', 'failures': failures,
                       'completed_tasks': len(datasets), 'total_tasks': len(plan['tasks'])})
            for failure in failures:
                progress('slurm-prepare', f'Incomplete task {failure["index"]}: {failure["sample_id"]}; '
                         f'{failure["reason"]}; retry this index and collect again')
            return 2
        destination = Path(plan['results'])
        final_receipt = path.parent / 'collection-result.json'
        if final_receipt.is_file():
            saved = _read(final_receipt)
            if saved['plan_sha256'] != file_sha256(path):
                raise DatasetError('Collected output belongs to another plan')
            result = _resolve(path.parent, saved['dataset'], 'collected dataset', external=True)
            if result != destination / 'dataset.json':
                raise DatasetError('Collected receipt points outside the declared final output')
            _validate_final(result.parent, datasets, path, plan)
            progress('slurm-prepare', f'Reused complete prepared dataset: {result}')
            return 0
        merged = {'records': [], 'parameters': {'parallel_prepare': {
            'plan_sha256': file_sha256(path), 'input_count': len(plan['tasks']),
            'execution': 'checked_per_assembly_results_no_collector_inference',
            'tasks': []}}, 'provenance': [], 'conflicts': [], 'retrieval': [],
                  'lineages': [], 'crosswalk': [], 'exclusions': []}
        catalogues = {}
        mlst_audits = []
        sources = path.parent / 'collected-sources'
        if sources.is_symlink():
            raise DatasetError('Collection evidence directory is a symlink')
        sources.mkdir(exist_ok=True)
        for index, (dataset, manifest) in enumerate(datasets):
            row = dict(dataset.samples[0])
            if row.get('assembly_reference'):
                row['assembly_reference'] = str((manifest.parent / row['assembly_reference']).resolve())
            for matrix in dataset.profiles:
                catalogue = matrix.catalogue
                if row['sample_id'] in matrix.sample_ids:
                    row.update(cgmlst_scheme=catalogue.scheme_id,
                               cgmlst_scheme_version=catalogue.scheme_version,
                               cgmlst_database_version=catalogue.database_version,
                               cgmlst_database_sha256=catalogue.database_sha256,
                               cgmlst_profile={k: v for k, v in matrix.profile(row['sample_id']).items()
                                               if v is not None},
                               cgmlst_loci=list(catalogue.loci),
                               cgmlst_locus_universe_complete=catalogue.complete)
            merged['records'].append(row)
            for key in ('provenance', 'conflicts', 'retrieval', 'lineages', 'crosswalk', 'exclusions'):
                merged[key].extend(getattr(dataset, key))
            for matrix in dataset.profiles:
                value = asdict(matrix.catalogue)
                catalogues[content_hash(value)] = value
            merged['parameters']['parallel_prepare']['tasks'].append({
                'sample_id': row['sample_id'], 'dataset_sha256': file_sha256(manifest),
                'original_result': str(manifest),
                'retained_evidence': f'sources/provider/task-{index:06d}'})
            mlst_audits.extend(dict(audit, retained_evidence_directory=
                                    f'sources/provider/task-{index:06d}/classical_mlst')
                               for audit in dataset.parameters.get('classical_mlst', {}).get('samples', []))
            # Keep original native calls, ST evidence and metadata receipts portable.
            target = sources / f'task-{index:06d}'
            _retain_evidence(manifest.parent / 'sources', target)
        merged['parameters']['classical_mlst'] = {'samples': mlst_audits,
            'policy': 'checked per-assembly classical MLST results; collector inference disabled'}
        source, catalogue_path = path.parent / 'collected-input.json', path.parent / 'catalogues.json'
        write_json(source, merged)
        write_json(catalogue_path, list(catalogues.values()))
        if destination.exists():
            _validate_final(destination, datasets, path, plan)
            result_path = destination / 'dataset.json'
        else:
            result = run_prepare(source, destination, input_kind='profiles', catalogues=catalogue_path,
                                 typing_config=plan['typing_config'], _provider_sources=sources, _infer_mlst=False)
            result_path = result.dataset_manifest
            _validate_final(destination, datasets, path, plan)
        write_json(final_receipt, {'plan_sha256': file_sha256(path),
                                  'dataset': reference(path.parent, result_path)})
        write_json(path.parent / 'collection-status.json', {'status': 'complete',
                   'completed_tasks': len(datasets), 'total_tasks': len(plan['tasks'])})
        progress('slurm-prepare', f'Collected all {len(datasets)} original inputs: {result_path}')
        return 0


def main():
    import argparse
    from chronoclade.progress import configure_progress
    configure_progress()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('operation', choices=['task', 'collect'])
    parser.add_argument('index', type=int, nargs='?')
    args = parser.parse_args()
    try:
        if args.operation == 'collect':
            return collect(args.plan)
        if args.index is None:
            parser.error('task requires an index')
        run_task(args.plan, args.index)
        return 0
    except (OSError, ValueError, WorkflowError) as error:
        print(f'Error: {error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
