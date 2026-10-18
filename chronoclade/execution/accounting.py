"""Optional bounded scheduler evidence; never authority for scientific completeness."""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

from chronoclade.artifacts import file_sha256

_STATES = frozenset({'PENDING', 'RUNNING', 'COMPLETED', 'CANCELLED', 'TIMEOUT',
    'OUT_OF_MEMORY', 'FAILED', 'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE',
    'REQUEUED', 'REVOKED', 'SUSPENDED'})
_ID = re.compile(r'[0-9]+(?:_[0-9]+)?\Z')


def task_accounting(path, plan, statuses, *, identity_key):
    """Return safe accounting for requested task indices, with explicit fallbacks."""
    path = Path(path)
    digest = file_sha256(path)
    arrays = []
    ledger = path.parent / 'submissions.tsv'
    if ledger.is_file() and ledger.stat().st_size <= 262144:
        # Legacy two-column records have no plan binding and cannot be trusted.
        for line in ledger.read_text().splitlines()[-1000:]:
            fields = line.split('\t')
            if len(fields) == 3 and fields[0] == 'array' and fields[2] == digest and fields[1].isdigit():
                arrays.append(fields[1])
    jobs, result = {}, {}
    for index, status in statuses.items():
        task = plan['tasks'][index]
        valid = (status.get('plan_sha256') == digest
                 and status.get(identity_key) == task[identity_key])
        job = str(status.get('slurm_job_id', '')) if valid else ''
        if not _ID.fullmatch(job) and arrays:
            job = f'{arrays[-1]}_{index}'
        if _ID.fullmatch(job):
            jobs[index] = job
            result[index] = {'availability': 'unavailable', 'job_id': job,
                             'reason': 'Accounting unavailable or delayed'}
        else:
            result[index] = {'availability': 'unavailable',
                             'reason': 'No plan-bound scheduler job ID'}
    if not jobs:
        return result
    # One bounded query, even if a large cohort loses all its workers.
    selected = list(jobs.items())[:1000]
    try:
        process = subprocess.run(['sacct', '-n', '-P', '-j', ','.join(dict.fromkeys(j for _, j in selected)),
            '--format=JobIDRaw,JobID,State,ExitCode'], capture_output=True, text=True, timeout=5,
            check=False)
    except (OSError, subprocess.TimeoutExpired):
        return result
    if process.returncode or len(process.stdout) > 262144:
        return result
    rows = {}
    for line in process.stdout.splitlines():
        fields = line.split('|')
        if len(fields) < 4:
            continue
        raw_job, display_job, raw_state, exit_code = (v.strip() for v in fields[:4])
        state = raw_state.split()[0].rstrip('+') if raw_state else ''
        if state in _STATES and re.fullmatch(r'[0-9]+:[0-9]+', exit_code):
            for job in (raw_job, display_job):
                if _ID.fullmatch(job):
                    rows[job] = {'availability': 'available', 'job_id': job, 'state': state,
                                 'exit_code': exit_code}
    for index, job in selected:
        if job in rows:
            result[index] = rows[job]
    return result


def accounting_reason(value):
    if value['availability'] == 'available':
        return f"Scheduler {value['state']} (exit {value['exit_code']}; job {value['job_id']})"
    return value['reason']
