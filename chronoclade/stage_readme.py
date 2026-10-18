"""Small shared navigation guides for prepared, tree and time output roots."""

from pathlib import Path
import json
import shlex
from urllib.parse import quote


def _ensemble_report_links(directory, stage):
    """Find completed reports through the saved manifests, not assumed layouts."""
    directory = Path(directory).resolve()
    try:
        ensemble = json.loads((directory / f'{stage}_ensemble.json').read_text())
    except (OSError, ValueError):
        return []
    links = []
    for run in ensemble.get('runs', []):
        raw = (run.get('result') or {}).get('path')
        if run.get('status') != 'complete' or not raw:
            continue
        manifest = (directory / raw).resolve()
        if not manifest.is_relative_to(directory):
            continue
        try:
            result = json.loads(manifest.read_text())
            report = (manifest.parent / result['report']['path']).resolve()
        except (OSError, ValueError, KeyError, TypeError):
            continue
        if report.is_relative_to(directory) and report.is_file():
            href = quote(report.relative_to(directory).as_posix(), safe='/')
            label = str(run.get('selection_id') or 'Completed selection')
            links.append(f'- [{label} report]({href})')
    return links


def write_stage_readme(directory, stage, *, report=None, ensemble=False, readiness=None,
                       published_directory=None, completed_runs=None, readiness_summary=None,
                       executor_plan=None):
    directory = Path(directory)
    published = Path(published_directory or directory).resolve()
    manifest = f'{stage}_ensemble.json' if ensemble else f'{stage}.json'
    lines = [f'# ChronoClade {stage} output', '']
    if report:
        lines += [f'Open [the report]({report}) to review this stage.', '']
    if ensemble:
        reports = _ensemble_report_links(directory, stage)
        if reports:
            lines += ['Completed selection reports:', '', *reports, '']
    lines += [f'`{manifest}` records the saved stage result and referenced evidence.',
              'Keep referenced files together; do not edit checksummed scientific evidence.', '']
    if stage == 'prepare':
        lines += ['- `dataset.json`: portable profiles, typing and metadata for cgMLST.',
                  '- `readiness.json`: per-sample profile, lineage and metadata readiness.',
                  '- Inspection CSVs: readable copies of the saved evidence.',
                  '- `sources/`: frozen input/provider records with provenance.', '',
                  f'Readiness: **{readiness or "see readiness.json"}**. A completed preparation',
                  'can retain incomplete typing; inspect the reasons before analysis.', '']
        if readiness_summary:
            lines += [readiness_summary, '']
        command = ['chronoclade', 'cgmlst', str(published / 'dataset.json'),
                   '--out', str(published.parent / 'cgmlst')]
        lines += ['Next, analyse the prepared profiles. Add `--fetch-context` when public',
                  'comparison profiles must be retrieved; assemblies are acquired later.', '']
    elif stage == 'tree':
        evidence_directory = '`runs/*/jobs/`' if ensemble else '`jobs/`'
        lines += [f'- {evidence_directory}: fingerprinted inference, alignment, recombination and temporal evidence.',
                  '- Assembly acquisition retains the exact selected identities; failures do not choose replacements.',
                  '- Inspect the temporal assessment before dating. Unsupported and not-assessed',
                  '  are different from computational failure.', '']
        if ensemble:
            lines += ['- `runs/`: separate tree results for each saved selection.',
                      '- Inspect failed selections in the ensemble manifest before interpreting agreement.', '']
        if ensemble and completed_runs == 0:
            command = None
            lines += ['No tree selections completed. Fix the failures before running time.',
                      'Open `tree_ensemble.json` and inspect each failed selection reason and',
                      ('its linked audit or log, then use the saved Slurm submission below.'
                       if executor_plan else
                       'its linked audit or log, then rerun the original tree command.'), '']
        else:
            command = ['chronoclade', 'time', str(published / manifest),
                       '--out', str(published.parent / 'time')]
            lines += ['Next, date the saved tree or ensemble when its assessment supports dating.', '']
    elif stage == 'time':
        lines += ['- Dated trees and node intervals are linked from the manifest when supported.',
                  '- Unsupported dating produces an assessment report without an endorsed dated tree.',
                  '- Node-date intervals are conditional model uncertainty. Ranges between',
                  '  selections describe sampling sensitivity, not confidence intervals.', '']
        if ensemble:
            lines += ['- `runs/`: individual time results.',
                      '- `comparisons/`: comparison report and exact shared-anchor evidence.', '']
        command = None
    else:
        raise ValueError(f'Unsupported stage guide: {stage}')
    if command:
        lines += ['```bash', shlex.join(command), '```', '',
                  'Run in the configured environment (or prefix with `pixi run`). These paths',
                  'refer to the current output location; update them if you move the folder.', '']
    recovery = (
        'Preparation never replaces an existing output. Choose a new output for revised inputs.'
        if stage == 'prepare' else
        'Repeat the original command to reuse compatible saved results. Changed sources or parameters are validated; use a new output for a distinct analysis.'
    )
    if executor_plan:
        plan = Path(executor_plan).resolve()
        recovery = (
            f'This ensemble was collected from the saved Slurm plan `{plan}`. '
            f'Inspect task and collector logs in `{plan.parent / "logs"}` and failed runs in the ensemble manifest. '
            'Resubmit the saved array and dependent collector to reuse compatible completed tasks:'
        )
    lines += ['## Reuse and recovery', '', recovery]
    if executor_plan:
        lines += ['', '```bash', shlex.join(['bash', str(plan.parent / 'submit.sh')]), '```', '',
                  'If all task jobs have finished and only collection failed, collect the saved task results:', '',
                  '```bash', shlex.join(['sbatch', str(plan.parent / 'collect.slurm')]), '```', '',
                  'To retry selected failed indices, replace `0,2` with the zero-based indices from the saved plan. '
                  'Wait for existing jobs to finish before retrying the same tasks:', '',
                  '```bash',
                  'retry_job=$(' + shlex.join(['sbatch', '--parsable', '--array=0,2',
                                             str(plan.parent / 'array.slurm')]) + ')',
                  'retry_job=${retry_job%%;*}',
                  'sbatch --dependency="afterany:$retry_job" ' + shlex.quote(str(plan.parent / 'collect.slurm')),
                  '```', '',
                  'Keep the saved plan, task receipts and result folders together. '
                  'Use a fresh job and result directory for changed inputs or options.']
    lines += [
              'Inspect saved readiness/assessment and command logs for missing references,',
              'assemblies or tools. A successful command does not establish biological validity.', '']
    (directory / 'README.md').write_text('\n'.join(lines), encoding='utf-8')
