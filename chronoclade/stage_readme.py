"""Small shared navigation guides for prepared, tree and time output roots."""

from pathlib import Path
import shlex


def write_stage_readme(directory, stage, *, report=None, ensemble=False, readiness=None,
                       published_directory=None, completed_runs=None):
    directory = Path(directory)
    published = Path(published_directory or directory).resolve()
    manifest = f'{stage}_ensemble.json' if ensemble else f'{stage}.json'
    lines = [f'# ChronoClade {stage} output', '']
    if report:
        lines += [f'Open [the report]({report}) to review this stage.', '']
    lines += [f'`{manifest}` records the saved stage result and referenced evidence.',
              'Keep referenced files together; do not edit checksummed scientific evidence.', '']
    if stage == 'prepare':
        lines += ['- `dataset.json`: portable profiles, typing and metadata for cgMLST.',
                  '- `readiness.json`: per-sample profile, lineage and metadata readiness.',
                  '- Inspection CSVs: readable copies of the saved evidence.',
                  '- `sources/`: frozen input/provider records with provenance.', '',
                  f'Readiness: **{readiness or "see readiness.json"}**. A completed preparation',
                  'can retain incomplete typing; inspect the reasons before analysis.', '']
        command = ['chronoclade', 'cgmlst', str(published / 'dataset.json'),
                   '--out', str(published.parent / 'cgmlst')]
        lines += ['Next, analyse the prepared profiles. Add `--fetch-context` when public',
                  'comparison profiles must be retrieved; assemblies are acquired later.', '']
    elif stage == 'tree':
        lines += ['- `jobs/`: fingerprinted inference, alignment, recombination and temporal evidence.',
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
                      'its linked audit or log, then rerun the original tree command.', '']
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
    lines += ['## Reuse and recovery', '',
              ('Preparation never replaces an existing output. Choose a new output for revised inputs.'
               if stage == 'prepare' else
               'Repeat the original command to reuse compatible saved results. Changed sources or parameters are validated; use a new output for a distinct analysis.'),
              'Inspect saved readiness/assessment and command logs for missing references,',
              'assemblies or tools. A successful command does not establish biological validity.', '']
    (directory / 'README.md').write_text('\n'.join(lines), encoding='utf-8')
