"""Refresh profile clock diagnostics from immutable saved numerical evidence."""

import json
import time
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import load_dataset
from chronoclade.errors import WorkflowError
from chronoclade.profile_report import write_profile_report
from chronoclade.report_components.root_comparison import add_root_comparisons
from chronoclade.selection_manifest import _resolve
from chronoclade.cgmlst.workflow import _relative_artifacts, _restore_artifact_paths


def refresh_clock_reports(manifest_path):
    """Retain profiles, membership, trees and selections; recompute diagnostics only."""
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('status') != 'complete' or not isinstance(manifest.get('blocks'), list):
        raise WorkflowError('A complete cgMLST manifest is required')
    root = manifest_path.parent
    dataset_path = _resolve(root, manifest['dataset'], 'dataset', external=True)
    dataset = load_dataset(dataset_path)
    if dataset.dataset_id != manifest['dataset_id']:
        raise WorkflowError('Prepared dataset identity differs from the cgMLST run')
    records = {row['sample_id']: dict(row) for row in dataset.samples}
    audits = []
    for block in manifest['blocks']:
        directory = root / block['block_id']
        analysis_path = directory / 'profile_analysis.json'
        analysis = _restore_artifact_paths(json.loads(analysis_path.read_text()), directory)
        checkpoint = json.loads((directory / 'numerical_analysis.json').read_text())
        ids = [ident for cohort in analysis['cohorts'] for ident in cohort['sample_ids']]
        if len(ids) != len(set(ids)) or set(ids) - records.keys():
            raise WorkflowError('Saved cohorts must contain unique prepared sample identities')
        # Check each source NJ tree against the frozen numerical checkpoint.
        for cohort in analysis['cohorts']:
            if not cohort.get('tree_path'):
                continue
            path = Path(cohort['tree_path'])
            relative = path.relative_to(directory).as_posix()
            if checkpoint['artifacts'].get(relative) != file_sha256(path):
                raise WorkflowError('Saved NJ tree differs from its numerical checkpoint')
        before = file_sha256(analysis_path)
        start = time.perf_counter()
        add_root_comparisons(analysis, [records[ident] for ident in ids], directory)
        elapsed = time.perf_counter() - start
        def is_input(ident):
            row = records[ident]
            return row.get('role', row.get('origin')) in {'local', 'query', 'focal', 'input'}
        groups = {'queries': [records[i] for i in ids if is_input(i)],
                  'context': [records[i] for i in ids if not is_input(i)]}
        provenance = {
            'source_dataset': dataset.dataset_id,
            'input_partition': analysis.get('input_partition'),
            'profile_context': dataset.parameters.get('profile_context'),
            'bounded_public_context_count': block['context_count'],
            'coverage': {name: {'total': len(rows), 'profiles_available': len(rows)}
                         for name, rows in groups.items()},
        }
        report = write_profile_report(analysis, directory=directory, provenance=provenance,
            report_label='cgMLST profile report', group_index_href='../index.html')
        write_json(analysis_path, _relative_artifacts(analysis, directory))
        block['report_sha256'] = file_sha256(report)
        audits.append({'block_id': block['block_id'], 'sample_count': len(ids),
            'diagnostic_seconds': elapsed, 'source_analysis_sha256': before,
            'refreshed_analysis_sha256': file_sha256(analysis_path),
            'root_comparisons': [{'cohort_id': r['cohort_id'], 'audit': r['audit']}
                                 for r in analysis['root_comparisons']]})
    manifest['clock_diagnostic_refresh'] = {'policy': 'saved_NJ_continuous_edge_MSE_diagnostic',
        'blocks': audits}
    write_json(manifest_path, manifest)
    return manifest


def refresh_command(
    manifest: Annotated[Path, typer.Argument(help='Completed cgmlst.json; refreshes reports in place')],
):
    """Recompute clock-root diagnostics without resampling or rebuilding NJ trees."""
    try:
        result = refresh_clock_reports(manifest)
    except (OSError, ValueError, KeyError, WorkflowError) as error:
        Console().print(f'Error: {error}', style='bold red', markup=False)
        raise typer.Exit(2) from error
    for row in result['clock_diagnostic_refresh']['blocks']:
        Console().print(f"{row['block_id']}: {row['sample_count']} profiles; diagnostics {row['diagnostic_seconds']:.2f}s")
