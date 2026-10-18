"""Refresh report views and clock diagnostics from saved numerical evidence."""

import json
import sys
import time
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import load_dataset
from chronoclade.errors import WorkflowError
from chronoclade.execution_provenance import execution_provenance
from chronoclade.profile_report import write_profile_report
from chronoclade.report_components.root_comparison import add_root_comparisons
from chronoclade.context_partitions import profile_records
from chronoclade.report_components.lineage_groups import render_lineage_groups
from chronoclade.report_components.date_diagnostic import redraw_saved_date_figures
from chronoclade.report_components.output_guide import write_output_guide
from chronoclade.selection_manifest import _resolve
from chronoclade.cgmlst.workflow import _relative_artifacts, _restore_artifact_paths
from chronoclade.cgmlst.report import write_partition_index


def refresh_clock_reports(manifest_path, *, recompute_diagnostics=False):
    """Repackage saved evidence; diagnostic recomputation is an explicit separate action."""
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('status') != 'complete' or not isinstance(manifest.get('blocks'), list):
        raise WorkflowError('A complete cgMLST manifest is required')
    root = manifest_path.parent
    dataset_path = _resolve(root, manifest['dataset'], 'dataset', external=True)
    dataset = load_dataset(dataset_path)
    if dataset.dataset_id != manifest['dataset_id']:
        raise WorkflowError('Prepared dataset identity differs from the cgMLST run')
    records = {row['sample_id']: row for row in profile_records(dataset, source_root=dataset_path)}
    rendered_with = execution_provenance("cgmlst-report", "saved-evidence-render-v1")
    diagnostics_with = execution_provenance("cgmlst-clock-diagnostics", "continuous-edge-MSE-v1") if recompute_diagnostics else None
    partition_path = _resolve(root, manifest['partitions'], 'partitions')
    partition = json.loads(partition_path.read_text())
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
            for key in ('tree_path', 'pcoa_csv', 'direct_pcoa_csv'):
                if not cohort.get(key):
                    continue
                path = Path(cohort[key])
                relative = path.relative_to(directory).as_posix()
                if checkpoint['artifacts'].get(relative) != file_sha256(path):
                    raise WorkflowError(f'Saved {key} differs from its numerical checkpoint')
        before = file_sha256(analysis_path)
        start = time.perf_counter()
        sys.setrecursionlimit(max(sys.getrecursionlimit(), 4 * len(ids) + 100))
        if recompute_diagnostics:
            add_root_comparisons(analysis, [records[ident] for ident in ids], directory)
            analysis["diagnostics_computed_with"] = diagnostics_with
        diagnostic_seconds = time.perf_counter() - start if recompute_diagnostics else 0.0
        date_figures = redraw_saved_date_figures(analysis, directory)
        analysis["rendered_with"] = rendered_with
        render_lineage_groups(analysis, [records[ident] for ident in ids], directory,
                              reuse_coordinates=True,
                              level=manifest["parameters"].get(
                                  "display_lin_level", analysis.get("display_group_level", 7)))
        elapsed = time.perf_counter() - start
        def is_input(ident):
            row = records[ident]
            return row.get('role', row.get('origin')) in {'local', 'query', 'focal', 'input'}
        groups = {'queries': [records[i] for i in ids if is_input(i)],
                  'context': [records[i] for i in ids if not is_input(i)]}
        provenance = {
            'relatedness_parameters': manifest['parameters'],
            'relatedness_audit': partition,
            'source_dataset': dataset.dataset_id,
            'input_partition': analysis.get('input_partition'),
            'profile_context': dataset.parameters.get('profile_context'),
            'bounded_public_context_count': block['context_count'],
            'coverage': {name: {'total': len(rows), 'profiles_available': len(rows)}
                         for name, rows in groups.items()},
        }
        report = write_profile_report(analysis, directory=directory, provenance=provenance,
            report_label='cgMLST profile report', group_index_href='../index.html',
            relatedness_records=groups['queries'] + groups['context'])
        write_json(analysis_path, _relative_artifacts(analysis, directory))
        block['report_sha256'] = file_sha256(report)
        audits.append({'block_id': block['block_id'], 'sample_count': len(ids),
            'diagnostic_seconds': diagnostic_seconds, 'refresh_seconds': elapsed, 'source_analysis_sha256': before,
            'date_figure_presentation': date_figures,
            'refreshed_analysis_sha256': file_sha256(analysis_path),
            'root_comparisons': [{'cohort_id': r['cohort_id'], 'audit': r['audit']}
                                 for r in analysis.get('root_comparisons', [])]})
    manifest['rendered_with'] = rendered_with
    manifest['clock_diagnostic_refresh'] = {
        'policy': 'saved_NJ_continuous_edge_MSE_diagnostic' if recompute_diagnostics else 'saved_evidence_presentation',
        'recomputed': recompute_diagnostics, 'computed_with': diagnostics_with, 'blocks': audits}

    index = write_partition_index(manifest['blocks'], partition,
        records_count=len(records), lin_level=manifest['parameters']['lin_level'],
        hiercc_level=manifest['parameters'].get('hiercc_level'), status=manifest['status'], output=root,
        display_lin_level=manifest['parameters'].get('display_lin_level', 7),
        profile_context=dataset.parameters.get('profile_context'))
    manifest['report']['sha256'] = file_sha256(index)
    write_json(manifest_path, manifest)
    write_output_guide(root, blocks=manifest['blocks'])
    return manifest


def refresh_command(
    manifest: Annotated[Path, typer.Argument(help='Completed cgmlst.json; refreshes reports in place')],
    recompute_diagnostics: Annotated[bool, typer.Option(help="Explicitly recompute root/date diagnostics; default preserves saved numerical results")] = False,
):
    """Repackage saved report evidence, with optional diagnostic recomputation."""
    try:
        result = refresh_clock_reports(manifest, recompute_diagnostics=recompute_diagnostics)
    except (OSError, ValueError, KeyError, WorkflowError) as error:
        Console().print(f'Error: {error}', style='bold red', markup=False)
        raise typer.Exit(2) from error
    for row in result['clock_diagnostic_refresh']['blocks']:
        Console().print(f"{row['block_id']}: {row['sample_count']} profiles; report refreshed; diagnostics {'recomputed' if recompute_diagnostics else 'preserved'}")
