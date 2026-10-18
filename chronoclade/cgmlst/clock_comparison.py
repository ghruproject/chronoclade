"""Checksum-bound, matched-panel cgMLST versus corrected-SNP root diagnostics."""
from __future__ import annotations

import csv
import json
from html import escape
from pathlib import Path
from typing import Annotated

import numpy as np
import typer
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from chronoclade.artifacts import write_json
from chronoclade.clock_rooting import optimise_clock_root
from chronoclade.datasets import load_dataset
from chronoclade.errors import WorkflowError
from chronoclade.matrix_distances import BinaryDistanceEvidence
from chronoclade.report_components.root_comparison import _nonnegative_diagnostic_copy
from chronoclade.report_components.stage import write_report
from chronoclade.selection_manifest import _read, _resolve, load_selection_ensemble
from chronoclade.cgmlst.clock_evidence import (
    _dense_distance_matrix, _distance_sample_ids, _read_clock_distance_evidence,
)
from chronoclade.selections import DistanceEvidence
from chronoclade.stage_artifacts import header, load_result, reference
from chronoclade.temporal_diagnostics import (
    allele_unit_diagnostic, date_regression, date_distance_point, diagnostic_display_units,
)
from chronoclade.cgmlst.scale_tree import IndexedTree

NOTE = ('Each selection uses identical dated genomes in both panels. The clock-compatible root '
        'minimises ordinary least-squares residual MSE on a separate nonnegative diagnostic copy. '
        'Dates choose the root; improved R² is descriptive and does not establish temporal signal. '
        'Slopes have different units and must not be compared as mutation rates. No date permutations '
        'or significance tests are performed here; original SNP files remain unchanged.')



def _comparison_computed_with():
    from chronoclade.execution_provenance import execution_provenance

    return execution_provenance(
        'matched-clock-comparison', 'matched-root-diagnostics-v1',
        software_version=header('matched-clock-comparison')['software']['version'],
    )


def _report_slope_units(fitted):
    return fitted.get('display_slope_units') or fitted.get('slope_units') or 'Unavailable'


def exact_panel_nj(evidence, ids, *, expected_loci=None):
    """Reconstruct NJ from this exact selected submatrix, never prune the pool NJ."""
    if isinstance(evidence, DistanceEvidence):
        if expected_loci is None:
            raise WorkflowError('Dense distance evidence requires its prepared locus count')
        evidence = _dense_distance_matrix(evidence, ids, expected_loci)
    if len(ids) != len(set(ids)) or not set(ids).issubset(evidence.sample_ids):
        raise WorkflowError('Selected identifiers must be unique and present in full-pool distances')
    position = {name: i for i, name in enumerate(evidence.sample_ids)}
    indices = [position[name] for name in ids]
    matrix = np.asarray(evidence.arrays['distance'][np.ix_(indices, indices)])
    if not np.all(np.isfinite(matrix)):
        raise WorkflowError('Exact selected cgMLST panel contains unavailable distance pairs')
    loci = np.unique(evidence.arrays['scheme_loci'][np.ix_(indices, indices)])
    if len(loci) != 1 or loci[0] < 1:
        raise WorkflowError('Exact selected cgMLST panel needs one recorded locus universe')
    if expected_loci is not None and int(loci[0]) != expected_loci:
        raise WorkflowError('Exact selected panel locus count disagrees with its prepared catalogue')
    tree = DistanceTreeConstructor().nj(DistanceMatrix(list(ids),
        [matrix[i, :i + 1].tolist() for i in range(len(ids))]))
    return tree, int(loci[0])


def matched_dates(rows, ids):
    """Declare one date panel before either fit, using the shared >=1900 policy."""
    if [row['sample_id'] for row in rows] != ids:
        raise WorkflowError('Metadata does not preserve the exact ordered selected identifiers')
    eligibility = date_regression([date_distance_point(row, 0) for row in rows],
        distance_field='root_to_tip', distance_units='', label='Matched date eligibility')
    names = {point['sample_id'] for point in eligibility['points']}
    return [row for row in rows if row['sample_id'] in names], eligibility


def _fit(tree, rows, units, loci=None, *, catalogue_complete=True):
    tree = IndexedTree.from_tree(tree)
    value = date_regression([date_distance_point(row, tree.distance(row['sample_id']))
        for row in rows],
        distance_field='root_to_tip', distance_units=units, label='Matched root-to-tip regression')
    if loci:
        value = allele_unit_diagnostic(value, loci, fixed_callable=False, catalogue_complete=catalogue_complete)
    residuals = [point['residual'] for point in value['points'] if point.get('residual') is not None]
    value['residual_mse'] = float(np.mean(np.square(residuals))) if residuals else None
    return value


def compare_panel(tree, rows, output, name, units, loci=None, *, catalogue_complete=True):
    diagnostic, policy = _nonnegative_diagnostic_copy(tree)
    rooted, audit = optimise_clock_root(IndexedTree.from_tree(diagnostic), rows)
    value = {'branch_length_policy': policy, 'root_audit': audit,
        'original_root': _fit(tree, rows, units, loci, catalogue_complete=catalogue_complete),
        'nonnegative_original_root': _fit(diagnostic, rows, units, loci, catalogue_complete=catalogue_complete),
        'clock_compatible': _fit(rooted, rows, units, loci, catalogue_complete=catalogue_complete) if rooted is not None else None}
    Phylo.write(tree, output / f'{name}_original.nwk', 'newick', format_branch_length='%.17g')
    Phylo.write(diagnostic, output / f'{name}_nonnegative.nwk', 'newick', format_branch_length='%.17g')
    if rooted is not None:
        Phylo.write(rooted, output / f'{name}_clock_compatible.nwk', 'newick', format_branch_length='%.17g')
    return value


def _figure(run, path):
    from chronoclade.plotting import headless_pyplot
    plt = headless_pyplot()
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for row_axes, name in zip(axes, ('cgmlst', 'corrected_snp'), strict=True):
        panel = run[name]
        for ax, key, title in zip(row_axes, ('nonnegative_original_root', 'clock_compatible'),
                                 ('Existing root, nonnegative copy', 'Clock-compatible root'), strict=True):
            fit = panel[key]
            if fit is None or fit['slope'] is None:
                ax.text(.5, .5, 'Fit unavailable', ha='center', transform=ax.transAxes)
                continue
            points = fit['points']
            x = [p['midpoint_year'] for p in points]
            ax.scatter(x, [p['root_to_tip'] for p in points], s=8, alpha=.45, color='#645485')
            ends = np.array([min(x), max(x)])
            ax.plot(ends, fit['centered_intercept'] + fit['slope'] * (ends-fit['reference_year']), color='#27212e')
            r2 = 'undefined' if fit['r_squared'] is None else f"{fit['r_squared']:.3f}"
            ax.set(title=f'{name.replace("_", " ")} · {title}\nR² {r2} · n={fit["n_points"]}',
                xlabel='Collection-date midpoint year', ylabel=diagnostic_display_units(fit))
            ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _check_binding(tree_result, selection_entry, ids, rows, dataset_rows, dataset_sha):
    if (tree_result['selection_id'] != selection_entry['selection_id']
            or tree_result['source']['selection']['sha256'] != selection_entry['sha256']
            or tree_result['source']['dataset']['sha256'] != dataset_sha
            or tree_result['selected_sample_ids'] != ids):
        raise WorkflowError('Tree does not bind the exact saved selection and prepared dataset')
    if rows != [dataset_rows[ident] for ident in ids]:
        raise WorkflowError('Selected metadata differs from the checksum-bound prepared dataset')


def compare_clocks(tree_ensemble, *, output):
    source = Path(tree_ensemble).expanduser().resolve()
    output = Path(output).expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise WorkflowError('Use an empty comparison output directory to preserve saved evidence')
    index = _read(source)
    if index.get('schema') != 'chronoclade.tree-ensemble' or index.get('schema_version') != 1:
        raise WorkflowError('Expected a versioned tree ensemble')
    ensemble_path = _resolve(source.parent, index['source']['ensemble'], 'source.ensemble', external=True)
    ensemble = load_selection_ensemble(ensemble_path)
    if [run.get('selection_id') for run in index['runs']] != [s['selection_id'] for s in ensemble['selections']]:
        raise WorkflowError('Tree ensemble must retain exact saved selections once in order')
    dataset_path = _resolve(ensemble_path.parent, ensemble['source']['dataset'], 'source.dataset', external=True)
    dataset = load_dataset(dataset_path)
    dataset_rows = {row['sample_id']: dict(row) for row in dataset.samples}
    evidence_path = _resolve(ensemble_path.parent, ensemble['source']['distance_evidence'], 'distance_evidence', external=True)
    evidence = _read_clock_distance_evidence(evidence_path, ensemble["source"])
    output.mkdir(parents=True, exist_ok=True)
    result = {**header('matched-clock-comparison'), 'status': 'running',
        'computed_with': _comparison_computed_with(), 'source': {'tree_ensemble': reference(output, source),
        'selection_ensemble': reference(output, ensemble_path), 'dataset': reference(output, dataset_path),
        'distance_evidence': reference(output, evidence_path)}, 'note': NOTE, 'runs': [],
        'parameters': {'minimum_collection_year': 1900, 'date_randomisations': 0,
            'profile_tree_method': 'fresh Biopython NJ on exact selected full-pool distance submatrix',
            'snp_tree_method': 'saved corrected SNP tree with existing biological root',
            'root_objective': 'OLS residual MSE', 'slope_comparison_across_panels': False}}
    for run, entry in zip(index['runs'], ensemble['selections'], strict=True):
        if run['status'] != 'complete':
            result['runs'].append({'selection_id': run['selection_id'], 'status': run['status'],
                'reason': run.get('reason', 'Source corrected tree unavailable')})
            continue
        selection_path = _resolve(ensemble_path.parent, entry, 'selection')
        selection = _read(selection_path)
        tree_path = _resolve(source.parent, run['result'], 'tree result', external=True)
        tree_result = load_result(tree_path, 'tree')
        ids = selection['selected_sample_ids']
        rows = json.loads(_resolve(tree_path.parent, tree_result['metadata'], 'metadata').read_text())
        _check_binding(tree_result, entry, ids, rows, dataset_rows, ensemble['source']['dataset']['sha256'])
        if tree_result.get('branch_units') != 'substitutions/site':
            raise WorkflowError('Corrected SNP comparison requires declared substitutions/site branch units')
        snp = Phylo.read(_resolve(tree_path.parent, tree_result['rooted_tree'], 'rooted_tree'), 'newick')
        for field in ('corrected_tree', 'rooted_tree'):
            checked = Phylo.read(_resolve(tree_path.parent, tree_result[field], field), 'newick')
            tips = [tip.name for tip in checked.get_terminals()]
            if len(tips) != len(ids) or set(tips) != set(ids):
                raise WorkflowError('Corrected SNP tree must preserve exact selected identifiers')
        eligible_rows, date_audit = matched_dates(rows, ids)
        directory = output / run['selection_id']
        directory.mkdir()
        evidence_ids = _distance_sample_ids(evidence)
        if not evidence_ids.issubset(dataset_rows):
            raise WorkflowError('Distance evidence identifiers are outside the prepared dataset')
        if isinstance(evidence, BinaryDistanceEvidence):
            nj, loci = exact_panel_nj(evidence, ids)
            catalogues = [matrix.catalogue for matrix in dataset.profiles
                if set(ids).issubset(matrix.sample_ids) and len(matrix.catalogue.loci) == loci]
        else:
            catalogues = [matrix.catalogue for matrix in dataset.profiles
                if set(ids).issubset(matrix.sample_ids)]
            if len(catalogues) != 1:
                raise WorkflowError('Exact selected panel must bind one prepared profile locus catalogue')
            nj, loci = exact_panel_nj(evidence, ids, expected_loci=len(catalogues[0].loci))
        if len(catalogues) != 1:
            raise WorkflowError('Exact selected panel must bind one prepared profile locus catalogue')
        catalogue_complete = catalogues[0].complete
        # A fixed, date-independent baseline makes each fresh NJ root explicit.
        nj.root_with_outgroup(ids[0])
        item = {'selection_id': run['selection_id'], 'status': 'complete',
            'source_tree': reference(output, tree_path), 'source_selection': reference(output, selection_path),
            'selected_sample_ids': ids, 'eligible_sample_ids': [r['sample_id'] for r in eligible_rows],
            'eligible_tip_count': len(eligible_rows), 'date_filter': date_audit['date_filter'],
            'date_exclusions': date_audit['exclusions'], 'normalization_loci': loci,
            'profile_baseline_root_sample_id': ids[0],
            'profile_catalogue_complete': catalogue_complete,
            'cgmlst': compare_panel(nj, eligible_rows, directory, 'cgmlst', 'cgMLST mismatch fraction', loci,
                catalogue_complete=catalogue_complete),
            'corrected_snp': compare_panel(snp, eligible_rows, directory, 'corrected_snp', 'substitutions/site')}
        write_json(directory / 'comparison.json', item)
        _figure(item, directory / 'root_comparison.svg')
        result['runs'].append(item)
        write_json(output / 'comparison.json', result)
    flat = []
    html = ['<p>' + escape(NOTE) + '</p><table><thead><tr><th>Selection</th><th>Panel</th><th>Dated genomes</th><th>R² existing</th><th>R² nonnegative existing</th><th>R² clock-compatible</th><th>Clock-compatible slope and units</th><th>Negative edges (total magnitude)</th></tr></thead><tbody>']
    for run in result['runs']:
        if run['status'] != 'complete':
            html.append(f'<tr><td>{escape(run["selection_id"])}</td><td colspan="7">{escape(run["reason"])}</td></tr>')
            continue
        for name in ('cgmlst', 'corrected_snp'):
            panel = run[name]
            fitted = panel['clock_compatible'] or {}
            policy = panel['branch_length_policy']
            row = {'selection': run['selection_id'], 'panel': name, 'n_dates': run['eligible_tip_count'],
                'r2_existing': panel['original_root']['r_squared'],
                'r2_nonnegative_existing': panel['nonnegative_original_root']['r_squared'],
                'r2_clock_compatible': fitted.get('r_squared'), 'slope': fitted.get('slope'),
                'slope_units': fitted.get('slope_units'), 'negative_edge_count': policy['negative_edge_count'],
                'negative_edge_total_magnitude': policy['negative_edge_total_magnitude']}
            flat.append(row)
            cells = [row['selection'], 'cgMLST' if name == 'cgmlst' else 'Corrected SNP', row['n_dates']]
            cells += [f'{row[k]:.4f}' if row[k] is not None else 'Unavailable'
                      for k in ('r2_existing', 'r2_nonnegative_existing', 'r2_clock_compatible')]
            units = _report_slope_units(fitted)
            slope = f"{row['slope']:.4g}" if row['slope'] is not None else 'Unavailable'
            cells += [f"{slope} {units}", f"{row['negative_edge_count']} ({row['negative_edge_total_magnitude']:.5g})"]
            html.append('<tr>' + ''.join('<td>' + escape(str(v)) + '</td>' for v in cells) + '</tr>')
    html.append('</tbody></table>')
    if flat:
        with (output / 'summary.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(flat[0]))
            writer.writeheader()
            writer.writerows(flat)
    result['status'] = 'complete'
    write_json(output / 'comparison.json', result)
    report_result = dict(result, selected_sample_ids=[], temporal_assessment={
        'code': 'sensitivity', 'supported': False, 'reason': 'Matched exploratory root sensitivity; temporal significance is not assessed.'})
    write_report(output / 'report.html', title='Matched cgMLST and corrected SNP root-to-tip comparison',
        result=report_result, extra=''.join(html), report_label='MATCHED CLOCK DIAGNOSTICS',
        figure_heading='Matched root-to-tip plots',
        overview=f'{len(flat)//2} saved selections · exact matched date panels · separate distance units',
        figures=[(f'{r["selection_id"]}/root_comparison.svg', r['selection_id'] + ' matched root diagnostics')
            for r in result['runs'] if r['status'] == 'complete'])
    return result


def compare_clocks_command(
    tree_ensemble: Annotated[Path, typer.Argument(help='Checked native tree_ensemble.json')],
    output: Annotated[Path, typer.Option('--out', '--output', '-o')] = Path('matched_clock_comparison'),
):
    """Compare matched cgMLST and corrected SNP panels under a shared root objective."""
    try:
        compare_clocks(tree_ensemble, output=output)
    except (WorkflowError, ValueError, OSError, KeyError) as error:
        from chronoclade.cli import _fail
        _fail(error)
    typer.echo(f'Matched root diagnostics: {output.resolve() / "report.html"}')
