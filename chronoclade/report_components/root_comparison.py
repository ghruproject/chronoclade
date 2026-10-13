"""Date-optimised root sensitivity; never changes the biological display root."""
import json
from pathlib import Path

from Bio import Phylo
import numpy as np

from chronoclade.temporal_diagnostics import date_regression, allele_unit_diagnostic


NOTE = ('The clock-compatible root minimises regression residual mean squared using '
        'collection-date interval midpoints. Dates are used to choose it, so an improved fit '
        'does not establish temporal signal. The main tree and state network retain their '
        'existing root. Both fits use the same dated records, topology and branch lengths.')
SOURCE = 'https://beast.community/tempest_tutorial'


def add_root_comparisons(analysis, records, directory):
    """Save a paired diagnostic for every saved NJ cohort, including SNP panels."""
    from chronoclade.clock_rooting import optimise_clock_root
    directory = Path(directory)
    by_id = {row['sample_id']: row for row in records}
    comparisons = []
    for cohort in analysis.get('cohorts', []):
        if not cohort.get('tree_path'):
            continue
        path = Path(cohort['tree_path'])
        tree = Phylo.read(path if path.is_absolute() else directory / path, 'newick')
        rows = [by_id[name] for name in cohort['sample_ids']]
        if len(rows) > 1500:
            optimised, audit = None, {'status': 'skipped',
                'reason': 'Exact root search is limited to 1,500 genomes; use a representative tree.'}
        else:
            optimised, audit = optimise_clock_root(tree, rows)
        result = {'cohort_id': cohort['cohort_id'], 'audit': audit, 'note': NOTE,
                  'criterion_source': SOURCE, 'implementation': 'ChronoClade; not a TempEst run'}
        if optimised is not None:
            previous = next((d for d in analysis.get('root_to_tip', [])
                             if d.get('cohort_id') == cohort['cohort_id']), {})
            def diagnostic(current):
                value = date_regression(
                    [dict(row, root_to_tip=current.distance(row['sample_id'])) for row in rows],
                    distance_field='root_to_tip', distance_units=analysis.get('distance_units', ''),
                    label='Root comparison')
                factor = previous.get('normalization_loci') or previous.get('scale_loci')
                if factor:
                    value = allele_unit_diagnostic(value, int(factor),
                        fixed_callable=previous.get('callable_denominators') == 'fixed full locus count',
                        catalogue_complete=not any(row.get('cgmlst_locus_universe_complete') is False
                                                   for row in rows))
                residuals = [p['residual'] for p in value['points'] if p.get('residual') is not None]
                value['residual_mse'] = float(np.mean(np.square(residuals))) if residuals else None
                return value
            baseline, fitted = diagnostic(tree), diagnostic(optimised)
            result.update(baseline_root=previous.get('root', 'Existing NJ root'),
                          baseline=baseline, clock_compatible=fitted,
                          slope_change=fitted['slope'] - baseline['slope'])
            tree_path = directory / f"{cohort['cohort_id']}_clock_compatible.nwk"
            Phylo.write(optimised, tree_path, 'newick', format_branch_length='%.17g')
            result['tree_path'] = str(tree_path)
            analysis.setdefault('paths', {})[f"{cohort['cohort_id']}_clock_compatible_tree"] = str(tree_path)
            result['figure_path'] = _figure(baseline, fitted,
                directory / f"{cohort['cohort_id']}_root_comparison.svg")
        evidence = directory / f"{cohort['cohort_id']}_root_comparison.json"
        evidence.write_text(json.dumps(result, indent=2) + '\n')
        result['audit_path'] = str(evidence)
        analysis.setdefault('paths', {})[f"{cohort['cohort_id']}_root_comparison"] = str(evidence)
        comparisons.append(result)
    analysis['root_comparisons'] = comparisons


def _figure(baseline, fitted, path):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True, sharey=True)
    for ax, value, title in zip(axes, (baseline, fitted), ('Existing root', 'Clock-compatible root')):
        for p in value['points']:
            x = p['midpoint_year']
            ax.errorbar(x, p['root_to_tip'],
                xerr=[[x-p['year_min']], [p['year_max']-x]], fmt='o',
                markersize=3, alpha=.45, color='#5c6570', linewidth=.6)
        years = [min(p['year_min'] for p in value['points']),
                 max(p['year_max'] for p in value['points'])]
        ax.plot(years, [value['centered_intercept'] + value['slope'] *
                       (x-value['reference_year']) for x in years], color='#17191e')
        r_squared = 'undefined' if value['r_squared'] is None else f"{value['r_squared']:.3g}"
        ax.set(title=f"{title}\nSlope {value['slope']:.3g} · R² {r_squared}",
               xlabel='Collection year')
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('Estimated allele differences from root' if baseline.get('normalization_loci')
                       else baseline['distance_units'])
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return str(path)


def root_comparison_html(comparisons, directory, asset):
    """Use the report's existing asset embedding and typography."""
    from html import escape
    if not comparisons:
        return ''
    content = ['<h4>How much does the root affect the slope?</h4>']
    for comparison in comparisons:
        content.append('<p>' + escape(comparison['note']) + '</p>')
        content.append('<p class="muted">ChronoClade implements the minimum-residual criterion '
                       'described in the <a href="' + SOURCE + '">TempEst tutorial</a>.</p>')
        if not comparison.get('clock_compatible'):
            content.append('<p>Comparison unavailable: ' + escape(str(
                comparison['audit'].get('reason', 'insufficient dated data'))) + '</p>')
            continue
        baseline, fitted = comparison['baseline'], comparison['clock_compatible']
        unit = ('estimated allele differences/year' if baseline.get('normalization_loci')
                else baseline.get('slope_units', 'distance units/year'))
        content.append('<div class="table-wrap"><table><thead><tr><th>Root</th><th>Slope (' +
            escape(unit) + ')</th><th>R²</th><th>Residual mean squared</th></tr></thead><tbody>')
        for label, diagnostic in ((comparison['baseline_root'], baseline),
                                   ('Clock-compatible (date-optimised)', fitted)):
            content.append('<tr><td>' + escape(str(label)) + '</td>' + ''.join(
                '<td>' + ('undefined' if diagnostic[key] is None else f'{diagnostic[key]:.4g}') + '</td>'
                for key in ('slope', 'r_squared', 'residual_mse')) + '</tr>')
        content.append(f"</tbody></table></div><p>Slope change: {comparison['slope_change']:+.4g} " +
                       escape(unit) + '.</p>')
        content.append(asset(directory, {'comparison': comparison['figure_path']},
                             ('comparison',), 'Root sensitivity comparison'))
    return ''.join(content)
