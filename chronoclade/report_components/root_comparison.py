"""Date-optimised root sensitivity; never changes the biological display root."""
import copy
import json
import math
from pathlib import Path
import sys

from Bio import Phylo
import numpy as np

from chronoclade.temporal_diagnostics import (
    date_regression, allele_unit_diagnostic, date_distance_point,
)


NOTE = ('The clock-compatible root minimises regression residual mean squared using '
        'collection-date interval midpoints. Dates are used to choose it, so an improved fit '
        'does not establish temporal signal. The main tree and state network retain their '
        'existing root. The primary comparison uses the same dated records, topology and '
        'nonnegative diagnostic branch lengths; signed original-branch sensitivity is reported separately.')
SOURCE = 'https://beast.community/tempest_tutorial'
CLOCKOR2_SOURCE = 'https://github.com/clockor2/clockor2/blob/main/src/features/engine/bestFittingRoot.ts'


def _nonnegative_diagnostic_copy(tree):
    """Clamp negative finite branches for a separate, auditable rooting diagnostic."""
    diagnostic = copy.copy(tree)
    clones = {}
    stack = [tree.root]
    while stack:
        original = stack.pop()
        cloned = copy.copy(original)
        cloned.clades = []
        for key, value in vars(original).items():
            if key != 'clades':
                setattr(cloned, key, copy.deepcopy(value))
        clones[original] = cloned
        stack.extend(original.clades)
    for original, cloned in clones.items():
        cloned.clades = [clones[child] for child in original.clades]
    diagnostic.root = clones[tree.root]
    negative = []
    stack = [diagnostic.root]
    while stack:
        parent = stack.pop()
        for child in parent.clades:
            length = child.branch_length
            if length is not None and math.isfinite(float(length)) and length < 0:
                negative.append(float(length))
                child.branch_length = 0.0
        stack.extend(parent.clades)
    return diagnostic, {
        'policy': 'negative finite branches clamped to zero in a diagnostic copy',
        'negative_edge_count': len(negative),
        'negative_edge_total_magnitude': float(sum(abs(value) for value in negative)),
        'negative_edge_max_magnitude': float(max(map(abs, negative), default=0.0)),
        'original_tree_modified': False,
    }


def _write_newick_deep_tree(tree, path):
    """Write with enough scoped recursion headroom for Bio.Phylo's renderer."""
    node_count = 0
    stack = [tree.root]
    while stack:
        node = stack.pop()
        node_count += 1
        stack.extend(node.clades)
    previous_limit = sys.getrecursionlimit()
    try:
        sys.setrecursionlimit(max(previous_limit, 4 * node_count + 1000))
        Phylo.write(tree, path, 'newick', format_branch_length='%.17g')
    finally:
        sys.setrecursionlimit(previous_limit)


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
        rows = [by_id[name] for name in cohort['sample_ids'] if name in by_id]
        eligible = date_regression(
            [date_distance_point(row, 0) for row in rows],
            distance_field='root_to_tip', distance_units='', label='Date eligibility')
        dated_ids = {point['sample_id'] for point in eligible['points']}
        eligible_rows = [row for row in rows if row['sample_id'] in dated_ids]
        diagnostic_tree, branch_policy = _nonnegative_diagnostic_copy(tree)
        from chronoclade.cgmlst.scale_tree import IndexedTree
        indexed_tree = IndexedTree.from_tree(tree)
        indexed_diagnostic_tree = IndexedTree.from_tree(diagnostic_tree)
        optimised, audit = optimise_clock_root(indexed_diagnostic_tree, eligible_rows)
        audit.update(
            diagnostic='minimum dated-tip OLS residual MSE over fixed-tree edges',
            eligible_sample_ids=sorted(dated_ids),
            eligible_tip_count=len(eligible_rows),
            date_filter=eligible['date_filter'],
            branch_length_policy=branch_policy,
            original_root_and_network_preserved=True,
            pairwise_distances_preserved_within_diagnostic_branch_policy=(
                audit.get('pairwise_distances_preserved') if optimised is not None else None),
            estimates_are_descriptive_distance_slopes=True,
            confidence_intervals_reported=False,
        )
        result = {'cohort_id': cohort['cohort_id'], 'audit': audit, 'note': NOTE,
                  'criterion_source': SOURCE, 'criterion_sources': [SOURCE, CLOCKOR2_SOURCE],
                  'implementation': 'ChronoClade analytical continuous-edge rooting; inspired by Clockor2',
                  'eligible_tip_count': len(eligible_rows),
                  'date_filter': eligible['date_filter'],
                  'branch_length_policy': branch_policy,
                  'rooting_diagnostic': 'Clock-compatible root minimising dated-tip residual MSE; exploratory only'}
        if optimised is not None:
            optimised = IndexedTree.from_tree(optimised)
            previous = next((d for d in analysis.get('root_to_tip', [])
                             if d.get('cohort_id') == cohort['cohort_id']), {})
            def diagnostic(current):
                value = date_regression(
                    [date_distance_point(row, current.distance(row['sample_id']))
                     for row in eligible_rows],
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
            baseline, fitted = diagnostic(indexed_diagnostic_tree), diagnostic(optimised)
            signed_baseline = diagnostic(indexed_tree) if branch_policy['negative_edge_count'] else None
            result.update(baseline_root=previous.get('root', 'Existing NJ root'),
                          baseline=baseline, clock_compatible=fitted,
                          signed_original_baseline=signed_baseline,
                          slope_change=(fitted['slope'] - baseline['slope']
                              if fitted['slope'] is not None and baseline['slope'] is not None else None))
            tree_path = directory / f"{cohort['cohort_id']}_clock_compatible.nwk"
            _write_newick_deep_tree(optimised, tree_path)
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
    from chronoclade.plotting import headless_pyplot

    plt = headless_pyplot()
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
                       'described in the <a href="' + SOURCE + '">TempEst tutorial</a>, with its '
                       'continuous-edge analytical rooting approach inspired by '
                       '<a href="' + CLOCKOR2_SOURCE + '">Clockor2</a>. ChronoClade performs this '
                       'analysis locally; it is not a TempEst or Clockor2 run.</p>')
        policy = comparison.get('branch_length_policy', {})
        content.append('<p class="muted">Rooting diagnostic: ' + escape(str(
            comparison.get('rooting_diagnostic', 'Exploratory root sensitivity'))) + '. Eligible dated tips: ' +
            escape(str(comparison.get('eligible_tip_count', 0))) + '. Branch policy: ' + escape(str(
                policy.get('policy', 'original branch lengths retained'))) +
            f"; {int(policy.get('negative_edge_count', 0))} negative edges, total magnitude " +
            f"{float(policy.get('negative_edge_total_magnitude', 0.0)):.4g}. Original tree unchanged.</p>")
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
        signed = comparison.get('signed_original_baseline')
        if signed:
            content.append('<tr><td>Signed original-branch sensitivity at existing root</td>' + ''.join(
                '<td>' + ('undefined' if signed[key] is None else f'{signed[key]:.4g}') + '</td>'
                for key in ('slope', 'r_squared', 'residual_mse')) + '</tr>')
        if comparison['slope_change'] is not None:
            content.append(f"</tbody></table></div><p>Slope change: {comparison['slope_change']:+.4g} " +
                           escape(unit) + '.</p>')
        content.append(asset(directory, {'comparison': comparison['figure_path']},
                             ('comparison',), 'Root sensitivity comparison'))
    return ''.join(content)
