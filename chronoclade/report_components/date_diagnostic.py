"""Draw saved regression evidence without estimating a root, slope or dates."""

import math
from pathlib import Path

from chronoclade.temporal_diagnostics import diagnostic_display_units


def _saved_title(diagnostic):
    if 'arbitrary' in str(diagnostic.get('interpretation', '')).casefold():
        return 'Exploratory date diagnostic · arbitrary tip root'
    if diagnostic.get('cohort_root'):
        return 'Exploratory full NJ date diagnostic · declared input tip root'
    return 'Exploratory NJ date diagnostic · saved root'


def write_saved_date_figure(diagnostic, path, *, title=None):
    """Render stored points/coefficients; missing points leave existing assets intact."""
    points = diagnostic.get('points')
    if not isinstance(points, list):
        return None
    field = diagnostic.get('distance_field', 'root_to_tip')
    for point in points:
        values = [point.get(key) for key in ('year_min', 'year_max', field)]
        if not all(isinstance(v, (int, float)) and math.isfinite(v) for v in values):
            return None
        if point['year_min'] > point['year_max']:
            return None
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    try:
        for point in points:
            midpoint = (point['year_min'] + point['year_max']) / 2
            ax.errorbar(midpoint, point[field],
                        xerr=[[midpoint-point['year_min']], [point['year_max']-midpoint]],
                        fmt='o', alpha=0.4, markersize=3)
        if points and diagnostic.get('slope') is not None:
            years = [min(p['year_min'] for p in points), max(p['year_max'] for p in points)]
            slope = diagnostic['slope']
            if diagnostic.get('centered_intercept') is not None and diagnostic.get('reference_year') is not None:
                fitted = [diagnostic['centered_intercept'] + slope*(year-diagnostic['reference_year']) for year in years]
            elif diagnostic.get('intercept') is not None:
                fitted = [diagnostic['intercept'] + slope*year for year in years]
            else:
                fitted = None
            if fitted is not None:
                ax.plot(years, fitted, color='#444', label='Saved regression')
        if not points:
            ax.text(0.5, 0.5, 'No valid collection dates available', transform=ax.transAxes, ha='center')
        ax.set(xlabel='Collection year (date interval shown)',
               ylabel=diagnostic_display_units(diagnostic) or 'Genetic distance',
               title=title or _saved_title(diagnostic))
        fig.tight_layout()
        fig.savefig(path)
    finally:
        plt.close(fig)
    return str(path)


def redraw_saved_date_figures(analysis, directory):
    """Refresh only existing report figures, retaining every numerical field."""
    directory = Path(directory).resolve()
    diagnostics = {row.get('cohort_id'): row for row in analysis.get('root_to_tip', [])}
    audit = []
    for cohort in analysis.get('cohorts', []):
        ident = cohort.get('cohort_id')
        row = diagnostics.get(ident)
        raw = cohort.get('root_to_tip_figure')
        if row is None or not raw:
            audit.append({'cohort_id': ident, 'status': 'skipped_no_saved_diagnostic_or_figure'})
            continue
        path = Path(raw)
        path = (path if path.is_absolute() else directory/path).resolve()
        if not path.is_relative_to(directory):
            raise ValueError('Saved date figure must remain within its report block')
        if path.suffix.lower() not in {'.svg', '.png'}:
            audit.append({'cohort_id': ident, 'status': 'skipped_unsupported_figure_format'})
            continue
        written = write_saved_date_figure(row, path)
        audit.append({'cohort_id': ident, 'status': 'redrawn_saved_evidence' if written else 'skipped_missing_or_invalid_saved_points',
                      'figure': path.relative_to(directory).as_posix()})
    return audit
