"""Country counts shared by dense and large profile reports."""
from collections import Counter
from pathlib import Path
import re


def write_country_figure(output, cohort_id, records):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if not re.fullmatch(r'[A-Za-z0-9_-]+', cohort_id):
        raise ValueError('Country figure requires a safe cohort identifier')
    records = list(records)
    counts = Counter(str(row.get('country') or 'Unknown') for row in records)
    inputs = Counter(str(row.get('country') or 'Unknown') for row in records
                     if row.get('role', row.get('origin')) in {'input', 'local', 'query', 'focal'})
    fig, ax = plt.subplots(figsize=(10, max(2.8, len(counts) * .34 + 1.1)))
    keys = sorted(counts, key=lambda key: (-counts[key], key))
    positions = list(range(len(keys)))
    local = [inputs[key] for key in keys]
    public = [counts[key] - inputs[key] for key in keys]
    ax.barh(positions, local, height=.56, color='#2166ac', label='Input genomes')
    ax.barh(positions, public, left=local, height=.56, color='#666666',
            label='Public comparisons')
    ax.set_yticks(positions, keys, fontsize=10, color='#17191f')
    total = sum(counts.values())
    largest = max(counts.values(), default=1)
    for index, key in enumerate(keys):
        ax.text(counts[key] + largest * .015, index,
                f'{counts[key]:,}  ({counts[key] / total:.1%})',
                va='center', fontsize=9, color='#5f6470')
    ax.invert_yaxis()
    ax.set_xlim(0, largest * 1.24)
    ax.set_xlabel('Genomes', fontsize=10, color='#5f6470', labelpad=10)
    ax.tick_params(axis='y', length=0, pad=9)
    ax.tick_params(axis='x', length=0, labelsize=9, colors='#5f6470')
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color='#e5e7ec', linewidth=.6)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if keys:
        ax.legend(loc='lower left', bbox_to_anchor=(0, 1.02), ncol=2,
                  frameon=False, fontsize=9, handlelength=1.2, columnspacing=2)
    else:
        ax.text(.5, .5, 'No country records available', transform=ax.transAxes,
                ha='center', color='#5f6470')
    fig.tight_layout(pad=1.5)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    path = output / f'{cohort_id}_countries.svg'
    fig.savefig(path)
    plt.close(fig)
    return str(path)
