"""Navigable output guides without relocating checksum-bound evidence files."""

from html import escape
from pathlib import Path
from urllib.parse import quote

from .tables import write_report_html
from .styles import report_styles


CATEGORIES = {
    'tables': ('Tables', 'Sample metadata, nearest neighbours, group membership and exclusions.'),
    'trees': ('Trees', 'Saved tree topology and interactive tree viewers.'),
    'figures': ('Figures and ordination', 'Distance projections, country summaries and date diagnostics.'),
    'networks': ('State-change networks', 'Interactive location or host networks and their supporting counts.'),
    'subsamples': ('Subsamples', 'Selected genome sets and manifests for the tree and time modules.'),
    'methods': ('Methods and numerical evidence', 'Analysis manifests, diagnostics and reproducible distance evidence.'),
}


def _category(path):
    name = path.name.lower()
    if 'selection_bundle' in path.parts or 'subsample' in name:
        return 'subsamples'
    if 'network' in name or 'transition' in name:
        return 'networks'
    if path.suffix in {'.nwk', '.newick', '.nex'} or 'tree' in name and path.suffix == '.html':
        return 'trees'
    if path.suffix in {'.svg', '.png', '.pdf'} or 'ordination' in name or 'pcoa' in name:
        return 'figures'
    if path.suffix in {'.csv', '.tsv'}:
        return 'tables'
    return 'methods'


def _page(title, body):
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>ChronoClade · {escape(title)}</title><style>{report_styles()}'
            'li{margin:8px 0;overflow-wrap:anywhere}small{color:#606675}'
            '</style></head><body><main class="shell"><header class="identity">'
            '<div class="identity-mark">ChronoClade · OUTPUT FILES</div>'
            f'<div class="identity-copy"><h1>{escape(title)}</h1></div></header>'
            f'<section class="stage"><div class="stage-body">{body}</div></section></main></body></html>')


def write_output_guide(directory, *, blocks=None):
    """Add README and categorized link folders; never copy or move run evidence."""
    directory = Path(directory)
    files = directory / 'files'
    files.mkdir(exist_ok=True)
    if blocks is not None:
        links = []
        lines = ['# ChronoClade cgMLST analysis', '',
                 'Start with [the group summary](index.html). Each LIN block has its own report.', '']
        for block in blocks:
            label = str(block.get('label', block['block_id']))
            report = block['report']
            href = quote(str(report), safe='/')
            guide = quote(f"{block['block_id']}/files/index.html", safe='/')
            links.append(f'<li><a href="../{href}">{escape(label)}</a> · '
                         f'<a href="../{guide}">Output files</a></li>')
            lines.append(f'- [{label}]({href}) — {block.get("input_count", "?")} inputs, '
                         f'{block.get("context_count", "?")} public comparisons.')
        body = '<p><a href="../index.html">Back to group summary</a></p><ul>' + ''.join(links) + '</ul>'
        body += '<p>The cgmlst.json manifest records dataset identity, parameters and block outputs. Partitions describe how inputs were split.</p>'
        lines += ['', '## Run evidence', '',
                  '- `cgmlst.json`: run manifest, parameters and links to block outputs.',
                  '- `partitions.json`: input partitioning and unresolved inputs.',
                  '- Each block: report, categorized file browser and preserved evidence files.']
    else:
        grouped = {key: [] for key in CATEGORIES}
        for path in sorted(directory.rglob('*')):
            relative = path.relative_to(directory)
            if (not path.is_file() or 'files' in relative.parts or path.name == 'README.md'
                    or path.name in {'profile_report.html', 'index.html'}):
                continue
            # Block guides only cover their own evidence, never sibling reports.
            grouped[_category(relative)].append((relative, path.stat().st_size))
        body = '<p><a href="../profile_report.html">Back to cgMLST report</a> · <a href="../README.md">Folder guide</a></p><ul>'
        lines = ['# ChronoClade cgMLST block', '',
                 'Open [profile_report.html](profile_report.html) for the results.',
                 'Browse [output files by type](files/index.html).', '',
                 'The report combines your input genomes and public comparisons in this LIN block.',
                 'It shows dataset composition, NJ relationships and ordination, exploratory date',
                 'diagnostics, state-change networks, and subsamples for SNP-tree analysis.', '',
                 '## Output folders', '']
        for key, (title, description) in CATEGORIES.items():
            target = files / key
            target.mkdir(exist_ok=True)
            entries = []
            for relative, size in grouped[key]:
                href = quote(relative.as_posix(), safe='/')
                entries.append(f'<li><a href="../../{href}">{escape(relative.as_posix())}</a> '
                               f'<small>({size:,} bytes)</small></li>')
            write_report_html(target / 'index.html', _page(title,
                f'<p><a href="../index.html">All output types</a> · <a href="../../profile_report.html">Report</a></p>'
                f'<p>{escape(description)}</p><ul>{"".join(entries)}</ul>'))
            body += f'<li><a href="{key}/index.html">{escape(title)}</a> ({len(entries)} files)<br>{escape(description)}</li>'
            lines.append(f'- [{title}](files/{key}/index.html): {description}')
        body += '</ul><p>These folders contain navigation pages. Original evidence files stay in their recorded paths so saved manifests and subsequent stages remain valid.</p>'
        lines += ['', '## Reuse and sharing', '',
                  '- `selection_bundle/ensemble.json`: input for the tree module; includes alternative representative sets.',
                  '- `binary_distances/`: full distance evidence; large files are stored once.',
                  '- `profile_analysis.json`: report data and diagnostics.',
                  '- `numerical_analysis.json`: frozen numerical checkpoint and evidence checksums.',
                  '- `sample_labels.csv`: public display labels mapped to stable internal identities.', '',
                  'The HTML report embeds its figures and interactive views and can be opened without a server.',
                  'Keep the output folder for downloads, navigation and subsequent analysis stages.',
                  'The files/ subfolders provide navigation; evidence remains at its recorded path.',
                  'NJ branches describe profile distances. Date regressions are exploratory;',
                  'they do not establish a mutation rate or transmission.']
    write_report_html(files / 'index.html', _page('Browse output files', body))
    (directory / 'README.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
