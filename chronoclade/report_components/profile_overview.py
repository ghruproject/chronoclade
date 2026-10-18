"""Reader-facing scope, relatives and group summaries; analytical IDs stay intact."""
import re
from collections import Counter, defaultdict
from pathlib import Path

from chronoclade.location_network.reconstruction import observed_host
from chronoclade.metadata_dates import sample_date_interval
from chronoclade.sample_labels import display_label

from .profile_summary import _dataset_metadata, _origin, _read_json, _reported, _selection_path, _table


def _bundle(directory):
    bundle = Path(directory) / 'selection_bundle'
    return bundle, _read_json(bundle / 'ensemble.json') or {}


def input_scope(data, provenance):
    """Explain the page's input denominator using explicit partition provenance."""
    scope = provenance.get('input_partition') or {}
    total, current, blocks = (scope.get(k) for k in ('total_inputs', 'block_inputs', 'blocks'))
    if not all(type(value) is int and value > 0 for value in (total, current, blocks)) or current > total:
        return ''
    if blocks == 1 and current == total:
        return f'All {total} input genomes are in this analysis block.'
    remaining = total - current
    kind = data.get('partition', {}).get('kind')
    name = ('cgLIN' if kind == 'cglin' else 'HierCC' if kind == 'hiercc'
            else 'species/ST profile' if kind == 'st'
            else 'species and compatible profile-scheme' if kind == 'profiles' else 'lineage')
    elsewhere = ('The remaining inputs are outside this block; check All groups for their assignment and status.'
                 if scope.get('assigned_inputs', total) < total else
                 f'The remaining {remaining} inputs are analysed in the other block reports, available through All groups.')
    return (f'Your original input dataset contains {total} genomes, partitioned into {blocks} {name} blocks. '
            f'This report covers {current} of those inputs. {elsewhere}')


def selection_handoff(directory, data):
    """An early count-only signpost; selection methods remain in their own section."""
    from html import escape
    bundle, ensemble = _bundle(directory)
    runs = ensemble.get('selections') or []
    if not runs:
        return ''
    path = _selection_path(bundle, runs[0])
    selection = _read_json(path) if path else None
    if not selection:
        return ''
    inputs = set(selection.get('query_ids') or [])
    contexts = set(selection.get('selected_context_ids') or [])
    pool = len(inputs | set(selection.get('available_context_ids') or []))
    selected = len(inputs | contexts)
    retained = 'all ' if data.get('coverage', {}).get('total_inputs') == len(inputs) else ''
    text = (f'{pool} genomes in this block → {selected} selected for SNP-tree analysis: '
            f'{retained}{len(inputs)} inputs + {len(contexts)} public comparisons.')
    if selected < 3:
        text += ' Tree stage unavailable: at least three selected genomes are required.'
    return (f'<p class="handoff"><strong>{escape(text)}</strong> '
            '<a href="#subsamples">See selection and downloads</a>.</p>')


def nearest_summary(data):
    """One row per input, keeping all equal minimum raw-distance matches."""
    matches = defaultdict(list)
    evidence = data.get('nearest_neighbours') or []
    if not isinstance(evidence, list):
        return ''
    for row in evidence:
        if isinstance(row, dict) and row.get('query_id'):
            matches[row['query_id']].append(row)
    labels = data.get('sample_labels') or {}
    direct = data.get('distance_method') == 'direct-dna'
    difference = 'snp_differences' if direct else 'allele_differences'
    denominator = 'shared_called_sites' if direct else 'shared_called_loci'
    units = 'SNP differences' if direct else 'allele differences'
    compared = 'sites' if direct else 'loci'
    rows = []
    for query, candidates in sorted(matches.items()):
        valid = [r for r in candidates if r.get('status') == 'matched'
                 and isinstance(r.get(difference), (int, float))]
        item = {'input': display_label(labels.get(query, query))}
        if not valid:
            item.update(distance='No comparable public profile', relatives='—', countries='—')
        else:
            minimum = min(r[difference] for r in valid)
            tied = [r for r in valid if r[difference] == minimum]
            calls = sorted({r[denominator] for r in tied if isinstance(r.get(denominator), (int, float))})
            coverage = str(calls[0]) if len(calls) == 1 else f'{calls[0]}–{calls[-1]}' if calls else 'unknown'
            accessions = sorted({display_label(labels.get(r['context_id'], r['context_id']))
                                 for r in tied if r.get('context_id')})
            shown = ', '.join(accessions[:2])
            if len(accessions) > 2:
                shown += f' (+{len(accessions) - 2} other accessions)'
            unit = units.removesuffix('s') if minimum == 1 else units
            noun = 'record' if len(tied) == 1 else 'records'
            item.update(distance=f'{minimum:g} {unit}; {coverage} compared {compared}',
                        relatives=f'{len(tied)} closest public {noun}: {shown}',
                        countries=', '.join(sorted({str(r.get('country') or 'Unknown') for r in tied})))
        rows.append(item)
    if not rows:
        return ''
    return ('<h3>Closest public relatives at a glance</h3>'
            + _table(rows, [('input', 'Input'), ('distance', 'Closest difference / coverage'),
                            ('relatives', 'Closest public records'), ('countries', 'Countries')], empty='')
            + '<p class="muted">Closest means the smallest raw difference among comparable public profiles in this block. '
              'All minimum-distance ties are retained. Expand Closest relatives for every accession and comparison.</p>')


def genetic_group_summary(directory, data):
    """Join group membership to source metadata without merging accession aliases."""
    bundle, ensemble = _bundle(directory)
    metadata = _dataset_metadata(bundle, ensemble, data)
    groups = data.get('lineage_groups', data.get('genetic_groups')) or []
    if not groups:
        return ''

    analysed_ids = {ident for group in groups if isinstance(group, dict)
                    for ident in group.get('sample_ids', [])}
    optional_fields = [field for field in ('region', 'city', 'host')
                       if any(_reported(metadata.get(ident, {}).get(field)) for ident in analysed_ids)]

    def composition(members, field):
        counts = Counter((observed_host(row.get(field)) if field == 'host' else str(row.get(field)).strip())
                         if _reported(row.get(field)) else 'Unknown' for row in members)
        return '; '.join(f'{value} {count}' for value, count in sorted(counts.items(), key=lambda p: (-p[1], p[0])))

    rows = []
    for group in groups:
        if not isinstance(group, dict):
            continue
        ids = group.get('sample_ids') or []
        members = [metadata.get(ident, {'sample_id': ident}) for ident in ids]
        roles = Counter(_origin(row) for row in members)
        years, missing = [], 0
        for row in members:
            date = sample_date_interval(row)
            if date.get('status') == 'valid':
                years.extend([int(date['start'][:4]), int(date['end'][:4])])
            else:
                missing += 1
        match = re.search(r'_group_(\d+)$', str(group.get('group_id', '')))
        label = f'Group {match.group(1)}' if match else str(group.get('group_id', 'Unlabelled'))
        period = (str(min(years)) if min(years) == max(years) else f'{min(years)}–{max(years)}') if years else 'Unknown'
        if missing:
            period += f'; {missing} undated'
        rows.append({'group': label, 'inputs': roles['input'] if not roles['other'] else 'Unknown',
                     'public': roles['context'] if not roles['other'] else 'Unknown', 'total': len(ids),
                     'countries': composition(members, 'country'),
                     **{field: composition(members, field) for field in optional_fields},
                     'years': period})
    basis = data.get('display_group_basis')
    if basis == 'reported_cglin':
        title = 'cgLIN groups at a glance'
        note = ('These groups share the reported cgLIN prefix at the displayed level. '
                'No additional distance cutoff defines membership.')
    elif basis == 'reported_hiercc':
        title = 'HierCC groups at a glance'
        note = (f'These groups use the reported HierCC assignment at {data["display_group_level"]}, '
                'within this analysis block. No profile-distance cutoff defines membership; '
                'missing or unresolved assignments remain unassigned.')
    elif basis == 'unassigned_lineage':
        title = 'Lineage status at a glance'
        note = ('All profiles share one unassigned-lineage display category because this ST-only analysis '
                'does not group by lineage. This category is not an inferred lineage or genetic cluster.')
    else:
        title = 'Genetic groups at a glance'
        note = ('Groups use the recorded direct profile-distance cutoff, not their separation on the two-dimensional plot. '
                'They are exploratory genetic groups, distinct from the selected analysis partition.')
    columns = [('group', 'Group'), ('inputs', 'Inputs'), ('public', 'Public'), ('total', 'Total'),
               ('countries', 'Countries (counts)')]
    headings = {'region': 'Regions (counts)', 'city': 'Cities (counts)', 'host': 'Hosts (counts)'}
    columns.extend((field, headings[field]) for field in optional_fields)
    columns.append(('years', 'Collection years'))
    composition_note = 'location and host' if 'host' in optional_fields else 'location'
    return (f'<h3>{title}</h3>'
            '<p>Group names match the ordination colour key. Counts combine inputs and public comparisons; '
            f'{composition_note} counts include missing metadata as Unknown.</p>'
            + _table(rows, columns, empty='')
            + f'<p class="muted">{note}</p>')


def grouping_rule(data):
    """Expose the saved clustering method and cutoff beside both linked views."""
    from html import escape
    if data.get('display_group_basis') == 'reported_cglin':
        level = data['display_group_level']
        return ('<p class="group-rule"><strong>What defines the groups?</strong> '
                f'Colours show the reported cgLIN prefixes at level {level}, within this analysis block. '
                'Genomes belong together when their first '
                f'{level} cgLIN positions match, using the same typing scheme. '
                'The codes come from the prepared typing results; no 2% distance cutoff or MDS clustering is used. '
                'Unassigned codes remain unassigned. Select a cgLIN group in either figure to highlight its genomes in both.</p>')
    if data.get('display_group_basis') == 'reported_hiercc':
        level = escape(str(data['display_group_level']))
        return ('<p class="group-rule"><strong>What defines the groups?</strong> '
                f'Colours show reported HierCC assignments at {level}, within this analysis block. '
                'Codes come from the prepared typing results; no profile-distance cutoff or MDS clustering defines membership. '
                'Missing or unresolved assignments remain unassigned. Select a HierCC group in either figure to highlight its genomes in both.</p>')
    if data.get('display_group_basis') == 'unassigned_lineage':
        return ('<p class="group-rule"><strong>What defines the groups?</strong> '
                'This ST-only analysis uses one unassigned-lineage display category. No lineage or genetic clusters '
                'are inferred from profile distances.</p>')
    cutoff = data.get('distance_threshold')
    groups = data.get('genetic_groups') or []
    if not groups or not isinstance(cutoff, (int, float)):
        return ''
    denominator = ('jointly compared nucleotide sites' if data.get('distance_method') == 'direct-dna'
                   else 'jointly called cgMLST loci')
    rule = (f'The {len(groups)} coloured genetic groups use complete-linkage clustering of direct profile distances. '
            f'Every pair within a group differs at no more than {cutoff:.6%} of {denominator}. ')
    # Percentages such as 2% read better than 2.000000%; retain small DNA cutoffs.
    rule = rule.replace(f'{cutoff:.6%}', f'{cutoff * 100:.6g}%')
    return ('<p class="group-rule"><strong>What defines the groups?</strong> ' + escape(rule)
            + 'Group numbers are local labels, not Pathogenwatch LIN codes. They are not inferred from the MDS plot '
              'and need not form single tree clades. Select a group in either figure to highlight its genomes in both.</p>')


def reported_lineage_annotations(data):
    """Show supplied nomenclature as annotations, separate from distance groups."""
    evidence = data.get('reported_lineage_annotations') or []
    if not evidence:
        return ''
    labels = data.get('sample_labels') or {}
    rows = []
    for row in evidence:
        ident = row.get('sample_id', '')
        hiercc = row.get('hiercc_codes') or {}
        codes = '; '.join(f'{level}: {code}' for level, code in sorted(hiercc.items())) if isinstance(hiercc, dict) else str(hiercc)
        rows.append({
            'sample': display_label(labels.get(ident, ident)),
            'cglin': row.get('cglin_raw') or 'Unavailable',
            'cglin_scope': row.get('cglin_scheme') or 'Unreported',
            'cglin_status': row.get('cglin_status') or 'Unreported',
            'hiercc': codes or 'Unavailable',
            'hiercc_scope': row.get('hiercc_scheme') or 'Unreported',
            'hiercc_status': row.get('hiercc_status') or 'Unreported',
        })
    return (
        '<details class="coverage-details"><summary>Optional reported cgLIN/HierCC labels</summary>'
        '<p>These labels were supplied with the profiles. They do not define the distance groups or '
        'establish cluster membership across incompatible schemes or reference versions. '
        'Unreported, provisional or unresolved labels remain unverified; complete scope and provenance '
        'are retained in the prepared dataset.</p>'
        + _table(rows, [('sample', 'Genome'), ('cglin', 'Reported cgLIN'),
                        ('cglin_scope', 'cgLIN scheme'), ('cglin_status', 'cgLIN status'),
                        ('hiercc', 'Reported HierCC'), ('hiercc_scope', 'HierCC scheme'),
                        ('hiercc_status', 'HierCC status')], empty='')
        + '</details>'
    )
