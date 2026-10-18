"""Presentation of saved context scope and reported-prefix relationships."""
from html import escape

from chronoclade.cglin import resolved_prefix
from chronoclade.typing_scopes import typing_scope


def scope_summary(parameters, audit, *, provider=None, records=(), analysis=None, block_context_count=None):
    """Render explicit scope/counting units without selecting or assigning genomes."""
    provider, analysis = provider or {}, analysis or {}
    partition = parameters.get('partition_by', audit.get('partition_by', 'lineage'))
    st_only = partition == 'st' or audit.get('policy') == 'supplied-inputs-species-st-profile-scope'
    hc = parameters.get('hiercc_level', audit.get('hiercc_level'))
    level = parameters.get('lin_level', 5)
    context = parameters.get('context_lin_level', audit.get('context_lin_level'))
    display = parameters.get('display_lin_level', 7)
    policy = provider.get('policy') or {}
    if st_only:
        scopes = 'Search: not run · Analysis: supplied profiles by species/ST · Colours: lineage unassigned.'
    elif hc:
        scopes = f'Saved search: {hc}' if provider.get('pools') else 'Search: not run in this stage'
        scopes += f' · Analysis: {hc} · Colours: {hc}.'
    else:
        if provider.get('pools'):
            scopes = f"Saved search: species/ST → LIN{policy.get('lin_level', context or level)}"
        else:
            scopes = 'Search: not run in this stage'
        scopes += f' · Context: LIN{context} input-prefix union' if context is not None else f' · Context: LIN{level} matching background'
        scopes += f' · Analysis: LIN{level} · Colours: LIN{display}.'
    counts = audit.get('context_scope_counts') or {}
    chain = []
    if provider.get('pools'):
        discovered = sum(p.get('discovered', 0) for p in provider['pools'])
        chain.append(f'{discovered:,} ST candidate records')
        if provider.get('accessible_matches') is not None:
            chain.append(f"{provider['accessible_matches']:,} prefix-matching source records")
    stages = [(counts.get('prepared_context'), 'prepared public comparisons'),
              (counts.get('scope_context'), 'within context scope')]
    if records:
        public = {r['sample_id'] for r in records if r.get('role', r.get('origin')) == 'context'}
        stages.append((len(public) if block_context_count is None else block_context_count,
                       'public comparisons in this block'))
        if analysis.get('cohorts') is not None:
            comparable = set().union(*(set(c['sample_ids']) for c in analysis['cohorts'])) if analysis['cohorts'] else set()
            stages.append((len(public & comparable), 'in NJ comparable cohorts'))
    stages = [(value, label) for value, label in stages if value is not None]
    equal = len(stages) > 1 and len({value for value, _ in stages}) == 1
    if equal and counts.get('prepared_context') is not None:
        count = stages[0][0]
        noun = 'unique public comparisons' if provider.get('context_deduplication') else 'prepared public comparisons'
        chain.append(f'{count:,} {noun}')
        included = (f'All {count:,} are included and comparable in this block.'
                    if records and analysis.get('cohorts') is not None else
                    f'All {count:,} are included in this block.' if records else
                    f'All {count:,} are within context scope.')
        count_text = ' → '.join(chain) + '. ' + included
    else:
        grouped = {}
        for value, label in stages:
            grouped.setdefault(value, []).append(label)
        separate = [f"{value:,} " + ' / '.join(labels) for value, labels in grouped.items()]
        count_text = (' → '.join(chain) + '; ' if chain else '') + '; '.join(separate)
        count_text += '.' if count_text else ''
    from chronoclade.provider_hiercc import reported_scope, DISCLOSURE
    reported = bool(parameters.get("provider_reported_hiercc")) or bool(provider.get("provider_reported_hiercc")) or any(reported_scope(r) for r in records)
    disclosure = ("<p>" + escape(DISCLOSURE) +
                  " Labels were retrieved, not assigned locally; nearest-neighbour distances are calculated separately.</p>") if reported else ""
    caveat = ('This is a background comparison scope. Prefix membership alone does not establish '
              'a tight pairwise cluster; nearest neighbours use actual distances within this scope.')
    if st_only:
        caveat = 'Only supplied profiles are analysed; no lineage assignment required. No public search or lineage groups are inferred.'
    return ('<div class="comparison-scope"><p><strong>Comparison scope.</strong> ' + escape(scopes) + '</p>'
            + ('<p>' + escape(count_text) + '</p>' if count_text else '')
            + disclosure + '<p class="muted">' + escape(caveat) + '</p></div>')


def _scope(row):
    scope = typing_scope(row, 'cglin')
    if scope is None or row.get('cglin_status') in {'failed', 'unassigned', 'missing', 'unsupported', 'conflict', 'malformed'}:
        return None
    return scope, row.get('cglin_database_version'), row.get('cglin_database_sha256')


def prefix_relationship(left, right, level):
    """Describe resolved reported prefixes; unknown suffixes are never filled."""
    a_scope, b_scope = _scope(left), _scope(right)
    if a_scope is None or b_scope is None:
        return 'Unresolved typing scope', 'Not assessed'
    if a_scope[0] != b_scope[0]:
        return 'Incompatible typing scopes', 'Not comparable'
    for a, b in zip(a_scope[1:], b_scope[1:]):
        if a and b and a != b:
            return 'Incompatible reference versions', 'Not comparable'
        if bool(a) != bool(b):
            return 'Reference compatibility unverified', 'Not assessed'
    a = next((p for d in range(10, 0, -1) if (p := resolved_prefix(left, d)) is not None), ())
    b = next((p for d in range(10, 0, -1) if (p := resolved_prefix(right, d)) is not None), ())
    if not a or not b:
        return f'Unresolved LIN{level} prefix', 'Not assessed'
    shared = 0
    for x, y in zip(a, b):
        if x != y:
            break
        shared += 1
    depth = f'≥{shared}' if shared == min(len(a), len(b)) and shared < 10 else str(shared)
    if min(len(a), len(b)) < level:
        relation = f'Unresolved LIN{level} prefix'
    else:
        relation = f'Same LIN{level} prefix' if a[:level] == b[:level] else f'Different LIN{level} prefixes'
    return relation, depth


def annotate_neighbours(neighbours, records, *, level, kind):
    """Copy display rows, leaving ranking, IDs and tie membership untouched."""
    if kind != 'cglin' or not records:
        return neighbours
    by_id = {r['sample_id']: r for r in records}
    result = []
    for row in neighbours:
        item = dict(row)
        if row.get('status') == 'matched':
            relation, depth = prefix_relationship(by_id.get(row.get('query_id'), {}), by_id.get(row.get('context_id'), {}), level)
            wording = (f'at least {depth[1:]} shared levels' if depth.startswith('≥')
                       else f'{depth} shared levels' if depth.isdigit() else '')
            item.update(display_prefix_relation=relation, shared_resolved_prefix_depth=depth,
                        lin_relationship=relation + (' · ' + wording if wording else ''))
        result.append(item)
    return result
