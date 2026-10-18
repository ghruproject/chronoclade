"""Explicit supplied-input analysis when no lineage partition is requested."""
import hashlib
import json

from chronoclade.context_partitions import mlst_comparison_records


def partition_profiles(records):
    """Compare species and exact profile scopes; ST/lineage labels are metadata.

    Missing ST never blocks local comparison. Conflicting or unknown reference
    namespaces are kept separate instead of translating allele numbers by name.
    Context is admitted only to blocks containing inputs.
    """
    from chronoclade.typing_scopes import typing_scope

    groups, unresolved, unmatched = {}, [], []
    for row in records:
        scope = typing_scope(row, 'cgmlst')
        species = ' '.join((row.get('species') or '').replace('_', ' ').casefold().split())
        calls = dict(row.get('cgmlst_profile') or {}, **(row.get('cgmlst_novel_alleles') or {}))
        reason = ('missing_species' if not species else
                  'missing_profile_or_profile_scope' if not scope or not any(
                      value not in (None, '', 0, '0', '-', '?') for value in calls.values()) else None)
        if reason:
            if row['role'] == 'input':
                unresolved.append({'sample_id': row['sample_id'], 'reason': reason})
            else:
                unmatched.append(row['sample_id'])
            continue
        key = (species, *scope, row.get('cgmlst_database_version') or None,
               row.get('cgmlst_database_sha256') or None,
               tuple(sorted(row.get('cgmlst_loci') or calls)))
        groups.setdefault(key, []).append(row)
    blocks = []
    for key, rows in sorted(groups.items(), key=lambda item: repr(item[0])):
        if not any(row['role'] == 'input' for row in rows):
            unmatched.extend(row['sample_id'] for row in rows)
            continue
        digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:12]
        sts = sorted({str(row['mlst_st']) for row in rows if row.get('mlst_st')})
        blocks.append({
            'block_id': 'profiles_' + digest, 'kind': 'profiles', 'level': None, 'prefix': [],
            'species': rows[0]['species'].replace('_', ' '),
            'mlst_st': sts[0] if len(sts) == 1 else 'multiple' if sts else 'unknown',
            'mlst_st_values': sts, 'mlst_scheme': None,
            'scheme': key[1], 'scheme_version': key[2],
            'database_version': key[3], 'database_sha256': key[4],
            'lineage_assignment': 'optional annotation; not used for comparison or selection',
            'records': sorted(rows, key=lambda row: row['sample_id']),
        })
    context_count = sum(row['role'] == 'context' for row in records)
    return blocks, {
        'policy': 'species-compatible-profile-scope', 'partition_by': 'profiles',
        'context_scope': 'All supplied compatible profiles; public discovery uses verified species/ST without mandatory lineage codes',
        'lin_level': None, 'context_lin_level': None, 'hiercc_level': None,
        'context_scope_counts': {'prepared_context': context_count,
                                'scope_context': context_count - len(unmatched),
                                'outside_scope_context': len(unmatched)},
        'outside_context_scope_ids': sorted(unmatched), 'unresolved_context_scope_inputs': [],
        'unresolved_inputs': unresolved, 'available_level_counts': [],
        'unmatched_context_ids': sorted(unmatched), 'ambiguous_mlst_namespace_ids': [],
        'live_context_retrieval': 'not_requested',
        'selection_group_basis': 'cgMLST distance-derived groups; optional lineage labels do not affect selection',
    }


def partition_input_profiles(records):
    """Keep species, MLST namespaces and exact profile scopes separate.

    This route does not discover public context or assign lineage codes. It is
    opt-in and rejects context records rather than silently dropping them.
    """
    if any(row['role'] != 'input' for row in records):
        raise ValueError('ST profile-only partition requires inputs only; remove prepared context or use lineage partitioning')
    prepared = mlst_comparison_records(records)
    groups, unresolved = {}, []
    for row in prepared:
        reason = None
        if row.get('_ambiguous_mlst_namespace'):
            reason = 'ambiguous_mlst_namespace'
        elif not row.get('species') or not row.get('mlst_st'):
            reason = 'missing_species_or_mlst_st'
        elif not row.get('cgmlst_scheme') or not row.get('cgmlst_scheme_version') or not any(
                value not in (None, '', 0, '0', '-', '?') for value in (row.get('cgmlst_profile') or {}).values()):
            reason = 'missing_profile_or_profile_scope'
        if reason:
            unresolved.append({'sample_id':row['sample_id'], 'reason':reason})
            continue
        key = (' '.join(row['species'].replace('_', ' ').casefold().split()),
               str(row['mlst_st']), row.get('_comparison_mlst_scheme'),
               row['cgmlst_scheme'], row['cgmlst_scheme_version'],
               row.get('cgmlst_database_version'), row.get('cgmlst_database_sha256'),
               tuple(row.get('cgmlst_loci') or []))
        groups.setdefault(key, []).append(row)
    blocks = []
    for key, rows in sorted(groups.items(), key=lambda item:repr(item[0])):
        digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:12]
        blocks.append({'block_id':'st_'+digest, 'kind':'st', 'level':None, 'prefix':[],
                       'species':rows[0]['species'].replace('_', ' '), 'mlst_st':key[1],
                       'mlst_scheme':key[2], 'scheme':key[3], 'scheme_version':key[4],
                       'database_version':key[5], 'database_sha256':key[6],
                       'lineage_assignment':'not_used; supplied-input profile analysis',
                       'records':sorted(rows, key=lambda row:row['sample_id'])})
    return blocks, {
        'policy':'supplied-inputs-species-st-profile-scope', 'partition_by':'st',
        'context_scope':'Supplied inputs only; no public context retrieval or lineage refinement',
        'mlst_namespace_policy':'Declared namespaces remain separate; missing metadata joins only an unambiguous declared namespace.',
        'ambiguous_mlst_namespace_ids':sorted(row['sample_id'] for row in prepared if row.get('_ambiguous_mlst_namespace')),
        'lin_level':None, 'context_lin_level':None, 'hiercc_level':None,
        'context_scope_counts':{'prepared_context':0, 'scope_context':0, 'outside_scope_context':0},
        'outside_context_scope_ids':[], 'unresolved_context_scope_inputs':[],
        'unresolved_inputs':unresolved, 'available_level_counts':[], 'unmatched_context_ids':[],
        'live_context_retrieval':'not_requested',
    }
