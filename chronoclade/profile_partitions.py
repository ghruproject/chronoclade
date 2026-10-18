"""Explicit supplied-input analysis when no lineage partition is requested."""
import hashlib
import json

from chronoclade.context_partitions import mlst_comparison_records


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
