"""Public count-only planning before detailed metadata or profile acquisition."""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from html import escape
from pathlib import Path

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import DatasetError, load_dataset
from chronoclade.mlst_typing import numeric_st
from chronoclade.pathogenwatch import PathogenwatchClient, content_hash, load_catalogue
from chronoclade.progress import progress


def plan_context(dataset_manifest, *, client=None, max_candidates=None, cache_root=None):
    """Preview every primary species/ST pool without genome details or exports.

    Counts are source records before profile compatibility and identity deduplication.
    A valid completed catalogue can supply its historical count without a new search.
    """
    from chronoclade.profile_context_provider import dataset_records

    if max_candidates is not None and (type(max_candidates) is not int or max_candidates < 1):
        raise DatasetError('Maximum context candidates must be a positive integer')
    source = Path(dataset_manifest).expanduser().resolve()
    if source.is_dir():
        source = source/'dataset.json'
    dataset = load_dataset(source, mmap_mode=None)
    rows = [r for r in dataset_records(dataset, source_root=source) if r['role'] == 'input']
    client = client or PathogenwatchClient()
    scopes = sorted({(r['species'], str(r.get('mlst_st') or '')) for r in rows})
    species_inputs = Counter(r['species'] for r in rows)
    species_candidates = Counter()
    organisms, pools = None, []
    for index, (species, st) in enumerate(scopes):
        pool = {'index': index, 'species': species, 'st': st, 'search_field': 'mlst',
                'candidate_count': None, 'status': 'unavailable'}
        pools.append(pool)
        if not numeric_st(st):
            pool['reason'] = 'missing_primary_mlst_st'
            continue
        relevant = [r for r in rows if r['species'] == species and str(r.get('mlst_st') or '') == st]
        if species.casefold() == 'escherichia coli' and any(
            str(r.get('mlst_scheme') or '').casefold() in
            {'ecoli', 'ecoli_2', 'ecoli_pasteur', 'pasteur', 'pubmlst:ecoli',
             'escherichia.pasteur8genemlst'} for r in relevant
        ):
            pool['reason'] = 'alternative_mlst_namespace_not_primary_achtman'
            continue
        cached = Path(cache_root)/f'public_{index}.json' if cache_root else None
        if cached and cached.exists():
            envelope = load_catalogue(cached)
            provenance = envelope['provenance']
            query = provenance['query']
            organism = provenance['organism_capabilities']
            if (provenance.get('species', '').casefold() != species.casefold()
                    or str((query.get('mlst') or [None])[0]) != st
                    or str(query.get('organismId')) != str(organism.get('organismId'))
                    or provenance.get('base_url', '').rstrip('/') != client.base_url.rstrip('/')):
                raise DatasetError(f'Frozen context catalogue {index} does not match its scope')
            count = provenance['expected_count']
            pool.update(count_source='completed_frozen_catalogue',
                        captured_at=provenance['retrieved_at'],
                        catalogue_sha256=file_sha256(cached))
        else:
            if organisms is None:
                organisms = client.supported_organisms()
            organism = next((r for r in organisms
                             if r.get('fullName', '').casefold() == species.casefold()), None)
            if organism is None or ('typing' in organism and 'MLST' not in (organism.get('typing') or [])):
                pool['reason'] = 'unsupported_provider_species_or_mlst'
                continue
            preview = client.request_json(
                'POST', '/api/search/genomes',
                body={'organismId': str(organism['organismId']), 'mlst': [st], 'qc': [True, False]},
                params={'limit': 1, 'sort': 'id'},
            )
            count = preview.get('meta', {}).get('count') if isinstance(preview, dict) else None
            if type(count) is not int or count < 0:
                raise DatasetError('Public context candidate count is unavailable; cannot estimate retrieval scope')
            pool.update(count_source='live_public_count_preview',
                        captured_at=datetime.now(timezone.utc).isoformat())
        pool.update(organism_id=str(organism['organismId']), candidate_count=count,
                    status='over_limit' if max_candidates is not None and count > max_candidates else 'ready',
                    cgmlst_export_batches_per_analysis_upper_bound=(count+99)//100)
        species_candidates[species] += count
        progress('context', f'{species} ST{st}: {count} public candidates before profile compatibility; '
                 'the context-size selection budget does not limit retrieval')
    upper = {name: count+species_inputs[name] for name, count in species_candidates.items()}
    plan = {
        'schema': 'chronoclade.context-preflight', 'schema_version': 1,
        'source_dataset_sha256': file_sha256(source),
        'provider_base_url': client.base_url.rstrip('/'),
        'created_at': datetime.now(timezone.utc).isoformat(),
        'primary_search_field': 'mlst', 'max_candidates_per_pool': max_candidates,
        'download_credentials_available': bool(getattr(client, 'api_key', None)),
        'candidate_limit_policy': 'refuse whole oversized pools; never truncate',
        'pools': pools,
        'status': 'ready' if all(p['status'] == 'ready' for p in pools) and pools else 'blocked',
        'cost': {
            'public_candidate_records': sum(species_candidates.values()),
            'metadata_detail_requests_upper_bound': sum(species_candidates.values()),
            'cgmlst_export_batches_per_analysis_upper_bound': sum(p.get('cgmlst_export_batches_per_analysis_upper_bound', 0) for p in pools),
            'guide_records_upper_bound_by_species': upper,
            'pairwise_comparisons_upper_bound': sum(n*(n-1)//2 for n in upper.values()),
            'one_float64_distance_matrix_bytes_upper_bound': sum(n*n*8 for n in upper.values()),
            'estimate_scope': 'Before deduplication/profile compatibility; one matrix only. Actual peak memory includes additional matrices, profiles and working arrays. Export estimates apply per advertised cgMLST analysis; they exclude optional lineage exports and polling. Frozen counts are historical, not a fresh snapshot.',
        },
        'selection_policy': 'Retrieve complete compatible pools before nearest-neighbour selection; preserve all inputs and mandatory ties. Context-size only limits downstream optional context selection.',
        'temporal_policy': 'Sampling dates must independently support dating; counts and context selection do not override the temporal gate.',
    }
    plan['plan_sha256'] = content_hash(plan)
    return plan


def require_allowed_pools(plan):
    if not plan['pools']:
        raise DatasetError('Public context preflight requires at least one input genome')
    oversized = [p for p in plan['pools'] if p['status'] == 'over_limit']
    if oversized:
        descriptions = '; '.join(f"{p['species']} ST{p['st']} has {p['candidate_count']} public candidates" for p in oversized)
        raise DatasetError(f'{descriptions}, exceeding --context-max-candidates {plan["max_candidates_per_pool"]}. '
                           'No pool was truncated; increase the limit deliberately or run local prepared-profile analysis.')
    unavailable = [p for p in plan['pools'] if p['status'] == 'unavailable']
    if unavailable:
        raise DatasetError('Public context preflight is incomplete: '+ '; '.join(
            f"{p['species']} ST{p['st'] or 'unknown'}: {p['reason']}" for p in unavailable))


def write_context_plan(plan, output, *, allow_existing=False):
    directory = Path(output).expanduser().resolve()
    if directory.exists() and any(directory.iterdir()) and not allow_existing:
        raise DatasetError(f'Context preflight output must be empty: {directory}')
    directory.mkdir(parents=True, exist_ok=True)
    write_json(directory/'context_preflight.json', plan)
    reasons = {
        'missing_primary_mlst_st': 'Primary seven-locus ST is missing. Run preparation with assembly typing before public discovery.',
        'alternative_mlst_namespace_not_primary_achtman': 'This E. coli ST uses an alternative MLST scheme; public discovery requires the primary Achtman scheme.',
        'unsupported_provider_species_or_mlst': 'The provider does not advertise primary MLST discovery for this species.',
    }
    limit = plan.get('max_candidates_per_pool')

    def pool_labels(pool):
        if pool['status'] == 'over_limit':
            return ('Exceeds candidate limit', f"Limit: {limit:,} records per species/ST pool. No records will be truncated.")
        if pool['status'] == 'unavailable':
            return ('Discovery unavailable', reasons.get(pool.get('reason'), 'Discovery scope could not be established.'))
        return ('Count checked', f'Limit: {limit:,} records per species/ST pool.' if limit is not None else 'No candidate limit configured.')

    headers = '<tr><th>Species</th><th>ST</th><th>Public candidates</th><th>Retrieval status</th><th>Details</th></tr>'
    body = ''.join('<tr>'+''.join('<td>'+escape(str(value))+'</td>' for value in
                   [p['species'], p['st'] or 'Unknown', f"{p['candidate_count']:,}" if p['candidate_count'] is not None else 'Unknown', *pool_labels(p)])+'</tr>'
                   for p in plan['pools'])
    cost = plan['cost']
    matrix_bytes = cost['one_float64_distance_matrix_bytes_upper_bound']
    matrix_size = f'{matrix_bytes:,} bytes'
    for unit, divisor in [('TiB', 1024**4), ('GiB', 1024**3), ('MiB', 1024**2), ('KiB', 1024)]:
        if matrix_bytes >= divisor:
            matrix_size = f'{matrix_bytes/divisor:,.1f} {unit} ({matrix_bytes:,} bytes)'
            break
    status = 'Counts checked' if plan['status'] == 'ready' else 'Needs attention'
    credentials = 'Configured' if plan['download_credentials_available'] else 'Missing; protected credentials are required to download profiles'
    document = ('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
                '<title>Public context preflight</title><style>body{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;line-height:1.5;color:#243445}'
                'table{border-collapse:collapse;width:100%;font-size:.95rem}th,td{padding:.7rem;text-align:left;border-bottom:1px solid #d5dde5;vertical-align:top}th{background:#edf3f7}'
                'a{color:#155e85}p{max-width:95ch}</style></head>'
                '<body><h1>Public context preflight</h1><p>Status: '+status+'</p><table>'+headers+body+'</table>'
                '<p>Total advertised source records: '+f"{cost['public_candidate_records']:,}"+
                '; profile export batches per analysis (upper bound): '+f"{cost['cgmlst_export_batches_per_analysis_upper_bound']:,}"+
                '; one distance matrix (upper bound): '+matrix_size+'.</p>'
                '<p>'+escape(cost['estimate_scope'])+'</p><p>'+escape(plan['selection_policy'])+'</p>'
                '<p>Protected profile-download credentials: '+credentials+'. '
                'This preview performs no genome-detail, profile or assembly downloads.</p>'
                '<a href="context_preflight.json" download>Download preflight JSON</a></body></html>')
    from chronoclade.report_components.downloads import enhance_report_downloads
    (directory/'index.html').write_text(enhance_report_downloads(document, directory))
    return directory/'context_preflight.json'
