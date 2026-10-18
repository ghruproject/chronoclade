"""Count all discovery pools before expensive retrieval, without exposing credentials."""
import json

import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.context_preflight import plan_context, require_allowed_pools
from chronoclade.datasets import LocusCatalogue, from_profile_records, write_dataset
from chronoclade.profile_context_provider import discover_profile_context


def test_preflight_report_explains_cost_and_limits_without_changing_plan(tmp_path):
    from chronoclade.context_preflight import write_context_plan

    plan = plan_context(source(tmp_path), client=Client({'1': 6890, '149': 23}), max_candidates=300)
    plan['cost']['one_float64_distance_matrix_bytes_upper_bound'] = 30_000_000_000
    output = tmp_path / 'cost-report'
    write_context_plan(plan, output)
    document = (output / 'index.html').read_text()
    assert 'Needs attention' in document
    assert 'Exceeds candidate limit' in document
    assert '6,890' in document and 'Limit: 300 records per species/ST pool' in document
    assert '27.9 GiB (30,000,000,000 bytes)' in document
    assert 'No records will be truncated' in document
    assert 'Missing; protected credentials are required to download profiles' in document
    assert 'over_limit' not in document
    assert json.loads((output / 'context_preflight.json').read_text()) == plan


def test_preflight_report_explains_missing_typing_and_escapes_provider_values(tmp_path):
    from chronoclade.context_preflight import write_context_plan

    plan = plan_context(source(tmp_path), client=Client(), max_candidates=30)
    plan['pools'][0].update(species='<script>provider</script>', st='', candidate_count=None,
                            status='unavailable', reason='missing_primary_mlst_st')
    plan.update(status='blocked', download_credentials_available=True)
    output = tmp_path / 'missing-typing-report'
    write_context_plan(plan, output)
    document = (output / 'index.html').read_text()
    assert 'Discovery unavailable' in document
    assert 'Primary seven-locus ST is missing' in document
    assert 'Run preparation with assembly typing' in document
    assert '&lt;script&gt;provider&lt;/script&gt;' in document
    assert '<script>provider</script>' not in document
    assert 'missing_primary_mlst_st' not in document
    assert 'Protected profile-download credentials: Configured' in document
    assert json.loads((output / 'context_preflight.json').read_text()) == plan


def source(tmp_path):
    rows = [dict(sample_id=f'q{st}', role='input', species='Staphylococcus aureus',
                 mlst_st=st, cgmlst_scheme='fixture', cgmlst_scheme_version='1',
                 cgmlst_profile={'a': '1', 'b': '2'}, cgmlst_status='resolved')
            for st in ['1', '149']]
    return write_dataset(from_profile_records(rows, [LocusCatalogue('fixture', '1', ('a', 'b'), 'fixture')]),
                         tmp_path/'prepared')


class Client:
    api_key = ''
    base_url = 'https://pathogen.watch'

    def __init__(self, counts=None):
        self.counts = counts or {'1': 12, '149': 23}
        self.calls = []

    def supported_organisms(self):
        return [{'fullName': 'Staphylococcus aureus', 'organismId': '1280', 'typing': ['MLST']}]

    def request_json(self, method, path, *, body=None, params=None):
        assert method == 'POST' and path == '/api/search/genomes'
        assert params == {'limit': 1, 'sort': 'id'}
        assert body['qc'] == [True, False]
        self.calls.append(body['mlst'][0])
        return {'meta': {'count': self.counts[body['mlst'][0]]}, 'genomes': []}

    def freeze_catalogue(self, *args, **kwargs):
        pytest.fail('Count planning and guards must precede every detailed retrieval')


def test_all_pool_counts_and_cost_are_available_without_credentials(tmp_path):
    client = Client()
    plan = plan_context(source(tmp_path), client=client, max_candidates=20)
    assert client.calls == ['1', '149']
    assert plan['status'] == 'blocked'
    assert plan['download_credentials_available'] is False
    assert [p['status'] for p in plan['pools']] == ['ready', 'over_limit']
    assert plan['cost']['public_candidate_records'] == 35
    assert plan['cost']['cgmlst_export_batches_per_analysis_upper_bound'] == 2
    assert plan['cost']['one_float64_distance_matrix_bytes_upper_bound'] == 37*37*8
    with pytest.raises(ValueError, match='No pool was truncated'):
        require_allowed_pools(plan)


def test_normal_fetch_preflights_all_pools_before_any_details(tmp_path):
    client = Client()
    client.api_key = 'fixture-secret-never-publish'
    with pytest.raises(ValueError, match='ST149 has 23 public candidates'):
        discover_profile_context(source(tmp_path), tmp_path/'context', client=client,
                                 lineage_filter=False, max_candidates=20)
    assert client.calls == ['1', '149']
    assert not (tmp_path/'context/dataset.json').exists()
    plan = tmp_path/'context_retrieval/preflight/context_preflight.json'
    assert json.loads(plan.read_text())['cost']['public_candidate_records'] == 35
    assert client.api_key not in ''.join(p.read_text() for p in plan.parent.iterdir())


def test_count_only_cli_saves_report_without_analysis_or_downloads(tmp_path, monkeypatch):
    import chronoclade.context_preflight as preflight
    client = Client()
    monkeypatch.setattr(preflight, 'PathogenwatchClient', lambda: client)
    output = tmp_path/'plan'
    result = CliRunner().invoke(app, ['cgmlst', str(source(tmp_path)), '--context-preflight',
                                     '--context-max-candidates', '30', '--out', str(output)])
    assert result.exit_code == 0, result.output
    assert '35 public candidate records' in result.output
    assert 'credentials are unavailable' in result.output
    assert (output/'context_preflight.json').is_file() and (output/'index.html').is_file()
    assert not (output/'cgmlst.json').exists()
    assert 'no genome-detail, profile or assembly downloads' in (output/'index.html').read_text()
    import base64
    import re
    document = (output/'index.html').read_text()
    payload = re.search(r'data-cc-download="([^"]+)"', document)[1]
    assert base64.b64decode(payload) == (output/'context_preflight.json').read_bytes()
    assert 'id="cc-offline-downloads"' in document


def test_count_only_cli_refusal_still_preserves_complete_plan(tmp_path, monkeypatch):
    import chronoclade.context_preflight as preflight
    client = Client()
    monkeypatch.setattr(preflight, 'PathogenwatchClient', lambda: client)
    output = tmp_path/'plan'
    result = CliRunner().invoke(app, ['cgmlst', str(source(tmp_path)), '--context-preflight',
                                     '--context-max-candidates', '20', '--out', str(output)])
    assert result.exit_code == 2, result.output
    plan = json.loads((output/'context_preflight.json').read_text())
    assert len(plan['pools']) == 2 and plan['status'] == 'blocked'
    assert not (output/'cgmlst.json').exists()


def test_every_oversized_pool_is_reported_after_all_counts_are_collected(tmp_path):
    client = Client()
    plan = plan_context(source(tmp_path), client=client, max_candidates=10)
    assert client.calls == ['1', '149']
    with pytest.raises(ValueError) as error:
        require_allowed_pools(plan)
    assert 'ST1 has 12 public candidates' in str(error.value)
    assert 'ST149 has 23 public candidates' in str(error.value)


def test_increasing_refusal_guard_on_resume_keeps_the_same_full_scopes(tmp_path):
    class RetryingClient(Client):
        api_key = 'fixture-key'

        def freeze_catalogue(self, path, **kwargs):
            assert kwargs['resume'] is True
            raise RuntimeError('Reached complete-pool retrieval after raised guard')
    client = RetryingClient()
    manifest = source(tmp_path)
    target = tmp_path/'context'
    with pytest.raises(ValueError, match='No pool was truncated'):
        discover_profile_context(manifest, target, client=client, lineage_filter=False,
                                 max_candidates=20)
    with pytest.raises(RuntimeError, match='complete-pool retrieval'):
        discover_profile_context(manifest, target, client=client, lineage_filter=False,
                                 max_candidates=30, resume=True)
    assert client.calls == ['1', '149', '1', '149']
    plan = json.loads((tmp_path/'context_retrieval/preflight/context_preflight.json').read_text())
    assert plan['cost']['public_candidate_records'] == 35
    assert [p['candidate_count'] for p in plan['pools']] == [12, 23]
    assert not (target/'dataset.json').exists()


def test_context_preparation_rejects_another_writer_before_retrieval(tmp_path):
    from chronoclade.prepare_recovery import DirectoryLock
    manifest = source(tmp_path)
    with DirectoryLock(tmp_path/'.context.context-retrieval.lock'):
        with pytest.raises(ValueError, match='Another process owns recovery lock'):
            discover_profile_context(manifest, tmp_path/'context', client=Client(), lineage_filter=False)
    assert not (tmp_path/'context_retrieval').exists()


def test_preflight_provider_error_is_a_concise_cli_failure(tmp_path, monkeypatch):
    import chronoclade.context_preflight as preflight
    from chronoclade.pathogenwatch import PathogenwatchError
    def interrupted(*args, **kwargs):
        raise PathogenwatchError('Pathogenwatch /api/search/genomes transport failed after retries')
    monkeypatch.setattr(preflight, 'plan_context', interrupted)
    output = tmp_path/'plan'
    result = CliRunner().invoke(app, ['cgmlst', str(source(tmp_path)), '--context-preflight',
                                     '--out', str(output)])
    assert result.exit_code == 2, result.output + repr(result.exception)
    assert 'Error: Pathogenwatch /api/search/genomes transport failed after retries' in ' '.join(result.output.split())
    assert 'Traceback' not in result.output
    assert not (output/'cgmlst.json').exists()
    assert not (output/'context_preflight.json').exists()


@pytest.mark.parametrize('reason', [
    'Pathogenwatch /api/genomes/details transport failed after retries',
    'Search count changed during pagination',
    'Catalogue capture receipt content hash mismatch',
])
def test_fetch_provider_error_retains_cache_and_reports_cli_recovery(tmp_path, monkeypatch, reason):
    import chronoclade.profile_context_provider as provider
    from chronoclade.pathogenwatch import PathogenwatchError
    output = tmp_path/'analysis'
    cache = tmp_path/'analysis_context_retrieval'
    cache.mkdir()
    accepted = cache/'accepted-request.raw'
    accepted.write_bytes(b'previous accepted retrieval evidence')
    def interrupted(*args, **kwargs):
        raise PathogenwatchError(reason)
    monkeypatch.setattr(provider, 'discover_profile_context', interrupted)
    result = CliRunner().invoke(app, ['cgmlst', str(source(tmp_path)), '--fetch-context',
                                     '--out', str(output)])
    assert result.exit_code == 2, result.output + repr(result.exception)
    assert reason in ' '.join(result.output.split()) and 'Error:' in result.output
    assert 'Traceback' not in result.output
    assert '--resume' in result.output and 'cache' in result.output.lower()
    assert accepted.read_bytes() == b'previous accepted retrieval evidence'
    assert not (output/'cgmlst.json').exists()
