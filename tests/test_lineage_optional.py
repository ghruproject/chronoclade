"""Profiles drive comparison and selection; nomenclature is optional evidence."""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.cgmlst.workflow import run_cgmlst
from chronoclade.context_partitions import profile_records
from chronoclade.datasets import LocusCatalogue, from_profile_records, load_dataset, write_dataset
from chronoclade.profile_context_provider import discover_profile_context
from chronoclade.profile_partitions import partition_profiles


def record(ident, *, st=None, role='input', **kwargs):
    return dict(sample_id=ident, species='Staphylococcus aureus', role=role,
                mlst_st=st, cgmlst_scheme='fixture', cgmlst_scheme_version='1',
                cgmlst_profile={'a': '1', 'b': '1'}, cgmlst_status='resolved', **kwargs)


def dataset(tmp_path, rows):
    ds = from_profile_records(rows, [LocusCatalogue('fixture', '1', ('a', 'b'), 'fixture')])
    return write_dataset(ds, tmp_path / 'prepared')


def test_default_profiles_do_not_require_st_or_lineages_and_keep_all_inputs(tmp_path):
    source = dataset(tmp_path, [record('q1'), record('q2', st='1'), record('q3', st='45')])
    result = run_cgmlst(source, tmp_path / 'analysis', bootstrap_replicates=0)
    manifest = json.loads(result.manifest_path.read_text())
    assert manifest['status'] == 'complete'
    assert manifest['parameters']['partition_by'] == 'profiles'
    assert len(manifest['blocks']) == 1
    block = manifest['blocks'][0]
    ensemble_path = tmp_path / 'analysis' / block['ensemble']
    ensemble = json.loads(ensemble_path.read_text())
    selected = json.loads((ensemble_path.parent / ensemble['selections'][0]['path']).read_text())
    assert set(selected['selected_sample_ids']) == {'q1', 'q2', 'q3'}
    analysis = json.loads((tmp_path / 'analysis' / block['block_id'] / 'profile_analysis.json').read_text())
    assert analysis['display_group_basis'] == 'exploratory_profile_distance'


def test_profile_partition_ignores_label_and_st_names_but_separates_scheme_identity(tmp_path):
    source = dataset(tmp_path, [record('a'), record('b'), record('c')])
    rows = profile_records(load_dataset(source))
    rows[1].update(mlst_st='99', cglin_raw='1.2.3.4.5', cglin_scheme='unrelated',
                   cglin_scheme_version='other', cglin_status='complete')
    blocks, _ = partition_profiles(rows)
    assert len(blocks) == 1
    rows[2]['cgmlst_database_sha256'] = 'a' * 64
    blocks, audit = partition_profiles(rows)
    assert len(blocks) == 2
    assert not audit['unresolved_inputs']


def test_generic_public_context_does_not_gate_on_lineage_and_preserves_nearest_ties(tmp_path, monkeypatch):
    from chronoclade import profile_context_provider as provider
    source = dataset(tmp_path, [record('q', st='1')])
    public = [dict(record('a', st='1'), source_genome_id='source-a', numeric_source_id=1),
              dict(record('b', st='1'), source_genome_id='source-b', numeric_source_id=2)]

    class Client:
        api_key = 'fixture-key'
        base_url = 'https://pathogen.watch'

        def supported_organisms(self):
            return [{'fullName': 'Staphylococcus aureus', 'organismId': '1280'}]

        def freeze_catalogue(self, path, **kwargs):
            value = {'rows': public, 'provenance': {'complete': True, 'expected_count': 2}}
            path.write_text(json.dumps(value))
            return value

    requests = []
    def exports(rows, client, output, provenance, download_names=None, resume=False):
        requests.append(download_names)
        assert download_names == {'cgmlst'}
        return rows
    monkeypatch.setattr(provider, '_grouped_analysis_exports', exports)
    prepared = discover_profile_context(source, tmp_path / 'context', client=Client(), lineage_filter=False)
    combined = load_dataset(prepared.dataset_manifest)
    assert len(combined.samples) == 3
    assert requests == [{'cgmlst'}]
    result = run_cgmlst(prepared.dataset_manifest, tmp_path / 'analysis', context_size=0,
                        nearest_per_query=1, bootstrap_replicates=0)
    manifest = json.loads(result.manifest_path.read_text())
    ensemble_path = result.manifest_path.parent / manifest['blocks'][0]['ensemble']
    ensemble = json.loads(ensemble_path.read_text())
    selected = json.loads((ensemble_path.parent / ensemble['selections'][0]['path']).read_text())
    assert len(selected['selected_sample_ids']) == 3  # Both tied neighbours are mandatory.


def test_optional_lineage_export_failure_does_not_block_profiles(tmp_path, monkeypatch):
    from chronoclade import profile_context_provider as provider
    rows = [dict(record('q', st='1'), species='Klebsiella pneumoniae')]
    source = dataset(tmp_path, rows)
    public = [dict(rows[0], sample_id='a', source_genome_id='source-a', numeric_source_id=1)]

    class Client:
        api_key = 'fixture-key'
        base_url = 'https://pathogen.watch'
        def supported_organisms(self):
            return [{'fullName': 'Klebsiella pneumoniae', 'organismId': '573'}]
        def freeze_catalogue(self, path, **kwargs):
            value = {'rows': public, 'provenance': {'complete': True, 'expected_count': 1}}
            path.write_text(json.dumps(value))
            return value
    def exports(rows, client, output, provenance, download_names=None, resume=False):
        if download_names == {'klebsiella-lincodes'}:
            provenance['analysis_export_status'] = 'unavailable'
        return rows
    monkeypatch.setattr(provider, '_grouped_analysis_exports', exports)
    prepared = discover_profile_context(source, tmp_path / 'context', client=Client(), lineage_filter=False)
    assert len(load_dataset(prepared.dataset_manifest).samples) == 2


def test_missing_st_blocks_only_public_discovery_and_incomplete_cli_is_nonzero(tmp_path):
    source = dataset(tmp_path, [record('q')])
    with pytest.raises(ValueError, match='local cgMLST analysis does not'):
        discover_profile_context(source, tmp_path / 'context', lineage_filter=False)
    result = CliRunner().invoke(app, ['cgmlst', str(source), '--partition-by', 'lineage',
                                     '--out', str(tmp_path / 'unavailable')])
    assert result.exit_code == 2
    assert 'Analysis unavailable' in result.stdout
    assert Path(tmp_path / 'unavailable/cgmlst.json').is_file()


@pytest.mark.parametrize('scheme', ['ecoli', 'ecoli_pasteur', 'pasteur', 'ECOLI'])
def test_alternative_ecoli_mlst_namespace_is_rejected_before_public_retrieval(tmp_path, scheme):
    source = dataset(tmp_path, [dict(record('q', st='22'), species='Escherichia coli',
                                     mlst_scheme=scheme)])

    class Client:
        api_key = 'fixture-key'
        base_url = 'https://pathogen.watch'

        def supported_organisms(self):
            pytest.fail('Alternative MLST namespace must be rejected before public retrieval')

        def freeze_catalogue(self, *args, **kwargs):
            pytest.fail('Alternative ST must not be used to retrieve a primary-namespace pool')

    with pytest.raises(ValueError, match='Achtman MLST.*alternative Pasteur namespace'):
        discover_profile_context(source, tmp_path / 'context', client=Client(), lineage_filter=False)
    assert not (tmp_path / 'context/dataset.json').exists()


def test_candidate_guard_refuses_without_truncation_or_detail_fetch(tmp_path):
    source = dataset(tmp_path, [record('q', st='1')])
    class Client:
        api_key = 'fixture-key'
        base_url = 'https://pathogen.watch'
        def supported_organisms(self):
            return [{'fullName': 'Staphylococcus aureus', 'organismId': '1280'}]
        def request_json(self, method, path, *, body=None, params=None):
            assert path == '/api/search/genomes' and params['limit'] == 1
            return {'meta': {'count': 11}, 'genomes': []}
        def freeze_catalogue(self, *args, **kwargs):
            pytest.fail('Candidate guard must run before fetching details')
    with pytest.raises(ValueError, match='No pool was truncated'):
        discover_profile_context(source, tmp_path / 'context', client=Client(),
                                 lineage_filter=False, max_candidates=10)
    assert not (tmp_path / 'context/dataset.json').exists()


def test_generic_missing_credentials_fail_before_network(tmp_path):
    source = dataset(tmp_path, [record('q', st='1')])
    class Client:
        api_key = ''
        base_url = 'https://pathogen.watch'
        def supported_organisms(self):
            pytest.fail('Missing download credentials must be reported before retrieval')
    with pytest.raises(ValueError, match='protected Pathogenwatch credentials'):
        discover_profile_context(source, tmp_path / 'context', client=Client(), lineage_filter=False)
