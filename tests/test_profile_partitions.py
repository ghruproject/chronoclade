import json

import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.context_partitions import profile_records
from chronoclade.datasets import LocusCatalogue, from_profile_records, write_dataset
from chronoclade.profile_partitions import partition_input_profiles
from chronoclade.cgmlst import workflow
from chronoclade.selection_manifest import load_selection_manifest


def source(tmp_path):
    records = [dict(sample_id=f'e{i}', origin='local', species='Escherichia coli',
                    mlst_st='95', mlst_scheme='Achtman', country='Australia',
                    collection_date=str(2010+i), cgmlst_scheme='fixture',
                    cgmlst_scheme_version='1', cgmlst_status='resolved',
                    cgmlst_profile={'a':str(i+1), 'b':'1', 'c':'1'}) for i in range(3)]
    dataset = from_profile_records(records, [LocusCatalogue('fixture', '1', ('a','b','c'), 'synthetic')])
    return dataset, write_dataset(dataset, tmp_path/'prepared')


def test_st_partition_scope_separation_and_context_rejection(tmp_path):
    dataset, _ = source(tmp_path)
    records = profile_records(dataset)
    records[1]['mlst_scheme'] = 'Pasteur'
    records[2]['mlst_scheme'] = None
    blocks, audit = partition_input_profiles(records)
    assert len(blocks) == 2
    assert audit['unresolved_inputs'] == [{'sample_id':'e2', 'reason':'ambiguous_mlst_namespace'}]
    assert all(block['kind'] == 'st' and block['prefix'] == [] for block in blocks)
    records[0]['role'] = 'context'
    with pytest.raises(ValueError, match='inputs only'):
        partition_input_profiles(records)


def test_opt_in_normal_workflow_preserves_exact_selection_and_unknown_lineage(tmp_path):
    _, manifest = source(tmp_path)
    default = workflow.run_cgmlst(manifest, tmp_path/'default', bootstrap_replicates=0)
    assert json.loads(default.manifest_path.read_text())['status'] == 'unavailable'
    result = workflow.run_cgmlst(manifest, tmp_path/'st', partition_by='st', bootstrap_replicates=0)
    saved = json.loads(result.manifest_path.read_text())
    assert saved['status'] == 'complete'
    assert saved['parameters']['partition_by'] == 'st'
    assert len(saved['blocks']) == 1
    block = saved['blocks'][0]
    assert 'no lineage partition' in block['label']
    ensemble = json.loads((tmp_path/'st'/block['ensemble']).read_text())
    assert ensemble
    selection_file = next((tmp_path/'st').rglob('selection-*.json'))
    selection = load_selection_manifest(selection_file)
    assert set(selection["selected_sample_ids"]) == {'e0','e1','e2'}
    partitions = json.loads((tmp_path/'st/partitions.json').read_text())
    assert partitions['blocks'][0]['sample_ids'] == ['e0','e1','e2']
    assert 'no public context' in partitions['context_scope']


def test_partition_mode_is_bound_to_interrupted_resume(monkeypatch, tmp_path):
    _, manifest = source(tmp_path)
    def fail(*args, **kwargs):
        raise RuntimeError('packaging interrupted')
    monkeypatch.setattr(workflow, 'write_profile_report', fail)
    with pytest.raises(RuntimeError, match='packaging interrupted'):
        workflow.run_cgmlst(manifest, tmp_path/'run', partition_by='st', bootstrap_replicates=0)
    with pytest.raises(ValueError, match='parameters differ'):
        workflow.run_cgmlst(manifest, tmp_path/'run', partition_by='lineage', bootstrap_replicates=0, resume=True)


def test_cli_blocks_context_before_any_retrieval(tmp_path):
    _, manifest = source(tmp_path)
    result = CliRunner().invoke(app, ['cgmlst', str(manifest), '--partition-by', 'st', '--fetch-context', '--out', str(tmp_path/'run')])
    assert result.exit_code == 2
    assert 'cannot fetch context' in result.stdout
    assert not (tmp_path/'run').exists()
