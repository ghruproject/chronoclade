"""Verticall reference coordinates, assembly repair and cohort-safe masking."""

import csv
import json
from pathlib import Path
import random
import shutil

from Bio import SeqIO
import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.errors import WorkflowError
from chronoclade.metadata import Sample
from chronoclade.native_pipeline import LineageFiles
from chronoclade.verticall_filter import read_regions, repair_assembly, run_verticall_mask


FIELDS = ['assembly_a', 'assembly_b', 'assembly_a_vertical_regions',
          'assembly_a_horizontal_regions', 'assembly_a_unaligned_regions']


def write_regions(path, values):
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(FIELDS)
        writer.writerows(values)


def test_repair_splits_ambiguities_and_preserves_called_bases(tmp_path):
    source = tmp_path / 'original.fasta'
    source.write_text('>duplicate\nAAANNccRtt\n>duplicate\nGG\n')
    target = tmp_path / 'repaired.fasta'
    receipt = repair_assembly(source, target)
    seqs = list(SeqIO.parse(target, 'fasta'))
    assert [str(r.seq) for r in seqs] == ['AAA', 'CC', 'TT', 'GG']
    assert len({r.id for r in seqs}) == 4
    assert receipt['ambiguous_bases_removed'] == 3 and receipt['called_bases'] == 9
    assert source.read_text() == '>duplicate\nAAANNccRtt\n>duplicate\nGG\n'


def test_regions_use_first_result_and_require_complete_classification(tmp_path):
    path = tmp_path / 'pairwise.tsv'
    write_regions(path, [['reference', 'sample000000', 'reference:0-20,reference:30-100', 'reference:20-30', ''],
                         ['reference', 'sample000000', 'reference:0-100', '', '']])
    result, extra = read_regions(path, {'sample000000': 'q'}, 100)
    assert result['sample000000']['horizontal'] == [(20, 30)] and extra == 1
    with pytest.raises(WorkflowError, match='omitted'):
        read_regions(path, {'sample000000': 'q', 'sample000001': 'c'}, 100)


@pytest.mark.parametrize('region', ['reference:0-101', 'reference:1-100', 'reference:0-50,reference:49-100', 'other:0-100'])
def test_invalid_reference_regions_are_refused(tmp_path, region):
    path = tmp_path / 'pairwise.tsv'
    write_regions(path, [['reference', 'sample000000', region, '', '']])
    with pytest.raises(WorkflowError):
        read_regions(path, {'sample000000': 'q'}, 100)


def test_two_reference_records_and_native_mask_are_reconciled(tmp_path, monkeypatch):
    files = LineageFiles.in_directory(tmp_path)
    members = []
    raw = 'A' * 500 + 'N' + 'A' * 1099
    for ident in ['q', 'c', 'd']:
        assembly = tmp_path / f'{ident}.fasta'
        assembly.write_text('>contig\n' + raw + '\n')
        members.append(Sample(ident, assembly, '', '', 'sp', 'group', 'local'))
    files.alignment.write_text(''.join(f'>{m.sample_id}\n{raw}\n' for m in members))
    reference = tmp_path / 'reference.fasta'
    reference.write_text('>chromosome\n' + 'A' * 1000 + '\n>plasmid\n' + 'A' * 600 + '\n')

    def native(command, **kwargs):
        if command[1] == 'pairwise':
            ref = Path(command[command.index('-r') + 1])
            size = len(next(SeqIO.parse(ref, 'fasta')))
            path = Path(command[command.index('-o') + 1])
            write_regions(path, [['reference', 'sample000000', f'reference:0-10,reference:20-{size}', 'reference:10-20', ''],
                                 ['reference', 'sample000001', f'reference:0-30,reference:40-{size}', '', 'reference:30-40'],
                                 ['reference', 'sample000002', f'reference:0-{size}', '', '']])
        elif command[1] == 'mask':
            rows = list(SeqIO.parse(command[command.index('-a') + 1], 'fasta'))
            assert rows[0].id == 'reference' and set(rows[0].seq) == {'A'}
            output = Path(command[command.index('-o') + 1])
            with output.open('w') as stream:
                for row in rows[1:]:
                    seq = str(row.seq)
                    if row.id == 'sample000000':
                        seq = seq[:10] + 'N' * 10 + seq[20:]
                    if row.id == 'sample000001':
                        seq = seq[:30] + '-' * 10 + seq[40:]
                    stream.write(f'>{row.id}\n{seq}\n')
        else:
            pytest.fail('Unexpected native command')

    monkeypatch.setattr('chronoclade.native_pipeline._run_command', native)
    filtered, summary = run_verticall_mask(files, members, reference, threads=2, force=False)
    assert summary['recombination_removed_sites'] == 20
    assert summary['unaligned_only_removed_sites'] == 20
    assert summary['incomplete_only_removed_sites'] == 1
    assert summary['retained_clonal_sites'] == 1559
    assert [r.id for r in SeqIO.parse(filtered, 'fasta')] == ['q', 'c', 'd']
    assert all(str(r.seq) == 'A' * 1559 for r in SeqIO.parse(filtered, 'fasta'))
    assert [r['record'] for r in summary['reference_records']] == ['chromosome', 'plasmid']
    assert summary['repaired_assemblies'] == 3
    repairs = json.loads((tmp_path / 'verticall/assembly_repair.json').read_text())
    assert all(r['ambiguous_bases_removed'] == 1 for r in repairs)


def test_cli_and_slurm_accept_verticall(tmp_path, monkeypatch):
    from test_tree_stage import selection_fixture
    from chronoclade.execution.jobs import SlurmResources
    from chronoclade.execution.slurm import write_slurm_ensemble

    selection, ensemble = selection_fixture(tmp_path)
    captured = []

    def run(*args, **kwargs):
        captured.append(kwargs)
        return {'temporal_assessment': {'code': 'not_assessed'}}

    monkeypatch.setattr('chronoclade.tree_stage.workflow.run_tree', run)
    result = CliRunner().invoke(app, ['tree', str(selection), '--recombination-method', 'verticall'])
    assert result.exit_code == 0, result.output
    assert captured[0]['recombination_method'] == 'verticall'
    path = write_slurm_ensemble('tree', ensemble, output=tmp_path / 'jobs', results=tmp_path / 'results',
                               resources=SlurmResources(2, 4, '00:10:00'), options={'recombination_method': 'verticall'})
    assert json.loads(path.read_text())['options']['recombination_method'] == 'verticall'


def test_workflow_records_verticall_and_retains_report_evidence(tmp_path, monkeypatch):
    from test_tree_stage import selection_fixture, native_stubs
    from chronoclade.tree_stage.workflow import run_tree
    from chronoclade.stage_artifacts import load_tree_result
    import chronoclade.native_pipeline as native

    native_stubs.__wrapped__(monkeypatch)
    original = native._run_core_phylogeny

    def core(*args, **kwargs):
        assert kwargs['recombination_method'] == 'verticall'
        tree, alignment, _ = original(*args, **kwargs)
        return tree, alignment, dict(method='verticall', samples=3,
                                     reference_records=[{'record': 'chr'}],
                                     recombination_removed_sites=0, retained_clonal_sites=2000,
                                     alignment_policy='Synthetic complete-site policy',
                                     limitation='Synthetic native-tool fixture')

    monkeypatch.setattr(native, '_run_core_phylogeny', core)
    selection, _ = selection_fixture(tmp_path, dates=('', '', ''))
    result = run_tree(selection, output=tmp_path / 'tree', recombination_method='verticall')
    assert result['parameters']['recombination_method'] == 'verticall'
    assert result['recombination']['method'] == 'verticall'
    assert 'verticall' in result['tools'] and 'minimap2' in result['tools']
    assert 'ClonalFrameML' not in result['tools'] and 'Profile' not in result['tools']
    assert load_tree_result(tmp_path / 'tree/tree.json') == result
    assert 'Recombination filtering: Verticall' in (tmp_path / 'tree' / result['report']['path']).read_text()


def test_preflight_only_requires_selected_verticall_tools(monkeypatch):
    monkeypatch.setattr('chronoclade.native_tools.native_platform_supported', lambda: True)
    available = {name: '/bin/' + name for name in ['ska', 'iqtree', 'verticall', 'minimap2', 'treetime']}
    monkeypatch.setattr('chronoclade.native_tools.tool_status', lambda: available)
    result = CliRunner().invoke(app, ['preflight', '--stage', 'tree', '--recombination-method', 'verticall'])
    assert result.exit_code == 0
    assert 'ClonalFrameML' not in result.output and 'Profile' not in result.output
    available.pop('minimap2')
    missing = CliRunner().invoke(app, ['preflight', '--stage', 'tree', '--recombination-method', 'verticall'])
    assert missing.exit_code == 2 and 'minimap2' in missing.output


@pytest.mark.integration
@pytest.mark.skipif(not shutil.which('verticall') or not shutil.which('minimap2'), reason='Verticall/minimap2 not installed')
def test_installed_verticall_retains_exact_genomes_and_native_mask(tmp_path):
    rng = random.Random(42)
    original = ''.join(rng.choices('ACGT', k=50000))
    altered = list(original)
    for pos in range(100, 50000, 1000):
        altered[pos] = next(c for c in 'ACGT' if c != altered[pos])
    members = []
    for ident, seq in [('q', original), ('c', ''.join(altered))]:
        path = tmp_path / f'{ident}.fasta'
        path.write_text('>original\n' + seq + '\n')
        members.append(Sample(ident, path, '', '', 'sp', 'group', 'local'))
    files = LineageFiles.in_directory(tmp_path)
    files.alignment.write_text('>q\n' + original + '\n>c\n' + ''.join(altered) + '\n')
    filtered, summary = run_verticall_mask(files, members, members[0].assembly, threads=1, force=False)
    assert summary['samples'] == 2
    assert summary['retained_clonal_sites'] > 0
    assert {r.id for r in SeqIO.parse(filtered, 'fasta')} == {'q', 'c'}
    assert all(set(r.seq) <= set('ACGT') for r in SeqIO.parse(filtered, 'fasta'))
