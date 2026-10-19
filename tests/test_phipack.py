"""PHIPack masks require valid tests and exact physical-coordinate accounting."""

import csv
import json
import random
import shutil

from Bio import SeqIO
import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.errors import WorkflowError
from chronoclade.native_pipeline import LineageFiles
from chronoclade.phipack import adjust_fdr, read_profile, run_profile_mask


def test_uninitialised_profile_values_are_never_used(tmp_path):
    path = tmp_path / 'Profile.csv'
    path.write_text('500,0\n1000,4.6e-310\n')
    log = tmp_path / 'profile.log'
    log.write_text('Checking for recombination at 500 between 0 and 1000\n'
                   'Too few informative sites to use normal approximation.\nPHI (Normal): --\n'
                   'Checking for recombination at 1000 between 500 and 1500\nPHI (Normal): -nan\n')
    result = read_profile(path, log, starts=[0, 500], window=1000)
    assert [r['p_value'] for r in result] == [None, None]


@pytest.mark.parametrize('value', ['0.2', 'nan', '-0.1'])
def test_profile_probability_disagreement_refused(tmp_path, value):
    path = tmp_path / 'Profile.csv'
    path.write_text(f'500,{value}\n')
    log = tmp_path / 'profile.log'
    log.write_text('Checking for recombination at 500 between 0 and 1000\nPHI (Normal): 1.00e-02\n')
    with pytest.raises(WorkflowError, match='probabilities'):
        read_profile(path, log, starts=[0], window=1000)


def test_fdr_uses_all_tested_windows_without_treating_missing_as_zero():
    rows = [{'p_value': p} for p in [None, 0.5, 0.04, 0.001]]
    adjust_fdr(rows)
    assert [r['q_value'] for r in rows] == [None, 0.5, 0.06, 0.003]


def test_mask_boundaries_final_tail_and_incomplete_columns(tmp_path, monkeypatch):
    files = LineageFiles.in_directory(tmp_path)
    seq = 'A' * 550 + 'N' + 'A' * 499
    files.alignment.write_text('>a\n' + seq + '\n>b\n' + 'A' * 1050 + '\n')
    reference = tmp_path / 'reference.fasta'
    reference.write_text('>chr\n' + 'A' * 600 + '\n>plasmid\n' + 'A' * 450 + '\n')

    def native(command, *, cwd, log, **kwargs):
        width = int(command[command.index('-n') + 1])
        step = int(command[command.index('-m') + 1])
        size = len(next(SeqIO.parse(command[command.index('-f') + 1], 'fasta')))
        values, messages = [], []
        for begin in range(0, size - width + 1, step):
            center = begin + width // 2
            p = 0.001 if cwd.name == 'block-00000' and begin == 0 else 0.8
            values.append(f'{center},{p}\n')
            messages.append(f'Checking for recombination at {center} between {begin} and {begin+width}\nPHI (Normal): {p}\n')
        (cwd / 'Profile.csv').write_text(''.join(values))
        log.write_text(''.join(messages))

    monkeypatch.setattr('chronoclade.native_pipeline._run_command', native)
    filtered, summary = run_profile_mask(files, reference, threads=2, force=False, window=400, step=200)
    assert summary['recombination_removed_sites'] == 400
    assert summary['incomplete_only_removed_sites'] == 1
    assert summary['retained_clonal_sites'] == 649
    assert summary['intervals'] == [{'start': 1, 'end': 400, 'reference_record': 'chr'}]
    assert all(len(r) == 649 for r in SeqIO.parse(filtered, 'fasta'))
    with (tmp_path / 'phipack/windows.csv').open() as stream:
        windows = list(csv.DictReader(stream))
    assert [(r['reference_record'], r['record_start_1based'], r['record_end_1based']) for r in windows] == [
        ('chr', '1', '400'), ('chr', '201', '600'), ('plasmid', '1', '400'), ('plasmid', '51', '450')]


@pytest.mark.integration
@pytest.mark.skipif(not shutil.which('Profile'), reason='PHIPack Profile not installed')
@pytest.mark.parametrize('recombined', [False, True])
def test_installed_profile_filters_synthetic_boundary_and_retains_untestable_alignment(tmp_path, recombined):
    rng = random.Random(10)
    sequences = [[] for _ in range(8)]
    for pos in range(20000):
        first = rng.choice('ACGT')
        other = rng.choice([c for c in 'ACGT' if c != first])
        pattern = [0, 0, 0, 0, 1, 1, 1, 1] if pos < 10000 else [0, 0, 1, 1, 0, 0, 1, 1]
        for i, group in enumerate(pattern):
            sequences[i].append(other if recombined and pos % 8 == 0 and group else first)
    files = LineageFiles.in_directory(tmp_path)
    files.alignment.write_text(''.join(f'>sample{i}\n' + ''.join(seq) + '\n' for i, seq in enumerate(sequences)))
    reference = tmp_path / 'reference.fasta'
    reference.write_text('>chromosome\n' + ''.join(sequences[0]) + '\n')
    filtered, summary = run_profile_mask(files, reference, threads=1, force=False)
    if recombined:
        assert summary['intervals'] == [{'start': 5001, 'end': 15000, 'reference_record': 'chromosome'}]
        assert summary['significant_windows'] == 1
        assert {len(r) for r in SeqIO.parse(filtered, 'fasta')} == {10000}
    else:
        assert summary['tested_windows'] == 0
        assert summary['untested_windows'] == 3
        assert summary['recombination_removed_sites'] == 0
        assert filtered.read_text() == files.alignment.read_text()


def test_cli_defaults_to_phipack_and_passes_explicit_options(tmp_path, monkeypatch):
    from test_tree_stage import selection_fixture

    selection, _ = selection_fixture(tmp_path)
    received = []

    def run(*args, **kwargs):
        received.append(kwargs)
        return {'temporal_assessment': {'code': 'not_assessed'}}

    monkeypatch.setattr('chronoclade.tree_stage.workflow.run_tree', run)
    for options in [[], ['--recombination-method', 'clonalframeml'],
                    ['--phipack-window', '20000', '--phipack-step', '10000', '--phipack-threshold', '0.01']]:
        result = CliRunner().invoke(app, ['tree', str(selection), *options])
        assert result.exit_code == 0, result.output
    assert [r['recombination_method'] for r in received] == ['phipack', 'clonalframeml', 'phipack']
    assert received[2]['phipack_window'] == 20000
    assert received[2]['phipack_step'] == 10000
    assert received[2]['phipack_threshold'] == 0.01


def test_slurm_saves_method_and_rejects_invalid_mask_settings(tmp_path):
    from test_tree_stage import selection_fixture
    from chronoclade.execution.jobs import SlurmResources
    from chronoclade.execution.slurm import write_slurm_ensemble

    _, ensemble = selection_fixture(tmp_path)
    kwargs = dict(results=tmp_path / 'results', resources=SlurmResources(2, 4, '00:10:00'))
    with pytest.raises(WorkflowError, match='step'):
        write_slurm_ensemble('tree', ensemble, output=tmp_path / 'bad', options={'phipack_window': 400, 'phipack_step': 500}, **kwargs)
    assert not (tmp_path / 'bad').exists()
    plan = write_slurm_ensemble('tree', ensemble, output=tmp_path / 'jobs', **kwargs)
    assert json.loads(plan.read_text())['options']['recombination_method'] == 'phipack'


def test_filtered_tree_is_reestimated_and_clonalframeml_is_not_called(tmp_path, monkeypatch):
    from chronoclade.metadata import Sample
    from chronoclade.native_pipeline import _run_core_phylogeny

    files = LineageFiles.in_directory(tmp_path)
    members = [Sample(i, tmp_path / f'{i}.fasta', '', '', 'sp', 'group', 'local') for i in 'abc']
    alignment = ''.join(f'>{i}\n' + 'A' * 1000 + '\n' for i in 'abc')
    filtered = tmp_path / 'phipack.filtered.fasta'
    commands = []

    def native(command, **kwargs):
        commands.append(command)
        if command[:2] == ['ska', 'map']:
            files.alignment.write_text(alignment)
        elif command[0] == 'iqtree':
            prefix = command[command.index('-pre') + 1]
            from pathlib import Path
            Path(prefix + '.treefile').write_text('(a:0.1,b:0.1,c:0.1);')
        elif command[0] != 'ska':
            pytest.fail(f'Unexpected native program: {command[0]}')

    def mask(*args, **kwargs):
        filtered.write_text(''.join(f'>{i}\n' + 'A' * 600 + '\n' for i in 'abc'))
        return filtered, {'method': 'phipack'}

    monkeypatch.setattr('chronoclade.native_pipeline._run_command', native)
    monkeypatch.setattr('chronoclade.phipack.run_profile_mask', mask)
    tree, used_alignment, summary = _run_core_phylogeny(files, members, members[0], 2, False,
                                                       recombination_method='phipack')
    inferred = [c for c in commands if c[0] == 'iqtree']
    assert len(inferred) == 2
    assert inferred[0][inferred[0].index('-s') + 1] == str(files.alignment)
    assert inferred[1][inferred[1].index('-s') + 1] == str(filtered)
    assert tree.name == 'phipack_tree.treefile'
    assert used_alignment == filtered and summary['method'] == 'phipack'
    assert not files.importations.exists() and not files.tree.exists()


def test_method_change_invalidates_reuse_and_is_saved(tmp_path, monkeypatch):
    from test_tree_stage import selection_fixture, native_stubs
    from chronoclade.tree_stage.workflow import run_tree

    calls = native_stubs.__wrapped__(monkeypatch)
    selection, _ = selection_fixture(tmp_path, dates=('', '', ''))
    first = run_tree(selection, output=tmp_path / 'tree')
    second = run_tree(selection, output=tmp_path / 'tree', recombination_method='clonalframeml')
    assert first['parameters']['recombination_method'] == 'phipack'
    assert second['recombination']['method'] == 'clonalframeml'
    assert first['fingerprint'] != second['fingerprint']
    assert len(calls) == 2
