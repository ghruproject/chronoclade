"""Branch support remains a separate, identity-checked inference."""

from io import StringIO
from pathlib import Path

from Bio import Phylo
import pytest

from chronoclade.errors import WorkflowError
from chronoclade.tree_stage.support import run_branch_support, splits, validate_replicates


@pytest.mark.parametrize('value', [-1, 1, 999, True, 1000.0])
def test_invalid_replicates(value):
    with pytest.raises(WorkflowError, match='at least 1000'):
        validate_replicates(value)


def test_split_keys_ignore_root_and_preserve_support():
    _, first = splits(Phylo.read(StringIO('((a,b)95,(c,d));'), 'newick'))
    _, second = splits(Phylo.read(StringIO('(a,(b,(c,d)95));'), 'newick'))
    assert first == second == {('a', 'b'): 95}


def test_support_inference_records_shared_and_absent_splits(tmp_path, monkeypatch):
    corrected = tmp_path/'corrected.nwk'
    corrected.write_text('((a,b),(c,(d,e)));')
    alignment = tmp_path/'filtered.fasta'
    alignment.write_text(''.join(f'>{name}\nACGT\n' for name in 'abcde'))
    original = corrected.read_bytes()
    commands = []

    def native(command, **kwargs):
        commands.append(command)
        tree, report = kwargs['expected']
        tree.write_text('((a,b)98,(d,(c,e)82));')
        report.write_text('UFBoot2')

    monkeypatch.setattr('chronoclade.lineage._run_command', native)
    result = run_branch_support(alignment=alignment, corrected_tree=corrected,
        directory=tmp_path, output_root=tmp_path, replicates=1000, threads=1, seed=42)
    assert corrected.read_bytes() == original
    assert result['shared_splits'] == 1
    assert result['corrected_splits'] == 2
    assert result['shared_splits_with_support'] == 1
    assert '--bnni' in commands[0]
    assert commands[0][commands[0].index('-seed')+1] == '42'
    assert 'False' in (tmp_path/'filtered_support_splits.csv').read_text()
    assert result['source_alignment']['sha256']


def test_wrong_tip_set_is_rejected(tmp_path, monkeypatch):
    corrected = tmp_path/'corrected.nwk'
    corrected.write_text('((a,b),(c,d));')
    alignment = tmp_path/'alignment.fasta'
    alignment.write_text('>a\nACGT\n')

    def native(command, **kwargs):
        tree, report = kwargs['expected']
        tree.write_text('((a,b),(c,wrong));')
        report.write_text('IQ-TREE')

    monkeypatch.setattr('chronoclade.lineage._run_command', native)
    with pytest.raises(WorkflowError, match='exact selected tip set'):
        run_branch_support(alignment=alignment, corrected_tree=corrected,
            directory=tmp_path, output_root=tmp_path, replicates=1000, threads=1, seed=42)


def test_disabled_needs_no_files(tmp_path):
    assert run_branch_support(alignment=Path('absent'), corrected_tree=Path('absent'),
        directory=tmp_path, output_root=tmp_path, replicates=0, threads=1, seed=42)['status'] == 'not_assessed'
