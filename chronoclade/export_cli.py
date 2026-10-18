"""Checked local viewer exports from saved module manifests."""

from enum import Enum
import json
from pathlib import Path
import re
import tempfile
from typing import Annotated

from Bio import Phylo
import typer
from rich.console import Console

from chronoclade.artifacts import file_sha256, write_json
from chronoclade.context_partitions import profile_records
from chronoclade.datasets import load_dataset
from chronoclade.errors import WorkflowError
from chronoclade.interoperability import write_grapetree_exports, write_microreact_exports
from chronoclade.selection_manifest import _read, _resolve


class Viewer(str, Enum):
    grapetree = "grapetree"
    microreact = "microreact"


def _choose(items, value, key, option):
    if not isinstance(items, list) or any(not isinstance(item, dict) for item in items):
        raise WorkflowError(f'Manifest {option} entries must be a list of objects')
    matches = [item for item in items if value is None or item.get(key) == value]
    if len(matches) != 1:
        available = ', '.join(str(item.get(key)) for item in items)
        raise WorkflowError(f"Choose {option}; available values: {available or 'none'}")
    return matches[0]


def _profile_source(path, manifest, block, cohort):
    schema = manifest.get('schema')
    if schema == 'chronoclade.prepared-dataset':
        if block is not None or cohort is not None:
            raise WorkflowError('--block/--cohort apply only to a cgmlst.json manifest')
        dataset = load_dataset(path)
        return dataset, profile_records(dataset), None, {}
    if (schema != 'chronoclade.cgmlst.analysis' or manifest.get('schema_version') != 1
            or manifest.get('status') != 'complete'):
        raise WorkflowError('Profile export requires a complete dataset.json or cgmlst.json')
    selected = _choose(manifest.get('blocks', []), block, 'block_id', '--block')
    block_id = selected['block_id']
    if not isinstance(block_id, str) or not re.fullmatch(r'[A-Za-z0-9_-]+', block_id):
        raise WorkflowError('Invalid block directory identity')
    directory = path.parent / block_id
    checkpoint_path = directory / 'numerical_analysis.json'
    checkpoint = _read(checkpoint_path)
    if not isinstance(checkpoint.get('analysis'), dict) or not isinstance(
        checkpoint.get('artifacts'), dict
    ):
        raise WorkflowError('Numerical checkpoint must contain analysis and artifact objects')
    selected_cohort = _choose(checkpoint['analysis']['cohorts'], cohort, 'cohort_id', '--cohort')
    ids = selected_cohort['sample_ids']
    if not isinstance(ids, list) or any(not isinstance(ident, str) for ident in ids):
        raise WorkflowError('Cohort sample identifiers must be a list of strings')
    dataset = load_dataset(_resolve(path.parent, manifest['dataset'], 'dataset', external=True))
    if dataset.dataset_id != manifest['dataset_id']:
        raise WorkflowError('Prepared dataset identity differs from cgMLST analysis')
    records = {row['sample_id']: row for row in profile_records(dataset)}
    if len(set(ids)) != len(ids) or set(ids) - records.keys():
        raise WorkflowError('Cohort identities differ from prepared data')
    relative = selected_cohort.get('tree_path')
    if not relative:
        raise WorkflowError('The selected cohort has no usable NJ tree')
    tree = _resolve(directory, {'path': relative,
        'sha256': checkpoint['artifacts'].get(relative)}, 'saved NJ tree')
    return dataset, [records[ident] for ident in ids], tree, {
        'block_id': block_id, 'cohort_id': selected_cohort['cohort_id'],
        'checkpoint_sha256': file_sha256(checkpoint_path),
    }


def export_saved(source, output, viewer, *, block=None, cohort=None):
    """Publish a new export directory after all identity and source checks pass."""
    path = Path(source).expanduser().resolve()
    output = Path(output).expanduser().resolve()
    viewer = Viewer(viewer)
    if output.exists():
        raise WorkflowError('Export output already exists; choose a new directory')
    manifest = _read(path)
    digest = file_sha256(path)
    schema = manifest.get('schema')
    audit, tree, dataset = {}, None, None
    if schema in {'chronoclade.tree', 'chronoclade.time'}:
        if viewer != Viewer.microreact or block is not None or cohort is not None:
            raise WorkflowError('Tree/time manifests support Microreact; profile options do not apply')
        from chronoclade.stage_artifacts import load_result, load_tree_result

        stage = schema.rsplit('.', 1)[-1]
        result = load_tree_result(path) if stage == 'tree' else load_result(path, 'time')
        reference = result.get('corrected_tree') if stage == 'tree' else result.get('dated_tree')
        if not reference:
            raise WorkflowError('No dated tree is available; export the saved tree-stage result instead')
        tree = _resolve(path.parent, reference, 'viewer tree')
        metadata = result.get('metadata') or result['artifacts'].get('selected_metadata.json')
        rows = json.loads(_resolve(path.parent, metadata, 'metadata').read_text())
        audit['branch_units'] = result['branch_units']
    else:
        dataset, rows, tree, audit = _profile_source(path, manifest, block, cohort)
        if viewer == Viewer.microreact and tree is None:
            raise WorkflowError('Microreact needs a tree; supply cgmlst.json, tree.json or time.json')
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.chronoclade-export-', dir=output.parent) as temporary:
        working = Path(temporary)
        target = working / 'export'
        if viewer == Viewer.grapetree:
            ids = {row['sample_id'] for row in rows}
            matrices = [matrix for matrix in dataset.profiles if ids <= set(matrix.sample_ids)]
            if len(matrices) != 1:
                raise WorkflowError('GrapeTree requires one profile catalogue covering the selected records')
            catalogue = matrices[0].catalogue
            audit.update(locus_universe_complete=catalogue.complete,
                         cgmlst_scheme=catalogue.scheme_id, locus_count=len(catalogue.loci))
            paths = write_grapetree_exports(rows, target, loci=catalogue.loci)
        else:
            if schema == 'chronoclade.time':
                from chronoclade.stage_artifacts import read_dated_tree

                converted = working / 'dated.newick'
                Phylo.write(read_dated_tree(tree), converted, 'newick', format_branch_length='%.17g')
                tree = converted
            paths = write_microreact_exports(tree, rows, target)
        if file_sha256(path) != digest:
            raise WorkflowError('Source manifest changed during export')
        evidence = {'schema': 'chronoclade.viewer-export', 'schema_version': 1,
            'viewer': viewer.value, 'sample_count': len(rows), 'source_sha256': digest,
            'source_schema': schema, 'scope': audit,
            'artifacts': {key: {'path': value.relative_to(target).as_posix(),
                               'sha256': file_sha256(value)} for key, value in paths.items()}}
        write_json(target / 'export.json', evidence)
        target.rename(output)
    return evidence


def export_command(
    source: Annotated[Path, typer.Argument(help='Saved dataset, cgMLST, tree or time manifest')],
    viewer: Annotated[Viewer, typer.Option(help='Target viewer')],
    output: Annotated[Path, typer.Option('--out', '--output', '-o')],
    block: Annotated[str | None, typer.Option(help='Explicit block when cgMLST has several')] = None,
    cohort: Annotated[str | None, typer.Option(help='Explicit cohort when a block has several')] = None,
):
    """Export checked local files for GrapeTree or Microreact without uploading."""
    try:
        evidence = export_saved(source, output, viewer, block=block, cohort=cohort)
    except (OSError, ValueError, KeyError, WorkflowError) as error:
        Console().print(f'Error: {error}', style='bold red', markup=False)
        raise typer.Exit(2) from error
    Console().print(f'{viewer.value}: {evidence["sample_count"]} genomes; files in {output.resolve()}')
    Console().print(f'Open {output.resolve() / "README.txt"} for viewer import instructions.')
