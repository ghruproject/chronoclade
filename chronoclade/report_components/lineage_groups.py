"""Authoritative typing groups for display, separate from exploratory distance cuts."""
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio import Phylo

from chronoclade.profile_types import _lineage
from chronoclade.provider_hiercc import reported_metadata, reported_scope


def assign_lineage_groups(analysis, records, *, level=7):
    """Assign reported lineage labels, or one neutral category in ST-only mode."""
    if isinstance(level, bool) or not isinstance(level, int) or not 5 <= level <= 10:
        raise ValueError("Display LIN level must be an integer between 5 and 10")
    partition = analysis.get('partition') or {}
    kind = partition.get('kind')
    if kind == 'profiles':
        # Numerical groups are descriptive distance cuts, never nomenclature
        # assignments. Provider codes remain scoped optional annotations.
        assignments = {}
        groups = []
        for group in analysis.get('genetic_groups', []):
            ids = group.get('sample_ids', [])
            label = 'Distance group ' + str(group['group_id'])
            assignments.update({ident: label for ident in ids})
            groups.append({**group, 'group_id': label, 'kind': 'profile_distance',
                           'method': 'exploratory complete-linkage profile distance group'})
        missing = [row['sample_id'] for row in records if row['sample_id'] not in assignments]
        if missing:
            label = 'Distance group unassigned'
            assignments.update({ident: label for ident in missing})
            groups.append({'group_id': label, 'sample_ids': missing, 'kind': 'profile_distance',
                           'method': 'no comparable distance group'})
        analysis['lineage_groups'] = groups
        analysis['display_group_basis'] = 'exploratory_profile_distance'
        analysis['display_group_level'] = None
        analysis['reported_lineage_annotations'] = [
            {'sample_id': row['sample_id'], 'role': row.get('role'),
             **{key: value for key, value in row.items()
                if key.startswith(('cglin_', 'hiercc_'))}}
            for row in records if row.get('cglin_raw') or row.get('hiercc_codes')
        ]
        return assignments
    if kind == 'hiercc':
        hiercc_level = partition.get('level')
        if not isinstance(hiercc_level, str) or not hiercc_level.startswith('HC') or not hiercc_level[2:].isdigit():
            raise ValueError('HierCC display requires the explicit HC level used for partitioning')
        members = defaultdict(list)
        namespaces = {}
        for row in records:
            assignment = _lineage(row, 'hiercc', hiercc_level)
            label = (f'HierCC {hiercc_level}: ' + str(assignment[2])
                     if assignment is not None else f'HierCC {hiercc_level}: unassigned')
            if assignment is not None:
                namespaces[label] = {'scheme': assignment[0],
                                     'scheme_version': row.get('hiercc_scheme_version') if reported_scope(row) else assignment[1],
                                     'code': assignment[2], **reported_metadata(row)}
            members[label].append(row['sample_id'])
        analysis['lineage_groups'] = [
            {'group_id': label, 'sample_ids': ids, 'method': 'reported HierCC assignment',
             'kind': 'hiercc', 'level': hiercc_level, **namespaces.get(label, {})}
            for label, ids in sorted(members.items())]
        analysis['display_group_basis'] = 'reported_hiercc'
        analysis['display_group_level'] = hiercc_level
        return {ident: label for label, ids in members.items() for ident in ids}
    if kind == 'st':
        label = 'Lineage unassigned'
        ids = [row['sample_id'] for row in records]
        analysis['lineage_groups'] = [
            {'group_id': label, 'sample_ids': ids, 'method': 'lineage not used in ST-only mode',
             'kind': 'unassigned'}
        ] if ids else []
        analysis['display_group_basis'] = 'unassigned_lineage'
        analysis['display_group_level'] = None
        return {ident: label for ident in ids}
    if kind != 'cglin':
        return None
    members = defaultdict(list)
    namespaces = {}
    for row in records:
        assignment = _lineage(row, 'cglin', level)
        if assignment is None:
            label = f'cgLIN {level}: unassigned'
        else:
            label = f'cgLIN {level}: ' + '.'.join(str(x) for x in assignment[2])
            namespaces[label] = {'scheme': assignment[0], 'scheme_version': assignment[1],
                                 'prefix': list(assignment[2])}
        members[label].append(row['sample_id'])
    analysis['lineage_groups'] = [
        {'group_id': label, 'sample_ids': ids, 'method': 'reported cgLIN prefix',
         'kind': 'cglin', 'level': level, **namespaces.get(label, {})}
        for label, ids in sorted(members.items())]
    analysis['display_group_basis'] = 'reported_cglin'
    analysis['display_group_level'] = level
    return {ident: label for label, ids in members.items() for ident in ids}


def render_lineage_groups(analysis, records, directory, *, reuse_coordinates=False, level=7):
    """Refresh presentation only from saved tree and distance-coordinate evidence."""
    from .group_tree import write_group_tree
    from .ordination import write_ordination_views
    from .tree_ordination import ordination_bases
    assignments = assign_lineage_groups(analysis, records, level=level)
    if assignments is None:
        return
    directory = Path(directory)
    by_id = {row['sample_id']: row for row in records}
    for cohort in analysis.get('cohorts') or []:
        if not cohort.get('tree_path'):
            continue
        rows = [by_id[ident] for ident in cohort['sample_ids']]
        def local_path(value):
            value = Path(value)
            return value if value.is_absolute() else directory / value
        with local_path(cohort.get('direct_pcoa_csv', cohort['pcoa_csv'])).open(newline='') as handle:
            direct = {row['sample_id']: [float(row['axis_1']), float(row['axis_2'])]
                      for row in csv.DictReader(handle)}
        tree = Phylo.read(local_path(cohort['tree_path']), 'newick')
        direct_coords = np.array([direct[row['sample_id']] for row in rows])
        if reuse_coordinates:
            with local_path(cohort['pcoa_csv']).open(newline='') as handle:
                saved = {row['sample_id']: [float(row['axis_1']), float(row['axis_2'])]
                         for row in csv.DictReader(handle)}
            coords = np.array([saved[row['sample_id']] for row in rows])
            bases = {
                'nj_tree': {'coordinates': coords, 'label': 'NJ-tree distances',
                            'axis_labels': ('MDS axis 1', 'MDS axis 2'),
                            'note': 'Saved NJ-tree distance projection. Negative branches were set to zero for ordination only; the saved tree is unchanged. The two-dimensional projection is approximate.'},
                'direct': {'coordinates': direct_coords, 'label': 'Direct profile distances',
                           'axis_labels': ('PCoA axis 1', 'PCoA axis 2'),
                           'note': 'Saved direct pairwise profile distance projection, before NJ fitting. Changing projection does not redefine the displayed groups.'},
            }
        else:
            coords, _, bases = ordination_bases(tree, rows, direct_coords)
        stem = cohort['cohort_id']
        cohort['tree_groups_html'] = write_group_tree(tree, rows, analysis['sample_labels'], assignments,
            directory / f'{stem}_tree_groups.html', distance_units=analysis.get('distance_units'),
            group_basis=analysis.get('display_group_basis'),
            group_threshold=analysis.get('distance_threshold'))
        cohort.update(write_ordination_views(directory, stem, coords, rows, groups=assignments,
            title='NJ-tree distance ordination', coordinate_views=bases, default_basis='nj_tree'))
