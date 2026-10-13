"""Authoritative typing groups for display, separate from exploratory distance cuts."""
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio import Phylo

from chronoclade.context_refinement import _lineage


def assign_lineage_groups(analysis, records):
    """Use reported LIN prefixes; never infer a code from distances or neighbours."""
    partition = analysis.get('partition') or {}
    if partition.get('kind') != 'cglin':
        return None
    level = 7
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


def render_lineage_groups(analysis, records, directory):
    """Refresh presentation only from saved tree and distance-coordinate evidence."""
    from .group_tree import write_group_tree
    from .ordination import write_ordination_views
    from .tree_ordination import ordination_bases
    assignments = assign_lineage_groups(analysis, records)
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
        coords, _, bases = ordination_bases(tree, rows, np.array([direct[row['sample_id']] for row in rows]))
        stem = cohort['cohort_id']
        cohort['tree_groups_html'] = write_group_tree(tree, rows, analysis['sample_labels'], assignments,
            directory / f'{stem}_tree_groups.html', distance_units=analysis.get('distance_units'))
        cohort.update(write_ordination_views(directory, stem, coords, rows, groups=assignments,
            title='NJ-tree distance ordination', coordinate_views=bases, default_basis='nj_tree'))
