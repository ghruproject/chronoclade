"""Ordinate NJ patristic distances with an explicit negative-branch policy."""
from pathlib import Path

import numpy as np


def nj_distances(tree, sample_ids, *, output=None, batch_size=256):
    """Sum nonnegative branch lengths between tips, without mutating the tree.

    Cross-child pairs are assigned once at their lowest common ancestor. A
    memory-mapped output keeps large blocks from needing another resident matrix.
    """
    names = list(sample_ids)
    terminals = tree.get_terminals()
    if (len(set(names)) != len(names) or len(terminals) != len(names)
            or {tip.name for tip in terminals} != set(names)):
        raise ValueError('NJ ordination requires unique IDs matching every tree tip')
    n = len(names)
    if output is None:
        matrix = np.zeros((n, n), dtype=float)
    else:
        matrix = np.lib.format.open_memmap(Path(output), mode='w+', dtype='float64', shape=(n, n))
        matrix[:] = 0
    indexes = {name: i for i, name in enumerate(names)}
    depths = {tree.root: 0.0}
    negative = []
    for node in tree.find_clades(order='preorder'):
        for child in node.clades:
            length = float(child.branch_length or 0)
            if not np.isfinite(length):
                raise ValueError('NJ ordination requires finite branch lengths')
            if length < 0:
                negative.append(length)
            depths[child] = depths[node] + max(0.0, length)
    tip_depths = np.array([depths[tip] for tip in sorted(terminals, key=lambda t: indexes[t.name])])
    descendants = {}
    for node in tree.find_clades(order='postorder'):
        if not node.clades:
            descendants[node] = np.array([indexes[node.name]], dtype=int)
            continue
        seen = np.array([], dtype=int)
        for child in node.clades:
            current = descendants.pop(child)
            for start in range(0, len(seen), batch_size):
                left = seen[start:start + batch_size]
                distances = np.maximum(0, tip_depths[left, None] + tip_depths[current] - 2 * depths[node])
                matrix[np.ix_(left, current)] = distances
                matrix[np.ix_(current, left)] = distances.T
            seen = np.concatenate((seen, current))
        descendants[node] = seen
    if isinstance(matrix, np.memmap):
        matrix.flush()
    return matrix, {
        'basis': 'nj_tree',
        'negative_branch_policy': 'set negative branches to zero for ordination only; saved NJ tree unchanged',
        'negative_branch_count': len(negative),
        'negative_branch_total_magnitude': -sum(negative),
    }


def pcoa(matrix, *, seed=42):
    """Classical MDS; use the existing iterative backend for large blocks."""
    if len(matrix) > 1500:
        from chronoclade.cgmlst.scale_tree import leading_pcoa
        return leading_pcoa(matrix, seed=seed)
    n = len(matrix)
    centered = np.eye(n) - np.ones((n, n)) / n
    values, vectors = np.linalg.eigh(-.5 * centered @ (matrix ** 2) @ centered)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    positive = np.maximum(values, 0)
    coordinates = np.zeros((n, 2))
    count = min(2, n)
    coordinates[:, :count] = vectors[:, :count] * np.sqrt(positive[:count])
    for column in range(count):
        anchor = np.argmax(np.abs(coordinates[:, column]))
        if coordinates[anchor, column] < 0:
            coordinates[:, column] *= -1
    return coordinates, {
        'pcoa_eigenvalues': values.tolist(),
        'pcoa_positive_axis_fraction': (positive[:2] / positive.sum()).tolist() if positive.sum() else [0, 0],
        'pcoa_full_spectrum_computed': True,
        'pcoa_solver': 'dense symmetric eigendecomposition of centered squared NJ-tree distances',
        'pcoa_inertia_note': 'Axis fractions use the total positive eigenvalue sum; negative eigenvalues are retained in evidence.',
    }


def ordination_bases(tree, records, direct_coordinates, *, output=None, seed=42):
    distances, audit = nj_distances(tree, [row['sample_id'] for row in records], output=output)
    coordinates, diagnostics = pcoa(distances, seed=seed)
    audit.update(diagnostics)
    count = audit['negative_branch_count']
    note = ('Distances sum branches between tips of the NJ tree. '
            f'{count} negative branches were set to zero for this ordination only; the saved NJ tree is unchanged. '
            'The two-dimensional projection is approximate. Changing projection does not redefine the displayed groups.')
    return coordinates, audit, {
        'nj_tree': {'coordinates': coordinates, 'label': 'NJ-tree distances',
                    'axis_labels': ('MDS axis 1', 'MDS axis 2'), 'note': note},
        'direct': {'coordinates': direct_coordinates, 'label': 'Direct profile distances',
                   'axis_labels': ('PCoA axis 1', 'PCoA axis 2'),
                   'note': 'Direct pairwise profile distances, before NJ fitting. This view can differ from the tree projection.'},
    }
