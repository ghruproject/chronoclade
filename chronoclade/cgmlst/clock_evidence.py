"""Load and validate dense or binary distance evidence for matched clocks."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from chronoclade.errors import WorkflowError
from chronoclade.matrix_distances import BinaryDistanceEvidence
from chronoclade.selection_manifest import _read
from chronoclade.selections import DistanceEvidence


def _read_clock_distance_evidence(path, source):
    """Load the validated dense pair records or binary matrix named by an ensemble."""
    data = _read(path)
    if not isinstance(data, dict):
        raise WorkflowError('Distance evidence must be a JSON object')
    if (data.get('method') != source.get('method')
            or data.get('distance_definition') != source.get('distance_definition')):
        raise WorkflowError('Distance evidence method and definition disagree with the ensemble')
    schema = data.get('schema')
    if schema == 'chronoclade.binary_distances':
        evidence = BinaryDistanceEvidence(path)
        if evidence.descriptor != data:
            raise WorkflowError('Binary distance descriptor disagrees with its manifest')
        return evidence
    if schema not in (None, ''):
        raise WorkflowError(f'Unsupported distance evidence schema: {schema!r}')
    records = data.get('records')
    if not isinstance(records, list) or any(not isinstance(row, dict) for row in records):
        raise WorkflowError('Dense distance evidence requires pair records')
    return DistanceEvidence(
        tuple(dict(row) for row in records),
        method=data['method'],
        distance_definition=data['distance_definition'],
    )


def _distance_sample_ids(evidence):
    if isinstance(evidence, BinaryDistanceEvidence):
        return set(evidence.sample_ids)
    if not isinstance(evidence, DistanceEvidence):
        raise WorkflowError('Unsupported matched-clock distance evidence adapter')
    identifiers = []
    for row in evidence.records:
        left, right = row.get('sample_id_1'), row.get('sample_id_2')
        if any(not isinstance(name, str) or not name for name in (left, right)):
            raise WorkflowError('Dense distance evidence requires nonempty pair identifiers')
        identifiers.extend((left, right))
    return set(identifiers)


def _dense_distance_matrix(evidence, ids, expected_loci):
    """Adapt validated dense pairs to the exact matrix interface used by NJ."""
    if (evidence.method not in {'alleles', 'cgmlst'}
            or evidence.distance_definition != 'fraction of mismatching jointly called cgMLST loci'):
        raise WorkflowError('Unsupported dense cgMLST distance method or definition')
    if type(expected_loci) is not int or expected_loci < 1:
        raise WorkflowError('Dense distance evidence requires one prepared locus catalogue')
    sample_ids = []
    for row in evidence.records:
        for field in ('sample_id_1', 'sample_id_2'):
            name = row.get(field)
            if not isinstance(name, str) or not name:
                raise WorkflowError('Dense distance evidence requires nonempty pair identifiers')
            if name not in sample_ids:
                sample_ids.append(name)
    if len(ids) != len(set(ids)) or not set(ids).issubset(sample_ids):
        raise WorkflowError('Selected identifiers must be unique and present in full-pool distances')
    position = {name: index for index, name in enumerate(sample_ids)}
    distance = np.full((len(sample_ids), len(sample_ids)), np.nan, dtype=np.float64)
    np.fill_diagonal(distance, 0)
    scheme_loci = np.full((len(sample_ids), len(sample_ids)), expected_loci, dtype=np.uint32)
    seen = set()

    def integer(row, field):
        value = row.get(field)
        if isinstance(value, bool):
            raise WorkflowError(f'Dense distance evidence has an invalid {field}')
        try:
            number = float(value)
        except (TypeError, ValueError):
            raise WorkflowError(f'Dense distance evidence has an invalid {field}') from None
        if not np.isfinite(number) or number < 0 or not number.is_integer():
            raise WorkflowError(f'Dense distance evidence has an invalid {field}')
        return int(number)

    for row in evidence.records:
        left, right = row['sample_id_1'], row['sample_id_2']
        if left == right:
            raise WorkflowError('Dense distance evidence must not contain self-pairs')
        key = tuple(sorted((left, right)))
        if key in seen:
            raise WorkflowError(f'Duplicate dense distance evidence for pair {key!r}')
        seen.add(key)
        raw = integer(row, 'allele_differences')
        shared = integer(row, 'shared_called_loci')
        recorded_loci = row.get('scheme_loci')
        if recorded_loci not in (None, '') and integer(row, 'scheme_loci') != expected_loci:
            raise WorkflowError('Dense distance evidence disagrees with the prepared locus catalogue')
        if raw > shared or shared > expected_loci:
            raise WorkflowError('Dense distance counts exceed the prepared locus catalogue')
        overlap = row.get('call_overlap')
        if overlap not in (None, ''):
            try:
                overlap = float(overlap)
            except (TypeError, ValueError):
                raise WorkflowError('Dense distance evidence has an invalid call_overlap') from None
            if (not np.isfinite(overlap) or not 0 <= overlap <= 1
                    or not np.isclose(overlap, shared / expected_loci, rtol=0, atol=1e-12)):
                raise WorkflowError('Dense call_overlap disagrees with the prepared locus catalogue')
        raw_distance = row.get('distance')
        if raw_distance in (None, ''):
            value = np.nan
        else:
            try:
                value = float(raw_distance)
            except (TypeError, ValueError):
                raise WorkflowError('Dense distance evidence has an invalid distance') from None
        if np.isinf(value) or (np.isfinite(value) and not 0 <= value <= 1):
            raise WorkflowError('Dense distance values must be between zero and one or unavailable')
        if np.isfinite(value):
            if shared == 0 or not np.isclose(value, raw / shared, rtol=0, atol=1e-12):
                raise WorkflowError('Dense distances disagree with mismatch/callable counts')
        i, j = position[left], position[right]
        distance[i, j] = distance[j, i] = value
    return SimpleNamespace(
        sample_ids=tuple(sample_ids),
        arrays={'distance': distance, 'scheme_loci': scheme_loci},
    )
