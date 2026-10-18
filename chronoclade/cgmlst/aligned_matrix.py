"""Memory-bounded all-pairs evidence from per-locus nucleotide alignments.

Allele vocabulary indices are categorical. Missing calls and unavailable DNA map
into a zero-evidence category, never a matching allele. Substitutions are counted
only at jointly A/C/G/T positions; aligned gaps are a separate measurement.
"""

from dataclasses import dataclass
import math
import mmap
from pathlib import Path
import time

import numpy as np

from chronoclade.cgmlst.locus_alignment import file_hash, save_json


@dataclass
class AlignedMatrix:
    """On-disk square evidence arrays, including measured diagonal evidence."""

    sample_ids: tuple
    snps: np.ndarray
    sites: np.ndarray
    joint_loci: np.ndarray
    gap_columns: np.ndarray
    resolved_alleles: np.ndarray
    provenance: dict
    eligibility: np.ndarray | None = None

    def eligible(self, min_loci=None, rows=slice(None)):
        """Return evidence eligibility; zero jointly callable sites never qualify."""
        minimum = self.provenance["minimum_compared_loci"] if min_loci is None else min_loci
        if not 1 <= minimum <= self.provenance["canonical_loci"]:
            raise ValueError("Minimum compared loci lies outside the canonical scheme")
        return (self.joint_loci[rows] >= minimum) & (self.sites[rows] > 0)

    def mismatch_fraction(self, rows=slice(None)):
        """Compute requested rows lazily; no-evidence cells are NaN."""
        sites = self.sites[rows]
        result = np.full(sites.shape, np.nan, dtype=np.float64)
        np.divide(self.snps[rows], sites, out=result, where=sites > 0)
        return result


def _locus_table(matrix, column, locus, width, aligned_alleles):
    """Collapse identical DNA, then measure each aligned allele pair exactly."""
    codes = matrix.codes[:, column]
    present = np.unique(codes)
    mapping = np.zeros(len(matrix.categories) + 1, dtype=np.uint32)
    sequences, dna_codes = [], {}
    for code in present:
        if not code:
            continue
        sequence = aligned_alleles.get((locus, matrix.categories[int(code) - 1]))
        if sequence is None:
            continue
        if len(sequence) != width or set(sequence) - set("ACGTN-"):
            raise ValueError(f"Invalid aligned DNA for locus {locus}")
        if sequence not in dna_codes:
            sequences.append(sequence)
            dna_codes[sequence] = len(sequences)
        mapping[code] = dna_codes[sequence]
    encoded = mapping[codes]
    count = len(sequences) + 1
    tables = [np.zeros((count, count), dtype=dtype)
              for dtype in (np.uint32, np.uint32, np.uint16, np.uint32)]
    if sequences:
        vectors = np.array([list(seq.encode()) for seq in sequences], dtype=np.uint8)
        bases = (vectors == ord("A")) | (vectors == ord("C")) | (vectors == ord("G")) | (vectors == ord("T"))
        gaps = vectors == ord("-")
        for index, (sequence, called, gap) in enumerate(zip(vectors, bases, gaps, strict=True), 1):
            joint = called & bases
            sites = np.count_nonzero(joint, axis=1)
            tables[0][index, 1:] = np.count_nonzero((sequence != vectors) & joint, axis=1)
            tables[1][index, 1:] = sites
            tables[2][index, 1:] = sites > 0
            tables[3][index, 1:] = np.count_nonzero((gap & bases) | (called & gaps), axis=1)
    return encoded, [(index, table) for index, table in enumerate(tables) if np.any(table)]


def build_aligned_matrix(matrix, alignment, output, *, block_size=128, min_loci=None,
                         progress=None):
    """Accumulate exact per-locus lookup tables in bounded row blocks.

    The five persisted evidence arrays use 16 bytes per sample pair. Working arrays
    use approximately ``16 * block_size * n`` bytes, plus one lookup temporary
    and small allele tables. Iterating blocks outside loci avoids rewriting the
    full on-disk matrices for each locus. No concatenation or float NxN array is
    required. Symmetric evidence is computed in the upper block triangle and
    mirrored. ``progress(phase, completed, total)`` is optional.
    """
    matrix.validate()
    if type(block_size) is not int or block_size < 1:
        raise ValueError("Block size must be a positive integer")
    ids, loci = tuple(matrix.sample_ids), tuple(matrix.catalogue.loci)
    alignment_ids = tuple(alignment.sample_ids or alignment.provenance.get("sample_ids", ())
                          or alignment.samples)
    if ids != alignment_ids or loci != tuple(p["locus"] for p in alignment.partitions):
        raise ValueError("Allele matrix and alignment must have identical sample and locus order")
    minimum = math.ceil(0.95 * len(loci)) if min_loci is None else min_loci
    if type(minimum) is not int or not 1 <= minimum <= len(loci):
        raise ValueError("Minimum compared loci lies outside the canonical scheme")
    if len(loci) > np.iinfo(np.uint16).max:
        raise ValueError("Joint locus counts exceed uint16 capacity")
    widths = [p["width"] for p in alignment.partitions]
    if any(type(width) is not int or width < 0 for width in widths):
        raise ValueError("Alignment partition widths must be nonnegative integers")
    if sum(widths) > np.iinfo(np.uint32).max:
        raise ValueError("Alignment counts exceed uint32 capacity")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "matrix.json").unlink(missing_ok=True)
    started = time.perf_counter()
    tables = []
    for column, (locus, width) in enumerate(zip(loci, widths, strict=True)):
        if width:
            encoded, measures = _locus_table(matrix, column, locus, width, alignment.aligned_alleles)
            tables.append((column, encoded, measures))
        if progress and ((column + 1) % 25 == 0 or column + 1 == len(loci)):
            progress("locus_tables", column + 1, len(loci))
    names = ("snps", "sites", "joint_loci", "gap_columns", "resolved_alleles")
    dtypes = (np.uint32, np.uint32, np.uint16, np.uint32, np.uint16)
    arrays = [np.lib.format.open_memmap(output / f"{name}.npy", mode="w+", dtype=dtype,
                                      shape=(len(ids), len(ids)))
              for name, dtype in zip(names, dtypes, strict=True)]
    eligibility = np.lib.format.open_memmap(output / "eligible.npy", mode="w+", dtype=bool,
                                           shape=(len(ids), len(ids)))
    for start in range(0, len(ids), block_size):
        stop = min(start + block_size, len(ids))
        accumulators = [np.zeros((stop - start, len(ids) - start), dtype=dtype) for dtype in dtypes]
        for column, encoded, measures in tables:
            left, right = encoded[start:stop, None], encoded[None, start:]
            for index, measure in measures:
                measured = measure[left, right]
                accumulators[index] += measured
                if index == 2:
                    different = matrix.codes[start:stop, column, None] != matrix.codes[None, start:, column]
                    accumulators[4] += different & (measured > 0)
        for array, accumulator in zip(arrays, accumulators, strict=True):
            array[start:stop, start:] = accumulator
            array[start:, start:stop] = accumulator.T
        eligible = (accumulators[2] >= minimum) & (accumulators[1] > 0)
        eligibility[start:stop, start:] = eligible
        eligibility[start:, start:stop] = eligible.T
        # Persist before evicting shared mappings. Large outputs need not remain
        # resident merely because all their row blocks were written in this run.
        for array in [*arrays, eligibility]:
            array.flush()
            if hasattr(array._mmap, "madvise") and hasattr(mmap, "MADV_DONTNEED"):
                array._mmap.madvise(mmap.MADV_DONTNEED)
        if progress:
            progress("matrix_rows", stop, len(ids))
    for array in [*arrays, eligibility]:
        array.flush()
    artifacts = {}
    for index, (name, dtype) in enumerate(zip((*names, "eligible"), (*dtypes, bool), strict=True), 1):
        path = output / f"{name}.npy"
        artifacts[name] = {"path": path.name, "dtype": np.dtype(dtype).name,
                           "sha256": file_hash(path)}
        if progress:
            progress("matrix_hashes", index, len(names) + 1)
    metadata = {"schema": "chronoclade.cgmlst.aligned-matrix", "schema_version": 1,
                "status": "complete", "sample_ids": list(ids), "canonical_loci": len(loci),
                "alignment_sites": sum(widths), "minimum_compared_loci": minimum,
                "block_size": block_size, "elapsed_seconds": time.perf_counter() - started,
                "implementation_sha256": file_hash(__file__),
                "alignment_provenance": alignment.provenance,
                "policy": "Jointly ACGT substitutions and sites; loci require at least one joint ACGT site; missing calls and unavailable sequences contribute no evidence; gaps counted separately; resolved_alleles counts original categorical differences only for loci with joint ACGT evidence",
                "arrays": artifacts}
    save_json(output / "matrix.json", metadata)
    return AlignedMatrix(ids, *arrays, metadata, eligibility)
