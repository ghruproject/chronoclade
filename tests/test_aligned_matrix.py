"""Exact scalar agreement and bounded reconstruction of aligned evidence."""

import json
import sys

import numpy as np
import pytest

from chronoclade.cgmlst.aligned_matrix import build_aligned_matrix
from chronoclade.cgmlst.locus_alignment import LocusAlignment, align_sequences, file_hash, load_alignment, read_fasta
from chronoclade.datasets.model import AlleleMatrix, LocusCatalogue


def panel():
    loci = ("variable", "fixed", "unavailable", "absent")
    profiles = {"a": {"variable": "1", "fixed": "1", "unavailable": "8"},
                "b": {"variable": "2", "fixed": "2", "unavailable": "9"},
                "c": {"variable": "3", "fixed": "1"}, "d": {}}
    matrix = AlleleMatrix.from_profiles(LocusCatalogue("test", "1", loci, "fixture"),
                                       tuple(profiles), profiles)
    alleles = {("variable", "1"): "AC--GTN", ("variable", "2"): "ATAC-TN",
               ("variable", "3"): "NN--NNN", ("fixed", "1"): "GG",
               ("fixed", "2"): "GG"}
    partitions = [{"locus": locus, "width": width} for locus, width in zip(loci, (7, 2, 3, 0))]
    alignment = LocusAlignment({}, partitions, {"sample_ids": list(profiles)}, alleles,
                               tuple(profiles))
    return matrix, alignment


def scalar(matrix, alignment, i, j):
    """Direct character oracle, independent of category tables and NumPy."""
    totals = [0, 0, 0, 0, 0]
    left, right = matrix.profile(matrix.sample_ids[i]), matrix.profile(matrix.sample_ids[j])
    for partition in alignment.partitions:
        locus = partition["locus"]
        a = alignment.aligned_alleles.get((locus, left[locus]))
        b = alignment.aligned_alleles.get((locus, right[locus]))
        if a is None or b is None:
            continue
        sites = 0
        for x, y in zip(a, b, strict=True):
            if x in "ACGT" and y in "ACGT":
                sites += 1
                totals[0] += x != y
            elif (x == "-" and y in "ACGT") or (y == "-" and x in "ACGT"):
                totals[3] += 1
        totals[1] += sites
        totals[2] += sites > 0
        totals[4] += sites > 0 and left[locus] != right[locus]
    return tuple(totals)


@pytest.mark.parametrize("block_size", [1, 3, 128])
def test_exact_scalar_evidence_with_missing_n_indels_and_identical_dna(tmp_path, block_size):
    matrix, alignment = panel()
    progress = []
    result = build_aligned_matrix(matrix, alignment, tmp_path, block_size=block_size,
                                  min_loci=2, progress=lambda *event: progress.append(event))
    for i in range(4):
        for j in range(4):
            assert tuple(int(array[i, j]) for array in
                         (result.snps, result.sites, result.joint_loci, result.gap_columns, result.resolved_alleles)) == scalar(matrix, alignment, i, j)
    assert result.resolved_alleles[0, 1] == 2  # Fixed locus keeps distinct allele labels.
    assert result.snps[0, 1] == 1
    assert result.sites[0, 1] == 5
    assert result.gap_columns[0, 1] == 3
    assert result.joint_loci[0, 1] == 2
    assert result.joint_loci[2, 2] == 1  # Ns and both-gap columns supply no joint base.
    assert result.eligibility[0, 1]
    assert not result.eligibility[2, 2]
    assert not result.eligibility[3].any()
    assert np.isnan(result.mismatch_fraction(rows=3)).all()
    assert result.mismatch_fraction(rows=0)[1] == 0.2
    assert np.array_equal(result.eligibility, result.eligible())
    assert isinstance(result.snps, np.memmap)
    assert np.array_equal(np.load(tmp_path / "snps.npy"), result.snps)
    assert json.loads((tmp_path / "matrix.json").read_text())["status"] == "complete"
    assert ("matrix_rows", 4, 4) in progress
    assert progress[-1] == ("matrix_hashes", 6, 6)
    for artifact in result.provenance["arrays"].values():
        assert file_hash(tmp_path / artifact["path"]) == artifact["sha256"]


def test_randomized_scalar_agreement(tmp_path):
    rng = np.random.default_rng(19)
    loci, ids = tuple(f"l{i}" for i in range(7)), tuple(f"s{i}" for i in range(19))
    profiles = {ident: {locus: str(rng.integers(0, 5)) for locus in loci} for ident in ids}
    matrix = AlleleMatrix.from_profiles(LocusCatalogue("test", "1", loci, "fixture"), ids, profiles)
    alleles = {(locus, str(code)): "".join(rng.choice(list("ACGTN-"), size=23))
               for locus in loci for code in range(1, 5)}
    alignment = LocusAlignment({}, [{"locus": locus, "width": 23} for locus in loci],
                               {}, alleles, ids)
    result = build_aligned_matrix(matrix, alignment, tmp_path, block_size=4)
    for i in range(len(ids)):
        for j in range(len(ids)):
            assert tuple(int(array[i, j]) for array in
                         (result.snps, result.sites, result.joint_loci, result.gap_columns, result.resolved_alleles)) == scalar(matrix, alignment, i, j)


@pytest.mark.parametrize("failure", ["order", "width", "dna", "block", "minimum"])
def test_invalid_alignment_evidence_fails_closed(tmp_path, failure):
    matrix, alignment = panel()
    kwargs = {}
    if failure == "order":
        alignment.sample_ids = tuple(reversed(alignment.sample_ids))
    elif failure == "width":
        alignment.partitions[0]["width"] = -1
    elif failure == "dna":
        alignment.aligned_alleles["variable", "1"] = "XXXXXXX"
    elif failure == "block":
        kwargs["block_size"] = 0
    else:
        kwargs["min_loci"] = 5
    with pytest.raises(ValueError):
        build_aligned_matrix(matrix, alignment, tmp_path, **kwargs)
    assert not (tmp_path / "matrix.json").exists()


@pytest.mark.parametrize("export", [True, False])
def test_streaming_alignment_matches_materialized_and_keeps_source_provenance(tmp_path, export):
    tool = tmp_path / "mafft"
    tool.write_text(f"#!{sys.executable}\nprint('fixture 1')\n")
    tool.chmod(0o755)
    rows = [{"sample_id": "a", "cgmlst_loci": ["x", "y"],
             "cgmlst_profile": {"x": "1"}},
            {"sample_id": "b", "cgmlst_loci": ["x", "y"],
             "cgmlst_profile": {}}]
    kwargs = {"sequences": {("x", "1"): "ACGT"}, "records": rows, "mafft": tool,
              "scope": {"source_database_version": None},
              "provenance": {"audit": "source evidence supplied by caller"}}
    ordinary = align_sequences(output=tmp_path / "ordinary", **kwargs)
    streamed = align_sequences(output=tmp_path / "stream", materialize_samples=False,
                               export_concatenation=export, **kwargs)
    reloaded = load_alignment(tmp_path / "stream", rows)
    assert reloaded.aligned_alleles == streamed.aligned_alleles
    assert reloaded.samples == {}
    assert load_alignment(tmp_path / "stream", rows, materialize_samples=True).samples == ordinary.samples
    rows[0]["cgmlst_profile"]["x"] = "changed"
    with pytest.raises(ValueError, match="sample profiles differ"):
        load_alignment(tmp_path / "stream", rows)
    rows[0]["cgmlst_profile"]["x"] = "1"
    assert streamed.samples == {}
    assert streamed.sample_ids == ("a", "b")
    assert streamed.aligned_alleles == ordinary.aligned_alleles
    assert streamed.provenance["scope"] == kwargs["scope"]
    assert streamed.provenance["catalogue"] == kwargs["provenance"]
    fasta = tmp_path / "stream" / "concatenated.fasta"
    if export:
        assert read_fasta(fasta.read_text()) == ordinary.samples
    else:
        assert not fasta.exists()
        assert streamed.provenance["concatenation_sha256"] is None
    with pytest.raises(ValueError, match="materialized"):
        next(streamed.pairs(min_loci=1))


@pytest.mark.parametrize("artifact", ["partitions.json", "missing-calls.json", "loci/0001/aligned.fasta",
                                      "loci/0001/allele-map.json", "concatenated.fasta"])
def test_reload_rejects_corrupted_persisted_evidence(tmp_path, artifact):
    tool = tmp_path / "mafft"
    tool.write_text(f"#!{sys.executable}\nprint('fixture 1')\n")
    tool.chmod(0o755)
    rows = [{"sample_id": "a", "cgmlst_loci": ["x"], "cgmlst_profile": {"x": "1"}}]
    output = tmp_path / "alignment"
    align_sequences({("x", "1"): "ACGT"}, rows, output, tool, scope={}, provenance={"fixture": True},
                    materialize_samples=False)
    with (output / artifact).open("a") as stream:
        stream.write("corrupt")
    with pytest.raises(ValueError, match="hash mismatch"):
        load_alignment(output, rows)
