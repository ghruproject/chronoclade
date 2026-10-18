"""Synthetic tests only; no labelled study outcomes are accessed."""

import csv
import json
import subprocess
import sys

import pytest

from chronoclade.datasets import AlleleMatrix, LocusCatalogue, PreparedDataset, write_dataset
from research.klebsiella_grouping.evaluate import evaluate, read_labels, time_bin


def fixture_dataset():
    ids = ("a", "b", "c", "d")
    loci = tuple(f"L{i}" for i in range(629))
    catalogue = LocusCatalogue("synthetic629", "v1", loci, "synthetic", database_sha256="a" * 64)
    calls = {ident: dict.fromkeys(loci, "1") for ident in ids[:3]}
    calls["b"]["L0"] = "2"
    calls["c"]["L0"] = calls["c"]["L1"] = "2"
    matrix = AlleleMatrix.from_profiles(catalogue, ids, calls)
    samples = tuple(
        {"sample_id": ident, "label": ident, "role": "input", "collection_date": "2020-01-01"}
        for ident in ids
    )
    lineages = tuple(
        {
            "sample_id": ident,
            "kind": "cglin",
            "scheme_id": "syntheticLIN",
            "scheme_version": "v1",
            "database_version": "snapshot1",
            "assignment": "0_0_0_0_0_0_0_0_0_0",
            "assignment_method": "synthetic",
            "resolution_status": "resolved",
            "evidence": {},
        }
        for ident in ids[:3]
    )
    return PreparedDataset(samples=samples, profiles=(matrix,), lineages=lineages)


def fixture_labels():
    return {
        ident: {
            "sample_id": ident,
            "study_id": "study1",
            "event_id": "event1" if ident in "ab" else "",
            "epi_status": "outbreak" if ident in "ab" else "control",
            "control_for_event": "" if ident in "ab" else "event1",
            "source": "synthetic independent records",
        }
        for ident in "abcd"
    }


def test_metrics_missingness_and_unknowns():
    result = evaluate(fixture_dataset(), fixture_labels())
    pairs = {(p["sample_a"], p["sample_b"]): p for p in result["pairs"]}
    assert pairs["a", "b"]["epi_label"] == "positive"
    assert pairs["a", "c"]["epi_label"] == "negative"
    assert pairs["c", "d"]["epi_label"] == "unknown"
    assert pairs["a", "b"]["allele_mismatches"] == 1
    assert pairs["a", "d"]["allele_mismatches"] is None
    assert pairs["a", "d"]["alleles_le_0"] is None
    assert not pairs["a", "d"]["eligible"]
    assert result["samples"][3]["called_loci"] == 0
    assert not result["samples"][3]["LIN5_assigned"]
    query = next(
        r for r in result["queries"] if r["sample_id"] == "a" and r["candidate"] == "alleles_le_1"
    )
    assert query["recovery"] == 1 and query["burden"] == 1
    lin = next(r for r in result["queries"] if r["sample_id"] == "a" and r["candidate"] == "LIN5")
    assert lin["burden"] == 2 and lin["negative_retrieved"] == 1
    nearest = next(
        r
        for r in result["queries"]
        if r["sample_id"] == "b" and r["candidate"] == "nearest_all_ties"
    )
    assert nearest["burden"] == nearest["nearest_ties_retained"] == 2
    assert json.loads(nearest["retrieved_ids"]) == ["a", "c"]
    assert not any(r["candidate"] == "nearest_all_ties" for r in result["candidate_counts"])
    labels = fixture_labels()
    labels["c"]["control_for_event"] = ""
    assert evaluate(fixture_dataset(), labels)["pairs"][1]["epi_label"] == "unknown"


def test_coverage_boundary_and_snapshot_abstention():
    dataset = fixture_dataset()
    codes = dataset.profiles[0].codes.copy()
    codes[1, :30] = 0
    matrix = AlleleMatrix(
        dataset.profiles[0].catalogue, dataset.sample_ids, dataset.profiles[0].categories, codes
    )
    altered = PreparedDataset(
        samples=dataset.samples, profiles=(matrix,), lineages=dataset.lineages
    )
    assert evaluate(altered, fixture_labels())["pairs"][0]["eligible"]
    codes[1, 30] = 0
    assert not evaluate(altered, fixture_labels())["pairs"][0]["eligible"]
    assert evaluate(altered, fixture_labels(), min_shared_loci=598)["pairs"][0]["eligible"]
    sensitivity = evaluate(altered, fixture_labels(), min_shared_loci=567)
    assert sensitivity["summary"]["coverage_policy"]["min_shared_loci"] == 567
    assert not evaluate(altered, fixture_labels(), min_shared_loci=629)["pairs"][0]["eligible"]
    dataset.lineages[0]["database_version"] = None
    pair = evaluate(dataset, fixture_labels())["pairs"][0]
    assert pair["LIN5"] is None
    assert time_bin("2020", "2020-01-01") == "unknown"
    for row in dataset.lineages:
        row["database_version"] = None
        row["evidence"]["cglin_database_sha256"] = "b" * 64
    assert evaluate(dataset, fixture_labels())["pairs"][0]["LIN5"] is True
    assert time_bin("2020-01-01", "2020-02-01") == "31-90"
    assert time_bin("2020-01-01", "2020-05-01") == ">90"


def test_cli_prepared_bundle(tmp_path):
    dataset = fixture_dataset()
    manifest = write_dataset(dataset, tmp_path / "prepared")
    labels_path = tmp_path / "labels.csv"
    rows = list(fixture_labels().values())
    with labels_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    assert read_labels(labels_path, dataset.sample_ids) == fixture_labels()
    output = tmp_path / "results"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "research.klebsiella_grouping.evaluate",
            "--dataset",
            str(manifest),
            "--labels",
            str(labels_path),
            "--output",
            str(output),
        ],
        check=True,
    )
    summary = json.loads((output / "summary.json").read_text())
    assert summary["eligible_pair_counts"] == {"positive": 1, "negative": 2}
    assert summary["individuals"] == 4 and summary["documented_events"] == 1
    assert summary["coverage_policy"]["min_shared_loci"] == 599
    with pytest.raises(ValueError, match="exactly"):
        read_labels(labels_path, ("a", "b", "c", "d", "e"))
