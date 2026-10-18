"""The native clock must receive intervals for imprecise collection dates."""

import csv
from pathlib import Path

import pytest


def _sample(ident, value):
    from chronoclade.metadata import Sample

    return Sample(ident, Path("assembly.fa"), value, "France", "Klebsiella pneumoniae", "Unassigned", "local")


@pytest.mark.parametrize(
    ("raw", "serialized"),
    [("1980", "1980-XX-XX"), ("1980-05", "1980-05-XX"),
     ("1980-05-17", "1980-05-17"), ("1980-XX-XX", "1980-XX-XX"),
     ("[1980.0:1981.0]", "[1980.0:1981.0]"), ("1980.5", "1980.5")],
)
def test_clock_serialization_preserves_declared_precision(raw, serialized):
    from chronoclade.metadata_dates import treetime_date

    assert treetime_date(raw) == serialized


def test_installed_treetime_reads_year_month_intervals_and_exact_day(tmp_path):
    from treetime.utils import parse_dates
    from chronoclade.temporal import _write_dates

    samples = [_sample("year", "1980"), _sample("month", "1980-05"), _sample("day", "1980-05-17")]
    path = tmp_path / "dates.csv"
    _write_dates(path, samples)
    parsed = parse_dates(str(path), name_col="sample_id", date_col="collection_date")
    assert len(parsed["year"]) == 2
    assert 1980 <= parsed["year"][0] < 1980.01
    assert 1980.99 < parsed["year"][1] < 1981
    assert len(parsed["month"]) == 2
    assert 0.07 < parsed["month"][1] - parsed["month"][0] < 0.09
    assert isinstance(parsed["day"], float)
    assert parsed["month"][0] < parsed["day"] < parsed["month"][1]
    assert [s.collection_date for s in samples] == ["1980", "1980-05", "1980-05-17"]


def test_permuted_value_precision_moves_with_assigned_date(tmp_path):
    from chronoclade.temporal import _write_dates
    from chronoclade.tree_stage.temporal import date_clusters
    from Bio import Phylo
    from io import StringIO

    samples = [_sample("a", "1980"), _sample("b", "1981-05"), _sample("c", "1980")]
    tree = Phylo.read(StringIO("((a:1,c:1):1,b:2);"), "newick")
    clusters_before = date_clusters(tree, samples)
    path = tmp_path / "permuted.csv"
    _write_dates(path, samples, ["1981-05", "1980", "1981-05"])
    rows = list(csv.DictReader(path.open()))
    assert {r["sample_id"]: r["collection_date"] for r in rows} == {
        "a": "1981-05-XX", "b": "1980-XX-XX", "c": "1981-05-XX",
    }
    assert date_clusters(tree, samples) == clusters_before
    assert {c["collection_date"] for c in clusters_before} == {"1980", "1981-05"}


@pytest.mark.parametrize("use_explicit_clock_file", [False, True])
def test_observed_and_dated_native_routes_use_intervals_without_metadata_mutation(tmp_path, monkeypatch, use_explicit_clock_file):
    from chronoclade import native_pipeline as native

    files = native.LineageFiles.in_directory(tmp_path)
    files.metadata.write_text("sample_id,collection_date,location\na,1980,France\nb,1981-05,France\n")
    metadata = files.metadata.read_bytes()
    if use_explicit_clock_file:
        (tmp_path / "clock_dates.csv").write_text("sample_id,collection_date\na,1980\nb,1981-05\n")
    captured = []

    def run(command, **kwargs):
        source = Path(command[command.index("--dates") + 1])
        assert source != files.metadata
        captured.append(list(csv.DictReader(source.open())))

    monkeypatch.setattr(native, "_run_command", run)
    native._run_observed_clock(files, files.tree, 1000, False)
    native._run_dated_tree(files, 1000, False)
    assert len(captured) == 2
    assert all({r["sample_id"]: r["collection_date"] for r in rows} == {
        "a": "1980-XX-XX", "b": "1981-05-XX",
    } for rows in captured)
    assert files.metadata.read_bytes() == metadata
