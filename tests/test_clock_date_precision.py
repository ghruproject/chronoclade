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


@pytest.mark.parametrize("value", ["1980-XX-XX", "1980-05-XX", "1980-XX-17", "1980.5", "[1980:1981.5]"])
def test_native_eligibility_accepts_valid_clock_dates_without_integer_split(value):
    from chronoclade.native_pipeline import _dated_members
    from chronoclade.metadata_dates import clock_date_interval

    sample = _sample("input", value)
    assert clock_date_interval(value)["status"] == "valid"
    assert _dated_members([sample]) == [sample]
    assert sample.collection_date == value


@pytest.mark.parametrize("value", ["[1982:1980]", "nan", "inf", "1980-13-XX", "1980-02-30", "0000.5", "before 2009", "2006-2009"])
def test_native_eligibility_rejects_invalid_or_undeclared_calibrations_without_crash(value):
    from chronoclade.native_pipeline import _dated_members

    sample = _sample("input", value)
    assert _dated_members([sample]) == []
    assert sample.collection_date == value


def test_native_future_interval_eligibility_uses_lower_bound_without_clipping():
    from datetime import date
    from chronoclade.native_pipeline import _dated_members
    from chronoclade.metadata_dates import clock_date_interval

    current = date.today().year
    spanning = _sample("spanning", f"[{current - 1}:{current + 1}]")
    future = _sample("future", f"[{current + 1}:{current + 2}]")
    assert _dated_members([spanning, future]) == [spanning]
    assert clock_date_interval(spanning.collection_date)["year_max"] == current + 1


def test_initial_published_cohort_retains_97_inputs_with_92_calibration_dates():
    from chronoclade.native_pipeline import _dated_members

    known = [_sample(f"known{i}", str(1932 + i % 84)) for i in range(92)]
    uncertain = [_sample(f"uncertain{i}", value) for i, value in enumerate(
        ["2006-2009", "before 2009", "1996-2001", "before 2009", "2008-2013"])]
    cohort = known + uncertain
    original = [sample.collection_date for sample in cohort]
    assert _dated_members(cohort) == known
    assert len(cohort) == 97 and [sample.collection_date for sample in cohort] == original


def test_legacy_clock_retains_uncertain_tips_without_inventing_date_bounds(tmp_path):
    import math
    from chronoclade.native_pipeline import LineageFiles, _clock_dates_path
    from treetime.utils import parse_dates

    files = LineageFiles.in_directory(tmp_path)
    files.metadata.write_text("sample_id,collection_date\nknown,1980\nrange,[1981:1982]\ncensored,before 2009\ninformal,2006-2009\n")
    original = files.metadata.read_bytes()
    target = _clock_dates_path(files)
    rows = list(csv.DictReader(target.open()))
    assert len(rows) == 4
    assert {r["sample_id"]: r["collection_date"] for r in rows} == {
        "known": "1980-XX-XX", "range": "[1981:1982]", "censored": "", "informal": "",
    }
    parsed = parse_dates(str(target), name_col="sample_id", date_col="collection_date")
    assert parsed["range"] == [1981.0, 1982.0]
    # This installed parser represents blank constraints as NaN; ClockTree
    # explicitly treats those as undated, retaining the structural tips.
    assert math.isnan(parsed["censored"]) and math.isnan(parsed["informal"])
    assert files.metadata.read_bytes() == original
