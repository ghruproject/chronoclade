import csv
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from beyondmlst.context import (
    DEFAULT_CONTEXT_METADATA,
    ContextCandidate,
    attach_downloaded_assemblies,
    context_metadata_provenance,
    filter_candidates,
    load_same_st_candidates,
    normalise_collection_date,
    parse_ska_distances,
    select_context,
    stratified_candidate_pool,
    write_candidate_table,
    write_combined_metadata,
)
from beyondmlst.metadata import Sample, read_metadata


def candidate(sample_id: str, country: str, year: str) -> ContextCandidate:
    return ContextCandidate(
        sample_id=sample_id,
        species="Escherichia coli",
        lineage="ST131",
        mlst_scheme="ecoli_achtman_4",
        mlst_st="131",
        collection_date=year,
        country=country,
        hq_filter="PASS",
    )


def test_normalises_only_supported_collection_dates() -> None:
    assert normalise_collection_date("2024") == "2024"
    assert normalise_collection_date("2024-03") == "2024-03"
    assert normalise_collection_date("2024-03-14T12:00:00") == "2024-03-14"
    assert normalise_collection_date("2024-13") == ""
    assert normalise_collection_date("not collected") == ""


def write_test_snapshot(path: Path) -> Path:
    rows = [
        {
            "sample_id": "SAMN1",
            "species": "Escherichia coli",
            "mlst_scheme": "ecoli_achtman_4",
            "mlst_st": "131",
            "mlst_status": "PERFECT",
            "collection_date": "2023-04",
            "country": "United Kingdom",
            "host": "Homo sapiens",
            "isolation_source": "blood",
            "hq_filter": "PASS",
            "completeness": 99.0,
            "contamination": 0.2,
            "genome_size": 5_000_000,
            "contig_n50": 100_000,
            "aws_url": "https://example/SAMN1.fa.gz",
        },
        {
            "sample_id": "SAMN2",
            "species": "Escherichia coli",
            "mlst_scheme": "ecoli_achtman_4",
            "mlst_st": "131",
            "mlst_status": "PERFECT",
            "collection_date": "",
            "country": "Philippines",
            "host": "Homo sapiens",
            "isolation_source": "urine",
            "hq_filter": "PASS",
            "completeness": 99.0,
            "contamination": 0.2,
            "genome_size": 5_000_000,
            "contig_n50": 100_000,
            "aws_url": "https://example/SAMN2.fa.gz",
        },
    ]
    pq.write_table(pa.Table.from_pylist(rows), path)
    return path


def test_compact_snapshot_returns_all_and_dated_same_st_candidates(tmp_path: Path) -> None:
    snapshot = write_test_snapshot(tmp_path / "context.parquet")

    accessions, candidates = load_same_st_candidates(
        snapshot,
        species="Escherichia coli",
        lineage="ST131",
        scheme="ecoli_achtman_4",
        st="131",
    )

    assert accessions == ["SAMN1", "SAMN2"]
    assert [candidate.sample_id for candidate in candidates] == ["SAMN1"]
    assert candidates[0].collection_date == "2023-04"
    assert candidates[0].completeness == "99.0"


def test_compact_snapshot_matches_underscore_separated_species(tmp_path: Path) -> None:
    snapshot = write_test_snapshot(tmp_path / "context.parquet")

    accessions, _ = load_same_st_candidates(
        snapshot,
        species="Escherichia_coli",
        lineage="ST131",
        scheme="ecoli_achtman_4",
        st="131",
    )

    assert accessions == ["SAMN1", "SAMN2"]


def test_snapshot_provenance_includes_checksum_and_manifest(tmp_path: Path) -> None:
    snapshot = write_test_snapshot(tmp_path / "context.parquet")
    snapshot.with_suffix(".json").write_text('{"atb_release": "test-release"}\n', encoding="utf-8")

    provenance = context_metadata_provenance(snapshot)

    assert provenance["release"] == "test-release"
    assert provenance["rows"] == 2
    assert len(str(provenance["sha256"])) == 64


def test_bundled_snapshot_contains_expected_st131_records() -> None:
    accessions, candidates = load_same_st_candidates(
        DEFAULT_CONTEXT_METADATA,
        species="Escherichia coli",
        lineage="ST131",
        scheme="ecoli_achtman_4",
        st="131",
    )

    assert len(accessions) == 13_579
    assert len(candidates) == 8_953
    assert DEFAULT_CONTEXT_METADATA.stat().st_size < 10 * 1024 * 1024
    assert "mlst_status" not in pq.ParquetFile(DEFAULT_CONTEXT_METADATA).schema.names


def test_bundled_snapshot_contains_expected_klebsiella_st258_records() -> None:
    accessions, candidates = load_same_st_candidates(
        DEFAULT_CONTEXT_METADATA,
        species="Klebsiella pneumoniae",
        lineage="ST258",
        scheme="klebsiella",
        st="258",
    )

    assert len(accessions) == 3_036
    assert len(candidates) == 2_644


def test_stratified_pool_is_reproducible_and_spans_strata() -> None:
    candidates = [
        candidate("UK1", "United Kingdom", "2021"),
        candidate("UK2", "United Kingdom", "2021"),
        candidate("PH1", "Philippines", "2022"),
        candidate("PH2", "Philippines", "2022"),
        candidate("IN1", "India", "2023"),
    ]

    first = stratified_candidate_pool(candidates, limit=3, seed=17)
    second = stratified_candidate_pool(candidates, limit=3, seed=17)

    assert [item.sample_id for item in first] == [item.sample_id for item in second]
    assert {item.country for item in first} == {"United Kingdom", "Philippines", "India"}


def test_filter_candidates_applies_metadata_constraints() -> None:
    candidates = [
        candidate("A", "United Kingdom:England", "2021"),
        candidate("B", "Philippines", "2024"),
    ]
    candidates[0].host = "Homo sapiens"
    candidates[0].isolation_source = "blood"

    result = filter_candidates(
        candidates,
        countries=["United Kingdom"],
        year_from=2020,
        year_to=2022,
        host="sapiens",
        isolation_source="blood",
    )

    assert [item.sample_id for item in result] == ["A"]


def test_parse_and_select_context_prioritises_focal_neighbours(tmp_path: Path) -> None:
    distance_file = tmp_path / "distances.tsv"
    distance_file.write_text(
        "Sample1\tSample2\tDistance\tMismatches (proportion)\tMatch count\tMismatch count\n"
        "F1\tC1\t2.00\t0.001\t100\t1\n"
        "F1\tC2\t30.00\t0.030\t100\t3\n"
        "F1\tC3\t40.00\t0.040\t100\t4\n"
        "F2\tC1\t35.00\t0.035\t100\t3\n"
        "F2\tC2\t3.00\t0.003\t100\t1\n"
        "F2\tC3\t20.00\t0.020\t100\t2\n",
        encoding="utf-8",
    )
    distances = parse_ska_distances(distance_file)
    candidates = [
        candidate("C1", "United Kingdom", "2021"),
        candidate("C2", "Philippines", "2022"),
        candidate("C3", "India", "2023"),
    ]

    selected, audit = select_context(
        candidates,
        ["F1", "F2"],
        distances,
        max_context=3,
        nearest_per_focal=1,
        seed=2,
    )

    assert {item.sample_id for item in selected} == {"C1", "C2", "C3"}
    assert next(item for item in selected if item.sample_id == "C1").selection_reason == (
        "nearest_to=F1"
    )
    assert next(item for item in selected if item.sample_id == "C2").selection_reason == (
        "nearest_to=F2"
    )
    assert audit["focal_neighbour_coverage"] == 1.0


def test_zero_distance_is_ranked_before_nonzero_distance() -> None:
    candidates = [
        candidate("IDENTICAL", "United Kingdom", "2021"),
        candidate("DIFFERENT", "Philippines", "2022"),
    ]
    distances = {
        ("IDENTICAL", "F1"): (0.0, 0.0),
        ("DIFFERENT", "F1"): (1.0, 0.001),
    }

    selected, _ = select_context(
        candidates,
        ["F1"],
        distances,
        max_context=2,
        nearest_per_focal=1,
        seed=2,
    )

    assert [item.sample_id for item in selected] == ["IDENTICAL", "DIFFERENT"]


def test_candidate_manifest_is_tab_separated(tmp_path: Path) -> None:
    manifest = write_candidate_table(
        tmp_path / "context_manifest.tsv",
        [candidate("SAMN1", "Philippines", "2023")],
    )

    header = manifest.read_text(encoding="utf-8").splitlines()[0]
    assert "sample_id\tspecies\tlineage" in header


def test_download_attachment_and_combined_metadata_round_trip(tmp_path: Path) -> None:
    assembly = tmp_path / "focal.fasta"
    assembly.write_text(">focal\nAAAA\n", encoding="utf-8")
    download_dir = tmp_path / "downloads"
    download_dir.mkdir()
    context_assembly = download_dir / "SAMN1.fa.gz"
    context_assembly.write_bytes(b"placeholder")
    context = candidate("SAMN1", "Philippines", "2023-02")
    attached, missing = attach_downloaded_assemblies([context], download_dir)
    assert not missing
    assert attached[0].assembly == str(context_assembly.resolve())

    focal = [
        Sample(
            "F1",
            assembly,
            "2024-01",
            "KIMS",
            "Escherichia coli",
            "ST131",
            "local",
        )
    ]
    combined = write_combined_metadata(tmp_path / "combined.csv", focal, attached)
    parsed = read_metadata(combined)

    assert [sample.origin for sample in parsed] == ["local", "context"]
    with combined.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[1]["location"] == "Philippines"
