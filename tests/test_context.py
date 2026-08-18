import csv
import sqlite3
import subprocess
from pathlib import Path

import beyondmlst.context as context_module
from beyondmlst.context import (
    ContextCandidate,
    attach_downloaded_assemblies,
    cache_official_atb_mlst,
    discover_same_st,
    discover_same_st_resilient,
    filter_candidates,
    lookup_atb_metadata,
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


def test_caches_official_mlst_parquet_atomically(tmp_path: Path, monkeypatch) -> None:
    source = tmp_path / "source.parquet"
    source.write_bytes(b"PAR1mock-parquet-content")
    monkeypatch.setattr(context_module, "OFFICIAL_ATB_MLST_URL", source.as_uri())

    result = cache_official_atb_mlst(tmp_path / "cache")

    assert result.read_bytes() == source.read_bytes()
    assert not result.with_suffix(".parquet.part").exists()
    provenance = result.with_suffix(".parquet.provenance.json")
    assert provenance.is_file()
    assert "official_allthebacteria_osf" in provenance.read_text(encoding="utf-8")


def test_resilient_discovery_repairs_missing_atbfetcher_cache(tmp_path: Path, monkeypatch) -> None:
    calls = 0

    def fake_discovery(**_kwargs) -> list[str]:
        nonlocal calls
        calls += 1
        if calls == 1:
            raise context_module.ContextError("upstream returned 401")
        return ["SAMN2", "SAMN1"]

    cached = tmp_path / "mlst.parquet"

    def fake_cache(_cache_dir: Path) -> Path:
        cached.write_bytes(b"PAR1fallback")
        return cached

    monkeypatch.setattr(context_module, "discover_same_st", fake_discovery)
    monkeypatch.setattr(context_module, "cache_official_atb_mlst", fake_cache)

    accessions, audit = discover_same_st_resilient(
        species="Escherichia coli",
        scheme="ecoli_achtman_4",
        st="131",
        cache_dir=tmp_path,
    )

    assert accessions == ["SAMN2", "SAMN1"]
    assert calls == 2
    assert audit["method"] == "atbfetcher_with_official_atb_mlst_fallback"


def test_discovery_ignores_atbfetcher_logging_before_tsv_header(
    tmp_path: Path, monkeypatch
) -> None:
    output = (
        "[20:01:05] INFO Loading cached MLST data\n"
        "sample\tmlst_scheme\tmlst_st\n"
        "SAMN2\tecoli_achtman_4\t131\n"
        "SAMN1\tecoli_achtman_4\t131\n"
    )
    completed = subprocess.CompletedProcess(["atbfetcher"], 0, output, "")
    monkeypatch.setattr(context_module, "_run_capture", lambda *_args, **_kwargs: completed)

    result = discover_same_st(
        species="Escherichia coli",
        scheme="ecoli_achtman_4",
        st="131",
        cache_dir=tmp_path,
    )

    assert result == ["SAMN1", "SAMN2"]


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


def test_lookup_atb_metadata_uses_dated_hq_assemblies(tmp_path: Path) -> None:
    database = tmp_path / "atb.metadata.test.sqlite"
    connection = sqlite3.connect(database)
    connection.executescript(
        """
        CREATE TABLE assembly (
          sample_accession TEXT, sylph_species TEXT, hq_filter TEXT,
          aws_url TEXT, asm_fasta_on_osf INTEGER
        );
        CREATE TABLE run (sample_accession TEXT, run_accession TEXT);
        CREATE TABLE ena_202505_used (
          run_accession TEXT, country TEXT, collection_date TEXT,
          host TEXT, isolation_source TEXT
        );
        CREATE TABLE checkm2 (
          sample_accession TEXT, Completeness_Specific REAL,
          Contamination REAL, Genome_Size INTEGER, Contig_N50 INTEGER
        );
        INSERT INTO assembly VALUES ('SAMN1','Escherichia coli','PASS','https://a',1);
        INSERT INTO assembly VALUES ('SAMN2','Escherichia coli','PASS','https://b',1);
        INSERT INTO run VALUES ('SAMN1','ERR1');
        INSERT INTO run VALUES ('SAMN2','ERR2');
        INSERT INTO ena_202505_used VALUES ('ERR1','United Kingdom','2023-04','Homo sapiens','blood');
        INSERT INTO ena_202505_used VALUES ('ERR2','Philippines','','Homo sapiens','blood');
        INSERT INTO checkm2 VALUES ('SAMN1',99.0,0.2,5000000,100000);
        INSERT INTO checkm2 VALUES ('SAMN2',99.0,0.2,5000000,100000);
        """
    )
    connection.commit()
    connection.close()

    result = lookup_atb_metadata(
        database,
        ["SAMN1", "SAMN2"],
        species="Escherichia coli",
        lineage="ST131",
        scheme="ecoli_achtman_4",
        st="131",
    )

    assert [item.sample_id for item in result] == ["SAMN1"]
    assert result[0].collection_date == "2023-04"
    assert result[0].completeness == "99.0"


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
