#!/usr/bin/env python3
"""Build the small AllTheBacteria metadata snapshot bundled with beyondMLST."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sqlite3
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from beyondmlst.context import normalise_collection_date


SPECIES_SCHEMES = {
    "Escherichia coli": "ecoli_achtman_4",
    "Klebsiella pneumoniae": "klebsiella",
}
ATB_RELEASE = "2025-05"
ATB_SQLITE_URL = "https://osf.io/download/my56u/"
ATB_MLST_URL = "https://osf.io/download/69c66d33fa3d973d94254f46/"
ATB_DOI = "https://doi.org/10.1101/2024.03.08.584059"


def normalise_species(value: str) -> str:
    return re.sub(r"_[A-Z]\b", "", value).strip()


def latest_ena_table(connection: sqlite3.Connection) -> str:
    names = sorted(
        (
            str(row[0])
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'ena_%'"
            )
            if re.fullmatch(r"ena_\d{8}", str(row[0]))
        ),
        reverse=True,
    )
    if not names:
        raise RuntimeError("No dated ENA metadata table was found in the ATB SQLite snapshot")
    return names[0]


def build_snapshot(sqlite_path: Path, mlst_path: Path) -> tuple[pd.DataFrame, str]:
    connection = sqlite3.connect(sqlite_path)
    assembly = pd.read_sql_query(
        """
        SELECT a.sample_accession AS sample_id,
               a.sylph_species AS species,
               a.hq_filter,
               a.aws_url,
               c.Completeness_Specific AS completeness,
               c.Contamination AS contamination,
               c.Genome_Size AS genome_size,
               c.Contig_N50 AS contig_n50
          FROM assembly a
     LEFT JOIN checkm2 c ON a.sample_accession = c.sample_accession
         WHERE a.hq_filter = 'PASS'
           AND a.asm_fasta_on_osf = 1
           AND (a.sylph_species LIKE 'Escherichia coli%'
                OR a.sylph_species LIKE 'Klebsiella pneumoniae%')
        """,
        connection,
    )
    assembly["normalised_species"] = assembly["species"].map(normalise_species)
    assembly = assembly[assembly["normalised_species"].isin(SPECIES_SCHEMES)].copy()
    assembly["mlst_scheme"] = assembly["normalised_species"].map(SPECIES_SCHEMES)

    ena_table = latest_ena_table(connection)
    connection.execute("CREATE TEMP TABLE wanted (sample_accession TEXT PRIMARY KEY)")
    connection.executemany(
        "INSERT INTO wanted VALUES (?)",
        ((sample_id,) for sample_id in assembly["sample_id"]),
    )
    metadata = pd.read_sql_query(
        f"""
        SELECT e.sample_accession AS sample_id, e.run_accession,
               e.collection_date, e.country, e.host, e.isolation_source
          FROM {ena_table} e
          JOIN wanted w ON e.sample_accession = w.sample_accession
        """,
        connection,
    )
    connection.close()
    metadata["has_date"] = metadata["collection_date"].fillna("").ne("")
    metadata = (
        metadata.sort_values(
            ["sample_id", "has_date", "run_accession"],
            ascending=[True, False, True],
        )
        .drop_duplicates("sample_id", keep="first")
        .drop(columns=["has_date", "run_accession"])
    )

    mlst = pd.read_parquet(
        mlst_path,
        columns=["sample", "mlst_scheme", "mlst_st", "mlst_status"],
    ).rename(columns={"sample": "sample_id"})
    mlst = mlst[mlst["mlst_scheme"].isin(SPECIES_SCHEMES.values())].copy()
    mlst = mlst[mlst["mlst_status"] == "PERFECT"]
    mlst["mlst_st"] = mlst["mlst_st"].astype(str)
    mlst = mlst[mlst["mlst_st"] != "-"]
    mlst = mlst.drop_duplicates(["sample_id", "mlst_scheme"], keep="first")
    mlst = mlst.drop(columns="mlst_status")

    result = assembly.merge(mlst, on=["sample_id", "mlst_scheme"], how="inner")
    result = result.merge(metadata, on="sample_id", how="left")
    result["species"] = result.pop("normalised_species")
    columns = [
        "sample_id",
        "species",
        "mlst_scheme",
        "mlst_st",
        "collection_date",
        "country",
        "host",
        "isolation_source",
        "hq_filter",
        "completeness",
        "contamination",
        "genome_size",
        "contig_n50",
        "aws_url",
    ]
    result = result[columns].sort_values(["species", "mlst_scheme", "mlst_st", "sample_id"])
    return result, ena_table


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sqlite", type=Path, required=True, help="ATB 2025-05 SQLite file")
    parser.add_argument("--mlst", type=Path, required=True, help="Official ATB mlst.parquet")
    parser.add_argument("--output", type=Path, required=True, help="Output Parquet path")
    parser.add_argument("--manifest", type=Path, required=True, help="Output provenance JSON")
    args = parser.parse_args()

    result, ena_table = build_snapshot(args.sqlite, args.mlst)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".part")
    table = pa.Table.from_pandas(result, preserve_index=False)
    table = table.replace_schema_metadata(
        {
            b"atb_release": ATB_RELEASE.encode(),
            b"atb_sqlite_url": ATB_SQLITE_URL.encode(),
            b"atb_mlst_url": ATB_MLST_URL.encode(),
            b"atb_citation": ATB_DOI.encode(),
            b"licence": b"MIT",
        }
    )
    pq.write_table(table, temporary, compression="zstd")
    temporary.replace(args.output)

    manifest = {
        "name": "beyondMLST AllTheBacteria context metadata",
        "atb_release": ATB_RELEASE,
        "ena_table": ena_table,
        "species_schemes": SPECIES_SCHEMES,
        "rows": len(result),
        "usable_collection_date_rows": int(
            result["collection_date"]
            .fillna("")
            .map(lambda value: bool(normalise_collection_date(str(value))))
            .sum()
        ),
        "columns": list(result.columns),
        "file": args.output.name,
        "size_bytes": args.output.stat().st_size,
        "sha256": sha256(args.output),
        "generated_at": datetime.now(UTC).isoformat(),
        "sources": {
            "sqlite": ATB_SQLITE_URL,
            "mlst": ATB_MLST_URL,
            "documentation": "https://allthebacteria.org/docs/sample_metadata/",
            "citation": ATB_DOI,
        },
        "licence": "MIT",
        "selection": (
            "HQ downloadable Escherichia coli/ecoli_achtman_4 and "
            "Klebsiella pneumoniae/klebsiella records with assigned STs"
        ),
    }
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(
        f"Wrote {len(result):,} rows to {args.output} "
        f"({args.output.stat().st_size / 1024 / 1024:.1f} MiB)"
    )


if __name__ == "__main__":
    main()
