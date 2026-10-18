"""Portable profile, tree and metadata exports for downstream viewers.

Exports use deterministic ChronoClade IDs in viewer files and retain a separate
crosswalk to source sample IDs. GrapeTree allele values are recoded as numeric
categories per locus; missing calls are always ``0``.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from Bio import Phylo

from chronoclade.datasets.model import _allele
from chronoclade.sample_labels import sample_label


METADATA_FIELDS = (
    "species",
    "mlst_st",
    "collection_date",
    "date_precision",
    "country",
    "region",
    "city",
    "host",
    "isolation_source",
    "cglin_raw",
    "hiercc_codes",
    "cgmlst_scheme",
    "cgmlst_scheme_version",
)


def _source_ids(records: Sequence[Mapping]) -> list[str]:
    ids = [row.get("sample_id") for row in records]
    if any(not isinstance(value, str) or not value for value in ids):
        raise ValueError("Every export record requires a nonempty sample_id")
    if len(set(ids)) != len(ids):
        raise ValueError("Export records contain duplicate sample_id values")
    return ids


def _export_ids(ids: Iterable[str]) -> dict[str, str]:
    result = {
        ident: "CC_" + hashlib.sha256(ident.encode("utf-8")).hexdigest()[:16]
        for ident in ids
    }
    if len(set(result.values())) != len(result):
        raise ValueError("Source IDs collided after safe export-ID generation")
    return result


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
                          allow_nan=False)
    return str(value)


def _metadata_records(records: Sequence[Mapping], ids: Mapping[str, str], fields):
    rows = []
    for row in records:
        source_id = row["sample_id"]
        label = sample_label(row) or row.get("strain_name") or source_id
        item = {
            "ID": ids[source_id],
            "chronoclade_sample_id": source_id,
            "sample_label": label,
            "source_genome_id": row.get("source_genome_id"),
            "accession": row.get("accession") or row.get("run_accession") or row.get("biosample"),
        }
        item.update({field: row.get(field) for field in fields})
        rows.append(item)
    return rows


def _write_delimited(path: Path, rows: Sequence[Mapping], fields, delimiter: str):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter=delimiter,
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _cell(row.get(key)) for key in fields})


def _new_output_dir(output_dir: str | Path) -> Path:
    output = Path(output_dir).expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Export output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    return output


def _normalized_call(value):
    try:
        return _allele(value)
    except ValueError as error:
        raise ValueError(f"Invalid categorical allele call: {error}") from error


def write_grapetree_exports(
    records: Iterable[Mapping],
    output_dir: str | Path,
    *,
    loci: Sequence[str] | None = None,
    metadata_fields: Sequence[str] = METADATA_FIELDS,
) -> dict[str, Path]:
    """Write a GrapeTree profile TSV, metadata TSV and auditable allele crosswalk.

    ``records`` are canonical rows with ``sample_id`` and ``cgmlst_profile``.
    If ``loci`` is omitted, the sorted union of observed profile keys is used and
    the output README states that it is not a claim of complete scheme coverage.
    """
    rows = [dict(row) for row in records]
    _source_ids(rows)
    if not rows:
        raise ValueError("At least one record is required for a GrapeTree export")
    unknown_fields = set(metadata_fields) - set(METADATA_FIELDS)
    if unknown_fields:
        raise ValueError(f"Unsupported export metadata fields: {sorted(unknown_fields)}")

    profiled = [row for row in rows if isinstance(row.get("cgmlst_profile"), Mapping)
                and row["cgmlst_profile"]]
    if not profiled:
        raise ValueError("No records contain a usable cgMLST profile")
    scopes = {
        (
            row.get("cgmlst_scheme"),
            row.get("cgmlst_scheme_version"),
            row.get("cgmlst_database_version"),
            row.get("cgmlst_database_sha256"),
            row.get("cgmlst_locus_universe_complete"),
            tuple(sorted(set(row.get("cgmlst_loci") or row["cgmlst_profile"]))),
        )
        for row in profiled
    }
    if len(scopes) != 1:
        raise ValueError("GrapeTree export requires one compatible cgMLST scheme/database/locus scope")
    missing_profiles = [row["sample_id"] for row in rows
                        if not isinstance(row.get("cgmlst_profile"), Mapping)
                        or not row["cgmlst_profile"]]
    if missing_profiles:
        raise ValueError(f"Records lack a cgMLST profile: {missing_profiles[:10]}")
    normalized = {
        row["sample_id"]: {locus: _normalized_call(value)
                           for locus, value in row["cgmlst_profile"].items()}
        for row in rows
    }
    all_missing = [row["sample_id"] for row in rows
                   if not any(normalized[row["sample_id"]].values())]
    if all_missing:
        raise ValueError(f"Records have no called cgMLST loci: {all_missing[:10]}")
    inferred_loci = loci is None
    if inferred_loci:
        loci = sorted({locus for row in profiled for locus in row["cgmlst_profile"]})
    else:
        loci = list(loci)
    if not loci or any(not isinstance(locus, str) or not locus for locus in loci):
        raise ValueError("GrapeTree export requires nonempty locus names")
    if len(set(loci)) != len(loci):
        raise ValueError("GrapeTree locus list contains duplicates")
    if any(set(row["cgmlst_profile"]) - set(loci) for row in profiled):
        raise ValueError("A profile contains loci outside the declared export locus list")

    ids = _export_ids(row["sample_id"] for row in profiled)
    categories: dict[str, dict[str, int]] = {}
    for locus in loci:
        values = sorted({normalized[row["sample_id"]][locus]
        for row in rows
                         if normalized[row["sample_id"]].get(locus) is not None})
        categories[locus] = {value: index + 1 for index, value in enumerate(values)}

    output = _new_output_dir(output_dir)
    profile_path, metadata_path, key_path = (
        output / "profiles.tsv", output / "metadata.tsv", output / "allele_key.csv"
    )
    with profile_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["#Strain", *loci])
        for row in rows:
            calls = []
            for locus in loci:
                raw = normalized[row["sample_id"]].get(locus)
                calls.append("0" if raw is None else str(categories[locus][raw]))
            writer.writerow([ids[row["sample_id"]], *calls])
    metadata_rows = _metadata_records(rows, ids, metadata_fields)
    metadata_columns = ["ID", "chronoclade_sample_id", "sample_label", "source_genome_id",
                        "accession", *metadata_fields]
    _write_delimited(metadata_path, metadata_rows, metadata_columns, "\t")
    key_rows = [{"locus": locus, "export_allele": code, "source_allele": allele}
                for locus, mapping in categories.items()
                for allele, code in mapping.items()]
    _write_delimited(key_path, key_rows, ("locus", "export_allele", "source_allele"), ",")
    readme = output / "README.txt"
    readme.write_text(
        "GrapeTree input files: profiles.tsv and metadata.tsv. IDs are shared.\n"
        "Use metadata.tsv ID as the metadata key and sample_label for readable labels.\n"
        "Alleles are recoded as numeric categories per locus; 0 means missing.\n"
        "allele_key.csv maps each export category back to its source allele.\n"
        + ("The locus list was inferred from observed profiles; it does not establish complete scheme coverage.\n"
           if inferred_loci else "The supplied locus list defines the exported loci.\n"),
        encoding="utf-8",
    )
    return {"profiles": profile_path, "metadata": metadata_path,
            "allele_key": key_path, "readme": readme}


def write_microreact_exports(
    tree_path: str | Path,
    records: Iterable[Mapping],
    output_dir: str | Path,
    *,
    metadata_fields: Sequence[str] = METADATA_FIELDS,
) -> dict[str, Path]:
    """Write a Microreact Newick/CSV pair with matching safe leaf and row IDs."""
    rows = [dict(row) for row in records]
    _source_ids(rows)
    by_id = {row["sample_id"]: row for row in rows}
    unknown_fields = set(metadata_fields) - set(METADATA_FIELDS)
    if unknown_fields:
        raise ValueError(f"Unsupported export metadata fields: {sorted(unknown_fields)}")
    source_tree = Path(tree_path).expanduser().resolve()
    tree = Phylo.read(source_tree, "newick")
    terminals = [tip.name for tip in tree.get_terminals()]
    if any(not isinstance(name, str) or not name for name in terminals):
        raise ValueError("Microreact tree requires a name for every terminal")
    if len(set(terminals)) != len(terminals):
        raise ValueError("Microreact tree contains duplicate terminal IDs")
    missing = sorted(set(terminals) - set(by_id))
    if missing:
        raise ValueError(f"Tree tips have no matching metadata records: {missing[:5]}")
    selected = [by_id[name] for name in terminals]
    ids = _export_ids(terminals)
    for tip in tree.get_terminals():
        tip.name = ids[tip.name]

    output = _new_output_dir(output_dir)
    safe_tree = output / "tree.newick"
    metadata = output / "metadata.csv"
    id_map = output / "id_map.csv"
    Phylo.write(tree, safe_tree, "newick", format_branch_length="%.17g")
    metadata_rows = _metadata_records(selected, ids, metadata_fields)
    columns = ["ID", "chronoclade_sample_id", "sample_label", "source_genome_id",
               "accession", *metadata_fields]
    _write_delimited(metadata, metadata_rows, columns, ",")
    _write_delimited(id_map, [{"export_id": ids[ident], "sample_id": ident}
                              for ident in terminals], ("export_id", "sample_id"), ",")
    readme = output / "README.txt"
    readme.write_text(
        "Upload tree.newick and metadata.csv to Microreact. Select ID as the shared key.\n"
        "Tree leaf labels and metadata IDs are deterministic safe IDs.\n"
        "id_map.csv retains the source sample IDs; internal Newick support labels are retained.\n",
        encoding="utf-8",
    )
    return {"tree": safe_tree, "metadata": metadata, "id_map": id_map, "readme": readme}
