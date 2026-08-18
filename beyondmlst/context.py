"""Reproducible acquisition and down-selection of public context genomes."""

from __future__ import annotations

import csv
import hashlib
import json
import random
import re
import shutil
import sqlite3
import subprocess
import urllib.error
import urllib.request
from collections import defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

from beyondmlst.metadata import Sample


OFFICIAL_ATB_MLST_URL = "https://osf.io/download/69c66d33fa3d973d94254f46/"


class ContextError(RuntimeError):
    """Raised when contextual-genome preparation cannot be completed safely."""


@dataclass
class ContextCandidate:
    """One same-ST public genome considered for contextual analysis."""

    sample_id: str
    species: str
    lineage: str
    mlst_scheme: str
    mlst_st: str
    collection_date: str = ""
    country: str = ""
    host: str = ""
    isolation_source: str = ""
    hq_filter: str = ""
    aws_url: str = ""
    completeness: str = ""
    contamination: str = ""
    genome_size: str = ""
    contig_n50: str = ""
    assembly: str = ""
    nearest_focal: str = ""
    min_ska_distance: float | None = None
    mismatch_proportion: float | None = None
    selection_reason: str = ""

    @property
    def year(self) -> str:
        return self.collection_date[:4] if len(self.collection_date) >= 4 else "Unknown"

    @property
    def stratum(self) -> tuple[str, str]:
        return (self.country or "Unknown", self.year)


MANIFEST_FIELDS = [
    "sample_id",
    "species",
    "lineage",
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
    "assembly",
    "nearest_focal",
    "min_ska_distance",
    "mismatch_proportion",
    "selection_reason",
    "aws_url",
]


def _run_capture(
    command: list[str], *, log: Path | None = None
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if log is not None:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(
            "COMMAND\n"
            + " ".join(command)
            + "\n\nSTDOUT\n"
            + completed.stdout
            + "\nSTDERR\n"
            + completed.stderr,
            encoding="utf-8",
        )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip() or "no diagnostic output"
        raise ContextError(f"Command failed ({command[0]}): {detail}")
    return completed


def atbfetcher_version(executable: str = "atbfetcher") -> str:
    """Return the installed atbfetcher version string."""

    completed = _run_capture([executable, "--version"])
    return completed.stdout.strip() or completed.stderr.strip() or "unknown"


def discover_same_st(
    *,
    species: str,
    scheme: str,
    st: str,
    cache_dir: Path,
    executable: str = "atbfetcher",
    log: Path | None = None,
) -> list[str]:
    """Use atbfetcher to return all high-quality samples assigned to one ST."""

    command = [
        executable,
        "mlst-query",
        "--species",
        species,
        "--scheme",
        scheme,
        "--st",
        st,
        "--format",
        "tsv",
        "--cache-dir",
        str(cache_dir),
    ]
    completed = _run_capture(command, log=log)
    lines = completed.stdout.splitlines()
    header_index = next(
        (
            index
            for index, line in enumerate(lines)
            if "sample" in {field.strip() for field in line.split("\t")}
        ),
        None,
    )
    if header_index is None:
        raise ContextError("atbfetcher MLST output did not contain a 'sample' column")
    reader = csv.DictReader(lines[header_index:], delimiter="\t")
    if reader.fieldnames is None or "sample" not in reader.fieldnames:
        raise ContextError("atbfetcher MLST output did not contain a 'sample' column")
    accessions = sorted({(row.get("sample") or "").strip() for row in reader})
    accessions = [accession for accession in accessions if accession]
    if not accessions:
        raise ContextError(f"No high-quality {species} ST{st} candidates were found")
    return accessions


def cache_official_atb_mlst(cache_dir: Path) -> Path:
    """Cache the current official ATB MLST Parquet used by atbfetcher.

    atbfetcher's pinned release still points to a retired R2 object. The current
    ATB CLI publishes the same MLST schema from OSF, so placing that file at the
    cache location expected by atbfetcher lets its normal query path continue.
    """

    cache_dir = cache_dir.expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    destination = cache_dir / "mlst.parquet"
    provenance = cache_dir / "mlst.parquet.provenance.json"
    if destination.is_file() and destination.stat().st_size > 8:
        with destination.open("rb") as handle:
            if handle.read(4) == b"PAR1":
                if not provenance.is_file():
                    _write_mlst_provenance(destination, provenance)
                return destination

    temporary = cache_dir / "mlst.parquet.part"
    request = urllib.request.Request(
        OFFICIAL_ATB_MLST_URL,
        headers={"User-Agent": "beyondmlst/0.1 (+https://github.com/ghruproject/beyondmlst)"},
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:  # noqa: S310
            with temporary.open("wb") as handle:
                shutil.copyfileobj(response, handle)
        if temporary.stat().st_size <= 8:
            raise ContextError("The official ATB MLST download was unexpectedly empty")
        with temporary.open("rb") as handle:
            if handle.read(4) != b"PAR1":
                raise ContextError("The official ATB MLST download was not a Parquet file")
        temporary.replace(destination)
        _write_mlst_provenance(destination, provenance)
    except (OSError, urllib.error.URLError) as error:
        temporary.unlink(missing_ok=True)
        if isinstance(error, ContextError):
            raise
        raise ContextError(f"Could not download the official ATB MLST table: {error}") from error
    return destination


def _write_mlst_provenance(table: Path, destination: Path) -> None:
    with table.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    destination.write_text(
        json.dumps(
            {
                "source": "official_allthebacteria_osf",
                "url": OFFICIAL_ATB_MLST_URL,
                "path": str(table),
                "size_bytes": table.stat().st_size,
                "sha256": digest,
                "cached_at": datetime.now().astimezone().isoformat(),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def discover_same_st_resilient(
    *,
    species: str,
    scheme: str,
    st: str,
    cache_dir: Path,
    executable: str = "atbfetcher",
    logs: Path | None = None,
) -> tuple[list[str], dict[str, str]]:
    """Query atbfetcher, repairing only its known missing MLST cache path."""

    initial_log = logs / "atbfetcher_mlst_query_initial.log" if logs else None
    try:
        accessions = discover_same_st(
            species=species,
            scheme=scheme,
            st=st,
            cache_dir=cache_dir,
            executable=executable,
            log=initial_log,
        )
        provenance = cache_dir.expanduser().resolve() / "mlst.parquet.provenance.json"
        if provenance.is_file():
            try:
                details = json.loads(provenance.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                details = {}
            return accessions, {
                "method": "atbfetcher_with_cached_official_atb_mlst",
                "mlst_table": str(cache_dir.expanduser().resolve() / "mlst.parquet"),
                "mlst_table_url": str(details.get("url", OFFICIAL_ATB_MLST_URL)),
                "mlst_table_sha256": str(details.get("sha256", "")),
            }
        return accessions, {"method": "atbfetcher", "mlst_table": "atbfetcher cache"}
    except ContextError as initial_error:
        cached_table = cache_dir.expanduser().resolve() / "mlst.parquet"
        if cached_table.is_file():
            raise
        cached_table = cache_official_atb_mlst(cache_dir)
        retry_log = logs / "atbfetcher_mlst_query_retry.log" if logs else None
        accessions = discover_same_st(
            species=species,
            scheme=scheme,
            st=st,
            cache_dir=cache_dir,
            executable=executable,
            log=retry_log,
        )
        return accessions, {
            "method": "atbfetcher_with_official_atb_mlst_fallback",
            "mlst_table": str(cached_table),
            "mlst_table_url": OFFICIAL_ATB_MLST_URL,
            "fallback_reason": str(initial_error),
        }


def normalise_collection_date(value: str) -> str:
    """Return a beyondMLST-compatible public collection date or an empty string."""

    value = value.strip()
    if not value:
        return ""
    value = value.split("T", 1)[0].split(" ", 1)[0]
    patterns = (
        (r"^(\d{4})$", "%Y"),
        (r"^(\d{4})-(\d{2})$", "%Y-%m"),
        (r"^(\d{4})-(\d{2})-(\d{2})$", "%Y-%m-%d"),
    )
    for pattern, date_format in patterns:
        if not re.fullmatch(pattern, value):
            continue
        try:
            parsed = datetime.strptime(value, date_format)
        except ValueError:
            return ""
        if 1800 <= parsed.year <= 2200:
            return value
    return ""


def _normalise_species(value: str) -> str:
    """Normalise GTDB letter suffixes for cross-table species matching."""

    return re.sub(r"_[A-Z]\b", "", value).casefold().strip()


def _chunks(values: list[str], size: int = 400) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def lookup_atb_metadata(
    db_path: Path,
    accessions: list[str],
    *,
    species: str,
    lineage: str,
    scheme: str,
    st: str,
) -> list[ContextCandidate]:
    """Look up ATB/ENA metadata for accessions discovered by atbfetcher.

    The current atbfetcher CLI intentionally prints only accessions for metadata
    queries, so this adapter reads the same downloaded ATB SQLite snapshot to
    preserve the metadata needed for reproducible selection and reporting.
    """

    db_path = db_path.expanduser().resolve()
    if not db_path.is_file():
        raise ContextError(
            f"ATB metadata database not found: {db_path}. Run `atbfetcher download-db`."
        )
    accession_set = set(accessions)
    assembly_by_sample: dict[str, sqlite3.Row] = {}
    ena_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    connection = sqlite3.connect(str(db_path))
    connection.row_factory = sqlite3.Row
    try:
        for chunk in _chunks(accessions):
            placeholders = ",".join("?" for _ in chunk)
            sql = f"""
                SELECT a.sample_accession AS sample, a.sylph_species AS species,
                       a.hq_filter, a.aws_url, a.asm_fasta_on_osf,
                       c.Completeness_Specific AS completeness,
                       c.Contamination AS contamination,
                       c.Genome_Size AS genome_size,
                       c.Contig_N50 AS contig_n50
                  FROM assembly a
             LEFT JOIN checkm2 c ON a.sample_accession = c.sample_accession
                 WHERE a.sample_accession IN ({placeholders})
            """
            for row in connection.execute(sql, chunk):
                assembly_by_sample[str(row["sample"])] = row

        ena_tables = sorted(
            (
                str(row["name"])
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type = 'table' AND name LIKE 'ena_%'"
                )
                if re.fullmatch(r"ena_\d{8}", str(row["name"]))
            ),
            reverse=True,
        )
        if not ena_tables:
            ena_tables = ["ena_202505_used"]
        ena_table = ena_tables[0]
        ena_columns = {
            str(row["name"]) for row in connection.execute(f"PRAGMA table_info({ena_table})")
        }
        run_to_sample: dict[str, str] = {}
        if "sample_accession" not in ena_columns:
            for row in connection.execute("SELECT run_accession, sample_accession FROM run"):
                sample_id = str(row["sample_accession"] or "")
                if sample_id in accession_set:
                    run_to_sample[str(row["run_accession"])] = sample_id

        ena_sql = "SELECT e.run_accession, e.country, e.collection_date, e.host, e.isolation_source"
        if "sample_accession" in ena_columns:
            connection.execute(
                "CREATE TEMP TABLE requested_accession (sample_accession TEXT PRIMARY KEY)"
            )
            connection.executemany(
                "INSERT INTO requested_accession VALUES (?)",
                ((accession,) for accession in accessions),
            )
            ena_sql += (
                f", e.sample_accession FROM {ena_table} e "
                "JOIN requested_accession q ON e.sample_accession = q.sample_accession"
            )
        else:
            ena_sql += f" FROM {ena_table} e"
        for row in connection.execute(ena_sql):
            if "sample_accession" in ena_columns:
                sample_id = str(row["sample_accession"] or "")
            else:
                sample_id = run_to_sample.get(str(row["run_accession"]), "")
            if sample_id not in accession_set:
                continue
            ena_by_sample[sample_id].append(
                {
                    "run_accession": str(row["run_accession"] or ""),
                    "country": str(row["country"] or ""),
                    "collection_date": str(row["collection_date"] or ""),
                    "host": str(row["host"] or ""),
                    "isolation_source": str(row["isolation_source"] or ""),
                }
            )
    except sqlite3.Error as error:
        raise ContextError(f"Could not query the ATB metadata snapshot: {error}") from error
    finally:
        connection.close()

    candidates: list[ContextCandidate] = []
    for accession in accessions:
        assembly_row = assembly_by_sample.get(accession)
        if assembly_row is None:
            continue
        metadata_rows = ena_by_sample.get(accession, [])
        metadata_rows.sort(
            key=lambda row: (
                not bool(normalise_collection_date(str(row["collection_date"] or ""))),
                str(row["run_accession"] or ""),
            )
        )
        if not metadata_rows:
            continue
        metadata_row = metadata_rows[0]
        if str(assembly_row["hq_filter"] or "") != "PASS" or not bool(
            assembly_row["asm_fasta_on_osf"]
        ):
            continue
        if _normalise_species(str(assembly_row["species"] or "")) != _normalise_species(species):
            continue
        collection_date = normalise_collection_date(str(metadata_row["collection_date"] or ""))
        if not collection_date:
            continue
        candidates.append(
            ContextCandidate(
                sample_id=accession,
                species=species,
                lineage=lineage,
                mlst_scheme=scheme,
                mlst_st=st,
                collection_date=collection_date,
                country=metadata_row["country"],
                host=metadata_row["host"],
                isolation_source=metadata_row["isolation_source"],
                hq_filter=str(assembly_row["hq_filter"] or ""),
                aws_url=str(assembly_row["aws_url"] or ""),
                completeness=str(assembly_row["completeness"] or ""),
                contamination=str(assembly_row["contamination"] or ""),
                genome_size=str(assembly_row["genome_size"] or ""),
                contig_n50=str(assembly_row["contig_n50"] or ""),
            )
        )
    return candidates


def filter_candidates(
    candidates: list[ContextCandidate],
    *,
    countries: list[str] | None = None,
    year_from: int | None = None,
    year_to: int | None = None,
    host: str | None = None,
    isolation_source: str | None = None,
) -> list[ContextCandidate]:
    """Apply explicit epidemiological filters to same-ST candidates."""

    country_values = [value.casefold() for value in countries or []]
    host_value = host.casefold() if host else ""
    source_value = isolation_source.casefold() if isolation_source else ""
    selected: list[ContextCandidate] = []
    for candidate in candidates:
        year = int(candidate.year)
        if country_values and not any(
            candidate.country.casefold().startswith(value) for value in country_values
        ):
            continue
        if year_from is not None and year < year_from:
            continue
        if year_to is not None and year > year_to:
            continue
        if host_value and host_value not in candidate.host.casefold():
            continue
        if source_value and source_value not in candidate.isolation_source.casefold():
            continue
        selected.append(candidate)
    return selected


def stratified_candidate_pool(
    candidates: list[ContextCandidate], *, limit: int, seed: int
) -> list[ContextCandidate]:
    """Select a reproducible country/year-balanced candidate pool before download."""

    if limit < 1:
        raise ValueError("candidate pool limit must be at least 1")
    if len(candidates) <= limit:
        return sorted(candidates, key=lambda candidate: candidate.sample_id)

    rng = random.Random(seed)
    strata: dict[tuple[str, str], list[ContextCandidate]] = defaultdict(list)
    for candidate in candidates:
        strata[candidate.stratum].append(candidate)
    for members in strata.values():
        members.sort(key=lambda candidate: candidate.sample_id)
        rng.shuffle(members)
    keys = sorted(strata)
    rng.shuffle(keys)

    result: list[ContextCandidate] = []
    while len(result) < limit:
        progressed = False
        for key in keys:
            if strata[key]:
                result.append(strata[key].pop())
                progressed = True
                if len(result) == limit:
                    break
        if not progressed:
            break
    return result


def write_accessions(path: Path, candidates: list[ContextCandidate]) -> Path:
    path.write_text(
        "".join(f"{candidate.sample_id}\n" for candidate in candidates), encoding="utf-8"
    )
    return path


def download_candidates(
    candidates: list[ContextCandidate],
    *,
    accessions_path: Path,
    output: Path,
    cache_dir: Path,
    source: str,
    threads: int,
    executable: str = "atbfetcher",
    log: Path | None = None,
) -> None:
    """Fetch a frozen candidate accession list with atbfetcher."""

    if not candidates:
        raise ContextError("No context candidates remain to download")
    output.mkdir(parents=True, exist_ok=True)
    command = [
        executable,
        "accessions",
        str(accessions_path),
        "--output",
        str(output),
        "--cache-dir",
        str(cache_dir),
        "--source",
        source,
        "--threads",
        str(threads),
    ]
    _run_capture(command, log=log)


def _assembly_identifier(path: Path) -> str:
    name = path.name
    if name.endswith(".gz"):
        name = name[:-3]
    for extension in (".fasta", ".fna", ".fa"):
        if name.endswith(extension):
            return name[: -len(extension)]
    return ""


def attach_downloaded_assemblies(
    candidates: list[ContextCandidate], output: Path
) -> tuple[list[ContextCandidate], list[str]]:
    """Attach downloaded FASTA paths and report missing accessions."""

    paths: dict[str, Path] = {}
    for path in output.rglob("*"):
        if not path.is_file():
            continue
        identifier = _assembly_identifier(path)
        if identifier:
            paths.setdefault(identifier, path.resolve())
    found: list[ContextCandidate] = []
    missing: list[str] = []
    for candidate in candidates:
        assembly = paths.get(candidate.sample_id)
        if assembly is None:
            missing.append(candidate.sample_id)
            continue
        candidate.assembly = str(assembly)
        found.append(candidate)
    return found, missing


def write_ska_inputs(path: Path, focal: list[Sample], candidates: list[ContextCandidate]) -> Path:
    rows = [f"{sample.sample_id}\t{sample.assembly}\n" for sample in focal]
    rows.extend(f"{candidate.sample_id}\t{candidate.assembly}\n" for candidate in candidates)
    path.write_text("".join(rows), encoding="utf-8")
    return path


def run_ska_screen(
    *, inputs: Path, output_dir: Path, threads: int, executable: str = "ska"
) -> Path:
    """Build one SKA index and calculate pairwise screening distances."""

    prefix = output_dir / "context_screen"
    skf = output_dir / "context_screen.skf"
    distances = output_dir / "context_distances.tsv"
    logs = output_dir / "logs"
    _run_capture(
        [executable, "build", "-f", str(inputs), "-o", str(prefix), "--threads", str(threads)],
        log=logs / "ska_build.log",
    )
    if not skf.is_file():
        raise ContextError(f"SKA did not create its expected index: {skf}")
    _run_capture(
        [
            executable,
            "distance",
            "--threads",
            str(threads),
            "-o",
            str(distances),
            str(skf),
        ],
        log=logs / "ska_distance.log",
    )
    if not distances.is_file():
        raise ContextError(f"SKA did not create its distance table: {distances}")
    return distances


def parse_ska_distances(path: Path) -> dict[tuple[str, str], tuple[float, float]]:
    """Parse SKA2 pairwise SNP and mismatch-proportion output."""

    result: dict[tuple[str, str], tuple[float, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Sample1", "Sample2", "Distance", "Mismatches (proportion)"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ContextError(f"Unrecognised SKA distance format: {path}")
        for row in reader:
            first = str(row["Sample1"])
            second = str(row["Sample2"])
            value = (float(row["Distance"]), float(row["Mismatches (proportion)"]))
            result[(first, second)] = value
            result[(second, first)] = value
    return result


def annotate_screening_distances(
    candidates: list[ContextCandidate],
    focal_ids: list[str],
    distances: dict[tuple[str, str], tuple[float, float]],
) -> list[ContextCandidate]:
    """Annotate every candidate with its nearest focal SKA screening distance."""

    annotated: list[ContextCandidate] = []
    for candidate in candidates:
        comparisons = [
            (distances[(candidate.sample_id, focal_id)], focal_id)
            for focal_id in focal_ids
            if (candidate.sample_id, focal_id) in distances
        ]
        if not comparisons:
            continue
        (distance, mismatch), focal_id = min(
            comparisons, key=lambda value: (value[0][0], value[0][1], value[1])
        )
        candidate.min_ska_distance = distance
        candidate.mismatch_proportion = mismatch
        candidate.nearest_focal = focal_id
        annotated.append(candidate)
    return annotated


def select_context(
    candidates: list[ContextCandidate],
    focal_ids: list[str],
    distances: dict[tuple[str, str], tuple[float, float]],
    *,
    max_context: int,
    nearest_per_focal: int,
    seed: int,
) -> tuple[list[ContextCandidate], dict[str, object]]:
    """Retain nearest neighbours, then fill with a balanced background."""

    if max_context < 1 or nearest_per_focal < 1:
        raise ValueError("max_context and nearest_per_focal must be at least 1")
    candidates = annotate_screening_distances(candidates, focal_ids, distances)
    by_id = {candidate.sample_id: candidate for candidate in candidates}
    nearest_for: dict[str, list[str]] = defaultdict(list)
    support: dict[str, int] = defaultdict(int)
    for focal_id in focal_ids:
        ranked = sorted(
            (
                (distances[(candidate.sample_id, focal_id)][0], candidate.sample_id)
                for candidate in candidates
                if (candidate.sample_id, focal_id) in distances
            ),
            key=lambda value: (value[0], value[1]),
        )[:nearest_per_focal]
        for _, candidate_id in ranked:
            nearest_for[candidate_id].append(focal_id)
            support[candidate_id] += 1

    selected_ids = sorted(
        nearest_for,
        key=lambda candidate_id: (
            -support[candidate_id],
            (
                by_id[candidate_id].min_ska_distance
                if by_id[candidate_id].min_ska_distance is not None
                else float("inf")
            ),
            candidate_id,
        ),
    )[:max_context]
    selected_set = set(selected_ids)
    for candidate_id in selected_ids:
        by_id[candidate_id].selection_reason = "nearest_to=" + "|".join(
            sorted(nearest_for[candidate_id])
        )

    remaining: dict[tuple[str, str], list[ContextCandidate]] = defaultdict(list)
    for candidate in candidates:
        if candidate.sample_id not in selected_set:
            remaining[candidate.stratum].append(candidate)
    for members in remaining.values():
        members.sort(
            key=lambda candidate: (
                (
                    candidate.min_ska_distance
                    if candidate.min_ska_distance is not None
                    else float("inf")
                ),
                candidate.sample_id,
            ),
            reverse=True,
        )
    keys = sorted(remaining)
    random.Random(seed).shuffle(keys)
    while len(selected_ids) < max_context:
        progressed = False
        for key in keys:
            if remaining[key]:
                candidate = remaining[key].pop()
                candidate.selection_reason = (
                    f"stratified_background={candidate.country or 'Unknown'}|{candidate.year}"
                )
                selected_ids.append(candidate.sample_id)
                selected_set.add(candidate.sample_id)
                progressed = True
                if len(selected_ids) == max_context:
                    break
        if not progressed:
            break

    selected = [by_id[candidate_id] for candidate_id in selected_ids]
    selected.sort(
        key=lambda candidate: (
            (
                candidate.min_ska_distance
                if candidate.min_ska_distance is not None
                else float("inf")
            ),
            candidate.sample_id,
        )
    )
    covered = {
        focal_id for candidate in selected for focal_id in nearest_for.get(candidate.sample_id, [])
    }
    audit = {
        "screened_candidates": len(candidates),
        "selected_contexts": len(selected),
        "max_context": max_context,
        "nearest_per_focal": nearest_per_focal,
        "focal_samples": len(focal_ids),
        "focal_samples_with_selected_neighbour": len(covered),
        "focal_neighbour_coverage": len(covered) / len(focal_ids) if focal_ids else 0.0,
    }
    return selected, audit


def write_candidate_table(path: Path, candidates: list[ContextCandidate]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MANIFEST_FIELDS,
            extrasaction="ignore",
            delimiter="\t",
        )
        writer.writeheader()
        for candidate in candidates:
            writer.writerow(asdict(candidate))
    return path


def write_combined_metadata(
    path: Path, focal: list[Sample], selected: list[ContextCandidate]
) -> Path:
    fields = [
        "sample_id",
        "assembly",
        "collection_date",
        "location",
        "species",
        "lineage",
        "origin",
        "is_reference",
        "patient_id",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for sample in focal:
            writer.writerow(
                {
                    "sample_id": sample.sample_id,
                    "assembly": sample.assembly,
                    "collection_date": sample.collection_date,
                    "location": sample.location,
                    "species": sample.species,
                    "lineage": sample.lineage,
                    "origin": sample.origin,
                    "is_reference": str(sample.is_reference).lower(),
                    "patient_id": sample.patient_id,
                }
            )
        for candidate in selected:
            writer.writerow(
                {
                    "sample_id": candidate.sample_id,
                    "assembly": candidate.assembly,
                    "collection_date": candidate.collection_date,
                    "location": candidate.country or "Public_context",
                    "species": candidate.species,
                    "lineage": candidate.lineage,
                    "origin": "context",
                    "is_reference": "false",
                    "patient_id": "",
                }
            )
    return path


def write_audit(path: Path, audit: dict[str, object]) -> Path:
    path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return path


def find_atb_database(cache_dir: Path, requested: Path | None = None) -> Path:
    """Resolve an explicit or cached ATB SQLite metadata snapshot."""

    if requested is not None:
        result = requested.expanduser().resolve()
        if result.is_file():
            return result
        raise ContextError(f"ATB metadata database not found: {result}")
    cache_dir = cache_dir.expanduser().resolve()
    candidates = sorted(cache_dir.glob("atb.metadata.*.sqlite"), reverse=True)
    if not candidates:
        raise ContextError(
            "ATB metadata database not found. Run `atbfetcher download-db` or pass --db-path."
        )
    return candidates[0]


def prepare_context(
    focal: list[Sample],
    *,
    species: str,
    lineage: str,
    scheme: str,
    st: str,
    output: Path,
    cache_dir: Path,
    db_path: Path | None,
    countries: list[str] | None,
    year_from: int | None,
    year_to: int | None,
    host: str | None,
    isolation_source: str | None,
    candidate_pool: int,
    max_context: int,
    nearest_per_focal: int,
    seed: int,
    threads: int,
    source: str,
    dry_run: bool,
    atbfetcher_executable: str = "atbfetcher",
    ska_executable: str = "ska",
) -> dict[str, object]:
    """Prepare a frozen, distance-screened public context set for one lineage."""

    if not focal:
        raise ContextError("No focal samples were supplied for context preparation")
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    logs = output / "logs"
    cache_dir = cache_dir.expanduser().resolve()
    resolved_db = find_atb_database(cache_dir, db_path)
    tool_version = atbfetcher_version(atbfetcher_executable)

    raw_accessions, discovery = discover_same_st_resilient(
        species=species,
        scheme=scheme,
        st=st,
        cache_dir=cache_dir,
        executable=atbfetcher_executable,
        logs=logs,
    )
    focal_ids = {sample.sample_id for sample in focal}
    raw_accessions = [accession for accession in raw_accessions if accession not in focal_ids]
    (output / "same_st_accessions.txt").write_text(
        "".join(f"{accession}\n" for accession in raw_accessions), encoding="utf-8"
    )

    metadata_candidates = lookup_atb_metadata(
        resolved_db,
        raw_accessions,
        species=species,
        lineage=lineage,
        scheme=scheme,
        st=st,
    )
    filtered = filter_candidates(
        metadata_candidates,
        countries=countries,
        year_from=year_from,
        year_to=year_to,
        host=host,
        isolation_source=isolation_source,
    )
    if not filtered:
        raise ContextError("No dated same-ST candidates passed the requested metadata filters")
    pool = stratified_candidate_pool(filtered, limit=candidate_pool, seed=seed)
    write_candidate_table(output / "candidate_pool.tsv", pool)
    accessions_path = write_accessions(output / "candidate_accessions.txt", pool)

    audit: dict[str, object] = {
        "species": species,
        "lineage": lineage,
        "mlst_scheme": scheme,
        "mlst_st": st,
        "atbfetcher_version": tool_version,
        "mlst_discovery": discovery,
        "atb_metadata_snapshot": {
            "path": str(resolved_db),
            "name": resolved_db.name,
            "size_bytes": resolved_db.stat().st_size,
            "modified": datetime.fromtimestamp(resolved_db.stat().st_mtime).isoformat(),
        },
        "filters": {
            "countries": countries or [],
            "year_from": year_from,
            "year_to": year_to,
            "host": host,
            "isolation_source": isolation_source,
        },
        "seed": seed,
        "same_st_accessions": len(raw_accessions),
        "dated_hq_metadata_candidates": len(metadata_candidates),
        "metadata_filtered_candidates": len(filtered),
        "candidate_pool_limit": candidate_pool,
        "candidate_pool": len(pool),
        "dry_run": dry_run,
    }
    if dry_run:
        audit["next_step"] = "Rerun without --dry-run to fetch, screen and select assemblies."
        write_audit(output / "context_selection.json", audit)
        return audit

    assemblies = output / "assemblies"
    download_candidates(
        pool,
        accessions_path=accessions_path,
        output=assemblies,
        cache_dir=cache_dir,
        source=source,
        threads=threads,
        executable=atbfetcher_executable,
        log=logs / "atbfetcher_download.log",
    )
    downloaded, missing = attach_downloaded_assemblies(pool, assemblies)
    if not downloaded:
        raise ContextError("atbfetcher completed but no candidate FASTA assemblies were found")

    ska_inputs = write_ska_inputs(output / "ska_inputs.tsv", focal, downloaded)
    distance_path = run_ska_screen(
        inputs=ska_inputs,
        output_dir=output,
        threads=threads,
        executable=ska_executable,
    )
    distance_values = parse_ska_distances(distance_path)
    selected, selection_audit = select_context(
        downloaded,
        [sample.sample_id for sample in focal],
        distance_values,
        max_context=max_context,
        nearest_per_focal=nearest_per_focal,
        seed=seed,
    )
    if not selected:
        raise ContextError("SKA screening did not yield any usable contextual genomes")

    screened = annotate_screening_distances(
        downloaded, [sample.sample_id for sample in focal], distance_values
    )
    write_candidate_table(output / "screened_candidates.tsv", screened)
    manifest = write_candidate_table(output / "context_manifest.tsv", selected)
    combined = write_combined_metadata(output / "combined_metadata.csv", focal, selected)
    audit.update(selection_audit)
    audit.update(
        {
            "downloaded_candidates": len(downloaded),
            "missing_downloads": missing,
            "source": source,
            "outputs": {
                "manifest": str(manifest),
                "combined_metadata": str(combined),
                "ska_distances": str(distance_path),
            },
        }
    )
    write_audit(output / "context_selection.json", audit)
    return audit
