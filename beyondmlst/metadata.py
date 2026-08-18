"""Metadata parsing and validation for beyondMLST."""

from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable, TextIO

REQUIRED_COLUMNS = {
    "sample_id",
    "assembly",
    "collection_date",
    "location",
    "species",
    "lineage",
    "origin",
}
SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9_.-]+$")
YEAR = re.compile(r"^(\d{4})$")
YEAR_MONTH = re.compile(r"^(\d{4})-(\d{2})$")
FULL_DATE = re.compile(r"^(\d{4})-(\d{2})-(\d{2})$")
AMBIGUOUS_DATE = re.compile(r"^\d{4}-(?:\d{2}|XX)-(?:\d{2}|XX)$")
NUMERIC_DATE = re.compile(r"^\d{4}(?:\.\d+)?$")
RANGE_DATE = re.compile(r"^\[\d{4}(?:\.\d+)?:\d{4}(?:\.\d+)?\]$")
TRUE_VALUES = {"1", "true", "yes", "y"}
FALSE_VALUES = {"", "0", "false", "no", "n"}


class MetadataError(ValueError):
    """Raised when input metadata are unsafe or incomplete."""


@dataclass(frozen=True)
class Sample:
    """One genome and its analysis metadata."""

    sample_id: str
    assembly: Path
    collection_date: str
    location: str
    species: str
    lineage: str
    origin: str
    is_reference: bool = False
    patient_id: str = ""

    @property
    def lineage_key(self) -> tuple[str, str]:
        return (self.species, self.lineage)


def _parse_boolean(value: str, *, row_number: int) -> bool:
    normalised = value.strip().lower()
    if normalised in TRUE_VALUES:
        return True
    if normalised in FALSE_VALUES:
        return False
    raise MetadataError(f"Row {row_number}: is_reference must be true/false, not {value!r}")


def _validate_date(value: str, *, row_number: int) -> None:
    if YEAR.fullmatch(value) or NUMERIC_DATE.fullmatch(value) or RANGE_DATE.fullmatch(value):
        return
    match = YEAR_MONTH.fullmatch(value)
    if match:
        year, month = map(int, match.groups())
        if 1 <= month <= 12 and 1800 <= year <= 2200:
            return
    match = FULL_DATE.fullmatch(value)
    if match:
        try:
            date(*map(int, match.groups()))
            return
        except ValueError:
            pass
    if AMBIGUOUS_DATE.fullmatch(value):
        return
    raise MetadataError(
        f"Row {row_number}: unsupported collection_date {value!r}; use YYYY, "
        "YYYY-MM, YYYY-MM-DD, YYYY-XX-XX, a numeric date, or [start:end]"
    )


def read_metadata(path: Path) -> list[Sample]:
    """Read and validate the study metadata CSV."""

    path = path.expanduser().resolve()
    if not path.is_file():
        raise MetadataError(f"Metadata file does not exist: {path}")

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise MetadataError("Metadata CSV has no header")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise MetadataError(f"Metadata CSV is missing columns: {', '.join(sorted(missing))}")

        samples: list[Sample] = []
        seen: set[str] = set()
        for row_number, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            if not SAFE_IDENTIFIER.fullmatch(sample_id):
                raise MetadataError(
                    f"Row {row_number}: sample_id {sample_id!r} must contain only letters, "
                    "numbers, '.', '_' or '-'"
                )
            if sample_id in seen:
                raise MetadataError(f"Row {row_number}: duplicate sample_id {sample_id!r}")
            seen.add(sample_id)

            assembly_text = (row.get("assembly") or "").strip()
            assembly = Path(assembly_text).expanduser()
            if not assembly.is_absolute():
                assembly = path.parent / assembly
            assembly = assembly.resolve()
            if not assembly.is_file():
                raise MetadataError(
                    f"Row {row_number}: assembly for {sample_id!r} does not exist: {assembly}"
                )

            collection_date = (row.get("collection_date") or "").strip()
            _validate_date(collection_date, row_number=row_number)

            values = {
                key: (row.get(key) or "").strip()
                for key in ("location", "species", "lineage", "origin")
            }
            empty = [key for key, value in values.items() if not value]
            if empty:
                raise MetadataError(f"Row {row_number}: empty values for {', '.join(empty)}")

            samples.append(
                Sample(
                    sample_id=sample_id,
                    assembly=assembly,
                    collection_date=collection_date,
                    location=values["location"],
                    species=values["species"],
                    lineage=values["lineage"],
                    origin=values["origin"],
                    is_reference=_parse_boolean(
                        row.get("is_reference") or "", row_number=row_number
                    ),
                    patient_id=(row.get("patient_id") or "").strip(),
                )
            )

    if not samples:
        raise MetadataError("Metadata CSV contains no samples")
    return samples


def group_samples(samples: Iterable[Sample]) -> dict[tuple[str, str], list[Sample]]:
    """Group samples into independently analysed species/lineage sets."""

    groups: dict[tuple[str, str], list[Sample]] = {}
    for sample in samples:
        groups.setdefault(sample.lineage_key, []).append(sample)
    return groups


def _open_fasta(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def fasta_lengths(path: Path) -> list[int]:
    """Return contig lengths from a plain or gzipped FASTA assembly."""

    lengths: list[int] = []
    current = 0
    with _open_fasta(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                if not lengths and current == 0 and line.startswith("@"):  # obvious FASTQ
                    raise MetadataError(
                        f"Expected an assembly FASTA, found FASTQ-like data: {path}"
                    )
                current += len(line)
    if current:
        lengths.append(current)
    if not lengths:
        raise MetadataError(f"Assembly contains no FASTA sequences: {path}")
    return lengths


def assembly_n50(path: Path) -> int:
    """Calculate assembly N50 for deterministic fallback reference selection."""

    lengths = sorted(fasta_lengths(path), reverse=True)
    threshold = sum(lengths) / 2
    cumulative = 0
    for length in lengths:
        cumulative += length
        if cumulative >= threshold:
            return length
    raise AssertionError("unreachable")


def assembly_length(path: Path) -> int:
    return sum(fasta_lengths(path))


def select_reference(samples: list[Sample]) -> Sample:
    """Use an explicit reference or select the highest-N50 assembly."""

    explicit = [sample for sample in samples if sample.is_reference]
    if len(explicit) > 1:
        names = ", ".join(sample.sample_id for sample in explicit)
        raise MetadataError(f"Multiple references selected in one lineage: {names}")
    if explicit:
        return explicit[0]
    return max(samples, key=lambda sample: (assembly_n50(sample.assembly), sample.sample_id))


def slugify_lineage(species: str, lineage: str) -> str:
    value = f"{species}__{lineage}"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
