from pathlib import Path

from chronoclade.evidence import (
    pairwise_distances,
    read_alignment,
)
from chronoclade.metadata import Sample


def samples(tmp_path: Path, *, include_third_local: bool = False) -> list[Sample]:
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig\nACGT\n", encoding="utf-8")
    result = [
        Sample("L1", assembly, "2023-01", "KIMS", "E_coli", "ST1", "local", patient_id="P1"),
        Sample("L2", assembly, "2024-02", "KIMS", "E_coli", "ST1", "local", patient_id="P1"),
    ]
    if include_third_local:
        result.append(
            Sample("L3", assembly, "2024-06", "KIMS", "E_coli", "ST1", "local", patient_id="P2")
        )
    result.append(Sample("C1", assembly, "2022", "Elsewhere", "E_coli", "ST1", "context"))
    return result


def write_alignment(path: Path, *, include_third_local: bool = False) -> Path:
    records = [">L1\nACGTNN\n", ">L2\nATGTAA\n"]
    if include_third_local:
        records.append(">L3\nATGCAA\n")
    records.append(">C1\nGTTTAA\n")
    path.write_text("".join(records), encoding="utf-8")
    return path


def temporal_assessment() -> dict[str, object]:
    return {
        "code": "not_supported",
        "supported": False,
        "reason": "Sampling dates did not outperform randomised dates.",
    }


def test_pairwise_clonal_snps_use_pair_specific_callable_sites(tmp_path: Path) -> None:
    alignment = write_alignment(tmp_path / "filtered.fasta")

    records = pairwise_distances(read_alignment(alignment), samples(tmp_path))
    local_pair = next(record for record in records if record["comparison"] == "focal_focal")

    assert local_pair["clonal_snps"] == 1
    assert local_pair["callable_sites"] == 4
    assert local_pair["period_comparison"] == "between_years"
    assert local_pair["patient_comparison"] == "same_patient"


def test_unknown_months_are_compared_at_year_precision(tmp_path: Path) -> None:
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig\nACGT\n", encoding="utf-8")
    uncertain = [
        Sample("A", assembly, "2024-XX-XX", "KIMS", "E_coli", "ST1", "local"),
        Sample("B", assembly, "2024-XX-XX", "KIMS", "E_coli", "ST1", "local"),
    ]

    records = pairwise_distances({"A": "ACGT", "B": "ATGT"}, uncertain)

    assert records[0]["period_comparison"] == "same_year"
