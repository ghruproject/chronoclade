from pathlib import Path

from beyondmlst.evidence import (
    build_public_health_evidence,
    pairwise_distances,
    patient_sensitivity,
    read_alignment,
)
from beyondmlst.metadata import Sample


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


def test_public_health_evidence_writes_matrices_heatmap_and_persistent_scenario(
    tmp_path: Path,
) -> None:
    alignment = write_alignment(tmp_path / "filtered.fasta")
    tree = tmp_path / "tree.newick"
    tree.write_text("((L1:0.1,L2:0.1)NODE:0.1,C1:0.2)ROOT;\n", encoding="utf-8")

    result = build_public_health_evidence(
        filtered_alignment=alignment,
        rooted_tree=tree,
        samples=samples(tmp_path),
        output=tmp_path,
        temporal_assessment=temporal_assessment(),
    )

    assert result["scenario"]["code"] == "persistent_local_lineage"
    assert result["scenario"]["confidence"] == "low"
    assert "One focal-only topology group" in result["scenario"]["decision_rule"]
    assert result["topology"]["candidate_local_group_count"] == 1
    assert result["patient_sensitivity"]["excluded_repeated_patient_samples"] == ["L2"]
    assert (tmp_path / "clonal_pairwise_distances.tsv").is_file()
    assert (tmp_path / "clonal_snp_matrix.tsv").is_file()
    assert (tmp_path / "pairwise_callable_sites.tsv").is_file()
    assert (tmp_path / "clonal_snp_heatmap.svg").is_file()
    assert (tmp_path / "public_health_evidence.json").is_file()


def test_separated_local_group_plus_longitudinal_group_is_mixed(tmp_path: Path) -> None:
    alignment = write_alignment(tmp_path / "filtered.fasta", include_third_local=True)
    tree = tmp_path / "tree.newick"
    tree.write_text("((L1,L2),C1,L3);\n", encoding="utf-8")

    result = build_public_health_evidence(
        filtered_alignment=alignment,
        rooted_tree=tree,
        samples=samples(tmp_path, include_third_local=True),
        output=tmp_path,
        temporal_assessment=temporal_assessment(),
    )

    assert result["scenario"]["code"] == "mixed"
    assert "At least two separated" in result["scenario"]["decision_rule"]
    assert result["topology"]["candidate_local_group_count"] == 2
    assert result["distance_summary"]["categories"]["within_candidate_groups"]["comparisons"] == 1
    assert result["distance_summary"]["categories"]["between_candidate_groups"]["comparisons"] == 2


def test_patient_sensitivity_retains_earliest_isolate_per_patient(tmp_path: Path) -> None:
    lineage_samples = samples(tmp_path)
    records = pairwise_distances(
        read_alignment(write_alignment(tmp_path / "filtered.fasta")), lineage_samples
    )

    sensitivity = patient_sensitivity(lineage_samples, records)

    assert sensitivity["representative_samples"] == ["L1"]
    assert sensitivity["excluded_repeated_patient_samples"] == ["L2"]
