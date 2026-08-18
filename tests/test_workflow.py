from pathlib import Path

from beyondmlst.metadata import Sample
from beyondmlst.workflow import plan


def test_plan_separates_species_and_lineages(tmp_path: Path) -> None:
    assembly = tmp_path / "sample.fasta"
    assembly.write_text(">contig\nAAAA\n", encoding="utf-8")
    samples = [
        Sample("E1", assembly, "2020", "KIMS", "E_coli", "ST131", "local"),
        Sample("E2", assembly, "2021", "KIMS", "E_coli", "ST131", "local"),
        Sample("K1", assembly, "2022", "KIMS", "K_pneumoniae", "ST15", "local"),
    ]

    items = plan(samples, min_samples=2)

    assert len(items) == 2
    assert items[0]["status"] == "ready"
    assert items[1]["status"] == "too_few_samples"
