import threading
import time
from pathlib import Path

import pytest

from beyondmlst import workflow
from beyondmlst.metadata import Sample
from beyondmlst.workflow import (
    allocate_resources,
    complete_alignment_sites,
    context_evidence,
    plan,
    read_context_manifest,
)


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


def test_resource_plan_respects_global_thread_budget() -> None:
    resources = allocate_resources(
        total_threads=16,
        lineage_jobs=3,
        randomisation_jobs=8,
        ready_lineages=5,
    )

    assert resources.lineage_workers == 3
    assert resources.threads_per_lineage == 5
    assert resources.randomisation_workers_per_lineage == 5
    assert resources.lineage_workers * resources.threads_per_lineage <= 16


def test_resource_plan_uses_all_threads_for_one_lineage() -> None:
    resources = allocate_resources(
        total_threads=8,
        lineage_jobs=2,
        randomisation_jobs=4,
        ready_lineages=1,
    )

    assert resources.lineage_workers == 1
    assert resources.threads_per_lineage == 8
    assert resources.randomisation_workers_per_lineage == 4


def test_complete_alignment_sites_excludes_ambiguous_columns(tmp_path: Path) -> None:
    alignment = tmp_path / "alignment.fasta"
    alignment.write_text(">A\nACGTNA\n>B\nACGT-A\n", encoding="utf-8")

    assert complete_alignment_sites(alignment) == 5


def test_context_manifest_is_summarised_for_the_lineage(tmp_path: Path) -> None:
    assembly = tmp_path / "sample.fasta"
    assembly.write_text(">contig\nAAAA\n", encoding="utf-8")
    members = [
        Sample("F1", assembly, "2024", "KIMS", "E_coli", "ST131", "local"),
        Sample("C1", assembly, "2023", "Philippines", "E_coli", "ST131", "context"),
    ]
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample_id\tspecies\tlineage\tcountry\tcollection_date\tnearest_focal\tmin_ska_distance\tselection_reason\n"
        "C1\tE_coli\tST131\tPhilippines\t2023\tF1\t3\tnearest_to=F1\n",
        encoding="utf-8",
    )

    evidence = context_evidence(
        members, read_context_manifest(manifest), directory=tmp_path / "lineage"
    )

    assert evidence["local_samples"] == 1
    assert evidence["context_samples"] == 1
    assert evidence["nearest_screening_contexts"][0]["min_ska_distance"] == 3.0
    assert (tmp_path / "lineage" / "context_manifest.tsv").is_file()


def test_workflow_runs_lineages_with_bounded_concurrency(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assembly = tmp_path / "sample.fasta"
    assembly.write_text(">contig\nAAAA\n", encoding="utf-8")
    samples = [
        Sample(
            f"S{lineage}{index}",
            assembly,
            str(2020 + index),
            "KIMS",
            "E_coli",
            f"ST{lineage}",
            "local",
        )
        for lineage in range(3)
        for index in range(2)
    ]
    active = 0
    maximum_active = 0
    assigned_threads: list[int] = []
    lock = threading.Lock()

    def fake_run_lineage(
        item: dict[str, object], _: list[Sample], **kwargs: object
    ) -> dict[str, object]:
        nonlocal active, maximum_active
        with lock:
            active += 1
            maximum_active = max(maximum_active, active)
            assigned_threads.append(int(kwargs["threads"]))
        time.sleep(0.02)
        with lock:
            active -= 1
        return item

    monkeypatch.setattr(workflow, "check_tools", lambda: None)
    monkeypatch.setattr(workflow, "_run_lineage", fake_run_lineage)
    summary = workflow.run_workflow(
        samples,
        output=tmp_path / "results",
        threads=4,
        lineage_jobs=2,
        randomisation_jobs=4,
        randomisations=10,
        temporal_p_value=0.05,
        min_samples=2,
        seed=1,
        force=False,
    )

    assert len(summary["lineages"]) == 3
    assert maximum_active == 2
    assert assigned_threads == [2, 2, 2]
