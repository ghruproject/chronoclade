from pathlib import Path
from random import Random

import pytest

from chronoclade.metadata import Sample
from chronoclade.workflow import native_platform_supported, run_workflow, tool_status


def write_variant_assembly(path: Path, mutations: dict[int, str]) -> None:
    rng = Random(42)
    sequence = [rng.choice("ACGT") for _ in range(3000)]
    for position, base in mutations.items():
        sequence[position] = base
    path.write_text(">chromosome\n" + "".join(sequence) + "\n", encoding="utf-8")


@pytest.mark.integration
def test_small_workflow_reaches_time_scaled_and_location_trees(tmp_path: Path) -> None:
    if not native_platform_supported():
        pytest.skip("the compiled workflow is validated natively on Linux")
    if any(path is None for path in tool_status().values()):
        pytest.skip("bioinformatics tools are not installed")

    mutation_sets = [
        {100: "C", 500: "G", 900: "T"},
        {100: "C", 500: "G", 1300: "T"},
        {100: "C", 1700: "G", 2100: "T"},
        {2500: "C", 2600: "G", 2700: "T"},
    ]
    samples = []
    for index, mutations in enumerate(mutation_sets):
        assembly = tmp_path / f"sample_{index}.fasta"
        write_variant_assembly(assembly, mutations)
        samples.append(
            Sample(
                sample_id=f"S{index}",
                assembly=assembly,
                collection_date=str(2020 + index),
                location="KIMS" if index < 3 else "Context",
                species="E_coli",
                lineage="ST_test",
                origin="local" if index < 3 else "context",
                is_reference=index == 0,
            )
        )

    output = tmp_path / "results"
    summary = run_workflow(
        samples,
        output=output,
        threads=1,
        lineage_jobs=1,
        randomisation_jobs=1,
        randomisations=1,
        # This plumbing test deliberately uses a permissive gate so CI executes
        # TreeTime's confidence-aware dated-tree branch with a tiny dataset.
        temporal_p_value=1.0,
        min_samples=4,
        seed=1,
        force=False,
    )

    assert len(summary["lineages"]) == 1
    assert (output / "E_coli__ST_test" / "clonalframeml.labelled_tree.newick").is_file()
    assert (output / "E_coli__ST_test" / "timetree" / "timetree.nexus").is_file()
    assert (output / "E_coli__ST_test" / "timetree" / "timetree.svg").is_file()
    assert (output / "E_coli__ST_test" / "timetree.png").is_file()
    assert (output / "E_coli__ST_test" / "timetree_with_confidence.svg").is_file()
    assert (output / "E_coli__ST_test" / "location" / "annotated_tree.nexus").is_file()
    assert (output / "E_coli__ST_test" / "clonal_pairwise_distances.tsv").is_file()
    assert (output / "E_coli__ST_test" / "clonal_snp_matrix.tsv").is_file()
    assert (output / "E_coli__ST_test" / "pairwise_callable_sites.tsv").is_file()
    assert (output / "E_coli__ST_test" / "clonal_snp_heatmap.svg").is_file()
    assert (output / "E_coli__ST_test" / "clonal_snp_heatmap.png").is_file()
    assert (output / "E_coli__ST_test" / "root_to_tip.png").is_file()
    assert (output / "E_coli__ST_test" / "date_randomisation.csv").is_file()
    assert (output / "E_coli__ST_test" / "date_randomisation.png").is_file()
    assert (output / "E_coli__ST_test" / "timetree_confidence.csv").is_file()
    assert (output / "E_coli__ST_test" / "node_dates.csv").is_file()
    assert (output / "E_coli__ST_test" / "supporting_results.zip").is_file()
    assert (output / "E_coli__ST_test" / "public_health_evidence.json").is_file()
    scenario = summary["lineages"][0]["public_health"]["scenario"]
    assert scenario["code"] in {
        "persistent_local_lineage",
        "multiple_introductions",
        "mixed",
        "indeterminate",
    }
