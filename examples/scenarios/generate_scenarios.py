#!/usr/bin/env python3
"""Generate controlled reports for all four public-health interpretations.

These small alignments and supplied trees are synthetic test fixtures. They are
designed to exercise report behaviour, not to estimate biological performance.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

from beyondmlst.evidence import build_public_health_evidence
from beyondmlst.metadata import Sample
from beyondmlst.report import assess_temporal_signal, write_lineage_report, write_summary_report

ALIGNMENT_LENGTH = 2_000


@dataclass(frozen=True)
class Scenario:
    slug: str
    lineage: str
    expected_code: str
    tree: str
    samples: tuple[tuple[str, str, str, str, str], ...]
    mutations: dict[str, set[int]]
    temporal_signal: bool = False


def _sequence(positions: set[int]) -> str:
    bases = ["A"] * ALIGNMENT_LENGTH
    for position in positions:
        bases[position] = "C"
    return "".join(bases)


def _temporal_result(supported: bool) -> dict[str, object]:
    randomised = [
        {"rate": 2e-7 + index * 1e-8, "r_squared": 0.04 + index * 0.01}
        for index in range(20)
    ]
    return {
        "observed": {
            "rate": 1.1e-6 if supported else 3.5e-7,
            "r_squared": 0.72 if supported else 0.18,
        },
        "requested_randomisations": 20,
        "successful_randomisations": 20,
        "p_value_r_squared": 0.047 if supported else 0.38,
        "randomised": randomised,
    }


def _scenarios() -> tuple[Scenario, ...]:
    contexts = (
        ("C1", "2021-04", "context", "Regional_context", ""),
        ("C2", "2022-10", "context", "International_context", ""),
        ("C3", "2024-02", "context", "Regional_context", ""),
    )
    context_mutations = {
        "C1": set(range(100, 116)),
        "C2": set(range(200, 219)),
        "C3": set(range(300, 322)),
    }
    return (
        Scenario(
            slug="persistent_local_lineage",
            lineage="SYNTH_PERSISTENT",
            expected_code="persistent_local_lineage",
            tree="(((F1:0.1,F2:0.1):0.1,F3:0.1):0.2,(C1:0.2,C2:0.2,C3:0.2):0.2)ROOT;",
            samples=(
                ("F1", "2022-01", "local", "KIMS", "P001"),
                ("F2", "2023-05", "local", "KIMS", "P002"),
                ("F3", "2024-08", "local", "KIMS", "P003"),
                *contexts,
            ),
            mutations={
                "F1": {1, 2, 3},
                "F2": {1, 2, 3, 4},
                "F3": {1, 2, 3, 4, 5},
                **context_mutations,
            },
            temporal_signal=True,
        ),
        Scenario(
            slug="multiple_introductions",
            lineage="SYNTH_MULTIPLE",
            expected_code="multiple_introductions",
            tree="((F1:0.1,C1:0.1):0.2,(F2:0.1,C2:0.1):0.2,(F3:0.1,C3:0.1):0.2)ROOT;",
            samples=(
                ("F1", "2022-01", "local", "KIMS", "P011"),
                ("F2", "2023-05", "local", "KIMS", "P012"),
                ("F3", "2024-08", "local", "KIMS", "P013"),
                *contexts,
            ),
            mutations={
                "F1": set(range(100, 115)),
                "F2": set(range(200, 218)),
                "F3": set(range(300, 321)),
                **context_mutations,
            },
        ),
        Scenario(
            slug="mixed",
            lineage="SYNTH_MIXED",
            expected_code="mixed",
            tree="(((F1:0.1,F2:0.1):0.1,C1:0.2):0.2,(F3:0.1,C2:0.1):0.2,(F4:0.1,C3:0.1):0.2)ROOT;",
            samples=(
                ("F1", "2022-01", "local", "KIMS", "P021"),
                ("F2", "2024-03", "local", "KIMS", "P022"),
                ("F3", "2023-05", "local", "KIMS", "P023"),
                ("F4", "2024-08", "local", "KIMS", "P024"),
                *contexts,
            ),
            mutations={
                "F1": {1, 2, 3},
                "F2": {1, 2, 3, 4},
                "F3": set(range(200, 218)),
                "F4": set(range(300, 321)),
                **context_mutations,
            },
        ),
        Scenario(
            slug="indeterminate",
            lineage="SYNTH_INDETERMINATE",
            expected_code="indeterminate",
            tree="((F1:0.1,C1:0.1):0.2,(C2:0.1,C3:0.1):0.2)ROOT;",
            samples=(
                ("F1", "2024-03", "local", "KIMS", "P031"),
                *contexts,
            ),
            mutations={"F1": {1, 2, 3}, **context_mutations},
        ),
    )


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text(
        "".join(f">{sample_id}\n{sequence}\n" for sample_id, sequence in records),
        encoding="utf-8",
    )


def _build_scenario(scenario: Scenario, output: Path) -> dict[str, object]:
    directory = output / scenario.slug
    assemblies = directory / "assemblies"
    assemblies.mkdir(parents=True, exist_ok=True)
    samples: list[Sample] = []
    alignment_records: list[tuple[str, str]] = []
    for sample_id, date, origin, location, patient_id in scenario.samples:
        sequence = _sequence(scenario.mutations[sample_id])
        assembly = assemblies / f"{sample_id}.fasta"
        _write_fasta(assembly, [(sample_id, sequence)])
        alignment_records.append((sample_id, sequence))
        samples.append(
            Sample(
                sample_id,
                assembly,
                date,
                location,
                "Escherichia_coli",
                scenario.lineage,
                origin,
                patient_id=patient_id,
            )
        )

    alignment = directory / "clonalframeml.filtered.fasta"
    tree = directory / "clonalframeml.labelled_tree.newick"
    _write_fasta(alignment, alignment_records)
    tree.write_text(scenario.tree + "\n", encoding="utf-8")
    temporal = _temporal_result(scenario.temporal_signal)
    assessment = assess_temporal_signal(temporal, p_value_threshold=0.05)
    public_health = build_public_health_evidence(
        filtered_alignment=alignment,
        rooted_tree=tree,
        samples=samples,
        output=directory,
        temporal_assessment=assessment,
    )
    actual_code = str(public_health["scenario"]["code"])
    if actual_code != scenario.expected_code:
        raise RuntimeError(
            f"{scenario.slug}: expected {scenario.expected_code}, observed {actual_code}"
        )

    focal = [sample for sample in samples if sample.origin != "context"]
    contexts = [sample for sample in samples if sample.origin == "context"]
    report = {
        "species": "Escherichia_coli",
        "lineage": scenario.lineage,
        "sample_count": len(samples),
        "distinct_dates": len({sample.collection_date for sample in focal}),
        "temporal_signal": temporal,
        "context": {
            "local_samples": len(focal),
            "context_samples": len(contexts),
            "context_locations": sorted({sample.location for sample in contexts}),
            "manifest_available": False,
            "interpretation": "Synthetic contextual genomes were supplied for report testing.",
        },
        "public_health": public_health,
    }
    (directory / "temporal_signal.json").write_text(
        json.dumps(temporal, indent=2) + "\n", encoding="utf-8"
    )
    (directory / "report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    write_lineage_report(report, directory=directory, p_value_threshold=0.05)
    return {
        "species": report["species"],
        "lineage": report["lineage"],
        "slug": scenario.slug,
        "sample_count": report["sample_count"],
        "temporal_status": assessment["code"],
        "temporal_signal": temporal,
        "context": report["context"],
        "public_health": public_health,
    }


def generate(output: Path) -> Path:
    output.mkdir(parents=True, exist_ok=True)
    lineages = [_build_scenario(scenario, output) for scenario in _scenarios()]
    (output / "README.txt").write_text(
        "Synthetic fixtures for interface and scenario-logic testing only.\n"
        "They are not biological validation datasets and must not be presented as such.\n",
        encoding="utf-8",
    )
    write_summary_report({"lineages": lineages}, output=output)
    return output / "index.html"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("scenario_reports"),
        help="Directory for the four generated scenario reports",
    )
    args = parser.parse_args()
    print(f"Generated synthetic scenario reports: {generate(args.output.resolve())}")


if __name__ == "__main__":
    main()
