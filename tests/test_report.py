from pathlib import Path

from beyondmlst.report import (
    assess_temporal_signal,
    write_lineage_report,
    write_summary_report,
)


def temporal_result(*, rate: float = 1e-6, p_value: float = 0.01) -> dict[str, object]:
    return {
        "observed": {"rate": rate, "r_squared": 0.65},
        "requested_randomisations": 4,
        "successful_randomisations": 4,
        "p_value_r_squared": p_value,
        "randomised": [
            {"rate": 2e-7, "r_squared": 0.10},
            {"rate": 4e-7, "r_squared": 0.15},
            {"rate": 6e-7, "r_squared": 0.20},
            {"rate": 8e-7, "r_squared": 0.25},
        ],
    }


def test_temporal_assessment_reports_supported_and_unsupported_results() -> None:
    supported = assess_temporal_signal(temporal_result(), p_value_threshold=0.05)
    unsupported = assess_temporal_signal(temporal_result(p_value=0.25), p_value_threshold=0.05)

    assert supported["code"] == "supported"
    assert supported["supported"] is True
    assert unsupported["code"] == "not_supported"
    assert unsupported["supported"] is False


def test_lineage_report_contains_visuals_verdict_and_guardrail(tmp_path: Path) -> None:
    (tmp_path / "clock").mkdir()
    (tmp_path / "timetree").mkdir()
    (tmp_path / "clock" / "root_to_tip_regression.svg").write_text(
        "<svg xmlns='http://www.w3.org/2000/svg'/>", encoding="utf-8"
    )
    (tmp_path / "timetree" / "timetree.svg").write_text(
        "<svg xmlns='http://www.w3.org/2000/svg'/>", encoding="utf-8"
    )
    report = {
        "species": "E_coli",
        "lineage": "ST131",
        "sample_count": 42,
        "distinct_dates": 12,
        "temporal_signal": temporal_result(),
        "context": {
            "local_samples": 30,
            "context_samples": 12,
            "context_locations": ["Philippines", "United Kingdom"],
            "manifest_available": True,
            "interpretation": "Contextual genomes are available.",
            "nearest_screening_contexts": [
                {
                    "sample_id": "SAMN1",
                    "country": "Philippines",
                    "collection_date": "2022",
                    "nearest_focal": "F1",
                    "min_ska_distance": 4.0,
                }
            ],
        },
    }

    output = write_lineage_report(report, directory=tmp_path, p_value_threshold=0.05)
    text = output.read_text(encoding="utf-8")

    assert "Temporal signal supported" in text
    assert "clock/root_to_tip_regression.svg" in text
    assert "date_randomisation.svg" in text
    assert "Time-scaled phylogeny" in text
    assert "Public contextual genomes" in text
    assert "SAMN1" in text
    assert "SKA distances" in text
    assert "does not prove direct transmission" in text
    assert (tmp_path / "date_randomisation.svg").is_file()


def test_unsupported_report_omits_dated_tree_visual(tmp_path: Path) -> None:
    report = {
        "species": "Klebsiella_pneumoniae",
        "lineage": "ST15",
        "sample_count": 20,
        "distinct_dates": 5,
        "temporal_signal": temporal_result(p_value=0.25),
    }

    output = write_lineage_report(report, directory=tmp_path, p_value_threshold=0.05)
    text = output.read_text(encoding="utf-8")

    assert "Temporal signal not supported" in text
    assert "Time-scaled phylogeny" not in text
    assert "cannot distinguish" in text


def test_summary_report_links_lineage_reports(tmp_path: Path) -> None:
    metrics = temporal_result()
    summary = {
        "lineages": [
            {
                "species": "E_coli",
                "lineage": "ST131",
                "slug": "E_coli__ST131",
                "sample_count": 42,
                "temporal_status": "supported",
                "temporal_signal": metrics,
            }
        ]
    }

    output = write_summary_report(summary, output=tmp_path)
    text = output.read_text(encoding="utf-8")

    assert "E_coli__ST131/report.html" in text
    assert "Supported" in text
