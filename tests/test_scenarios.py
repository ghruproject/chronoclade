import json
from pathlib import Path

from examples.scenarios.generate_scenarios import generate


def test_controlled_scenarios_generate_all_four_expected_interpretations(
    tmp_path: Path,
) -> None:
    index = generate(tmp_path)

    expected = {
        "persistent_local_lineage": "persistent_local_lineage",
        "multiple_introductions": "multiple_introductions",
        "mixed": "mixed",
        "indeterminate": "indeterminate",
    }
    observed = {}
    for directory, code in expected.items():
        report = json.loads((tmp_path / directory / "report.json").read_text(encoding="utf-8"))
        observed[directory] = report["public_health"]["scenario"]["code"]
        html = (tmp_path / directory / "report.html").read_text(encoding="utf-8")
        assert 'id="root-to-tip"' in html
        assert '<section class="stage" id="randomisation"' in html
        assert '<section class="stage" id="interpretation"' in html
        assert "What does the genomic evidence support?" in html
        assert code in html

    assert observed == expected
    assert index.is_file()
    assert "not biological validation datasets" in (tmp_path / "README.txt").read_text(
        encoding="utf-8"
    )
