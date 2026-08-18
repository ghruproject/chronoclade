from pathlib import Path

import pytest

from beyondmlst.temporal import TemporalError, parse_clock


def test_parse_treetime_clock(tmp_path: Path) -> None:
    output = tmp_path / "molecular_clock.txt"
    output.write_text(
        "Root-Tip-Regression:\n --rate:\t1.234e-06\n --r^2:  \t0.67\n",
        encoding="utf-8",
    )

    assert parse_clock(output) == {"rate": 1.234e-06, "r_squared": 0.67}


def test_rejects_unrecognised_clock_output(tmp_path: Path) -> None:
    output = tmp_path / "molecular_clock.txt"
    output.write_text("not a clock\n", encoding="utf-8")

    with pytest.raises(TemporalError):
        parse_clock(output)
