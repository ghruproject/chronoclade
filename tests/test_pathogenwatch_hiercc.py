import csv
import gzip
import hashlib
import io
from pathlib import Path

import pytest

from chronoclade.pathogenwatch_hiercc import HierCCExportError, parse_hiercc_export

FIXTURE = Path(__file__).parent / "fixtures" / "pathogenwatch-hclink-562-4-id108048.csv.gz"
IDENT = "iLGx7naVnHXTvwPF8DEgAo"


def csv_gzip(headers, rows):
    text = io.StringIO(newline="")
    writer = csv.writer(text)
    writer.writerow(headers)
    writer.writerows(rows)
    return gzip.compress(text.getvalue().encode())


def fixture_result():
    raw = FIXTURE.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == "c5ccab9dabf1f1ab4a831ea385c2819bac8db9ba5d187ab63556c9a29a5fdd76"
    return raw, parse_hiercc_export(raw, [IDENT], job="hclink-562-4")[0]


def test_real_public_single_genome_response_preserves_partial_levels_and_unknown_scope():
    raw, row = fixture_result()
    assert row["source_genome_id"] == IDENT
    assert row["hiercc_codes"] == {
        "HC50": "7", "HC100": "7", "HC200": "7", "HC400": "7",
        "HC1100": "7", "HC1500": "7", "HC2000": "7", "HC2350": "1",
    }
    assert row["hiercc_blank_levels"] == ["HC0", "HC2", "HC5", "HC10", "HC20"]
    assert all(row["hiercc_raw_levels"][level] == "" for level in row["hiercc_blank_levels"])
    assert row["hiercc_closest_st"] == "246645"
    assert row["hiercc_distance_to_closest_st"] == 35.7
    assert row["hiercc_distance_to_closest_st_raw"] == "35.7"
    assert row["hiercc_export_sha256"] == hashlib.sha256(raw).hexdigest()
    assert row["hiercc_provider_job"] == "hclink-562-4"
    assert row["hiercc_scheme_version"] is None
    assert row["hiercc_database_version"] is None
    assert row["hiercc_database_sha256"] is None
    assert row["hiercc_status"] == "resolved"


def test_explicit_zero_cluster_is_assigned_but_blank_remains_missing():
    raw = csv_gzip(
        ["Genome ID", "Genome Name", "Closest ST", "Distance to closest ST", "HC0", "HC2"],
        [["g1", "sample", "", "", "0", ""]],
    )
    row = parse_hiercc_export(raw, ["g1"], job="hclink-562-4")[0]
    assert row["hiercc_codes"] == {"HC0": "0"}
    assert row["hiercc_blank_levels"] == ["HC2"]
    assert row["hiercc_status"] == "resolved"


def test_all_blank_levels_are_retained_as_unassigned_not_zero():
    raw = csv_gzip(
        ["Genome ID", "Genome Name", "Closest ST", "Distance to closest ST", "HC0"],
        [["g1", "sample", "", "", ""]],
    )
    row = parse_hiercc_export(raw, ["g1"], job="hclink-562-4")[0]
    assert row["hiercc_codes"] == {}
    assert row["hiercc_blank_levels"] == ["HC0"]
    assert row["hiercc_status"] == "unassigned"


@pytest.mark.parametrize(
    "raw,requested,job",
    [
        (b"not gzip or CSV", [IDENT], "hclink-562-4"),
        (csv_gzip(["Genome ID", "Genome Name", "Closest ST", "Distance to closest ST", "HC0"], [["other", "x", "", "", "1"]]), [IDENT], "hclink-562-4"),
        (csv_gzip(["Genome ID", "Genome Name", "Closest ST", "Distance to closest ST", "HC0"], [[IDENT, "x", "", "", "1"], [IDENT, "x", "", "", "1"]]), [IDENT], "hclink-562-4"),
        (csv_gzip(["Genome ID", "Genome Name", "Closest ST", "Distance to closest ST", "HC0"], [[IDENT, "x", "", "", "1"]]), [IDENT, "missing"], "hclink-562-4"),
        (csv_gzip(["Genome ID", "Genome Name", "Closest ST", "Distance to closest ST", "HC0"], [[IDENT, "x", "", "", "1"]]), [IDENT], "bad/job"),
    ],
)
def test_identity_and_input_validation_rejects_unusable_export(raw, requested, job):
    with pytest.raises(HierCCExportError):
        parse_hiercc_export(raw, requested, job=job)


@pytest.mark.parametrize("value", ["-1", "1.2", "NaN", "abc"])
def test_invalid_cluster_values_are_rejected(value):
    raw = csv_gzip(
        ["Genome ID", "Genome Name", "Closest ST", "Distance to closest ST", "HC0"],
        [["g1", "sample", "", "", value]],
    )
    with pytest.raises(HierCCExportError, match="cluster IDs"):
        parse_hiercc_export(raw, ["g1"], job="hclink-562-4")


@pytest.mark.parametrize("distance", ["-0.1", "NaN", "Infinity", "1e999", "nope"])
def test_invalid_closest_distances_are_rejected(distance):
    raw = csv_gzip(
        ["Genome ID", "Genome Name", "Closest ST", "Distance to closest ST", "HC0"],
        [["g1", "sample", "ST1", distance, "1"]],
    )
    with pytest.raises(HierCCExportError, match="distance"):
        parse_hiercc_export(raw, ["g1"], job="hclink-562-4")


def test_duplicate_normalized_level_headers_rejected():
    raw = csv_gzip(
        ["Genome ID", "Genome Name", "Closest ST", "Distance to closest ST", "HC0", "HC00"],
        [["g1", "sample", "", "", "1", "1"]],
    )
    with pytest.raises(HierCCExportError, match="normalized HierCC level"):
        parse_hiercc_export(raw, ["g1"], job="hclink-562-4")


def test_ragged_and_malformed_csv_are_rejected():
    headers = ["Genome ID", "Genome Name", "Closest ST", "Distance to closest ST", "HC0"]
    with pytest.raises(HierCCExportError, match="extra unnamed"):
        parse_hiercc_export(
            b"Genome ID,Genome Name,Closest ST,Distance to closest ST,HC0\ng1,x,,,1,extra\n",
            ["g1"], job="hclink-562-4",
        )
    with pytest.raises(HierCCExportError, match="malformed CSV"):
        parse_hiercc_export(
            b'Genome ID,Genome Name,Closest ST,Distance to closest ST,HC0\ng1,"unterminated,,,1\n',
            ["g1"], job="hclink-562-4",
        )
    with pytest.raises(HierCCExportError, match="required column values"):
        parse_hiercc_export(
            csv_gzip(headers, [["g1", "sample"]]), ["g1"], job="hclink-562-4"
        )
    with pytest.raises(HierCCExportError, match="malformed CSV header"):
        parse_hiercc_export(
            b'"unterminated\nGenome ID,Genome Name,Closest ST,Distance to closest ST,HC0\ng1,x,,,1\n',
            ["g1"], job="hclink-562-4",
        )
