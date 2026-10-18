import csv
import gzip
import hashlib
import io
from pathlib import Path

import pytest

from chronoclade.pathogenwatch import PathogenwatchClient
from chronoclade.profile_inputs import (
    ProfileInputError,
    _canonical,
    _frozen_cglin,
    collection_id,
    materialise_assemblies,
    parse_cgmlst_export,
)

UUID = "jX5cwsoyJ1KDquMqssUzAD"


def test_canonical_adapter_preserves_public_assembly_identity_metadata():
    checksum = "a" * 40
    row = _canonical(
        {
            "source_genome_id": "public-id",
            "id": 42,
            "mlst_st": "147",
            "checksum": checksum,
            "length": 1234,
            "qc": True,
        },
        origin="context",
        species="Klebsiella pneumoniae",
    )
    assert row["numeric_source_id"] == 42
    assert row["source_checksum"] == checksum
    assert row["source_length"] == 1234
    assert row["qc_pass"] is True


def test_canonical_adapter_preserves_explicit_city_without_inference_from_nuts2():
    row = _canonical(
        {
            "source_genome_id": "public-id",
            "metadata/User defined metadata/City": "Oxford",
            "metadata/User defined metadata/Region": "Oxfordshire",
            "metadata/User defined metadata/NUTS2": "UKJ1",
        },
        origin="context",
        species="Klebsiella pneumoniae",
    )
    assert row["city"] == "Oxford"
    assert row["region"] == "Oxfordshire"
    assert row["nuts2"] == "UKJ1"

    nuts_only = _canonical(
        {"source_genome_id": "other-id", "nuts2": "UKJ1"},
        origin="context",
        species="Klebsiella pneumoniae",
    )
    assert nuts_only["city"] == ""


@pytest.mark.parametrize(
    "text",
    [
        UUID,
        f"https://pathogen.watch/collections/{UUID}-readable-slug",
        f"https://next.pathogen.watch/collections/{UUID}",
    ],
)
def test_collection_identifier(text):
    assert collection_id(text) == UUID


@pytest.mark.parametrize(
    "text",
    [
        "https://evil.example/collections/" + UUID,
        "https://pathogen.watch/genomes/" + UUID,
        "https://user:pass@pathogen.watch/collections/" + UUID,
        "bad",
    ],
)
def test_collection_identifier_rejects_untrusted_url(text):
    with pytest.raises(ProfileInputError):
        collection_id(text)


def export(rows):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=["Genome ID", "Gene", "Allele ID"])
    writer.writeheader()
    writer.writerows(rows)
    return gzip.compress(stream.getvalue().encode())


def test_verified_export_shape_exact_ids_and_ambiguous_loci():
    raw = export(
        [
            {"Genome ID": "A", "Gene": "g1", "Allele ID": "2"},
            {"Genome ID": "A", "Gene": "g2", "Allele ID": "3"},
            {"Genome ID": "A", "Gene": "g2", "Allele ID": "4"},
            {"Genome ID": "A", "Gene": "g3", "Allele ID": "0"},
        ]
    )
    row = parse_cgmlst_export(raw, ["A"], job="cgmlst-1280-2")[0]
    assert row["cgmlst_profile"] == {"g1": "2"}
    assert row["cgmlst_ambiguous_loci"] == ["g2"]
    assert row["cgmlst_loci"] == ["g1", "g2", "g3"]
    assert row["cgmlst_locus_universe_complete"] is False
    assert "cgmlst_called_fraction" not in row
    with pytest.raises(ProfileInputError, match="unrequested"):
        parse_cgmlst_export(raw, ["B"], job="cgmlst-1280-2")


def test_deferred_assembly_acquisition_only_selected_rows(tmp_path, monkeypatch):
    seen = []

    def download(rows, **kwargs):
        seen.extend(row["source_genome_id"] for row in rows)
        path = tmp_path / "a.fa"
        path.write_text(">a\nACGT\n")
        return {"A": path}, [{"status": "downloaded"}]

    monkeypatch.setattr("chronoclade.profile_inputs.download_assemblies", download)
    client = PathogenwatchClient(api_key="fixture")
    records = [
        {
            "sample_id": "PW_A",
            "source_genome_id": "A",
            "origin": "context",
            "species": "E. coli",
            "lineage": "ST1",
        }
    ]
    result = materialise_assemblies(records, output=tmp_path, client=client)
    assert seen == ["A"]
    assert result[0].origin == "context"
    assert "assembly" not in records[0]


def test_explicit_frozen_cglin_unmatched_ids_cannot_keep_live_namespace():
    from chronoclade.cglin import normalise_assignment

    live = normalise_assignment(
        {
            "source_genome_id": "unmatched",
            "cglin_raw": "1_2_3_4_5_6_7_8_9_10",
            "cglin_scheme_version": "live-version",
            "cgst": "42",
        }
    )
    live.update(cgmlst_profile={"g1": "1"}, cglin_database_sha256="live-fingerprint")
    frozen = normalise_assignment(
        {
            "source_genome_id": "another",
            "cglin_raw": "1_2_3_4_5_6_7_8_9_10",
            "cglin_scheme_version": "unknown",
        }
    )
    row = _frozen_cglin([live], [frozen])[0]
    assert row["cglin_status"] == "missing"
    assert row["cglin_group_5"] == ""
    assert row["cglin_scheme_version"] == "unknown"
    assert row["cglin_scope_source"] == "user_frozen_export"
    assert row["cglin_export_record_count"] == 0
    assert row["cgst"] == ""
    assert "cglin_database_sha256" not in row
    assert row["cgmlst_profile"] == {"g1": "1"}


def test_cgmlst_identity_mismatch_retries_exact_batch_without_accepting_foreign_calls(
    tmp_path, monkeypatch
):
    from chronoclade.profile_inputs import _analysis_exports

    wrong = export([{"Genome ID": "foreign", "Gene": "g1", "Allele ID": "9"}])
    correct = export([{"Genome ID": "A", "Gene": "g1", "Allele ID": "1"}])
    requests, waits = [], []

    def fetch(route, **kwargs):
        requests.append((route, kwargs["body"]))
        return wrong if len(requests) == 1 else correct

    monkeypatch.setattr("chronoclade.profile_inputs.request_download", fetch)
    monkeypatch.setattr("chronoclade.profile_inputs.time.sleep", waits.append)
    provenance = {}
    rows = _analysis_exports(
        [{"source_genome_id": "A", "numeric_source_id": 17}],
        [{"name": "cgmlst", "job": "cgmlst-573-2"}],
        PathogenwatchClient(api_key="fixture"),
        tmp_path,
        provenance,
    )
    assert requests == [
        ("/api/downloads/cgmlst?job=cgmlst-573-2", {"ids": "17"}),
        ("/api/downloads/cgmlst?job=cgmlst-573-2", {"ids": "17"}),
    ]
    rejected = provenance["analysis_export_identity_retries"]
    assert len(rejected) == 1
    assert rejected[0]["requested_source_ids"] == ["A"]
    assert rejected[0]["rejected_export_sha256"] == hashlib.sha256(wrong).hexdigest()
    assert Path(rejected[0]["rejected_export"]).read_bytes() == wrong
    assert rows[0]["cgmlst_profile"] == {"g1": "1"}
    assert waits == [1]
    assert provenance["analysis_exports"][-1]["status"] == "downloaded"


def test_cgmlst_identity_mismatch_exhaustion_keeps_profile_missing(tmp_path, monkeypatch):
    from chronoclade.profile_inputs import _analysis_exports

    wrong = export([{"Genome ID": "foreign", "Gene": "g1", "Allele ID": "9"}])
    requests, waits = [], []

    def fetch(route, **kwargs):
        requests.append(route)
        return wrong

    monkeypatch.setattr("chronoclade.profile_inputs.request_download", fetch)
    monkeypatch.setattr("chronoclade.profile_inputs.time.sleep", waits.append)
    provenance = {}
    rows = _analysis_exports(
        [{"source_genome_id": "A", "numeric_source_id": 17}],
        [{"name": "cgmlst", "job": "cgmlst-573-2"}],
        PathogenwatchClient(api_key="fixture"),
        tmp_path,
        provenance,
    )
    assert len(requests) == 3
    assert len(provenance["analysis_export_identity_retries"]) == 3
    assert not rows[0].get("cgmlst_profile")
    assert waits == [1, 2]
    assert provenance["analysis_exports"][-1]["status"] == "unavailable"


def test_cgmlst_malformed_locus_is_not_retried_as_transient_identity(tmp_path, monkeypatch):
    from chronoclade.profile_inputs import _analysis_exports

    requests = []

    def fetch(*args, **kwargs):
        requests.append(kwargs["body"])
        return export([{"Genome ID": "A", "Gene": "", "Allele ID": "1"}])

    monkeypatch.setattr("chronoclade.profile_inputs.request_download", fetch)
    provenance = {}
    result = _analysis_exports(
        [{"source_genome_id": "A", "numeric_source_id": 17}],
        [{"name": "cgmlst", "job": "cgmlst-573-2"}],
        PathogenwatchClient(api_key="fixture"),
        tmp_path,
        provenance,
    )
    assert requests == [{"ids": "17"}]
    assert not result[0].get("cgmlst_profile")
    assert "analysis_export_identity_retries" not in provenance
    assert provenance["analysis_exports"][-1]["reason"] == "cgMLST export contains an empty locus"
