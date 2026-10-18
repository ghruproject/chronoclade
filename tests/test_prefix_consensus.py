"""Prefix-specific resolution preserves actual returned native assignment evidence."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest
import chronoclade
from chronoclade.cglin import CGLINError, normalise_assignment, resolved_prefix
from chronoclade.context_partitions import _group_key, partition_records, profile_records
from chronoclade.datasets import LocusCatalogue, from_profile_records, load_dataset, write_dataset
from chronoclade.profile_context_provider import _matching
from chronoclade.query_typing import _assignment_result

SANDBOX = Path(__file__).resolve().parents[1]
assert Path(chronoclade.__file__).resolve().parent == SANDBOX / "chronoclade"
CONFIG = {"kind": "plincer", "scheme": "scgMLST629_S", "scheme_version": "frozen-test"}
PREFIX = (0, 0, 105, 6, 0, 0, 0, 5)


def raw_assignment():
    return json.loads((SANDBOX / "tests/fixtures/lost-plincer-assignment.json").read_text())


def row(ident, role, assignment):
    return dict(
        assignment,
        sample_id=ident,
        role=role,
        origin="local" if role == "input" else "context",
        species="Klebsiella pneumoniae",
        mlst_st="258",
        mlst_scheme="mlst",
        cgmlst_scheme="test",
        cgmlst_scheme_version="1",
        cgmlst_profile={"a": "1", "b": "1"},
        cgmlst_status="resolved",
        cglin_database_sha256="same-frozen-reference",
    )


def query():
    return row(
        "q",
        "input",
        _assignment_result({"cgST": "134", "LINcode": [str(n) for n in (*PREFIX, 1, 0)]}, CONFIG),
    )


def test_saved_zero_distance_neighbour_keeps_known_eight_and_unknown_suffix():
    raw = raw_assignment()
    before = deepcopy(raw)
    candidate = _assignment_result(raw, CONFIG)
    assert raw == before
    assert candidate["cglin_raw"] == "0_0_105_6_0_0_0_5_*_*"
    assert candidate["cglin_provisional"] is True
    assert candidate["cglin_status"] == "provisional"
    assert (
        candidate["cglin_prefix_consensus"]["raw_assignment_sha256"]
        == hashlib.sha256(
            json.dumps(raw, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
    )
    assert resolved_prefix(candidate, 8) == PREFIX
    assert resolved_prefix(candidate, 9) is None
    assert resolved_prefix(candidate, 10) is None
    assert (
        candidate["cglin_export_row"]["cglin_prefix_consensus"]
        == candidate["cglin_prefix_consensus"]
    )
    richer = normalise_assignment(candidate, depths=(5, 8, 9, 10))
    assert richer["cglin_status_8"] == "resolved"
    assert richer["cglin_status_9"] == "partial"
    context = row("c", "context", candidate)
    assert _matching(query(), context, "cglin", 8)
    assert not _matching(query(), context, "cglin", 9)
    blocks, audit = partition_records([query(), context], lin_level=5, context_lin_level=8)
    assert [r["sample_id"] for r in blocks[0]["records"]] == ["c", "q"]
    assert audit["outside_context_scope_ids"] == []


@pytest.mark.parametrize("code", ["0_0_105_6_0_0_0_5_1_0", "0_0_105_6_0_0_0_5_*_*"])
def test_external_provisional_without_consensus_is_never_mistaken_for_proven_prefix(code):
    candidate = normalise_assignment(
        {"LINcode": code, "provisional": True}, scheme_version="frozen-test"
    )
    assert candidate["cglin_status_5"] == "provisional"
    context = row("c", "context", candidate)
    assert resolved_prefix(context, 5) is None
    assert not _matching(query(), context, "cglin", 8)
    assert _group_key(context, 5, None) is None


def test_unknown_requested_prefix_refuses_without_copying_or_fallback():
    raw = {
        "cgST": "new",
        "LINcode": [str(n) for n in (*PREFIX, 1, 0)],
        "matches": [{"LINcode": ["0", "0", "105", "6", "0", "0", "99", "1", "0", "0"]}],
    }
    unresolved = row("q", "input", _assignment_result(raw, CONFIG))
    assert unresolved["cglin_raw"] == "0_0_105_6_0_0_*_*_*_*"
    assert resolved_prefix(unresolved, 5) == PREFIX[:5]
    assert resolved_prefix(unresolved, 8) is None
    assert not _matching(unresolved, row("c", "context", query()), "cglin", 8)
    blocks, audit = partition_records([unresolved], lin_level=5, context_lin_level=8)
    assert len(blocks) == 1 and len(blocks[0]["records"]) == 1
    assert audit["unresolved_context_scope_inputs"] == [
        {"sample_id": "q", "reason": "missing_or_unresolved_context_scope_typing"}
    ]


def test_tampered_consensus_and_true_different_known_prefix_refuse():
    candidate = _assignment_result(raw_assignment(), CONFIG)
    candidate["cglin_prefix_consensus"]["known_prefix"][-1] = 99
    assert resolved_prefix(candidate, 8) is None
    with pytest.raises(CGLINError, match="disagrees"):
        normalise_assignment(candidate)
    outside = row(
        "c",
        "context",
        _assignment_result(
            {"cgST": "known", "LINcode": [str(n) for n in (*PREFIX[:7], 50, 1, 0)]}, CONFIG
        ),
    )
    assert not _matching(query(), outside, "cglin", 8)
    blocks, audit = partition_records([query(), outside], lin_level=5, context_lin_level=8)
    assert audit["outside_context_scope_ids"] == ["c"]
    assert len(blocks[0]["records"]) == 1


def test_proof_and_raw_assignment_survive_prepared_bundle_roundtrip(tmp_path):
    candidate = row("c", "context", _assignment_result(raw_assignment(), CONFIG))
    dataset = from_profile_records(
        [query(), candidate], [LocusCatalogue("test", "1", ("a", "b"), "frozen")]
    )
    manifest = write_dataset(dataset, tmp_path / "prepared")
    restored = {r["sample_id"]: r for r in profile_records(load_dataset(manifest))}
    assert restored["c"]["cglin_prefix_consensus"] == candidate["cglin_prefix_consensus"]
    assert restored["c"]["cglin_raw"] == candidate["cglin_raw"]
    assert restored["c"]["cglin_provisional"] is True
    assert _matching(restored["q"], restored["c"], "cglin", 8)
    blocks, _ = partition_records(list(restored.values()), lin_level=5, context_lin_level=8)
    assert len(blocks[0]["records"]) == 2


def test_unresolved_requested_scope_has_conditional_human_readable_report(tmp_path):
    from chronoclade.cgmlst.report import write_partition_index

    unresolved = row(
        "q",
        "input",
        _assignment_result(
            {
                "cgST": "new",
                "LINcode": ["0", "0", "105", "6", "0", "0", "0", "5", "1", "0"],
                "matches": [{"LINcode": ["0", "0", "105", "6", "0", "0", "99", "0", "0", "0"]}],
            },
            CONFIG,
        ),
    )
    _, audit = partition_records([unresolved], lin_level=5, context_lin_level=9)
    path = write_partition_index(
        [],
        audit,
        records_count=1,
        lin_level=5,
        hiercc_level=None,
        status="partial",
        output=tmp_path,
    )
    text = path.read_text()
    assert "requested LIN9 context cannot be defined" in text
    assert "Requested LIN9 public context is unresolved" in text
    assert "missing_or_unresolved_context_scope_typing" not in text
    assert "zero retained contexts does not mean no relatives exist" in text
    assert "Status: partial" in text
    _, resolved_audit = partition_records([query()], lin_level=5, context_lin_level=8)
    path = write_partition_index(
        [],
        resolved_audit,
        records_count=1,
        lin_level=5,
        hiercc_level=None,
        status="complete",
        output=tmp_path,
    )
    text = path.read_text()
    assert "context cannot be defined" not in text
    assert "zero retained contexts" not in text


def test_native_unknown_hash_suffix_is_preserved_and_not_cast_to_allele_number():
    raw = {"cgST": "*novel", "LINcode": ["0", "0", "105", "6", "0", "*novel"]}
    assignment = _assignment_result(raw, CONFIG)
    assert assignment["cglin_raw"].startswith("0_0_105_6_0_*novel")
    assert resolved_prefix(assignment, 5) == PREFIX[:5]
    assert resolved_prefix(assignment, 6) is None
