"""Native query progress identifies coverage without changing typing evidence."""

import json
import logging

import pytest

from chronoclade import query_typing
from chronoclade.metadata import Sample
from chronoclade.progress import configure_progress


@pytest.mark.parametrize("profiles_only", [False, True])
def test_two_queries_show_progress_and_quiet_preserves_results(
    tmp_path, monkeypatch, capsys, profiles_only
):
    caller = {
        "command": ["caller"],
        "scheme": "klebsiella_1",
        "scheme_version": "fixture-v1",
        "index_dir": tmp_path,
        "index_loci": ["a", "b", "c", "d"],
    }
    assigner = {
        "kind": "plincer",
        "command": ["plincer"],
        "canonical_loci": ["a", "b", "c", "d"],
        "scheme_file": "scheme",
        "profiles_file": "profiles",
        "alleles_db": "alleles",
    }
    monkeypatch.setattr(
        query_typing,
        "load_typing_config",
        lambda *a, **k: {
            "organisms": {"Klebsiella pneumoniae": {"cgmlst": caller, "assignment": assigner}}
        },
    )
    monkeypatch.setattr(
        query_typing, "execution_provenance", lambda *a: {"identity_sha256": "fixture"}
    )
    monkeypatch.setattr(
        query_typing, "_provenance", lambda *a: {"databases": {"index_dir": {"sha256": "a" * 64}}}
    )
    monkeypatch.setattr(query_typing, "_assignment_result", lambda *a: {"cglin_resolved_depth": 10})
    calls = []

    def native(command, payload, destination):
        calls.append(command[0])
        if command[0] == "plincer":
            return {"fixture": "assignment"}
        first = b"first" in payload
        return {
            "scheme": "klebsiella_1",
            "schemeSize": 4,
            "genes": ["a", "b", "c", "d"],
            "code": "1__" + "b" * 40 + "_3" if first else "1_2_3_4",
            "alleles": {"d": [{"id": 3}, {"id": 4}]} if first else {},
        }

    monkeypatch.setattr(query_typing, "_run", native)
    samples = []
    for ident in ("first", "second"):
        assembly = tmp_path / f"{ident}.fa"
        assembly.write_text(f">{ident}\nACGT\n")
        samples.append(
            Sample(ident, assembly, "2020", "Greece", "Klebsiella pneumoniae", "ST258", "local")
        )
    (tmp_path / "config").write_text("{}")
    logger = logging.getLogger("chronoclade.progress")
    previous = (logger.handlers[:], logger.level, logger.propagate)
    try:
        configure_progress(quiet=False)
        normal = query_typing.type_query_assemblies(
            samples, tmp_path / "config", tmp_path / "normal", profiles_only=profiles_only
        )
        messages = capsys.readouterr().err
        assert "Typing sample 1/2: first" in messages
        assert "Typing sample 2/2: second" in messages
        assert (
            "Completed typing sample 1/2: first; 2/4 loci called; 1 missing; 1 ambiguous"
            in messages
        )
        assert (
            "Completed typing sample 2/2: second; 4/4 loci called; 0 missing; 0 ambiguous"
            in messages
        )
        assert "b" * 40 not in messages
        configure_progress(quiet=True)
        quiet = query_typing.type_query_assemblies(
            samples, tmp_path / "config", tmp_path / "quiet", profiles_only=profiles_only
        )
        assert capsys.readouterr().err == ""
    finally:
        logger.handlers, logger.level, logger.propagate = previous
    assert normal == quiet
    assert [row["sample_id"] for row in normal] == ["first", "second"]
    assert normal[0]["cgmlst_profile"] == {"a": "1"}
    assert normal[0]["cgmlst_novel_alleles"] == {"c": "b" * 40}
    assert normal[0]["cgmlst_missing_loci"] == ["b"]
    assert normal[0]["cgmlst_ambiguous_loci"] == ["d"]
    assert len(calls) == (4 if profiles_only else 8)
    for name in ("normal", "quiet"):
        saved = json.loads((tmp_path / name / "query_typing.json").read_text())
        assert saved["assignments"] == normal
