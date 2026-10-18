"""Tree/time host panels reuse the existing tree and offline network widget."""

import json
import re
from io import StringIO

from Bio import Phylo
import pytest

from chronoclade.location_network.reconstruction import (
    build_location_network,
    build_state_network,
)
from chronoclade.location_network.state_selector import state_network_selector
from chronoclade.profile_report import _state_network_views
from chronoclade.report_components.stage import write_report
from chronoclade.time_stage.workflow import run_time
from chronoclade.tree_stage.workflow import run_tree
from test_tree_stage import selection_fixture, native_stubs as _native_stubs


@pytest.fixture
def native_stubs(monkeypatch):
    return _native_stubs.__wrapped__(monkeypatch)


def report_result():
    return {
        "selected_sample_ids": ["a", "b", "c", "d"],
        "dating_status": "dated",
        "temporal_assessment": {"code": "not_assessed", "supported": False, "reason": "fixture"},
    }


def state_fixture():
    tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
    rows = [
        {"sample_id": ident, "country": "UK", "host": host, "role": role}
        for ident, host, role in (
            ("a", "Human", "input"),
            ("b", "Human", "context"),
            ("c", "Chicken", "context"),
            ("d", "Unknown", "context"),
        )
    ]
    return tree, rows


def test_country_default_and_country_model_are_unchanged_when_host_is_available(tmp_path):
    tree, rows = state_fixture()
    country_before = build_location_network(tree, rows, tree_basis="sequence")
    country = build_location_network(tree, rows, tree_basis="sequence")
    host = build_state_network(tree, rows, state_field="host", tree_basis="sequence")
    path = write_report(
        tmp_path / "report.html",
        title="Fixture",
        result=report_result(),
        rows=rows,
        network=country,
        host_network=host,
    )
    html = path.read_text()

    assert country == country_before
    assert country["optimum_changes"] == 0
    assert host["optimum_changes"] == 1
    assert 'data-state-panel="country">' in html
    assert html.index('data-state-panel="country">') < html.index('data-state-panel="host" hidden')
    assert 'value="country">Location (country)</option>' in html
    assert 'data-widget="host-network"' in html
    assert "Host colour key" in html and "Arrows show host-state changes" in html
    payloads = [json.loads(raw) for raw in re.findall(r"data-network-data>(.*?)</script>", html, re.S)]
    host_payload = next(payload for payload in payloads if payload["stateLabel"] == "host")
    assert {node["id"] for node in host_payload["nodes"]} == {"Human", "Chicken"}
    assert all("country" not in node for node in host_payload["nodes"])
    assert "data-network-action=\"export-png\"" in html
    assert "data-network-action=\"export-html\"" in html
    profile_panel = _state_network_views(
        {"state_networks": {"country": [country], "host": [host]},
         "metadata_summary": {"host": {"reported": 3}}},
        tmp_path,
    )
    assert "cc-state-network-selector" in profile_panel


def test_missing_host_panel_is_omitted_and_single_host_state_has_no_links(tmp_path):
    tree, rows = state_fixture()
    country = build_location_network(tree, rows, tree_basis="sequence")
    assert state_network_selector("country only") == "country only"

    no_host_rows = [{**row, "host": "Unknown"} for row in rows]
    no_host = build_state_network(tree, no_host_rows, state_field="host", tree_basis="sequence")
    assert no_host["known_state_count"] == 0
    no_host_path = write_report(
        tmp_path / "no-host.html", title="No host", result=report_result(),
        rows=no_host_rows, network=country,
    )
    no_host_html = no_host_path.read_text()
    assert "cc-state-network-selector" not in no_host_html
    assert "host-network" not in no_host_html

    one_host_rows = [{**row, "host": "Human"} for row in rows]
    one_host = build_state_network(tree, one_host_rows, state_field="host", tree_basis="sequence")
    assert one_host["status"] == "insufficient_metadata"
    assert one_host["known_state_count"] == 1
    assert one_host["directed_edges"] == []
    one_path = write_report(
        tmp_path / "one-host.html", title="One host", result=report_result(),
        rows=one_host_rows, network=country, host_network=one_host,
    )
    one_html = one_path.read_text()
    assert "1 distinct known host state" in one_html
    assert "At least two distinct known hosts are required" in one_html
    assert 'data-widget="host-network"' not in one_html


def test_tree_writes_optional_host_audit_on_same_tree_and_suppresses_absent_host(tmp_path, native_stubs):
    (tmp_path / "host").mkdir()
    selection, _ = selection_fixture(
        tmp_path / "host", dates=("", "", ""), hosts=("Human", "Chicken", "Unknown")
    )
    host_result = run_tree(selection, output=tmp_path / "tree-host", assess_temporal=False)
    host_path = tmp_path / "tree-host" / host_result["host_state_network"]["path"]
    host = json.loads(host_path.read_text())
    country = json.loads((tmp_path / "tree-host" / "jobs" / host_result["fingerprint"] / "location_network.json").read_text())
    report = (tmp_path / "tree-host" / host_result["report"]["path"]).read_text()
    assert host["tree_sha256"] == country["tree_sha256"]
    assert host["root_policy"] == country["root_policy"]
    assert host["reconstruction"]["root_id"] == country["reconstruction"]["root_id"]
    assert host["reconstruction"]["newick"] == country["reconstruction"]["newick"]
    assert host["status"] == "available"
    assert "cc-state-network-selector" in report and 'data-widget="host-network"' in report

    (tmp_path / "all-states").mkdir()
    selection, _ = selection_fixture(
        tmp_path / "all-states", dates=("", "", ""),
        hosts=("Human", "Chicken", "Swine"),
        countries=("UK", "UK", "France"),
        regions=("England", "Scotland", "Bretagne"),
        cities=("Oxford", "Oxford", "Rennes"),
    )
    all_output = tmp_path / "tree-all-states"
    all_result = run_tree(selection, output=all_output, assess_temporal=False)
    refs = all_result["state_networks"]
    assert set(refs) == {"country", "host", "region", "city"}
    assert refs["country"]["path"].endswith("location_network.json")
    assert all_result["host_state_network"] == refs["host"]
    country = json.loads((all_output / refs["country"]["path"]).read_text())
    for field in ("host", "region", "city"):
        ref = refs[field]
        assert ref["path"].endswith(f"{field}_state_network.json")
        assert ref["sha256"] == __import__("hashlib").sha256(
            (all_output / ref["path"]).read_bytes()
        ).hexdigest()
        state = json.loads((all_output / ref["path"]).read_text())
        assert state["state_field"] == field
        assert state["tree_sha256"] == country["tree_sha256"]
        assert state["branch_units"] == country["branch_units"]
        assert state["root_policy"] == country["root_policy"]
        assert state["reconstruction"]["root_id"] == country["reconstruction"]["root_id"]
        assert state["reconstruction"]["newick"] == country["reconstruction"]["newick"]
        if field == "region":
            assert {row["state"] for row in state["state_definitions"]} == {
                "UK / England", "UK / Scotland", "France / Bretagne"
            }
        if field == "city":
            assert {row["state"] for row in state["state_definitions"]} == {
                "UK / England / Oxford", "UK / Scotland / Oxford",
                "France / Bretagne / Rennes",
            }
    report = (all_output / all_result["report"]["path"]).read_text()
    assert all(f'data-state-panel="{field}"' in report for field in ("country", "host", "region", "city"))
    assert 'value="region">Region</option>' in report
    assert 'value="city">City</option>' in report

    (tmp_path / "no-host").mkdir()
    no_host_selection, _ = selection_fixture(
        tmp_path / "no-host", dates=("", "", "")
    )
    no_host_result = run_tree(
        no_host_selection, output=tmp_path / "tree-no-host", assess_temporal=False
    )
    no_host_job = tmp_path / "tree-no-host" / "jobs" / no_host_result["fingerprint"]
    no_host_report = (tmp_path / "tree-no-host" / no_host_result["report"]["path"]).read_text()
    assert "host_state_network" not in no_host_result
    assert set(no_host_result["state_networks"]) == {"country"}
    assert not (no_host_job / "host_state_network.json").exists()
    assert "cc-state-network-selector" not in no_host_report


def test_dated_time_report_uses_same_host_topology_without_running_native_tools(
    tmp_path, native_stubs, monkeypatch
):
    (tmp_path / "dated").mkdir()
    selection, _ = selection_fixture(
        tmp_path / "dated", dates=("2018", "2019", "2020"),
        hosts=("Human", "Chicken", "Unknown"),
        regions=("England", "Scotland", "Bretagne"),
        cities=("Oxford", "Glasgow", "Rennes"),
    )
    run_tree(selection, output=tmp_path / "dated-tree", assess_temporal=False)

    def fake_dated_tree(files, *args, **kwargs):
        source = Phylo.read(files.clock_dir / "rerooted.newick", "newick")
        files.time_tree.parent.mkdir(parents=True, exist_ok=True)
        Phylo.write(source, files.time_tree, "nexus")

    monkeypatch.setattr("chronoclade.time_stage.workflow.require_tools", lambda names: None)
    monkeypatch.setattr("chronoclade.time_stage.workflow.tools", lambda names: {name: {"path": None} for name in names})
    monkeypatch.setattr("chronoclade.lineage._run_dated_tree", fake_dated_tree)
    monkeypatch.setattr("chronoclade.temporal_report.write_timetree_confidence", lambda job: {"status": "fixture"})
    monkeypatch.setattr("chronoclade.temporal_report.write_timetree_figures", lambda *args, **kwargs: None)

    time_result = run_time(
        tmp_path / "dated-tree/tree.json", output=tmp_path / "dated-time", allow_unsupported=True
    )
    host_path = tmp_path / "dated-time" / time_result["host_state_network"]["path"]
    host = json.loads(host_path.read_text())
    country_path = tmp_path / "dated-time" / "jobs" / time_result["fingerprint"] / "location_network.json"
    country = json.loads(country_path.read_text())
    report = (tmp_path / "dated-time" / time_result["report"]["path"]).read_text()
    assert time_result["dating_status"] == "dated"
    assert host["tree_sha256"] == __import__("hashlib").sha256(
        (tmp_path / "dated-time" / time_result["dated_tree"]["path"]).read_bytes()
    ).hexdigest()
    assert host["root_policy"] == country["root_policy"]
    assert host["reconstruction"]["root_id"] == country["reconstruction"]["root_id"]
    assert host["reconstruction"]["newick"] == country["reconstruction"]["newick"]
    assert "cc-state-network-selector" in report and 'data-widget="host-network"' in report
    assert set(time_result["state_networks"]) == {"country", "host", "region", "city"}
    assert time_result["host_state_network"] == time_result["state_networks"]["host"]
    for field in ("host", "region", "city"):
        state_ref = time_result["state_networks"][field]
        state = json.loads((tmp_path / "dated-time" / state_ref["path"]).read_text())
        assert state["state_field"] == field
        assert state["tree_sha256"] == country["tree_sha256"]
        assert state["branch_units"] == country["branch_units"] == "years"
        assert state["root_policy"] == country["root_policy"]
        assert state["reconstruction"]["root_id"] == country["reconstruction"]["root_id"]
        assert state["reconstruction"]["newick"] == country["reconstruction"]["newick"]
        if field == "region":
            assert {row["state"] for row in state["state_definitions"]} == {
                "Greece / England", "Greece / Scotland", "UK / Bretagne"
            }
        if field == "city":
            assert {row["state"] for row in state["state_definitions"]} == {
                "Greece / England / Oxford", "Greece / Scotland / Glasgow",
                "UK / Bretagne / Rennes",
            }
    assert 'data-state-panel="region" hidden' in report
    assert 'data-state-panel="city" hidden' in report
