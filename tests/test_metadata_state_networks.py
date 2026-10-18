"""Categorical reconstructions share a tree, but never borrow other observations."""

import copy
import json
import re
from io import StringIO

import pytest
from Bio import Phylo

from chronoclade.location_network.reconstruction import (
    build_location_network, build_optional_state_networks, build_state_network, observed_state,
)
from chronoclade.location_network.state_selector import state_network_panel, state_network_selector
from chronoclade.location_network.widget import interactive_network_html


def example():
    tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
    rows = [dict(sample_id=name, country=country, region=region, city=city, host=host,
                 role="input" if name in {"a", "c"} else "context")
            for name, country, region, city, host in [
                ("a", "Greece", "Attica", "Athens", "Human"),
                ("b", "Greece", "Attica", "Athens", "Human"),
                ("c", "Italy", "Lazio", "Rome", "Chicken"),
                ("d", "Italy", "Lazio", "Rome", None),
            ]]
    return tree, rows


@pytest.mark.parametrize("field", ["country", "host", "region", "city"])
def test_each_field_reconstructs_its_own_observations_on_same_tree(field):
    tree, rows = example()
    before = copy.deepcopy(rows)
    country = build_location_network(tree, rows, tree_basis="sequence")
    model = build_state_network(tree, rows, state_field=field, tree_basis="sequence")
    assert model["optimum_changes"] == 1
    assert model["reconstruction"]["newick"] == country["reconstruction"]["newick"]
    assert rows == before
    assert build_location_network(tree, rows, tree_basis="sequence") == country


def test_homonymous_places_keep_parent_geography_and_missing_targets_stay_unknown():
    tree, rows = example()
    rows = [dict(row, region="Central", city="Springfield") for row in rows]
    for field in ["region", "city"]:
        model = build_state_network(tree, rows, state_field=field, tree_basis="sequence")
        states = {node["state"] for node in model["nodes"]}
        assert len(states) == 2
        assert any(state.startswith("Greece / Central") for state in states)
        assert any(state.startswith("Italy / Central") for state in states)
    assert observed_state({"country": "Greece", "region": "Attica", "city": None}, "city") is None
    assert observed_state({"country": "Greece", "region": None, "city": "Athens"}, "region") is None


def test_omitted_fields_and_outside_tree_records_cannot_add_states():
    tree, rows = example()
    rows = [dict(row, host=None, region="N/A", city=None) for row in rows]
    rows.append(dict(sample_id="outside", country="France", host="Bird", region="Paris", city="Paris"))
    assert build_optional_state_networks(tree, rows, tree_basis="sequence") == {}


def test_generic_categorical_field_and_city_widget_use_state_values():
    tree, rows = example()
    custom = build_state_network(tree, [dict(r, source="clinical" if r["sample_id"] in {"a", "b"}
                                            else "environmental") for r in rows],
                                 state_field="source", tree_basis="sequence")
    assert custom["optimum_changes"] == 1
    city = build_state_network(tree, rows, state_field="city", tree_basis="sequence")
    markup = interactive_network_html(city, "cities")
    assert 'data-widget="city-network"' in markup and "All cities" in markup
    payload = json.loads(re.search(r"data-network-data>(.*?)</script>", markup, re.S)[1])
    assert payload["stateLabel"] == "city"
    assert {node["state"] for node in payload["nodes"]} == {
        "Greece / Attica / Athens", "Italy / Lazio / Rome"
    }
    assert all("country" not in node for node in payload["nodes"])


def test_snp_nearest_comparisons_use_the_selected_state_and_callable_sites():
    tree, rows = example()
    networks = build_optional_state_networks(
        tree, rows, tree_basis="sequence", nearest_neighbours=[
            dict(focal_sample="a", context_sample="d", clonal_snps=3, callable_sites=1000)
        ]
    )
    edge = networks["city"]["nearest_edges"][0]
    assert (edge["source"], edge["target"]) == ("Greece / Attica / Athens", "Italy / Lazio / Rome")
    assert edge["snp_differences_min"] == 3 and edge["shared_sites_min"] == 1000


def test_shared_selector_omits_absent_fields_and_explains_single_state():
    tree, rows = example()
    optional = build_optional_state_networks(tree, rows, tree_basis="sequence")
    panels = {field: state_network_panel(model, basis="Test") for field, model in optional.items()}
    markup = state_network_selector("Country panel", panels=panels)
    for field in ["country", "host", "region", "city"]:
        assert f'value="{field}"' in markup
    one = build_state_network(tree, [dict(r, city="Athens", region="Attica", country="Greece")
                                    for r in rows], state_field="city", tree_basis="sequence")
    panel = state_network_panel(one)
    assert "At least two distinct known cities" in panel
    assert "data-network-canvas" not in panel
    assert state_network_selector(None, panels={"host": "Host panel"}) == "Host panel"
