import io
import re
from pathlib import Path

from Bio import Phylo

from chronoclade.report_components.group_tree import (
    linked_group_highlights,
    tree_viewer,
    write_group_tree,
)


def topology(tree):
    return [
        (node.name, node.branch_length, tuple(child.name for child in node.clades))
        for node in tree.find_clades(order="preorder")
    ]


def test_full_tree_preserves_topology_and_negative_lengths_and_uses_accessions(tmp_path):
    tree = Phylo.read(io.StringIO("((PW_123:-2,b:1):3,c:4);"), "newick")
    before = topology(tree)
    records = [
        {"sample_id": "PW_123", "country": "UK", "host": "human", "collection_date": "2020"},
        {"sample_id": "b", "country": "France", "host": "human", "collection_date": "2021"},
        {"sample_id": "c", "country": "Japan", "collection_date": "2022"},
    ]
    path = write_group_tree(
        tree, records, {"PW_123": "SRR1 (PW_123)", "b": "SAMN2", "c": "GCA_3.1"},
        {"PW_123": "cohort_1_group_2", "b": "cohort_1_group_1", "c": "cohort_1_group_2"},
        tmp_path / "groups.html", distance_units="mismatch fraction",
        group_basis="distance_threshold", group_threshold=0.02,
    )
    html = Path(path).read_text(encoding="utf-8")
    assert topology(tree) == before
    assert 'class="branch negative"' in html
    assert 'data-length="-2"' in html
    assert "SRR1" in html and "SAMN2" in html and "GCA_3.1" in html
    assert "PW_123" in html  # The private/raw ID remains data, but never appears in tip metadata text.
    assert "| ID: ${t.id}" not in html
    assert "Unknown host" in html
    assert "Groups use the declared direct-profile distance threshold (mismatch fraction 0.02), not tree clades." in html
    assert "mismatch fraction" in html
    assert "https://" not in html


def test_group_colors_follow_sorted_display_names_and_metadata_json_is_safe(tmp_path):
    tree = Phylo.read(io.StringIO("(a:1,b:1,c:1);"), "newick")
    records = [
        {"sample_id": "a", "origin": "input", "country": "<script>alert(1)</script>"},
        {"sample_id": "b", "origin": "context"},
        {"sample_id": "c", "origin": "context"},
    ]
    path = write_group_tree(
        tree, records, {"a": "Input A", "b": "Public B", "c": "Public C"},
        {"a": "cohort_1_group_2", "b": "cohort_1_group_1", "c": "cohort_1_group_2"},
        tmp_path / "safe.html",
    )
    html = Path(path).read_text(encoding="utf-8")
    # Lexical display order is Group 1, Group 2, mapped to the prescribed palette.
    assert re.search(r'"Group 1":"#2166ac","Group 2":"#b2182b"', html)
    assert "\\u003cscript\\u003e" in html
    assert "<script>alert(1)</script>" not in html
    assert 'class="tip-marker input"' in html
    assert 'class="tip-marker context"' in html
    assert "postMessage" in html


def test_tree_viewer_embeds_contained_file_and_parent_sync_routes_only_known_sources(tmp_path):
    (tmp_path / "cohort_a.html").write_text("<!doctype html>", encoding="utf-8")
    outside = tmp_path.parent / "outside.html"
    outside.write_text("<!doctype html>", encoding="utf-8")
    cohorts = [
        {"cohort_id": "cohort-a", "tree_groups_html": "cohort_a.html"},
        {"cohort_id": "missing", "tree_groups_html": "missing.html"},
        {"cohort_id": "outside", "tree_groups_html": outside},
    ]
    markup = tree_viewer(cohorts, tmp_path)
    assert 'data-group-cohort="cohort-a"' in markup
    assert 'src="cohort_a.html"' in markup
    assert 'title="NJ tree genetic groups"' in markup
    assert "cohort_1" not in markup
    assert "missing.html" not in markup and "outside.html" not in markup
    assert 'height:700px' in markup
    assert "Open the complete interactive tree" in markup
    assert tree_viewer(None, tmp_path) == ""
    script = linked_group_highlights()
    assert "frame.contentWindow===event.source" in script
    assert "frame.dataset.groupCohort===cohort" in script
    assert "typeof group!=='string'" in script and 'group.length>200' in script


def tree_payload(html):
    import json
    match = re.search(r'<script type="application/json" id="tree-data">(.*?)</script>', html)
    return json.loads(match.group(1))


def test_large_group_selector_is_bounded_and_full_tree_is_navigable(tmp_path):
    from Bio.Phylo.BaseTree import Clade, Tree

    names = [f"sample_{index}" for index in range(512)]
    tree = Tree(root=Clade(clades=[Clade(name=name, branch_length=index / 100)
                                  for index, name in enumerate(names)]))
    groups = {name: f"cohort_group_{index}" for index, name in enumerate(names)}
    records = [{"sample_id": name, "origin": "input" if index == 511 else "context"}
               for index, name in enumerate(names)]
    path = write_group_tree(tree, records, {}, groups, tmp_path / "large.html")
    html = Path(path).read_text()
    data = tree_payload(html)
    assert len(data["tips"]) == 512
    assert len(data["groups"]) == 512
    assert len(data["selector"]["compactGroups"]) == 21
    assert "cohort_group_511" in data["selector"]["compactGroups"]
    assert html.count('data-length="') == 512
    assert 'max-height:140px;overflow-y:auto' in html
    assert 'Search all groups' in html and 'Show all groups' in html
    assert 'height="400"' in html and 'height="10276"' not in html
    assert "let selected=null,located=null,mode='overview'" in html
    assert "function fitAll()" in html and 'w:data.width,h:data.height' in html
    assert "function focusTips(indexes)" in html
    assert "function readLabels()" in html
    assert "data.rowSpace*sy>=14" in html
    assert "matrix(${sx} 0 0 ${sy}" in html
    assert "translate(${x} ${y})" in html  # Unscaled text and marker layer.
    assert "scale(1 ${space/20})" not in html
    for control in ('zoom-in', 'zoom-out', 'zoom-inputs', 'zoom-group', 'locate-tip'):
        assert f'id="{control}"' in html
    assert "addEventListener('wheel'" in html
    assert "addEventListener('pointermove'" in html
    assert "event.shiftKey?event.deltaY:event.deltaX" in html
    assert "event.source!==window.parent" in html
    assert "data.tips[i].group_id===group" in html


def test_group_counts_rank_before_names_and_aliases_keep_identity(tmp_path):
    from Bio.Phylo.BaseTree import Clade, Tree
    from chronoclade.report_components.group_controls import (
        group_colors, group_label_map, group_selector_config,
    )

    ids = [f"cohort_group_{i}" for i in range(30)]
    labels = group_label_map(ids)
    rows = [{"group_id": name, "origin": "context"} for name in ids]
    rows += [{"group_id": ids[-1], "origin": "context"}] * 3
    config = group_selector_config(rows, group_key="group_id", labels=labels)
    assert config["compactGroups"][0] == ids[-1]
    assert len(config["compactGroups"]) == 20
    assert config["compactGroups"][1] == "cohort_group_0"
    raw = ["first_group_1", "second_group_1"]
    names = group_label_map(raw)
    assert names == {key: f"Group 1 [{key}]" for key in raw}
    tree = Tree(root=Clade(clades=[Clade(name="a", branch_length=1),
                                  Clade(name="b", branch_length=2)]))
    path = write_group_tree(tree, [{"sample_id": "b"}, {"sample_id": "a"}], {},
                            dict(zip(["a", "b"], raw)), tmp_path / "aliases.html")
    data = tree_payload(Path(path).read_text())
    assert [row["id"] for row in data["tips"]] == ["a", "b"]
    assert data["groups"] == raw
    assert data["selector"]["labels"] == names
    assert data["colors"] == group_colors(raw, labels=names)


def test_tree_optional_metadata_hides_absent_host_and_preserves_explicit_city(tmp_path):
    import json
    tree = Phylo.read(io.StringIO('(a:1,b:2);'), 'newick')
    before = topology(tree)
    records = [{'sample_id':'a','city':'Athens','host':'unknown','nuts2':'EL30'},
               {'sample_id':'b','city':None,'host':'not provided'}]
    path = write_group_tree(tree, records, {}, {}, tmp_path/'optional.html')
    html = Path(path).read_text()
    payload = json.loads(re.search(r'id="tree-data">(.*?)</script>', html).group(1))
    assert payload['metadata_fields'] == ['city']
    assert all('host' not in row and 'region' not in row for row in payload['tips'])
    assert payload['tips'][0]['city'] == 'Athens'
    assert payload['tips'][1]['city'] == 'Unknown city'
    assert 'EL30' not in html
    assert topology(tree) == before
