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
    )
    html = Path(path).read_text(encoding="utf-8")
    assert topology(tree) == before
    assert 'class="branch negative"' in html
    assert 'data-length="-2"' in html
    assert "SRR1" in html and "SAMN2" in html and "GCA_3.1" in html
    assert "PW_123" in html  # The private/raw ID remains data, but never appears in tip metadata text.
    assert "| ID: ${t.id}" not in html
    assert "Unknown host" in html
    assert "Groups defined from direct profile distances, not tree clades." in html
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
