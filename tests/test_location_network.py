"""Shared tree/location boundaries, independent of cgMLST subgroup selection."""

from io import StringIO

import pytest
from Bio import Phylo

from chronoclade.datasets import PreparedDataset, samples_from_records
from chronoclade.location_network.reconstruction import build_location_network
from chronoclade.location_network.tree import draw_country_tree
from chronoclade.profile_network import build_profile_network


def biological_tree():
    return Phylo.read(StringIO("((q:0.1,c:0.1):0.2,u:0.3);"), "newick")


def records():
    return [
        dict(sample_id="q", country="Greece", origin="query"),
        dict(sample_id="c", country="Italy", origin="context"),
        dict(sample_id="u", country=None, origin="context"),
    ]


def test_tree_provenance_is_required_and_unknown_distances_are_rejected():
    with pytest.raises(TypeError, match="tree_basis"):
        build_location_network(biological_tree(), records())
    for basis in ("unknown", "distance_matrix", None):
        with pytest.raises(ValueError, match="declared biological tree basis"):
            build_location_network(biological_tree(), records(), tree_basis=basis)
    with pytest.raises(TypeError, match="Bio.Phylo biological tree"):
        build_location_network([[0, 1], [1, 0]], records(), tree_basis="sequence")


def test_sequence_and_independently_provided_trees_share_exact_location_calculation():
    cgmlst = build_location_network(biological_tree(), records(), tree_basis="cgmlst")
    for basis in ("sequence", "provided_phylogeny"):
        shared = build_location_network(biological_tree(), records(), tree_basis=basis)
        assert shared["tree_basis"] == basis
        assert "NJ" not in shared["direction_scope"]
        # Provenance wording changes; the same supplied tree yields exactly the
        # same history, uncertainty, sample counts, palette and graph metrics.
        shared["direction_scope"] = cgmlst["direction_scope"]
        shared["tree_basis"] = cgmlst["tree_basis"]
        assert shared == cgmlst
        assert "views" not in shared  # LIN group ownership stays with cgMLST.


def test_prepared_roles_preserve_input_counts_roots_and_tree_marker(tmp_path):
    legacy = records()
    dataset = PreparedDataset(samples=samples_from_records(legacy))
    prepared = dataset.samples
    assert all("origin" not in row for row in prepared)
    before = build_location_network(biological_tree(), legacy, tree_basis="sequence")
    after = build_location_network(biological_tree(), prepared, tree_basis="sequence")
    assert before == after
    assert after["tested_roots"][0] == "q"
    assert after["input_ids"] == ["q"]
    # Render both role boundaries to ensure input outlining and bold tip labels
    # remain identical in addition to the scientific counts.
    for rows, name in ((legacy, "legacy.png"), (prepared, "prepared.png")):
        draw_country_tree(after, rows, tmp_path / name)
    assert (tmp_path / "legacy.png").read_bytes() == (tmp_path / "prepared.png").read_bytes()


def test_prepared_role_controls_lin_subgroups_and_overrides_source_origin():
    rows = [
        dict(
            row,
            role="input" if row["origin"] == "query" else "context",
            cglin_raw="0,0,0,0,1",
            cglin_scheme="scgMLST629_S",
            cglin_scheme_version="v1",
        )
        for row in records()
    ]
    rows[1]["origin"] = "query"  # Provenance cannot override the prepared role.
    data = build_profile_network(biological_tree(), rows)
    assert data["input_ids"] == ["q"]
    assert data["views"][0]["input_ids"] == ["q"]
    assert data["nodes"][1]["input_count"] == 0


def test_shared_tree_renderer_accepts_method_specific_scale_and_dating_text(tmp_path):
    network = build_location_network(biological_tree(), records(), tree_basis="sequence")
    path = tmp_path / "sequence-tree.svg"
    draw_country_tree(
        network,
        records(),
        path,
        title="Country-coloured sequence tree",
        distance_label="Substitutions per site",
        dating_note="",
    )
    svg = path.read_text()
    assert "Country-coloured sequence tree" in svg
    assert "Substitutions per site" in svg
    assert "Showing 3 of 3 genomes" in svg
    assert "cgMLST" not in svg
    assert "profiles" not in svg
    assert "NJ tree" not in svg
    assert "This tree is not dated" not in svg


@pytest.mark.parametrize("plot_kind", ["country", "root_comparison"])
def test_worker_plots_override_gui_backend_in_fresh_process(tmp_path, plot_kind):
    """Fresh workers must not inherit a GUI backend from macOS/user configuration."""
    import os
    from pathlib import Path
    import subprocess
    import sys
    import textwrap

    program = textwrap.dedent('''
        from concurrent.futures import ThreadPoolExecutor
        from io import StringIO
        from pathlib import Path
        import sys
        from threading import Barrier
        from Bio import Phylo
        from chronoclade.location_network.reconstruction import build_location_network
        from chronoclade.location_network.tree import draw_country_tree
        from chronoclade.report_components.root_comparison import _figure
        assert 'matplotlib.pyplot' not in sys.modules
        output = Path(sys.argv[1])
        plot_kind = sys.argv[2]
        gate = Barrier(2)
        def plot(index):
            rows = [dict(sample_id='q', country='UK', role='input'),
                    dict(sample_id='c', country='Italy', role='context'),
                    dict(sample_id='u', country='Greece', role='context')]
            tree = Phylo.read(StringIO('((q:0.1,c:0.1):0.2,u:0.3);'), 'newick')
            network = build_location_network(tree, rows, tree_basis='sequence')
            gate.wait()
            path = output / f'tree-{index}.svg'
            if plot_kind == "country":
                draw_country_tree(network, rows, path)
            else:
                diagnostic = dict(points=[dict(midpoint_year=2020, year_min=2019, year_max=2021, root_to_tip=0.1)], centered_intercept=0.1, slope=0.01, reference_year=2020, r_squared=0.5, distance_units="substitutions/site")
                _figure(diagnostic, diagnostic, path)
            assert '<svg' in path.read_text()
            return path
        with ThreadPoolExecutor(max_workers=2) as workers:
            paths = list(workers.map(plot, range(2)))
        import matplotlib
        assert matplotlib.get_backend().lower() == 'agg'
        assert len(paths) == 2
    ''')
    environment = {**os.environ, "MPLBACKEND": "MacOSX", "MPLCONFIGDIR": str(tmp_path / "mpl")}
    completed = subprocess.run(
        [sys.executable, "-c", program, str(tmp_path), plot_kind],
        cwd=Path(__file__).resolve().parents[1], env=environment,
        capture_output=True, text=True, timeout=60,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert all((tmp_path / f"tree-{index}.svg").stat().st_size > 1000 for index in range(2))
