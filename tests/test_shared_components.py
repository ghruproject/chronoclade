"""Shared foundations must remain usable independently of analysis stages."""

import subprocess
import sys

import pytest

from chronoclade.errors import WorkflowError
from chronoclade.report_components.tables import enhance_report_tables
from chronoclade.selections import selected_context_records


def test_shared_style_and_selection_imports_do_not_load_scientific_stages():
    subprocess.run(
        [sys.executable, "-c", """
import sys
from chronoclade.report_components.styles import report_styles
from chronoclade.selections import selected_context_records
from chronoclade.datasets import PreparedDataset, load_dataset
assert "@font-face" in report_styles()
for name in (
    "chronoclade.temporal_report", "chronoclade.profile_analysis",
    "chronoclade.cgmlst.workflow", "chronoclade.native_pipeline", "chronoclade.pathogenwatch",
    "chronoclade.profile_types",
):
    assert name not in sys.modules, name
"""],
        check=True,
    )


def test_saved_selection_preserves_identity_and_does_not_mutate_pool():
    queries = [{"sample_id": "query"}]
    context = [{"sample_id": "a", "accession": "SAMN1"}, {"sample_id": "b"}]
    selection = dict(
        available_context_ids=["a", "b"], query_ids=["query"],
        selected_context_ids=["a"], selected_sample_ids=["query", "a"],
        decisions=[dict(sample_id="a", reason="user_requested")],
    )
    assert selected_context_records(selection, queries, context) == [
        dict(sample_id="a", accession="SAMN1", selection_reason="user_requested")
    ]
    assert "selection_reason" not in context[0]
    selection["decisions"][0].pop("reason")
    with pytest.raises(WorkflowError, match="decisions"):
        selected_context_records(selection, queries, context)


def test_large_report_tables_default_to_ten_rows_and_keep_page_size_choices():
    html = enhance_report_tables("<!doctype html><html><head></head><body></body></html>")

    assert "page = 0, size = 10" in html
    assert "[10,25,50,100,'All']" in html
    assert "select.value='10'" in html
