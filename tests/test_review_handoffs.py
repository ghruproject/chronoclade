"""Review regressions: tell readers which next stages can actually run."""

import json

from chronoclade.cgmlst.report import _tree_eligibility
from chronoclade.prepare_stage import _next_stage_summary, _parallel_prepare_report, _report
from chronoclade.profile_report import write_profile_report
from chronoclade.report_components.profile_summary import subsample_summary
from chronoclade.report_components.profile_overview import reported_lineage_annotations
from chronoclade.report_components.relatedness import scope_summary
from chronoclade.stage_readme import write_stage_readme


def _selection(tmp_path, count):
    bundle = tmp_path / "selection_bundle"
    bundle.mkdir()
    ids = [f"sample-{i}" for i in range(count)]
    (bundle / "selection.json").write_text(json.dumps({
        "query_ids": ids, "selected_sample_ids": ids, "selection_id": "test",
    }))
    (bundle / "ensemble.json").write_text(json.dumps({
        "selections": [{"path": "selection.json", "selection_id": "test"}],
    }))


def test_small_selection_report_keeps_evidence_without_unrunnable_commands(tmp_path):
    _selection(tmp_path, 2)
    document = subsample_summary(tmp_path, {})
    assert "requires at least three selected genomes" in document
    assert "chronoclade tree" not in document
    assert 'href="selection_bundle/selection.json" download' in document
    # Eligibility uses selected genomes, not the total available profile pool.
    status = _tree_eligibility({"ensemble": "selection_bundle/ensemble.json",
                                "context_count": 100}, tmp_path)
    assert "Blocked: 1/1" in status


def test_eligible_selection_report_has_tree_handoff(tmp_path):
    _selection(tmp_path, 3)
    assert "chronoclade tree" in subsample_summary(tmp_path, {})
    assert "Eligible: 1/1" in _tree_eligibility(
        {"ensemble": "selection_bundle/ensemble.json"}, tmp_path)


def test_download_attributes_do_not_replace_interactive_or_figure_navigation(tmp_path):
    (tmp_path / "tree.nwk").write_text("(a:1,b:1);")
    (tmp_path / "figure.svg").write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
    document = write_profile_report({
        "records_count": 2,
        "cohorts": [{"cohort_id": "c", "sample_ids": ["a", "b"],
                     "tree_path": "tree.nwk", "tree_figure": "figure.svg"}],
    }, directory=tmp_path).read_text()
    assert 'href="tree.nwk" download' in document
    assert 'data-source-asset="figure.svg"' in document
    assert 'href="data:image/svg+xml' in document


def test_readiness_explains_incomplete_catalogue_separately_from_missing_calls():
    row = {"role": "input", "sample_id": "q", "profile_ready": True,
           "profiles": [{"catalogue_complete": False}]}
    text = _next_stage_summary({"samples": [row]})
    assert "Continue with cgmlst" in text
    assert "catalogue is incomplete" in text
    assert "Missing cgLIN/HierCC labels do not prevent" in text
    row["profile_ready"] = False
    assert "lack called cgMLST profiles" in _next_stage_summary({"samples": [row]})


def test_slurm_readme_recovery_keeps_executor_and_ensemble_layout(tmp_path):
    plan = tmp_path / "slurm jobs" / "array.json"
    write_stage_readme(tmp_path, "tree", ensemble=True, completed_runs=0, executor_plan=plan)
    text = (tmp_path / "README.md").read_text()
    assert "`runs/*/jobs/`" in text
    assert "submit.sh" in text and str(plan) in text
    assert "Repeat the original command" not in text
    assert "rerun the original tree command" not in text
    assert "chronoclade time" not in text
    write_stage_readme(tmp_path, "tree")
    assert "`jobs/`" in (tmp_path / "README.md").read_text()


def test_profile_scope_does_not_claim_mandatory_lineage_or_prefix_matching():
    text = scope_summary({"partition_by": "profiles"}, {}, provider={
        "pools": [{"discovered": 10}], "accessible_matches": 8,
    })
    assert "compatible cgMLST scheme" in text
    assert "8 compatible source records" in text
    assert "codes and ST assignments are not required" in text
    assert "LIN5 matching" not in text and "prefix-matching" not in text


def test_ensemble_readme_links_reports_using_saved_result_paths(tmp_path):
    run = tmp_path / 'runs' / 'saved selection'
    (run / 'jobs' / 'fingerprint').mkdir(parents=True)
    (run / 'jobs' / 'fingerprint' / 'report.html').write_text('report')
    (run / 'tree.json').write_text(json.dumps({
        'report': {'path': 'jobs/fingerprint/report.html'},
    }))
    (tmp_path / 'tree_ensemble.json').write_text(json.dumps({'runs': [{
        'status': 'complete', 'selection_id': 'selection-001',
        'result': {'path': 'runs/saved selection/tree.json'},
    }]}))
    write_stage_readme(tmp_path, 'tree', ensemble=True, completed_runs=1)
    assert '[selection-001 report](runs/saved%20selection/jobs/fingerprint/report.html)' in (
        tmp_path / 'README.md').read_text()


def test_reported_codes_are_escaped_and_distinguished_from_distance_group_assignments():
    text = reported_lineage_annotations({
        'sample_labels': {'q': 'ERR123'},
        'reported_lineage_annotations': [{
            'sample_id': 'q', 'cglin_raw': '<script>', 'cglin_scheme': 'scheme-a',
            'hiercc_codes': {'HC1100': 12}, 'hiercc_status': 'resolved',
        }],
    })
    assert 'ERR123' in text and '&lt;script&gt;' in text and '<script>' not in text
    assert 'HC1100: 12' in text and 'scheme-a' in text
    assert 'do not define the distance groups' in text


def test_prepare_report_separates_available_profiles_from_optional_lineage_completeness():
    audit = {
        'samples': [{'sample_id': 'q', 'label': 'ERR123', 'role': 'input',
                     'profile_ready': True, 'mlst_st_ready': True, 'mlst_st': '149',
                     'mlst_scheme': 'saureus', 'mlst_status': 'resolved',
                     'profiles': [{'called_loci': 2195, 'loci': 2208,
                                   'scheme_id': 'saureus_1', 'catalogue_complete': True}],
                     'readiness': 'incomplete',
                     'reasons': ['Lineage readiness is not assessed for this species'],
                     'missing_metadata': [], 'assembly_status': 'bundled'}],
        'operations': {'provider_retrieval': 'not_run_local_assemblies',
                       'typing': 'run_missing_profile_inputs',
                       'reference_database': 'ready_validated_configuration'},
    }
    original = json.dumps(audit, sort_keys=True)
    document = _report(audit)
    assert 'Called missing input cgMLST profiles' in document
    assert 'run_missing_profile_inputs' not in document
    assert '<th scope="col">Profile availability</th>' in document
    assert '<td>Available</td>' in document and '<td>incomplete</td>' not in document
    assert 'Optional lineage annotations' in document
    assert 'Lineage readiness is not assessed for this species' in document
    assert json.dumps(audit, sort_keys=True) == original
    audit['samples'][0]['profile_ready'] = False
    assert '<td>Unavailable</td>' in _report(audit)


def test_parallel_prepare_summary_links_retained_worker_evidence_without_claiming_collector_calls(tmp_path):
    directory = tmp_path / 'sources/provider/task-000000'
    for relative in ['provider/native_typing/query_typing.json', 'classical_mlst/audit.json']:
        path = directory / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('{}')
    audit = {'samples': [{'sample_id': 'a', 'label': 'ERR123', 'role': 'input',
                          'profile_ready': True, 'mlst_st_ready': True},
                         {'sample_id': 'b', 'label': 'ERR456', 'role': 'input',
                          'profile_ready': False, 'mlst_st_ready': True}],
             'operations': {'typing': 'not_run_frozen_import'}}
    original = json.dumps(audit, sort_keys=True)
    summary = _parallel_prepare_report(audit, {'tasks': [
        {'sample_id': 'a', 'retained_evidence': 'sources/provider/task-000000'},
        {'sample_id': 'b', 'retained_evidence': 'sources/provider/missing-worker'},
    ]}, tmp_path)
    assert '2 checked per-assembly worker results' in summary
    assert '1 / 2 inputs have cgMLST profiles and 2 / 2 have classical MLST STs' in summary
    assert 'collector imported validated worker results; it did not perform typing' in summary
    assert 'href="sources/provider/task-000000/provider/native_typing/query_typing.json"' in summary
    assert 'href="sources/provider/task-000000/classical_mlst/audit.json"' in summary
    assert 'missing-worker/' not in summary
    assert json.dumps(audit, sort_keys=True) == original
