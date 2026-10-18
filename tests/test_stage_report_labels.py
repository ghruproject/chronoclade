from chronoclade.report_components.stage import write_report
from chronoclade.cgmlst.report import write_partition_index


def test_shared_report_uses_public_labels_in_table_and_neighbours(tmp_path):
    path = tmp_path/'report.html'
    upstream = tmp_path/'cgMLST & context.html'
    upstream.write_text('Saved upstream report')
    write_report(path, title='Labels', rows=[
        {'sample_id':'PW_a', 'label':'PW_a', 'accession':'ERR123456'},
        {'sample_id':'PW_b', 'label':'ERR234567 (PW_b)'}], result={
            'selected_sample_ids':['PW_a','PW_b'],
            'temporal_assessment':{'code':'not_assessed','supported':False,'reason':'fixture'},
            'nearest_contexts_by_corrected_snps':[{'focal_sample':'PW_a', 'context_sample':'PW_b',
                                                 'clonal_snps':2, 'callable_sites':1000}]},
                 navigation=(("cgMLST report", upstream), ("Missing report", tmp_path/'missing.html')))
    html = path.read_text()
    assert '>ERR123456</td>' in html
    assert '>ERR234567</td>' in html
    assert '>PW_a</td>' not in html
    assert '(PW_b)</td>' not in html
    assert 'href="cgMLST%20%26%20context.html"' in html
    assert 'Missing report' not in html
    assert 'href="#evidence"' in html and 'id="evidence"' in html


def test_partition_index_names_the_active_partition_without_cross_lineage_labels(tmp_path):
    audit = {'unresolved_inputs': [], 'available_level_counts': []}
    path = write_partition_index([], audit, records_count=2, lin_level=5,
        hiercc_level='HC10', status='complete', output=tmp_path)
    html = path.read_text()
    assert 'Analysis: HC10' in html
    assert 'LIN depth' not in html and 'cgLIN' not in html

    st_audit = {**audit, 'partition_by': 'st',
                'policy': 'supplied-inputs-species-st-profile-scope',
                'context_scope': 'Supplied inputs only; no public context retrieval or lineage refinement'}
    st_path = write_partition_index([], st_audit, records_count=2, lin_level=5,
        hiercc_level=None, status='complete', output=tmp_path)
    st_html = st_path.read_text()
    assert 'supplied profiles by species/ST' in st_html
    assert 'no lineage assignment required' in st_html
    assert 'Supplied inputs only; no public context retrieval or lineage refinement' in st_html
    assert 'LIN depth' not in st_html and 'HierCC' not in st_html

    lineage_path = write_partition_index([], audit, records_count=2, lin_level=7,
        hiercc_level=None, status='complete', output=tmp_path)
    assert 'Analysis: LIN7' in lineage_path.read_text()
