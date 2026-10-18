import json

from chronoclade.report_components.profile_summary import (
    combined_dataset_summary,
    subsample_summary,
)


def test_combined_dataset_summary_puts_country_totals_first_and_collapses_subregions():
    html = combined_dataset_summary({
        "metadata_geography": [
            {"origin": "local", "country": "Kenya", "region": "Nairobi", "city": "Nairobi", "count": 2},
            {"origin": "context", "country": "Kenya", "region": "Nairobi", "city": "Nairobi", "count": 3},
            {"origin": "local", "country": "Uganda", "region": "Central", "city": "Kampala", "count": 1},
        ],
        "metadata_summary": {
            "host": {"total": 6, "reported": 5, "missing": 1},
            "collection_date": {"total": 6, "reported": 4, "missing": 2},
        },
    })

    assert html.index("Country totals") < html.index("Location detail")
    assert '<details class="coverage-details"><summary>Location detail</summary>' in html
    assert "Kenya" in html and "5" in html and "3" in html
    assert "Host and collection-date coverage" in html
    assert "<script>" not in combined_dataset_summary({
        "metadata_geography": [{"origin": "local", "country": "<script>", "count": 1}]
    })


def test_nuts2_is_not_used_as_city_or_as_an_empty_location_detail():
    html = combined_dataset_summary({"metadata_geography": [
        {"origin": "local", "country": "Unknown", "region": "Unknown", "nuts2": "Unknown", "count": 2},
        {"origin": "context", "country": "Greece", "region": "Unknown", "nuts2": "EL54", "count": 3},
    ]})
    assert "Unknown" in html
    assert "Location detail" not in html
    assert "EL54" not in html
    assert "NUTS2" not in html


def test_city_only_metadata_has_country_and_city_without_an_empty_region_column():
    html = combined_dataset_summary({"metadata_geography": [
        {"origin": "local", "country": "Greece", "city": "Athens", "count": 2},
        {"origin": "context", "country": "Greece", "city": "Athens", "count": 3},
    ]})
    detail = html.split("Location detail", 1)[1]
    assert 'scope="col">Country' in detail and 'scope="col">City' in detail
    assert 'scope="col">Region' not in detail
    assert "Athens" in detail and ">5</td>" in detail


def test_unknown_only_host_and_location_fields_are_omitted():
    html = combined_dataset_summary({
        "metadata_geography": [{"country": "Greece", "region": "not provided", "city": "Unknown city", "count": 3}],
        "metadata_summary": {"host": {"total": 3, "reported": 0, "missing": 3},
                             "collection_date": {"total": 3, "reported": 0, "missing": 3}},
        "metadata_composition": [{"field": "Host", "value": "Unknown", "total": 3}],
    })
    assert "Location detail" not in html
    assert "Host" not in html
    assert "Collection-date coverage" in html


def test_partial_host_metadata_keeps_reported_hosts_and_missing_counts():
    html = combined_dataset_summary({
        "metadata_composition": [{"field": "Host", "value": "Human", "total": 2},
                                 {"field": "Host", "value": "Unknown", "total": 1}],
        "metadata_summary": {"host": {"total": 3, "reported": 2, "missing": 1}},
    })
    assert "Host breakdown" in html and "Human" in html and "Unknown" in html


def test_region_detail_is_hidden_when_only_unknown_locations_are_reported():
    html = combined_dataset_summary({"metadata_geography": [
        {"origin": "local", "country": "Kenya", "region": "Unknown", "nuts2": "", "count": 1},
    ]})

    assert "Country totals" in html
    assert "Kenya" in html
    assert "Location detail" not in html


def test_subsample_summary_reports_exact_counts_links_and_collapsed_metadata(tmp_path):
    bundle = tmp_path / "selection_bundle"
    (bundle / "selections").mkdir(parents=True)
    selection = {
        "selection_id": "selection-001",
        "query_ids": ["q1"],
        "available_context_ids": ["c1", "c2", "c3"],
        "selected_context_ids": ["c1", "c2"],
        "selected_sample_ids": ["q1", "c1", "c2"],
        "required_nearest_context_ids": ["c1"],
        "nearest_ranking": "raw_snp_differences_when_available",
        "decisions": [
            {"sample_id": "c1", "reasons": ["snp_neighbour:q1:rank=1"]},
            {"sample_id": "c2", "reasons": ["genetic_group_time_region_representative"]},
        ],
    }
    (bundle / "selections" / "selection-001.json").write_text(json.dumps(selection))
    ensemble = {
        "selections": [{"selection_id": "selection-001", "path": "selections/selection-001.json"}],
        "pairwise_overlap": [],
    }
    (bundle / "ensemble.json").write_text(json.dumps(ensemble))

    html = subsample_summary(tmp_path, {
        "sample_metadata": [
            {"sample_id": "q1", "country": "Kenya", "collection_date": "2020", "host": "human"},
            {"sample_id": "c1", "country": "Kenya", "collection_date": "2021", "host": "human"},
            {"sample_id": "c2", "country": "Uganda", "host": "cattle"},
            {"sample_id": "c3", "country": "Uganda", "collection_date": "2019", "host": "cattle"},
        ]
    })

    assert "1 input genomes + 2 context genomes = 3 genomes" in html
    assert 'href="selection_bundle/selections/selection-001.json"' in html
    assert "Mandatory nearest" in html and "1" in html
    assert "raw differences" in html and "boundary ties" in html
    assert "Distinct alternatives" not in html and "Overlap" not in html
    assert "Country, collection-date and host comparison" not in html
    assert '<details class="coverage-details"><summary>selection-001 metadata comparison</summary>' in html
    assert "2 / 3 reported" in html and "3 / 4 reported" in html


def test_subsample_summary_explains_absent_bundle(tmp_path):
    html = subsample_summary(tmp_path, {})
    assert "No subsamples were prepared for this dataset." in html



def test_multiple_subsamples_show_distinctness_and_actual_overlap(tmp_path):
    bundle = tmp_path / "selection_bundle"
    (bundle / "selections").mkdir(parents=True)
    runs = []
    for n, contexts in enumerate((["c1", "c2"], ["c1", "c3"]), 1):
        path = f"selections/selection-{n:03d}.json"
        (bundle / path).write_text(json.dumps({"selection_id": str(n), "query_ids": ["q"],
            "available_context_ids": ["c1", "c2", "c3"], "selected_context_ids": contexts,
            "selected_sample_ids": ["q"] + contexts}))
        runs.append({"path": path})
    (bundle / "ensemble.json").write_text(json.dumps({"selections": runs,
        "pairwise_overlap": [{"jaccard": 1/3}]}))
    html = subsample_summary(tmp_path, {})
    assert "2 / 2" in html
    assert "3 unique genomes" in html
    assert "0.33–0.33" in html
