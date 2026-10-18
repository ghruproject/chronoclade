from chronoclade.temporal_diagnostics import date_regression, allele_unit_diagnostic


def test_allele_unit_scale_keeps_fraction_evidence_and_residuals():
    diagnostic = date_regression([dict(sample_id=str(i), collection_date=str(2000+i), d=i*.01)
                                  for i in range(3)], distance_field='d',
                                 distance_units='fraction', label='test')
    scaled = allele_unit_diagnostic(diagnostic, 100)
    assert scaled['slope'] == diagnostic['slope'] * 100
    assert scaled['fraction_slope'] == diagnostic['slope']
    assert scaled['points'][1]['d'] == 1
    assert diagnostic['points'][1]['d'] == .01
    assert 'scheme-equivalent' in scaled['slope_units']
    fixed = allele_unit_diagnostic(diagnostic, 100, fixed_callable=True)
    assert 'fixed callable' in fixed['distance_units']


def test_reader_labels_are_short_without_discarding_units_or_changing_values():
    from chronoclade.profile_report import _root_to_tip

    original = date_regression(
        [dict(sample_id=str(i), collection_date=str(2000+i), d=i*.01)
         for i in range(3)], distance_field="d", distance_units="fraction", label="test"
    )
    for fixed in (False, True):
        scaled = allele_unit_diagnostic(original, 629, fixed_callable=fixed)
        assert scaled["display_distance_units"] == "Estimated allele differences"
        assert scaled["display_slope_units"] == "Estimated allele differences/year"
        assert "629" in scaled["slope_units"]
        assert scaled["slope"] == original["slope"] * 629
        html = _root_to_tip([{**scaled, "cohort_id": "cohort_1",
                              "midpoint_date_slope": scaled["slope"]}])
        assert "Estimated allele differences/year" in html
        assert scaled["slope_units"] not in html
    assert "display_slope_units" not in original


def test_unscaled_dna_or_legacy_diagnostic_keeps_declared_display_units():
    from chronoclade.profile_report import _root_to_tip

    html = _root_to_tip([dict(cohort_id="cohort_1", slope_units="SNP differences per year",
                              midpoint_date_slope=1.2)])
    assert "SNP differences per year" in html
    assert "Estimated allele differences/year" not in html
