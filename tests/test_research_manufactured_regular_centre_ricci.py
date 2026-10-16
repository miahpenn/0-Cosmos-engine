"""Regression controls for the independent analytic regular-centre Ricci oracle."""
from __future__ import annotations

from tools.audit_manufactured_regular_centre_ricci import analyze


def test_manufactured_regular_centre_error_classes_and_wrong_parity_control():
    report = analyze((80, 120, 160, 240, 320))
    assert report["diagnostic_contract_passed"] is True
    assert report["scientific_admission"] is False
    assert report["production_evolution_steps"] == 0
    orders = report["descriptive_orders_on_manufactured_profile"]
    assert all(3.5 <= row["order"] <= 4.5 for row in orders["bssn_cell0_error"])
    assert all(1.7 <= row["order"] <= 2.3 for row in orders["metric_cell0_error"])
    assert all(1.7 <= row["order"] <= 2.3 for row in orders["route_gap_cell0"])
    assert all(3.5 <= row["order"] <= 4.5 for row in orders["bssn_annulus_error"])
    assert all(3.5 <= row["order"] <= 4.5 for row in orders["metric_annulus_error"])
    assert all(
        row["wrong_parity_cell0_error_abs"] > 1.0e3 * row["metric_first2_max_abs_error"]
        for row in report["rows"]
    )
    assert all(
        row["determinant_max_abs"] < 2.0e-14 and row["all_finite"]
        for row in report["rows"]
    )
