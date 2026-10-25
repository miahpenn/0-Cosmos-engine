import numpy as np

from tools.audit_exact_profile_momentum_tangent import _apply_exact_stage_terms, _summary


def test_apply_exact_stage_terms_replaces_native_explicit_rates_and_keeps_rhs_terms():
    rates = {
        "geometry": {
            "a": np.array([1.0, 1.0]),
            "b": np.array([2.0, 2.0]),
            "X": np.array([3.0, 3.0]),
            "Aa": np.array([4.0, 4.0]),
            "K": np.array([5.0, 5.0]),
            "Lambda": np.array([6.0, 6.0]),
        },
        "scalars": {"S": np.array([7.0, 7.0])},
        "matter": {},
    }
    terms = {
        "explicit": {
            "a": np.array([11.0, 11.0]),
            "b": np.array([12.0, 12.0]),
            "X": np.array([13.0, 13.0]),
        },
        "primary_l2": {"Aa": np.array([14.0, 14.0]), "K": np.array([15.0, 15.0])},
        "primary_l3": {"Aa": np.array([16.0, 16.0]), "K": np.array([17.0, 17.0])},
        "lambda_l2": np.array([18.0, 18.0]),
        "lambda_l3": np.array([19.0, 19.0]),
    }
    out = _apply_exact_stage_terms(rates, terms)
    assert np.array_equal(out["geometry"]["a"], [11.0, 11.0])
    assert np.array_equal(out["geometry"]["b"], [12.0, 12.0])
    assert np.array_equal(out["geometry"]["X"], [13.0, 13.0])
    assert np.array_equal(out["geometry"]["Aa"], [30.0, 30.0])
    assert np.array_equal(out["geometry"]["K"], [32.0, 32.0])
    assert np.array_equal(out["geometry"]["Lambda"], [37.0, 37.0])
    assert np.array_equal(out["scalars"]["S"], [7.0, 7.0])
    assert np.array_equal(rates["geometry"]["a"], [1.0, 1.0])
    assert out["geometry_stage_terms"] is terms


def test_summary_keeps_origin_cells_and_physical_annulus_separate():
    r = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 11.5])
    values = np.array([-1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0])
    record = _summary(values, r)
    assert record["signed_cells_0_1_2_6"] == {
        "0": -1.0, "1": 2.0, "2": -3.0, "6": -7.0
    }
    assert record["max_abs_full_grid"] == 7.0
    assert record["max_abs_r_le_10"] == 6.0
    assert record["max_abs_fixed_annulus_r5_to10"] == 6.0
