from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from audit_centre_ricci_origin_series import fit_even_series


def test_even_series_fit_recovers_known_coefficients_on_fixed_radius_window():
    r = (np.arange(64, dtype=float) + 0.5) * 0.125
    expected = (4.5e-4, -2.0e-5, 1.25e-6)
    values = expected[0] + expected[1] * r**2 + expected[2] * r**4
    result = fit_even_series(r, values, 3.0)
    np.testing.assert_allclose(
        [result["coefficients"]["c0"], result["coefficients"]["c2"], result["coefficients"]["c4"]],
        expected, rtol=2e-11, atol=2e-14,
    )
    assert result["max_abs_residual"] < 2e-15


def test_excluding_first_cell_estimates_regular_continuation_and_detects_outlier():
    r = (np.arange(32, dtype=float) + 0.5) * 0.25
    expected = (4.5e-4, -2.0e-5, 1.25e-6)
    values = expected[0] + expected[1] * r**2 + expected[2] * r**4
    values[0] += 8.0e-4
    all_cells = fit_even_series(r, values, 3.0)
    nonfirst = fit_even_series(r, values, 3.0, exclude_first_cell=True)
    np.testing.assert_allclose(
        [nonfirst["coefficients"]["c0"], nonfirst["coefficients"]["c2"], nonfirst["coefficients"]["c4"]],
        expected, rtol=2e-11, atol=2e-14,
    )
    assert all_cells["max_abs_residual"] > 5e-5
    np.testing.assert_allclose(
        nonfirst["first_cell_extrapolation"]["signed_observed_minus_fitted"],
        8.0e-4, rtol=2e-11, atol=2e-14,
    )


def test_even_series_fit_uses_fixed_physical_radius_not_fixed_cell_count():
    coarse_r = (np.arange(40, dtype=float) + 0.5) * 0.5
    fine_r = (np.arange(80, dtype=float) + 0.5) * 0.25
    profile = lambda r: 4.5e-4 - 2.0e-5*r**2 + 1.25e-6*r**4
    coarse = fit_even_series(coarse_r, profile(coarse_r), 3.0)
    fine = fit_even_series(fine_r, profile(fine_r), 3.0)
    assert coarse["fit_cells"] == 6
    assert fine["fit_cells"] == 12
    np.testing.assert_allclose(
        [coarse["coefficients"]["c0"], coarse["coefficients"]["c2"], coarse["coefficients"]["c4"]],
        [fine["coefficients"]["c0"], fine["coefficients"]["c2"], fine["coefficients"]["c4"]],
        rtol=2e-11, atol=2e-14,
    )
