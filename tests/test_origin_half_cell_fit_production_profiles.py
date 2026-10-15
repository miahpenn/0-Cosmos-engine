"""Unit tests for the saved-production-profile half-cell sensitivity audit."""
import math
import numpy as np
import pytest

from tools.audit_origin_half_cell_fit_production_profiles import (
    analyze_artifact,
    fit_half_cell_moment_from_samples,
)


def polynomial_source(r):
    r = np.asarray(r, dtype=float)
    return 1.0 + 0.3 * r**2 - 0.04 * r**4 + 0.004 * r**6


@pytest.mark.parametrize("dr", [0.5, 1.0 / 3.0, 0.25])
def test_degree_three_even_polynomial_half_cell_moment_is_exact(dr):
    n = int(round(80.0 / dr))
    r = (np.arange(n) + 0.5) * dr
    source = polynomial_source(r)
    fit = fit_half_cell_moment_from_samples(source, dr, 2.5, 3)
    r0 = dr / 2.0
    exact = r0 + 0.3 * r0**3 / 3.0 - 0.04 * r0**5 / 5.0 + 0.004 * r0**7 / 7.0
    assert math.isclose(fit["fitted_half_cell_moment"], exact, rel_tol=0.0, abs_tol=2.0e-13)
    assert fit["sample_count"] >= 5


def test_local_six_cell_window_uses_six_samples_and_preserves_finite_diagnostics():
    dr = 0.25
    n = 320
    r = (np.arange(n) + 0.5) * dr
    source = polynomial_source(r) + 2.0e-4 * np.exp(-r**2)
    fit = fit_half_cell_moment_from_samples(source, dr, 6.0 * dr, 3)
    assert fit["sample_count"] == 6
    assert np.isfinite(fit["signed_correction_vs_point_sample"])
    assert np.isfinite(fit["max_abs_loo_source_prediction_residual"])


def test_saturated_fit_is_rejected_instead_of_presented_as_validated():
    source = np.ones(20)
    with pytest.raises(ValueError, match=r"degree\+2"):
        fit_half_cell_moment_from_samples(source, dr=0.5, fit_radius=2.0, degree_in_r_squared=3)


def test_wrong_artifact_digest_fails_closed(tmp_path):
    wrong = tmp_path / "wrong.zip"
    wrong.write_bytes(b"not the pinned production artifact")
    with pytest.raises(RuntimeError, match="artifact ZIP SHA-256 mismatch"):
        analyze_artifact(wrong)
