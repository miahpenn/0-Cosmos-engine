import numpy as np

from tools.audit_centre_ricci_fixed_s2 import (
    _fixed_s2_fit,
    _fitted_even_polynomial,
    _moment_from_fit,
    _parity_composite_simpson_integral,
)


def test_fixed_s2_fit_recovers_higher_even_coefficients_without_refitting_s0_s2():
    r = (np.arange(32, dtype=float) + 0.5) * 0.125
    s0 = 1.0
    s2 = -4.0 * np.pi * 0.01**2
    s4, s6, s8 = 2.0e-4, -1.0e-5, 4.0e-7
    source = s0 + s2*r**2 + s4*r**4 + s6*r**6 + s8*r**8
    fit = _fixed_s2_fit(source, r, 3.5)
    np.testing.assert_allclose(
        [fit["s0_fixed"], fit["s2_fixed"], fit["s4_fitted"], fit["s6_fitted"], fit["s8_fitted"]],
        [s0, s2, s4, s6, s8],
        rtol=3.0e-7,
        atol=2.0e-12,
    )
    assert fit["fit_max_abs_source_residual"] < 5.0e-11
    assert fit["coefficients_finite"] is True


def test_fixed_s2_half_cell_integral_matches_exact_even_polynomial():
    r = (np.arange(32, dtype=float) + 0.5) * 0.125
    r0 = float(r[0])
    s0 = 1.0
    s2 = -4.0 * np.pi * 0.01**2
    s4, s6, s8 = 2.0e-4, -1.0e-5, 4.0e-7
    source = s0 + s2*r**2 + s4*r**4 + s6*r**6 + s8*r**8
    fit = _fixed_s2_fit(source, r, 3.5)
    moment, _ = _moment_from_fit(fit, r0, "fixed_s2")
    exact = (
        s0*r0 + s2*r0**3/3.0 + s4*r0**5/5.0
        + s6*r0**7/7.0 + s8*r0**9/9.0
    )
    assert abs(moment-exact) < 2.0e-16


def test_fixed_s2_fit_requires_sufficient_cells():
    r = (np.arange(4, dtype=float) + 0.5) * 1.0
    source = 1.0 - 4.0*np.pi*0.01**2*r*r
    try:
        _fixed_s2_fit(source, r, 1.0)
    except ValueError:
        pass
    else:
        raise AssertionError("under-sampled constrained fit must fail closed")


def test_degree5_comparison_retains_s10_and_integrates_all_six_even_coefficients():
    r = (np.arange(64, dtype=float) + 0.5) * 0.125
    r_fit = 3.5
    coefficients = np.array([1.0, -2.0e-4, 1.0e-5, -3.0e-7, 1.0e-8, -2.0e-7])
    source = sum(value * r ** (2 * power) for power, value in enumerate(coefficients))

    fit = _fitted_even_polynomial(source, r, float(r[0]), r_fit, degree=5)
    assert fit["coefficient_names"] == ["s0", "s2", "s4", "s6", "s8", "s10"]
    assert len(fit["coefficients_even_powers"]) == 6
    np.testing.assert_allclose(fit["coefficients_even_powers"], coefficients, rtol=1.0e-5, atol=5.0e-12)

    r0 = 1.0
    moment, meta = _moment_from_fit(fit, r0, "fitted")
    exact = sum(value * r0 ** (2 * power + 1) / (2 * power + 1)
                for power, value in enumerate(coefficients))
    assert abs(moment - exact) < 1.0e-11
    assert meta["integrated_coefficient_names"][-1] == "s10"
    assert len(meta["integrated_coefficients_even_powers"]) == 6



def test_parity_composite_simpson_integrates_even_quadratic_source():
    for dr in (0.5, 1.0/3.0, 0.25):
        r = (np.arange(32, dtype=float) + 0.5)*dr
        s0 = 1.0
        s2 = -4.0*np.pi*0.01**2
        source = s0 + s2*r**2
        integral = _parity_composite_simpson_integral(
            source, r, dr, origin_value=s0
        )
        B = integral/r
        np.testing.assert_allclose(
            B[:8], s0+(s2/3.0)*r[:8]**2, rtol=0.0, atol=3.0e-14
        )
