"""Regression tests for strict-JSON failure evidence in the PIRK L2 observer."""
import json

import numpy as np
import pytest

from tools.audit_tracefree_ricci_proper_radial import _write_strict_json_report, _independent_metric_crosscheck, _d4


def test_nonfinite_diagnostic_failure_still_writes_strict_json(tmp_path):
    report = {
        "status": "diagnostic_crosscheck_failure",
        "errors": ["N=80: nonfinite observer input"],
        "resolution_reports": [
            {
                "resolution": 80,
                "all_checks_pass": False,
                "stage_summaries": [
                    {
                        "stage": "final_implicit_l2_rhs",
                        "finite": False,
                        "independent_metric_crosscheck": {
                            "ricci_gap_max_abs_full": float("nan"),
                            "ricci_gap_peak_full": float("inf"),
                            "direct_d2_hessian_gap_max_abs_full": np.float64("-inf"),
                            "centre_witness": np.array([0.125, np.nan, -np.inf]),
                        },
                    }
                ],
            }
        ],
    }
    output = tmp_path / "report.json"

    _write_strict_json_report(output, report)

    raw = output.read_text(encoding="utf-8")
    decoded = json.loads(
        raw,
        parse_constant=lambda token: pytest.fail(f"non-standard JSON constant: {token}"),
    )
    stage = decoded["resolution_reports"][0]["stage_summaries"][0]
    assert decoded["status"] == "diagnostic_crosscheck_failure"
    assert decoded["errors"] == ["N=80: nonfinite observer input"]
    assert stage["finite"] is False
    metric = stage["independent_metric_crosscheck"]
    assert metric["ricci_gap_max_abs_full"] is None
    assert metric["ricci_gap_peak_full"] is None
    assert metric["direct_d2_hessian_gap_max_abs_full"] is None
    assert metric["centre_witness"] == [0.125, None, None]
    assert "NaN" not in raw
    assert "Infinity" not in raw


def test_manufactured_constant_curvature_regular_profile_is_negative_control():
    from engine.v55_pirk_adapter import vendor_modules
    grid_ops, vacuum, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(160, 8.0)
    r = np.asarray(grid.centers, dtype=float)
    k = 1.0e-3
    f = 1.0 / np.sqrt(1.0-k*r*r)
    geometry = vacuum.flat_state(grid)
    geometry.a[:] = f**(4.0/3.0)
    geometry.b[:] = f**(-2.0/3.0)
    geometry.X[:] = f**(-1.0/3.0)
    geometry.alpha[:] = 1.0
    geometry.beta[:] = 0.0
    d1 = lambda values, parity: grid.cell_derivative_fourth(values, parity=parity)
    geometry.Lambda[:] = d1(geometry.a,1)/(2.0*geometry.a*geometry.a)-d1(geometry.b,1)/(geometry.a*geometry.b)+2.0/r*(1.0/geometry.b-1.0/geometry.a)
    terms = vacuum.geometry_terms(grid,geometry)
    capture = {"r":r,"a":geometry.a,"b":geometry.b,"X":geometry.X,"alpha":geometry.alpha,"Lambda":geometry.Lambda,"R":terms["R"],"Rr_r":terms["Rr_r"],"DrDr_alpha":terms["DrDr_alpha"],"lap_alpha":terms["lap_alpha"]}
    arrays, summary = _independent_metric_crosscheck(capture,float(grid.dr))
    interior = r <= 2.0
    assert np.max(np.abs(arrays["ricci_tracefree_metric"][interior])) < 5.0e-6
    assert np.max(np.abs(arrays["ricci_tracefree_vendor"][interior])) < 5.0e-6
    assert np.max(np.abs(arrays["ricci_tracefree_gap"][interior])) < 5.0e-6
    assert summary["ricci_gap_max_abs_centre_i0_i1"] == pytest.approx(np.max(np.abs(arrays["ricci_tracefree_gap"][:2])))
    annulus = summary["ricci_gap_fixed_areal_annulus_R5_R8"]
    assert annulus["count"] > 0 and annulus["max_abs"] is not None

def test_manufactured_nonconstant_regular_profile_matches_nonzero_ricci_oracle():
    """Validate both curvature routes against a smooth profile with nonzero TF Ricci."""
    from engine.v55_pirk_adapter import vendor_modules

    grid_ops, vacuum, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(160, 8.0)
    r = np.asarray(grid.centers, dtype=float)

    # Use proper distance ell=r and areal radius R(r)=r(1+c3*r^2).
    # The closed-form radial trace-free Ricci eigenvalue is:
    #   (2/3) * ((R'^2 - 1)/R^2 - R''/R).
    # Unlike the constant-curvature control above, this is nonzero for r>0.
    c3 = 5.0e-2
    areal_exact = r * (1.0 + c3 * r * r)
    areal_prime_exact = 1.0 + 3.0 * c3 * r * r
    areal_second_exact = 6.0 * c3 * r
    tf_ricci_from_derivatives = (2.0 / 3.0) * (
        (areal_prime_exact**2 - 1.0) / areal_exact**2
        - areal_second_exact / areal_exact
    )
    # Algebraically simplify the same oracle to avoid cancellation near r=0:
    # TF Ricci = 2*c3^2*r^2/(1+c3*r^2)^2, hence TF Ricci/r^2 -> 2*c3^2.
    tf_ricci_exact = 2.0 * c3**2 * r**2 / (1.0 + c3 * r**2) ** 2
    tf_ricci_over_r2_exact = 2.0 * c3**2 / (1.0 + c3 * r**2) ** 2
    # The unsimplified expression subtracts O(c3) terms to leave O(c3^2*r^2)
    # at the first cell; allow its ~1e-13 absolute cancellation roundoff. The
    # centre checks below use the stable expression, not this cancellation path.
    assert np.allclose(tf_ricci_from_derivatives, tf_ricci_exact, rtol=1.0e-11, atol=2.0e-13)

    # Represent the same physical metric with the production BSSN
    # determinant-one conformal parameterization:
    # q=R/r, a=q^(-4/3), b=q^(2/3), X=q^(-2/3).
    # This preserves sqrt(a)/X=1, R_areal=r*sqrt(b)/X, and a*b^2=1.
    q_exact = areal_exact / r
    geometry = vacuum.flat_state(grid)
    geometry.X[:] = q_exact ** (-2.0 / 3.0)
    geometry.a[:] = q_exact ** (-4.0 / 3.0)
    geometry.b[:] = q_exact ** (2.0 / 3.0)
    geometry.alpha[:] = 1.0
    geometry.beta[:] = 0.0
    assert np.allclose(geometry.a * geometry.b**2, 1.0, rtol=0.0, atol=1.0e-13)
    assert np.allclose(geometry.a / geometry.X**2, 1.0, rtol=0.0, atol=1.0e-13)
    assert np.allclose(r * np.sqrt(geometry.b) / geometry.X, areal_exact, rtol=1.0e-13, atol=1.0e-13)

    # Parity control: areal radius is odd, while the radial metric factor is
    # even. The declared odd parity differentiates the analytic cubic exactly
    # in the first two cells; deliberately breaking parity must not pass.
    h = float(grid.dr)
    proper_parity_areal_prime = _d4(areal_exact, h, -1, 1)
    broken_parity_areal_prime = _d4(areal_exact, h, +1, 1)
    centre_first_two = np.arange(r.size) < 2
    assert np.max(np.abs(proper_parity_areal_prime[centre_first_two] - areal_prime_exact[centre_first_two])) < 1.0e-10
    assert np.max(np.abs(broken_parity_areal_prime[centre_first_two] - areal_prime_exact[centre_first_two])) > 1.0e-2

    d1 = lambda values, parity: grid.cell_derivative_fourth(values, parity=parity)
    geometry.Lambda[:] = (
        d1(geometry.a, 1) / (2.0 * geometry.a * geometry.a)
        - d1(geometry.b, 1) / (geometry.a * geometry.b)
        + 2.0 / r * (1.0 / geometry.b - 1.0 / geometry.a)
    )
    terms = vacuum.geometry_terms(grid, geometry)
    capture = {
        "r": r,
        "a": geometry.a,
        "b": geometry.b,
        "X": geometry.X,
        "alpha": geometry.alpha,
        "Lambda": geometry.Lambda,
        "R": terms["R"],
        "Rr_r": terms["Rr_r"],
        "DrDr_alpha": terms["DrDr_alpha"],
        "lap_alpha": terms["lap_alpha"],
    }
    arrays, summary = _independent_metric_crosscheck(capture, float(grid.dr))

    # Keep the centre witness distinct: the existing production issue is
    # localized there, whereas this oracle check is for a well-resolved annulus.
    centre = r < 0.2
    annulus = (r >= 1.0) & (r <= 2.0)
    assert np.any(centre) and np.any(annulus)
    assert np.all(np.isfinite(tf_ricci_exact))
    assert np.all(np.isfinite(tf_ricci_over_r2_exact))
    assert np.all(np.isfinite(arrays["ricci_tracefree_metric"]))
    assert np.all(np.isfinite(arrays["ricci_tracefree_vendor"]))
    assert np.max(np.abs(tf_ricci_exact[annulus])) > 1.0e-3

    # Sign control: this profile's radial trace-free Ricci is positive. A
    # globally sign-flipped implementation must fail, not pass on |Ricci|.
    assert np.all(tf_ricci_exact[annulus] > 0.0)
    assert np.all(arrays["ricci_tracefree_metric"][annulus] > 0.0)
    assert np.all(arrays["ricci_tracefree_vendor"][annulus] > 0.0)

    # Explicitly test the first two centre cells against the analytic oracle
    # using TF Ricci/r^2, which has the finite nonzero limit 2*c3^2 at r=0.
    # Keep the signed checks; norms alone would hide a sign error.
    for key in ("ricci_tracefree_metric", "ricci_tracefree_vendor"):
        scaled = arrays[key] / r**2
        centre_error = scaled[centre_first_two] - tf_ricci_over_r2_exact[centre_first_two]
        assert np.all(arrays[key][centre_first_two] > 0.0), f"{key}: wrong centre sign"
        assert np.max(np.abs(centre_error)) < 5.0e-4, f"{key}: centre oracle mismatch {centre_error}"

    # Fourth-order stencils should reproduce this smooth, analytic profile to
    # substantially better than 5% on the annulus.
    atol = 2.0e-4
    metric_error = arrays["ricci_tracefree_metric"][annulus] - tf_ricci_exact[annulus]
    vendor_error = arrays["ricci_tracefree_vendor"][annulus] - tf_ricci_exact[annulus]
    assert np.max(np.abs(metric_error)) < atol
    assert np.max(np.abs(vendor_error)) < atol
    assert summary["ricci_gap_max_abs_centre_i0_i1"] == pytest.approx(
        np.max(np.abs(arrays["ricci_tracefree_gap"][:2]))
    )

def test_manufactured_regular_centre_series_validates_c5_oracle_and_witnesses():
    """Check the r^5 centre-series term with a stable analytic Ricci oracle."""
    from engine.v55_pirk_adapter import vendor_modules

    grid_ops, vacuum, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(160, 8.0)
    r = np.asarray(grid.centers, dtype=float)
    x = r**2

    # Proper distance ell=r and areal radius
    # R(ell)=ell+c3*ell^3+c5*ell^5. For
    # Ric_TF^r_r=(2/3)*[(R'^2-1)/R^2 - R''/R], expansion gives
    # (2*c3^2 - 20*c5/3)*r^2 + O(r^4).
    c3 = 5.0e-2
    c5 = -2.0e-4
    q_exact = 1.0 + c3 * x + c5 * x**2
    areal_exact = r * q_exact
    areal_prime_exact = 1.0 + 3.0 * c3 * x + 5.0 * c5 * x**2
    areal_second_exact = r * (6.0 * c3 + 20.0 * c5 * x)

    # Cancellation-stable exact form: factor out r^2 before evaluating.
    centre_polynomial = (
        (3.0 * c3**2 - 10.0 * c5)
        + 4.0 * c3 * c5 * x
        + 5.0 * c5**2 * x**2
    )
    tf_ricci_exact = (2.0 / 3.0) * x * centre_polynomial / q_exact**2
    tf_ricci_over_r2_exact = (2.0 / 3.0) * centre_polynomial / q_exact**2
    centre_coefficient = 2.0 * c3**2 - (20.0 / 3.0) * c5
    assert tf_ricci_over_r2_exact[0] == pytest.approx(centre_coefficient, abs=1.0e-5)
    assert abs(centre_coefficient - 2.0 * c3**2) > 1.0e-3  # c5 is actually exercised.

    # Verify the stable closed form against the unsimplified geometric identity.
    # The latter can lose precision close to the origin, so this tolerance is
    # limited to the analytic identity and is not reused for stencil validation.
    tf_from_derivatives = (2.0 / 3.0) * (
        (areal_prime_exact**2 - 1.0) / areal_exact**2
        - areal_second_exact / areal_exact
    )
    assert np.allclose(tf_from_derivatives, tf_ricci_exact, rtol=1.0e-10, atol=3.0e-13)

    # Same physical geometry in the determinant-one BSSN representation:
    # a*b^2=1, a/X^2=1 (so proper distance is ell=r), and R=r*sqrt(b)/X.
    geometry = vacuum.flat_state(grid)
    geometry.X[:] = q_exact ** (-2.0 / 3.0)
    geometry.a[:] = q_exact ** (-4.0 / 3.0)
    geometry.b[:] = q_exact ** (2.0 / 3.0)
    geometry.alpha[:] = 1.0
    geometry.beta[:] = 0.0
    assert np.allclose(geometry.a * geometry.b**2, 1.0, rtol=0.0, atol=1.0e-13)
    assert np.allclose(geometry.a / geometry.X**2, 1.0, rtol=0.0, atol=1.0e-13)
    assert np.allclose(r * np.sqrt(geometry.b) / geometry.X, areal_exact, rtol=1.0e-13, atol=1.0e-13)

    # Retain an explicit parity-breaking negative control. The r^5 term is
    # smooth and odd; the proper odd-parity stencil should recover R' at the
    # first cells while the deliberately even assignment should not.
    h = float(grid.dr)
    proper_parity_prime = _d4(areal_exact, h, -1, 1)
    broken_parity_prime = _d4(areal_exact, h, +1, 1)
    centre_first_two = np.arange(r.size) < 2
    assert np.max(np.abs(proper_parity_prime[centre_first_two] - areal_prime_exact[centre_first_two])) < 1.0e-8
    assert np.max(np.abs(broken_parity_prime[centre_first_two] - areal_prime_exact[centre_first_two])) > 1.0e-2

    d1 = lambda values, parity: grid.cell_derivative_fourth(values, parity=parity)
    geometry.Lambda[:] = (
        d1(geometry.a, 1) / (2.0 * geometry.a * geometry.a)
        - d1(geometry.b, 1) / (geometry.a * geometry.b)
        + 2.0 / r * (1.0 / geometry.b - 1.0 / geometry.a)
    )
    terms = vacuum.geometry_terms(grid, geometry)
    capture = {
        "r": r,
        "a": geometry.a,
        "b": geometry.b,
        "X": geometry.X,
        "alpha": geometry.alpha,
        "Lambda": geometry.Lambda,
        "R": terms["R"],
        "Rr_r": terms["Rr_r"],
        "DrDr_alpha": terms["DrDr_alpha"],
        "lap_alpha": terms["lap_alpha"],
    }
    arrays, summary = _independent_metric_crosscheck(capture, float(grid.dr))
    annulus = (r >= 1.0) & (r <= 2.0)
    assert np.any(annulus)
    assert np.all(np.isfinite(tf_ricci_exact))
    assert np.all(np.isfinite(tf_ricci_over_r2_exact))
    assert np.max(np.abs(tf_ricci_exact[annulus])) > 1.0e-3

    # The five-point first derivative has a known exact truncation term on
    # this quintic: D4[R] = R' - h^4*R^(5)/30 = R' - 4*c5*h^4.
    # The five-point second derivative is exact for degree-five R because
    # its leading error uses R^(6)=0. Predict the discrete metric-route result
    # independently, rather than incorrectly comparing it to the continuum
    # oracle alone in the first cells where Ricci/r^2 amplifies this term.
    delta_prime = -4.0 * c5 * h**4
    tf_metric_discrete_expected = tf_ricci_exact + (2.0 / 3.0) * (
        2.0 * areal_prime_exact * delta_prime + delta_prime**2
    ) / areal_exact**2
    assert np.allclose(
        arrays["ricci_tracefree_metric"][centre_first_two],
        tf_metric_discrete_expected[centre_first_two],
        rtol=1.0e-8,
        atol=1.0e-10,
    )

    # The centre-scaled continuum oracle remains separately inspectable. Signed
    # centre values are kept as witnesses; the vendor route is not forced to
    # equal a continuum value when its discrete conversion is still under study.
    for key in ("ricci_tracefree_metric", "ricci_tracefree_vendor"):
        assert np.all(np.isfinite(arrays[key][centre_first_two])), f"{key}: nonfinite centre witness"
        assert np.all(np.isfinite(arrays[key][annulus])), f"{key}: nonfinite annulus witness"
        assert np.all(arrays[key][annulus] > 0.0), f"{key}: wrong annulus sign"
        assert np.all(arrays[key][centre_first_two] > 0.0), f"{key}: wrong centre sign"
        scaled = arrays[key] / r**2
        centre_error = scaled[centre_first_two] - tf_ricci_over_r2_exact[centre_first_two]
        assert np.all(np.isfinite(centre_error)), f"{key}: nonfinite c5 centre comparison"
        annulus_error = arrays[key][annulus] - tf_ricci_exact[annulus]
        assert np.max(np.abs(annulus_error)) < 2.0e-4, f"{key}: c5 annulus mismatch"

    assert summary["ricci_gap_max_abs_centre_i0_i1"] == pytest.approx(
        np.max(np.abs(arrays["ricci_tracefree_gap"][:2]))
    )

