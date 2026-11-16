"""GEAR-67 S-D quartic: lambda > 0 is consistent everywhere.

Re-baselined for the hand-floor removal step 2 (Nov 16 2026): production now runs lambda = 1.57e6 (a ruler matched
today, engine/model_constants.py); lambda = 0 is the labelled hand-floor reference path. The uncoupled-quartic checks
below run on that reference path with lambda set by hand, as they did on the old production path; the coupled
production path is tested in tests/test_hand_floor_removal_step2.py."""
import math

import numpy as np
import pytest

import engine.scalar_system as ss


def test_production_lambda_is_the_today_matched_ruler_and_reference_is_zero():
    from engine import model_constants as mc

    assert ss.SD_QUARTIC_LAMBDA == mc.SD_QUARTIC_LAMBDA
    assert ss.SD_QUARTIC_LAMBDA == pytest.approx(1.56889e6, rel=1e-5)
    with ss.hand_floor_reference():
        assert ss.SD_QUARTIC_LAMBDA == 0.0
    assert ss.SD_QUARTIC_LAMBDA == mc.SD_QUARTIC_LAMBDA


def test_lambda_zero_gradient_is_exactly_the_quadratic_force():
    rng = np.random.default_rng(1)
    S, D = rng.normal(size=50), rng.normal(size=50)
    dS, dD = ss.sd_potential_gradient(S, D, 0.0)
    np.testing.assert_array_equal(dS, S)
    np.testing.assert_array_equal(dD, -D)


@pytest.mark.parametrize("lam", [0.1, 1.0, 7.5])
def test_gradient_is_the_derivative_of_the_gear67_potential(lam):
    rng = np.random.default_rng(2)
    S, D = rng.normal(size=40), rng.normal(size=40)
    dS, dD = ss.sd_potential_gradient(S, D, lam)
    h = 1e-6
    fdS = (ss.sd_potential(S + h, D, lam) - ss.sd_potential(S - h, D, lam)) / (2 * h)
    fdD = (ss.sd_potential(S, D + h, lam) - ss.sd_potential(S, D - h, lam)) / (2 * h)
    np.testing.assert_allclose(dS, fdS, rtol=1e-7, atol=1e-8)
    np.testing.assert_allclose(dD, fdD, rtol=1e-7, atol=1e-8)


def test_gear67_reduction_from_two_elements():
    """V_SD equals two elements with on-site lambda x^4/4 and bond -x1 x2 (g = 1, kappa = 0)."""
    rng = np.random.default_rng(3)
    lam = 0.7
    S, D = rng.normal(size=30), rng.normal(size=30)
    x1, x2 = (D + S) / math.sqrt(2.0), (D - S) / math.sqrt(2.0)   # D ejected (sum), S surviving (difference)
    two_element = -x1 * x2 + (lam / 4.0) * (x1**4 + x2**4)
    np.testing.assert_allclose(ss.sd_potential(S, D, lam), two_element, rtol=1e-12, atol=1e-12)


def test_force_on_D_holds_it_where_matter_clears_the_threshold():
    """D-dot-dot = (1 - 3 lambda S^2 / 2) D - (lambda/2) D^3: held iff S^2 > 2/(3 lambda) = D*^2 / 3."""
    lam = 0.5
    D = np.array([1e-6])
    threshold = 2.0 / (3.0 * lam)
    for S2, held in ((0.5 * threshold, False), (2.0 * threshold, True)):
        _, dD = ss.sd_potential_gradient(np.sqrt(np.array([S2])), D, lam)
        assert bool(-dD[0] / D[0] < 0.0) == held
    assert threshold == pytest.approx((2.0 / lam) / 3.0)


def test_engine_rhs_and_stress_use_the_quartic_when_lambda_is_set():
    from engine.production_kernel import V55ProductionKernel

    with ss.hand_floor_reference():           # uncoupled path; lambda set by hand below, restored on exit
        kernel = V55ProductionKernel()
        state = kernel.initialize(resolution=24, r_max=40.0)
        grid, geometry = state.grid, state.geometry
        r = np.asarray(grid.centers)
        import dataclasses
        fields = dataclasses.replace(
            state.scalars, S=0.3 * np.exp(-(r / 6) ** 2), D=0.2 * np.exp(-(r / 5) ** 2)
        )
        base_rhs = ss.scalar_rhs_arrays(grid, geometry, fields)
        base_e = ss.scalar_projection(grid, geometry, fields)[0]
        ss.SD_QUARTIC_LAMBDA = 2.0
        rhs = ss.scalar_rhs_arrays(grid, geometry, fields)
        e = ss.scalar_projection(grid, geometry, fields)[0]
    dS, dD = ss.sd_potential_gradient(fields.S, fields.D, 2.0)
    np.testing.assert_allclose(rhs.PS - base_rhs.PS, -geometry.alpha * (dS - fields.S), rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(rhs.PD - base_rhs.PD, -geometry.alpha * (dD + fields.D), rtol=1e-12, atol=1e-15)
    quartic = 0.25 * (fields.S**4 + 6 * fields.S**2 * fields.D**2 + fields.D**4)
    np.testing.assert_allclose(e - base_e, quartic, rtol=1e-12, atol=1e-15)


def test_initial_hamiltonian_constraint_still_converges_with_lambda():
    from engine import v55_pirk_adapter as adapter
    from engine.production_kernel import V55ProductionKernel

    errs = []
    for n in (60, 120):
        kernel = V55ProductionKernel()
        with ss.hand_floor_reference():       # uncoupled quartic on the reference pocket, lambda set by hand
            ss.SD_QUARTIC_LAMBDA = 50.0
            state = kernel.initialize(resolution=n, r_max=40.0, amplitude=0.02, D_amplitude=0.02)
            _, vac, _ = adapter.vendor_modules()
            c = vac.constraints(state.grid, state.geometry)
            p = kernel.diagnostics(state, profiles=True)
        p = p.get("profiles", p)
        H = np.asarray(c["hamiltonian"]) - 16 * math.pi * np.asarray(p["total_rho"])
        r = np.asarray(state.grid.centers)
        errs.append(float(np.max(np.abs(H[(r >= 2) & (r <= 15)]))))
    assert errs[0] / errs[1] >= 8.0


def test_reduced_local_lane_refuses_a_nonzero_lambda():
    from engine.local import LocalParams, local_rhs

    y = (0.1, 0.0, 1e-3, 0.0)
    with ss.hand_floor_reference():
        assert len(local_rhs(y, LocalParams())) == 4
    # production (lambda = 1.57e6 since step 2) is refused by the provenance lane
    with pytest.raises(NotImplementedError, match="lambda = 0"):
        local_rhs(y, LocalParams())
