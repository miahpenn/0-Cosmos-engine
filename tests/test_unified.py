from math import isclose, pi

from engine.unified import (
    UnifiedUnits,
    ScalarStress,
    scalar_stress,
    sum_stress,
    momentum_curvature,
    radial_B_rhs,
    misner_sharp_mass,
    trapping_indicator,
    misner_sharp_current_residual,
)


def test_archive_normalization_bridge():
    u = UnifiedUnits()
    assert isclose(u.kappa, 8.0 * pi)
    assert isclose(u.rho_phys(8.0 * pi), 1.0)


def test_scalar_stress_and_sum():
    s = scalar_stress(phi_r=2.0, Pi=3.0, potential=0.5, a=1.0)
    assert s == ScalarStress(rho=7.0, p_r=6.0, p_t=0.0, j_r=-6.0)
    total = sum_stress(s, s)
    assert total.rho == 14.0
    assert total.j_r == -12.0


def test_constraints_and_mass_relations():
    Kr = momentum_curvature(0.2, 0.1, 2.0, 0.3)
    assert isclose(Kr, 0.2 + 2.0 * 0.1 + 4.0 * pi * 2.0 * 0.3)

    rhs = radial_B_rhs(2.0, 0.01, Kr, 0.2)
    assert rhs == radial_B_rhs(2.0, 0.01, Kr, 0.2)

    assert misner_sharp_mass(2.0, 0.9, 0.1) == 0.11
    assert trapping_indicator(0.9, 0.1) == 0.89


def test_current_gate_zero_for_matching_data():
    flux = 0.7
    r = 3.0
    dm_dt = 4.0 * pi * r * r * flux
    assert isclose(misner_sharp_current_residual(dm_dt, r, flux), 0.0)
