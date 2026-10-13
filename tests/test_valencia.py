from math import isclose

from engine.valencia import (
    FluidPrimitive,
    SphericalMetric,
    mixed_flux,
    primitive_to_conserved,
    spherical_metric_from_bssn,
)


def flat_metric():
    return SphericalMetric(
        alpha=1.0,
        beta=0.0,
        gamma_rr=1.0,
        gamma_rr_inv=1.0,
        gamma_thth=4.0,
        gamma_thth_inv=0.25,
        sqrt_gamma=2.0,
    )


def test_dust_rest_state():
    m = flat_metric()
    q = FluidPrimitive(rho=3.0, pressure=0.0, v_r=0.0, gamma_rr=1.0)
    u = primitive_to_conserved(m, q)
    assert isclose(u.rest, 6.0)
    assert isclose(u.momentum_r, 0.0)
    assert isclose(u.energy_t, -6.0)


def test_radiation_flux_in_flat_limit():
    m = flat_metric()
    q = FluidPrimitive(rho=3.0, pressure=1.0, v_r=0.0, gamma_rr=1.0)
    u = mixed_flux(m, q)
    assert u.rest == 0.0
    assert isclose(u.energy_t, 0.0)
    assert isclose(u.momentum_r, 2.0)


def test_spherical_metric_mapping():
    m = spherical_metric_from_bssn(
        r=2.0, a=1.0, b=1.0, X=1.0, alpha=1.0, beta=0.0
    )
    assert isclose(m.gamma_rr, 1.0)
    assert isclose(m.gamma_thth, 4.0)
    assert isclose(m.sqrt_gamma, 4.0)
