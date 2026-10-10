import math
from math import isclose

import numpy as np

from engine.matter_system import (
    BSSNMetricSlice,
    MetricDerivativeSet,
    Species,
    evolve_species,
    initialize_radiation,
)
from engine.valencia import (
    FluidConserved,
    FluidPrimitive,
    SphericalMetric,
    mixed_flux,
    primitive_to_conserved,
    recover_radiation,
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


def test_radiation_reconstruction_preserves_equation_of_state():
    from engine.matter_system import _reconstructed_primitive

    prim = [
        FluidPrimitive(rho=3.0, pressure=1.0, v_r=0.10, gamma_rr=1.0),
        FluidPrimitive(rho=6.0, pressure=2.0, v_r=0.20, gamma_rr=1.0),
        FluidPrimitive(rho=9.0, pressure=3.0, v_r=0.30, gamma_rr=1.0),
    ]

    left = _reconstructed_primitive(
        prim, 1, "left", 1.0, Species.RADIATION
    )
    right = _reconstructed_primitive(
        prim, 1, "right", 1.0, Species.RADIATION
    )

    assert isclose(left.pressure, left.rho / 3.0)
    assert isclose(right.pressure, right.rho / 3.0)


def test_radiation_transport_preserves_zero_rest_component():
    n = 4
    r = np.arange(n, dtype=float) + 0.5
    zeros = np.zeros(n, dtype=float)
    metric = BSSNMetricSlice(
        r=r,
        a=np.ones(n),
        b=np.ones(n),
        X=np.ones(n),
        alpha=np.ones(n),
        beta=zeros,
        Aa=zeros,
        K=zeros,
        Lambda=zeros,
        B=zeros,
    )
    derivatives = MetricDerivativeSet(
        time={key: zeros.copy() for key in ("tt", "tr", "rr", "thth")},
        radial={
            key: zeros.copy()
            for key in ("tt", "tr", "rr", "thth", "alpha", "beta")
        },
    )
    metrics = [
        spherical_metric_from_bssn(
            float(ri), 1.0, 1.0, 1.0, 1.0, 0.0
        )
        for ri in r
    ]
    state = initialize_radiation(
        metrics, np.full(n, 1.0e-6), np.zeros(n)
    )

    advanced = evolve_species(
        metric, derivatives, state, Species.RADIATION, 1.0e-3
    )

    assert np.allclose(advanced.rest, 0.0, atol=1.0e-30)



def test_low_density_near_null_radiation_recovery_is_relative_accuracy():
    """Tiny physical energy must not turn the inversion tolerance into an absolute 1e-12."""
    metric = SphericalMetric(
        alpha=1.0,
        beta=0.0,
        gamma_rr=0.6938513563030541,
        gamma_rr_inv=1.4412308787982382,
        gamma_thth=5000.0,
        gamma_thth_inv=1.0 / 5000.0,
        sqrt_gamma=7579.229307576932,
    )
    energy = 3.5729677738253134e-8
    ratio = 0.999571710404479
    radial_momentum = ratio * energy / math.sqrt(metric.gamma_rr_inv)
    state = FluidConserved(
        rest=0.0,
        energy_t=metric.sqrt_gamma * energy,
        momentum_r=metric.sqrt_gamma * radial_momentum,
    )

    recovered = recover_radiation(metric, state)
    W = recovered.lorentz()
    h = recovered.rho + recovered.pressure
    energy_reconstructed = h * W * W - recovered.pressure
    momentum_reconstructed = h * W * W * metric.gamma_rr * recovered.v_r

    assert recovered.pressure > 0.0
    assert abs(energy_reconstructed - energy) / energy < 1.0e-10
    assert abs(momentum_reconstructed - radial_momentum) / abs(radial_momentum) < 1.0e-10
