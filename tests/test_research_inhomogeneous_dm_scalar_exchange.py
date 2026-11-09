"""Research test: DM/scalar exchange conserves total energy and momentum on an
INHOMOGENEOUS state, through the production kernel's own right-hand side.

Diagnostic-only. No production edits. scientific_admission is unchanged.

What it checks
--------------
The existing exchange tests are pointwise on hand-built arrays, or return a
pair that cancels by construction (tests/test_exchange.py). This test instead
runs V55ProductionKernel._rhs on a state that is non-uniform in phi, Pi, the
dark-matter density and velocity, and in the metric (lapse and a non-zero
shift), and isolates the exchange contribution by differencing the RHS at the
locked BETA_DM against the RHS at beta = 0. It then compares, cell by cell,

  dust:    d_t(sqrt(gamma) E)|exch, d_t(sqrt(gamma) S_r)|exch     (conserved vars)
  scalar:  sqrt(gamma) Pi dPi/dt|exch,  -sqrt(gamma) phi_r dPi/dt|exch
           (the Pi-dependent parts of scalar_projection's e and j)

with both divided by the same KAPPA, and requires dust + scalar = 0.

Sign convention under test (mostly plus; nabla_mu T_dm^{mu nu} = +Q^nu,
Q_nu = beta rho d_nu phi; nabla_mu T_phi^{mu nu} = -Q^nu). This is the "+"
branch of section 13A of the V5.5 gate document, the second of the two
statements the document makes.
"""
import dataclasses

import numpy as np
import pytest

import engine.production_kernel as pk
import engine.v55_pirk_adapter as adapter
from engine.matter_system import Species, _metric_arrays, initialize_dust
from engine.scalar_system import KAPPA
from engine.v55_matter import metric_slice_from_q


def _vendor_present() -> bool:
    return adapter.VENDOR.exists() and any(adapter.VENDOR.iterdir())


pytestmark = pytest.mark.skipif(
    not _vendor_present(), reason="pinned vendor kernel not initialised"
)


def _inhomogeneous_state(n=60, r_max=40.0):
    kernel = pk.V55ProductionKernel()
    st = kernel.initialize(resolution=n, r_max=r_max)
    r = np.asarray(st.grid.centers)
    # Smooth, regular-at-centre profiles (even in r), all non-uniform.
    bump = np.exp(-(r / 9.0) ** 2)
    st.scalars = dataclasses.replace(
        st.scalars,
        phi=(st.scalars.phi + 0.05 * bump * (1.0 + 0.6 * np.cos(r / 3.0))
             + 0.004 * r * np.sin(r / 2.5) * np.exp(-(r / 20.0) ** 2)),
        Pi=st.scalars.Pi + 0.004 * bump * (1.0 - 0.4 * np.sin(r / 7.0)),
    )
    st.geometry.alpha = st.geometry.alpha * (1.0 + 0.08 * bump)
    # Non-zero shift (odd in r) so the coordinate/normal conversion is live.
    st.geometry.beta = 0.03 * (r / 10.0) * np.exp(-(r / 12.0) ** 2)
    st.geometry.a = st.geometry.a * (1.0 + 0.05 * bump)
    metric = metric_slice_from_q(st.grid, st.geometry)
    rho = 3.0e-5 * (1.0 + 0.8 * bump + 0.2 * np.cos(r / 4.0))
    vr = 0.12 * np.tanh(r / 6.0) * np.exp(-(r / 15.0) ** 2)
    st.matter.dark_matter = initialize_dust(_metric_arrays(metric), rho, vr)
    return kernel, st, metric


def _exchange_parts(kernel, st, beta_dm):
    saved = pk.BETA_DM
    try:
        pk.BETA_DM = beta_dm
        srhs, mrhs, _ = kernel._rhs(st)
    finally:
        pk.BETA_DM = saved
    return srhs, mrhs["dark_matter"]


def _transfer(kernel, st, metric):
    srhs_b, dm_b = _exchange_parts(kernel, st, pk.BETA_DM)
    srhs_0, dm_0 = _exchange_parts(kernel, st, 0.0)
    dPi = srhs_b.Pi - srhs_0.Pi
    d_energy = dm_b.energy_t - dm_0.energy_t
    d_mom = dm_b.momentum_r - dm_0.momentum_r
    d_rest = dm_b.rest - dm_0.rest

    sqrtg = np.array([m.sqrt_gamma for m in _metric_arrays(metric)])
    phi_r = kernel._d1(st.grid, st.scalars.phi, 1)
    Pi = st.scalars.Pi
    scal_energy = sqrtg * Pi * dPi / KAPPA
    scal_mom = -sqrtg * phi_r * dPi / KAPPA
    return dict(
        dust_energy=d_energy / KAPPA,
        dust_mom=d_mom / KAPPA,
        dust_rest=d_rest,
        scal_energy=scal_energy,
        scal_mom=scal_mom,
        dPi=dPi,
    )


def test_exchange_is_equal_and_opposite_on_inhomogeneous_state():
    kernel, st, metric = _inhomogeneous_state()
    t = _transfer(kernel, st, metric)

    # Not degenerate: the transfer is non-zero, non-uniform, and changes sign.
    for key in ("dust_energy", "dust_mom", "scal_energy", "scal_mom"):
        assert np.max(np.abs(t[key])) > 1e-12, key
    assert np.ptp(t["dust_energy"]) > 0.1 * np.max(np.abs(t["dust_energy"]))
    assert (t["dust_mom"].min() < 0.0 < t["dust_mom"].max())

    scale_e = np.max(np.abs(t["dust_energy"]))
    scale_p = np.max(np.abs(t["dust_mom"]))
    np.testing.assert_allclose(
        t["dust_energy"] + t["scal_energy"], 0.0, atol=1e-12 * scale_e
    )
    np.testing.assert_allclose(
        t["dust_mom"] + t["scal_mom"], 0.0, atol=1e-12 * scale_p
    )


def test_exchange_negative_controls_are_observable():
    kernel, st, metric = _inhomogeneous_state()
    t = _transfer(kernel, st, metric)
    scale_e = np.max(np.abs(t["dust_energy"]))
    scale_p = np.max(np.abs(t["dust_mom"]))
    # Reversed scalar sign (the other branch of section 13A) must fail.
    assert np.max(np.abs(t["dust_energy"] - t["scal_energy"])) > 0.5 * scale_e
    assert np.max(np.abs(t["dust_mom"] - t["scal_mom"])) > 0.5 * scale_p
    # A missing lapse factor on the scalar source must fail.
    alpha = np.asarray(st.geometry.alpha)
    assert np.max(np.abs(
        t["dust_energy"] + t["scal_energy"] / alpha
    )) > 1e-3 * scale_e


def test_rest_mass_source_matches_variable_mass_law():
    """d_t(sqrt(gamma) D)|exch = -alpha sqrt(gamma) beta rho u.d phi."""
    kernel, st, metric = _inhomogeneous_state()
    t = _transfer(kernel, st, metric)
    metrics = _metric_arrays(metric)
    from engine.matter_system import primitives
    prim = primitives(metrics, st.matter.dark_matter, Species.DARK_MATTER)
    dphi_r = kernel._d1(st.grid, st.scalars.phi, 1)
    srhs0, _ = _exchange_parts(kernel, st, 0.0)
    dphi_t = srhs0.phi
    expected = np.empty(len(metrics))
    for i, (m, p) in enumerate(zip(metrics, prim)):
        W = p.lorentz()
        u_dphi = W / m.alpha * (dphi_t[i] + (m.alpha * p.v_r - m.beta) * dphi_r[i])
        expected[i] = -m.alpha * m.sqrt_gamma * pk.BETA_DM * p.rho * u_dphi
    np.testing.assert_allclose(
        t["dust_rest"], expected, rtol=1e-12, atol=1e-14 * np.max(np.abs(expected))
    )
