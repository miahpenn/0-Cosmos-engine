"""Research test: GEAR-91 coupling (research/gear91_coupling_v1.py). No production change.

(a) switches off -> bit-identical to production; (b) the S/D-phi exchange balances cell by cell through the kernel's
own RHS against the stress projection's potential, with negative controls (missing 8 pi, reversed sign); (c) the
universal switch keeps Omega_S^2 / Gamma_D^2 = 1/alpha0^2 for every phi; plus the potential's own identities.
"""
import dataclasses

import numpy as np
import pytest

import engine.production_kernel as pk
import engine.v55_pirk_adapter as adapter
import engine.v55_matter as vm
from engine.matter_system import _metric_arrays
from engine.v55_matter import metric_slice_from_q
from research import gear91_coupling_v1 as g91

KAPPA = g91.KAPPA


def _vendor_present() -> bool:
    return adapter.VENDOR.exists() and any(adapter.VENDOR.iterdir())


needs_vendor = pytest.mark.skipif(not _vendor_present(), reason="pinned vendor kernel not initialised")


# ---------------------------------------------------------------- the potential itself
def test_saddle_adds_nothing_for_any_phi_and_v_is_even_in_d():
    phi = np.linspace(-3, 3, 13)
    z = np.zeros_like(phi)
    for bS, bD in ((0.0, g91.BETA), (g91.BETA, g91.BETA), (0.3, -0.2)):
        assert np.all(g91.potential(z, z, phi, bS, bD, 1.0, 50.0) == 0.0)
        S, D = 0.1 + 0 * phi, 0.07 + 0 * phi
        np.testing.assert_array_equal(g91.potential(S, D, phi, bS, bD, 1.0, 50.0),
                                      g91.potential(S, -D, phi, bS, bD, 1.0, 50.0))


@pytest.mark.parametrize("lam", [0.0, 3.0])
def test_gradient_is_the_derivative_of_the_potential(lam):
    rng = np.random.default_rng(91)
    S, D, phi = rng.normal(size=30), rng.normal(size=30), rng.normal(size=30)
    args = (0.07, -0.04, 1.3, lam)
    dS, dD, dphi = g91.gradient(S, D, phi, *args)
    h = 1e-6
    V = lambda s, d, p: g91.potential(s, d, p, *args)  # noqa: E731
    np.testing.assert_allclose(dS, (V(S + h, D, phi) - V(S - h, D, phi)) / (2 * h), rtol=1e-7, atol=1e-8)
    np.testing.assert_allclose(dD, (V(S, D + h, phi) - V(S, D - h, phi)) / (2 * h), rtol=1e-7, atol=1e-8)
    np.testing.assert_allclose(dphi, (V(S, D, phi + h) - V(S, D, phi - h)) / (2 * h), rtol=1e-7, atol=1e-8)


def test_universal_switch_keeps_the_saddle_ratio_for_every_phi():
    phi = np.linspace(-5, 5, 11)
    np.testing.assert_allclose(g91.ratio_check(phi, alpha0=1.0), 1.0, rtol=1e-6)
    np.testing.assert_allclose(g91.ratio_check(phi, alpha0=0.8), 1 / 0.64, rtol=1e-6)


def test_g_kappa_log_with_s_uncoupled():
    phi = np.array([0.0, 0.179055, 1.0])
    g, kap, ratio = g91.g_kappa(phi, beta_D=g91.BETA, beta_S=0.0)
    KD = np.exp(-2 * g91.BETA * phi)
    np.testing.assert_allclose(g, (1 + KD) / 2)
    np.testing.assert_allclose(kap, (1 - KD) / 2)
    assert kap[0] == 0.0 and np.isinf(ratio[0])
    assert np.all(kap[1:] < 0.0)          # beta_D phi < 0 gives kappa < 0 (meaning not derived)


# ---------------------------------------------------------------- through the kernel
def _state(n=48, r_max=40.0):
    kernel = pk.V55ProductionKernel()
    st = kernel.initialize(resolution=n, r_max=r_max)
    r = np.asarray(st.grid.centers)
    bump = np.exp(-(r / 9.0) ** 2)
    sc = st.scalars
    st = dataclasses.replace(st, scalars=dataclasses.replace(
        sc,
        S=0.02 * np.exp(-(r / 6.0) ** 2) * (1 + 0.3 * np.cos(r / 3)),
        PS=0.004 * bump,
        D=0.015 * np.exp(-(r / 8.0) ** 2) * (1 - 0.2 * np.sin(r / 4)),
        PD=-0.003 * bump,
        phi=sc.phi + 0.3 * bump * (1.0 + 0.5 * np.cos(r / 3.0)),
        Pi=sc.Pi + 0.01 * bump * (1.0 - 0.4 * np.sin(r / 7.0)),
    ))
    st.geometry.alpha = st.geometry.alpha * (1.0 + 0.08 * bump)
    return kernel, st


@needs_vendor
def test_switches_off_is_bit_identical_to_production():
    kernel, st = _state()
    ref = kernel.step(st, 0.01)
    with g91.coupled(beta_S=0.0, beta_D=0.0):
        got = kernel.step(st, 0.01)
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        np.testing.assert_array_equal(getattr(got.scalars, name), getattr(ref.scalars, name))
    for name in ("a", "b", "X", "K", "Aa", "Lambda", "alpha"):
        np.testing.assert_array_equal(getattr(got.geometry, name), getattr(ref.geometry, name))


def _balance(kernel, st, **kw):
    """phi side from the kernel RHS (coupled minus coupled-without-phi-source); S/D side from the projection's
    potential (d e / d phi at fixed S, D, P, coupled minus uncoupled)."""
    with g91.coupled(**kw):
        s_on, _, _ = kernel._rhs(st)
    with g91.coupled(**{**kw, "phi_source": False}):
        s_off, _, _ = kernel._rhs(st)
    dPi = s_on.Pi - s_off.Pi
    h = 1e-6

    def e_of(dphi, on):
        f = dataclasses.replace(st.scalars, phi=st.scalars.phi + dphi)
        if on:
            with g91.coupled(beta_S=kw.get("beta_S", 0.0), beta_D=kw.get("beta_D", g91.BETA)):
                return vm.scalar_projection(st.grid, st.geometry, f)[0]
        return vm.scalar_projection(st.grid, st.geometry, f)[0]
    de_dphi = ((e_of(h, True) - e_of(-h, True)) - (e_of(h, False) - e_of(-h, False))) / (2 * h)
    metric = metric_slice_from_q(st.grid, st.geometry)
    sqrtg = np.array([m.sqrt_gamma for m in _metric_arrays(metric)])
    alpha = np.asarray(st.geometry.alpha)
    Pi = st.scalars.Pi
    phi_r = kernel._d1(st.grid, st.scalars.phi, 1)
    return dict(
        sd_energy=alpha * sqrtg * de_dphi * Pi,
        phi_energy=sqrtg * Pi * dPi / KAPPA,
        sd_mom=-alpha * sqrtg * de_dphi * phi_r,
        phi_mom=-sqrtg * phi_r * dPi / KAPPA,
    )


@needs_vendor
@pytest.mark.parametrize("bS,bD", [(0.0, g91.BETA), (g91.BETA, 0.0), (g91.BETA, g91.BETA)])
def test_exchange_balances_cell_by_cell(bS, bD):
    kernel, st = _state()
    t = _balance(kernel, st, beta_S=bS, beta_D=bD)
    for key in ("sd_energy", "phi_energy", "sd_mom", "phi_mom"):
        assert np.max(np.abs(t[key])) > 1e-14, key
    se, sp = np.max(np.abs(t["sd_energy"])), np.max(np.abs(t["sd_mom"]))
    np.testing.assert_allclose(t["sd_energy"] + t["phi_energy"], 0.0, atol=1e-6 * se)
    np.testing.assert_allclose(t["sd_mom"] + t["phi_mom"], 0.0, atol=1e-6 * sp)


@needs_vendor
def test_exchange_negative_controls_fail():
    kernel, st = _state()
    for bad in ({"kappa_factor": 1.0}, {"phi_source_sign": -1.0}):
        t = _balance(kernel, st, beta_S=0.0, beta_D=g91.BETA, **bad)
        se = np.max(np.abs(t["sd_energy"]))
        assert np.max(np.abs(t["sd_energy"] + t["phi_energy"])) > 0.5 * se, bad


@needs_vendor
def test_sd_rhs_uses_the_projected_potential():
    """dPS/dt and dPD/dt change by -alpha d(Delta V)/dS, -alpha d(Delta V)/dD of the projection's Delta V."""
    kernel, st = _state()
    s0, _, _ = kernel._rhs(st)
    with g91.coupled(beta_S=g91.BETA, beta_D=g91.BETA):
        s1, _, _ = kernel._rhs(st)
    h = 1e-7
    alpha = np.asarray(st.geometry.alpha)

    def dV(field, eps):
        f = dataclasses.replace(st.scalars, **{field: getattr(st.scalars, field) + eps})
        with g91.coupled(beta_S=g91.BETA, beta_D=g91.BETA):
            on = vm.scalar_projection(st.grid, st.geometry, f)[0]
        return on - vm.scalar_projection(st.grid, st.geometry, f)[0]
    for field, P in (("S", "PS"), ("D", "PD")):
        expect = -alpha * (dV(field, h) - dV(field, -h)) / (2 * h)
        got = getattr(s1, P) - getattr(s0, P)
        np.testing.assert_allclose(got, expect, rtol=1e-5, atol=1e-6 * np.max(np.abs(expect)))
