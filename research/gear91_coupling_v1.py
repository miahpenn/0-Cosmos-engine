"""RESEARCH ONLY. GEAR-91 coupling of the S-D pair to phi (HANDOFF-to-projects Nov 10, Task 2).

    V_SD(S, D; phi) = e^(-2 beta_S phi) S^2/2 - e^(-2 beta_D phi) alpha0^2 D^2/2 + (lambda/8)(S^4 + 6 S^2 D^2 + D^4)

GEAR-91, AMENDMENT-GEAR-91 (saddle vacuum 0*, supplied by phi; phi_ref = 0 shared with dark matter's GEAR-78 convention;
the 8 pi in the exchange), AMENDMENT-2-GEAR-91, AMENDMENT-3-GEAR-91 (separate switches; S uncoupled by default).

- Separate switches: beta_S (default 0: S is tangible matter; a coupled S meets the Cassini bound ~100x over) and
  beta_D (default beta = -0.04, the locked dark-matter coupling). Universal coupling is beta_S = beta_D = beta, a
  switch, not a default; the physics choice is Jeremiah's.
- One action. S, D are normalised at G_N = 1, phi at 8 pi G_N = 1 (its stress is divided by KAPPA = 8 pi), so
      L = (1/8pi)[-(d phi)^2/2 - V_phi] - (dS)^2/2 - (dD)^2/2 - V_SD(S, D; phi)
  and phi's equation gains 8 pi dV_SD/dphi:  Box phi = V_phi' + 8 pi dV_SD/dphi.  S and D see the phi-dependent masses
  through dV_SD/dS, dV_SD/dD. The stress of S/D carries the full V_SD.
- S/D add nothing at the saddle for any phi (V_SD(0, 0; phi) = 0): the saddle's vacuum is 0*, supplied by phi.

No production file changes. The coupling is applied by `coupled(...)`, a context manager that wraps the two
production entry points the kernel uses (engine.production_kernel.scalar_rhs_arrays and
engine.v55_matter.scalar_projection) and adds the difference between the coupled and the production potential.
With beta_S = beta_D = 0 nothing is added (bit-identical to production at the same lambda).

Unpaid inputs to name with every result: lambda (GEAR-67 3(e)); V0 and phi_today (GEAR-92) whenever phi is on
("at production V0"); beta_S / beta_D choice.
"""
from __future__ import annotations

import contextlib

import numpy as np

import engine.production_kernel as _pk
import engine.scalar_system as _ss
import engine.v55_matter as _vm

KAPPA = _ss.KAPPA
BETA = -0.04          # locked dark-matter coupling (GEAR-78)

_CFG = {"beta_S": 0.0, "beta_D": BETA, "alpha0": 1.0, "phi_source": True, "kappa_factor": KAPPA,
        "phi_source_sign": 1.0}


def _lam():
    return _ss.SD_QUARTIC_LAMBDA


def potential(S, D, phi, beta_S=None, beta_D=None, alpha0=None, lam=None):
    bS = _CFG["beta_S"] if beta_S is None else beta_S
    bD = _CFG["beta_D"] if beta_D is None else beta_D
    a0 = _CFG["alpha0"] if alpha0 is None else alpha0
    lam = _lam() if lam is None else lam
    return (np.exp(-2 * bS * phi) * 0.5 * S**2 - np.exp(-2 * bD * phi) * 0.5 * a0**2 * D**2
            + (lam / 8.0) * (S**4 + 6.0 * S**2 * D**2 + D**4))


def gradient(S, D, phi, beta_S=None, beta_D=None, alpha0=None, lam=None):
    """(dV/dS, dV/dD, dV/dphi)."""
    bS = _CFG["beta_S"] if beta_S is None else beta_S
    bD = _CFG["beta_D"] if beta_D is None else beta_D
    a0 = _CFG["alpha0"] if alpha0 is None else alpha0
    lam = _lam() if lam is None else lam
    eS, eD = np.exp(-2 * bS * phi), np.exp(-2 * bD * phi)
    dS = eS * S + (lam / 2.0) * S**3 + 1.5 * lam * S * D**2
    dD = -eD * a0**2 * D + (lam / 2.0) * D**3 + 1.5 * lam * S**2 * D
    dphi = -bS * eS * S**2 + bD * eD * a0**2 * D**2
    return dS, dD, dphi


def g_kappa(phi, beta_D=None, beta_S=None):
    """GEAR-03 element couplings read off the quadratic part (log only; meaning not derived, AMENDMENT-3-GEAR-91 §2):
    K_S = e^(-2 beta_S phi), K_D = e^(-2 beta_D phi); g = (K_S + K_D)/2, kappa = (K_S - K_D)/2, and g/kappa."""
    bS = _CFG["beta_S"] if beta_S is None else beta_S
    bD = _CFG["beta_D"] if beta_D is None else beta_D
    KS, KD = np.exp(-2 * bS * np.asarray(phi)), np.exp(-2 * bD * np.asarray(phi))
    g, kap = 0.5 * (KS + KD), 0.5 * (KS - KD)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(kap != 0.0, g / kap, np.inf)
    return g, kap, ratio


def _active():
    return _CFG["beta_S"] != 0.0 or _CFG["beta_D"] != 0.0 or _CFG["alpha0"] != 1.0


def _base_gradient(S, D):
    return _ss.sd_potential_gradient(S, D)


def make_rhs(orig):
    def scalar_rhs_arrays(grid, geometry, fields, beta_dm=-0.04, rho_dm=None):
        out = orig(grid, geometry, fields, beta_dm=beta_dm, rho_dm=rho_dm)
        if not _active():
            return out
        bS0, bD0 = _base_gradient(fields.S, fields.D)
        dS, dD, dphi = gradient(fields.S, fields.D, fields.phi)
        a = geometry.alpha
        PS = out.PS - a * (dS - bS0)
        PD = out.PD - a * (dD - bD0)
        Pi = out.Pi
        if _CFG["phi_source"]:
            Pi = Pi - _CFG["phi_source_sign"] * a * _CFG["kappa_factor"] * dphi
        return _ss.ScalarFields(out.S, PS, out.D, PD, out.phi, Pi)
    return scalar_rhs_arrays


def make_projection(orig):
    def scalar_projection(grid, geometry, fields):
        e, pr, pt, j = orig(grid, geometry, fields)
        if not _active():
            return e, pr, pt, j
        dV = potential(fields.S, fields.D, fields.phi) - _ss.sd_potential(fields.S, fields.D)
        return e + dV, pr - dV, pt - dV, j
    return scalar_projection


@contextlib.contextmanager
def coupled(beta_S=0.0, beta_D=BETA, alpha0=1.0, phi_source=True, kappa_factor=KAPPA, phi_source_sign=1.0):
    """Apply the GEAR-91 coupling to the production kernel for the duration of the block.

    kappa_factor and phi_source_sign exist for negative controls only (a missing 8 pi, a reversed sign)."""
    saved_cfg = dict(_CFG)
    rhs0, proj0 = _pk.scalar_rhs_arrays, _vm.scalar_projection
    _CFG.update(beta_S=float(beta_S), beta_D=float(beta_D), alpha0=float(alpha0), phi_source=bool(phi_source),
                kappa_factor=float(kappa_factor), phi_source_sign=float(phi_source_sign))
    _pk.scalar_rhs_arrays = make_rhs(rhs0)
    _vm.scalar_projection = make_projection(proj0)
    try:
        yield _CFG
    finally:
        _pk.scalar_rhs_arrays, _vm.scalar_projection = rhs0, proj0
        _CFG.clear(); _CFG.update(saved_cfg)


def ratio_check(phi, alpha0=1.0, beta=BETA):
    """Universal switch: Omega_S^2 / Gamma_D^2 at the saddle = 1/alpha0^2 for every phi."""
    eps = 1e-5
    S0 = np.zeros_like(np.asarray(phi, dtype=float))
    d2S = (potential(S0 + eps, S0, phi, beta, beta, alpha0, 0.0) - 2 * potential(S0, S0, phi, beta, beta, alpha0, 0.0)
           + potential(S0 - eps, S0, phi, beta, beta, alpha0, 0.0)) / eps**2
    d2D = (potential(S0, S0 + eps, phi, beta, beta, alpha0, 0.0) - 2 * potential(S0, S0, phi, beta, beta, alpha0, 0.0)
           + potential(S0, S0 - eps, phi, beta, beta, alpha0, 0.0)) / eps**2
    return d2S / (-d2D)


__all__ = ["potential", "gradient", "g_kappa", "coupled", "ratio_check", "BETA", "KAPPA"]
