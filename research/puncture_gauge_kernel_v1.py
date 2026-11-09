"""RESEARCH ONLY. V5.5 kernel with the damped moving-puncture gauge instead of CMC.

``PunctureGaugeKernel`` subclasses ``V55ProductionKernel`` and replaces only its time step. The lapse
and shift are evolved (advective 1+log lapse, beta_t = B, damped Gamma-driver B_t = 3/4 Lambda_t
- eta B), following the vendor's ``damped_moving_puncture_pirk2_step`` stage for stage, while the
scalar fields, the fluids and their stress-energy are the production kernel's own (``_rhs``,
``adapter.geometry_stage_terms`` with matter sources in the Aa, K and Lambda blocks). No production
file is changed and CMC is not involved.

Stage structure (PIRK2, as the vendor step):
  stage 1: u = (a, b, X, alpha, beta) explicit Euler; (Aa, K) with 1/2 L2(u0) + 1/2 L2(u1) + L3(u0);
           Lambda with 1/2 LL2(u0) + 1/2 LL2(stage1) + LL3(u0); B from the Gamma-driver;
           scalars and fluids explicit Euler.
  final:   u trapezoid; (Aa, K) with 1/2 [L2(u0) + L2(u_final) + L3(u0) + L3(stage1)];
           Lambda likewise; B trapezoid; scalars and fluids trapezoid.
Optional Kreiss-Oliger dissipation (vendor filter) on the geometry, as in the vendor's puncture runs.
No outer-boundary reconstruction is applied (the production one is written for CMC); use a large
r_max so the boundary stays causally disconnected over the run.

Not covered: the radial Maxwell sector (electromagnetic state must be None).
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

import engine.production_kernel as pk
from engine import v55_pirk_adapter as adapter
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState
from engine.matter_system import ConservedSpecies

U_NAMES = ("a", "b", "X", "alpha", "beta")
V_NAMES = ("Aa", "K")
SPECIES = ("dark_matter", "baryons", "radiation")


class PunctureGaugeKernel(pk.V55ProductionKernel):
    def __init__(self, *, eta: float = 2.0, ko_epsilon: float = 0.0, **kwargs):
        super().__init__(**kwargs)
        if not (np.isfinite(eta) and eta >= 0.0):
            raise ValueError("eta must be finite and nonnegative")
        if not (np.isfinite(ko_epsilon) and ko_epsilon >= 0.0):
            raise ValueError("ko_epsilon must be finite and nonnegative")
        self.eta = float(eta)
        self.ko_epsilon = float(ko_epsilon)

    # ------------------------------------------------------------------ stage terms
    def _terms(self, grid, geometry, scalars, matter):
        _, vacuum, moving = adapter.vendor_modules()
        explicit = moving.moving_puncture_explicit_rhs(grid, geometry, ko_epsilon=self.ko_epsilon)
        t = adapter.geometry_stage_terms(
            grid, geometry, scalars, matter, lambda_m=pk.LAMBDA_M, explicit_rhs_override=explicit
        )
        if self.ko_epsilon:
            ko = vacuum.ko_dissipation
            t = dict(t)
            t["primary_l3"] = {k: np.asarray(v) + ko(grid, getattr(geometry, k), 1, epsilon=self.ko_epsilon)
                               if k in V_NAMES else v for k, v in t["primary_l3"].items()}
            t["lambda_l3"] = np.asarray(t["lambda_l3"]) + ko(grid, geometry.Lambda, -1, epsilon=self.ko_epsilon)
        return t

    def _l2(self, grid, geometry):
        _, vacuum, _ = adapter.vendor_modules()
        return vacuum.primary_l2_rhs(grid, geometry)

    def _ll2(self, grid, geometry):
        _, vacuum, _ = adapter.vendor_modules()
        return vacuum.lambda_l2_rhs(grid, geometry, lambda_m=pk.LAMBDA_M)

    @staticmethod
    def _matter_combo(m0, rhs_a, rhs_b, dt, w_a, w_b):
        out = []
        for name in SPECIES:
            s0 = getattr(m0, name)
            ra, rb = rhs_a[name], rhs_b[name] if rhs_b is not None else None

            def comb(field):
                v = getattr(s0, field) + dt * w_a * getattr(ra, field)
                if rb is not None:
                    v = v + dt * w_b * getattr(rb, field)
                return v
            out.append(ConservedSpecies(comb("rest"), comb("energy_t"), comb("momentum_r")))
        return V55MatterState(*out, charged_dust=None)

    # ------------------------------------------------------------------ step
    def step(self, state: pk.ProductionState, dt: float) -> pk.ProductionState:
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt must be finite and positive")
        if state.electromagnetic is not None or state.matter.charged_dust is not None:
            raise NotImplementedError("the puncture-gauge research kernel has no radial Maxwell sector")
        grid = state.grid
        g0, s0, m0 = state.geometry.copy(), state.scalars, state.matter
        st0 = replace(state, geometry=g0)

        t0 = self._terms(grid, g0, s0, m0)
        srhs0, mrhs0, _ = self._rhs(st0)

        # ---- stage 1
        g_u1 = g0.copy()
        for name in U_NAMES:
            setattr(g_u1, name, getattr(g0, name) + dt * t0["explicit"][name])
        s1 = ScalarFields(*(getattr(s0, n) + dt * getattr(srhs0, n) for n in ("S", "PS", "D", "PD", "phi", "Pi")))
        m1 = self._matter_combo(m0, mrhs0, None, dt, 1.0, 0.0)
        l2_u1 = self._l2(grid, g_u1)
        g1 = g_u1.copy()
        for name in V_NAMES:
            setattr(g1, name, getattr(g0, name) + dt * (
                0.5 * t0["primary_l2"][name] + 0.5 * l2_u1[name] + t0["primary_l3"][name]))
        g1.Lambda = g0.Lambda + dt * (0.5 * t0["lambda_l2"] + 0.5 * self._ll2(grid, g1) + t0["lambda_l3"])
        g1.B = g0.B + 0.75 * (g1.Lambda - g0.Lambda) - dt * self.eta * g0.B
        g1.assert_finite_positive()
        st1 = replace(state, geometry=g1, scalars=s1, matter=m1)
        t1 = self._terms(grid, g1, s1, m1)
        srhs1, mrhs1, _ = self._rhs(st1)

        # ---- final
        gnew = g1.copy()
        for name in U_NAMES:
            setattr(gnew, name, getattr(g0, name) + 0.5 * dt * (t0["explicit"][name] + t1["explicit"][name]))
        l2_final = self._l2(grid, gnew)
        for name in V_NAMES:
            setattr(gnew, name, getattr(g0, name) + 0.5 * dt * (
                t0["primary_l2"][name] + l2_final[name] + t0["primary_l3"][name] + t1["primary_l3"][name]))
        gnew.Lambda = g0.Lambda + 0.5 * dt * (
            t0["lambda_l2"] + self._ll2(grid, gnew) + t0["lambda_l3"] + t1["lambda_l3"])
        gnew.B = g0.B + 0.75 * (gnew.Lambda - g0.Lambda) - 0.5 * dt * self.eta * (g0.B + g1.B)
        gnew.assert_finite_positive()
        snew = ScalarFields(*(getattr(s0, n) + 0.5 * dt * (getattr(srhs0, n) + getattr(srhs1, n))
                              for n in ("S", "PS", "D", "PD", "phi", "Pi")))
        mnew = self._matter_combo(m0, mrhs0, mrhs1, dt, 0.5, 0.5)

        tau = state.tau + 0.5 * dt * (float(g0.alpha[0]) + float(gnew.alpha[0]))
        return pk.ProductionState(
            grid=grid, geometry=gnew, scalars=snew, matter=mnew, electromagnetic=None,
            t=state.t + dt, tau=tau, e_folds=state.e_folds,
            cycle=state.cycle, history=list(state.history), handoffs=list(state.handoffs),
        )
