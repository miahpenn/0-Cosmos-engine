"""Research-only: signed-H COSMOS cycle through the full contraction, in e-fold time.

Brought into the repo from the Drive archive module ``cosmos_efold.py`` (Oct 7, 2026),
with one change: the end of the contraction is a DENSITY condition, not an e-fold count.

Integration (unchanged from the archive module)
----------------------------------------------
The homogeneous COSMOS lane (``engine.cosmos``) is a signed-H dynamical system: H is an
evolved state, so the turnaround is the equations' own H = 0 crossing. In cosmic time the
contracting branch collapses in finite time and the integrator stalls about 10 e-folds below
the turnaround, so:

    Stage 1 (cosmic time): today -> turnaround -> ``handoff`` e-folds below it, using
        ``engine.cosmos.cosmos_rhs`` unchanged.
    Stage 2 (e-fold time): N = ln(a / a_turn) is the independent variable; the SAME
        equations divided by H. H stays evolved (as ln|H|, sign inherited), so the
        Friedmann residual is a witness, not enforced. Densities are carried as logs.

No bounce, reset, branch flip, stop law, fitted coefficient or clock conversion is added.

End of contraction (the change)
-------------------------------
GEAR 57's depth (47.203 = 123 ln10 / 6) and CORRECTION-GEAR-50's audit (47.12 = ln(ratio)/6)
both convert one physical condition, "the total density reaches the Planck density", into an
e-fold count under w = 1 scaling from TODAY's density. The integrated contraction starts at
the turnaround, below today's density, and w briefly exceeds 1 after it, so neither count is
the depth of this trajectory. Here the observer stops where

    rho_total / rho_today = PLANCK_TO_CRITICAL_DENSITY_RATIO

(rho_today = 3 H_today^2 in the lane's units, i.e. today's critical density), and the depth
is REPORTED: e-folds from the turnaround, e-folds below today, and the phi span. The ratio
is CORRECTION-GEAR-50's audited value, recomputed below from the same inputs.

``CosmosParams.n_contract`` is not read here. An explicit ``depth`` (e-folds below the
turnaround) is still accepted as an observer-only stopping point for comparisons.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import log, pi as PI

import numpy as np
from scipy.integrate import solve_ivp

from engine.cosmos import (
    CosmosParams,
    construct_present_day,
    cosmos_rhs,
    cosmos_rho_p,
)

# CORRECTION-GEAR-50 (Sep 6, 2026): rho_Planck / rho_crit with H_0 = 2.20e-18 s^-1 (h = 0.678).
_C = 299_792_458.0            # m s^-1
_HBAR = 1.054_571_817e-34     # J s
_G_SI = 6.674_30e-11          # m^3 kg^-1 s^-2
_H0_SI = 2.20e-18             # s^-1, the value the audit used
RHO_PLANCK_SI = _C**5 / (_HBAR * _G_SI**2)
RHO_CRIT_TODAY_SI = 3.0 * _H0_SI**2 / (8.0 * PI * _G_SI)
PLANCK_TO_CRITICAL_DENSITY_RATIO = RHO_PLANCK_SI / RHO_CRIT_TODAY_SI   # 5.954e122
# The Planck scale is defined only up to factors like (8 pi)^2: the reduced Planck density
# c^5 / (hbar (8 pi G)^2) is 632x lower (10^119.97). The ratio rho_total/rho_today is the same
# in any unit system, so the lane's 8 pi G = 1 units do not choose between them; the choice
# moves the endpoint by ln(64 pi^2)/6 = 1.07 e-folds at w = 1. Both are reported.
REDUCED_PLANCK_TO_CRITICAL_DENSITY_RATIO = PLANCK_TO_CRITICAL_DENSITY_RATIO / (8.0 * PI) ** 2

# Power-law exponents: rho_i a^n_i = const for sources with no exchange.
_POWER = {"rho_b": 3.0, "rho_r": 4.0, "rho_shear": 6.0}


@dataclass
class CycleTrace:
    t: np.ndarray            # cosmic time since the present-day state
    ln_a_today: np.ndarray   # ln(a / a_today)
    N_turn: np.ndarray       # ln(a / a_turn) (0 at the turnaround, < 0 after)
    H: np.ndarray
    phi: np.ndarray
    pi_phi: np.ndarray
    rho_dm: np.ndarray
    rho_b: np.ndarray
    rho_r: np.ndarray
    rho_shear: np.ndarray
    stage: np.ndarray        # 1 = cosmic-time stage, 2 = e-fold-time stage
    t_turn: float
    a_turn_over_a_today: float
    rho_today: float
    endpoint: str            # "planck_density" or "depth"
    status: str


def _rho_p(H, phi, pi, dm, b, r, s, p):
    return cosmos_rho_p((1.0, H, phi, pi, dm, b, r, s), p)


def integrate_cycle(
    p: CosmosParams = CosmosParams(),
    depth: float | None = None,
    density_ratio: float = PLANCK_TO_CRITICAL_DENSITY_RATIO,
    handoff: float = 1.0,
    rtol1: float = 1e-12,
    atol1: float = 1e-18,
    rtol2: float = 1e-11,
    t_max: float = 5000.0,
    n_limit: float = 200.0,
) -> CycleTrace:
    """Integrate today -> turnaround -> the end of contraction.

    By default the observer stops where rho_total / rho_today = ``density_ratio``. With
    ``depth`` it stops ``depth`` e-folds below the turnaround instead.
    """
    s0 = construct_present_day(p)
    y0 = np.array([s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear])
    rho_today = 3.0 * s0.H**2

    def f(t, y):
        return np.array(cosmos_rhs(tuple(y), p))

    # ---- Stage 1a: today -> turnaround (first downward H crossing) ----
    def ev_turn(t, y):
        return y[1]
    ev_turn.terminal = True
    ev_turn.direction = -1
    A = solve_ivp(f, (0.0, t_max), y0, events=ev_turn, rtol=rtol1, atol=atol1, method="DOP853")
    if A.status != 1:
        raise RuntimeError(f"no turnaround found: {A.message}")
    y_turn = A.y[:, -1].copy()
    t_turn = float(A.t[-1])
    a_turn = float(y_turn[0])

    # ---- Stage 1b: turnaround -> handoff e-folds below it ----
    def ev_hand(t, y):
        return np.log(y[0] / a_turn) + handoff
    ev_hand.terminal = True
    ev_hand.direction = -1
    B = solve_ivp(f, (t_turn, t_turn + t_max), y_turn, events=ev_hand,
                  rtol=rtol1, atol=atol1, method="DOP853")
    if B.status != 1:
        raise RuntimeError(f"handoff not reached: {B.message}")
    if B.y[1, -1] >= 0.0:
        raise RuntimeError("H did not go negative after the turnaround")

    t1 = np.concatenate([A.t, B.t[1:]])
    y1 = np.hstack([A.y, B.y[:, 1:]])
    stage1 = np.full(t1.shape, 1)

    # ---- Stage 2: e-fold time, H evolved as ln|H| ----
    yb = B.y[:, -1]
    sign = np.sign(yb[1])
    N0 = float(np.log(yb[0] / a_turn))

    def g(N, z):
        phi, pi, ldm, lb, lr, ls, u, _t = z
        dm, b, r, s = np.exp([ldm, lb, lr, ls])
        H = sign * np.exp(u)
        rho, pr = _rho_p(H, phi, pi, dm, b, r, s, p)
        dVd = -p.V0 * p.lam * np.exp(-p.lam * phi)
        return np.array([
            pi / H,
            -3.0 * pi - (dVd - p.beta * dm) / H,
            -3.0 - p.beta * pi / H,
            -3.0,
            -4.0,
            -6.0,
            -0.5 * (rho + pr) / (H * H),
            1.0 / H,
        ])

    if depth is None:
        endpoint = "planck_density"
        target = log(density_ratio * rho_today)

        def ev_end(N, z):
            phi, pi, ldm, lb, lr, ls, u, _t = z
            dm, b, r, s = np.exp([ldm, lb, lr, ls])
            rho, _ = _rho_p(sign * np.exp(u), phi, pi, dm, b, r, s, p)
            return np.log(rho) - target
        ev_end.terminal = True
        ev_end.direction = 1
        N_end, events = -float(n_limit), [ev_end]
    else:
        endpoint = "depth"
        N_end, events = -float(depth), None

    z0 = np.array([
        yb[2], yb[3],
        np.log(yb[4]), np.log(yb[5]), np.log(yb[6]), np.log(yb[7]),
        np.log(abs(yb[1])), 0.0,
    ])
    atol = np.array([1e-14, 1e-30, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-12])
    C = solve_ivp(g, (N0, N_end), z0, rtol=rtol2, atol=atol, method="DOP853", events=events)
    if endpoint == "planck_density":
        status = "ok" if C.status == 1 else f"density endpoint not reached: {C.message}"
    else:
        status = "ok" if C.status == 0 else f"stage2 stopped: {C.message}"

    NC = C.t[1:]
    phi, pi = C.y[0, 1:], C.y[1, 1:]
    dm, b, r, s = np.exp(C.y[2:6, 1:])
    H2 = sign * np.exp(C.y[6, 1:])
    t2 = t1[-1] + C.y[7, 1:]
    a2 = a_turn * np.exp(NC)

    a = np.concatenate([y1[0], a2])
    return CycleTrace(
        t=np.concatenate([t1, t2]),
        ln_a_today=np.log(a / s0.a),
        N_turn=np.log(a / a_turn),
        H=np.concatenate([y1[1], H2]),
        phi=np.concatenate([y1[2], phi]),
        pi_phi=np.concatenate([y1[3], pi]),
        rho_dm=np.concatenate([y1[4], dm]),
        rho_b=np.concatenate([y1[5], b]),
        rho_r=np.concatenate([y1[6], r]),
        rho_shear=np.concatenate([y1[7], s]),
        stage=np.concatenate([stage1, np.full(NC.shape, 2)]),
        t_turn=t_turn,
        a_turn_over_a_today=a_turn / s0.a,
        rho_today=rho_today,
        endpoint=endpoint,
        status=status,
    )


def _crossing(trace: CycleTrace, rho: np.ndarray, ratio: float) -> dict:
    """Depth where rho_total/rho_today first reaches ``ratio`` after the turnaround."""
    i_turn = int(np.argmin(np.abs(trace.N_turn)))
    # rho -> 0 at the turnaround itself; a floor only there keeps the log defined.
    x = np.log(np.maximum(rho[i_turn:], 1e-300) / trace.rho_today) - np.log(ratio)
    j = np.where((x[:-1] < 0.0) & (x[1:] >= 0.0))[0]
    if j.size == 0:
        return {"reached": False}
    k = i_turn + int(j[0])
    w = -x[k - i_turn] / (x[k - i_turn + 1] - x[k - i_turn])
    interp = lambda arr: float(arr[k] + w * (arr[k + 1] - arr[k]))  # noqa: E731
    return {
        "reached": True,
        "efolds_from_turnaround": -interp(trace.N_turn),
        "efolds_below_today": -interp(trace.ln_a_today),
        "phi_span_from_turnaround": abs(interp(trace.phi) - float(trace.phi[i_turn])),
    }


def contraction_summary(trace: CycleTrace) -> dict:
    """Depth of the contraction, reported rather than assumed."""
    i_turn = int(np.argmin(np.abs(trace.N_turn)))
    rho, pr = total_density(trace)
    return {
        "planck_convention_crossings": {
            "full_planck_density_c5_over_hbar_G2": _crossing(
                trace, rho, PLANCK_TO_CRITICAL_DENSITY_RATIO),
            "reduced_planck_density_over_8pi_squared": _crossing(
                trace, rho, REDUCED_PLANCK_TO_CRITICAL_DENSITY_RATIO),
        },
        "endpoint": trace.endpoint,
        "a_turn_over_a_today": float(trace.a_turn_over_a_today),
        "efolds_from_turnaround": float(-trace.N_turn[-1]),
        "efolds_below_today": float(-trace.ln_a_today[-1]),
        "phi_span_from_turnaround": float(abs(trace.phi[-1] - trace.phi[i_turn])),
        "density_ratio_at_end": float(rho[-1] / trace.rho_today),
        "w_at_end": float(pr[-1] / rho[-1]),
        # w is singular at the turnaround itself (rho -> 0), so measure from 1 e-fold in.
        "w_max_from_one_efold_after_turnaround": float(np.max((pr / rho)[trace.N_turn <= -1.0])),
    }


def total_density(trace: CycleTrace, p: CosmosParams = CosmosParams()):
    """rho_total and p_total along the trace (vectorised)."""
    V = p.V0 * (np.exp(-p.lam * trace.phi) - p.D)
    kin = 0.5 * trace.pi_phi ** 2
    rho = kin + V + trace.rho_dm + trace.rho_b + trace.rho_r + trace.rho_shear
    pr = kin - V + trace.rho_r / 3.0 + trace.rho_shear
    return rho, pr


def friedmann_residual(trace: CycleTrace, p: CosmosParams = CosmosParams()):
    """(3 H^2 - rho) / S, with S = sum of |component densities|; a witness, never enforced."""
    rho, _ = total_density(trace, p)
    V = p.V0 * (np.exp(-p.lam * trace.phi) - p.D)
    scale = (0.5 * trace.pi_phi ** 2 + np.abs(V) + trace.rho_dm + trace.rho_b
             + trace.rho_r + trace.rho_shear)
    return (3.0 * trace.H ** 2 - rho) / np.maximum(scale, 1e-300)


def power_law_invariants(trace: CycleTrace):
    """max |d ln(rho_i a^n_i)| for baryons, radiation and shear."""
    out = {}
    for name, n in _POWER.items():
        inv = np.log(getattr(trace, name)) + n * trace.ln_a_today
        out[name] = float(np.max(np.abs(inv - inv[0])))
    return out


def radiation_gauge(trace: CycleTrace, p: CosmosParams = CosmosParams()):
    """Log-space radiation diagnostic, regular through the turnaround (see archive module)."""
    rho, pr = total_density(trace, p)
    kin = 0.5 * trace.pi_phi ** 2
    V = p.V0 * (np.exp(-p.lam * trace.phi) - p.D)
    scale = kin + np.abs(V) + trace.rho_dm + trace.rho_b + trace.rho_r + trace.rho_shear
    inv = np.log(trace.rho_r) + 4.0 * trace.ln_a_today
    comps = np.vstack([kin, np.maximum(V, 0.0), trace.rho_dm, trace.rho_b,
                       trace.rho_r, trace.rho_shear])
    return {
        "log10_r_over_scale": np.log10(trace.rho_r / np.maximum(scale, 1e-300)),
        "log10_omega_r": np.log10(trace.rho_r / np.maximum(np.abs(rho), 1e-300)),
        "log10_r_over_kin": np.log10(trace.rho_r / np.maximum(kin, 1e-300)),
        "invariant_drift": inv - inv[0],
        "Hdot": -0.5 * (rho + pr),
        "dominant": np.argmax(comps, axis=0),
    }


def crunch_time_estimate(trace: CycleTrace):
    """Remaining cosmic time to the crunch, t_c - t = 1/(3|H|), exact only where w = 1."""
    H = trace.H
    est = np.full(H.shape, np.nan)
    m = H < 0.0
    est[m] = 1.0 / (3.0 * np.abs(H[m]))
    return est


if __name__ == "__main__":
    import json

    print(json.dumps(contraction_summary(integrate_cycle()), indent=2))
