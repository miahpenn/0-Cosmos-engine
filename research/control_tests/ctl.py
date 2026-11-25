"""Known-physics control tests (Lumen² design 1BK-Xj3r..., Amendments 1-3). Research only; no engine file changes.

Runtime overrides (Jeremiah, Nov 25 04:46 UTC, "Set V0 and lambda to 0"):
scalar_system.V0 = 0 and scalar_system.SD_QUARTIC_LAMBDA = 0, set after initialize() and before any state is built.
"""
import json
import math
import time
from pathlib import Path

import numpy as np

from engine import scalar_system as ss, v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel, ProductionState
from engine.stress_energy import assemble_total_stress_energy

R_MAX = 12.0
SAVE_EVERY_T = 0.3
DUST_DENSITY = 1e-150  # D*D must not underflow in recover_dust
T_END = 3.0
FIELDS_G = ("a", "b", "X", "Aa", "K", "Lambda")
FIELDS_S = ("S", "PS", "D", "PD", "phi", "Pi")
_, VAC, _ = ad.vendor_modules()
_ORIG = (ss.V0, ss.SD_QUARTIC_LAMBDA)


def make(dr, r_max=R_MAX, *, eight_node_areal=False):
    n = round(r_max / dr)
    kernel = V55ProductionKernel()
    ss.V0, ss.SD_QUARTIC_LAMBDA = _ORIG   # initialize() builds the production background, then the overrides apply
    state = kernel.initialize(resolution=n, r_max=r_max, outer_edge="open_edge", initial_structure="none",
                              eight_node_areal=eight_node_areal)
    ss.V0 = 0.0
    ss.SD_QUARTIC_LAMBDA = 0.0
    g = state.geometry
    for k in VAC.ALL_FIELDS:
        getattr(g, k)[:] = 0.0
    for k in ("a", "b", "X", "alpha"):
        getattr(g, k)[:] = 1.0
    state.scalars = ss.ScalarFields(*(np.zeros(n) for _ in FIELDS_S))
    # Dust cannot be zero (valencia.recover_dust needs rest > 0): keep the background dust at rest, scaled to
    # DUST_DENSITY (rest, energy_t, momentum_r are linear in density at rest). Radiation is zero.
    for sp in ("dark_matter", "baryons"):
        s = getattr(state.matter, sp)
        scale = DUST_DENSITY / float(np.max(np.abs(s.rest)))
        for k in ("rest", "energy_t", "momentum_r"):
            getattr(s, k)[:] *= scale
    for k in ("rest", "energy_t", "momentum_r"):
        getattr(state.matter.radiation, k)[:] = 0.0
    return kernel, state


def constraints(state):
    st = assemble_total_stress_energy(state.grid, state.geometry, state.scalars, state.matter)
    raw = VAC.constraints(state.grid, state.geometry)
    return raw["hamiltonian"] - 16 * np.pi * st.rho, raw["momentum"] - 8 * np.pi * st.j


def semi(kernel, state):
    """Engine's semi-discrete RHS at the stored state (stage 0 of step): solved alpha, kdot, RHS of every field."""
    grid = state.grid
    g0 = state.geometry.copy()
    alpha, kdot = kernel._solve_lapse_for_state(grid, g0, state.scalars, state.matter, state.electromagnetic)
    g0.alpha = alpha
    g0.beta.fill(0.0); g0.B.fill(0.0)
    st0 = ProductionState(grid=grid, geometry=g0, scalars=state.scalars, matter=state.matter,
                          electromagnetic=state.electromagnetic, t=state.t, tau=state.tau, e_folds=state.e_folds)
    srhs, mrhs, _ = kernel._rhs(st0)
    gt = kernel._geometry_stage_terms(grid, g0, state.scalars, state.matter, state.electromagnetic)
    R = {k: np.asarray(gt["explicit"][k]) for k in ("a", "b", "X")}
    for k in ("Aa", "K"):
        R[k] = np.asarray(gt["primary_l2"][k]) + np.asarray(gt["primary_l3"][k])
    R["Lambda"] = np.asarray(gt["lambda_l2"]) + np.asarray(gt["lambda_l3"])
    for k in FIELDS_S:
        R["s_" + k] = np.asarray(getattr(srhs, k))
    return g0, float(kdot), R, mrhs


def perturbed(state, R, h):
    g = state.geometry.copy()
    for k in FIELDS_G:
        setattr(g, k, getattr(g, k) + h * R[k])
    s = ss.ScalarFields(*(getattr(state.scalars, k) + h * R["s_" + k] for k in FIELDS_S))
    return ProductionState(grid=state.grid, geometry=g, scalars=s, matter=state.matter,
                           electromagnetic=state.electromagnetic, t=state.t, tau=state.tau, e_folds=state.e_folds)


def cdot(state, R, h):
    Hp, Mp = constraints(perturbed(state, R, h))
    Hm, Mm = constraints(perturbed(state, R, -h))
    return (Hp - Hm) / (2 * h), (Mp - Mm) / (2 * h)


def delta2_minus_delta4(grid, g):
    """(Delta_2 band operator - Delta_4) alpha, both on g (with g.alpha the solved alpha). Edge cell: nan."""
    a, b, X, al = (np.asarray(getattr(g, k)) for k in ("a", "b", "X", "alpha"))
    r = np.asarray(grid.centers); h = grid.dr
    d1 = lambda v: grid.cell_derivative_fourth(v, parity=1)
    c = -0.5 * d1(a) / a + d1(b) / b - d1(X) / X + 2.0 / r
    inv = X * X / a
    ext = np.concatenate(([al[0]], al))  # even ghost alpha[-1] = alpha[0]
    lap2 = np.full_like(al, np.nan)
    lo, mid, hi = ext[:-2], ext[1:-1], ext[2:]
    lap2[:-1] = inv[:-1] * ((lo - 2 * mid + hi) / h**2 + c[:-1] * (hi - lo) / (2 * h))
    lap4 = VAC.geometry_terms(grid, g)["lap_alpha"]
    return lap2 - lap4


def fixed_points(v, dr):
    """Cell 0, then the dr .05 band cells r .075,.125,.175,.225 (r^2-weighted averages of finer cells)."""
    f = round(0.05 / dr)
    r = (np.arange(len(v)) + 0.5) * dr
    out = [float(v[0])]
    for j in (1, 2, 3, 4):
        sl = slice(j * f, j * f + f)
        out.append(float(np.sum(r[sl]**2 * v[sl]) / np.sum(r[sl]**2)))
    return out


def save_arrays(state, g0, kdot, R, extra):
    d = {f"g_{k}": np.asarray(getattr(state.geometry, k)) for k in VAC.ALL_FIELDS}
    d.update({f"s_{k}": np.asarray(getattr(state.scalars, k)) for k in FIELDS_S})
    d["alpha_rhs"] = np.asarray(g0.alpha); d["kdot"] = np.array(kdot)
    d.update({f"R_{k}": v for k, v in R.items()})
    d.update(extra)
    return d


def run(name, kernel, state, dt, out_dir, *, with_cdot, t_end=T_END, on_save=None, max_steps=None, on_step=None):
    """on_step(n, state, per_save) -> dict or None, called on every step's state before the step (Test 2 near-save
    readings); its non-None rows go to rec['near_save']. Default None: unchanged from stage 1."""
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    dr = state.grid.dr
    steps = round(t_end / dt); per_save = round(SAVE_EVERY_T / dt)
    assert abs(steps * dt - t_end) < 1e-12 and abs(per_save * dt - SAVE_EVERY_T) < 1e-12 and per_save % 2 == 0
    rec = {"name": name, "dr": dr, "dt": dt, "n": state.grid.n, "r_max": state.grid.r_max, "saves": [], "per_step": []}
    t0 = time.time(); arrays = {}
    for n in range(steps + 1):
        do_save = n % per_save == 0
        if with_cdot or do_save:
            g0, kdot, R, _ = semi(kernel, state)
            H, M = constraints(state)
        if with_cdot:
            cH, cM = cdot(state, R, dt)
            rec["per_step"].append([n * dt] + fixed_points(H, dr) + fixed_points(M, dr)
                                   + fixed_points(cH, dr) + fixed_points(cM, dr))
        if on_step is not None:
            x = on_step(n, state, per_save)
            if x is not None:
                rec.setdefault("near_save", []).append(x)
        if do_save:
            row = {"step": n, "t": n * dt, "H": fixed_points(H, dr), "M": fixed_points(M, dr),
                   "K_rhs_nonuniform": fixed_points(R["K"] - kdot, dr),
                   "identity_d2_minus_d4": fixed_points(delta2_minus_delta4(state.grid, g0), dr),
                   "kdot": kdot, "cell0": {k: float(getattr(state.geometry, k)[0]) for k in ("a", "b", "X", "Aa", "K", "alpha")},
                   "alpha_rhs_cell0": float(g0.alpha[0]), "wall": time.time() - t0}
            if on_save is not None:
                row.update(on_save(state, g0, kdot, R))
            rec["saves"].append(row)
            k = n // per_save
            for key, v in save_arrays(state, g0, kdot, R, {"H": H, "M": M}).items():
                arrays[f"t{k}_{key}"] = v
            np.savez_compressed(out_dir / f"{name}.npz", **arrays)
            (out_dir / f"{name}.json").write_text(json.dumps(rec))
        if n == steps or (max_steps is not None and n >= max_steps):
            break
        state = kernel.step(state, dt)
    rec["wall_total"] = time.time() - t0; rec["steps_done"] = n
    (out_dir / f"{name}.json").write_text(json.dumps(rec))
    return state, rec
