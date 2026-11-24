"""Reanalysis of stored f_ring A6 data (#591): where |M| sits in r at each state save, and the shape of the stored
per-step M near the origin. No time step is taken. M is re-evaluated on the saved states with the same snapshot()
the runs used; the stored nine sample cells must match it exactly before anything is read.

Usage (repo root): python -m research.w1_fring_a6.m_profile OUT.json
"""
import json
import sys
from pathlib import Path

import numpy as np

from engine import scalar_system as ss, v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel
from research.w1_fring_smoke import origin_trace as ot

ARMS = Path("results/research/w1_fring_a6/arms")
STATIC = Path("results/research/w1_static_slice/native_static.npz")
KEYS = ("0.015", "0.8235", "1.6455", "3.291")
PERIOD = 3.29


def load_state(state, z, vac):
    for k in vac.ALL_FIELDS:
        np.copyto(getattr(state.geometry, k), z[f"geometry_{k}"])
    state.scalars = ss.ScalarFields(*(np.array(z[f"scalars_{k}"]) for k in ("S", "PS", "D", "PD", "phi", "Pi")))
    for species in ("dark_matter", "baryons", "radiation"):
        s = getattr(state.matter, species)
        for k in ("rest", "energy_t", "momentum_r"):
            np.copyto(getattr(s, k), z[f"matter_{species}_{k}"])


def extent(r, v):
    a = abs(v - v[-1])
    return {f"{th:g}_of_peak_out_to_r": float(r[np.where(a > th*a.max())[0][-1]]) for th in (.5, .1, .01, .001)}


def profile(r, M):
    a = abs(M)
    i = int(a.argmax())
    out = {"peak_cell": i, "peak_r": float(r[i]), "peak_M": float(M[i]),
           "half_peak_out_to_r": float(r[np.where(a > .5*a[i])[0][-1]]),
           "half_peak_from_r": float(r[np.where(a > .5*a[i])[0][0]])}
    for lo, hi in ((0, .5), (.5, 2), (2, 5), (5, 10), (10, 30), (30, 60)):
        w = (r >= lo) & (r < hi)
        out[f"max_abs_M_r{lo:g}-{hi:g}"] = float(a[w].max())
    return out


def series_shape(t, m):
    t, m = np.asarray(t), np.asarray(m)
    y = np.log(m)
    dy = np.diff(y)
    q = {f"t{x:g}": float(np.interp(x, t, m)) for x in (.015, PERIOD/4, PERIOD/2, 3*PERIOD/4, 3.291)}
    sign = np.sign(dy)
    turns = int(np.sum(sign[1:] != sign[:-1]))
    smooth = np.convolve(y, np.ones(21)/21, mode="valid")
    ts = t[10:-10]
    ds = np.diff(smooth)
    s = np.sign(ds)
    return {"values_at": q, "steps": int(len(t)), "rising_steps_fraction": float(np.mean(dy > 0)),
            "sign_changes_of_step_to_step_change": turns,
            "smoothed_21_step_local_maxima_t": [float(ts[i+1]) for i in range(len(s)-1) if s[i] > 0 > s[i+1]],
            "smoothed_21_step_local_minima_t": [float(ts[i+1]) for i in range(len(s)-1) if s[i] < 0 < s[i+1]]}


def main(out):
    go, vac, _ = ad.vendor_modules()
    grid = go.SphericalCellGrid(1200, 60.)
    r = np.asarray(grid.centers)
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=1200, r_max=60., outer_edge="open_edge", eight_node_areal=True,
                              initial_structure="none")
    with np.load(STATIC) as z:
        result = {"note": "reanalysis of stored A6 saves; no time step", "period": PERIOD,
                  "initial_field_extent": {k: extent(r, z["dr0.05_"+k]) for k in ("S", "D")}, "arms": {}}
    for name in ("A1", "A2", "B1", "B2"):
        rec = json.loads((ARMS/f"{name}.json").read_text())
        arm = {"saves": {}}
        for key in KEYS:
            with np.load(ARMS/f"{name}-t{key}.npz") as z:
                load_state(state, z, vac)
            snap = ot.snapshot(grid, state.geometry, state.scalars, state.matter, vac)
            stored = rec["saves"][key]["row"]["M"]
            if snap["cells"]["M"] != stored:
                raise SystemExit(f"{name} t{key}: recomputed M differs from the stored sample cells")
            gap = max(abs(a - b) for a, b in zip(snap["cells"]["M"], stored))
            raw = vac.constraints(grid, state.geometry)
            from engine.stress_energy import assemble_total_stress_energy
            M = raw["momentum"] - 8*np.pi*assemble_total_stress_energy(grid, state.geometry, state.scalars, state.matter).j
            arm["saves"][key] = {"stored_cells_match": gap == 0, **profile(r, M),
                                 "norm_l2_dr": float(np.sqrt(np.sum(M**2)*.05))}
        ps = [p for p in rec["per_step"] if p["t"] >= .015 - 1e-12]
        arm["per_step_cells_0_5"] = series_shape([p["t"] for p in ps], [max(abs(x) for x in p["M"]) for p in ps])
        result["arms"][name] = arm
    Path(out).write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
