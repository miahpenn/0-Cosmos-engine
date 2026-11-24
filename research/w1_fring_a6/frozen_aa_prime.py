"""Frozen-slice D_r Aa comparison, regularized_aa_prime vs the current even stencil, on stored states
(method fixed in results/research/w1_fring_a6_frozen_aa_prime/PRERUN.md). No time step.
Usage (repo root): python -m research.w1_fring_a6.frozen_aa_prime OUT.json"""
import json
import sys
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from research.w1_fring_a6.m_profile import load_state

ARMS = Path("results/research/w1_fring_a6/arms")
RUNS = Path("results/research/w1_fring_spatial_rerun/runs")
KEYS = ("0.015", "0.8235", "1.6455", "3.291")
EPS = np.finfo(float).eps
BAND = (.05, .25)


def absd(grid, u):
    e = np.empty(grid.n + 4)
    e[2:-2] = np.abs(u)
    e[1], e[0] = e[2], e[3]
    e[-2:] = e[-3]
    return (e[:-4] + 8*e[1:-3] + 8*e[3:-1] + e[4:])/(12*grid.dr)


class Reader:
    def __init__(self, go, vac):
        self.go, self.vac, self.grids = go, vac, {}

    def grid(self, n):
        if n not in self.grids:
            self.grids[n] = self.go.SphericalCellGrid(n, 60.)
        return self.grids[n]

    def dM(self, grid, Aa):
        return self.vac.regularized_aa_prime(grid, Aa) - self.vac._d(grid, Aa, 1)

    def rnd(self, grid, Aa):
        r = np.asarray(grid.centers)
        Al = 1.5*Aa/r**2
        return EPS*(absd(grid, Aa) + (4/3)*r*np.abs(Al) + (2/3)*r**2*absd(grid, Al))

    def read(self, grid, Aa, M_current=None):
        vac, r, dr = self.vac, np.asarray(grid.centers), grid.dr
        Q = vac._quadratic_extrapolate_first(grid, Aa/r**2)
        delta = float(Aa[0] - r[0]**2*Q)
        Aa_reg = Aa.copy()
        Aa_reg[0] = r[0]**2*Q
        dM, dMr, rnd = self.dM(grid, Aa), self.dM(grid, Aa_reg), self.rnd(grid, Aa)
        k = np.zeros(grid.n)
        k[:3] = np.array([4, -14/3, 2])/dr
        band = np.where((r >= BAND[0]) & (r <= BAND[1]))[0]
        q1 = np.abs(dMr)/rnd
        i1 = int(np.argmax(q1))
        out = {"delta": delta, "dM_cells_0_5": dM[:6].tolist(), "predicted_cells_0_2": (delta*k[:3]).tolist(),
               "rnd_cells_0_5": rnd[:6].tolist(),
               "linear_split_gap_over_rnd": float(np.max(np.abs(dM - dMr - delta*k)/rnd)),
               "R1": {"max_ratio_to_rnd": float(q1[i1]), "cell": i1, "dM_reg_there": float(dMr[i1]),
                      "band_max_ratio_to_rnd": float(q1[band].max()), "pass": bool(q1.max() <= 10)}}
        resolved = abs(4*delta/dr) > 10*rnd[0]
        beyond = bool(np.all(np.abs(dM[3:]) <= 10*rnd[3:]))
        q = np.abs(dM[3:])/rnd[3:]
        r2 = {"resolved": bool(resolved), "beyond_cell2_max_ratio_to_rnd": float(q.max()),
              "beyond_cell2_max_cell": int(np.argmax(q)) + 3}
        if resolved:
            tol = np.maximum(.01*np.abs(delta*k[:3]), 10*rnd[:3])
            r2.update(ok_cells_0_2=(np.abs(dM[:3] - delta*k[:3]) <= tol).tolist(),
                      rel_gap_cells_0_2=(np.abs(dM[:3] - delta*k[:3])/np.abs(delta*k[:3])).tolist(),
                      ratio_1_0=float(dM[1]/dM[0]), ratio_2_0=float(dM[2]/dM[0]))
            r2["pass"] = bool(all(r2["ok_cells_0_2"]) and beyond)
        else:
            r2["pass"] = "unresolved"
        if M_current is not None:
            r2["current_M0"] = float(M_current[0])
            r2["spike_share_of_current_M0"] = float((16/3)*delta/dr/M_current[0])
        out["R2"] = r2
        qb = [(int(i), float(abs(dM[i])/rnd[i])) for i in band if i >= 3]
        x = (np.arange(grid.n) + .5)*dr
        j0 = max(0, int(np.searchsorted(x, BAND[0])) - 2)
        nodes = np.arange(j0, j0 + 4)
        w = [np.prod([(BAND[0] - x[m])/(x[l] - x[m]) for m in nodes if m != l]) for l in nodes]
        interp = float(sum(wl*dM[l] for wl, l in zip(w, nodes)))
        w2 = dict(zip(nodes.tolist(), w)).get(2, 0.)
        out["R3"] = {"band_cells": [int(band[0]), int(band[-1])],
                     "max_ratio_to_rnd_band_cells_ge3": max(v for _, v in qb),
                     "pass": bool(all(v <= 10 for _, v in qb)),
                     "interp_dM_at_R0.05": interp, "predicted_2delta_over_dr_times_w2": float(2*delta/dr*w2)}
        q4 = np.abs(dM)/rnd
        out["R4_all_cells_max_ratio_to_rnd"] = float(q4.max())
        out["R4_all_cells_at_roundoff"] = bool(q4.max() <= 10)
        cur = np.abs(self.vac._d(grid, Aa, 1))
        a = np.abs(dMr)
        out["regular_part_profile (added after reading, not gated)"] = {
            "max_abs_dM_reg_cells_0_5": float(a[:6].max()), "max_abs_dM_reg_band": float(a[band].max()),
            "max_abs_dM_reg": float(a.max()), "cell_of_max": int(a.argmax()),
            "dM_reg_cells_0_5": dMr[:6].tolist(), "delta_k_cells_0_2": (delta*k[:3]).tolist(),
            "band_max_dM_reg_over_current_Aa_prime": float(np.max(a[band]/cur[band])),
            "last_r_above_1e-3_of_max": float(r[np.where(a > 1e-3*a.max())[0][-1]]),
            "last_r_above_1e-1_of_max": float(r[np.where(a > 1e-1*a.max())[0][-1]])}
        return out


def synthetic(rd, grid):
    r, dr = np.asarray(grid.centers), grid.dr
    res = {}
    for name, Aa in (("r2", 1e-7*r**2), ("r4", 1e-7*r**4)):
        res[name] = rd.read(grid, Aa)["R1"]
    A6 = 1e-7
    Aa = A6*r**6
    pred = 24*A6*dr**4*r
    dM = rd.dM(grid, Aa)
    res["r6"] = {"max_rel_gap_to_24A6dr4r_cells_0_99": float(np.max(np.abs(dM[:100] - pred[:100])/pred[:100])),
                 "R1": rd.read(grid, Aa)["R1"]}
    base = 1e-7*r**2
    for name, cell in (("delta_cell0", 0), ("delta_cell1", 1)):
        Aa = base.copy()
        Aa[cell] += 4.5e-13
        rr = rd.read(grid, Aa)
        res[name] = {"R2": rr["R2"], "R1": rr["R1"]["pass"]}
    ok = (res["r2"]["pass"] and res["r4"]["pass"] and res["r6"]["max_rel_gap_to_24A6dr4r_cells_0_99"] < 1e-6
          and res["delta_cell0"]["R2"]["pass"] is True and res["delta_cell1"]["R2"]["pass"] is not True)
    res["discriminates"] = bool(ok)
    return res


def main(out):
    go, vac, _ = ad.vendor_modules()
    rd = Reader(go, vac)
    grid = rd.grid(1200)
    res = {"synthetic": synthetic(rd, grid)}
    if not res["synthetic"]["discriminates"]:
        raise SystemExit("STOP: synthetic cases do not discriminate")
    state = V55ProductionKernel().initialize(resolution=1200, r_max=60., outer_edge="open_edge",
                                             eight_node_areal=True, initial_structure="none")
    r = np.asarray(grid.centers)
    a6, checks = {}, {}
    for arm in ("A1", "A2", "B1", "B2"):
        rec = json.loads((ARMS/f"{arm}.json").read_text())
        for key in KEYS:
            with np.load(ARMS/f"{arm}-t{key}.npz") as z:
                load_state(state, z, vac)
            g = state.geometry
            j = assemble_total_stress_energy(grid, g, state.scalars, state.matter).j
            M = vac.constraints(grid, g)["momentum"] - 8*np.pi*j
            stored = np.array(rec["saves"][key]["row"]["M"][:6])
            checks[f"{arm}-{key}"] = bool(M[:6].tobytes() == stored.tobytes())
            if not checks[f"{arm}-{key}"]:
                raise SystemExit(f"STOP: {arm} {key} does not reproduce stored M")
            src = vac.momentum_constraint
            orig = vac._d
            vac.momentum_constraint = lambda gr, st: (src(gr, st) - orig(gr, st.Aa, 1)
                                                      + vac.regularized_aa_prime(gr, st.Aa))
            M_patched = vac.constraints(grid, g)["momentum"] - 8*np.pi*j
            vac.momentum_constraint = src
            rr = rd.read(grid, np.array(g.Aa), M)
            rr["patched_minus_current_vs_dM_max_abs"] = float(np.max(np.abs((M_patched - M) - rd.dM(grid, np.array(g.Aa)))))
            Aa_proj = np.array(g.Aa)
            Aa_proj[0] = r[0]**2*vac._quadratic_extrapolate_first(grid, Aa_proj/r**2)
            rr["R4_on_projected_A6"] = rd.read(grid, Aa_proj)["R4_all_cells_max_ratio_to_rnd"]
            a6[f"{arm}-t{key}"] = rr
    res["input_check_stored_M_bitwise"] = checks
    res["a6"] = a6
    fine = {}
    for f in sorted(RUNS.glob("*-split.npz")):
        with np.load(f) as z:
            rr_, three = np.array(z["r"]), np.array(z["three_Aa_over_r"])
        gr = rd.grid(len(rr_))
        fine[f.name.replace("-split.npz", "")] = rd.read(gr, rr_*three/3)
    res["spatial_rerun_t0.015"] = fine
    Path(out).write_text(json.dumps(res, indent=1))
    print("done", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
