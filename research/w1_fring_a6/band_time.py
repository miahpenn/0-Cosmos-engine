"""Band M (R 0.05-0.25) against time on stored data: A6 ring arms A1/A2 (dr 0.05, to t 3.291) and the dr 0.0125
replays (t 0.015). No time step.
- Cell centres: max |M| over cells with r in [0.05, 0.25] (dr 0.05: cells 1-4), every recorded step.
- Aster's definition: max |cubic spline of M through areal R < 0.4| at the 17 radii 0.05..0.25. Exact at the A6 saves
  (full state); every step only approximately, from cells 0-5 in r (labelled).
Usage (repo root): python -m research.w1_fring_a6.band_time ASTER_RUNS_DIR OUT.json [M|H]"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from engine import v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from research.w1_fring_a6.m_profile import load_state

ARMS = Path("results/research/w1_fring_a6/arms")
RUNS = Path("results/research/w1_fring_spatial_rerun/runs")
RADII = np.linspace(.05, .25, 17)
KEYS = ("0.015", "0.8235", "1.6455", "3.291")


def summary(t, v):
    t, v = np.asarray(t), np.asarray(v)
    i = int(np.argmax(v))
    edges = np.arange(0, t[-1] + 1e-12, .25)
    env = [{"t_from": float(a), "t_to": float(b), "max": float(v[(t >= a) & (t < b)].max())}
           for a, b in zip(edges[:-1], edges[1:]) if np.any((t >= a) & (t < b))]
    late = v[t >= edges[-1]] if np.any(t >= edges[-1]) else []
    if len(late):
        env.append({"t_from": float(edges[-1]), "t_to": float(t[-1]), "max": float(np.max(late))})
    return {"peak_t": float(t[i]), "peak": float(v[i]), "last_t": float(t[-1]), "last": float(v[-1]),
            "peak_over_last": float(v[i]/v[-1]), "envelope_0.25": env}


def main(aster_runs, out, Q="M"):
    res = {"quantity": Q, "a6": {}, "dr0.0125": {}, "signs_and_dt_ratio": {}}
    go, vac, _ = ad.vendor_modules()
    grid = go.SphericalCellGrid(1200, 60.)
    r = np.asarray(grid.centers)
    state = V55ProductionKernel().initialize(resolution=1200, r_max=60., outer_edge="open_edge",
                                             eight_node_areal=True, initial_structure="none")
    for arm in ("A1", "A2"):
        rec = json.loads((ARMS/f"{arm}.json").read_text())
        rows = rec["per_step"]
        t = [row["t"] for row in rows]
        cc = [float(np.max(np.abs(row[Q][1:5]))) for row in rows]
        ap = [float(np.max(np.abs(CubicSpline(r[:6], row[Q][:6])(RADII)))) for row in rows]
        saves = {}
        for key in KEYS:
            with np.load(ARMS/f"{arm}-t{key}.npz") as z:
                load_state(state, z, vac)
            g = state.geometry
            st = assemble_total_stress_energy(grid, g, state.scalars, state.matter)
            raw = vac.constraints(grid, g)
            M = raw["momentum"] - 8*np.pi*st.j if Q == "M" else raw["hamiltonian"] - 16*np.pi*st.rho
            stored = np.array(rec["saves"][key]["row"][Q][:6])
            if M[:6].tobytes() != stored.tobytes():
                raise SystemExit(f"STOP: {arm} {key} {Q} does not reproduce the stored save")
            R = r*np.sqrt(g.b)/g.X
            m = R < .4
            saves[key] = {"cells_1_4": M[1:5].tolist(),
                          "sign_flips_of_3": int(np.sum(np.sign(M[2:5]) != np.sign(M[1:4]))),
                          "cell_centres": float(np.max(np.abs(M[1:5]))),
                          "aster_17_radii": float(np.max(np.abs(CubicSpline(R[m], M[m])(RADII))))}
        res["a6"][arm] = {"dt": rec["dt"], "cell_centres_every_step": summary(t, cc),
                          "aster_approx_cells0_5_every_step": summary(t, ap), "saves_exact": saves}
    for key in KEYS:
        s1, s2 = res["a6"]["A1"]["saves_exact"][key], res["a6"]["A2"]["saves_exact"][key]
        res["signs_and_dt_ratio"][key] = {k: s2[k]/s1[k] for k in ("cell_centres", "aster_17_radii")}
    for tag in (("A-dr0.0125-dt0.000375", "A-dr0.0125-dt0.0001875") if Q == "M" else ()):
        rows = json.loads((RUNS/f"{tag}.json").read_text())["rows"]
        res["dr0.0125"][tag + " (our slice, every step, cell centres)"] = summary(
            [x["t"] for x in rows], [x["band_max_abs_M"] for x in rows])
    c = json.loads((Path(aster_runs)/"spatial_comparison.json").read_text())
    for dt in ("0.000375", "0.0001875"):
        rows = c["runs"]["0.0125"][dt]["rows"]
        with np.load(Path(aster_runs)/f"sp/0.0125-{dt}.npz") as z:
            Ms, rr = z[Q], z["r"]
        w = (rr >= .05) & (rr <= .25)
        sg = np.sign(Ms[-1][w])
        res["dr0.0125"][f"Aster input dt {dt} (5 samples)"] = {
            "t": [x["t"] for x in rows], "aster_17_radii": [x["band_max"][Q] for x in rows],
            "cell_centres": [float(np.max(np.abs(m[w]))) for m in Ms],
            "sign_flips_t0.015": int(np.sum(sg[1:] != sg[:-1])), "of": int(w.sum()) - 1}
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:4])
