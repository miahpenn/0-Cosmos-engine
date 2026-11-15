"""RESEARCH ONLY. A short lattice-cell run that also stores the fluid primitives at each snapshot, so the lane-follow
test (results/channel-companions-2026-11-14) can be made independent of the assembled stress totals.
Approved as a short run, not a campaign (Nov 14 2026). Same cell, rule and driver as research/lattice_cell_run_v1.py.
Usage: python research/lattice_cell_run_primitives_v1.py CASE x_c N RMAX T out.npz [closure=ii] [start=cold]
"""
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

import research.lattice_cell_v1 as lc  # noqa: E402
from engine.matter_system import Species, primitives  # noqa: E402
from engine.v55_matter import _metrics_from_slice, metric_slice_from_q  # noqa: E402
from research.lattice_cell_run_v1 import KEYS, snapshot  # noqa: E402

PRIM = ("p_rho_dm", "p_rho_b", "p_rho_r", "p_W_dm", "p_W_b", "p_W_r")


def prims(st):
    mets = _metrics_from_slice(metric_slice_from_q(st.grid, st.geometry))
    out = {}
    for tag, name, sp in (("dm", "dark_matter", Species.DARK_MATTER), ("b", "baryons", Species.BARYON), ("r", "radiation", Species.RADIATION)):
        q = primitives(mets, getattr(st.matter, name), sp)
        out[f"p_rho_{tag}"] = np.array([x.rho for x in q]); out[f"p_W_{tag}"] = np.array([x.lorentz() for x in q])
    return out


def run(case, x_c, N, RMAX, T, out, closure="ii", start="cold", cfl=0.2, every=0.5):
    k, st, r, H0 = lc.cell_initial(case, x_c, 0.0, N, RMAX, closure=closure, start=start)
    dt = cfl * float(st.grid.dr); n_every = max(1, int(round(every / dt))); nsteps = int(round(T / dt))
    rec = {kk: [] for kk in KEYS + PRIM}
    for n in range(nsteps + 1):
        if n % n_every == 0:
            s = snapshot(k, st); s.update(prims(st))
            for kk in KEYS + PRIM:
                rec[kk].append(np.array(s[kk], dtype=float).copy())
            if n % (n_every * 8) == 0:
                print(f"t={st.t:6.2f} alpha0={st.geometry.alpha[0]:.4f} K0={st.geometry.K[0]:+.5f}", flush=True)
        if n < nsteps:
            st = k.step(st, dt)
    np.savez_compressed(out, r=r, H0=H0, status="ok", case=case, x_c=x_c, closure=closure, start=start, cfl=cfl,
                        **{kk: np.array(v) for kk, v in rec.items()})


if __name__ == "__main__":
    a = sys.argv
    run(a[1], float(a[2]), int(a[3]), float(a[4]), float(a[5]), a[6], closure=a[7] if len(a) > 7 else "ii",
        start=a[8] if len(a) > 8 else "cold")
    print("done")
