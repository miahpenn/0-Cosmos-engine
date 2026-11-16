"""RESEARCH ONLY. Lattice-cell run with the GEAR-91 coupling on (beta_S = 0, beta_D given): the floor run approved Nov 15
2026. The one-budget cell already has the hand floor off and D's ejected vacuum as the floor in geometry; with
beta_D != 0, phi also feels that floor (8 pi dV_SD/dphi). Same driver and records as research/lattice_cell_run_v1.py.
Usage: python research/lattice_cell_run_coupled_v1.py CASE x_c N RMAX T out.npz beta_D [closure=ii] [every=0.5]
"""
import sys

sys.path.insert(0, ".")

import research.gear91_coupling_v1 as g91  # noqa: E402
import research.lattice_cell_run_v1 as drv  # noqa: E402
import research.lattice_cell_v1 as lc  # noqa: E402

if __name__ == "__main__":
    a = sys.argv
    beta_D = float(a[7]); closure = a[8] if len(a) > 8 else "ii"; every = float(a[9]) if len(a) > 9 else 0.5
    base = lc.cell_initial
    lc.cell_initial = lambda case, x_c, _b, N, RMAX, **kw: base(case, x_c, beta_D, N, RMAX, **kw)  # noqa: E731
    print(f"[coupled cell {a[1]} x_c={a[2]} N={a[3]} r_max={a[4]} T={a[5]} beta_D={beta_D} closure={closure} every={every}]", flush=True)
    with g91.coupled(beta_S=0.0, beta_D=beta_D):
        print(drv.run(a[1], float(a[2]), int(a[3]), float(a[4]), float(a[5]), a[6], closure=closure, every=every))
