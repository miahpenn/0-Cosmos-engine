"""RESEARCH ONLY. Acceptance test T1 of the lattice-cell wall (HANDOFF Nov 14 2026 (4)): the cell with no pocket.
Reads lattice_cell_run_v1 output for CASE = nomatter and compares with research/cell_background_ode_v1.
Usage: python research/lattice_cell_t1_v1.py RUN.npz [RUN2.npz ...]
"""
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402

import research.dm_one_budget_v1 as dm  # noqa: E402
from engine import scalar_system as ss  # noqa: E402
from engine.production_kernel import BETA_DM  # noqa: E402
from research.cell_background_ode_v1 import background  # noqa: E402

if __name__ == "__main__":
    print("T1: no pocket. Spread = (max - min)/|mean| over the cell, largest over all snapshots.")
    for p in sys.argv[1:]:
        d = np.load(p, allow_pickle=True)
        t = d["t"]; r = d["r"]
        sp = lambda f: float(np.max(np.ptp(f, axis=1) / np.maximum(np.abs(f.mean(1)), 1e-300)))  # noqa: E731
        ss.V0 = dm.C2 * dm.V0_PROD; ss.D_POT = 0.0
        dm0, b0, rr0 = dm.C2 * dm.vm.RHO_DM_PRESENT, dm.C2 * dm.vm.RHO_B_PRESENT, dm.C2 * dm.vm.RHO_R_PRESENT
        VD = 8 * math.pi * (-1.0 / (2 * dm.LAM))
        bg = background(float(d["H0"]), float(d["phi"][0, 0]), float(d["Pi"][0, 0]), dm0, b0, rr0, ss.V0, ss.LAMBDA, 0.0, VD,
                        BETA_DM, float(t[-1]), n_out=len(t))
        H = -d["K"].mean(1) / 3
        i_turn = int(np.argmax(H < 0)) if np.any(H < 0) else -1
        relH = np.abs(H - bg["H"]) / float(d["H0"])
        print(f"\n{p.split('/')[-1]}: N = {len(r)}, r_max = {r[-1] + 0.5*(r[1]-r[0]):g}, closure {d['closure']}, cfl {float(d['cfl'])}, t_end = {t[-1]:.1f}, status {d['status']}")
        print(f"  spreads: K {sp(d['K']):.1e}, X {sp(d['X']):.1e}, alpha {sp(d['alpha']):.1e}, phi {sp(d['phi']):.1e}, Pi {sp(d['Pi']):.1e},"
              f" D {sp(d['D']):.1e}, dust {sp(d['rho_dm']):.1e}; max |a/b - 1| {np.abs(d['a']/d['b'] - 1).max():.1e};"
              f" max |Aa| {np.abs(d['Aa']).max():.1e}; max |alpha - 1| {np.abs(d['alpha'] - 1).max():.1e}")
        print(f"  turnaround (H = 0): cell t = {t[i_turn]:.2f} (first snapshot with H < 0), ODE t = "
              f"{bg['t'][int(np.argmax(bg['H'] < 0))]:.2f}; e-folds at turnaround (ODE) {bg['N'][int(np.argmax(bg['H'] < 0))]:.4f}")
        for tt in (10, 20, 40, t[-1]):
            i = int(np.argmin(np.abs(t - tt)))
            print(f"  t = {t[i]:6.2f}: H cell {H[i]:+.10e}, ODE {bg['H'][i]:+.10e}, |diff|/H0 {relH[i]:.2e};"
                  f" phi diff {abs(d['phi'][i].mean() - bg['phi'][i]):.1e}; dust ratio - 1 {d['rho_dm'][i].mean()/bg['dm'][i] - 1:+.1e}")
        print(f"  max |H_cell - H_ODE|/H0 over the run: {relH.max():.2e};  max constraint residual |H| {np.abs(d['Hc']).max():.1e}, |M| {np.abs(d['Mc']).max():.1e}")
