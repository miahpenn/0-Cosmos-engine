"""RESEARCH ONLY. Step 1 of the hand-floor removal: the production kernel's runs (research/floor_swap_run_v1.py) against
the lanes (research/floor_swap_lane_v1.py; engine.cosmos for the old path). Proper-volume means over r < 20 (away from
the production outer edge). Note: results/floor-swap-2026-11-15/README.md.
Usage: python research/floor_swap_compare_v1.py
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

from engine.cosmos import HAND_FLOOR_REFERENCE_PARAMS, CosmosParams, construct_present_day, cosmos_rhs  # noqa: E402
from research.floor_swap_lane_v1 import BD, LAM_PHI0, LAM_TODAY, make  # noqa: E402

SCR = "/tmp/claude-0/-home-claude/a388b6ac-1f6f-5a55-8963-a6c8cbfba238/scratchpad/cell/"
P = HAND_FLOOR_REFERENCE_PARAMS  # hand-floor lane (labelled reference since step 2, Nov 16 2026)


def means(d, rcut=20.0):
    r = d["r"]; m = r < rcut; out = {}
    for k in ("phi", "Pi", "D"):
        out[k] = np.array([float(np.sum(d[k][i][m] * d["dV"][i][m]) / np.sum(d["dV"][i][m])) for i in range(len(d["t"]))])
    out["H"] = np.array([-float(np.mean(d["K"][i][m])) / 3 for i in range(len(d["t"]))])
    out["spread_phi"] = np.array([float(np.ptp(d["phi"][i][m])) for i in range(len(d["t"]))])
    return out


if __name__ == "__main__":
    s0 = construct_present_day(P); H0 = s0.H
    ya = [s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear]
    A = solve_ivp(lambda t, y: list(cosmos_rhs(tuple(y), P)), (0, 1200), ya, method="DOP853", rtol=1e-11, atol=1e-18, dense_output=True)
    lanes = {"fs_old": (lambda t: A.sol(t), None)}
    for f, lam in (("fs_l152", LAM_PHI0), ("fs_l157", LAM_TODAY)):
        rhs, y0f = make(lam, True)
        s = solve_ivp(rhs, (0, 1200), y0f(), method="DOP853", rtol=1e-11, atol=1e-20, dense_output=True, max_step=0.5)
        lanes[f] = (s.sol, lam)
    print("production kernel (N 20, r_max 40, r < 20 means) minus lane; t in lane units (Hubble time 176.4)")
    print("  run      t      phi kernel   d phi       Pi kernel    dPi/Pi      H kernel     dH/H        D kernel / D_rest(phi)   phi spread")
    for f, lab in (("fs_old", "old path: hand floor, lam 0 (engine lane)"), ("fs_l152", "swap, lam 1.52e6"), ("fs_l157", "swap, lam 1.57e6")):
        d = np.load(SCR + f + ".npz", allow_pickle=True); t = d["t"]; M = means(d); sol, lam = lanes[f]
        print(f"  {lab}  (status: run to t = {t[-1]:.0f})")
        for tq in (0, 25, 50, 100, 150, 200, 300, 400, 500, 700, 900, 1050):
            if tq > t[-1]:
                continue
            i = int(np.argmin(np.abs(t - tq))); y = sol(t[i])
            Drest = math.sqrt(2 * math.exp(-2 * BD * M["phi"][i]) / lam) if lam else float("nan")
            print(f"    {t[i]:7.1f}  {M['phi'][i]:.6f}  {M['phi'][i]-y[2]:+.1e}   {M['Pi'][i]:+.6f}  {(M['Pi'][i]-y[3])/abs(y[3]):+.1e}   {M['H'][i]:+.6f}  {(M['H'][i]-y[1])/abs(y[1]):+.1e}"
                  f"   {M['D'][i]/Drest if lam else 0.0:.6f}                {M['spread_phi'][i]:.1e}")
    # The production kernel drifts from its own lane on the old path too (open outer edge over Hubble times: gap #3).
    # So the swap is judged differentially: (swap - old) in the kernel against (swap - old) in the lanes.
    print("\ndifferential: (swap - old path) in the kernel vs in the lanes")
    print("  run              t      d phi kernel   d phi lane    ratio     d H kernel    d H lane     ratio")
    d0 = np.load(SCR + "fs_old.npz", allow_pickle=True); M0 = means(d0); t0 = d0["t"]
    for f in ("fs_l152", "fs_l157"):
        d = np.load(SCR + f + ".npz", allow_pickle=True); M = means(d); t = d["t"]; sol, _ = lanes[f]; sol0 = lanes["fs_old"][0]
        for tq in (25, 50, 100, 150, 200, 250, 300):
            i = int(np.argmin(np.abs(t - tq))); j = int(np.argmin(np.abs(t0 - t[i])))
            if abs(t0[j] - t[i]) > 1e-6 or i >= len(t):
                continue
            dk_p, dk_H = M["phi"][i] - M0["phi"][j], M["H"][i] - M0["H"][j]
            y, y0 = sol(t[i]), sol0(t[i]); dl_p, dl_H = y[2] - y0[2], y[1] - y0[1]
            print(f"  {f}  {t[i]:6.1f}   {dk_p:+.4e}    {dl_p:+.4e}   {dk_p/dl_p:.3f}    {dk_H:+.3e}    {dl_H:+.3e}   {dk_H/dl_H:.3f}")
