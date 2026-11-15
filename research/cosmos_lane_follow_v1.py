"""RESEARCH ONLY. Does a cosmos lane that carries the S-D species follow the cell over time?
(Proposed handoff item 2, done on STORED trajectories: no new lattice run. Note: COMPANIONS_NOTE.md.)

The lane is integrated on its own, in the dust clock, from the cell's averages at t0, and compared with the cell's
stored averages afterwards. It is told nothing after t0.

    dH/dtau = -H^2 - (rho + 3p)/6                 Raychaudhuri form (B4)
    rho + 3p = 2 pi^2 - 2 V(phi) + rho_dm + rho_b + 2 rho_r - 2 rho_vac + (1 + 3w) rho_SD
    dphi = pi;  dpi = -3H pi - V'(phi) + beta rho_dm;  drho_dm = -3H rho_dm - beta rho_dm pi
    drho_b = -3H rho_b;  drho_r = -4H rho_r;  drho_SD = -3H (1 + w) rho_SD
Variants: "none" (the lane as it is, no S-D species), "w=0", "w=read" (w from the cell's own cycle mean after settling),
and the control "friedmann" (the lane's -(rho+p)/2 form with the S-D species at w = 0).
rho_b and rho_r at t0 are solved from the stored totals (velocities neglected); V0 and the vacuum constant are the
one-budget values. The companion terms (a, phi, pi, rho) are NOT fed in: their size is what the residual measures.
Usage: python research/cosmos_lane_follow_v1.py RUN.npz [RUN.npz ...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

from research.b6_frozen_test_v1 import KAP, LAM, P_FULL, P_NOD, stored, v_sd  # noqa: E402
from research.channel_companions_v1 import BETA, slice_terms  # noqa: E402


def averages(d):
    t = d["t"]; rows = []
    for i in range(len(t)):
        s = slice_terms(d, i); f = stored(d, i); dV = d["dV"][i]
        av = lambda q: float(np.sum(q * dV) / np.sum(dV))                                 # noqa: E731
        qfl, rfl = av(f["q_fl"]), av(f["rho_fl"])
        rho_b = 2 * ((rfl - s["rho"]) - (qfl - 0.5 * s["rho"])); rho_r = rfl - s["rho"] - rho_b
        rows.append((t[i], s["N"], -s["K"] / 3, s["phi"], s["Pi"], s["rho"], rho_b, rho_r, KAP * s["e"], s["w"] / max(s["e"], 1e-300),
                     math.log(np.sum(dV)) / 3))
    A = np.array(rows)
    tau = np.concatenate([[0], np.cumsum(0.5 * (A[1:, 1] + A[:-1, 1]) * np.diff(A[:, 0]))])
    return A, tau


def lane(y0, taus, P, rho_vac, w, form):
    def rhs(_, y):
        lna, H, phi, pi, dm, b, r, sd = y
        V = P.V0 * (math.exp(-P.lam * phi) - P.D); Vp = -P.V0 * P.lam * math.exp(-P.lam * phi)
        rho = 0.5 * pi * pi + V + dm + b + r + rho_vac + sd
        r3p = 2 * pi * pi - 2 * V + dm + b + 2 * r - 2 * rho_vac + (1 + 3 * w) * sd
        dH = -H * H - r3p / 6 if form == "raychaudhuri" else -0.5 * (rho + (r3p - rho) / 3)
        return [H, dH, pi, -3 * H * pi - Vp + BETA * dm, -3 * H * dm - BETA * dm * pi, -3 * H * b, -4 * H * r, -3 * H * (1 + w) * sd]
    s = solve_ivp(rhs, (taus[0], taus[-1]), y0, t_eval=taus, method="DOP853", rtol=1e-11, atol=1e-18)
    return s.y.T


if __name__ == "__main__":
    for path in sys.argv[1:]:
        d = np.load(path, allow_pickle=True)
        A, tau = averages(d); t = A[:, 0]
        has_D = bool(np.any(d["D"][0] != 0)); P = P_FULL if has_D else P_NOD
        rho_vac = KAP * v_sd(0.0, math.sqrt(2 / LAM)) if has_D else 0.0
        print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, {d['start']}, r_max {d['r'][-1]+0.05:g}, closure {d['closure']};"
              f" turnaround at t = {t[np.argmin(np.abs(A[:,2]))]:.1f}; rho_SD/rho_dm at t = 12: {A[24,8]/A[24,5]:.4f}")
        for t0 in (12.0, 0.0):
            i0 = int(np.argmin(np.abs(t - t0)))
            w_read = float(np.mean(A[(t >= 12) & (t < 54), 9]))
            y0 = [A[i0, 10], A[i0, 2], A[i0, 3], A[i0, 4], A[i0, 5], A[i0, 6], A[i0, 7], A[i0, 8]]
            print(f"  lane started from the cell's averages at t = {t0:g} (rho_b/rho_dm {A[i0,6]/A[i0,5]:.4f}, rho_r/rho_dm {A[i0,7]/A[i0,5]:.4f}); w read = {w_read:+.4f}")
            print("    variant        at t:   dH/H(t0)     d ln a      d phi/phi-change   d pi/|pi|   d rho_dm/rho_dm   d rho_SD/rho_SD")
            for name, w, form, sd in (("none", 0.0, "raychaudhuri", 0.0), ("w=0", 0.0, "raychaudhuri", 1.0), ("w=read", w_read, "raychaudhuri", 1.0),
                                      ("friedmann", 0.0, "friedmann", 1.0)):
                y = list(y0); y[7] *= sd
                Y = lane(y, tau[i0:], P, rho_vac, w, form)
                for tq in (24.0, 39.0, 54.0, 66.0):
                    if tq <= t0:
                        continue
                    i = int(np.argmin(np.abs(t - tq))); k = i - i0
                    sdres = (Y[k, 7] / A[i, 8] - 1) if sd else float("nan")
                    print(f"    {name:10s}    {t[i]:5.1f}:  {(Y[k,1]-A[i,2])/abs(A[i0,2]):+.2e}   {Y[k,0]-A[i,10]:+.2e}   {(Y[k,2]-A[i,3])/abs(A[i,3]-A[i0,3]):+.2e}"
                          f"          {(Y[k,3]-A[i,4])/abs(A[i,4]):+.2e}   {Y[k,4]/A[i,5]-1:+.2e}         {sdres:+.2e}")
