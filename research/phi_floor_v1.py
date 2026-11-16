"""RESEARCH ONLY. The floor run (approved Nov 15 2026): GEAR-91 on (beta_S = 0, beta_D = -0.04), hand floor off, D's
ejected vacuum the only floor, lambda = 5000 one-budget cell. Compares with beta_D = 0 (phi blind to the floor) and with
the derived background prediction. Note: results/phi-floor-2026-11-15/README.md.

Derived background (homogeneous, archive normalisation as the lane: H^2 = rho/3):
  V(phi) = V0_1B e^(-phi);  vacuum rho_v(phi) = 8 pi V_SD at the sea's rest = -8 pi e^(-4 bD phi)/(2 lam);  p_v = -rho_v
  Pi' = -3 H Pi - V' - rho_v'(phi) + beta_dm rho_dm;  rho_dm' = -3 H rho_dm - beta_dm rho_dm Pi;  baryons a^-3, radiation a^-4
  H' = -H^2 - (rho + 3p)/6.  With bD = 0, rho_v is the constant floor and phi is blind to it.
Checks: (1) background cell (no matter) vs this prediction, for bD = 0 and -0.04: phi, Pi, H; geometry counts the floor
once (H agrees with the prediction that has the floor once). (2) D's sea follows phi: D_s = D* e^(-bD phi).
(3) pocket cell: what the coupling changes in the reset (held step, release) and in phi's spread (phi feels no floor where
D is held near 0: dV_SD/dphi = bD e^(-2 bD phi) D^2).
Usage: python research/phi_floor_v1.py
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

from research.b6_frozen_test_v1 import LAM, V0_1B  # noqa: E402
from research.channel_companions_v1 import BETA  # noqa: E402

SCR = "/tmp/claude-0/-home-claude/a388b6ac-1f6f-5a55-8963-a6c8cbfba238/scratchpad/cell/"


def pv(d, i, q):
    return float(np.sum(d[q][i] * d["dV"][i]) / np.sum(d["dV"][i]))


def predict(d, bD):
    from research.cosmos_lane_follow_v1 import averages
    A, tau = averages(d)
    y0 = [A[0, 10], A[0, 2], A[0, 3], A[0, 4], A[0, 5], A[0, 6], A[0, 7]]

    def rv(phi):
        return -8 * math.pi * math.exp(-4 * bD * phi) / (2 * LAM)

    def rhs(_, y):
        lna, H, phi, pi, dm, b, r = y
        V = V0_1B * math.exp(-phi); Vp = -V
        rvp = (16 * math.pi * bD / LAM) * math.exp(-4 * bD * phi)                 # d rho_v / d phi
        r3p = 2 * pi * pi - 2 * V + dm + b + 2 * r - 2 * rv(phi)
        return [H, -H * H - r3p / 6, pi, -3 * H * pi - Vp - rvp + BETA * dm, -3 * H * dm - BETA * dm * pi, -3 * H * b, -4 * H * r]
    s = solve_ivp(rhs, (tau[0], tau[-1]), y0, t_eval=tau, method="DOP853", rtol=1e-11, atol=1e-18)
    return A, s.y.T


if __name__ == "__main__":
    print("Derived: extra force on phi from D's floor / phi's own slope at phi_today:",
          f"{(16*math.pi*0.04/LAM)*math.exp(-4*-0.04*0.179055) / (V0_1B*math.exp(-0.179055)):.4f}")
    runs = {"bD = 0": ("g_nm_b0", 0.0), "bD = -0.04": ("g_nm_bD", -0.04)}
    print("\n1. Background cells (no matter) against the derived background; differences cell - prediction")
    print("   run          t      phi - phi0 (cell)   (pred)      dphi      Pi cell     dPi/|Pi|     H cell      dH/H        D/D* cell   e^(-bD phi)")
    store = {}
    for lab, (f, bD) in runs.items():
        d = np.load(SCR + f + ".npz", allow_pickle=True); t = d["t"]
        A, Y = predict(d, bD); store[lab] = (d, A, Y)
        for tq in (10, 20, 30, 40, 50, 60, 70):
            i = int(np.argmin(np.abs(t - tq)))
            if i >= len(Y):
                continue
            Dm = pv(d, i, "D") / math.sqrt(2 / LAM)
            print(f"   {lab:11s} {t[i]:5.1f}   {A[i,3]-A[0,3]:+.6f}        {Y[i,2]-Y[0,2]:+.6f}   {A[i,3]-Y[i,2]:+.1e}   {A[i,4]:+.5f}   {(A[i,4]-Y[i,3])/abs(A[i,4]):+.1e}   {A[i,2]:.6f}   {(A[i,2]-Y[i,1])/A[i,2]:+.1e}   {Dm:.6f}    {math.exp(-bD*A[i,3]):.6f}")
    (d0, A0, _), (d1, A1, _) = store["bD = 0"], store["bD = -0.04"]
    t = d0["t"]
    print("\n   coupling's effect on the background (bD = -0.04 minus bD = 0):")
    for tq in (20, 30, 40, 50, 56):
        i = int(np.argmin(np.abs(t - tq)))
        if i >= min(len(A0), len(A1)):
            continue
        print(f"     t {t[i]:4.0f}: d(phi - phi0) = {(A1[i,3]-A1[0,3])-(A0[i,3]-A0[0,3]):+.3e} ({((A1[i,3]-A1[0,3])/(A0[i,3]-A0[0,3])-1)*100:+.2f}% of phi's roll),"
              f"  dH/H = {(A1[i,2]-A0[i,2])/A0[i,2]:+.2e}")
    print("\n3. Pocket cells, x_c 2.0 r_max 24: bD = 0 (u_x20_cold) vs bD = -0.04 (g_x20_bD)")
    from research.handoff_contract_v1 import RULES, handoff_time, pieces
    for lab, f in (("bD = 0", "u_x20_cold"), ("bD = -0.04", "g_x20_bD")):
        d = np.load(SCR + f + ".npz", allow_pickle=True); t = d["t"]; r = d["r"]
        Pr = {rr: pieces(d, rr) for rr in RULES}; th, _ = handoff_time(d, Pr)
        held = Pr["box2"][:, 2]; ipk = int(np.nanargmax(held))
        print(f"   {lab:11s} status {d['status']}; handoff t_h = {th.get('box2')}; raw held / S-D energy: peak {held[ipk]:.4f} at t = {t[ipk]:.1f};"
              f" at t = 8, 12, 24: {held[16]:.4f} {held[24]:.4f} {held[48]:.4f}")
        sp = []
        for tq in (10, 30, 50, 70):
            i = int(np.argmin(np.abs(t - tq)))
            if i >= len(t):
                continue
            ph = d["phi"][i]; pin = float(np.mean(ph[r < 6])); pout = float(np.mean(ph[r > 16]))
            sp.append(f"t={t[i]:.0f}: {pin - pout:+.2e}")
        print("      phi(pocket r<6) - phi(sea r>16): " + "  ".join(sp))
        print("      the collapse: t, H, x in the core (2<s^2>, r<4), D/D_sea in the core, held/S-D energy (raw)")
        for tq in (40, 44, 48, 52, 54, 56, 58):
            i = int(np.argmin(np.abs(t - tq)))
            if abs(t[i] - tq) > 0.6:
                continue
            m = r < 4; dV = d["dV"][i]
            x = 2 * LAM * float(np.sum(d["S"][i]**2 * dV * m) / np.sum(dV * m))
            Dr = float(np.sum(np.abs(d["D"][i]) * dV * m) / np.sum(dV * m)) / float(np.mean(np.abs(d["D"][i][r > 20])))
            print(f"        t {t[i]:4.1f}  H {-float(d['K'][i].mean())/3:+.4f}  x_core {x:6.3f}  D_core/D_sea {Dr:.3f}  held {held[i]:.4f}")


def held_coupled(path):
    """Raw held step with the sea's rest set by phi: h = (eD^2/2)(1 - Ds^2/(2 eD))^2, eD = e^(-2 bD phi) (bD = -0.04),
    over the S-D energy above the local vacuum -eD^2/2 (scaled). The fixed-vacuum measure of the handoff note reads
    the phi-driven shift of D_s as held step (0.24 at t = 12, 0.8 at t = 24 in the coupled cell): an artefact."""
    from research.reset_handoff_v1 import SQ, d_slow
    d = np.load(SCR + path + ".npz", allow_pickle=True); t = d["t"]; dV = d["dV"]; dr = float(d["r"][1] - d["r"][0])
    Ds = d_slow(d, "box2"); eD = np.exp(0.08 * d["phi"])
    h = 0.5 * eD**2 * (1 - Ds**2 / (2 * eD))**2
    out = []
    for i in range(len(t)):
        s, dd, ps, pd = (SQ * d[k][i] for k in ("S", "D", "PS", "PD"))
        g2 = (np.gradient(s, dr)**2 + np.gradient(dd, dr)**2) / d["grr"][i]
        V = 0.5 * s**2 - 0.5 * eD[i] * dd**2 + (s**4 + 6 * s**2 * dd**2 + dd**4) / 8 + 0.5 * eD[i]**2
        e = 0.5 * (ps**2 + pd**2) + 0.5 * g2 + V
        out.append(np.nansum(h[i] * dV[i]) / np.sum(e * dV[i]) if not np.isnan(Ds[i]).any() else np.nan)
    return t, np.array(out)
