"""RESEARCH ONLY. Test of the bound-pocket virial in an expanding background (chat derivation, Drive:
"DERIVATION-bound-pocket-virial-in-an-expanding-background ... Nov15-2026"). Stored states only.
Note: results/handoff-contract-2026-11-15/CORE_VIRIAL.md.

Identity under test (S-D species only, core r < r_b in the cell's coordinate, proper distance l from the centre):
    Sum T^ii dV  =  (1/2) I''  -  (a''/a) I  -  W_x  +  S_b                                   (1)
  T^ii = 3 p_field (the field's full stress, flow included); I = Sum e l^2 dV; E = Sum e dV (e: S-D energy above vacuum,
  scaled); W_x = -Sum e l g_d dV, g_d = G dM_d / R_areal^2 (G = 1 units, dM_d from the total overdensity against the
  cell mean, per steradian x 4 pi), times KAP/LAM to bring e to the cell's units (cancels in ratios except for W_x: see
  W_UNITS); S_b = l_b T^rr(l_b) R_areal(l_b)^2 (the surface term; T^rr = radial stress, which for a field is not the
  isotropic p). Derivatives in the core's proper time (label time / mean core lapse); a''/a from the cell's proper
  volume in the fluids' time (wall lapse ~ 1).
Rest-frame pressure (2nd form of the derivation): 3 Sum p_rest = Sum T^ii - Sum (e + p) v^2, v = (R'/R) l, R^2 = I/E.
Rulers: core radius r_b = 4, 6, 8; smoothing (boxcar pi, applied twice / three times); time windows 13-30, 30-54.
Usage: python research/core_virial_v1.py RUN.npz [...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from research.b6_frozen_test_v1 import KAP, LAM  # noqa: E402
from research.reset_handoff_v1 import SQ, boxcar  # noqa: E402

TS = math.pi


def per_snapshot(d, i, r_b):
    r = d["r"]; dr = float(r[1] - r[0]); dV = d["dV"][i]; grr = d["grr"][i]
    s, dd, ps, pd = (SQ * d[k][i] for k in ("S", "D", "PS", "PD"))
    sr, dr_ = np.gradient(s, dr), np.gradient(dd, dr)
    Kk = 0.5 * (ps**2 + pd**2); G = 0.5 * (sr**2 + dr_**2) / grr
    V = 0.5 * (s**2 - dd**2) + (s**4 + 6 * s**2 * dd**2 + dd**4) / 8 + 0.5
    e = Kk + G + V; p = Kk - G / 3 - V; Trr = Kk + G - V
    ell = np.concatenate([[0.5 * dr * math.sqrt(grr[0])], 0.5 * dr * math.sqrt(grr[0]) + np.cumsum(0.5 * (np.sqrt(grr[1:]) + np.sqrt(grr[:-1])) * dr)])
    Ra = d["R"][i]
    rho = d["rho"][i]; rbar = float(np.sum(rho * dV) / np.sum(dV))
    dM = 4 * math.pi * np.cumsum((rho - rbar) * dV)                    # overdensity's mass inside each shell, G = 1
    g_d = dM / np.maximum(Ra, 1e-12)**2
    m = r < r_b; jb = int(np.searchsorted(r, r_b)) - 1
    conv = KAP / LAM                                                    # scaled S-D energy -> the cell's rho units
    return dict(E=float(np.sum(e * dV * m)), I=float(np.sum(e * ell**2 * dV * m)), TII=float(np.sum(3 * p * dV * m)),
                W=float(-np.sum(e * ell * g_d * dV * m)) / 1.0, Sb=float(ell[jb] * Trr[jb] * Ra[jb]**2),
                ep=float(np.sum((e + p) * ell**2 * dV * m)), al=float(np.sum(d["alpha"][i] * e * dV * m) / np.sum(e * dV * m)),
                lnV=math.log(float(np.sum(dV))) / 3, rho_sd_conv=conv, rho_c=float(rho[0]), e_c=float(e[0]))


def smooth(t, y, n):
    for _ in range(n):
        y = boxcar(t, y, TS)
    return y


def analyse(path, r_bs=(4.0, 6.0, 8.0), n_s=2):
    d = np.load(path, allow_pickle=True); t = d["t"]
    out = {}
    for r_b in r_bs:
        rows = [per_snapshot(d, i, r_b) for i in range(len(t))]
        A = {k: np.array([o[k] for o in rows]) for k in rows[0]}
        al = float(np.mean(A["al"][(t > 12) & (t < 54)]))
        S = {k: smooth(t, A[k], n_s) for k in ("E", "I", "TII", "W", "Sb", "ep", "lnV")}
        tau = t * al                                                    # core proper time (mean lapse)
        Idd = np.gradient(np.gradient(S["I"], tau), tau)
        lna = S["lnV"]; add = np.gradient(np.gradient(np.exp(lna), t), t) / np.exp(lna)   # a''/a, fluids' time
        Rdot_R = 0.5 * np.gradient(np.log(S["I"] / S["E"]), tau)
        flow = S["ep"] * Rdot_R**2
        out[r_b] = dict(t=t, E=S["E"], TII=S["TII"], half_Idd=0.5 * Idd, tidal=-add * S["I"], W=S["W"], Sb=S["Sb"], flow=flow, al=al)
    return d, out


if __name__ == "__main__":
    # W_x's units: e is scaled S-D energy; g_d is in the cell's G = 1 units, so W_x as written is already per unit of e
    # (the acceleration does not care what is accelerated). Same for the tidal term. No conversion needed.
    for path in sys.argv[1:]:
        d, out = analyse(path)
        print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, r_max {d['r'][-1]+0.05:g}, closure {d['closure']}")
        print("  terms of (1), each / (3E), time-averaged over the window; LHS = measured sum T^ii/(3E) = w with the flow included")
        print("  r_b  window   LHS (w)    (1/2)I''   -(a''/a)I   -W_x       S_b        RHS sum    LHS-RHS  |  flow/(3E)  w_rest measured  w_rest predicted  | l_b H")
        for r_b, o in out.items():
            t = o["t"]
            for a, b in ((13, 30), (30, 54), (13, 54)):
                m = (t >= a) & (t < b) & ~np.isnan(o["half_Idd"])
                f = lambda q: float(np.mean(q[m] / (3 * o["E"][m])))   # noqa: E731
                lhs, hI, td, W, Sb, fl = f(o["TII"]), f(o["half_Idd"]), f(o["tidal"]), f(-o["W"]), f(o["Sb"]), f(o["flow"])
                rhs = hI + td + W + Sb
                i = int(np.argmin(np.abs(t - 0.5 * (a + b)))); grr = d["grr"][i]; dr = float(d["r"][1] - d["r"][0])
                lb = float(np.sum(np.sqrt(grr[d["r"] < r_b])) * dr); H = -float(d["K"][i].mean()) / 3
                print(f"  {r_b:3.0f}  {a}-{b}   {lhs:+.4f}   {hI:+.4f}    {td:+.4f}     {W:+.4f}    {Sb:+.4f}    {rhs:+.4f}    {lhs-rhs:+.4f}  |  {fl:+.4f}     {lhs-fl:+.4f}          {rhs-fl:+.4f}          | {lb*H:.2f}")


# ---- second form: the identity in comoving terms (no large cancellation) ----
# The cell's slice observers move with the Hubble flow (CMC, zero shift), so the measured field stress is already the
# pressure relative to the flow. Writing l = a chi and I = a^2 J (J = Sum e chi^2 dV), the identity's Hubble pieces
# cancel analytically:
#   (1/2) I'' - (a''/a) I - Sum (e + p) H^2 l^2  =  (1/2) a^2 J'' + 2 a^2 H J' - w_c H^2 I
# so the prediction for the core's pressure relative to the flow is
#   3 Sum p = (1/2) a^2 J'' + 2 a^2 H J' - W_x + S_b - w_c H^2 I                                   (2)
# (w_c H^2 I is second order and is evaluated with the measured w).

def analyse2(path, r_bs=(4.0, 6.0, 8.0), n_s=2, clock="core"):
    d = np.load(path, allow_pickle=True); t = d["t"]
    out = {}
    for r_b in r_bs:
        rows = [per_snapshot(d, i, r_b) for i in range(len(t))]
        A = {k: np.array([o[k] for o in rows]) for k in rows[0]}
        al = float(np.mean(A["al"][(t > 12) & (t < 54)])) if clock == "core" else 1.0
        a = np.exp(A["lnV"] - A["lnV"][0]); J = A["I"] / a**2
        S = {k: smooth(t, q, n_s) for k, q in (("E", A["E"]), ("J", J), ("TII", A["TII"]), ("W", A["W"]), ("Sb", A["Sb"]),
                                                 ("I", A["I"]), ("lna", np.log(a)))}
        tau = t * al; a_s = np.exp(S["lna"]); H = np.gradient(S["lna"], tau)
        Jd = np.gradient(S["J"], tau); Jdd = np.gradient(Jd, tau)
        w_c = S["TII"] / (3 * S["E"])
        out[r_b] = dict(t=t, E=S["E"], TII=S["TII"], j2=0.5 * a_s**2 * Jdd, j1=2 * a_s**2 * H * Jd, W=-S["W"], Sb=S["Sb"],
                        hw=-w_c * H**2 * S["I"], al=al)
    return d, out


def report2(path, clock="core"):
    d, out = analyse2(path, clock=clock)
    print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, r_max {d['r'][-1]+0.05:g}; form (2), clock: {clock}; each term / (3E), window means")
    print("  r_b  window   measured w   (1/2)a^2 J''  2a^2 H J'   -W_x      S_b       -wH^2I     predicted w   meas - pred")
    for r_b, o in out.items():
        t = o["t"]
        for a, b in ((13, 30), (30, 54), (13, 54)):
            m = (t >= a) & (t < b) & ~np.isnan(o["j2"])
            f = lambda q: float(np.mean(q[m] / (3 * o["E"][m])))       # noqa: E731
            lhs, j2, j1, W, Sb, hw = (f(o[k]) for k in ("TII", "j2", "j1", "W", "Sb", "hw"))
            pr = j2 + j1 + W + Sb + hw
            print(f"  {r_b:3.0f}  {a}-{b}   {lhs:+.4f}      {j2:+.4f}       {j1:+.4f}    {W:+.4f}   {Sb:+.4f}   {hw:+.4f}    {pr:+.4f}       {lhs-pr:+.4f}")


def flux_check(path, r_bs=(4.0, 6.0, 8.0), n_s=2):
    """Size of the neglected boundary flux. Energy through the comoving boundary, from the core's energy balance:
    F = dE/dtau + 3 H (1 + w) E (> 0: into the core). Its contribution to J' is F chi_b^2; compared with the measured J'."""
    d = np.load(path, allow_pickle=True); t = d["t"]; r = d["r"]
    for r_b in r_bs:
        rows = [per_snapshot(d, i, r_b) for i in range(len(t))]
        A = {k: np.array([o[k] for o in rows]) for k in rows[0]}
        al = float(np.mean(A["al"][(t > 12) & (t < 54)]))
        a = np.exp(A["lnV"] - A["lnV"][0]); J = A["I"] / a**2
        jb = int(np.searchsorted(r, r_b)) - 1
        chib2 = np.array([(np.sum(np.sqrt(d["grr"][i][:jb + 1])) * float(r[1] - r[0]) / a[i])**2 for i in range(len(t))])
        E, Js, lna, TII = (smooth(t, q, n_s) for q in (A["E"], J, np.log(a), A["TII"]))
        tau = t * al; H = np.gradient(lna, tau); w = TII / (3 * E)
        F = np.gradient(E, tau) + 3 * H * (1 + w) * E
        Jd = np.gradient(Js, tau)
        for a_, b_ in ((13, 30), (30, 54)):
            m = (t >= a_) & (t < b_) & ~np.isnan(F)
            print(f"  r_b {r_b:3.0f}  {a_}-{b_}: F/(H E) = {np.mean(F[m]/(H[m]*E[m])):+.3f}   J'/(H J) = {np.mean(Jd[m]/(H[m]*Js[m])):+.3f}"
                  f"   flux share of J' (F chi_b^2 / J') = {np.mean(F[m]*smooth(t, chib2, n_s)[m])/np.mean(Jd[m]):+.2f}")


# ---- third form (corrects the second): the identity in the central observer's frame, built from the slice ----
# Slice observers move outward at H l relative to the central (Fermi) frame. To first order in that boost and with the
# field's momentum density j (radial, slice frame):
#   T^ii_F = T^ii_slice + 2 H Sum j l dV + H^2 Sum (e + p) l^2 dV
# Identity (1) in the central observer's proper time tau, with H = -K/3 (CMC: uniform per proper time):
#   T^ii_F = (1/2) d^2I/dtau^2 - (dH/dtau + H^2) I - W_x + S_b                                  (3)
# The sign of j is fixed by the core's energy balance: dE/dtau = -(flux out) + work, checked below.

def analyse3(path, r_b=6.0, n_s=2, jsign=1.0):
    d = np.load(path, allow_pickle=True); t = d["t"]; r = d["r"]; dr = float(r[1] - r[0])
    rows = []
    for i in range(len(t)):
        o = per_snapshot(d, i, r_b)
        grr = d["grr"][i]; dV = d["dV"][i]; m = r < r_b
        s, dd, ps, pd = (SQ * d[k][i] for k in ("S", "D", "PS", "PD"))
        j = jsign * -(ps * np.gradient(s, dr) + pd * np.gradient(dd, dr)) / np.sqrt(grr)
        ell = np.concatenate([[0.5 * dr * math.sqrt(grr[0])], 0.5 * dr * math.sqrt(grr[0]) + np.cumsum(0.5 * (np.sqrt(grr[1:]) + np.sqrt(grr[:-1])) * dr)])
        jb = int(np.searchsorted(r, r_b)) - 1
        o.update(Jl=float(np.sum(j * ell * dV * m)), jb=float(j[jb] * d["R"][i][jb]**2 * d["alpha"][i][jb]), H=-float(d["K"][i].mean()) / 3)
        rows.append(o)
    A = {k: np.array([o[k] for o in rows]) for k in rows[0]}
    al = float(np.mean(A["al"][(t > 12) & (t < 54)])); tau = t * al
    S = {k: smooth(t, A[k], n_s) for k in ("E", "I", "TII", "W", "Sb", "ep", "Jl", "jb", "H")}
    H = S["H"]
    lhs = S["TII"] + 2 * H * S["Jl"] + H**2 * S["ep"]
    rhs_terms = dict(hI=0.5 * np.gradient(np.gradient(S["I"], tau), tau), tid=-(np.gradient(H, tau) + H**2) * S["I"], W=-S["W"], Sb=S["Sb"])
    flux = S["jb"] / al                                                # energy out through r_b per unit tau (per sr)
    return t, S, lhs, rhs_terms, np.gradient(S["E"], tau), flux


def report3(path, r_bs=(4.0, 6.0, 8.0)):
    print(f"\n=== {path.split('/')[-1]}: form (3); each / (3E), window means")
    for r_b in r_bs:
        best = None
        for js in (1.0, -1.0):
            t, S, lhs, R, Edot, flux = analyse3(path, r_b, jsign=js)
            m = (t >= 13) & (t < 54) & ~np.isnan(R["hI"])
            c = float(np.corrcoef(Edot[m], -flux[m])[0, 1])
            if best is None or c > best[0]:
                best = (c, js, t, S, lhs, R, Edot, flux)
        c, js, t, S, lhs, R, Edot, flux = best
        print(f"  r_b {r_b:g}: sign of j fixed by dE/dtau vs -flux (correlation {c:+.2f}, sign {js:+g}); flux out / (H E) = {np.nanmean((flux/(S['H']*S['E']))[(t>13)&(t<54)]):+.3f}")
        print("     window   measured w(slice)  boost terms   T^ii_F    (1/2)I''   -(H'+H^2)I   -W_x     S_b      RHS      resid (T^ii_F - RHS)")
        for a, b in ((13, 30), (30, 54), (13, 54)):
            m = (t >= a) & (t < b) & ~np.isnan(R["hI"])
            f = lambda q: float(np.mean(q[m] / (3 * S["E"][m])))     # noqa: E731
            w_s = f(S["TII"]); L = f(lhs); terms = [f(R[k]) for k in ("hI", "tid", "W", "Sb")]; rhs = sum(terms)
            print(f"     {a}-{b}    {w_s:+.4f}             {L-w_s:+.4f}       {L:+.4f}   " + "   ".join(f"{x:+.4f}" for x in terms) + f"   {rhs:+.4f}   {L-rhs:+.4f}")


# ---- fourth form: exact boundary terms (comoving boundary) and chat's Amendment 1 (the large terms cancelled on paper) ----
# Exact moment identity for a region whose boundary moves at v_b = H l_b (comoving), central frame:
#   Sum T^ii_F = (1/2) I'' - (H' + H^2) I - W_x + S_b - B1 - B2
#   B1 = l_b A j_F,b v_b                       (momentum carried across the moving boundary)
#   S_b must then be the central-frame radial stress, T^rr_F = T^rr_slice + 2 H l j + (e + p) H^2 l^2; the boost part is
#   folded into B1 here (S_b stays slice-frame), leaving B1 - boost = -l_b^2 A H j_slice.
#   B2 = (1/2) d/dtau [ l_b^2 A (e_b v_b - j_F,b) ]   (energy carried across it, weighted by l_b^2)
# with j_F = j_slice + (e + p) H l (first-order boost), A = R_areal^2 (per steradian).
# Amendment 1 (chat): with R_c^2 = I/E, nu = d ln R_c / d ln a, nu' = d nu / d ln a, for a self-similar core of steady
# energy, the rest-frame pressure is
#   w_rest = (R_c H)^2 [nu' + (nu - 1)(nu + 1 + H'/H^2)] / 3 + (-W_x + S_b) / (3E)
# ((nu - 1)(nu - 1/2 - 3 w_bg / 2) in a flat background; H'/H^2 is measured here, the cell being closed).
# The slice observers move at H l, the core at nu H l, so the slice-frame pressure carries the relative flow:
#   w_slice = w_rest + (1 + w)(nu - 1)^2 (R_c H)^2 / 3

def analyse4(path, r_b, n_s=2):
    d = np.load(path, allow_pickle=True); t = d["t"]; r = d["r"]; dr = float(r[1] - r[0])
    rows = []
    for i in range(len(t)):
        o = per_snapshot(d, i, r_b)
        grr = d["grr"][i]; dV = d["dV"][i]; m = r < r_b; H = -float(d["K"][i].mean()) / 3
        s, dd, ps, pd = (SQ * d[k][i] for k in ("S", "D", "PS", "PD"))
        sr, dr_ = np.gradient(s, dr), np.gradient(dd, dr)
        j = -(ps * sr + pd * dr_) / np.sqrt(grr)
        Kk = 0.5 * (ps**2 + pd**2); G = 0.5 * (sr**2 + dr_**2) / grr
        V = 0.5 * (s**2 - dd**2) + (s**4 + 6 * s**2 * dd**2 + dd**4) / 8 + 0.5
        e = Kk + G + V; p = Kk - G / 3 - V
        ell = np.concatenate([[0.5 * dr * math.sqrt(grr[0])], 0.5 * dr * math.sqrt(grr[0]) + np.cumsum(0.5 * (np.sqrt(grr[1:]) + np.sqrt(grr[:-1])) * dr)])
        jb = int(np.searchsorted(r, r_b)) - 1; A = d["R"][i][jb]**2; lb = ell[jb]
        jF = j[jb] + (e[jb] + p[jb]) * H * lb; vb = H * lb
        o.update(H=H, Jl=float(np.sum(j * ell * dV * m)), B1=lb * A * jF * vb - lb * A * (2 * H * lb * j[jb] + (e[jb] + p[jb]) * (H * lb)**2),
                 B2in=lb**2 * A * (e[jb] * vb - jF))
        rows.append(o)
    A_ = {k: np.array([o[k] for o in rows]) for k in rows[0]}
    al = float(np.mean(A_["al"][(t > 12) & (t < 54)])); tau = t * al
    S = {k: smooth(t, A_[k], n_s) for k in ("E", "I", "TII", "W", "Sb", "ep", "Jl", "H", "B1", "B2in")}
    H = S["H"]; Hd = np.gradient(H, tau)
    TF = S["TII"] + 2 * H * S["Jl"] + H**2 * S["ep"]
    hI = 0.5 * np.gradient(np.gradient(S["I"], tau), tau); tid = -(Hd + H**2) * S["I"]
    B2 = 0.5 * np.gradient(S["B2in"], tau)
    # Amendment 1
    Rc = np.sqrt(S["I"] / S["E"]); lna = np.cumsum(np.concatenate([[0], 0.5 * (H[1:] + H[:-1]) * np.diff(tau)]))
    nu = np.gradient(np.log(Rc), tau) / H; nup = np.gradient(nu, tau) / H
    w_s = S["TII"] / (3 * S["E"])
    tug = (Rc * H)**2 * (nup + (nu - 1) * (nu + 1 + Hd / H**2)) / 3
    frame = (1 + w_s) * (nu - 1)**2 * (Rc * H)**2 / 3
    return dict(t=t, E=S["E"], w_s=w_s, TF=TF, hI=hI, tid=tid, W=-S["W"], Sb=S["Sb"], B1=S["B1"], B2=B2,
                tug=tug, frame=frame, nu=nu, nup=nup, RcH=Rc * H, q=Hd / H**2)


def report4(path, r_bs=(6.0, 8.0, 10.0)):
    print(f"\n=== {path.split('/')[-1]}")
    for r_b in r_bs:
        o = analyse4(path, r_b); t = o["t"]
        print(f"  r_b {r_b:g}")
        print("     window  R_cH   nu     nu'    H'/H^2 | exact identity /3E: T^ii_F   RHS(no flux)  -B1      -B2      RHS exact  resid | "
              "Amendment 1: tug    -W_x/3E   S_b/3E   w_rest pred  + frame   = w_slice pred   measured   meas-pred")
        for a, b in ((13, 30), (30, 54), (13, 54)):
            m = (t >= a) & (t < b) & ~np.isnan(o["hI"]) & ~np.isnan(o["nup"])
            f = lambda q: float(np.mean(q[m] / (3 * o["E"][m])))       # noqa: E731
            g = lambda q: float(np.mean(q[m]))                          # noqa: E731
            TF, hI, tid, W, Sb, B1, B2 = (f(o[k]) for k in ("TF", "hI", "tid", "W", "Sb", "B1", "B2"))
            rhs0 = hI + tid + W + Sb; rhs = rhs0 - B1 - B2
            tug, fr, ws = g(o["tug"]), g(o["frame"]), g(o["w_s"])
            wr = tug + W + Sb
            print(f"     {a}-{b}   {g(o['RcH']):.3f}  {g(o['nu']):.3f}  {g(o['nup']):+.3f}  {g(o['q']):+.3f} |       {TF:+.4f}   {rhs0:+.4f}      {-B1:+.4f}  {-B2:+.4f}  {rhs:+.4f}   {TF-rhs:+.4f} |"
                  f"              {tug:+.4f}   {W:+.4f}   {Sb:+.4f}   {wr:+.4f}      {fr:+.4f}    {wr+fr:+.4f}          {ws:+.4f}   {ws-wr-fr:+.4f}")
