"""RESEARCH ONLY. Two-component closure in a research copy of the cosmos lane (HANDOFF Nov 15 2026 (2), Task 1).
Note: results/two-component-lane-2026-11-15/README.md. Stored states and the short fresh runs of PR #485; no new runs.

Lane (Raychaudhuri form, the fluids' proper time):
    dH/dtau = -H^2 - (rho + 3p)/6,  rho + 3p = 2 pi^2 - 2V + rho_dm + rho_b + 2 rho_r - 2 rho_vac + rho_m - 2 rho_h
    rho_m: S-D matter, w = 0;  rho_h: the held vacuum step, w = -1;  exchange Q = -d rho_h/dtau (into rho_m):
    d rho_m/dtau = -3H rho_m - d rho_h/dtau.
Release law (Derived, EOS_NOTE section 2): where S holds D_s out, D's slow part obeys (scaled, mean field over S's
period, s^2 averaged to x/2)
    d'' + 3H d' = d (1 - (3/4) x) - d^3/2,     x = 2 lambda <S^2> = the local concentration, diluting as a^-3,
held while x > 4/3 (the averaged band edge; Floquet gives 1.28), released when x falls below it: D grows at
gamma = sqrt(1 - 3x/4) from its value then, and the held energy h(d) = (1/2)(1 - d^2/2)^2 goes into S-D matter.
Inputs read from the cell at t0, shell by shell: D, its rate, x, and the shells' volume weights. Nothing fitted.
Usage: python research/two_component_lane_v1.py law RUN.npz ... | lane RUN.npz ...
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

from research.b6_frozen_test_v1 import KAP, LAM, P_FULL, v_sd  # noqa: E402
from research.channel_companions_v1 import BETA  # noqa: E402

RHO_VAC = KAP * v_sd(0.0, math.sqrt(2 / LAM))
SQ = math.sqrt(LAM)


def h_of(d):
    return 0.5 * (1 - 0.5 * d * d)**2


def shells(d, i0):
    """Per-shell inputs at snapshot i0 (scaled): d, d', x, proper-volume weight."""
    s, ps = SQ * d["S"][i0], SQ * d["PS"][i0]
    dd, pd = SQ * d["D"][i0], SQ * d["PD"][i0]
    e_sad = 0.5 * ps**2 + 0.5 * s**2 + s**4 / 8                       # S's energy if D were held out
    xA = -2 + 2 * np.sqrt(1 + 2 * e_sad)
    w = d["dV"][i0] / d["dV"][i0].sum()
    return dd, pd * d["alpha"][i0], 0.97 * xA, w                       # x = 2<s^2> = 0.97 x_A in the band


def held_measured(d):
    """Cell-mean held energy (scaled), from the slow part of D (Hann over 4.5 time units, edge-normalised)."""
    D = SQ * d["D"]; k = 9; wv = np.hanning(k + 2)[1:-1]
    num = np.array([np.convolve(D[:, j], wv, mode="same") for j in range(D.shape[1])]).T
    den = np.convolve(np.ones(D.shape[0]), wv, mode="same")[:, None]
    Db = num / den
    return np.array([np.sum(h_of(Db[i]) * d["dV"][i]) / np.sum(d["dV"][i]) for i in range(len(d["t"]))])


def averages(d, fresh):
    if fresh:
        from research.lane_follow_independent_v1 import rows
        R = rows(d)
        A = np.array([[x["t"], x["N_all"], x["H"], x["phi"], x["Pi"], x["dm"], x["b"], x["r"], x["sd"], 0.0, x["lna"]] for x in R])
    else:
        from research.cosmos_lane_follow_v1 import averages as av
        A, _ = av(d)
    tau = np.concatenate([[0], np.cumsum(0.5 * (A[1:, 1] + A[:-1, 1]) * np.diff(A[:, 0]))])
    return A, tau


def lane(y0, sh, taus, variant):
    """y0 = [lna, H, phi, pi, dm, b, r, rho_m]; sh = shell inputs (d, d', x, w) or None. variant: main | no_step |
    held_constant | dust_throughout | none."""
    P = P_FULL
    d0, v0, x0, w = sh
    act = h_of(d0) > 1e-6                                                # shells that hold any step
    dA, vA, xA, wA = d0[act], v0[act], x0[act], w[act]
    n = int(act.sum()); lna0 = y0[0]
    h0 = KAP / LAM * float(np.sum(wA * h_of(dA)))
    # S-D matter starts as the cell's S-D energy less the held step; "dust_throughout" carries the step as dust;
    # "no_step" leaves the step's energy out altogether.
    m0 = y0[7] if variant == "dust_throughout" else y0[7] - h0
    def rhs(_, y):
        lna, H, phi, pi, dm, b, r, rm = y[:8]
        dd, vv, hm = y[8:8 + n], y[8 + n:8 + 2 * n], y[8 + 2 * n:8 + 3 * n]
        if variant == "main":
            x = xA * np.exp(-3 * (lna - lna0))
            acc = dd * (1 - 0.75 * x) - 0.5 * dd**3 - 3 * H * vv
            dh = np.minimum(0.0, -(1 - 0.5 * dd * dd) * dd * vv) * (hm > 0)   # one-way: the held step only fades
        else:
            acc = np.zeros(n); vv = np.zeros(n); dh = np.zeros(n)
        rh = {"main": KAP / LAM * float(np.sum(wA * hm)), "held_constant": h0}.get(variant, 0.0)
        V = P.V0 * (math.exp(-P.lam * phi) - P.D); Vp = -P.V0 * P.lam * math.exp(-P.lam * phi)
        r3p = 2 * pi * pi - 2 * V + dm + b + 2 * r - 2 * RHO_VAC + rm - 2 * rh
        dH = -H * H - r3p / 6
        Q = -KAP / LAM * float(np.sum(wA * dh))                            # = -d rho_h/dtau >= 0, into S-D matter
        return [H, dH, pi, -3 * H * pi - Vp + BETA * dm, -3 * H * dm - BETA * dm * pi, -3 * H * b, -4 * H * r,
                -3 * H * rm + Q] + list(vv) + list(acc) + list(dh)
    yy = list(y0[:7]) + [m0] + list(dA) + list(vA) + list(h_of(dA))
    s = solve_ivp(rhs, (taus[0], taus[-1]), yy, t_eval=taus, method="DOP853", rtol=1e-10, atol=1e-16)
    Y = s.y.T
    rh = np.array([KAP / LAM * float(np.sum(wA * np.maximum(Y[k, 8 + 2 * n:8 + 3 * n], 0))) if variant == "main"
                   else (h0 if variant == "held_constant" else 0.0) for k in range(len(Y))])
    return Y[:, :8], rh, h0


def law_grad(d, A, tau, i0=0):
    """The same release law with D's gradient term kept: shells on the areal radius at t0, Laplacian scaled by
    (a0/a)^2 as the lane expands (a sub-grid pocket model, not a cosmos-lane closure)."""
    d0, v0, x0, w = shells(d, i0)
    R = d["R"][i0]; n = len(R)
    Hs = np.interp(tau, tau, A[:, 2]); lnas = A[:, 10]
    def lap(f):
        fr = np.gradient(f, R); return np.gradient(fr, R) + 2 * fr / np.maximum(R, 1e-6)
    def rhs(tt, y):
        dd, vv = y[:n], y[n:]
        H = float(np.interp(tt, tau, A[:, 2])); lna = float(np.interp(tt, tau, lnas))
        sc = math.exp(-2 * (lna - lnas[i0])); x = x0 * np.exp(-3 * (lna - lnas[i0]))
        acc = sc * lap(dd) + dd * (1 - 0.75 * x) - 0.5 * dd**3 - 3 * H * vv
        acc[-1] = 0.0; vv = vv.copy(); vv[-1] = 0.0
        return np.concatenate([vv, acc])
    sol = solve_ivp(rhs, (tau[i0], tau[-1]), np.concatenate([d0, v0]), t_eval=tau[i0:], method="RK45", rtol=1e-7, atol=1e-10, max_step=0.05)
    del Hs
    return np.array([KAP / LAM * float(np.sum(w * h_of(sol.y[:n, k]))) for k in range(sol.y.shape[1])])


def law(path):
    d = np.load(path, allow_pickle=True); fresh = "p_rho_dm" in d.files
    A, tau = averages(d, fresh); t = A[:, 0]
    meas = KAP / LAM * held_measured(d)
    rg = law_grad(d, A, tau)
    sh = shells(d, 0)
    Y, rh, h0 = lane([A[0, 10], A[0, 2], A[0, 3], A[0, 4], A[0, 5], A[0, 6], A[0, 7], A[0, 8]], sh, tau, "main")
    print(f"\n=== {path.split('/')[-1]}: release law from t = 0 (shell inputs read at t = 0; lane's own a and H)")
    print("     t     held step / S-D energy: cell (slow D)   law (no gradient)   law with D's gradient | held / its t=0.5 value: cell   law   law with gradient")
    for tq in (0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15):
        i = int(np.argmin(np.abs(t - tq)))
        print(f"   {t[i]:5.1f}       {meas[i]/A[i,8]:.4f}                       {rh[i]/A[i,8]:.4f}              {rg[i]/A[i,8]:.4f}                 |  {meas[i]/meas[1]:.4f}   {rh[i]/rh[1]:.4f}   {rg[i]/rg[1]:.4f}")
    half = lambda y: float(np.interp(0.5, (y / y[1])[::-1], t[::-1]))  # noqa: E731
    print(f"   half-release time (from t = 0.5): cell {half(meas):.2f}, law {half(rh):.2f}, law with gradient {half(rg):.2f};  step's share of S-D energy at t = 0.5: cell {meas[1]/A[1,8]:.3f}, law {rh[1]/A[1,8]:.3f}")


def lane_test(path):
    d = np.load(path, allow_pickle=True); fresh = "p_rho_dm" in d.files
    A, tau = averages(d, fresh); t = A[:, 0]
    tag = "fresh (densities read)" if fresh else "stored (rho_b, rho_r solved from totals)"
    print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, r_max {d['r'][-1]+0.05:g}, {tag}")
    for t0 in (0.0, 12.0):
        i0 = int(np.argmin(np.abs(t - t0)))
        sh = shells(d, i0)
        y0 = [A[i0, 10], A[i0, 2], A[i0, 3], A[i0, 4], A[i0, 5], A[i0, 6], A[i0, 7], A[i0, 8]]
        print(f"  started at t = {t0:g}: held step {KAP/LAM*float(np.sum(sh[3]*h_of(sh[0])))/A[i0,8]:.3f} of the S-D energy")
        print("    variant            at t:    dH/H(t0)    d ln a      d rho_dm/rho_dm   d rho_SD/rho_SD (matter + step)")
        for var in ("main", "no_step", "held_constant", "dust_throughout", "none"):
            yv = list(y0)
            if var == "none":
                yv[7] = 0.0
            Y, rh, _ = lane(yv, sh, tau[i0:], var if var != "none" else "no_step")
            if var == "none":
                Y, rh, _ = lane(yv, sh, tau[i0:], "dust_throughout")
            for tq in (6.0, 12.0, 24.0, 39.0, 54.0, 66.0):
                if tq <= t0 or tq > t[-1] + 1e-9:
                    continue
                i = int(np.argmin(np.abs(t - tq))); k = i - i0
                sd = (Y[k, 7] + rh[k]) / A[i, 8] - 1 if var != "none" else float("nan")
                print(f"    {var:16s}  {t[i]:5.1f}:  {(Y[k,1]-A[i,2])/abs(A[i0,2]):+.2e}   {Y[k,0]-A[i,10]:+.2e}   {Y[k,4]/A[i,5]-1:+.2e}         {sd:+.2e}")


if __name__ == "__main__":
    fn = {"law": law, "lane": lane_test}[sys.argv[1]]
    for p in sys.argv[2:]:
        fn(p)
