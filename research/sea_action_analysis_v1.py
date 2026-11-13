"""RESEARCH ONLY. Analysis for HANDOFF to projects Nov 13 2026 (Tasks 1-3): S's action, and D split into its settled
and travelling parts. Reads the snapshots of research/sd_sea_v1.py (flat) and research/dm_one_budget_v1.py (gravity).

S's action (the new standing measure). Per cell, over a window of W time units:
    j = sqrt(<S_tau^2> <S^2>)          (time averages over the window; S_tau = derivative in local proper time)
For a harmonic oscillation S = A cos(w tau) this is A^2 w / 2 = E / w, with E the oscillation's full energy INCLUDING
whatever sets its frequency (its share of the S-D interaction), so no convention for splitting the interaction energy
enters. J = sum of j over the pocket's proper volume. For an anharmonic oscillation j differs from the exact action
integral by a few per cent at fixed shape (it is the same functional throughout, so ratios in time are what is read).
Constant J: retuning. Falling J: drain. Rising J: real gain.

D's split (Task 1). delta = D - <D> (D's oscillation about its window mean at that cell: about +D* or -D* where D has
settled, about its excluded value inside dense matter). delta and delta_t are split by spatial wavenumber (FFT of r delta, extended oddly through
the origin, in spherical geometry; of delta in the row) at k_c:
    (a) settled   : k < k_c   (group velocity k / sqrt(2 + k^2) < 0.33 at k_c = 0.5)
    (b) travelling: k >= k_c
k_c = 0.5 is the reading; k_c = 0.3 and 1.0 are reported as rulers. Energy and pressure of each part to quadratic order
about the window mean, with the local curvature m^2 = -1 + (3 lambda/2)(<D>^2 + <S^2>) (= 2, D's settled frequency^2,
at D* away from matter):
    e = delta_t^2/2 + delta'^2/2 + m^2 delta^2/2,   p = delta_t^2/2 - delta'^2/6 - m^2 delta^2/2   (isotropic average)
w = <p>/<e> over the window and region (dust: w = 0; radiation: 1/3). Anharmonic corrections are neglected.
Usage: python research/sea_action_analysis_v1.py rates475 FULL.npz NOD.npz NOMATTER.npz  (PR #475 table, proper gradient)
       python research/sea_action_analysis_v1.py flat DIR        (Tasks 1, 2 on research/sd_sea_v1.py output)
       python research/sea_action_analysis_v1.py grav DIR [tmax] (Tasks 0, 1, 3 on research/dm_one_budget_v1.py output)
"""
import glob
import math
import os
import sys

import numpy as np

LAM = 5000.0
DSTAR = math.sqrt(2 / LAM)


def windows(t, W, t0=0.0):
    edges = np.arange(t0, t[-1] - 0.75 * W + 1e-9, W)
    return [(a, a + W, (t >= a) & (t < a + W)) for a in edges if ((t >= a) & (t < a + W)).sum() >= 4]


def action(S, St, vol, mask, sel):
    """J and the window-mean S energy proxy over the cells in mask, window sel."""
    s2 = (S[sel] ** 2).mean(0); p2 = (St[sel] ** 2).mean(0)
    return float(np.sum((np.sqrt(s2 * p2) * vol)[mask]))


def lowpass(f, r, sph, kc, dr):
    """Sharp spatial low-pass at k_c. Spherical: on u = r f, extended oddly through the origin (u(-r) = -u(r)), so the
    filter keeps u(0) = 0 and does not leak into small r."""
    n = f.shape[-1]
    if sph:
        u = f * r
        ue = np.concatenate([-u[..., ::-1], u], axis=-1)
    else:
        ue = f
    F = np.fft.rfft(ue, axis=-1)
    k = 2 * math.pi * np.fft.rfftfreq(ue.shape[-1], dr)
    F[..., k >= kc] = 0.0
    out = np.fft.irfft(F, n=ue.shape[-1], axis=-1)
    return out[..., n:] / r if sph else out


def split_D(D, Dt, r, sph, dr, sel, kc, S=None, grad=None):
    """Window-averaged energy and pressure densities of the settled (a) and travelling (b) parts of D's oscillation."""
    Dm = D[sel].mean(0)
    s2 = (S[sel] ** 2).mean(0) if S is not None else 0.0
    m2 = np.clip(-1.0 + 1.5 * LAM * (Dm**2 + s2), 0.0, None)      # local curvature of V in D at the window mean
    dl = D[sel] - Dm; dlt = Dt[sel]
    a, at = lowpass(dl, r, sph, kc, dr), lowpass(dlt, r, sph, kc, dr)
    b, bt = dl - a, dlt - at
    g = grad if grad is not None else (lambda f: np.gradient(f, dr, axis=-1))
    out = {}
    for name, (x, xt) in (("a", (a, at)), ("b", (b, bt))):
        xr = g(x)
        e = 0.5 * xt**2 + 0.5 * xr**2 + 0.5 * m2 * x**2
        p = 0.5 * xt**2 - xr**2 / 6 - 0.5 * m2 * x**2
        out[name] = (e.mean(0), p.mean(0))
    return out


# ------------------------------------------------------------------------------------------------ flat runs (Tasks 1, 2)
def flat(dirn):
    files = sorted(glob.glob(os.path.join(dirn, "*.npz")))
    files = [f for f in files if "test" not in f]
    print("Flat space, lambda = 5000, cold start, R = 8, seed 1e-10; windows of 10 time units; pocket = r < R + 3.")
    print("J/J0 = S's action relative to the first window; E/E0 = S's own energy (excluding interaction) likewise.")
    for f in files:
        d = np.load(f, allow_pickle=True)
        if str(d["status"]) != "ok":
            print(f"\n{os.path.basename(f)}: {d['status']}")
            continue
        geom, bc, x_in, m_sea, dr = str(d["geom"]), str(d["bc"]), float(d["x_in"]), float(d["m_sea"]), float(d["dr"])
        t, r, dist, vol = d["t"], d["r"], d["dist"], d["wvol"]
        sph = geom == "sph"
        dr_s = float(r[1] - r[0])
        S, St, D, Dt, chi, chit = (d[k].astype(float) for k in ("S", "St", "D", "Dt", "chi", "chit"))
        sea = m_sea > 0
        pocket = dist < float(d["R"]) + 3.0
        far = (dist > 22) & (dist < 32)
        tot = d["tot"]
        Svals = []
        print(f"\n{geom} {bc:8s} x_in = {x_in}  sea m = {('%.4f' % m_sea) if sea else 'none'}  dr = {dr}"
              f"   energy: start {tot[0,1]:.4e}, end in box {tot[-1,1]:.4e}, absorbed {tot[-1,2]:.4e}, "
              f"drift {(tot[-1,1] + tot[-1,2] - tot[0,1]) / tot[0,1]:+.1e}")
        hdr = "   window   J/J0    E/E0  | settled D (k<0.5): E_pocket   E_all    w_all  | travelling: E_all  w_all | rulers k_c=0.3/1.0: settled E_all"
        if sea:
            hdr = ("   window   J/J0    E/E0  | pocket KE+GE: D       chi    | far (r 22-32): D       chi   | far uniform rms, D / chi (in D*)"
                   " | far normal modes rms: slow / fast (in D*)")
            m2 = m_sea**2; c3 = 1.5 * LAM * (float(d["Dvac"]) + float(d["chivac"]))**2
            wv, vv = np.linalg.eigh(np.array([[-1 + c3, c3], [c3, m2 + c3]]))
            print(f"  vacuum (derived): D = {float(d['Dvac']):.5f}, chi = {float(d['chivac']):.5f}, V = {float(d['Vvac']):.3e}"
                  f" ({float(d['Vvac'])*2*LAM:+.3f} x 1/(2 lambda)); normal-mode frequencies {np.sqrt(wv[0]):.4f}, {np.sqrt(wv[1]):.4f};"
                  f" slow mode (D, chi) = ({vv[0,0]:+.3f}, {vv[1,0]:+.3f})")
        print(hdr)
        J0 = E0 = None
        for a, b, sel in windows(t, 10.0):
            J = action(S, St, vol, pocket, sel)
            eS = 0.5 * St[sel]**2 + 0.5 * np.gradient(S[sel], dr_s, axis=-1)**2 + 0.5 * S[sel]**2 + LAM / 8 * S[sel]**4
            ES = float(np.sum((eS.mean(0) * vol)[pocket]))
            if J0 is None:
                J0, E0 = J, ES
            Svals.append((a, J / J0, ES / E0))
            if not sea:
                sp = split_D(D, Dt, r, sph, dr_s, sel, 0.5, S=S)
                ea, pa = sp["a"]; eb, pb = sp["b"]
                rul = [float(np.sum(split_D(D, Dt, r, sph, dr_s, sel, kc, S=S)["a"][0] * vol)) for kc in (0.3, 1.0)]
                print(f"  {a:5.0f}-{b:<4.0f} {J/J0:6.3f}  {ES/E0:6.3f}  |  {np.sum((ea*vol)[pocket]):.3e}  {np.sum(ea*vol):.3e}"
                      f"  {np.sum(pa*vol)/max(np.sum(ea*vol),1e-300):+.3f} |  {np.sum(eb*vol):.3e}  {np.sum(pb*vol)/max(np.sum(eb*vol),1e-300):+.3f}"
                      f" |  {rul[0]:.3e} / {rul[1]:.3e}")
            else:
                def kg(F, Ft, msk):
                    return float(np.sum(((0.5 * Ft[sel]**2 + 0.5 * np.gradient(F[sel], dr_s, axis=-1)**2).mean(0) * vol)[msk]))
                dD = D[sel][:, far].mean(1) - float(d["Dvac"]); dc = chi[sel][:, far].mean(1) - float(d["chivac"])
                q = np.vstack([dD, dc]).T @ vv          # normal-mode coordinates (orthonormal eigenvectors)
                rms = lambda y: float(np.sqrt(np.mean((y - y.mean())**2))) / DSTAR  # noqa: E731
                print(f"  {a:5.0f}-{b:<4.0f} {J/J0:6.3f}  {ES/E0:6.3f}  |  {kg(D, Dt, pocket):.3e}  {kg(chi, chit, pocket):.3e} |"
                      f"  {kg(D, Dt, far):.3e}  {kg(chi, chit, far):.3e} |  {rms(dD):.4f} / {rms(dc):.4f}"
                      f" |  {rms(q[:, 0]):.4f} / {rms(q[:, 1]):.4f}")
        if not sea:
            # Task 1(a): does energy or action flow from the settled part back into S? correlation of window changes
            Ea = []
            for a, b, sel in windows(t, 10.0):
                Ea.append(float(np.sum(split_D(D, Dt, r, sph, dr_s, sel, 0.5, S=S)["a"][0] * vol)))
            Sv = np.array(Svals); Ea = np.array(Ea)
            dJ, dE, dA = np.diff(Sv[:, 1]), np.diff(Sv[:, 2]), np.diff(Ea)
            cJ = np.corrcoef(dJ, dA)[0, 1] if len(dJ) > 2 else float("nan")
            cE = np.corrcoef(dE, dA)[0, 1] if len(dJ) > 2 else float("nan")
            up = dJ > 0
            print(f"  exchange test: corr(dJ_S, dE_settled) = {cJ:+.2f}, corr(dE_S, dE_settled) = {cE:+.2f}; windows where J_S rises:"
                  f" {int(up.sum())}/{len(dJ)}, in those the settled energy falls in {int((dA[up] < 0).sum())}")
            # Task 3 (flat): settled D's radial profile against S's
            sel = t >= t[-1] - 40
            ea = split_D(D, Dt, r, sph, dr_s, sel, 0.5, S=S)["a"][0]
            eSp = (0.5 * St[sel]**2 + 0.5 * S[sel]**2 + LAM / 8 * S[sel]**4).mean(0)
            pts = [0.5, 4, 8, 11, 15, 20, 30]
            print("  last 40 units, profiles (energy density): " + ", ".join(
                f"r={p}: settled {ea[np.argmin(abs(dist-p))]:.2e} / S {eSp[np.argmin(abs(dist-p))]:.2e}" for p in pts))


# ------------------------------------------------------------------------------------------- gravity runs (Tasks 0, 1, 3)
def load_grav(path, tmax=None):
    d = np.load(path)
    r = d["r"]; dr = float(r[1] - r[0]); t = d["t"]
    n = len(t) if tmax is None else int(np.searchsorted(t, tmax + 1e-9))
    d = {k: (d[k][:n] if d[k].ndim and len(d[k]) == len(t) else d[k]) for k in d.files}
    t = d["t"]
    grr = d["grr"] if "grr" in d else np.ones_like(d["S"])
    return dict(r=r, dr=dr, t=t, S=d["S"], PS=d["PS"], D=d["D"], PD=d["PD"], dV=d["dV"], grr=grr, rho=d["rho_dm"],
                alpha=d["alpha"], name=os.path.basename(path)[:-4])


def grav(dirn, W=6.0, tmax=None):
    runs = {os.path.basename(f)[:-4]: load_grav(f, tmax) for f in sorted(glob.glob(os.path.join(dirn, "*.npz")))
            if not os.path.basename(f).startswith("cfl")}
    print("Gravity on, one-budget rule, lambda = 5000 (c^2 = 304.92, H0 = 0.099); pocket = coordinate r < 11, proper volume.")
    print(f"Windows of {W:g} coordinate time units. S_tau = PS (normal-observer proper-time derivative); proper gradient"
          " S'^2 / gamma_rr (the metric correction).")
    for name, R in runs.items():
        t, r, dr = R["t"], R["r"], R["dr"]
        pocket = r < 11.0
        vol = R["dV"]
        print(f"\n{name}: t_end = {t[-1]:.1f}")
        if "nomatter" in name:
            continue
        print("   window    J/J0    E_S/E_S0  E_S+int/E_S0 | settled D (k<0.5): E_pocket   E(r<20)   w(r<20) | travelling E(r<20) w | settled E k_c 0.3/1.0")
        J0 = E0 = None
        for a, b, sel in windows(t, W):
            vm = vol[sel].mean(0)
            J = action(R["S"], R["PS"], vm, pocket, sel)
            Sp = np.gradient(R["S"][sel], dr, axis=-1)
            eS = 0.5 * R["PS"][sel]**2 + 0.5 * Sp**2 / R["grr"][sel] + 0.5 * R["S"][sel]**2 + LAM / 8 * R["S"][sel]**4
            eI = 0.75 * LAM * R["S"][sel]**2 * R["D"][sel]**2
            ES = float(np.sum((eS.mean(0) * vm)[pocket])); EI = float(np.sum((eI.mean(0) * vm)[pocket]))
            if J0 is None:
                J0, E0 = J, ES
            line = f"  {a:5.0f}-{b:<4.0f} {J/J0:7.3f}   {ES/E0:7.3f}    {(ES+EI)/E0:7.3f}    |"
            if "noD" not in name:
                inner = r < 20
                g = (lambda f: np.gradient(f, dr, axis=-1) / np.sqrt(R["grr"][sel]))
                sp = split_D(R["D"], R["PD"], r, True, dr, sel, 0.5, S=R["S"], grad=g)
                ea, pa = sp["a"]; eb, pb = sp["b"]
                rul = [float(np.sum((split_D(R["D"], R["PD"], r, True, dr, sel, kc, S=R["S"], grad=g)["a"][0] * vm)[inner])) for kc in (0.3, 1.0)]
                line += (f"  {np.sum((ea*vm)[pocket]):.3e}  {np.sum((ea*vm)[inner]):.3e}  {np.sum((pa*vm)[inner])/np.sum((ea*vm)[inner]):+.3f} |"
                         f"  {np.sum((eb*vm)[inner]):.3e} {np.sum((pb*vm)[inner])/np.sum((eb*vm)[inner]):+.3f} | {rul[0]:.3e} / {rul[1]:.3e}")
            print(line)
        # Task 3: profiles over the last 16 units: settled D, S, and dark matter relative to the no-matter control
        ctrl = None
        for c in runs:
            if "nomatter" in c and c.split("_")[0] == name.split("_")[0] and c.split("_")[-1] == name.split("_")[-1]:
                ctrl = runs[c]
        sel = t >= t[-1] - 16
        if "noD" not in name:
            g = (lambda f: np.gradient(f, dr, axis=-1) / np.sqrt(R["grr"][sel]))
            ea = split_D(R["D"], R["PD"], r, True, dr, sel, 0.5, S=R["S"], grad=g)["a"][0]
            eSp = (0.5 * R["PS"][sel]**2 + 0.5 * R["S"][sel]**2 + LAM / 8 * R["S"][sel]**4).mean(0)
            dm = None
            if ctrl is not None and len(ctrl["t"]) >= len(t):
                dm = (R["rho"][sel] / ctrl["rho"][:len(t)][sel]).mean(0) - 1
            pts = [0.5, 4, 8, 11, 14, 18, 22]
            print("  last 16 units, profiles: " + "; ".join(
                f"r={p}: settled D {ea[np.argmin(abs(r-p))]:.2e}, S {eSp[np.argmin(abs(r-p))]:.2e}"
                + (f", DM excess {dm[np.argmin(abs(r-p))]:+.3f}" if dm is not None else "") for p in pts))


if __name__ == "__main__":
    if sys.argv[1] == "rates475":
        pass
    elif sys.argv[1] == "grav":
        grav(sys.argv[2], tmax=float(sys.argv[3]) if len(sys.argv) > 3 else None)
    else:
        flat(sys.argv[2])


def rates475(full, noD, ctrl):
    """PR #475's rate table (b) recomputed with the proper gradient (S'^2 / gamma_rr). Instantaneous E_S in r < 11."""
    def es(p):
        d = np.load(p); r = d["r"]; dr = float(r[1] - r[0]); ins = r < 11.0
        g = d["grr"] if "grr" in d.files else np.ones_like(d["S"])
        e = (0.5 * d["PS"]**2 + 0.5 * np.gradient(d["S"], dr, axis=-1)**2 / g + 0.5 * d["S"]**2 + LAM / 8 * d["S"]**4)
        return d["t"], (e * d["dV"])[:, ins].sum(1), d
    t, Ef, df = es(full); _, En, _ = es(noD)
    dc = np.load(ctrl)
    n = min(len(t), len(dc["t"]), len(En))
    ins = df["r"] < 11.0
    edm = (df["dm_E"][:n][:, ins].sum(1) - dc["dm_E"][:n][:, ins].sum(1)) / (8 * math.pi)
    print(f"{os.path.basename(full)}: E_S0 = {Ef[0]:.4e}; E_S/E_S0 at t = 8/16/24: "
          + "/".join(f"{Ef[np.argmin(abs(t - x))] / Ef[0]:.3f}" for x in (8, 16, 24))
          + "; no-D: " + "/".join(f"{En[np.argmin(abs(t - x))] / En[0]:.3f}" for x in (8, 16, 24)))
    for a, b in ((0, 4), (4, 12), (12, 24)):
        m = (t[:n] >= a) & (t[:n] <= b)
        dex = np.polyfit(t[:n][m], edm[m], 1)[0]
        dS = np.polyfit(t[:n][m], (Ef[:n] - En[:n])[m], 1)[0]
        print(f"  t = {a}-{b}: DM energy in {dex:+.2e} per unit t; d/dt[E_S(full) - E_S(noD)] {dS:+.2e}; ratio {abs(dex / dS):.1f}")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "rates475":
    rates475(*sys.argv[2:5])
