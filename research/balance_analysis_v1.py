"""RESEARCH ONLY. Analysis for HANDOFF to projects Nov 14 2026, Tasks 1 and 2 (GEAR-102 §3-§5, ledger amendment 5).

Notation (ledger amendment 5), per cell over a window of W time units:
  D at rest  : <D> (window mean), reported against D*
  slow ring  : the oscillation delta = D - <D> at spatial wavenumber k < k_c (k_c = 0.5; rulers 0.3 and 1.0)
  D_r        : the same oscillation at k >= k_c (outgoing / travelling waves)
Energies of the two ring parts to quadratic order with the local curvature (research/sea_action_analysis_v1.py):
  e = delta_t^2/2 + delta'^2/2 + m^2 delta^2/2,   m^2 = -1 + (3 lambda/2)(<D>^2 + <S^2>)
Actions: J_S = sum sqrt(<S_t^2><S^2>) dV (PR #476); J_ring = sum sqrt(<delta_t^2><delta^2>) dV, and the same for each
ring part separately.

Task 1, the balance ratio: (slow ring + D_r) energy / S's energy, and J_ring / J_S, at the plateau, over the whole
domain and in the pocket (|x - centre| < R + 3). Direction of flow: W_SD (work done on D by S through the interaction
term; > 0: S -> D) and W_DS (> 0: D -> S) per window, in total and as the fraction of the pocket's cells where it is
positive.

Task 2, the resting-point shift. Averaging D's equation over the ring at fixed S (small amplitude, D = D* + mu + delta):
  2 mu = -(3 lambda/2) D* (<delta^2> + <S^2>)    so    mu = -(3 lambda/8) D* a^2  -  (3 lambda/4) D* <S^2>,
with a^2 = 2 <delta^2> the local ring amplitude squared. The first term is GEAR-102 §3's prediction (toward the
saddle); the second is the shift matter itself causes. Reported: the measured mu = <D> - D* against both, per cell,
with a least-squares slope (measured on predicted) in cells where the prediction is small (|pred| < 0.2 D*) and D is
on its +D* or -D* branch. Implied vacuum-energy shift at the shifted resting point: V(D* + mu) - V(D*) = mu^2 (to
second order; V'' = 2 at D*), summed and set against the ring energy.
Usage: python research/balance_analysis_v1.py DIR [t_plateau]       (all *.npz from research/sd_balance_v1.py or
       research/sd_sea_v1.py without the sea)
"""
import glob
import math
import os
import sys

import numpy as np

sys.path.insert(0, ".")
from research.sea_action_analysis_v1 import lowpass, windows  # noqa: E402


def ring_parts(D, Dt, S, r, sph, dr, sel, lam, kc):
    Dm = D[sel].mean(0)
    s2 = (S[sel] ** 2).mean(0)
    m2 = np.clip(-1.0 + 1.5 * lam * (Dm**2 + s2), 0.0, None)
    dl = D[sel] - Dm; dlt = Dt[sel]
    a, at = lowpass(dl, r, sph, kc, dr), lowpass(dlt, r, sph, kc, dr)
    b, bt = dl - a, dlt - at
    g = lambda f: np.gradient(f, dr, axis=-1)  # noqa: E731
    out = {}
    for nm, (x, xt) in (("slow", (a, at)), ("Dr", (b, bt)), ("ring", (dl, dlt))):
        e = (0.5 * xt**2 + 0.5 * g(x)**2 + 0.5 * m2 * x**2).mean(0)
        j = np.sqrt((xt**2).mean(0) * (x**2).mean(0))
        out[nm] = (e, j)
    return Dm, s2, (dl**2).mean(0), out


def analyse(path, t_plat):
    d = np.load(path, allow_pickle=True)
    lam = float(d["lam"]) if "lam" in d.files else 5000.0
    Dstar = math.sqrt(2 / lam)
    geom, bc, x_in, dr0 = str(d["geom"]), str(d["bc"]), float(d["x_in"]), float(d["dr"])
    L = float(d["L"]) if "L" in d.files else (40.0 if geom == "row" else float("nan"))
    t, r, dist, vol = d["t"], d["r"], d["dist"], d["wvol"]
    sph = geom == "sph"; dr = float(r[1] - r[0])
    S, St, D, Dt = (d[k].astype(float) for k in ("S", "St", "D", "Dt"))
    pocket = dist < float(d["R"]) + 3.0
    far = dist > float(d["R"]) + 6.0
    hasW = "WSD" in d.files
    rows = []
    shift = []
    for a, b, sel in windows(t, 10.0):
        eS = (0.5 * St[sel]**2 + 0.5 * np.gradient(S[sel], dr, axis=-1)**2 + 0.5 * S[sel]**2 + lam / 8 * S[sel]**4).mean(0)
        JS = np.sqrt((St[sel]**2).mean(0) * (S[sel]**2).mean(0))
        Dm, s2, d2, parts = ring_parts(D, Dt, S, r, sph, dr, sel, lam, 0.5)
        row = dict(t=a, ES=np.sum(eS * vol), ESp=np.sum((eS * vol)[pocket]), JS=np.sum(JS * vol), JSp=np.sum((JS * vol)[pocket]))
        for nm, (e, j) in parts.items():
            row["E_" + nm] = np.sum(e * vol); row["Ep_" + nm] = np.sum((e * vol)[pocket])
            row["J_" + nm] = np.sum(j * vol); row["Jp_" + nm] = np.sum((j * vol)[pocket])
        for kc in (0.3, 1.0):
            _, _, _, pk = ring_parts(D, Dt, S, r, sph, dr, sel, lam, kc)
            row[f"E_slow_{kc}"] = np.sum(pk["slow"][0] * vol)
        if hasW:
            WSD, WDS = d["WSD"][sel].astype(float).sum(0), d["WDS"][sel].astype(float).sum(0)
            row["WSD"] = np.sum(WSD * vol) / (b - a); row["WDS"] = np.sum(WDS * vol) / (b - a)
            row["fSD"] = float(np.sum(vol[pocket] * (WSD[pocket] > 0)) / np.sum(vol[pocket]))
            row["fDS"] = float(np.sum(vol[pocket] * (WDS[pocket] > 0)) / np.sum(vol[pocket]))
        rows.append(row)
        if a >= t_plat:
            branch = np.abs(Dm) > 0.5 * Dstar
            sgn = np.sign(Dm)
            mu = (np.abs(Dm) - Dstar) / Dstar                      # measured shift toward the saddle is negative
            p_ring = -(3 * lam / 8) * 2 * d2                         # -(3 lambda/8) D* a^2, in units of D*  (a^2 = 2<delta^2>)
            p_S = -(3 * lam / 4) * s2                                 # -(3 lambda/4) <S^2>, in D*
            shift.append((mu, p_ring, p_S, branch, sgn, vol, far, pocket, d2, Dm))
    return dict(name=os.path.basename(path)[:-4], lam=lam, geom=geom, bc=bc, x_in=x_in, L=L, dr=dr0, rows=rows,
                shift=shift, Dstar=Dstar, t_end=float(t[-1]))


def report(res, t_plat):
    rows = [r for r in res["rows"] if r["t"] >= t_plat]
    print(f"\n{res['name']}: {res['geom']} {res['bc']} x_in = {res['x_in']} lambda = {res['lam']:g} L = {res['L']:g} dr = {res['dr']}"
          f"  (plateau windows t >= {t_plat:g}: {len(rows)})")
    if not rows:
        return None
    J0 = res["rows"][0]["JS"]; E0 = res["rows"][0]["ES"]
    f = lambda k: np.array([r[k] for r in rows])  # noqa: E731
    JS, ES = f("JS"), f("ES")
    ring_E = f("E_ring"); slowE, DrE = f("E_slow"), f("E_Dr")
    ringJ = f("J_ring")
    ms = lambda y: f"{y.mean():.3f} ± {y.std():.3f}"  # noqa: E731
    Jp0 = res["rows"][0]["JSp"]
    print(f"  S: J/J0 {ms(JS / J0)} (whole domain), {ms(f('JSp') / Jp0)} (pocket); E/E0 {ms(ES / E0)}")
    print(f"  ring energy / S energy (whole domain): {ms(ring_E / ES)}   [slow ring {ms(slowE / ES)}, D_r {ms(DrE / ES)}]")
    print(f"  ring action / S action (whole domain): {ms(ringJ / JS)}   [slow ring {ms(f('J_slow') / JS)}, D_r {ms(f('J_Dr') / JS)}]")
    print(f"  in the pocket: ring E / S E {ms(f('Ep_ring') / f('ESp'))}, ring J / S J {ms(f('Jp_ring') / f('JSp'))}")
    print(f"  slow-ring energy at k_c = 0.3 / 0.5 / 1.0 (relative to S's): {ms(f('E_slow_0.3') / ES)} / {ms(slowE / ES)} / {ms(f('E_slow_1.0') / ES)}")
    out = dict(JSp=(f("JSp") / Jp0).mean(), JS=(JS / J0).mean(), ratioE=(ring_E / ES).mean(), ratioJ=(ringJ / JS).mean())
    if "WSD" in rows[0]:
        WSD, WDS = f("WSD"), f("WDS")
        print(f"  flow per unit t, in units of S's energy (whole domain): W_SD (S->D) {(WSD / ES).mean():+.2e} ± {(WSD / ES).std():.2e},"
              f" W_DS (D->S) {(WDS / ES).mean():+.2e} ± {(WDS / ES).std():.2e};"
              f" windows with W_SD > 0: {int((WSD > 0).sum())}/{len(WSD)}, with W_DS > 0: {int((WDS > 0).sum())}/{len(WDS)}")
        print(f"  pocket cells with net S->D in a window: {ms(f('fSD'))}; with net D->S: {ms(f('fDS'))}"
              f"; mean |W_SD| / (S energy) per unit t: {np.mean(np.abs(WSD) / ES):.2e}")
        out.update(WSD=WSD.mean(), WDS=WDS.mean(), fSD=f("fSD").mean())
    return out


def report_shift(res):
    if not res["shift"]:
        return
    mu = np.concatenate([s[0] for s in res["shift"]]); pr = np.concatenate([s[1] for s in res["shift"]])
    pS = np.concatenate([s[2] for s in res["shift"]]); br = np.concatenate([s[3] for s in res["shift"]])
    vol = np.concatenate([s[5] for s in res["shift"]]); far = np.concatenate([s[6] for s in res["shift"]])
    pocket = np.concatenate([s[7] for s in res["shift"]]); d2 = np.concatenate([s[8] for s in res["shift"]])
    ok = br & (np.abs(pr + pS) < 0.2)
    for name, m in (("far from matter (S ~ 0)", ok & far), ("pocket", ok & pocket), ("all cells", ok)):
        if m.sum() < 5:
            print(f"  shift, {name}: too few cells on a branch with a small prediction")
            continue
        P = (pr + pS)[m]; Pr = pr[m]; M = mu[m]; w = vol[m]
        slope = np.sum(w * P * M) / np.sum(w * P * P)
        slope_r = np.sum(w * Pr * M) / np.sum(w * Pr * Pr) if np.any(Pr != 0) else float("nan")
        resid = np.sqrt(np.sum(w * (M - P) ** 2) / np.sum(w))
        print(f"  shift, {name}: cells {int(m.sum())}; mean measured mu {np.average(M, weights=w):+.2e} D*,"
              f" predicted ring term {np.average(Pr, weights=w):+.2e}, S term {np.average(pS[m], weights=w):+.2e};"
              f" slope measured/predicted (ring + S) {slope:.3f}, (ring only) {slope_r:.3f}; rms residual {resid:.1e} D*")
    # implied vacuum-energy shift: V(D* + mu) - V(D*) = mu^2 (units of D*^2 -> energy D*^2 mu^2), against the ring energy
    m = br & far
    if m.sum() > 5:
        Ds = res["Dstar"]
        dV = np.sum((Ds * mu[m]) ** 2 * vol[m]); ringE = np.sum(2 * d2[m] * vol[m])   # ring energy ~ m^2 <delta^2> = 2<delta^2>
        print(f"  implied vacuum-energy shift far from matter: sum mu^2 D*^2 dV = {dV:.3e} against ring energy {ringE:.3e}"
              f" (ratio {dV / ringE:.2e}); per volume {dV / np.sum(vol[m]):.2e} against the vacuum depth 1/(2 lambda) = {1/(2*res['lam']):.2e}")


if __name__ == "__main__":
    dirn = sys.argv[1]
    t_plat = float(sys.argv[2]) if len(sys.argv) > 2 else 200.0
    print("Balance ratio and resting-point shift (flat, no phi / no vacuum); windows of 10 time units.")
    for f in sorted(glob.glob(os.path.join(dirn, "*.npz"))):
        if "test" in os.path.basename(f):
            continue
        d = np.load(f, allow_pickle=True)
        if str(d["status"]) != "ok" or float(d["m_sea"]) > 0:
            continue
        res = analyse(f, t_plat)
        report(res, t_plat)
        report_shift(res)
