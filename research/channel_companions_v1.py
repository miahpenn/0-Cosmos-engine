"""RESEARCH ONLY. The companion terms of the local-to-cosmos channel: what the cell's inhomogeneity feeds into each
cosmos-lane equation, and the averaged equation of the S-D species the lane does not yet carry.
Derivation: results/channel-companions-2026-11-14/COMPANIONS_NOTE.md. Stored frozen states only; no new runs.

Averaging rule (CMC, zero shift, proper volume): d<f>/dt = <d_t f> - K (<alpha f> - <alpha><f>).
Dust clock: d tau = N dt, N = <alpha rho_dm>/<rho_dm> (dark-matter rest mass; W = 1 on the stored states).

    a:      d ln a/d tau = H <alpha>/N                                     lane: H
    phi:    d phi/d tau  = [<alpha Pi> - K Cov(alpha, phi)]/N              lane: pi
    pi:     d Pi/d tau   = [K <alpha> Pi - <alpha V'(phi)> + beta <alpha rho_dm>]/N     lane: K pi - V'(phi) + beta rho_dm
    rho_dm: d rho/d tau  = [K <alpha> rho - beta <alpha rho_dm Pi>]/N      lane: K rho - beta rho pi
    S-D:    d e/d tau    = [K <alpha> e + (K/3) <alpha S> + <alpha Aa (pr - pt)> - <j . grad alpha>]/N
(bars = proper-volume averages; K = -3H). Each identity is also checked against the centred time difference of the
stored averages.
Usage: python research/channel_companions_v1.py RUN.npz [RUN.npz ...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from research.b6_frozen_test_v1 import KAP, LAM, P_FULL, P_NOD, v_sd  # noqa: E402

BETA = -0.04


def slice_terms(d, i):
    dV, al, K = d["dV"][i], d["alpha"][i], float(d["K"][i].mean())
    av = lambda q: float(np.sum(q * dV) / np.sum(dV))                                     # noqa: E731
    P = P_FULL if np.any(d["D"][i] != 0) else P_NOD
    dr = float(d["r"][1] - d["r"][0]); grr = d["grr"][i]
    gr = lambda q: np.gradient(q, dr)                                                     # noqa: E731
    phi, Pi, rho = d["phi"][i], d["Pi"][i], d["rho_dm"][i]
    S, PS, D, PD = d["S"][i], d["PS"][i], d["D"][i], d["PD"][i]
    N = av(al * rho) / av(rho); am = av(al)
    Vp = lambda x: -P.V0 * P.lam * np.exp(-P.lam * x)                                     # noqa: E731
    vac = v_sd(0.0, math.sqrt(2 / LAM)) if np.any(D != 0) else 0.0
    g2 = (gr(S)**2 + gr(D)**2) / grr
    e = 0.5 * PS**2 + 0.5 * PD**2 + 0.5 * g2 + v_sd(S, D) - vac
    tr = 1.5 * (PS**2 + PD**2) - 0.5 * g2 - 3 * (v_sd(S, D) - vac)                       # pr + 2 pt
    jda = -(PS * gr(S) + PD * gr(D)) * gr(al) / grr
    o = dict(t=float(d["t"][i]), K=K, N=N, am=am, phi=av(phi), Pi=av(Pi), rho=av(rho), e=av(e), F=0.0)
    # cell's averaged right-hand sides per unit dust time, and the lane's at the averaged fields
    o["cell_phi"] = (av(al * Pi) - K * (av(al * phi) - am * av(phi))) / N
    o["lane_phi"] = av(Pi)
    o["phi_parts"] = dict(lapse_weight=av(al * Pi) / N - av(Pi), cov=-K * (av(al * phi) - am * av(phi)) / N)
    o["cell_Pi"] = (K * am * av(Pi) - av(al * Vp(phi)) + BETA * av(al * rho)) / N
    o["lane_Pi"] = K * av(Pi) - float(Vp(av(phi))) + BETA * av(rho)
    o["Pi_parts"] = dict(clock=K * av(Pi) * (am / N - 1), potential=-(av(al * Vp(phi)) / N - float(Vp(av(phi)))),
                         beta=BETA * (av(al * rho) / N - av(rho)))
    o["cell_rho"] = (K * am * av(rho) - BETA * av(al * rho * Pi)) / N
    o["lane_rho"] = K * av(rho) - BETA * av(rho) * av(Pi)
    o["rho_parts"] = dict(clock=K * av(rho) * (am / N - 1), beta=-BETA * (av(al * rho * Pi) / N - av(rho) * av(Pi)))
    o["a_part"] = am / N - 1
    o["cell_e"] = (K * am * av(e) + K / 3 * av(al * tr) + av(al * d["Aa"][i] * g2) - av(jda)) / N
    o["e_parts"] = dict(dilution=K * am * av(e) / N, pressure=K / 3 * av(al * tr) / N, shear=av(al * d["Aa"][i] * g2) / N,
                        lapse=-av(jda) / N)
    o["w"] = av(tr) / 3
    o["grad_share"] = av(0.5 * g2) / max(av(e), 1e-300)
    o["raw"] = dict(phi=av(phi), Pi=av(Pi), rho=av(rho), e=av(e))
    return o


if __name__ == "__main__":
    for path in sys.argv[1:]:
        d = np.load(path, allow_pickle=True); t = d["t"]
        T = [slice_terms(d, i) for i in range(len(t))]
        print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, {d['start']}, r_max {d['r'][-1]+0.05:g}, closure {d['closure']}")
        # identity check against centred differences of the stored averages (snapshots 0.5 apart)
        chk = {}
        for key, cell in (("phi", "cell_phi"), ("Pi", "cell_Pi"), ("rho", "cell_rho"), ("e", "cell_e")):
            v = np.array([x["raw"][key] for x in T]); c = np.array([x[cell] * x["N"] for x in T])
            fd = (v[2:] - v[:-2]) / (t[2:] - t[:-2])
            sm = (c[2:] + 4 * c[1:-1] + c[:-2]) / 6                           # Simpson mean of the rate over the same interval
            chk[key] = float(np.median(np.abs(fd - sm) / np.maximum(np.abs(fd), 1e-300)))
        print("  identities against time differences of the stored averages (median |difference|/|rate|): "
              + ", ".join(f"{k} {v:.1e}" for k, v in chk.items()))
        print("  companion terms per unit dust time, as fractions of the lane's own right-hand side (window means over 6)")
        print("      t   | a: <alpha>/N-1 | phi: lapse weight   cov(alpha,phi) | pi: clock      potential   beta      | rho_dm: clock   beta      |"
              " S-D species: rho_SD_E/rho_dm   w = p/rho   pressure/dilution   shear work/dil   lapse work/dil   gradient share")
        for tq in (3, 9, 15, 27, 39, 51, 63):
            sel = [x for x in T if tq - 3 <= x["t"] < tq + 3]
            m = lambda f: float(np.mean([f(x) for x in sel]))                             # noqa: E731
            lp, lP, lr = m(lambda x: abs(x["lane_phi"])), m(lambda x: abs(x["lane_Pi"])), m(lambda x: abs(x["lane_rho"]))
            dil = m(lambda x: abs(x["e_parts"]["dilution"]))
            print(f"    {tq:4.0f}  |  {m(lambda x: x['a_part']):+.2e}     |  {m(lambda x: x['phi_parts']['lapse_weight'])/lp:+.2e}       {m(lambda x: x['phi_parts']['cov'])/lp:+.2e}      |"
                  f"  {m(lambda x: x['Pi_parts']['clock'])/lP:+.2e}  {m(lambda x: x['Pi_parts']['potential'])/lP:+.2e}  {m(lambda x: x['Pi_parts']['beta'])/lP:+.1e} |"
                  f"  {m(lambda x: x['rho_parts']['clock'])/lr:+.2e}     {m(lambda x: x['rho_parts']['beta'])/lr:+.2e} |"
                  f"     {KAP*m(lambda x: x['e'])/m(lambda x: x['rho']):.4f}            {m(lambda x: x['w'])/m(lambda x: x['e']):+.4f}      {m(lambda x: x['e_parts']['pressure'])/dil:+.4f}"
                  f"             {m(lambda x: x['e_parts']['shear'])/dil:+.2e}        {m(lambda x: x['e_parts']['lapse'])/dil:+.2e}        {m(lambda x: x['grad_share']):.3f}")
