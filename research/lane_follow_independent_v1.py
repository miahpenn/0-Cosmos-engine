"""RESEARCH ONLY. The lane-follow test and the evolved B6 build made independent: baryon and radiation densities, the
Lorentz factors and the clock weight are read from stored fluid primitives (research/lattice_cell_run_primitives_v1.py),
not solved from the assembled totals. Note: results/lane-follow-independent-2026-11-14/README.md.
Usage: python research/lane_follow_independent_v1.py RUN_WITH_PRIMITIVES.npz [...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from research.b6_frozen_test_v1 import KAP, LAM, P_FULL, evaluate, stored, v_sd  # noqa: E402
from research.channel_companions_v1 import slice_terms  # noqa: E402
from research.cosmos_lane_follow_v1 import lane  # noqa: E402


def rows(d):
    t = d["t"]; out = []
    for i in range(len(t)):
        dV, al = d["dV"][i], d["alpha"][i]
        av = lambda q: float(np.sum(q * dV) / np.sum(dV))                                 # noqa: E731
        s = slice_terms(d, i); f = stored(d, i)
        dm, b, r = d["p_rho_dm"][i], d["p_rho_b"][i], d["p_rho_r"][i]
        Wd, Wb, Wr = d["p_W_dm"][i], d["p_W_b"][i], d["p_W_r"][i]
        rest = dm * Wd + b * Wb                                                           # conserved rest density / sqrt(gamma)
        N_all = av(al / np.where(dm * Wd + b * Wb > 0, (dm * Wd * Wd + b * Wb * Wb) / rest, 1.0) * rest) / av(rest)
        N_dm = av(al * dm) / av(dm)
        v2 = lambda W: 1 - 1 / W**2                                                       # noqa: E731
        q_fl = 0.5 * dm * Wd**2 * (1 + v2(Wd)) + 0.5 * b * Wb**2 * (1 + v2(Wb)) + (2.0 / 3.0) * r * Wr**2 * (1 + v2(Wr)) + r / 3
        # inferred from totals, as in cosmos_lane_follow_v1 (the non-independent route)
        qf, rf = av(f["q_fl"]), av(f["rho_fl"])
        b_inf = 2 * ((rf - av(dm)) - (qf - 0.5 * av(dm))); r_inf = rf - av(dm) - b_inf
        g = dict(f); g.update(q_fl=q_fl, rho_b=b, rho_r=r, rest_dm=rest, W_dm=np.where(rest > 0, (dm * Wd * Wd + b * Wb * Wb) / rest, 1.0))
        g.pop("fluid_from_totals", None)
        e = evaluate(g)
        out.append(dict(t=float(t[i]), N_all=N_all, N_dm=N_dm, H=-s["K"] / 3, phi=s["phi"], Pi=s["Pi"], dm=av(dm), b=av(b), r=av(r),
                        sd=KAP * s["e"], lna=math.log(np.sum(dV)) / 3, b_inf=b_inf, r_inf=r_inf,
                        vel=av(q_fl) / (0.5 * av(dm) + 0.5 * av(b) + av(r)) - 1, build=e["build"],
                        miss=(e["ray"] - e["cell"]) / abs(e["cell"]), missC=(e["with_C"] - e["cell"]) / abs(e["cell"]),
                        vmax=float(np.sqrt(np.max(v2(Wd))))))
    return out


if __name__ == "__main__":
    rho_vac = KAP * v_sd(0.0, math.sqrt(2 / LAM))
    for path in sys.argv[1:]:
        d = np.load(path, allow_pickle=True); R = rows(d); t = np.array([x["t"] for x in R])
        print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, {d['start']}, r_max {d['r'][-1]+0.05:g}, closure {d['closure']}")
        print("  evolved frozen states, every part of C built from fields and fluid primitives (B6 made independent):")
        print("      t     independent build |<alpha Q> - Kdot|/Kdot   miss without channel   with channel   largest dust speed   fluid velocity part of rho+3p   N(all rest mass) - N(dark matter only)")
        for tq in (3, 6, 12, 18, 24):
            x = R[int(np.argmin(np.abs(t - tq)))]
            print(f"    {x['t']:5.1f}      {x['build']:.1e}                                 {x['miss']:+.3e}             {x['missC']:+.1e}       {x['vmax']:.2e}             {x['vel']:+.2e}"
                  f"                      {x['N_all']-x['N_dm']:+.2e}")
        i0 = int(np.argmin(np.abs(t - 12.0))); x0 = R[i0]
        print(f"  at t = 12: rho_b/rho_dm read {x0['b']/x0['dm']:.6f}, solved from totals {x0['b_inf']/x0['dm']:.6f};  rho_r/rho_dm read {x0['r']/x0['dm']:.3e}, solved {x0['r_inf']/x0['dm']:.3e}")
        for clock in ("N_all", "N_dm"):
            Nn = np.array([x[clock] for x in R])
            tau = np.concatenate([[0], np.cumsum(0.5 * (Nn[1:] + Nn[:-1]) * np.diff(t))])
            print(f"  lane integrated alone from t = 12, clock = {'all rest mass (dark matter + baryons), alpha/W' if clock == 'N_all' else 'dark matter only'}; initial rho_b, rho_r READ from primitives")
            print("    variant        at t:   dH/H(t0)     d ln a      d phi/phi-change   d pi/|pi|   d rho_dm/rho_dm   d rho_b/rho_b   d rho_r/rho_r   d rho_SD/rho_SD")
            for name, form, sdf, src in (("none", "raychaudhuri", 0.0, "read"), ("w=0", "raychaudhuri", 1.0, "read"),
                                         ("friedmann", "friedmann", 1.0, "read"), ("w=0, solved", "raychaudhuri", 1.0, "solved")):
                b0, r0 = (x0["b"], x0["r"]) if src == "read" else (x0["b_inf"], x0["r_inf"])
                Y = lane([x0["lna"], x0["H"], x0["phi"], x0["Pi"], x0["dm"], b0, r0, x0["sd"] * sdf], tau[i0:], P_FULL, rho_vac, 0.0, form)
                for tq in (18.0, 24.0):
                    i = int(np.argmin(np.abs(t - tq))); k = i - i0; x = R[i]
                    sd = (Y[k, 7] / x["sd"] - 1) if sdf else float("nan")
                    print(f"    {name:12s}  {x['t']:5.1f}:  {(Y[k,1]-x['H'])/abs(x0['H']):+.2e}   {Y[k,0]-x['lna']:+.2e}   {(Y[k,2]-x['phi'])/abs(x['phi']-x0['phi']):+.2e}"
                          f"          {(Y[k,3]-x['Pi'])/abs(x['Pi']):+.2e}   {Y[k,4]/x['dm']-1:+.2e}         {Y[k,5]/x['b']-1:+.2e}       {Y[k,6]/x['r']-1:+.2e}       {sd:+.2e}")
