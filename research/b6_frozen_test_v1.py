"""RESEARCH ONLY. B6 of HANDOFF Nov 14 2026 (5): the frozen-state test of the derived local-to-cosmos channel, on the
dust clock (decision Nov 14: the cosmos lane's t is the proper time of its own fluids, rest-mass weighted).
No time stepping. Nothing in engine/ is changed; the cosmos lane's own functions are called as they are.

    cell:    dH/dtau = -(dK/dt)/(3 N),  dK/dt from the cell's lapse solve,  N = <alpha/W> weighted by dust rest mass
    cosmos:  dH/dtau = -H^2 - (rho + 3p)/6            (Raychaudhuri form), rho, p from engine.cosmos.cosmos_rho_p at
             the proper-volume-averaged fields (+ D at rest in its vacuum, one-budget cells)
    channel: - C/3,  C = [ (3/2)<alpha Aa^2> + 8 pi <alpha (PS^2 + PD^2 - V_SD)> - vacuum + <alpha (Pi^2 - V_c)>
             + <alpha q_fluid> + <alpha> K^2/3 ] / N - Q(averaged fields)

Two kinds of frozen state:
  fresh   a fresh initial slice of the cell. Every part of C is built from the fields and the fluid primitives
          (rho, v), not from the engine's assembled rho + S. Fully independent.
  stored  a stored evolved snapshot (Nov 14 (4) runs). The snapshots hold the fields but not the baryon and radiation
          primitives, so the fluid part is taken as (assembled total) - (scalar parts). Independent for the scalars,
          the shear, the lapse weighting and the clock only. Labelled as such.
Negative controls: (1) S's action as the averaging weight; (2) the cosmos lane's -(rho + p)/2 form with -C/3 added;
(3) a pocket without D.
Usage: python research/b6_frozen_test_v1.py fresh | stored RUN.npz [RUN.npz ...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from engine.cosmos import CosmosParams, cosmos_rho_p  # noqa: E402

LAM = 5000.0
KAP = 8 * math.pi
V0_1B = KAP / (2 * LAM * 8.2425e-6) * 8.242522415500654e-5   # one-budget rule: c^2 V0, c^2 = 304.92 (as research/dm_one_budget_v1.py)
# The cells' own physics, stated (Nov 16 2026): no coupled floor (beta_D = 0 in these cells; the vacuum enters as rho_vac
# below) and imported, fully coupled dust (f_D = 1). Since #503 the defaults had silently added production's coupled
# floor (lambda = 1.57e6) inside cosmos_rho_p; since the dark-matter replacement an eight-slot call raises.
P_FULL = CosmosParams(V0=V0_1B, D=0.0, coupled_floor=False, f_D=1.0)    # hand floor off; D's ejected vacuum is the floor
P_NOD = CosmosParams(V0=V0_1B, D=0.10, coupled_floor=False, f_D=1.0)    # control without D: the hand floor restored


def v_sd(S, D):
    return 0.5 * S**2 - 0.5 * D**2 + LAM / 8 * (S**4 + 6 * S**2 * D**2 + D**4)


def evaluate(f, weight=None):
    """f: dict of arrays on one slice. Returns the test's numbers. weight: averaging weight for the cosmos lane's
    fields (None = proper volume, the derived measure)."""
    dV, al, K = f["dV"], f["alpha"], float(np.mean(f["K"]))
    av = lambda q: float(np.sum(q * dV) / np.sum(dV))                                     # noqa: E731
    w = dV if weight is None else weight * dV
    bar = lambda q: float(np.sum(q * w) / np.sum(w))                                      # noqa: E731
    has_D = bool(np.any(f["D"] != 0))
    P1B = P_FULL if has_D else P_NOD
    rho_vac = KAP * v_sd(0.0, math.sqrt(2 / LAM)) if has_D else 0.0                       # Einstein-normalised, p = -rho
    q_phi = f["Pi"]**2 - P1B.V0 * (np.exp(-P1B.lam * f["phi"]) - P1B.D)
    q_sd = KAP * (f["PS"]**2 + f["PD"]**2 - v_sd(f["S"], f["D"]))
    shear = 1.5 * f["Aa"]**2
    aQ = av(al * (shear + q_sd + q_phi + f["q_fl"])) + av(al) * K * K / 3
    wd = f["rest_dm"] * dV
    N = float(np.sum(al / f["W_dm"] * wd) / np.sum(wd))                                   # dust clock per unit label
    H = -K / 3
    # cosmos lane at the averaged fields, through its own function
    y = (1.0, H, bar(f["phi"]), bar(f["Pi"]), bar(f["rho_dm"]), bar(f["rho_b"]), bar(f["rho_r"]), 0.0)
    rho, p = cosmos_rho_p(y, P1B)
    rho += rho_vac; p -= rho_vac
    if f.get("fluid_from_totals"):                 # stored states: rho_b, rho_r are not individually known
        Vb = P1B.V0 * (math.exp(-P1B.lam * bar(f["phi"])) - P1B.D)
        src = bar(f["Pi"])**2 - Vb + bar(f["q_fl"]) - rho_vac
        rho = 0.5 * bar(f["Pi"])**2 + Vb + bar(f["rho_fl"]) + rho_vac
    else:
        src = 0.5 * (rho + 3 * p)
    C = aQ / N - (K * K / 3 + src)
    cell = -f["kdot"] / (3 * N)
    ray = -H * H - src / 3
    form2 = -0.5 * (rho + (2 * src - rho) / 3) - C / 3          # -(rho + p)/2 - C/3, with p = (2 src - rho)/3
    return dict(H=H, N=N, cell=cell, ray=ray, C=C, with_C=ray - C / 3, form2=form2, F=av(f["q_fl"]),
                build=abs(aQ - f["kdot"]) / abs(f["kdot"]), friedmann_gap=H * H - rho / 3)


def fresh(case, x_c, N, RMAX, start="cold"):
    import research.lattice_cell_v1 as lc
    from engine.matter_system import Species, primitives
    from engine.v55_matter import _metrics_from_slice, metric_slice_from_q
    k, st, r, H0 = lc.cell_initial(case, x_c, 0.0, N, RMAX, closure="ii", start=start)
    g, sc = st.geometry, st.scalars
    al, kdot = k._solve_lapse_for_state(st.grid, g, sc, st.matter, None)
    mets = _metrics_from_slice(metric_slice_from_q(st.grid, g))
    pr = {n: primitives(mets, getattr(st.matter, n), sp) for n, sp in
          (("dark_matter", Species.DARK_MATTER), ("baryons", Species.BARYON), ("radiation", Species.RADIATION))}
    rho = {n: np.array([q.rho for q in v]) for n, v in pr.items()}
    W = {n: np.array([q.lorentz() for q in v]) for n, v in pr.items()}
    v2 = {n: 1.0 - 1.0 / W[n]**2 for n in W}
    q_fl = (0.5 * rho["dark_matter"] * W["dark_matter"]**2 * (1 + v2["dark_matter"])
            + 0.5 * rho["baryons"] * W["baryons"]**2 * (1 + v2["baryons"])
            + (4.0 / 3.0) * rho["radiation"] * W["radiation"]**2 * (1 + v2["radiation"]) / 2 + rho["radiation"] / 3)
    vols = np.asarray(st.grid.volumes, float)
    return dict(dV=np.sqrt(g.a) * g.b / g.X**3 * vols, alpha=np.asarray(al), K=np.asarray(g.K), Aa=np.asarray(g.Aa),
                S=np.asarray(sc.S), PS=np.asarray(sc.PS), D=np.asarray(sc.D), PD=np.asarray(sc.PD), phi=np.asarray(sc.phi),
                Pi=np.asarray(sc.Pi), rho_dm=rho["dark_matter"], rho_b=rho["baryons"], rho_r=rho["radiation"], q_fl=q_fl,
                rest_dm=rho["dark_matter"] * W["dark_matter"], W_dm=W["dark_matter"], kdot=float(kdot))


def stored(d, i):
    P1B = P_FULL if np.any(d["D"][i] != 0) else P_NOD
    Vc = P1B.V0 * (np.exp(-d["phi"][i]) - P1B.D)
    q_phi = d["Pi"][i]**2 - Vc
    q_sd = KAP * (d["PS"][i]**2 + d["PD"][i]**2 - v_sd(d["S"][i], d["D"][i]))
    tot = 4 * math.pi * (d["rho"][i] + d["pr"][i] + 2 * d["pt"][i])
    dr = float(d["r"][1] - d["r"][0])
    grad = lambda q: np.gradient(q, dr)**2 / d["grr"][i]                                  # noqa: E731
    e_sc = (0.5 * d["Pi"][i]**2 + 0.5 * grad(d["phi"][i]) + Vc
            + KAP * (0.5 * d["PS"][i]**2 + 0.5 * d["PD"][i]**2 + 0.5 * grad(d["S"][i]) + 0.5 * grad(d["D"][i])
                     + v_sd(d["S"][i], d["D"][i])))
    return dict(dV=d["dV"][i], alpha=d["alpha"][i], K=d["K"][i], Aa=d["Aa"][i], S=d["S"][i], PS=d["PS"][i], D=d["D"][i],
                PD=d["PD"][i], phi=d["phi"][i], Pi=d["Pi"][i], rho_dm=d["rho_dm"][i], rho_b=0 * d["rho_dm"][i],
                rho_r=0 * d["rho_dm"][i], q_fl=tot - q_phi - q_sd, rho_fl=KAP * d["rho"][i] - e_sc,
                rest_dm=d["rho_dm"][i], W_dm=np.ones_like(d["rho_dm"][i]), kdot=float(d["kdot"][i]), fluid_from_totals=True)


HEAD = ("   dH/dtau: cell        cosmos (Raychaudhuri)   miss without channel   C/F       with channel: miss   independent build"
        "   | control 1: J weight, miss   | control 2: -(rho+p)/2 - C/3, miss   (Friedmann gap H^2 - rho/3)")


def line(tag, f, jw):
    e = evaluate(f); ej = evaluate(f, weight=jw)
    s = abs(e["cell"])
    ej = dict(ej, with_C=ej["ray"] - e["C"] / 3)        # the derived channel handed to a cosmos lane fed J-weighted fields
    return (f"  {tag:30s} {e['cell']:+.6e}  {e['ray']:+.6e}   {(e['ray']-e['cell'])/s:+.3e}           {e['C']/e['F']:+.4f}   "
            f"{(e['with_C']-e['cell'])/s:+.1e}           {e['build']:.1e}        |  {(ej['with_C']-e['cell'])/s:+.3e} (C_J/F {ej['C']/e['F']:+.3f})"
            f"   |  {(e['form2']-e['cell'])/s:+.3e}   ({e['friedmann_gap']:+.2e}, N = {e['N']:.5f})")


if __name__ == "__main__":
    print("B6 frozen-state test on the dust clock. Misses are (cosmos - cell)/|cell|. F = <(rho_dm + rho_b)/2 + rho_r>.")
    if sys.argv[1] == "fresh":
        print("FRESH INITIAL SLICES (fully independent build; S at its turning point, so the S-D part is its potential alone)")
        print(" " * 30 + HEAD)
        for case, xc, n, R, start in (("full", 2.0, 160, 16, "cold"), ("full", 2.0, 240, 24, "cold"), ("full", 2.0, 360, 36, "cold"),
                                      ("full", 1.5, 240, 24, "cold"), ("full", 2.0, 240, 24, "hot"), ("noD", 2.0, 240, 24, "cold"),
                                      ("nomatter", 2.0, 120, 12, "cold")):
            f = fresh(case, xc, n, R, start)
            jw = f["S"]**2 + 1e-300
            print(line(f"{case} x_c {xc} r_max {R} {start}", f, jw))
    else:
        for path in sys.argv[2:]:
            d = np.load(path, allow_pickle=True); t = d["t"]
            print(f"\nSTORED {path.split('/')[-1]} (x_c {float(d['x_c']):g}, {d['start']}, r_max {d['r'][-1]+0.05:g}, closure {d['closure']}); fluid part from the assembled totals")
            print(" " * 30 + HEAD)
            for tq in (3, 9, 15, 27, 39, 51, 63):
                sel = np.where((t >= tq - 3) & (t < tq + 3))[0]
                jw = np.sqrt((d["PS"][sel]**2).mean(0) * (d["S"][sel]**2).mean(0)) + 1e-300
                for i in (int(np.argmin(np.abs(t - tq))),):
                    print(line(f"t = {t[i]:.1f}", stored(d, i), jw))
            # window means of the misses (the S-D part alternates in sign within a cycle)
            print("   window means over 6 time units: miss without channel / with channel / control 1 / control 2 ; C/F")
            for tq in (3, 15, 39, 63):
                sel = np.where((t >= tq - 3) & (t < tq + 3))[0]
                jw = np.sqrt((d["PS"][sel]**2).mean(0) * (d["S"][sel]**2).mean(0)) + 1e-300
                E = [evaluate(stored(d, i)) for i in sel]; EJ = [evaluate(stored(d, i), weight=jw) for i in sel]
                EJ = [dict(b_, with_C=b_["ray"] - a_["C"] / 3) for a_, b_ in zip(E, EJ)]
                m = lambda g: float(np.mean([g(a, b) for a, b in zip(E, EJ)]))            # noqa: E731
                print(f"     t = {tq:2d}:  {m(lambda a, b: (a['ray']-a['cell'])/abs(a['cell'])):+.3e} / {m(lambda a, b: (a['with_C']-a['cell'])/abs(a['cell'])):+.1e}"
                      f" / {m(lambda a, b: (b['with_C']-a['cell'])/abs(a['cell'])):+.3e} / {m(lambda a, b: (a['form2']-a['cell'])/abs(a['cell'])):+.3e} ;  {m(lambda a, b: a['C']/a['F']):+.4f}")
