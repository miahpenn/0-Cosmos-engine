"""B7 W1 outer-cell C = Lambda - Lambda_metric trace (research only, no engine change).

Gates fixed before the run (Lumen2 DERIVATION 1SN10MQ6l0K4DaWt0eFRx7QJsRuUboaPJvsrGQTZbFYA with Amendments
1ju8hVnEDtIIrzAxH0lQu1QvsFGQJ1gLMEaMMz7dDU8M and 1tJGspPs8L4RF6sv5LWws-4sx6QvViRALAc00jSf_5PY; go from
Jeremiah 2026-11-23 22:47 UTC). Runs on Aster's draft initializer branch aster/w1-profile (cfae43c), the
production kernel step is unchanged. Run from that checkout's root:

    OPENBLAS_NUM_THREADS=1 python -m research.b7_w1_outer_cell_trace OUTDIR

References (all 60-digit mpmath, outer window of the last 8 cells):
  dep  : departure side-copy. a-1, b-1, X-1 start from the initializer's own source departure integrated with
         exact eight-node weights, then take each stage's increment from the engine's binary64 rates. It never
         feeds back. Lambda is the engine's stored Lambda (not near 1, no storage cancellation).
  ctl  : arithmetic-only control. Exact-weight hp arithmetic on the stored binary64 a, b.
"""
from __future__ import annotations

import json
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np

from engine import v55_initial as vi, v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy

mp.mp.dps = 60
EPS = 2.0 ** -52
R_MAX = 60.0
import os
HORIZON = float(os.environ.get("TRACE_HORIZON", "0.003"))
W = 8  # outer window
G0_EXPECTED = {  # Aster, production-kernel replay tables, peak |M| over last three cells, units 1e-12
    ("W1", 1200, 0.00075): 0.522, ("W1", 1200, 0.000375): 0.577,
    ("PC2", 1200, 0.00075): 0.141, ("PC2", 1200, 0.000375): 1.087,
    ("W1", 2400, 0.00075): 1.507, ("W1", 2400, 0.000375): 1.296,
    ("PC2", 2400, 0.00075): 0.986, ("PC2", 2400, 0.000375): 1.977,
}
S_OUT = Fr(32, 3)
W_OUT = [Fr(c, 12) for c in (3, -16, 36, -48, 25)]      # cell n-1, nodes n-5..n-1
W_PEN = [Fr(c, 12) for c in (-1, 6, -18, 10, 3)]        # cell n-2, nodes n-5..n-1
W_CEN = [Fr(c, 12) for c in (1, -8, 0, 8, -1)]          # cell i, nodes i-2..i+2

go, vacuum, _ = ad.vendor_modules()
import bssn_characteristic_boundary as bcb  # noqa: E402  (vendor path set by vendor_modules)


def F_floor(dr):
    return float(Fr(3, 2) * S_OUT) * EPS / dr


def mpf(x):
    return mp.mpf(float(x))


# ---------------------------------------------------------------- hp window operators
def D1(vals, k, dr):
    """Exact-weight first derivative at window index k (window = last W cells, k=W-1 is the outer cell)."""
    if k == W - 1:
        w, nodes = W_OUT, range(W - 5, W)
    elif k == W - 2:
        w, nodes = W_PEN, range(W - 5, W)
    else:
        w, nodes = W_CEN, range(k - 2, k + 3)
    s = mp.mpf(0)
    for c, j in zip(w, nodes):
        s += mp.mpf(c.numerator) / c.denominator * vals[j]
    return s / dr


def lam_metric(qa, qb, r, k, dr):
    a, b = 1 + qa[k], 1 + qb[k]
    return D1(qa, k, dr) / (2 * a * a) - D1(qb, k, dr) / (a * b) + 2 * (qa[k] - qb[k]) / (a * b) / r[k]


def hp_write(qa, qb, qX, g, rho, j, r, dr):
    """Vendor apply_light_constraint_boundary in hp, metric from departures, other operands binary64."""
    a = [1 + q for q in qa]; b = [1 + q for q in qb]; X = [1 + q for q in qX]
    Aa = [mpf(v) for v in g.Aa[-W:]]; K = [mpf(v) for v in g.K[-W:]]; L = [mpf(v) for v in g.Lambda[-W:]]

    def dd(k):
        return (D1(qa, k, dr) / (2 * a[k]), D1(qb, k, dr) / (2 * b[k]), -D1(qX, k, dr) / (2 * X[k]))

    def omega(k):
        da, db, chi = dd(k)
        return mp.sqrt(a[k]) / X[k] * (Aa[k] - 2 * K[k] / 3) - mp.mpf(2) / 3 * (da - db - a[k] * L[k] + 2 * chi)

    def source(k):
        da, db, chi = dd(k)
        rr = r[k]; sa = mp.sqrt(a[k]); pi = mp.pi
        return (a[k] * Aa[k] ** 2 / (4 * X[k] ** 2) - a[k] * K[k] ** 2 / (9 * X[k] ** 2)
                + sa / X[k] * (Aa[k] * (da - 3 * db - 4 * chi - 3 / rr) - K[k] * (mp.mpf(2) / 3 * da + mp.mpf(4) / 3 * chi))
                + a[k] * L[k] * da / 3
                + 2 * a[k] * da / (b[k] * rr) - 8 * a[k] * db / (3 * b[k] * rr) - 5 * a[k] / (3 * b[k] * rr ** 2)
                + mp.mpf(4) / 3 * chi ** 2 - mp.mpf(4) / 3 * chi * da + mp.mpf(8) / 3 * chi * db + 8 * chi / (3 * rr)
                - da ** 2 / 3 + db ** 2 + 2 * db / rr + 5 / (3 * rr ** 2)
                + 8 * pi * sa * mpf(j[k - W]) / X[k] + 8 * pi * a[k] * mpf(rho[k - W]) / (3 * X[k] ** 2))

    target = omega(W - 3) + 2 * dr * source(W - 2)
    k = W - 1
    da, db, chi = dd(k)
    return 2 * K[k] / 3 + X[k] / mp.sqrt(a[k]) * (target + mp.mpf(2) / 3 * (da - db - a[k] * L[k] + 2 * chi))


# ---------------------------------------------------------------- setup departure (eight-node, exact weights)
_EIGHT = {}


def eight_weights(offset, left):
    key = (offset, left)
    if key not in _EIGHT:
        nodes = [Fr(i - offset) + Fr(1, 2) for i in range(8)]
        mom = [(Fr(1, 2) ** p - Fr(left) ** p) / p for p in range(1, 9)]
        M = [[nodes[i] ** p for i in range(8)] + [mom[p]] for p in range(8)]
        for c in range(8):
            piv = next(rw for rw in range(c, 8) if M[rw][c] != 0)
            M[c], M[piv] = M[piv], M[c]
            for rw in range(8):
                if rw != c and M[rw][c] != 0:
                    f = M[rw][c] / M[c][c]
                    M[rw] = [x - f * y for x, y in zip(M[rw], M[c])]
        sol = [M[i][8] / M[i][i] for i in range(8)]
        _EIGHT[key] = [mp.mpf(x.numerator) / x.denominator for x in sol]
    return _EIGHT[key]


def setup_departures(ops, n):
    r = ops["r"]; dr = mpf(ops["dr"])
    rr = [mpf(x) for x in r]
    dep = [-8 * mp.pi * rr[i] ** 2 * mpf(ops["rho"][i]) + 2 * rr[i] ** 2 * mpf(ops["Kr"][i]) * mpf(ops["Kt"][i])
           + rr[i] ** 2 * mpf(ops["Kt"][i]) ** 2 for i in range(n)]
    cells = np.arange(n); starts = np.minimum(cells - 4, n - 8)
    acc = mp.mpf(0); qB = {}
    for i in range(n):
        idx = starts[i] + np.arange(8)
        idx = np.where(idx < 0, -1 - idx, idx)
        w = eight_weights(int(i - starts[i]), 0 if i == 0 else Fr(-1, 2))
        acc += dr * mp.fsum(wk * dep[jk] for wk, jk in zip(w, idx))
        if i >= n - W:
            qB[i] = acc / rr[i]
    out = {}
    for name, p in (("a", mp.mpf(-2) / 3), ("b", mp.mpf(1) / 3), ("X", mp.mpf(1) / 6)):
        out[name] = [mp.expm1(p * mp.log1p(qB[i])) for i in range(n - W, n)]
    out["qB"] = [qB[i] for i in range(n - W, n)]
    return out


# ---------------------------------------------------------------- instrumented run
def run(fixture, n, dt):
    K = V55ProductionKernel()
    grid = go.SphericalCellGrid(n, R_MAX)
    dr64 = float(grid.dr); dr = mpf(dr64)
    r = [mpf(x) for x in grid.centers[-W:]]
    FF = F_floor(dr64)

    ops = {}
    orig_int = vi._eight_node_origin_integral

    def hook(grid_, r_, dr_, f):
        loc = sys._getframe(1).f_locals
        ops.update(r=np.array(r_), dr=dr_, rho=np.array(loc["rho"]), Kr=np.array(loc["Kr"]), Kt=np.array(loc["Kt"]),
                   f=np.array(f))
        return orig_int(grid_, r_, dr_, f)

    vi._eight_node_origin_integral = hook
    try:
        if fixture == "W1":
            d = np.load(Path("results/research/b7_w1_initializer_draft/native.npz"))
            tag = f"dr{60 / n:g}"
            p = vi.W1ArealProfile(grid.centers, **{k: d[tag + "_" + k].copy() for k in ("S", "PS", "D", "PD")})
            state = K.initialize(resolution=n, r_max=R_MAX, outer_edge="open_edge", initial_structure="w1_research",
                                 w1_profile=p, eight_node_areal=True)
        else:
            state = K.initialize(resolution=n, r_max=R_MAX, outer_edge="open_edge",
                                 initial_structure="rescaled_pocket", pocket_x=.5, eight_node_areal=True)
    finally:
        vi._eight_node_origin_integral = orig_int
    assert ops, "eight-node hook did not fire"
    dep0 = setup_departures(ops, n)
    g = state.geometry
    grid = state.grid
    sanity = dict(
        a_minus_dep=max(abs(float(mpf(g.a[n - W + k]) - 1 - dep0["a"][k])) for k in range(W)),
        b_minus_dep=max(abs(float(mpf(g.b[n - W + k]) - 1 - dep0["b"][k])) for k in range(W)),
        X_minus_dep=max(abs(float(mpf(g.X[n - W + k]) - 1 - dep0["X"][k])) for k in range(W)),
    )

    def C64(geo):
        from engine.native_geometry.geometry import metric_connection_lambda
        return geo.Lambda[-3:] - metric_connection_lambda(grid, geo)[-3:]

    def Cdep(geo, qa, qb):
        return [mpf(geo.Lambda[n - W + k]) - lam_metric(qa, qb, r, k, dr) for k in range(W - 3, W)]

    def Cctl(geo):
        qa = [mpf(v) - 1 for v in geo.a[-W:]]; qb = [mpf(v) - 1 for v in geo.b[-W:]]
        return [mpf(geo.Lambda[n - W + k]) - lam_metric(qa, qb, r, k, dr) for k in range(W - 3, W)]

    def mom(geo, j):
        return vacuum.constraints(grid, geo)["momentum"] - 8 * np.pi * j

    def mom_state(s):
        tot = assemble_total_stress_energy(s.grid, s.geometry, s.scalars, s.matter)
        return mom(s.geometry, tot.j)

    rec = dict(fixture=fixture, n=n, dt=dt, dr=dr64, F=FF, sanity=sanity, stages=[], steps=[])
    # setup
    c64 = C64(g); cd = Cdep(g, dep0["a"], dep0["b"]); cc = Cctl(g)
    rec["setup"] = dict(C64=[float(x) for x in c64], Cdep=[float(x) for x in cd], Cctl=[float(x) for x in cc],
                        C64_minus_Cdep=[float(mpf(c64[i]) - cd[i]) for i in range(3)],
                        C64_minus_Cctl=[float(mpf(c64[i]) - cc[i]) for i in range(3)])
    peak = float(np.max(np.abs(mom_state(state)[-3:])))

    # instrumentation
    cap = {}
    orig_gst = K._geometry_stage_terms
    orig_l2 = vacuum.primary_l2_rhs
    orig_ll2 = vacuum.lambda_l2_rhs
    orig_write = bcb.apply_light_constraint_boundary

    def gst(*a, **k):
        out = orig_gst(*a, **k)
        cap.setdefault("gt", []).append({
            "E": {nm: np.array(out["explicit"][nm]) for nm in ("a", "b", "X")},
            "ll2": np.array(out["lambda_l2"]), "ll3": np.array(out["lambda_l3"])})
        return out

    def l2(grid_, st):
        out = orig_l2(grid_, st)
        if sys._getframe(1).f_code.co_name == "step":
            cap["l2_last"] = (grid_, st.copy(), {k: np.array(out[k]) for k in ("Aa", "K")})
        return out

    def ll2(grid_, st, **kw):
        out = orig_ll2(grid_, st, **kw)
        if sys._getframe(1).f_code.co_name == "step":
            cap["ll2_last"] = np.array(out)
        return out

    def write(grid_, st, rho, j):
        pre = st.copy()
        cap.setdefault("writes", []).append(dict(pre=pre, rho=np.array(rho[-W:]), j=np.array(j[-W:]),
                                                 l2_final=cap.get("l2_last"), ll2_final=cap.get("ll2_last")))
        orig_write(grid_, st, rho, j)
        cap["writes"][-1]["post"] = st.copy()
        cap["writes"][-1]["mom_post"] = float(mom(st, j)[-1])

    K._geometry_stage_terms = gst
    vacuum.primary_l2_rhs = l2
    vacuum.lambda_l2_rhs = ll2
    bcb.apply_light_constraint_boundary = write

    qa, qb, qX = dep0["a"], dep0["b"], dep0["X"]
    acc_piece = {"Aa": [mp.mpf(0)] * 3, "K": [mp.mpf(0)] * 3}
    dtm = mpf(dt)
    try:
        for step in range(round(HORIZON / dt)):
            cap.clear()
            g0 = state.geometry.copy()
            state = K.step(state, dt)
            gt0, gtp, gt1 = cap["gt"][0], cap["gt"][1], cap["gt"][2]
            w1, w2 = cap["writes"]
            # departure copy: predictor and corrector metric
            sl = slice(n - W, n)
            qa1 = [qa[k] + dtm * mpf(gt0["E"]["a"][sl][k]) for k in range(W)]
            qb1 = [qb[k] + dtm * mpf(gt0["E"]["b"][sl][k]) for k in range(W)]
            qX1 = [qX[k] + dtm * mpf(gt0["E"]["X"][sl][k]) for k in range(W)]
            h = dtm / 2
            qa2 = [qa[k] + h * (mpf(gt0["E"]["a"][sl][k]) + mpf(gt1["E"]["a"][sl][k])) for k in range(W)]
            qb2 = [qb[k] + h * (mpf(gt0["E"]["b"][sl][k]) + mpf(gt1["E"]["b"][sl][k])) for k in range(W)]
            qX2 = [qX[k] + h * (mpf(gt0["E"]["X"][sl][k]) + mpf(gt1["E"]["X"][sl][k])) for k in range(W)]

            # ---- T1: stored update equals hp rebuild from the operands used
            def t1(stored, base, inc):
                worst = 0.0
                for k in range(W - 3, W):
                    i = n - W + k
                    exact = mpf(base[i]) + inc(i)
                    tol = float(np.spacing(abs(stored[i]))) + 4 * float(np.spacing(abs(float(inc(i))) or 1e-300))
                    worst = max(worst, float(abs(mpf(stored[i]) - exact)) / tol)
                return worst

            p1, p2 = w1["pre"], w2["pre"]
            t1r = {}
            for nm in ("a", "b", "X"):
                t1r["pred_" + nm] = t1(getattr(p1, nm), getattr(g0, nm), lambda i, nm=nm: dtm * mpf(gt0["E"][nm][i]))
                t1r["corr_" + nm] = t1(getattr(p2, nm), getattr(g0, nm),
                                       lambda i, nm=nm: h * (mpf(gt0["E"][nm][i]) + mpf(gt1["E"][nm][i])))
            t1r["pred_Lambda"] = t1(p1.Lambda, g0.Lambda, lambda i: dtm * (mpf(gt0["ll2"][i]) / 2 + mpf(gtp["ll2"][i]) / 2
                                                                           + mpf(gt0["ll3"][i])))
            t1r["corr_Lambda"] = t1(p2.Lambda, g0.Lambda, lambda i: h * (mpf(gt0["ll2"][i]) + mpf(w2["ll2_final"][i])
                                                                         + mpf(gt0["ll3"][i]) + mpf(gt1["ll3"][i])))
            write_only = True
            for w_ in (w1, w2):
                for nm in ("a", "b", "X", "K", "Lambda", "alpha", "beta", "B"):
                    write_only &= bool(np.array_equal(getattr(w_["pre"], nm), getattr(w_["post"], nm)))
                write_only &= bool(np.array_equal(w_["pre"].Aa[:-1], w_["post"].Aa[:-1]))

            # ---- T3 / control: stage increments of C relative to the step base g0
            cb64 = C64(g0); cbd = Cdep(g0, qa, qb); cbc = Cctl(g0)
            for stage, w_, qa_s, qb_s, qX_s in (("pred", w1, qa1, qb1, qX1), ("corr", w2, qa2, qb2, qX2)):
                ca64 = C64(w_["pre"]); cad = Cdep(w_["pre"], qa_s, qb_s); cac = Cctl(w_["pre"])
                d64 = [mpf(ca64[i]) - mpf(cb64[i]) for i in range(3)]
                ddp = [cad[i] - cbd[i] for i in range(3)]
                dct = [cac[i] - cbc[i] for i in range(3)]
                aa_dep = hp_write(qa_s, qb_s, qX_s, w_["pre"], w_["rho"], w_["j"], r, dr)
                rec["stages"].append(dict(
                    step=step, stage=stage,
                    dC64=[float(x) for x in d64], dCdep=[float(x) for x in ddp],
                    diff_dep=[float(d64[i] - ddp[i]) for i in range(3)],
                    diff_ctl=[float(d64[i] - dct[i]) for i in range(3)],
                    C64=[float(x) for x in ca64], Cdep=[float(x) for x in cad],
                    Aa_written=float(w_["post"].Aa[-1]), Aa_dep=float(aa_dep),
                    Aa_diff=float(mpf(w_["post"].Aa[-1]) - aa_dep), mom_post=w_["mom_post"]))

            # ---- T2: mixed-stage piece in the corrector L2 (predictor Lambda on the final metric vs carried)
            grid_, fstate, l2eng = w2["l2_final"]
            lm_pred = [lam_metric(qa1, qb1, r, k, dr) for k in range(W)]
            lm_fin = [lam_metric(qa2, qb2, r, k, dr) for k in range(W)]
            carry = fstate.copy()
            for k in range(2, W):  # centred metric connection needs k>=2 in the window
                carry.Lambda[n - W + k] = float(mpf(fstate.Lambda[n - W + k]) + lm_fin[k] - lm_pred[k])
            l2c = orig_l2(grid_, carry)
            piece = {}
            for nm in ("Aa", "K"):
                eng = l2eng[nm]; car = np.array(l2c[nm])
                vals = [h * (mpf(eng[n - 3 + i]) - mpf(car[n - 3 + i])) for i in range(3)]
                acc_piece[nm] = [acc_piece[nm][i] + vals[i] for i in range(3)]
                piece[nm] = [float(v) for v in vals]
            rec["steps"].append(dict(step=step, t1_worst_over_tol=t1r, write_only_Aa_last=write_only, piece=piece))
            qa, qb, qX = qa2, qb2, qX2
            peak = max(peak, float(np.max(np.abs(mom_state(state)[-3:]))))
    finally:
        K._geometry_stage_terms = orig_gst
        vacuum.primary_l2_rhs = orig_l2
        vacuum.lambda_l2_rhs = orig_ll2
        bcb.apply_light_constraint_boundary = orig_write
    rec["G0_peak_e12"] = peak * 1e12
    rec["acc_piece"] = {nm: [float(x) for x in acc_piece[nm]] for nm in acc_piece}
    return rec


def rms(x):
    return float(np.sqrt(np.mean(np.square(x)))) if len(x) else float("nan")


def gates(runs):
    out = {}
    key = lambda r: (r["fixture"], r["n"], r["dt"])
    R = {key(r): r for r in runs}
    g0 = {f"{k[0]} {k[1]} {k[2]}": dict(measured=round(R[k]["G0_peak_e12"], 3), expected=v,
                                         exact=R[k]["G0_peak_e12"], match=round(R[k]["G0_peak_e12"], 3) == v)
          for k, v in G0_EXPECTED.items()}
    out["G0"] = dict(rows=g0, PASS=all(v["match"] for v in g0.values()))
    t1 = max(max(v for v in s["t1_worst_over_tol"].values()) for r in runs for s in r["steps"])
    wo = all(s["write_only_Aa_last"] for r in runs for s in r["steps"])
    out["T1"] = dict(worst_residual_over_tolerance=t1, write_touches_only_last_Aa=wo, PASS=bool(t1 <= 1.0 and wo))
    # T2
    t2 = {}
    for fx in ("W1", "PC2"):
        for n in (1200, 2400):
            a, b = R[(fx, n, 0.00075)], R[(fx, n, 0.000375)]
            FF = a["F"]
            def P(rr):
                return max(abs(rr["acc_piece"]["Aa"][0]), abs(rr["acc_piece"]["Aa"][1]),
                           *[abs(x) for x in rr["acc_piece"]["K"]])
            pa, pb = P(a), P(b)
            ratio = pa / pb if pb else float("inf")
            if pa < FF and pb < FF:
                verdict = "BELOW F at both dt: mixed-stage term does not matter at the edge"
            else:
                verdict = "PASS" if 3.5 <= ratio <= 4.5 else "FAIL"
            t2[f"{fx} {n}"] = dict(P_full=pa, P_half=pb, ratio=ratio, F=FF, verdict=verdict,
                                   acc_full=a["acc_piece"], acc_half=b["acc_piece"])
    out["T2"] = t2
    # T3 + setup + control + T4
    t3 = {"per_run": {}, "ratios": {}}
    worst3 = 0.0; worst4 = 0.0; mom_ok = True; setup_ok = True
    for r in runs:
        FF = r["F"]
        dd = [s["diff_dep"][2] for s in r["stages"]]
        dc = [s["diff_ctl"][2] for s in r["stages"]]
        w = max(abs(x) for x in dd) / FF
        worst3 = max(worst3, w)
        su = abs(r["setup"]["C64_minus_Cdep"][2]) / FF
        setup_ok &= su <= 1.0
        w4 = max(abs(s["Aa_diff"]) for s in r["stages"]) / FF
        worst4 = max(worst4, w4)
        lim = 3e-12 if r["n"] == 1200 else 8e-12
        mp_ = max(abs(s["mom_post"]) for s in r["stages"])
        mom_ok &= mp_ <= lim
        t3["per_run"][f"{r['fixture']} {r['n']} {r['dt']}"] = dict(
            max_diff_dep_over_F=w, rms_diff_dep=rms(dd), rms_diff_ctl=rms(dc),
            ctl_over_dep_rms=rms(dc) / rms(dd) if rms(dd) else float("nan"),
            setup_C64_minus_Cdep_over_F=su, setup_C64_minus_Cctl=r["setup"]["C64_minus_Cctl"][2],
            T4_max_Aa_diff_over_F=w4, T4_max_mom_post=mp_, T4_mom_limit=lim,
            cells_n3_n2_max_diff_dep=[max(abs(s["diff_dep"][i]) for s in r["stages"]) for i in (0, 1)],
            sign_lag1_corr=float(np.corrcoef(dd[:-1], dd[1:])[0, 1]) if len(dd) > 2 else None)
    rat_ok = True
    for fx in ("W1", "PC2"):
        for dt in (0.00075, 0.000375):
            a, b = R[(fx, 1200, dt)], R[(fx, 2400, dt)]
            ra = rms([s["diff_dep"][2] for s in b["stages"]]) / rms([s["diff_dep"][2] for s in a["stages"]])
            t3["ratios"][f"{fx} {dt}"] = ra
            rat_ok &= 1.3 <= ra <= 3.0
    out["setup"] = dict(PASS=bool(setup_ok))
    out["T3"] = dict(worst_stage_diff_over_F=worst3, refinement_ratios=t3["ratios"],
                     PASS=bool(worst3 <= 2.0 and rat_ok), detail=t3["per_run"],
                     control="report-only; predicted 1/4 to 1/2 of the departure difference")
    out["T4"] = dict(worst_Aa_diff_over_F=worst4, limit_over_F=1.7, momentum_within_ceiling=bool(mom_ok),
                     PASS=bool(worst4 <= 1.7 and mom_ok))
    return out


def main(outdir):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    runs = []
    for fx in ("W1", "PC2"):
        for n in (1200, 2400):
            for dt in (0.00075, 0.000375):
                rr = run(fx, n, dt)
                runs.append(rr)
                print(fx, n, dt, "G0 peak %.4f e-12" % rr["G0_peak_e12"], flush=True)
    res = dict(gates=gates(runs), runs=runs)
    (outdir / "report.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps(res["gates"], indent=1, default=float)[:6000])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "results/research/b7_w1_outer_cell_trace")
