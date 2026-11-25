"""60-digit readback of the stage-1 outputs (python3 -m research.control_tests.mp_readback OUT). Research only.

Runs under plain python3 with mpmath and numpy (no engine import). Reads only the files that test1.py wrote in OUT
and writes OUT/readback.json. The gate values are the predeclared ones (A6 B, A9, A12, A13 section 1, A15); nothing
here is fitted to or tuned on the data it reads.

  inputs:  A13 section 1 (a, b, X, Aa, K, Lambda RHS: |stored - 60 digits| <= gamma_m Sum|s_k|), B1, B3 (i)-(iii),
           B4, B5, B6 (G1b.1 identity at t = 0), and the independence check (60 digits rounded vs the engine
           within 100 eps x the largest term) on every input slice.
  arms:    G1b.1 identity at each save of the 1b arm (dr .05 only; the order part is withdrawn, A15 section 1);
           at each recorded step, the GK.1 prediction at the primary-predictor and accepted stages (A6 A1, A15
           section 2), the K/Aa RHS evaluation rounding passed into M (A9 section 7, A12), and the K' term in
           60 digits beside its binary64 value.
Edge cells (last two) are readings only (A6 B0, A9 section 1) and are left out of every pass/fail.
"""
import json
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

from research.control_tests import mp60 as M

EPS = 2.0 ** -52
U = 2.0 ** -53
EDGE = 2
FIELDS = ("a", "b", "X", "Aa", "K", "Lambda", "alpha")


def fixed_points(v, dr):
    """Same as ctl.fixed_points: cell 0, then r^2-weighted averages onto r .075, .125, .175, .225."""
    v = np.asarray(v, dtype=float)
    f = round(0.05 / dr)
    r = (np.arange(len(v)) + 0.5) * dr
    out = [float(v[0])]
    for j in (1, 2, 3, 4):
        sl = slice(j * f, j * f + f)
        out.append(float(np.sum(r[sl] ** 2 * v[sl]) / np.sum(r[sl] ** 2)))
    return out


def state(z, prefix, alpha_key=None):
    s = {k: M.arr(z[prefix + k]) for k in FIELDS}
    if alpha_key is not None:
        s["alpha"] = M.arr(z[alpha_key])
    return s


def matter(z, prefix):
    m = {k: M.arr(z[f"{prefix}{k}"]) for k in ("rho", "pr", "pt")}
    m["j"] = [x + y for x, y in zip(M.arr(z[f"{prefix}j_scalar"]), M.arr(z[f"{prefix}j_fluid"]))]
    return m


def rhs_all(g, s, m):
    ex = M.explicit_rhs(g, s)
    l2, l3 = M.rhs_primary(g, s, m)
    lam = [x + y for x, y in zip(M.lambda_l2(g, s), M.lambda_l3(g, s, m))]
    return {"a": ex["a"], "b": ex["b"], "X": ex["X"],
            "Aa": [x + y for x, y in zip(l2["Aa"], l3["Aa"])], "K": [x + y for x, y in zip(l2["K"], l3["K"])],
            "Lambda": lam}, l2, l3


def compare(stored, xs, n_int):
    """A13 section 1: worst |stored - 60 digits| / (gamma_m Sum|s_k|) over interior cells, and the independence
    check |stored - round(60 digits)| <= 100 eps max_k a_k."""
    st = np.asarray(stored, dtype=float)
    v, b = M.vals(xs), M.bound(xs)
    diff = np.array([abs(mp.mpf(float(st[i])) - xs[i].v) for i in range(len(xs))], dtype=float)
    ratio = diff[:n_int] / np.maximum(b[:n_int], 1e-300)
    amax = max(float(x.a) for x in xs[:n_int])
    indep = float(np.max(np.abs(st[:n_int] - v[:n_int])))
    return {"worst_ratio": float(np.max(ratio)), "cell": int(np.argmax(ratio)), "pass": bool(np.max(ratio) <= 1.0),
            "edge_ratio": [float(x) for x in diff[n_int:] / np.maximum(b[n_int:], 1e-300)],
            "independence": {"max_diff": indep, "limit": 100 * EPS * amax, "pass": indep <= 100 * EPS * amax}}


def b_g(g, s, dr):
    """B6: B_G = gamma_8 (Sum|w_D4| + Sum|w_D2|) max|alpha| / dr^2, per cell (c dr from the solve's own c)."""
    a, b, X = s["a"], s["b"], s["X"]
    ap, bp, Xp = g.d(a, 1, 1), g.d(b, 1, 1), g.d(X, 1, 1)
    amax = max(abs(x.v) for x in s["alpha"])
    out = []
    for i in range(g.n):
        c = -0.5 * ap[i].v / a[i].v + bp[i].v / b[i].v - Xp[i].v / X[i].v + 2 / g.r[i].v
        cdr = abs(c) * g.dr.v
        w4 = mp.mpf(64) / 12 + cdr * mp.mpf(18) / 12   # 4th-order d2 (1,16,30,16,1)/12 and c d1 (1,8,8,1)/12
        w2 = 4 + cdr                                   # 3-point d2 (1,2,1) and c d1 (1,1)/2
        out.append(float(M.gamma(8) * (w4 + w2) * amax / g.dr.v ** 2))
    return np.array(out)


def identity(z, prefix, g, s, dr):
    """G1b.1: |K RHS_nonuniform - (Delta_2 - Delta_4) alpha_RHS| <= 4 B_G per fixed point; counts if signal >= 10 B_G."""
    d24 = M.delta2_minus_delta4(g, s)
    sig = np.array([np.nan if x is None else float(x.v) for x in d24])
    nonu = np.asarray(z[prefix + "R_K"], dtype=float) - float(z[prefix + "kdot"])
    bg = b_g(g, s, dr)
    n_int = g.n - EDGE
    fp = lambda v: fixed_points(np.where(np.isnan(v[:n_int]), 0.0, v[:n_int]), dr)
    res, sg, bgp = fp(nonu - sig), fp(sig), fp(bg)
    pts = [{"residual": r_, "signal": s_, "B_G": b_, "pass": abs(r_) <= 4 * b_, "counts": abs(s_) >= 10 * b_}
           for r_, s_, b_ in zip(res, sg, bgp)]
    return {"points": pts, "pass": all(p["pass"] for p in pts), "n_counting": sum(p["counts"] for p in pts)}


def h_largest_term(g, s, rho):
    """A2 P4's 'largest term in H', read as the largest top-level summand of the vendor H (the nine summands of the
    R bracket times -X^2/a, the Aa and K terms, and 16 pi rho), per cell. Reading only; see t2_checks."""
    a, b, X, Lam, r = s["a"], s["b"], s["X"], s["Lambda"], g.r
    ap, bp, app, bpp = g.d(a, 1, 1), g.d(b, 1, 1), g.d(a, 1, 2), g.d(b, 1, 2)
    Xp, Xpp, Lp = g.d(X, 1, 1), g.d(X, 1, 2), g.d(Lam, -1, 1)
    out = []
    for i in range(g.n):
        A, B, x, L, rr = a[i].v, b[i].v, X[i].v, Lam[i].v, r[i].v
        cp = -Xp[i].v / (2 * x)
        cpp = -(Xpp[i].v / x - (Xp[i].v / x) ** 2) / 2
        ml = (1 - A / B) / rr ** 2
        t = [app[i].v / (2 * A), bpp[i].v / B, -A * Lp[i].v, -(ap[i].v / A) ** 2, (bp[i].v / B) ** 2 / 2,
             2 * (3 - A / B) * bp[i].v / (rr * B), 4 * ml, 8 * (cpp + cp ** 2),
             -8 * cp * (ap[i].v / (2 * A) - bp[i].v / B - 2 / rr)]
        inv = x * x / A
        Aa, K = s["Aa"][i].v, s["K"][i].v
        out.append(max([abs(inv * y) for y in t] + [abs(Aa * Aa * mp.mpf("1.5")), abs(2 * K * K / 3),
                                                     abs(16 * mp.pi * rho[i].v)]))
    return out


def ps_rhs(g, s, S, P):
    """scalar_rhs_arrays' P_S RHS, engine operation order, beta = 0, V' = S (lambda = 0, BETA_S = 0, D = phi = 0)."""
    a, b, X, al, K, r = s["a"], s["b"], s["X"], s["alpha"], s["K"], g.r
    ap, bp, Xp, alp = (g.d(v, 1, 1) for v in (a, b, X, al))
    fp, fpp, pp = g.d(S, 1, 1), g.d(S, 1, 2), g.d(P, 1, 1)
    out = []
    for i in range(g.n):
        inv = (X[i] ** 2) / a[i]
        c = ((0.5 * ap[i]) / a[i] - bp[i] / b[i]) + Xp[i] / X[i] - 2.0 / r[i]
        lap = inv * (fpp[i] - fp[i] * c)
        out.append(M.E(0) * pp[i] + al[i] * lap + (al[i] * K[i]) * P[i] + (inv * alp[i]) * fp[i] - al[i] * S[i])
    return out


def t2_checks(z, g, s, s_in, m, n_int, meta):
    """A2 section 2, Test 2 (with A13 section 1 for the a, b, X, Aa, K, Lambda RHS, done in input_checks)."""
    out = {}
    # Newton residual (A2 P4 as fixed by A22 item 2): the gate is the B5 form, |H| <= 10 eps x Sum|s_k| per non-edge
    # cell, summed down to the stencil weights. The literal reading (largest top-level summand) is supplementary.
    H64 = np.asarray(z["H"], dtype=float)
    Hmp = [h - 16.0 * M.PI * rho for h, rho in zip(M.hamiltonian_vac(g, s_in), m["rho"])]
    big = h_largest_term(g, s_in, m["rho"])
    lit = np.array([abs(H64[i]) / float(10 * EPS * big[i]) for i in range(n_int)])
    b5 = np.array([abs(H64[i]) / float(10 * EPS * Hmp[i].a) for i in range(n_int)])
    out["newton_residual"] = {
        "literal_worst_ratio": float(np.max(lit)), "literal_cell": int(np.argmax(lit)),
        "literal_pass": bool(np.max(lit) <= 1), "literal_cells_failing": int(np.sum(lit > 1)),
        "B5_form_worst_ratio": float(np.max(b5)), "B5_form_pass": bool(np.max(b5) <= 1),
        "pass": bool(np.max(b5) <= 1), "gate": "B5 form (A22 item 2); literal_* supplementary",
        "max_abs_H_non_edge": float(np.max(np.abs(H64[:n_int]))),
        "max_abs_H60_at_stored_X": float(max(abs(h.v) for h in Hmp[:n_int])),
        "newton_log": meta.get("newton"),
        "note": "A22 item 2: B5 form is the gate"}
    # M(0) = 0 bitwise (P_S = 0, so j = 0)
    Mz = np.asarray(z["M"], dtype=float)
    out["M0_bitwise"] = {"nonzero_cells": int(np.count_nonzero(Mz)), "pass": bool(np.count_nonzero(Mz) == 0)}
    # S RHS = alpha P_S (= 0), bitwise
    al = np.asarray(z["alpha_rhs"], dtype=float)
    want = al * np.asarray(z["s_PS"], dtype=float)
    out["S_rhs_pi_term"] = {"equal_bitwise": bool(np.array_equal(np.asarray(z["R_s_S"], dtype=float), want)),
                            "max_abs": float(np.max(np.abs(z["R_s_S"])))}
    out["S_rhs_pi_term"]["pass"] = out["S_rhs_pi_term"]["equal_bitwise"]
    # P_S RHS: (i) the stored binary64 against 60 digits on the same inputs (A13 section 1 form);
    # (ii) the stored binary64 against 60 digits on the exact Gaussian samples (analytic (Laplacian - 1) plus the
    #      stencil truncation, predicted in 60 digits), within the evaluation bound plus the input rounding of S.
    S, P = M.arr(z["s_S"]), M.arr(z["s_PS"])
    out["PS_rhs_vs_60_digits"] = compare(z["R_s_PS"], ps_rhs(g, s, S, P), n_int)
    amp, sig = mp.mpf(meta["amp"]), mp.mpf(meta.get("sigma", 1.0))
    Sx = [M.E(amp * mp.exp(-(x.v / sig) ** 2)) for x in g.r]
    px = ps_rhs(g, s, Sx, P)
    st = np.asarray(z["R_s_PS"], dtype=float)
    dif = np.array([float(abs(mp.mpf(float(st[i])) - px[i].v)) for i in range(n_int)])
    bnd = np.array([float((M.gamma(px[i].n) + M.U) * px[i].a) for i in range(n_int)])
    # analytic: alpha (X^2/a)(S'' - S' c) + (X^2/a) alpha' S' - alpha S, analytic S', S'', 60-digit geometry stencils
    a_, b_, X_, al_ = s["a"], s["b"], s["X"], s["alpha"]
    ap, bp, Xp, alp = (g.d(v, 1, 1) for v in (a_, b_, X_, al_))
    an = []
    for i in range(g.n):
        rr, Sv = g.r[i].v, Sx[i].v
        S1, S2 = -2 * rr * Sv / sig ** 2, (4 * rr ** 2 / sig ** 4 - 2 / sig ** 2) * Sv
        inv = X_[i].v ** 2 / a_[i].v
        c = ap[i].v / (2 * a_[i].v) - bp[i].v / b_[i].v + Xp[i].v / X_[i].v - 2 / rr
        an.append(al_[i].v * inv * (S2 - S1 * c) + inv * alp[i].v * S1 - al_[i].v * Sv)
    trunc = np.array([float(px[i].v - an[i]) for i in range(g.n)])
    meas = np.array([float(mp.mpf(float(st[i])) - an[i]) for i in range(g.n)])
    out["PS_rhs_vs_analytic"] = {
        "worst_ratio_to_bound": float(np.max(dif / np.maximum(bnd, 1e-300))), "cell": int(np.argmax(dif / np.maximum(bnd, 1e-300))),
        "pass": bool(np.all(dif <= bnd)),
        "stored_minus_analytic_fixed_points": fixed_points(meas[:n_int], float(g.dr.v)),
        "predicted_truncation_fixed_points": fixed_points(trunc[:n_int], float(g.dr.v)),
        "max_abs_predicted_truncation": float(np.max(np.abs(trunc[:n_int]))),
        "note": "A22 item 3: gate is |stored - 60 digits on exact samples| <= (gamma_n + u) Sum|s_k|, on the same stored "
                "X, alpha (solved), a, b as the engine call; the *_analytic and truncation entries are readings"}
    return out


def input_checks(z, dr, kind, meta=None):
    n = len(z["g_a"])
    g = M.Grid(n, dr)
    n_int = n - EDGE
    s = state(z, "g_", alpha_key="alpha_rhs")       # the RHS slice: the solved alpha
    s_in = state(z, "g_")                             # the solve's input slice: the stored alpha
    m = matter(z, "m_")
    rhs, _, _ = rhs_all(g, s, m)
    out = {"rhs": {k: compare(z[f"R_{k}"], rhs[k], n_int) for k in rhs}}
    # B1 / B5: H against the dust source, |H_stored - H_60| <= gamma_m Sum|s_k|
    H = [h - 16.0 * M.PI * rho for h, rho in zip(M.hamiltonian_vac(g, s_in), m["rho"])]
    out["B5_H"] = compare(z["H"], H, n_int)
    out["B1"] = {"H": fixed_points(z["H"], dr), "16pi_rho_dust": fixed_points(16 * math.pi * z["m_rho"], dr)}
    # B4: |M(0)| <= (2/3) gamma_4 1.5 |K0| / dr
    K0 = abs(float(z["g_K"][0]))
    lim = float((mp.mpf(2) / 3) * M.gamma(4) * mp.mpf("1.5") * K0 / g.dr.v)
    out["B4"] = {"M0": float(z["M"][0]), "limit": lim, "pass": abs(float(z["M"][0])) <= lim}
    # B3 (ii): Kdot against the 60-digit projection on the input slice
    kv, kb = M.kdot_projection(g, s_in, m)
    kd = float(z["kdot"])
    out["B3_ii_kdot"] = {"stored": kd, "mp": float(kv), "bound": float(kb), "pass": abs(mp.mpf(kd) - kv) <= kb}
    if kind == "1a":
        al = np.asarray(z["alpha_rhs"], dtype=float)
        lim1 = 6.5 * n * n * U
        out["B3_i"] = {"max_abs_alpha_minus_1": float(np.max(np.abs(al - 1))), "limit": lim1,
                       "pass": float(np.max(np.abs(al - 1))) <= lim1}
        nz = int(np.count_nonzero(np.asarray(z["R_Aa"], dtype=float)[:n_int]))
        # reading only: the RHS slice carries the solved alpha (|alpha - 1| up to B3 (i)), so its Aa RHS need not be 0.
        # The gated B3 (iii) is the alpha = 1 repeat (A9 section 2), from b3iii.py, merged in main().
        out["solved_alpha_Aa_rhs_nonzero_interior_cells"] = {"count": nz}
    elif kind == "1b":
        out["B6_G1b1"] = identity(z, "", g, s, dr)
    elif kind == "T2":
        out.update(t2_checks(z, g, s, s_in, m, n_int, meta or {}))
    return out


def rms(v):
    v = np.asarray(v, dtype=float)
    return float(np.sqrt(np.mean(v * v)))


def recorded_step(z, dr, dt):
    """GK.1 prediction, rounding shares and the 60-digit K' term at one recorded step."""
    n = len(z["start_a"])
    g = M.Grid(n, dr)
    n_int = n - EDGE
    S = {name: (state(z, f"{name}_"), matter(z, f"{name}_m_"))
         for name in ("start", "stage0", "explicit_predictor", "primary", "final_l2_call", "accepted")}
    calls = {}
    for name in ("stage0", "explicit_predictor", "primary"):
        s, m = S[name]
        l2, l3 = M.rhs_primary(g, s, m)
        calls[name] = {"l2": l2, "l3": l3, "D": M.delta2_minus_delta4(g, s), "kdot": float(z[f"call_{name}_kdot"])}
    s, m = S["final_l2_call"]
    calls["final_l2"] = {"l2": M.rhs_primary(g, s, m)[0]}
    v = lambda xs: np.array([np.nan if x is None else float(x.v) for x in xs])
    h = mp.mpf(dt)
    half = lambda xs, ys: [None if (x is None or y is None) else (x.v + y.v) * h / 2 for x, y in zip(xs, ys)]
    c0, ce, c1, cf = calls["stage0"], calls["explicit_predictor"], calls["primary"], calls["final_l2"]
    K0 = np.asarray(z["start_K"], dtype=float)
    # primary predictor: K1 - K0 - dt(kdot0 + kdot_e)/2 = dt/2 (D0 + De) + dt/2 (l3_0 - l3_e)
    pred_p = [None if d is None else d + (a.v - b.v) * h / 2
              for d, a, b in zip(half(c0["D"], ce["D"]), c0["l3"]["K"], ce["l3"]["K"])]
    meas_p = np.asarray(z["primary_K"], dtype=float) - K0 - dt * (c0["kdot"] + ce["kdot"]) / 2
    # accepted: Knew - K0 - dt(kdot0 + kdot1)/2 = dt/2 (D0 + D1) + dt/2 (l2_final - l2_1)
    pred_a = [None if d is None else d + (a.v - b.v) * h / 2
              for d, a, b in zip(half(c0["D"], c1["D"]), cf["l2"]["K"], c1["l2"]["K"])]
    meas_a = np.asarray(z["accepted_K"], dtype=float) - K0 - dt * (c0["kdot"] + c1["kdot"]) / 2
    gk = {}
    for nm, pred, meas in (("primary", pred_p, meas_p), ("accepted", pred_a, meas_a)):
        pv = np.array([np.nan if x is None else float(x) for x in pred])
        gk[nm] = {"pred": fixed_points(np.nan_to_num(pv[:n_int]), dr), "measured": fixed_points(meas[:n_int], dr),
                  "residual": fixed_points(np.nan_to_num(meas[:n_int] - pv[:n_int]), dr)}
    # rounding of each K/Aa RHS call (stored - 60 digits), through the accepted-step weights dt/2, into M
    stored = lambda c, blk, k: np.asarray(z[f"call_{c}_{blk}_{k}"], dtype=float)
    used = (("stage0", "primary_l2", c0["l2"]), ("final_l2", "primary_l2", cf["l2"]),
            ("stage0", "primary_l3", c0["l3"]), ("primary", "primary_l3", c1["l3"]))
    dK = np.zeros(n); dA = np.zeros(n)
    per_call = {}
    for c, blk, mpblk in used:
        eK = stored(c, blk, "K") - v(mpblk["K"]); eA = stored(c, blk, "Aa") - v(mpblk["Aa"])
        dK += 0.5 * dt * eK; dA += 0.5 * dt * eA
        per_call[f"{c}_{blk}"] = {"K_rms_interior": rms(eK[:n_int]), "Aa_rms_interior": rms(eA[:n_int])}
    sA = S["start"][0]
    gt = M.geometry_terms(g, sA)
    dAE, dKE = M.arr(dA), M.arr(dK)
    mt = M.momentum_terms(g, dict(sA, Aa=dAE, K=dKE), gt)
    Mr = v(mt["M_vac"])
    dM = v(M.momentum_terms(g, S["accepted"][0])["M_vac"]) - v(M.momentum_terms(g, sA, gt)["M_vac"])
    band = lambda x: fixed_points(np.nan_to_num(x[:n_int]), dr)[1:]
    share = {"M_from_rounding_band": band(Mr), "dM_step_band": band(dM),
             "band_rms_ratio": rms(band(Mr)) / max(rms(band(dM)), 1e-300),
             "native_sign_change_fraction": float(np.mean(np.sign(Mr[1:n_int]) != np.sign(Mr[:n_int - 1]))),
             "per_call": per_call}
    kp = {}
    for nm in ("start", "stage0", "explicit_predictor", "primary", "final_l2_call", "accepted"):
        s = S[nm][0]
        kp[nm] = {"mp60": fixed_points(np.nan_to_num(v(M.k_prime(g, s))[:n_int]), dr)}
    return {"GK1": gk, "rounding_share": share, "K_prime": kp}


def main(out):
    out = Path(out)
    res = {"inputs": {}, "arms": {}}
    inp = out / "inputs"
    if (inp / "inputs.json").exists():
        meta = json.loads((inp / "inputs.json").read_text())
        for name, p in meta.items():
            with np.load(inp / f"{name}.npz") as z:
                res["inputs"][name] = input_checks(z, p["dr"], name[:2], p)
        b3 = inp / "b3iii.json"                       # A9 section 2 alpha = 1 repeat (research.control_tests.b3iii)
        if b3.exists():
            for name, v in json.loads(b3.read_text()).items():
                res["inputs"].setdefault(name, {})["B3_iii_alpha1_repeat"] = v
    for js in sorted(out.glob("*_rstage.json")):
        arm = js.name[: -len("_rstage.json")]
        rec = json.loads((out / f"{arm}.json").read_text())
        dr, dt = rec["dr"], rec["dt"]
        a = {"saves": [], "steps": []}
        if arm.startswith("1b") and abs(dr - 0.05) < 1e-12:
            with np.load(out / f"{arm}.npz") as z:
                k = 0
                while f"t{k}_g_a" in z:
                    n = len(z[f"t{k}_g_a"])
                    g = M.Grid(n, dr)
                    s = state(z, f"t{k}_g_", alpha_key=f"t{k}_alpha_rhs")
                    a["saves"].append(identity(z, f"t{k}_", g, s, dr))
                    k += 1
        for i, _ in enumerate(json.loads(js.read_text())):
            with np.load(out / f"{arm}_rstage{i}.npz") as z:
                a["steps"].append(recorded_step(z, dr, dt))
        res["arms"][arm] = a
    (out / "readback.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1])
