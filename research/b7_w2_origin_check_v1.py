"""RESEARCH ONLY. B7 W2 zero-step origin check (Aster's review, Drive 1nEdruskrM23fooJmQTu2jo75hevhiNDUtMFVRWj2sm8).

Pre-run design, gates and pass/fail criteria: coordination log entry 2026-11-23 ~00:45 UTC (fixed before any run).
Go: Jeremiah, main chat, 2026-11-23 00:31:38 UTC, "Go on the origin check". No time evolution. Engine, vendor and the
W2 clone (research/b7_ring_constraint_clone_v1.py) are used read-only; nothing in engine/ or vendor/ changes.

Usage: python research/b7_w2_origin_check_v1.py RING.npz OUTDIR   (RING.npz = source C t300 member, sha f3675f71...)

Steps: 1 reproduce/localize; 2a quadrature unit test; 2b manufactured smooth slice with an independent ODE integral;
2c quadrature swap on the W2 samples; 3a even (spline) reconstruction of the archived arrays; 3b u = 1.0144 UNTESTED (W1).
"""
import json
import math
import os
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402
from scipy.interpolate import CubicSpline  # noqa: E402

from engine import scalar_system as _ss  # noqa: E402
from engine import v55_pirk_adapter as adapter  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402
from engine.scalar_system import cosmos_potential, sd_potential  # noqa: E402
from engine.stress_energy import assemble_total_stress_energy  # noqa: E402
from engine.v55_initial import CorrectedInitialData, _cumulative_origin_integral, _d1  # noqa: E402
from engine.v55_matter import RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT, initialize_from_archive  # noqa: E402
from research import b7_ring_constraint_clone_run_v1 as w2  # noqa: E402
from research import b7_ring_constraint_clone_v1 as clone  # noqa: E402
from research import b7_ring_source_v1 as src  # noqa: E402

# ---- fixed before the run (log entry 2026-11-23 ~00:45 UTC) ----
R_MAX = 60.0
DRS = (0.2, 0.1, 0.05)
G4C_FACTOR = 10.0
REPRO_RTOL = 1e-12
QUAD_EXACT_RTOL = 1e-13
ODE_RTOL, ODE_RTOL_CHECK, ODE_AGREE = 1e-13, 1e-11, 1e-10
STENCIL_KEEP, QUAD_DROP = 0.5, 0.1          # 2b and 2c readings
SAMPLING_DROP, NO_EFFECT = 10.0, 2.0        # step 3 readings
ANNULUS = (2.0, 50.0)
FIXED_R = 0.3
MAN_W2 = 25.0                               # manufactured profile width^2
MAN_FAR_R = 40.0
FIELDS = ("S", "PS", "D", "PD")
TOY_KEYS = {"S": "t300_s", "PS": "t300_ps", "D": "t300_d", "PD": "t300_pd"}
LAM = _ss.SD_QUARTIC_LAMBDA


# ---------------------------------------------------------------- reconstruction variants
def spline_even_integral(r, f):
    """Independent reference: exact integral of a cubic spline through the even-mirrored samples."""
    x = np.r_[-r[::-1], r]
    cs = CubicSpline(x, np.r_[f[::-1], f])
    P = cs.antiderivative()
    return P(r) - P(0.0)


def even6_integral(r, f, dr):
    """Junction test (research only): one uniform rule at every cell. Each interval [r_{i-1}, r_i], and [0, r_0] first,
    is integrated exactly over the 6-point Lagrange interpolant through the nearest cell centres, with even ghosts
    f(-r_k) = f(r_k) (cell-centred: r_{-k-1} = -r_k) at the origin and a one-sided window only at the outer edge."""
    n = len(r)
    f = np.asarray(f, dtype=float)
    pieces = np.empty(n)
    for i in range(n):
        lo = min(max(i - 3, -3), n - 6)
        idx = np.arange(lo, lo + 6)
        x = idx + 0.5                                   # positions in units of dr
        xa, xb = (0.0, 0.5) if i == 0 else (i - 0.5, i + 0.5)
        c = 0.5 * (xa + xb)
        V = np.vander(x - c, 6, increasing=True)
        p = np.arange(6)
        m = ((xb - c) ** (p + 1) - (xa - c) ** (p + 1)) / (p + 1)
        w = np.linalg.solve(V.T, m)
        vals = np.where(idx >= 0, f[np.abs(idx)], f[np.clip(-idx - 1, 0, n - 1)])
        pieces[i] = dr * float(w @ vals)
    return np.cumsum(pieces)


def _rho(S, PS, D, PD, Sp, Dp, phi, Pi, B, rho_r):
    rhoS = 0.5 * (PS ** 2 + B * Sp ** 2)
    rhoD = 0.5 * (PD ** 2 + B * Dp ** 2) + sd_potential(S, D, phi=phi)
    rhoPhi = (0.5 * Pi ** 2 + cosmos_potential(phi)) / (8.0 * math.pi)
    rhoMatter = (RHO_DM_PRESENT + RHO_B_PRESENT + rho_r) / (8.0 * math.pi)
    return rhoS + rhoD + rhoPhi + rhoMatter


def first_q4(r, f):
    """The engine's first piece, same arithmetic as _cumulative_origin_integral."""
    A = np.vstack([np.ones(3), r[:3] ** 2, r[:3] ** 4]).T
    c = np.linalg.solve(A, f[:3])
    return c[0] * r[0] + c[1] * r[0] ** 3 / 3.0 + c[2] * r[0] ** 5 / 5.0


def first_s6(r, f, dr):
    """First-piece test arm FP6: even sextic through cells 0-3, fitted in x = r/dr, integrated over [0, r_0]."""
    x = r[:4] / dr
    A = np.vstack([np.ones(4), x ** 2, x ** 4, x ** 6]).T
    c = np.linalg.solve(A, f[:4])
    x0 = x[0]
    return dr * (c[0] * x0 + c[1] * x0 ** 3 / 3.0 + c[2] * x0 ** 5 / 5.0 + c[3] * x0 ** 7 / 7.0)


def reconstruct(grid, vac_factory, S, PS, D, PD, phi, Pi, *, quadrature="engine", src_derivs=None, B_override=None,
                n_iter=6, first_override=None):
    """The W2 clone's arithmetic (Hubble K_t, PS terms on, uniform dust), with exactly one of three knobs changed:
    quadrature ('engine' | 'spline'), src_derivs (analytic (Sp, Dp) used in the Hamiltonian source and its j only),
    B_override (B supplied from outside). K, j for the geometry, conversion and evaluator are always the clone's."""
    assert _ss.coupling_active()
    r = np.asarray(grid.centers)
    dr = float(grid.dr)
    rho_r = RHO_R_PRESENT
    H0 = clone.background_H0(True, True)
    Kt = np.full(len(r), -H0)
    Sp, Dp, Ktp = _d1(grid, S, 1), _d1(grid, D, 1), _d1(grid, Kt, 1)
    j = -(PD * Dp) - PS * Sp
    sSp, sDp = (Sp, Dp) if src_derivs is None else src_derivs
    sj = -(PD * sDp) - PS * sSp
    quad = {"engine": lambda f: _cumulative_origin_integral(grid, r, dr, f),
            "spline": lambda f: spline_even_integral(r, f),
            "even6": lambda f: even6_integral(r, f, dr),
            "fp6": lambda f: _cumulative_origin_integral(grid, r, dr, f) - first_q4(r, f) + first_s6(r, f, dr),
            "fpo": lambda f: _cumulative_origin_integral(grid, r, dr, f) - first_q4(r, f) + first_override}[quadrature]
    B = np.ones_like(r)
    aux = {}
    if B_override is None:
        for _ in range(n_iter):
            rho = _rho(S, PS, D, PD, sSp, sDp, phi, Pi, B, rho_r)
            Kr = Kt + r * Ktp + 4.0 * math.pi * r * sj
            srcv = 1.0 - 8.0 * math.pi * r * r * rho + 2.0 * r * r * Kr * Kt + r * r * Kt * Kt
            Bn = quad(srcv) / r
            aux["last_dB_absmax"] = float(np.max(np.abs(Bn - B)))
            B = Bn
        aux["src_frozen"] = srcv
    else:
        B = np.asarray(B_override, dtype=float)
    A = 1.0 / np.sqrt(B)
    Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
    K = Kr + 2.0 * Kt
    a, b, X = A ** (4.0 / 3.0), A ** (-2.0 / 3.0), A ** (-1.0 / 3.0)
    Lam = _d1(grid, a, 1) / (2.0 * a * a) - _d1(grid, b, 1) / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)
    st = vac_factory(grid)
    st.a[:], st.b[:], st.X[:] = a, b, X
    st.alpha[:], st.beta[:] = 1.0, 0.0
    st.Aa[:], st.K[:], st.Lambda[:], st.B[:] = Kr - K / 3.0, K, Lam, 0.75 * Lam
    scal = _ss.ScalarFields(S, PS, D, PD, phi, Pi)
    mat = initialize_from_archive(grid, st, include_radiation=True, charged_dust_density=None)
    aux.update(B=B, Kr=Kr, Kt=Kt)
    return CorrectedInitialData(geometry=st, scalars=scal, matter=mat, H0=H0, include_radiation=True), aux


def exact_B(grid, fn, rtol):
    """Independent integral: (rB)' = src(r, B) with analytic sources, y = rB linear in y. fn(r) -> fields and derivs."""
    H0 = clone.background_H0(True, True)
    Kt = -H0

    def rhs(t, y):
        f = fn(np.array([t]))
        r0 = _rho(f["S"], f["PS"], f["D"], f["PD"], f["Sp"], f["Dp"], f["phi"], f["Pi"], 0.0, RHO_R_PRESENT)[0]
        r1 = 0.5 * (f["Sp"][0] ** 2 + f["Dp"][0] ** 2)
        jj = -(f["PD"][0] * f["Dp"][0]) - f["PS"][0] * f["Sp"][0]
        Kr = Kt + 4.0 * math.pi * t * jj
        return [1.0 - 8.0 * math.pi * t * t * r0 - 8.0 * math.pi * t * r1 * y[0] + 2.0 * t * t * Kr * Kt + t * t * Kt * Kt]

    r = np.asarray(grid.centers)
    sol = solve_ivp(rhs, (0.0, float(r[-1])), [0.0], method="DOP853", rtol=rtol, atol=1e-15, t_eval=r)
    return sol.y[0] / r


# ---------------------------------------------------------------- evaluator
def evaluate(kernel, grid, init):
    """Per-cell signed H and M exactly as kernel.diagnostics forms them; before and after the lapse solve."""
    _, vac, _ = adapter.vendor_modules()
    out = {}
    for tag in ("pre_lapse", "post_lapse"):
        if tag == "post_lapse":
            st = clone.research_state(kernel, grid, init)
            geom, scal, mat = st.geometry, st.scalars, st.matter
        else:
            geom, scal, mat = init.geometry, init.scalars, init.matter
        raw = vac.constraints(grid, geom)
        tot = assemble_total_stress_energy(grid, geom, scal, mat)
        out[tag] = {"H": raw["hamiltonian"] - 16.0 * math.pi * tot.rho, "M": raw["momentum"] - 8.0 * math.pi * tot.j,
                    "H_geom": np.asarray(raw["hamiltonian"]), "16pi_rho": 16.0 * math.pi * np.asarray(tot.rho),
                    "M_geom": np.asarray(raw["momentum"]), "8pi_j": 8.0 * math.pi * np.asarray(tot.j)}
    out["kernel_diag"] = kernel.diagnostics(st)
    return out


def summary(grid, ev, aux=None, ncell=6):
    r = np.asarray(grid.centers)
    H = ev["post_lapse"]["H"]
    ann = (r >= ANNULUS[0]) & (r <= ANNULUS[1])
    i = int(np.argmax(np.abs(H)))
    s = {"H01": float(np.max(np.abs(H[:2]))), "H_full_absmax": float(np.max(np.abs(H))),
         "H_annulus_absmax": float(np.max(np.abs(H[ann]))), "argmax_idx": i, "argmax_r": float(r[i]),
         "H_at_r0.3": float(np.interp(FIXED_R, r, H)),
         "M_full_absmax": float(np.max(np.abs(ev["post_lapse"]["M"]))),
         "H01_pre_lapse": float(np.max(np.abs(ev["pre_lapse"]["H"][:2]))),
         "kernel_H_full_absmax": float(ev["kernel_diag"]["hamiltonian_full_absmax"]),
         "cells": {k: [float(x) for x in ev["post_lapse"][k][:ncell]] for k in ("H", "H_geom", "16pi_rho", "M", "M_geom", "8pi_j")}}
    if aux is not None:
        s["cells"].update({k: [float(x) for x in aux[k][:ncell]] for k in ("B", "Kr") if k in aux})
        s["Kt"] = float(aux["Kt"][0])
        if "last_dB_absmax" in aux:
            s["last_dB_absmax"] = aux["last_dB_absmax"]
    return s


def even_quartic(r, f):
    A = np.vstack([np.ones(3), r[:3] ** 2, r[:3] ** 4]).T
    return [float(c) for c in np.linalg.solve(A, f[:3])]


# ---------------------------------------------------------------- inputs
def _ring(z, grid, how):
    r_t, rc = z["r"], np.asarray(grid.centers)
    out = {}
    for k in FIELDS:
        y = z[TOY_KEYS[k]]
        if how == "linear":
            v = np.interp(rc, r_t, y)                                     # W2 exactly
        else:
            v = CubicSpline(np.r_[-r_t[::-1], r_t], np.r_[y[::-1], y])(rc)  # even reconstruction
        out[k] = v / math.sqrt(LAM)
    out["phi"], out["Pi"] = np.full(len(rc), clone.PHI0), np.full(len(rc), clone.PI0)
    return out


def manufactured_params(z):
    """Rule fixed pre-run: A = ring value at cell 0 (dr 0.05), central curvature c2 from the even quartic fit on cells
    0-2, f_far = ring value at r = 40; f = f_far + (A - f_far + c2' r^2) exp(-r^2/25), c2' = c2 + (A - f_far)/25."""
    r_t = z["r"]
    p = {}
    for k in FIELDS:
        y = z[TOY_KEYS[k]] / math.sqrt(LAM)
        A, c2 = y[0], even_quartic(r_t, y)[1]
        far = float(np.interp(MAN_FAR_R, r_t, y))
        p[k] = {"A": float(A), "c2": float(c2), "far": far, "c2p": float(c2 + (A - far) / MAN_W2)}
    return p


def manufactured_fn(p):
    def fn(r):
        r = np.asarray(r, dtype=float)
        e = np.exp(-r * r / MAN_W2)
        out = {}
        for k in FIELDS:
            q = p[k]
            amp = q["A"] - q["far"] + q["c2p"] * r * r
            out[k] = q["far"] + amp * e
            out[k + "_r"] = (2.0 * q["c2p"] * r - 2.0 * r / MAN_W2 * amp) * e
        out["Sp"], out["Dp"] = out["S_r"], out["D_r"]
        out["phi"], out["Pi"] = np.full(r.shape, clone.PHI0), np.full(r.shape, clone.PI0)
        return out
    return fn


# ---------------------------------------------------------------- run
def reading(ratios, keep, drop):
    if all(x >= keep for x in ratios):
        return "KEEPS"
    if all(x <= drop for x in ratios):
        return "DROPS"
    return "MIXED"


def main(npz, outdir):
    os.makedirs(outdir, exist_ok=True)
    kernel = V55ProductionKernel()
    go, vac, _ = adapter.vendor_modules()
    z = np.load(npz)
    r0 = json.load(open("results/research/b7_ring_constraint_clone_r0/report.json"))
    rep = {"member_sha256": src.sha256_file(npz), "runtime": src.environment_report(), "lambda": LAM,
           "held_fixed": "uniform dust dm0; G4c = ring H <= 10x PC2 (same grid); Hubble K_t; PS terms on; n_iter 6; r_max 60",
           "steps": {}}
    # 2a quadrature unit test (grid only)
    s2a = {}
    for dr in DRS:
        g = go.SphericalCellGrid(int(round(R_MAX / dr)), R_MAX)
        r = np.asarray(g.centers)
        c = (1.0, -0.3, 0.02, -1e-3)
        f4 = c[0] + c[1] * r ** 2 + c[2] * r ** 4
        F4 = c[0] * r + c[1] * r ** 3 / 3 + c[2] * r ** 5 / 5
        f6 = f4 + c[3] * r ** 6
        F6 = F4 + c[3] * r ** 7 / 7
        q4 = _cumulative_origin_integral(g, r, dr, f4)
        q6 = _cumulative_origin_integral(g, r, dr, f6)
        s2a["dr%g" % dr] = {"first_centre_rel_err_r4": float(abs(q4[0] - F4[0]) / abs(F4[0])),
                            "first_centre_err_r6": float(abs(q6[0] - F6[0])),
                            "cells0_5_abs_err_r6": [float(x) for x in np.abs(q6 - F6)[:6]],
                            "spline_ref_first_centre_err_r6": float(abs(spline_even_integral(r, f6)[0] - F6[0]))}
    s2a["PASS_first_centre_exact_through_r4"] = all(v["first_centre_rel_err_r4"] <= QUAD_EXACT_RTOL for v in s2a.values())
    e6 = [s2a["dr%g" % d]["first_centre_err_r6"] for d in DRS]
    s2a["r6_first_centre_order"] = [float(math.log(e6[i] / e6[i + 1], 2)) for i in range(2)]
    rep["steps"]["2a"] = s2a
    mp = manufactured_params(z)
    fn = manufactured_fn(mp)
    rep["manufactured_params"] = mp
    S1, S2b, S2c, S3a = {}, {}, {}, {}
    for dr in DRS:
        key = "dr%g" % dr
        g = go.SphericalCellGrid(int(round(R_MAX / dr)), R_MAX)
        r = np.asarray(g.centers)
        e2 = w2.engine_slice(g, vac, "rescaled_pocket", pocket_x=0.5)
        pc2 = summary(g, evaluate(kernel, g, e2))
        pc2_ref = pc2["kernel_H_full_absmax"]
        # step 1
        f = _ring(z, g, "linear")
        a = [f[k] for k in ("S", "PS", "D", "PD", "phi", "Pi")]
        base, aux = reconstruct(g, vac.flat_state, *a)
        sb = summary(g, evaluate(kernel, g, base), aux)
        arch = r0["grids"][key]["ring"]["hubble"]["H_full_absmax"]
        cl, _ = clone.reconstruct_slice(g, vac.flat_state, *a)
        kz, _ = clone.reconstruct_slice(g, vac.flat_state, *a, Kt_mode="zero")
        S1[key] = {"ring": sb, "PC2": pc2, "archived_H_full_absmax": arch,
                   "repro_rel": abs(sb["kernel_H_full_absmax"] - arch) / arch,
                   "variant_equals_clone_bitwise": all(np.array_equal(getattr(cl.geometry, q), getattr(base.geometry, q)) for q in w2.GEOM_FIELDS),
                   "PC2_archived": r0["grids"][key]["ring"]["hubble"]["PC2_ref"]["H_full_absmax"],
                   "Kt_zero_control": summary(g, evaluate(kernel, g, kz)),
                   "ring_over_PC2": sb["H_full_absmax"] / pc2_ref, "G4c": sb["H_full_absmax"] <= G4C_FACTOR * pc2_ref}
        # 2b manufactured
        m = fn(r)
        ma = [m[k] for k in ("S", "PS", "D", "PD", "phi", "Pi")]
        ME, auxE = reconstruct(g, vac.flat_state, *ma)
        MQ, auxQ = reconstruct(g, vac.flat_state, *ma, src_derivs=(m["Sp"], m["Dp"]))
        MQ40, _ = reconstruct(g, vac.flat_state, *ma, src_derivs=(m["Sp"], m["Dp"]), n_iter=40)
        Bx = exact_B(g, fn, ODE_RTOL)
        Bx2 = exact_B(g, fn, ODE_RTOL_CHECK)
        MX, auxX = reconstruct(g, vac.flat_state, *ma, B_override=Bx)
        sE, sQ, sX = (summary(g, evaluate(kernel, g, x), y) for x, y in ((ME, auxE), (MQ, auxQ), (MX, auxX)))
        sQ40 = summary(g, evaluate(kernel, g, MQ40))
        S2b[key] = {"M_E": sE, "M_Q": sQ, "M_X": sX, "M_Q_niter40_H01": sQ40["H01"],
                    "ode_agree_rel": float(np.max(np.abs(Bx - Bx2) / np.abs(Bx))),
                    "B_X_minus_B_E_cells0_2": [float(x) for x in (Bx - auxE["B"])[:3]],
                    "ratio_X_over_E": sX["H01"] / sE["H01"], "ratio_Q_over_E": sQ["H01"] / sE["H01"],
                    "M_E_over_PC2": sE["H_full_absmax"] / pc2_ref, "M_X_over_PC2": sX["H_full_absmax"] / pc2_ref,
                    "parity_check": {"S_r_cell0_over_r": float(m["Sp"][0] / r[0]), "B_X_minus_1_over_r2_cells0_2": [float(x) for x in ((Bx - 1) / r ** 2)[:3]]}}
        # 2c quadrature swap on the W2 samples
        sw, auxS = reconstruct(g, vac.flat_state, *a, quadrature="spline")
        ss_ = summary(g, evaluate(kernel, g, sw), auxS)
        frozen = aux["src_frozen"]
        dF = spline_even_integral(r, frozen) - _cumulative_origin_integral(g, r, float(g.dr), frozen)
        S2c[key] = {"swap": ss_, "ratio_swap_over_base": ss_["H01"] / sb["H01"], "dF_frozen_cells0_2": [float(x) for x in dF[:3]],
                    "dB_rebuilt_cells0_2": [float(x) for x in (auxS["B"] - aux["B"])[:3]],
                    "G4c": ss_["H_full_absmax"] <= G4C_FACTOR * pc2_ref}
        # 3a even reconstruction of the archived arrays
        fe = _ring(z, g, "spline")
        ae = [fe[k] for k in ("S", "PS", "D", "PD", "phi", "Pi")]
        ev3, aux3 = reconstruct(g, vac.flat_state, *ae)
        s3 = summary(g, evaluate(kernel, g, ev3), aux3)
        coincide = bool(np.all(np.min(np.abs(r[:, None] - z["r"][None, :]), axis=1) <= 1e-12))
        S3a[key] = {"even_spline": s3, "ratio_spline_over_linear": s3["H01"] / sb["H01"],
                    "G4c": s3["H_full_absmax"] <= G4C_FACTOR * pc2_ref, "ring_over_PC2": s3["H_full_absmax"] / pc2_ref,
                    "direct_sampling_applicable": coincide,
                    "near_origin_fit": {k: {"linear": even_quartic(r, f[k]), "spline": even_quartic(r, fe[k])} for k in FIELDS},
                    "interp_diff_absmax_cells0_5": {k: float(np.max(np.abs(fe[k] - f[k])[:6])) for k in FIELDS},
                    "interp_diff_absmax_all": {k: float(np.max(np.abs(fe[k] - f[k]))) for k in FIELDS}}
    # gates and readings (criteria from the pre-run entry)
    rep["steps"]["1"] = S1
    rep["steps"]["1"]["REPRO_PASS"] = all(S1["dr%g" % d]["repro_rel"] <= REPRO_RTOL for d in DRS)
    ode_ok = all(S2b["dr%g" % d]["ode_agree_rel"] <= ODE_AGREE for d in DRS)
    S2b["ODE_VALID"] = ode_ok
    S2b["reading"] = reading([S2b["dr%g" % d]["ratio_X_over_E"] for d in DRS], STENCIL_KEEP, QUAD_DROP) if ode_ok else "UNTESTED"
    S2c["reading"] = reading([S2c["dr%g" % d]["ratio_swap_over_base"] for d in DRS], STENCIL_KEEP, QUAD_DROP)
    r3 = [S3a["dr%g" % d]["ratio_spline_over_linear"] for d in DRS]
    if all(x <= 1 / SAMPLING_DROP for x in r3) or all(S3a["dr%g" % d]["G4c"] for d in DRS):
        S3a["reading"] = "ACCOUNTS"
    elif all(1 / NO_EFFECT <= x <= NO_EFFECT for x in r3):
        S3a["reading"] = "NO_EFFECT"
    else:
        S3a["reading"] = "PARTIAL"
    rep["steps"]["2b"], rep["steps"]["2c"], rep["steps"]["3a"] = S2b, S2c, S3a
    rep["steps"]["3b"] = {"result": "UNTESTED", "why": "chart-consistent u=1.0144 slice needs W1 (Aster's spec); no proxy"}
    stencil = S2b["reading"] == "KEEPS" and S2c["reading"] == "KEEPS"
    quad = S2c["reading"] == "DROPS" and S2b["reading"] == "DROPS"
    samp = S3a["reading"] == "ACCOUNTS"
    named = [n for n, ok in (("stencil", stencil), ("quadrature", quad), ("sampling", samp)) if ok]
    controls_ok = rep["steps"]["1"]["REPRO_PASS"] and s2a["PASS_first_centre_exact_through_r4"]
    rep["verdict"] = {"controls_ok": controls_ok, "cause": named[0] if (controls_ok and len(named) == 1) else "OPEN",
                      "readings": {"2b": S2b["reading"], "2c": S2c["reading"], "3a": S3a["reading"], "3b": "UNTESTED"}}
    json.dump(rep, open(os.path.join(outdir, "report.json"), "w"), indent=1, default=float)
    return rep


# ---------------------------------------------------------------- junction test (pre-run log entry ~01:00 UTC Nov 23)
V_EXACT_RTOL = 1e-12
V_ORDER_MIN = 5.5
V_JUMP_MAX = 10.0
J_DROP, J_KEEP, J_MATCH_X = 0.1, 0.5, 2.0


def gate_v(go, order_rmax=None):
    """order_rmax None = r0 gate (whole domain); 5.0 = r1 gate (pre-run entry ~01:05 UTC Nov 23)."""
    out = {}
    for dr in DRS:
        g = go.SphericalCellGrid(int(round(R_MAX / dr)), R_MAX)
        r = np.asarray(g.centers)
        d = {}
        for p in (0, 2, 4):
            F = r ** (p + 1) / (p + 1)
            d["max_rel_err_r%d" % p] = float(np.max(np.abs(even6_integral(r, r ** p, dr) - F) / np.abs(F)))
        F6 = r ** 7 / 7
        e = np.abs(even6_integral(r, r ** 6, dr) - F6)
        ee = np.abs(_cumulative_origin_integral(g, r, dr, r ** 6) - F6)
        sel = np.ones_like(r, dtype=bool) if order_rmax is None else (r <= order_rmax)
        d.update(r6_max_err=float(e[sel].max()), order_rmax=order_rmax, r6_cells0_5=[float(x) for x in e[:6]],
                 engine_r6_max_err=float(ee.max()), engine_r6_cells0_5=[float(x) for x in ee[:6]])
        out["dr%g" % dr] = d
    ex = all(out[k]["max_rel_err_r%d" % p] <= V_EXACT_RTOL for k in out for p in (0, 2, 4))
    m = [out["dr%g" % d]["r6_max_err"] for d in DRS]
    order = [float(math.log(m[i] / m[i + 1], 2)) for i in range(2)]
    c = out["dr0.05"]["r6_cells0_5"]
    jump = c[1] / c[0] if c[0] > 0 else float("inf")
    out.update(exact_through_r4=ex, r6_order=order, cell1_over_cell0_dr005=jump,
               PASS=bool(ex and min(order) >= V_ORDER_MIN and jump <= V_JUMP_MAX))
    return out


def main_junction(npz, outdir, order_rmax=None):
    os.makedirs(outdir, exist_ok=True)
    kernel = V55ProductionKernel()
    go, vac, _ = adapter.vendor_modules()
    z = np.load(npz)
    rep = {"member_sha256": src.sha256_file(npz), "runtime": src.environment_report(), "gate_V": gate_v(go, order_rmax), "grids": {}}
    if not rep["gate_V"]["PASS"]:
        rep["verdict"] = "STOPPED: gate V failed, no swap run"
        json.dump(rep, open(os.path.join(outdir, "report.json"), "w"), indent=1, default=float)
        return rep
    fn = manufactured_fn(manufactured_params(z))
    for dr in DRS:
        g = go.SphericalCellGrid(int(round(R_MAX / dr)), R_MAX)
        r = np.asarray(g.centers)
        pc2 = summary(g, evaluate(kernel, g, w2.engine_slice(g, vac, "rescaled_pocket", pocket_x=0.5)))["kernel_H_full_absmax"]
        f = _ring(z, g, "linear")
        a = [f[k] for k in ("S", "PS", "D", "PD", "phi", "Pi")]
        base, ab = reconstruct(g, vac.flat_state, *a)
        q6, aq = reconstruct(g, vac.flat_state, *a, quadrature="even6")
        sb, sq = summary(g, evaluate(kernel, g, base), ab), summary(g, evaluate(kernel, g, q6), aq)
        m = fn(r)
        ma = [m[k] for k in ("S", "PS", "D", "PD", "phi", "Pi")]
        mE, aE = reconstruct(g, vac.flat_state, *ma)
        mQ, aQ = reconstruct(g, vac.flat_state, *ma, quadrature="even6")
        mX, aX = reconstruct(g, vac.flat_state, *ma, B_override=exact_B(g, fn, ODE_RTOL))
        sE, sQ6, sX = (summary(g, evaluate(kernel, g, x), y) for x, y in ((mE, aE), (mQ, aQ), (mX, aX)))
        e2 = w2.engine_slice(g, vac, "rescaled_pocket", pocket_x=0.5).scalars
        p6, _ = reconstruct(g, vac.flat_state, e2.S, e2.PS, e2.D, e2.PD, e2.phi, e2.Pi, quadrature="even6")
        rep["grids"]["dr%g" % dr] = {
            "PC2_engine_H": pc2, "PC2_rebuilt_Q6_H": summary(g, evaluate(kernel, g, p6))["kernel_H_full_absmax"],
            "ring_engine": sb, "ring_Q6": sq, "ring_ratio_Q6_over_engine": sq["H01"] / sb["H01"],
            "ring_engine_over_PC2": sb["H_full_absmax"] / pc2, "ring_Q6_over_PC2": sq["H_full_absmax"] / pc2,
            "G4c_Q6": sq["H_full_absmax"] <= G4C_FACTOR * pc2,
            "dB_Q6_minus_engine_cells0_2": [float(x) for x in (aq["B"] - ab["B"])[:3]],
            "man_E_H01": sE["H01"], "man_Q6_H01": sQ6["H01"], "man_X_H01": sX["H01"],
            "man_Q6_over_X": sQ6["H01"] / sX["H01"], "man_Q6_cells": sQ6["cells"]["H"], "man_X_cells": sX["cells"]["H"]}
    G = [rep["grids"]["dr%g" % d] for d in DRS]
    rr = [x["ring_ratio_Q6_over_engine"] for x in G]
    if all(x <= J_DROP for x in rr) and all(x["man_Q6_over_X"] <= J_MATCH_X for x in G):
        v = "JUNCTION IMPLICATED"
    elif all(x >= J_KEEP for x in rr):
        v = "NOT THE JUNCTION"
    else:
        v = "MIXED"
    rep["verdict"] = {"reading": v, "G4c_Q6_all_grids": all(x["G4c_Q6"] for x in G)}
    json.dump(rep, open(os.path.join(outdir, "report.json"), "w"), indent=1, default=float)
    return rep


# ---------------------------------------------------------------- first-piece test (pre-run entry ~01:20 UTC Nov 23)
EPS = np.finfo(float).eps
G1_RTOL, G2_ORDER_MIN, G2_FLOOR, G3_FLOOR, G4_FLOOR = 1e-12, 8.5, 100.0, 100.0, 10.0
R1_MATCH_X, R1_KEEP, R2_DROP, R2_KEEP = 2.0, 0.5, 0.1, 0.5


def first_piece_gates(go):
    from scipy.special import erf
    g1 = {}
    for dr in DRS:
        r = (np.arange(8) + 0.5) * dr
        g1["dr%g" % dr] = {"r%d" % p: float(abs(first_s6(r, r ** p, dr) - r[0] ** (p + 1) / (p + 1)) / (r[0] ** (p + 1) / (p + 1)))
                           for p in (0, 2, 4, 6)}
    g1_pass = all(v <= G1_RTOL for d in g1.values() for v in d.values())
    g2 = {}
    for dr in (0.4, 0.2, 0.1):
        r = (np.arange(8) + 0.5) * dr
        f = np.exp(-4.0 * r * r)
        F0 = math.sqrt(math.pi) / 4.0 * erf(2.0 * r[0])
        e6, e4 = abs(first_s6(r, f, dr) - F0), abs(first_q4(r, f) - F0)
        g2["dr%g" % dr] = {"err_fp6": e6, "err_engine_q4": e4, "fp6_over_floor": e6 / (EPS * abs(F0)),
                           "discriminating": e6 >= G2_FLOOR * EPS * abs(F0)}
    keep = [d for d in (0.4, 0.2, 0.1) if g2["dr%g" % d]["discriminating"]]
    o6 = [math.log(g2["dr%g" % keep[i]]["err_fp6"] / g2["dr%g" % keep[i + 1]]["err_fp6"], keep[i] / keep[i + 1]) for i in range(len(keep) - 1)]
    o4 = [math.log(g2["dr%g" % a]["err_engine_q4"] / g2["dr%g" % b]["err_engine_q4"], a / b) for a, b in ((0.4, 0.2), (0.2, 0.1))]
    g2_state = "UNTESTED" if len(o6) < 2 else ("PASS" if min(o6) >= G2_ORDER_MIN else "FAIL")
    return {"G1": g1, "G1_PASS": g1_pass, "G2": g2, "G2_orders_fp6": o6, "G2_orders_engine_q4": o4, "G2_state": g2_state}


def exact_first_piece(fn, r0, B0):
    """Oracle: quad of the continuum source over [0, r_0]; B(r) = 1 + (B0 - 1)(r/r_0)^2 (enters only via 8 pi r^2 rho_1 B)."""
    from scipy.integrate import quad
    H0 = clone.background_H0(True, True)
    Kt = -H0

    def srcf(t):
        f = fn(np.array([t]))
        B = 1.0 + (B0 - 1.0) * (t / r0) ** 2
        rho = _rho(f["S"], f["PS"], f["D"], f["PD"], f["Sp"], f["Dp"], f["phi"], f["Pi"], B, RHO_R_PRESENT)[0]
        jj = -(f["PD"][0] * f["Dp"][0]) - f["PS"][0] * f["Sp"][0]
        Kr = Kt + 4.0 * math.pi * t * jj
        return 1.0 - 8.0 * math.pi * t * t * rho + 2.0 * t * t * Kr * Kt + t * t * Kt * Kt

    val, err = quad(srcf, 0.0, r0, epsabs=1e-30, epsrel=1.2e-14, limit=200)
    return val, err


def main_first_piece(npz, outdir):
    os.makedirs(outdir, exist_ok=True)
    kernel = V55ProductionKernel()
    go, vac, _ = adapter.vendor_modules()
    z = np.load(npz)
    rep = {"member_sha256": src.sha256_file(npz), "runtime": src.environment_report(), "gates": first_piece_gates(go), "grids": {}}
    run_fp6 = rep["gates"]["G1_PASS"] and rep["gates"]["G2_state"] == "PASS"
    fn = manufactured_fn(manufactured_params(z))
    for dr in DRS:
        key = "dr%g" % dr
        g = go.SphericalCellGrid(int(round(R_MAX / dr)), R_MAX)
        r = np.asarray(g.centers)
        e2 = w2.engine_slice(g, vac, "rescaled_pocket", pocket_x=0.5)
        pc2 = summary(g, evaluate(kernel, g, e2))["kernel_H_full_absmax"]
        s2 = e2.scalars
        p6, _ = reconstruct(g, vac.flat_state, s2.S, s2.PS, s2.D, s2.PD, s2.phi, s2.Pi, quadrature="even6")
        floor = summary(g, evaluate(kernel, g, p6))["kernel_H_full_absmax"]
        m = fn(r)
        ma = [m[k] for k in ("S", "PS", "D", "PD", "phi", "Pi")]
        mE, aE = reconstruct(g, vac.flat_state, *ma)
        Bx = exact_B(g, fn, ODE_RTOL)
        mX, aX = reconstruct(g, vac.flat_state, *ma, B_override=Bx)
        fx, fx_err = exact_first_piece(fn, float(r[0]), float(Bx[0]))
        fE = float(first_q4(r, aE["src_frozen"]))
        mO, aO = reconstruct(g, vac.flat_state, *ma, quadrature="fpo", first_override=fx)
        sE, sX, sO = (summary(g, evaluate(kernel, g, x), y) for x, y in ((mE, aE), (mX, aX), (mO, aO)))
        G = {"PC2_engine_H": pc2, "floor_PC2_Q6_H": floor,
             "G3_first_exact": fx, "G3_quad_err_est": fx_err, "G3_first_engine": fE, "G3_diff": fx - fE,
             "G3_discriminating": abs(fx - fE) >= G3_FLOOR * EPS * abs(fx) and fx_err < 0.01 * abs(fx - fE),
             "G4_discriminating": sX["H01"] >= G4_FLOOR * floor,
             "man_E": sE, "man_X": sX, "man_FPO": sO, "man_FPO_over_X": sO["H01"] / sX["H01"], "man_FPO_over_E": sO["H01"] / sE["H01"],
             "dB_FPO_minus_X_cells0_2": [float(x) for x in (aO["B"] - Bx)[:3]], "dB_E_minus_X_cells0_2": [float(x) for x in (aE["B"] - Bx)[:3]]}
        if run_fp6:
            f = _ring(z, g, "linear")
            a = [f[k] for k in ("S", "PS", "D", "PD", "phi", "Pi")]
            rb, ab = reconstruct(g, vac.flat_state, *a)
            r6, a6 = reconstruct(g, vac.flat_state, *a, quadrature="fp6")
            m6, am6 = reconstruct(g, vac.flat_state, *ma, quadrature="fp6")
            sb, s6, sm6 = summary(g, evaluate(kernel, g, rb), ab), summary(g, evaluate(kernel, g, r6), a6), summary(g, evaluate(kernel, g, m6), am6)
            G.update(ring_E=sb, ring_FP6=s6, ring_FP6_over_E=s6["H01"] / sb["H01"], ring_FP6_over_PC2=s6["H_full_absmax"] / pc2,
                     G4c_FP6=s6["H_full_absmax"] <= G4C_FACTOR * pc2, man_FP6=sm6, man_FP6_over_X=sm6["H01"] / sX["H01"],
                     dB_ring_FP6_minus_E_cells0_2=[float(x) for x in (a6["B"] - ab["B"])[:3]])
        rep["grids"][key] = G
    Gs = [rep["grids"]["dr%g" % d] for d in DRS]
    disc = all(x["G3_discriminating"] and x["G4_discriminating"] for x in Gs)
    if not disc:
        r1 = "UNTESTED"
    elif all(x["man_FPO"]["H01"] <= R1_MATCH_X * x["man_X"]["H01"] for x in Gs):
        r1 = "FIRST PIECE ACCOUNTS"
    elif all(x["man_FPO_over_E"] >= R1_KEEP for x in Gs):
        r1 = "DOES NOT"
    else:
        r1 = "MIXED"
    if not run_fp6:
        r2 = "UNTESTED (FP6 gates)"
    else:
        rr = [x["ring_FP6_over_E"] for x in Gs]
        r2 = "FIX WORKS" if all(v <= R2_DROP for v in rr) else ("NO EFFECT" if all(v >= R2_KEEP for v in rr) else "PARTIAL")
    rep["verdict"] = {"R1": r1, "R2": r2, "cause": "FIRST PIECE" if r1 == "FIRST PIECE ACCOUNTS" else "OPEN",
                      "G4c_FP6_all_grids": all(x.get("G4c_FP6", False) for x in Gs) if run_fp6 else None}
    json.dump(rep, open(os.path.join(outdir, "report.json"), "w"), indent=1, default=float)
    return rep


if __name__ == "__main__":
    if len(sys.argv) > 3 and sys.argv[3] == "first_piece":
        main_first_piece(*sys.argv[1:3])
        sys.exit(0)
    if len(sys.argv) > 3 and sys.argv[3] == "junction_r1":
        main_junction(*sys.argv[1:3], order_rmax=5.0)
    elif len(sys.argv) > 3 and sys.argv[3] == "junction":
        main_junction(*sys.argv[1:3])
    else:
        main(*sys.argv[1:3])
