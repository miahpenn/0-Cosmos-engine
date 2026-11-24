"""C9 W1 reruns (gates: results/research/w1_c9_reruns/PRERUN.md). Zero-step slices only.
Usage (repo root, OPENBLAS_NUM_THREADS=1): python3 research/w1_c9_reruns/run.py native.npz members.json OUTDIR
"""
import hashlib
import json
import math
import os
import sys
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.getcwd())
import numpy as np  # noqa: E402

from engine import scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402
from engine.stress_energy import assemble_total_stress_energy  # noqa: E402
from research import b7_w2_origin_check_v1 as oc  # noqa: E402

npz_path, members_path, outdir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
outdir.mkdir(parents=True, exist_ok=True)
z = np.load(npz_path)
members = json.loads(members_path.read_text())
go, vac, _ = ad.vendor_modules()
kernel = V55ProductionKernel()
PHI0, PI0 = oc.clone.PHI0, oc.clone.PI0
GRIDS = ((0.05, 1200), (0.025, 2400))
E_RING, KAPPA = 1.1525e-4, 0.630519
W = 1.0 / KAPPA
F_REF = {0.025: 2.906e-12, 0.05: 8.268e-13}
EC_REF = {0.025: 1.4102e-10, 0.05: 2.2306e-9}
K1_TOL = {0.025: 2e-7, 0.05: 3.2e-6}
K3_TOL = {0.025: 0.03, 0.05: 8.268e-13 / 2.2306e-9 + 0.009}
JS574 = {0.05: 1.552569083079575e-09, 0.025: 1.549082904208777e-09}
INTEGRAL = {"B": "_eight_node_origin_integral", "A": "_cumulative_origin_integral"}
rep = {"native_sha256": hashlib.sha256(npz_path.read_bytes()).hexdigest(), "grids": {}}


def sech2(M):
    return lambda r: 3 * KAPPA**3 * M / math.pi**3 / np.cosh(KAPPA * r) ** 2


def gauss(M):
    return lambda r: M / (math.pi**1.5 * W**3) * np.exp(-(r / W) ** 2)


def build(grid, run, structure, inj=None, count=None, **kw):
    name = INTEGRAL[run]
    real = getattr(vi, name)

    def wrapped(g, r, h, f):
        if inj is not None:
            f = f - 8.0 * math.pi * r * r * inj(r)
        return real(g, r, h, f)
    patches = [mock.patch.object(vi, name, wrapped)]
    if count is not None:
        realp = ss.scalar_momentum_projection

        def proj(*a):
            count.append(1)
            return realp(*a)
        patches.append(mock.patch.object(vi._scalar_system, "scalar_momentum_projection", proj))
    with patches[0]:
        if len(patches) > 1:
            with patches[1]:
                return vi.build_initial_data(grid, vac.flat_state, initial_structure=structure,
                                             eight_node_areal=(run == "B"), **kw)
        return vi.build_initial_data(grid, vac.flat_state, initial_structure=structure,
                                     eight_node_areal=(run == "B"), **kw)


def H_of(grid, init):
    return np.asarray(oc.evaluate(kernel, grid, init)["post_lapse"]["H"], dtype=float)


def first(v, h, n):
    N = len(v)
    i = np.arange(N)
    start = np.minimum(i - n // 2, N - n)
    ix = start[:, None] + np.arange(n)
    vals = v[np.where(ix < 0, -ix - 1, ix)] - v[:, None]
    d = np.empty(N)
    for off in np.unique(i - start):
        sel = i - start == off
        wts = np.linalg.solve(np.vander(np.arange(n) - off, n, increasing=True).T, np.eye(n)[1])
        d[sel] = vals[sel] @ wts / h
    return d


def direct(grid, init, h):
    # As Aster's full.py: geometry of the initial slice, 16 pi rho from the post-lapse readout.
    e = oc.evaluate(kernel, grid, init)["post_lapse"]
    geom = init.geometry
    r = np.asarray(grid.centers)
    B = geom.X ** 2 / geom.a
    terms = -1.5 * geom.Aa ** 2 + 2 / 3 * geom.K ** 2 - e["16pi_rho"]
    return {str(n): float(np.max(np.abs(2 * ((1 - B) / r - first(B, h, n)) / r + terms))) for n in (7, 9)}


u = (members["dr0.05"]["drhat"] / 0.05) ** 2
_r05 = np.asarray(go.SphericalCellGrid(1200, 60.0).centers)
MAN_P = {}
for _k in oc.FIELDS:
    _y = np.asarray(z["dr0.05_" + _k], dtype=float)
    _A, _c2 = _y[0], oc.even_quartic(_r05, _y)[1]
    _far = float(np.interp(oc.MAN_FAR_R, _r05, _y))
    MAN_P[_k] = {"A": float(_A), "c2": float(_c2), "far": _far, "c2p": float(_c2 + (_A - _far) / oc.MAN_W2)}
fnm = oc.manufactured_fn(MAN_P)
for dr, n in GRIDS:
    key = "dr%g" % dr
    grid = go.SphericalCellGrid(n, 60.0)
    r = np.asarray(grid.centers)
    G = {}
    none_B = build(grid, "B", "none")
    H0 = float(none_B.H0)
    t = key + "_"
    S, PS, D, PD = (np.array(z[t + k], dtype=np.float64) for k in ("S", "PS", "D", "PD"))
    phi, Pi = np.full(n, PHI0), np.full(n, PI0)
    m = members[key]
    label = (f"Aster B7-W1-grid-aligned-comparison.zip native.npz sha256 {rep['native_sha256']}; u={u:.6f}; theta0=0; "
             f"snapshot phase 0; matching: harmonic-balance M2 Newton solve on engine-aligned nodes (seed: even cubic of "
             f"archived harmonics), Hubble chart map; toy dr={m['drhat']!r}; toy matched outer radius={m['Rhat']!r}; "
             f"H0_chart={H0!r}")
    sl = vi.SuppliedSlice(ss.ScalarFields(S, PS, D, PD, phi, Pi), np.array(r), H0, 1e-7, label)

    # ---- G0
    Drest = float(ss.sd_rest_value(PHI0))
    edge = {}
    for nm, arr in (("S", S), ("PS", PS), ("PD", PD), ("D-D_rest", D - Drest)):
        top = float(np.max(np.abs(arr)))
        edge[nm] = {"edge_max": float(np.max(np.abs(arr[-5:]))), "grid_max": top,
                    "ratio": float(np.max(np.abs(arr[-5:]))) / top}
    G["G0"] = {"G0.1_label": label, "G0.1_PASS": all(s in label for s in ("u=", "theta0", "phase", "matching", "toy dr",
                                                                          "outer radius", "H0_chart")),
               "G0.2_PASS": bool(np.array_equal(sl.r, r)), "G0.2_note": "producer nodes within 2e-14 of centres (Aster, Read)",
               "G0.3_H0_chart": H0, "G0.3_PASS": H0 != 0.0 and abs(sl.H0_chart - H0) <= 1e-12 * H0,
               "G0.4_PASS": bool(np.all(phi == PHI0) and np.all(Pi == PI0)),
               "G0.5": edge, "G0.5_PASS": all(v["ratio"] <= 1e-6 for v in edge.values())}

    runs = {}
    for run in ("B", "A"):
        R_ = {}
        cnt = []
        ring = build(grid, run, "supplied_slice", count=cnt, supplied_slice=sl)
        none = none_B if run == "B" else build(grid, run, "none")
        Hr, Hn = H_of(grid, ring), H_of(grid, none)
        order = np.argsort(-np.abs(Hr))
        R, F = float(abs(Hr[order[0]])), float(np.max(np.abs(Hn)))
        R_["ring"] = {"R": R, "argmax": int(order[0]), "second": int(order[1]), "second_val": float(abs(Hr[order[1]])),
                      "areal_last_step": float(ring.areal_last_step), "momentum_projection_calls": len(cnt)}
        R_["F"] = F
        K = {}
        for tag, prof in (("sech2_k10", sech2(10 * E_RING)), ("sech2_k100", sech2(100 * E_RING)),
                          ("gauss_k10", gauss(10 * E_RING))):
            init = build(grid, run, "none", inj=prof)
            H = H_of(grid, init)
            Hc = 16 * math.pi * prof(r)
            err = np.abs(H - Hc)
            K[tag] = {"e_c": float(err.max()), "e_c_cell": int(np.argmax(err)), "argmax_H": int(np.argmax(np.abs(H))),
                      "Hc_max": float(Hc.max()), "rel": float(err.max() / Hc.max()),
                      "e_at_ring_argmax": float(err[order[0]])}
        R_["K"] = K
        e10, e100, eg = K["sech2_k10"]["e_c"], K["sech2_k100"]["e_c"], K["gauss_k10"]["e_c"]
        C, Cg = e10 / 10, eg / 10
        chk = {"K0": abs(F / F_REF[dr] - 1) <= 0.10,
               "K1": K["sech2_k10"]["argmax_H"] == 0 and K["sech2_k10"]["rel"] <= K1_TOL[dr],
               "K2": e10 >= 20 * F, "K3_value": abs(10 * e10 / e100 - 1), "K3_tol": K3_TOL[dr],
               "K5": abs(e10 / EC_REF[dr] - 1) <= 0.02}
        chk["K3"] = chk["K3_value"] <= K3_TOL[dr]
        R_["checks"] = chk
        R_["C"], R_["C_gauss"] = C, Cg
        R_["R_over_C"], R_["R_over_C_gauss"], R_["R_over_F"] = R / C, R / Cg, R / F
        R_["sep_e10_over_F"] = e10 / F
        if R + F <= 10 * C:
            verdict = "PASS"
        elif R - F > 10 * C:
            verdict = "FAIL"
        else:
            verdict = "UNDETERMINED"
        kfail = [k for k in ("K0", "K1", "K2", "K3", "K5") if not chk[k]]
        R_["G4c_rule"] = verdict
        R_["G4c_reading"] = verdict if not kfail else f"UNDETERMINED (K failed: {','.join(kfail)}; rule alone {verdict})"
        # old control and direct curvature, recorded only
        pocket = build(grid, run, "rescaled_pocket", pocket_x=0.5)
        P = float(np.max(np.abs(H_of(grid, pocket))))
        dr_ring, dr_pocket = direct(grid, ring, dr), direct(grid, pocket, dr)
        R_["old_control"] = {"pocket_max": P, "ring_over_pocket": R / P,
                             "direct_ring_over_pocket": {k: dr_ring[k] / dr_pocket[k] for k in dr_ring}}
        # G1 (Run B and A both recorded; gates read on Run B)
        st = kernel.initialize(resolution=n, r_max=60.0, initial_structure="supplied_slice", supplied_slice=sl,
                               eight_node_areal=(run == "B"))
        mom = vac.constraints(st.grid, st.geometry)["momentum"]
        tot = assemble_total_stress_energy(st.grid, st.geometry, st.scalars, st.matter)
        res = float(np.max(np.abs(mom - 8 * np.pi * tot.j)))
        Sp = vi._d1(st.grid, st.scalars.S, 1)
        jS = -st.scalars.PS * Sp
        res_noS = float(np.max(np.abs(mom - 8 * np.pi * (tot.j - jS))))
        mjS = float(np.max(np.abs(jS)))
        R_["G1"] = {"G1.1_residual": res, "G1.1_threshold": 1e-6 * 8 * np.pi * mjS,
                    "G1.1_PASS": res <= 1e-6 * 8 * np.pi * mjS,
                    "G1.2_ratio": res_noS / res if res else float("inf"), "G1.2_PASS": res_noS >= 1e3 * res,
                    "G1.3_max_jS": mjS, "G1.3_rel": abs(mjS / JS574[dr] - 1), "G1.3_PASS": abs(mjS / JS574[dr] - 1) <= 1e-12,
                    "G1.4_last_dB": float(ring.areal_last_step), "G1.4_PASS": ring.areal_last_step <= 1e-13,
                    "G1.5_projection_calls": len(cnt), "G1.5_PASS": len(cnt) == 1}
        # G2 (recorded on both runs; read on Run B)
        rho_r = assemble_total_stress_energy(grid, ring.geometry, ring.scalars, ring.matter).rho
        rho_n = assemble_total_stress_energy(grid, none.geometry, none.scalars, none.matter).rho
        Bg = ring.geometry.X ** 2 / ring.geometry.a
        dE = 4 * math.pi * r * r * (np.asarray(rho_r) - np.asarray(rho_n)) * dr
        inner = r < 25.0
        E = {"flat_r<25": float(dE[inner].sum()), "flat_box": float(dE.sum()),
             "areal_r<25": float((dE / np.sqrt(Bg))[inner].sum()), "areal_box": float((dE / np.sqrt(Bg)).sum())}
        R_["G2"] = {"E_ring": E, "E_ring_rel_to_1.1525e-4": {k: v / E_RING - 1 for k, v in E.items()},
                    "E_ring_PASS_flat_r<25": abs(E["flat_r<25"] / E_RING - 1) <= 0.01,
                    "max_abs_PD_field": float(np.max(np.abs(PD))), "PD_note": "field amplitude vs toy radiated power: not comparable",
                    "f_ring": "NOT COMPUTED (definition not available)"}
        runs[run] = R_
        print(key, run, verdict, R_["G4c_reading"], "R/C=%.4g" % (R / C), "F=%.4g e10=%.4g" % (F, e10), flush=True)
    G["runs"] = runs
    # G4 coherent B (#568 method on the aligned member), Run A and Run B
    fields = dict(S=S, PS=PS, D=D, PD=PD, phi=phi, Pi=Pi)
    a = [fields[k] for k in ("S", "PS", "D", "PD", "phi", "Pi")]
    B5 = oc.exact_B(grid, oc.spline_fields_fn(r, fields, 5), oc.CB_ODE[0])
    B5b = oc.exact_B(grid, oc.spline_fields_fn(r, fields, 5), oc.CB_ODE[1])
    B3 = oc.exact_B(grid, oc.spline_fields_fn(r, fields, 3), oc.CB_ODE[0])
    uB = np.abs(B5 - B3) + np.abs(B5 - B5b)
    G4 = {}
    for run in ("A", "B"):
        real = getattr(vi, INTEGRAL[run])
        with mock.patch.object(oc, "_cumulative_origin_integral", real), mock.patch.object(vi, "_cumulative_origin_integral", real):
            def arm(B=None):
                init, aux = oc.reconstruct(grid, vac.flat_state, *a, B_override=B)
                return H_of(grid, init)[:2], aux
            mm = fnm(r)
            ma = [mm[k] for k in ("S", "PS", "D", "PD", "phi", "Pi")]

            def marm(B=None):
                init, _ = oc.reconstruct(grid, vac.flat_state, *ma, B_override=B)
                return H_of(grid, init)[:2]
            hXm = marm(oc.exact_B(grid, fnm, oc.CB_ODE[0]))
            hQm = marm(oc.exact_B(grid, oc.spline_fields_fn(r, mm, 5), oc.CB_ODE[0]))
            hEm = marm()
            V_ref, V_eng = float(np.linalg.norm(hQm - hXm)), float(np.linalg.norm(hEm - hXm))
            hE, aux = arm()
            BE = aux["B"]
            eng = build(grid, run, "supplied_slice", supplied_slice=sl)
            Beng = eng.geometry.X ** 2 / eng.geometry.a
            hF, hF3, hFb = arm(B5)[0], arm(B3)[0], arm(B5b)[0]
            dEF = float(np.linalg.norm(hE - hF))
            u_h = float(np.linalg.norm(hF - hF3) + np.linalg.norm(hF - hFb))
            sw = {}
            for tag, cells in (("S0", [0]), ("SO", list(range(1, 6))), ("SB", list(range(6)))):
                hx = arm(oc.swap(BE, B5, cells))[0]
                sw[tag] = float(1.0 - np.linalg.norm(hx - hF) / dEF) if dEF else float("nan")
            c40 = slice(0, 41)
            dB = np.abs(BE - B5)
            ring_imprint = np.abs(Beng - (build(grid, run, "none").geometry.X ** 2 / build(grid, run, "none").geometry.a))
            G4[run] = {"gate_V": V_ref <= oc.CB_V * V_eng, "V_ref_err": V_ref, "V_engine_err": V_eng,
                       "research_vs_engine_B_maxdiff": float(np.max(np.abs(BE - Beng))),
                       "dEF": dEF, "u_h": u_h, "gate_U": u_h <= oc.CB_U * dEF, "U_ratio": u_h / dEF if dEF else None,
                       "share_S0": sw["S0"], "share_cells1_5": sw["SO"], "share_SB": sw["SB"],
                       "max_dB_0_40": float(dB[c40].max()), "max_uB_0_40": float(uB[c40].max()),
                       "cells_dB_over_uB_0_40": int(np.sum(dB[c40] > uB[c40])),
                       "G4.4_PASS": float(dB[c40].max()) <= float(uB[c40].max()),
                       "G4.3_defect_over_imprint_1_40": float(np.max(dB[1:41] / np.maximum(ring_imprint[1:41], 1e-300))),
                       "dB_cells0_5": [float(x) for x in (BE - B5)[:6]], "uB_cells0_5": [float(x) for x in uB[:6]]}
            print(key, "G4", run, {k: G4[run][k] for k in ("max_dB_0_40", "max_uB_0_40", "share_cells1_5", "U_ratio")}, flush=True)
    G["G4"] = G4
    rep["grids"][key] = G

rep["G4.2_RunA_share_change"] = abs(rep["grids"]["dr0.05"]["G4"]["A"]["share_cells1_5"] - rep["grids"]["dr0.025"]["G4"]["A"]["share_cells1_5"])
rep["G4.2_RunA_ASYMPTOTIC"] = rep["G4.2_RunA_share_change"] <= 0.1 and all(rep["grids"][k]["G4"]["A"]["gate_U"] for k in ("dr0.05", "dr0.025"))
rep["K4_ratio"] = {run: rep["grids"]["dr0.05"]["runs"][run]["K"]["sech2_k10"]["e_c"] / rep["grids"]["dr0.025"]["runs"][run]["K"]["sech2_k10"]["e_c"] for run in ("A", "B")}
rep["K4_PASS"] = {run: 12 <= v <= 20 for run, v in rep["K4_ratio"].items()}
(outdir / "report.json").write_text(json.dumps(rep, indent=1, default=float))
print("K4", rep["K4_ratio"], "G4.2", rep["G4.2_RunA_share_change"], rep["G4.2_RunA_ASYMPTOTIC"])
