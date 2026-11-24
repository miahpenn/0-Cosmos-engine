"""C9 W1 input-check rerun per AMENDMENT 2 (Drive 1bpHE2NiBAoYmbMiOOakOqYhBrFARJLd1PsaP6baWmb8).

Scope (Amendment 2, 5.1): input checks only, both grids, evaluated first: G0.1-G0.4 (incl. the engine's
supplied_slice entry with edge_tolerance per Lumen²'s ruling) and G0.5'a-d; G0.5'e (direction, AMENDMENT 3
Drive 1E6gfRSm6XvoDd8bbZqdLwNfFNnFaAiocBHQLYUGNYCk, section 9) for runs declared after Nov 24 02:34 UTC; #577 ran without it. If native.npz (sha256) and the engine
build are unchanged from #576 and every input check passes, the #576 readings stand (5.2); they are not recomputed.
Corrected (static-time) slice, W1 derivation AMENDMENTS 3 and 4 (Drive 1JU2EhMUs7dSio23VVu0WIPMt33iAY6u9gbD9EKE64aY,
1J_dGcUVzAy9jceOwxByyvuQ1pyix9hHiC9owyoRNItM), with --static BUILD_JSON from research/w1_static_slice/build.py:
G0.1 label carries the A3 3.7 fields, G0.3 checks H0_phase, 'd uses the A3 4.5 PD bound, and G0.5'f (static-frame
flux ratio in [0.99, 1.01], A3 4.3, adopted Nov 24 03:26 UTC) is gated. Without --static, 'f is recorded only.
Usage: python3 research/w1_c9_g05p/run.py native.npz members.json OUTDIR [--static build.json]
"""
import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.getcwd())
import numpy as np  # noqa: E402

from engine import scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad  # noqa: E402
from research import b7_w2_origin_check_v1 as oc  # noqa: E402

npz_path, members_path, outdir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
STATIC = json.loads(Path(sys.argv[5]).read_text()) if len(sys.argv) > 5 and sys.argv[4] == "--static" else None
outdir.mkdir(parents=True, exist_ok=True)
z = np.load(npz_path)
members = json.loads(members_path.read_text())
go, vac, _ = ad.vendor_modules()
PHI0, PI0 = oc.clone.PHI0, oc.clone.PI0

# ---- declared before the run (PRERUN.md)
NATIVE_SHA_576 = "16151c445c9dcfcc34cf61e28e268e991ab2bfea629beaa3aa52fadc4c5a98ea"
ENGINE_BLOBS_576 = {"engine/v55_initial.py": "d613ab44efa6919c6de208a45ac0c2ee8559d25a",
                    "engine/production_kernel.py": "849014d34c532b1f142def2b9b30ae02809e81d7"}
GRIDS = ((0.025, 2400), (0.05, 1200))
EDGE_TOL = {0.025: 1.943e-7, 0.05: 1.949e-7}   # Lumen² ruling: 2x model PD envelope; set after fault-test baseline
OMEGA_U = 3.818758 / 2
OM2, K2 = 2 * OMEGA_U, 3.543171
OM3, K3 = 3 * OMEGA_U, 5.363607
P_D_REF = 4.207e-10
A_D_PRED = 2.2245e-6
ENV = {0.025: (3.708e-8, 9.693e-8), 0.05: (3.709e-8, 9.696e-8)}
if STATIC:   # A3 4.5: last-cell PD envelope factor (1 - v k/Omega + v^2)
    ENV = {0.025: (3.708e-8, 1.1331e-7), 0.05: (3.709e-8, 1.1333e-7)}
F_BAND = (0.99, 1.01)                    # A3 4.3
E_EXPECT_STATIC = {"D": 0.756, "S": 0.713}  # A3 4.4, recorded beside 'e
U, E_RING = 1.0144275, 1.1525e-4


def fit(r, a2, k):
    M = np.stack([np.ones_like(r), np.cos(2 * k * r), np.sin(2 * k * r)], 1)
    A, c, s = np.linalg.lstsq(M, a2, rcond=None)[0]
    return float(A), float(math.hypot(c, s) / A), float(math.atan2(s, c))


rep = {"native_sha256": hashlib.sha256(npz_path.read_bytes()).hexdigest(),
       "engine_blobs": {p: subprocess.check_output(["git", "hash-object", p], text=True).strip() for p in ENGINE_BLOBS_576},
       "edge_tolerance": {str(k): v for k, v in EDGE_TOL.items()}, "grids": {}}
rep["unchanged_from_576"] = rep["native_sha256"] == NATIVE_SHA_576 and rep["engine_blobs"] == ENGINE_BLOBS_576
rep["slice"] = "static (W1 A3/A4)" if STATIC else "admitted"

for dr, n in GRIDS:
    key = "dr%g" % dr
    grid = go.SphericalCellGrid(n, 60.0)
    r = np.asarray(grid.centers)
    H0 = float(vi.build_initial_data(grid, vac.flat_state, initial_structure="none", eight_node_areal=True).H0)
    S, PS, D, PD = (np.array(z[key + "_" + k], dtype=np.float64) for k in ("S", "PS", "D", "PD"))
    phi, Pi = np.full(n, PHI0), np.full(n, PI0)
    Ds = float(ss.sd_rest_value(PHI0))
    m = members[key]
    u_lab = (m["drhat"] / dr) ** 2
    label = (f"Aster B7-W1-grid-aligned-comparison.zip native.npz sha256 {rep['native_sha256']}; u={u_lab:.6f}; "
             f"theta0=0; snapshot phase 0; matching: harmonic-balance M2 Newton solve on engine-aligned nodes; "
             f"Hubble chart map; toy dr={m['drhat']!r}; toy matched outer radius={m['Rhat']!r}; H0_chart={H0!r}")
    if STATIC:
        gsha = hashlib.sha256(b"".join(np.ascontiguousarray(a).tobytes() for a in (S, PS, D, PD))).hexdigest()
        H0_phase = float(STATIC["settings"]["H0"])
        label = (f"W1 corrected slice (A3/A4) file sha256 {rep['native_sha256']}; grid arrays sha256 {gsha}; "
                 f"slice_time=static (T = t + H0 R^2/2); H0_phase={H0_phase!r}; u={u_lab:.6f}; theta0=0; "
                 f"snapshot phase omega_u H0 r^2/2 per solver node; matching: Aster's M2 member, stored HB vectors "
                 f"x through the producer's mapped(); harmonics S {STATIC['harmonics']['S']} D {STATIC['harmonics']['D']}; "
                 f"control theta=0 bitwise vs native.npz {STATIC['native_npz_sha256']}: "
                 f"{'PASS' if STATIC['control_PASS'] else 'FAIL'}; toy dr={m['drhat']!r}; "
                 f"toy matched outer radius={m['Rhat']!r}; H0_chart={H0!r}")
    sl = vi.SuppliedSlice(ss.ScalarFields(S, PS, D, PD, phi, Pi), np.array(r), H0, EDGE_TOL[dr], label)
    G = {}
    # ---- G0.1-G0.4 and the engine entry
    edge = {nm: float(np.max(np.abs(a[-5:]))) for nm, a in (("D-D_s", D - Ds), ("S", S), ("PS", PS), ("PD", PD))}
    over = [nm for nm, v in edge.items() if v > EDGE_TOL[dr]]
    try:
        vi._supplied_arrays(sl, r, PHI0, PI0, H0)
        entry = "ACCEPTED"
    except ValueError as e:
        entry = "REFUSED: " + str(e)
    G["G0.1_PASS"] = all(s in label for s in ("u=", "theta0", "phase", "matching", "toy dr", "outer radius", "H0_chart"))
    if STATIC:
        G["G0.1_label"] = label
        G["G0.1_PASS"] = (G["G0.1_PASS"] and STATIC["control_PASS"]
                          and all(s in label for s in ("slice_time=static", "H0_phase", "harmonics", "control", "sha256")))
    G["G0.2_PASS"] = bool(np.array_equal(sl.r, r))
    G["G0.3_H0"] = H0
    G["G0.3_PASS"] = H0 != 0.0 and abs(sl.H0_chart - H0) <= 1e-12 * H0
    if STATIC:
        G["G0.3_H0_phase"] = H0_phase
        G["G0.3_PASS"] = G["G0.3_PASS"] and abs(H0_phase - H0) <= 1e-12 * H0
    G["G0.4_PASS"] = bool(np.all(phi == PHI0) and np.all(Pi == PI0))
    G["entry"] = {"last5_max": edge, "edge_tolerance": EDGE_TOL[dr], "over": over, "engine": entry}
    if any(nm in over for nm in ("D-D_s", "PD")):
        G["entry_reading"] = "FAIL (D or PD over the model-based tolerance)"
    elif over:
        G["entry_reading"] = "UNDETERMINED (S/PS over the scalar; no model basis for S/PS) - stop"
    else:
        G["entry_reading"] = "PASS"
    # ---- G0.5'
    dD, dS = grid.cell_derivative_fourth(D, parity=1), grid.cell_derivative_fourth(S, parity=1)
    PDr, PSr = PD - H0 * r * dD, PS - H0 * r * dS
    aD2 = r**2 * ((D - Ds) ** 2 + (PDr / OM2) ** 2)
    aS2 = r**2 * (S**2 + (PSr / OM3) ** 2)
    wD, wS = (r >= 20) & (r <= 60 - 5 * dr), (r >= 35) & (r <= 60 - 5 * dr)
    AD, mD, phD = fit(r[wD], aD2[wD], K2)
    AS, mS, phS = fit(r[wS], aS2[wS], K3)
    P_tail = 2 * math.pi * OM2 * K2 * AD
    cD, cS = float(np.max(np.abs(aD2[-5:] / AD - 1))), float(np.max(np.abs(aS2[-5:] / AS - 1)))
    sc = math.sqrt(AD) / A_D_PRED * 1.03
    G["G0.5p"] = {
        "D": {"A": AD, "m": mD, "phase": phD, "last5_dev": cD}, "S": {"A": AS, "m": mS, "phase": phS, "last5_dev": cS},
        "P_tail": P_tail, "P_tail_rel": P_tail / P_D_REF - 1,
        "a_PASS": abs(P_tail / P_D_REF - 1) <= 0.03,
        "b_PASS": mD <= 0.01 and mS <= 0.01,
        "c_PASS": cD <= 0.03 + mD and cS <= 0.03 + mS,
        "d": {"D_last5": edge["D-D_s"], "PD_last5": edge["PD"], "D_bound": ENV[dr][0] * sc, "PD_bound": ENV[dr][1] * sc,
              "D_pred": ENV[dr][0], "PD_pred": ENV[dr][1]},
    }
    G["G0.5p"]["d_PASS"] = edge["D-D_s"] <= ENV[dr][0] * sc and edge["PD"] <= ENV[dr][1] * sc
    # G0.5'e direction (AMENDMENT 3, 9.5/9.10; adopted Nov 24 02:34 UTC, runs declared after that only):
    # window mean of the rest-frame flux 4 pi r^2 (-P_rest f_r) > 0 for D and S. Ratio to P_tail recorded, no band.
    FD = float(np.mean((4 * math.pi * r**2 * (-PDr * dD))[wD]))
    FS = float(np.mean((4 * math.pi * r**2 * (-PSr * dS))[wS]))
    PS_tail = 2 * math.pi * OM3 * K3 * AS
    G["G0.5p"]["e"] = {"F_win_D": FD, "F_win_S": FS, "F_win_D/P_tail": FD / P_tail, "F_win_S/P_S_tail": FS / PS_tail}
    G["G0.5p"]["e_PASS"] = FD > 0 and FS > 0
    if STATIC:
        G["G0.5p"]["e"]["expected_ratio_static"] = E_EXPECT_STATIC
    # G0.5'f static-frame flux (W1 A3 4.3; gated for the corrected slice, recorded otherwise):
    # F_st = window mean of 4 pi r^2 (-P_st G_st), P_st = P - H0 r f_r, G_st = f_r - H0 r P_st (engine derivatives)
    GDst, GSst = dD - H0 * r * PDr, dS - H0 * r * PSr
    FstD = float(np.mean((4 * math.pi * r**2 * (-PDr * GDst))[wD]))
    FstS = float(np.mean((4 * math.pi * r**2 * (-PSr * GSst))[wS]))
    G["G0.5p"]["f"] = {"F_st_D/P_tail": FstD / P_tail, "F_st_S/P_S_tail": FstS / PS_tail, "band": F_BAND,
                       "gated": bool(STATIC)}
    G["G0.5p"]["f_PASS"] = all(F_BAND[0] <= q <= F_BAND[1] for q in (FstD / P_tail, FstS / PS_tail))
    # ---- recorded only (5.3)
    w = 4 * math.pi * r * r * dr
    G["recorded"] = {
        "P_S_tail": 2 * math.pi * OM3 * K3 * AS,
        "f_snap0": {str(rc): float(-(U / E_RING) * np.sum(((D**2 - Ds**2) * w)[r <= rc])) for rc in (25.0, 60.0)},
        "f_snap0_note": "E_ring = 1.1525e-4 (reference); r_cut 25 and 60",
        "peaks": {nm: {"max": float(np.max(np.abs(a))), "cell": int(np.argmax(np.abs(a)))}
                  for nm, a in (("S", S), ("PS", PS), ("D-D_s", D - Ds), ("PD", PD))},
        "old_G0.5_ratios": {nm: float(np.max(np.abs(a[-5:])) / np.max(np.abs(a)))
                            for nm, a in (("S", S), ("PS", PS), ("D-D_s", D - Ds), ("PD", PD))},
    }
    g5 = G["G0.5p"]
    G["input_checks_PASS"] = (G["G0.1_PASS"] and G["G0.2_PASS"] and G["G0.3_PASS"] and G["G0.4_PASS"]
                              and G["entry_reading"] == "PASS" and entry == "ACCEPTED"
                              and g5["a_PASS"] and g5["b_PASS"] and g5["c_PASS"] and g5["d_PASS"]
                              and g5["e_PASS"] and (g5["f_PASS"] or not STATIC))
    rep["grids"][key] = G

allp = all(g["input_checks_PASS"] for g in rep["grids"].values())
rep["outcome"] = ("INPUT CHECKS PASS on both grids (corrected slice); new bytes: the zero-step C9 set must rerun "
                  "on it before it is used (W1 A3 6)" if allp and STATIC else
                  "INPUT CHECKS PASS on both grids; #576 readings stand (same bytes, same engine; Amendment 2, 5.2)"
                  if allp and rep["unchanged_from_576"] else
                  "input checks pass but inputs/engine changed: full C9 set must rerun" if allp else
                  "REFUSED/UNDETERMINED at input checks: nothing after them is read")
(outdir / "report.json").write_text(json.dumps(rep, indent=1, default=float))
print(rep["outcome"])
