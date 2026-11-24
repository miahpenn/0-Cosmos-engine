"""f_ring Stage 0 timing smoke (declarations: results/research/w1_fring_smoke/PRERUN.md).

DESIGN-W1-f_ring-homeostasis with Amendments 1-4: 100 steps per arm, arms A (corrected W1 slice) and B (background
"none"), same grid, edge (open), eight-node start, dt = cfl dr (cfl 0.03). dr 0.025 gated; dr 0.05 pair recorded only.
Stops S0a, S0b (A2 2.2 for H; A4 1 for M), S0c, S0e; S0d recorded (A3 3). Full-run guard G_run evaluated and recorded.
Before any step: implementation check of the static-frame flux estimator (cross term G_st = f_R - v P_st) and fault
cases for every stop rule; if any fault check misbehaves, nothing is stepped.

Usage (repo root):
  python3 research/w1_fring_smoke/run.py native_static.npz OUTDIR            (driver)
  python3 research/w1_fring_smoke/run.py native_static.npz OUTDIR --arm A025  (one arm; the driver spawns these)
"""
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, os.getcwd())
import numpy as np  # noqa: E402

from engine import lattice_wall as _wall, scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad  # noqa: E402
from engine.adaptive_step import advance_with_radiation_admissibility_retries  # noqa: E402
from engine.cosmology_observables import effective_hubble  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402
from engine.run_production import assert_finite_state  # noqa: E402
from engine.stress_energy import assemble_total_stress_energy  # noqa: E402
from engine.worldtube import areal_radius, areal_radius_derivative  # noqa: E402
from research import b7_w2_origin_check_v1 as oc  # noqa: E402

npz_path, outdir = Path(sys.argv[1]), Path(sys.argv[2])
ARM = sys.argv[4] if len(sys.argv) > 4 and sys.argv[3] == "--arm" else None
outdir.mkdir(parents=True, exist_ok=True)
go, vac, _ = ad.vendor_modules()
kernel = V55ProductionKernel()
PHI0, PI0 = oc.clone.PHI0, oc.clone.PI0

# ---- declared in PRERUN.md before any step
CFL, NSTEPS, R_MAX = 0.03, 100, 60.0
EPS = 2.2e-16
EDGE_TOL = {0.025: 1.943e-7, 0.05: 1.949e-7}
GROWTH = math.exp(3 * 1.909379 * 0.075)          # 1.54 (A1 4.2 S0b), evaluated at t_100 of dr 0.025
GROWTH = round(GROWTH, 2)
H_MULT, M_MULT = 100, 1000                      # A2 2.2 (H); A4 1 (M)
S0C = 1.01
GRUN_FRAC = 1e-2
R_F = 20.0
STEPS_FULL = {0.025: 52650, 0.05: 26325}        # A1 2.1 (t_end 39.49 / dt)
CEILING_H = 24.0                                # S0e: serial arm-hours for the full run (PRERUN)
IMPL_TOL = 1e-3                                 # flux estimator vs input-check estimator, step 0 (PRERUN)
WIN = {"D": (20.0, 59.85), "S": (35.0, 59.60)}  # 'f windows at dr 0.025 (#581)
ARMS = {"A025": (0.025, 2400, True), "B025": (0.025, 2400, False),
        "A05": (0.05, 1200, True), "B05": (0.05, 1200, False)}
LABELS = json.loads(Path("results/research/w1_c9_static/input_checks/report.json").read_text())


def build(arm):
    dr, n, ring = ARMS[arm]
    if not ring:
        return kernel.initialize(resolution=n, r_max=R_MAX, initial_structure="none",
                                 outer_edge="open_edge", eight_node_areal=True)
    z = np.load(npz_path)
    key = "dr%g" % dr
    grid = go.SphericalCellGrid(n, R_MAX)
    r = np.asarray(grid.centers)
    H0 = float(vi.build_initial_data(grid, vac.flat_state, initial_structure="none", eight_node_areal=True).H0)
    S, PS, D, PD = (np.array(z[key + "_" + k], dtype=np.float64) for k in ("S", "PS", "D", "PD"))
    sl = vi.SuppliedSlice(ss.ScalarFields(S, PS, D, PD, np.full(n, PHI0), np.full(n, PI0)), r.copy(), H0,
                          EDGE_TOL[dr], LABELS["grids"][key]["G0.1_label"])
    return kernel.initialize(resolution=n, r_max=R_MAX, initial_structure="supplied_slice", supplied_slice=sl,
                             outer_edge="open_edge", eight_node_areal=True)


def residuals(state):
    grid = _wall.bind(state.grid, state.geometry)
    geom = state.geometry
    tot = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    raw = vac.constraints(grid, geom)
    H = raw["hamiltonian"] - 16.0 * math.pi * tot.rho
    M = raw["momentum"] - 8.0 * math.pi * tot.j
    return grid, geom, tot, raw, H, M


def finite(state, H, M):
    try:
        assert_finite_state(state)
    except Exception:  # noqa: BLE001
        return False
    return bool(np.all(np.isfinite(H)) and np.all(np.isfinite(M)))


def flux_density(grid, geom, f, P, Hub, cross=True):
    """Static-frame flux density 4 pi R^2 (-P_st G_st) (f_ring A2 1.7): v = H(t) R, f_R = f_r / R_r."""
    R = areal_radius(grid, geom)
    fR = grid.cell_derivative_fourth(f, parity=1) / areal_radius_derivative(grid, geom)
    v = Hub * R
    Pst = P - v * fR
    Gst = fR - v * Pst if cross else fR
    return R, 4 * math.pi * R**2 * (-Pst * Gst)


def measure(state):
    grid, geom, tot, raw, H, M = residuals(state)
    sc = state.scalars
    R = areal_radius(grid, geom)
    Hub = effective_hubble(geom, grid.volumes)
    out = {"t": float(state.t), "finite": finite(state, H, M),
           "H_max": float(np.max(np.abs(H))), "H_argR": float(R[np.argmax(np.abs(H))]),
           "M_max": float(np.max(np.abs(M))), "M_argR": float(R[np.argmax(np.abs(M))]),
           "S_max": float(np.max(np.abs(sc.S))), "H_t": float(Hub),
           "D_nonuni": float(np.max(np.abs((sc.D - sc.D[0])[R < 55.0])) / ss.sd_rest_value(PHI0))}
    for nm, f, P in (("D", sc.D, sc.PD), ("S", sc.S, sc.PS)):
        for cross in (True, False):
            Rr, g = flux_density(grid, geom, f, P, Hub, cross)
            out[f"flux{nm}_{'st' if cross else 'nocross'}_Rf"] = float(np.interp(R_F, Rr, g))
    return out, (grid, geom, tot, raw, H, M)


def step0_floors(state):
    m, (grid, geom, tot, raw, H, M) = measure(state)
    Aa = geom.Aa
    terms = {"Ricci": float(np.max(np.abs(vac.geometry_terms(grid, geom)["R"]))),
             "A2": float(np.max(np.abs(Aa**2 + 2 * (0.5 * Aa) ** 2))),
             "K2": float(np.max(np.abs(2 / 3 * geom.K**2))), "16pi_rho": float(np.max(np.abs(16 * math.pi * tot.rho)))}
    dr = float(grid.dr)
    rfH = max(m["H_max"], EPS * max(terms.values()))
    rfM = max(m["M_max"], EPS * float(np.max(np.abs(geom.K))) / dr)
    return m, {"H_terms_max": terms, "max_abs_K": float(np.max(np.abs(geom.K))), "r0_H": m["H_max"],
               "r0_M": m["M_max"], "r_floor_H": rfH, "r_floor_M": rfM}, (grid, geom, tot)


def s0b(r0, rf, r100, mult):
    return r100 > GROWTH * r0 + mult * rf


def s0c(s0, smax_series):
    return max(smax_series) > S0C * s0


def grun(res, C, rf, mult):
    return res > max(GRUN_FRAC * C, mult * rf)


def s0e(t_step):
    return sum(t_step[a] * STEPS_FULL[ARMS[a][0]] for a in t_step) / 3600.0


# ------------------------------------------------------------------ one arm
if ARM:
    st = build(ARM)
    dt = CFL * float(st.grid.dr)
    m0, fl, _ = step0_floors(st)
    rows = [m0]
    walls = []
    for k in range(NSTEPS):
        t0 = time.perf_counter()
        st = advance_with_radiation_admissibility_retries(kernel, st, dt).state
        walls.append(time.perf_counter() - t0)
        m, _ = measure(st)
        rows.append(m)
        if not m["finite"]:
            break
    (outdir / f"arm_{ARM}.json").write_text(json.dumps({"arm": ARM, "dt": dt, "floors": fl, "rows": rows,
                                                        "step_wall": walls}, indent=0))
    sys.exit(0)

# ------------------------------------------------------------------ driver: prechecks before any step
pre = {}
s = {a: build(a) for a in ("A025", "B025")}
fl = {}
for a in s:
    m0, fl[a], ctx = step0_floors(s[a])
    fl[a]["m0"] = m0
gA, geomA, totA = step0_floors(s["A025"])[2]
_, geomB, totB = step0_floors(s["B025"])[2]
drho = np.asarray(totA.rho) - np.asarray(totB.rho)
C_H, C_M = 16 * math.pi * float(np.max(drho)), 8 * math.pi * float(np.max(drho))
pre["floors"], pre["C_H"], pre["C_M"] = fl, C_H, C_M

# Implementation check: evolved-state static flux at step 0 vs the input-check estimator (H0 r, engine operator)
z = np.load(npz_path)
r = np.asarray(gA.centers)
H0 = float(vi.build_initial_data(gA, vac.flat_state, initial_structure="none", eight_node_areal=True).H0)
HubA = effective_hubble(geomA, gA.volumes)
impl = {"H_t0": float(HubA), "H0": H0}
for nm, f, P in (("D", s["A025"].scalars.D, s["A025"].scalars.PD), ("S", s["A025"].scalars.S, s["A025"].scalars.PS)):
    F, Pz = np.asarray(z["dr0.025_" + nm]), np.asarray(z["dr0.025_P" + nm])
    dF = gA.cell_derivative_fourth(F, parity=1)
    Pst = Pz - H0 * r * dF
    g_in = 4 * math.pi * r**2 * (-Pst * (dF - H0 * r * Pst))
    w = (r > WIN[nm][0]) & (r < WIN[nm][1])
    top = float(np.max(np.abs(g_in[w])))
    var = {}
    for tag, kw in (("static", dict(Hub=HubA, cross=True)), ("fault_no_cross", dict(Hub=HubA, cross=False)),
                    ("fault_H_zero", dict(Hub=0.0, cross=True))):
        Rr, g = flux_density(gA, geomA, f, P, kw["Hub"], kw["cross"])
        var[tag] = float(np.max(np.abs(np.interp(r[w], Rr, g) - g_in[w])) / top)
    impl[nm] = {**var, "static_PASS": var["static"] <= IMPL_TOL,
                "faults_fail": var["fault_no_cross"] > IMPL_TOL and var["fault_H_zero"] > IMPL_TOL}
pre["impl_check"] = impl

# Fault cases for every stop rule (must trip on the fault and pass the clean case)
fc = {}
st_nan = build("A025")
st_nan.scalars.PD[1200] = np.nan
_, (_, _, _, _, Hn, Mn) = measure(st_nan)
fc["S0a"] = {"clean_finite": fl["A025"]["m0"]["finite"], "nan_cell_detected": not finite(st_nan, Hn, Mn)}
for a in ("A025", "B025"):
    f_ = fl[a]
    seq_ok_H, seq_ok_M = 1.5 * f_["r0_H"], 1.5 * f_["r0_M"]
    inst_H = f_["r_floor_H"] * 1.06 ** NSTEPS + f_["r0_H"]     # > 1.05 per step (A2 2.2)
    inst_M = f_["r_floor_M"] * 1.08 ** NSTEPS + f_["r0_M"]     # > 1.072 per step (A4 1)
    fc[f"S0b_{a}"] = {"healthy_H_pass": not s0b(f_["r0_H"], f_["r_floor_H"], seq_ok_H, H_MULT),
                      "unstable_H_trip": s0b(f_["r0_H"], f_["r_floor_H"], inst_H, H_MULT),
                      "healthy_M_pass": not s0b(f_["r0_M"], f_["r_floor_M"], seq_ok_M, M_MULT),
                      "unstable_M_trip": s0b(f_["r0_M"], f_["r_floor_M"], inst_M, M_MULT)}
# the residual measure itself responds to a grid-scale metric defect (checkerboard 1e-10 in K)
st_k = build("A025")
st_k.geometry.K = st_k.geometry.K + 1e-10 * (-1.0) ** np.arange(st_k.grid.n)
mk, _ = measure(st_k)
fc["S0b_measure"] = {"M_with_K_checkerboard": mk["M_max"], "trips_vs_floor":
                     s0b(fl["A025"]["r0_M"], fl["A025"]["r_floor_M"], mk["M_max"], M_MULT)}
s0 = fl["A025"]["m0"]["S_max"]
fc["S0c"] = {"x1.005_pass": not s0c(s0, [s0, 1.005 * s0]), "x1.02_trip": s0c(s0, [s0, 1.02 * s0])}
fc["S0e"] = {"slow_trip": s0e({a: 1.0 for a in ARMS}) > CEILING_H, "fast_pass": s0e({a: 0.1 for a in ARMS}) <= CEILING_H}
fc["G_run"] = {"H_2e-2_trip": grun(2e-2 * C_H, C_H, fl["A025"]["r_floor_H"], H_MULT),
               "H_5e-3_pass": not grun(5e-3 * C_H, C_H, fl["A025"]["r_floor_H"], H_MULT),
               "M_2e-2_trip": grun(2e-2 * C_M, C_M, fl["A025"]["r_floor_M"], M_MULT),
               "M_5e-3_pass": not grun(5e-3 * C_M, C_M, fl["A025"]["r_floor_M"], M_MULT)}
pre["fault_checks"] = fc
flat = [v for d in fc.values() for k, v in d.items() if isinstance(v, bool)]
pre["prechecks_PASS"] = (all(flat) and all(impl[nm]["static_PASS"] and impl[nm]["faults_fail"] for nm in ("D", "S")))
(outdir / "prechecks.json").write_text(json.dumps(pre, indent=1, default=float))
print("prechecks", pre["prechecks_PASS"], json.dumps(impl), flush=True)
if not pre["prechecks_PASS"]:
    sys.exit("STOP: a stop rule or the flux estimator does not discriminate; nothing stepped")
del s, st_nan, st_k

# ------------------------------------------------------------------ the smoke: four arms in parallel
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
procs = [subprocess.Popen([sys.executable, __file__, str(npz_path), str(outdir), "--arm", a], env=env) for a in ARMS]
codes = [p.wait() for p in procs]
res = {a: json.loads((outdir / f"arm_{a}.json").read_text()) for a in ARMS}
rep = {"prechecks": pre, "arms": {}, "stops": {}}
t_step = {}
for a, d in res.items():
    rows, f_ = d["rows"], d["floors"]
    t_step[a] = float(np.median(d["step_wall"]))
    last = rows[-1]
    arm = {"dt": d["dt"], "steps": len(rows) - 1, "t_end": last["t"], "step_wall_median": t_step[a],
           "step_wall_mean": float(np.mean(d["step_wall"])), "floors": f_,
           "S0a_nonfinite": not all(x["finite"] for x in rows),
           "S0b_H": {"r0": f_["r0_H"], "r100": last["H_max"], "limit": GROWTH * f_["r0_H"] + H_MULT * f_["r_floor_H"],
                     "trip": s0b(f_["r0_H"], f_["r_floor_H"], last["H_max"], H_MULT), "argR": last["H_argR"]},
           "S0b_M": {"r0": f_["r0_M"], "r100": last["M_max"], "limit": GROWTH * f_["r0_M"] + M_MULT * f_["r_floor_M"],
                     "trip": s0b(f_["r0_M"], f_["r_floor_M"], last["M_max"], M_MULT), "argR": last["M_argR"]},
           "H_series": [x["H_max"] for x in rows], "M_series": [x["M_max"] for x in rows],
           "H_t": [rows[0]["H_t"], last["H_t"]]}
    if ARMS[a][2]:
        arm["S0c"] = {"S0": rows[0]["S_max"], "max": max(x["S_max"] for x in rows), "trip": s0c(rows[0]["S_max"],
                                                                                             [x["S_max"] for x in rows])}
    else:
        arm["S0d_recorded"] = {"step1": rows[1]["D_nonuni"], "step10": rows[10]["D_nonuni"],
                               "step100": rows[-1]["D_nonuni"]} if len(rows) > 100 else "incomplete"
    if a == "A025":
        arm["G_run"] = {"H_trip_any": any(grun(x["H_max"], C_H, f_["r_floor_H"], H_MULT) for x in rows),
                        "M_trip_any": any(grun(x["M_max"], C_M, f_["r_floor_M"], M_MULT) for x in rows),
                        "H_max_over_C_H": max(x["H_max"] for x in rows) / C_H,
                        "M_max_over_C_M": max(x["M_max"] for x in rows) / C_M}
    rep["arms"][a] = arm
# static-frame flux at R_f = 20, ring minus control, recorded per step (dr 0.025)
for pair, (A_, B_) in {"dr0.025": ("A025", "B025"), "dr0.05": ("A05", "B05")}.items():
    ra, rb = res[A_]["rows"], res[B_]["rows"]
    rep["flux_Rf20_" + pair] = {k: [x[k] - y[k] for x, y in zip(ra, rb)]
                                for k in ("fluxD_st_Rf", "fluxD_nocross_Rf", "fluxS_st_Rf", "fluxS_nocross_Rf")}
cost = s0e({a: t_step[a] for a in ("A025", "B025")})
cost_all = s0e(t_step)
rep["S0e"] = {"step_wall_median": t_step, "full_run_serial_hours_dr025_pair": cost,
              "full_run_serial_hours_all_four": cost_all, "ceiling_hours": CEILING_H, "trip": cost > CEILING_H}
g = rep["arms"]
trips = [f"S0a {a}" for a in ("A025", "B025") if g[a]["S0a_nonfinite"]]
trips += [f"S0b {c} {a}" for a in ("A025", "B025") for c in ("H", "M") if g[a][f"S0b_{c}"]["trip"]]
trips += ["S0c A025"] if g["A025"]["S0c"]["trip"] else []
trips += ["S0e"] if rep["S0e"]["trip"] else []
rec = [f"S0a {a}" for a in ("A05", "B05") if g[a]["S0a_nonfinite"]]
rec += [f"S0b {c} {a}" for a in ("A05", "B05") for c in ("H", "M") if g[a][f"S0b_{c}"]["trip"]]
rec += ["S0c A05"] if g["A05"]["S0c"]["trip"] else []
rep["exit_codes"] = codes
rep["stops_tripped"], rep["record_only_dr005_would_trip"] = trips, rec
rep["outcome"] = ("STAGE 0 PASS: no stop holds (dr 0.025 pair); full run needs its own go" if not trips and all(c == 0 for c in codes)
                  else "STOP: " + ", ".join(trips or ["arm exit code"]))
(outdir / "report.json").write_text(json.dumps(rep, indent=1, default=float))
print(rep["outcome"], json.dumps(rep["S0e"]))
