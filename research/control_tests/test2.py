"""Test 2 runner: weak S pulse in flat vacuum (design 1BK-Xj3r... section 3; A2 P4, P9-P11; A3; A5 arm list; A21).

Research only, engine venv. Code only until Jeremiah's go for stage 2. Nothing here reads a gate: the arms write
JSON/npz readings, and the 60-digit readback (mp_readback, input checks) and score2 (stage-2 gates) score them
afterwards, against the predeclared gates.

  python -m research.control_tests.test2 inputs OUT                  t = 0 input-check slices (zero steps)
  python -m research.control_tests.test2 arm NAME OUT                one stage-2 arm (STAGE2)
  python -m research.control_tests.test2 arm NAME OUT --stage3       a stage-3 arm (needs its own go)
  python -m research.control_tests.test2 cdot_floor STAGE1_ARMS OUT.json
        reanalysis of the stored stage-1 1a saves (no step): Cdot at each save, the floor for GS.2 (A2 P11)

Initial data (design 3.2, A2 P4): S = amp exp(-r^2/sigma^2), sigma = 1, P_S = D = P_D = phi = Pi = 0; a = b = 1,
Aa = K = Lambda = 0 (Lambda_metric = 0 at a = b = 1); alpha = 1 as the first solve's input; dust 1e-150 as in Test 1
(ctl.make); V0 = SD_QUARTIC_LAMBDA = 0 (ctl.make). X from a Newton solve of the engine's own discrete
H(X) - 16 pi rho(X) = 0, H from vendor constraints and rho from assemble_total_stress_energy (ctl.constraints).

Newton normalisation (confirmed by Lumen^2, A22 item 1; exact for the flat exterior psi = 1 + m/2r): constant X is a null mode of the discrete H at
rho = 0, so the last row (an edge cell, a reading only) is replaced by r (psi - 1) = const across the last two cells,
psi = X^(-1/2): the exterior conformal factor 1 + m/(2r). Stop: the step no longer halves the interior residual
(the rounding floor), at most MAX_NEWTON iterations. The A2 P4 stop rule, in its B5 form (A22 item 2:
|H| <= 10 eps x Sum|s_k| per cell), is scored in mp_readback, not here.

Amplitudes: T2-4 at 3e-2 (A21, Jeremiah's "Yes to both", 08:07:50 UTC); every other arm as A5.
"""
import json
import sys
from pathlib import Path

import numpy as np

from engine import scalar_system as ss
from research.control_tests import ctl, observer

DT0, DT1 = 3.75e-4, 1.875e-4
T_END = 2.1
SIGMA = 1.0
EDGE = 2
MAX_NEWTON = 25
H_FD = 2.0 ** -20

# name: (amp, dr, dt, r_max, per-step Cdot, R-stage)
ARMS = {
    "T2-3": (1e-2, 0.05, DT0, 8.0, True, True),
    "T2-4": (3e-2, 0.025, DT0, 8.0, False, True),
    "T2-6": (5e-3, 0.05, DT0, 8.0, False, False),
    "T2-E": (1e-2, 0.05, DT0, 12.0, False, False),
    "T2-1": (1e-5, 0.05, DT0, 8.0, False, False),
    "T2-2": (1e-5, 0.025, DT0, 8.0, False, False),
    "T2-5": (1e-2, 0.05, DT1, 8.0, True, True),
}
STAGE2 = ("T2-3", "T2-4", "T2-6", "T2-E")


def gaussian(r, amp):
    return amp * np.exp(-(r * r) / (SIGMA * SIGMA))


def newton_X(state):
    """X solving the engine's discrete H = 16 pi rho on every row but the last (the normalisation row)."""
    g, r = state.geometry, np.asarray(state.grid.centers)
    n = len(r)

    def F(X):
        g.X[:] = X
        H = ctl.constraints(state)[0]
        out = np.array(H, dtype=float)
        psi = X ** -0.5
        out[-1] = r[-1] * (psi[-1] - 1.0) - r[-2] * (psi[-2] - 1.0)
        return out

    def size(f):
        return float(np.max(np.abs(f[:n - EDGE])))

    X = np.ones(n)
    f = F(X)
    best, log = (size(f), X.copy()), [{"iter": 0, "max_abs_H_non_edge": size(f), "norm_row": float(f[-1])}]
    for it in range(1, MAX_NEWTON + 1):
        J = np.empty((n, n))
        for j in range(n):
            Xp, Xm = X.copy(), X.copy()
            Xp[j] += H_FD; Xm[j] -= H_FD
            J[:, j] = (F(Xp) - F(Xm)) / (2 * H_FD)
        X = X + np.linalg.solve(J, -F(X))
        if not np.all(X > 0):
            raise RuntimeError("Newton left X > 0")
        f = F(X)
        log.append({"iter": it, "max_abs_H_non_edge": size(f), "norm_row": float(f[-1])})
        if size(f) < best[0]:
            best = (size(f), X.copy())
        if size(f) > 0.5 * log[-2]["max_abs_H_non_edge"]:
            break
    g.X[:] = best[1]
    return best[1], log


def build_t2(dr, amp, r_max):
    kernel, state = ctl.make(dr, r_max)
    n = state.grid.n
    S = gaussian(np.asarray(state.grid.centers, dtype=float), amp)
    state.scalars = ss.ScalarFields(S, *(np.zeros(n) for _ in range(5)))
    _, log = newton_X(state)
    return kernel, state, {"amp": amp, "sigma": SIGMA, "newton": log}


def input_name(amp, dr, r_max):
    return f"T2_amp{amp:g}_dr{dr}_r{r_max:g}"


def inputs(out):
    """Zero-step t = 0 slices for the A2 section 2 Test 2 input checks, one per distinct (amp, dr, r_max)."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    meta = {}
    for amp, dr, rm in sorted({(a[0], a[1], a[3]) for a in ARMS.values()}):
        name = input_name(amp, dr, rm)
        k, s, p = build_t2(dr, amp, rm)
        g0, kdot, R, _ = ctl.semi(k, s)
        H, M = ctl.constraints(s)
        tot = observer.matter_totals(s.grid, s.geometry, s.scalars, s.matter)
        extra = {"H": H, "M": M, **{f"m_{x}": v for x, v in tot.items()}}
        np.savez_compressed(out / f"{name}.npz", **ctl.save_arrays(s, g0, kdot, R, extra))
        meta[name] = {"dr": s.grid.dr, "n": s.grid.n, "r_max": s.grid.r_max, **p}
    (out / "inputs.json").write_text(json.dumps(meta, indent=1))


def matter_arrays(st):
    from engine import v55_matter as vm
    tot = vm.total_matter_projection(st.grid, st.geometry, st.scalars, st.matter)
    return {"m_rho": tot["rho"], "m_pr": tot["pr"], "m_pt": tot["pt"], "m_j_scalar": tot["scalar"]["j"],
            "m_j_fluid": tot["fluid"]["j"], "m_rho_scalar": tot["scalar"]["rho"]}


def arm(name, out, stage3=False):
    if name not in STAGE2 and not stage3:
        raise SystemExit(f"{name} is a stage-3 arm (A5 section 4): pass --stage3 only on its own go")
    amp, dr, dt, r_max, with_cdot, rstage_on = ARMS[name]
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    kernel, state, p = build_t2(dr, amp, r_max)
    rec_kernel = observer.Recorder()
    rstage, checked, nsave = [], [], []

    def on_save(st, g0, kdot, R):
        k = len(nsave); nsave.append(k)
        cH, cM = ctl.cdot(st, R, dt)
        arrays = {"cH": cH, "cM": cM, **{x: np.asarray(v, dtype=float) for x, v in matter_arrays(st).items()}}
        np.savez_compressed(out / f"{name}_save{k}.npz", **arrays)
        row = {"cH": ctl.fixed_points(cH, dr), "cM": ctl.fixed_points(cM, dr)}
        if rstage_on:
            if not checked:
                observer.check_bitwise(kernel, rec_kernel, st, dt)
                checked.append(True)
            new, log = rec_kernel.recorded_step(st, dt)
            reading, ra = observer.r_stage(st.grid, st, new, log, dr)
            np.savez_compressed(out / f"{name}_rstage{k}.npz", **ra)
            rstage.append({"save": k, "t": st.t, "terms": reading})
            (out / f"{name}_rstage.json").write_text(json.dumps(rstage))
            row["rstage_file"] = f"{name}_rstage{k}.npz"
        return row

    def on_step(n, st, per_save):
        # H, M at the steps next to each save, for d_t H at the save (G2.4, the #602 rule)
        if n % per_save in (1, per_save - 1, per_save - 2):
            H, M = ctl.constraints(st)
            return {"step": n, "t": n * dt, "H": ctl.fixed_points(H, dr), "M": ctl.fixed_points(M, dr)}
        return None

    _, rec = ctl.run(name, kernel, state, dt, out, with_cdot=with_cdot, t_end=T_END, on_save=on_save,
                     on_step=on_step)
    rec["params"] = p
    rec["amp"] = amp
    (out / f"{name}.json").write_text(json.dumps(rec))


def cdot_floor(arms_dir, out):
    """Cdot (A2 P11, h = dt) on each stored stage-1 1a save, at the fixed points. Zero steps."""
    arms_dir = Path(arms_dir)
    res = {}
    for js in sorted(arms_dir.glob("1a_*.json")):
        if js.name.endswith("_rstage.json"):
            continue
        rec = json.loads(js.read_text())
        dr, dt, rm = rec["dr"], rec["dt"], rec["r_max"]
        rows = []
        with np.load(arms_dir / f"{js.stem}.npz") as z:
            for k, sv in enumerate(rec["saves"]):
                _, state = ctl.make(dr, rm)
                g = state.geometry
                for f in ctl.VAC.ALL_FIELDS:
                    getattr(g, f)[:] = z[f"t{k}_g_{f}"]
                state.scalars = ss.ScalarFields(*(np.array(z[f"t{k}_s_{f}"]) for f in ctl.FIELDS_S))
                R = {f: np.array(z[f"t{k}_R_{f}"]) for f in ctl.FIELDS_G}
                R.update({f"s_{f}": np.array(z[f"t{k}_R_s_{f}"]) for f in ctl.FIELDS_S})
                cH, cM = ctl.cdot(state, R, dt)
                rows.append({"t": sv["t"], "cH": ctl.fixed_points(cH, dr), "cM": ctl.fixed_points(cM, dr)})
        res[js.stem] = {"dr": dr, "dt": dt, "r_max": rm, "saves": rows}
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "inputs":
        inputs(sys.argv[2])
    elif sys.argv[1] == "arm":
        arm(sys.argv[2], sys.argv[3], stage3="--stage3" in sys.argv[4:])
    elif sys.argv[1] == "cdot_floor":
        cdot_floor(sys.argv[2], sys.argv[3])
    else:
        raise SystemExit(__doc__)
