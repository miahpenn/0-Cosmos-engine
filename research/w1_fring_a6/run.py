"""f_ring A6 runs with the A7 in-run write ledger (Jeremiah's go, Nov 24 2026 07:11 UTC). Choices: PRERUN.md.

Arms at dr 0.05, open edge, eight-node start, KO off, engine unchanged (origin_trace checks it against #584):
  A1 ring, lambda 0.03 (dt 1.5e-3, 2194 steps)   A2 ring, lambda 0.015 (dt 7.5e-4, 4388 steps)
  B1 control, lambda 0.03                        B2 control, lambda 0.015
each to t 3.291, with state saves at t 0.015, 0.8235, 1.6455, 3.291.
Usage (repo root): python -m research.w1_fring_a6.run OUT --arm A1   (one arm), then python -m research.w1_fring_a6.run OUT
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from engine.native_geometry.evolution import native_cmc_explicit_metric_rhs
from engine.production_kernel import NativeSphericalCellGrid
from engine.v55_pirk_adapter import vendor_modules
from research.w1_fring_a6 import gates
from research.w1_fring_smoke import areal_trace, metric_balance, origin_trace

WIDE = metric_balance.WIDE
DR, END = 0.05, 3.291
SAVES = (0.015, 0.8235, 1.6455, 3.291)
ARMS = {"A1": ("A", 1), "A2": ("A", 2), "B1": ("B", 1), "B2": ("B", 2)}
TRACE = Path("results/research/w1_fring_origin_trace")
X_TERMS = ("tangent", "increment_arithmetic", "write")
AB_TERMS = ("tangent", "integration", "increment_arithmetic", "write")
CELLS = slice(0, 6)


def hashes():
    paths = [Path(__file__), Path(gates.__file__), Path(origin_trace.__file__), Path(areal_trace.__file__),
             Path(metric_balance.__file__), TRACE / "A-0.05.json", TRACE / "B-0.05.json"]
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def coords(log_fields):
    a, b = log_fields
    return np.stack((a + 2*b, (b - a)/3))  # D = ln a + 2 ln b, T = (ln b - ln a)/3


def numbers(values):
    return metric_balance.numbers(values)


def save_state(path, state, vacuum, extra):
    arrays = {f"geometry_{k}": np.asarray(getattr(state.geometry, k)) for k in vacuum.ALL_FIELDS}
    arrays.update({f"scalars_{k}": np.asarray(getattr(state.scalars, k)) for k in ("S", "PS", "D", "PD", "phi", "Pi")})
    for species in ("dark_matter", "baryons", "radiation"):
        s = getattr(state.matter, species)
        arrays.update({f"matter_{species}_{k}": np.asarray(getattr(s, k)) for k in ("rest", "energy_t", "momentum_r")})
    if state.matter.charged_dust is not None or state.electromagnetic is not None:
        raise ValueError("expected the uncharged W1 state")
    arrays["clock"] = np.array([state.t, state.tau, state.e_folds])
    np.savez(path, **arrays, **extra)


def run_arm(name, out):
    arm, divisor = ARMS[name]
    source = hashes()
    _, vacuum, _ = vendor_modules()
    n = round(60/DR)
    native = NativeSphericalCellGrid(n, 60.)
    dt = .03*DR/divisor
    steps = round(END/dt)
    save_steps = {round(t/dt): t for t in SAVES}
    trace = json.loads((TRACE / f"{arm}-0.05.json").read_text())["rows"] if divisor == 1 else None
    x_tot = {k: np.zeros(n, dtype=WIDE) for k in X_TERMS}
    ab_tot = {k: np.zeros((2, n), dtype=WIDE) for k in AB_TERMS}
    stages, carry = {}, {"previous": None, "initial": None, "count": 0}
    per_step, saves, walls = [], {}, [time.perf_counter()]

    def stage(stage_name, grid, geometry):
        if stage_name == "after_initial_lapse_refresh":
            stages.clear()
        stages[stage_name] = geometry
        if stage_name != "corrected_after_boundary":
            return
        start, primary = stages["after_initial_lapse_refresh"], stages["primary_predictor_after_lapse"]
        predictor, stored = stages["explicit_predictor"], stages["corrected_before_boundary"]
        fields = np.stack([start.a, start.b, start.X])
        if carry["previous"] is None:
            carry["initial"] = fields.astype(WIDE)
        else:
            np.testing.assert_array_equal(fields, carry["previous"])
        end = np.stack([stored.a, stored.b, stored.X])
        np.testing.assert_array_equal(end, np.stack([geometry.a, geometry.b, geometry.X]))
        r0, r1 = native_cmc_explicit_metric_rhs(native, start), native_cmc_explicit_metric_rhs(native, primary)
        keys = ("a", "b", "X")
        np.testing.assert_array_equal([getattr(start, k) + dt*r0[k] for k in keys],
                                      [getattr(predictor, k) for k in keys])
        increment = np.stack([(dt/2)*(r0[k] + r1[k]) for k in keys])
        np.testing.assert_array_equal(fields + increment, end)

        # D, T: conformal_ledger arithmetic (#587/#588).
        wide_start, wide_end = fields[:2].astype(WIDE), end[:2].astype(WIDE)
        primary_ab = np.stack([primary.a, primary.b]).astype(WIDE)
        f0 = start.alpha.astype(WIDE)*start.Aa.astype(WIDE)
        f1 = primary.alpha.astype(WIDE)*primary.Aa.astype(WIDE)
        rates = np.array([[-2], [1]], dtype=WIDE)
        tangent = (WIDE(dt)/2)*rates*(f0 + f1)
        ideal = (WIDE(dt)/2)*rates*(wide_start*f0 + primary_ab*f1)
        ideal_log = np.log1p(ideal/wide_start)
        inc = increment[:2].astype(WIDE)
        write = wide_end - wide_start - inc
        parts = {"tangent": tangent, "integration": ideal_log - tangent,
                 "increment_arithmetic": np.log1p((inc - ideal)/(wide_start + ideal)),
                 "write": np.log1p(write/(wide_start + inc))}
        for k, v in parts.items():
            ab_tot[k] += coords(v)
        # X: x_update_ledger arithmetic, rate kept whole (no factor split).
        g0 = start.X.astype(WIDE)*start.alpha.astype(WIDE)*start.K.astype(WIDE)
        g1 = primary.X.astype(WIDE)*primary.alpha.astype(WIDE)*primary.K.astype(WIDE)
        x_rate = (WIDE(dt)/2)*(g0 + g1)/3
        x_tot["tangent"] += x_rate
        x_tot["increment_arithmetic"] += increment[2].astype(WIDE) - x_rate
        x_tot["write"] += end[2].astype(WIDE) - fields[2].astype(WIDE) - increment[2].astype(WIDE)
        carry["previous"] = end
        carry["count"] += 1

    def accepted(state, vac, row):
        i = len(per_step)
        if trace is not None and i < len(trace) and row != trace[i]:
            raise ValueError(f"accepted slice disagrees with the pinned origin trace at step {i}")
        H, M = row["cells"]["H"][:6], row["cells"]["M"][:6]
        finite = bool(np.all(np.isfinite(H)) and np.all(np.isfinite(M)))
        per_step.append({"t": state.t, "H": H, "M": M, "finite": finite})
        walls.append(time.perf_counter())
        if i not in save_steps:
            return
        if carry["count"] != i:
            raise ValueError("ledger and accepted steps out of step")
        t = save_steps[i]
        obs = areal_trace.observe(state, vac, row)
        endpoint = np.stack(metric_balance.metric_arrays(state.geometry, WIDE))
        x_identity = float(np.max(abs(endpoint[2] - carry["initial"][2] - sum(x_tot.values()))))
        ab_identity = float(np.max(abs(coords(np.log1p((endpoint[:2] - carry["initial"][:2])/carry["initial"][:2]))
                                       - sum(ab_tot.values()))))
        if x_identity > 1e-12 or ab_identity > 1e-12:
            raise ValueError(f"ledger does not telescope at t {t}: {x_identity} {ab_identity}")
        cells = row["cells"]
        saves[f"{t:g}"] = {
            "step": i, "t": state.t, "observe": obs["cells"],
            "row": cells, "ledger_identity_gap": {"X": x_identity, "DT": ab_identity},
            "B_minus_one": numbers((state.geometry.X**2/state.geometry.a - 1)[CELLS]),
            "U_minus_one": numbers((np.sqrt(state.geometry.b)/state.geometry.X - 1)[CELLS]),
        }
        save_state(out / f"{name}-t{t:g}.npz", state, vac, {
            "endpoint": endpoint, "initial": carry["initial"],
            **{f"X_{k}": v for k, v in x_tot.items()}, **{f"DT_{k}": v for k, v in ab_tot.items()}})
        (out / f"{name}.json").write_text(json.dumps(
            {"arm": name, "ring": arm == "A", "lambda": .03/divisor, "dt": dt, "steps": steps, "sha256": source,
             "saves": saves, "per_step": per_step, "step_wall": float(np.mean(np.diff(walls)))},
            allow_nan=False))

    def quiet_record(self, stage_name, grid, geometry, scalars, matter):
        self.stage_observer(stage_name, grid, geometry.copy())
    original = origin_trace.Observer.record
    origin_trace.Observer.record = quiet_record
    try:
        result = origin_trace.run(DR, arm, divisor, end_time=END, accepted_observer=accepted,
                                  stage_observer=stage, stage_step=None)
    finally:
        origin_trace.Observer.record = original
    if result["attempts"] != steps or len(saves) != len(SAVES) or hashes() != source:
        raise ValueError("unexpected attempts, saves or changed sources")
    print(name, "done", steps, per_step[-1]["H"][0], flush=True)


def pair_ledger(grid, left, right):
    """Curvature responses of one dt pair at one save on frozen common backgrounds (x_endpoint, conformal_ledger)."""
    lv, rv = left["endpoint"], right["endpoint"]
    dx = {k: left[f"X_{k}"] - right[f"X_{k}"] for k in X_TERMS}
    dab = {k: left[f"DT_{k}"] - right[f"DT_{k}"] for k in AB_TERMS}
    np.testing.assert_array_equal(left["initial"], right["initial"])
    out = {}
    for width in areal_trace.WIDTHS:
        curv = lambda arrays: metric_balance.curvature(grid, arrays, width)  # noqa: E731
        full = curv(tuple(lv)) - curv(tuple(rv))
        a, b = (lv[0] + rv[0])/2, (lv[1] + rv[1])/2
        xbase = curv((a, b, rv[2]))
        x_only = curv((a, b, lv[2])) - xbase
        x_parts = {k: curv((a, b, rv[2] + v)) - xbase for k, v in dx.items()}
        x_closure = x_only - sum(x_parts.values())
        xm = (lv[2] + rv[2])/2
        abase = curv((rv[0], rv[1], xm))
        ab_total = curv((lv[0], lv[1], xm)) - abase

        def response(det, tf):
            return curv((rv[0]*np.exp(det/3 - 2*tf), rv[1]*np.exp(det/3 + tf), xm)) - abase
        ab_parts = {f"{label}_{k}": response(v[0], 0) if i == 0 else response(0, v[1])
                    for k, v in dab.items() for i, label in enumerate(("D", "T"))}
        ab_closure = ab_total - sum(ab_parts.values())
        interaction = full - x_only - ab_total
        write = x_parts["write"] + ab_parts["D_write"] + ab_parts["T_write"]
        increment = (x_parts["increment_arithmetic"] + ab_parts["D_increment_arithmetic"]
                     + ab_parts["T_increment_arithmetic"])
        out[width] = {"full": full, "X_only": x_only, "X_parts": x_parts, "X_closure": x_closure,
                      "DT_total": ab_total, "DT_parts": ab_parts, "DT_closure": ab_closure,
                      "interaction": interaction, "write": write, "increment_arithmetic": increment,
                      "closure_magnitude": abs(x_closure) + abs(ab_closure) + abs(interaction)}
    return out


def combine(out):
    source = hashes()
    arms = {}
    for name in ARMS:
        arms[name] = json.loads((out / f"{name}.json").read_text())
        if arms[name]["sha256"] != source:
            raise ValueError("sources changed since the run")
    report = {"source": origin_trace.SOURCE, "vendor": origin_trace.VENDOR, "input_sha256": origin_trace.DIGEST,
              "sha256": source, "lambda": {k: v["lambda"] for k, v in arms.items()},
              "step_wall_s": {k: v["step_wall"] for k, v in arms.items()}}

    # First check, before anything else is read.
    first = gates.first_check(arms["A1"]["saves"]["0.015"]["observe"]["H_metric_9"][0])
    report["first_check"] = first
    if not first["pass"]:
        print(json.dumps(report, indent=1, allow_nan=False))
        raise SystemExit("first check failed: nothing else is read")

    per_arm = {}
    for name, data in arms.items():
        rows = {}
        for t, s in data["saves"].items():
            o, c = s["observe"], s["row"]
            H_terms = [o["R_production"][:6], o["H_source"][:6], o["H_kinematic"][:6]]
            M_terms = [c["M_Aa"][:6], c["M_K"][:6], c["M_other_geometry"][:6], c["M_source"][:6]]
            rows[t] = {
                "H_production": o["H_production"][:6], "H_direct_9": o["H_metric_9"][:6],
                "H_direct_7": o["H_metric_7"][:6], "connection_term": o["R_connection"][:6],
                "H_adjusted_9": o["H_adjusted_9"][:6], "production_minus_adjusted_9": o["R_gap_9"][:6],
                "B_minus_one": s["B_minus_one"], "U_minus_one": s["U_minus_one"],
                "H_terms": {"Ricci": H_terms[0], "16pi_rho": H_terms[1], "kinematic": H_terms[2]},
                "M": c["M"][:6], "M_terms": {"M_Aa": M_terms[0], "M_K": M_terms[1],
                                             "M_other_geometry": M_terms[2], "M_source": M_terms[3]},
                "3i_H": gates.relative_size(o["H_production"][:6], H_terms),
                "3i_M": gates.relative_size(c["M"][:6], M_terms),
                "ledger_identity_gap": s["ledger_identity_gap"],
            }
        steps = [p for p in data["per_step"] if p["t"] >= 0.015 - 1e-12]
        hmax = [max(abs(h) for h in p["H"]) for p in steps]
        per_arm[name] = {"saves": rows, "finite": all(p["finite"] for p in data["per_step"]),
                         "growth": gates.growth([p["t"] for p in steps], hmax),
                         "3i_all_saves": all(r["3i_H"]["pass"] and r["3i_M"]["pass"] for r in rows.values())}
    report["arms"] = per_arm

    grid = vendor_modules()[0].SphericalCellGrid(round(60/DR), 60.)
    pairs = {}
    for coarse, fine in (("A1", "A2"), ("B1", "B2")):
        pair = {}
        for t in SAVES:
            key = f"{t:g}"
            with np.load(out / f"{coarse}-t{key}.npz") as L, np.load(out / f"{fine}-t{key}.npz") as R:
                ledger = pair_ledger(grid, dict(L), dict(R))
            obs_c, obs_f = arms[coarse]["saves"][key]["observe"], arms[fine]["saves"][key]["observe"]
            entry = {}
            for width, led in ledger.items():
                dH = obs_c[f"H_metric_{width}"][0] - obs_f[f"H_metric_{width}"][0]
                reading = gates.order_reading(dH, float(led["write"][0]), float(led["increment_arithmetic"][0]),
                                              float(led["closure_magnitude"][0]), obs_f[f"H_metric_{width}"][0])
                entry[width] = {
                    "order": reading, "dH_direct": dH,
                    "dR_direct_full": numbers(led["full"]),
                    "production_H_difference_incl_writes": obs_c["H_production"][0] - obs_f["H_production"][0],
                    "X_only": numbers(led["X_only"]), "X_parts": {k: numbers(v) for k, v in led["X_parts"].items()},
                    "X_closure": numbers(led["X_closure"]), "DT_total": numbers(led["DT_total"]),
                    "DT_parts": {k: numbers(v) for k, v in led["DT_parts"].items()},
                    "DT_closure": numbers(led["DT_closure"]), "interaction": numbers(led["interaction"]),
                    "write": numbers(led["write"]), "increment_arithmetic": numbers(led["increment_arithmetic"]),
                }
            pair[key] = entry
        pairs[f"{coarse}-{fine}"] = pair
    report["pairs"] = pairs
    for name, (arm, _) in ARMS.items():
        pair = pairs["A1-A2" if arm == "A" else "B1-B2"]
        fail_end = pair[f"{SAVES[-1]:g}"][9]["order"]["reading"] == "fail"
        per_arm[name]["3iii"] = gates.instability(per_arm[name]["finite"], per_arm[name]["growth"], fail_end)
    meeting = [ARMS[n][1] for n in ("A1", "A2") if per_arm[n]["3i_all_saves"]]
    report["largest_lambda_meeting_3i_ring"] = .03/min(meeting) if meeting else None
    print(json.dumps(report, indent=1, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--arm", choices=tuple(ARMS))
    args = parser.parse_args()
    if np.finfo(WIDE).nmant < 63:
        raise RuntimeError("at least 64 significand bits required")
    args.out.mkdir(parents=True, exist_ok=True)
    if args.arm:
        run_arm(args.arm, args.out)
    else:
        combine(args.out)
