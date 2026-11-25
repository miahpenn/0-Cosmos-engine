"""W1 PIRK corrector pairing: fixed-grid timestep refinement series (research only).

Design DESIGN-W1-PIRK-pairing-dt-refinement-Nov25-2026 (Drive 1MZrugyHkzJWhH3OeTl_VoW9LJjwaZEv7FAS3iRsw8GI) with
Amendment 1 (1u-9OTudwOuI74yEXIs8rjxIcYZrVYCwCA4F-nxjYVAo), Amendment 2 (1sgCtsP0--_qMa-fsjUDf8xF1_9ylNyAO1uIh0JVk2Y8) and
Amendment 3 (1AWx0eRcRlHsxqcyeo10-nAAeiUC86XY-HMfyo_qt-7c). Code only until Jeremiah's go. Nothing here changes the engine; arms run production arithmetic and production _d2.

  python -m research.w1_fring_smoke.pairing_dt_series runs OUT            command lines, longest first (for xargs -P 4)
  python -m research.w1_fring_smoke.pairing_dt_series run CASE ARM LEVEL OUT [--floor]
  python -m research.w1_fring_smoke.pairing_dt_series score OUT           input checks and gates G0-G2 -> OUT/score.json

Arms: P is production pairing; Q replaces the corrector's right-hand L2 for Aa and K with the primary-stage L2 (the
endpoint swap of corrector_pairing.py, A1.3 narrow form). Levels 2dt0, dt0, dt0/2, dt0/4 (20/40/80/160 steps) to
T = 0.015 (A2.1); t is summed per step (A1.7: the RHS does not read t; adaptive_step checks t advanced by dt).
Floor runs: a, b, X perturbed by +-1 ulp per cell after the start checks, signs from SEED (section 2).

Input checks, per run: C0 (P at dt0, step 1: archived start digest, stage ledgers and endpoint digest); C3 (Q, every
step: L2 recomputed on the primary predictor equals the substituted value bitwise, else the run stops, A1.3). C1
(step-1 primary predictor of P and Q, every field, bitwise, A1.2) and C2 (max|X_pred - X_expl| on r .05-.25 at the last
step, per level pair, A1.5, A2.3) are scored from the saved arrays. Observables (A1.4): shell_rms (8 shells,
r .05-.25) of M (unnormalised), Aa and K at accepted_end of the last step.
"""
import argparse
from copy import deepcopy
from dataclasses import fields
import hashlib
import json
import math
from pathlib import Path
from unittest.mock import patch

import numpy as np

from engine.stress_energy import assemble_total_stress_energy
from research.w1_fring_smoke import live_lapse, stage_momentum as replay
from research.w1_fring_smoke import capture_lapse_stages as capture
from research.w1_fring_smoke.corrector_pairing import PairEndpoint
from research.w1_fring_smoke.frozen_centre_derivative import shell_rms
from research.w1_fring_smoke.ricci_history import initial_state

DT0 = 0.000375
LEVELS = {"2dt0": (2 * DT0, 20), "dt0": (DT0, 40), "dt0_2": (DT0 / 2, 80), "dt0_4": (DT0 / 4, 160)}
ORDER = ("2dt0", "dt0", "dt0_2", "dt0_4")
ARMS = {"P": (), "Q": ("Aa", "K")}
CASES = ("0.025", "0.0125")
SEED = 20261125
PERTURBED = ("a", "b", "X")
FIELDS = ("M", "Aa", "K")
BAND = (0.05, 0.25)
HIGH, LOW, VALID = 1.7, 1.3, 10.0
REFERENCE = Path("results/research/w1_fring_live_ricci_precision")


class SeriesPair(PairEndpoint):
    """PairEndpoint over many steps: the swap counter and captures reset every step (A1.1)."""

    def __init__(self, vacuum, fields_):
        super().__init__(vacuum, fields_)
        self.capture_matter = False
        self.primary_state = None
        self.primary_grid = None

    def _geometry_stage_terms(self, grid, geometry, scalars, matter, electromagnetic):
        if self.geometry_calls == 2:
            self.primary_grid = grid
            if self.capture_matter:
                self.primary_state = (deepcopy(scalars), deepcopy(matter))
        return super()._geometry_stage_terms(grid, geometry, scalars, matter, electromagnetic)

    def step(self, state, dt):
        self.swaps = 0
        self.primary_geometry = None
        self.right_metric = None
        self.primary_grid = None
        return super().step(state, dt)


def group_arrays(prefix, group):
    return {f"{prefix}.{f.name}": np.asarray(getattr(group, f.name)) for f in fields(group)
            if isinstance(getattr(group, f.name), np.ndarray)}


def primary_arrays(kernel):
    scalars, matter = kernel.primary_state
    out = {**group_arrays("geometry", kernel.primary_geometry), **group_arrays("scalars", scalars)}
    for species in ("dark_matter", "baryons", "radiation"):
        out.update(group_arrays(f"matter.{species}", getattr(matter, species)))
    return out


def perturb(state):
    rng = np.random.default_rng(SEED)
    digest = hashlib.sha256()
    for name in PERTURBED:
        values = getattr(state.geometry, name)
        signs = rng.choice((-1.0, 1.0), size=values.shape)
        values[:] = np.nextafter(values, np.where(signs > 0, np.inf, -np.inf))
        digest.update(name.encode()); digest.update(signs.astype("<f8").tobytes())
    return {"seed": SEED, "fields": list(PERTURBED), "signs_sha256": digest.hexdigest()}


def band_mask(r):
    return (r >= BAND[0]) & (r < BAND[1])


def end_fields(vacuum, state):
    g, grid = state.geometry, state.grid
    total = assemble_total_stress_energy(grid, g, state.scalars, state.matter)
    M = vacuum.momentum_constraint(grid, g) - 8 * np.pi * total.j
    return {"M": np.asarray(M), "Aa": np.asarray(g.Aa), "K": np.asarray(g.K)}


def c0(step_ledger, end_state, reference):
    out = {"start_checked": True}
    expected = reference["production"]
    mismatched = [name for name, stage in step_ledger["stages"].items()
                  if stage != {k: v for k, v in expected["stages"][name].items() if k != "M_statistics"}]
    out["stage_mismatches"] = mismatched
    out["end_sha256"] = live_lapse.state_digest(end_state)
    out["end_matches"] = out["end_sha256"] == expected["end_sha256"]
    out["pass"] = not mismatched and out["end_matches"]
    return out


def run(case, arm, level, out, floor=False):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    name = f"{case}_{arm}_{level}" + ("_floor" if floor else "")
    dt, steps = LEVELS[level]
    meta = replay.provenance()
    go, vacuum, _ = replay.pins.ad.vendor_modules()
    reference = json.loads((REFERENCE / f"{case}.json").read_text())
    state, archived_dt = initial_state(case, go, vacuum, meta)
    if archived_dt != DT0:
        raise ValueError("archived first-step dt is not dt0")
    start = live_lapse.state_digest(state)
    if start != reference["start_sha256"]:
        raise ValueError("C0: start state differs from the archived initializer")
    record = {"case": case, "arm": arm, "level": level, "dt": dt, "steps": steps, "floor": floor,
              "start_sha256": start, "engine": replay.pins.SOURCE, "vendor": replay.pins.VENDOR,
              "script_sha256": replay.sha(Path(__file__)), "status": "running"}
    if floor:
        record["perturbation"] = perturb(state)
    kernel = SeriesPair(vacuum, ARMS[arm])
    r = np.asarray(state.grid.centers)
    mask = band_mask(r)
    shared_differs = 0
    history = []
    try:
        for n in range(1, steps + 1):
            kernel.capture_matter = n == 1 and not floor
            state, ledger = replay.advance(kernel, state, dt, True)
            if kernel.swaps != 1:
                raise ValueError("corrector L2 not isolated in this step")
            if ARMS[arm]:   # C3, every step (A1.3)
                fresh = vacuum.primary_l2_rhs(kernel.primary_grid, kernel.primary_geometry)
                for key in ARMS[arm]:
                    if not np.array_equal(np.asarray(fresh[key]), np.asarray(kernel.primary_l2[key])):
                        raise ValueError(f"C3 failed at step {n}: recomputed primary L2 {key} differs")
            shared_differs += any(not np.array_equal(getattr(kernel.primary_geometry, k), getattr(kernel.right_metric, k))
                                  for k in ("alpha", "beta", "Lambda", "Aa", "K"))
            if n == 1 and not floor:
                np.savez_compressed(out / f"{name}_primary1.npz", **primary_arrays(kernel))
                if arm == "P" and level == "dt0" and not floor:
                    record["C0"] = c0(ledger, state, reference)
                    if not record["C0"]["pass"]:
                        raise ValueError("C0 failed: step 1 does not reproduce the archive")
            now = end_fields(vacuum, state)
            if not all(np.all(np.isfinite(v)) for v in now.values()):
                raise FloatingPointError(f"non-finite field at step {n}")
            history.append({f: float(np.max(np.abs(now[f][mask]))) for f in FIELDS})   # A2.4 reading only
        gap = np.abs(np.asarray(kernel.primary_geometry.X) - np.asarray(kernel.right_metric.X))
        values = end_fields(vacuum, state)
        np.savez_compressed(out / f"{name}_end.npz", r=r, **values)
        record.update({
            "status": "complete", "t_end": state.t, "end_sha256": live_lapse.state_digest(state),
            "C3_checked_steps": steps if ARMS[arm] else 0,
            "steps_with_operands_differing_outside_abX": int(shared_differs),
            "C2_gap_band_max": float(np.max(gap[mask])), "C2_gap_all_max": float(np.max(gap)),
            "shell_rms": {f: shell_rms(r[:replay.LOCAL_CELLS], values[f][:replay.LOCAL_CELLS]) for f in FIELDS},
            "band_max_abs": {f: float(np.max(np.abs(values[f][mask]))) for f in FIELDS},
            "band_max_abs_by_step": history,
        })
    except Exception as exc:   # recorded; a 2dt0 failure reads UNDETERMINED (A2.4), any other failure is reported
        record.update({"status": "failed", "error": f"{type(exc).__name__}: {exc}", "band_max_abs_by_step": history})
        (out / f"{name}.json").write_text(json.dumps(record, indent=1, default=float))
        raise
    (out / f"{name}.json").write_text(json.dumps(record, indent=1, default=float))


def runs(out):
    lines = []
    for level in reversed(ORDER):       # longest first: fine grid, then most steps
        for case in reversed(CASES):
            for arm in ARMS:
                for floor in (False, True):
                    lines.append(f"run {case} {arm} {level} {out}" + (" --floor" if floor else ""))
    print("\n".join(lines))


# ---------- scoring (gates fixed before any data, design section 6, A1, A2) ----------

def load(out, case, arm, level, floor=False):
    name = f"{case}_{arm}_{level}" + ("_floor" if floor else "")
    js = Path(out) / f"{name}.json"
    if not js.exists():
        return None, None
    rec = json.loads(js.read_text())
    if rec["status"] != "complete":
        return rec, None
    with np.load(Path(out) / f"{name}_end.npz") as z:
        return rec, {k: np.array(z[k]) for k in ("r", *FIELDS)}


def norm(r, values):
    n = replay.LOCAL_CELLS
    return shell_rms(r[:n], values[:n])["rms"]


def ratio_verdicts(ps, high="SECOND_ORDER", low="FIRST_ORDER_SIGNAL"):
    counted = [p for p in ps if p["counted"]]
    if not counted:
        return "UNDETERMINED"
    if all(p["p"] >= HIGH for p in counted):
        return high
    if all(p["p"] <= LOW for p in counted):
        return low
    return "UNDETERMINED"


def blow_up(rec2, rec1):
    """A3.4: at each matched time (2dt0 step k = dt0 step 2k), any band max |F| above 2x the dt0 value is a blow-up.
    An error or non-finite value (a failed run) also counts; a missing run cannot be checked and reads as blown up."""
    if rec2 is None or rec1 is None or rec2.get("status") != "complete" or rec1.get("status") != "complete":
        return {"blown_up": True, "reason": "2dt0 or dt0 run failed or missing"}
    h2, h1 = rec2["band_max_abs_by_step"], rec1["band_max_abs_by_step"]
    worst = {f: max(h2[k][f] / h1[2 * k + 1][f] if h1[2 * k + 1][f] > 0 else math.inf for k in range(len(h2)))
             for f in FIELDS}
    exceeded = sorted(f for f, x in worst.items() if x > 2.0)
    return {"blown_up": bool(exceeded), "worst_ratio_to_dt0": worst, "exceeded": exceeded}


def score_case(out, case):
    res = {"input_checks": {}, "C2": {}, "floors": {}, "G1": {}, "G2": {}, "readings": {}}
    data, recs = {}, {}
    for arm in ARMS:
        for level in ORDER:
            for floor in (False, True):
                rec, arr = load(out, case, arm, level, floor)
                recs[(arm, level, floor)], data[(arm, level, floor)] = rec, arr
    # C0
    c0rec = recs[("P", "dt0", False)]
    res["input_checks"]["C0"] = (c0rec or {}).get("C0", {"pass": False, "note": "missing"})
    # C1: step-1 primary predictor, P vs Q, every field bitwise, at each level (unperturbed runs)
    c1 = {}
    for level in ORDER:
        p, q = (Path(out) / f"{case}_{a}_{level}_primary1.npz" for a in ARMS)
        if p.exists() and q.exists():
            with np.load(p) as zp, np.load(q) as zq:
                diff = sorted(k for k in set(zp.files) | set(zq.files)
                              if k not in zp.files or k not in zq.files or not np.array_equal(zp[k], zq[k]))
            c1[level] = {"fields_compared": len(set(zp.files) | set(zq.files)), "differing": diff, "pass": not diff}
        else:
            c1[level] = {"pass": False, "note": "missing"}
    res["input_checks"]["C1"] = c1
    # C3 and run status
    res["input_checks"]["runs"] = {f"{a}_{l}" + ("_floor" if f else ""): (rec or {}).get("status", "missing")
                                   for (a, l, f), rec in recs.items()}
    invalid = (not res["input_checks"]["C0"].get("pass") or not all(v["pass"] for v in c1.values() if v.get("note") != "missing")
               or any(rec is not None and rec["status"] == "failed" and l != "2dt0" for (a, l, f), rec in recs.items()))
    res["input_checks"]["series_valid"] = not invalid
    ok = {lv: all(data[(a, lv, f)] is not None for a in ARMS for f in (False, True)) for lv in ORDER}
    # C2 (P arm; Q recorded as a reading): log2 gap ratio per level pair >= 1.7
    for arm in ARMS:
        gaps = {lv: recs[(arm, lv, False)]["C2_gap_band_max"] for lv in ORDER if data[(arm, lv, False)] is not None}
        pairs = {}
        for a, b in zip(ORDER[:-1], ORDER[1:]):
            if a in gaps and b in gaps and gaps[a] > 0 and gaps[b] > 0:
                x = math.log2(gaps[a] / gaps[b])
                pairs[f"{a}/{b}"] = {"log2": x, "pass": x >= HIGH}
            else:
                pairs[f"{a}/{b}"] = {"pass": False, "note": "missing"}
        res["C2"][arm] = {"gap_band_max": gaps, "pairs": pairs}
    c2 = res["C2"]["P"]["pairs"]
    asymptotic = c2["dt0/dt0_2"]["pass"] and c2["dt0_2/dt0_4"]["pass"]   # original C2
    use_2dt0 = c2["2dt0/dt0"]["pass"] and ok["2dt0"]                       # A2.3, A2.4
    res["input_checks"]["C2_asymptotic"] = asymptotic
    res["input_checks"]["C2_2dt0_counts"] = use_2dt0
    blown = {arm: blow_up(recs[(arm, "2dt0", False)], recs[(arm, "dt0", False)]) for arm in ARMS}   # A3.4
    res["input_checks"]["blow_up_2dt0"] = blown
    for field in FIELDS:
        fl = {}
        for arm in ARMS:
            for lv in ORDER:
                u, v = data[(arm, lv, False)], data[(arm, lv, True)]
                if u is not None and v is not None:
                    fl[(arm, lv)] = norm(u["r"], u[field] - v[field])
        res["floors"][field] = {f"{a}_{l}": x for (a, l), x in fl.items()}

        def ratio(num, den, levels, floor_arms, label):
            r = data[("P", "dt0", False)]["r"] if data[("P", "dt0", False)] is not None else None
            if num is None or den is None or any((a, l) not in fl for a in floor_arms for l in levels):
                return {"ratio": label, "counted": False,
                        "note": "2dt0 run failed or missing (A2.4)" if "2dt0" in levels else "missing run"}
            nn, nd = norm(r, num), norm(r, den)
            threshold = VALID * max(fl[(a, l)] for a in floor_arms for l in levels)
            pre = "2dt0" in levels and (not use_2dt0 or any(blown[a]["blown_up"] for a in floor_arms))
            counted = nn >= threshold and nd >= threshold and asymptotic and not pre
            return {"ratio": label, "num": nn, "den": nd, "p": math.log2(nn / nd) if nd > 0 and nn > 0 else None,
                    "G0_threshold": threshold, "counted": bool(counted and nd > 0 and nn > 0),
                    **({"note": "2dt0 pre-asymptotic, failed or blown up (A2.3, A2.4, A3.4)"} if pre else {})}

        def F(arm, lv):
            d = data[(arm, lv, False)]
            return None if d is None else d[field]

        def sub(x, y):
            return None if x is None or y is None else x - y

        g1 = {}
        for arm in ARMS:
            E = {ORDER[i]: sub(F(arm, ORDER[i]), F(arm, ORDER[i + 1])) for i in range(3)}
            ps = [ratio(E[ORDER[i]], E[ORDER[i + 1]], ORDER[i:i + 3], (arm,), f"E({ORDER[i]})/E({ORDER[i+1]})")
                  for i in range(2)]
            g1[arm] = {"ratios": ps, "verdict": ratio_verdicts(ps)}
        res["G1"][field] = {"verdict": g1["P"]["verdict"], "P": g1["P"], "Q_reading": g1["Q"]}
        D = {lv: sub(F("Q", lv), F("P", lv)) for lv in ORDER}
        pd = [ratio(D[ORDER[i]], D[ORDER[i + 1]], ORDER[i:i + 2], ("P", "Q"), f"Dq({ORDER[i]})/Dq({ORDER[i+1]})")
              for i in range(3)]
        res["G2"][field] = {"ratios": pd, "verdict": ratio_verdicts(pd, "CONSISTENT", "DEFECT")}
        res["readings"][field] = {f"Dq_{lv}_band_max_abs": (float(np.max(np.abs(D[lv][band_mask(data[('P', lv, False)]['r'])])))
                                                         if D[lv] is not None else None) for lv in ORDER}
    if not res["input_checks"]["series_valid"]:
        res["note"] = "an input check failed or a non-2dt0 run failed: gates are not read"
    res["primary"] = {f: res["G2"][f]["verdict"] for f in ("M", "Aa")} if res["input_checks"]["series_valid"] else "INVALID"
    return res


def score(out):
    result = {case: score_case(out, case) for case in CASES
              if (Path(out) / f"{case}_P_dt0.json").exists()}
    (Path(out) / "score.json").write_text(json.dumps(result, indent=1, default=float))
    print(json.dumps({c: {"valid": v["input_checks"]["series_valid"], "C2_asymptotic": v["input_checks"]["C2_asymptotic"],
                          "G1": {f: v["G1"][f]["verdict"] for f in FIELDS},
                          "G2": {f: v["G2"][f]["verdict"] for f in FIELDS}} for c, v in result.items()}, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("run"); a.add_argument("case", choices=CASES); a.add_argument("arm", choices=tuple(ARMS))
    a.add_argument("level", choices=ORDER); a.add_argument("out"); a.add_argument("--floor", action="store_true")
    sub.add_parser("runs").add_argument("out")
    sub.add_parser("score").add_argument("out")
    args = p.parse_args()
    if args.cmd == "run":
        run(args.case, args.arm, args.level, args.out, args.floor)
    elif args.cmd == "runs":
        runs(args.out)
    else:
        score(args.out)


if __name__ == "__main__":
    main()
