"""Test a common primary endpoint for selected W1 corrector L2/L3 terms."""
import argparse
from copy import deepcopy
from decimal import localcontext
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

from research.w1_fring_smoke import live_lapse, matched_lapse, stage_momentum as replay
from research.w1_fring_smoke import capture_lapse_stages as capture
from research.w1_fring_smoke.live_ricci_precision import stats
from research.w1_fring_smoke.live_lapse_audit import EDGES
from research.w1_fring_smoke.first_step_rate_precision import diff_at, l2_cell
from research.w1_fring_smoke.momentum_precision import wide
from research.w1_fring_smoke.ricci_history import initial_state

ARMS = {"production": (), "shear": ("Aa",), "K": ("K",), "both": ("Aa", "K")}
CELLS = replay.LOCAL_CELLS


class PairEndpoint(replay.Observe):
    """Exchange only the selected right-hand terms in the live corrector."""
    def __init__(self, vacuum, fields):
        super().__init__(vacuum)
        self.fields = fields
        self.primary_l2 = None
        self.swaps = 0
        self.terms_before = None

    def _geometry_stage_terms(self, grid, geometry, scalars, matter, electromagnetic):
        number = self.geometry_calls
        term = super()._geometry_stage_terms(grid, geometry, scalars, matter, electromagnetic)
        if number == 2:
            self.primary_l2 = term["primary_l2"]
            self.primary_geometry = geometry.copy()
        return term

    def step(self, state, dt):
        self.primary_l2 = None
        original = self.vacuum.primary_l2_rhs

        def l2(grid, geometry):
            rate = original(grid, geometry)
            if self.primary_l2 is not None and self.boundary_calls == 1:
                self.swaps += 1
                self.right_metric = geometry.copy()
                self.terms_before = {key: rate[key].copy() for key in capture.PRIMARY_FIELDS}
                for key in self.fields:
                    rate[key] = self.primary_l2[key]
            return rate

        with patch.object(self.vacuum, "primary_l2_rhs", l2):
            result = super().step(state, dt)
        if self.swaps != 1:
            raise ValueError("did not isolate a single corrector evaluation")
        return result


def momentum_response(grid, vacuum, endpoint_metric, changes):
    r = grid.centers[:CELLS]
    a = changes["Aa"]
    k = changes["K"]
    derivative = lambda values: vacuum._d(grid, values, 1)[:CELLS]
    factor = (3/r - 3*derivative(endpoint_metric.X)/endpoint_metric.X[:CELLS]
              + 1.5*derivative(endpoint_metric.b)/endpoint_metric.b[:CELLS])
    return {"Aa": derivative(a) + factor*a[:CELLS], "K": -2/3*derivative(k)}


def profile(r, values):
    return {**stats(r, values), "signed": values.tolist()}


def pairing_attribution(grid, vacuum, primary, corrected, primary_l2, corrected_l2, dt):
    """Read the two operands with the same widened operator on both."""
    radius = grid.centers[:CELLS]
    n = CELLS - 2
    fields = ("a", "b", "X", "Lambda", "alpha")
    parity = (1, 1, 1, -1, 1)
    for field in ("Lambda", "Aa", "K", "alpha", "beta"):
        np.testing.assert_array_equal(getattr(primary,field),getattr(corrected,field))
    with localcontext() as ctx:
        ctx.prec = 60
        dr = wide(grid.dr)

        def evaluate(geometry):
            values = [list(map(wide, getattr(geometry,key)[:CELLS+2])) for key in fields]
            for i in range(CELLS):
                d1 = [diff_at(v,i,dr,parity=p) for v,p in zip(values,parity,strict=True)]
                d2 = [diff_at(values[j],i,dr,second=True) for j in (0,1,2,4)]
                yield l2_cell(wide(radius[i]),[v[i] for v in values],d1,d2)

        weight = wide(dt)/2
        delta = [[weight*(a-b) for a,b in zip(right,left,strict=True)]
                 for right,left in zip(evaluate(primary),evaluate(corrected),strict=True)]
        k, shear_hessian, shear_ricci = zip(*delta,strict=True)
        dx, db = (vacuum._d(grid,v,1)[:n] for v in (corrected.X,corrected.b))
        coefficient = 3/radius[:n] - 3*dx/corrected.X[:n] + 1.5*db/corrected.b[:n]

        def shear_response(values):
            return np.array([float(diff_at(values,i,dr) + wide(coefficient[i])*values[i])
                             for i in range(n)])

        widened = {"K": np.array([float(-2*diff_at(k,i,dr)/3) for i in range(n)]),
                   "Aa_hessian": shear_response(shear_hessian),
                   "Aa_ricci": shear_response(shear_ricci)}
    binary = momentum_response(grid,vacuum,corrected,{
        key: .5*dt*(primary_l2[key]-corrected_l2[key]) for key in capture.PRIMARY_FIELDS
    })
    total = sum(widened.values())
    return {
        "widened": {key: profile(radius[:n],v) for key,v in widened.items()},
        "widened_sum": profile(radius[:n],total),
        "binary64_minus_widened": profile(radius[:n],sum(binary.values())[:n]-total),
        "metric_operand_changes": {
            key: {"max_absolute": float(np.max(abs(getattr(primary,key)[:CELLS]-getattr(corrected,key)[:CELLS]))),
                  "changed_cells": int(np.count_nonzero(getattr(primary,key)[:CELLS]!=getattr(corrected,key)[:CELLS]))}
            for key in ("a","b","X")
        },
    }


def first_step(case, go, vacuum, meta, archived):
    state, dt = initial_state(case, go, vacuum, meta)
    grid = state.grid
    r = grid.centers[:CELLS]
    digest = live_lapse.state_digest(state)
    if digest != archived["start_sha256"]:
        raise ValueError("first step has a different initializer")
    baseline, step = replay.advance(replay.Observe(vacuum), deepcopy(state), dt, True)
    if live_lapse.state_digest(baseline) != archived["production"]["end_sha256"]:
        raise ValueError("unmodified branch differs from the archived endpoint")
    m0 = np.asarray(step["stages"]["accepted_end"]["rows"]["M"])
    reference = np.asarray(step["stages"]["initial_lapse"]["rows"]["M"])
    pre_corrector = ("accepted_start", "initial_lapse", "explicit_predictor",
                     "primary_before_boundary", "primary_after_lapse")
    output = {
        "t_start": state.t, "t_end": baseline.t, "start_sha256": digest,
        "production": {"M": profile(r,m0), "step_delta_M": profile(r,m0-reference),
                       "constraints": live_lapse.end_constraints(vacuum,baseline)},
    }
    for arm, fields in ARMS.items():
        if not fields:
            continue
        observer = PairEndpoint(vacuum,fields)
        branch, result = replay.advance(observer,deepcopy(state),dt,True)
        if any(result["stages"][key] != step["stages"][key] for key in pre_corrector):
            raise ValueError("the pairing intervention changed a predictor")
        changes = live_lapse.compare_stages(r, step["stages"], result["stages"])
        m1 = np.asarray(result["stages"]["accepted_end"]["rows"]["M"])
        increments = {
            key: .5*dt*(observer.primary_l2[key]-observer.terms_before[key])
                 if key in fields else np.zeros(grid.n)
            for key in capture.PRIMARY_FIELDS
        }
        readback = momentum_response(grid,vacuum,baseline.geometry,increments)
        linear = sum(readback.values())
        output[arm] = {
            "M": profile(r,m1), "delta_M": profile(r,m1-m0),
            "step_delta_M": profile(r,m1-reference),
            "prewrite_response": profile(r,linear),
            "live_minus_prewrite_response": profile(r,m1-m0-linear),
            "constraints": live_lapse.end_constraints(vacuum,branch),
            "corrector_comparison": changes["accepted_end"],
        }
    output["operator_attribution"] = pairing_attribution(
        grid,vacuum,observer.primary_geometry,observer.right_metric,
        observer.primary_l2,observer.terms_before,dt,
    )
    if live_lapse.state_digest(state) != digest:
        raise ValueError("shared initializer changed")
    return output


def late_step(snapshot, grid, vacuum):
    dt = float(snapshot["dt"])
    r = grid.centers[:CELLS]
    metric = matched_lapse.fields_from(snapshot,"corrector_copy",vacuum)
    primary = matched_lapse.fields_from(snapshot,"primary",vacuum)
    np.testing.assert_array_equal(metric.alpha,primary.alpha)
    right = vacuum.primary_l2_rhs(grid,primary)
    original = vacuum.primary_l2_rhs(grid,metric)
    names = ("start_primary_l2","start_primary_l3","primary_primary_l3")

    def written(key, endpoint):
        first, l3_initial, l3_primary = (snapshot[f"{name}_{key}"] for name in names)
        return snapshot[f"start_{key}"] + .5*dt*(first+endpoint+l3_initial+l3_primary)

    before = metric.copy()
    for key in capture.PRIMARY_FIELDS:
        np.testing.assert_array_equal(original[key],snapshot[f"corrector_copy_primary_l2_{key}"])
        np.testing.assert_array_equal(right[key],snapshot[f"primary_primary_l2_{key}"])
        value = written(key,original[key])
        np.testing.assert_array_equal(value[:CELLS],snapshot[f"corrector_before_boundary_{key}"])
        setattr(before,key,value)
    baseline = vacuum.momentum_constraint(grid,before)[:CELLS]
    m0 = snapshot["corrector_before_boundary_M"]
    reference = snapshot["initial_lapse_M"]
    output = {"t_start": float(snapshot["t_start"]), "t_end": float(snapshot["t_end"]),
              "production": {"M": profile(r,m0), "step_delta_M": profile(r,m0-reference)},
              "operator_attribution": pairing_attribution(grid,vacuum,primary,metric,right,original,dt)}
    for arm, fields in ARMS.items():
        if not fields:
            continue
        after = before.copy()
        for key in fields:
            setattr(after,key,written(key,right[key]))
        direct = vacuum.momentum_constraint(grid,after)[:CELLS]-baseline
        increments = {
            key: .5*dt*(right[key]-original[key]) if key in fields else np.zeros(grid.n)
            for key in capture.PRIMARY_FIELDS
        }
        readback = momentum_response(grid,vacuum,before,increments)
        output[arm] = {
            "M": profile(r,m0+direct), "delta_M": profile(r,direct),
            "step_delta_M": profile(r,m0+direct-reference),
            "prewrite_response": profile(r,sum(readback.values())),
            "written_minus_prewrite_response": profile(r,direct-sum(readback.values())),
        }
    return output


def compare_timestep(first, timing):
    t = first["t_end"]
    if first["start_sha256"] != timing["start_sha256"]:
        raise ValueError("timestep comparator initializer differs")
    old, half = (timing["result"][f"production_dt_over_{i}"]["rows"][-1] for i in (1,2))
    if old["t_end"] != t or half["t_end"] != t:
        raise ValueError("timing and pairing endpoints differ")
    full = old["stages"]["accepted_end"]["M"]
    alternative = half["stages"]["accepted_end"]["M"]
    np.testing.assert_array_equal(first["production"]["M"]["signed"],full["signed"])
    full_shells = np.asarray(full["signs"]["shell_means"])
    half_shells = np.asarray(alternative["signs"]["shell_means"])
    timing_delta = half_shells-full_shells
    pairing_delta = np.asarray(first["both"]["delta_M"]["signs"]["shell_means"])
    w = np.diff(EDGES**3)

    def inner(a, b):
        return np.average(a*b,weights=w)

    timing_rms = np.sqrt(inner(timing_delta,timing_delta))
    pairing_rms = np.sqrt(inner(pairing_delta,pairing_delta))
    remainder = pairing_delta-timing_delta
    return {
        "timing_delta_shells": timing_delta.tolist(),
        "pairing_delta_shells": pairing_delta.tolist(),
        "timing_delta_rms": float(timing_rms),
        "pairing_delta_rms": float(pairing_rms),
        "pairing_minus_timing_delta_rms": float(np.sqrt(inner(remainder,remainder))),
        "signed_cosine": float(inner(pairing_delta,timing_delta)/(pairing_rms*timing_rms)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",choices=("0.025","0.0125"))
    parser.add_argument("late_snapshot",type=Path)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    meta = replay.provenance()
    go,vacuum,_ = replay.pins.ad.vendor_modules()
    root = Path("results/research")
    first_path = root/"w1_fring_live_ricci_precision"/f"{args.case}.json"
    first = json.loads(first_path.read_text())
    late_path = root/"w1_fring_live_lapse"/f"{args.case}.json"
    late = json.loads(late_path.read_text())
    with np.load(args.late_snapshot,allow_pickle=False) as snapshot:
        for key,value in (("engine",replay.pins.SOURCE),("vendor",replay.pins.VENDOR),
                          ("input_sha256",first["input_sha256"]),("dr",first["dr"]),("dt",first["dt"])):
            if snapshot[key].item() != value:
                raise ValueError(f"late snapshot differs in {key}")
        if snapshot["capture_sha256"].item() != replay.sha(Path(capture.__file__)):
            raise ValueError("late snapshot capture source differs")
        end = late["production"]["stages"]["corrector_before_boundary"]["rows"]
        np.testing.assert_array_equal(snapshot["corrector_before_boundary_M"],end["M"])
        grid = go.SphericalCellGrid(len(snapshot["start_K"]),60.)
        late_comparison = late_step(snapshot,grid,vacuum)
    first_comparison = first_step(args.case,go,vacuum,meta,first)
    timing_path = root/"w1_fring_ricci_history"/f"{args.case}-timing.json"
    timing = compare_timestep(first_comparison,json.loads(timing_path.read_text()))
    report = {
        "dr": first["dr"], "dt": first["dt"], "input_sha256": first["input_sha256"],
        "r": first["r"], "engine": replay.pins.SOURCE, "vendor": replay.pins.VENDOR,
        "script_sha256": replay.sha(Path(__file__)), "frozen_sha256": replay.sha(args.late_snapshot),
        "first_reference_sha256": replay.sha(first_path), "late_reference_sha256": replay.sha(late_path),
        "first": first_comparison, "late": late_comparison,
        "timing_reference_sha256": replay.sha(timing_path), "timing_comparison": timing,
    }
    args.output.write_text(json.dumps(report,allow_nan=False,separators=(",", ":"))+"\n")


if __name__ == "__main__":
    main()
