"""Full R80 evolution using the archived R160 global CMC Kdot history.

Diagnostic counterfactual: preserves the production matter/geometry equations,
uses the same R80 matter-face history and R160 outer-alpha history as the
boundary-matched run, but prescribes the interpolated R160 full-domain Kdot
target inside each CMC lapse solve. This is a gauge-target experiment only;
it does not edit production source and must not be read as new physics.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import traceback
import hashlib

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import engine.cmc_gauge as cmc
from engine.production_kernel import V55ProductionKernel

HELPER_PATH = ROOT / "diagnostics" / "curvature_consistency" / "boundary-matched-cmc-r80-evolution.py"
_spec = importlib.util.spec_from_file_location("boundary_matched_cmc_helper", HELPER_PATH)
if _spec is None or _spec.loader is None:
    raise RuntimeError(f"cannot load shared matched-boundary diagnostic helpers: {HELPER_PATH}")
helper = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(helper)

RESOLUTION = 160
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-4
INCLUDE_RADIATION = True
CFL = 0.0075
FINAL_TIME = 12.0
SAMPLE_INTERVAL = 0.25
TARGET_TIMES = tuple(float(x) for x in np.arange(
    SAMPLE_INTERVAL, FINAL_TIME + 0.5 * SAMPLE_INTERVAL, SAMPLE_INTERVAL
))
OUT_DIR = ROOT / "runs" / "cmc-r80-r160-kdot-counterfactual"


class KdotSeries:
    def __init__(self, rows, sample_stride: int = 1):
        if not isinstance(sample_stride, int) or sample_stride < 1:
            raise ValueError("sample_stride must be a positive integer")
        if len(rows) < 2:
            raise ValueError("R160 Kdot series needs at least two archived rows")
        retained = list(rows[::sample_stride])
        # Always preserve the final archived endpoint, even when N-1 is not
        # divisible by sample_stride. This keeps the t=12 boundary anchored.
        if retained[-1]["t"] != rows[-1]["t"]:
            retained.append(rows[-1])
        self.source_sample_count = int(len(rows))
        self.sample_stride = int(sample_stride)
        self.t = np.asarray([row["t"] for row in retained], dtype=float)
        self.kdot = np.asarray([row["target_Kdot_full_domain"] for row in retained], dtype=float)
        if len(self.t) < 2 or not np.all(np.isfinite(self.t)) or not np.all(np.isfinite(self.kdot)):
            raise ValueError("R160 Kdot series is missing or non-finite")
        if not np.all(np.diff(self.t) > 0.0):
            raise ValueError("R160 Kdot reference times are not strictly increasing")

    def __call__(self, t):
        return float(np.interp(float(t), self.t, self.kdot))

    def describe(self):
        return {
            "source_sample_count": self.source_sample_count,
            "sample_count": int(len(self.t)),
            "sample_stride": self.sample_stride,
            "t_range": [float(self.t[0]), float(self.t[-1])],
            "kdot_range": [float(np.min(self.kdot)), float(np.max(self.kdot))],
            "kdot_t0": float(self.kdot[0]),
            "kdot_t12": float(self.kdot[-1]),
            "interpolation": "linear in time through retained archived R160 full-domain target_Kdot_full_domain samples; first and final endpoints preserved",
        }


class R160TargetKernel(V55ProductionKernel):
    def __init__(self, alpha_series, face_series, kdot_series):
        self._active_stage_time = 0.0
        self._kdot_series = kdot_series
        self.target_call_log = []
        self.face_call_log = []

        def boundary(t):
            self.face_call_log.append(float(t))
            return face_series(float(t))

        super().__init__(
            diagnostic_outer_alpha=alpha_series,
            diagnostic_matter_boundary=boundary,
        )

    def _outer_alpha_at(self, t: float):
        # The inherited step evaluates this immediately before each lapse solve;
        # that stage time selects the matching archived global Kdot.
        self._active_stage_time = float(t)
        return super()._outer_alpha_at(t)

    def _solve_lapse(self, grid, geometry, scalars, matter, outer_alpha=None):
        stage_t = float(self._active_stage_time)
        prescribed_kdot = self._kdot_series(stage_t)
        original = cmc.target_kdot

        def target_override(_grid, _geometry, _scalars, _matter, outer_frac=0.20):
            return prescribed_kdot

        cmc.target_kdot = target_override
        try:
            alpha, used = cmc.solve_cmc_lapse(
                grid, geometry, scalars, matter, outer_alpha=outer_alpha
            )
        finally:
            cmc.target_kdot = original
        self.target_call_log.append({
            "stage_t": stage_t,
            "target_Kdot_R160": prescribed_kdot,
            "solver_returned_target": float(used),
            "outer_alpha": None if outer_alpha is None else float(outer_alpha),
        })
        return alpha, used


def write_progress(state, records, accepted_steps, status, failure=None):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "cmc_r80_r160_kdot_counterfactual_progress_v1",
        "diagnostic_only": True,
        "status": status,
        "accepted_steps": int(accepted_steps),
        "sample_count": len(records),
        "sample_times": [float(row["t"]) for row in records],
        "current_t": float(state.t),
        "current_tau": float(state.tau),
        "failure": failure,
    }
    (OUT_DIR / "campaign_progress.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )


def main():
    if len(sys.argv) not in (3, 4):
        raise SystemExit(
            "usage: python diagnostics/curvature_consistency/cmc-r80-r160-kdot-counterfactual.py <r80-baseline-profiles.jsonl> <r160-reference-profiles.jsonl> [target_sample_stride]"
        )
    r80_path, r160_path = map(Path, sys.argv[1:3])
    try:
        target_sample_stride = int(sys.argv[3]) if len(sys.argv) == 4 else 1
    except ValueError as exc:
        raise SystemExit("target_sample_stride must be a positive integer") from exc
    if target_sample_stride < 1:
        raise SystemExit("target_sample_stride must be a positive integer")
    r80_rows = helper.read_jsonl(str(r80_path))
    r160_rows = helper.read_jsonl(str(r160_path))
    if not np.allclose(
        [row["t"] for row in r80_rows],
        [row["t"] for row in r160_rows],
        rtol=0.0,
        atol=1.0e-12,
    ):
        raise ValueError("R80 and R160 target-history sample times do not match")
    if int(r160_rows[0]["full_N"]) != 320 or abs(float(r160_rows[0]["full_r_max"]) - 160.0) > 1.0e-12:
        raise ValueError("unexpected R160 full-domain reference")

    alpha_series = helper.BoundaryAlphaSeries()
    face_series = helper.R80BoundaryFaceSeries(r80_path)
    kdot_series = KdotSeries(r160_rows, sample_stride=target_sample_stride)
    kernel = R160TargetKernel(alpha_series, face_series, kdot_series)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    profiles_path = OUT_DIR / "cmc_r80_r160_kdot_counterfactual_profiles.jsonl"
    profiles_path.write_text("")

    state = kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        amplitude=AMPLITUDE,
        width=WIDTH,
        D_amplitude=D_AMPLITUDE,
        include_radiation=INCLUDE_RADIATION,
    )
    if abs(float(state.grid.dr) - 0.5) > 1.0e-13:
        raise ValueError(f"unexpected grid spacing {state.grid.dr}")
    initial_boundary_error = abs(
        float(state.geometry.alpha[-1]) - alpha_series(0.0)
    )
    if initial_boundary_error > 1.0e-12:
        raise AssertionError(
            f"initial CMC solve did not impose reference outer alpha: {initial_boundary_error}"
        )

    first = helper.snapshot_row(state, 0.0, "t=0")
    records = [first]
    with profiles_path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(first, sort_keys=True, allow_nan=False) + "\n")
        stream.flush()
    write_progress(state, records, 0, "in_progress")

    dt_nominal = CFL * float(state.grid.dr)
    accepted_steps = 0
    target_index = 0
    failure = None

    while state.t < FINAL_TIME - 1.0e-13:
        dt = min(dt_nominal, FINAL_TIME - state.t)
        try:
            next_state = kernel.step(state, dt)
            next_state.geometry.assert_finite_positive()
            state = next_state
            accepted_steps += 1
        except Exception as exc:
            failure = {
                "stage": "evolution_step",
                "step_attempt": accepted_steps + 1,
                "t_before": float(state.t),
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(limit=16),
            }
            break

        if accepted_steps % 128 == 0:
            print("[R160_KDOT_CMC_R80] HEARTBEAT " + json.dumps({
                "t": float(state.t), "accepted_steps": accepted_steps,
                "tau": float(state.tau),
            }, sort_keys=True), flush=True)

        while target_index < len(TARGET_TIMES) and (
            state.t + 0.5 * dt >= TARGET_TIMES[target_index]
            or abs(state.t - TARGET_TIMES[target_index]) < 1.0e-12
        ):
            target = TARGET_TIMES[target_index]
            row = helper.snapshot_row(state, target, f"target={target:g}")
            records.append(row)
            with profiles_path.open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
                stream.flush()
            target_index += 1
            write_progress(state, records, accepted_steps, "in_progress", failure)
            print("[R160_KDOT_CMC_R80] CHECKPOINT " + json.dumps({
                "t": float(state.t), "accepted_steps": accepted_steps,
                "sample_count": len(records), "tau": float(state.tau),
                "center_alpha": float(state.geometry.alpha[0]),
                "target_Kdot": kdot_series(float(state.t)),
            }, sort_keys=True), flush=True)

    complete = (
        failure is None
        and state.t >= FINAL_TIME - 1.0e-10
        and len(records) == len(TARGET_TIMES) + 1
    )
    status = "completed" if complete else "numerical_or_sampling_failure"

    comparison = None
    if len(records) == len(r80_rows) == len(r160_rows):
        comparison = helper.compare_records(records, r80_rows, r160_rows)

    stage_times = np.asarray([x["stage_t"] for x in kernel.target_call_log], dtype=float)
    prescribed = np.asarray([x["target_Kdot_R160"] for x in kernel.target_call_log], dtype=float)
    returned = np.asarray([x["solver_returned_target"] for x in kernel.target_call_log], dtype=float)
    target_call_error = float(np.max(np.abs(prescribed - returned))) if len(prescribed) else None
    report = {
        "schema": "cmc_r80_r160_kdot_counterfactual_v1",
        "status": status,
        "diagnostic_only": True,
        "code_commit": __import__("os").environ.get("GITHUB_SHA", "local"),
        "production_physics_changed": False,
        "production_kernel_source_changed": False,
        "settings": {
            "N": RESOLUTION,
            "r_max": R_MAX,
            "dr": float(state.grid.dr),
            "CFL": CFL,
            "final_time": FINAL_TIME,
            "sample_interval": SAMPLE_INTERVAL,
            "amplitude": AMPLITUDE,
            "width": WIDTH,
            "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "evolution_matter_boundary": "archived baseline R80 outer-cell a,b,X,alpha,beta, linearly interpolated",
            "CMC_outer_alpha": "archived R160 alpha at r=79.75, linearly interpolated",
            "CMC_target_Kdot": f"archived R160 full-domain target_Kdot_full_domain, every-{target_sample_stride}th archived target sample (final endpoint preserved), linearly interpolated and prescribed at every CMC lapse solve",
        },
        "reference_sources": {
            "R80_profile": str(r80_path),
            "R80_profile_sha256": hashlib.sha256(r80_path.read_bytes()).hexdigest(),
            "R160_profile": str(r160_path),
            "R160_profile_sha256": hashlib.sha256(r160_path.read_bytes()).hexdigest(),
            "R160_target_Kdot_series": kdot_series.describe(),
            "R160_alpha_boundary_series": alpha_series.describe(),
            "R80_matter_boundary_series": {
                "times": int(len(face_series.t)),
                "t_range": [float(face_series.t[0]), float(face_series.t[-1])],
                "source_profile_sha256": hashlib.sha256(r80_path.read_bytes()).hexdigest(),
            },
        },
        "accepted_steps": accepted_steps,
        "sample_count": len(records),
        "final_t": float(state.t),
        "final_tau": float(state.tau),
        "final_center_alpha": float(state.geometry.alpha[0]),
        "initial_outer_alpha_error": initial_boundary_error,
        "lapse_solve_target_calls": int(len(kernel.target_call_log)),
        "lapse_solve_stage_time_range": (
            [float(np.min(stage_times)), float(np.max(stage_times))]
            if len(stage_times) else None
        ),
        "max_target_value_delivery_error": target_call_error,
        "prescribed_target_values_all_finite": bool(
            len(prescribed) > 0 and np.all(np.isfinite(prescribed))
        ),
        "matter_boundary_provider_calls": int(len(kernel.face_call_log)),
        "failure": failure,
        "comparison": comparison,
        "interpretation_limit": (
            "This is a full evolving gauge-target counterfactual. It holds an interpolated "
            "version of the archived R160 Kdot history and does not constitute a self-consistent "
            "R160 solution on the R80 domain. Any changed tau is causal evidence for sensitivity "
            "to this prescribed target history within the counterfactual, not proof that it is the "
            "unique physical gauge. The report records retained/source target samples so interpolation "
            "resolution can be compared directly across runs."
        ),
    }
    report_path = OUT_DIR / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    if comparison is not None:
        (OUT_DIR / "comparison.json").write_text(
            json.dumps(comparison, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
    write_progress(state, records, accepted_steps, status, failure)
    print("[R160_KDOT_CMC_R80] FINAL " + json.dumps({
        "status": status,
        "accepted_steps": accepted_steps,
        "sample_count": len(records),
        "final_t": float(state.t),
        "final_tau": float(state.tau),
        "final_tau_minus_R160": (
            None if not r160_rows else float(state.tau - float(r160_rows[-1]["tau"]))
        ),
        "target_call_count": len(kernel.target_call_log),
        "max_target_delivery_error": target_call_error,
        "report": str(report_path),
    }, sort_keys=True), flush=True)
    if status != "completed":
        raise RuntimeError("R160-target CMC R80 trajectory did not complete")


if __name__ == "__main__":
    main()
