"""Separated boundary-matched R=80 evolution and archived-state comparison.

Run only after the admission workflow succeeds. The normal mode integrates one
R=80 trajectory to t=12 and compares all saved slices to the archived R=80
baseline and the archived R=160 trajectory cropped to the shared grid.
--smoke-only uses the same real boundary histories but advances just one step.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pathlib
import sys
import traceback

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine.production_kernel import V55ProductionKernel
from diagnostics.curvature_consistency.boundary_alpha_series import (
    BoundaryAlphaSeries,
    R80BoundaryFaceSeries,
)
from diagnostics.curvature_consistency.strong_D_evolved_state_r80_capture import (
    GEOMETRY_ARRAYS,
    SCALAR_ARRAYS,
    MATTER_NAMES,
    CONSERVED_ARRAYS,
    state_snapshot,
)

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
OUT_DIR = ROOT / "runs" / "boundary-matched-cmc-r80"


def read_jsonl(path):
    with open(path) as stream:
        rows = [json.loads(line) for line in stream if line.strip()]
    if len(rows) != 49:
        raise ValueError(f"expected 49 archived samples in {path}, found {len(rows)}")
    times = np.asarray([row["t"] for row in rows], dtype=float)
    if not np.all(np.isfinite(times)) or not np.all(np.diff(times) > 0):
        raise ValueError(f"invalid time series in {path}")
    return rows


def snapshot_row(state, requested_time, label):
    snap = state_snapshot(state)
    for group, values in snap.items():
        entries = values.items()
        for name, series in entries:
            if group == "matter":
                for field, arr in series.items():
                    if not np.all(np.isfinite(np.asarray(arr, dtype=float))):
                        raise FloatingPointError(f"non-finite {group}.{name}.{field}")
            elif not np.all(np.isfinite(np.asarray(series, dtype=float))):
                raise FloatingPointError(f"non-finite {group}.{name}")
    return {
        "schema": "boundary_matched_cmc_r80_v1",
        "diagnostic_only": True,
        "label": label,
        "requested_time": float(requested_time),
        "t": float(state.t),
        "tau": float(state.tau),
        "e_folds": float(state.e_folds),
        "domain": {"N": int(state.grid.n), "r_max": float(R_MAX), "dr": float(state.grid.dr)},
        "radius": [float(x) for x in np.asarray(state.grid.centers, dtype=float)],
        "state_snapshot": snap,
    }


def write_progress(state, records, accepted_steps, status, failure=None):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "boundary_matched_cmc_r80_progress_v1",
        "status": status,
        "diagnostic_only": True,
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


def make_kernel(alpha_series, face_series, face_call_log=None):
    def boundary(t):
        values = face_series(float(t))
        if face_call_log is not None:
            face_call_log.append(float(t))
        return values
    return V55ProductionKernel(
        diagnostic_outer_alpha=alpha_series,
        diagnostic_matter_boundary=boundary,
    )


def smoke(r80_path, r160_path):
    baseline_rows = read_jsonl(r80_path)
    r160_rows = read_jsonl(r160_path)
    if not np.allclose(
        [x["t"] for x in baseline_rows], [x["t"] for x in r160_rows],
        rtol=0.0, atol=1.0e-12,
    ):
        raise ValueError("R80/R160 reference sample times do not match")
    alpha_series = BoundaryAlphaSeries()
    face_series = R80BoundaryFaceSeries(r80_path)
    call_log = []

    baseline_kernel = V55ProductionKernel()
    baseline = baseline_kernel.initialize(
        resolution=RESOLUTION, r_max=R_MAX, amplitude=AMPLITUDE, width=WIDTH,
        D_amplitude=D_AMPLITUDE, include_radiation=INCLUDE_RADIATION,
    )
    matched_kernel = make_kernel(alpha_series, face_series, call_log)
    matched = matched_kernel.initialize(
        resolution=RESOLUTION, r_max=R_MAX, amplitude=AMPLITUDE, width=WIDTH,
        D_amplitude=D_AMPLITUDE, include_radiation=INCLUDE_RADIATION,
    )
    prescribed_initial = alpha_series(0.0)
    if abs(float(matched.geometry.alpha[-1]) - prescribed_initial) > 1.0e-14:
        raise AssertionError("initial CMC solve did not impose reference outer alpha")

    dt = CFL * float(matched.grid.dr)
    base_next = baseline_kernel.step(baseline, dt)
    matched_next = matched_kernel.step(matched, dt)
    base_next.geometry.assert_finite_positive()
    matched_next.geometry.assert_finite_positive()
    prescribed_next = alpha_series(float(matched_next.t))
    if abs(float(matched_next.geometry.alpha[-1]) - prescribed_next) > 1.0e-14:
        raise AssertionError("final stage CMC solve did not impose the stage-time reference alpha")
    center_alpha_delta = float(matched_next.geometry.alpha[0] - base_next.geometry.alpha[0])
    if abs(center_alpha_delta) <= 1.0e-8:
        raise AssertionError("override did not measurably move the centre lapse in one step")
    if len(call_log) < 2 or not any(abs(t - dt) < 1.0e-12 for t in call_log):
        raise AssertionError("matter-boundary provider was not called at the end-of-step stage")
    result = {
        "status": "passed",
        "resolution": RESOLUTION,
        "dt": float(dt),
        "baseline_t": float(base_next.t),
        "matched_outer_alpha_t0": float(matched.geometry.alpha[-1]),
        "matched_outer_alpha_t1": float(matched_next.geometry.alpha[-1]),
        "prescribed_outer_alpha_t1": float(prescribed_next),
        "center_alpha_delta_after_one_step": center_alpha_delta,
        "matter_boundary_provider_calls": len(call_log),
        "matter_boundary_provider_time_min": float(min(call_log)),
        "matter_boundary_provider_time_max": float(max(call_log)),
        "interpretation": "one-step implementation gate only; not a physical result",
    }
    print("[BOUNDARY_MATCHED_CMC_SMOKE] " + json.dumps(result, sort_keys=True), flush=True)
    return result


def flatten_snapshot(snapshot):
    flat = {}
    for group in ("geometry", "scalars"):
        for name, values in snapshot[group].items():
            flat[f"{group}.{name}"] = np.asarray(values, dtype=float)
    for species, fields in snapshot["matter"].items():
        for name, values in fields.items():
            flat[f"matter.{species}.{name}"] = np.asarray(values, dtype=float)
    return flat


def compare_records(matched_rows, r80_rows, r160_rows):
    if not (len(matched_rows) == len(r80_rows) == len(r160_rows)):
        raise ValueError("matched/R80/R160 sample counts differ")
    for i, (m, b, f) in enumerate(zip(matched_rows, r80_rows, r160_rows)):
        if max(abs(float(m["t"]) - float(b["t"])), abs(float(m["t"]) - float(f["t"]))) > 1.0e-12:
            raise ValueError(f"sample-time mismatch at sample {i}")
    radii = np.asarray(r80_rows[-1]["radius"], dtype=float)
    r160_radii = np.asarray(r160_rows[-1]["radius_full"], dtype=float)[:RESOLUTION]
    if not np.array_equal(radii, r160_radii):
        raise ValueError("R80 and cropped R160 cell-centre radii are not exactly identical")
    core = radii <= 20.0
    final_index = len(matched_rows) - 1

    baseline_tau_delta = np.asarray([b["tau"] - f["tau"] for b, f in zip(r80_rows, r160_rows)])
    matched_tau_delta = np.asarray([m["tau"] - f["tau"] for m, f in zip(matched_rows, r160_rows)])
    fields = flatten_snapshot(matched_rows[0]["state_snapshot"]).keys()
    result = {}
    for path in fields:
        base_diffs = []
        matched_diffs = []
        baseline_final = matched_final = None
        for i, (m, b, f) in enumerate(zip(matched_rows, r80_rows, r160_rows)):
            bm = flatten_snapshot(b["state_snapshot"])[path]
            mm = flatten_snapshot(m["state_snapshot"])[path]
            fm = flatten_snapshot(f["state_snapshot"])[path][:RESOLUTION]
            if bm.shape != mm.shape or bm.shape != fm.shape:
                raise ValueError(f"field shape mismatch at {path}, sample {i}")
            base_diffs.append(bm - fm)
            matched_diffs.append(mm - fm)
            if i == final_index:
                baseline_final, matched_final = bm - fm, mm - fm
        base_stack = np.concatenate(base_diffs)
        matched_stack = np.concatenate(matched_diffs)
        base_abs = np.abs(baseline_final)
        matched_abs = np.abs(matched_final)
        base_idx = int(np.argmax(base_abs))
        matched_idx = int(np.argmax(matched_abs))
        base_max = float(base_abs.max())
        matched_max = float(matched_abs.max())
        result[path] = {
            "all_saved_slices_max_abs_baseline_minus_R160": float(np.max(np.abs(base_stack))),
            "all_saved_slices_max_abs_matched_minus_R160": float(np.max(np.abs(matched_stack))),
            "all_saved_slices_rms_baseline_minus_R160": float(np.sqrt(np.mean(base_stack * base_stack))),
            "all_saved_slices_rms_matched_minus_R160": float(np.sqrt(np.mean(matched_stack * matched_stack))),
            "t12_max_abs_baseline_minus_R160": base_max,
            "t12_max_abs_matched_minus_R160": matched_max,
            "t12_rms_baseline_minus_R160": float(np.sqrt(np.mean(baseline_final * baseline_final))),
            "t12_rms_matched_minus_R160": float(np.sqrt(np.mean(matched_final * matched_final))),
            "t12_core_max_abs_baseline_minus_R160": float(base_abs[core].max()),
            "t12_core_max_abs_matched_minus_R160": float(matched_abs[core].max()),
            "t12_radius_of_baseline_max": float(radii[base_idx]),
            "t12_radius_of_matched_max": float(radii[matched_idx]),
            "t12_max_error_reduction_fraction": (
                float((base_max - matched_max) / base_max) if base_max > 0 else None
            ),
            "t12_core_error_reduction_fraction": (
                float((base_abs[core].max() - matched_abs[core].max()) / base_abs[core].max())
                if base_abs[core].max() > 0 else None
            ),
        }
    return {
        "field_comparison": result,
        "tau": {
            "baseline_minus_R160_final": float(baseline_tau_delta[-1]),
            "matched_minus_R160_final": float(matched_tau_delta[-1]),
            "baseline_max_abs_saved_slice_difference": float(np.max(np.abs(baseline_tau_delta))),
            "matched_max_abs_saved_slice_difference": float(np.max(np.abs(matched_tau_delta))),
        },
    }


def render_report(report, comparison):
    chosen = (
        "geometry.alpha", "geometry.a", "geometry.b", "geometry.X",
        "geometry.Aa", "geometry.K", "geometry.Lambda",
        "scalars.D", "scalars.PD", "scalars.phi", "scalars.PS",
        "matter.dark_matter.momentum_r", "matter.dark_matter.energy_t",
        "matter.baryons.momentum_r", "matter.radiation.momentum_r",
    )
    lines = [
        "# Boundary-matched CMC R=80 experiment — results",
        "",
        f"- Status: **{report['status']}**",
        f"- Commit under test: `{report['code_commit']}`",
        f"- Steps: {report['accepted_steps']}; saved states: {report['sample_count']}; final t: {report['final_t']:.17g}",
        f"- CMC-alpha reference SHA-256: `{report['alpha_reference_sha256']}`",
        "- Matter boundary inputs: archived R=80 last-cell a,b,X,alpha,beta, interpolated in time and frozen for matter RHS evaluation.",
        "- CMC Dirichlet: R=160 alpha(r=79.75), interpolated in time.",
        "- Interpolation caveat: source histories are sampled about every 0.25 time units, much coarser than dt=0.00375; stage-level interpolation error is not rigorously bounded.",
        "",
        "## t=12 mismatch versus archived R=160 state",
        "",
        "| Field | R80 baseline max error | Matched R80 max error | Max-error reduction | Baseline core max | Matched core max |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for field in chosen:
        d = comparison["field_comparison"].get(field)
        if d is None:
            continue
        reduction = d["t12_max_error_reduction_fraction"]
        reduction_text = "n/a" if reduction is None else f"{100.0 * reduction:.2f}%"
        lines.append(
            f"| {field} | {d['t12_max_abs_baseline_minus_R160']:.6g} | "
            f"{d['t12_max_abs_matched_minus_R160']:.6g} | {reduction_text} | "
            f"{d['t12_core_max_abs_baseline_minus_R160']:.6g} | "
            f"{d['t12_core_max_abs_matched_minus_R160']:.6g} |"
        )
    lines += [
        "",
        "## Proper-time accumulator (retain the engine's exact label)",
        "",
        f"- Baseline R80 minus R160 at final saved slice: {comparison['tau']['baseline_minus_R160_final']:.12g}",
        f"- Matched R80 minus R160 at final saved slice: {comparison['tau']['matched_minus_R160_final']:.12g}",
        "",
        "## Interpretation limits",
        "",
        "This is a controlled numerical counterfactual, not proof that a gauge choice is physically irrelevant. Positive error reduction means the boundary-matched trajectory is numerically closer to the archived R=160 trajectory for that field and metric; it does not by itself prove causation. Negative reduction means the mismatch increased. The outer-boundary histories are interpolated from quarter-time snapshots, so small changes must be interpreted against that limitation. No physical front, propagation speed, bounce, or turnaround is inferred here.",
        "",
    ]
    return "\n".join(lines)


def run_experiment(r80_path, r160_path):
    baseline_rows = read_jsonl(r80_path)
    r160_rows = read_jsonl(r160_path)
    if not np.allclose([x["t"] for x in baseline_rows], [x["t"] for x in r160_rows], rtol=0.0, atol=1.0e-12):
        raise ValueError("R80 and R160 reference sample times differ")
    if len(baseline_rows) != len(TARGET_TIMES) + 1:
        raise ValueError("unexpected baseline R80 sample count")
    if int(r160_rows[0]["full_N"]) != 320 or float(r160_rows[0]["full_r_max"]) != 160.0:
        raise ValueError("unexpected R160 reference domain")
    alpha_series = BoundaryAlphaSeries()
    face_series = R80BoundaryFaceSeries(r80_path)
    face_call_log = []
    kernel = make_kernel(alpha_series, face_series, face_call_log)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    profiles = OUT_DIR / "boundary_matched_r80_profiles.jsonl"
    profiles.write_text("")

    state = kernel.initialize(
        resolution=RESOLUTION, r_max=R_MAX, amplitude=AMPLITUDE, width=WIDTH,
        D_amplitude=D_AMPLITUDE, include_radiation=INCLUDE_RADIATION,
    )
    dt_nominal = CFL * float(state.grid.dr)
    if abs(float(state.grid.dr) - 0.5) > 1.0e-13:
        raise ValueError(f"unexpected dr: {state.grid.dr}")
    if abs(float(state.geometry.alpha[-1]) - alpha_series(0.0)) > 1.0e-14:
        raise AssertionError("initial matched CMC alpha does not equal reference series value")

    records = [snapshot_row(state, 0.0, "t=0")]
    with profiles.open("a") as stream:
        stream.write(json.dumps(records[0], sort_keys=True, allow_nan=False) + "\n")
        stream.flush()
    write_progress(state, records, 0, "in_progress")
    target_index, accepted_steps, failure = 0, 0, None

    while state.t < FINAL_TIME - 1.0e-13:
        dt = min(dt_nominal, FINAL_TIME - state.t)
        try:
            next_state = kernel.step(state, dt)
            next_state.geometry.assert_finite_positive()
            if abs(float(next_state.geometry.alpha[-1]) - alpha_series(float(next_state.t))) > 1.0e-14:
                raise AssertionError("accepted CMC boundary does not match reference alpha at accepted time")
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
            print("[BOUNDARY_MATCHED_CMC] HEARTBEAT " + json.dumps({
                "t": float(state.t), "accepted_steps": accepted_steps,
            }, sort_keys=True), flush=True)

        while target_index < len(TARGET_TIMES) and (
            state.t + 0.5 * dt >= TARGET_TIMES[target_index]
            or abs(state.t - TARGET_TIMES[target_index]) < 1.0e-12
        ):
            target = TARGET_TIMES[target_index]
            row = snapshot_row(state, target, f"target={target:g}")
            records.append(row)
            with profiles.open("a") as stream:
                stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
                stream.flush()
            target_index += 1
            write_progress(state, records, accepted_steps, "in_progress", failure)
            print("[BOUNDARY_MATCHED_CMC] CHECKPOINT " + json.dumps({
                "t": float(state.t), "accepted_steps": accepted_steps,
                "sample_count": len(records),
                "center_alpha": float(state.geometry.alpha[0]),
                "outer_alpha": float(state.geometry.alpha[-1]),
                "tau": float(state.tau),
            }, sort_keys=True), flush=True)

    complete = failure is None and state.t >= FINAL_TIME - 1.0e-10 and len(records) == len(TARGET_TIMES) + 1
    status = "completed" if complete else "numerical_or_sampling_failure"
    comparison = None
    if len(records) == len(baseline_rows) == len(r160_rows):
        comparison = compare_records(records, baseline_rows, r160_rows)

    series_path = pathlib.Path(__file__).resolve().parent / "r160_alpha_boundary_r79_75.json"
    report = {
        "schema": "boundary_matched_cmc_r80_report_v1",
        "status": status,
        "diagnostic_only": True,
        "code_commit": os.environ.get("GITHUB_SHA", "local"),
        "accepted_steps": accepted_steps,
        "sample_count": len(records),
        "final_t": float(state.t),
        "final_tau": float(state.tau),
        "final_e_folds": float(state.e_folds),
        "failure": failure,
        "settings": {
            "N": RESOLUTION, "r_max": R_MAX, "dr": float(state.grid.dr),
            "CFL": CFL, "final_time": FINAL_TIME, "sample_interval": SAMPLE_INTERVAL,
            "amplitude": AMPLITUDE, "width": WIDTH, "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "production_physics_changed": False,
            "CMC_operator_coefficients_changed": False,
            "matter_equations_changed": False,
            "CMC_boundary_override_enabled": True,
            "matter_boundary_face_values_frozen": True,
        },
        "alpha_reference_sha256": hashlib.sha256(series_path.read_bytes()).hexdigest(),
        "matter_boundary_provider_calls": len(face_call_log),
        "matter_boundary_provider_time_min": float(min(face_call_log)) if face_call_log else None,
        "matter_boundary_provider_time_max": float(max(face_call_log)) if face_call_log else None,
        "comparison": comparison,
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    if comparison is not None:
        (OUT_DIR / "comparison.json").write_text(json.dumps(comparison, indent=2, sort_keys=True, allow_nan=False) + "\n")
        (OUT_DIR / "report.md").write_text(render_report(report, comparison) + "\n")
    write_progress(state, records, accepted_steps, status, failure)
    print("[BOUNDARY_MATCHED_CMC] FINAL " + json.dumps({
        "status": status, "accepted_steps": accepted_steps, "sample_count": len(records),
        "final_t": float(state.t), "final_tau": float(state.tau),
        "matter_boundary_provider_calls": len(face_call_log),
        "report": str(OUT_DIR / "report.json"),
    }, sort_keys=True), flush=True)
    if status != "completed":
        raise RuntimeError("boundary-matched CMC trajectory did not complete cleanly")
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r80-profile", required=True)
    parser.add_argument("--r160-profile", required=True)
    parser.add_argument("--smoke-only", action="store_true")
    args = parser.parse_args()
    if args.smoke_only:
        smoke(args.r80_profile, args.r160_profile)
    else:
        run_experiment(args.r80_profile, args.r160_profile)


if __name__ == "__main__":
    main()
