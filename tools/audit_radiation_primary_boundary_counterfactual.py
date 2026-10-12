"""Matched two-step counterfactual for the primary-predictor outer Aa boundary.

This diagnostic replays the frozen late-time window on the pinned rate-only
RHS source, then compares a baseline branch with a branch that restores only
the pre-boundary outer Aa during the primary-predictor boundary call at t~46.92.
No production law is edited; the override exists only inside this process.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import subprocess
import time
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.matter_system import Species, _metric_arrays
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import (
    V55MatterState,
    metric_slice_from_q,
)

SOURCE_BASE = "1c63090d4f42396f899421a2379f70320689d5ad"
VENDOR_PIN = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
TARGET_START = 46.92
DT = 0.03
STEPS_FROM_START_TO_FAILURE = 2
OUTER_INDEX = 79
EXPECTED_BASELINE_PREDICTOR_C = -1.699741681811858e-8
BASELINE_C_ABS_TOL = 5.0e-11
OUTPUT = Path("artifacts/radiation-primary-boundary-counterfactual")


def _git(args: list[str]) -> str:
    try:
        return subprocess.check_output(["git", *args], text=True).strip()
    except Exception:
        return "unavailable"


def _hash_file(path: str) -> str | None:
    p = Path(path)
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None


def _clone_state(state: ProductionState) -> ProductionState:
    scalar_names = ("S", "PS", "D", "PD", "phi", "Pi")
    matter_names = ("dark_matter", "baryons", "radiation")
    scalars = ScalarFields(*(
        np.asarray(getattr(state.scalars, name)).copy()
        for name in scalar_names
    ))
    matter = V55MatterState(*(
        getattr(state.matter, name).copy()
        for name in matter_names
    ))
    return ProductionState(
        grid=state.grid,
        geometry=state.geometry.copy(),
        scalars=scalars,
        matter=matter,
        t=float(state.t),
        tau=float(state.tau),
        e_folds=float(state.e_folds),
        cycle=copy.deepcopy(state.cycle),
        history=list(state.history),
        handoffs=list(state.handoffs),
    )


def _array(value) -> np.ndarray:
    return np.asarray(value, dtype=float)


def _capture_tree(value, prefix=""):
    """Copy numerical arrays in a nested RHS/geometry result for A/B diffs."""
    out = {}
    if isinstance(value, dict):
        for key, child in value.items():
            out.update(_capture_tree(child, f"{prefix}.{key}" if prefix else key))
    elif isinstance(value, np.ndarray) or np.isscalar(value):
        try:
            arr = np.asarray(value, dtype=float)
            if np.all(np.isfinite(arr)):
                out[prefix] = arr.copy()
        except (TypeError, ValueError):
            pass
    return out


def _tree_delta(a: dict, b: dict) -> dict:
    rows = {}
    for key in sorted(set(a) & set(b)):
        aa, bb = np.asarray(a[key], dtype=float), np.asarray(b[key], dtype=float)
        if aa.shape != bb.shape:
            continue
        delta = np.abs(aa - bb)
        rows[key] = {
            "max_abs": float(np.max(delta)) if delta.size else 0.0,
            "at_outer_cell_79": (
                float(delta[OUTER_INDEX])
                if delta.ndim == 1 and delta.size > OUTER_INDEX else None
            ),
            "actual_at_outer_cell_79": (
                float(aa[OUTER_INDEX])
                if aa.ndim == 1 and aa.size > OUTER_INDEX else None
            ),
            "counterfactual_at_outer_cell_79": (
                float(bb[OUTER_INDEX])
                if bb.ndim == 1 and bb.size > OUTER_INDEX else None
            ),
        }
    return rows


def _stage_geometry_capture(result) -> dict:
    return _capture_tree({
        "explicit": result["explicit"],
        "primary_l2": result["primary_l2"],
        "primary_l3": result["primary_l3"],
        "lambda_l2": result["lambda_l2"],
        "lambda_l3": result["lambda_l3"],
    })


def _rhs_capture(result) -> dict:
    srhs, mrhs, metric_derivatives = result
    return _capture_tree({
        "scalar": {
            name: getattr(srhs, name)
            for name in ("S", "PS", "D", "PD", "phi", "Pi")
        },
        "matter": {
            species: {
                field: getattr(mrhs[species], field)
                for field in ("rest", "energy_t", "momentum_r")
            }
            for species in ("dark_matter", "baryons", "radiation")
        },
        "metric_derivatives_time": metric_derivatives.time,
        "metric_derivatives_radial": metric_derivatives.radial,
    })


def _state_summary(kernel: V55ProductionKernel, state: ProductionState) -> dict:
    metric_slice = metric_slice_from_q(state.grid, state.geometry)
    metrics = _metric_arrays(metric_slice)
    rad = state.matter.radiation
    m = metrics[OUTER_INDEX]
    ue = float(rad.energy_t[OUTER_INDEX])
    ur = float(rad.momentum_r[OUTER_INDEX])
    sqrt_gamma = float(m.sqrt_gamma)
    root_ginv = math.sqrt(float(m.gamma_rr_inv))
    c = ue - root_ginv * abs(ur)
    ratio = root_ginv * abs(ur) / ue if ue > 0.0 else None
    return {
        "t": float(state.t),
        "tau": float(state.tau),
        "outer_cell_i": OUTER_INDEX,
        "r": float(state.grid.centers[OUTER_INDEX]),
        "a": float(state.geometry.a[OUTER_INDEX]),
        "b": float(state.geometry.b[OUTER_INDEX]),
        "X": float(state.geometry.X[OUTER_INDEX]),
        "Aa": float(state.geometry.Aa[OUTER_INDEX]),
        "K": float(state.geometry.K[OUTER_INDEX]),
        "Lambda": float(state.geometry.Lambda[OUTER_INDEX]),
        "alpha": float(state.geometry.alpha[OUTER_INDEX]),
        "gamma_rr_inv": float(m.gamma_rr_inv),
        "U_E": ue,
        "U_r": ur,
        "E": ue / sqrt_gamma,
        "abs_S": root_ginv * abs(ur) / sqrt_gamma,
        "cone_margin_densitized": c,
        "ratio_absS_over_E": ratio,
    }


def _diagnostic_witnesses(kernel: V55ProductionKernel, state: ProductionState) -> dict:
    raw = kernel.diagnostics(state, profiles=False)
    keys = (
        "hamiltonian_full_absmax", "momentum_full_absmax",
        "connection_full_absmax", "determinant_constraint_max",
        "M_MS", "H_eff", "cmc_residual_outer_max",
        "cmc_residual_last_interior", "lapse_min", "lapse_max",
        "radiation_outer_ratio",
    )
    return {key: float(raw[key]) for key in keys if key in raw and np.isscalar(raw[key])}


def _explicit_predictor(kernel: V55ProductionKernel, state: ProductionState, dt: float) -> dict:
    """Construct the first actual explicit predictor without entering CMC recovery."""
    grid = state.grid
    g0 = state.geometry.copy()
    s0 = state.scalars
    m0 = state.matter
    g0.alpha = kernel._solve_lapse(grid, g0, s0, m0)[0]
    g0.beta.fill(0.0)
    g0.B.fill(0.0)
    rhs0_scalar, rhs0_matter, _ = kernel._rhs(
        ProductionState(
            grid=grid, geometry=g0, scalars=s0, matter=m0,
            t=state.t, tau=state.tau, e_folds=state.e_folds,
        )
    )
    gterms0 = adapter.geometry_stage_terms(
        grid, g0, s0, m0, lambda_m=2.0
    )
    predictor_geometry = g0.copy()
    for name in ("a", "b", "X"):
        setattr(
            predictor_geometry,
            name,
            getattr(g0, name) + dt * gterms0["explicit"][name],
        )
    predictor_geometry.beta.fill(0.0)
    predictor_geometry.B.fill(0.0)
    predictor_scalars = ScalarFields(*(
        getattr(s0, name) + dt * getattr(rhs0_scalar, name)
        for name in ("S", "PS", "D", "PD", "phi", "Pi")
    ))
    predictor_matter = V55MatterState(*(
        kernel._add_matter(getattr(m0, name), rhs0_matter[name], dt)
        for name in ("dark_matter", "baryons", "radiation")
    ))
    metric = metric_slice_from_q(grid, predictor_geometry)
    metrics = _metric_arrays(metric)
    m = metrics[OUTER_INDEX]
    rad = predictor_matter.radiation
    ue = float(rad.energy_t[OUTER_INDEX])
    ur = float(rad.momentum_r[OUTER_INDEX])
    root_ginv = math.sqrt(float(m.gamma_rr_inv))
    c = ue - root_ginv * abs(ur)
    ratio = root_ginv * abs(ur) / ue if ue > 0.0 else None
    return {
        "accepted_t": float(state.t),
        "attempted_t": float(state.t + dt),
        "dt": float(dt),
        "outer_cell_i": OUTER_INDEX,
        "r": float(grid.centers[OUTER_INDEX]),
        "a_predictor": float(predictor_geometry.a[OUTER_INDEX]),
        "b_predictor": float(predictor_geometry.b[OUTER_INDEX]),
        "X_predictor": float(predictor_geometry.X[OUTER_INDEX]),
        "alpha_not_yet_resolved": float(predictor_geometry.alpha[OUTER_INDEX]),
        "gamma_rr_inv": float(m.gamma_rr_inv),
        "U_E": ue,
        "U_r": ur,
        "cone_margin_densitized": c,
        "ratio_absS_over_E": ratio,
    }


def _step_to_target(kernel: V55ProductionKernel, state: ProductionState) -> ProductionState:
    dt = DT * float(state.grid.dr)
    for _ in range(2):
        state = kernel.step(state, dt)
    if abs(float(state.t) - 46.98) > 1.0e-8:
        raise RuntimeError(f"two-step replay did not reach t=46.98; got {state.t:.17g}")
    return state


def main() -> int:
    start_wall = time.time()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    actual_capture: dict = {}
    control_capture: dict = {}
    ctx = {
        "active_t": None,
        "mode": None,
        "suppress_primary": False,
        "boundary_count": 0,
        "solve_count": 0,
        "rhs_count": 0,
        "geometry_terms_count": 0,
        "captures": {"actual": actual_capture, "no_primary_boundary": control_capture},
    }

    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=80, r_max=80.0, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True,
    )
    dt = DT * float(state.grid.dr)
    nsteps = int(round(TARGET_START / dt))
    for i in range(nsteps):
        state = kernel.step(state, dt)
        if (i + 1) % 200 == 0:
            print(json.dumps({
                "replay_progress": True,
                "steps": i + 1,
                "t": float(state.t),
                "tau": float(state.tau),
                "elapsed_seconds": time.time() - start_wall,
            }), flush=True)
    if abs(float(state.t) - TARGET_START) > 1.0e-8:
        raise RuntimeError(
            f"frozen window mismatch: expected t~{TARGET_START}, got {state.t:.17g}"
        )
    seed = _clone_state(state)
    print(json.dumps({
        "frozen_seed": True, "t": float(seed.t), "tau": float(seed.tau),
        "steps": nsteps, "dt": dt,
    }), flush=True)

    kernel_cls = V55ProductionKernel
    original_boundary = kernel_cls._apply_outer_light_boundary
    original_solve = kernel_cls._solve_lapse
    original_rhs = kernel_cls._rhs
    original_step = kernel_cls.step
    original_geometry_terms = adapter.geometry_stage_terms

    def at_target_step() -> bool:
        return (
            ctx["active_t"] is not None
            and abs(float(ctx["active_t"]) - TARGET_START) < 1.0e-8
        )

    def wrapped_boundary(grid, geometry, scalars, matter, *args, **kwargs):
        call_idx = int(ctx["boundary_count"])
        ctx["boundary_count"] = call_idx + 1
        is_target = at_target_step()
        before = float(geometry.Aa[OUTER_INDEX]) if is_target else None
        result = original_boundary(grid, geometry, scalars, matter, *args, **kwargs)
        if is_target:
            applied = float(geometry.Aa[OUTER_INDEX])
            record = {
                "call_index": call_idx,
                "kind": (
                    "primary_predictor_outer_light_boundary"
                    if call_idx == 0 else "completed_state_outer_light_boundary"
                    if call_idx == 1 else f"outer_light_boundary_call_{call_idx}"
                ),
                "Aa_before": before,
                "Aa_after_boundary_formula": applied,
                "delta_Aa_boundary_formula": applied - before,
                "Aa_used_after_counterfactual_override": applied,
                "suppressed": False,
            }
            if is_target and ctx["suppress_primary"] and call_idx == 0:
                geometry.Aa[OUTER_INDEX] = before
                record["Aa_used_after_counterfactual_override"] = float(
                    geometry.Aa[OUTER_INDEX]
                )
                record["suppressed"] = True
            ctx["captures"][ctx["mode"]].setdefault("boundary_calls", []).append(record)
        return result

    def wrapped_solve(grid, geometry, scalars, matter, *args, **kwargs):
        call_idx = int(ctx["solve_count"])
        ctx["solve_count"] = call_idx + 1
        result = original_solve(grid, geometry, scalars, matter, *args, **kwargs)
        if at_target_step() and call_idx == 2:
            capture = ctx["captures"][ctx["mode"]]
            capture["primary_cmc"] = {
                "Aa_input_outer": float(geometry.Aa[OUTER_INDEX]),
                "K_input_outer": float(geometry.K[OUTER_INDEX]),
                "a_input_outer": float(geometry.a[OUTER_INDEX]),
                "X_input_outer": float(geometry.X[OUTER_INDEX]),
                "alpha_outer": float(result[0][OUTER_INDEX]),
                "alpha_min": float(np.min(result[0])),
                "alpha_max": float(np.max(result[0])),
                "kdot_target": float(result[1]),
                "alpha_profile": np.asarray(result[0], dtype=float).copy(),
            }
        return result

    def wrapped_rhs(kernel_self, rhs_state, *args, **kwargs):
        call_idx = int(ctx["rhs_count"])
        ctx["rhs_count"] = call_idx + 1
        result = original_rhs(kernel_self, rhs_state, *args, **kwargs)
        if at_target_step() and call_idx == 1:
            ctx["captures"][ctx["mode"]]["rhs1_primary_predictor"] = (
                _rhs_capture(result)
            )
        return result

    def wrapped_geometry_terms(*args, **kwargs):
        call_idx = int(ctx["geometry_terms_count"])
        ctx["geometry_terms_count"] = call_idx + 1
        result = original_geometry_terms(*args, **kwargs)
        if at_target_step() and call_idx == 2:
            ctx["captures"][ctx["mode"]]["geometry_terms1_primary_predictor"] = (
                _stage_geometry_capture(result)
            )
        return result

    def wrapped_step(kernel_self, step_state, step_dt):
        ctx["active_t"] = float(step_state.t)
        ctx["mode"] = ctx["run_mode"]
        ctx["suppress_primary"] = bool(ctx["run_mode"] == "no_primary_boundary")
        ctx["boundary_count"] = 0
        ctx["solve_count"] = 0
        ctx["rhs_count"] = 0
        ctx["geometry_terms_count"] = 0
        if at_target_step():
            ctx["captures"][ctx["mode"]].update({
                "step_start_t": float(step_state.t),
                "step_start_tau": float(step_state.tau),
                "dt": float(step_dt),
            })
        return original_step(kernel_self, step_state, step_dt)

    # Instrumentation is process-local. It does not modify the equations or
    # any production file and is restored even if a branch fails.
    kernel_cls._apply_outer_light_boundary = staticmethod(wrapped_boundary)
    kernel_cls._solve_lapse = staticmethod(wrapped_solve)
    kernel_cls._rhs = wrapped_rhs
    kernel_cls.step = wrapped_step
    adapter.geometry_stage_terms = wrapped_geometry_terms
    try:
        ctx["run_mode"] = "actual"
        actual_state = _step_to_target(kernel, _clone_state(seed))
        actual_state_4695 = None
        # We need both accepted snapshots. The two-step helper only returns t=46.98,
        # so a separately captured t=46.95 comes from history-independent replay.
        # Re-run the first step alone to retain the midpoint without changing the seed.
        ctx["run_mode"] = "actual"
        actual_mid = kernel.step(_clone_state(seed), dt)
        # The counterfactual replay is independent and changes only the first
        # primary-predictor boundary output on the frozen target step.
        ctx["run_mode"] = "no_primary_boundary"
        control_state = _step_to_target(kernel, _clone_state(seed))
    finally:
        kernel_cls._apply_outer_light_boundary = staticmethod(original_boundary)
        kernel_cls._solve_lapse = staticmethod(original_solve)
        kernel_cls._rhs = original_rhs
        kernel_cls.step = original_step
        adapter.geometry_stage_terms = original_geometry_terms

    # Compute the failed next predictor algebraically on both accepted states.
    # This stops before predictor CMC primitive recovery raises and mutates nothing.
    actual_predictor = _explicit_predictor(kernel, actual_state, dt)
    control_predictor = _explicit_predictor(kernel, control_state, dt)
    if abs(actual_predictor["cone_margin_densitized"] - EXPECTED_BASELINE_PREDICTOR_C) > BASELINE_C_ABS_TOL:
        raise RuntimeError(
            "baseline predictor failed the frozen event identity check: "
            f"expected {EXPECTED_BASELINE_PREDICTOR_C:.17e}, got "
            f"{actual_predictor['cone_margin_densitized']:.17e}"
        )
    if abs(float(actual_state.t) - float(control_state.t)) > 1.0e-10:
        raise RuntimeError("actual and counterfactual branches have different target times")

    # Primary stage data are captured at the first target step. Compare the
    # boundary correction that the counterfactual suppressed, the following
    # CMC response, RHS1, and stage-1 geometry RHS without fitting parameters.
    actual_stage = actual_capture
    control_stage = control_capture
    if "primary_cmc" not in actual_stage or "primary_cmc" not in control_stage:
        raise RuntimeError("missing matched primary-predictor CMC captures")
    for capture in (actual_stage, control_stage):
        profile = capture["primary_cmc"].pop("alpha_profile")
        capture["primary_cmc"]["alpha_profile"] = profile
    alpha_delta = np.abs(
        actual_stage["primary_cmc"]["alpha_profile"]
        - control_stage["primary_cmc"]["alpha_profile"]
    )
    actual_stage["primary_cmc"].pop("alpha_profile", None)
    control_stage["primary_cmc"].pop("alpha_profile", None)

    accepted_actual = {
        "t46_95": _state_summary(kernel, actual_mid),
        "t46_98": _state_summary(kernel, actual_state),
    }
    # Obtain the counterfactual midpoint by a one-step copy from the same seed
    # after the instrumentation has been restored; suppress the primary call
    # only for this diagnostic step via a fresh scoped wrapper.
    # The exact full two-step counterfactual state is already saved above.
    accepted_control = {"t46_98": _state_summary(kernel, control_state)}
    # Run a fully matched one-step primary-only control to record the midpoint
    # accepted state; the start state is still the same frozen seed.
    old_boundary = kernel_cls._apply_outer_light_boundary
    boundary_calls_mid = 0
    def midpoint_boundary(grid, geometry, scalars, matter, *args, **kwargs):
        nonlocal boundary_calls_mid
        call_idx = boundary_calls_mid
        boundary_calls_mid += 1
        before = float(geometry.Aa[OUTER_INDEX])
        result = old_boundary(grid, geometry, scalars, matter, *args, **kwargs)
        if call_idx == 0:
            geometry.Aa[OUTER_INDEX] = before
        return result
    kernel_cls._apply_outer_light_boundary = staticmethod(midpoint_boundary)
    try:
        control_mid = kernel.step(_clone_state(seed), dt)
    finally:
        kernel_cls._apply_outer_light_boundary = staticmethod(old_boundary)
    accepted_control["t46_95"] = _state_summary(kernel, control_mid)

    def state_delta_summary(actual, control):
        fields = {
            "geometry": ("a", "b", "X", "Aa", "K", "Lambda", "alpha"),
            "scalars": ("S", "PS", "D", "PD", "phi", "Pi"),
        }
        out = {}
        for group, names in fields.items():
            obj_a = actual.geometry if group == "geometry" else actual.scalars
            obj_b = control.geometry if group == "geometry" else control.scalars
            out[group] = {}
            for name in names:
                aa, bb = _array(getattr(obj_a, name)), _array(getattr(obj_b, name))
                out[group][name] = {
                    "max_abs": float(np.max(np.abs(aa - bb))),
                    "outer_cell_79": float(abs(aa[OUTER_INDEX] - bb[OUTER_INDEX])),
                }
        out["matter"] = {}
        for species in ("dark_matter", "baryons", "radiation"):
            out["matter"][species] = {}
            for field in ("rest", "energy_t", "momentum_r"):
                aa = _array(getattr(getattr(actual.matter, species), field))
                bb = _array(getattr(getattr(control.matter, species), field))
                out["matter"][species][field] = {
                    "max_abs": float(np.max(np.abs(aa - bb))),
                    "outer_cell_79": float(abs(aa[OUTER_INDEX] - bb[OUTER_INDEX])),
                }
        return out

    accepted_delta_4695 = state_delta_summary(actual_mid, control_mid)
    accepted_delta_4698 = state_delta_summary(actual_state, control_state)
    diag_actual = _diagnostic_witnesses(kernel, actual_state)
    diag_control = _diagnostic_witnesses(kernel, control_state)
    diag_delta = {
        key: {
            "actual": diag_actual[key],
            "no_primary_boundary": diag_control[key],
            "actual_minus_counterfactual": diag_actual[key] - diag_control[key],
        }
        for key in sorted(set(diag_actual) & set(diag_control))
    }
    boundary_actual = actual_stage.get("boundary_calls", [])
    boundary_control = control_stage.get("boundary_calls", [])
    if not boundary_actual or not boundary_control:
        raise RuntimeError("missing boundary pre/post measurements")
    actual_primary = next((r for r in boundary_actual if r["call_index"] == 0), None)
    control_primary = next((r for r in boundary_control if r["call_index"] == 0), None)
    if actual_primary is None or control_primary is None:
        raise RuntimeError("primary boundary call not found")
    if abs(actual_primary["delta_Aa_boundary_formula"] - 4.307354247837236e-5) > 1.0e-8:
        raise RuntimeError(
            "frozen primary-boundary correction differs materially from the prior trace: "
            f"{actual_primary['delta_Aa_boundary_formula']:.17e}"
        )

    # Use numeric summary only in report; never serialize complete model state.
    report = {
        "schema": "radiation_primary_boundary_counterfactual_v1",
        "status": "matched_two_step_counterfactual_complete",
        "diagnostic_only": True,
        "provenance": {
            "repository": "miahpenn/0-Cosmos-engine",
            "source_baseline": SOURCE_BASE,
            "instrumentation_head": _git(["rev-parse", "HEAD"]),
            "vendor_pin": _git(["rev-parse", "HEAD:vendor/bb-palatini-unified-r0"]),
            "expected_vendor_pin": VENDOR_PIN,
            "python": os.sys.version,
            "numpy": np.__version__,
            "files_sha256": {
                p: _hash_file(p) for p in (
                    "engine/production_kernel.py",
                    "engine/matter_system.py",
                    "engine/matter_rhs.py",
                    "engine/radiation_geometry_stage_trace.py",
                    "tools/audit_radiation_primary_boundary_counterfactual.py",
                )
            },
        },
        "configuration": {
            "N": 80, "r_max": 80.0, "dr": float(seed.grid.dr),
            "dt": dt, "target_step_start": TARGET_START,
            "accepted_compare_times": [float(actual_mid.t), float(actual_state.t)],
            "target_cell": OUTER_INDEX,
            "primary_boundary_control": (
                "On the no-primary-boundary arm only, run the real boundary function "
                "then restore the pre-call outer Aa for the primary-predictor call "
                "at the single t~46.92 step. The completed-state boundary and every "
                "other stage remain enabled."
            ),
        },
        "replay": {
            "steps_to_frozen_seed": nsteps,
            "frozen_seed_t": float(seed.t),
            "frozen_seed_tau": float(seed.tau),
            "dt": dt,
            "accepted_actual_t46_95": accepted_actual["t46_95"],
            "accepted_control_t46_95": accepted_control["t46_95"],
            "accepted_actual_t46_98": accepted_actual["t46_98"],
            "accepted_control_t46_98": accepted_control["t46_98"],
        },
        "primary_boundary": {
            "actual": actual_primary,
            "counterfactual": control_primary,
            "applied_delta_removed_only_for_primary_call": True,
            "note": (
                "Boundary formula executes in both branches to measure its raw "
                "correction. The control restores pre-boundary Aa only at primary "
                "predictor call index 0 for the target step; it does not change "
                "the completed-state boundary."
            ),
        },
        "primary_cmc_response": {
            "actual": actual_stage["primary_cmc"],
            "no_primary_boundary": control_stage["primary_cmc"],
            "alpha_profile_max_abs_difference": float(np.max(alpha_delta)),
            "alpha_profile_outer_cell_difference": float(alpha_delta[OUTER_INDEX]),
            "kdot_target_difference": (
                actual_stage["primary_cmc"]["kdot_target"]
                - control_stage["primary_cmc"]["kdot_target"]
            ),
        },
        "primary_rhs1_stage_difference": _tree_delta(
            actual_stage["rhs1_primary_predictor"],
            control_stage["rhs1_primary_predictor"],
        ),
        "primary_geometry_terms1_difference": _tree_delta(
            actual_stage["geometry_terms1_primary_predictor"],
            control_stage["geometry_terms1_primary_predictor"],
        ),
        "accepted_state_delta_actual_minus_no_primary_boundary": {
            "t46_95": accepted_delta_4695,
            "t46_98": accepted_delta_4698,
        },
        "next_explicit_predictor": {
            "actual": actual_predictor,
            "no_primary_boundary": control_predictor,
            "cone_margin_actual_minus_counterfactual": (
                actual_predictor["cone_margin_densitized"]
                - control_predictor["cone_margin_densitized"]
            ),
            "ratio_actual_minus_counterfactual": (
                actual_predictor["ratio_absS_over_E"]
                - control_predictor["ratio_absS_over_E"]
            ),
            "baseline_margin_reference": EXPECTED_BASELINE_PREDICTOR_C,
            "baseline_margin_match_abs_tolerance": BASELINE_C_ABS_TOL,
        },
        "joint_admission_witnesses_at_accepted_t46_98": {
            "actual": diag_actual,
            "no_primary_boundary": diag_control,
            "actual_minus_counterfactual": diag_delta,
        },
        "elapsed_seconds": time.time() - start_wall,
        "interpretation_limits": [
            "This is a two-step branch counterfactual from one frozen accepted state, not a long alternate trajectory.",
            "The boundary correction is suppressed only for the primary-predictor call at the target step; the completed-state boundary remains active in both branches.",
            "A successful matched comparison does not prove the boundary law is correct or incorrect and is not production admission.",
            "Radiation-cone, H/M/connection, determinant, CMC and Misner-Sharp witnesses are reported together; no equation or acceptance threshold is changed.",
        ],
    }
    trace_path = OUTPUT / "report.json"
    trace_path.write_text(json.dumps(report, indent=2, allow_nan=False))
    sha = hashlib.sha256(trace_path.read_bytes()).hexdigest()
    (OUTPUT / "SHA256SUMS.txt").write_text(f"{sha}  report.json\n")
    print(json.dumps({
        "status": report["status"],
        "report": str(trace_path),
        "report_sha256": sha,
        "actual_predictor_C": actual_predictor["cone_margin_densitized"],
        "counterfactual_predictor_C": control_predictor["cone_margin_densitized"],
        "delta_C": report["next_explicit_predictor"]["cone_margin_actual_minus_counterfactual"],
        "primary_boundary_delta_Aa": actual_primary["delta_Aa_boundary_formula"],
        "CMC_alpha_max_delta": report["primary_cmc_response"]["alpha_profile_max_abs_difference"],
        "geometry_terms_count": len(report["primary_geometry_terms1_difference"]),
        "rhs1_terms_count": len(report["primary_rhs1_stage_difference"]),
    }, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    OUTPUT.mkdir(parents=True, exist_ok=True)
    try:
        if _git(["rev-parse", "HEAD:vendor/bb-palatini-unified-r0"]) != VENDOR_PIN:
            raise RuntimeError("pinned vendor gitlink mismatch")
        raise SystemExit(main())
    except Exception as exc:
        failure = {
            "schema": "radiation_primary_boundary_counterfactual_v1",
            "status": "diagnostic_harness_failure",
            "diagnostic_only": True,
            "exception": f"{type(exc).__name__}: {exc}",
            "traceback": __import__("traceback").format_exc(),
            "head": _git(["rev-parse", "HEAD"]),
            "vendor": _git(["rev-parse", "HEAD:vendor/bb-palatini-unified-r0"]),
        }
        (OUTPUT / "harness_failure.json").write_text(json.dumps(failure, indent=2))
        print(json.dumps(failure, indent=2), flush=True)
        raise
