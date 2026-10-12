#!/usr/bin/env python3
"""Replay an archived accepted radiation-failure state at four dt values.

This is a bounded one-step diagnostic. Every arm starts from the same exact
stored numerical arrays. It does not run the long lead-in trajectory and does
not alter the equations, cone condition, boundary law, or retry policy.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from engine import matter_rhs as matter_rhs_module
from engine import matter_system as matter_system_module
from engine import production_kernel as production_kernel_module
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.matter_system import Species

SOURCE_BASE = "cdcec8282bff0443f2fe4674b15929822ff38736"
VENDOR_PIN = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_STATE_SHA256 = "ccd1fa0fbabe65c51c3dd7ac2c9227d6d25aa680c169eef1afd3fca4326294cf"
DT_ARMS = (0.03, 0.015, 0.0075, 0.00375)
OUTER_CELLS = 5
SPECIES_FIELDS = ("rest", "energy_t", "momentum_r")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
GEOMETRY_FIELDS = ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def git_value(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], text=True).strip()
    except Exception:
        return "unavailable"


def load_state(npz_path: Path, include_radiation: bool = True):
    data = np.load(npz_path, allow_pickle=False)
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=80, r_max=80.0, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=include_radiation,
    )
    if not np.array_equal(np.asarray(state.grid.centers), data["accepted_r"]):
        raise RuntimeError("frozen-state radial grid differs from the expected N=80/rmax=80 grid")
    for name in GEOMETRY_FIELDS:
        key = f"accepted_geometry_{name}"
        if key not in data:
            raise RuntimeError(f"accepted-state snapshot lacks {key}")
        setattr(state.geometry, name, np.asarray(data[key], dtype=float).copy())
    for name in SCALAR_FIELDS:
        setattr(
            state.scalars, name,
            np.asarray(data[f"accepted_scalar_{name}"], dtype=float).copy(),
        )
    for species_name in ("dark_matter", "baryons", "radiation"):
        species = getattr(state.matter, species_name)
        for field in SPECIES_FIELDS:
            setattr(
                species, field,
                np.asarray(data[f"accepted_matter_{species_name}_{field}"], dtype=float).copy(),
            )
    state.t = float(data["accepted_t"])
    state.tau = float(data["accepted_tau"])
    state.e_folds = float(data["accepted_e_folds"])
    if abs(state.t - 47.130000000001026) > 1.0e-10:
        raise RuntimeError(f"unexpected frozen-state time {state.t:.17g}")
    return kernel, state


def cone(geom, rad, grid):
    metrics = matter_system_module._metric_arrays(metric_slice_from_q(grid, geom))
    gamma_inv = np.asarray([m.gamma_rr_inv for m in metrics], dtype=float)
    sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics], dtype=float)
    ue = np.asarray(rad.energy_t, dtype=float)
    ur = np.asarray(rad.momentum_r, dtype=float)
    C = ue - np.sqrt(gamma_inv) * np.abs(ur)
    ratio = np.full_like(ue, np.nan)
    good = ue > 0
    ratio[good] = np.sqrt(gamma_inv[good]) * np.abs(ur[good]) / ue[good]
    return {
        "gamma_rr_inv": gamma_inv,
        "sqrt_gamma": sqrt_gamma,
        "U_E": ue,
        "U_r": ur,
        "C": C,
        "ratio": ratio,
    }


def max_abs(value):
    array = np.asarray(value, dtype=float)
    finite = np.abs(array[np.isfinite(array)])
    return float(np.max(finite)) if finite.size else None


def geometry_only_witnesses(grid, geom):
    """Raw vacuum/geometric constraints plus Misner–Sharp geometry; no matter recovery."""
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(grid, geom)
    raw_summary = {}
    for key, value in raw.items():
        maximum = max_abs(value)
        if maximum is not None:
            raw_summary[key] = maximum

    r = np.asarray(grid.centers, dtype=float)
    areal_radius = r * np.sqrt(np.asarray(geom.b)) / np.asarray(geom.X)
    dR = grid.cell_derivative_fourth(areal_radius, parity=1)
    Ktheta = np.asarray(geom.K) / 3.0 - np.asarray(geom.Aa) / 2.0
    normal_dR = -areal_radius * Ktheta
    chi = production_kernel_module.trapping_indicator_from_areal_radius(
        dR, normal_dR, np.asarray(geom.X) ** 2 / np.asarray(geom.a)
    )
    mass = production_kernel_module.misner_sharp_mass_from_chi(areal_radius, chi)
    return {
        "raw_geometric_constraint_maxima": raw_summary,
        "areal_radius_outer": float(areal_radius[-1]),
        "chi_outer": float(chi[-1]),
        "M_MS_outer": float(mass[-1]),
        "M_MS_max_abs": max_abs(mass),
    }


def clone_scalar_rhs(srhs):
    return ScalarFields(*(
        np.asarray(getattr(srhs, name), dtype=float).copy()
        for name in SCALAR_FIELDS
    ))


def scalar_add(base, rates, dt):
    return ScalarFields(*(
        np.asarray(getattr(base, name), dtype=float)
        + float(dt) * np.asarray(getattr(rates, name), dtype=float)
        for name in SCALAR_FIELDS
    ))


def capture_rhs_and_predictor(kernel, base_state, dt):
    """Run one real step with transparent observers; never intercepts a failure."""
    ctx = {
        "rhs_count": 0,
        "geometry_count": 0,
        "sources": {},
        "active_radiation": False,
        "source_index": 0,
        "rhs0_geometry": None,
        "scalar_rhs0": None,
        "matter_rhs0": None,
        "metric_predictor": None,
        "radiation_predictor": None,
        "validation_stage": None,
    }

    original_rhs = kernel._rhs
    original_geometry_terms = adapter.geometry_stage_terms
    original_validator = production_kernel_module.validate_species_admissibility
    original_source = matter_system_module._valencia_source
    original_evolve = matter_rhs_module.evolve_species

    def observed_evolve(*args, **kwargs):
        species = kwargs.get("species")
        if species is None and len(args) >= 4:
            species = args[3]
        is_radiation = species is Species.RADIATION
        if is_radiation:
            ctx["active_radiation"] = True
            ctx["source_index"] = 0
            ctx["sources"] = {}
        try:
            return original_evolve(*args, **kwargs)
        finally:
            if is_radiation:
                ctx["active_radiation"] = False

    def observed_source(*args, **kwargs):
        result = original_source(*args, **kwargs)
        if ctx["active_radiation"]:
            i = int(ctx["source_index"])
            ctx["sources"][i] = (float(result[0]), float(result[1]))
            ctx["source_index"] = i + 1
        return result

    def observed_rhs(rhs_state, *args, **kwargs):
        result = original_rhs(rhs_state, *args, **kwargs)
        index = int(ctx["rhs_count"])
        ctx["rhs_count"] = index + 1
        if index == 0:
            srhs, mrhs, _ = result
            ctx["rhs0_geometry"] = rhs_state.geometry.copy()
            ctx["scalar_rhs0"] = clone_scalar_rhs(srhs)
            ctx["matter_rhs0"] = {
                name: rate.copy() for name, rate in mrhs.items()
            }
            ctx["rhs0_rad"] = rhs_state.matter.radiation.copy()
        return result

    def observed_geometry_terms(*args, **kwargs):
        result = original_geometry_terms(*args, **kwargs)
        index = int(ctx["geometry_count"])
        ctx["geometry_count"] = index + 1
        if index == 0:
            geometry = args[1]
            ctx["geometry_rhs0_input"] = geometry.copy()
            ctx["geometry_rhs0_explicit"] = {
                name: np.asarray(result["explicit"][name], dtype=float).copy()
                for name in ("a", "b", "X")
            }
        return result

    def observed_validator(metric, state, species, *, stage):
        if stage == "explicit-predictor" and ctx["metric_predictor"] is None:
            ctx["metric_predictor"] = {
                name: np.asarray(getattr(metric, name), dtype=float).copy()
                for name in ("r",) + GEOMETRY_FIELDS
            }
            ctx["radiation_predictor"] = state.copy()
            ctx["validation_stage"] = stage
        return original_validator(metric, state, species, stage=stage)

    # Patch only the observer hooks and restore all of them after exactly one attempt.
    kernel._rhs = observed_rhs
    adapter.geometry_stage_terms = observed_geometry_terms
    production_kernel_module.validate_species_admissibility = observed_validator
    matter_system_module._valencia_source = observed_source
    matter_rhs_module.evolve_species = observed_evolve
    start = time.time()
    candidate = None
    failure = None
    try:
        candidate = kernel.step(base_state, float(dt))
    except ValueError as exc:
        message = str(exc)
        if "stage='explicit-predictor'" not in message:
            raise
        failure = {
            "class": type(exc).__name__,
            "message": message,
            "validator_stage": ctx["validation_stage"],
        }
    finally:
        kernel._rhs = original_rhs
        adapter.geometry_stage_terms = original_geometry_terms
        production_kernel_module.validate_species_admissibility = original_validator
        matter_system_module._valencia_source = original_source
        matter_rhs_module.evolve_species = original_evolve

    if ctx["rhs0_geometry"] is None or ctx["matter_rhs0"] is None:
        raise RuntimeError("RHS0 capture missing")
    if ctx["geometry_rhs0_input"] is None or ctx["geometry_rhs0_explicit"] is None:
        raise RuntimeError("actual explicit geometry RHS0 capture missing")
    if ctx["metric_predictor"] is None or ctx["radiation_predictor"] is None:
        raise RuntimeError("the explicit-predictor validator was not reached")
    if len(ctx["sources"]) != base_state.grid.n:
        raise RuntimeError(
            f"radiation source capture incomplete: {len(ctx['sources'])}/{base_state.grid.n}"
        )

    # Reconstruct the full explicit predictor state from the same RHS0 used by the
    # real kernel, then overwrite the captured metric/radiation fields with values
    # observed immediately before the real stage validator.
    predictor_geometry = ctx["rhs0_geometry"].copy()
    for name in GEOMETRY_FIELDS:
        setattr(predictor_geometry, name, ctx["metric_predictor"][name].copy())
    predictor_matter = base_state.matter.copy()
    for species_name in ("dark_matter", "baryons", "radiation"):
        updated = kernel._add_matter(
            getattr(base_state.matter, species_name),
            ctx["matter_rhs0"][species_name],
            float(dt),
        )
        setattr(predictor_matter, species_name, updated)
    predictor_matter.radiation = ctx["radiation_predictor"].copy()
    predictor_scalars = scalar_add(
        base_state.scalars, ctx["scalar_rhs0"], float(dt)
    )

    rad_budget = radiation_margin_budget(
        base_state, predictor_geometry, predictor_matter.radiation,
        ctx["rhs0_rad"], ctx["matter_rhs0"]["radiation"],
        ctx["sources"], float(dt),
    )
    metric_budget = metric_predictor_budget(
        base_state.grid, ctx["rhs0_geometry"], predictor_geometry,
        ctx["geometry_rhs0_explicit"], float(dt),
    )
    state_geom_witness = geometry_only_witnesses(base_state.grid, predictor_geometry)

    # Full sourced constraints are safe to evaluate only at an admissible completed
    # state. For a rejected predictor report raw geometry constraints separately.
    final_diagnostics = None
    if candidate is not None:
        full = kernel.diagnostics(candidate, profiles=False)
        keys = (
            "hamiltonian_full_absmax", "momentum_full_absmax",
            "connection_full_absmax", "determinant_constraint_max",
            "M_MS", "H_eff", "cmc_residual_outer_max",
            "cmc_residual_last_interior", "radiation_outer_ratio",
        )
        final_diagnostics = {
            key: float(full[key])
            for key in keys if key in full and np.isscalar(full[key])
        }

    return {
        "dt": float(dt),
        "dt_over_nominal": float(dt / DT_ARMS[0]),
        "status": "expected_explicit_predictor_rejection" if failure else "completed_one_step",
        "final_time": None if candidate is None else float(candidate.t),
        "failure": failure,
        "predictor_budget": rad_budget,
        "geometry_predictor_budget": metric_budget,
        "predictor_geometry_witnesses": state_geom_witness,
        "completed_state_full_diagnostics": final_diagnostics,
        "predictor_state_reconstruction": {
            "radiation_energy_t_max_abs_difference_from_validator": max_abs(
                predictor_matter.radiation.energy_t - ctx["radiation_predictor"].energy_t
            ),
            "radiation_momentum_r_max_abs_difference_from_validator": max_abs(
                predictor_matter.radiation.momentum_r - ctx["radiation_predictor"].momentum_r
            ),
            "predictor_scalar_fields_reconstructed": True,
            "cycle_history_restored_from_npz": False,
        },
        "rhs0_call_count": int(ctx["rhs_count"]),
        "geometry_rhs0_call_count": int(ctx["geometry_count"]),
        "radiation_source_cells_captured": len(ctx["sources"]),
        "wall_seconds": time.time() - start,
    }


def radiation_margin_budget(base_state, predictor_geometry, predictor_rad,
                            rhs0_rad, radiation_rhs, sources, dt):
    g0 = base_state.geometry
    c0 = cone(g0, rhs0_rad, base_state.grid)
    cp = cone(predictor_geometry, predictor_rad, base_state.grid)
    total_e = np.asarray(radiation_rhs.energy_t, dtype=float)
    total_s = np.asarray(radiation_rhs.momentum_r, dtype=float)
    src_e = np.asarray([sources[i][0] for i in range(base_state.grid.n)], dtype=float)
    src_s = np.asarray([sources[i][1] for i in range(base_state.grid.n)], dtype=float)
    flux_e, flux_s = total_e - src_e, total_s - src_s
    ue0, ur0 = c0["U_E"], c0["U_r"]
    root0 = np.sqrt(c0["gamma_rr_inv"])
    ue_flux, ur_flux = ue0 + dt * flux_e, ur0 + dt * flux_s
    c_flux = ue_flux - root0 * np.abs(ur_flux)
    ue_source, ur_source = ue_flux + dt * src_e, ur_flux + dt * src_s
    c_source = ue_source - root0 * np.abs(ur_source)
    c_pred_same_conserved = ue_source - np.sqrt(cp["gamma_rr_inv"]) * np.abs(ur_source)
    c_pred_actual = cp["C"]
    delta_flux = c_flux - c0["C"]
    delta_source = c_source - c_flux
    delta_metric = c_pred_same_conserved - c_source
    # The production predictor's conservative state and the recomputed RHS split
    # may differ by roundoff because species_rhs currently uses unit-step subtract-back.
    reconstruction_e = predictor_rad.energy_t - (ue0 + dt * total_e)
    reconstruction_s = predictor_rad.momentum_r - (ur0 + dt * total_s)
    n = base_state.grid.n
    lo = max(0, n - OUTER_CELLS)
    rows = []
    for i in range(lo, n):
        rows.append({
            "i": int(i), "r": float(base_state.grid.centers[i]),
            "C_accepted_metric": float(c0["C"][i]),
            "C_after_flux_on_accepted_metric": float(c_flux[i]),
            "C_after_sources_on_accepted_metric": float(c_source[i]),
            "C_same_conserved_on_predictor_metric": float(c_pred_same_conserved[i]),
            "C_actual_predictor": float(c_pred_actual[i]),
            "delta_C_flux": float(delta_flux[i]),
            "delta_C_sources": float(delta_source[i]),
            "delta_C_metric_fixed_conserved": float(delta_metric[i]),
            "delta_C_net": float(c_pred_actual[i] - c0["C"][i]),
            "closure_error_split": float(
                delta_flux[i] + delta_source[i] + (c_pred_actual[i]-c_source[i])
                - (c_pred_actual[i]-c0["C"][i])
            ),
            "U_E_accepted": float(ue0[i]), "U_r_accepted": float(ur0[i]),
            "radiation_source_energy_rate": float(src_e[i]),
            "radiation_source_momentum_rate": float(src_s[i]),
            "radiation_flux_energy_rate": float(flux_e[i]),
            "radiation_flux_momentum_rate": float(flux_s[i]),
            "reconstruction_error_U_E": float(reconstruction_e[i]),
            "reconstruction_error_U_r": float(reconstruction_s[i]),
            "gamma_rr_inv_accepted": float(c0["gamma_rr_inv"][i]),
            "gamma_rr_inv_predictor": float(cp["gamma_rr_inv"][i]),
            "a_accepted": float(g0.a[i]), "a_predictor": float(predictor_geometry.a[i]),
            "X_accepted": float(g0.X[i]), "X_predictor": float(predictor_geometry.X[i]),
            "ratio_predictor": float(cp["ratio"][i]),
        })
    finite_ratio = cp["ratio"][np.isfinite(cp["ratio"])]
    max_i = int(np.nanargmax(cp["ratio"]))
    return {
        "budget_order": [
            "accepted radiation on accepted metric",
            "flux-only conserved increment evaluated on accepted metric",
            "source increment evaluated on accepted metric",
            "same reconstructed conservative state evaluated on predictor metric",
            "actual predictor radiation on predictor metric",
        ],
        "outer_cells": rows,
        "whole_grid_min_C_predictor": float(np.min(cp["C"])),
        "whole_grid_min_C_index": int(np.argmin(cp["C"])),
        "whole_grid_max_ratio_predictor": float(np.max(finite_ratio)) if finite_ratio.size else None,
        "whole_grid_max_ratio_index": max_i,
        "source_cells_missing": 0,
        "max_abs_reconstruction_error_U_E": max_abs(reconstruction_e),
        "max_abs_reconstruction_error_U_r": max_abs(reconstruction_s),
        "max_abs_split_closure_error": max(abs(row["closure_error_split"]) for row in rows),
    }


def metric_predictor_budget(grid, accepted_geometry, predictor_geometry,
                            explicit_rhs, dt):
    a0, b0, X0 = map(
        lambda name: np.asarray(getattr(accepted_geometry, name), dtype=float),
        ("a", "b", "X"),
    )
    ap, bp, Xp = map(
        lambda name: np.asarray(getattr(predictor_geometry, name), dtype=float),
        ("a", "b", "X"),
    )
    rhs = {name: np.asarray(explicit_rhs[name], dtype=float) for name in ("a", "b", "X")}
    a_raw, b_raw, X_raw = a0 + dt * rhs["a"], b0 + dt * rhs["b"], X0 + dt * rhs["X"]
    gamma0, gamma_a = X0**2 / a0, X0**2 / a_raw
    gamma_raw, gamma_pred = X_raw**2 / a_raw, Xp**2 / ap
    delta_a, delta_X, delta_projection = gamma_a-gamma0, gamma_raw-gamma_a, gamma_pred-gamma_raw
    explicit_a_formula = -2.0 * np.asarray(accepted_geometry.alpha) * a0 * np.asarray(accepted_geometry.Aa)
    explicit_X_formula = (
        np.asarray(accepted_geometry.alpha) * X0 * np.asarray(accepted_geometry.K) / 3.0
    )
    lo = max(0, grid.n - OUTER_CELLS)
    rows = []
    for i in range(lo, grid.n):
        rows.append({
            "i": int(i), "r": float(grid.centers[i]),
            "a_accepted": float(a0[i]), "a_rhs_euler": float(a_raw[i]),
            "a_predictor": float(ap[i]),
            "b_accepted": float(b0[i]), "b_rhs_euler": float(b_raw[i]),
            "b_predictor": float(bp[i]),
            "X_accepted": float(X0[i]), "X_rhs_euler": float(X_raw[i]),
            "X_predictor": float(Xp[i]),
            "gamma_rr_inv_accepted": float(gamma0[i]),
            "gamma_rr_inv_after_a": float(gamma_a[i]),
            "gamma_rr_inv_after_X": float(gamma_raw[i]),
            "gamma_rr_inv_predictor": float(gamma_pred[i]),
            "delta_gamma_rr_inv_a": float(delta_a[i]),
            "delta_gamma_rr_inv_X": float(delta_X[i]),
            "delta_gamma_rr_inv_projection_or_other": float(delta_projection[i]),
            "delta_gamma_rr_inv_total": float(gamma_pred[i]-gamma0[i]),
            "gamma_identity_closure": float(
                delta_a[i] + delta_X[i] + delta_projection[i] - (gamma_pred[i]-gamma0[i])
            ),
            "rhs_a_identity_residual": float(rhs["a"][i]-explicit_a_formula[i]),
            "rhs_X_identity_residual": float(rhs["X"][i]-explicit_X_formula[i]),
        })
    return {
        "outer_cells": rows,
        "max_abs_gamma_identity_closure": max(abs(r["gamma_identity_closure"]) for r in rows),
        "max_abs_rhs_a_identity_residual": max(abs(r["rhs_a_identity_residual"]) for r in rows),
        "max_abs_rhs_X_identity_residual": max(abs(r["rhs_X_identity_residual"]) for r in rows),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("runs/source-artifact"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/radiation-frozen-state-stage-budget"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    state_path = args.input_dir / "boundary_failure_accepted_state.npz"
    if not state_path.is_file():
        raise FileNotFoundError(f"missing frozen-state artifact: {state_path}")
    state_hash = sha256(state_path)
    if state_hash != EXPECTED_STATE_SHA256:
        raise RuntimeError(f"frozen-state SHA-256 mismatch: {state_hash}")
    metadata_path = args.input_dir / "boundary_failure_accepted_state.json"
    metadata = json.loads(metadata_path.read_text())
    if metadata.get("rate_fix_source_commit") != SOURCE_BASE:
        raise RuntimeError("frozen-state source commit does not match the pinned stacked source")
    if metadata.get("pinned_vendor_commit") != VENDOR_PIN:
        raise RuntimeError("frozen-state vendor pin mismatch")

    started = time.time()
    kernel, seed = load_state(state_path)
    baseline_full = kernel.diagnostics(seed, profiles=False)
    baseline_keys = (
        "hamiltonian_full_absmax", "momentum_full_absmax",
        "connection_full_absmax", "determinant_constraint_max",
        "M_MS", "H_eff", "cmc_residual_outer_max",
        "cmc_residual_last_interior", "radiation_outer_ratio",
    )
    baseline = {key: float(baseline_full[key]) for key in baseline_keys
                if key in baseline_full and np.isscalar(baseline_full[key])}

    arms = []
    for dt in DT_ARMS:
        result = capture_rhs_and_predictor(
            kernel, copy.deepcopy(seed), dt
        )
        arms.append(result)
        print(json.dumps({
            "arm_dt": dt,
            "status": result["status"],
            "max_ratio": result["predictor_budget"]["whole_grid_max_ratio_predictor"],
            "min_C": result["predictor_budget"]["whole_grid_min_C_predictor"],
            "constraint_witnesses": result["predictor_geometry_witnesses"]["raw_geometric_constraint_maxima"],
        }), flush=True)

    report = {
        "schema": "radiation_frozen_state_dt_stage_budget_v1",
        "status": "all_four_arms_captured",
        "diagnostic_only": True,
        "provenance": {
            "repository": "miahpenn/0-Cosmos-engine",
            "source_base": SOURCE_BASE,
            "instrumentation_branch_head": git_value("rev-parse", "HEAD"),
            "vendor_pin": git_value("rev-parse", "HEAD:vendor/bb-palatini-unified-r0"),
            "expected_vendor_pin": VENDOR_PIN,
            "frozen_state_npz_sha256": state_hash,
            "frozen_state_metadata_sha256": sha256(metadata_path),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "source_scope": (
                "Only engine/radiation_geometry_stage_trace.py may differ from the pinned "
                "stacked source under this gate; this standalone audit does not alter equations."
            ),
        },
        "frozen_input": {
            "t": float(seed.t), "tau": float(seed.tau),
            "dt_nominal_rejected": 0.03, "N": int(seed.grid.n),
            "r_max": float(seed.grid.r_max),
            "exact_state_arrays_restored": True,
            "cycle_history_restored_from_npz": False,
            "initial_full_sourced_constraint_witnesses": baseline,
            "expected_rejection": metadata.get("failure", {}),
        },
        "arms": arms,
        "interpretation_limits": [
            "Every arm starts from the same exact archived numerical arrays at t=47.13; this does not replay the preceding 1571-step lead-in.",
            "The nominal dt arm is expected to reject at the existing explicit-predictor guard. Smaller dt arms are independent one-step candidates, not a multi-step adaptive evolution.",
            "Raw geometric constraints on the explicit predictor are computed without recovering matter primitives. Full sourced H/M and completed-state diagnostics are available only for admissible completed arms.",
            "The accepted-state JSON does not include the cycle/history ledger objects; those are freshly initialized. They do not feed the explicit predictor equations or the constraint/metric-margin budgets reported here.",
            "The unit-step/subtract-back RHS contract may leave tiny finite precision reconstruction differences; these are measured rather than suppressed.",
            "No clipping, floor, damping, source change, accepted-metric recovery, cone relaxation, or retry policy is used.",
        ],
        "elapsed_seconds": time.time() - started,
    }
    report_path = args.output_dir / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False))
    digest = sha256(report_path)
    (args.output_dir / "SHA256SUMS.txt").write_text(f"{digest}  report.json\n")
    print(json.dumps({
        "status": report["status"],
        "report": str(report_path),
        "report_sha256": digest,
        "elapsed_seconds": time.time() - started,
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
