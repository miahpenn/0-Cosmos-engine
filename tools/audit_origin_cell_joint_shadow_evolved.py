#!/usr/bin/env python3
"""Full accepted-state diagnostic for a copied regular-origin initializer closure.

No production source is edited. Four arms are evolved through the same pinned
V55ProductionKernel.step(): production baseline, half-cell-integral candidate,
regular-origin CMC-row only, and joint candidate+CMC row. The initializer copy
differs only in the first half-cell integral; all later radial intervals retain
the source-pinned trapezoid recurrence. The candidate's matter is reinitialized
on its rebuilt metric before the production CMC solve. Scientific admission is
never inferred from this diagnostic.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import cosmos_potential
from engine.v55_matter import (
    RHO_DM_PRESENT,
    RHO_B_PRESENT,
    RHO_R_PRESENT,
    initialize_from_archive,
)
from engine import v55_pirk_adapter as adapter
from tools.audit_continuous_origin_half_cell_oracle import continuous_reference
from tools.audit_origin_cell_joint_shadow_one_step import (
    VENDOR_COMMIT,
    SOURCE_COMMIT,
    apply_continuous_metric_shadow,
    build_continuous_metric_B,
    constraint_measure,
    metric_lambda_from_geometry,
    regular_origin_cmc_lapse,
    state_digest,
    verify_source_pins,
    _install_regular_row,
)

RESOLUTION = int(os.environ.get("ORIGIN_EVOLUTION_RESOLUTION", "120"))
R_MAX = 80.0
TARGET_TIME = 0.60
CFL = 0.03
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
INCLUDE_RADIATION = True
DEFAULT_ARMS = (
    "production_baseline",
    "half_cell_integral_only",
    "regular_cmc_row_only",
    "joint_integral_and_cmc_row",
)
SUPPORTED_ARMS = frozenset((*DEFAULT_ARMS, "full_continuous_metric_profile"))
ARMS = tuple(
    name.strip()
    for name in os.environ.get("ORIGIN_EVOLUTION_ARMS", ",".join(DEFAULT_ARMS)).split(",")
    if name.strip()
)
if not ARMS or len(set(ARMS)) != len(ARMS) or any(name not in SUPPORTED_ARMS for name in ARMS):
    raise ValueError(f"invalid ORIGIN_EVOLUTION_ARMS selection: {ARMS}")
OUT = Path(os.environ.get(
    "ORIGIN_EVOLUTION_OUT",
    f"results/origin-cell-joint-shadow-evolved-N{RESOLUTION}-t0.6-2026-10-24",
))
BASELINE_DIR = Path(os.environ.get(
    "BASELINE_ARTIFACT_DIR",
    f"baseline-reference/N{RESOLUTION}",
))
EXPECTED_CONTINUOUS_ORACLE_BLOB = "6e6f70cc4e57db67815f3cab41558d14a7c63797"

def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()

def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

def verify_candidate_inputs() -> dict:
    pin = verify_source_pins()
    source_blob = _git(
        "rev-parse", f"{SOURCE_COMMIT}:tools/audit_continuous_origin_half_cell_oracle.py"
    )
    runtime_blob = _git(
        "rev-parse", "HEAD:tools/audit_continuous_origin_half_cell_oracle.py"
    )
    if source_blob != EXPECTED_CONTINUOUS_ORACLE_BLOB or runtime_blob != source_blob:
        raise RuntimeError(
            f"continuous-origin oracle source drift: {source_blob} / {runtime_blob}"
        )
    manifest_path = BASELINE_DIR / "SHA256SUMS.json"
    report_path = BASELINE_DIR / "report.json"
    state_path = BASELINE_DIR / "final_state.npz"
    for path in (manifest_path, report_path, state_path):
        if not path.exists():
            raise FileNotFoundError(f"required baseline artifact file is absent: {path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for name, expected in manifest.items():
        artifact_path = BASELINE_DIR / name
        if not artifact_path.exists() or _sha256(artifact_path) != expected:
            raise RuntimeError(f"baseline artifact hash mismatch: {artifact_path}")
    source_report = json.loads(report_path.read_text(encoding="utf-8"))
    if source_report.get("source_pins", {}).get("baseline_commit") != SOURCE_COMMIT:
        raise RuntimeError("baseline artifact source commit mismatch")
    if source_report.get("source_pins", {}).get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError("baseline artifact vendor commit mismatch")
    cfg = source_report.get("configuration", {})
    if (
        int(cfg.get("resolution", -1)) != RESOLUTION
        or not math.isclose(float(cfg.get("r_max", float("nan"))), R_MAX, abs_tol=1e-12)
        or not math.isclose(float(cfg.get("target_time", float("nan"))), TARGET_TIME, abs_tol=1e-12)
        or not math.isclose(float(cfg.get("dt", float("nan"))), CFL * R_MAX / RESOLUTION, abs_tol=1e-12)
    ):
        raise RuntimeError(f"baseline artifact configuration does not match: {cfg}")
    return {
        "source_pin_check": pin,
        "continuous_oracle_blob": source_blob,
        "baseline_artifact": {
            "directory": str(BASELINE_DIR),
            "manifest_sha256": _sha256(manifest_path),
            "report_sha256": _sha256(report_path),
            "final_state_sha256": _sha256(state_path),
            "manifest": manifest,
            "configuration": cfg,
            "initial_state_sha256": source_report.get("initial_state_sha256"),
            "final_state_sha256_reported": source_report.get("final_state_sha256"),
        },
    }

def _candidate_half_cell_initializer(state) -> dict:
    """Replace only I(0->r0) in a copied initializer recurrence; reproject matter."""
    grid = state.grid
    r = np.asarray(grid.centers, dtype=np.float64)
    dr = float(grid.dr)
    r0 = 0.5 * dr
    n = int(grid.n)
    geom = state.geometry.copy()
    scalars = state.scalars
    S = np.asarray(scalars.S, dtype=np.float64)
    PS = np.asarray(scalars.PS, dtype=np.float64)
    D = np.asarray(scalars.D, dtype=np.float64)
    PD = np.asarray(scalars.PD, dtype=np.float64)
    phi = np.asarray(scalars.phi, dtype=np.float64)
    Pi = np.asarray(scalars.Pi, dtype=np.float64)

    rho_r = RHO_R_PRESENT if INCLUDE_RADIATION else 0.0
    phi0 = 0.179055
    pi0 = 0.00341913
    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += RHO_DM_PRESENT + RHO_B_PRESENT + rho_r
    H0 = math.sqrt(max(rho0 / 3.0, 0.0))
    Kt = np.full(n, -H0, dtype=np.float64)
    Sp = grid.cell_derivative_fourth(S, parity=1)
    Dp = grid.cell_derivative_fourth(D, parity=1)
    Ktp = grid.cell_derivative_fourth(Kt, parity=1)
    j = -(PD * Dp)
    rho_matter = (RHO_DM_PRESENT + RHO_B_PRESENT + rho_r) / (8.0 * math.pi)

    # The same analytic continuous-source half-cell integral used by the
    # archived origin oracle; it changes only the initial interval [0,r0].
    continuum = continuous_reference(r0)
    I0 = float(continuum["reference_moment"])
    continuous_B0 = float(continuum["reference_B"])
    baseline_B0 = float(np.asarray(state.geometry.X, dtype=np.float64)[0] ** 6)

    B = np.ones_like(r)
    source = np.empty_like(r)
    last_iteration_max_abs_change = None
    for _ in range(6):
        rhoS = 0.5 * (PS**2 + B * Sp**2) + 0.5 * S**2
        rhoD = 0.5 * (PD**2 + B * Dp**2) - 0.5 * D**2
        rho_phi = (0.5 * Pi**2 + cosmos_potential(phi)) / (8.0 * math.pi)
        rho = rhoS + rhoD + rho_phi + rho_matter
        Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
        source = (
            1.0 - 8.0 * math.pi * r * r * rho
            + 2.0 * r * r * Kr * Kt
            + r * r * Kt * Kt
        )
        integ = np.empty_like(r)
        integ[0] = I0
        if n > 1:
            integ[1:] = I0 + np.cumsum(
                0.5 * dr * (source[1:] + source[:-1])
            )
        B_new = integ / r
        last_iteration_max_abs_change = float(np.max(np.abs(B_new - B)))
        B = B_new

    if not np.all(np.isfinite(B)) or np.any(B <= 0.0):
        raise FloatingPointError("half-cell candidate generated invalid metric quotient")
    A = 1.0 / np.sqrt(B)
    Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
    K = Kr + 2.0 * Kt
    Aa = Kr - K / 3.0
    geom.a[:] = A ** (4.0 / 3.0)
    geom.b[:] = A ** (-2.0 / 3.0)
    geom.X[:] = A ** (-1.0 / 3.0)
    geom.K[:] = K
    geom.Aa[:] = Aa
    geom.Lambda[:] = metric_lambda_from_geometry(grid, geom)
    geom.alpha.fill(1.0)
    geom.beta.fill(0.0)
    geom.B[:] = 0.75 * geom.Lambda
    # Matter conserved variables are rebuilt from the candidate metric, as in
    # the actual initializer; unlike the prior metric-only shadow, they are not
    # frozen from the baseline geometry.
    candidate_matter = initialize_from_archive(
        grid, geom, include_radiation=INCLUDE_RADIATION
    )
    state.geometry = geom
    state.matter = candidate_matter
    return {
        "half_cell_integral": I0,
        "continuous_B0": continuous_B0,
        "baseline_discrete_B0": baseline_B0,
        "candidate_B0": float(B[0]),
        "baseline_minus_continuous_B0": float(baseline_B0 - continuous_B0),
        "candidate_minus_continuous_B0": float(B[0] - continuous_B0),
        "source_iterations": 6,
        "last_iteration_max_abs_change": last_iteration_max_abs_change,
        "all_B_finite_positive": True,
        "matter_reinitialized_on_candidate_metric": True,
    }

def _state_arrays(state) -> dict[str, np.ndarray]:
    arrays = {
        "t": np.asarray(state.t, dtype=np.float64),
        "tau": np.asarray(state.tau, dtype=np.float64),
        "e_folds": np.asarray(state.e_folds, dtype=np.float64),
        "r": np.asarray(state.grid.centers, dtype=np.float64),
    }
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        arrays[name] = np.asarray(getattr(state.geometry, name), dtype=np.float64)
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        arrays[name] = np.asarray(getattr(state.scalars, name), dtype=np.float64)
    for short, long in (
        ("dm", "dark_matter"),
        ("baryon", "baryons"),
        ("radiation", "radiation"),
    ):
        item = getattr(state.matter, long)
        for name in ("rest", "energy_t", "momentum_r"):
            arrays[f"{short}_{name}"] = np.asarray(getattr(item, name), dtype=np.float64)
    return arrays

def _compare_baseline_initial_geometry(state) -> dict:
    path = BASELINE_DIR / "stage_profiles.npz"
    actual = _state_arrays(state)
    keys = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda")
    result = {}
    with np.load(path, allow_pickle=False) as z:
        for name in keys:
            archive_key = f"0000_accepted_000__{name}"
            if archive_key not in z or name not in actual:
                result[name] = {"status": "missing_field"}
                continue
            expected = np.asarray(z[archive_key], dtype=np.float64)
            value = np.asarray(actual[name], dtype=np.float64)
            result[name] = {
                "status": "compared" if expected.shape == value.shape else "shape_mismatch",
                "max_abs_delta": float(np.max(np.abs(value - expected)))
                    if expected.shape == value.shape else None,
            }
    vals = [v["max_abs_delta"] for v in result.values()
            if v.get("status") == "compared" and v.get("max_abs_delta") is not None]
    return {
        "all_geometry_fields_compared": all(v.get("status") == "compared" for v in result.values()),
        "max_abs_delta_any_initial_geometry_field": max(vals) if vals else None,
        "field_max_abs_delta": result,
    }

def _compare_baseline_to_archive(state) -> dict:
    path = BASELINE_DIR / "final_state.npz"
    max_by_field = {}
    with np.load(path, allow_pickle=False) as z:
        actual = _state_arrays(state)
        for key in z.files:
            if key not in actual:
                max_by_field[key] = {"max_abs_delta": None, "status": "field_not_mapped"}
                continue
            expected = np.asarray(z[key], dtype=np.float64)
            value = np.asarray(actual[key], dtype=np.float64)
            if expected.shape != value.shape:
                max_by_field[key] = {
                    "max_abs_delta": None, "status": "shape_mismatch",
                    "expected_shape": list(expected.shape), "actual_shape": list(value.shape),
                }
            else:
                max_by_field[key] = {
                    "max_abs_delta": float(np.max(np.abs(value - expected))) if value.size else 0.0,
                    "status": "compared",
                }
    compared = [x["max_abs_delta"] for x in max_by_field.values()
                if x.get("status") == "compared" and x.get("max_abs_delta") is not None]
    return {
        "all_fields_compared": all(v.get("status") == "compared" for v in max_by_field.values()),
        "max_abs_delta_any_saved_field": max(compared) if compared else None,
        "field_max_abs_delta": max_by_field,
    }

def _metrics_at_state(state) -> tuple[dict, dict]:
    metrics, profiles = constraint_measure(state)
    r = np.asarray(state.grid.centers, dtype=np.float64)
    annulus = (r >= 2.0) & (r <= 8.0)
    H = profiles["H"]
    M = profiles["M"]
    C = profiles["C_lambda"]
    Aa = profiles["Aa"]
    r0 = 0.5 * float(state.grid.dr)
    metric = {
        **metrics,
        "Aa_over_r0_squared_cell0": float(Aa[0] / (r0 * r0)),
        "H_cell_0_5": [float(x) for x in H[:6]],
        "M_cell_0_5": [float(x) for x in M[:6]],
        "C_lambda_cell_0_5": [float(x) for x in C[:6]],
        "H_max_abs_2_le_r_le_8": float(np.max(np.abs(H[annulus]))),
        "M_max_abs_2_le_r_le_8": float(np.max(np.abs(M[annulus]))),
        "C_lambda_max_abs_2_le_r_le_8": float(np.max(np.abs(C[annulus]))),
        "Aa_at_shared_radii_1_to_5": {
            str(float(q)): float(np.interp(q, r, Aa)) for q in (1., 2., 3., 4., 5.)
        },
        "H_at_shared_radii_1_to_5": {
            str(float(q)): float(np.interp(q, r, H)) for q in (1., 2., 3., 4., 5.)
        },
        "M_at_shared_radii_1_to_5": {
            str(float(q)): float(np.interp(q, r, M)) for q in (1., 2., 3., 4., 5.)
        },
        "C_lambda_at_shared_radii_1_to_5": {
            str(float(q)): float(np.interp(q, r, C)) for q in (1., 2., 3., 4., 5.)
        },
    }
    return metric, profiles

def _run_arm(arm: str) -> tuple[dict, dict]:
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=RESOLUTION, r_max=R_MAX,
        amplitude=AMPLITUDE, width=WIDTH,
        D_amplitude=D_AMPLITUDE, include_radiation=INCLUDE_RADIATION,
    )
    candidate_info = None
    if arm == "full_continuous_metric_profile":
        B_continuous = build_continuous_metric_B(
            state.grid, subdivisions=2048
        )
        apply_continuous_metric_shadow(state, B_continuous)
        candidate_info = {
            "method": "rebuild the full continuous-source radial metric profile over all radii using 2048 substeps per coarse cell",
            "conformal_metric_profile": "a, b, X rebuilt from B=X^6; Lambda recomputed from that profile",
            "held_fixed": [
                "A_a and K from the source-pinned initializer",
                "all scalar fields and momenta",
                "all conserved matter arrays",
            ],
            "production_initializer_file_modified": False,
            "matter_reinitialized_on_candidate_metric": False,
        }
    elif arm in ("half_cell_integral_only", "joint_integral_and_cmc_row"):
        candidate_info = _candidate_half_cell_initializer(state)
    has_regular_cmc = arm in (
        "regular_cmc_row_only",
        "joint_integral_and_cmc_row",
        "full_continuous_metric_profile",
    )
    if has_regular_cmc:
        _install_regular_row(kernel)
        state.geometry.alpha = regular_origin_cmc_lapse(
            state.grid, state.geometry, state.scalars, state.matter
        )[0].copy()
    elif candidate_info is not None:
        state.geometry.alpha = kernel._solve_lapse_for_state(
            state.grid, state.geometry, state.scalars, state.matter,
            state.electromagnetic,
        )[0].copy()

    initial_digest = state_digest(state)
    initial_metric, initial_profiles = _metrics_at_state(state)
    time_history = [{
        "step": 0, "time": float(state.t), **{
            key: initial_metric[key] for key in (
                "Aa_cell0", "Aa_over_r0_squared_cell0",
                "H_max_abs_r_le_3", "M_max_abs_r_le_3", "C_lambda_max_abs_r_le_3",
                "H_max_abs_2_le_r_le_8", "M_max_abs_2_le_r_le_8",
                "C_lambda_max_abs_2_le_r_le_8",
            )
        }
    }]
    dt = CFL * float(state.grid.dr)
    steps = int(round(TARGET_TIME / dt))
    baseline_replay = None
    if arm == "production_baseline":
        # Compare the initial state digest and later the final fields to the
        # same pinned saved trajectory; no tolerance retargeting is applied.
        with open(BASELINE_DIR / "report.json", "r", encoding="utf-8") as f:
            baseline_report = json.load(f)
        baseline_replay = {
            "initial_state_sha256_recreated": initial_digest,
            "initial_geometry_comparison": _compare_baseline_initial_geometry(state),
            "note": "The local state-digest helper is not assumed identical to the archived report's digest function; initial geometry is compared field-by-field, and all saved final-state fields are compared below.",
        }
    for step_index in range(1, steps + 1):
        before = state_digest(state)
        next_state = kernel.step(state, dt)
        if state_digest(state) != before:
            raise RuntimeError(f"{arm}: kernel.step mutated the supplied input state at step {step_index}")
        state = next_state
        metric, _ = _metrics_at_state(state)
        time_history.append({
            "step": step_index, "time": float(state.t),
            "Aa_cell0": metric["Aa_cell0"],
            "Aa_over_r0_squared_cell0": metric["Aa_over_r0_squared_cell0"],
            "H_max_abs_r_le_3": metric["H_max_abs_r_le_3"],
            "M_max_abs_r_le_3": metric["M_max_abs_r_le_3"],
            "C_lambda_max_abs_r_le_3": metric["C_lambda_max_abs_r_le_3"],
            "H_max_abs_2_le_r_le_8": metric["H_max_abs_2_le_r_le_8"],
            "M_max_abs_2_le_r_le_8": metric["M_max_abs_2_le_r_le_8"],
            "C_lambda_max_abs_2_le_r_le_8": metric["C_lambda_max_abs_2_le_r_le_8"],
        })
    final_metric, final_profiles = _metrics_at_state(state)
    if not all(np.all(np.isfinite(v)) for v in _state_arrays(state).values()):
        raise FloatingPointError(f"{arm}: non-finite final state array")
    output_arrays = {
        f"{arm}_final_{key}": value.copy()
        for key, value in _state_arrays(state).items()
    }
    # Preserve cellwise final residuals for independent annulus/centre analysis.
    for key, value in final_profiles.items():
        output_arrays[f"{arm}_final_profile_{key}"] = np.asarray(value).copy()
    if arm == "production_baseline":
        baseline_replay["final_state_comparison"] = _compare_baseline_to_archive(state)
    if abs(float(state.t) - TARGET_TIME) > 1.0e-12:
        raise RuntimeError(f"{arm}: accepted time {state.t} != target {TARGET_TIME}")
    return ({
        "arm": arm,
        "resolution": RESOLUTION,
        "r_max": R_MAX,
        "dr": float(state.grid.dr),
        "dt": dt,
        "accepted_steps": steps,
        "final_time": float(state.t),
        "initial_state_sha256": initial_digest,
        "initializer_candidate": candidate_info,
        "regular_origin_cmc_row_used_at_initialization_and_each_stage": has_regular_cmc,
        "initial": initial_metric,
        "accepted": final_metric,
        "baseline_replay": baseline_replay,
        "time_history": time_history,
        "all_final_arrays_finite": True,
        "step_inputs_immutable": True,
    }, output_arrays)

def main() -> int:
    pin = verify_candidate_inputs()
    OUT.mkdir(parents=True, exist_ok=True)
    runs = []
    profiles = {}
    for arm in ARMS:
        row, arrays = _run_arm(arm)
        runs.append(row)
        profiles.update({f"N{RESOLUTION}_{key}": value for key, value in arrays.items()})
    profile_path = OUT / "profiles.npz"
    np.savez_compressed(profile_path, **profiles)
    report = {
        "schema": "origin-cell-joint-shadow-evolved/v1",
        "status": "completed_diagnostic_only",
        "source_pins": pin,
        "runtime": {
            "python": sys.version, "numpy": np.__version__,
            "scipy": scipy.__version__, "platform": platform.platform(),
        },
        "configuration": {
            "resolution": RESOLUTION, "r_max": R_MAX,
            "target_time": TARGET_TIME, "CFL": CFL,
            "dt_rule": "dt=0.03*dr", "amplitude": AMPLITUDE,
            "width": WIDTH, "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "arms": list(ARMS),
        },
        "initializer_candidate": {
            "scope": "candidate-specific details are recorded in each run's initializer_candidate field",
            "production_initializer_file_modified": False,
            "production_or_vendor_files_modified": False,
        },
        "runs": runs,
        "acceptance": {
            "all_final_arrays_finite": all(row["all_final_arrays_finite"] for row in runs),
            "all_step_inputs_immutable": all(row["step_inputs_immutable"] for row in runs),
            "production_source_or_vendor_changed": False,
            "baseline_match_is_a_diagnostic_not_a_science_gate": True,
            "scientific_admission": False,
        },
        "limitations": [
            "The half-cell initializer candidate is executed in a copied diagnostic recurrence; engine/v55_initial.py is not edited.",
            "The alternate CMC centre row is a diagnostic control, not an authorized production gauge change.",
            "Even the joint arm is a research candidate and is not production validation.",
            "A full t=0.6 trajectory does not by itself close the independent matter-coupled momentum-propagation identity or total stress-energy conservation gates.",
        ],
        "scientific_admission": False,
    }
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    checksum_path = OUT / "SHA256SUMS.txt"
    checksum_path.write_text(
        f"{_sha256(report_path)}  report.json\n{_sha256(profile_path)}  profiles.npz\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "resolution": RESOLUTION,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "baseline_final_state_max_abs_delta": next(
            (run.get("baseline_replay", {}).get("final_state_comparison", {}).get("max_abs_delta_any_saved_field")
             for run in runs if run["arm"] == "production_baseline"), None
        ),
        "output_dir": str(OUT),
        "report_sha256": _sha256(report_path),
        "profiles_sha256": _sha256(profile_path),
        "scientific_admission": False,
    }, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
