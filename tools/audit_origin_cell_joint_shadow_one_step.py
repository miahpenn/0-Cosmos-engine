#!/usr/bin/env python3
"""One-step production-kernel comparison from a joint regular-origin shadow.

Diagnostic-only: production and vendor source remain immutable. Four arms separate
the baseline, metric-only shadow, CMC-row-only shadow, and joint metric+CMC shadow.
All arms use the same V55ProductionKernel.step() PIRK ordering, RHS, and boundary
code; only the explicitly named regular-origin CMC-row arms substitute the central
row in the lapse solve. Scientific admission is never inferred from this run.
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

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

import numpy as np
import scipy
from scipy.linalg import solve_banded

from engine.cmc_gauge import _proper_volume_weights, _project_kdot_from_weights
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import total_matter_projection
from engine import v55_pirk_adapter as adapter
from tools.audit_continuous_origin_half_cell_oracle import continuous_coefficients

SOURCE_COMMIT = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (60, 120, 240)
R_MAX = 80.0
CFL = 0.03
METRIC_SUBDIVISIONS = 2048
OUTPUT_DIR = Path(os.environ.get(
    "ORIGIN_STEP_OUT",
    "results/origin-cell-joint-shadow-one-step-2026-10-24",
))
PROTECTED_PATHS = (
    "engine/production_kernel.py",
    "engine/v55_initial.py",
    "engine/cmc_gauge.py",
    "engine/v55_pirk_adapter.py",
    "engine/v55_matter.py",
    "engine/matter_system.py",
    "engine/matter_rhs.py",
    "engine/scalar_system.py",
    "engine/stress_energy.py",
    "engine/valencia.py",
    "engine/worldtube.py",
)

def git(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=REPOSITORY_ROOT, text=True
    ).strip()

def verify_source_pins() -> dict:
    expected_source = os.environ.get("SOURCE_COMMIT", SOURCE_COMMIT)
    expected_vendor = os.environ.get("VENDOR_COMMIT", VENDOR_COMMIT)
    if expected_source != SOURCE_COMMIT:
        raise RuntimeError(f"source pin env mismatch: {expected_source}")
    if expected_vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor pin env mismatch: {expected_vendor}")
    git("fetch", "origin", SOURCE_COMMIT)
    vendor = git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor gitlink mismatch: {vendor}")
    vendor_checkout = git(
        "-C", "vendor/bb-palatini-unified-r0", "rev-parse", "HEAD"
    )
    if vendor_checkout != VENDOR_COMMIT:
        raise RuntimeError(f"vendor checkout mismatch: {vendor_checkout}")
    protected_diff = git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--",
        "engine", "vendor", ".gitmodules", "requirements.txt",
    )
    if protected_diff:
        raise RuntimeError(
            "protected engine/vendor/dependency paths differ from source pin: "
            + protected_diff
        )
    blobs = {}
    for path in PROTECTED_PATHS:
        source_blob = git("rev-parse", f"{SOURCE_COMMIT}:{path}")
        head_blob = git("rev-parse", f"HEAD:{path}")
        if source_blob != head_blob:
            raise RuntimeError(
                f"protected blob mismatch for {path}: {source_blob} != {head_blob}"
            )
        blobs[path] = {"source_blob": source_blob, "runtime_blob": head_blob}
    return {
        "source_commit": SOURCE_COMMIT,
        "diagnostic_head": git("rev-parse", "HEAD"),
        "vendor_commit": vendor,
        "vendor_checkout": vendor_checkout,
        "protected_engine_vendor_dependency_diff_empty": True,
        "protected_source_blobs": blobs,
    }

def state_digest(state) -> str:
    digest = hashlib.sha256()
    def update(name: str, values) -> None:
        array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
        digest.update(name.encode("utf-8") + b"\0")
        digest.update(str(array.shape).encode("ascii") + b"\0")
        digest.update(array.tobytes())
    digest.update(f"t={state.t:.17e};tau={state.tau:.17e};efolds={state.e_folds:.17e}".encode())
    update("grid.centers", state.grid.centers)
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        update("geometry." + name, getattr(state.geometry, name))
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        update("scalars." + name, getattr(state.scalars, name))
    for species_name in ("dark_matter", "baryons", "radiation", "charged_dust"):
        item = getattr(state.matter, species_name, None)
        if item is None:
            digest.update(("matter." + species_name + "=None").encode())
            continue
        for name in ("rest", "energy_t", "momentum_r"):
            update(f"matter.{species_name}.{name}", getattr(item, name))
    return digest.hexdigest()

def build_continuous_metric_B(grid, subdivisions: int = METRIC_SUBDIVISIONS):
    """Sample the regular continuum initializer ODE at the production centres."""
    n = int(grid.n)
    dr = float(grid.dr)
    fine_n = n * int(subdivisions)
    h = dr / float(subdivisions)
    radii = np.arange(fine_n + 1, dtype=np.float64) * h
    coeff_A, coeff_q = continuous_coefficients(radii)
    integral = np.zeros(fine_n + 1, dtype=np.float64)
    B = np.ones(fine_n + 1, dtype=np.float64)
    source = np.empty(fine_n + 1, dtype=np.float64)
    source[0] = 1.0
    for i in range(1, fine_n + 1):
        denom = 1.0 + 0.5 * h * coeff_q[i] / radii[i]
        integral[i] = (
            integral[i - 1] + 0.5 * h * (source[i - 1] + coeff_A[i])
        ) / denom
        B[i] = integral[i] / radii[i]
        source[i] = coeff_A[i] - coeff_q[i] * B[i]
    sample_indices = np.arange(n, dtype=np.int64) * subdivisions + subdivisions // 2
    sampled = B[sample_indices]
    if sampled.shape != (n,) or not np.all(np.isfinite(sampled)) or np.any(sampled <= 0):
        raise FloatingPointError("continuous-origin metric B sample is invalid")
    return sampled

def metric_lambda_from_geometry(grid, geometry) -> np.ndarray:
    r = np.asarray(grid.centers, dtype=np.float64)
    a = np.asarray(geometry.a, dtype=np.float64)
    b = np.asarray(geometry.b, dtype=np.float64)
    ap = grid.cell_derivative_fourth(a, parity=1)
    bp = grid.cell_derivative_fourth(b, parity=1)
    return ap / (2.0 * a * a) - bp / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)

def apply_continuous_metric_shadow(state, B_continuous: np.ndarray) -> None:
    """Replace only conformal metric profiles; hold Aa/K/matter conserved data fixed."""
    geometry = state.geometry
    B = np.asarray(B_continuous, dtype=np.float64)
    # Same spherical conformal family as build_initial_data: det(gamma)=1 and b=X^2.
    X = np.power(B, 1.0 / 6.0)
    geometry.X = X.copy()
    geometry.b = X * X
    geometry.a = 1.0 / (geometry.b * geometry.b)
    geometry.Lambda = metric_lambda_from_geometry(state.grid, geometry)
    geometry.B = 0.75 * geometry.Lambda
    geometry.beta = np.zeros_like(geometry.beta)
    if not all(np.all(np.isfinite(np.asarray(getattr(geometry, name))))
               for name in ("a", "b", "X", "Lambda", "B")):
        raise FloatingPointError("non-finite metric shadow profile")
    if np.any(geometry.a <= 0.0) or np.any(geometry.b <= 0.0) or np.any(geometry.X <= 0.0):
        raise FloatingPointError("non-positive metric shadow profile")

def regular_origin_cmc_lapse(grid, geometry, scalars, matter):
    """Copied CMC solve; change only the regular-centre central row.

    Interior and outer rows, proper-volume Kdot projection, matter projection,
    coefficient definitions, and banded solve match engine/cmc_gauge.py.
    Central limit: alpha'' + c alpha' -> 6 alpha_2, with
    alpha_2 ~= (alpha[1]-alpha[0])/(2 h^2).
    """
    weights, total_weight = _proper_volume_weights(grid, geometry)
    total = total_matter_projection(grid, geometry, scalars, matter)
    a = np.asarray(geometry.a, dtype=np.float64)
    b = np.asarray(geometry.b, dtype=np.float64)
    X = np.asarray(geometry.X, dtype=np.float64)
    r = np.asarray(grid.centers, dtype=np.float64)
    inv = X * X / a
    ap = grid.cell_derivative_fourth(a, parity=1)
    bp = grid.cell_derivative_fourth(b, parity=1)
    Xp = grid.cell_derivative_fourth(X, parity=1)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r
    Q = (
        1.5 * np.asarray(geometry.Aa, dtype=np.float64) ** 2
        + np.asarray(geometry.K, dtype=np.float64) ** 2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total["rho"], dtype=np.float64)
            + np.asarray(total["pr"], dtype=np.float64)
            + 2.0 * np.asarray(total["pt"], dtype=np.float64)
        )
    )
    kdot = _project_kdot_from_weights(
        grid, geometry, scalars, matter, weights, total_weight
    )
    m = Q / inv
    rhs = -float(kdot) / inv
    n = int(grid.n)
    h = float(grid.dr)
    band = np.zeros((3, n), dtype=np.float64)
    bvec = np.zeros(n, dtype=np.float64)
    band[1, 0] = -3.0 / h**2 - m[0]
    band[0, 1] = 3.0 / h**2
    bvec[0] = rhs[0]
    for i in range(1, n - 1):
        ci = c[i] / (2.0 * h)
        band[2, i - 1] = 1.0 / h**2 - ci
        band[1, i] = -2.0 / h**2 - m[i]
        band[0, i + 1] = 1.0 / h**2 + ci
        bvec[i] = rhs[i]
    band[1, n - 1] = 1.0
    bvec[n - 1] = 1.0
    alpha = solve_banded((1, 1), band, bvec, check_finite=False)
    if not np.all(np.isfinite(alpha)) or float(np.min(alpha)) <= 0.0:
        raise FloatingPointError("regular-origin CMC shadow returned invalid lapse")
    return alpha, float(kdot)

def _install_regular_row(kernel) -> None:
    def solve_lapse_for_state(grid, geometry, scalars, matter, electromagnetic):
        if electromagnetic is not None:
            raise RuntimeError("this diagnostic intentionally excludes the optional Maxwell arm")
        return regular_origin_cmc_lapse(grid, geometry, scalars, matter)
    # Instance attribute is called with the same argument list as the kernel method.
    kernel._solve_lapse_for_state = solve_lapse_for_state

def constraint_measure(state) -> dict:
    grid, geom = state.grid, state.geometry
    _, vacuum, _ = adapter.vendor_modules()
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    raw = vacuum.constraints(grid, geom)
    H = np.asarray(raw["hamiltonian"], dtype=np.float64) - 16.0 * math.pi * np.asarray(total.rho, dtype=np.float64)
    M = np.asarray(raw["momentum"], dtype=np.float64) - 8.0 * math.pi * np.asarray(total.j, dtype=np.float64)
    lambda_metric = metric_lambda_from_geometry(grid, geom)
    C = np.asarray(geom.Lambda, dtype=np.float64) - lambda_metric
    r = np.asarray(grid.centers, dtype=np.float64)
    central = r <= 3.0
    arrays = {
        "Aa": np.asarray(geom.Aa, dtype=np.float64),
        "H": H,
        "M": M,
        "C_lambda": C,
        "Lambda": np.asarray(geom.Lambda, dtype=np.float64),
        "Lambda_metric": lambda_metric,
    }
    r0 = 0.5 * float(grid.dr)
    result = {
        "t": float(state.t),
        "r0": r0,
        "Aa_cell0": float(geom.Aa[0]),
        "Aa_over_r0_squared_cell0": float(geom.Aa[0] / (r0 * r0)),
        "H_cell_0_5": [float(x) for x in H[:6]],
        "M_cell_0_5": [float(x) for x in M[:6]],
        "C_lambda_cell_0_5": [float(x) for x in C[:6]],
        "Lambda_cell_0_5": [float(x) for x in np.asarray(geom.Lambda)[:6]],
        "Lambda_metric_cell_0_5": [float(x) for x in lambda_metric[:6]],
        "H_max_abs_r_le_3": float(np.max(np.abs(H[central]))),
        "M_max_abs_r_le_3": float(np.max(np.abs(M[central]))),
        "C_lambda_max_abs_r_le_3": float(np.max(np.abs(C[central]))),
        "H_max_abs_full_grid": float(np.max(np.abs(H))),
        "M_max_abs_full_grid": float(np.max(np.abs(M))),
        "C_lambda_max_abs_full_grid": float(np.max(np.abs(C))),
        "finite": all(np.all(np.isfinite(values)) for values in arrays.values()),
    }
    return result, arrays

def run_arm(resolution: int, arm: str, B_continuous: np.ndarray | None) -> tuple[dict, dict]:
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution, r_max=R_MAX, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True,
    )
    if B_continuous is not None:
        apply_continuous_metric_shadow(state, B_continuous)
    joint_cmc = arm in ("cmc_row_only", "joint_metric_and_cmc_row")
    if joint_cmc:
        _install_regular_row(kernel)
        state.geometry.alpha = regular_origin_cmc_lapse(
            state.grid, state.geometry, state.scalars, state.matter
        )[0].copy()
    elif arm == "metric_only":
        state.geometry.alpha = kernel._solve_lapse_for_state(
            state.grid, state.geometry, state.scalars, state.matter, state.electromagnetic
        )[0].copy()
    # Baseline state already carries the initialized production CMC lapse.

    initial_digest = state_digest(state)
    initial_metrics, initial_arrays = constraint_measure(state)
    before_step_digest = state_digest(state)
    dt = CFL * float(state.grid.dr)
    candidate = kernel.step(state, dt)
    after_step_input_digest = state_digest(state)
    if before_step_digest != after_step_input_digest:
        raise RuntimeError(f"{arm}/N={resolution}: production step mutated its input state")
    accepted_metrics, accepted_arrays = constraint_measure(candidate)
    if not initial_metrics["finite"] or not accepted_metrics["finite"]:
        raise FloatingPointError(f"{arm}/N={resolution}: non-finite constraint witness")
    if abs(candidate.t - dt) > 64.0 * np.finfo(float).eps * max(1.0, abs(dt)):
        raise RuntimeError(f"{arm}/N={resolution}: accepted time mismatch")
    profile = {}
    for stage, arrays in (("initial", initial_arrays), ("accepted", accepted_arrays)):
        for name, values in arrays.items():
            profile[f"N{resolution}_{arm}_{stage}_{name}"] = np.asarray(values).copy()
    row = {
        "arm": arm,
        "resolution": int(resolution),
        "r_max": R_MAX,
        "dr": float(state.grid.dr),
        "dt": dt,
        "target_t": dt,
        "initial_state_sha256": initial_digest,
        "input_state_unchanged_by_step": True,
        "metric_shadow": B_continuous is not None,
        "regular_origin_cmc_row_used_at_every_cmc_stage": bool(joint_cmc),
        "matter_conserved_variables_held_fixed_during_metric_rebuild": bool(B_continuous is not None),
        "initial": initial_metrics,
        "accepted": accepted_metrics,
    }
    return row, profile

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def main() -> int:
    pins = verify_source_pins()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    reports = []
    all_profiles = {}
    for resolution in RESOLUTIONS:
        # Generate the same regular continuum metric reference once, then reuse
        # it across the four copied one-step arms at this resolution.
        kernel = V55ProductionKernel()
        reference_state = kernel.initialize(
            resolution=resolution, r_max=R_MAX, amplitude=0.01, width=7.0,
            D_amplitude=1.0e-10, include_radiation=True,
        )
        B_continuous = build_continuous_metric_B(reference_state.grid)
        arms = []
        for arm in (
            "baseline",
            "metric_only",
            "cmc_row_only",
            "joint_metric_and_cmc_row",
        ):
            metric_shadow = B_continuous if arm in (
                "metric_only", "joint_metric_and_cmc_row"
            ) else None
            row, profiles = run_arm(resolution, arm, metric_shadow)
            arms.append(row)
            all_profiles.update(profiles)
        reports.append({
            "resolution": resolution,
            "dt": CFL * R_MAX / resolution,
            "subdivisions_per_coarse_dr": METRIC_SUBDIVISIONS,
            "arms": arms,
        })
    profile_path = OUTPUT_DIR / "profiles.npz"
    np.savez_compressed(profile_path, **all_profiles)
    report = {
        "schema": "origin-cell-joint-shadow-one-step/v1",
        "status": "completed_diagnostic_only",
        "source_pins": pins,
        "runtime": {
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
        },
        "configuration": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "amplitude": 0.01,
            "width": 7.0,
            "D_amplitude": 1.0e-10,
            "include_radiation": True,
            "CFL": CFL,
            "dt_rule": "dt = 0.03 * dr",
            "metric_reference": "continuous initializer radial-source ODE sampled at production centres; 2048 substeps per coarse dr",
            "shadow_limit": "a/b/X/Lambda are rebuilt from the continuous metric; Aa, K, scalar fields, and conserved matter variables are held fixed. The total stress projection is recomputed on each geometry.",
        },
        "arm_definitions": {
            "baseline": "Unmodified production initializer and pinned production CMC solve.",
            "metric_only": "Continuous-source metric shadow; CMC solve and every production stage retain the pinned central row.",
            "cmc_row_only": "Production metric; copied regular-origin CMC central row is used at initialization and each CMC stage.",
            "joint_metric_and_cmc_row": "Continuous-source metric plus copied regular-origin CMC central row at initialization and each CMC stage.",
            "shared": "All arms call the pinned V55ProductionKernel.step(); RHS, PIRK stage ordering, outer-boundary calls, matter/scalar updates, and accepted-state construction are otherwise unchanged.",
        },
        "runs": reports,
        "acceptance": {
            "all_initial_and_accepted_profiles_finite": True,
            "step_input_immutable": True,
            "production_source_or_vendor_changed": False,
            "scientific_admission": False,
        },
        "limitations": [
            "This is one accepted step only, not a converged continuum/evolution study.",
            "The metric shadow holds Aa/K and conserved matter variables fixed rather than performing a full metric-matter fixed-point initialization.",
            "The regular-origin CMC central row is a diagnostic sensitivity control, not an authorized production gauge change.",
            "A finite or improved one-step residual does not close the full matter-coupled momentum-propagation or total stress-energy conservation gates.",
        ],
        "scientific_admission": False,
    }
    report_path = OUTPUT_DIR / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    checksum_path = OUTPUT_DIR / "SHA256SUMS.txt"
    checksum_path.write_text(
        f"{sha256(report_path)}  report.json\n{sha256(profile_path)}  profiles.npz\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "resolutions": list(RESOLUTIONS),
        "arms_per_resolution": 4,
        "output_dir": str(OUTPUT_DIR),
        "report_sha256": sha256(report_path),
        "profiles_sha256": sha256(profile_path),
        "scientific_admission": False,
    }, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
