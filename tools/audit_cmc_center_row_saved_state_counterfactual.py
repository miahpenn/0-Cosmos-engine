#!/usr/bin/env python3
"""Offline CMC-centre-row counterfactual on saved origin-half-cell A/B states.

The observer verifies source-pinned artifact bundles, reconstructs the saved
accepted state at t=0.6, and compares the production alpha[0]=alpha[1] row with
a copied CMC solve using the exact pinned parity-aware D1/D2 operators on all
PDE rows. It advances no states and never edits production source.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import solve_cmc_lapse
from engine.matter_system import ConservedSpecies
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState, total_matter_projection
from tools import audit_evolved_origin_aa_regularity as source_audit

SOURCE_COMMIT = source_audit.SOURCE_COMMIT
VENDOR_COMMIT = source_audit.VENDOR_COMMIT
MEMBERS = ("report.json", "stage_profiles.npz", "final_state.npz", "SHA256SUMS.json")
EPS = float(np.finfo(np.float64).eps)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def _read_member(path: Path, name: str) -> bytes:
    if path.is_dir():
        target = path / name
        if not target.is_file():
            raise FileNotFoundError(f"{path} is missing {name}")
        return target.read_bytes()
    if not path.is_file() or not zipfile.is_zipfile(path):
        raise ValueError(f"expected ZIP or extracted artifact directory: {path}")
    with zipfile.ZipFile(path) as archive:
        try:
            return archive.read(name)
        except KeyError as exc:
            raise FileNotFoundError(f"{path} is missing {name}") from exc


def _load_bundle(path_text: str, resolution: int, arm: str) -> dict:
    path = Path(path_text).resolve()
    blobs = {name: _read_member(path, name) for name in MEMBERS}
    manifest = json.loads(blobs["SHA256SUMS.json"])
    member_hashes = {}
    for name in ("report.json", "stage_profiles.npz", "final_state.npz"):
        actual = _sha(blobs[name])
        if manifest.get(name) != actual:
            raise RuntimeError(
                f"{arm} N={resolution}: checksum mismatch for {name}; "
                f"expected={manifest.get(name)!r}, actual={actual}"
            )
        member_hashes[name] = actual

    report = json.loads(blobs["report.json"])
    if report.get("schema") != "source_pinned_evolved_origin_Aa_regularity_v1":
        raise RuntimeError(f"{arm}: unexpected report schema {report.get('schema')!r}")
    if report.get("scientific_admission") is not False:
        raise RuntimeError(f"{arm}: archive does not preserve scientific_admission=false")
    if report.get("production_or_vendor_files_modified") is not False:
        raise RuntimeError(f"{arm}: source archive indicates production/vendor edits")
    pins = report.get("source_pins", {})
    if (
        pins.get("baseline_commit") != SOURCE_COMMIT
        or pins.get("vendor_commit") != VENDOR_COMMIT
        or pins.get("production_vendor_dependency_diff_empty") is not True
    ):
        raise RuntimeError(f"{arm}: source/vendor/no-diff pins failed")
    for rel_path, row in pins.get("blobs", {}).items():
        expected = source_audit.SOURCE_BLOBS.get(rel_path)
        if (
            expected is None
            or row.get("expected_blob") != expected
            or row.get("base_blob") != expected
            or row.get("head_blob") != expected
            or row.get("matches") is not True
        ):
            raise RuntimeError(f"{arm}: source blob pin failed for {rel_path}")
    if set(pins.get("blobs", {})) != set(source_audit.SOURCE_BLOBS):
        raise RuntimeError(f"{arm}: source pin table is incomplete")

    config = report.get("configuration", {})
    for key, expected, atol in (
        ("resolution", resolution, 0.0),
        ("r_max", 80.0, 1e-12),
        ("target_time", 0.6, 1e-12),
        ("cfl", 0.03, 1e-12),
    ):
        value = float(config.get(key, math.nan))
        if not math.isfinite(value) or not math.isclose(value, expected, rel_tol=0.0, abs_tol=atol):
            raise RuntimeError(f"{arm}: configuration mismatch for {key}: {value}")
    if arm == "candidate" and report.get("initialization_arm") != "origin_half_cell_candidate":
        raise RuntimeError("candidate bundle does not identify the origin-half-cell initializer")

    checkpoint = np.load(io.BytesIO(blobs["final_state.npz"]), allow_pickle=False)
    arrays = {name: np.asarray(checkpoint[name]).copy() for name in checkpoint.files}
    expected_r = (np.arange(resolution, dtype=float) + 0.5) * (80.0 / resolution)
    if (
        arrays.get("r", np.empty(0)).shape != (resolution,)
        or not np.all(np.isfinite(arrays["r"]))
        or not np.allclose(arrays["r"], expected_r, rtol=0.0, atol=1e-12)
    ):
        raise RuntimeError(f"{arm}: checkpoint grid does not match N={resolution}, r_max=80")
    if not math.isclose(float(arrays["t"]), 0.6, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError(f"{arm}: final checkpoint is not at t=0.6")
    checkpoint_digest = _sha(blobs["final_state.npz"])
    manifest_digest = _sha(blobs["SHA256SUMS.json"])
    return {
        "path": str(path),
        "report": report,
        "arrays": arrays,
        "member_hashes": member_hashes,
        "checkpoint_sha256": checkpoint_digest,
        "manifest_sha256": manifest_digest,
    }


def _restore_state(kernel: V55ProductionKernel, bundle: dict, resolution: int) -> ProductionState:
    arrays = bundle["arrays"]
    # Initialize only to obtain the pinned grid and geometry container types;
    # every saved evolution field below is then restored from the verified NPZ.
    template = kernel.initialize(resolution=resolution, r_max=80.0, include_radiation=True)
    geom = template.geometry.copy()
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        values = np.asarray(arrays[name], dtype=float)
        if values.shape != template.grid.centers.shape or not np.all(np.isfinite(values)):
            raise ValueError(f"checkpoint geometry field {name} is malformed")
        setattr(geom, name, values.copy())

    scalars = ScalarFields(*(np.asarray(arrays[name], dtype=float).copy()
                             for name in ("S", "PS", "D", "PD", "phi", "Pi")))
    matter = V55MatterState(
        ConservedSpecies(*(np.asarray(arrays[prefix + suffix], dtype=float).copy()
                           for suffix in ("_rest", "_energy_t", "_momentum_r")))
        for prefix in ()
    )
    # Build the three conserved sectors explicitly to keep their order clear.
    matter = V55MatterState(
        dark_matter=ConservedSpecies(
            arrays["dm_rest"].copy(), arrays["dm_energy_t"].copy(), arrays["dm_momentum_r"].copy()
        ),
        baryons=ConservedSpecies(
            arrays["baryon_rest"].copy(), arrays["baryon_energy_t"].copy(), arrays["baryon_momentum_r"].copy()
        ),
        radiation=ConservedSpecies(
            arrays["radiation_rest"].copy(), arrays["radiation_energy_t"].copy(), arrays["radiation_momentum_r"].copy()
        ),
    )
    state = ProductionState(
        grid=template.grid,
        geometry=geom,
        scalars=scalars,
        matter=matter,
        t=float(arrays["t"]),
        tau=float(arrays["tau"]),
        e_folds=float(arrays["e_folds"]),
        electromagnetic=None,
    )
    saved_digest = source_audit._state_digest(state)
    expected_digest = bundle["report"].get("final_state_sha256")
    if saved_digest != expected_digest:
        raise RuntimeError(
            f"restored {bundle['path']} checkpoint digest mismatch: "
            f"expected={expected_digest}, actual={saved_digest}"
        )
    return state


def _d_matrices(grid) -> tuple[np.ndarray, np.ndarray]:
    n = grid.n
    d1 = np.zeros((n, n), dtype=float)
    d2 = np.zeros((n, n), dtype=float)
    for j in range(n):
        basis = np.zeros(n, dtype=float)
        basis[j] = 1.0
        d1[:, j] = grid.cell_derivative_fourth(basis, parity=1)
        d2[:, j] = grid.cell_second_derivative_fourth(basis, parity=1)
    return d1, d2


def _operator_compatible_alpha(state: ProductionState) -> tuple[np.ndarray, dict]:
    """Copy the CMC solve, changing only its centre row to the production D1/D2 PDE."""
    grid, geom = state.grid, state.geometry
    n = grid.n
    if n < 5:
        raise ValueError("operator-compatible CMC row requires at least five cells")
    if float(np.max(np.abs(np.asarray(geom.beta, dtype=float)))) != 0.0:
        raise ValueError("the copied CMC row audit expects the zero-shift branch")

    _, target_kdot = solve_cmc_lapse(grid, geom.copy(), state.scalars, state.matter)
    a = np.asarray(geom.a, dtype=float)
    b = np.asarray(geom.b, dtype=float)
    X = np.asarray(geom.X, dtype=float)
    Aa = np.asarray(geom.Aa, dtype=float)
    K = np.asarray(geom.K, dtype=float)
    inv = X * X / a
    if not all(np.all(np.isfinite(v)) for v in (a, b, X, inv)) or np.any(inv <= 0.0):
        raise FloatingPointError("invalid inverse radial metric in CMC counterfactual")

    d1_mat, d2_mat = _d_matrices(grid)
    d1a = d1_mat @ a
    d1b = d1_mat @ b
    d1x = d1_mat @ X
    c = -0.5 * d1a / a + d1b / b - d1x / X + 2.0 / grid.centers
    total = total_matter_projection(grid, geom, state.scalars, state.matter)
    q = (
        1.5 * Aa * Aa + K * K / 3.0
        + 4.0 * math.pi * (
            np.asarray(total["rho"], dtype=float)
            + np.asarray(total["pr"], dtype=float)
            + 2.0 * np.asarray(total["pt"], dtype=float)
        )
    )
    m = q / inv

    _, vacuum, _ = adapter.vendor_modules()
    l3_vacuum = np.asarray(vacuum.primary_l3_rhs(grid, geom)["K"], dtype=float)
    quadratic = Aa * Aa + 2.0 * (-0.5 * Aa) ** 2 + K * K / 3.0
    ko_k = l3_vacuum - np.asarray(geom.alpha, dtype=float) * quadratic
    rhs = (ko_k - float(target_kdot)) / inv
    matrix = d2_mat + c[:, None] * d1_mat - np.diag(m)
    matrix[-1, :] = 0.0
    matrix[-1, -1] = 1.0
    rhs[-1] = 1.0

    alpha = np.linalg.solve(matrix, rhs)
    if not np.all(np.isfinite(alpha)) or float(np.min(alpha)) <= 0.0:
        raise FloatingPointError("operator-compatible CMC lapse is not finite and positive")

    residual = matrix @ alpha - rhs
    row_scale = np.abs(d2_mat) @ np.abs(alpha) + np.abs(c) * (np.abs(d1_mat) @ np.abs(alpha))
    row_scale += np.abs(m * alpha) + np.abs(rhs)
    row_scale[-1] = abs(alpha[-1]) + 1.0
    row_bound = 512.0 * EPS * row_scale
    operator_residual = d2_mat @ alpha + c * (d1_mat @ alpha) - m * alpha - (ko_k - float(target_kdot)) / inv
    operator_residual[-1] = alpha[-1] - 1.0
    if (
        not np.all(np.abs(residual) <= row_bound)
        or not np.all(np.abs(operator_residual) <= row_bound)
    ):
        raise FloatingPointError("operator-compatible CMC row residual exceeds roundoff envelope")

    return alpha, {
        "target_kdot_from_existing_CMC_projection": float(target_kdot),
        "matrix_max_abs_residual": float(np.max(np.abs(residual))),
        "operator_max_abs_residual": float(np.max(np.abs(operator_residual))),
        "raw_operand_scaled_rows_pass": True,
        "pde_rows": int(n - 1),
        "outer_normalization": float(alpha[-1]),
    }


def _fit_even(r: np.ndarray, values: np.ndarray, *, start: int, stop: int, degree: int) -> np.ndarray:
    x = r[start:stop] ** 2
    y = np.asarray(values[start:stop], dtype=float)
    actual_degree = min(degree, x.size - 1)
    return np.polynomial.polynomial.polyfit(x, y, deg=actual_degree)


def _channel_snapshot(state: ProductionState) -> dict:
    grid, geom = state.grid, state.geometry
    _, vacuum, _ = adapter.vendor_modules()
    terms = vacuum.geometry_terms(grid, geom)
    l2 = vacuum.primary_l2_rhs(grid, geom)
    l3 = adapter.primary_l3_with_matter(grid, geom, state.scalars, state.matter)
    r = np.asarray(grid.centers, dtype=float)
    inv = np.asarray(geom.X, dtype=float) ** 2 / np.asarray(geom.a, dtype=float)
    c2 = float(_fit_even(r, np.asarray(geom.alpha), start=1, stop=4, degree=2)[1])
    inv0 = float(_fit_even(r, inv, start=1, stop=4, degree=2)[0])
    k_l2_regular = -6.0 * inv0 * c2

    alpha_fit = _fit_even(r, np.asarray(geom.alpha), start=1, stop=min(8, r.size), degree=3)
    k_fit = _fit_even(r, np.asarray(l2["K"]), start=1, stop=min(8, r.size), degree=3)
    aa_lapse = -(np.asarray(terms["DrDr_alpha"]) - np.asarray(terms["lap_alpha"]) / 3.0)
    aa_ricci = np.asarray(geom.alpha) * (
        np.asarray(terms["Rr_r"]) - np.asarray(terms["R"]) / 3.0
    )
    constraints = vacuum.constraints(grid, geom)
    return {
        "alpha0": float(geom.alpha[0]),
        "alpha1": float(geom.alpha[1]),
        "alpha0_minus_even_fit": float(geom.alpha[0] - alpha_fit[0]),
        "alpha_c2_from_cells_1_to_3": c2,
        "inverse_metric_center_estimate": inv0,
        "K_L2_cell0": float(l2["K"][0]),
        "K_L2_even_fit_intercept_cells_1_to_7": float(k_fit[0]),
        "K_L2_cell0_minus_fit": float(l2["K"][0] - k_fit[0]),
        "K_L2_regular_even_center_estimate": float(k_l2_regular),
        "K_L2_cell0_minus_regular_even_estimate": float(l2["K"][0] - k_l2_regular),
        "K_L2_cells_0_to_3": [float(x) for x in np.asarray(l2["K"][:4])],
        "K_L3_cells_0_to_3": [float(x) for x in np.asarray(l3["K"][:4])],
        "K_total_cells_0_to_3": [float(x) for x in (np.asarray(l2["K"]) + np.asarray(l3["K"]))[:4]],
        "Aa_L2_cells_0_to_3": [float(x) for x in np.asarray(l2["Aa"][:4])],
        "Aa_L2_lapse_hessian_cells_0_to_3": [float(x) for x in aa_lapse[:4]],
        "Aa_L2_tracefree_ricci_cells_0_to_3": [float(x) for x in aa_ricci[:4]],
        "Aa_L2_split_max_abs_closure": float(np.max(np.abs(
            np.asarray(l2["Aa"]) - aa_lapse - aa_ricci
        ))),
        "constraints": {name: [float(v) for v in np.asarray(values[:4])]
                        for name, values in constraints.items()},
        "constraints_max_abs": {name: float(np.max(np.abs(values)))
                                for name, values in constraints.items()},
    }


def _replace_only_lapse(state: ProductionState, alpha: np.ndarray) -> ProductionState:
    geom = state.geometry.copy()
    geom.alpha = np.asarray(alpha, dtype=float).copy()
    return ProductionState(
        grid=state.grid,
        geometry=geom,
        scalars=state.scalars,
        matter=state.matter,
        t=state.t,
        tau=state.tau,
        e_folds=state.e_folds,
        cycle=state.cycle,
        history=list(state.history),
        handoffs=list(state.handoffs),
        electromagnetic=None,
    )


def _audit_pair(kernel: V55ProductionKernel, n: int, baseline: dict, candidate: dict) -> dict:
    outputs = {}
    for arm, bundle in (("baseline", baseline), ("candidate", candidate)):
        state = _restore_state(kernel, bundle, n)
        before_digest = source_audit._state_digest(state)
        production = _channel_snapshot(state)
        alpha_cf, solve_report = _operator_compatible_alpha(state)
        cf_state = _replace_only_lapse(state, alpha_cf)
        counterfactual = _channel_snapshot(cf_state)
        _, vacuum, _ = adapter.vendor_modules()
        l2_cf = vacuum.primary_l2_rhs(none_grid := state.grid, cf_state.geometry)
        l3_cf = adapter.primary_l3_with_matter(none_grid, cf_state.geometry, cf_state.scalars, cf_state.matter)
        target = solve_report["target_kdot_from_existing_CMC_projection"]
        kdot_cf = np.asarray(l2_cf["K"], dtype=float) + np.asarray(l3_cf["K"], dtype=float)
        error = kdot_cf - target
        scale = np.abs(np.asarray(cf_state.geometry.X)**2 / np.asarray(cf_state.geometry.a)) * (
            np.abs(np.asarray(l2_cf["K"])) + np.abs(np.asarray(l3_cf["K"]))
        ) + abs(target)
        bound = 1024.0 * EPS * np.maximum(scale, np.finfo(float).tiny)
        if not np.all(np.abs(error[:-1]) <= bound[:-1]):
            raise FloatingPointError(f"{arm} N={n}: production K RHS does not close to fixed CMC target")
        constraint_delta = {
            key: float(np.max(np.abs(
                np.asarray(counterfactual["constraints"][key], dtype=float)
                - np.asarray(production["constraints"][key], dtype=float)
            )))
            for key in production["constraints"]
        }
        # The state itself remains frozen; only the diagnostic copy's alpha changes.
        after_digest = source_audit._state_digest(state)
        if before_digest != after_digest:
            raise RuntimeError(f"{arm} N={n}: counterfactual mutated the saved state")

        outputs[arm] = {
            "saved_time": float(state.t),
            "checkpoint_sha256": bundle["checkpoint_sha256"],
            "checkpoint_member_sha256": bundle["member_hashes"]["final_state.npz"],
            "saved_state_sha256": before_digest,
            "production_lapse": production,
            "operator_compatible_lapse": counterfactual,
            "operator_compatible_solve": solve_report,
            "lapse_change_max_abs": float(np.max(np.abs(alpha_cf - np.asarray(state.geometry.alpha)))),
            "lapse_change_rms": float(np.sqrt(np.mean((alpha_cf - np.asarray(state.geometry.alpha)) ** 2))),
            "K_rhs_minus_fixed_target_max_abs_in_pde_rows": float(np.max(np.abs(error[:-1]))),
            "K_rhs_raw_operand_scaled_closure_pass": True,
            "static_constraint_max_abs_delta": constraint_delta,
            "saved_input_unchanged": before_digest == after_digest,
            "interpretation_limit": "This is a frozen saved-state lapse/operator counterfactual; it does not evolve the candidate or establish continuum convergence.",
        }

    return {
        "resolution": n,
        "time": 0.6,
        "baseline_run_id": 40127716287,
        "candidate_run_id": 40133940265 if n == 60 else 40128647367,
        "arms": outputs,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=int, choices=(60, 120, 240), required=True)
    parser.add_argument("--baseline", required=True, help="saved baseline artifact ZIP or directory")
    parser.add_argument("--candidate", required=True, help="saved origin-half-cell candidate artifact ZIP or directory")
    parser.add_argument("--out", required=True, help="output JSON path")
    args = parser.parse_args()

    provenance = source_audit.verify_source_pins()
    if provenance["vendor_commit"] != VENDOR_COMMIT:
        raise RuntimeError("runtime vendor pin mismatch")
    n = args.resolution
    baseline = _load_bundle(args.baseline, n, "baseline")
    candidate = _load_bundle(args.candidate, n, "candidate")
    kernel = V55ProductionKernel()
    resolution_result = _audit_pair(kernel, n, baseline, candidate)
    report = {
        "schema": "current_source_cmc_center_row_saved_state_counterfactual_v1",
        "status": "completed_offline_saved_state_counterfactual",
        "source_pins": provenance,
        "configuration": {
            "resolution": n,
            "r_max": 80.0,
            "time": 0.6,
            "baseline_source_commit": SOURCE_COMMIT,
            "vendor_commit": VENDOR_COMMIT,
            "baseline_reused_without_reevolution": True,
            "candidate_reused_without_reevolution": True,
            "production_evolution_steps_advanced": 0,
        },
        "resolution_result": resolution_result,
        "production_or_vendor_files_modified": False,
        "scientific_admission": False,
        "interpretation_limits": [
            "The copied PDE-row solve shares the pinned derivative operators with the production K RHS; it tests operator compatibility, not an independent discretization.",
            "The CMC target is held at the value projected from the saved production lapse so the row change is the isolated counterfactual; this is not a joint nonlinear CMC solve.",
            "Static H/M/connection/determinant witnesses are expected to remain unchanged because only alpha is replaced; this test reports RHS and lapse response, not evolved constraint propagation.",
            "No production equation, initializer, coefficient, gauge, boundary, dependency, vendor file, threshold, or admission flag is changed.",
        ],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "resolution": n,
        "out": str(out),
        "baseline_lapse_change_max_abs": resolution_result["arms"]["baseline"]["lapse_change_max_abs"],
        "candidate_lapse_change_max_abs": resolution_result["arms"]["candidate"]["lapse_change_max_abs"],
        "baseline_K_rhs_closure": resolution_result["arms"]["baseline"]["K_rhs_minus_fixed_target_max_abs_in_pde_rows"],
        "candidate_K_rhs_closure": resolution_result["arms"]["candidate"]["K_rhs_minus_fixed_target_max_abs_in_pde_rows"],
        "scientific_admission": False,
    }, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
