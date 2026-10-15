"""One-step coupled evolution/constraint-propagation control for off-constraint Lambda modes.

Research-only: compare pinned baseline and copied initial states differing only by
a manufactured smooth odd Lambda perturbation. Two dt values test whether the actual
Aa response approaches the frozen L2 prediction. Production sources are untouched.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_pirk_adapter import vendor_modules
from tools.audit_current_main_centre_proper_radial_ricci import (
    _finite, _max_abs, _state_digest, _verify_source_pins, proper_radial_ricci,
)

PI = float(np.pi)
SOURCE_COMMIT = "285a2ee0c1a3b97cddeb4519d8e749109b115aab"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (80, 160)
WAVELENGTHS = (16.0, 8.0, 4.0)
DT_OVER_DR = (0.01, 0.005)
AMPLITUDE = 1.0e-4
R_MAX = 80.0
OUT = Path(os.environ.get(
    "COUPLED_LAMBDA_OUT",
    "results/coupled-lambda-constraint-propagation-2026-10-15",
))


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _profile_summary(values, r, mask, centre_cells=6) -> dict:
    values = np.asarray(values, dtype=float)
    r = np.asarray(r, dtype=float)
    if values.shape != r.shape or not _finite(values):
        raise ValueError("profile must be finite and match the radial grid")
    imax = int(np.argmax(np.abs(values)))
    n = min(centre_cells, values.size)
    ic = int(np.argmax(np.abs(values[:n])))
    return {
        "cell0_signed": float(values[0]),
        "cell1_signed": float(values[1]),
        "first_six_signed": [float(x) for x in values[:n]],
        "first_six_absmax": float(abs(values[ic])),
        "first_six_signed_at_absmax": float(values[ic]),
        "first_six_r_at_absmax": float(r[ic]),
        "full_absmax": float(abs(values[imax])),
        "full_signed_at_absmax": float(values[imax]),
        "full_r_at_absmax": float(r[imax]),
        "full_rms": float(np.sqrt(np.mean(values * values))),
        "fixed_interior_absmax": _max_abs(values[mask]),
        "fixed_interior_rms": float(np.sqrt(np.mean(values[mask] ** 2))),
    }


def _constraint_profiles(state, vacuum) -> dict[str, np.ndarray]:
    """Return BSSN and independently reconstructed metric-only constraint witnesses."""
    grid, geom = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=float)
    raw = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    # The BSSN Hamiltonian intentionally uses the independently evolved Lambda
    # inside its Ricci scalar. Keep that Lambda-dependent witness explicit.
    H_bssn = (
        np.asarray(raw["hamiltonian"], dtype=float)
        - 16.0 * PI * np.asarray(total.rho, dtype=float)
    )
    M = np.asarray(raw["momentum"], dtype=float) - 8.0 * PI * np.asarray(total.j, dtype=float)
    d1 = lambda v, parity=1: grid.cell_derivative_fourth(v, parity=parity)
    a, b = np.asarray(geom.a, dtype=float), np.asarray(geom.b, dtype=float)
    X = np.asarray(geom.X, dtype=float)
    lambda_metric = d1(a, 1) / (2.0 * a * a) - d1(b, 1) / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)

    # Independent metric-only Hamiltonian: reconstruct the spatial Ricci scalar
    # from areal radius and proper radial distance, then add the same algebraic
    # extrinsic-curvature and matter terms. This channel does not use geom.Lambda.
    metric_ricci = proper_radial_ricci(
        grid, r * np.sqrt(b) / X, np.sqrt(a) / X
    )
    R_metric = (
        np.asarray(metric_ricci["Ric_l_l"], dtype=float)
        + 2.0 * np.asarray(metric_ricci["Ric_theta_theta"], dtype=float)
    )
    Aa = np.asarray(geom.Aa, dtype=float)
    Ab = -0.5 * Aa
    K = np.asarray(geom.K, dtype=float)
    H_metric_only = (
        R_metric - (Aa * Aa + 2.0 * Ab * Ab)
        + (2.0 / 3.0) * K * K
        - 16.0 * PI * np.asarray(total.rho, dtype=float)
    )

    result = {
        "hamiltonian": H_bssn,
        "hamiltonian_metric_only": H_metric_only,
        "momentum": M,
        "connection": np.asarray(geom.Lambda, dtype=float) - lambda_metric,
        "determinant": a * b * b - 1.0,
    }
    if not all(_finite(v) for v in result.values()):
        raise FloatingPointError("non-finite constraint witness after coupled step")
    return result


def _relative_rms(numerator, denominator, mask) -> float:
    n = np.asarray(numerator, dtype=float)[mask]
    d = np.asarray(denominator, dtype=float)[mask]
    scale = max(float(np.sqrt(np.mean(d*d))), np.finfo(float).tiny)
    return float(np.sqrt(np.mean(n*n)) / scale)


def _run_pair(kernel, base_state, vacuum, wavelength: float, dt: float) -> tuple[dict, dict]:
    grid = base_state.grid
    r = np.asarray(grid.centers, dtype=float)
    dr = float(grid.dr)
    geom = base_state.geometry
    k = 2.0 * PI / float(wavelength)
    delta_lambda = AMPLITUDE * np.sin(k * r)
    exact_delta_lambda_prime = AMPLITUDE * k * np.cos(k * r)
    d1 = lambda v, parity=1: grid.cell_derivative_fourth(v, parity=parity)
    d4 = lambda v: grid.cell_derivative_fourth(v, parity=-1)
    a = np.asarray(geom.a, dtype=float)
    X = np.asarray(geom.X, dtype=float)
    predicted_tf_response = (
        (2.0 / 3.0) * X * X * d4(delta_lambda)
        + 0.5 * X * X * (d1(a, 1) / a) * delta_lambda
    )
    predicted_aa_rate = np.asarray(geom.alpha, dtype=float) * predicted_tf_response
    interior = (r >= 5.0) & (r <= R_MAX - 5.0)
    if np.count_nonzero(interior) < 10:
        raise ValueError("fixed interior region has too few cells")

    base_before = _state_digest(base_state)
    base_input = copy.deepcopy(base_state)
    off_input = copy.deepcopy(base_state)
    off_input.geometry.Lambda = np.asarray(off_input.geometry.Lambda, dtype=float) + delta_lambda
    off_before = _state_digest(off_input)
    base_new = kernel.step(base_input, dt)
    off_new = kernel.step(off_input, dt)
    if _state_digest(base_state) != base_before:
        raise RuntimeError("step observer mutated original initialized state")
    if _state_digest(base_input) != base_before:
        raise RuntimeError("baseline input state changed during evolution")
    if _state_digest(off_input) != off_before:
        raise RuntimeError("step observer mutated perturbed input state")

    def garr(state, name):
        return np.asarray(getattr(state.geometry, name), dtype=float)

    metric_delta = {
        name: garr(off_new, name) - garr(base_new, name)
        for name in ("a", "b", "X", "Aa", "K", "Lambda", "alpha")
    }
    actual_aa_rate = metric_delta["Aa"] / dt
    aa_rate_error = actual_aa_rate - predicted_aa_rate
    constraints_base = _constraint_profiles(base_new, vacuum)
    constraints_off = _constraint_profiles(off_new, vacuum)
    constraint_delta = {
        key: constraints_off[key] - constraints_base[key]
        for key in constraints_base
    }
    initial_constraints_base = _constraint_profiles(base_state, vacuum)
    initial_constraints_off = _constraint_profiles(off_input, vacuum)

    # At fixed metric, the pinned BSSN scalar-Ricci bracket gives the signed
    # identity delta(H_BSSN) = X^2 D4(delta Lambda). The independent metric-only
    # Hamiltonian must remain unchanged; do not incorrectly require BSSN H to.
    initial_bssn_h_delta = (
        initial_constraints_off["hamiltonian"] - initial_constraints_base["hamiltonian"]
    )
    predicted_initial_bssn_h_delta = X * X * d4(delta_lambda)
    initial_bssn_h_closure = (
        initial_bssn_h_delta - predicted_initial_bssn_h_delta
    )
    initial_bssn_h_scale = max(
        _max_abs(initial_bssn_h_delta),
        _max_abs(predicted_initial_bssn_h_delta),
        np.finfo(float).tiny,
    )
    initial_bssn_h_tolerance = max(
        5.0e-12, 2048.0 * np.finfo(float).eps * initial_bssn_h_scale
    )
    if _max_abs(initial_bssn_h_closure) > initial_bssn_h_tolerance:
        raise RuntimeError(
            "fixed-metric BSSN Hamiltonian/Lambda response identity failed: "
            f"{_max_abs(initial_bssn_h_closure):.6e}"
        )

    error_rel = _relative_rms(aa_rate_error, predicted_aa_rate, interior)
    actual_rms = float(np.sqrt(np.mean(actual_aa_rate[interior] ** 2)))
    predicted_rms = float(np.sqrt(np.mean(predicted_aa_rate[interior] ** 2)))
    all_arrays = [
        delta_lambda, exact_delta_lambda_prime, predicted_tf_response,
        predicted_aa_rate, actual_aa_rate, aa_rate_error,
        *metric_delta.values(), *constraint_delta.values(),
        initial_bssn_h_delta, predicted_initial_bssn_h_delta, initial_bssn_h_closure,
    ]
    if not all(_finite(v) for v in all_arrays):
        raise FloatingPointError("non-finite coupled-mode response")

    report = {
        "resolution": int(grid.n), "dr": dr,
        "wavelength": float(wavelength), "points_per_wavelength": float(wavelength / dr),
        "wavenumber": float(k), "dt": float(dt), "dt_over_dr": float(dt / dr),
        "perturbation_amplitude": AMPLITUDE,
        "original_baseline_state_sha256": base_before,
        "off_mode_input_state_sha256": off_before,
        "original_baseline_state_unchanged": _state_digest(base_state) == base_before,
        "baseline_input_state_unchanged": _state_digest(base_input) == base_before,
        "off_mode_input_state_unchanged": _state_digest(off_input) == off_before,
        "one_accepted_step_per_arm": True,
        "actual_Aa_rate_delta_vs_frozen_L2_prediction_relative_rms_error": error_rel,
        "actual_Aa_rate_delta_rms_fixed_interior": actual_rms,
        "predicted_Aa_rate_delta_rms_fixed_interior": predicted_rms,
        "actual_Aa_rate_delta_max_abs_fixed_interior": _max_abs(actual_aa_rate[interior]),
        "predicted_Aa_rate_delta_max_abs_fixed_interior": _max_abs(predicted_aa_rate[interior]),
        "delta_Aa_final_max_abs": _max_abs(metric_delta["Aa"]),
        "delta_K_final_max_abs": _max_abs(metric_delta["K"]),
        "delta_a_final_max_abs": _max_abs(metric_delta["a"]),
        "delta_b_final_max_abs": _max_abs(metric_delta["b"]),
        "delta_X_final_max_abs": _max_abs(metric_delta["X"]),
        "delta_Lambda_final_max_abs": _max_abs(metric_delta["Lambda"]),
        "lambda_mode_input_max_abs": _max_abs(delta_lambda),
        "lambda_mode_retention_ratio_final": float(
            np.sqrt(np.mean(metric_delta["Lambda"][interior] ** 2))
            / max(np.sqrt(np.mean(delta_lambda[interior] ** 2)), np.finfo(float).tiny)
        ),
        "initial_constraint_witnesses_baseline": {
            key: _profile_summary(val, r, interior) for key, val in initial_constraints_base.items()
        },
        "initial_constraint_witnesses_off_mode": {
            key: _profile_summary(val, r, interior) for key, val in initial_constraints_off.items()
        },
        "initial_lambda_only_bssn_hamiltonian_response": {
            "identity": "delta_H_BSSN = X^2 * D4(delta_Lambda) at fixed metric and fixed algebraic fields",
            "observed_delta": _profile_summary(initial_bssn_h_delta, r, interior),
            "predicted_delta": _profile_summary(predicted_initial_bssn_h_delta, r, interior),
            "closure": _profile_summary(initial_bssn_h_closure, r, interior),
            "closure_max_abs": _max_abs(initial_bssn_h_closure),
            "metric_only_hamiltonian_delta_max_abs": _max_abs(
                initial_constraints_off["hamiltonian_metric_only"]
                - initial_constraints_base["hamiltonian_metric_only"]
            ),
        },
        "final_constraint_witnesses_baseline": {
            key: _profile_summary(val, r, interior) for key, val in constraints_base.items()
        },
        "final_constraint_witnesses_off_mode": {
            key: _profile_summary(val, r, interior) for key, val in constraints_off.items()
        },
        "final_constraint_delta_off_minus_baseline": {
            key: _profile_summary(val, r, interior) for key, val in constraint_delta.items()
        },
        "fixed_interior_constraint_delta_absmax": {
            key: _max_abs(val[interior]) for key, val in constraint_delta.items()
        },
        "finite": True,
    }
    profiles = {
        "r": r,
        "delta_lambda_input": delta_lambda,
        "exact_delta_lambda_prime": exact_delta_lambda_prime,
        "predicted_tf_ricci_response": predicted_tf_response,
        "predicted_Aa_L2_rate_delta": predicted_aa_rate,
        "actual_Aa_rate_delta": actual_aa_rate,
        "Aa_rate_closure_error": aa_rate_error,
        **{f"delta_{key}_after_step": val for key, val in metric_delta.items()},
        **{f"constraint_delta_{key}": val for key, val in constraint_delta.items()},
        **{f"initial_constraint_baseline_{key}": val for key, val in initial_constraints_base.items()},
        **{f"initial_constraint_off_mode_{key}": val for key, val in initial_constraints_off.items()},
        "initial_bssn_hamiltonian_delta_off_minus_baseline": initial_bssn_h_delta,
        "predicted_initial_bssn_hamiltonian_delta_from_lambda": predicted_initial_bssn_h_delta,
        "initial_bssn_hamiltonian_response_closure": initial_bssn_h_closure,
        **{f"final_constraint_baseline_{key}": val for key, val in constraints_base.items()},
        **{f"final_constraint_off_mode_{key}": val for key, val in constraints_off.items()},
    }
    return report, profiles


def run_resolution(resolution: int, vacuum) -> tuple[dict, dict]:
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution, r_max=R_MAX, amplitude=0.01,
        width=7.0, D_amplitude=1.0e-10, include_radiation=True,
    )
    before = _state_digest(state)
    reports, profiles = [], {}
    for dt_index, factor in enumerate(DT_OVER_DR):
        dt = factor * float(state.grid.dr)
        for wavelength in WAVELENGTHS:
            row, arrays = _run_pair(kernel, state, vacuum, wavelength, dt)
            reports.append(row)
            key = f"N{resolution}_dt{dt_index}_lambda{wavelength:g}".replace(".", "p")
            profiles.update({f"{key}_{name}": val for name, val in arrays.items()})
    if _state_digest(state) != before:
        raise RuntimeError(f"N={resolution}: initialized state changed during coupled observer")
    return {
        "resolution": int(resolution), "dr": float(state.grid.dr),
        "initial_state_sha256": before, "initialized_state_unchanged": _state_digest(state) == before,
        "cases": reports,
        "all_cases_finite": all(case["finite"] for case in reports),
        "all_input_states_unchanged": all(
            case["original_baseline_state_unchanged"]
            and case["baseline_input_state_unchanged"]
            and case["off_mode_input_state_unchanged"] for case in reports
        ),
    }, profiles


def main() -> int:
    provenance = _verify_source_pins()
    if provenance.get("source_commit") != SOURCE_COMMIT:
        raise RuntimeError("engine source anchor differs from declared pin")
    if provenance.get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError("vendor pin differs from declared pin")
    if not provenance.get("engine_vendor_dependency_diff_empty"):
        raise RuntimeError("production engine/vendor/dependency diff is not empty")
    if not all(item.get("matches") for item in provenance.get("source_blobs", {}).values()):
        raise RuntimeError("one or more pinned engine source blobs failed")
    _, vacuum, _ = vendor_modules()
    rows, all_profiles = [], {}
    for resolution in RESOLUTIONS:
        row, profiles = run_resolution(resolution, vacuum)
        rows.append(row)
        all_profiles.update(profiles)
    if not all(row["initialized_state_unchanged"] for row in rows):
        raise RuntimeError("one or more initialized state digests changed")
    if not all(row["all_cases_finite"] and row["all_input_states_unchanged"] for row in rows):
        raise FloatingPointError("one or more coupled propagation cases failed integrity checks")
    by_key = {}
    for row in rows:
        for case in row["cases"]:
            by_key.setdefault((row["resolution"], case["wavelength"]), {})[case["dt_over_dr"]] = case
    temporal = []
    for (resolution, wavelength), pair in sorted(by_key.items()):
        coarse_dt, fine_dt = sorted(pair.keys(), reverse=True)
        e0 = pair[coarse_dt]["actual_Aa_rate_delta_vs_frozen_L2_prediction_relative_rms_error"]
        e1 = pair[fine_dt]["actual_Aa_rate_delta_vs_frozen_L2_prediction_relative_rms_error"]
        temporal.append({
            "resolution": resolution, "wavelength": wavelength,
            "coarse_dt_over_dr": coarse_dt, "fine_dt_over_dr": fine_dt,
            "coarse_relative_rms_error": e0, "fine_relative_rms_error": e1,
            "fine_error_le_coarse_with_25pct_tolerance": bool(e1 <= 1.25 * e0 + 1.0e-8),
        })
    report = {
        "schema": "v55_coupled_lambda_off_constraint_step_propagation_v1",
        "status": "completed_research_only_one_step_off_constraint_propagation_control",
        "scientific_admission": False, "production_or_vendor_files_modified": False,
        "evolution_steps_per_arm": 1, "source_pins": provenance,
        "declared_source_pins": {"engine_source": SOURCE_COMMIT, "vendor": VENDOR_COMMIT},
        "configuration": {
            "resolutions": list(RESOLUTIONS), "r_max": R_MAX,
            "wavelengths": list(WAVELENGTHS), "dt_over_dr": list(DT_OVER_DR),
            "perturbation_amplitude": AMPLITUDE,
            "perturbation": "delta_Lambda = epsilon*sin(2*pi*r/wavelength), all other initial fields copied unchanged",
        },
        "resolutions": rows, "temporal_refinement": temporal,
        "all_initialized_states_unchanged": all(row["initialized_state_unchanged"] for row in rows),
        "all_profiles_finite": all(row["all_cases_finite"] for row in rows),
        "all_input_states_unchanged": all(row["all_input_states_unchanged"] for row in rows),
        "interpretation_limits": [
            "One accepted step per arm over two bounded dt values; not a full-horizon stability/convergence campaign.",
            "The odd Lambda perturbation is a manufactured off-constraint mode, not a constrained physical initial datum.",
            "The predicted Aa response is the frozen-state Ricci-bearing L2 response; finite-step differences include predictor/corrector and lapse-refresh effects.",
            "Constraint residual deltas are signed comparisons of perturbed and baseline outputs, not by themselves proof of stability or violation.",
            "No production operator change is made and no chain-rule substitution is tested in the actual step.",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    csv_path = OUT / "coupled_step_witnesses.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "N", "dr", "wavelength", "points_per_wavelength", "dt", "dt_over_dr",
            "Aa_rate_relative_RMS_error_vs_frozen_L2", "actual_Aa_rate_RMS",
            "predicted_Aa_rate_RMS", "delta_Aa_absmax", "delta_K_absmax",
            "delta_a_absmax", "delta_b_absmax", "delta_X_absmax", "delta_Lambda_absmax",
            "Lambda_mode_retention_ratio", "delta_H_absmax", "delta_H_metric_only_absmax",
            "delta_M_absmax", "delta_connection_absmax", "delta_determinant_absmax",
            "final_H_off_absmax", "final_M_off_absmax", "final_connection_off_absmax",
        ])
        for row in rows:
            for case in row["cases"]:
                writer.writerow([
                    row["resolution"], row["dr"], case["wavelength"], case["points_per_wavelength"],
                    case["dt"], case["dt_over_dr"],
                    case["actual_Aa_rate_delta_vs_frozen_L2_prediction_relative_rms_error"],
                    case["actual_Aa_rate_delta_rms_fixed_interior"],
                    case["predicted_Aa_rate_delta_rms_fixed_interior"],
                    case["delta_Aa_final_max_abs"], case["delta_K_final_max_abs"],
                    case["delta_a_final_max_abs"], case["delta_b_final_max_abs"],
                    case["delta_X_final_max_abs"], case["delta_Lambda_final_max_abs"],
                    case["lambda_mode_retention_ratio_final"],
                    case["fixed_interior_constraint_delta_absmax"]["hamiltonian"],
                    case["fixed_interior_constraint_delta_absmax"]["hamiltonian_metric_only"],
                    case["fixed_interior_constraint_delta_absmax"]["momentum"],
                    case["fixed_interior_constraint_delta_absmax"]["connection"],
                    case["fixed_interior_constraint_delta_absmax"]["determinant"],
                    case["final_constraint_witnesses_off_mode"]["hamiltonian"]["full_absmax"],
                    case["final_constraint_witnesses_off_mode"]["momentum"]["full_absmax"],
                    case["final_constraint_witnesses_off_mode"]["connection"]["full_absmax"],
                ])
    npz_path = OUT / "profiles.npz"
    np.savez_compressed(npz_path, **all_profiles)
    (OUT / "SHA256SUMS.txt").write_text(
        f"{_sha256(report_path)}  report.json\n"
        f"{_sha256(csv_path)}  coupled_step_witnesses.csv\n"
        f"{_sha256(npz_path)}  profiles.npz\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "audit_head": provenance.get("audit_head"),
        "cases": [{
            "N": row["resolution"], "wavelength": case["wavelength"],
            "dt_over_dr": case["dt_over_dr"],
            "Aa_rate_error": case["actual_Aa_rate_delta_vs_frozen_L2_prediction_relative_rms_error"],
            "delta_Aa": case["delta_Aa_final_max_abs"],
            "delta_K": case["delta_K_final_max_abs"],
            "delta_Lambda": case["delta_Lambda_final_max_abs"],
            "delta_H": case["fixed_interior_constraint_delta_absmax"]["hamiltonian"],
            "delta_H_metric_only": case["fixed_interior_constraint_delta_absmax"]["hamiltonian_metric_only"],
            "delta_M": case["fixed_interior_constraint_delta_absmax"]["momentum"],
            "delta_connection": case["fixed_interior_constraint_delta_absmax"]["connection"],
        } for row in rows for case in row["cases"]],
        "initialized_state_unchanged": report["all_initialized_states_unchanged"],
        "all_profiles_finite": report["all_profiles_finite"],
        "all_input_states_unchanged": report["all_input_states_unchanged"],
        "errors": [],
    }, indent=2, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
