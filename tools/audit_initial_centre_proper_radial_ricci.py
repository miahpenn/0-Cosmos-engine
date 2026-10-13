#!/usr/bin/env python3
"""No-evolution audit of regular-centre geometry on production initial slices.

This observer compares the pinned BSSN Ricci source with an independent
spherical warped-product reconstruction from areal radius and proper radial
distance. It intentionally reports discrepancies without changing the engine
or imposing a scientific-admission threshold.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

SOURCE_COMMIT = "29d828b4d72b017c12a20b6dc90549f2511f025c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_BLOBS = {
    "engine/production_kernel.py": "2954e7a79656849b5f088053ea0f4ac97d4bb443",
    "engine/v55_initial.py": "4cfbfe4be6ac3ce6ed492f4e435186997d94116a",
    "engine/v55_pirk_adapter.py": "f671b838d9aefd8926702d9c2235209f52cd6815",
    "engine/cmc_gauge.py": "b1e6b3e95ff37702509e1b38ec6401d10e2482c3",
    "engine/v55_matter.py": "0f71164494e534ca4226c56c10058845d9a14e1d",
    "engine/stress_energy.py": "c31c647f74fab6893254b3dea1366eada15b3110",
}
RESOLUTIONS = (60, 80, 120)
R_MAX = 80.0
CENTRE_CELLS = 6
OUT = Path("results/initial-centre-proper-radial-ricci")


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def _arr(values) -> np.ndarray:
    return np.asarray(values, dtype=float).copy()


def _max_abs(values) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.max(np.abs(values))) if values.size else 0.0


def _finite(values) -> bool:
    return bool(np.all(np.isfinite(np.asarray(values, dtype=float))))


def _state_digest(state) -> str:
    digest = hashlib.sha256()
    arrays = {}
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        arrays["geometry." + name] = getattr(state.geometry, name)
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        arrays["scalars." + name] = getattr(state.scalars, name)
    for sector in ("dark_matter", "baryons", "radiation"):
        item = getattr(state.matter, sector)
        for name in ("rest", "energy_t", "momentum_r"):
            arrays["matter." + sector + "." + name] = getattr(item, name)
    for name in sorted(arrays):
        values = np.ascontiguousarray(np.asarray(arrays[name]))
        digest.update(name.encode("utf-8"))
        digest.update(str(values.dtype).encode("ascii"))
        digest.update(repr(values.shape).encode("ascii"))
        digest.update(values.tobytes())
    return digest.hexdigest()


def proper_radial_ricci(grid, areal_radius, proper_factor):
    """Compute mixed Ricci eigenvalues from dl=f(r)dr and R_areal(r).

    For ds^2 = dl^2 + R(l)^2 dOmega^2:
      Ric^l_l = -2 R_ll/R
      Ric^theta_theta = (1 - R_l^2 - R R_ll)/R^2.
    The construction does not use the evolved conformal connection Lambda or
    the vendor's BSSN Ricci brackets. Derivatives use the pinned grid operator
    with odd parity for R and even parity for dl/dr.
    """
    R = _arr(areal_radius)
    f = _arr(proper_factor)
    if R.shape != np.asarray(grid.centers).shape or f.shape != R.shape:
        raise ValueError("areal radius and proper factor must match the grid")
    if np.any(R <= 0.0) or np.any(f <= 0.0):
        raise ValueError("areal radius and proper factor must be positive")
    R_r = grid.cell_derivative_fourth(R, parity=-1)
    R_rr = grid.cell_second_derivative_fourth(R, parity=-1)
    f_r = grid.cell_derivative_fourth(f, parity=1)
    R_l = R_r / f
    R_ll = (R_rr - R_r * f_r / f) / (f * f)
    Ric_l_l = -2.0 * R_ll / R
    Ric_theta_theta = (1.0 - R_l * R_l - R * R_ll) / (R * R)
    tf = (2.0 / 3.0) * (Ric_l_l - Ric_theta_theta)
    return {
        "R_r": R_r,
        "R_rr": R_rr,
        "f_r": f_r,
        "R_l": R_l,
        "R_ll": R_ll,
        "Ric_l_l": Ric_l_l,
        "Ric_theta_theta": Ric_theta_theta,
        "tf_ricci": tf,
    }


def even_taylor_fit(r, values, count=CENTRE_CELLS):
    """Fit c0 + c2*r^2 + c4*r^4 on the first cell-centred samples."""
    r = _arr(r)
    values = _arr(values)
    if r.shape != values.shape or count < 3 or r.size < count:
        raise ValueError("invalid centre profile for even Taylor fit")
    rr = r[:count]
    yy = values[:count]
    scale = float(rr[-1])
    x = (rr / scale) ** 2
    design = np.column_stack((np.ones(count), x, x * x))
    coeff_scaled, *_ = np.linalg.lstsq(design, yy, rcond=None)
    fitted = design @ coeff_scaled
    residual = yy - fitted
    physical = [
        float(coeff_scaled[0]),
        float(coeff_scaled[1] / (scale * scale)),
        float(coeff_scaled[2] / (scale**4)),
    ]
    return {
        "basis": "c0 + c2*r^2 + c4*r^4",
        "fit_cells": int(count),
        "c0": physical[0],
        "c2": physical[1],
        "c4": physical[2],
        "max_abs_residual": _max_abs(residual),
        "max_abs_profile": _max_abs(yy),
        "signed_residuals": [float(x) for x in residual],
        "fitted_values": [float(x) for x in fitted],
    }


def _verify_source_pins() -> dict:
    source_head = _git("rev-parse", "HEAD")
    source_vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if source_vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor pin mismatch: expected {VENDOR_COMMIT}, found {source_vendor}")
    pins = {}
    for path, expected in SOURCE_BLOBS.items():
        baseline_blob = _git("rev-parse", f"{SOURCE_COMMIT}:{path}")
        head_blob = _git("rev-parse", f"HEAD:{path}")
        pins[path] = {
            "expected_blob": expected,
            "baseline_blob": baseline_blob,
            "audit_head_blob": head_blob,
            "matches": baseline_blob == expected and head_blob == expected,
        }
        if not pins[path]["matches"]:
            raise RuntimeError(f"source blob mismatch for {path}: {pins[path]}")
    changed = _git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--",
        "engine", "vendor", ".gitmodules", "requirements.txt",
    )
    if changed:
        raise RuntimeError("engine/vendor/dependency paths changed from source pin:\n" + changed)
    return {
        "audit_head": source_head,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": source_vendor,
        "source_blobs": pins,
        "engine_vendor_dependency_diff_empty": True,
    }


def _profile_row(resolution, index, r, arrays):
    row = {"resolution": resolution, "i": index, "r": float(r[index])}
    for name, values in arrays.items():
        row[name] = float(values[index])
    return row


def main() -> int:
    provenance = _verify_source_pins()

    from engine import v55_pirk_adapter as adapter
    from engine.production_kernel import V55ProductionKernel

    _, vacuum, _ = adapter.vendor_modules()
    OUT.mkdir(parents=True, exist_ok=True)
    reports = []
    csv_rows = []
    profiles: dict[str, np.ndarray] = {}
    problems = []

    for resolution in RESOLUTIONS:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=resolution,
            r_max=R_MAX,
            amplitude=0.01,
            width=7.0,
            D_amplitude=1.0e-10,
            include_radiation=True,
        )
        grid = state.grid
        g = state.geometry
        r = _arr(grid.centers)
        before = _state_digest(state)

        a, b, X = _arr(g.a), _arr(g.b), _arr(g.X)
        alpha, Lambda = _arr(g.alpha), _arr(g.Lambda)
        areal = r * np.sqrt(b) / X
        proper_factor = np.sqrt(a) / X
        direct = proper_radial_ricci(grid, areal, proper_factor)
        vendor_terms = vacuum.geometry_terms(grid, g)
        vendor_tf = _arr(vendor_terms["Rr_r"] - vendor_terms["R"] / 3.0)
        vendor_l2_tf = alpha * vendor_tf
        direct_l2_tf = alpha * direct["tf_ricci"]
        tf_delta = vendor_tf - direct["tf_ricci"]
        l2_delta = vendor_l2_tf - direct_l2_tf

        d1 = lambda v, parity=1: grid.cell_derivative_fourth(v, parity=parity)
        lambda_from_metric = (
            d1(a, 1) / (2.0 * a * a)
            - d1(b, 1) / (a * b)
            + 2.0 / r * (1.0 / b - 1.0 / a)
        )
        a_minus_b_over_r2 = (a - b) / (r * r)
        lambda_over_r = Lambda / r

        metric_arrays = {
            "a": a,
            "b": b,
            "X": X,
            "alpha": alpha,
            "Lambda": Lambda,
            "a_minus_b_over_r2": a_minus_b_over_r2,
            "Lambda_over_r": lambda_over_r,
            "areal_radius": areal,
            "proper_factor": proper_factor,
            "vendor_tf_ricci": vendor_tf,
            "direct_tf_ricci": direct["tf_ricci"],
            "tf_ricci_difference": tf_delta,
            "vendor_l2_tf_ricci": vendor_l2_tf,
            "direct_l2_tf_ricci": direct_l2_tf,
            "l2_tf_ricci_difference": l2_delta,
            "R_r": direct["R_r"],
            "R_rr": direct["R_rr"],
            "f_r": direct["f_r"],
            "R_l": direct["R_l"],
            "R_ll": direct["R_ll"],
            "Ric_l_l": direct["Ric_l_l"],
            "Ric_theta_theta": direct["Ric_theta_theta"],
            "a_prime_over_r": d1(a, 1) / r,
            "b_prime_over_r": d1(b, 1) / r,
            "X_prime_over_r": d1(X, 1) / r,
            "alpha_prime_over_r": d1(alpha, 1) / r,
            "lambda_metric_identity_residual": Lambda - lambda_from_metric,
            "areal_radius_minus_coordinate_radius": areal - r,
        }

        taylor = {}
        for name in ("a", "b", "X", "alpha", "a_minus_b_over_r2", "Lambda_over_r"):
            taylor[name] = even_taylor_fit(r, metric_arrays[name])

        identities = {
            "a_times_b_squared_minus_one_max_abs": _max_abs(a * b * b - 1.0),
            "b_minus_X_squared_max_abs": _max_abs(b - X * X),
            "a_minus_X_minus_four_max_abs": _max_abs(a - X ** (-4.0)),
            "Lambda_metric_reconstruction_max_abs": _max_abs(Lambda - lambda_from_metric),
            "areal_radius_minus_r_max_abs": _max_abs(areal - r),
        }

        finite_names = list(metric_arrays)
        all_finite = all(_finite(metric_arrays[name]) for name in finite_names)
        after = _state_digest(state)
        unchanged = before == after
        if not all_finite:
            problems.append(f"N={resolution}: non-finite profile in initialized-state audit")
        if not unchanged:
            problems.append(f"N={resolution}: input state mutated by observer")

        first = {
            "i": 0,
            "r": float(r[0]),
            "vendor_tf_ricci": float(vendor_tf[0]),
            "direct_tf_ricci": float(direct["tf_ricci"][0]),
            "signed_difference": float(tf_delta[0]),
            "vendor_l2_tf_ricci": float(vendor_l2_tf[0]),
            "direct_l2_tf_ricci": float(direct_l2_tf[0]),
            "l2_signed_difference": float(l2_delta[0]),
            "Lambda_over_r": float(lambda_over_r[0]),
            "a_minus_b_over_r2": float(a_minus_b_over_r2[0]),
        }
        centre = []
        for i in range(min(CENTRE_CELLS, len(r))):
            centre.append({
                "i": i,
                "r": float(r[i]),
                "a_minus_b_over_r2": float(a_minus_b_over_r2[i]),
                "Lambda_over_r": float(lambda_over_r[i]),
                "a_prime_over_r": float(metric_arrays["a_prime_over_r"][i]),
                "b_prime_over_r": float(metric_arrays["b_prime_over_r"][i]),
                "X_prime_over_r": float(metric_arrays["X_prime_over_r"][i]),
                "alpha_prime_over_r": float(metric_arrays["alpha_prime_over_r"][i]),
                "vendor_tf_ricci": float(vendor_tf[i]),
                "direct_tf_ricci": float(direct["tf_ricci"][i]),
                "signed_difference": float(tf_delta[i]),
                "vendor_l2_tf_ricci": float(vendor_l2_tf[i]),
                "direct_l2_tf_ricci": float(direct_l2_tf[i]),
                "l2_signed_difference": float(l2_delta[i]),
                "Ric_l_l": float(direct["Ric_l_l"][i]),
                "Ric_theta_theta": float(direct["Ric_theta_theta"][i]),
            })

        denom = max(
            _max_abs(vendor_tf[:CENTRE_CELLS]),
            _max_abs(direct["tf_ricci"][:CENTRE_CELLS]),
            np.finfo(float).tiny,
        )
        report = {
            "resolution": resolution,
            "r_max": R_MAX,
            "dr": float(grid.dr),
            "evolution_steps": 0,
            "coordinate_time": float(state.t),
            "proper_time": float(state.tau),
            "radiation_enabled": True,
            "input_state_sha256_before": before,
            "input_state_sha256_after": after,
            "input_state_unchanged": unchanged,
            "all_profiles_finite": all_finite,
            "algebraic_and_connection_identities": identities,
            "first_cell": first,
            "first_six_cells": centre,
            "ricci_comparison_first_six": {
                "vendor_tf_max_abs": _max_abs(vendor_tf[:CENTRE_CELLS]),
                "direct_tf_max_abs": _max_abs(direct["tf_ricci"][:CENTRE_CELLS]),
                "signed_difference_max_abs": _max_abs(tf_delta[:CENTRE_CELLS]),
                "relative_difference_to_larger_signal": float(
                    _max_abs(tf_delta[:CENTRE_CELLS]) / denom
                ),
                "vendor_l2_tf_max_abs": _max_abs(vendor_l2_tf[:CENTRE_CELLS]),
                "direct_l2_tf_max_abs": _max_abs(direct_l2_tf[:CENTRE_CELLS]),
                "l2_signed_difference_max_abs": _max_abs(l2_delta[:CENTRE_CELLS]),
            },
            "even_taylor_fits_first_six_cells": taylor,
            "interpretation_limits": [
                "A successful observer workflow is not scientific admission.",
                "The direct warped-product reconstruction is independent of evolved Lambda and the vendor Ricci brackets, but shares the pinned radial grid derivative operator.",
                "Taylor fits summarize the first six cell-centred values; fit residuals are evidence, not a post-hoc pass/fail tolerance.",
                "A discrepancy localizes a representation/operator compatibility question but does not alone prove which equation is wrong.",
            ],
        }
        reports.append(report)

        prefix = f"N{resolution}_"
        for name, values in metric_arrays.items():
            profiles[prefix + name] = _arr(values)
        for name, result in taylor.items():
            profiles[prefix + name + "_taylor_residual_first6"] = _arr(result["signed_residuals"])
            for i, r_value in enumerate(r):
                csv_rows.append(_profile_row(resolution, i, r_value, metric_arrays))

    status = "completed_diagnostic_only" if not problems else "observer_integrity_failure"
    result = {
        "schema": "v55_initial_centre_proper_radial_ricci_v1",
        "status": status,
        "source_pins": provenance,
        "resolutions": reports,
        "errors": problems,
        "no_evolution": True,
        "production_or_vendor_files_modified": False,
        "scientific_admission": False,
    }

    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    npz_path = OUT / "profiles.npz"
    np.savez_compressed(npz_path, **profiles)
    csv_path = OUT / "profiles.csv"
    fieldnames = ["resolution", "i", "r"] + list(reports[0]["resolutions"][0].keys()) if False else None
    columns = ["resolution", "i", "r"] + list(next(iter(reports))["first_six_cells"][0].keys()) if False else None
    # Export the full signed profiles in a deterministic, flat table.
    all_names = list(reports and [
        "a", "b", "X", "alpha", "Lambda", "a_minus_b_over_r2", "Lambda_over_r",
        "areal_radius", "proper_factor", "vendor_tf_ricci", "direct_tf_ricci",
        "tf_ricci_difference", "vendor_l2_tf_ricci", "direct_l2_tf_ricci",
        "l2_tf_ricci_difference", "R_r", "R_rr", "f_r", "R_l", "R_ll",
        "Ric_l_l", "Ric_theta_theta", "a_prime_over_r", "b_prime_over_r",
        "X_prime_over_r", "alpha_prime_over_r", "lambda_metric_identity_residual",
        "areal_radius_minus_coordinate_radius",
    ])
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["resolution", "i", "r"] + all_names)
        for n in RESOLUTIONS:
            source = {key[len(f"N{n}_"):]: value for key, value in profiles.items()
                      if key.startswith(f"N{n}_") and not key.endswith("_taylor_residual_first6")}
            for i in range(source["a"].size):
                writer.writerow([n, i, float(i * (R_MAX / n) + 0.5 * (R_MAX / n))] +
                                [float(source[name][i]) for name in all_names])

    checksum_paths = (report_path, npz_path, csv_path)
    (OUT / "SHA256SUMS.txt").write_text(
        "\n".join(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}"
                  for path in checksum_paths) + "\n",
        encoding="utf-8",
    )
    (OUT / "SOURCE_COMMIT.txt").write_text(SOURCE_COMMIT + "\n", encoding="utf-8")
    (OUT / "VENDOR_COMMIT.txt").write_text(VENDOR_COMMIT + "\n", encoding="utf-8")

    print(json.dumps({
        "status": status,
        "head": provenance["audit_head"],
        "resolutions": [
            {
                "N": item["resolution"],
                "centre_vendor_tf": item["first_cell"]["vendor_tf_ricci"],
                "centre_direct_tf": item["first_cell"]["direct_tf_ricci"],
                "centre_difference": item["first_cell"]["signed_difference"],
                "first_six_relative_difference": item["ricci_comparison_first_six"]["relative_difference_to_larger_signal"],
                "state_unchanged": item["input_state_unchanged"],
            }
            for item in reports
        ],
        "errors": problems,
    }, indent=2, allow_nan=False), flush=True)
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
