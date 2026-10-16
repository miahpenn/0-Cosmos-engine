#!/usr/bin/env python3
"""Read-only audit of the proposed scalar-DM covariant-source correction.

This checks PR #280's exact candidate snapshot against independently assembled
Valencia/scalar exchange identities. It also records the finite-volume
face-flux + geometric-source + exchange-source budget for one frozen profile.
No production time step is taken and no production source is changed.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np

SOURCE_COMMIT = "1f142b0b2f3bc832a3030027f139d0d6c0f87095"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_BLOBS = {
    "engine/matter_system.py": "cd70b715a4b779adb0e933f4badccc57277a3dc8",
    "engine/scalar_system.py": "4c3a13b4e58fabbe06d8bb1316cef73d2db8e1f4",
    "engine/matter_rhs.py": "0153c232ee15b5a9d3dbac37ea858e4d2d730011",
    "engine/production_kernel.py": "b5de14ee1d9b7ffb954b8a08ff888c8953d722c2",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
}
BETA_DM = -0.04
KAPPA = 8.0 * math.pi
TOL = 5.0e-12


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def arr(name: str, x) -> np.ndarray:
    out = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError(f"{name} has a non-finite value")
    return out


def rms(x) -> float:
    x = arr("RMS profile", x)
    return float(np.sqrt(np.mean(x * x)))


def profile(x, r) -> dict:
    x = arr("summary profile", x)
    i = int(np.argmax(np.abs(x)))
    return {
        "cells": int(x.size), "rms": rms(x),
        "max_abs": float(abs(x[i])), "max_index": i,
        "max_radius": float(r[i]), "signed_at_max_abs": float(x[i]),
    }


def verify_provenance(root: Path) -> dict:
    commit = git(root, "rev-parse", "HEAD")
    vendor = git(root / "vendor" / "bb-palatini-unified-r0", "rev-parse", "HEAD")
    blobs = {p: git(root, "rev-parse", f"HEAD:{p}") for p in SOURCE_BLOBS}
    out = {
        "source_commit": commit, "declared_source_commit": SOURCE_COMMIT,
        "vendor_commit": vendor, "declared_vendor_commit": VENDOR_COMMIT,
        "source_blobs": blobs, "declared_source_blobs": SOURCE_BLOBS,
        "source_commit_matches": commit == SOURCE_COMMIT,
        "vendor_commit_matches": vendor == VENDOR_COMMIT,
        "source_blobs_match": blobs == SOURCE_BLOBS,
    }
    if not all((out["source_commit_matches"], out["vendor_commit_matches"], out["source_blobs_match"])):
        raise RuntimeError("candidate source provenance failed: " + json.dumps(out))
    return out


class AuditGrid:
    """Derivative-protocol adapter; paired RHS differencing cancels its stencil."""
    def __init__(self, centers):
        self.centers = arr("grid centers", centers)
        drs = np.diff(self.centers)
        if drs.size == 0 or not np.allclose(drs, drs[0], rtol=2e-14, atol=0.0):
            raise ValueError("AuditGrid requires a uniform grid")
        self.dr = float(drs[0])

    def cell_derivative_fourth(self, values, parity=1):
        # Source-term coefficient checks use differences of two identical RHS
        # calls; this derivative protocol adapter is not a convergence stencil.
        values = arr("AuditGrid D1", values)
        return np.gradient(values, self.dr, edge_order=2)

    def cell_second_derivative_fourth(self, values, parity=1):
        values = arr("AuditGrid D2", values)
        return np.gradient(np.gradient(values, self.dr, edge_order=2), self.dr, edge_order=2)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=Path("results/candidate-source-closure/report.json"))
    ap.add_argument("--resolution", type=int, default=160)
    ap.add_argument("--r-max", type=float, default=24.0)
    args = ap.parse_args()

    source_raw = os.environ.get("COSMOS_CANDIDATE_SOURCE_ROOT")
    if not source_raw:
        raise RuntimeError("COSMOS_CANDIDATE_SOURCE_ROOT must point to the exact PR #280 checkout")
    source_root = Path(source_raw).resolve()
    if not source_root.is_dir():
        raise RuntimeError(f"candidate source checkout not found: {source_root}")
    sys.path.insert(0, str(source_root))
    provenance = verify_provenance(source_root)

    from engine.matter_system import (
        BSSNMetricSlice, MetricDerivativeSet, Species, _metric_arrays,
        _species_rhs_components, initialize_dust, primitives,
    )
    from engine.scalar_system import ScalarFields, scalar_rhs_arrays

    n, r_max = int(args.resolution), float(args.r_max)
    if n < 32 or not math.isfinite(r_max) or r_max <= 0:
        raise ValueError("resolution must be >=32 and r_max positive/finite")
    dr = r_max / n
    r = (np.arange(n, dtype=float) + 0.5) * dr
    e = np.exp(-(r * r) / 100.0)
    alpha = 0.65 + 0.20 * (r * r) / (r * r + 9.0)
    beta_shift = 0.025 * r / (1.0 + r * r / 16.0)
    metric = BSSNMetricSlice(
        r=r,
        a=1.0 + 0.04 * e,
        b=1.0 + 0.012 * e,
        X=1.0 + 0.009 * e,
        alpha=alpha,
        beta=beta_shift,
        Aa=2.0e-4 * r * r * np.exp(-(r * r) / 80.0),
        K=-0.025 + 0.003 * e,
        Lambda=np.zeros_like(r),
        B=np.zeros_like(r),
    )
    grid = AuditGrid(r)
    geometry = SimpleNamespace(
        a=metric.a, b=metric.b, X=metric.X,
        alpha=metric.alpha, beta=metric.beta, K=metric.K,
    )
    fields = ScalarFields(
        S=0.010 * np.exp(-(r * r) / 28.0),
        PS=0.003 * np.exp(-(r * r) / 32.0),
        D=0.015 * np.exp(-(r * r) / 24.0),
        PD=0.004 * np.exp(-(r * r) / 30.0),
        phi=0.12 * np.exp(-(r * r) / 70.0),
        Pi=0.04 * np.exp(-(r * r) / 40.0),
    )
    scalar_base = scalar_rhs_arrays(grid, geometry, fields, beta_dm=BETA_DM, rho_dm=None)
    rho_seed = 0.008 + 0.012 * np.exp(-(r * r) / 36.0)
    v_seed = 0.018 * r / (1.0 + r * r / 16.0)
    metrics = _metric_arrays(metric)
    state = initialize_dust(metrics, rho_seed, v_seed)
    primitive = primitives(metrics, state, Species.DARK_MATTER)
    rho = arr("recovered DM primitive rho", [q.rho for q in primitive])
    scalar_full = scalar_rhs_arrays(grid, geometry, fields, beta_dm=BETA_DM, rho_dm=rho)
    dphi_t = arr("dphi_t=scalar phi RHS", scalar_base.phi)
    dphi_r = arr("D_r(phi)", grid.cell_derivative_fourth(fields.phi, parity=1))

    zeros = np.zeros(n, dtype=float)
    metric_derivatives = MetricDerivativeSet(
        time={},
        radial={"alpha": zeros.copy(), "beta": zeros.copy(), "rr": zeros.copy(), "thth": zeros.copy()},
    )
    components = _species_rhs_components(
        metric, metric_derivatives, state, Species.DARK_MATTER,
        dphi_t=dphi_t, dphi_r=dphi_r, beta_dm=BETA_DM,
    )
    face_flux, source_e, source_s, q_e, q_s, q_rest = components
    sqrt_gamma = arr("sqrt gamma", [m.sqrt_gamma for m in metrics])
    Q_t = arr("Q_t", BETA_DM * rho * dphi_t)
    Q_r = arr("Q_r", BETA_DM * rho * dphi_r)

    expected_q_e = arr("expected qE", -sqrt_gamma * (Q_t - beta_shift * Q_r))
    expected_q_s = arr("expected qSr", alpha * sqrt_gamma * Q_r)
    pi_normal = arr("normal scalar derivative reconstructed from dphi_t", (dphi_t - beta_shift * dphi_r) / alpha)
    lorentz = arr("DM Lorentz factors", [p.lorentz() for p in primitive])
    v = arr("recovered DM velocities", [p.v_r for p in primitive])
    expected_q_rest = arr(
        "expected variable-mass rest-current source",
        -alpha * sqrt_gamma * BETA_DM * rho * lorentz * (pi_normal + v * dphi_r),
    )
    delta_Pi = arr("scalar Pi coupling increment", scalar_full.Pi - scalar_base.Pi)
    expected_delta_Pi = arr("expected covariant scalar source increment", alpha * BETA_DM * rho)

    resid_e = arr("qE residual", q_e - expected_q_e)
    resid_s = arr("qSr residual", q_s - expected_q_s)
    resid_rest = arr("qD residual", q_rest - expected_q_rest)
    resid_pi = arr("scalar Pi residual", delta_Pi - expected_delta_Pi)

    scalar_energy_rate = arr("scalar exchange energy rate", fields.Pi * delta_Pi / KAPPA)
    fluid_energy_rate = arr("fluid exchange energy rate", q_e / sqrt_gamma / KAPPA)
    total_energy_rate = arr("scalar+fluid source energy rate", scalar_energy_rate + fluid_energy_rate)
    scalar_momentum_rate = arr("scalar exchange radial-current rate", -dphi_r * delta_Pi / KAPPA)
    fluid_momentum_rate = arr("fluid exchange radial-current rate", q_s / sqrt_gamma / KAPPA)
    total_momentum_rate = arr("scalar+fluid source radial-current rate", scalar_momentum_rate + fluid_momentum_rate)

    # Full frozen-state conservative RHS as face flux divergence + metric source
    # + exchange source, keeping each contribution separately auditable.
    rhs_rest = arr("Valencia rest RHS", -np.diff(face_flux[:, 0]) / dr + q_rest)
    rhs_energy = arr("Valencia energy RHS", -np.diff(face_flux[:, 1]) / dr + source_e + q_e)
    rhs_momentum = arr("Valencia momentum RHS", -np.diff(face_flux[:, 2]) / dr + source_s + q_s)
    fv_rows = {}
    for name, k, rhs, sources in (
        ("rest", 0, rhs_rest, (q_rest,)),
        ("energy", 1, rhs_energy, (source_e, q_e)),
        ("momentum", 2, rhs_momentum, (source_s, q_s)),
    ):
        boundary = -(float(face_flux[-1, k]) - float(face_flux[0, k]))
        source_integral = dr * sum(float(np.sum(item)) for item in sources)
        lhs = dr * float(np.sum(rhs))
        fv_rows[name] = {
            "integrated_rhs": lhs,
            "outer_minus_inner_flux_with_sign": boundary,
            "integrated_local_sources": source_integral,
            "telescoping_identity_residual": lhs - boundary - source_integral,
            "rhs_rms": rms(rhs),
            "flux_divergence_rms": rms(-np.diff(face_flux[:, k]) / dr),
            "source_rms": [rms(item) for item in sources],
            "inner_face_flux": float(face_flux[0, k]),
            "outer_face_flux": float(face_flux[-1, k]),
        }

    def finite_error(name, x, scale):
        x = arr(name, x)
        denominator = max(rms(scale), np.finfo(float).tiny)
        return {
            "residual": profile(x, r),
            "relative_rms": rms(x) / denominator,
            "closure_passed": bool(float(np.max(np.abs(x))) <= TOL),
        }

    report = {
        "schema": "issue4_candidate_scalar_dm_exchange_source_closure_v1",
        "status": "completed_diagnostic_only",
        "scientific_admission": False,
        "total_time_evolved_bianchi_gate_closed": False,
        "production_steps_advanced": 0,
        "production_or_vendor_files_modified": False,
        "candidate_provenance": provenance,
        "configuration": {
            "resolution": n, "r_max": r_max, "dr": dr, "beta_dm": BETA_DM,
            "alpha_min": float(np.min(alpha)), "alpha_max": float(np.max(alpha)),
            "nonzero_shift": True, "moving_dust": True, "nonzero_scalar_gradient": True,
        },
        "exchange_source_checks": {
            "energy_qE": finite_error("qE residual", resid_e, expected_q_e),
            "radial_momentum_qSr": finite_error("qSr residual", resid_s, expected_q_s),
            "variable_mass_rest_current_qD": finite_error("qD residual", resid_rest, expected_q_rest),
            "scalar_pi_covariant_source": finite_error("Pi residual", resid_pi, expected_delta_Pi),
        },
        "combined_local_exchange_budget": {
            "scalar_energy_rate": profile(scalar_energy_rate, r),
            "fluid_energy_rate": profile(fluid_energy_rate, r),
            "total_energy_rate": profile(total_energy_rate, r),
            "energy_exchange_closure_passed": bool(float(np.max(np.abs(total_energy_rate))) <= TOL),
            "scalar_radial_current_rate": profile(scalar_momentum_rate, r),
            "fluid_radial_current_rate": profile(fluid_momentum_rate, r),
            "total_radial_current_rate": profile(total_momentum_rate, r),
            "radial_current_exchange_closure_passed": bool(float(np.max(np.abs(total_momentum_rate))) <= TOL),
        },
        "finite_volume_budget": {
            "interpretation": "This is exact discrete accounting of the source candidate's own face flux, metric source and exchange arrays, not a standalone continuum-conservation proof.",
            "channels": fv_rows,
            "all_telescoping_identities_close": all(
                abs(row["telescoping_identity_residual"]) <= 2.0e-10 for row in fv_rows.values()
            ),
        },
        "gate_result": {
            "all_four_source_channels_close": False,
            "scalar_fluid_local_energy_exchange_closes": bool(float(np.max(np.abs(total_energy_rate))) <= TOL),
            "scalar_fluid_local_radial_current_exchange_closes": bool(float(np.max(np.abs(total_momentum_rate))) <= TOL),
            "finite_volume_telescope_closes": all(
                abs(row["telescoping_identity_residual"]) <= 2.0e-10 for row in fv_rows.values()
            ),
            "scientific_admission": False,
        },
    }
    checks = report["exchange_source_checks"]
    report["gate_result"]["all_four_source_channels_close"] = all(v["closure_passed"] for v in checks.values())
    report["gate_result"]["all_diagnostic_controls_pass"] = bool(
        report["gate_result"]["all_four_source_channels_close"]
        and report["gate_result"]["scalar_fluid_local_energy_exchange_closes"]
        and report["gate_result"]["scalar_fluid_local_radial_current_exchange_closes"]
        and report["gate_result"]["finite_volume_telescope_closes"]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "source_pins_match": True,
        "all_source_channels_close": report["gate_result"]["all_four_source_channels_close"],
        "local_energy_exchange_closes": report["gate_result"]["scalar_fluid_local_energy_exchange_closes"],
        "local_momentum_exchange_closes": report["gate_result"]["scalar_fluid_local_radial_current_exchange_closes"],
        "finite_volume_telescope_closes": report["gate_result"]["finite_volume_telescope_closes"],
        "scientific_admission": False,
        "report": str(args.output),
    }, indent=2))
    if not report["gate_result"]["all_diagnostic_controls_pass"]:
        raise RuntimeError("candidate-source closure check failed; inspect the retained report")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
