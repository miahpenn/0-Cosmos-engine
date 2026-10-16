#!/usr/bin/env python3
"""Bounded, read-only Valencia stage budget for Issue #4 Gate C.

This audit freezes a smooth spherical profile and evaluates the pinned
production matter RHS. It reports signed face-flux divergences, metric sources,
and exchange sources separately, independently reconstructs the Valencia
metric-source contraction, and tests the conditional covariant exchange map.

It does NOT evolve a production trajectory and does NOT claim that the
finite-volume telescope is an independent proof of div(T_total)=0. A passing
run is an implementation/source discriminator only; total stress conservation,
the Brown/Bianchi propagation oracle, physical S/D normalization, and
scientific admission remain separate gates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
SOURCE_ANCHOR = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
VENDOR_ANCHOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_BLOBS = {
    "engine/production_kernel.py": "9685e76bb5a6f26d6bf383e240de24f37523f265",
    "engine/matter_system.py": "e52a117690253ddad7b24d116145a65474a351d2",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
}
BETA_DM = -0.04
ABS_TOL = 5.0e-12
REL_TOL = 5.0e-12


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def finite(name: str, values) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return out


def rms(values) -> float:
    values = finite("RMS input", values)
    return float(np.sqrt(np.mean(values * values)))


def summarize(values, radii) -> dict:
    arr = finite("summary profile", values)
    r = finite("summary radii", radii)
    i = int(np.argmax(np.abs(arr)))
    return {
        "rms": rms(arr),
        "max_abs": float(abs(arr[i])),
        "signed_at_max_abs": float(arr[i]),
        "max_index": i,
        "max_radius": float(r[i]),
        "signed_first_four_cells": [float(x) for x in arr[:4]],
        "signed_at_r_near_1": float(arr[int(np.argmin(abs(r - 1.0)))]),
        "signed_at_r_near_3": float(arr[int(np.argmin(abs(r - 3.0)))]),
        "signed_at_r_near_6": float(arr[int(np.argmin(abs(r - 6.0)))]),
    }


def provenance() -> dict:
    actual_blobs = {
        path: git("rev-parse", f"HEAD:{path}") for path in SOURCE_BLOBS
    }
    vendor = git(
        "-C", "vendor/bb-palatini-unified-r0", "rev-parse", "HEAD"
    )
    mismatches = {
        path: {"expected": expected, "actual": actual_blobs[path]}
        for path, expected in SOURCE_BLOBS.items()
        if actual_blobs[path] != expected
    }
    if mismatches:
        raise RuntimeError(
            "production source blobs differ from the declared baseline: "
            + json.dumps(mismatches, sort_keys=True)
        )
    if vendor != VENDOR_ANCHOR:
        raise RuntimeError(
            f"vendor pin mismatch: expected {VENDOR_ANCHOR}, got {vendor}"
        )
    return {
        "repository_commit": git("rev-parse", "HEAD"),
        "production_source_anchor": SOURCE_ANCHOR,
        "vendor_commit": vendor,
        "source_blobs": actual_blobs,
        "source_blobs_match": True,
        "vendor_commit_matches": True,
    }


def independent_metric_source(metric, q, K: float, Aa: float,
                              dr_alpha: float, dr_beta: float,
                              dr_grr: float, dr_gthth: float) -> tuple[float, float]:
    """Reconstruct the covariant Valencia metric-source contraction.

    This is independently expressed in the audit rather than obtained from
    production's _valencia_source helper. Formula follows the declared
    Eulerian energy/momentum source convention; it does not validate the
    physical normalization of every matter sector.
    """
    W = q.lorentz()
    h = q.rho + q.pressure
    E = h * W * W - q.pressure
    S_r = h * W * W * metric.gamma_rr * q.v_r
    S_up = metric.gamma_rr_inv * S_r
    Srr = (
        h * W * W * q.v2 + q.pressure
    ) * metric.gamma_rr_inv
    Sthth = q.pressure * metric.gamma_thth_inv
    Krr = metric.gamma_rr * (K / 3.0 + Aa)
    Kthth = metric.gamma_thth * (K / 3.0 - Aa / 2.0)
    energy = metric.sqrt_gamma * (
        metric.alpha * (Krr * Srr + 2.0 * Kthth * Sthth)
        - S_up * dr_alpha
    )
    momentum = metric.sqrt_gamma * (
        -E * dr_alpha
        + S_r * dr_beta
        + 0.5 * metric.alpha * (
            Srr * dr_grr + 2.0 * Sthth * dr_gthth
        )
    )
    return float(energy), float(momentum)


def build_profile(n: int, r_max: float):
    from engine.production_kernel import V55ProductionKernel
    from engine.matter_system import BSSNMetricSlice

    # Use the production grid object for its actual radial/parity derivative
    # operators. Initialization only; no RHS step or trajectory is advanced.
    kernel = V55ProductionKernel()
    seed = kernel.initialize(
        resolution=n, r_max=r_max, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True,
    )
    grid = seed.grid
    r = finite("grid centers", grid.centers)
    e = np.exp(-(r * r) / 100.0)
    metric = BSSNMetricSlice(
        r=r.copy(),
        a=1.0 + 0.04 * e,
        b=1.0 + 0.012 * e,
        X=1.0 + 0.009 * e,
        alpha=0.65 + 0.20 * (r * r) / (r * r + 9.0),
        beta=0.025 * r / (1.0 + r * r / 16.0),
        Aa=2.0e-4 * r * r * np.exp(-(r * r) / 80.0),
        K=-0.025 + 0.003 * e,
        Lambda=np.zeros_like(r),
        B=np.zeros_like(r),
    )
    geometry = type("GeometryView", (), {})()
    for name in ("a", "b", "X", "alpha", "beta", "K"):
        setattr(geometry, name, getattr(metric, name))
    return grid, metric, geometry


def rhs_as_arrays(rhs) -> dict[str, np.ndarray]:
    return {
        "rest": finite("RHS rest", rhs.rest),
        "energy_t": finite("RHS energy", rhs.energy_t),
        "momentum_r": finite("RHS momentum", rhs.momentum_r),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=160)
    parser.add_argument("--r-max", type=float, default=24.0)
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/issue4-valencia-total-stage-budget/report.json"),
    )
    args = parser.parse_args()
    if args.resolution < 64 or not math.isfinite(args.r_max) or args.r_max <= 0:
        raise ValueError("resolution must be >=64 and r-max positive and finite")

    source_provenance = provenance()
    from engine.matter_system import (
        Species, _metric_arrays, _species_rhs_components,
        geometry_metric_derivatives, initialize_dust, initialize_radiation,
        primitives,
    )
    from engine.matter_rhs import species_rhs
    from engine.scalar_system import ScalarFields, scalar_rhs_arrays

    grid, metric, geometry = build_profile(args.resolution, args.r_max)
    r = finite("radii", metric.r)
    dr = float(r[1] - r[0])
    # Production derivative operator and its declared origin parity.
    derivs = geometry_metric_derivatives(metric, grid.cell_derivative_fourth)
    metric_objects = _metric_arrays(metric)
    n = r.size

    fields = ScalarFields(
        S=0.010 * np.exp(-(r * r) / 28.0),
        PS=0.003 * np.exp(-(r * r) / 32.0),
        D=0.015 * np.exp(-(r * r) / 24.0),
        PD=0.004 * np.exp(-(r * r) / 30.0),
        phi=0.12 * np.exp(-(r * r) / 70.0),
        Pi=0.04 * np.exp(-(r * r) / 40.0),
    )
    scalar_uncoupled = scalar_rhs_arrays(
        grid, geometry, fields, beta_dm=BETA_DM, rho_dm=None
    )
    dphi_t = finite("coordinate d_t(phi)", scalar_uncoupled.phi)
    dphi_r = finite(
        "D_r(phi)", grid.cell_derivative_fourth(fields.phi, parity=1)
    )

    rho_dm_seed = 0.008 + 0.012 * np.exp(-(r * r) / 36.0)
    v_dm = 0.018 * r / (1.0 + r * r / 16.0)
    rho_b_seed = 0.004 + 0.003 * np.exp(-(r * r) / 20.0)
    v_b = -0.012 * r / (1.0 + r * r / 12.0)
    rho_rad_seed = 0.0007 + 0.0002 * np.exp(-(r * r) / 18.0)
    v_rad = 0.006 * r / (1.0 + r * r / 15.0)
    initial_states = {
        "dark_matter": (
            Species.DARK_MATTER, initialize_dust(metric_objects, rho_dm_seed, v_dm)
        ),
        "baryons": (
            Species.BARYON, initialize_dust(metric_objects, rho_b_seed, v_b)
        ),
        "radiation": (
            Species.RADIATION, initialize_radiation(metric_objects, rho_rad_seed, v_rad)
        ),
    }

    species_reports = {}
    signed_profiles = {"r": [float(x) for x in r]}
    controls = {
        "all_profiles_finite": True,
        "independent_metric_sources_match": True,
        "production_rhs_matches_signed_component_assembly": True,
        "finite_volume_telescopes_close": True,
    }
    max_metric_source_abs = 0.0
    max_component_rhs_abs = 0.0
    max_telescoping_abs = 0.0

    for name, (species, state) in initial_states.items():
        prim = primitives(metric_objects, state, species)
        components = _species_rhs_components(
            metric, derivs, state, species,
            dphi_t=dphi_t if species is Species.DARK_MATTER else None,
            dphi_r=dphi_r if species is Species.DARK_MATTER else None,
            beta_dm=BETA_DM,
        )
        flux, source_e, source_s, q_e, q_s = components
        flux = finite(f"{name} face flux", flux)
        source_e = finite(f"{name} energy metric source", source_e)
        source_s = finite(f"{name} momentum metric source", source_s)
        q_e = finite(f"{name} energy exchange", q_e)
        q_s = finite(f"{name} momentum exchange", q_s)

        expected_e = np.empty(n, dtype=float)
        expected_s = np.empty(n, dtype=float)
        for i, (m, q) in enumerate(zip(metric_objects, prim)):
            expected_e[i], expected_s[i] = independent_metric_source(
                m, q, float(metric.K[i]), float(metric.Aa[i]),
                float(derivs.radial["alpha"][i]),
                float(derivs.radial["beta"][i]),
                float(derivs.radial["rr"][i]),
                float(derivs.radial["thth"][i]),
            )
        metric_e_residual = source_e - expected_e
        metric_s_residual = source_s - expected_s
        metric_scale = max(
            1.0, float(np.max(np.abs(expected_e))),
            float(np.max(np.abs(expected_s))),
        )
        metric_closed = (
            float(np.max(np.abs(metric_e_residual))) <= ABS_TOL * metric_scale
            and float(np.max(np.abs(metric_s_residual))) <= ABS_TOL * metric_scale
        )
        controls["independent_metric_sources_match"] &= bool(metric_closed)
        max_metric_source_abs = max(
            max_metric_source_abs,
            float(np.max(np.abs(metric_e_residual))),
            float(np.max(np.abs(metric_s_residual))),
        )

        flux_div = {
            "rest": finite(f"{name} rest flux divergence", -np.diff(flux[:, 0]) / dr),
            "energy_t": finite(f"{name} energy flux divergence", -np.diff(flux[:, 1]) / dr),
            "momentum_r": finite(f"{name} momentum flux divergence", -np.diff(flux[:, 2]) / dr),
        }
        source_terms = {
            "rest": np.zeros(n, dtype=float),
            "energy_t": source_e + q_e,
            "momentum_r": source_s + q_s,
        }
        rhs_expected = {
            channel: finite(f"{name} assembled {channel}", flux_div[channel] + source_terms[channel])
            for channel in ("rest", "energy_t", "momentum_r")
        }
        rhs_direct = rhs_as_arrays(species_rhs(
            metric, derivs, state, species,
            dphi_t=dphi_t if species is Species.DARK_MATTER else None,
            dphi_r=dphi_r if species is Species.DARK_MATTER else None,
            beta_dm=BETA_DM,
        ))
        rhs_residuals = {
            channel: finite(f"{name} direct-vs-assembled {channel}", rhs_direct[channel] - rhs_expected[channel])
            for channel in rhs_expected
        }
        rhs_max = max(float(np.max(np.abs(x))) for x in rhs_residuals.values())
        controls["production_rhs_matches_signed_component_assembly"] &= (
            rhs_max <= ABS_TOL * max(
                1.0, *(float(np.max(np.abs(x))) for x in rhs_direct.values())
            )
        )
        max_component_rhs_abs = max(max_component_rhs_abs, rhs_max)

        fv_balance = {}
        for channel, column in (("rest", 0), ("energy_t", 1), ("momentum_r", 2)):
            local = source_terms[channel]
            rhs = rhs_expected[channel]
            lhs_int = dr * float(np.sum(rhs))
            boundary = -(float(flux[-1, column]) - float(flux[0, column]))
            source_int = dr * float(np.sum(local))
            residual = lhs_int - boundary - source_int
            fv_balance[channel] = {
                "integrated_rhs": lhs_int,
                "negative_outer_plus_inner_face_flux": boundary,
                "integrated_metric_plus_exchange_source": source_int,
                "signed_telescope_residual": residual,
                "telescope_closed": abs(residual) <= 5.0e-11 * max(
                    1.0, abs(lhs_int), abs(boundary), abs(source_int)
                ),
            }
            max_telescoping_abs = max(max_telescoping_abs, abs(residual))
            controls["finite_volume_telescopes_close"] &= bool(
                fv_balance[channel]["telescope_closed"]
            )

        # Store signed per-cell data, retaining source/flux separation for audit.
        signed_profiles[name] = {
            "face_flux_rest_energy_momentum": flux.tolist(),
            "flux_divergence": {k: v.tolist() for k, v in flux_div.items()},
            "metric_source_energy": source_e.tolist(),
            "metric_source_momentum": source_s.tolist(),
            "exchange_source_energy": q_e.tolist(),
            "exchange_source_momentum": q_s.tolist(),
            "assembled_rhs": {k: v.tolist() for k, v in rhs_expected.items()},
            "metric_source_residual_energy": metric_e_residual.tolist(),
            "metric_source_residual_momentum": metric_s_residual.tolist(),
        }

        species_reports[name] = {
            "metric_source_reconstruction": {
                "energy_residual": summarize(metric_e_residual, r),
                "momentum_residual": summarize(metric_s_residual, r),
                "max_abs_residual": max(
                    float(np.max(np.abs(metric_e_residual))),
                    float(np.max(np.abs(metric_s_residual))),
                ),
                "passed": bool(metric_closed),
                "interpretation": "Independent formula implementation of the declared Eulerian Valencia metric-source contraction; not a full continuum stress-divergence test.",
            },
            "production_rhs_component_accounting": {
                "channel_residuals": {
                    k: summarize(v, r) for k, v in rhs_residuals.items()
                },
                "max_abs_residual": rhs_max,
                "passed": bool(controls["production_rhs_matches_signed_component_assembly"]),
                "interpretation": "Checks direct RHS wiring against the helper's signed component arrays; both paths share the production component helper, so this is not independent physical conservation evidence.",
            },
            "finite_volume_signed_accounting": fv_balance,
        }

    # Independently compare the live baseline exchange factors with the locked,
    # conditional covariant-to-coordinate-time projections.
    dm_prim = primitives(metric_objects, initial_states["dark_matter"][1], Species.DARK_MATTER)
    rho_dm = finite("DM primitive density", [p.rho for p in dm_prim])
    W_dm = finite("DM Lorentz factor", [p.lorentz() for p in dm_prim])
    v_dm_actual = finite("DM recovered velocity", [p.v_r for p in dm_prim])
    sqrt_gamma = finite("sqrt gamma", [m.sqrt_gamma for m in metric_objects])
    Q_t = finite("Q_t", BETA_DM * rho_dm * dphi_t)
    Q_r = finite("Q_r", BETA_DM * rho_dm * dphi_r)
    dm_report = species_reports["dark_matter"]
    q_e_actual = np.asarray(signed_profiles["dark_matter"]["exchange_source_energy"])
    q_s_actual = np.asarray(signed_profiles["dark_matter"]["exchange_source_momentum"])
    q_e_expected = -sqrt_gamma * (Q_t - metric.beta * Q_r)
    q_s_expected = metric.alpha * sqrt_gamma * Q_r
    q_s_residual = q_s_actual - q_s_expected
    q_s_predicted_defect = (1.0 - metric.alpha) * sqrt_gamma * Q_r
    q_s_prediction_residual = q_s_residual - q_s_predicted_defect
    q_e_residual = q_e_actual - q_e_expected

    pi_normal = (dphi_t - metric.beta * dphi_r) / metric.alpha
    q_rest_expected = (
        -metric.alpha * sqrt_gamma * BETA_DM * rho_dm * W_dm
        * (pi_normal + v_dm_actual * dphi_r)
    )
    scalar_coupled = scalar_rhs_arrays(
        grid, geometry, fields, beta_dm=BETA_DM, rho_dm=rho_dm
    )
    scalar_delta_pi = finite(
        "scalar Pi coupling increment",
        scalar_coupled.Pi - scalar_uncoupled.Pi,
    )
    scalar_delta_pi_expected = metric.alpha * BETA_DM * rho_dm
    scalar_pi_residual = scalar_delta_pi - scalar_delta_pi_expected

    q_s_rel = rms(q_s_residual) / max(rms(q_s_expected), np.finfo(float).tiny)
    q_d_rel = rms(q_rest_expected) / max(rms(q_rest_expected), np.finfo(float).tiny)
    pi_rel = rms(scalar_pi_residual) / max(rms(scalar_delta_pi_expected), np.finfo(float).tiny)
    exchange_controls = {
        "energy_qE_matches_covector_projection": bool(np.max(np.abs(q_e_residual)) <= ABS_TOL),
        "radial_momentum_mismatch_matches_predicted_lapse_factor": bool(
            np.allclose(q_s_residual, q_s_predicted_defect, rtol=REL_TOL, atol=ABS_TOL)
        ),
        "radial_momentum_mismatch_detected_for_nonunit_lapse": bool(
            rms(q_s_expected) > 1.0e-12 and q_s_rel > 1.0e-2
        ),
        "variable_mass_rest_current_source_absent_in_live_rhs": bool(
            float(np.max(np.abs(q_rest_expected))) > 1.0e-12
        ),
        "scalar_pi_lapse_mismatch_detected": bool(
            rms(scalar_delta_pi) > 1.0e-12 and pi_rel > 1.0e-2
        ),
        "scalar_unit_lapse_controls_deferred": True,
    }
    controls["all_profiles_finite"] = True
    for arr in (
        q_e_residual, q_s_residual, q_s_prediction_residual,
        q_rest_expected, scalar_delta_pi, scalar_pi_residual
    ):
        controls["all_profiles_finite"] &= bool(np.all(np.isfinite(arr)))

    signed_profiles["conditional_exchange_witness"] = {
        "Q_t": Q_t.tolist(),
        "Q_r": Q_r.tolist(),
        "qE_actual_minus_expected": q_e_residual.tolist(),
        "qSr_actual_minus_expected": q_s_residual.tolist(),
        "qSr_residual_minus_predicted_lapse_defect": q_s_prediction_residual.tolist(),
        "expected_qD_source_missing_from_live_RHS": q_rest_expected.tolist(),
        "scalar_delta_Pi_actual": scalar_delta_pi.tolist(),
        "scalar_delta_Pi_expected_alpha_weighted": scalar_delta_pi_expected.tolist(),
        "scalar_delta_Pi_residual": scalar_pi_residual.tolist(),
    }
    exchange_report = {
        "conditional_law": "Q_nu = beta_DM*rho_DM*partial_nu(phi), scalar receives the opposite exchange under the declared model",
        "beta_dm": BETA_DM,
        "qE_projection": {
            "expected_formula": "-sqrt(gamma)*(Q_t-beta^r*Q_r)",
            "residual": summarize(q_e_residual, r),
            "passed": exchange_controls["energy_qE_matches_covector_projection"],
        },
        "qSr_projection": {
            "expected_formula": "alpha*sqrt(gamma)*Q_r",
            "actual_legacy_formula": "sqrt(gamma)*Q_r",
            "relative_rms_residual": q_s_rel,
            "residual": summarize(q_s_residual, r),
            "predicted_lapse_defect_closes": exchange_controls["radial_momentum_mismatch_matches_predicted_lapse_factor"],
            "mismatch_detected": exchange_controls["radial_momentum_mismatch_detected_for_nonunit_lapse"],
        },
        "qD_rest_current": {
            "expected_formula_under_variable_mass_dust_premise": "-alpha*sqrt(gamma)*beta*rho*W*(Pi+v^r*D_r(phi))",
            "max_abs_expected_source": float(np.max(np.abs(q_rest_expected))),
            "baseline_rhs_contains_this_source": False,
            "interpretation": "Conditional diagnostic discriminator, not an accepted production-law change.",
        },
        "scalar_Pi_exchange": {
            "actual_increment_formula": "beta_DM*rho_DM",
            "expected_increment_under_covariant_coordinate_time_source": "alpha*beta_DM*rho_DM",
            "relative_rms_residual": pi_rel,
            "residual": summarize(scalar_pi_residual, r),
            "mismatch_detected": exchange_controls["scalar_pi_lapse_mismatch_detected"],
        },
        "controls": exchange_controls,
    }

    input_digest = hashlib.sha256()
    for name, values in (
        ("r", r), ("a", metric.a), ("b", metric.b), ("X", metric.X),
        ("alpha", metric.alpha), ("beta", metric.beta), ("Aa", metric.Aa),
        ("K", metric.K), ("phi", fields.phi), ("Pi", fields.Pi),
        ("dphi_t", dphi_t), ("dphi_r", dphi_r),
    ):
        input_digest.update(name.encode("utf-8"))
        input_digest.update(np.ascontiguousarray(values, dtype=np.float64).tobytes())

    result = {
        "schema": "issue4_total_valencia_stage_budget_v1",
        "status": "completed_diagnostic_only",
        "source_provenance": source_provenance,
        "input_state_sha256": input_digest.hexdigest(),
        "configuration": {
            "resolution": int(r.size), "r_max": float(args.r_max),
            "dr": dr, "beta_dm": BETA_DM,
            "alpha_min": float(np.min(metric.alpha)),
            "alpha_max": float(np.max(metric.alpha)),
            "nonzero_shift": bool(np.max(np.abs(metric.beta)) > 0.0),
            "moving_dark_matter": bool(np.max(np.abs(v_dm_actual)) > 0.0),
            "nonzero_scalar_gradient": bool(np.max(np.abs(dphi_r)) > 0.0),
            "production_steps_advanced": 0,
        },
        "species": species_reports,
        "conditional_exchange_source_audit": exchange_report,
        "controls": controls,
        "maxima": {
            "independent_metric_source_abs_residual": max_metric_source_abs,
            "direct_rhs_component_accounting_abs_residual": max_component_rhs_abs,
            "finite_volume_telescope_abs_residual": max_telescoping_abs,
        },
        "gate_result": {
            "source_provenance_matches": True,
            "independent_metric_source_formula_matches": bool(controls["independent_metric_sources_match"]),
            "production_rhs_component_accounting_closes": bool(controls["production_rhs_matches_signed_component_assembly"]),
            "finite_volume_accounting_closes": bool(controls["finite_volume_telescopes_close"]),
            "baseline_covariant_exchange_mismatch_reproduced": bool(
                exchange_controls["energy_qE_matches_covector_projection"]
                and exchange_controls["radial_momentum_mismatch_matches_predicted_lapse_factor"]
                and exchange_controls["radial_momentum_mismatch_detected_for_nonunit_lapse"]
                and exchange_controls["variable_mass_rest_current_source_absent_in_live_rhs"]
                and exchange_controls["scalar_pi_lapse_mismatch_detected"]
            ),
            "independent_total_stress_energy_divergence_closed": False,
            "brown_bianchi_propagation_gate_closed": False,
            "archive_SD_physical_normalization_closed": False,
            "scientific_admission": False,
        },
        "interpretation": (
            "Finite-volume telescoping and direct RHS/component agreement are "
            "accounting/wiring checks, not independent div(T_total)=0 evidence. "
            "Metric-source reconstruction is an independent formula implementation "
            "for the declared Eulerian Valencia convention. Exchange mismatches "
            "are conditional on the stated variable-mass covariant law. No production "
            "step or equation change was performed."
        ),
        "signed_profiles": signed_profiles,
    }
    result["gate_result"]["all_bounded_diagnostic_checks_pass"] = bool(
        result["gate_result"]["source_provenance_matches"]
        and result["gate_result"]["independent_metric_source_formula_matches"]
        and result["gate_result"]["production_rhs_component_accounting_closes"]
        and result["gate_result"]["finite_volume_accounting_closes"]
        and result["gate_result"]["baseline_covariant_exchange_mismatch_reproduced"]
        and controls["all_profiles_finite"]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": result["status"],
        "resolution": n,
        "source_pins_match": True,
        "independent_metric_sources_match": result["gate_result"]["independent_metric_source_formula_matches"],
        "finite_volume_accounting_closes": result["gate_result"]["finite_volume_accounting_closes"],
        "baseline_covariant_exchange_mismatch_reproduced": result["gate_result"]["baseline_covariant_exchange_mismatch_reproduced"],
        "independent_total_stress_energy_divergence_closed": False,
        "scientific_admission": False,
        "artifact": str(args.output),
    }, sort_keys=True))
    return 0 if result["gate_result"]["all_bounded_diagnostic_checks_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
