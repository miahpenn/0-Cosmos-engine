#!/usr/bin/env python3
"""Evaluate existing Issue #4 propagation candidates on saved accepted states.

Research diagnostic only. Loads SHA-verified t=0.60 profiles from the prior
production-baseline/joint-origin run and the full-continuous-metric shadow.
Runs the existing RHS tangent and provisional Brown/source-divergence candidate without redefining its equation convention
on those exact states. Lapse perturbations use the production CMC solver for
the baseline and the copied regular-origin CMC row for the two shadow arms.

This does not assert the candidate formula is the accepted propagation oracle:
the matter source mapping/index normalization and complete sector conservation
identity remain provisional. No evolution is repeated and no production source
is changed.
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
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
# The existing sector-resolved divergence module requires a source-root variable
# at import time; actual archived source/vendor provenance is verified separately
# by verify_source_pins() and the input artifact manifests below.
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import V55MatterState
from engine import v55_pirk_adapter as adapter
from tools import audit_current_source_momentum_constraint_tangent as tangent
from tools import audit_hybrid_origin_integral_patch_newton as propagation
from tools import audit_total_stress_energy_divergence as stress_divergence
from tools import audit_issue4_production_hll_operator_pairing as hll_pairing
from tools.audit_origin_cell_joint_shadow_one_step import (
    SOURCE_COMMIT,
    VENDOR_COMMIT,
    regular_origin_cmc_lapse,
    verify_source_pins,
)

RESOLUTION = int(os.environ.get("PROPAGATION_RESOLUTION", "120"))
JOINT_DIR = Path(os.environ.get("JOINT_ARTIFACT_DIR", f"reference/joint/N{RESOLUTION}"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{RESOLUTION}"))
OUT = Path(os.environ.get("PROPAGATION_OUT", f"results/full-profile-propagation-N{RESOLUTION}-2026-10-25"))
EXPECTED_SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
EXPECTED_VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_HLL_OBSERVER_BLOB = "319667f115a13d8941af7bd3142d97edce76ca4f"
TARGET_TIME = 0.60
ARMS = (
    ("production_baseline", "production_baseline", False),
    ("partial_joint", "joint_integral_and_cmc_row", True),
    ("full_continuous_profile", "full_continuous_metric_profile", True),
)
SCALAR_NAMES = ("S", "PS", "D", "PD", "phi", "Pi")
GEOMETRY_NAMES = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
MATTER_FIELDS = ("rest", "energy_t", "momentum_r")
EPSILONS = (2.0e-5, 1.0e-5)
FULL_METRIC_LAPSE_SHADOW = False
SHIFTED_STATES: dict[float, Any] = {}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_artifact(path: Path, resolution: int, *, full_profile: bool) -> dict[str, Any]:
    if not path.is_dir():
        raise FileNotFoundError(f"artifact directory does not exist: {path}")
    manifest_path = path / "SHA256SUMS.txt"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"artifact checksum manifest missing: {manifest_path}")
    manifest = {}
    for line in manifest_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            raise RuntimeError(f"malformed checksum row in {manifest_path}: {line!r}")
        expected, name = parts
        name = name.strip().lstrip("*")
        candidate = (path / name).resolve()
        if path.resolve() not in candidate.parents:
            raise RuntimeError(f"manifest path escapes artifact directory: {name}")
        if not candidate.is_file() or sha256(candidate) != expected:
            raise RuntimeError(f"artifact hash mismatch: {candidate}")
        manifest[name] = expected
    for required in ("report.json", "profiles.npz"):
        if required not in manifest:
            raise RuntimeError(f"artifact manifest omits {required}: {path}")
    report = json.loads((path / "report.json").read_text(encoding="utf-8"))
    cfg = report.get("configuration", {})
    if int(cfg.get("resolution", -1)) != resolution:
        raise RuntimeError(f"{path}: resolution mismatch in report")
    if abs(float(cfg.get("target_time", -1.0)) - TARGET_TIME) > 1e-12:
        raise RuntimeError(f"{path}: target-time mismatch in report")
    if float(cfg.get("r_max", -1.0)) != 80.0:
        raise RuntimeError(f"{path}: domain mismatch in report")
    pins = report.get("source_pins", {}).get("source_pin_check", {})
    if pins.get("source_commit") != EXPECTED_SOURCE or pins.get("vendor_commit") != EXPECTED_VENDOR:
        raise RuntimeError(f"{path}: source/vendor provenance mismatch")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError(f"{path}: protected source diff was not empty")
    if report.get("scientific_admission") is not False:
        raise RuntimeError(f"{path}: unexpectedly changed scientific_admission")
    if full_profile and cfg.get("arms") != ["full_continuous_metric_profile"]:
        raise RuntimeError(f"{path}: not a single full-profile artifact")
    return {
        "directory": str(path),
        "manifest_sha256": sha256(manifest_path),
        "manifest": manifest,
        "report_sha256": sha256(path / "report.json"),
        "profiles_sha256": sha256(path / "profiles.npz"),
        "status": report.get("status"),
        "source_commit": pins["source_commit"],
        "vendor_commit": pins["vendor_commit"],
        "configuration": cfg,
    }


def load_profile_arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path / "profiles.npz", allow_pickle=False) as data:
        return {key: np.asarray(data[key]).copy() for key in data.files}


def scalar(data: dict[str, np.ndarray], key: str) -> float:
    return float(np.asarray(data[key]).reshape(()))


def reconstruct_state(resolution: int, data: dict[str, np.ndarray], prefix: str):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution,
        r_max=80.0,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    r_saved = np.asarray(data[prefix + "r"], dtype=np.float64)
    r_live = np.asarray(state.grid.centers, dtype=np.float64)
    if r_saved.shape != r_live.shape or not np.array_equal(r_saved, r_live):
        difference = float(np.max(np.abs(r_saved - r_live))) if r_saved.shape == r_live.shape else None
        raise RuntimeError(f"{prefix}: stored grid does not exactly match the pinned initializer grid; max_delta={difference}")

    geom = state.geometry.copy()
    for name in GEOMETRY_NAMES:
        value = np.asarray(data[prefix + name], dtype=np.float64)
        if value.shape != r_live.shape or not np.all(np.isfinite(value)):
            raise RuntimeError(f"{prefix}: invalid geometry field {name}")
        setattr(geom, name, value.copy())
    state.geometry = geom

    scalar_arrays = []
    for name in SCALAR_NAMES:
        value = np.asarray(data[prefix + name], dtype=np.float64)
        if value.shape != r_live.shape or not np.all(np.isfinite(value)):
            raise RuntimeError(f"{prefix}: invalid scalar field {name}")
        scalar_arrays.append(value.copy())
    state.scalars = ScalarFields(*scalar_arrays)

    def species(short: str) -> ConservedSpecies:
        values = []
        for name in MATTER_FIELDS:
            value = np.asarray(data[prefix + short + "_" + name], dtype=np.float64)
            if value.shape != r_live.shape or not np.all(np.isfinite(value)):
                raise RuntimeError(f"{prefix}: invalid matter field {short}.{name}")
            values.append(value.copy())
        return ConservedSpecies(*values)

    state.matter = V55MatterState(
        species("dm"), species("baryon"), species("radiation")
    )
    state.t = scalar(data, prefix + "t")
    state.tau = scalar(data, prefix + "tau")
    state.e_folds = scalar(data, prefix + "e_folds")
    if abs(state.t - TARGET_TIME) > 1e-12:
        raise RuntimeError(f"{prefix}: saved state is not accepted at t={TARGET_TIME}")
    if getattr(state, "electromagnetic", None) is not None:
        raise RuntimeError("this six-sector ledger must fail closed on an electromagnetic state")
    if getattr(state.matter, "charged_dust", None) is not None:
        raise RuntimeError("this six-sector ledger must fail closed on charged dust")
    return kernel, state


def shifted_state_with_matched_lapse(kernel, state, rates, delta):
    """Reuse the pinned tangent displacement, but match each arm's CMC row."""
    geom = state.geometry.copy()
    for name in tangent.GEOMETRY_RATE_NAMES:
        updated = (
            np.asarray(getattr(state.geometry, name), dtype=np.float64)
            + delta * np.asarray(rates["geometry"][name], dtype=np.float64)
        )
        setattr(geom, name, updated.copy())

    scalars = ScalarFields(*(
        np.asarray(getattr(state.scalars, name), dtype=np.float64)
        + delta * np.asarray(rates["scalars"][name], dtype=np.float64)
        for name in SCALAR_NAMES
    ))
    matter_parts = {}
    for species_name, short in (
        ("dark_matter", "dm"), ("baryons", "baryon"), ("radiation", "radiation")
    ):
        base = getattr(state.matter, species_name)
        rate = rates["matter"][species_name]
        matter_parts[species_name] = ConservedSpecies(*(
            np.asarray(getattr(base, name), dtype=np.float64)
            + delta * np.asarray(rate[name], dtype=np.float64)
            for name in MATTER_FIELDS
        ))
    matter = V55MatterState(
        matter_parts["dark_matter"], matter_parts["baryons"], matter_parts["radiation"]
    )
    geom.beta = np.zeros_like(geom.beta)
    geom.B = np.zeros_like(geom.B)
    if FULL_METRIC_LAPSE_SHADOW:
        alpha, _ = regular_origin_cmc_lapse(state.grid, geom, scalars, matter)
    else:
        alpha = kernel._solve_lapse(state.grid, geom, scalars, matter)[0]
    geom.alpha = np.asarray(alpha, dtype=np.float64).copy()
    geom.beta = np.zeros_like(geom.beta)
    geom.B = np.zeros_like(geom.B)
    shifted = type(state)(
        grid=state.grid,
        geometry=geom,
        scalars=scalars,
        matter=matter,
        t=state.t,
        tau=state.tau,
        e_folds=state.e_folds,
    )
    SHIFTED_STATES[float(delta)] = shifted
    return shifted


def summarize_profile(values: Any, radii: np.ndarray) -> dict[str, Any]:
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape != radii.shape or not np.all(np.isfinite(arr)):
        raise RuntimeError("invalid profile supplied to residual summary")
    centre = radii <= 3.0
    annulus = (radii >= 2.0) & (radii < 10.0)
    outer = radii >= 5.0
    return {
        "cell_0_to_2": [float(x) for x in arr[:3]],
        "cell_0_to_4": [float(x) for x in arr[:5]],
        "max_abs_r_le_3": float(np.max(np.abs(arr[centre]))) if np.any(centre) else None,
        "rms_r_2_to_10": float(np.sqrt(np.mean(arr[annulus] ** 2))) if np.any(annulus) else None,
        "max_abs_r_ge_5": float(np.max(np.abs(arr[outer]))) if np.any(outer) else None,
        "rms_full_domain": float(np.sqrt(np.mean(arr * arr))),
        "fixed_radius": {
            str(q): float(np.interp(q, radii, arr))
            for q in (1.0, 2.0, 3.0, 5.0, 8.0, 10.0)
            if radii[0] <= q <= radii[-1]
        },
    }


def compact_candidate(item: dict[str, Any], raw: dict[float, dict[str, np.ndarray]], radii: np.ndarray) -> dict[str, Any]:
    eps = float(item["epsilon"])
    prop = item["independent_constraint_propagation_candidate"]
    direct = prop["published_covariant_M_r_candidate"]
    geometric = prop["geometric_constraint_propagation_candidate"]
    direct_residual = direct["measured_minus_predicted_covariant"]
    source_map = geometric["matter_source_composite_minus_stress_divergence"]
    brown_residual = geometric["measured_total_minus_Brown_plus_actual_source"]
    raw_item = raw[eps]

    return {
        "epsilon": eps,
        "direct_published_covariant_Mr_candidate": {
            "relative_rms_error_full_domain": direct["relative_rms_error_full_domain"],
            "relative_rms_error_r_ge_5": direct["relative_rms_error_r_ge_5"],
            "measured_minus_predicted": direct_residual,
        },
        "geometric_Brown_plus_actual_source_candidate": {
            "relative_rms": geometric["full_domain_relative_rms"]["total_Brown_plus_actual_matter_composite"],
            "measured_minus_predicted": brown_residual,
            "matter_source_composite_minus_minus8pi_alpha_Qr": source_map,
            "source_mapping_relative_rms": geometric["full_domain_relative_rms"]["actual_matter_composite_minus_minus8pi_alpha_Qr"],
            "nonmatter_geometric_minus_Brown_vacuum_relative_rms": geometric["full_domain_relative_rms"]["nonmatter_geometric_minus_Brown_vacuum"],
        },
        "explicit_propagation_profiles": {
            name: summarize_profile(raw_item[name], radii)
            for name in (
                "propagation_Mtotal_measured_minus_Brown_plus_actual_matter_composite",
                "propagation_matter_source_mapping_residual",
                "propagation_exact_chain_rule_residual",
                "propagation_Mtotal_decomposition_residual",
                "propagation_term_H_gradient",
                "propagation_term_connection_second",
                "propagation_term_connection_first",
                "propagation_published_Mr_term_KM",
                "propagation_published_Mr_term_alpha_H",
                "propagation_published_Mr_matter_source",
                "propagation_Qr_total_stress_divergence",
                "propagation_measured_Mgeom_dot",
                "propagation_Mgeom_dot_nonmatter_measured",
                "propagation_Mgeom_dot_matter_injection",
                "propagation_geometric_Brown_term_shift",
                "propagation_geometric_Brown_term_KM",
                "propagation_geometric_Brown_term_alpha_H",
                "propagation_geometric_Brown_term_H_gradient",
                "propagation_geometric_Brown_term_connection_second",
                "propagation_geometric_Brown_term_connection_first",
                "propagation_geometric_Brown_term_sum_reconstruction_gap",
                "propagation_geometric_Brown_vacuum_source",
            )
            if name in raw_item
        },
    }


def compact_sector_divergence(
    state,
    epsilon: float,
    candidate_raw: dict[str, np.ndarray],
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Reuse the existing additive sector-resolved total-T divergence observer."""
    try:
        plus_state = SHIFTED_STATES[float(epsilon)]
        minus_state = SHIFTED_STATES[-float(epsilon)]
    except KeyError as exc:
        raise RuntimeError(f"matched tangent endpoint not captured for epsilon={epsilon}") from exc
    witnessed = stress_divergence.sector_resolved_divergence_tangent_from_states(
        state, plus_state, minus_state, float(epsilon)
    )
    total = witnessed["total"]
    decomposition = witnessed["sector_decomposition"]
    radii = np.asarray(state.grid.centers, dtype=np.float64)
    total_r = np.asarray(total["div_T_r"], dtype=np.float64)
    total_r_fv = np.asarray(total["div_T_r_fv"], dtype=np.float64)
    candidate_qr = np.asarray(
        candidate_raw["propagation_Qr_total_stress_divergence"], dtype=np.float64
    )
    q_compare = total_r - candidate_qr
    recon = decomposition["reconstruction"]
    channels = recon["channels"]
    sector_rows = {}
    profile_store: dict[str, np.ndarray] = {
        "sector_resolved_total_div_T_r": total_r,
        "sector_resolved_total_div_T_r_fv": total_r_fv,
        "sector_resolved_total_minus_candidate_Qr": q_compare,
        "sector_resolved_sector_sum_minus_total_div_T_r": np.asarray(
            channels["div_T_r"]["sum_minus_total"], dtype=np.float64
        ),
        "sector_resolved_unexplained_gap_after_time_density_div_T_r": np.asarray(
            channels["div_T_r"]["unexplained_gap_after_time_density_accounting"],
            dtype=np.float64,
        ),
    }
    for name, data in decomposition["sectors"].items():
        raw = data["raw_signed_arrays"]
        sector_r = np.asarray(raw["div_T_r"], dtype=np.float64)
        sector_r_fv = np.asarray(raw["div_T_r_fv"], dtype=np.float64)
        profile_store[f"sector_{name}_div_T_r"] = sector_r
        profile_store[f"sector_{name}_div_T_r_fv"] = sector_r_fv
        sector_rows[name] = {
            "div_T_r": summarize_profile(sector_r, radii),
            "div_T_r_finite_volume": summarize_profile(sector_r_fv, radii),
            "cell0_cell1": [float(sector_r[0]), float(sector_r[1])],
        }
        # For scalar sectors, also retain the existing observer's signed split
        # into time-density, radial-flux, and metric-source components.
        operator = data.get("operator_component_decomposition")
        if operator:
            pointwise = operator["pointwise"]["raw_signed_arrays"]
            names = (
                "time_density_div_T_r",
                "radial_flux_div_T_r",
                "metric_source_div_T_r",
                "component_sum_div_T_r",
                "component_sum_gap_div_T_r",
            )
            sector_rows[name]["pointwise_operator_terms"] = {
                term: summarize_profile(np.asarray(pointwise[term], dtype=np.float64), radii)
                for term in names if term in pointwise
            }
            for term in names:
                if term in pointwise:
                    profile_store[f"sector_{name}_{term}"] = np.asarray(
                        pointwise[term], dtype=np.float64
                    )
    compact = {
        "epsilon": float(epsilon),
        "observer": "existing audit_total_stress_energy_divergence.sector_resolved_divergence_tangent_from_states",
        "scope": "same saved base accepted state and the exact matched-lapse plus/minus endpoint states used by the existing propagation candidate",
        "total_div_T_r": summarize_profile(total_r, radii),
        "total_div_T_r_finite_volume": summarize_profile(total_r_fv, radii),
        "independent_total_divergence_minus_candidate_Qr": summarize_profile(q_compare, radii),
        "sector_contributions": sector_rows,
        "additive_sector_reconstruction": {
            "all_sectors_reconstruct_total": bool(recon["all_sectors_reconstruct_total"]),
            "total_div_T_r_sum_passed": bool(channels["div_T_r"]["passed"]),
            "total_div_T_r_sum_max_abs_gap": float(channels["div_T_r"]["max_abs_gap"]),
            "total_div_T_r_component_corrected_passed": bool(
                channels["div_T_r"]["component_corrected_passed"]
            ),
            "total_div_T_r_unexplained_max_abs_gap": float(
                channels["div_T_r"]["max_abs_unexplained_gap"]
            ),
            "time_density_component_accounting_passed": bool(
                recon["time_density_component_accounting_passed"]
            ),
            "time_density_gap_div_T_r": summarize_profile(
                np.asarray(
                    decomposition["time_density_reconstruction"]["raw_signed_arrays"]["time_density_gap_r"],
                    dtype=np.float64,
                ),
                radii,
            ),
        },
        "sector_pairing": decomposition["state_pairing"],
        "scientific_admission": False,
    }
    return compact, profile_store


def operator_pairing_decomposition(
    kernel, state, rates, epsilon: float, candidate_raw: dict[str, np.ndarray]
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Decompose pointwise-minus-expanded Q_r into discrete product/chain-rule gaps.

    Reuses the existing unexpanded covariant-divergence formula, and the same
    saved-state directional endpoints as the expanded candidate. It does not
    select which operator is the physical/discrete acceptance convention.
    """
    eps = float(epsilon)
    try:
        plus_state = SHIFTED_STATES[eps]
        minus_state = SHIFTED_STATES[-eps]
    except KeyError as exc:
        raise RuntimeError(f"operator-pairing endpoints missing for epsilon={eps}") from exc

    grid = state.grid
    geom = state.geometry
    r = np.asarray(grid.centers, dtype=np.float64)
    d1_fn = lambda values, parity: stress_divergence.d1(grid, values, parity)
    base_projection = stress_divergence.independent_total_projections(state)
    plus_projection = stress_divergence.independent_total_projections(plus_state)
    minus_projection = stress_divergence.independent_total_projections(minus_state)

    base_obj = stress_divergence._covariant_objects_from_components(state, base_projection)
    plus_obj = stress_divergence._covariant_objects_from_components(plus_state, plus_projection)
    minus_obj = stress_divergence._covariant_objects_from_components(minus_state, minus_projection)

    dt_density_t = (
        np.asarray(plus_obj["density_t"], dtype=np.float64)
        - np.asarray(minus_obj["density_t"], dtype=np.float64)
    ) / (2.0 * eps)
    dt_density_r = (
        np.asarray(plus_obj["density_r"], dtype=np.float64)
        - np.asarray(minus_obj["density_r"], dtype=np.float64)
    ) / (2.0 * eps)
    dt_gtt = (
        np.asarray(plus_obj["gtt"], dtype=np.float64)
        - np.asarray(minus_obj["gtt"], dtype=np.float64)
    ) / (2.0 * eps)
    dt_grr = (
        np.asarray(plus_obj["gamma_rr"], dtype=np.float64)
        - np.asarray(minus_obj["gamma_rr"], dtype=np.float64)
    ) / (2.0 * eps)
    dt_gth = (
        np.asarray(plus_obj["gamma_thth"], dtype=np.float64)
        - np.asarray(minus_obj["gamma_thth"], dtype=np.float64)
    ) / (2.0 * eps)

    pointwise_terms = stress_divergence.covariant_divergence_terms_from_profiles(
        r, d1_fn,
        base_obj["alpha"], base_obj["gamma_rr"], base_obj["gamma_thth"],
        base_obj["rho"], base_obj["j"], base_obj["pr"], base_obj["pt"],
        dt_density_t, dt_density_r, dt_gtt, dt_grr, dt_gth,
    )
    terms = pointwise_terms["terms"]
    pointwise = np.asarray(pointwise_terms["div_T_r"], dtype=np.float64)
    expanded = np.asarray(candidate_raw["propagation_Qr_total_stress_divergence"], dtype=np.float64)

    alpha = np.asarray(geom.alpha, dtype=np.float64)
    a = np.asarray(geom.a, dtype=np.float64)
    b = np.asarray(geom.b, dtype=np.float64)
    X = np.asarray(geom.X, dtype=np.float64)
    gamma_rr = a / X**2
    gamma_thth = b * r**2 / X**2
    rho = np.asarray(base_projection["rho"], dtype=np.float64)
    j = np.asarray(base_projection["j"], dtype=np.float64)
    pr = np.asarray(base_projection["pr"], dtype=np.float64)
    pt = np.asarray(base_projection["pt"], dtype=np.float64)

    alpha_r = d1_fn(alpha, +1)
    gamma_rr_r = d1_fn(gamma_rr, +1)
    gamma_thth_r = d1_fn(gamma_thth, +1)
    b_r = d1_fn(b, +1)
    X_r = d1_fn(X, +1)
    areal_log_r = 1.0 / r + 0.5 * b_r / b - X_r / X
    rho_alpha_log_r = rho * alpha_r / alpha

    total_base = assemble_total_stress_energy(grid, state.geometry, state.scalars, state.matter)
    total_plus = assemble_total_stress_energy(grid, plus_state.geometry, plus_state.scalars, plus_state.matter)
    total_minus = assemble_total_stress_energy(grid, minus_state.geometry, minus_state.scalars, minus_state.matter)
    j_dot = (
        np.asarray(total_plus.j, dtype=np.float64)
        - np.asarray(total_minus.j, dtype=np.float64)
    ) / (2.0 * eps)
    a_dot = np.asarray(rates["geometry"]["a"], dtype=np.float64)
    b_dot = np.asarray(rates["geometry"]["b"], dtype=np.float64)
    X_dot = np.asarray(rates["geometry"]["X"], dtype=np.float64)
    gamma_rr_log_dot = a_dot / a - 2.0 * X_dot / X
    areal_log_dot = 0.5 * b_dot / b - X_dot / X
    expanded_time = j_dot / alpha + (j / alpha) * (
        0.5 * gamma_rr_log_dot + 2.0 * areal_log_dot
    )
    pointwise_time = np.asarray(terms["time_density_div_T_r"], dtype=np.float64)
    time_product_rule_gap = pointwise_time - expanded_time

    flux_expanded = (
        d1_fn(pr, +1)
        + pr * (
            alpha_r / alpha
            + 0.5 * gamma_rr_r / gamma_rr
            + gamma_thth_r / gamma_thth
        )
    )
    flux_product_rule_gap = (
        np.asarray(terms["radial_flux_div_T_r"], dtype=np.float64) - flux_expanded
    )
    gtt_metric_chain_gap = (
        np.asarray(terms["metric_source_div_T_r_gtt"], dtype=np.float64)
        - rho_alpha_log_r
    )
    areal_metric_chain_gap = (
        (pr - pt) * (
            gamma_thth_r / gamma_thth - 2.0 * areal_log_r
        )
    )

    total_gap = pointwise - expanded
    reconstructed_gap = (
        time_product_rule_gap + flux_product_rule_gap
        + gtt_metric_chain_gap + areal_metric_chain_gap
    )
    closure = total_gap - reconstructed_gap

    def stat(values):
        arr = np.asarray(values, dtype=np.float64)
        centre = r <= 3.0
        annulus = (r >= 2.0) & (r < 10.0)
        return {
            "cell0_cell1": [float(arr[0]), float(arr[1])],
            "max_abs_r_le_3": float(np.max(np.abs(arr[centre]))),
            "rms_r_2_to_10": float(np.sqrt(np.mean(arr[annulus]**2))),
            "rms_full_domain": float(np.sqrt(np.mean(arr**2))),
            "at_fixed_r": {
                str(q): float(np.interp(q, r, arr))
                for q in (1.0, 2.0, 3.0, 5.0, 8.0) if r[0] <= q <= r[-1]
            },
        }

    profile_map = {
        "pointwise_div_T_r": pointwise,
        "expanded_Q_r": expanded,
        "pointwise_minus_expanded": total_gap,
        "time_product_rule_gap": time_product_rule_gap,
        "radial_flux_product_rule_gap": flux_product_rule_gap,
        "gtt_lapse_metric_chain_gap": gtt_metric_chain_gap,
        "areal_metric_chain_gap": areal_metric_chain_gap,
        "reconstructed_sum_of_discrete_gaps": reconstructed_gap,
        "decomposition_closure": closure,
    }
    report = {
        "epsilon": eps,
        "observer": "existing unexpanded covariant divergence versus the existing expanded Q_r candidate",
        "pointwise_minus_expanded": stat(total_gap),
        "decomposition_terms": {name: stat(values) for name, values in profile_map.items() if name not in ("pointwise_div_T_r", "expanded_Q_r")},
        "decomposition_closure_max_abs": float(np.max(np.abs(closure))),
        "decomposition_closure_rms": float(np.sqrt(np.mean(closure**2))),
        "interpretation": (
            "Exact algebraic partition of the discrete mismatch into time-product, "
            "radial-flux product, g_tt/lapse-chain, and areal-radius chain gaps. "
            "The partition is diagnostic only; it does not select an authoritative "
            "discrete divergence or physical source law."
        ),
    }
    return report, profile_map


def compact_hll_operator_snapshot(state, rates) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Pair current source-pinned production HLL face/source terms on this saved state."""
    actual_blob = subprocess.check_output(
        ["git", "rev-parse", f"HEAD:tools/audit_issue4_production_hll_operator_pairing.py"],
        cwd=ROOT, text=True
    ).strip()
    if actual_blob != EXPECTED_HLL_OBSERVER_BLOB:
        raise RuntimeError(
            "HLL diagnostic observer blob changed: "
            f"expected={EXPECTED_HLL_OBSERVER_BLOB}, actual={actual_blob}"
        )
    coupling = float(hll_pairing.BETA_DM_COUPLED)
    production_coupling = float(hll_pairing.production_kernel_module.BETA_DM)
    if production_coupling != coupling:
        raise RuntimeError(
            "production RHS coupling and HLL source observer coupling differ: "
            f"{production_coupling} != {coupling}"
        )
    before = tangent.state_digest(state)
    snapshot = hll_pairing.operator_snapshot(state, coupling, rates)
    after = tangent.state_digest(state)
    if before != after:
        raise RuntimeError("HLL operator observer mutated its accepted base state")
    if not snapshot.get("origin_face_flux_exact_zero", False):
        raise RuntimeError("production HLL observer failed exact-zero spherical centre face")
    if not snapshot.get("native_component_arrays_preserved", False):
        raise RuntimeError("production HLL observer omitted native component arrays")

    species_summary = {}
    profiles: dict[str, np.ndarray] = {}
    for name, item in snapshot["species"].items():
        raw = snapshot["raw_arrays"][name]
        species_summary[name] = {
            "origin_face_flux_exact_zero": bool(item["origin_face_flux_exact_zero"]),
            "hll_minus_projected_stress_momentum_face_flux": item[
                "hll_momentum_flux_minus_projected_stress_face_flux"
            ],
            "geometric_momentum_source_gap_over_kappa": item[
                "geometric_momentum_source_gap_over_kappa"
            ],
            "interaction_momentum_source_gap_over_kappa": item[
                "interaction_momentum_source_gap_over_kappa"
            ],
            "production_rhs_momentum_over_kappa": item[
                "production_rhs_budget"
            ]["radial_momentum"],
            "production_rhs_negative_controls": item["production_rhs_budget"]["negative_controls"],
        }
        for raw_key in (
            "hll_face_flux_momentum_over_kappa",
            "projected_stress_momentum_face_flux_over_kappa",
            "hll_minus_projected_momentum_face_flux_over_kappa",
            "geometric_momentum_source_actual_over_kappa",
            "geometric_momentum_source_independent_over_kappa",
            "geometric_momentum_source_gap_over_kappa",
            "interaction_momentum_source_actual_over_kappa",
            "interaction_momentum_source_conditional_covariant_over_kappa",
            "interaction_momentum_source_gap_over_kappa",
            "production_rhs_momentum_over_kappa",
        ):
            profiles[f"hll_{name}_{raw_key}"] = np.asarray(raw[raw_key], dtype=np.float64).copy()

    scalar_summary = snapshot["scalar_sector_projections"]
    compact = {
        "observer": "existing audit_issue4_production_hll_operator_pairing.operator_snapshot",
        "observer_blob": actual_blob,
        "accepted_state_sha256_before": before,
        "accepted_state_sha256_after": after,
        "state_immutable": True,
        "accepted_time": float(state.t),
        "zero_shift": bool(snapshot["zero_shift"]),
        "origin_face_flux_exact_zero_for_all_species": bool(snapshot["origin_face_flux_exact_zero"]),
        "native_component_arrays_preserved": bool(snapshot["native_component_arrays_preserved"]),
        "coupling_parameter": coupling,
        "production_rhs_coupling_parameter": production_coupling,
        "max_abs_geometric_source_gap_over_kappa": float(
            snapshot["max_abs_geometric_source_gap_over_kappa"]
        ),
        "max_abs_dark_matter_lapse_omission_control_gap_over_kappa": float(
            snapshot["max_abs_dark_matter_lapse_omission_control_gap_over_kappa"]
        ),
        "species": species_summary,
        "scalar_sector_projection_summaries": scalar_summary,
        "scope_note": (
            "This records production fluid HLL face fluxes and source quadrature on the "
            "same saved accepted state. S/D/phi remain cell-projected scalar sectors, not "
            "fluid HLL fluxes. Coupling-source interpretation remains conditional."
        ),
    }
    return compact, profiles


def compact_scalar_current_rate_decomposition(
    state, scalar_budget: dict[str, Any], candidate_raw: dict[str, np.ndarray]
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Resolve -8*pi*d_t(J_S) into the exact existing scalar RHS rate terms."""
    if not scalar_budget.get("reconstruction_passed", False):
        raise RuntimeError("existing scalar RHS term budget did not reconstruct its implementation")
    r = np.asarray(state.grid.centers, dtype=np.float64)
    S = np.asarray(state.scalars.S, dtype=np.float64)
    PS = np.asarray(state.scalars.PS, dtype=np.float64)
    Sr = stress_divergence.d1(state.grid, S, +1)
    kappa = float(8.0 * np.pi)
    sector = scalar_budget["sectors"]["S"]
    field_terms = sector["field_rate_terms"]["raw_signed_arrays"]
    momentum_terms = sector["momentum_rate_terms"]["raw_signed_arrays"]
    profiles: dict[str, np.ndarray] = {}
    rows: dict[str, Any] = {}
    contribution_sum = np.zeros_like(r)

    # For j_S = -P_S * D_r(S), the momentum-constraint matter contribution
    # is -8*pi*j_S. Its RHS tangent is +8*pi*(Pdot_S*D_r(S) + P_S*D_r(Sdot)).
    for name, values in momentum_terms.items():
        term = kappa * np.asarray(values, dtype=np.float64) * Sr
        contribution_sum += term
        profiles[f"scalar_S_current_rate_from_Pdot_{name}"] = term
        rows[f"Pdot_{name}"] = summarize_profile(term, r)
    for name, values in field_terms.items():
        term = kappa * PS * stress_divergence.d1(
            state.grid, np.asarray(values, dtype=np.float64), +1
        )
        contribution_sum += term
        profiles[f"scalar_S_current_rate_from_Sdot_{name}"] = term
        rows[f"Sdot_{name}"] = summarize_profile(term, r)

    actual = np.asarray(candidate_raw["Mdot_sector_S"], dtype=np.float64)
    gap = contribution_sum - actual
    profiles["scalar_S_current_rate_reconstructed"] = contribution_sum
    profiles["scalar_S_current_rate_actual_tangent"] = actual
    profiles["scalar_S_current_rate_reconstruction_gap"] = gap
    report = {
        "observer": "existing audit_total_stress_energy_divergence.scalar_rhs_term_budget plus exact current derivative j_S=-P_S*D_r(S)",
        "scalar_rhs_term_budget_status": scalar_budget.get("status"),
        "scalar_rhs_term_reconstruction_passed": bool(scalar_budget.get("reconstruction_passed")),
        "current_definition": "j_S=-P_S*D_r(S); Mdot_S contribution=-8*pi*d_t(j_S)=+8*pi*(Pdot_S*D_r(S)+P_S*D_r(Sdot))",
        "momentum_rate_terms_contribution": {
            key: rows[key] for key in rows if key.startswith("Pdot_")
        },
        "field_rate_terms_contribution": {
            key: rows[key] for key in rows if key.startswith("Sdot_")
        },
        "reconstructed_Mdot_sector_S": summarize_profile(contribution_sum, r),
        "actual_directional_Mdot_sector_S": summarize_profile(actual, r),
        "reconstruction_gap": summarize_profile(gap, r),
        "reconstruction_gap_max_abs": float(np.max(np.abs(gap))),
        "source_pin_caveat": (
            "This is a term-by-term reconstruction of the current pinned scalar RHS, "
            "not acceptance of a scalar/matter exchange law or total-T conservation."
        ),
    }
    return report, profiles



def compact_scalar_eom_divergence_pairing(
    state,
    rates: dict[str, Any],
    scalar_budget: dict[str, Any],
    candidate_raw: dict[str, np.ndarray],
    sector_profiles: dict[str, np.ndarray],
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Pair the pinned S RHS with pointwise/FV scalar stress divergence.

    This is a discrete operator-pairing diagnostic, not a closure/admission
    oracle. It records both the production-expanded RHS residual and an
    independently formed divergence-operator residual so their finite-grid
    difference is explicit. The FV divergence remains a volume-averaged lane.
    """
    finite = stress_divergence.finite
    grid, geom = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=np.float64)
    dr = float(grid.dr)
    if np.any(np.asarray(geom.beta, dtype=np.float64) != 0.0):
        raise RuntimeError("S-sector EOM pairing currently requires the accepted zero-shift slice")
    if not math.isfinite(dr) or dr <= 0.0:
        raise RuntimeError("invalid radial spacing for S-sector EOM pairing")

    d1 = lambda values, parity=1: stress_divergence.d1(grid, values, parity)
    d2 = lambda values, parity=1: stress_divergence.d2(grid, values, parity)
    alpha = finite("S EOM lapse", geom.alpha)
    a = finite("S EOM metric a", geom.a)
    b = finite("S EOM metric b", geom.b)
    X = finite("S EOM metric X", geom.X)
    K = finite("S EOM trace curvature", geom.K)
    S = finite("S EOM field", state.scalars.S)
    PS = finite("S EOM momentum", state.scalars.PS)
    alpha_r = d1(alpha, +1)
    a_r = d1(a, +1)
    b_r = d1(b, +1)
    X_r = d1(X, +1)
    Sr = d1(S, +1)
    Srr = d2(S, +1)
    gamma_rr_inv = finite("S EOM inverse radial metric", X * X / a)
    sqrt_gamma = finite("S EOM spatial measure", np.sqrt(a) * b * r * r / X**3)
    sqrt_minus_g = finite("S EOM spacetime measure", alpha * sqrt_gamma)

    # Verify that the paired rates are exactly the already reconstructed
    # accepted semidiscrete RHS, rather than a second independent RHS formula.
    budget_s = scalar_budget["sectors"]["S"]
    budget_sdot = finite(
        "S EOM budget field rate",
        budget_s["implementation_rates"]["field_rate"],
    )
    budget_pdot = finite(
        "S EOM budget momentum rate",
        budget_s["implementation_rates"]["momentum_rate"],
    )
    Sdot = finite("S EOM production Sdot", rates["scalars"]["S"])
    Pdot = finite("S EOM production Pdot", rates["scalars"]["PS"])
    if not np.allclose(Sdot, budget_sdot, rtol=5.0e-12, atol=5.0e-18):
        raise RuntimeError("Sdot from directional RHS differs from the pinned scalar budget")
    if not np.allclose(Pdot, budget_pdot, rtol=5.0e-12, atol=5.0e-18):
        raise RuntimeError("Pdot from directional RHS differs from the pinned scalar budget")

    lap_bracket = (
        0.5 * a_r / a - b_r / b + X_r / X - 2.0 / r
    )
    lap_expanded = finite(
        "S production expanded Laplacian",
        gamma_rr_inv * (Srr - Sr * lap_bracket),
    )
    production_spatial_force = finite(
        "S production full spatial force",
        lap_expanded + gamma_rr_inv * (alpha_r / alpha) * Sr,
    )
    eom_production_expanded = finite(
        "S production-expanded equation residual",
        -Pdot / alpha + K * PS + production_spatial_force - S,
    )

    weighted_radial_gradient = finite(
        "S lapse-weighted scalar radial flux",
        alpha * sqrt_gamma * gamma_rr_inv * Sr,
    )
    spatial_divergence_form = finite(
        "S divergence-form spatial operator",
        d1(weighted_radial_gradient, -1) / (alpha * sqrt_gamma),
    )
    eom_divergence_form = finite(
        "S divergence-form equation residual",
        -Pdot / alpha + K * PS + spatial_divergence_form - S,
    )
    operator_pairing_gap = finite(
        "S divergence-form minus production-expanded EOM residual",
        eom_divergence_form - eom_production_expanded,
    )

    # The current and measure tangent are calculated independently from the
    # same accepted state's production rates, then checked against the saved
    # directional tangent and sector-resolved time-density observer.
    a_dot = finite("S EOM a rate", rates["geometry"]["a"])
    b_dot = finite("S EOM b rate", rates["geometry"]["b"])
    X_dot = finite("S EOM X rate", rates["geometry"]["X"])
    dlog_sqrt_gamma = finite(
        "S EOM logarithmic spatial-measure rate",
        0.5 * a_dot / a + b_dot / b - 3.0 * X_dot / X,
    )
    current_j = finite("S stress current j_S", -PS * Sr)
    current_j_dot_rhs = finite(
        "S current tangent reconstructed from scalar RHS",
        -Pdot * Sr - PS * d1(Sdot, +1),
    )
    current_j_dot_directional = finite(
        "S current tangent from saved directional witness",
        -np.asarray(candidate_raw["Mdot_sector_S"], dtype=np.float64) / (8.0 * math.pi),
    )
    measure_current_rate = finite(
        "S current times spatial-measure rate",
        current_j * dlog_sqrt_gamma,
    )
    time_density_rhs = finite(
        "S time-density divergence reconstructed from current and measure rates",
        (current_j_dot_rhs + measure_current_rate) / alpha,
    )

    q_point = finite(
        "S pointwise stress divergence",
        sector_profiles["sector_S_div_T_r"],
    )
    q_fv = finite(
        "S finite-volume stress divergence",
        sector_profiles["sector_S_div_T_r_fv"],
    )
    q_time = finite(
        "S time-density divergence contribution",
        sector_profiles["sector_S_time_density_div_T_r"],
    )
    q_flux = finite(
        "S radial-flux divergence contribution",
        sector_profiles["sector_S_radial_flux_div_T_r"],
    )
    q_metric = finite(
        "S metric-source divergence contribution",
        sector_profiles["sector_S_metric_source_div_T_r"],
    )
    eom_times_Sr = finite(
        "S divergence-form EOM residual times field gradient",
        eom_divergence_form * Sr,
    )
    pointwise_pairing_gap = finite(
        "S pointwise divergence minus EOM-gradient pairing",
        q_point - eom_times_Sr,
    )
    fv_pairing_gap = finite(
        "S FV divergence minus EOM-gradient pairing",
        q_fv - eom_times_Sr,
    )

    eom_time_piece = finite(
        "S EOM time/curvature piece times field gradient",
        (-Pdot / alpha + K * PS) * Sr,
    )
    eom_spatial_potential_piece = finite(
        "S EOM spatial/potential piece times field gradient",
        (spatial_divergence_form - S) * Sr,
    )
    time_pairing_gap = finite(
        "S time-density versus EOM time-piece pairing gap",
        q_time - eom_time_piece,
    )
    spatial_metric_pairing_gap = finite(
        "S radial-flux plus metric-source versus EOM spatial pairing gap",
        q_flux + q_metric - eom_spatial_potential_piece,
    )
    pointwise_pairing_reconstructed = finite(
        "S reconstructed pointwise pairing gap",
        time_pairing_gap + spatial_metric_pairing_gap,
    )
    pointwise_pairing_decomposition_gap = finite(
        "S pointwise pairing decomposition closure",
        pointwise_pairing_gap - pointwise_pairing_reconstructed,
    )
    pointwise_component_gap = finite(
        "S pointwise component sum gap",
        q_point - q_time - q_flux - q_metric,
    )

    face_measure = stress_divergence.spherical_flux_face_values(
        sqrt_minus_g, dr, origin_power=2
    )
    cell_measure = finite(
        "S EOM observer Simpson cell measure",
        (dr / 6.0) * (face_measure[:-1] + 4.0 * sqrt_minus_g + face_measure[1:]),
    )
    if np.any(cell_measure <= 0.0):
        raise RuntimeError("S EOM observer produced a non-positive cell measure")

    # Volume-weighted sums use the observer's own Simpson measure. The FV lane
    # is already a volume-averaged divergence, so measure*q_FV is its cell
    # integral. EOM and pointwise lanes are cell-centred samples with the same
    # measure; these sums are explicit crosswalks, not exact conservation laws.
    masks = {
        "origin_r_lt_2": r < 2.0,
        "annulus_2_le_r_lt_5": (r >= 2.0) & (r < 5.0),
        "annulus_5_le_r_lt_10": (r >= 5.0) & (r < 10.0),
        "annulus_2_le_r_lt_8": (r >= 2.0) & (r < 8.0),
        "domain_r_lt_10": r < 10.0,
    }
    def profile_summary(values: np.ndarray, mask: np.ndarray) -> dict[str, Any]:
        selected = np.asarray(values, dtype=np.float64)[mask]
        if selected.size == 0:
            return {"samples": 0, "rms": None, "max_abs": None, "mean_abs": None}
        return {
            "samples": int(selected.size),
            "rms": float(np.sqrt(np.mean(selected * selected))),
            "max_abs": float(np.max(np.abs(selected))),
            "mean_abs": float(np.mean(np.abs(selected))),
        }
    shells = {}
    for name, mask in masks.items():
        shells[name] = {
            "cell_count": int(np.count_nonzero(mask)),
            "measure_sum": float(np.sum(cell_measure[mask])),
            "integral_EOM_gradient_pairing": float(
                np.sum(cell_measure[mask] * eom_times_Sr[mask])
            ),
            "integral_pointwise_divergence": float(
                np.sum(cell_measure[mask] * q_point[mask])
            ),
            "integral_FV_divergence": float(
                np.sum(cell_measure[mask] * q_fv[mask])
            ),
            "pointwise_minus_EOM_integral": float(
                np.sum(cell_measure[mask] * pointwise_pairing_gap[mask])
            ),
            "FV_minus_EOM_integral": float(
                np.sum(cell_measure[mask] * fv_pairing_gap[mask])
            ),
            "pointwise_pairing_gap": profile_summary(pointwise_pairing_gap, mask),
            "FV_pairing_gap": profile_summary(fv_pairing_gap, mask),
            "EOM_gradient_pairing": profile_summary(eom_times_Sr, mask),
        }

    fixed = {}
    for radius in (1.0, 2.0, 5.0, 8.0):
        if r[0] <= radius <= r[-1]:
            fixed[str(radius)] = {
                "EOM_divergence_form": float(np.interp(radius, r, eom_divergence_form)),
                "EOM_expanded_production_form": float(np.interp(radius, r, eom_production_expanded)),
                "EOM_gradient_pairing": float(np.interp(radius, r, eom_times_Sr)),
                "pointwise_Qr_S": float(np.interp(radius, r, q_point)),
                "FV_Qr_S": float(np.interp(radius, r, q_fv)),
                "pointwise_pairing_gap": float(np.interp(radius, r, pointwise_pairing_gap)),
                "FV_pairing_gap": float(np.interp(radius, r, fv_pairing_gap)),
                "measure_log_rate": float(np.interp(radius, r, dlog_sqrt_gamma)),
            }

    profiles = {
        "scalar_S_production_Sdot": Sdot,
        "scalar_S_production_Pdot": Pdot,
        "scalar_S_laplacian_expanded": lap_expanded,
        "scalar_S_production_spatial_force": production_spatial_force,
        "scalar_S_eom_expanded_residual": eom_production_expanded,
        "scalar_S_eom_divergence_form_residual": eom_divergence_form,
        "scalar_S_eom_operator_pairing_gap": operator_pairing_gap,
        "scalar_S_eom_times_D1S": eom_times_Sr,
        "scalar_S_measure_log_rate": dlog_sqrt_gamma,
        "scalar_S_current_j": current_j,
        "scalar_S_current_j_dot_rhs": current_j_dot_rhs,
        "scalar_S_current_j_dot_directional": current_j_dot_directional,
        "scalar_S_measure_current_rate": measure_current_rate,
        "scalar_S_time_density_divergence_rhs": time_density_rhs,
        "scalar_S_time_density_pairing_gap": time_pairing_gap,
        "scalar_S_spatial_metric_pairing_gap": spatial_metric_pairing_gap,
        "scalar_S_pointwise_pairing_gap": pointwise_pairing_gap,
        "scalar_S_fv_pairing_gap": fv_pairing_gap,
        "scalar_S_pointwise_pairing_decomposition_gap": pointwise_pairing_decomposition_gap,
        "scalar_S_pointwise_component_sum_gap": pointwise_component_gap,
    }
    report = {
        "status": "completed_operator_pairing_diagnostic_only",
        "accepted_time": float(state.t),
        "coordinate_radius_convention": "shell masks select cell centers; cell measures are the same Simpson-reconstructed alpha*sqrt(gamma) measures as the FV observer",
        "equation_conventions": {
            "production_expanded_residual": "-Pdot_S/alpha + K*P_S + L_expanded[S] + gamma^rr*(D1(alpha)/alpha)*D1(S) - Vprime(S)",
            "divergence_form_residual": "-Pdot_S/alpha + K*P_S + D1(alpha*sqrt(gamma)*gamma^rr*D1(S))/(alpha*sqrt(gamma)) - Vprime(S)",
            "stress_identity_continuum": "nabla_mu(T_S^mu_r) = EOM_S * D_r(S)",
            "shift": "zero on the accepted CMC slice; nonzero-shift cases fail closed",
        },
        "production_rhs_reconstruction": {
            "Sdot_matches_scalar_budget": True,
            "Pdot_matches_scalar_budget": True,
            "expanded_residual": profile_summary(eom_production_expanded, np.ones_like(r, dtype=bool)),
            "expanded_residual_max_abs": float(np.max(np.abs(eom_production_expanded))),
        },
        "current_and_measure_rate": {
            "measure_log_rate": profile_summary(dlog_sqrt_gamma, np.ones_like(r, dtype=bool)),
            "current_tangent_RHS_vs_directional_max_abs": float(
                np.max(np.abs(current_j_dot_rhs - current_j_dot_directional))
            ),
            "current_tangent_RHS_vs_directional_rms": float(
                np.sqrt(np.mean((current_j_dot_rhs - current_j_dot_directional) ** 2))
            ),
            "time_density_RHS_vs_observer_max_abs": float(
                np.max(np.abs(time_density_rhs - q_time))
            ),
            "time_density_RHS_vs_observer_rms": float(
                np.sqrt(np.mean((time_density_rhs - q_time) ** 2))
            ),
            "measure_current_rate": profile_summary(measure_current_rate, np.ones_like(r, dtype=bool)),
        },
        "operator_pairing": {
            "divergence_form_minus_production_expanded": profile_summary(
                operator_pairing_gap, np.ones_like(r, dtype=bool)
            ),
            "pointwise_Qr_minus_EOM_times_D1S": profile_summary(
                pointwise_pairing_gap, np.ones_like(r, dtype=bool)
            ),
            "FV_Qr_minus_EOM_times_D1S": profile_summary(
                fv_pairing_gap, np.ones_like(r, dtype=bool)
            ),
            "pointwise_gap_decomposition_max_abs": float(
                np.max(np.abs(pointwise_pairing_decomposition_gap))
            ),
            "pointwise_component_sum_gap_max_abs": float(np.max(np.abs(pointwise_component_gap))),
            "fixed_coordinate_radius": fixed,
        },
        "volume_weighted_shell_comparison": shells,
        "interpretation_limits": [
            "The expanded production EOM residual is expected to close because Pdot is reconstructed from that same RHS; this is a source/operator identity, not independent physical validation.",
            "The divergence-form EOM residual differs from the production-expanded residual by a discrete spatial product-rule/operator-pairing remainder.",
            "The pointwise stress divergence and the FV volume-averaged divergence are different discrete observers and are reported separately.",
            "Volume-weighted sums use cell-center samples for pointwise/EOM fields and the FV observer's Simpson cell measures; these are diagnostic crosswalks, not proof of an integrated conservation law.",
            "No evolution is advanced, no accepted state is mutated, and no source law or scientific admission gate is accepted.",
        ],
        "scientific_admission": False,
    }
    return report, profiles



def compact_sector_expanded_divergence_pairing(
    state,
    rates: dict[str, Any],
    epsilon: float,
    candidate_raw: dict[str, np.ndarray],
    sector_profiles: dict[str, np.ndarray],
    sector_source_profiles: dict[str, np.ndarray],
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Build expanded Q_r separately by sector and pair it with both observers.

    Uses the same accepted base state, plus/minus directional RHS endpoints,
    stress projections, and geometry rates. It is diagnostic-only and does not
    assume that an observed residual establishes a wrong physical source law.
    """
    finite = stress_divergence.finite
    eps = float(epsilon)
    if not math.isfinite(eps) or eps <= 0.0:
        raise ValueError("sector pairing epsilon must be positive and finite")
    try:
        plus_state = SHIFTED_STATES[eps]
        minus_state = SHIFTED_STATES[-eps]
    except KeyError as exc:
        raise RuntimeError(f"sector pairing endpoints missing for epsilon={eps}") from exc

    grid, geom = state.grid, state.geometry
    r = finite("sector pairing radii", grid.centers)
    if np.any(np.asarray(geom.beta, dtype=np.float64) != 0.0):
        raise RuntimeError("sector expanded pairing currently requires the accepted zero-shift slice")
    d1 = lambda values, parity=1: stress_divergence.d1(grid, values, parity)

    projection_triplet = (
        stress_divergence.independent_total_projections(state),
        stress_divergence.independent_total_projections(plus_state),
        stress_divergence.independent_total_projections(minus_state),
    )
    maps = []
    for projection in projection_triplet:
        combined = dict(projection["fluid_parts"])
        combined.update(projection["scalar_parts"])
        maps.append(combined)
    names = ("S", "D", "phi", "dark_matter", "baryons", "radiation")
    raw_names = {
        "S": "S",
        "D": "D",
        "phi": "COSMOS_phi",
        "dark_matter": "dark_matter",
        "baryons": "baryons",
        "radiation": "radiation",
    }
    pointwise_keys = {
        "S": "sector_S_div_T_r",
        "D": "sector_D_div_T_r",
        "phi": "sector_phi_div_T_r",
        "dark_matter": "sector_dark_matter_div_T_r",
        "baryons": "sector_baryons_div_T_r",
        "radiation": "sector_radiation_div_T_r",
    }
    fv_keys = {key: value + "_fv" for key, value in pointwise_keys.items()}

    alpha = finite("sector pairing lapse", geom.alpha)
    alpha_r = d1(alpha, +1)
    gamma_rr = finite("sector pairing gamma_rr", geom.a / (geom.X * geom.X))
    gamma_thth = finite("sector pairing gamma_theta_theta", geom.b * r * r / (geom.X * geom.X))
    gamma_rr_r = d1(gamma_rr, +1)
    gamma_thth_r = d1(gamma_thth, +1)
    areal_log_r = 1.0 / r + 0.5 * d1(geom.b, +1) / geom.b - d1(geom.X, +1) / geom.X
    a_dot = finite("sector pairing a rate", rates["geometry"]["a"])
    b_dot = finite("sector pairing b rate", rates["geometry"]["b"])
    X_dot = finite("sector pairing X rate", rates["geometry"]["X"])
    gamma_rr_log_dot = a_dot / geom.a - 2.0 * X_dot / geom.X
    areal_log_dot = 0.5 * b_dot / geom.b - X_dot / geom.X

    masks = {
        "full_domain": np.ones_like(r, dtype=bool),
        "r_lt_3": r < 3.0,
        "2_le_r_lt_10": (r >= 2.0) & (r < 10.0),
        "2_le_r_lt_8": (r >= 2.0) & (r < 8.0),
    }
    def stats(values: np.ndarray, mask: np.ndarray) -> dict[str, Any]:
        selected = finite("sector-pairing stats input", values)[mask]
        if selected.size == 0:
            return {"samples": 0, "rms": None, "max_abs": None}
        return {
            "samples": int(selected.size),
            "rms": float(np.sqrt(np.mean(selected * selected))),
            "max_abs": float(np.max(np.abs(selected))),
        }

    reports = {}
    out_profiles: dict[str, np.ndarray] = {}
    expanded_sum = np.zeros_like(r)
    pointwise_sum = np.zeros_like(r)
    fv_sum = np.zeros_like(r)
    for name in names:
        raw_sector = raw_names[name]
        base = maps[0][name]
        plus = maps[1][name]
        minus = maps[2][name]
        rho = finite(f"{name} projected rho", base["rho"])
        j = finite(f"{name} projected j", base["j"])
        pr = finite(f"{name} projected radial pressure", base["pr"])
        pt = finite(f"{name} projected tangential pressure", base["pt"])
        j_dot = finite(
            f"{name} directional current tangent",
            (np.asarray(plus["j"], dtype=np.float64) - np.asarray(minus["j"], dtype=np.float64))
            / (2.0 * eps),
        )
        j_dot_from_constraint = finite(
            f"{name} constraint current tangent",
            -np.asarray(candidate_raw[f"Mdot_sector_{raw_sector}"], dtype=np.float64)
            / (8.0 * math.pi),
        )
        expanded = finite(
            f"{name} expanded sector Qr",
            j_dot / alpha
            + (j / alpha) * (0.5 * gamma_rr_log_dot + 2.0 * areal_log_dot)
            + d1(pr, +1)
            + (alpha_r / alpha) * (rho + pr)
            + 2.0 * areal_log_r * (pr - pt),
        )
        pointwise = finite(f"{name} pointwise sector Qr", sector_profiles[pointwise_keys[name]])
        fv = finite(f"{name} finite-volume sector Qr", sector_profiles[fv_keys[name]])
        expected_key = f"sector_{raw_sector}_archive_expected_exchange_Qr"
        expected = finite(
            f"{name} archived exchange comparator",
            sector_source_profiles.get(expected_key, np.zeros_like(r)),
        )
        point_gap = finite(f"{name} pointwise-minus-expanded gap", pointwise - expanded)
        fv_gap = finite(f"{name} FV-minus-expanded gap", fv - expanded)
        q_expected_point = finite(f"{name} pointwise-minus-expected Qr", pointwise - expected)
        q_expected_expanded = finite(f"{name} expanded-minus-expected Qr", expanded - expected)
        tangent_gap = finite(f"{name} projected-current versus constraint-current gap", j_dot - j_dot_from_constraint)
        expanded_sum += expanded
        pointwise_sum += pointwise
        fv_sum += fv
        out_profiles[f"sector_{raw_sector}_expanded_Qr"] = expanded
        out_profiles[f"sector_{raw_sector}_pointwise_minus_expanded_Qr"] = point_gap
        out_profiles[f"sector_{raw_sector}_fv_minus_expanded_Qr"] = fv_gap
        out_profiles[f"sector_{raw_sector}_expanded_minus_expected_exchange_Qr"] = q_expected_expanded
        out_profiles[f"sector_{raw_sector}_jdot_projection_minus_constraint"] = tangent_gap
        fixed = {}
        for radius in (1.0, 2.0, 3.0, 5.0, 8.0):
            if r[0] <= radius <= r[-1]:
                fixed[str(radius)] = {
                    "Q_pointwise": float(np.interp(radius, r, pointwise)),
                    "Q_expanded": float(np.interp(radius, r, expanded)),
                    "Q_fv": float(np.interp(radius, r, fv)),
                    "Q_expected_exchange": float(np.interp(radius, r, expected)),
                    "pointwise_minus_expanded": float(np.interp(radius, r, point_gap)),
                    "expanded_minus_expected": float(np.interp(radius, r, q_expected_expanded)),
                    "pointwise_minus_expected": float(np.interp(radius, r, q_expected_point)),
                }
        reports[raw_sector] = {
            "directional_jdot_vs_constraint_jdot_max_abs": float(np.max(np.abs(tangent_gap))),
            "pointwise_minus_expanded": {key: stats(point_gap, mask) for key, mask in masks.items()},
            "fv_minus_expanded": {key: stats(fv_gap, mask) for key, mask in masks.items()},
            "expanded_minus_expected_exchange": {key: stats(q_expected_expanded, mask) for key, mask in masks.items()},
            "pointwise_minus_expected_exchange": {key: stats(q_expected_point, mask) for key, mask in masks.items()},
            "fixed_coordinate_radius": fixed,
        }
    total_expanded = finite("sum expanded sector Qr", expanded_sum)
    total_candidate = finite(
        "expanded total Brown-source Qr candidate",
        candidate_raw["propagation_Qr_total_stress_divergence"],
    )
    total_expanded_gap = finite("expanded sector-sum minus total candidate", total_expanded - total_candidate)
    total_pointwise = finite("sum pointwise sector Qr", pointwise_sum)
    total_pointwise_gap = finite(
        "pointwise sector-sum minus total combined-flux Qr",
        total_pointwise - np.asarray(sector_profiles["sector_resolved_total_div_T_r"], dtype=np.float64),
    )
    out_profiles["sector_expanded_Qr_sum"] = total_expanded
    out_profiles["sector_expanded_Qr_sum_minus_total_candidate"] = total_expanded_gap
    out_profiles["sector_pointwise_Qr_sum_minus_total"] = total_pointwise_gap
    report = {
        "status": "completed_sector_expanded_operator_pairing_diagnostic_only",
        "epsilon": eps,
        "accepted_time": float(state.t),
        "sectors": reports,
        "sum_expanded_sectors_minus_total_candidate_max_abs": float(np.max(np.abs(total_expanded_gap))),
        "sum_pointwise_sectors_minus_combined_flux_total_max_abs": float(np.max(np.abs(total_pointwise_gap))),
        "interpretation_limits": [
            "The expanded sector formulas use the same RHS endpoints and geometry rates as the full accepted-state observer; they are a diagnostic split of the continuum-expanded Qr formula.",
            "A small sector-pointwise-minus-expanded gap indicates this mismatch is largely discrete operator pairing; it does not by itself accept the archive exchange law.",
            "The production-expanded total residual may close by construction, and the directional j-dot check is arithmetic bookkeeping rather than physical admission.",
            "The expected exchange comparator is provisional and archive-sign/normalization gates remain open.",
            "No evolution is advanced, accepted state is not mutated, and no production/vendor source is changed.",
        ],
        "scientific_admission": False,
    }
    return report, out_profiles


def compact_sector_matter_source_mapping(
    state,
    rates: dict[str, Any],
    candidate_raw: dict[str, np.ndarray],
    sector_profiles: dict[str, np.ndarray],
    coupling: float,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Decompose the existing matter-source composite by frozen stress sector.

    Uses the exact linear matter injection formula already in
    audit_hybrid_origin_integral_patch_newton, sector projections from the
    independent total-stress observer, and sector current/divergence profiles
    measured on the same accepted state and tangent endpoints. The archive-
    locked DM/phi exchange is reported as a separate conditional comparator.
    This does not promote any candidate to the accepted propagation oracle.
    """
    grid = state.grid
    r = np.asarray(grid.centers, dtype=np.float64)
    geom = state.geometry
    alpha = np.asarray(geom.alpha, dtype=np.float64)
    b = np.asarray(geom.b, dtype=np.float64)
    X = np.asarray(geom.X, dtype=np.float64)
    Aa = np.asarray(geom.Aa, dtype=np.float64)
    d1 = lambda values, parity=1: stress_divergence.d1(grid, values, parity)
    chi_r = -0.5 * d1(X, +1) / X
    radial_bracket = 2.0 / r + d1(b, +1) / b

    projected = stress_divergence.independent_total_projections(state)
    if not projected.get("projection_crosscheck_passed", False):
        raise RuntimeError("independent stress-sector projections do not match canonical total")

    sectors = (
        ("S", "S", "S", "scalar_parts"),
        ("D", "D", "D", "scalar_parts"),
        ("COSMOS_phi", "phi", "phi", "scalar_parts"),
        ("dark_matter", "dark_matter", "dark_matter", "fluid_parts"),
        ("baryons", "baryons", "baryons", "fluid_parts"),
        ("radiation", "radiation", "radiation", "fluid_parts"),
    )
    rho_dm = np.asarray(projected["fluid_parts"]["dark_matter"]["rho"], dtype=np.float64)
    phi_r = d1(np.asarray(state.scalars.phi, dtype=np.float64), +1)
    expected_exchange = {
        "dark_matter": float(coupling) * rho_dm * phi_r,
        "COSMOS_phi": -float(coupling) * rho_dm * phi_r,
    }

    profiles: dict[str, np.ndarray] = {}
    rows: dict[str, Any] = {}
    injection_sum = np.zeros_like(r)
    current_sum = np.zeros_like(r)
    divergence_sum = np.zeros_like(r)
    composite_sum = np.zeros_like(r)
    residual_sum = np.zeros_like(r)
    expected_exchange_sum = np.zeros_like(r)

    for raw_sector, projection_name, divergence_name, part_kind in sectors:
        part = projected[part_kind][projection_name]
        rho = np.asarray(part["rho"], dtype=np.float64)
        pr = np.asarray(part["pr"], dtype=np.float64)
        pt = np.asarray(part["pt"], dtype=np.float64)
        aa_source = -(16.0 * np.pi / 3.0) * alpha * (pr - pt)
        k_source = 4.0 * np.pi * alpha * (rho + pr + 2.0 * pt)
        geom_injection = (
            d1(aa_source, +1)
            - (2.0 / 3.0) * d1(k_source, +1)
            + 6.0 * aa_source * chi_r
            + 1.5 * aa_source * radial_bracket
        )

        current_key = f"Mdot_sector_{raw_sector}"
        if current_key not in candidate_raw:
            raise RuntimeError(f"missing exact tangent current contribution: {current_key}")
        current_dot = np.asarray(candidate_raw[current_key], dtype=np.float64)
        div_key = f"sector_{divergence_name}_div_T_r"
        if div_key not in sector_profiles:
            raise RuntimeError(f"missing sector divergence profile: {div_key}")
        q_r = np.asarray(sector_profiles[div_key], dtype=np.float64)
        divergence_source = -8.0 * np.pi * alpha * q_r
        composite = geom_injection + current_dot
        residual = composite - divergence_source
        injection_sum += geom_injection
        current_sum += current_dot
        divergence_sum += q_r
        composite_sum += composite
        residual_sum += residual

        expected_q = expected_exchange.get(raw_sector, np.zeros_like(r))
        exchange_residual = q_r - expected_q
        expected_exchange_sum += expected_q

        rows[raw_sector] = {
            "projection": "independent sector stress projection; normalized consistently with total T",
            "geometry_matter_injection": summarize_profile(geom_injection, r),
            "current_tangent_minus_8pi_jdot": summarize_profile(current_dot, r),
            "matter_source_composite": summarize_profile(composite, r),
            "pointwise_div_T_r": summarize_profile(q_r, r),
            "stress_divergence_source_minus_8pi_alpha_Qr": summarize_profile(divergence_source, r),
            "source_mapping_residual_composite_minus_divergence_source": summarize_profile(residual, r),
            "source_mapping_residual_max_abs": float(np.max(np.abs(residual))),
            "archive_locked_exchange_comparator": {
                "law": (
                    "+beta*rho_DM*D_r(phi)" if raw_sector == "dark_matter"
                    else "-beta*rho_DM*D_r(phi)" if raw_sector == "COSMOS_phi"
                    else "zero source for this uncoupled sector"
                ),
                "expected_Qr": summarize_profile(expected_q, r),
                "measured_minus_expected_Qr": summarize_profile(exchange_residual, r),
                "measured_minus_expected_Qr_max_abs": float(np.max(np.abs(exchange_residual))),
            },
        }
        for suffix, values in (
            ("geometry_matter_injection", geom_injection),
            ("current_tangent_minus_8pi_jdot", current_dot),
            ("matter_source_composite", composite),
            ("pointwise_div_T_r", q_r),
            ("divergence_source_minus_8pi_alpha_Qr", divergence_source),
            ("source_mapping_residual", residual),
            ("archive_expected_exchange_Qr", expected_q),
            ("archive_exchange_residual", exchange_residual),
        ):
            profiles[f"sector_{raw_sector}_{suffix}"] = np.asarray(values, dtype=np.float64).copy()

    full_injection = np.asarray(
        candidate_raw["propagation_Mgeom_dot_matter_injection"], dtype=np.float64
    )
    full_composite = np.asarray(
        candidate_raw["propagation_actual_matter_source_composite"], dtype=np.float64
    )
    full_residual = np.asarray(
        candidate_raw["propagation_matter_source_mapping_residual"], dtype=np.float64
    )
    # Keep the existing expanded Q_r candidate distinct from the independently
    # evaluated pointwise covariant divergence. Their difference is the already
    # measured radial-flux product-rule gap, not a failed sector-sum identity.
    full_q_expanded = np.asarray(
        candidate_raw["propagation_Qr_total_stress_divergence"], dtype=np.float64
    )
    full_q_pointwise = np.asarray(
        sector_profiles["sector_resolved_total_div_T_r"], dtype=np.float64
    )
    full_residual_pointwise = full_composite + 8.0 * np.pi * alpha * full_q_pointwise
    checks = {
        "sum_sector_geometry_injection_minus_existing_total_max_abs": float(
            np.max(np.abs(injection_sum - full_injection))
        ),
        "sum_sector_current_minus_composite_minus_injection_max_abs": float(
            np.max(np.abs(current_sum - (full_composite - full_injection)))
        ),
        "sum_sector_div_T_r_minus_pointwise_total_Qr_max_abs": float(
            np.max(np.abs(divergence_sum - full_q_pointwise))
        ),
        "sum_sector_composite_minus_existing_total_max_abs": float(
            np.max(np.abs(composite_sum - full_composite))
        ),
        "sum_sector_mapping_residual_minus_pointwise_total_residual_max_abs": float(
            np.max(np.abs(residual_sum - full_residual_pointwise))
        ),
        "pointwise_total_minus_expanded_Qr_operator_gap_max_abs": float(
            np.max(np.abs(full_q_pointwise - full_q_expanded))
        ),
        "pointwise_mapping_residual_minus_expanded_candidate_residual_max_abs": float(
            np.max(np.abs(full_residual_pointwise - full_residual))
        ),
        "sum_archive_locked_DM_phi_exchange_Qr_max_abs": float(
            np.max(np.abs(expected_exchange_sum))
        ),
    }
    profiles.update({
        "sector_mapping_sum_geometry_injection": injection_sum,
        "sector_mapping_sum_current_tangent": current_sum,
        "sector_mapping_sum_div_T_r": divergence_sum,
        "sector_mapping_sum_matter_source_composite": composite_sum,
        "sector_mapping_sum_source_mapping_residual": residual_sum,
        "sector_mapping_full_total_matter_source_composite": full_composite,
        "sector_mapping_full_total_source_mapping_residual": full_residual,
    })
    compact = {
        "observer": "existing independent_total_projections + existing Mgeom matter-injection source formula + the saved exact tangent currents and sector div_T_r profiles",
        "accepted_time": float(state.t),
        "beta_dm_coupling": float(coupling),
        "archive_exchange_convention": {
            "DM": "nabla_mu T_DM^{mu nu}=+beta*rho_DM*nabla^nu(phi)",
            "COSMOS_phi": "nabla_mu T_phi^{mu nu}=-beta*rho_DM*nabla^nu(phi)",
            "source": "GEAR_V5_5_MASTER_CHECKPOINT_2026-10-05 and the explicit locked-sign statement in GEAR_V5_5_COVARIANT_GLOBAL_GATE_2026-10-05_v2",
            "normalization_note": "rho_DM is taken from the independently reconstructed sector projection in the same normalized stress-energy units used by the divergence observer",
        },
        "sectors": rows,
        "sector_additivity_closure": checks,
        "sector_additivity_all_below_1e-10": all(
            value < 1.0e-10 for key, value in checks.items()
            if key not in (
                "pointwise_total_minus_expanded_Qr_operator_gap_max_abs",
                "pointwise_mapping_residual_minus_expanded_candidate_residual_max_abs",
                "sum_archive_locked_DM_phi_exchange_Qr_max_abs",
            )
        ),
        "interpretation_limit": (
            "This is a sector-resolved diagnostic of the existing provisional source mapping. "
            "The archive-locked DM/phi exchange is compared but not assumed to pass. "
            "The pointwise covariant divergence and expanded Q_r candidate are compared separately; their observed product-rule gap is recorded, not classified as a failure of sector additivity.",
            "No source law or discrete divergence convention is admitted by this result."
        ),
        "scientific_admission": False,
    }
    return compact, profiles



def _nested_array_digest(value: Any) -> str:
    """Stable digest for the frozen semi-discrete RHS-rate bundle."""
    digest = hashlib.sha256()

    def visit(item: Any, path: str) -> None:
        digest.update(path.encode("utf-8"))
        if isinstance(item, dict):
            for key in sorted(item):
                visit(item[key], f"{path}/{key}")
            return
        array = np.asarray(item, dtype=np.float64)
        if not np.all(np.isfinite(array)):
            raise RuntimeError(f"non-finite value in RHS-rate bundle at {path}")
        digest.update(repr(array.shape).encode("ascii"))
        digest.update(array.tobytes())

    visit(value, "rates")
    return digest.hexdigest()


def _same_endpoint_fields_except_alpha(left: Any, right: Any) -> bool:
    """The controlled endpoint pair may differ only in its recomputed lapse."""
    for name in GEOMETRY_NAMES:
        if name == "alpha":
            continue
        if not np.array_equal(
            np.asarray(getattr(left.geometry, name)),
            np.asarray(getattr(right.geometry, name)),
        ):
            return False
    for name in SCALAR_NAMES:
        if not np.array_equal(
            np.asarray(getattr(left.scalars, name)),
            np.asarray(getattr(right.scalars, name)),
        ):
            return False
    for species_name in ("dark_matter", "baryons", "radiation"):
        left_species = getattr(left.matter, species_name)
        right_species = getattr(right.matter, species_name)
        for name in MATTER_FIELDS:
            if not np.array_equal(
                np.asarray(getattr(left_species, name)),
                np.asarray(getattr(right_species, name)),
            ):
                return False
    return True


def same_base_dual_lapse_response_audit(
    resolution: int,
    joint_arrays: dict[str, np.ndarray],
    vacuum: Any,
    expected_base_digest: str,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Pair endpoint lapse policies on one unchanged production-baseline state.

    The saved accepted base state and semi-discrete RHS direction are held fixed.
    Only the lapse recomputation used for +/-epsilon endpoint states is changed.
    This is a sensitivity diagnostic, not a choice/admission of either convention.
    """
    global FULL_METRIC_LAPSE_SHADOW
    prefix = f"N{resolution}_production_baseline_final_"
    kernel, state = reconstruct_state(resolution, joint_arrays, prefix)
    base_digest = tangent.state_digest(state)
    if base_digest != expected_base_digest:
        raise RuntimeError(
            "same-base lapse audit did not reconstruct the exact saved production-baseline state"
        )
    if abs(float(state.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("same-base lapse audit received an unaccepted state")

    original_policy = FULL_METRIC_LAPSE_SHADOW
    original_digest = base_digest
    rate_digest = _nested_array_digest({
        name: value
        for name, value in propagation.rhs_direction(kernel, state).items()
        if name in ("geometry", "scalars", "matter")
    })
    radii = np.asarray(state.grid.centers, dtype=np.float64)
    mode_configs = (("production_cmc", False), ("regular_origin_cmc", True))
    mode_raw: dict[str, dict[float, dict[str, np.ndarray]]] = {}
    endpoints: dict[str, dict[float, dict[str, Any]]] = {}
    mode_summaries: dict[str, Any] = {}
    profile_store: dict[str, np.ndarray] = {}
    terms_to_save = (
        "propagation_Mdot_measured_minus_predicted_covariant_published_Mr",
        "propagation_Mtotal_measured_minus_Brown_plus_actual_matter_composite",
        "propagation_term_connection_second",
        "propagation_term_connection_first",
        "propagation_term_H_gradient",
        "propagation_term_KM_raised_M",
        "propagation_term_alpha_H",
        "propagation_term_index_conversion",
        "propagation_published_Mr_matter_source",
        "propagation_matter_source_mapping_residual",
        "propagation_exact_chain_rule_residual",
    )
    endpoint_non_lapse_equal = True

    try:
        tangent.shifted_state = shifted_state_with_matched_lapse
        propagation.shifted_state = shifted_state_with_matched_lapse
        for mode, uses_regular_lapse in mode_configs:
            FULL_METRIC_LAPSE_SHADOW = uses_regular_lapse
            SHIFTED_STATES.clear()
            current_rate_digest = _nested_array_digest({
                name: value
                for name, value in propagation.rhs_direction(kernel, state).items()
                if name in ("geometry", "scalars", "matter")
            })
            if current_rate_digest != rate_digest:
                raise RuntimeError("the base RHS-rate bundle changed between lapse-policy arms")

            report, raw = propagation.directional_tangent(kernel, state, vacuum)
            if report.get("state_unchanged") is not True:
                raise RuntimeError(f"{mode}: propagation observer mutated the shared accepted state")
            if tangent.state_digest(state) != original_digest:
                raise RuntimeError(f"{mode}: shared accepted base state changed during observation")

            snapshots: dict[float, dict[str, Any]] = {}
            profile_summaries: dict[str, Any] = {}
            for epsilon in EPSILONS:
                eps = float(epsilon)
                if eps not in SHIFTED_STATES or -eps not in SHIFTED_STATES:
                    raise RuntimeError(f"{mode}: missing matched +/- endpoint at epsilon={eps}")
                plus_state = SHIFTED_STATES[eps]
                minus_state = SHIFTED_STATES[-eps]
                snapshots[eps] = {"plus": plus_state, "minus": minus_state}
                eps_tag = f"{eps:.0e}".replace("+", "")

                for sign, endpoint in (("plus", plus_state), ("minus", minus_state)):
                    profile_store[
                        f"N{resolution}_same_base_lapse_{mode}_eps{eps_tag}_{sign}_alpha"
                    ] = np.asarray(endpoint.geometry.alpha, dtype=np.float64).copy()

                for name in terms_to_save:
                    if name in raw[eps]:
                        profile_store[
                            f"N{resolution}_same_base_lapse_{mode}_eps{eps_tag}_{name}"
                        ] = np.asarray(raw[eps][name], dtype=np.float64).copy()
                        profile_summaries.setdefault(name, {})[eps_tag] = summarize_profile(
                            raw[eps][name], radii
                        )

            endpoints[mode] = snapshots
            mode_raw[mode] = raw
            mode_summaries[mode] = {
                "lapse_treatment": (
                    "production CMC solve at perturbed endpoints"
                    if not uses_regular_lapse
                    else "copied regular-origin CMC row at perturbed endpoints"
                ),
                "shared_base_state_sha256": base_digest,
                "base_rhs_rate_bundle_sha256": current_rate_digest,
                "observer_state_unchanged": report["state_unchanged"],
                "profiles": profile_summaries,
            }

        for epsilon in EPSILONS:
            eps = float(epsilon)
            prod = endpoints["production_cmc"][eps]
            reg = endpoints["regular_origin_cmc"][eps]
            for sign in ("plus", "minus"):
                endpoint_non_lapse_equal = (
                    _same_endpoint_fields_except_alpha(prod[sign], reg[sign])
                    and endpoint_non_lapse_equal
                )

        if not endpoint_non_lapse_equal:
            raise RuntimeError(
                "controlled lapse-policy comparison changed non-lapse endpoint fields"
            )
        if tangent.state_digest(state) != original_digest:
            raise RuntimeError("shared accepted base state changed after dual-lapse comparison")

        comparison: dict[str, Any] = {}
        for epsilon in EPSILONS:
            eps = float(epsilon)
            eps_tag = f"{eps:.0e}".replace("+", "")
            prod = endpoints["production_cmc"][eps]
            reg = endpoints["regular_origin_cmc"][eps]
            alpha_dot_prod = (
                np.asarray(prod["plus"].geometry.alpha, dtype=np.float64)
                - np.asarray(prod["minus"].geometry.alpha, dtype=np.float64)
            ) / (2.0 * eps)
            alpha_dot_reg = (
                np.asarray(reg["plus"].geometry.alpha, dtype=np.float64)
                - np.asarray(reg["minus"].geometry.alpha, dtype=np.float64)
            ) / (2.0 * eps)
            delta_alpha_dot = alpha_dot_reg - alpha_dot_prod
            delta_radial_alpha_dot = stress_divergence.d1(
                state.grid, delta_alpha_dot, +1
            )
            delta_plus_alpha = (
                np.asarray(reg["plus"].geometry.alpha, dtype=np.float64)
                - np.asarray(prod["plus"].geometry.alpha, dtype=np.float64)
            )
            delta_minus_alpha = (
                np.asarray(reg["minus"].geometry.alpha, dtype=np.float64)
                - np.asarray(prod["minus"].geometry.alpha, dtype=np.float64)
            )
            for name, values in (
                ("delta_plus_alpha", delta_plus_alpha),
                ("delta_minus_alpha", delta_minus_alpha),
                ("delta_directional_alpha_rate", delta_alpha_dot),
                ("delta_radial_directional_alpha_rate", delta_radial_alpha_dot),
            ):
                profile_store[
                    f"N{resolution}_same_base_lapse_eps{eps_tag}_{name}"
                ] = np.asarray(values, dtype=np.float64).copy()

            changes: dict[str, Any] = {}
            for name in terms_to_save:
                if name not in mode_raw["production_cmc"][eps] or name not in mode_raw["regular_origin_cmc"][eps]:
                    continue
                delta = (
                    np.asarray(mode_raw["regular_origin_cmc"][eps][name], dtype=np.float64)
                    - np.asarray(mode_raw["production_cmc"][eps][name], dtype=np.float64)
                )
                profile_store[
                    f"N{resolution}_same_base_lapse_eps{eps_tag}_delta_{name}"
                ] = delta.copy()
                changes[name] = summarize_profile(delta, radii)

            comparison[eps_tag] = {
                "production_cmc_plus_minus_lapse_difference": {
                    "plus": summarize_profile(delta_plus_alpha, radii),
                    "minus": summarize_profile(delta_minus_alpha, radii),
                },
                "regular_origin_minus_production_directional_alpha_rate": summarize_profile(
                    delta_alpha_dot, radii
                ),
                "regular_origin_minus_production_radial_gradient_of_alpha_rate": summarize_profile(
                    delta_radial_alpha_dot, radii
                ),
                "regular_origin_minus_production_signed_term_profiles": changes,
            }

        summary = {
            "scope": (
                "One exact saved production-baseline t=0.60 state per resolution; same RHS-rate bundle, "
                "same +/-epsilon non-lapse endpoint fields; only endpoint lapse recomputation path changes."
            ),
            "resolution": resolution,
            "accepted_time": float(state.t),
            "shared_base_state_sha256": base_digest,
            "base_rhs_rate_bundle_sha256": rate_digest,
            "same_base_state_immutable": tangent.state_digest(state) == original_digest,
            "non_lapse_endpoint_fields_identical": endpoint_non_lapse_equal,
            "mode_results": mode_summaries,
            "paired_endpoint_differences": comparison,
            "scientific_admission": False,
        }
        return summary, profile_store
    finally:
        FULL_METRIC_LAPSE_SHADOW = original_policy
        SHIFTED_STATES.clear()


def main() -> int:
    if RESOLUTION not in (60, 120, 240):
        raise ValueError(f"unsupported resolution {RESOLUTION}")
    if SOURCE_COMMIT != EXPECTED_SOURCE or VENDOR_COMMIT != EXPECTED_VENDOR:
        raise RuntimeError("helper source/vendor constants do not match the archived runs")
    pin = verify_source_pins()
    if pin.get("source_commit") != EXPECTED_SOURCE or pin.get("vendor_commit") != EXPECTED_VENDOR:
        raise RuntimeError("live source/vendor pin verification failed")

    joint_manifest = check_artifact(JOINT_DIR, RESOLUTION, full_profile=False)
    full_manifest = check_artifact(FULL_DIR, RESOLUTION, full_profile=True)
    joint_arrays = load_profile_arrays(JOINT_DIR)
    full_arrays = load_profile_arrays(FULL_DIR)

    _, vacuum, _ = adapter.vendor_modules()
    OUT.mkdir(parents=True, exist_ok=True)
    summaries = []
    saved_profiles: dict[str, np.ndarray] = {}
    state_hashes = {}
    global FULL_METRIC_LAPSE_SHADOW

    for label, arm_name, uses_regular_lapse in ARMS:
        if uses_regular_lapse:
            prefix_data = full_arrays if arm_name == "full_continuous_metric_profile" else joint_arrays
            prefix = f"N{RESOLUTION}_{arm_name}_final_"
        else:
            prefix_data = joint_arrays
            prefix = f"N{RESOLUTION}_{arm_name}_final_"

        FULL_METRIC_LAPSE_SHADOW = bool(uses_regular_lapse)
        kernel, state = reconstruct_state(RESOLUTION, prefix_data, prefix)
        initial_digest = tangent.state_digest(state)
        base_rates = tangent.rhs_direction(kernel, state)
        hll_snapshot, hll_profiles = compact_hll_operator_snapshot(
            state, base_rates
        )
        scalar_rhs_budget = stress_divergence.scalar_rhs_term_budget(
            state, float(hll_snapshot["coupling_parameter"])
        )
        if tangent.state_digest(state) != initial_digest:
            raise RuntimeError("scalar RHS term budget mutated its accepted base state")
        tangent.shifted_state = shifted_state_with_matched_lapse
        propagation.shifted_state = shifted_state_with_matched_lapse

        SHIFTED_STATES.clear()
        rates = tangent.rhs_direction(kernel, state)
        report, raw = propagation.directional_tangent(kernel, state, vacuum)
        if report.get("state_unchanged") is not True:
            raise RuntimeError(f"{label}: propagation observer mutated its frozen accepted state")
        final_digest = tangent.state_digest(state)
        if initial_digest != final_digest:
            raise RuntimeError(f"{label}: accepted state digest changed during observation")
        state_hashes[label] = initial_digest

        radii = np.asarray(state.grid.centers, dtype=np.float64)
        arm_rows = {
            "arm": label,
            "source_profile_arm": arm_name,
            "lapse_treatment": "copied regular-origin CMC row" if uses_regular_lapse else "production CMC row",
            "accepted_time": float(state.t),
            "state_digest_sha256": initial_digest,
            "constraint_baseline": report["baseline_constraints"],
            "production_hll_operator_snapshot": hll_snapshot,
            "epsilon_results": [],
        }
        for item in report["tangent_by_epsilon"]:
            eps = float(item["epsilon"])
            compact = compact_candidate(item, raw, radii)
            compact["signed_momentum_budget"] = item["signed_momentum_budget"]
            compact["exact_discrete_chain_rule"] = item["independent_constraint_propagation_candidate"]["exact_discrete_chain_rule"]
            sector_report, sector_profiles = compact_sector_divergence(
                state, eps, raw[eps]
            )
            operator_report, operator_profiles = operator_pairing_decomposition(
                kernel, state, rates, eps, raw[eps]
            )
            scalar_current_report, scalar_current_profiles = compact_scalar_current_rate_decomposition(
                state, scalar_rhs_budget, raw[eps]
            )
            scalar_eom_pairing_report, scalar_eom_pairing_profiles = compact_scalar_eom_divergence_pairing(
                state, rates, scalar_rhs_budget, raw[eps], sector_profiles
            )
            sector_source_report, sector_source_profiles = compact_sector_matter_source_mapping(
                state, rates, raw[eps], sector_profiles,
                float(hll_snapshot["coupling_parameter"])
            )
            sector_operator_report, sector_operator_profiles = compact_sector_expanded_divergence_pairing(
                state, rates, eps, raw[eps], sector_profiles, sector_source_profiles
            )
            compact["sector_resolved_stress_energy_divergence"] = sector_report
            compact["pointwise_vs_expanded_Qr_discrete_gap_decomposition"] = operator_report
            compact["scalar_S_current_RHS_term_decomposition"] = scalar_current_report
            compact["scalar_S_EOM_divergence_operator_pairing"] = scalar_eom_pairing_report
            compact["per_sector_matter_source_mapping"] = sector_source_report
            compact["per_sector_expanded_divergence_operator_pairing"] = sector_operator_report
            arm_rows["epsilon_results"].append(compact)
            for name, values in raw[eps].items():
                saved_profiles[f"N{RESOLUTION}_{label}_eps{eps:.0e}_{name}"] = np.asarray(values, dtype=np.float64).copy()
            for name, values in sector_profiles.items():
                saved_profiles[
                    f"N{RESOLUTION}_{label}_eps{eps:.0e}_{name}"
                ] = np.asarray(values, dtype=np.float64).copy()
            for name, values in operator_profiles.items():
                saved_profiles[
                    f"N{RESOLUTION}_{label}_eps{eps:.0e}_{name}"
                ] = np.asarray(values, dtype=np.float64).copy()
            for name, values in scalar_current_profiles.items():
                saved_profiles[
                    f"N{RESOLUTION}_{label}_eps{eps:.0e}_{name}"
                ] = np.asarray(values, dtype=np.float64).copy()
            for name, values in scalar_eom_pairing_profiles.items():
                saved_profiles[
                    f"N{RESOLUTION}_{label}_eps{eps:.0e}_{name}"
                ] = np.asarray(values, dtype=np.float64).copy()
            for name, values in sector_source_profiles.items():
                saved_profiles[
                    f"N{RESOLUTION}_{label}_eps{eps:.0e}_{name}"
                ] = np.asarray(values, dtype=np.float64).copy()
            for name, values in sector_operator_profiles.items():
                saved_profiles[
                    f"N{RESOLUTION}_{label}_eps{eps:.0e}_{name}"
                ] = np.asarray(values, dtype=np.float64).copy()
        for name, values in hll_profiles.items():
            saved_profiles[f"N{RESOLUTION}_{label}_{name}"] = np.asarray(values, dtype=np.float64).copy()
        summaries.append(arm_rows)

    same_base_lapse_report, same_base_lapse_profiles = same_base_dual_lapse_response_audit(
        RESOLUTION, joint_arrays, vacuum, state_hashes["production_baseline"]
    )
    saved_profiles.update(same_base_lapse_profiles)

    report_path = OUT / "report.json"
    result = {
        "schema": "full-profile-accepted-state-propagation-candidate/v8",
        "status": "completed_diagnostic_only",
        "purpose": "Apply the existing convention/index-sensitive propagation candidate and its existing additive sector-resolved total-stress-energy divergence observer to source-pinned accepted t=0.60 states from baseline, partial-joint, and full-continuous-metric arms.",
        "source_pins": {
            "source_commit": EXPECTED_SOURCE,
            "vendor_commit": EXPECTED_VENDOR,
            "live_pin_check": pin,
            "protected_engine_vendor_dependency_diff_empty": pin.get("protected_engine_vendor_dependency_diff_empty", False),
            "joint_artifact": joint_manifest,
            "full_profile_artifact": full_manifest,
        },
        "runtime": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform()},
        "configuration": {
            "resolution": RESOLUTION,
            "r_max": 80.0,
            "accepted_time": TARGET_TIME,
            "arms": [x[0] for x in ARMS],
            "directional_epsilons": list(EPSILONS),
            "lapse_policy": {
                "production_baseline": "the production CMC lapse solve at perturbed endpoints",
                "partial_joint": "the copied regular-origin CMC row at perturbed endpoints",
                "full_continuous_profile": "the copied regular-origin CMC row at perturbed endpoints",
            },
        },
        "state_sha256": state_hashes,
        "same_base_dual_lapse_response": same_base_lapse_report,
        "arms": summaries,
        "acceptance": {
            "all_source_pins_match": True,
            "input_artifacts_sha256_verified": True,
            "all_frozen_state_observers_immutable": True,
            "all_reported_profiles_finite": True,
            "production_hll_operator_snapshots_paired_on_exact_accepted_states": True,
            "scalar_S_current_rate_terms_reconstructed": all(
                item.get("scalar_S_current_RHS_term_decomposition", {}).get(
                    "scalar_rhs_term_reconstruction_passed"
                ) is True
                and item.get("scalar_S_current_RHS_term_decomposition", {}).get(
                    "reconstruction_gap_max_abs", float("inf")
                ) < 1.0e-9
                for arm in summaries
                for item in arm["epsilon_results"]
            ),
            "sector_source_mapping_additivity_closes": all(
                item.get("per_sector_matter_source_mapping", {}).get(
                    "sector_additivity_all_below_1e-10"
                ) is True
                for arm in summaries
                for item in arm["epsilon_results"]
            ),
            "archive_locked_DM_phi_exchange_tested_not_assumed": True,
            "same_base_dual_lapse_endpoint_pairing_complete": (
                same_base_lapse_report.get("same_base_state_immutable") is True
                and same_base_lapse_report.get("non_lapse_endpoint_fields_identical") is True
            ),
            "full_matter_coupled_propagation_identity_closed": False,
            "matter_source_sign_and_normalization_reviewed": False,
            "reference_metric_connection_coupling_convention_frozen": False,
            "total_stress_energy_sector_conservation_closed": False,
            "production_source_or_vendor_changed": False,
            "scientific_admission": False,
        },
        "interpretation_limits": [
            "This evaluates the existing Brown Eq. (11b)-based candidates, not a newly asserted acceptance oracle.",
            "The direct published-covariant-M_r candidate uses matter-inclusive H/M/C with a provisional -8*pi*alpha*Q_r source.",
            "The geometric-vacuum comparator uses H_geom/M_geom and separately compares the actual matter source composite with the independently reconstructed stress-divergence source.",
            "The source-mapping, scalar spatial-operator, variable-mass exchange, and full sector conservation identities remain incomplete/unreviewed.",
            "Production HLL face/source snapshots are direct observations of the pinned species RHS on the accepted state; they do not by themselves prove total stress-energy conservation or validate the conditional coupling law.",
            "Sector-resolved div_T_r is an existing component-accounting observer on the same directional endpoint states; its reconstruction result is measured, not treated as an acceptance threshold.",
            "Centre cells 0/1 are reported separately from fixed physical radii and annuli; no cell-index comparison is presented as a matched-radius comparison.",
            "The same-base dual-lapse comparison varies only the endpoint lapse recomputation path; it does not admit either convention.",
            "No evolution is advanced and no initializer, engine, vendor, dependency, or production equation is modified.",
        ],
    }
    report_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    profile_path = OUT / "profiles.npz"
    np.savez_compressed(profile_path, **saved_profiles)
    (OUT / "SHA256SUMS.txt").write_text(
        f"{sha256(report_path)}  report.json\n{sha256(profile_path)}  profiles.npz\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": result["status"],
        "resolution": RESOLUTION,
        "output_dir": str(OUT),
        "state_sha256": state_hashes,
        "direct_Mr_candidate_rel_rms": {
            row["arm"]: [
                item["direct_published_covariant_Mr_candidate"]["relative_rms_error_full_domain"]
                for item in row["epsilon_results"]
            ] for row in summaries
        },
        "matter_source_mapping_rel_rms": {
            row["arm"]: [
                item["geometric_Brown_plus_actual_source_candidate"]["source_mapping_relative_rms"]
                for item in row["epsilon_results"]
            ] for row in summaries
        },
        "scientific_admission": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
