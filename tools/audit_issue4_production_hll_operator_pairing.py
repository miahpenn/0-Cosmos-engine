#!/usr/bin/env python3
"""Diagnostic-only comparison to the pinned production per-species HLL operator.

This observer does not advance the production state. It calls the pinned
species RHS component function on a frozen state, records per-species HLL
face fluxes and source splits, and independently contracts each fluid's
contravariant stress with the radial metric derivative. A conditional
covector-source comparison is reported without accepting the physical premise.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import traceback
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
DEFAULT_EPSILONS = (2.0e-5, 1.0e-5)
BETA_DM_COUPLED = -0.04

# Load the already reviewed source-pinned state/RHS observer. It pins its own
# source blobs and inserts the production checkout before importing engine/.
AUDIT_PATH = ROOT / "tools" / "audit_total_stress_energy_divergence.py"
_AUDIT_SPEC = importlib.util.spec_from_file_location(
    "_issue4_pinned_total_t_divergence", AUDIT_PATH
)
if _AUDIT_SPEC is None or _AUDIT_SPEC.loader is None:
    raise ImportError(f"cannot load pinned total-T auditor: {AUDIT_PATH}")
AUDIT = importlib.util.module_from_spec(_AUDIT_SPEC)
sys.modules[_AUDIT_SPEC.name] = AUDIT
_AUDIT_SPEC.loader.exec_module(AUDIT)

from engine import production_kernel as production_kernel_module  # noqa: E402
from engine.matter_system import (  # noqa: E402
    Species,
    _metric_arrays,
    _species_rhs_components,
    geometry_metric_derivatives,
    primitives,
    species_projection,
    initialize_dust,
    initialize_radiation,
)
from engine.scalar_system import KAPPA  # noqa: E402
from engine.valencia import stress_tensor_contravariant  # noqa: E402


def finite(name: str, value: Any) -> np.ndarray:
    out = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return out


def summarize_faces(values: Any, dr: float) -> dict[str, Any]:
    x = finite("face profile summary", values)
    radii = np.arange(x.size, dtype=np.float64) * float(dr)
    mask = (radii >= 2.0) & (radii <= 10.0)
    i = int(np.argmax(np.abs(x)))
    return {
        "max_abs": float(abs(x[i])),
        "signed_at_max_abs": float(x[i]),
        "max_index": i,
        "max_radius": float(radii[i]),
        "rms": float(np.sqrt(np.mean(x * x))),
        "annulus_2_to_10_rms": float(np.sqrt(np.mean(x[mask] ** 2))) if np.any(mask) else None,
        "face_0": float(x[0]),
        "face_1": float(x[1]) if x.size > 1 else None,
        "face_2": float(x[2]) if x.size > 2 else None,
    }


def _metric_bundle(state):
    metric_slice = AUDIT.metric_slice_from_q(state.grid, state.geometry)
    d1_fn = lambda values, parity: AUDIT.d1(state.grid, values, parity)
    derivatives = geometry_metric_derivatives(metric_slice, d1_fn)
    metrics = _metric_arrays(metric_slice)
    return metric_slice, derivatives, metrics, d1_fn


def _species_states(state):
    result = [
        ("dark_matter", state.matter.dark_matter, Species.DARK_MATTER),
        ("baryons", state.matter.baryons, Species.BARYON),
        ("radiation", state.matter.radiation, Species.RADIATION),
    ]
    if getattr(state.matter, "charged_dust", None) is not None:
        result.append(("charged_dust", state.matter.charged_dust, Species.BARYON))
    return result


def operator_snapshot(state, coupling: float, rates: dict[str, Any]) -> dict[str, Any]:
    """Record pinned HLL fluxes and sources plus independent radial contractions."""
    if np.any(np.asarray(state.geometry.beta) != 0.0):
        raise ValueError("this production-source witness supports only the native zero-shift slice")
    grid = state.grid
    metric_slice, derivatives, metrics, _ = _metric_bundle(state)
    dr = float(grid.dr)
    r = finite("cell radii", grid.centers)
    alpha = finite("lapse", state.geometry.alpha)
    sqrt_gamma = finite("sqrt gamma", np.asarray([m.sqrt_gamma for m in metrics]))
    dphi_t = finite("partial_t(phi) from the pinned RHS", rates["scalars"]["phi"])
    dphi_r = AUDIT.d1(grid, state.scalars.phi, +1)
    projection_check = AUDIT.independent_total_projections(state)
    if not projection_check["projection_crosscheck_passed"]:
        raise RuntimeError("independent scalar+fluid projection disagrees with canonical assembly")
    scalar_sector_summary: dict[str, Any] = {}
    raw_scalar_arrays: dict[str, Any] = {}
    for scalar_name, scalar_channels in projection_check["scalar_parts"].items():
        scalar_sector_summary[scalar_name] = {}
        raw_scalar_arrays[scalar_name] = {}
        for channel in ("rho", "pr", "pt", "j"):
            values = finite(
                f"scalar sector {scalar_name} {channel}",
                scalar_channels[channel],
            )
            scalar_sector_summary[scalar_name][channel] = AUDIT.summarize(values, r)
            raw_scalar_arrays[scalar_name][channel] = values.tolist()

    out: dict[str, Any] = {}
    all_raw: dict[str, Any] = {}
    maximum_geom_gap = 0.0
    maximum_source_gap = 0.0

    for name, conserved, species in _species_states(state):
        prim = primitives(metrics, conserved, species)
        face_flux, source_e, source_s, q_e, q_s = _species_rhs_components(
            metric_slice,
            derivatives,
            conserved,
            species,
            dphi_t=dphi_t,
            dphi_r=dphi_r,
            beta_dm=float(coupling),
        )
        face_flux = finite(f"{name} production HLL face flux", face_flux)
        source_e = finite(f"{name} geometric energy source", source_e)
        source_s = finite(f"{name} geometric momentum source", source_s)
        q_e = finite(f"{name} interaction energy source", q_e)
        q_s = finite(f"{name} interaction momentum source", q_s)

        projection = species_projection(metrics, prim)
        projected_radial_flux_cells = finite(
            f"{name} projected alpha sqrt-gamma p_r",
            alpha * sqrt_gamma * np.asarray(projection["pr"]),
        )
        projected_radial_flux_faces = AUDIT.spherical_flux_face_values(
            projected_radial_flux_cells, dr, origin_power=2
        )
        actual_hll_radial_flux = finite(
            f"{name} normalized production HLL radial flux",
            face_flux[:, 2] / KAPPA,
        )
        hll_minus_projected = finite(
            f"{name} HLL minus reconstructed stress flux",
            actual_hll_radial_flux - projected_radial_flux_faces,
        )
        sign_flipped_hll_gap = finite(
            f"{name} sign-flipped HLL versus projected momentum flux",
            -actual_hll_radial_flux - projected_radial_flux_faces,
        )

        geometric_expected = np.empty_like(r)
        interaction_expected = np.zeros_like(r)
        q_r_profile = np.zeros_like(r)
        for i, (metric, q) in enumerate(zip(metrics, prim)):
            tensor = stress_tensor_contravariant(metric, q)
            dg = derivatives.radial
            geometric_expected[i] = 0.5 * metric.sqrt_minus_g * (
                tensor["tt"] * dg["tt"][i]
                + 2.0 * tensor["tr"] * dg["tr"][i]
                + tensor["rr"] * dg["rr"][i]
                + 2.0 * tensor["thth"] * dg["thth"][i]
            )
            if species is Species.DARK_MATTER and coupling != 0.0:
                q_r_profile[i] = float(coupling) * float(q.rho) * dphi_r[i]
                # Conditional covector identity: the coordinate source in the
                # conserved momentum equation is sqrt(-g) Q_r. The pinned
                # implementation is measured separately below; this is not
                # an assumption that the physical law is accepted.
                interaction_expected[i] = metric.sqrt_minus_g * q_r_profile[i]

        actual_geometric = finite(f"{name} normalized geometric momentum source", source_s / KAPPA)
        expected_geometric = finite(
            f"{name} independently contracted geometric momentum source",
            geometric_expected / KAPPA,
        )
        geom_gap = finite(
            f"{name} production minus independent geometric source",
            actual_geometric - expected_geometric,
        )
        actual_interaction = finite(
            f"{name} normalized production interaction momentum source", q_s / KAPPA
        )
        expected_interaction = finite(
            f"{name} conditional sqrt(-g) Q_r source", interaction_expected / KAPPA
        )
        source_gap = finite(
            f"{name} production interaction source minus conditional covariant source",
            actual_interaction - expected_interaction,
        )
        omitted_geometric_source_residual = finite(
            f"{name} residual if geometric momentum source is omitted",
            -actual_geometric,
        )
        omitted_interaction_source_residual = finite(
            f"{name} residual if interaction momentum source is omitted",
            -actual_interaction,
        )
        maximum_geom_gap = max(maximum_geom_gap, float(np.max(np.abs(geom_gap))))
        if name == "dark_matter":
            maximum_source_gap = max(maximum_source_gap, float(np.max(np.abs(source_gap))))

        out[name] = {
            "species": species.value,
            "face_flux_rms_by_channel": {
                "rest_over_kappa": summarize_faces(face_flux[:, 0] / KAPPA, dr),
                "energy_over_kappa": summarize_faces(face_flux[:, 1] / KAPPA, dr),
                "momentum_over_kappa": summarize_faces(actual_hll_radial_flux, dr),
            },
            "hll_momentum_flux_minus_projected_stress_face_flux": summarize_faces(
                hll_minus_projected, dr
            ),
            "negative_controls": {
                "sign_flipped_hll_minus_projected_momentum_flux": summarize_faces(
                    sign_flipped_hll_gap, dr
                ),
                "omitted_geometric_momentum_source_residual": AUDIT.summarize(
                    omitted_geometric_source_residual, r
                ),
                "omitted_interaction_momentum_source_residual": AUDIT.summarize(
                    omitted_interaction_source_residual, r
                ),
            },
            "geometric_energy_source_actual_over_kappa": AUDIT.summarize(
                finite(f"{name} normalized geometric energy source", source_e / KAPPA), r
            ),
            "interaction_energy_source_actual_over_kappa": AUDIT.summarize(
                finite(f"{name} normalized interaction energy source", q_e / KAPPA), r
            ),
            "geometric_momentum_source_actual_over_kappa": AUDIT.summarize(actual_geometric, r),
            "geometric_momentum_source_independent_over_kappa": AUDIT.summarize(expected_geometric, r),
            "geometric_momentum_source_gap_over_kappa": AUDIT.summarize(geom_gap, r),
            "interaction_momentum_source_actual_over_kappa": AUDIT.summarize(actual_interaction, r),
            "interaction_momentum_source_conditional_covariant_over_kappa": AUDIT.summarize(expected_interaction, r),
            "interaction_momentum_source_gap_over_kappa": AUDIT.summarize(source_gap, r),
            "interaction_source_identity_note": (
                "For dark matter this compares the pinned source with sqrt(-g) Q_r "
                "under the declared conditional covector premise. A mismatch is reported, "
                "not admitted as a production repair."
            ),
        }
        all_raw[name] = {
            "hll_face_flux_rest_over_kappa": (face_flux[:, 0] / KAPPA).tolist(),
            "hll_face_flux_energy_over_kappa": (face_flux[:, 1] / KAPPA).tolist(),
            "hll_face_flux_momentum_over_kappa": actual_hll_radial_flux.tolist(),
            "projected_stress_momentum_face_flux_over_kappa": projected_radial_flux_faces.tolist(),
            "hll_minus_projected_momentum_face_flux_over_kappa": hll_minus_projected.tolist(),
            "negative_control_sign_flipped_hll_minus_projected_momentum_face_flux_over_kappa": sign_flipped_hll_gap.tolist(),
            "negative_control_omitted_geometric_momentum_source_residual_over_kappa": omitted_geometric_source_residual.tolist(),
            "negative_control_omitted_interaction_momentum_source_residual_over_kappa": omitted_interaction_source_residual.tolist(),
            "geometric_energy_source_actual_over_kappa": (source_e / KAPPA).tolist(),
            "interaction_energy_source_actual_over_kappa": (q_e / KAPPA).tolist(),
            "geometric_momentum_source_actual_over_kappa": actual_geometric.tolist(),
            "geometric_momentum_source_independent_over_kappa": expected_geometric.tolist(),
            "geometric_momentum_source_gap_over_kappa": geom_gap.tolist(),
            "interaction_momentum_source_actual_over_kappa": actual_interaction.tolist(),
            "interaction_momentum_source_conditional_covariant_over_kappa": expected_interaction.tolist(),
            "interaction_momentum_source_gap_over_kappa": source_gap.tolist(),
            "scalar_dphi_r": dphi_r.tolist(),
        }

    return {
        "beta_dm": float(coupling),
        "zero_shift": True,
        "species": out,
        "raw_arrays": all_raw,
        "scalar_sector_projections": scalar_sector_summary,
        "raw_scalar_arrays": raw_scalar_arrays,
        "scalar_sector_note": (
            "S, D and phi scalar stress projections are recorded separately; the scalar spatial operator is not treated as a fluid HLL face flux."
        ),
        "max_abs_geometric_source_gap_over_kappa": maximum_geom_gap,
        "max_abs_dark_matter_conditional_interaction_source_gap_over_kappa": maximum_source_gap,
    }


def make_manufactured_operator_probe(kernel, n: int, r_max: float):
    """Return a copied synthetic state with nontrivial regular fluid/scalar profiles.

    This is an operator-discrimination fixture, not an evolved production state.
    The native gauge remains zero shift; the profile perturbations exist only in
    this disposable copy and are hash-checked for immutability.
    """
    state = kernel.initialize(resolution=n, r_max=r_max, include_radiation=True)
    r = finite("probe cell radii", state.grid.centers)
    state.geometry.alpha[:] = 0.75 + 0.20 * np.exp(-(r / 6.0) ** 2)
    state.scalars.phi[:] = 0.01 + 0.001 * r**2
    if np.any(np.asarray(state.geometry.beta) != 0.0):
        raise RuntimeError("manufactured operator probe must preserve native zero shift")

    metric_slice = AUDIT.metric_slice_from_q(state.grid, state.geometry)
    metrics = _metric_arrays(metric_slice)
    dm_prim = primitives(metrics, state.matter.dark_matter, Species.DARK_MATTER)
    b_prim = primitives(metrics, state.matter.baryons, Species.BARYON)
    r_prim = primitives(metrics, state.matter.radiation, Species.RADIATION)

    dm_rho = finite("probe dark-matter base rho", np.asarray([q.rho for q in dm_prim]))
    b_rho = finite("probe baryon base rho", np.asarray([q.rho for q in b_prim]))
    rad_rho = finite("probe radiation base rho", np.asarray([q.rho for q in r_prim]))
    v_dm = 0.04 * (r / 6.0) * np.exp(-(r / 6.0) ** 2)
    v_b = -0.025 * (r / 5.0) * np.exp(-(r / 5.0) ** 2)
    v_rad = 0.035 * (r / 7.0) * np.exp(-(r / 7.0) ** 2)

    state.matter.dark_matter = initialize_dust(
        metrics, dm_rho * (1.0 + 0.08 * np.exp(-(r / 5.0) ** 2)), v_dm
    )
    state.matter.baryons = initialize_dust(
        metrics, b_rho * (1.0 + 0.04 * np.exp(-(r / 8.0) ** 2)), v_b
    )
    state.matter.radiation = initialize_radiation(
        metrics, rad_rho * (1.0 + 0.06 * np.exp(-(r / 7.0) ** 2)), v_rad
    )
    return state


def summarize_probe(state, coupling: float) -> dict[str, Any]:
    n = len(state.grid.centers)
    state_hash_before = AUDIT.state_digest(state)
    zeros = np.zeros(n, dtype=np.float64)
    # Use explicit d(phi)/dt=0 for this spatially manufactured source control.
    snapshot = operator_snapshot(state, coupling, {"scalars": {"phi": zeros}})
    projection = AUDIT.independent_total_projections(state)
    state_unchanged = AUDIT.state_digest(state) == state_hash_before
    if not projection["projection_crosscheck_passed"]:
        raise RuntimeError("manufactured operator probe projection cross-check failed")
    if not state_unchanged:
        raise RuntimeError("manufactured operator probe mutated the copied state")

    max_hll_projection_gap = max(
        float(np.max(np.abs(np.asarray(piece["hll_minus_projected_momentum_face_flux_over_kappa"]))))
        for piece in snapshot["raw_arrays"].values()
    )
    dm_gap = float(np.max(np.abs(np.asarray(
        snapshot["raw_arrays"]["dark_matter"]["interaction_momentum_source_gap_over_kappa"]
    ))))
    max_sign_flip = max(
        float(np.max(np.abs(np.asarray(
            piece["negative_control_sign_flipped_hll_minus_projected_momentum_face_flux_over_kappa"]
        ))))
        for piece in snapshot["raw_arrays"].values()
    )
    max_omitted_geo = max(
        float(np.max(np.abs(np.asarray(
            piece["negative_control_omitted_geometric_momentum_source_residual_over_kappa"]
        ))))
        for piece in snapshot["raw_arrays"].values()
    )
    max_omitted_interaction = max(
        float(np.max(np.abs(np.asarray(
            piece["negative_control_omitted_interaction_momentum_source_residual_over_kappa"]
        ))))
        for piece in snapshot["raw_arrays"].values()
    )
    return {
        "synthetic_discriminator_only": True,
        "native_shift": 0.0,
        "initial_state_sha256": state_hash_before,
        "state_unchanged": state_unchanged,
        "projection_crosscheck_passed": True,
        "max_abs_hll_minus_projected_momentum_face_flux_over_kappa": max_hll_projection_gap,
        "max_abs_dark_matter_conditional_qr_source_gap_over_kappa": dm_gap,
        "max_abs_sign_flipped_hll_gap_over_kappa": max_sign_flip,
        "max_abs_omitted_geometric_source_residual_over_kappa": max_omitted_geo,
        "max_abs_omitted_interaction_source_residual_over_kappa": max_omitted_interaction,
        "nontrivial_hll_flux_discriminator": bool(max_hll_projection_gap > 1.0e-14),
        "nontrivial_conditional_qr_discriminator": bool(dm_gap > 1.0e-12) if coupling != 0.0 else True,
        "nontrivial_negative_controls": bool(
            max_sign_flip > 1.0e-14
            and max_omitted_geo > 1.0e-14
            and (coupling == 0.0 or max_omitted_interaction > 1.0e-12)
        ),
        "operator_snapshot": snapshot,
    }


def run_resolution(n: int, r_max: float, epsilons: tuple[float, ...]) -> dict[str, Any]:
    kernel = AUDIT.V55ProductionKernel()
    state = kernel.initialize(resolution=n, r_max=r_max, include_radiation=True)
    AUDIT.assert_supported_baseline_source_scope(state)
    if np.any(np.asarray(state.geometry.beta) != 0.0):
        raise RuntimeError("expected native initialized CMC shift beta=0")
    original_hash = AUDIT.state_digest(state)
    original_beta_dm = float(production_kernel_module.BETA_DM)
    arms = {}
    tangent_witnesses: dict[str, Any] = {}
    try:
        for label, coupling in (("uncoupled_control", 0.0), ("coupled_baseline", BETA_DM_COUPLED)):
            production_kernel_module.BETA_DM = float(coupling)
            base_rates = AUDIT.rhs_direction(kernel, state)
            base = operator_snapshot(state, coupling, base_rates)
            tangent_witnesses[label] = {}
            for epsilon in epsilons:
                plus = AUDIT.shifted_state(kernel, state, base_rates, +float(epsilon))
                minus = AUDIT.shifted_state(kernel, state, base_rates, -float(epsilon))
                AUDIT.assert_supported_baseline_source_scope(plus)
                AUDIT.assert_supported_baseline_source_scope(minus)
                if AUDIT.state_digest(state) != original_hash:
                    raise RuntimeError("directional shifted-state construction mutated the frozen base")
                plus_rates = AUDIT.rhs_direction(kernel, plus)
                minus_rates = AUDIT.rhs_direction(kernel, minus)
                plus_ops = operator_snapshot(plus, coupling, plus_rates)
                minus_ops = operator_snapshot(minus, coupling, minus_rates)
                species_tangent = {}
                for species_name in base["raw_arrays"]:
                    plus_arrays = plus_ops["raw_arrays"][species_name]
                    minus_arrays = minus_ops["raw_arrays"][species_name]
                    species_tangent[species_name] = {}
                    for key in (
                        "hll_face_flux_momentum_over_kappa",
                        "geometric_energy_source_actual_over_kappa",
                        "interaction_energy_source_actual_over_kappa",
                        "geometric_momentum_source_actual_over_kappa",
                        "interaction_momentum_source_actual_over_kappa",
                    ):
                        derivative = finite(
                            f"{species_name} {key} tangent",
                            (
                                np.asarray(plus_arrays[key], dtype=np.float64)
                                - np.asarray(minus_arrays[key], dtype=np.float64)
                            ) / (2.0 * float(epsilon)),
                        )
                        # Face arrays and cell source arrays have distinct
                        # radius coordinates; summarize with their own helper.
                        species_tangent[species_name][key] = (
                            summarize_faces(derivative, float(state.grid.dr))
                            if "hll_face_flux" in key
                            else AUDIT.summarize(derivative, state.grid.centers)
                        )
                scalar_tangent = {}
                for scalar_name, channels in base["raw_scalar_arrays"].items():
                    scalar_tangent[scalar_name] = {}
                    for channel in ("rho", "pr", "pt", "j"):
                        derivative = finite(
                            f"scalar {scalar_name} {channel} tangent",
                            (
                                np.asarray(plus_ops["raw_scalar_arrays"][scalar_name][channel], dtype=np.float64)
                                - np.asarray(minus_ops["raw_scalar_arrays"][scalar_name][channel], dtype=np.float64)
                            ) / (2.0 * float(epsilon)),
                        )
                        scalar_tangent[scalar_name][channel] = AUDIT.summarize(
                            derivative, state.grid.centers
                        )
                tangent_witnesses[label][f"epsilon_{epsilon:.8g}"] = {
                    "epsilon": float(epsilon),
                    "plus_state_sha256": AUDIT.state_digest(plus),
                    "minus_state_sha256": AUDIT.state_digest(minus),
                    "species_operator_tangent_summaries": species_tangent,
                    "scalar_sector_projection_tangent_summaries": scalar_tangent,
                }
            arms[label] = base
    finally:
        production_kernel_module.BETA_DM = original_beta_dm

    state_unchanged = AUDIT.state_digest(state) == original_hash
    if not state_unchanged:
        raise RuntimeError("production HLL diagnostic mutated its frozen base state")
    projection = AUDIT.independent_total_projections(state)
    if not projection["projection_crosscheck_passed"]:
        raise RuntimeError("independent total stress projection disagrees with canonical assembly")
    all_values_finite = True
    for arm in arms.values():
        for species_arrays in arm["raw_arrays"].values():
            for values in species_arrays.values():
                all_values_finite &= bool(np.all(np.isfinite(np.asarray(values, dtype=float))))
        for scalar_arrays in arm["raw_scalar_arrays"].values():
            for values in scalar_arrays.values():
                all_values_finite &= bool(np.all(np.isfinite(np.asarray(values, dtype=float))))
    # A separate disposable manufactured copy excites nonzero HLL currents and
    # a nonzero radial scalar gradient, avoiding a vacuous zero-profile baseline.
    probe_state = make_manufactured_operator_probe(kernel, n, r_max)
    probe_hash_before = AUDIT.state_digest(probe_state)
    probe_arms = {
        "uncoupled_control": summarize_probe(probe_state, 0.0),
        "coupled_conditional_source": summarize_probe(probe_state, BETA_DM_COUPLED),
    }
    probe_unchanged = AUDIT.state_digest(probe_state) == probe_hash_before
    if not probe_unchanged:
        raise RuntimeError("manufactured probe changed between source arms")
    probe_hll_nontrivial = all(
        arm["nontrivial_hll_flux_discriminator"] for arm in probe_arms.values()
    )
    probe_q_nontrivial = probe_arms["coupled_conditional_source"]["nontrivial_conditional_qr_discriminator"]
    if not probe_hll_nontrivial or not probe_q_nontrivial:
        raise RuntimeError(
            "manufactured operator control was vacuous: "
            f"hll={probe_hll_nontrivial}, conditional_qr={probe_q_nontrivial}"
        )

    return {
        "resolution": int(n),
        "r_max": float(r_max),
        "dr": float(state.grid.dr),
        "initial_state_sha256": original_hash,
        "base_state_unchanged": state_unchanged,
        "zero_shift_native_slice": True,
        "projection_crosscheck_passed": True,
        "all_hll_source_profiles_finite": bool(all_values_finite),
        "production_steps_advanced": 0,
        "arms": arms,
        "directional_epsilon_witnesses": tangent_witnesses,
        "manufactured_operator_probe": {
            "state_unchanged": probe_unchanged,
            "synthetic_discriminator_only": True,
            "profile": {
                "alpha": "0.75 + 0.20 exp(-(r/6)^2)",
                "phi": "0.01 + 0.001 r^2",
                "dm_rho": "base * (1 + 0.08 exp(-(r/5)^2))",
                "dm_v_r": "0.04 (r/6) exp(-(r/6)^2)",
                "baryon_rho": "base * (1 + 0.04 exp(-(r/8)^2)",
                "baryon_v_r": "-0.025 (r/5) exp(-(r/5)^2)",
                "radiation_rho": "base * (1 + 0.06 exp(-(r/7)^2))",
                "radiation_v_r": "0.035 (r/7) exp(-(r/7)^2)",
            },
            "arms": probe_arms,
        },
    }


def run(resolutions: tuple[int, ...], r_max: float, epsilons: tuple[float, ...]) -> dict[str, Any]:
    provenance = AUDIT.verify_provenance()
    rows = [run_resolution(n, r_max, epsilons) for n in resolutions]
    all_unchanged = all(row["base_state_unchanged"] for row in rows)
    all_finite = all(row["all_hll_source_profiles_finite"] for row in rows)
    probe_passed = all(
        row["manufactured_operator_probe"]["state_unchanged"]
        and all(
            arm["nontrivial_hll_flux_discriminator"]
            and arm["nontrivial_conditional_qr_discriminator"]
            and arm["nontrivial_negative_controls"]
            for arm in row["manufactured_operator_probe"]["arms"].values()
        )
        for row in rows
    )
    max_geom = max(
        arm["max_abs_geometric_source_gap_over_kappa"]
        for row in rows for arm in row["arms"].values()
    )
    max_q_gap = max(
        arm["max_abs_dark_matter_conditional_interaction_source_gap_over_kappa"]
        for row in rows for arm in row["arms"].values()
    )
    return {
        "status": "completed_diagnostic_only",
        "title": "Pinned production HLL face-flux and source operator pairing",
        "source_provenance": provenance,
        "configuration": {
            "production_source_commit": SOURCE_COMMIT,
            "vendor_commit": VENDOR_COMMIT,
            "resolutions": list(resolutions),
            "r_max": float(r_max),
            "epsilons": list(epsilons),
            "native_shift": 0.0,
            "couplings": {"uncoupled_control": 0.0, "coupled_baseline": BETA_DM_COUPLED},
        },
        "gate_result": {
            "source_pins_verified": True,
            "independent_total_projection_matches_canonical": True,
            "all_hll_source_profiles_finite": all_finite,
            "input_states_unchanged": all_unchanged,
            "native_zero_shift_scope_verified": True,
            "manufactured_operator_probe_nontrivial": probe_passed,
            "manufactured_probe_negative_controls_nontrivial": probe_passed,
            "scalar_sector_kept_separate": True,
            "production_steps_advanced": 0,
            "max_abs_independent_geometric_source_gap_over_kappa": max_geom,
            "max_abs_dark_matter_conditional_interaction_source_gap_over_kappa": max_q_gap,
            "conditional_exchange_law_accepted": False,
            "independent_total_stress_energy_divergence_closed": False,
            "scientific_admission": False,
        },
        "interpretation": {
            "geometric_source": "Compare the pinned production momentum geometric source with an independently computed 1/2 sqrt(-g) T^{ab} partial_r g_ab contraction.",
            "hll_flux": "Record actual pinned per-species HLL face fluxes and compare them with face reconstruction of cell-projected radial stress. Differences characterize the operator pair; no exact face equality is presumed for the dissipative HLL scheme.",
            "interaction_source": "Compare actual pinned dark-matter q_s with sqrt(-g) Q_r only under the explicitly conditional covector premise. Preserve any mismatch; do not change production source terms in this diagnostic.",
            "limits": [
                "This is a frozen-state RHS-component readout; no finite production timestep is advanced.",
                "The scalar-field spatial operator remains separate from the fluid HLL operator.",
                "It does not close the total stress-energy divergence or accept the physical exchange law.",
            ],
        },
        "resolution_results": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[60, 80, 120, 160])
    parser.add_argument("--r-max", type=float, default=24.0)
    parser.add_argument("--epsilons", nargs="+", type=float, default=list(DEFAULT_EPSILONS))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload: dict[str, Any]
    code = 0
    try:
        if len(args.epsilons) != 2 or any(not np.isfinite(x) or x <= 0.0 for x in args.epsilons):
            raise ValueError("exactly two finite positive directional epsilons are required")
        if len(args.resolutions) < 2 or any(n < 16 for n in args.resolutions):
            raise ValueError("at least two resolutions with n>=16 are required")
        payload = run(tuple(args.resolutions), args.r_max, tuple(args.epsilons))
        payload["error"] = None
    except Exception as exc:
        code = 1
        payload = {
            "status": "failed_closed",
            "scientific_admission": False,
            "gate_result": {
                "source_pins_verified": False,
                "all_hll_source_profiles_finite": False,
                "input_states_unchanged": False,
                "production_steps_advanced": 0,
                "independent_total_stress_energy_divergence_closed": False,
                "scientific_admission": False,
            },
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
        }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload.get("status"),
        "source_pins_verified": payload.get("gate_result", {}).get("source_pins_verified"),
        "input_states_unchanged": payload.get("gate_result", {}).get("input_states_unchanged"),
        "production_steps_advanced": payload.get("gate_result", {}).get("production_steps_advanced"),
        "independent_total_stress_energy_divergence_closed": False,
        "scientific_admission": False,
        "error": payload.get("error"),
    }, sort_keys=True))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
