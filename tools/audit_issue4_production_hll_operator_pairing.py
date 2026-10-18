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


def summarize_signed_budget(values: Any, radii: Any) -> dict[str, Any]:
    """Keep signed first-cell, second-cell, and fixed physical-annulus witnesses."""
    x = finite("signed discrete budget", values)
    r = finite("signed discrete budget cell radii", radii)
    if x.shape != r.shape or x.ndim != 1 or x.size < 2:
        raise ValueError("signed budget and cell-centre radii must be matching 1-D arrays")

    mask = (r >= 2.0) & (r <= 10.0)
    def cell(i: int) -> dict[str, Any]:
        return {"index": int(i), "radius": float(r[i]), "signed_value": float(x[i])}
    if np.any(mask):
        indices = np.flatnonzero(mask)
        local = int(indices[int(np.argmax(np.abs(x[mask])))])
        annulus = {
            "cell_count": int(indices.size),
            "r_min": float(np.min(r[mask])),
            "r_max": float(np.max(r[mask])),
            "signed_rms": float(np.sqrt(np.mean(x[mask] ** 2))),
            "max_abs": float(abs(x[local])),
            "signed_at_max_abs": float(x[local]),
            "index_at_max_abs": local,
            "radius_at_max_abs": float(r[local]),
        }
    else:
        annulus = {"cell_count": 0, "signed_rms": None, "max_abs": None}
    i_max = int(np.argmax(np.abs(x)))
    return {
        "cell_0": cell(0),
        "cell_1": cell(1),
        "annulus_2_to_10": annulus,
        "global": {
            "signed_rms": float(np.sqrt(np.mean(x * x))),
            "max_abs": float(abs(x[i_max])),
            "signed_at_max_abs": float(x[i_max]),
            "index_at_max_abs": i_max,
            "radius_at_max_abs": float(r[i_max]),
        },
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


def assert_supported_operator_scope(state) -> None:
    """Fail closed outside the declared DM/baryon/radiation, zero-shift scope.

    The pinned production kernel adds the optional Maxwell Lorentz source
    outside _species_rhs_components, while the paired total-T projection does
    not include optional EM stress by default. Capturing the fluid helper alone
    is therefore insufficient for either charged-dust or EM-enabled states.
    """
    matter = getattr(state, "matter", None)
    charged_dust = getattr(matter, "charged_dust", None)
    electromagnetic = getattr(state, "electromagnetic", None)
    if charged_dust is not None or electromagnetic is not None:
        active = []
        if charged_dust is not None:
            active.append("charged_dust")
        if electromagnetic is not None:
            active.append("radial Maxwell")
        raise ValueError(
            "optional charged-dust/Maxwell state is outside this observer's "
            "verified scope (" + ", ".join(active) + "); refusing partial RHS "
            "or total-stress accounting. A separate EM-on observer must capture "
            "the Maxwell stress projection and Lorentz source."
        )


def operator_snapshot(state, coupling: float, rates: dict[str, Any]) -> dict[str, Any]:
    """Record pinned HLL fluxes and sources plus independent radial contractions."""
    assert_supported_operator_scope(state)
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
        face_flux, source_e, source_s, q_e, q_s, q_rest = _species_rhs_components(
            metric_slice,
            derivatives,
            conserved,
            species,
            dphi_t=dphi_t,
            dphi_r=dphi_r,
            beta_dm=float(coupling),
        )
        face_flux = finite(f"{name} production HLL face flux", face_flux)
        if face_flux.shape != (r.size + 1, 3):
            raise ValueError(
                f"{name} face flux shape {face_flux.shape} != {(r.size + 1, 3)}"
            )
        origin_face_flux_exact_zero = bool(
            np.array_equal(face_flux[0], np.zeros(3, dtype=np.float64))
        )
        if not origin_face_flux_exact_zero:
            raise RuntimeError(
                f"{name} production HLL centre-face flux is not exactly zero"
            )
        source_e = finite(f"{name} geometric energy source", source_e)
        source_s = finite(f"{name} geometric momentum source", source_s)
        q_e = finite(f"{name} interaction energy source", q_e)
        q_s = finite(f"{name} interaction momentum source", q_s)
        for source_name, source_values in (
            ("source_e", source_e), ("source_s", source_s),
            ("q_e", q_e), ("q_s", q_s),
        ):
            if source_values.shape != r.shape:
                raise ValueError(
                    f"{name} {source_name} shape {source_values.shape} != {r.shape}"
                )

        # Exact semi-discrete convention used by pinned evolve_species: negative
        # face-flux difference over dr, then geometric and interaction sources.
        # Apply KAPPA normalization consistently to every term in each channel.
        rhs_rest_over_kappa = finite(
            f"{name} combined rest RHS over KAPPA",
            (-np.diff(face_flux[:, 0]) / dr) / KAPPA,
        )
        rhs_energy_over_kappa = finite(
            f"{name} combined energy RHS over KAPPA",
            (-np.diff(face_flux[:, 1]) / dr + source_e + q_e) / KAPPA,
        )
        rhs_momentum_over_kappa = finite(
            f"{name} combined radial momentum RHS over KAPPA",
            (-np.diff(face_flux[:, 2]) / dr + source_s + q_s) / KAPPA,
        )
        # Controls perturb one term at a time in this same signed combined RHS.
        rhs_momentum_sign_flipped_flux_over_kappa = finite(
            f"{name} sign-flipped-flux momentum RHS control",
            (np.diff(face_flux[:, 2]) / dr + source_s + q_s) / KAPPA,
        )
        rhs_momentum_omitted_geometric_source_over_kappa = finite(
            f"{name} omitted-geometric-source momentum RHS control",
            (-np.diff(face_flux[:, 2]) / dr + q_s) / KAPPA,
        )
        rhs_momentum_omitted_interaction_source_over_kappa = finite(
            f"{name} omitted-interaction-source momentum RHS control",
            (-np.diff(face_flux[:, 2]) / dr + source_s) / KAPPA,
        )

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
        interaction_lapse_omitted_expected = np.zeros_like(r)
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
                # Explicit negative control: omit the lapse from the coordinate
                # source. This is not the covariant target; it must be distinguishable
                # from the implemented alpha*sqrt(gamma)*Q_r source.
                interaction_lapse_omitted_expected[i] = metric.sqrt_gamma * q_r_profile[i]

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
        expected_lapse_omitted_interaction = finite(
            f"{name} lapse-omitted conditional interaction source",
            interaction_lapse_omitted_expected / KAPPA,
        )
        lapse_omission_gap = finite(
            f"{name} production interaction source minus lapse-omitted control",
            actual_interaction - expected_lapse_omitted_interaction,
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
            # The gate requires a non-vacuous contrast against an intentionally
            # lapse-omitted negative control, not a discrepancy from the covariant
            # source that the current implementation is expected to reproduce.
            maximum_source_gap = max(
                maximum_source_gap, float(np.max(np.abs(lapse_omission_gap)))
            )

        out[name] = {
            "species": species.value,
            "origin_face_flux_exact_zero": origin_face_flux_exact_zero,
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
            "interaction_momentum_source_lapse_omitted_control_over_kappa": AUDIT.summarize(expected_lapse_omitted_interaction, r),
            "interaction_momentum_source_lapse_omitted_control_gap_over_kappa": AUDIT.summarize(lapse_omission_gap, r),
            "interaction_source_identity_note": (
                "For dark matter this compares the pinned source with sqrt(-g) Q_r "
                "under the declared conditional covector premise. A mismatch is reported, "
                "not admitted as a production repair."
            ),
            "production_rhs_budget": {
                "discrete_convention": "rhs_i = -(F[i+1]-F[i])/dr + geometric_source_i + interaction_source_i; rest has no source term; all channels divided by KAPPA=8*pi.",
                "spatial_measure_note": "This is the native evolve_species semi-discrete update convention, not an extra cell-volume reweighting and not a total-T closure claim.",
                "rest": summarize_signed_budget(rhs_rest_over_kappa, r),
                "energy": summarize_signed_budget(rhs_energy_over_kappa, r),
                "radial_momentum": summarize_signed_budget(rhs_momentum_over_kappa, r),
                "negative_controls": {
                    "radial_momentum_sign_flipped_flux": summarize_signed_budget(rhs_momentum_sign_flipped_flux_over_kappa, r),
                    "radial_momentum_omitted_geometric_source": summarize_signed_budget(rhs_momentum_omitted_geometric_source_over_kappa, r),
                    "radial_momentum_omitted_interaction_source": summarize_signed_budget(rhs_momentum_omitted_interaction_source_over_kappa, r),
                },
            },
        }
        all_raw[name] = {
            # Preserve the exact five native _species_rhs_components outputs
            # before scaling, alongside the normalized comparison arrays.
            "native_hll_face_flux": face_flux.tolist(),
            "native_geometric_energy_source": source_e.tolist(),
            "native_geometric_momentum_source": source_s.tolist(),
            "native_interaction_energy_source": q_e.tolist(),
            "native_interaction_momentum_source": q_s.tolist(),
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
            "interaction_momentum_source_lapse_omitted_control_over_kappa": expected_lapse_omitted_interaction.tolist(),
            "interaction_momentum_source_lapse_omitted_control_gap_over_kappa": lapse_omission_gap.tolist(),
            "production_rhs_rest_over_kappa": rhs_rest_over_kappa.tolist(),
            "production_rhs_energy_over_kappa": rhs_energy_over_kappa.tolist(),
            "production_rhs_momentum_over_kappa": rhs_momentum_over_kappa.tolist(),
            "production_rhs_momentum_sign_flipped_flux_control_over_kappa": rhs_momentum_sign_flipped_flux_over_kappa.tolist(),
            "production_rhs_momentum_omitted_geometric_source_control_over_kappa": rhs_momentum_omitted_geometric_source_over_kappa.tolist(),
            "production_rhs_momentum_omitted_interaction_source_control_over_kappa": rhs_momentum_omitted_interaction_source_over_kappa.tolist(),
            "scalar_dphi_r": dphi_r.tolist(),
        }

    required_native_arrays = (
        "native_hll_face_flux",
        "native_geometric_energy_source",
        "native_geometric_momentum_source",
        "native_interaction_energy_source",
        "native_interaction_momentum_source",
    )
    return {
        "beta_dm": float(coupling),
        "zero_shift": True,
        "origin_face_flux_exact_zero": all(
            item["origin_face_flux_exact_zero"] for item in out.values()
        ),
        "native_component_arrays_preserved": all(
            all(key in item for key in required_native_arrays)
            for item in all_raw.values()
        ),
        "species": out,
        "raw_arrays": all_raw,
        "scalar_sector_projections": scalar_sector_summary,
        "raw_scalar_arrays": raw_scalar_arrays,
        "scalar_sector_note": (
            "S, D and phi scalar stress projections are recorded separately; the scalar spatial operator is not treated as a fluid HLL face flux."
        ),
        "max_abs_geometric_source_gap_over_kappa": maximum_geom_gap,
        "max_abs_dark_matter_lapse_omission_control_gap_over_kappa": maximum_source_gap,
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
        snapshot["raw_arrays"]["dark_matter"]["interaction_momentum_source_lapse_omitted_control_gap_over_kappa"]
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


def matched_total_t_witness(
    kernel,
    state,
    coupling: float,
    epsilons: tuple[float, ...],
    expected_state_hash: str,
) -> dict[str, Any]:
    """Run the existing total-T directional/control audit on the HLL observer's exact state.

    Separate workflow artifacts can share source pins and grid settings while
    still recording different bitwise initial-state hashes. This wrapper makes
    the relationship explicit by running both observers against this same
    in-memory state and fail-closing if the digest or epsilon contract differs.
    """
    before_hash = AUDIT.state_digest(state)
    if before_hash != expected_state_hash:
        raise RuntimeError("shared-state witness received a state different from the HLL base")

    witness = AUDIT.audit_case(kernel, state, float(coupling), epsilons)
    after_hash = AUDIT.state_digest(state)
    if after_hash != before_hash:
        raise RuntimeError("matched total-T witness mutated the shared frozen base state")
    if witness.get("initial_state_sha256") != before_hash:
        raise RuntimeError("total-T witness did not retain the shared HLL base-state digest")

    epsilon_rows = witness.get("epsilon_rows", [])
    if len(epsilon_rows) != len(epsilons) or any(
        not np.isclose(float(row.get("epsilon", float("nan"))), float(eps),
                       rtol=0.0, atol=0.0)
        for row, eps in zip(epsilon_rows, epsilons)
    ):
        raise RuntimeError("matched total-T witness epsilon rows do not match the HLL arms")

    raw_arrays = witness.get("raw_arrays", {})
    raw_arrays_finite = bool(raw_arrays) and all(
        np.all(np.isfinite(np.asarray(values, dtype=np.float64)))
        for values in raw_arrays.values()
    )
    if not witness.get("all_profiles_finite", False) or not raw_arrays_finite:
        raise FloatingPointError("matched total-T witness contains a non-finite profile")

    witness["matched_state_contract"] = {
        "same_in_memory_base_state_as_hll_operator_snapshot": True,
        "initial_state_sha256": before_hash,
        "state_digest_unchanged_after_total_t_audit": True,
        "matched_epsilon_values": [float(value) for value in epsilons],
        "all_total_t_directional_arrays_finite": raw_arrays_finite,
        "production_steps_advanced": 0,
        "conditional_exchange_law_accepted": False,
        "independent_total_stress_energy_divergence_closed": False,
        "scientific_admission": False,
    }
    return witness


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
                plus_state_sha256 = AUDIT.state_digest(plus)
                minus_state_sha256 = AUDIT.state_digest(minus)
                plus_rates = AUDIT.rhs_direction(kernel, plus)
                minus_rates = AUDIT.rhs_direction(kernel, minus)
                plus_ops = operator_snapshot(plus, coupling, plus_rates)
                minus_ops = operator_snapshot(minus, coupling, minus_rates)

                # Compute the existing covariant total-T witness directly from
                # these exact perturbed objects, not independently reconstructed
                # copies. That makes the species and total-T arrays a true paired
                # observation at this base state, coupling, and epsilon.
                paired_total_t = AUDIT.sector_resolved_divergence_tangent_from_states(
                    state, plus, minus, float(epsilon)
                )
                total_t = paired_total_t["total"]
                same_state_pair_verified = bool(
                    AUDIT.state_digest(state) == original_hash
                    and AUDIT.state_digest(plus) == plus_state_sha256
                    and AUDIT.state_digest(minus) == minus_state_sha256
                )
                if not same_state_pair_verified:
                    raise RuntimeError(
                        "shared HLL/total-T directional state pair changed during read-only observation"
                    )
                total_t_arrays = {
                    key: finite(f"paired total-T {key}", total_t[key])
                    for key in ("div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv")
                }
                if any(
                    values.shape != state.grid.centers.shape
                    for values in total_t_arrays.values()
                ):
                    raise ValueError("paired total-T profiles do not match the frozen cell-centre grid")
                same_state_total_t = {
                    "state_pairing_verified": same_state_pair_verified,
                    "pairing_method": (
                        "same in-memory plus/minus states passed to the HLL/source observer "
                        "and the existing total-T directional witness"
                    ),
                    "base_state_sha256": original_hash,
                    "plus_state_sha256": plus_state_sha256,
                    "minus_state_sha256": minus_state_sha256,
                    "epsilon": float(epsilon),
                    "raw_signed_arrays": {
                        key: values.tolist() for key, values in total_t_arrays.items()
                    },
                    "div_T_t_summary": total_t["div_T_t_summary"],
                    "div_T_r_summary": total_t["div_T_r_summary"],
                    "div_T_t_fv_summary": total_t["div_T_t_fv_summary"],
                    "div_T_r_fv_summary": total_t["div_T_r_fv_summary"],
                    "fv_face_flux_method": total_t["fv_face_flux_method"],
                    "projection_crosscheck_passed": bool(
                        total_t["projection_crosscheck_passed"]
                    ),
                    "sector_decomposition": paired_total_t["sector_decomposition"],
                }
                if not same_state_total_t["projection_crosscheck_passed"]:
                    raise RuntimeError("paired total-T projection cross-check failed")
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
                        "production_rhs_rest_over_kappa",
                        "production_rhs_energy_over_kappa",
                        "production_rhs_momentum_over_kappa",
                    ):
                        derivative = finite(
                            f"{species_name} {key} tangent",
                            (
                                np.asarray(plus_arrays[key], dtype=np.float64)
                                - np.asarray(minus_arrays[key], dtype=np.float64)
                            ) / (2.0 * float(epsilon)),
                        )
                        # Face arrays live at faces; combined RHS arrays live
                        # at cell centres. Keep 0, 1, and the fixed annulus in
                        # each of the two directional epsilon tangent arms.
                        if "hll_face_flux" in key:
                            species_tangent[species_name][key] = summarize_faces(
                                derivative, float(state.grid.dr)
                            )
                        elif key.startswith("production_rhs_"):
                            species_tangent[species_name][key] = summarize_signed_budget(
                                derivative, state.grid.centers
                            )
                        else:
                            species_tangent[species_name][key] = AUDIT.summarize(
                                derivative, state.grid.centers
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
                    "plus_state_sha256": plus_state_sha256,
                    "minus_state_sha256": minus_state_sha256,
                    "same_state_total_T_directional_witness": same_state_total_t,
                    "species_operator_tangent_summaries": species_tangent,
                    "scalar_sector_projection_tangent_summaries": scalar_tangent,
                    # Keep raw plus/minus source/flux arrays as well as their
                    # tangent summaries so directional evidence is reproducible.
                    "plus_species_operator_arrays": plus_ops["raw_arrays"],
                    "minus_species_operator_arrays": minus_ops["raw_arrays"],
                    "plus_scalar_sector_arrays": plus_ops["raw_scalar_arrays"],
                    "minus_scalar_sector_arrays": minus_ops["raw_scalar_arrays"],
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

    # Pair the total-T directional witness with the exact same in-memory state
    # used for the native HLL/source arrays above. This reuses the reviewed
    # audit_case and retains its conditional source controls.
    matched_total_t_witnesses = {
        "uncoupled_control": matched_total_t_witness(
            kernel, state, 0.0, epsilons, original_hash
        ),
        "coupled_baseline": matched_total_t_witness(
            kernel, state, BETA_DM_COUPLED, epsilons, original_hash
        ),
    }
    matched_total_t_state_hashes_match = all(
        item["matched_state_contract"]["initial_state_sha256"] == original_hash
        and item["initial_state_sha256"] == original_hash
        for item in matched_total_t_witnesses.values()
    )
    if not matched_total_t_state_hashes_match:
        raise RuntimeError("total-T and HLL/source observers did not retain one shared base-state hash")

    return {
        "resolution": int(n),
        "r_max": float(r_max),
        "dr": float(state.grid.dr),
        "initial_state_sha256": original_hash,
        "base_state_unchanged": state_unchanged,
        "zero_shift_native_slice": True,
        "origin_face_flux_exact_zero": all(
            arm["origin_face_flux_exact_zero"] for arm in arms.values()
        ),
        "native_component_arrays_preserved": all(
            arm["native_component_arrays_preserved"] for arm in arms.values()
        ),
        "directional_operator_arrays_preserved": all(
            "plus_species_operator_arrays" in witness
            and "minus_species_operator_arrays" in witness
            and "plus_scalar_sector_arrays" in witness
            and "minus_scalar_sector_arrays" in witness
            for arm_witnesses in tangent_witnesses.values()
            for witness in arm_witnesses.values()
        ),
        "projection_crosscheck_passed": True,
        "all_hll_source_profiles_finite": bool(all_values_finite),
        "production_steps_advanced": 0,
        "arms": arms,
        "directional_epsilon_witnesses": tangent_witnesses,
        "matched_total_t_witnesses": matched_total_t_witnesses,
        "matched_total_t_state_hashes_match": matched_total_t_state_hashes_match,
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
    budget_keys = (
        "production_rhs_rest_over_kappa",
        "production_rhs_energy_over_kappa",
        "production_rhs_momentum_over_kappa",
    )
    combined_budget_reported = all(
        all(key in species_arrays for key in budget_keys)
        and all(np.all(np.isfinite(np.asarray(species_arrays[key], dtype=np.float64))) for key in budget_keys)
        for row in rows
        for arm in row["arms"].values()
        for species_arrays in arm["raw_arrays"].values()
    )
    combined_budget_tangents_reported = all(
        all(key in species_tangent for key in budget_keys)
        for row in rows
        for epsilon_arms in row["directional_epsilon_witnesses"].values()
        for witness in epsilon_arms.values()
        for species_tangent in witness["species_operator_tangent_summaries"].values()
    )
    if not combined_budget_reported or not combined_budget_tangents_reported:
        raise RuntimeError("combined signed RHS budget or its epsilon tangents are missing")
    matched_total_t_witness_reported = all(
        set(row.get("matched_total_t_witnesses", {})) == {
            "uncoupled_control", "coupled_baseline"
        }
        and all(
            item.get("matched_state_contract", {}).get(
                "same_in_memory_base_state_as_hll_operator_snapshot"
            ) is True
            for item in row["matched_total_t_witnesses"].values()
        )
        for row in rows
    )
    matched_total_t_state_hashes_match = all(
        row.get("matched_total_t_state_hashes_match") is True
        for row in rows
    )
    matched_total_t_profiles_finite = all(
        item.get("matched_state_contract", {}).get(
            "all_total_t_directional_arrays_finite"
        ) is True
        for row in rows
        for item in row.get("matched_total_t_witnesses", {}).values()
    )
    same_state_epsilon_witness_reported = all(
        witness.get("same_state_total_T_directional_witness", {}).get(
            "state_pairing_verified"
        ) is True
        and witness["same_state_total_T_directional_witness"].get(
            "base_state_sha256"
        ) == row.get("initial_state_sha256")
        and witness["same_state_total_T_directional_witness"].get(
            "plus_state_sha256"
        ) == witness.get("plus_state_sha256")
        and witness["same_state_total_T_directional_witness"].get(
            "minus_state_sha256"
        ) == witness.get("minus_state_sha256")
        and witness["same_state_total_T_directional_witness"].get(
            "epsilon"
        ) == witness.get("epsilon")
        and witness["same_state_total_T_directional_witness"].get(
            "projection_crosscheck_passed"
        ) is True
        and all(
            len(witness["same_state_total_T_directional_witness"].get(
                "raw_signed_arrays", {}
            ).get(key, [])) == row["resolution"]
            and np.all(np.isfinite(np.asarray(
                witness["same_state_total_T_directional_witness"]["raw_signed_arrays"][key],
                dtype=np.float64,
            )))
            for key in ("div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv")
        )
        for row in rows
        for epsilon_arms in row["directional_epsilon_witnesses"].values()
        for witness in epsilon_arms.values()
    )
    if not (
        matched_total_t_witness_reported
        and matched_total_t_state_hashes_match
        and matched_total_t_profiles_finite
        and same_state_epsilon_witness_reported
    ):
        raise RuntimeError(
            "base-state total-T audit or exact plus/minus state-pair witness gate failed"
        )
    expected_sector_names = ["dark_matter", "baryons", "radiation", "S", "D", "phi"]
    sector_resolved_divergence_reported = all(
        witness.get("same_state_total_T_directional_witness", {}).get(
            "sector_decomposition", {}
        ).get("sector_names") == expected_sector_names
        and set(witness["same_state_total_T_directional_witness"][
            "sector_decomposition"
        ].get("sectors", {})) == set(expected_sector_names)
        and witness["same_state_total_T_directional_witness"][
            "sector_decomposition"
        ].get("state_pairing", {}).get("all_sector_arrays_finite") is True
        and all(
            set(sector.get("raw_signed_arrays", {})) == {
                "div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv"
            }
            and all(
                len(values) == row["resolution"]
                and np.all(np.isfinite(np.asarray(values, dtype=np.float64)))
                for values in sector["raw_signed_arrays"].values()
            )
            and all(
                "cell_0" in summary and "cell_1" in summary
                and "annulus_2_to_10_rms" in summary
                for summary in sector.get("summaries", {}).values()
            )
            for sector in witness["same_state_total_T_directional_witness"][
                "sector_decomposition"
            ]["sectors"].values()
        )
        for row in rows
        for epsilon_arms in row["directional_epsilon_witnesses"].values()
        for witness in epsilon_arms.values()
    )
    sector_raw_additive_reconstruction_passed = all(
        witness.get("same_state_total_T_directional_witness", {}).get(
            "sector_decomposition", {}
        ).get("reconstruction", {}).get("all_sectors_reconstruct_total") is True
        for row in rows
        for epsilon_arms in row["directional_epsilon_witnesses"].values()
        for witness in epsilon_arms.values()
    )
    sector_time_density_component_accounting_passed = all(
        witness.get("same_state_total_T_directional_witness", {}).get(
            "sector_decomposition", {}
        ).get("reconstruction", {}).get("time_density_component_accounting_passed") is True
        and all(
            channel.get("component_corrected_passed") is True
            for channel in witness["same_state_total_T_directional_witness"][
                "sector_decomposition"
            ]["reconstruction"]["channels"].values()
        )
        for row in rows
        for epsilon_arms in row["directional_epsilon_witnesses"].values()
        for witness in epsilon_arms.values()
    )
    sector_additive_reconstruction_passed = sector_time_density_component_accounting_passed
    max_geom = max(
        arm["max_abs_geometric_source_gap_over_kappa"]
        for row in rows for arm in row["arms"].values()
    )
    max_q_gap = max(
        arm["max_abs_dark_matter_lapse_omission_control_gap_over_kappa"]
        for row in rows for arm in row["arms"].values()
    )
    return {
        "status": "completed_diagnostic_only",
        "title": "Pinned production HLL face-flux and source operator pairing",
        "scope": {
            "supported_fluid_species": ["dark_matter", "baryons", "radiation"],
            "native_shift": 0.0,
            "optional_charged_dust_supported": False,
            "optional_radial_maxwell_supported": False,
            "unsupported_optional_lanes_fail_closed": True,
            "scalar_stress_sectors": ["S", "D", "phi"],
            "scalar_spatial_operator_budget_completed": False,
            "same_process_total_T_directional_witness_available": matched_total_t_witness_reported,
            "shared_initial_state_hash_verified": matched_total_t_state_hashes_match,
            "same_state_epsilon_total_T_witness_available": same_state_epsilon_witness_reported,
            "sector_resolved_total_T_decomposition_available": sector_resolved_divergence_reported,
            "sector_additive_reconstruction_passed": sector_additive_reconstruction_passed,
            "sector_raw_sum_strictly_matches_total": sector_raw_additive_reconstruction_passed,
            "sector_time_density_component_accounting_passed": sector_time_density_component_accounting_passed,
            "full_cross_sector_total_T_reconciliation_completed": False,
            "interpretation_boundary": (
                "Per-species HLL/source arrays and combined RHS budgets are "
                "constituent-level evolution accounting. S/D/phi stress projections "
                "remain separate. This observer does not assemble a complete scalar "
                "stress-divergence/source budget or a closed matter+scalar total-T "
                "identity; total divergence/source accounting remains open."
            ),
        },
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
            "origin_face_flux_exact_zero": all(row["origin_face_flux_exact_zero"] for row in rows),
            "native_component_arrays_preserved": all(row["native_component_arrays_preserved"] for row in rows),
            "directional_operator_arrays_preserved": all(row["directional_operator_arrays_preserved"] for row in rows),
            "manufactured_operator_probe_nontrivial": probe_passed,
            "manufactured_probe_negative_controls_nontrivial": probe_passed,
            "combined_discrete_budget_reported": combined_budget_reported,
            "combined_budget_tangents_reported": combined_budget_tangents_reported,
            "matched_total_t_witness_reported": matched_total_t_witness_reported,
            "matched_total_t_state_hashes_match": matched_total_t_state_hashes_match,
            "matched_total_t_profiles_finite": matched_total_t_profiles_finite,
            "same_state_epsilon_witness_reported": same_state_epsilon_witness_reported,
            "sector_resolved_divergence_reported": sector_resolved_divergence_reported,
            "sector_additive_reconstruction_passed": sector_additive_reconstruction_passed,
            "sector_raw_additive_reconstruction_passed": sector_raw_additive_reconstruction_passed,
            "sector_time_density_component_accounting_passed": sector_time_density_component_accounting_passed,
            "scalar_sector_kept_separate": True,
            "unsupported_optional_lanes_fail_closed": True,
            "optional_charged_dust_em_supported": False,
            "scalar_spatial_operator_budget_completed": False,
            "full_cross_sector_total_T_reconciliation_completed": False,
            "production_steps_advanced": 0,
            "max_abs_independent_geometric_source_gap_over_kappa": max_geom,
            "max_abs_dark_matter_lapse_omission_control_gap_over_kappa": max_q_gap,
            "conditional_exchange_law_accepted": False,
            "independent_total_stress_energy_divergence_closed": False,
            "scientific_admission": False,
        },
        "interpretation": {
            "combined_budget": "Reports the native per-species semi-discrete RHS with pinned face-flux difference and source signs, divided by KAPPA; retains signed cells 0 and 1 and the cell-centre annulus 2<=r<=10 for base arrays and both directional epsilon tangent arms. It is not a claim of total-T closure.",
            "geometric_source": "Compare the pinned production momentum geometric source with an independently computed 1/2 sqrt(-g) T^{ab} partial_r g_ab contraction.",
            "hll_flux": "Record actual pinned per-species HLL face fluxes and compare them with face reconstruction of cell-projected radial stress. Differences characterize the operator pair; no exact face equality is presumed for the dissipative HLL scheme.",
            "interaction_source": "Compare actual pinned dark-matter q_s with the conditional covariant target sqrt(-g) Q_r, and separately retain a lapse-omitted negative control. The implemented covariant-target gap is not required to be nonzero; the negative-control gap must be nontrivial. Do not change production source terms in this diagnostic.",
            "cross_sector_interpretation": "The total-T directional witness is decomposed into dark-matter, baryon, radiation, S, D and phi stress contributions on the exact same state pair. The raw sector-divergence sum is retained beside the total without hiding any gap. A separate time-density derivative ledger predicts the portion of that raw gap attributable to forming one total density difference versus summing six per-sector density differences; only the remaining unexplained gap is checked against the unchanged tolerance. This validates numerical accounting for that component, not physical stress-energy closure. The full scalar evolution RHS/source terms are not reconciled, so full cross-sector closure remains incomplete.",
            "optional_em_scope": "Optional charged-dust and radial-Maxwell states are rejected before projection/capture because the Lorentz source and EM stress projection are outside the current witness.",
            "limits": [
                "This is a frozen-state RHS-component readout; no finite production timestep is advanced.",
                "The scalar-field spatial operator remains separate from the fluid HLL operator; scalar spatial/source decomposition remains incomplete.",
                "The report keeps raw sector-sum gaps and separately records the independently computed time-density component that explains them. The component-accounted residual gate does not replace or overwrite raw values, and it does not reconcile each scalar evolution RHS/source with its stress tensor or establish total stress-energy closure or physical exchange-law acceptance.",
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
