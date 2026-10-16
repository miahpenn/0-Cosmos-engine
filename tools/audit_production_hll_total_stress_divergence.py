"""Audit the source-pinned production HLL matter operator against total-T FV fluxes.

This is a frozen-state diagnostic. It never calls a production step or changes
engine/vendor source. Fluid HLL face fluxes and source terms are preserved by
species; scalar stress fluxes remain a separate reconstructed contribution.
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_ROOT_RAW = os.environ.get("COSMOS_PINNED_SOURCE_ROOT")
if not SOURCE_ROOT_RAW:
    raise RuntimeError("COSMOS_PINNED_SOURCE_ROOT must point to the exact pinned source checkout")
SOURCE_ROOT = Path(SOURCE_ROOT_RAW).resolve()
if not SOURCE_ROOT.is_dir():
    raise RuntimeError(f"pinned production source directory does not exist: {SOURCE_ROOT}")

# Load the merged, source-pinned independent total-T observer. It supplies the
# provenance checks, immutable-state fingerprint, independent projection
# reconstruction, and established cell-centred/FV comparator.
_AUDIT_SPEC = importlib.util.spec_from_file_location(
    "_issue4_pinned_total_stress_energy_audit",
    ROOT / "tools" / "audit_total_stress_energy_divergence.py",
)
if _AUDIT_SPEC is None or _AUDIT_SPEC.loader is None:
    raise ImportError("cannot load the merged total-stress-energy divergence auditor")
AUDIT = importlib.util.module_from_spec(_AUDIT_SPEC)
sys.modules[_AUDIT_SPEC.name] = AUDIT
_AUDIT_SPEC.loader.exec_module(AUDIT)

# The merged auditor orders the pinned source ahead of the diagnostic checkout.
from engine import production_kernel as production_kernel_module  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402
from engine.matter_system import (  # noqa: E402
    Species,
    _species_rhs_components,
    geometry_metric_derivatives,
)
from engine.scalar_system import KAPPA, scalar_rhs_arrays  # noqa: E402
from engine.v55_matter import dm_density, metric_slice_from_q  # noqa: E402

SPECIES = (
    ("dark_matter", Species.DARK_MATTER),
    ("baryons", Species.BARYON),
    ("radiation", Species.RADIATION),
)


def _array(name: str, value: Any, shape: tuple[int, ...] | None = None) -> np.ndarray:
    out = np.asarray(value, dtype=np.float64)
    if shape is not None and out.shape != shape:
        raise ValueError(f"{name} shape {out.shape} does not match {shape}")
    if not np.all(np.isfinite(out)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return out


def _profile_summary(values: np.ndarray) -> dict[str, Any]:
    x = _array("profile summary input", values)
    if x.ndim != 1 or x.size == 0:
        raise ValueError("profile summary requires a nonempty one-dimensional array")
    i = int(np.argmax(np.abs(x)))
    return {
        "max_abs": float(abs(x[i])),
        "signed_at_max_abs": float(x[i]),
        "max_index": i,
        "rms": float(np.sqrt(np.mean(x * x))),
        "cell_0": float(x[0]),
        "cell_1": float(x[1]) if x.size > 1 else None,
    }


def rhs_from_native_components(
    face_flux: np.ndarray,
    source_e: np.ndarray,
    source_s: np.ndarray,
    q_e: np.ndarray,
    q_s: np.ndarray,
    dr: float,
) -> dict[str, np.ndarray]:
    """Reconstruct the semi-discrete RHS from the five native operator outputs."""
    face_flux = _array("native face flux", face_flux)
    if face_flux.ndim != 2 or face_flux.shape[1] != 3 or face_flux.shape[0] < 2:
        raise ValueError("native face flux must have shape (n+1, 3), n>=1")
    n = face_flux.shape[0] - 1
    source_e = _array("geometric energy source", source_e, (n,))
    source_s = _array("geometric momentum source", source_s, (n,))
    q_e = _array("interaction energy source", q_e, (n,))
    q_s = _array("interaction momentum source", q_s, (n,))
    if not np.isfinite(dr) or dr <= 0.0:
        raise ValueError("dr must be positive and finite")
    return {
        "rest": -np.diff(face_flux[:, 0]) / float(dr),
        "energy_t": -np.diff(face_flux[:, 1]) / float(dr) + source_e + q_e,
        "momentum_r": -np.diff(face_flux[:, 2]) / float(dr) + source_s + q_s,
    }


def negative_control_rhs(
    face_flux: np.ndarray,
    source_e: np.ndarray,
    source_s: np.ndarray,
    q_e: np.ndarray,
    q_s: np.ndarray,
    dr: float,
) -> dict[str, Any]:
    """Preserve deliberately wrong sign and omitted-geometric-source controls."""
    normal = rhs_from_native_components(face_flux, source_e, source_s, q_e, q_s, dr)
    flipped = rhs_from_native_components(-np.asarray(face_flux), source_e, source_s, q_e, q_s, dr)
    zeros = np.zeros_like(np.asarray(source_e, dtype=np.float64))
    omitted = rhs_from_native_components(face_flux, zeros, zeros, q_e, q_s, dr)
    return {
        "sign_flipped_face_flux": {
            "rhs": {k: v.tolist() for k, v in flipped.items()},
            "delta_vs_native": {k: _profile_summary(flipped[k] - normal[k]) for k in normal},
        },
        "omitted_geometric_source": {
            "rhs": {k: v.tolist() for k, v in omitted.items()},
            "delta_vs_native": {k: _profile_summary(omitted[k] - normal[k]) for k in normal},
        },
    }


def _serialize_rates(rate: dict[str, Any]) -> dict[str, list[float]]:
    return {
        name: _array(f"scalar RHS {name}", values).tolist()
        for name, values in rate.items()
    }


def _check_zero_shift(state: Any, label: str) -> None:
    if np.any(np.asarray(state.geometry.beta, dtype=np.float64) != 0.0):
        raise RuntimeError(f"{label}: source audit currently requires the zero-shift baseline")


def _capture_state(
    kernel: V55ProductionKernel,
    state: Any,
    rates: dict[str, Any],
    label: str,
) -> dict[str, Any]:
    """Capture all native matter operators on a frozen state and test their RHS."""
    _check_zero_shift(state, label)
    state_hash_before = AUDIT.state_digest(state)
    grid = state.grid
    metric = metric_slice_from_q(grid, state.geometry)
    derivatives = geometry_metric_derivatives(
        metric, lambda values, parity: kernel._d1(grid, values, parity)
    )
    rho_dm = dm_density(metric, state.matter)
    scalar_rates = scalar_rhs_arrays(
        grid,
        state.geometry,
        state.scalars,
        beta_dm=production_kernel_module.BETA_DM,
        rho_dm=rho_dm,
    )
    dphi_t = scalar_rates.phi
    dphi_r = kernel._d1(grid, state.scalars.phi, +1)

    projection = AUDIT.independent_total_projections(state)
    r = _array("cell-centre radii", grid.centers)
    if r.size < 4:
        raise ValueError("the finite-volume HLL audit requires at least four cells")
    spacing = np.diff(r)
    dr = float(spacing[0])
    if not np.allclose(spacing, dr, rtol=5.0e-13, atol=5.0e-15):
        raise ValueError("the finite-volume HLL audit requires a uniform grid")

    alpha = _array("lapse", state.geometry.alpha, r.shape)
    gamma_rr = _array("gamma_rr", state.geometry.a / state.geometry.X**2, r.shape)
    gamma_thth = _array("gamma_thth", state.geometry.b * r**2 / state.geometry.X**2, r.shape)
    sqrt_gamma = np.sqrt(gamma_rr) * gamma_thth
    alpha_faces = np.zeros(r.size + 1, dtype=np.float64)
    alpha_faces[1:-1] = 0.5 * (alpha[:-1] + alpha[1:])
    alpha_faces[-1] = alpha[-1]

    species_reports: dict[str, Any] = {}
    species_arrays: dict[str, dict[str, np.ndarray]] = {}
    all_rhs_matches = True
    all_profiles_finite = True
    wrong_flux_effects: list[float] = []
    omitted_source_effects: list[float] = []

    for name, kind in SPECIES:
        conserved = getattr(state.matter, name)
        face_flux, source_e, source_s, q_e, q_s = _species_rhs_components(
            metric,
            derivatives,
            conserved,
            kind,
            dphi_t=dphi_t if kind is Species.DARK_MATTER else None,
            dphi_r=dphi_r if kind is Species.DARK_MATTER else None,
            beta_dm=production_kernel_module.BETA_DM,
        )
        face_flux = _array(f"{name} face flux", face_flux, (r.size + 1, 3))
        source_e = _array(f"{name} geometric energy source", source_e, r.shape)
        source_s = _array(f"{name} geometric momentum source", source_s, r.shape)
        q_e = _array(f"{name} interaction energy source", q_e, r.shape)
        q_s = _array(f"{name} interaction momentum source", q_s, r.shape)
        components_rhs = rhs_from_native_components(face_flux, source_e, source_s, q_e, q_s, dr)
        kernel_rhs = rates["matter"][name]
        rhs_errors = {}
        rhs_matches = True
        for channel in ("rest", "energy_t", "momentum_r"):
            diff = components_rhs[channel] - _array(
                f"{name} kernel RHS {channel}", kernel_rhs[channel], r.shape
            )
            rhs_errors[channel] = {
                "max_abs": float(np.max(np.abs(diff))),
                "allclose": bool(np.allclose(
                    components_rhs[channel], kernel_rhs[channel],
                    rtol=5.0e-13, atol=5.0e-18,
                )),
            }
            rhs_matches = rhs_matches and rhs_errors[channel]["allclose"]
        all_rhs_matches = all_rhs_matches and rhs_matches

        fluid_projection = {
            k: _array(f"{name} projected {k}", v, r.shape)
            for k, v in projection["fluid_parts"][name].items()
        }
        # Projection arrays use Einstein-normalized stress units (divide the
        # native Valencia/HLL arrays by KAPPA before comparing.
        reference_flux_t_cells = (
            -(alpha**2) * sqrt_gamma * fluid_projection["j"] / gamma_rr
        )
        reference_flux_r_cells = alpha * sqrt_gamma * fluid_projection["pr"]
        reference_flux_t_faces = AUDIT.spherical_flux_face_values(
            reference_flux_t_cells, dr, origin_power=3
        )
        reference_flux_r_faces = AUDIT.spherical_flux_face_values(
            reference_flux_r_cells, dr, origin_power=2
        )
        native_flux_t_faces = -alpha_faces * face_flux[:, 1] / KAPPA
        native_flux_r_faces = face_flux[:, 2] / KAPPA
        gap_t = native_flux_t_faces - reference_flux_t_faces
        gap_r = native_flux_r_faces - reference_flux_r_faces

        controls = negative_control_rhs(face_flux, source_e, source_s, q_e, q_s, dr)
        wrong_t = controls["sign_flipped_face_flux"]["delta_vs_native"]["energy_t"]["max_abs"]
        wrong_s = controls["sign_flipped_face_flux"]["delta_vs_native"]["momentum_r"]["max_abs"]
        omitted_t = controls["omitted_geometric_source"]["delta_vs_native"]["energy_t"]["max_abs"]
        omitted_s = controls["omitted_geometric_source"]["delta_vs_native"]["momentum_r"]["max_abs"]
        wrong_flux_effects.extend([wrong_t, wrong_s])
        omitted_source_effects.extend([omitted_t, omitted_s])

        species_arrays[name] = {
            "native_face_flux": face_flux,
            "source_e": source_e,
            "source_s": source_s,
            "q_e": q_e,
            "q_s": q_s,
            "native_flux_t_faces": native_flux_t_faces,
            "native_flux_r_faces": native_flux_r_faces,
        }
        species_reports[name] = {
            "native_outputs": {
                "face_flux_rest_energy_momentum": face_flux.tolist(),
                "geometric_energy_source": source_e.tolist(),
                "geometric_radial_momentum_source": source_s.tolist(),
                "interaction_energy_source": q_e.tolist(),
                "interaction_radial_momentum_source": q_s.tolist(),
            },
            "rhs_from_native_components": {
                k: _array(f"{name} reconstructed RHS {k}", v, r.shape).tolist()
                for k, v in components_rhs.items()
            },
            "comparison_to_production_rhs_wrapper": {
                "matches": rhs_matches,
                "channels": rhs_errors,
            },
            "independent_fluid_stress_projection": {
                k: v.tolist() for k, v in fluid_projection.items()
            },
            "hll_vs_projected_flux_faces": {
                "normalization": "native HLL components divided by KAPPA=8*pi",
                "time_mixed_component_native_faces": native_flux_t_faces.tolist(),
                "time_mixed_component_projected_reconstruction": reference_flux_t_faces.tolist(),
                "time_mixed_component_native_minus_projected": gap_t.tolist(),
                "time_mixed_component_gap_summary": _profile_summary(gap_t),
                "radial_mixed_component_native_faces": native_flux_r_faces.tolist(),
                "radial_mixed_component_projected_reconstruction": reference_flux_r_faces.tolist(),
                "radial_mixed_component_native_minus_projected": gap_r.tolist(),
                "radial_mixed_component_gap_summary": _profile_summary(gap_r),
            },
            "negative_controls": controls,
        }
        all_profiles_finite = all_profiles_finite and all(
            np.all(np.isfinite(x))
            for x in (
                face_flux, source_e, source_s, q_e, q_s,
                *components_rhs.values(), *fluid_projection.values(),
                native_flux_t_faces, native_flux_r_faces,
                reference_flux_t_faces, reference_flux_r_faces, gap_t, gap_r,
            )
        )

    scalar_parts = {
        sector: {
            key: _array(f"scalar projection {sector}.{key}", value, r.shape)
            for key, value in fields.items()
        }
        for sector, fields in projection["scalar_parts"].items()
    }
    # Scalar sectors do not use the fluid HLL face operator; keep their stress
    # flux reconstruction separate and only combine the two operators afterward.
    scalar_flux_t_faces: dict[str, np.ndarray] = {}
    scalar_flux_r_faces: dict[str, np.ndarray] = {}
    for sector, fields in scalar_parts.items():
        cell_t = -(alpha**2) * sqrt_gamma * fields["j"] / gamma_rr
        cell_r = alpha * sqrt_gamma * fields["pr"]
        scalar_flux_t_faces[sector] = AUDIT.spherical_flux_face_values(
            cell_t, dr, origin_power=3
        )
        scalar_flux_r_faces[sector] = AUDIT.spherical_flux_face_values(
            cell_r, dr, origin_power=2
        )

    hll_fluid_t_faces = sum(
        (items["native_flux_t_faces"] for items in species_arrays.values()),
        start=np.zeros(r.size + 1, dtype=np.float64),
    )
    hll_fluid_r_faces = sum(
        (items["native_flux_r_faces"] for items in species_arrays.values()),
        start=np.zeros(r.size + 1, dtype=np.float64),
    )
    scalar_t_faces = sum(scalar_flux_t_faces.values(), start=np.zeros(r.size + 1))
    scalar_r_faces = sum(scalar_flux_r_faces.values(), start=np.zeros(r.size + 1))
    hll_plus_scalar_t_faces = hll_fluid_t_faces + scalar_t_faces
    hll_plus_scalar_r_faces = hll_fluid_r_faces + scalar_r_faces

    total_rho = _array("total projected rho", projection["rho"], r.shape)
    total_pr = _array("total projected pr", projection["pr"], r.shape)
    total_j = _array("total projected j", projection["j"], r.shape)
    reference_total_flux_t_cells = -(alpha**2) * sqrt_gamma * total_j / gamma_rr
    reference_total_flux_r_cells = alpha * sqrt_gamma * total_pr
    reference_total_flux_t_faces = AUDIT.spherical_flux_face_values(
        reference_total_flux_t_cells, dr, origin_power=3
    )
    reference_total_flux_r_faces = AUDIT.spherical_flux_face_values(
        reference_total_flux_r_cells, dr, origin_power=2
    )

    # The HLL path is a different fluid face discretization. A sign-flipped
    # dark-matter face contribution is retained as a negative control.
    dm_arrays = species_arrays["dark_matter"]
    sign_flip_t_faces = (
        hll_plus_scalar_t_faces - 2.0 * dm_arrays["native_flux_t_faces"]
    )
    sign_flip_r_faces = (
        hll_plus_scalar_r_faces - 2.0 * dm_arrays["native_flux_r_faces"]
    )

    if AUDIT.state_digest(state) != state_hash_before:
        raise RuntimeError(f"{label}: native operator capture mutated the frozen input state")

    return {
        "label": label,
        "state_sha256": state_hash_before,
        "zero_shift": True,
        "projection_crosscheck_passed": bool(projection["projection_crosscheck_passed"]),
        "projection_crosscheck": projection["projection_crosscheck"],
        "native_species_rhs_decomposition_matches_kernel_rhs": all_rhs_matches,
        "production_steps_advanced": 0,
        "optional_sector_scope": {
            "charged_dust": "absent; this baseline audit fails closed if the optional carrier is present",
            "radial_maxwell": "absent; no electromagnetic stress/source was silently folded into matter",
        },
        "species": species_reports,
        "scalar_operator": {
            "rhs": _serialize_rates(rates["scalars"]),
            "stress_projection_by_sector": {
                sector: {key: value.tolist() for key, value in fields.items()}
                for sector, fields in scalar_parts.items()
            },
            "projected_reconstructed_face_fluxes_by_sector": {
                sector: {
                    "time_mixed_component": scalar_flux_t_faces[sector].tolist(),
                    "radial_mixed_component": scalar_flux_r_faces[sector].tolist(),
                }
                for sector in scalar_parts
            },
            "combined_scalar_time_mixed_faces": scalar_t_faces.tolist(),
            "combined_scalar_radial_mixed_faces": scalar_r_faces.tolist(),
        },
        "combined_faces": {
            "native_species_hll_fluid_time_mixed_component": hll_fluid_t_faces.tolist(),
            "native_species_hll_fluid_radial_mixed_component": hll_fluid_r_faces.tolist(),
            "scalar_reconstructed_time_mixed_component": scalar_t_faces.tolist(),
            "scalar_reconstructed_radial_mixed_component": scalar_r_faces.tolist(),
            "hll_fluid_plus_scalar_time_mixed_component": hll_plus_scalar_t_faces.tolist(),
            "hll_fluid_plus_scalar_radial_mixed_component": hll_plus_scalar_r_faces.tolist(),
            "total_projection_reconstructed_time_mixed_component": reference_total_flux_t_faces.tolist(),
            "total_projection_reconstructed_radial_mixed_component": reference_total_flux_r_faces.tolist(),
            "hll_plus_scalar_minus_total_projection_time": (hll_plus_scalar_t_faces - reference_total_flux_t_faces).tolist(),
            "hll_plus_scalar_minus_total_projection_radial": (hll_plus_scalar_r_faces - reference_total_flux_r_faces).tolist(),
            "hll_plus_scalar_minus_total_projection_time_summary": _profile_summary(
                hll_plus_scalar_t_faces - reference_total_flux_t_faces
            ),
            "hll_plus_scalar_minus_total_projection_radial_summary": _profile_summary(
                hll_plus_scalar_r_faces - reference_total_flux_r_faces
            ),
            "sign_flipped_dark_matter_face_control_time": sign_flip_t_faces.tolist(),
            "sign_flipped_dark_matter_face_control_radial": sign_flip_r_faces.tolist(),
        },
        "negative_control_sensitivity": {
            "sign_flipped_face_flux_max_effect": float(max(wrong_flux_effects, default=0.0)),
            "omitted_geometric_source_max_effect": float(max(omitted_source_effects, default=0.0)),
        },
        "all_profiles_finite": bool(all_profiles_finite),
    }


def covariant_fv_with_custom_flux_faces(
    r: np.ndarray,
    d1_fn,
    alpha: np.ndarray,
    gamma_rr: np.ndarray,
    gamma_thth: np.ndarray,
    rho: np.ndarray,
    j: np.ndarray,
    pr: np.ndarray,
    pt: np.ndarray,
    dt_density_t: np.ndarray,
    dt_density_r: np.ndarray,
    dt_gtt: np.ndarray,
    dt_gamma_rr: np.ndarray,
    dt_gamma_thth: np.ndarray,
    flux_t_faces: np.ndarray,
    flux_r_faces: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Reuse the independent FV measure/source quadrature with caller-supplied flux faces."""
    r = _array("custom FV radii", r)
    if r.ndim != 1 or r.size < 4:
        raise ValueError("custom FV divergence requires n>=4")
    spacing = np.diff(r)
    if not np.all(np.isfinite(spacing)) or np.any(spacing <= 0.0):
        raise ValueError("custom FV radii must be increasing and finite")
    dr = float(spacing[0])
    if not np.allclose(spacing, dr, rtol=5.0e-13, atol=5.0e-15):
        raise ValueError("custom FV divergence requires a uniform grid")
    if not np.isclose(r[0], 0.5 * dr, rtol=5.0e-12, atol=5.0e-14):
        raise ValueError("custom FV spherical grid must have first center dr/2")

    alpha = _array("custom FV lapse", alpha, r.shape)
    grr = _array("custom FV gamma_rr", gamma_rr, r.shape)
    gth = _array("custom FV gamma_thth", gamma_thth, r.shape)
    E = _array("custom FV rho", rho, r.shape)
    j = _array("custom FV j", j, r.shape)
    pr = _array("custom FV pr", pr, r.shape)
    pt = _array("custom FV pt", pt, r.shape)
    dt_t = _array("custom FV dt density t", dt_density_t, r.shape)
    dt_r = _array("custom FV dt density r", dt_density_r, r.shape)
    gtt_t = _array("custom FV dt gtt", dt_gtt, r.shape)
    grr_t = _array("custom FV dt gamma_rr", dt_gamma_rr, r.shape)
    gth_t = _array("custom FV dt gamma_thth", dt_gamma_thth, r.shape)
    ft = _array("custom FV time flux faces", flux_t_faces, (r.size + 1,))
    fr = _array("custom FV radial flux faces", flux_r_faces, (r.size + 1,))
    sqrt_gamma = np.sqrt(grr) * gth
    sqrt_minus_g = alpha * sqrt_gamma
    if np.any(alpha <= 0.0) or np.any(grr <= 0.0) or np.any(gth <= 0.0):
        raise ValueError("custom FV divergence requires positive lapse and spatial metric")

    Ttt = E / alpha**2
    Trr = pr / grr
    Tthth = pt / gth
    source_t = -0.5 * (Ttt * gtt_t + Trr * grr_t + 2.0 * Tthth * gth_t)
    gtt_r = d1_fn(-(alpha**2), +1)
    grr_r = d1_fn(grr, +1)
    gth_r = d1_fn(gth, +1)
    source_r = -0.5 * (Ttt * gtt_r + Trr * grr_r + 2.0 * Tthth * gth_r)

    source_density_t = sqrt_minus_g * source_t
    source_density_r = sqrt_minus_g * source_r
    volume_faces = AUDIT.spherical_flux_face_values(sqrt_minus_g, dr, origin_power=2)
    dt_t_faces = AUDIT.spherical_flux_face_values(dt_t, dr, origin_power=2)
    dt_r_faces = AUDIT.spherical_flux_face_values(dt_r, dr, origin_power=3)
    source_t_faces = AUDIT.spherical_flux_face_values(source_density_t, dr, origin_power=2)
    source_r_faces = AUDIT.spherical_flux_face_values(source_density_r, dr, origin_power=1)

    def simpson_cells(centers: np.ndarray, faces: np.ndarray, name: str) -> np.ndarray:
        profile = _array(f"custom FV cell-center {name}", centers, r.shape)
        faces = _array(f"custom FV faces {name}", faces, (r.size + 1,))
        return (dr / 6.0) * (faces[:-1] + 4.0 * profile + faces[1:])

    volume = simpson_cells(sqrt_minus_g, volume_faces, "volume")
    if np.any(volume <= 0.0):
        raise ValueError("custom FV cell volumes must be positive")
    int_dt_t = simpson_cells(dt_t, dt_t_faces, "dt density t")
    int_dt_r = simpson_cells(dt_r, dt_r_faces, "dt density r")
    int_source_t = simpson_cells(source_density_t, source_t_faces, "source t")
    int_source_r = simpson_cells(source_density_r, source_r_faces, "source r")
    div_t = (int_dt_t + np.diff(ft) + int_source_t) / volume
    div_r = (int_dt_r + np.diff(fr) + int_source_r) / volume
    return _array("custom FV div t", div_t, r.shape), _array("custom FV div r", div_r, r.shape)


def _tangent_inputs(kernel, base, rates, epsilon: float):
    plus_state = AUDIT.shifted_state(kernel, base, rates, +epsilon)
    minus_state = AUDIT.shifted_state(kernel, base, rates, -epsilon)
    _check_zero_shift(plus_state, f"plus tangent epsilon={epsilon:.17g}")
    _check_zero_shift(minus_state, f"minus tangent epsilon={epsilon:.17g}")
    plus_obj = AUDIT.covariant_objects(plus_state)
    minus_obj = AUDIT.covariant_objects(minus_state)
    inv = 1.0 / (2.0 * epsilon)
    return plus_state, minus_state, {
        "dt_density_t": (plus_obj["density_t"] - minus_obj["density_t"]) * inv,
        "dt_density_r": (plus_obj["density_r"] - minus_obj["density_r"]) * inv,
        "dt_gtt": (plus_obj["gtt"] - minus_obj["gtt"]) * inv,
        "dt_gamma_rr": (plus_obj["gamma_rr"] - minus_obj["gamma_rr"]) * inv,
        "dt_gamma_thth": (plus_obj["gamma_thth"] - minus_obj["gamma_thth"]) * inv,
    }


def _audit_resolution(n: int, r_max: float, epsilons: tuple[float, float]) -> dict[str, Any]:
    kernel = V55ProductionKernel()
    original_state = kernel.initialize(resolution=n, r_max=r_max, include_radiation=True)
    AUDIT.assert_supported_baseline_source_scope(original_state)
    if getattr(original_state.matter, "charged_dust", None) is not None:
        raise RuntimeError("optional charged dust is outside this baseline diagnostic")
    if getattr(original_state, "electromagnetic", None) is not None:
        raise RuntimeError("optional radial Maxwell is outside this baseline diagnostic")

    original_hash = AUDIT.state_digest(original_state)
    state = copy.deepcopy(original_state)
    state_hash = AUDIT.state_digest(state)
    if state_hash != original_hash:
        raise RuntimeError("copied frozen state does not match initialized state hash")
    base_rates = AUDIT.rhs_direction(kernel, state)
    base_capture = _capture_state(kernel, state, base_rates, "base")
    _check_zero_shift(state, "base")
    tangent_reports: dict[str, Any] = {}
    all_finite = bool(base_capture["all_profiles_finite"])
    all_rhs_matches = bool(base_capture["native_species_rhs_decomposition_matches_kernel_rhs"])
    control_significant = base_capture["negative_control_sensitivity"]["sign_flipped_face_flux_max_effect"] > 0.0
    source_control_significant = base_capture["negative_control_sensitivity"]["omitted_geometric_source_max_effect"] > 0.0
    r = _array("grid centers", state.grid.centers)
    base_obj = AUDIT.covariant_objects(state)
    d1_fn = lambda values, parity: AUDIT.d1(state.grid, values, parity)

    for epsilon in epsilons:
        plus_state, minus_state, tangent = _tangent_inputs(kernel, state, base_rates, epsilon)
        plus_hash_before = AUDIT.state_digest(plus_state)
        minus_hash_before = AUDIT.state_digest(minus_state)
        plus_rates = AUDIT.rhs_direction(kernel, plus_state)
        minus_rates = AUDIT.rhs_direction(kernel, minus_state)
        plus_capture = _capture_state(kernel, plus_state, plus_rates, f"plus epsilon={epsilon:.17g}")
        minus_capture = _capture_state(kernel, minus_state, minus_rates, f"minus epsilon={epsilon:.17g}")
        if AUDIT.state_digest(plus_state) != plus_hash_before:
            raise RuntimeError("plus directional state mutated during observer capture")
        if AUDIT.state_digest(minus_state) != minus_hash_before:
            raise RuntimeError("minus directional state mutated during observer capture")

        projection_args = (
            r, d1_fn,
            base_obj["alpha"], base_obj["gamma_rr"], base_obj["gamma_thth"],
            base_obj["rho"], base_obj["j"], base_obj["pr"], base_obj["pt"],
            tangent["dt_density_t"], tangent["dt_density_r"],
            tangent["dt_gtt"], tangent["dt_gamma_rr"], tangent["dt_gamma_thth"],
        )
        # Established comparator: total-T face flux reconstructed from the
        # independent cell-centred projections.
        ref_d1_t, ref_d1_r = AUDIT.covariant_divergence_from_profiles(*projection_args)
        ref_fv_t, ref_fv_r = AUDIT.covariant_divergence_fv_from_profiles(*projection_args)
        custom_ref_t, custom_ref_r = covariant_fv_with_custom_flux_faces(
            *projection_args,
            base_obj["flux_t"] if False else AUDIT.spherical_flux_face_values(
                base_obj["flux_t"], float(np.diff(r)[0]), origin_power=3
            ),
            AUDIT.spherical_flux_face_values(
                base_obj["flux_r"], float(np.diff(r)[0]), origin_power=2
            ),
        )
        ref_control_close = bool(
            np.allclose(custom_ref_t, ref_fv_t, rtol=2.0e-12, atol=2.0e-18)
            and np.allclose(custom_ref_r, ref_fv_r, rtol=2.0e-12, atol=2.0e-18)
        )
        if not ref_control_close:
            raise RuntimeError("custom-flux FV wrapper disagrees with the established total-T FV comparator")

        hll_faces = base_capture["combined_faces"]
        hll_t = np.asarray(hll_faces["hll_fluid_plus_scalar_time_mixed_component"], dtype=np.float64)
        hll_r = np.asarray(hll_faces["hll_fluid_plus_scalar_radial_mixed_component"], dtype=np.float64)
        hll_d_t, hll_d_r = covariant_fv_with_custom_flux_faces(*projection_args, hll_t, hll_r)
        flip_t = np.asarray(hll_faces["sign_flipped_dark_matter_face_control_time"], dtype=np.float64)
        flip_r = np.asarray(hll_faces["sign_flipped_dark_matter_face_control_radial"], dtype=np.float64)
        flip_d_t, flip_d_r = covariant_fv_with_custom_flux_faces(*projection_args, flip_t, flip_r)

        tangent_reports[f"{epsilon:.17g}"] = {
            "epsilon": float(epsilon),
            "plus_state_sha256": plus_hash_before,
            "minus_state_sha256": minus_hash_before,
            "plus_state_capture": plus_capture,
            "minus_state_capture": minus_capture,
            "time_derivative_witnesses": {
                name: _profile_summary(values) for name, values in tangent.items()
            },
            "cell_centered_total_T_comparator": {
                "div_T_t": AUDIT.summarize(ref_d1_t, r),
                "div_T_r": AUDIT.summarize(ref_d1_r, r),
            },
            "projected_total_T_FV_comparator": {
                "div_T_t": AUDIT.summarize(ref_fv_t, r),
                "div_T_r": AUDIT.summarize(ref_fv_r, r),
                "custom_flux_reconstruction_crosscheck_passed": ref_control_close,
            },
            "production_HLL_fluid_plus_scalar_projected_FV": {
                "div_T_t": AUDIT.summarize(hll_d_t, r),
                "div_T_r": AUDIT.summarize(hll_d_r, r),
                "delta_vs_projected_total_T_FV_t": AUDIT.summarize(hll_d_t - ref_fv_t, r),
                "delta_vs_projected_total_T_FV_r": AUDIT.summarize(hll_d_r - ref_fv_r, r),
            },
            "negative_control_sign_flipped_dark_matter_HLL_faces": {
                "div_T_t": AUDIT.summarize(flip_d_t, r),
                "div_T_r": AUDIT.summarize(flip_d_r, r),
                "delta_vs_native_HLL_FV_t": AUDIT.summarize(flip_d_t - hll_d_t, r),
                "delta_vs_native_HLL_FV_r": AUDIT.summarize(flip_d_r - hll_d_r, r),
            },
            "all_profiles_finite": all(
                np.all(np.isfinite(x))
                for x in (
                    ref_d1_t, ref_d1_r, ref_fv_t, ref_fv_r, custom_ref_t,
                    custom_ref_r, hll_d_t, hll_d_r, flip_d_t, flip_d_r,
                    *tangent.values(),
                )
            ),
        }
        all_finite = all_finite and bool(
            plus_capture["all_profiles_finite"] and minus_capture["all_profiles_finite"]
            and tangent_reports[f"{epsilon:.17g}"]["all_profiles_finite"]
        )
        all_rhs_matches = all_rhs_matches and bool(
            plus_capture["native_species_rhs_decomposition_matches_kernel_rhs"]
            and minus_capture["native_species_rhs_decomposition_matches_kernel_rhs"]
        )
        control_significant = control_significant or bool(
            AUDIT.summarize(flip_d_t - hll_d_t, r)["max_abs"] > 0.0
            or AUDIT.summarize(flip_d_r - hll_d_r, r)["max_abs"] > 0.0
        )
        source_control_significant = source_control_significant or bool(
            plus_capture["negative_control_sensitivity"]["omitted_geometric_source_max_effect"] > 0.0
            or minus_capture["negative_control_sensitivity"]["omitted_geometric_source_max_effect"] > 0.0
        )

    if AUDIT.state_digest(state) != state_hash or AUDIT.state_digest(original_state) != original_hash:
        raise RuntimeError("base/original frozen state mutated during diagnostic")
    return {
        "resolution": int(n),
        "r_max": float(r_max),
        "base_state_sha256": state_hash,
        "original_state_sha256": original_hash,
        "base_capture": base_capture,
        "directional_tangent_checks": tangent_reports,
        "gate_result": {
            "input_states_unchanged": True,
            "production_steps_advanced": 0,
            "zero_shift_baseline_and_tangents": True,
            "independent_projection_matches_canonical": bool(base_capture["projection_crosscheck_passed"]),
            "native_species_rhs_decomposition_matches_kernel_rhs": bool(all_rhs_matches),
            "all_profiles_finite": bool(all_finite),
            "sign_flipped_face_control_changes_operator_balance": bool(control_significant),
            "omitted_geometric_source_control_changes_operator_rhs": bool(source_control_significant),
            "projected_total_T_crosscheck_matches_custom_FV": all(
                row["projected_total_T_FV_comparator"]["custom_flux_reconstruction_crosscheck_passed"]
                for row in tangent_reports.values()
            ),
            "production_HLL_total_stress_energy_divergence_closed": False,
            "scientific_admission": False,
        },
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    provenance = AUDIT.verify_provenance()
    results = [_audit_resolution(n, args.r_max, tuple(args.epsilons)) for n in args.resolutions]
    gate = {
        "source_pins_verified": provenance.get("source_blobs_match", False),
        "independent_projection_matches_canonical": all(
            item["gate_result"]["independent_projection_matches_canonical"] for item in results
        ),
        "native_species_rhs_decomposition_matches_kernel_rhs": all(
            item["gate_result"]["native_species_rhs_decomposition_matches_kernel_rhs"] for item in results
        ),
        "all_profiles_finite": all(item["gate_result"]["all_profiles_finite"] for item in results),
        "input_states_unchanged": all(item["gate_result"]["input_states_unchanged"] for item in results),
        "production_steps_advanced": 0,
        "sign_flipped_face_control_changes_operator_balance": all(
            item["gate_result"]["sign_flipped_face_control_changes_operator_balance"] for item in results
        ),
        "omitted_geometric_source_control_changes_operator_rhs": all(
            item["gate_result"]["omitted_geometric_source_control_changes_operator_rhs"] for item in results
        ),
        "projected_total_T_crosscheck_matches_custom_FV": all(
            item["gate_result"]["projected_total_T_crosscheck_matches_custom_FV"] for item in results
        ),
        "production_HLL_total_stress_energy_divergence_closed": False,
        "scientific_admission": False,
    }
    integrity = all(
        gate[key] is True for key in (
            "source_pins_verified", "independent_projection_matches_canonical",
            "native_species_rhs_decomposition_matches_kernel_rhs", "all_profiles_finite",
            "input_states_unchanged", "sign_flipped_face_control_changes_operator_balance",
            "omitted_geometric_source_control_changes_operator_rhs",
            "projected_total_T_crosscheck_matches_custom_FV",
        )
    )
    return {
        "schema": "issue4_production_hll_total_stress_divergence_v1",
        "status": "completed_diagnostic_only" if integrity else "completed_with_failed_integrity_gate",
        "source_provenance": provenance,
        "production_source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "diagnostic_commit": AUDIT.git(ROOT, "rev-parse", "HEAD"),
        "configuration": {
            "resolutions": list(args.resolutions),
            "r_max": float(args.r_max),
            "epsilons": list(args.epsilons),
            "coupling": float(production_kernel_module.BETA_DM),
            "baseline": "initialized frozen state; zero shift; no production step",
            "operator_capture": "engine.matter_system._species_rhs_components",
            "fluid_flux_normalization": "divide native Valencia/HLL fluxes by KAPPA=8*pi for Einstein-normalized T comparison",
            "scalar_operator": "kept separate; stress flux reconstructed from independent scalar projections, not assigned a fluid HLL flux",
        },
        "gate_result": gate,
        "resolution_results": results,
        "interpretation": (
            "The observer captures native per-species HLL face arrays and native geometric/interaction source arrays, "
            "checks their discrete RHS decomposition against the pinned production wrapper, and compares the actual "
            "HLL fluid face operator plus separately reconstructed scalar stress flux with the independent total-T FV "
            "comparator. Face-gap/residual differences are measured discretization evidence and are not presumed to "
            "close. Sign-flipped-face and omitted-geometric-source controls are retained. No production step or "
            "equation change occurs; the conditional exchange premise, full evolved-state conservation, and scientific "
            "admission remain unestablished."
        ),
        "scientific_admission": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[60, 80, 120, 160])
    parser.add_argument("--r-max", type=float, default=24.0)
    parser.add_argument("--epsilons", nargs=2, type=float, default=[2.0e-5, 1.0e-5])
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/issue4-production-hll-total-stress-divergence/report.json"),
    )
    args = parser.parse_args()
    if not args.resolutions or any(n < 60 for n in args.resolutions):
        parser.error("resolutions must be nonempty and n>=60 for the pinned audit")
    if not np.isfinite(args.r_max) or args.r_max <= 10.0:
        parser.error("r-max must exceed the fixed-annulus outer radius 10")
    if any(not np.isfinite(e) or e <= 0.0 for e in args.epsilons) or args.epsilons[0] == args.epsilons[1]:
        parser.error("epsilons must be distinct positive finite values")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        report = run(args)
        exit_code = 0 if report["status"] == "completed_diagnostic_only" else 1
    except Exception as exc:
        report = {
            "schema": "issue4_production_hll_total_stress_divergence_v1",
            "status": "failed_closed",
            "production_source_commit": SOURCE_COMMIT,
            "vendor_commit": VENDOR_COMMIT,
            "scientific_admission": False,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        exit_code = 1
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report.get("status"),
        "source_commit": report.get("production_source_commit"),
        "vendor_commit": report.get("vendor_commit"),
        "scientific_admission": False,
        "output": str(args.output),
        "error": report.get("error"),
    }, indent=2))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
