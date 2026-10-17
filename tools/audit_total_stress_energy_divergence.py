#!/usr/bin/env python3
"""Independent frozen-state covariant total-stress-energy divergence audit.

This is a diagnostic-only reconstruction of
    nabla_mu T^mu{}_nu =
      (1/sqrt(-g)) d_mu(sqrt(-g) T^mu{}_nu)
      - (1/2) T^{ab} d_nu g_ab
for spherical ADM data with beta^r=0 on the measured slice. Total stress
projections are independently rebuilt from scalar fields and fluid primitives,
then cross-checked against the canonical projection used by the production
consumers. Directional time derivatives are measured from two-sided
perturbations along the pinned semidiscrete RHS, including a CMC lapse
re-solve. We report the coupled result against a beta_DM=0 control; we do not
declare the physical exchange law accepted and never change production code.
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
import traceback
from typing import Any, Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_BLOBS = {
    "engine/production_kernel.py": "9685e76bb5a6f26d6bf383e240de24f37523f265",
    "engine/matter_system.py": "e52a117690253ddad7b24d116145a65474a351d2",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/v55_pirk_adapter.py": "eeee0492d80889d70e52b0946563eee5054c0f11",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
}
BETA_DM_COUPLED = -0.04

SOURCE_ROOT_RAW = os.environ.get("COSMOS_PINNED_SOURCE_ROOT")
if not SOURCE_ROOT_RAW:
    raise RuntimeError("COSMOS_PINNED_SOURCE_ROOT must point to the exact pinned source checkout")
SOURCE_ROOT = Path(SOURCE_ROOT_RAW).resolve()
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SOURCE_ROOT))

from engine import production_kernel as production_kernel_module  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402
from engine.scalar_system import (
    KAPPA, ScalarFields, cosmos_potential, cosmos_potential_prime, scalar_rhs_arrays,
)  # noqa: E402
from engine.matter_system import Species, _metric_arrays, primitives  # noqa: E402
from engine.stress_energy import assemble_total_stress_energy  # noqa: E402
from engine.v55_matter import dm_density, metric_slice_from_q  # noqa: E402
import importlib.util as importlib_util  # noqa: E402

_HELPER_PATH = ROOT / "tools" / "audit_current_source_momentum_constraint_tangent.py"
_HELPER_SPEC = importlib_util.spec_from_file_location(
    "_issue4_pinned_current_source_tangent_helper", _HELPER_PATH
)
if _HELPER_SPEC is None or _HELPER_SPEC.loader is None:
    raise ImportError(f"cannot load pinned tangent helper: {_HELPER_PATH}")
_HELPER_MODULE = importlib_util.module_from_spec(_HELPER_SPEC)
sys.modules[_HELPER_SPEC.name] = _HELPER_MODULE
_HELPER_SPEC.loader.exec_module(_HELPER_MODULE)
rhs_direction = _HELPER_MODULE.rhs_direction
shifted_state = _HELPER_MODULE.shifted_state
state_digest = _HELPER_MODULE.state_digest
assert_supported_baseline_source_scope = _HELPER_MODULE.assert_supported_baseline_source_scope


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def finite(name: str, value: Any) -> np.ndarray:
    out = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return out


def verify_provenance() -> dict[str, Any]:
    commit = git(SOURCE_ROOT, "rev-parse", "HEAD")
    vendor = git(SOURCE_ROOT / "vendor" / "bb-palatini-unified-r0", "rev-parse", "HEAD")
    blobs = {path: git(SOURCE_ROOT, "rev-parse", f"HEAD:{path}") for path in SOURCE_BLOBS}
    mismatches = {
        path: {"expected": sha, "actual": blobs[path]}
        for path, sha in SOURCE_BLOBS.items() if blobs[path] != sha
    }
    if commit != SOURCE_COMMIT or vendor != VENDOR_COMMIT or mismatches:
        raise RuntimeError(json.dumps({
            "source_commit": commit, "expected_source_commit": SOURCE_COMMIT,
            "vendor_commit": vendor, "expected_vendor_commit": VENDOR_COMMIT,
            "blob_mismatches": mismatches,
        }, indent=2))
    return {
        "source_commit": commit, "vendor_commit": vendor,
        "source_blobs": blobs, "source_blobs_match": True,
        "diagnostic_commit": git(ROOT, "rev-parse", "HEAD"),
    }


def d1(grid, values, parity: int) -> np.ndarray:
    return finite("radial derivative", grid.cell_derivative_fourth(values, parity=parity))


def d2(grid, values, parity: int) -> np.ndarray:
    """Use the same grid second-derivative interface as the pinned scalar RHS."""
    if hasattr(grid, "cell_second_derivative_fourth"):
        return finite(
            "radial second derivative",
            grid.cell_second_derivative_fourth(values, parity=parity),
        )
    return finite("radial second derivative", grid.d2(values, parity))


def independent_total_projections(state) -> dict[str, np.ndarray]:
    """Rebuild E, S_r, p_r and p_t from scalar fields/fluid primitives.

    This intentionally does not call scalar_projection, species_projection or
    total_matter_projection. It cross-checks its result against the canonical
    stress_energy assembly after constructing each contribution explicitly.
    """
    grid = state.grid
    geom = state.geometry
    r = finite("r", grid.centers)
    invr = finite("X^2/a", geom.X**2 / geom.a)
    rho = np.zeros_like(r)
    pr = np.zeros_like(r)
    pt = np.zeros_like(r)
    j = np.zeros_like(r)

    scalar_sectors = (
        (state.scalars.S, state.scalars.PS, 0.5 * state.scalars.S**2, 1.0, "S"),
        (state.scalars.D, state.scalars.PD, -0.5 * state.scalars.D**2, 1.0, "D"),
        (state.scalars.phi, state.scalars.Pi, cosmos_potential(state.scalars.phi), KAPPA, "phi"),
    )
    scalar_parts = {}
    for field, momentum, potential, norm, name in scalar_sectors:
        field = finite(f"{name} field", field)
        momentum = finite(f"{name} conjugate momentum", momentum)
        field_r = d1(grid, field, +1)
        grad2 = invr * field_r**2
        e_part = (0.5 * momentum**2 + 0.5 * grad2 + potential) / norm
        pr_part = (0.5 * momentum**2 + 0.5 * grad2 - potential) / norm
        pt_part = (0.5 * momentum**2 - 0.5 * grad2 - potential) / norm
        j_part = (-momentum * field_r) / norm
        rho += e_part
        pr += pr_part
        pt += pt_part
        j += j_part
        scalar_parts[name] = {
            "rho": finite(f"{name} rho", e_part),
            "pr": finite(f"{name} pr", pr_part),
            "pt": finite(f"{name} pt", pt_part),
            "j": finite(f"{name} j", j_part),
        }

    metric = metric_slice_from_q(grid, geom)
    metrics = _metric_arrays(metric)
    fluid_parts = {}
    fluid_states = (
        ("dark_matter", state.matter.dark_matter, Species.DARK_MATTER),
        ("baryons", state.matter.baryons, Species.BARYON),
        ("radiation", state.matter.radiation, Species.RADIATION),
    )
    for name, conserved, species in fluid_states:
        prims = primitives(metrics, conserved, species)
        arrays = {key: np.empty_like(r) for key in ("rho", "pr", "pt", "j")}
        for i, (m, q) in enumerate(zip(metrics, prims)):
            W = q.lorentz()
            h = q.rho + q.pressure
            eulerian_E = h * W * W - q.pressure
            momentum_r = h * W * W * m.gamma_rr * q.v_r
            arrays["rho"][i] = eulerian_E / KAPPA
            arrays["pr"][i] = (h * W * W * q.v2 + q.pressure) / KAPPA
            arrays["pt"][i] = q.pressure / KAPPA
            arrays["j"][i] = momentum_r / KAPPA
        for key in arrays:
            arrays[key] = finite(f"{name} {key}", arrays[key])
            if key == "rho": rho += arrays[key]
            elif key == "pr": pr += arrays[key]
            elif key == "pt": pt += arrays[key]
            else: j += arrays[key]
        fluid_parts[name] = arrays

    independent = {
        "rho": finite("independent total rho", rho),
        "pr": finite("independent total radial stress", pr),
        "pt": finite("independent total tangential stress", pt),
        "j": finite("independent total radial current", j),
    }
    canonical = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    comparisons = {}
    for key in independent:
        delta = independent[key] - np.asarray(getattr(canonical, key), dtype=np.float64)
        comparisons[key] = {
            "max_abs": float(np.max(np.abs(delta))),
            "rms": float(np.sqrt(np.mean(delta**2))),
            "matches_canonical": bool(np.allclose(
                independent[key], np.asarray(getattr(canonical, key)),
                rtol=5.0e-12, atol=2.0e-20,
            )),
        }
    return {
        **independent,
        "scalar_parts": scalar_parts,
        "fluid_parts": fluid_parts,
        "projection_crosscheck": comparisons,
        "projection_crosscheck_passed": all(x["matches_canonical"] for x in comparisons.values()),
    }


def spherical_flux_face_values(
    values: np.ndarray, dr: float, *, origin_power: int = 2
) -> np.ndarray:
    """Reconstruct regular spherical profiles to faces in the parity variable r^2.

    For a regular spherical profile F(r)=r^p*g(r^2), divide out its known
    origin power p, then interpolate g as a cubic polynomial in r^2 using
    four neighboring cell centers at every face. This honours regular parity
    at the origin and avoids the lower-order near-origin error of ordinary-r
    cubic interpolation. The outer face uses a one-sided cubic extrapolation.
    """
    values = finite("finite-volume cell-centered profile", values)
    if values.ndim != 1 or values.size < 4:
        raise ValueError("finite-volume face reconstruction requires n>=4")
    if not math.isfinite(float(dr)) or float(dr) <= 0.0:
        raise ValueError("finite-volume dr must be positive and finite")
    if origin_power not in (1, 2, 3):
        raise ValueError("regular spherical origin power must be 1, 2 or 3")
    n = values.size
    faces = np.zeros(n + 1, dtype=np.float64)
    centers = (np.arange(n, dtype=np.float64) + 0.5) * float(dr)
    regularized = finite(
        "finite-volume parity-regularized profile",
        values / (centers ** origin_power),
    )
    # Face 0 is exactly zero for all supported regular powers p>0. At every
    # positive face, use a 4-point Lagrange interpolant for g(r^2).
    for j in range(1, n + 1):
        start = min(max(j - 2, 0), n - 4)
        nodes = centers[start:start + 4] ** 2
        target = (float(j) * float(dr)) ** 2
        weights = np.ones(4, dtype=np.float64)
        for k in range(4):
            for ell in range(4):
                if ell != k:
                    weights[k] *= (target - nodes[ell]) / (nodes[k] - nodes[ell])
        faces[j] = (float(j) * float(dr)) ** origin_power * float(
            weights @ regularized[start:start + 4]
        )
    return finite("finite-volume reconstructed face profile", faces)


def finite_volume_radial_flux_derivative(
    values: np.ndarray, dr: float, *, origin_power: int = 2
) -> np.ndarray:
    """Cellwise derivative from telescoping reconstructed spherical face fluxes."""
    faces = spherical_flux_face_values(values, dr, origin_power=origin_power)
    return finite("finite-volume radial flux derivative", np.diff(faces) / float(dr))


def covariant_divergence_from_profiles(
    r: np.ndarray,
    d1_fn: Callable[[np.ndarray, int], np.ndarray],
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
) -> tuple[np.ndarray, np.ndarray]:
    """Return (nabla_mu T^mu_t, nabla_mu T^mu_r) for zero shift.

    Mixed-tensor divergence identity:
      D_nu = (d_t(sqrt(-g) T^t_nu)+d_r(sqrt(-g) T^r_nu))/sqrt(-g)
             - 1/2 T^{ab} d_nu g_ab.
    Angular terms are represented by the factor 2*T^{theta theta}*d g_thth.
    """
    r = finite("divergence radii", r)
    alpha = finite("divergence lapse", alpha)
    grr = finite("divergence gamma_rr", gamma_rr)
    gth = finite("divergence gamma_thth", gamma_thth)
    E = finite("divergence Eulerian rho", rho)
    j = finite("divergence radial current", j)
    pr = finite("divergence radial pressure", pr)
    pt = finite("divergence tangential pressure", pt)
    sqrt_gamma = np.sqrt(grr) * gth
    sqrt_minus_g = alpha * sqrt_gamma
    if np.any(alpha <= 0.0) or np.any(grr <= 0.0) or np.any(gth <= 0.0):
        raise ValueError("divergence requires positive lapse and spatial metric")

    Ttt = E / (alpha**2)
    Trr = pr / grr
    Tthth = pt / gth
    j_up = j / grr

    # Coordinate mixed components for beta^r=0:
    # T^t_t=-E, T^r_t=-alpha*S^r, T^t_r=j/alpha, T^r_r=p_r.
    flux_t = -(alpha**2) * sqrt_gamma * j_up
    flux_r = alpha * sqrt_gamma * pr
    radial_flux_t = finite("D_r(sqrt(-g) T^r_t)", d1_fn(flux_t, -1))
    radial_flux_r = finite("D_r(sqrt(-g) T^r_r)", d1_fn(flux_r, +1))

    source_t = -0.5 * (
        Ttt * finite("dt g_tt", dt_gtt)
        + Trr * finite("dt gamma_rr", dt_gamma_rr)
        + 2.0 * Tthth * finite("dt gamma_thth", dt_gamma_thth)
    )
    gtt_r = d1_fn(-(alpha**2), +1)
    grr_r = d1_fn(grr, +1)
    gth_r = d1_fn(gth, +1)
    source_r = -0.5 * (
        Ttt * gtt_r + Trr * grr_r + 2.0 * Tthth * gth_r
    )
    denom = finite("sqrt(-g) divergence denominator", sqrt_minus_g)
    div_t = finite(
        "covariant time-component divergence",
        (finite("dt sqrt(-g) T^t_t", dt_density_t) + radial_flux_t) / denom + source_t,
    )
    div_r = finite(
        "covariant radial-component divergence",
        (finite("dt sqrt(-g) T^t_r", dt_density_r) + radial_flux_r) / denom + source_r,
    )
    return div_t, div_r


def covariant_divergence_fv_from_profiles(
    r: np.ndarray,
    d1_fn: Callable[[np.ndarray, int], np.ndarray],
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
) -> tuple[np.ndarray, np.ndarray]:
    """Volume-integrated covariant divergence with telescoping face fluxes.

    Unlike the pointwise D1 observer, this companion integrates the time-density
    and metric-source densities over each spherical cell with Simpson quadrature,
    divides by the integrated sqrt(-g) measure, and uses a conservative
    difference of reconstructed total-T radial face fluxes. This is a diagnostic
    reconstruction from source-pinned cell-centered projections, not the
    production per-species HLL face flux.
    """
    r = finite("finite-volume divergence radii", r)
    if r.ndim != 1 or r.size < 4:
        raise ValueError("finite-volume divergence requires n>=4")
    spacing = np.diff(r)
    if not np.all(np.isfinite(spacing)) or np.any(spacing <= 0.0):
        raise ValueError("finite-volume radii must be increasing and finite")
    dr = float(spacing[0])
    if not np.allclose(spacing, dr, rtol=5.0e-13, atol=5.0e-15):
        raise ValueError("finite-volume divergence requires a uniform grid")
    if not math.isclose(float(r[0]), 0.5 * dr, rel_tol=5.0e-12, abs_tol=5.0e-14):
        raise ValueError("finite-volume spherical grid must have first center dr/2")

    alpha = finite("finite-volume lapse", alpha)
    grr = finite("finite-volume gamma_rr", gamma_rr)
    gth = finite("finite-volume gamma_thth", gamma_thth)
    E = finite("finite-volume Eulerian rho", rho)
    j = finite("finite-volume radial current", j)
    pr = finite("finite-volume radial pressure", pr)
    pt = finite("finite-volume tangential pressure", pt)
    sqrt_gamma = np.sqrt(grr) * gth
    sqrt_minus_g = finite("finite-volume sqrt(-g) measure", alpha * sqrt_gamma)
    if np.any(alpha <= 0.0) or np.any(grr <= 0.0) or np.any(gth <= 0.0):
        raise ValueError("finite-volume divergence requires positive lapse and spatial metric")

    Ttt = E / (alpha**2)
    Trr = pr / grr
    Tthth = pt / gth
    flux_t = -(alpha**2) * sqrt_gamma * (j / grr)
    flux_r = alpha * sqrt_gamma * pr
    flux_t_faces = spherical_flux_face_values(flux_t, dr, origin_power=3)
    flux_r_faces = spherical_flux_face_values(flux_r, dr, origin_power=2)

    source_t = -0.5 * (
        Ttt * finite("FV dt g_tt", dt_gtt)
        + Trr * finite("FV dt gamma_rr", dt_gamma_rr)
        + 2.0 * Tthth * finite("FV dt gamma_thth", dt_gamma_thth)
    )
    gtt_r = d1_fn(-(alpha**2), +1)
    grr_r = d1_fn(grr, +1)
    gth_r = d1_fn(gth, +1)
    source_r = -0.5 * (
        Ttt * gtt_r + Trr * grr_r + 2.0 * Tthth * gth_r
    )

    # Integrate each cell's conservative densities with Simpson's rule. The
    # source-density radial component scales as r near the origin because
    # T^theta_theta*d_r(gamma_theta_theta)~1/r while sqrt(-g)~r^2.
    dt_t = finite("FV dt density t", dt_density_t)
    dt_r = finite("FV dt density r", dt_density_r)
    source_density_t = finite("FV time metric-source density", sqrt_minus_g * source_t)
    source_density_r = finite("FV radial metric-source density", sqrt_minus_g * source_r)

    volume_faces = spherical_flux_face_values(sqrt_minus_g, dr, origin_power=2)
    dt_t_faces = spherical_flux_face_values(dt_t, dr, origin_power=2)
    dt_r_faces = spherical_flux_face_values(dt_r, dr, origin_power=3)
    source_t_faces = spherical_flux_face_values(source_density_t, dr, origin_power=2)
    source_r_faces = spherical_flux_face_values(source_density_r, dr, origin_power=1)

    def simpson_cells(centers: np.ndarray, faces: np.ndarray, name: str) -> np.ndarray:
        profile = finite(f"FV cell-center {name}", centers)
        face_profile = finite(f"FV face {name}", faces)
        if face_profile.size != profile.size + 1:
            raise ValueError(f"FV face/cell shape mismatch for {name}")
        return finite(f"FV cell integral {name}", (float(dr) / 6.0) * (
            face_profile[:-1] + 4.0 * profile + face_profile[1:]
        ))

    cell_volume = simpson_cells(sqrt_minus_g, volume_faces, "sqrt(-g) measure")
    int_dt_t = simpson_cells(dt_t, dt_t_faces, "time-component density rate")
    int_dt_r = simpson_cells(dt_r, dt_r_faces, "radial-component density rate")
    int_source_t = simpson_cells(source_density_t, source_t_faces, "time metric-source density")
    int_source_r = simpson_cells(source_density_r, source_r_faces, "radial metric-source density")
    if np.any(cell_volume <= 0.0):
        raise ValueError("finite-volume spherical cell measures must be positive")

    div_t = finite(
        "finite-volume averaged covariant time-component divergence",
        (int_dt_t + np.diff(flux_t_faces) + int_source_t) / cell_volume,
    )
    div_r = finite(
        "finite-volume averaged covariant radial-component divergence",
        (int_dt_r + np.diff(flux_r_faces) + int_source_r) / cell_volume,
    )
    return div_t, div_r


def _covariant_objects_from_components(state, components: dict[str, Any]) -> dict[str, Any]:
    """Build mixed-tensor density/flux components from one additive stress sector."""
    r = finite("covariant grid radii", state.grid.centers)
    alpha = finite("covariant lapse", state.geometry.alpha)
    gamma_rr = finite("gamma_rr", state.geometry.a / state.geometry.X**2)
    gamma_thth = finite("gamma_thth", state.geometry.b * r**2 / state.geometry.X**2)
    rho = finite("covariant sector rho", components["rho"])
    pr = finite("covariant sector radial pressure", components["pr"])
    pt = finite("covariant sector tangential pressure", components["pt"])
    j = finite("covariant sector radial current", components["j"])
    sqrt_gamma = np.sqrt(gamma_rr) * gamma_thth
    sqrt_minus_g = alpha * sqrt_gamma
    j_up = j / gamma_rr
    return {
        "rho": rho, "pr": pr, "pt": pt, "j": j,
        "alpha": alpha, "gamma_rr": gamma_rr, "gamma_thth": gamma_thth,
        "sqrt_gamma": sqrt_gamma, "sqrt_minus_g": sqrt_minus_g,
        "density_t": -sqrt_minus_g * rho,
        "density_r": sqrt_gamma * j,
        "flux_t": -(alpha**2) * sqrt_gamma * j_up,
        "flux_r": alpha * sqrt_gamma * pr,
        "gtt": -(alpha**2),
        "projection_crosscheck": components.get("projection_crosscheck", {}),
        "projection_crosscheck_passed": bool(components.get("projection_crosscheck_passed", True)),
    }


def covariant_objects(state) -> dict[str, Any]:
    projection = independent_total_projections(state)
    return _covariant_objects_from_components(state, projection)

def summarize(values: np.ndarray, r: np.ndarray) -> dict[str, Any]:
    x = finite("profile summary", values)
    i = int(np.argmax(np.abs(x)))
    annulus = (r >= 2.0) & (r <= 10.0)
    fixed = {
        str(target): float(np.interp(target, r, x))
        for target in (2.0, 3.0, 5.0, 8.0) if r[0] <= target <= r[-1]
    }
    return {
        "max_abs": float(np.abs(x[i])),
        "signed_at_max_abs": float(x[i]),
        "max_index": i,
        "max_radius": float(r[i]),
        "rms": float(np.sqrt(np.mean(x**2))),
        "annulus_2_to_10_rms": float(np.sqrt(np.mean(x[annulus]**2))) if np.any(annulus) else None,
        "cell_0": float(x[0]),
        "cell_1": float(x[1]) if x.size > 1 else None,
        "fixed_radius_samples": fixed,
    }


def scalar_rhs_term_budget(state, coupling: float) -> dict[str, Any]:
    """Decompose the pinned scalar RHS into signed, reconstructing term profiles.

    This is an implementation/source crosswalk: it records what the pinned RHS
    computes and does not accept the conditional phi--DM exchange premise.
    """
    coupling = float(coupling)
    if not math.isfinite(coupling):
        raise ValueError("scalar RHS term budget coupling must be finite")

    grid = state.grid
    geom = state.geometry
    fields = state.scalars
    r = finite("scalar RHS budget radii", grid.centers)
    alpha = finite("scalar RHS budget lapse", geom.alpha)
    beta = finite("scalar RHS budget shift", geom.beta)
    K = finite("scalar RHS budget trace curvature", geom.K)
    a = finite("scalar RHS budget a", geom.a)
    b = finite("scalar RHS budget b", geom.b)
    X = finite("scalar RHS budget X", geom.X)
    invr = finite("scalar RHS budget inverse radial metric", X**2 / a)
    ap = d1(grid, a, +1)
    bp = d1(grid, b, +1)
    Xp = d1(grid, X, +1)
    alphap = d1(grid, alpha, +1)
    metric_bracket = finite(
        "scalar RHS budget radial Laplacian bracket",
        0.5 * ap / a - bp / b + Xp / X - 2.0 / r,
    )
    metric_slice = metric_slice_from_q(grid, geom)
    rho_dm = finite(
        "scalar RHS budget dark-matter density",
        dm_density(metric_slice, state.matter),
    )
    rhs = scalar_rhs_arrays(
        grid, geom, fields, beta_dm=coupling, rho_dm=rho_dm
    )

    definitions = (
        ("S", fields.S, fields.PS, fields.S, "S", "PS"),
        ("D", fields.D, fields.PD, -fields.D, "D", "PD"),
        (
            "phi", fields.phi, fields.Pi,
            cosmos_potential_prime(fields.phi), "phi", "Pi"
        ),
    )
    sector_reports: dict[str, Any] = {}
    all_finite = True
    all_reconstructed = True
    closure_rtol = 5.0e-12
    closure_atol = 5.0e-18

    def pack_terms(terms: dict[str, np.ndarray], name: str) -> dict[str, Any]:
        checked = {
            key: finite(f"{name} term {key}", value)
            for key, value in terms.items()
        }
        return {
            "raw_signed_arrays": {key: value.tolist() for key, value in checked.items()},
            "summaries": {key: summarize(value, r) for key, value in checked.items()},
        }

    for sector, field, momentum, potential_prime, field_key, momentum_key in definitions:
        f = finite(f"scalar RHS {sector} field", field)
        p = finite(f"scalar RHS {sector} momentum", momentum)
        fp = d1(grid, f, +1)
        fpp = d2(grid, f, +1)
        pp = d1(grid, p, +1)
        zero = np.zeros_like(r, dtype=np.float64)

        field_terms = {
            "lapse_momentum": alpha * p,
            "shift_advection": beta * fp,
        }
        momentum_terms = {
            "shift_advection": beta * pp,
            "lapse_laplacian_second_derivative": alpha * invr * fpp,
            "lapse_metric_gradient": -alpha * invr * fp * metric_bracket,
            "trace_curvature": alpha * K * p,
            "lapse_gradient": invr * alphap * fp,
            "potential_derivative": -alpha * finite(
                f"{sector} potential derivative", potential_prime
            ),
            "conditional_dm_exchange": coupling * rho_dm if sector == "phi" else zero,
        }

        field_sum = np.zeros_like(r, dtype=np.float64)
        for value in field_terms.values():
            field_sum = field_sum + value
        momentum_sum = np.zeros_like(r, dtype=np.float64)
        for value in momentum_terms.values():
            momentum_sum = momentum_sum + value

        actual_field_rate = finite(
            f"{sector} implementation field rate", getattr(rhs, field_key)
        )
        actual_momentum_rate = finite(
            f"{sector} implementation momentum rate", getattr(rhs, momentum_key)
        )
        field_gap = finite(
            f"{sector} reconstructed-minus-implementation field rate",
            field_sum - actual_field_rate,
        )
        momentum_gap = finite(
            f"{sector} reconstructed-minus-implementation momentum rate",
            momentum_sum - actual_momentum_rate,
        )
        field_passed = bool(np.allclose(
            field_sum, actual_field_rate, rtol=closure_rtol, atol=closure_atol
        ))
        momentum_passed = bool(np.allclose(
            momentum_sum, actual_momentum_rate, rtol=closure_rtol, atol=closure_atol
        ))
        all_reconstructed = all_reconstructed and field_passed and momentum_passed
        sector_reports[sector] = {
            "field_rate_terms": pack_terms(field_terms, f"{sector} field rate"),
            "momentum_rate_terms": pack_terms(momentum_terms, f"{sector} momentum rate"),
            "implementation_rates": {
                "field_rate": actual_field_rate.tolist(),
                "momentum_rate": actual_momentum_rate.tolist(),
            },
            "reconstructed_rates": {
                "field_rate": field_sum.tolist(),
                "momentum_rate": momentum_sum.tolist(),
            },
            "reconstruction": {
                "field_rate_gap": field_gap.tolist(),
                "momentum_rate_gap": momentum_gap.tolist(),
                "field_rate_gap_summary": summarize(field_gap, r),
                "momentum_rate_gap_summary": summarize(momentum_gap, r),
                "field_rate_passed": field_passed,
                "momentum_rate_passed": momentum_passed,
                "rtol": closure_rtol,
                "atol": closure_atol,
            },
        }

    all_finite = all(
        np.all(np.isfinite(np.asarray(values, dtype=np.float64)))
        for item in sector_reports.values()
        for block in (
            item["field_rate_terms"]["raw_signed_arrays"],
            item["momentum_rate_terms"]["raw_signed_arrays"],
            item["implementation_rates"],
            item["reconstructed_rates"],
            item["reconstruction"],
        )
        for values in block.values()
        if isinstance(values, (list, tuple, np.ndarray))
    )
    return {
        "status": "scalar_rhs_term_crosswalk_only",
        "reported": True,
        "coupling_beta_dm": coupling,
        "state_sha256": state_digest(state),
        "conditional_exchange_model_law_accepted": False,
        "rho_dm_raw_signed_array": rho_dm.tolist(),
        "rho_dm_summary": summarize(rho_dm, r),
        "sectors": sector_reports,
        "all_profiles_finite": bool(all_finite),
        "reconstruction_passed": bool(all_reconstructed),
        "reconstruction_rtol": closure_rtol,
        "reconstruction_atol": closure_atol,
        "term_convention": {
            "field_rate": "alpha*p + beta^r*D_r(f)",
            "momentum_rate": (
                "beta^r*D_r(p) + alpha*(X^2/a)*D2(f) "
                "- alpha*(X^2/a)*D_r(f)*(0.5*D_r(a)/a-D_r(b)/b+D_r(X)/X-2/r) "
                "+ alpha*K*p + (X^2/a)*D_r(alpha)*D_r(f) - alpha*V'(f)"
            ),
            "conditional_dm_exchange": (
                "phi only: beta_DM*rho_DM as implemented in the pinned RHS; "
                "reported conditionally, not treated as accepted physics"
            ),
        },
    }


def scalar_rhs_term_coupling_increment(
    uncoupled: dict[str, Any], coupled: dict[str, Any], radii: np.ndarray
) -> dict[str, Any]:
    """Report coupled-minus-control term changes on the identical base state."""
    r = finite("scalar coupling-increment radii", radii)
    permitted_change = ("phi", "momentum_rate_terms", "conditional_dm_exchange")
    by_sector: dict[str, Any] = {}
    only_declared_exchange_changed = True
    for sector in ("S", "D", "phi"):
        by_sector[sector] = {}
        for category in ("field_rate_terms", "momentum_rate_terms"):
            unc_terms = uncoupled["sectors"][sector][category]["raw_signed_arrays"]
            cou_terms = coupled["sectors"][sector][category]["raw_signed_arrays"]
            if set(unc_terms) != set(cou_terms):
                raise RuntimeError(f"scalar RHS term names changed for {sector}/{category}")
            for term in unc_terms:
                delta = finite(
                    f"{sector} {category} {term} coupled-minus-control",
                    np.asarray(cou_terms[term], dtype=np.float64)
                    - np.asarray(unc_terms[term], dtype=np.float64),
                )
                allowed = (sector, category, term) == permitted_change
                if not allowed and np.any(delta != 0.0):
                    only_declared_exchange_changed = False
                by_sector[sector][f"{category}.{term}"] = {
                    "raw_signed_array": delta.tolist(),
                    "summary": summarize(delta, r),
                    "expected_to_change_only_under_declared_rhs_split": bool(allowed),
                }

    coupling_delta = float(coupled["coupling_beta_dm"] - uncoupled["coupling_beta_dm"])
    expected = finite(
        "expected conditional phi RHS coupling increment",
        coupling_delta * np.asarray(coupled["rho_dm_raw_signed_array"], dtype=np.float64),
    )
    measured = np.asarray(
        by_sector["phi"]["momentum_rate_terms.conditional_dm_exchange"]["raw_signed_array"],
        dtype=np.float64,
    )
    expected_match = bool(np.allclose(measured, expected, rtol=5.0e-12, atol=5.0e-18))
    return {
        "reported": True,
        "status": "same_base_scalar_rhs_term_difference",
        "uncoupled_state_sha256": uncoupled["state_sha256"],
        "coupled_state_sha256": coupled["state_sha256"],
        "same_base_state": uncoupled["state_sha256"] == coupled["state_sha256"],
        "conditional_exchange_model_law_accepted": False,
        "phi_exchange_increment_matches_declared_rhs_term": expected_match,
        "only_declared_phi_exchange_term_changed": bool(only_declared_exchange_changed),
        "sectors": by_sector,
        "conditional_phi_exchange_increment": {
            "formula": "(beta_DM,coupled - beta_DM,uncoupled)*rho_DM",
            "raw_signed_array": measured.tolist(),
            "expected_raw_signed_array": expected.tolist(),
            "summary": summarize(measured, r),
        },
    }


def _divergence_tangent_from_objects(
    state, base: dict[str, Any], plus: dict[str, Any], minus: dict[str, Any], epsilon: float
) -> dict[str, Any]:
    """Apply the existing mixed-tensor divergence to an already-projected state triplet."""
    if not math.isfinite(float(epsilon)) or float(epsilon) <= 0.0:
        raise ValueError("epsilon must be finite and positive")
    if np.any(np.asarray(state.geometry.beta) != 0.0):
        raise ValueError("this derivation currently assumes beta^r=0 at the frozen slice")
    inv = 1.0 / (2.0 * float(epsilon))
    d1_fn = lambda values, parity: d1(state.grid, values, parity)
    time_density_t = (plus["density_t"] - minus["density_t"]) * inv
    time_density_r = (plus["density_r"] - minus["density_r"]) * inv
    dt_gtt = (plus["gtt"] - minus["gtt"]) * inv
    dt_grr = (plus["gamma_rr"] - minus["gamma_rr"]) * inv
    dt_gth = (plus["gamma_thth"] - minus["gamma_thth"]) * inv
    div_t, div_r = covariant_divergence_from_profiles(
        state.grid.centers, d1_fn,
        base["alpha"], base["gamma_rr"], base["gamma_thth"],
        base["rho"], base["j"], base["pr"], base["pt"],
        time_density_t, time_density_r, dt_gtt, dt_grr, dt_gth,
    )
    div_t_fv, div_r_fv = covariant_divergence_fv_from_profiles(
        state.grid.centers, d1_fn,
        base["alpha"], base["gamma_rr"], base["gamma_thth"],
        base["rho"], base["j"], base["pr"], base["pt"],
        time_density_t, time_density_r, dt_gtt, dt_grr, dt_gth,
    )
    return {
        "epsilon": float(epsilon),
        "div_T_t": div_t, "div_T_r": div_r,
        "div_T_t_fv": div_t_fv, "div_T_r_fv": div_r_fv,
        "div_T_t_summary": summarize(div_t, state.grid.centers),
        "div_T_r_summary": summarize(div_r, state.grid.centers),
        "div_T_t_fv_summary": summarize(div_t_fv, state.grid.centers),
        "div_T_r_fv_summary": summarize(div_r_fv, state.grid.centers),
        "fv_minus_cell_centered_div_T_t_summary": summarize(div_t_fv - div_t, state.grid.centers),
        "fv_minus_cell_centered_div_T_r_summary": summarize(div_r_fv - div_r, state.grid.centers),
        "fv_face_flux_method": "four-point cubic cell-center-to-face reconstruction; zero inner spherical flux; outer cubic extrapolation",
        "projection_crosscheck": base.get("projection_crosscheck", {}),
        "projection_crosscheck_passed": bool(base.get("projection_crosscheck_passed", True)),
        "time_derivative_witnesses": {
            "d_t(sqrt(-g) T^t_t)": summarize(time_density_t, state.grid.centers),
            "d_t(sqrt(-g) T^t_r)": summarize(time_density_r, state.grid.centers),
            "d_t(g_tt)": summarize(dt_gtt, state.grid.centers),
            "d_t(gamma_rr)": summarize(dt_grr, state.grid.centers),
            "d_t(gamma_thth)": summarize(dt_gth, state.grid.centers),
        },
    }


def _divergence_tangent_from_projection_triplet(
    state, base_projection, plus_projection, minus_projection, epsilon: float
) -> dict[str, Any]:
    base = _covariant_objects_from_components(state, base_projection)
    plus = _covariant_objects_from_components(state, plus_projection)
    minus = _covariant_objects_from_components(state, minus_projection)
    return _divergence_tangent_from_objects(state, base, plus, minus, epsilon)


def divergence_tangent_from_states(
    state, plus_state, minus_state, epsilon: float
) -> dict[str, Any]:
    """Measure total-T divergence using a supplied plus/minus state pair."""
    if not math.isfinite(float(epsilon)) or float(epsilon) <= 0.0:
        raise ValueError("epsilon must be finite and positive")
    assert_supported_baseline_source_scope(plus_state)
    assert_supported_baseline_source_scope(minus_state)
    projections = (
        independent_total_projections(state),
        independent_total_projections(plus_state),
        independent_total_projections(minus_state),
    )
    return _divergence_tangent_from_projection_triplet(
        state, projections[0], projections[1], projections[2], float(epsilon)
    )


def sector_resolved_divergence_tangent_from_states(
    state, plus_state, minus_state, epsilon: float
) -> dict[str, Any]:
    """Resolve the same total-T directional witness into additive fluid/S/D/phi sectors.

    Each sector shares the exact frozen base, plus/minus objects, epsilon, and
    metric-derivative witnesses of the total projection. This is an accounting
    decomposition, not acceptance of a physical exchange law.
    """
    if not math.isfinite(float(epsilon)) or float(epsilon) <= 0.0:
        raise ValueError("epsilon must be finite and positive")
    assert_supported_baseline_source_scope(plus_state)
    assert_supported_baseline_source_scope(minus_state)
    if np.any(np.asarray(state.geometry.beta) != 0.0):
        raise ValueError("sector decomposition currently assumes beta^r=0 at the frozen slice")
    projections = (
        independent_total_projections(state),
        independent_total_projections(plus_state),
        independent_total_projections(minus_state),
    )
    object_triplet = tuple(
        _covariant_objects_from_components(s, p)
        for s, p in zip((state, plus_state, minus_state), projections)
    )
    base_all, plus_all, minus_all = object_triplet
    total = _divergence_tangent_from_objects(
        state, base_all, plus_all, minus_all, float(epsilon)
    )
    names = ("dark_matter", "baryons", "radiation", "S", "D", "phi")
    projection_maps = []
    for projection in projections:
        combined = dict(projection.get("fluid_parts", {}))
        combined.update(projection.get("scalar_parts", {}))
        projection_maps.append(combined)
    if any(set(piece) != set(names) for piece in projection_maps):
        raise RuntimeError("sector projection set differs from the declared DM/baryon/radiation/S/D/phi scope")
    inv = 1.0 / (2.0 * float(epsilon))
    d1_fn = lambda values, parity: d1(state.grid, values, parity)
    dt_gtt = (plus_all["gtt"] - minus_all["gtt"]) * inv
    dt_grr = (plus_all["gamma_rr"] - minus_all["gamma_rr"]) * inv
    dt_gth = (plus_all["gamma_thth"] - minus_all["gamma_thth"]) * inv
    channels = ("div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv")
    sector_results: dict[str, Any] = {}
    sector_dt_density_t_sum = np.zeros_like(base_all["rho"], dtype=np.float64)
    sector_dt_density_r_sum = np.zeros_like(base_all["rho"], dtype=np.float64)
    for name in names:
        b_obj = _covariant_objects_from_components(state, projection_maps[0][name])
        p_obj = _covariant_objects_from_components(plus_state, projection_maps[1][name])
        m_obj = _covariant_objects_from_components(minus_state, projection_maps[2][name])
        dt_rho_t = (p_obj["density_t"] - m_obj["density_t"]) * inv
        dt_rho_r = (p_obj["density_r"] - m_obj["density_r"]) * inv
        sector_dt_density_t_sum += dt_rho_t
        sector_dt_density_r_sum += dt_rho_r
        div_t, div_r = covariant_divergence_from_profiles(
            state.grid.centers, d1_fn,
            b_obj["alpha"], b_obj["gamma_rr"], b_obj["gamma_thth"],
            b_obj["rho"], b_obj["j"], b_obj["pr"], b_obj["pt"],
            dt_rho_t, dt_rho_r, dt_gtt, dt_grr, dt_gth,
        )
        div_t_fv, div_r_fv = covariant_divergence_fv_from_profiles(
            state.grid.centers, d1_fn,
            b_obj["alpha"], b_obj["gamma_rr"], b_obj["gamma_thth"],
            b_obj["rho"], b_obj["j"], b_obj["pr"], b_obj["pt"],
            dt_rho_t, dt_rho_r, dt_gtt, dt_grr, dt_gth,
        )
        arrays = {
            "div_T_t": finite(f"{name} sector div_T_t", div_t),
            "div_T_r": finite(f"{name} sector div_T_r", div_r),
            "div_T_t_fv": finite(f"{name} sector div_T_t FV", div_t_fv),
            "div_T_r_fv": finite(f"{name} sector div_T_r FV", div_r_fv),
        }
        sector_results[name] = {
            "raw_signed_arrays": {k: v.tolist() for k, v in arrays.items()},
            "summaries": {k: summarize(v, state.grid.centers) for k, v in arrays.items()},
        }
    # Isolate the directional density term before interpreting any full-divergence
    # sector-sum residual. This determines whether the failing time component is
    # already present in (d_t sqrt(-g) T^t_nu), before spatial flux/source terms.
    total_dt_density_t = (plus_all["density_t"] - minus_all["density_t"]) * inv
    total_dt_density_r = (plus_all["density_r"] - minus_all["density_r"]) * inv
    time_density_gap_t = finite(
        "time-density sector sum minus total for nu=t",
        sector_dt_density_t_sum - total_dt_density_t,
    )
    time_density_gap_r = finite(
        "time-density sector sum minus total for nu=r",
        sector_dt_density_r_sum - total_dt_density_r,
    )
    predicted_pointwise_gap_t = finite(
        "predicted pointwise divergence gap from nu=t time-density term",
        time_density_gap_t / base_all["sqrt_minus_g"],
    )
    predicted_pointwise_gap_r = finite(
        "predicted pointwise divergence gap from nu=r time-density term",
        time_density_gap_r / base_all["sqrt_minus_g"],
    )
    zero = np.zeros_like(state.grid.centers, dtype=np.float64)
    d1_fn_for_gap = lambda values, parity: d1(state.grid, values, parity)
    predicted_fv_gap_t, predicted_fv_gap_r = covariant_divergence_fv_from_profiles(
        state.grid.centers, d1_fn_for_gap,
        base_all["alpha"], base_all["gamma_rr"], base_all["gamma_thth"],
        zero, zero, zero, zero,
        time_density_gap_t, time_density_gap_r, zero, zero, zero,
    )
    time_density_gap_summary = {
        "term_definition": (
            "sum of per-sector finite-difference density rates minus the total-projection "
            "finite-difference density rate, before adding radial flux and metric source terms"
        ),
        "raw_signed_arrays": {
            "time_density_gap_t": time_density_gap_t.tolist(),
            "time_density_gap_r": time_density_gap_r.tolist(),
            "predicted_div_T_t_gap_pointwise": predicted_pointwise_gap_t.tolist(),
            "predicted_div_T_r_gap_pointwise": predicted_pointwise_gap_r.tolist(),
            "predicted_div_T_t_gap_fv": predicted_fv_gap_t.tolist(),
            "predicted_div_T_r_gap_fv": predicted_fv_gap_r.tolist(),
        },
        "summaries": {
            "time_density_gap_t": summarize(time_density_gap_t, state.grid.centers),
            "time_density_gap_r": summarize(time_density_gap_r, state.grid.centers),
            "predicted_div_T_t_gap_pointwise": summarize(predicted_pointwise_gap_t, state.grid.centers),
            "predicted_div_T_r_gap_pointwise": summarize(predicted_pointwise_gap_r, state.grid.centers),
            "predicted_div_T_t_gap_fv": summarize(predicted_fv_gap_t, state.grid.centers),
            "predicted_div_T_r_gap_fv": summarize(predicted_fv_gap_r, state.grid.centers),
        },
    }

    atol, rtol = 5.0e-18, 5.0e-10
    component_predicted_gaps = {
        "div_T_t": predicted_pointwise_gap_t,
        "div_T_r": predicted_pointwise_gap_r,
        "div_T_t_fv": predicted_fv_gap_t,
        "div_T_r_fv": predicted_fv_gap_r,
    }
    recon_channels: dict[str, Any] = {}
    for key in channels:
        sector_sum = np.sum(np.stack([
            np.asarray(sector_results[n]["raw_signed_arrays"][key], dtype=np.float64)
            for n in names
        ], axis=0), axis=0)
        gap = finite(f"{key} additive sector reconstruction gap",
                     sector_sum - np.asarray(total[key], dtype=np.float64))
        passed = bool(np.allclose(sector_sum, np.asarray(total[key], dtype=np.float64),
                                  rtol=rtol, atol=atol))
        predicted_gap = finite(
            f"{key} predicted time-density component gap", component_predicted_gaps[key]
        )
        unexplained_gap = finite(
            f"{key} residual after time-density component accounting", gap - predicted_gap
        )
        component_corrected_passed = bool(np.allclose(
            unexplained_gap, np.zeros_like(unexplained_gap), rtol=rtol, atol=atol
        ))
        recon_channels[key] = {
            "passed": passed,
            "sum_of_sectors": sector_sum.tolist(),
            "sum_minus_total": gap.tolist(),
            "gap_summary": summarize(gap, state.grid.centers),
            "max_abs_gap": float(np.max(np.abs(gap))),
            "time_density_component_predicted_gap": predicted_gap.tolist(),
            "time_density_component_predicted_gap_summary": summarize(predicted_gap, state.grid.centers),
            "unexplained_gap_after_time_density_accounting": unexplained_gap.tolist(),
            "unexplained_gap_summary": summarize(unexplained_gap, state.grid.centers),
            "max_abs_unexplained_gap": float(np.max(np.abs(unexplained_gap))),
            "component_corrected_passed": component_corrected_passed,
            "rtol": rtol, "atol": atol,
        }
    decomposition = {
        "status": "sector_resolved_diagnostic_only",
        "sector_names": list(names),
        "sectors": sector_results,
        "reconstruction": {
            "sector_names": list(names),
            "all_sectors_reconstruct_total": all(x["passed"] for x in recon_channels.values()),
            "all_sectors_reconstruct_total_after_time_density_accounting": all(
                x["component_corrected_passed"] for x in recon_channels.values()
            ),
            "time_density_component_accounting_passed": all(
                x["component_corrected_passed"] for x in recon_channels.values()
            ),
            "raw_sum_max_gap": max(x["max_abs_gap"] for x in recon_channels.values()),
            "unexplained_max_gap": max(x["max_abs_unexplained_gap"] for x in recon_channels.values()),
            "channels": recon_channels,
        },
        "time_density_reconstruction": time_density_gap_summary,
        "state_pairing": {
            "base_state_sha256": state_digest(state),
            "plus_state_sha256": state_digest(plus_state),
            "minus_state_sha256": state_digest(minus_state),
            "epsilon": float(epsilon),
            "all_sector_arrays_finite": all(
                np.all(np.isfinite(np.asarray(values, dtype=np.float64)))
                for sector in sector_results.values()
                for values in sector["raw_signed_arrays"].values()
            ),
        },
    }
    return {"total": total, "sector_decomposition": decomposition}


def divergence_tangent(kernel, state, rates, epsilon: float) -> dict[str, Any]:
    """Compatibility wrapper that constructs a two-sided pair, then measures it."""
    plus_state = shifted_state(kernel, state, rates, +epsilon)
    minus_state = shifted_state(kernel, state, rates, -epsilon)
    return divergence_tangent_from_states(state, plus_state, minus_state, epsilon)

def conditional_rest_current_source(
    state, rates: dict[str, Any], coupling: float
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Derive the rest-current source implied by the declared covector law.

    For pressureless dust T^mu_nu=rho*u^mu*u_nu and
    nabla_mu T^mu_nu=Q_nu=beta_dm*rho*partial_nu(phi), contraction with u^nu
    gives nabla_mu(rho*u^mu)=-beta_dm*rho*u^nu*partial_nu(phi).

    The Valencia conserved rest variable here is sqrt(gamma)*rho*W. Hence its
    conservative RHS source is
      S_rest = -sqrt(gamma)*beta_dm*rho*W*
               (partial_t(phi)+(alpha*v^r-beta^r)*partial_r(phi)).
    This is a conditional model-law discriminator only. It does not modify
    production rates or infer that this source law has been accepted.
    """
    grid = state.grid
    geom = state.geometry
    metric_slice = metric_slice_from_q(grid, geom)
    metrics = _metric_arrays(metric_slice)
    prims = primitives(metrics, state.matter.dark_matter, Species.DARK_MATTER)
    rho = finite("conditional rest-source rho", np.asarray([q.rho for q in prims]))
    W = finite("conditional rest-source Lorentz factor", np.asarray([q.lorentz() for q in prims]))
    v_r = finite("conditional rest-source radial velocity", np.asarray([q.v_r for q in prims]))
    sqrt_gamma = finite("conditional rest-source sqrt gamma", np.asarray([m.sqrt_gamma for m in metrics]))
    alpha = finite("conditional rest-source lapse", geom.alpha)
    beta_r = finite("conditional rest-source shift", geom.beta)
    dphi_t = finite("conditional rest-source d_t(phi)", rates["scalars"]["phi"])
    dphi_r = d1(grid, state.scalars.phi, +1)
    transport_phi_rate = finite(
        "conditional rest-source convective phi derivative",
        dphi_t + (alpha * v_r - beta_r) * dphi_r,
    )
    source = finite(
        "conditional conserved rest-current source",
        -sqrt_gamma * float(coupling) * rho * W * transport_phi_rate,
    )
    parts = {
        "rho": rho,
        "W": W,
        "v_r": v_r,
        "sqrt_gamma": sqrt_gamma,
        "dphi_t": dphi_t,
        "dphi_r": dphi_r,
        "transport_phi_rate": transport_phi_rate,
    }
    return source, parts


def copy_rates_with_rest_source(rates: dict[str, Any], source: np.ndarray) -> dict[str, Any]:
    """Copy only RHS containers needed to add a test-only DM rest source."""
    copied = dict(rates)
    copied["matter"] = {
        species: {
            name: np.array(value, dtype=np.float64, copy=True)
            for name, value in channels.items()
        }
        for species, channels in rates["matter"].items()
    }
    copied["matter"]["dark_matter"]["rest"] += np.asarray(source, dtype=np.float64)
    return copied


def audit_case(kernel, state, coupling: float, epsilons: tuple[float, float]) -> dict[str, Any]:
    original_digest = state_digest(state)
    old_beta = production_kernel_module.BETA_DM
    try:
        production_kernel_module.BETA_DM = float(coupling)
        rates = rhs_direction(kernel, state)
    finally:
        production_kernel_module.BETA_DM = old_beta
    tangent_rows = []
    sector_resolved_epsilon_rows = []
    for epsilon in epsilons:
        # Materialize one shared +/- state pair and reuse it for the total and
        # six-sector divergence witnesses; no parallel or independently built
        # state pair is introduced for the scalar budget.
        plus_state = shifted_state(kernel, state, rates, +float(epsilon))
        minus_state = shifted_state(kernel, state, rates, -float(epsilon))
        paired = sector_resolved_divergence_tangent_from_states(
            state, plus_state, minus_state, float(epsilon)
        )
        tangent_rows.append(paired["total"])
        sector_resolved_epsilon_rows.append(paired["sector_decomposition"])
    scalar_term_budget = scalar_rhs_term_budget(state, coupling)
    rest_source, rest_source_parts = conditional_rest_current_source(state, rates, coupling)
    rates_with_conditional_source = copy_rates_with_rest_source(rates, rest_source)
    conditional_source_tangent_rows = [
        divergence_tangent(kernel, state, rates_with_conditional_source, eps)
        for eps in epsilons
    ]

    # Independent companion predicted by the pinned-source exchange audit:
    # Pi_dot currently includes beta_dm*rho_dm, while the covariant wave-source
    # conversion under the same Q_nu premise requires alpha*beta_dm*rho_dm.
    # This is inserted only into another copied RHS direction for comparison.
    metric_slice = metric_slice_from_q(state.grid, state.geometry)
    rho_dm = finite("conditional Pi-lapse rho_dm", dm_density(metric_slice, state.matter))
    scalar_pi_lapse_correction = finite(
        "conditional Pi-dot lapse correction",
        (np.asarray(state.geometry.alpha) - 1.0) * float(coupling) * rho_dm,
    )
    rates_with_both_conditional_sources = copy_rates_with_rest_source(rates, rest_source)
    rates_with_both_conditional_sources["scalars"] = {
        name: np.array(value, dtype=np.float64, copy=True)
        for name, value in rates["scalars"].items()
    }
    rates_with_both_conditional_sources["scalars"]["Pi"] += scalar_pi_lapse_correction
    joint_source_tangent_rows = [
        divergence_tangent(kernel, state, rates_with_both_conditional_sources, eps)
        for eps in epsilons
    ]

    if state_digest(state) != original_digest:
        raise RuntimeError("frozen input state changed during covariant-divergence audit")
    all_tangent_rows = tangent_rows + conditional_source_tangent_rows + joint_source_tangent_rows
    projections_ok = all(row["projection_crosscheck_passed"] for row in all_tangent_rows)
    finite_ok = all(
        np.all(np.isfinite(row["div_T_t"]))
        and np.all(np.isfinite(row["div_T_r"]))
        and np.all(np.isfinite(row["div_T_t_fv"]))
        and np.all(np.isfinite(row["div_T_r_fv"]))
        for row in all_tangent_rows
    ) and bool(
        np.all(np.isfinite(rest_source))
        and np.all(np.isfinite(scalar_pi_lapse_correction))
        and scalar_term_budget["all_profiles_finite"]
        and all(row["state_pairing"]["all_sector_arrays_finite"]
                for row in sector_resolved_epsilon_rows)
    )
    r = np.asarray(state.grid.centers, dtype=np.float64)
    annulus = (r >= 2.0) & (r <= 10.0)
    return {
        "coupling_beta_dm": float(coupling),
        "initial_state_sha256": original_digest,
        "zero_shift_at_frozen_slice": bool(np.all(np.asarray(state.geometry.beta) == 0.0)),
        "projection_crosscheck_passed": bool(projections_ok),
        "all_profiles_finite": bool(finite_ok),
        "scalar_rhs_term_budget": scalar_term_budget,
        "sector_resolved_epsilon_rows": sector_resolved_epsilon_rows,
        "epsilon_rows": [
            {key: value for key, value in row.items() if key not in ("div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv")}
            for row in tangent_rows
        ],
        "conditional_rest_source_epsilon_rows": [
            {key: value for key, value in row.items() if key not in ("div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv")}
            for row in conditional_source_tangent_rows
        ],
        "joint_conditional_source_epsilon_rows": [
            {key: value for key, value in row.items() if key not in ("div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv")}
            for row in joint_source_tangent_rows
        ],
        "conditional_pi_dot_lapse_correction": {
            "formula": "(alpha-1)*beta_dm*rho_dm",
            "profile": summarize(scalar_pi_lapse_correction, np.asarray(state.grid.centers, dtype=np.float64)),
            "model_law_accepted": False,
        },
        "conditional_rest_current_source": {
            "premise": "nabla_mu T^mu_nu = Q_nu = beta_dm*rho*partial_nu(phi) for pressureless dark matter",
            "conservative_source_formula": "-sqrt(gamma)*beta_dm*rho*W*(partial_t(phi)+(alpha*v^r-beta^r)*partial_r(phi))",
            "profile": summarize(rest_source, r),
            "profile_parts": {name: summarize(value, r) for name, value in rest_source_parts.items()},
            "model_law_accepted": False,
        },
        "raw_arrays": {
            f"epsilon_{i}_div_T_t": row["div_T_t"].tolist()
            for i, row in enumerate(tangent_rows)
        } | {
            f"epsilon_{i}_div_T_r": row["div_T_r"].tolist()
            for i, row in enumerate(tangent_rows)
        } | {
            f"epsilon_{i}_conditional_rest_source_div_T_t": row["div_T_t"].tolist()
            for i, row in enumerate(conditional_source_tangent_rows)
        } | {
            f"epsilon_{i}_conditional_rest_source_div_T_r": row["div_T_r"].tolist()
            for i, row in enumerate(conditional_source_tangent_rows)
        } | {
            f"epsilon_{i}_joint_conditional_sources_div_T_t": row["div_T_t"].tolist()
            for i, row in enumerate(joint_source_tangent_rows)
        } | {
            f"epsilon_{i}_joint_conditional_sources_div_T_r": row["div_T_r"].tolist()
            for i, row in enumerate(joint_source_tangent_rows)
        } | {
            f"epsilon_{i}_div_T_t_fv": row["div_T_t_fv"].tolist()
            for i, row in enumerate(tangent_rows)
        } | {
            f"epsilon_{i}_div_T_r_fv": row["div_T_r_fv"].tolist()
            for i, row in enumerate(tangent_rows)
        } | {
            f"epsilon_{i}_div_T_t_fv_minus_cell_centered": (
                row["div_T_t_fv"] - row["div_T_t"]
            ).tolist()
            for i, row in enumerate(tangent_rows)
        } | {
            f"epsilon_{i}_div_T_r_fv_minus_cell_centered": (
                row["div_T_r_fv"] - row["div_T_r"]
            ).tolist()
            for i, row in enumerate(tangent_rows)
        } | {
            f"epsilon_{i}_conditional_rest_source_div_T_t_fv": row["div_T_t_fv"].tolist()
            for i, row in enumerate(conditional_source_tangent_rows)
        } | {
            f"epsilon_{i}_conditional_rest_source_div_T_r_fv": row["div_T_r_fv"].tolist()
            for i, row in enumerate(conditional_source_tangent_rows)
        } | {
            f"epsilon_{i}_conditional_rest_source_div_T_r_fv_minus_cell_centered": (
                row["div_T_r_fv"] - row["div_T_r"]
            ).tolist()
            for i, row in enumerate(conditional_source_tangent_rows)
        } | {
            f"epsilon_{i}_joint_conditional_sources_div_T_t_fv": row["div_T_t_fv"].tolist()
            for i, row in enumerate(joint_source_tangent_rows)
        } | {
            f"epsilon_{i}_joint_conditional_sources_div_T_r_fv": row["div_T_r_fv"].tolist()
            for i, row in enumerate(joint_source_tangent_rows)
        } | {
            f"epsilon_{i}_joint_conditional_sources_div_T_r_fv_minus_cell_centered": (
                row["div_T_r_fv"] - row["div_T_r"]
            ).tolist()
            for i, row in enumerate(joint_source_tangent_rows)
        } | {
            "conditional_rest_current_source": rest_source.tolist(),
            "conditional_pi_dot_lapse_correction": scalar_pi_lapse_correction.tolist(),
        },
        "epsilon_sensitivity": {
            "div_T_t_annulus_rms": float(np.sqrt(np.mean((
                tangent_rows[0]["div_T_t"] - tangent_rows[1]["div_T_t"]
            )[annulus]**2))),
            "div_T_r_annulus_rms": float(np.sqrt(np.mean((
                tangent_rows[0]["div_T_r"] - tangent_rows[1]["div_T_r"]
            )[annulus]**2))),
        },
        "conditional_rest_source_epsilon_sensitivity": {
            "div_T_t_annulus_rms": float(np.sqrt(np.mean((
                conditional_source_tangent_rows[0]["div_T_t"] - conditional_source_tangent_rows[1]["div_T_t"]
            )[annulus]**2))),
            "div_T_r_annulus_rms": float(np.sqrt(np.mean((
                conditional_source_tangent_rows[0]["div_T_r"] - conditional_source_tangent_rows[1]["div_T_r"]
            )[annulus]**2))),
        },
        "joint_conditional_source_epsilon_sensitivity": {
            "div_T_t_annulus_rms": float(np.sqrt(np.mean((
                joint_source_tangent_rows[0]["div_T_t"] - joint_source_tangent_rows[1]["div_T_t"]
            )[annulus]**2))),
            "div_T_r_annulus_rms": float(np.sqrt(np.mean((
                joint_source_tangent_rows[0]["div_T_r"] - joint_source_tangent_rows[1]["div_T_r"]
            )[annulus]**2))),
        },
    }


def sector_divergence_coupling_increment(
    uncoupled_rows: list[dict[str, Any]],
    coupled_rows: list[dict[str, Any]],
    radii: np.ndarray,
) -> dict[str, Any]:
    """Compare sector divergence profiles across coupling arms and same epsilons."""
    r = finite("sector divergence coupling-increment radii", radii)
    if len(uncoupled_rows) != len(coupled_rows):
        raise ValueError("coupled/control sector epsilon row counts differ")
    names = ("dark_matter", "baryons", "radiation", "S", "D", "phi")
    channels = ("div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv")
    epsilon_rows = []
    for unc_row, cou_row in zip(uncoupled_rows, coupled_rows):
        if float(unc_row["state_pairing"]["epsilon"]) != float(cou_row["state_pairing"]["epsilon"]):
            raise ValueError("coupled/control sector epsilon values differ")
        by_sector = {}
        for sector in names:
            unc_raw = unc_row["sectors"][sector]["raw_signed_arrays"]
            cou_raw = cou_row["sectors"][sector]["raw_signed_arrays"]
            by_sector[sector] = {}
            for channel in channels:
                delta = finite(
                    f"{sector} {channel} coupled-minus-control divergence",
                    np.asarray(cou_raw[channel], dtype=np.float64)
                    - np.asarray(unc_raw[channel], dtype=np.float64),
                )
                by_sector[sector][channel] = {
                    "raw_signed_array": delta.tolist(),
                    "summary": summarize(delta, r),
                }
        epsilon_rows.append({
            "epsilon": float(unc_row["state_pairing"]["epsilon"]),
            "uncoupled_state_pairing": unc_row["state_pairing"],
            "coupled_state_pairing": cou_row["state_pairing"],
            "sectors": by_sector,
        })
    return {
        "reported": True,
        "status": "same_state_pair_sector_coupling_difference",
        "sector_names": list(names),
        "channels": list(channels),
        "epsilon_rows": epsilon_rows,
        "production_steps_advanced": 0,
        "scientific_admission": False,
    }


def audit_resolution(n: int, r_max: float, epsilons: tuple[float, float]) -> dict[str, Any]:
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=n, r_max=r_max, include_radiation=True)
    assert_supported_baseline_source_scope(state)
    if np.any(np.asarray(state.geometry.beta) != 0.0):
        raise RuntimeError("expected the standard initialized gauge to have zero shift")
    input_hash_before = state_digest(state)
    uncoupled = audit_case(kernel, state, 0.0, epsilons)
    coupled = audit_case(kernel, state, BETA_DM_COUPLED, epsilons)
    if state_digest(state) != input_hash_before:
        raise RuntimeError("audit changed original initialized state")
    annulus = (state.grid.centers >= 2.0) & (state.grid.centers <= 10.0)
    # The difference isolates effects of the declared coupling on the same
    # frozen physical slice; it is not assumed to be a priori zero.
    delta_t = (
        np.asarray(coupled["raw_arrays"]["epsilon_0_div_T_t"])
        - np.asarray(uncoupled["raw_arrays"]["epsilon_0_div_T_t"])
    )
    delta_r = (
        np.asarray(coupled["raw_arrays"]["epsilon_0_div_T_r"])
        - np.asarray(uncoupled["raw_arrays"]["epsilon_0_div_T_r"])
    )
    delta_t_with_rest_source = (
        np.asarray(coupled["raw_arrays"]["epsilon_0_conditional_rest_source_div_T_t"])
        - np.asarray(uncoupled["raw_arrays"]["epsilon_0_conditional_rest_source_div_T_t"])
    )
    delta_r_with_rest_source = (
        np.asarray(coupled["raw_arrays"]["epsilon_0_conditional_rest_source_div_T_r"])
        - np.asarray(uncoupled["raw_arrays"]["epsilon_0_conditional_rest_source_div_T_r"])
    )
    delta_t_rest_correction = finite(
        "change in coupled D_t residual from conditional DM rest source",
        delta_t_with_rest_source - delta_t,
    )
    annulus_before_rest_rms = float(np.sqrt(np.mean(delta_t[annulus] ** 2)))
    annulus_after_rest_rms = float(np.sqrt(np.mean(delta_t_with_rest_source[annulus] ** 2)))
    # The source-pinned scalar RHS adds beta_dm*rho_dm to Pi_dot, but the
    # covariant Klein-Gordon conversion requires alpha*beta_dm*rho_dm.
    # The prior source audit gives the conservative local energy-rate defect
    # (1-alpha)*beta_dm*rho_dm*Pi/KAPPA. Since the diagnostic measures the
    # covariant mixed component nabla_mu T^mu_t (where T^t_t=-rho), the
    # predicted divergence residual has the opposite sign, (alpha-1)*... .
    # Test the signed profile, not only its RMS; this remains conditional.
    metric_slice = metric_slice_from_q(state.grid, state.geometry)
    raw_rho_dm = finite("lapse-discriminator DM density", dm_density(metric_slice, state.matter))
    predicted_lapse_defect = finite(
        "derived covariant D_t lapse-factor defect",
        (np.asarray(state.geometry.alpha) - 1.0)
        * BETA_DM_COUPLED * raw_rho_dm * np.asarray(state.scalars.Pi) / KAPPA,
    )
    lapse_defect_residual = finite(
        "measured coupling increment minus derived lapse defect",
        delta_t - predicted_lapse_defect,
    )
    joint_residual_same_sign = finite(
        "post-rest-source D_t increment minus lapse defect",
        delta_t_with_rest_source - predicted_lapse_defect,
    )
    joint_residual_opposite_sign = finite(
        "post-rest-source D_t increment plus lapse defect",
        delta_t_with_rest_source + predicted_lapse_defect,
    )
    predicted_lapse_rms = float(np.sqrt(np.mean(predicted_lapse_defect[annulus] ** 2)))
    joint_same_rms = float(np.sqrt(np.mean(joint_residual_same_sign[annulus] ** 2)))
    joint_opposite_rms = float(np.sqrt(np.mean(joint_residual_opposite_sign[annulus] ** 2)))
    lapse_defect_scale = max(
        float(np.sqrt(np.mean(predicted_lapse_defect[annulus] ** 2))),
        np.finfo(float).tiny,
    )
    scalar_term_increment = scalar_rhs_term_coupling_increment(
        uncoupled["scalar_rhs_term_budget"],
        coupled["scalar_rhs_term_budget"],
        state.grid.centers,
    )
    sector_divergence_increment = sector_divergence_coupling_increment(
        uncoupled["sector_resolved_epsilon_rows"],
        coupled["sector_resolved_epsilon_rows"],
        state.grid.centers,
    )
    return {
        "resolution": int(n),
        "r_max": float(r_max),
        "dr": float(state.grid.dr),
        "initial_state_sha256": input_hash_before,
        "uncoupled_control": uncoupled,
        "coupled_baseline": coupled,
        "coupling_increment_difference": {
            "div_T_t": summarize(delta_t, state.grid.centers),
            "div_T_r": summarize(delta_r, state.grid.centers),
        },
        "scalar_rhs_term_coupling_increment": scalar_term_increment,
        "sector_divergence_coupling_increment": sector_divergence_increment,
        "conditional_rest_current_source_discriminator": {
            "premise": "For pressureless dust with nabla_mu T^mu_nu=beta_dm*rho*partial_nu(phi), contraction with u^nu gives nabla_mu(rho*u^mu)=-beta_dm*rho*u^nu*partial_nu(phi).",
            "production_rest_rhs_has_interaction_source": False,
            "tested_change": "Add the derived source only to a copied RHS direction for the directional tangent; no production file or state is changed.",
            "conditional_coupling_increment_Dt_before": summarize(delta_t, state.grid.centers),
            "conditional_coupling_increment_Dt_after": summarize(delta_t_with_rest_source, state.grid.centers),
            "conditional_coupling_increment_Dr_before": summarize(delta_r, state.grid.centers),
            "conditional_coupling_increment_Dr_after": summarize(delta_r_with_rest_source, state.grid.centers),
            "change_due_to_test_source_Dt": summarize(delta_t_rest_correction, state.grid.centers),
            "fixed_annulus_2_to_10": {
                "dt_rms_before": annulus_before_rest_rms,
                "dt_rms_after": annulus_after_rest_rms,
                "dt_rms_reduction_factor": float(annulus_before_rest_rms / max(annulus_after_rest_rms, np.finfo(float).tiny)),
                "conditional_source_reduces_dt_residual": bool(annulus_after_rest_rms < annulus_before_rest_rms),
                "dr_rms_before": float(np.sqrt(np.mean(delta_r[annulus] ** 2))),
                "dr_rms_after": float(np.sqrt(np.mean(delta_r_with_rest_source[annulus] ** 2))),
            },
            "joint_conditional_lapse_remainder": {
                "derived_formula": "(alpha-1)*beta_dm*rho_dm*Pi/KAPPA (opposite sign to conservative energy-rate defect because T^t_t=-rho)",
                "predicted_signed_profile": summarize(predicted_lapse_defect, state.grid.centers),
                "post_rest_source_coupling_increment": summarize(delta_t_with_rest_source, state.grid.centers),
                "rms_residual_after_subtracting_prediction": joint_same_rms,
                "rms_residual_after_adding_prediction": joint_opposite_rms,
                "predicted_lapse_rms": predicted_lapse_rms,
                "post_rest_source_d_t_rms": annulus_after_rest_rms,
                "same_sign_closure_factor": float(predicted_lapse_rms / max(joint_same_rms, np.finfo(float).tiny)),
                "opposite_sign_closure_factor": float(predicted_lapse_rms / max(joint_opposite_rms, np.finfo(float).tiny)),
                "candidate_sign_match": bool(joint_same_rms < joint_opposite_rms),
                "closure_admitted": False,
            },
            "model_law_accepted": False,
            "match_admitted": False,
        },
        "lapse_factor_exchange_discriminator": {
            "derived_formula": "(alpha-1)*beta_dm*rho_dm*Pi/KAPPA (opposite sign to conservative energy-rate defect because T^t_t=-rho)",
            "predicted_Dt": summarize(predicted_lapse_defect, state.grid.centers),
            "measured_coupling_increment_Dt_before_rest_source": summarize(delta_t, state.grid.centers),
            "measured_coupling_increment_Dt_after_rest_source": summarize(delta_t_with_rest_source, state.grid.centers),
            "residual_after_derived_lapse_defect_before_rest_source": summarize(lapse_defect_residual, state.grid.centers),
            "fixed_annulus_2_to_10": {
                "predicted_rms": lapse_defect_scale,
                "measured_rms_before_rest_source": float(np.sqrt(np.mean(delta_t[annulus] ** 2))),
                "measured_rms_after_rest_source": annulus_after_rest_rms,
                "residual_rms_before_rest_source": float(np.sqrt(np.mean(lapse_defect_residual[annulus] ** 2))),
                "relative_rms_residual_to_prediction_before_rest_source": float(
                    np.sqrt(np.mean(lapse_defect_residual[annulus] ** 2)) / lapse_defect_scale
                ),
            },
            "match_admitted": False,
        },
        "raw_arrays": {
            "coupling_increment_Dt_before_rest_source": delta_t.tolist(),
            "coupling_increment_Dt_after_rest_source": delta_t_with_rest_source.tolist(),
            "derived_lapse_factor_defect_Dt": predicted_lapse_defect.tolist(),
            "joint_conditional_residual_same_sign": joint_residual_same_sign.tolist(),
            "joint_conditional_residual_opposite_sign": joint_residual_opposite_sign.tolist(),
        },
        "limits": [
            "Finite-difference time tangent includes the CMC lapse re-solve and the full frozen semi-discrete RHS direction.",
            "The FV comparator uses parity-aware total-T face reconstruction and Simpson-integrated time/source densities and spherical cell measure; it is diagnostic-only and is not the production per-species HLL flux.",
            "Any coupled residual is reported, not judged from a conditional interaction law; it requires independent source-law/model review.",
            "The audit has beta^r=0 at the frozen initial slice; nonzero-shift extension is a separate control.",
        ],
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    provenance = verify_provenance()
    if len(args.epsilons) != 2 or any(e <= 0.0 for e in args.epsilons):
        raise ValueError("provide exactly two positive directional epsilons")
    if not args.resolutions or any(n < 32 for n in args.resolutions):
        raise ValueError("resolutions must each be >=32")
    if not math.isfinite(args.r_max) or args.r_max <= 0.0:
        raise ValueError("r_max must be positive and finite")
    per_resolution = [
        audit_resolution(n, float(args.r_max), (float(args.epsilons[0]), float(args.epsilons[1])))
        for n in args.resolutions
    ]
    return {
        "schema": "issue4_independent_covariant_total_stress_energy_divergence_v1",
        "status": "completed_diagnostic_only",
        "scientific_admission": False,
        "source_provenance": provenance,
        "configuration": {
            "resolutions": list(args.resolutions),
            "r_max": float(args.r_max),
            "directional_epsilons": list(args.epsilons),
            "production_steps_advanced": 0,
            "nonzero_shift_at_frozen_slice": False,
        },
        "resolution_results": per_resolution,
        "gate_result": {
            "source_pins_verified": True,
            "independent_total_projection_matches_canonical": all(
                result["uncoupled_control"]["projection_crosscheck_passed"]
                and result["coupled_baseline"]["projection_crosscheck_passed"]
                for result in per_resolution
            ),
            "all_profiles_finite": all(
                result["uncoupled_control"]["all_profiles_finite"]
                and result["coupled_baseline"]["all_profiles_finite"]
                for result in per_resolution
            ),
            "scalar_rhs_term_budget_reported": all(
                result["uncoupled_control"]["scalar_rhs_term_budget"]["reported"]
                and result["coupled_baseline"]["scalar_rhs_term_budget"]["reported"]
                for result in per_resolution
            ),
            "scalar_rhs_term_reconstruction_passed": all(
                result["uncoupled_control"]["scalar_rhs_term_budget"]["reconstruction_passed"]
                and result["coupled_baseline"]["scalar_rhs_term_budget"]["reconstruction_passed"]
                for result in per_resolution
            ),
            "sector_resolved_scalar_divergence_reported": all(
                len(result["uncoupled_control"]["sector_resolved_epsilon_rows"]) == len(args.epsilons)
                and len(result["coupled_baseline"]["sector_resolved_epsilon_rows"]) == len(args.epsilons)
                and result["scalar_rhs_term_coupling_increment"]["reported"]
                and result["sector_divergence_coupling_increment"]["reported"]
                for result in per_resolution
            ),
            "scalar_rhs_coupling_increment_consistent": all(
                result["scalar_rhs_term_coupling_increment"]["same_base_state"]
                and result["scalar_rhs_term_coupling_increment"]["phi_exchange_increment_matches_declared_rhs_term"]
                and result["scalar_rhs_term_coupling_increment"]["only_declared_phi_exchange_term_changed"]
                for result in per_resolution
            ),
            "input_states_unchanged": True,
            "independent_total_stress_energy_divergence_closed": False,
            "scientific_admission": False,
        },
        "interpretation": (
            "This first independent covariant-divergence witness separates a zero-coupling "
            "control from the current beta_DM=-0.04 model on identical initialized states. "
            "It independently reconstructs total T projections and uses the mixed-tensor "
            "covariant divergence identity. Residuals are measured diagnostic evidence, "
            "not automatically attributed to physics: first separate cell-centered versus "
            "finite-volume spatial discretization, epsilon sensitivity, gauge response, and "
            "the conditional exchange-law premise. No production equation is modified."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[60, 80, 120])
    parser.add_argument("--r-max", type=float, default=24.0)
    parser.add_argument("--epsilons", nargs=2, type=float, default=[2.0e-5, 1.0e-5])
    parser.add_argument("--output", type=Path, default=Path(
        "results/issue4-independent-covariant-total-T-divergence/report.json"
    ))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        report = run(args)
        code = 0
    except Exception as exc:
        report = {
            "schema": "issue4_independent_covariant_total_stress_energy_divergence_v1",
            "status": "failed_closed",
            "scientific_admission": False,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        code = 1
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report.get("status"),
        "gate_result": report.get("gate_result"),
        "configuration": report.get("configuration"),
        "scientific_admission": False,
        "artifact": str(args.output),
    }, sort_keys=True), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
