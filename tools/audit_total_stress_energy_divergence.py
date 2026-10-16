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
from engine.scalar_system import KAPPA, ScalarFields, cosmos_potential  # noqa: E402
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


def spherical_flux_face_values(values: np.ndarray, dr: float) -> np.ndarray:
    """Reconstruct spherical integrated fluxes from cell centers to faces.

    Four cell-center values define a cubic interpolant at each face. The
    inner face is exactly zero because a regular spherical integrated radial
    flux carries the r^2 area factor. The outer face is a one-sided cubic
    extrapolation. This is a diagnostic reconstruction of the already
    projected total flux, not the production per-species HLL face flux.
    """
    flux = finite("finite-volume cell-centered flux", values)
    if flux.ndim != 1 or flux.size < 4:
        raise ValueError("finite-volume face reconstruction requires n>=4")
    if not math.isfinite(float(dr)) or float(dr) <= 0.0:
        raise ValueError("finite-volume dr must be positive and finite")
    n = flux.size
    faces = np.zeros(n + 1, dtype=np.float64)
    # Uniform spherical-cell centers are x_i=(i+1/2)dr; face j is x=j*dr.
    for j in range(1, n + 1):
        start = min(max(j - 2, 0), n - 4)
        nodes = np.arange(start, start + 4, dtype=np.float64) + 0.5
        target = float(j)
        weights = np.ones(4, dtype=np.float64)
        for k in range(4):
            for ell in range(4):
                if ell != k:
                    weights[k] *= (target - nodes[ell]) / (nodes[k] - nodes[ell])
        faces[j] = float(weights @ flux[start:start + 4])
    return finite("finite-volume face flux", faces)


def finite_volume_radial_flux_derivative(values: np.ndarray, dr: float) -> np.ndarray:
    """Cellwise derivative from telescoping reconstructed spherical face fluxes."""
    faces = spherical_flux_face_values(values, dr)
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
    """Covariant divergence with finite-volume radial flux differences.

    The time tangent and metric-source terms match the cell-centered
    derivative observer. Only the radial flux derivative changes: the same
    projected total-T flux is reconstructed to spherical faces with a cubic
    interpolant and differenced conservatively across each cell. This isolates
    spatial flux-discretization sensitivity; it is not the production HLL flux.
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
    sqrt_minus_g = alpha * sqrt_gamma
    if np.any(alpha <= 0.0) or np.any(grr <= 0.0) or np.any(gth <= 0.0):
        raise ValueError("finite-volume divergence requires positive lapse and spatial metric")

    Ttt = E / (alpha**2)
    Trr = pr / grr
    Tthth = pt / gth
    flux_t = -(alpha**2) * sqrt_gamma * (j / grr)
    flux_r = alpha * sqrt_gamma * pr
    radial_flux_t = finite_volume_radial_flux_derivative(flux_t, dr)
    radial_flux_r = finite_volume_radial_flux_derivative(flux_r, dr)

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
    denom = finite("finite-volume sqrt(-g) denominator", sqrt_minus_g)
    div_t = finite(
        "finite-volume covariant time-component divergence",
        (finite("FV dt sqrt(-g) T^t_t", dt_density_t) + radial_flux_t) / denom + source_t,
    )
    div_r = finite(
        "finite-volume covariant radial-component divergence",
        (finite("FV dt sqrt(-g) T^t_r", dt_density_r) + radial_flux_r) / denom + source_r,
    )
    return div_t, div_r


def covariant_objects(state) -> dict[str, np.ndarray]:
    projection = independent_total_projections(state)
    r = finite("covariant grid radii", state.grid.centers)
    alpha = finite("covariant lapse", state.geometry.alpha)
    gamma_rr = finite("gamma_rr", state.geometry.a / state.geometry.X**2)
    gamma_thth = finite("gamma_thth", state.geometry.b * r**2 / state.geometry.X**2)
    sqrt_gamma = np.sqrt(gamma_rr) * gamma_thth
    sqrt_minus_g = alpha * sqrt_gamma
    j_up = projection["j"] / gamma_rr
    return {
        "rho": projection["rho"],
        "pr": projection["pr"],
        "pt": projection["pt"],
        "j": projection["j"],
        "alpha": alpha,
        "gamma_rr": gamma_rr,
        "gamma_thth": gamma_thth,
        "sqrt_gamma": sqrt_gamma,
        "sqrt_minus_g": sqrt_minus_g,
        "density_t": -sqrt_minus_g * projection["rho"],  # sqrt(-g) T^t_t
        "density_r": sqrt_gamma * projection["j"],       # sqrt(-g) T^t_r
        "flux_t": -(alpha**2) * sqrt_gamma * j_up,       # sqrt(-g) T^r_t
        "flux_r": alpha * sqrt_gamma * projection["pr"], # sqrt(-g) T^r_r
        "gtt": -(alpha**2),
        "projection_crosscheck": projection["projection_crosscheck"],
        "projection_crosscheck_passed": projection["projection_crosscheck_passed"],
    }


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


def divergence_tangent(kernel, state, rates, epsilon: float) -> dict[str, Any]:
    plus_state = shifted_state(kernel, state, rates, +epsilon)
    minus_state = shifted_state(kernel, state, rates, -epsilon)
    assert_supported_baseline_source_scope(plus_state)
    assert_supported_baseline_source_scope(minus_state)
    if np.any(np.asarray(state.geometry.beta) != 0.0):
        raise ValueError("this derivation currently assumes beta^r=0 at the frozen slice")

    base = covariant_objects(state)
    plus = covariant_objects(plus_state)
    minus = covariant_objects(minus_state)
    inv = 1.0 / (2.0 * epsilon)
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
        "div_T_t": div_t,
        "div_T_r": div_r,
        "div_T_t_fv": div_t_fv,
        "div_T_r_fv": div_r_fv,
        "div_T_t_summary": summarize(div_t, state.grid.centers),
        "div_T_r_summary": summarize(div_r, state.grid.centers),
        "div_T_t_fv_summary": summarize(div_t_fv, state.grid.centers),
        "div_T_r_fv_summary": summarize(div_r_fv, state.grid.centers),
        "fv_minus_cell_centered_div_T_t_summary": summarize(div_t_fv - div_t, state.grid.centers),
        "fv_minus_cell_centered_div_T_r_summary": summarize(div_r_fv - div_r, state.grid.centers),
        "fv_face_flux_method": "four-point cubic cell-center-to-face reconstruction; zero inner spherical flux; outer cubic extrapolation",
        "projection_crosscheck": base["projection_crosscheck"],
        "projection_crosscheck_passed": base["projection_crosscheck_passed"],
        "time_derivative_witnesses": {
            "d_t(sqrt(-g) T^t_t)": summarize(time_density_t, state.grid.centers),
            "d_t(sqrt(-g) T^t_r)": summarize(time_density_r, state.grid.centers),
            "d_t(g_tt)": summarize(dt_gtt, state.grid.centers),
            "d_t(gamma_rr)": summarize(dt_grr, state.grid.centers),
            "d_t(gamma_thth)": summarize(dt_gth, state.grid.centers),
        },
    }


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
    tangent_rows = [divergence_tangent(kernel, state, rates, eps) for eps in epsilons]
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
    )
    r = np.asarray(state.grid.centers, dtype=np.float64)
    annulus = (r >= 2.0) & (r <= 10.0)
    return {
        "coupling_beta_dm": float(coupling),
        "initial_state_sha256": original_digest,
        "zero_shift_at_frozen_slice": bool(np.all(np.asarray(state.geometry.beta) == 0.0)),
        "projection_crosscheck_passed": bool(projections_ok),
        "all_profiles_finite": bool(finite_ok),
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
            "Spatial divergence uses the exact pinned parity-aware D1 but cell-centred projected T; finite-volume fluid reconstruction is not identical to that point-sampled covariant operator.",
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
