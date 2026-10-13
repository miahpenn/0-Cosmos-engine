#!/usr/bin/env python3
"""Source-pinned, research-only joined matter/EM/quasi-local ledger.

This tool observes a bounded, paired q/m=0 versus q/m=0.25 radial-Maxwell
campaign through the current production kernel. It does not alter evolution,
choose a physical carrier law, or assert a global conserved energy in a
dynamical spacetime. Full signed stage profiles go to a JSON artifact; compact
numbers go to stdout and the GitHub Actions step summary.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

# Running a file under tools/ places that directory first on sys.path. Add the
# repository root explicitly so "python tools/audit_em_total_energy_budget.py"
# resolves the in-tree engine package consistently in CI and locally.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from engine import matter_rhs as matter_rhs_module
from engine import v55_pirk_adapter as adapter
from engine.campaign import CampaignConfig
from engine.production_kernel import V55ProductionKernel
from engine.radial_maxwell import RadialMaxwellState
from engine.run_production import run_campaign
from engine.stress_energy import assemble_total_stress_energy
from engine.worldtube import areal_radius, areal_radius_derivative, invariant_chi, misner_sharp


FOUR_PI = 4.0 * math.pi
RESOLUTIONS = (24, 48, 96)
R_MAX = 24.0
FINAL_TIME = 0.03
CFL = 0.01
CHECKPOINT_INTERVAL = 0.01
TEMPORAL_REFINEMENT_RESOLUTION = 48
TIME_STEP_CFL_VALUES = (0.01, 0.005, 0.0025)
DENSITY_PROFILE = ((0.0, 2.0e-5), (R_MAX, 2.0e-5))
CHARGE_TO_MASS_ARMS = (0.0, 0.25)
FIXED_AREAL_RADIUS = 10.0

# Blob IDs from main 5c51bf509504365df394e6cca96bf1d117f4614d. The audit is
# deliberately tied to these operators; update the pin only after a new source
# review and a new exact-head result, never silently on a moving main branch.
SOURCE_PINS = {
    "engine/campaign.py": "119942d9debf63b99997cda96c541754e7e8aded",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/cosmology_observables.py": "4d6f39308685bae5e5bad2e17baa9c4043b67b4e",
    "engine/invariant_diagnostics.py": "4bbcc6174c978b07c4f2175d1dcfeccb2d0c12b1",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/matter_system.py": "05e23d2368120ea8a1a6aad3fa438625d5380750",
    "engine/production_contract.py": "ef1606eb5a01e0f5906a42de77afd35d3b9fffed",
    "engine/production_kernel.py": "4e0bdccc2ac614c48c34521663a6d53b404b6747",
    "engine/radial_maxwell.py": "81269f1d1334c46c4c2f44ba816d792722c0a38e",
    "engine/run_production.py": "726f3fd09d46142862574d348aa0755f83bc0612",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/v55_matter.py": "acc8f8ceae28d82511c9a5ad9d2ca670e09fb0f7",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
    "engine/v55_pirk_adapter.py": "eeee0492d80889d70e52b0946563eee5054c0f11",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
}
VENDOR_GITLINK = "vendor/bb-palatini-unified-r0"
VENDOR_PIN = "d6052d605673ce9d82cdc99b0df79855fa2ea215"


def verify_source_pins(root: Path) -> dict:
    """Fail closed unless every production source blob and vendor gitlink match."""
    observed = {}
    for path, expected in SOURCE_PINS.items():
        completed = subprocess.run(
            ["git", "rev-parse", f"HEAD:{path}"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
        actual = completed.stdout.strip()
        observed[path] = {"expected": expected, "observed": actual, "match": actual == expected}
        if actual != expected:
            raise RuntimeError(
                f"source pin mismatch for {path}: expected {expected}, got {actual}"
            )
    completed = subprocess.run(
        ["git", "rev-parse", f"HEAD:{VENDOR_GITLINK}"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    actual_vendor = completed.stdout.strip()
    observed[VENDOR_GITLINK] = {
        "expected": VENDOR_PIN,
        "observed": actual_vendor,
        "match": actual_vendor == VENDOR_PIN,
    }
    if actual_vendor != VENDOR_PIN:
        raise RuntimeError(
            f"vendor gitlink mismatch: expected {VENDOR_PIN}, got {actual_vendor}"
        )
    # Pull-request Actions jobs normally execute on a synthetic merge commit.
    # Keep the execution tree distinct from the PR's source branch head.
    execution_checkout_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    github_sha = os.environ.get("GITHUB_SHA", "").strip()
    if github_sha and github_sha != execution_checkout_sha:
        raise RuntimeError(
            "GitHub event SHA does not match the actual execution checkout: "
            f"GITHUB_SHA={github_sha}, checkout={execution_checkout_sha}"
        )
    pr_source_head_sha = os.environ.get("PR_SOURCE_HEAD_SHA", "").strip() or None
    source_head_sha = pr_source_head_sha or github_sha or execution_checkout_sha
    if len(source_head_sha) != 40 or any(ch not in "0123456789abcdefABCDEF" for ch in source_head_sha):
        raise RuntimeError(f"invalid source-head SHA for report provenance: {source_head_sha!r}")
    return {
        "execution_checkout_sha": execution_checkout_sha,
        "source_head_sha": source_head_sha,
        "source_head_kind": (
            "pull_request.head.sha" if pr_source_head_sha
            else ("GITHUB_SHA" if github_sha else "execution-checkout fallback (local run)")
        ),
        "production_blob_pins": observed,
    }


def field_energy_functional(dr: float, centers, geometry_arrays: dict, charge_faces) -> float:
    """Discrete total radial electric energy used by the live EM exchange test."""
    r = np.asarray(centers, dtype=float)
    q = np.asarray(charge_faces, dtype=float)
    if q.shape != (r.size + 1,):
        raise ValueError("charge_faces must have one more entry than centers")
    weight = (
        np.asarray(geometry_arrays["X"], dtype=float)
        * np.sqrt(np.asarray(geometry_arrays["a"], dtype=float))
        / (np.square(r) * np.asarray(geometry_arrays["b"], dtype=float))
    )
    q_cell = 0.5 * (q[:-1] + q[1:])
    return 0.5 * float(dr) * float(np.sum(weight * q_cell**2))


def interpolate_at_areal_radius(areal_radii, profile, target: float) -> float:
    """Interpolate a cell profile at one fixed areal radius, fail-closed."""
    R = np.asarray(areal_radii, dtype=float)
    values = np.asarray(profile, dtype=float)
    if R.ndim != 1 or values.shape != R.shape or R.size < 2:
        raise ValueError("areal radii and profile must be matching one-dimensional arrays")
    if not np.all(np.isfinite(R)) or not np.all(np.isfinite(values)):
        raise ValueError("areal-radius interpolation inputs must be finite")
    if not np.all(np.diff(R) > 0.0):
        raise ValueError("areal radius must be strictly increasing for tube interpolation")
    if not float(R[0]) <= float(target) <= float(R[-1]):
        raise ValueError(f"target areal radius {target} is outside the sampled domain")
    return float(np.interp(float(target), R, values))


def _array(value):
    return np.asarray(value, dtype=float).copy()


def _reconstruct_energy_rhs_with_roundoff_contract(observed, named_terms, *, label):
    """Reconstruct a cellwise RHS and fail closed outside a roundoff-scaled bound."""
    observed = np.asarray(observed, dtype=float)
    if observed.ndim != 1 or observed.size == 0:
        raise RuntimeError(f"{label}: observed RHS must be a nonempty one-dimensional array")
    if not np.all(np.isfinite(observed)):
        raise RuntimeError(f"{label}: observed RHS contains non-finite values")

    expected = np.zeros_like(observed)
    scale = np.abs(observed).copy()
    for term_name, raw_term in named_terms:
        term = np.asarray(raw_term, dtype=float)
        if term.shape != observed.shape:
            raise RuntimeError(
                f"{label}: {term_name} shape {term.shape} does not match RHS shape {observed.shape}"
            )
        if not np.all(np.isfinite(term)):
            raise RuntimeError(f"{label}: {term_name} contains non-finite values")
        expected += term
        scale += np.abs(term)

    if not np.all(np.isfinite(expected)) or not np.all(np.isfinite(scale)):
        raise RuntimeError(f"{label}: non-finite RHS reconstruction or operand scale")
    residual = observed - expected
    bound = 64.0 * np.finfo(np.float64).eps * np.maximum(
        scale, np.finfo(np.float64).tiny
    )
    if not np.all(np.isfinite(residual)) or not np.all(np.abs(residual) <= bound):
        ratio = np.abs(residual) / bound
        worst = int(np.argmax(ratio))
        raise RuntimeError(
            f"{label}: RHS reconstruction exceeded scale-aware roundoff bound "
            f"at cell {worst}: residual={residual[worst]:.17g}, "
            f"bound={bound[worst]:.17g}, ratio={ratio[worst]:.17g}"
        )
    return expected, residual, float(np.max(np.abs(residual) / bound))


def _state_energy_integrals(state) -> dict[str, float]:
    dr = float(state.grid.dr)
    out = {}
    for name in ("dark_matter", "baryons", "radiation", "charged_dust"):
        species = getattr(state.matter, name, None)
        if species is not None:
            out[name] = FOUR_PI * dr * float(np.sum(species.energy_t))
    return out


def _em_energy(state) -> float:
    if state.electromagnetic is None:
        return 0.0
    geom = state.geometry
    return field_energy_functional(
        state.grid.dr,
        state.grid.centers,
        {"X": geom.X, "a": geom.a, "b": geom.b},
        state.electromagnetic.enclosed_charge_faces,
    )


def _stress_projection(kernel, state):
    projection = kernel._electromagnetic_projection(
        state.grid, state.geometry, state.electromagnetic
    )
    if projection is None:
        return assemble_total_stress_energy(
            state.grid, state.geometry, state.scalars, state.matter
        )
    return assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter,
        electromagnetic_projection=projection,
    )


def _misner_sharp_chain_rule_rate_profiles(
    R, Rdot, Rr, Rdot_r, gamma_inv, gamma_inv_dot, Ktheta, Ktheta_dot
) -> dict[str, np.ndarray]:
    """Differentiate M_MS=R(1-chi)/2 at fixed coordinate from semidiscrete RHS.

    With chi = gamma^{rr} (D_r R)^2 - (R Ktheta)^2, the derivative uses
    the supplied live-stage rates for R, D_r R, gamma^{rr}, and Ktheta.
    The fixed-areal-radius tube sweep is applied separately by the caller.
    """
    raw = {
        "R": R,
        "Rdot": Rdot,
        "Rr": Rr,
        "Rdot_r": Rdot_r,
        "gamma_inv": gamma_inv,
        "gamma_inv_dot": gamma_inv_dot,
        "Ktheta": Ktheta,
        "Ktheta_dot": Ktheta_dot,
    }
    arrays = {name: np.asarray(value, dtype=float) for name, value in raw.items()}
    first = arrays["R"]
    if first.ndim != 1 or first.size == 0:
        raise RuntimeError("Misner-Sharp chain-rule R must be a nonempty 1-D array")
    for name, value in arrays.items():
        if value.shape != first.shape:
            raise RuntimeError(
                f"Misner-Sharp chain-rule {name} shape {value.shape} does not match {first.shape}"
            )
        if not np.all(np.isfinite(value)):
            raise RuntimeError(f"Misner-Sharp chain-rule {name} contains non-finite values")

    R = arrays["R"]
    Rdot = arrays["Rdot"]
    Rr = arrays["Rr"]
    Rdot_r = arrays["Rdot_r"]
    gamma_inv = arrays["gamma_inv"]
    gamma_inv_dot = arrays["gamma_inv_dot"]
    Ktheta = arrays["Ktheta"]
    Ktheta_dot = arrays["Ktheta_dot"]

    normal_dR = R * Ktheta
    normal_dR_dot = Rdot * Ktheta + R * Ktheta_dot
    chi = gamma_inv * Rr**2 - normal_dR**2
    chi_dot = (
        gamma_inv_dot * Rr**2
        + 2.0 * gamma_inv * Rr * Rdot_r
        - 2.0 * normal_dR * normal_dR_dot
    )
    mass = 0.5 * R * (1.0 - chi)
    mass_dot_coordinate = 0.5 * (
        Rdot * (1.0 - chi) - R * chi_dot
    )
    outputs = {
        "chi": chi,
        "chi_dot": chi_dot,
        "mass": mass,
        "mass_dot_coordinate": mass_dot_coordinate,
    }
    if not all(np.all(np.isfinite(value)) for value in outputs.values()):
        raise RuntimeError("Misner-Sharp chain-rule evaluation produced non-finite output")
    return outputs


def _fixed_target_interpolant_rate(
    R, Rdot, mass, mass_dot_coordinate, *, target: float
) -> dict:
    """Differentiate the exact piecewise-linear M(R,t) interpolant at fixed R.

    This is the derivative of the same moving-bracket interpolant used for the
    accepted-state tube observable, not a separately interpolated continuum
    chain rule. At an exact interior knot, use the right-hand bracket (and the
    left-hand bracket at the upper domain edge); the one-sided convention is
    reported through the active bracket and theta.
    """
    raw = {"R": R, "Rdot": Rdot, "mass": mass, "mass_dot_coordinate": mass_dot_coordinate}
    arrays = {name: np.asarray(value, dtype=float) for name, value in raw.items()}
    first = arrays["R"]
    if first.ndim != 1 or first.size < 2:
        raise RuntimeError("fixed-R interpolant rate R must be a 1-D array with at least two cells")
    for name, value in arrays.items():
        if value.shape != first.shape:
            raise RuntimeError(
                f"fixed-R interpolant rate {name} shape {value.shape} does not match {first.shape}"
            )
        if not np.all(np.isfinite(value)):
            raise RuntimeError(f"fixed-R interpolant rate {name} contains non-finite values")
    if not np.all(np.diff(first) > 0.0):
        raise RuntimeError("fixed-R interpolant rate R must be strictly increasing")
    target = float(target)
    if not math.isfinite(target):
        raise RuntimeError("fixed-R interpolant target must be finite")
    if target < float(first[0]) or target > float(first[-1]):
        raise RuntimeError("fixed-R interpolant target is outside the sampled domain")

    # Match np.interp's linear segment and retain the actual moving bracket.
    left = int(np.searchsorted(first, target, side="right") - 1)
    left = max(0, min(left, first.size - 2))
    right = left + 1
    span = float(first[right] - first[left])
    theta = float((target - first[left]) / span)
    blend = lambda values: float((1.0 - theta) * values[left] + theta * values[right])
    interp_mass = blend(arrays["mass"])
    interp_Rdot = blend(arrays["Rdot"])
    interp_mass_dot = blend(arrays["mass_dot_coordinate"])
    piecewise_gradient = float((arrays["mass"][right] - arrays["mass"][left]) / span)
    fixed_R_sweep = -interp_Rdot * piecewise_gradient
    fixed_R_rate = interp_mass_dot + fixed_R_sweep
    scalars = (
        interp_mass, interp_Rdot, interp_mass_dot, piecewise_gradient,
        fixed_R_sweep, fixed_R_rate, theta,
    )
    if not all(math.isfinite(value) for value in scalars):
        raise RuntimeError("fixed-R interpolant rate produced a non-finite result")
    return {
        "left_index": left,
        "right_index": right,
        "theta": theta,
        "interpolated_mass": interp_mass,
        "interpolated_Rdot": interp_Rdot,
        "interpolated_mass_dot_coordinate": interp_mass_dot,
        "piecewise_linear_mass_gradient": piecewise_gradient,
        "fixed_R_sweep": fixed_R_sweep,
        "mass_dot_fixed_R": fixed_R_rate,
    }


def _worldtube_sample(kernel, state, *, geometry_stage_capture=None) -> dict:
    """Evaluate native MS and EM-inclusive stress/flux/work on R=10."""
    grid, geom = state.grid, state.geometry
    R = np.asarray(areal_radius(grid, geom), dtype=float)
    if not np.all(np.diff(R) > 0.0):
        raise ValueError("native areal radius is not strictly increasing")
    Rr = np.asarray(areal_radius_derivative(grid, geom), dtype=float)
    chi = np.asarray(invariant_chi(grid, geom), dtype=float)
    mass = np.asarray(misner_sharp(grid, geom), dtype=float)
    Mr = np.asarray(grid.cell_derivative_fourth(mass, parity=-1), dtype=float)

    if geometry_stage_capture is None:
        _, _, moving = adapter.vendor_modules()
        explicit = moving.moving_puncture_explicit_rhs(grid, geom)
    else:
        explicit = geometry_stage_capture["explicit"]
    a_dot = np.asarray(explicit["a"], dtype=float)
    b_dot = np.asarray(explicit["b"], dtype=float)
    X_dot = np.asarray(explicit["X"], dtype=float)
    Rdot = R * (0.5 * b_dot / geom.b - X_dot / geom.X)

    total = _stress_projection(kernel, state)
    gamma_rr_inv = geom.X**2 / geom.a
    at = lambda profile: interpolate_at_areal_radius(R, profile, FIXED_AREAL_RADIUS)

    Rr_tube = at(Rr)
    Rdot_tube = at(Rdot)
    alpha_tube = at(geom.alpha)
    j_tube = at(total.j)
    pr_tube = at(total.pr)
    rho_tube = at(total.rho)
    gamma_inv_tube = at(gamma_rr_inv)
    mass_tube = at(mass)
    Mr_tube = at(Mr)
    tube_speed = -Rdot_tube / Rr_tube

    # Signs follow the existing production-kernel diagnostic; these quantities
    # are reported as an audited decomposition, not assumed to be a proven law.
    stress_flux = (
        -FOUR_PI * FIXED_AREAL_RADIUS**2
        * alpha_tube * gamma_inv_tube * j_tube * Rr_tube
    )
    pressure_work = (
        -FOUR_PI * FIXED_AREAL_RADIUS**2 * pr_tube * Rdot_tube
    )
    tube_sweep = Mr_tube * tube_speed
    predicted_fixed_R_rate = stress_flux + pressure_work + tube_sweep
    mass_gradient_rhs = FOUR_PI * FIXED_AREAL_RADIUS**2 * (
        rho_tube * Rr_tube + j_tube * Rdot_tube / alpha_tube
    )
    mass_gradient_residual = Mr_tube - mass_gradient_rhs

    chain_rule = None
    if geometry_stage_capture is not None:
        split = geometry_stage_capture
        Aa_dot = (
            np.asarray(split["primary_l2"]["Aa"], dtype=float)
            + np.asarray(split["primary_l3"]["Aa"], dtype=float)
        )
        K_dot = (
            np.asarray(split["primary_l2"]["K"], dtype=float)
            + np.asarray(split["primary_l3"]["K"], dtype=float)
        )
        gamma_inv_dot = gamma_rr_inv * (
            2.0 * X_dot / geom.X - a_dot / geom.a
        )
        Ktheta = geom.K / 3.0 - geom.Aa / 2.0
        Ktheta_dot = K_dot / 3.0 - Aa_dot / 2.0
        Rdot_r = np.asarray(
            grid.cell_derivative_fourth(Rdot, parity=-1), dtype=float
        )
        chain_profiles = _misner_sharp_chain_rule_rate_profiles(
            R, Rdot, Rr, Rdot_r, gamma_rr_inv, gamma_inv_dot,
            Ktheta, Ktheta_dot,
        )
        chain_Mr = np.asarray(
            grid.cell_derivative_fourth(chain_profiles["mass"], parity=-1),
            dtype=float,
        )
        chain_Mdot_coord = at(chain_profiles["mass_dot_coordinate"])
        chain_Mr_tube = at(chain_Mr)
        chain_sweep = chain_Mr_tube * tube_speed
        chain_Mdot_fixed_R = chain_Mdot_coord + chain_sweep
        chain_coordinate_residual = chain_Mdot_coord - (stress_flux + pressure_work)
        chain_fixed_R_residual = chain_Mdot_fixed_R - predicted_fixed_R_rate

        # Direct derivative of the same np.interp(R*=10, R, M) observable used
        # for accepted-state mass. Preserve all stage profiles and the exact
        # moving bracket so this can be audited independently of D_r M / D_r R.
        interpolant_rate = _fixed_target_interpolant_rate(
            R, Rdot, chain_profiles["mass"],
            chain_profiles["mass_dot_coordinate"], target=FIXED_AREAL_RADIUS,
        )
        if not math.isclose(
            interpolant_rate["interpolated_mass"], mass_tube,
            rel_tol=0.0, abs_tol=64.0 * np.finfo(float).eps * max(1.0, abs(mass_tube)),
        ):
            raise RuntimeError("direct fixed-R interpolant mass does not match accepted observable")
        interpolant_coordinate_residual = (
            interpolant_rate["interpolated_mass_dot_coordinate"] - (stress_flux + pressure_work)
        )
        interpolant_fixed_R_residual = (
            interpolant_rate["mass_dot_fixed_R"] - predicted_fixed_R_rate
        )
        chain_rule = {
            "tube_chain_rule_Rdot": at(Rdot),
            "tube_chain_rule_Rr": at(Rr),
            "tube_chain_rule_Rdot_r": at(Rdot_r),
            "tube_chain_rule_gamma_inv": at(gamma_rr_inv),
            "tube_chain_rule_gamma_inv_dot": at(gamma_inv_dot),
            "tube_chain_rule_Ktheta": at(Ktheta),
            "tube_chain_rule_Ktheta_dot": at(Ktheta_dot),
            "tube_chain_rule_Aa_dot": at(Aa_dot),
            "tube_chain_rule_K_dot": at(K_dot),
            "tube_chain_rule_chi_dot": at(chain_profiles["chi_dot"]),
            "tube_chain_rule_chi_dot_absmax": float(np.max(np.abs(chain_profiles["chi_dot"]))),
            "tube_chain_rule_mass_dot_coordinate": chain_Mdot_coord,
            "tube_chain_rule_mass_gradient": chain_Mr_tube,
            "tube_chain_rule_fixed_R_sweep": chain_sweep,
            "tube_chain_rule_mass_dot_fixed_R": chain_Mdot_fixed_R,
            "tube_chain_rule_coordinate_rate_minus_flux_pressure": chain_coordinate_residual,
            "tube_chain_rule_fixed_R_residual": chain_fixed_R_residual,
            "tube_chain_rule_gamma_inv_dot_absmax": float(np.max(np.abs(gamma_inv_dot))),
            "tube_chain_rule_Ktheta_dot_absmax": float(np.max(np.abs(Ktheta_dot))),
            "tube_interpolant_bracket_left_index": int(interpolant_rate["left_index"]),
            "tube_interpolant_bracket_right_index": int(interpolant_rate["right_index"]),
            "tube_interpolant_theta": interpolant_rate["theta"],
            "tube_interpolant_interpolated_mass": interpolant_rate["interpolated_mass"],
            "tube_interpolant_Rdot": interpolant_rate["interpolated_Rdot"],
            "tube_interpolant_mass_dot_coordinate": interpolant_rate["interpolated_mass_dot_coordinate"],
            "tube_interpolant_piecewise_dM_dR": interpolant_rate["piecewise_linear_mass_gradient"],
            "tube_interpolant_fixed_R_sweep": interpolant_rate["fixed_R_sweep"],
            "tube_interpolant_mass_dot_fixed_R": interpolant_rate["mass_dot_fixed_R"],
            "tube_interpolant_coordinate_rate_minus_flux_pressure": interpolant_coordinate_residual,
            "tube_interpolant_fixed_R_residual": interpolant_fixed_R_residual,
            "tube_interpolant_witness_profiles": {
                "R": np.asarray(R, dtype=float).copy(),
                "Rdot": np.asarray(Rdot, dtype=float).copy(),
                "mass": np.asarray(chain_profiles["mass"], dtype=float).copy(),
                "mass_dot_coordinate": np.asarray(chain_profiles["mass_dot_coordinate"], dtype=float).copy(),
            },
        }

    em_rho = total.electromagnetic_rho
    em_rho_absmax = (
        0.0 if em_rho is None else float(np.max(np.abs(em_rho)))
    )
    gauss_absmax = 0.0
    gauss_bound_ratio_max = 0.0
    if state.electromagnetic is not None:
        eta = float(state.electromagnetic.charge_to_rest_mass)
        rest = np.asarray(state.matter.charged_dust.rest, dtype=float)
        faces = np.asarray(state.electromagnetic.enclosed_charge_faces, dtype=float)
        shell = FOUR_PI * eta * float(grid.dr) * rest
        gauss = np.diff(faces) - shell
        operand_scale = np.maximum.reduce((
            np.abs(faces[:-1]), np.abs(faces[1:]), np.abs(np.diff(faces)),
            np.abs(shell),
        ))
        bound = 64.0 * np.finfo(np.float64).eps * np.maximum(
            operand_scale, np.finfo(np.float64).tiny
        )
        gauss_absmax = float(np.max(np.abs(gauss)))
        gauss_bound_ratio_max = float(np.max(np.abs(gauss) / bound))
        if not np.all(np.abs(gauss) <= bound):
            raise AssertionError(
                f"independently recomputed per-shell Gauss bound failed at t={state.t}"
            )

    history_row = state.history[-1] if state.history else {}
    constraint_witnesses = {}
    for key in (
        "hamiltonian_full_absmax", "momentum_full_absmax",
        "connection_full_absmax", "determinant_constraint_max",
        "hamiltonian_max", "momentum_max", "connection_max",
    ):
        if key in history_row:
            value = history_row[key]
            if isinstance(value, (int, float, np.number)) and math.isfinite(float(value)):
                constraint_witnesses[key] = float(value)

    return {
        "t": float(state.t),
        "tau": float(state.tau),
        "resolution": int(grid.n),
        "q_over_m": (
            None if state.electromagnetic is None
            else float(state.electromagnetic.charge_to_rest_mass)
        ),
        "shift_absmax": float(np.max(np.abs(geom.beta))),
        "matter_energy_integrals": _state_energy_integrals(state),
        "em_field_energy": _em_energy(state),
        "outer_enclosed_charge": (
            0.0 if state.electromagnetic is None
            else float(state.electromagnetic.enclosed_charge_faces[-1])
        ),
        "em_stress_rho_absmax": em_rho_absmax,
        "gauss_shell_residual_absmax_recomputed": gauss_absmax,
        "gauss_roundoff_bound_ratio_max": gauss_bound_ratio_max,
        "fixed_areal_radius": FIXED_AREAL_RADIUS,
        "tube_coordinate_r": at(grid.centers),
        "tube_mass_misner_sharp": mass_tube,
        "tube_Rr": Rr_tube,
        "tube_Rdot_at_fixed_coordinate": Rdot_tube,
        "tube_coordinate_speed_for_fixed_R": tube_speed,
        "tube_mass_gradient": Mr_tube,
        "tube_stress_flux": stress_flux,
        "tube_radial_pressure_work": pressure_work,
        "tube_coordinate_sweep": tube_sweep,
        "tube_predicted_fixed_R_mass_rate": predicted_fixed_R_rate,
        "tube_native_mass_gradient_identity_rhs": mass_gradient_rhs,
        "tube_native_mass_gradient_identity_residual": mass_gradient_residual,
        **(chain_rule or {}),
        "constraint_witnesses": constraint_witnesses,
    }


def _fixed_R_rate_snapshot(sample: dict) -> dict:
    """Keep signed fixed-R rate components for one actual PIRK RHS stage."""
    keys = (
        "tube_coordinate_r",
        "tube_mass_misner_sharp",
        "tube_stress_flux",
        "tube_radial_pressure_work",
        "tube_coordinate_sweep",
        "tube_predicted_fixed_R_mass_rate",
        "tube_native_mass_gradient_identity_residual",
        "tube_chain_rule_Rdot",
        "tube_chain_rule_Rr",
        "tube_chain_rule_Rdot_r",
        "tube_chain_rule_gamma_inv",
        "tube_chain_rule_gamma_inv_dot",
        "tube_chain_rule_Ktheta",
        "tube_chain_rule_Ktheta_dot",
        "tube_chain_rule_Aa_dot",
        "tube_chain_rule_K_dot",
        "tube_chain_rule_chi_dot",
        "tube_chain_rule_chi_dot_absmax",
        "tube_chain_rule_mass_dot_coordinate",
        "tube_chain_rule_mass_gradient",
        "tube_chain_rule_mass_dot_fixed_R",
        "tube_chain_rule_coordinate_rate_minus_flux_pressure",
        "tube_chain_rule_fixed_R_residual",
        "tube_chain_rule_fixed_R_sweep",
        "tube_interpolant_bracket_left_index",
        "tube_interpolant_bracket_right_index",
        "tube_interpolant_theta",
        "tube_interpolant_interpolated_mass",
        "tube_interpolant_Rdot",
        "tube_interpolant_mass_dot_coordinate",
        "tube_interpolant_piecewise_dM_dR",
        "tube_interpolant_fixed_R_sweep",
        "tube_interpolant_mass_dot_fixed_R",
        "tube_interpolant_coordinate_rate_minus_flux_pressure",
        "tube_interpolant_fixed_R_residual",
    )
    snapshot = {"t": float(sample["t"]), **{key: float(sample[key]) for key in keys}}
    if "tube_interpolant_witness_profiles" in sample:
        snapshot["tube_interpolant_witness_profiles"] = sample["tube_interpolant_witness_profiles"]
    return snapshot


def _species_name(stage_state, species_state) -> str:
    for name in ("dark_matter", "baryons", "radiation", "charged_dust"):
        current = getattr(stage_state.matter, name, None)
        if current is not None and current is species_state:
            return name
    raise RuntimeError("could not identify the species state at an active RHS stage")


class AuditedKernel(V55ProductionKernel):
    """Production kernel with passive observers around unchanged stage methods."""

    def __init__(self):
        self.audit_runs = []
        self._active_run = None
        self._active_step = None
        self._active_rhs = None

    def initialize(self, *args, **kwargs):
        state = super().initialize(*args, **kwargs)
        self._active_run = {
            "resolution": int(state.grid.n),
            "q_over_m": (
                None if state.electromagnetic is None
                else float(state.electromagnetic.charge_to_rest_mass)
            ),
            "states": [state],
            "steps": [],
        }
        self.audit_runs.append(self._active_run)
        return state

    def step(self, state, dt):
        record = {
            "t_before": float(state.t),
            "dt": float(dt),
            "state_before": state,
            "rhs_stages": [],
            "geometry_stages": [],
            "carrier_rest_fluxes": [],
        }
        self._active_step = record
        try:
            result = super().step(state, dt)
        except Exception:
            self._active_step = None
            self._active_rhs = None
            raise
        record["state_after"] = result
        record["t_after"] = float(result.t)
        self._active_run["steps"].append(record)
        self._active_run["states"].append(result)
        self._active_step = None
        self._active_rhs = None
        return result

    def _rhs(self, state):
        if self._active_step is None:
            return super()._rhs(state)
        record = {
            "stage_index": len(self._active_step["rhs_stages"]),
            "stage": "start" if len(self._active_step["rhs_stages"]) == 0 else "primary_predictor",
            "t": float(state.t),
            "stage_state": state,
            "components": {},
            "lorentz": None,
        }
        self._active_step["rhs_stages"].append(record)
        self._active_rhs = record
        try:
            scalar_rhs, matter_rhs, derivatives = super()._rhs(state)
        finally:
            self._active_rhs = None

        dr = float(state.grid.dr)
        for name, rhs in matter_rhs.items():
            components = record["components"].get(name)
            if components is None:
                raise RuntimeError(f"missing captured finite-volume components for {name}")
            if name == "charged_dust":
                if record["lorentz"] is None:
                    raise RuntimeError("charged-dust RHS stage missed actual Lorentz-source call")
                force_e = record["lorentz"]["energy"]
                force_s = record["lorentz"]["momentum"]
            else:
                force_e = np.zeros_like(components["source_energy"])
                force_s = np.zeros_like(components["source_momentum"])
            flux = components["face_flux"]
            flux_divergence = -np.diff(flux[:, 1]) / dr
            expected_energy_rhs, rhs_energy_residual, rhs_energy_bound_ratio_max = (
                _reconstruct_energy_rhs_with_roundoff_contract(
                    rhs.energy_t,
                    (
                        ("finite-volume energy-flux divergence", flux_divergence),
                        ("geometric energy source", components["source_energy"]),
                        ("scalar-exchange energy source", components["dm_exchange_energy"]),
                        ("Lorentz energy source", force_e),
                    ),
                    label=f"{name} {record['stage']} energy RHS",
                )
            )
            components["rhs_energy_profile"] = _array(rhs.energy_t)
            components["rhs_energy_from_named_terms_profile"] = _array(expected_energy_rhs)
            components["rhs_energy_component_residual_profile"] = _array(rhs_energy_residual)
            components["rhs_energy_component_roundoff_bound_ratio_max"] = rhs_energy_bound_ratio_max
            components["rhs_energy_integrated_rate"] = (
                FOUR_PI * dr * float(np.sum(rhs.energy_t))
            )
            components["boundary_energy_rate"] = (
                FOUR_PI * float(flux[0, 1] - flux[-1, 1])
            )
            components["geometric_energy_source_rate"] = (
                FOUR_PI * dr * float(np.sum(components["source_energy"]))
            )
            components["dm_scalar_exchange_energy_rate"] = (
                FOUR_PI * dr * float(np.sum(components["dm_exchange_energy"]))
            )
            components["lorentz_energy_rate"] = (
                FOUR_PI * dr * float(np.sum(force_e))
            )
            components["componentwise_energy_rate"] = (
                components["boundary_energy_rate"]
                + components["geometric_energy_source_rate"]
                + components["dm_scalar_exchange_energy_rate"]
                + components["lorentz_energy_rate"]
            )
            components["componentwise_energy_rate_residual"] = (
                components["rhs_energy_integrated_rate"]
                - components["componentwise_energy_rate"]
            )
        record["matter_rhs"] = matter_rhs
        return scalar_rhs, matter_rhs, derivatives

    def _geometry_stage_terms(self, grid, geometry, scalars, matter, electromagnetic):
        result = super()._geometry_stage_terms(
            grid, geometry, scalars, matter, electromagnetic
        )
        if self._active_step is not None:
            self._active_step["geometry_stages"].append({
                "geometry": {
                    name: _array(getattr(geometry, name))
                    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K")
                },
                "explicit": {
                    name: _array(result["explicit"][name])
                    for name in ("a", "b", "X")
                },
                "primary_l2": {
                    name: _array(result["primary_l2"][name])
                    for name in ("Aa", "K")
                },
                "primary_l3": {
                    name: _array(result["primary_l3"][name])
                    for name in ("Aa", "K")
                },
            })
        return result

    def _charged_dust_mass_flux_faces(self, state):
        values = super()._charged_dust_mass_flux_faces(state)
        if self._active_step is not None:
            self._active_step["carrier_rest_fluxes"].append(_array(values))
        return values


@contextmanager
def _instrument_rhs(kernel: AuditedKernel):
    """Observe component return values and actual Lorentz calls; restore on exit."""
    original_components = matter_rhs_module._species_rhs_components
    original_lorentz = RadialMaxwellState.lorentz_source

    def components_wrapper(*args, **kwargs):
        result = original_components(*args, **kwargs)
        active = kernel._active_rhs
        if active is not None:
            metric, metric_derivatives, species_state, species = args[:4]
            name = _species_name(active["stage_state"], species_state)
            face_flux, source_energy, source_momentum, dm_energy, dm_momentum = result
            active["components"][name] = {
                "species": str(getattr(species, "value", species)),
                "face_flux": _array(face_flux),
                "source_energy": _array(source_energy),
                "source_momentum": _array(source_momentum),
                "dm_exchange_energy": _array(dm_energy),
                "dm_exchange_momentum": _array(dm_momentum),
            }
        return result

    def lorentz_wrapper(self, *args, **kwargs):
        result = original_lorentz(self, *args, **kwargs)
        active = kernel._active_rhs
        if active is not None:
            force_energy, force_momentum, values = result
            active["lorentz"] = {
                "energy": _array(force_energy),
                "momentum": _array(force_momentum),
                "charge_density": _array(values["rho_e"]),
                "current_r": _array(values["j_e_r"]),
            }
        return result

    matter_rhs_module._species_rhs_components = components_wrapper
    RadialMaxwellState.lorentz_source = lorentz_wrapper
    try:
        yield
    finally:
        matter_rhs_module._species_rhs_components = original_components
        RadialMaxwellState.lorentz_source = original_lorentz


def _field_stage_powers(stage: dict, rest_mass_flux, explicit_geometry: dict) -> dict:
    """Use the same stage-local discretization/sign convention as #180."""
    state = stage["stage_state"]
    if state.electromagnetic is None:
        return {
            "charge_power": 0.0, "metric_power": 0.0,
            "matter_lorentz_power": 0.0, "outer_current_power": 0.0,
        }
    grid = state.grid
    geom = state.geometry
    em = state.electromagnetic
    r = np.asarray(grid.centers, dtype=float)
    q = np.asarray(em.enclosed_charge_faces, dtype=float)
    q_cell = 0.5 * (q[:-1] + q[1:])
    weight = (
        np.asarray(geom.X, dtype=float) * np.sqrt(np.asarray(geom.a, dtype=float))
        / (np.square(r) * np.asarray(geom.b, dtype=float))
    )
    current_faces = em.charge_flux_faces(rest_mass_flux)
    q_cell_dot = -0.5 * (current_faces[:-1] + current_faces[1:])
    charge_power = float(grid.dr) * float(np.sum(weight * q_cell * q_cell_dot))
    log_weight_dot = (
        explicit_geometry["X"] / geom.X
        + 0.5 * explicit_geometry["a"] / geom.a
        - explicit_geometry["b"] / geom.b
    )
    metric_power = 0.5 * float(grid.dr) * float(
        np.sum(weight * q_cell**2 * log_weight_dot)
    )
    lorentz = stage["lorentz"]
    if lorentz is None:
        raise RuntimeError("EM stage has no observed Lorentz source")
    matter_lorentz_power = FOUR_PI * float(grid.dr) * float(np.sum(lorentz["energy"]))
    outer_current_power = (
        -0.5 * float(grid.dr) * weight[-1] * q_cell[-1] * current_faces[-1]
    )
    return {
        "charge_power": charge_power,
        "metric_power": metric_power,
        "matter_lorentz_power": matter_lorentz_power,
        "outer_current_power": float(outer_current_power),
    }


def _stage_serializable(stage: dict) -> dict:
    state = stage["stage_state"]
    species_rows = {}
    for name, comp in stage["components"].items():
        species_rows[name] = {
            "species_kind": comp["species"],
            "face_flux_rest_mass": comp["face_flux"][:, 0].tolist(),
            "face_flux_energy": comp["face_flux"][:, 1].tolist(),
            "face_flux_radial_momentum": comp["face_flux"][:, 2].tolist(),
            "geometric_energy_source_profile": comp["source_energy"].tolist(),
            "geometric_momentum_source_profile": comp["source_momentum"].tolist(),
            "dm_scalar_exchange_energy_source_profile": comp["dm_exchange_energy"].tolist(),
            "dm_scalar_exchange_momentum_source_profile": comp["dm_exchange_momentum"].tolist(),
            "lorentz_energy_source_profile": (
                stage["lorentz"]["energy"].tolist() if name == "charged_dust"
                and stage["lorentz"] is not None else [0.0] * int(state.grid.n)
            ),
            "lorentz_momentum_source_profile": (
                stage["lorentz"]["momentum"].tolist() if name == "charged_dust"
                and stage["lorentz"] is not None else [0.0] * int(state.grid.n)
            ),
            "boundary_energy_rate": comp["boundary_energy_rate"],
            "geometric_energy_source_rate": comp["geometric_energy_source_rate"],
            "dm_scalar_exchange_energy_rate": comp["dm_scalar_exchange_energy_rate"],
            "lorentz_energy_rate": comp["lorentz_energy_rate"],
            "rhs_energy_integrated_rate": comp["rhs_energy_integrated_rate"],
            "componentwise_energy_rate": comp["componentwise_energy_rate"],
            "componentwise_energy_rate_residual": comp["componentwise_energy_rate_residual"],
            "rhs_energy_component_residual_absmax": float(
                np.max(np.abs(comp["rhs_energy_component_residual_profile"]))
            ),
            "rhs_energy_component_roundoff_bound_ratio_max": (
                comp["rhs_energy_component_roundoff_bound_ratio_max"]
            ),
        }
    return {
        "stage": stage["stage"],
        "t": stage["t"],
        "resolution": int(state.grid.n),
        "beta_absmax": float(np.max(np.abs(state.geometry.beta))),
        "species": species_rows,
        "lorentz_charge_density_absmax": (
            0.0 if stage["lorentz"] is None
            else float(np.max(np.abs(stage["lorentz"]["charge_density"])))
        ),
        "lorentz_current_absmax": (
            0.0 if stage["lorentz"] is None
            else float(np.max(np.abs(stage["lorentz"]["current_r"])))
        ),
    }


def _serialize_step(kernel, step: dict) -> dict:
    stages = step["rhs_stages"]
    geometry_stages = step["geometry_stages"]
    fluxes = step["carrier_rest_fluxes"]
    if len(stages) != 2 or len(geometry_stages) != 3 or len(fluxes) != 2:
        raise RuntimeError(
            "PIRK observer did not see expected start/predictor RHS, three geometry "
            "evaluations, and paired charged-dust rest-current fluxes"
        )
    stage_power_rows = []
    for stage, geom_index, flux in (
        (stages[0], 0, fluxes[0]),
        (stages[1], 2, fluxes[1]),
    ):
        stage_power_rows.append(_field_stage_powers(
            stage, flux, geometry_stages[geom_index]["explicit"]
        ))
    avg = lambda key: 0.5 * (
        stage_power_rows[0][key] + stage_power_rows[1][key]
    )

    before = step["state_before"]
    after = step["state_after"]
    dt = float(step["dt"])
    # Preserve rates evaluated on the actual stage states, not only accepted slices.
    fixed_R_stage_rates = {
        "start": _fixed_R_rate_snapshot(_worldtube_sample(
            kernel, stages[0]["stage_state"], geometry_stage_capture=geometry_stages[0]
        )),
        "primary_predictor": _fixed_R_rate_snapshot(_worldtube_sample(
            kernel, stages[1]["stage_state"], geometry_stage_capture=geometry_stages[2]
        )),
    }
    matter_before = _state_energy_integrals(before)
    matter_after = _state_energy_integrals(after)
    matter_step_rows = {}
    for name in matter_before:
        stage0 = stages[0]["components"][name]["rhs_energy_integrated_rate"]
        stage1 = stages[1]["components"][name]["rhs_energy_integrated_rate"]
        predicted = 0.5 * (stage0 + stage1)
        observed = (matter_after[name] - matter_before[name]) / dt
        terms = []
        for stage in stages:
            comp = stage["components"][name]
            terms.append({
                "boundary_energy_rate": comp["boundary_energy_rate"],
                "geometric_energy_source_rate": comp["geometric_energy_source_rate"],
                "dm_scalar_exchange_energy_rate": comp["dm_scalar_exchange_energy_rate"],
                "lorentz_energy_rate": comp["lorentz_energy_rate"],
                "rhs_energy_integrated_rate": comp["rhs_energy_integrated_rate"],
                "componentwise_energy_rate_residual": comp["componentwise_energy_rate_residual"],
            })
        matter_step_rows[name] = {
            "energy_before": matter_before[name],
            "energy_after": matter_after[name],
            "finite_difference_integrated_rate": observed,
            "trapezoid_rhs_integrated_rate": predicted,
            "temporal_update_remainder": observed - predicted,
            "stage_terms": terms,
        }

    em_before = _em_energy(before)
    em_after = _em_energy(after)
    em_rate = (em_after - em_before) / dt
    em_metric = avg("metric_power")
    em_lorentz = avg("matter_lorentz_power")
    em_outer = avg("outer_current_power")
    local_em_residual = em_rate - em_metric + em_lorentz - em_outer

    matter_boundary = 0.5 * sum(
        stages[k]["components"][name]["boundary_energy_rate"]
        for k in (0, 1)
        for name in stages[k]["components"]
    )
    matter_geometry_source = 0.5 * sum(
        stages[k]["components"][name]["geometric_energy_source_rate"]
        for k in (0, 1)
        for name in stages[k]["components"]
    )
    dm_exchange_source = 0.5 * sum(
        stages[k]["components"][name]["dm_scalar_exchange_energy_rate"]
        for k in (0, 1)
        for name in stages[k]["components"]
    )
    measured_matter_rate = sum(
        matter_step_rows[name]["finite_difference_integrated_rate"]
        for name in matter_step_rows
    )
    measured_combined_rate = measured_matter_rate + em_rate
    known_combined_terms = (
        matter_boundary + matter_geometry_source + dm_exchange_source
        + em_metric + em_outer
    )
    # This is intentionally named an unassigned remainder: it is not a global
    # conservation residual because scalar-field energy evolution, metric/gravity
    # accounting, and any Killing-energy premise are not fully assembled here.
    combined_remainder = measured_combined_rate - known_combined_terms

    return {
        "t_before": float(step["t_before"]),
        "t_after": float(step["t_after"]),
        "dt": dt,
        "fixed_R_stage_rates": fixed_R_stage_rates,
        "matter_species": matter_step_rows,
        "em_field_energy_before": em_before,
        "em_field_energy_after": em_after,
        "em_field_energy_rate_measured": em_rate,
        "em_charge_current_power_start_predictor_average": avg("charge_power"),
        "em_metric_power_start_predictor_average": em_metric,
        "em_matter_lorentz_power_start_predictor_average": em_lorentz,
        "em_outer_current_power_start_predictor_average": em_outer,
        "em_local_exchange_residual": local_em_residual,
        "em_local_exchange_relative_residual": abs(local_em_residual) / max(
            abs(avg("charge_power")) + abs(em_lorentz) + abs(em_metric) + abs(em_outer),
            np.finfo(float).tiny,
        ),
        "matter_lorentz_power_start_predictor_average": 0.5 * sum(
            sum(
                stages[k]["components"][name]["lorentz_energy_rate"]
                for name in stages[k]["components"]
            )
            for k in (0, 1)
        ),
        "matter_boundary_energy_rate_average": matter_boundary,
        "matter_geometric_energy_source_rate_average": matter_geometry_source,
        "dm_scalar_exchange_energy_rate_average": dm_exchange_source,
        "combined_matter_plus_em_energy_rate_measured": measured_combined_rate,
        "combined_rate_from_named_boundary_and_source_terms": known_combined_terms,
        "unassigned_matter_plus_em_rate_remainder_not_global_conservation_residual": combined_remainder,
        "source_stage_rows": [_stage_serializable(stage) for stage in stages],
    }


def _finish_run(kernel, run: dict, result, output_dir: Path) -> dict:
    if result.status != "completed" or not math.isclose(
        float(result.t_final), FINAL_TIME, rel_tol=0.0, abs_tol=1.0e-11
    ):
        raise RuntimeError(
            f"campaign did not complete N={run['resolution']}: "
            f"status={result.status}, t_final={result.t_final}"
        )
    retry_path = output_dir / f"N{run['resolution']}" / "radiation_step_retries.json"
    retries = json.loads(retry_path.read_text(encoding="utf-8"))
    if retries:
        raise RuntimeError(
            f"unexpected radiation retry(s) in bounded audit N={run['resolution']}: {retries}"
        )
    state_rows = [_worldtube_sample(kernel, state) for state in run["states"]]
    if len(state_rows) < 4:
        raise RuntimeError("fixed-areal-radius rate needs at least four accepted samples")
    times = np.asarray([row["t"] for row in state_rows], dtype=float)
    masses = np.asarray([row["tube_mass_misner_sharp"] for row in state_rows], dtype=float)
    fd_rates = np.gradient(masses, times, edge_order=2)
    for row, rate in zip(state_rows, fd_rates):
        row["tube_fixed_R_mass_rate_finite_difference"] = float(rate)
        row["tube_fixed_R_mass_rate_residual"] = (
            float(rate) - row["tube_predicted_fixed_R_mass_rate"]
        )

    steps = [_serialize_step(kernel, step) for step in run["steps"]]
    return {
        "resolution": int(run["resolution"]),
        "q_over_m": float(run["q_over_m"]),
        "steps": len(steps),
        "t_final": float(result.t_final),
        "source_line": "explicit opt-in separate charged dust; Gaussian c=1 research convention",
        "accepted_slice_rows": state_rows,
        "step_rows": steps,
    }


def _fixed_R_interval_balance(run: dict, *, label: str = "run") -> list[dict]:
    """Compare adjacent accepted-state M_MS rates with endpoint and live-stage rates."""
    states = run.get("accepted_slice_rows", [])
    steps = run.get("step_rows", [])
    if len(states) < 2:
        raise RuntimeError(f"{label}: need at least two accepted states")
    if len(steps) != len(states) - 1:
        raise RuntimeError(f"{label}: expected one recorded PIRK step per accepted-state interval")
    times = np.asarray([row["t"] for row in states], dtype=float)
    masses = np.asarray([row["tube_mass_misner_sharp"] for row in states], dtype=float)
    predicted = np.asarray([row["tube_predicted_fixed_R_mass_rate"] for row in states], dtype=float)
    if not (np.all(np.isfinite(times)) and np.all(np.isfinite(masses)) and np.all(np.isfinite(predicted))):
        raise RuntimeError(f"{label}: non-finite accepted-state fixed-R operand")
    dts = np.diff(times)
    if np.any(dts <= 0.0):
        raise RuntimeError(f"{label}: accepted-state times must be strictly increasing")

    rows = []
    for i, step in enumerate(steps):
        dt = float(dts[i])
        if not math.isclose(float(step["t_before"]), float(times[i]), rel_tol=0.0, abs_tol=1e-12):
            raise RuntimeError(f"{label}: step {i} start time does not match accepted state")
        if not math.isclose(float(step["t_after"]), float(times[i + 1]), rel_tol=0.0, abs_tol=1e-12):
            raise RuntimeError(f"{label}: step {i} end time does not match accepted state")
        if not math.isclose(float(step["dt"]), dt, rel_tol=1e-10, abs_tol=1e-14):
            raise RuntimeError(f"{label}: step {i} dt does not match accepted interval")
        stage_rates = step.get("fixed_R_stage_rates")
        if not isinstance(stage_rates, dict) or not {"start", "primary_predictor"} <= set(stage_rates):
            raise RuntimeError(f"{label}: step {i} missing start/predictor fixed-R stage rates")

        mass_rate = float((masses[i + 1] - masses[i]) / dt)
        endpoint_prediction = float(0.5 * (predicted[i] + predicted[i + 1]))
        start_rate = float(stage_rates["start"]["tube_predicted_fixed_R_mass_rate"])
        predictor_rate = float(stage_rates["primary_predictor"]["tube_predicted_fixed_R_mass_rate"])
        stage_pair_average = 0.5 * (start_rate + predictor_rate)
        stage_chain_fields = (
            "tube_chain_rule_Rdot",
            "tube_chain_rule_Rr",
            "tube_chain_rule_Rdot_r",
            "tube_chain_rule_gamma_inv_dot",
            "tube_chain_rule_Ktheta_dot",
            "tube_chain_rule_Aa_dot",
            "tube_chain_rule_K_dot",
            "tube_chain_rule_chi_dot",
            "tube_chain_rule_mass_dot_coordinate",
            "tube_chain_rule_mass_dot_fixed_R",
            "tube_chain_rule_coordinate_rate_minus_flux_pressure",
            "tube_chain_rule_fixed_R_residual",
            "tube_chain_rule_fixed_R_sweep",
            "tube_interpolant_bracket_left_index",
            "tube_interpolant_bracket_right_index",
            "tube_interpolant_theta",
            "tube_interpolant_interpolated_mass",
            "tube_interpolant_Rdot",
            "tube_interpolant_mass_dot_coordinate",
            "tube_interpolant_piecewise_dM_dR",
            "tube_interpolant_fixed_R_sweep",
            "tube_interpolant_mass_dot_fixed_R",
            "tube_interpolant_coordinate_rate_minus_flux_pressure",
            "tube_interpolant_fixed_R_residual",
        )
        for stage_name in ("start", "primary_predictor"):
            missing = [key for key in stage_chain_fields if key not in stage_rates[stage_name]]
            if missing:
                raise RuntimeError(
                    f"{label}: stage {stage_name} is missing chain-rule outputs: {missing}"
                )
        rows.append({
            "interval_index": i,
            "t_start": float(times[i]),
            "t_end": float(times[i + 1]),
            "dt": dt,
            "mass_rate_finite_difference": mass_rate,
            "accepted_state_endpoint_trapezoid_prediction": endpoint_prediction,
            "accepted_interval_residual": mass_rate - endpoint_prediction,
            "start_stage_predicted_rate": start_rate,
            "primary_predictor_stage_predicted_rate": predictor_rate,
            "stage_pair_average_prediction_diagnostic_only": stage_pair_average,
            "mass_rate_minus_stage_pair_average_diagnostic_only": mass_rate - stage_pair_average,
            "start_stage_components": stage_rates["start"],
            "primary_predictor_stage_components": stage_rates["primary_predictor"],
            **{
                f"start_stage_{key}": float(stage_rates["start"][key])
                for key in stage_chain_fields
            },
            **{
                f"primary_predictor_stage_{key}": float(stage_rates["primary_predictor"][key])
                for key in stage_chain_fields
            },
        })
    return rows


def _pair_fixed_R_interval_rows(charged: list[dict], neutral: list[dict], *, label: str) -> list[dict]:
    """Pair charged and neutral intervals only when their accepted time grids align."""
    if len(charged) != len(neutral):
        raise RuntimeError(f"{label}: charged/neutral interval counts do not match")
    result = []
    delta_fields = (
        "mass_rate_finite_difference",
        "accepted_state_endpoint_trapezoid_prediction",
        "accepted_interval_residual",
        "start_stage_predicted_rate",
        "primary_predictor_stage_predicted_rate",
        "stage_pair_average_prediction_diagnostic_only",
        "mass_rate_minus_stage_pair_average_diagnostic_only",
        "start_stage_tube_chain_rule_Rdot",
        "primary_predictor_stage_tube_chain_rule_Rdot",
        "start_stage_tube_chain_rule_Rr",
        "primary_predictor_stage_tube_chain_rule_Rr",
        "start_stage_tube_chain_rule_Rdot_r",
        "primary_predictor_stage_tube_chain_rule_Rdot_r",
        "start_stage_tube_chain_rule_gamma_inv_dot",
        "primary_predictor_stage_tube_chain_rule_gamma_inv_dot",
        "start_stage_tube_chain_rule_Ktheta_dot",
        "primary_predictor_stage_tube_chain_rule_Ktheta_dot",
        "start_stage_tube_chain_rule_Aa_dot",
        "primary_predictor_stage_tube_chain_rule_Aa_dot",
        "start_stage_tube_chain_rule_K_dot",
        "primary_predictor_stage_tube_chain_rule_K_dot",
        "start_stage_tube_chain_rule_chi_dot",
        "primary_predictor_stage_tube_chain_rule_chi_dot",
        "start_stage_tube_chain_rule_mass_dot_coordinate",
        "primary_predictor_stage_tube_chain_rule_mass_dot_coordinate",
        "start_stage_tube_chain_rule_mass_dot_fixed_R",
        "primary_predictor_stage_tube_chain_rule_mass_dot_fixed_R",
        "start_stage_tube_chain_rule_coordinate_rate_minus_flux_pressure",
        "primary_predictor_stage_tube_chain_rule_coordinate_rate_minus_flux_pressure",
        "start_stage_tube_chain_rule_fixed_R_residual",
        "primary_predictor_stage_tube_chain_rule_fixed_R_residual",
        "start_stage_tube_chain_rule_fixed_R_sweep",
        "primary_predictor_stage_tube_chain_rule_fixed_R_sweep",
        "start_stage_tube_interpolant_theta",
        "primary_predictor_stage_tube_interpolant_theta",
        "start_stage_tube_interpolant_mass_dot_coordinate",
        "primary_predictor_stage_tube_interpolant_mass_dot_coordinate",
        "start_stage_tube_interpolant_Rdot",
        "primary_predictor_stage_tube_interpolant_Rdot",
        "start_stage_tube_interpolant_piecewise_dM_dR",
        "primary_predictor_stage_tube_interpolant_piecewise_dM_dR",
        "start_stage_tube_interpolant_fixed_R_sweep",
        "primary_predictor_stage_tube_interpolant_fixed_R_sweep",
        "start_stage_tube_interpolant_mass_dot_fixed_R",
        "primary_predictor_stage_tube_interpolant_mass_dot_fixed_R",
        "start_stage_tube_interpolant_coordinate_rate_minus_flux_pressure",
        "primary_predictor_stage_tube_interpolant_coordinate_rate_minus_flux_pressure",
        "start_stage_tube_interpolant_fixed_R_residual",
        "primary_predictor_stage_tube_interpolant_fixed_R_residual",
    )
    for i, (c, n) in enumerate(zip(charged, neutral)):
        for time_key in ("t_start", "t_end", "dt"):
            if not math.isclose(float(c[time_key]), float(n[time_key]), rel_tol=0.0, abs_tol=1e-12):
                raise RuntimeError(f"{label}: charged/neutral {time_key} mismatch at interval {i}")
        row = {"interval_index": i, "t_start": c["t_start"], "t_end": c["t_end"], "dt": c["dt"]}
        for key in delta_fields:
            row[f"charged_minus_neutral_{key}"] = float(c[key] - n[key])
        result.append(row)
    return result


def _run_fixed_N_temporal_refinement(temp_root: Path) -> list[dict]:
    """Run a paired fixed-N CFL sweep; no production operators are changed."""
    temporal_runs = []
    for cfl in TIME_STEP_CFL_VALUES:
        for eta in CHARGE_TO_MASS_ARMS:
            arm_dir = temp_root / f"fixed_N_cfl_{cfl:.5f}" / f"q_over_m_{eta:.2f}"
            config = CampaignConfig(
                resolutions=(TEMPORAL_REFINEMENT_RESOLUTION,),
                r_max=R_MAX,
                cfl=cfl,
                final_time=FINAL_TIME,
                checkpoint_interval=CHECKPOINT_INTERVAL,
                output_dir=str(arm_dir),
                radial_maxwell_charge_to_mass=eta,
                charged_dust_density_profile=DENSITY_PROFILE,
            )
            kernel = AuditedKernel()
            with _instrument_rhs(kernel):
                results = run_campaign(kernel=kernel, config=config)
            if len(results) != 1 or len(kernel.audit_runs) != 1:
                raise RuntimeError("fixed-N temporal refinement did not collect exactly one run")
            run = _finish_run(kernel, kernel.audit_runs[0], results[0], arm_dir)
            run["cfl"] = float(cfl)
            run["nominal_dt"] = float(cfl * (R_MAX / TEMPORAL_REFINEMENT_RESOLUTION))
            run["accepted_interval_balance"] = _fixed_R_interval_balance(
                run, label=f"N={TEMPORAL_REFINEMENT_RESOLUTION}, CFL={cfl}, q/m={eta}"
            )
            temporal_runs.append(run)
    return temporal_runs


def _summarize_fixed_N_temporal_refinement(runs: list[dict]) -> dict:
    """Summarize matched charged/neutral runs without fitting a physical tolerance."""
    rows = []
    for cfl in TIME_STEP_CFL_VALUES:
        neutral = next(r for r in runs if r["cfl"] == cfl and r["q_over_m"] == 0.0)
        charged = next(r for r in runs if r["cfl"] == cfl and r["q_over_m"] == 0.25)
        n_intervals = neutral["accepted_interval_balance"]
        c_intervals = charged["accepted_interval_balance"]
        paired = _pair_fixed_R_interval_rows(c_intervals, n_intervals, label=f"fixed-N CFL={cfl}")
        if not paired:
            raise RuntimeError(f"fixed-N CFL={cfl}: no accepted intervals to summarize")
        final = paired[-1]
        rows.append({
            "resolution": TEMPORAL_REFINEMENT_RESOLUTION,
            "cfl": float(cfl),
            "nominal_dt": float(charged["nominal_dt"]),
            "neutral_steps": int(neutral["steps"]),
            "charged_steps": int(charged["steps"]),
            "charged_final_field_energy": float(charged["accepted_slice_rows"][-1]["em_field_energy"]),
            "charged_final_outer_charge": float(charged["accepted_slice_rows"][-1]["outer_enclosed_charge"]),
            "charged_final_gauss_shell_residual": float(charged["accepted_slice_rows"][-1]["gauss_shell_residual_absmax_recomputed"]),
            "max_abs_charged_interval_residual": max(abs(x["accepted_interval_residual"]) for x in c_intervals),
            "max_abs_neutral_interval_residual": max(abs(x["accepted_interval_residual"]) for x in n_intervals),
            "first_interval_charged_minus_neutral_residual": paired[0]["charged_minus_neutral_accepted_interval_residual"],
            "final_interval_charged_minus_neutral_residual": final["charged_minus_neutral_accepted_interval_residual"],
            "max_abs_charged_minus_neutral_interval_residual": max(
                abs(x["charged_minus_neutral_accepted_interval_residual"]) for x in paired
            ),
            "final_interval_charged_minus_neutral_mass_rate": final["charged_minus_neutral_mass_rate_finite_difference"],
            "final_interval_charged_minus_neutral_endpoint_prediction": final["charged_minus_neutral_accepted_state_endpoint_trapezoid_prediction"],
            "final_interval_charged_minus_neutral_start_stage_rate": final["charged_minus_neutral_start_stage_predicted_rate"],
            "final_interval_charged_minus_neutral_predictor_stage_rate": final["charged_minus_neutral_primary_predictor_stage_predicted_rate"],
            "final_interval_charged_minus_neutral_stage_pair_average_diagnostic": final["charged_minus_neutral_stage_pair_average_prediction_diagnostic_only"],
            "final_interval_charged_minus_neutral_mass_rate_minus_stage_pair_average": final["charged_minus_neutral_mass_rate_minus_stage_pair_average_diagnostic_only"],
            "final_interval_charged_minus_neutral_start_chain_fixed_R_residual": final["charged_minus_neutral_start_stage_tube_chain_rule_fixed_R_residual"],
            "final_interval_charged_minus_neutral_predictor_chain_fixed_R_residual": final["charged_minus_neutral_primary_predictor_stage_tube_chain_rule_fixed_R_residual"],
            "final_interval_charged_minus_neutral_predictor_interpolant_mass_rate": final["charged_minus_neutral_primary_predictor_stage_tube_interpolant_mass_dot_fixed_R"],
            "final_interval_charged_minus_neutral_predictor_interpolant_coordinate_residual": final["charged_minus_neutral_primary_predictor_stage_tube_interpolant_coordinate_rate_minus_flux_pressure"],
            "final_interval_charged_minus_neutral_predictor_interpolant_fixed_R_residual": final["charged_minus_neutral_primary_predictor_stage_tube_interpolant_fixed_R_residual"],
            "final_interval_charged_minus_neutral_start_chain_coordinate_rate": final["charged_minus_neutral_start_stage_tube_chain_rule_mass_dot_coordinate"],
            "final_interval_charged_minus_neutral_predictor_chain_coordinate_rate": final["charged_minus_neutral_primary_predictor_stage_tube_chain_rule_mass_dot_coordinate"],
            "paired_interval_rows": paired,
        })
    return {
        "resolution": TEMPORAL_REFINEMENT_RESOLUTION,
        "r_max": R_MAX,
        "final_time": FINAL_TIME,
        "cfl_values": list(TIME_STEP_CFL_VALUES),
        "q_over_m_arms": list(CHARGE_TO_MASS_ARMS),
        "interpretation": [
            "Spatial resolution is fixed while the nominal step decreases by factors of two.",
            "Accepted-interval secants are compared with trapezoidal endpoint rates and are kept distinct from actual PIRK start/predictor-stage rates.",
            "The average of start/predictor diagnostic rates is not asserted to be the production integrator's exact quadrature.",
            "The semidiscrete Misner-Sharp chain-rule witness differentiates the mass definition using actual captured explicit and primary L2/L3 stage rates; its stress/pressure remainder is a separate diagnostic, not a global conservation proof.",
            "A second stage witness differentiates the same piecewise-linear np.interp(R*=10, R, M) observable, preserving the active moving bracket, theta, and full R/Rdot/M/Mdot coordinate profiles; it is compared directly with the pointwise chain-rule derivative and stress/pressure terms.",
            "Charged-minus-neutral differences are only a discriminator; no physical/global residual threshold is imposed.",
        ],
        "rows": rows,
    }


def _summaries(runs: list[dict]) -> dict:
    rows = []
    for resolution in RESOLUTIONS:
        neutral = next(r for r in runs if r["resolution"] == resolution and r["q_over_m"] == 0.0)
        charged = next(r for r in runs if r["resolution"] == resolution and r["q_over_m"] == 0.25)
        n_last, c_last = neutral["accepted_slice_rows"][-1], charged["accepted_slice_rows"][-1]
        rows.append({
            "resolution": resolution,
            "neutral_steps": neutral["steps"],
            "charged_steps": charged["steps"],
            "neutral_outer_charge_final": neutral["accepted_slice_rows"][-1]["outer_enclosed_charge"],
            "charged_outer_charge_final": c_last["outer_enclosed_charge"],
            "charged_em_field_energy_final": c_last["em_field_energy"],
            "charged_em_stress_rho_absmax_final": c_last["em_stress_rho_absmax"],
            "charged_gauss_shell_residual_absmax_final": c_last["gauss_shell_residual_absmax_recomputed"],
            "charged_gauss_roundoff_bound_ratio_max": max(
                r["gauss_roundoff_bound_ratio_max"] for r in charged["accepted_slice_rows"]
            ),
            "charged_lorentz_energy_power_absmax": max(
                abs(s["em_matter_lorentz_power_start_predictor_average"])
                for s in charged["step_rows"]
            ),
            "neutral_em_energy_absmax": max(
                abs(s["em_field_energy_after"]) for s in neutral["step_rows"]
            ),
            "charged_minus_neutral_fixed_R_mass_final": (
                c_last["tube_mass_misner_sharp"] - n_last["tube_mass_misner_sharp"]
            ),
            "charged_fixed_R_mass_rate_residual_absmax": max(
                abs(r["tube_fixed_R_mass_rate_residual"]) for r in charged["accepted_slice_rows"]
            ),
            "charged_native_mass_gradient_identity_residual_absmax": max(
                abs(r["tube_native_mass_gradient_identity_residual"])
                for r in charged["accepted_slice_rows"]
            ),
            "max_constraint_witnesses_final": c_last["constraint_witnesses"],
        })
    return {
        "paired_neutral_vs_charged": rows,
        "interpretation": [
            "The full signed stage face flux and source profiles are retained in accepted JSON rows.",
            "Componentwise matter finite-volume terms and the locally measured EM exchange are kept separate.",
            "The combined matter-plus-EM remainder is not a global conservation residual; scalar and gravity energy accounting are not fully assembled here.",
            "The quasi-local fixed-areal-radius ledger is a discrete diagnostic. It uses interpolated R=10, the production stress-flux/pressure-work sign convention, and an explicit coordinate-tube sweep.",
            "The experiment is bounded to t=0.03 and zero shift; it does not establish nonzero-shift coupling, long-horizon stability, continuum convergence, physical normalization/carrier derivation, or scientific admission.",
        ],
    }


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"non-finite number in report: {value!r}")
    return value


def run_joined_audit(root: Path, output_path: Path) -> dict:
    provenance = verify_source_pins(root)
    all_runs = []
    with tempfile.TemporaryDirectory(prefix="em-budget-neutral-charged-") as tmp:
        temp_root = Path(tmp)
        for eta in CHARGE_TO_MASS_ARMS:
            arm_dir = temp_root / f"q_over_m_{eta:.2f}"
            config = CampaignConfig(
                resolutions=RESOLUTIONS,
                r_max=R_MAX,
                cfl=CFL,
                final_time=FINAL_TIME,
                checkpoint_interval=CHECKPOINT_INTERVAL,
                output_dir=str(arm_dir),
                radial_maxwell_charge_to_mass=eta,
                charged_dust_density_profile=DENSITY_PROFILE,
            )
            kernel = AuditedKernel()
            with _instrument_rhs(kernel):
                results = run_campaign(kernel=kernel, config=config)
            if len(results) != len(RESOLUTIONS) or len(kernel.audit_runs) != len(RESOLUTIONS):
                raise RuntimeError("observer did not collect every requested resolution")
            for run, result in zip(kernel.audit_runs, results):
                all_runs.append(_finish_run(kernel, run, result, arm_dir))
        temporal_runs = _run_fixed_N_temporal_refinement(temp_root)

    # Verify the explicit neutral control, but do not assert a fitted threshold
    # on any physical/global residual.
    for resolution in RESOLUTIONS:
        neutral = next(r for r in all_runs if r["resolution"] == resolution and r["q_over_m"] == 0.0)
        for row in neutral["accepted_slice_rows"]:
            if row["em_field_energy"] != 0.0 or row["outer_enclosed_charge"] != 0.0:
                raise AssertionError("q/m=0 neutral-control field/charge must be exactly zero")
    for resolution in RESOLUTIONS:
        charged = next(r for r in all_runs if r["resolution"] == resolution and r["q_over_m"] == 0.25)
        if not any(row["em_field_energy"] > 0.0 for row in charged["accepted_slice_rows"]):
            raise AssertionError("charged arm must produce nonzero EM field energy")

    report = {
        "classification": "source-pinned, test/report-only joined EM/matter/quasi-local ledger; research only",
        "source_anchor": "main source blobs from 5c51bf509504365df394e6cca96bf1d117f4614d",
        "provenance": provenance,
        "settings": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "final_time": FINAL_TIME,
            "cfl": CFL,
            "checkpoint_interval": CHECKPOINT_INTERVAL,
            "density_profile": [list(x) for x in DENSITY_PROFILE],
            "q_over_m_arms": list(CHARGE_TO_MASS_ARMS),
            "fixed_areal_radius": FIXED_AREAL_RADIUS,
            "vendor_gitlink": VENDOR_PIN,
        },
        "runs": all_runs,
        "summary": _summaries(all_runs),
        "fixed_N_temporal_refinement": _summarize_fixed_N_temporal_refinement(temporal_runs),
        "fixed_N_temporal_refinement_runs": temporal_runs,
    }
    report = _json_safe(report)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return report


def write_step_summary(report: dict) -> None:
    import os
    summary_raw = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary_raw:
        return
    summary_path = Path(summary_raw)
    summary = report["summary"]
    lines = [
        "## Joined radial-EM / matter / Misner–Sharp ledger",
        "",
        "Source-pinned, research-only paired run; q/m=0 versus 0.25; N=24/48/96; t=0.03; Gaussian c=1 convention.",
        "",
        f"Execution checkout SHA: `{report['provenance']['execution_checkout_sha']}`; source head SHA (`{report['provenance']['source_head_kind']}`): `{report['provenance']['source_head_sha']}`.",
        "",
        "| N | Neutral steps | Charged steps | Final charged U_EM | Final max Gauss residual | Max R=10 fixed-tube balance remainder |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary["paired_neutral_vs_charged"]:
        lines.append(
            f"| {row['resolution']} | {row['neutral_steps']} | {row['charged_steps']} | "
            f"{row['charged_em_field_energy_final']:.6e} | "
            f"{row['charged_gauss_shell_residual_absmax_final']:.6e} | "
            f"{row['charged_fixed_R_mass_rate_residual_absmax']:.6e} |"
        )
    temporal = report["fixed_N_temporal_refinement"]
    lines.extend([
        "",
        "### Fixed-N time-step discriminator",
        "",
        f"Resolution held at N={temporal['resolution']}; q/m=0 and 0.25 paired at CFL values {temporal['cfl_values']}.",
        "",
        "| CFL | nominal dt | steps | final Δ interval residual | Δ measured M rate | Δ endpoint prediction | Δ pointwise chain fixed-R remainder | Δ direct-interpolant fixed-R remainder |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for row in temporal["rows"]:
        lines.append(
            f"| {row['cfl']:.5g} | {row['nominal_dt']:.6e} | {row['charged_steps']} | "
            f"{row['final_interval_charged_minus_neutral_residual']:.6e} | "
            f"{row['final_interval_charged_minus_neutral_mass_rate']:.6e} | "
            f"{row['final_interval_charged_minus_neutral_endpoint_prediction']:.6e} | "
            f"{row['final_interval_charged_minus_neutral_predictor_chain_fixed_R_residual']:.6e} | "
            f"{row['final_interval_charged_minus_neutral_predictor_interpolant_fixed_R_residual']:.6e} |"
        )
    lines.extend([
        "",
        "The fixed-N sweep retains adjacent accepted-state interval secants, endpoint-rate trapezoids, actual start/predictor-stage rates, the pointwise chain-rule rate, and the derivative of the direct fixed-R piecewise-linear interpolant as separate witnesses. This does not identify a missing physics term or establish a global conserved energy.",
        "",
        "The full signed matter face-flux/source profiles and EM stage terms are in the uploaded JSON artifact. A nonzero fixed-tube remainder is reported as measured, not fitted away. The combined matter-plus-EM remainder is not a global conservation residual because scalar and gravity accounting are not fully assembled.",
        "",
    ])
    with summary_path.open("a", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    report = run_joined_audit(args.repo_root.resolve(), args.output.resolve())
    write_step_summary(report)
    concise = {
        "status": "completed",
        "classification": report["classification"],
        "execution_checkout_sha": report["provenance"]["execution_checkout_sha"],
        "source_head_sha": report["provenance"]["source_head_sha"],
        "source_head_kind": report["provenance"]["source_head_kind"],
        "source_pins_passed": all(
            value["match"]
            for value in report["provenance"]["production_blob_pins"].values()
        ),
        "paired_neutral_vs_charged": report["summary"]["paired_neutral_vs_charged"],
        "fixed_N_temporal_refinement": report["fixed_N_temporal_refinement"],
        "full_report_path": str(args.output.resolve()),
    }
    print("EM_TOTAL_ENERGY_LEDGER_REPORT=" + json.dumps(concise, sort_keys=True, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
