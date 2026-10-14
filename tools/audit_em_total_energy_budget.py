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

from tools.fixed_r_discriminators import (
    frozen_geometry_stress_balance,
    pure_sector_rates,
)


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


def _fixed_R_composite_flux_work_profiles(R, alpha, gamma_inv, current, Rr, pr, Rdot):
    """Build signed tube flux/work profiles before interpolation.

    This exposes the interpolation-order ambiguity: product-of-interpolated
    factors and interpolation-of-the-composite local terms are different
    discrete operators. Both are reported; neither is selected by its residual.
    """
    raw = {
        "R": R, "alpha": alpha, "gamma_inv": gamma_inv, "current": current,
        "Rr": Rr, "pr": pr, "Rdot": Rdot,
    }
    arrays = {name: np.asarray(value, dtype=float) for name, value in raw.items()}
    radii = arrays["R"]
    if radii.ndim != 1 or radii.size < 2:
        raise RuntimeError("composite fixed-R flux/work needs at least two radial cells")
    for name, values in arrays.items():
        if values.shape != radii.shape:
            raise RuntimeError(
                f"composite fixed-R {name} shape {values.shape} does not match {radii.shape}"
            )
        if not np.all(np.isfinite(values)):
            raise RuntimeError(f"composite fixed-R {name} contains non-finite values")
    if not np.all(np.diff(radii) > 0.0):
        raise RuntimeError("composite fixed-R areal radii must be strictly increasing")
    # Form the local physical spherical area before interpolation. Using the
    # target area at every cell would silently revert to product-of-interpolated
    # factors rather than interpolating the composite local flux/work profiles.
    area_factor = -FOUR_PI * np.square(radii)
    flux_profile = (
        area_factor * arrays["alpha"] * arrays["gamma_inv"]
        * arrays["current"] * arrays["Rr"]
    )
    pressure_work_profile = area_factor * arrays["pr"] * arrays["Rdot"]
    if not np.all(np.isfinite(flux_profile)) or not np.all(np.isfinite(pressure_work_profile)):
        raise RuntimeError("composite fixed-R flux/work profile is non-finite")
    return flux_profile, pressure_work_profile


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



def _sector_ktheta_source_rates(
    total, *, R, alpha, Ktheta, Rdot, target: float,
    live_primary_l3=None, vacuum_primary_l3=None,
) -> dict:
    """Audit live L3 matter sources and the signed fixed-R work identity.

    Report both interpolate-after-product and multiply-interpolated-factors
    conventions; these discrete operators need not agree. For zero shift,
    Rdot = -alpha*R*Ktheta and Ktheta_dot|matter = 4*pi*alpha*(rho/3+pr).
    With S=R**3*Ktheta*Ktheta_dot|matter, W=-4*pi*R**2*pr*Rdot,
    and D=(4*pi/3)*R**2*rho*Rdot, the consistent identity is
    S - W + D = 4*pi*R**2*(Rdot+alpha*R*Ktheta)*(rho/3+pr).

    This validates a source insertion, not a causal split of the complete
    Misner-Sharp tangent or a global energy-conservation law.
    """
    raw = {"R": R, "alpha": alpha, "Ktheta": Ktheta, "Rdot": Rdot}
    arrays = {name: np.asarray(value, dtype=float) for name, value in raw.items()}
    radii = arrays["R"]
    if radii.ndim != 1 or radii.size < 2:
        raise RuntimeError("Ktheta source audit R must be a 1-D array with at least two cells")
    for name, values in arrays.items():
        if values.shape != radii.shape:
            raise RuntimeError(
                f"Ktheta source audit {name} shape {values.shape} does not match {radii.shape}"
            )
        if not np.all(np.isfinite(values)):
            raise RuntimeError(f"Ktheta source audit {name} contains non-finite values")
    if not np.all(np.diff(radii) > 0.0):
        raise RuntimeError("Ktheta source audit R must be strictly increasing")
    target = float(target)
    if not math.isfinite(target) or target < radii[0] or target > radii[-1]:
        raise RuntimeError("Ktheta source audit target must be finite and inside the R domain")

    def _profile(name, value):
        out = np.asarray(value, dtype=float)
        if out.shape != radii.shape or not np.all(np.isfinite(out)):
            raise RuntimeError(f"Ktheta source audit {name} must be finite and match R")
        return out

    at = lambda profile: interpolate_at_areal_radius(radii, profile, target)
    alpha_values, Ktheta_values, Rdot_values = (
        arrays["alpha"], arrays["Ktheta"], arrays["Rdot"]
    )
    alpha_t, Ktheta_t, Rdot_t = at(alpha_values), at(Ktheta_values), at(Rdot_values)
    closure = Rdot_values + alpha_values * radii * Ktheta_values
    closure_t = at(closure)

    sector_rows = {}
    for label in ("scalar", "fluid", "electromagnetic"):
        raw_sector = {
            name: getattr(total, f"{label}_{name}", None)
            for name in ("rho", "pr", "pt")
        }
        if all(value is None for value in raw_sector.values()):
            present = False
            sector = {name: np.zeros_like(radii) for name in raw_sector}
        elif any(value is None for value in raw_sector.values()):
            raise RuntimeError(f"Ktheta source audit sector {label} has incomplete rho/pr/pt")
        else:
            present = True
            sector = {
                name: _profile(f"{label}_{name}", value)
                for name, value in raw_sector.items()
            }

        rho, pr, pt = sector["rho"], sector["pr"], sector["pt"]
        k_dot = FOUR_PI * alpha_values * (rho + pr + 2.0 * pt)
        aa_dot = -(16.0 * math.pi / 3.0) * alpha_values * (pr - pt)
        kt_dot = k_dot / 3.0 - aa_dot / 2.0
        source_profile = radii**3 * Ktheta_values * kt_dot
        work_profile = -FOUR_PI * radii**2 * pr * Rdot_values
        density_profile = (FOUR_PI / 3.0) * radii**2 * Rdot_values * rho
        identity_profile = source_profile - work_profile + density_profile
        identity_rhs_profile = FOUR_PI * radii**2 * closure * (rho / 3.0 + pr)

        rho_t, pr_t, pt_t = at(rho), at(pr), at(pt)
        factorized_kt_dot = FOUR_PI * alpha_t * (rho_t / 3.0 + pr_t)
        factorized_source = target**3 * Ktheta_t * factorized_kt_dot
        factorized_work = -FOUR_PI * target**2 * pr_t * Rdot_t
        factorized_density = (FOUR_PI / 3.0) * target**2 * Rdot_t * rho_t
        factorized_identity = factorized_source - factorized_work + factorized_density
        factorized_closure = Rdot_t + alpha_t * target * Ktheta_t
        factorized_identity_rhs = FOUR_PI * target**2 * factorized_closure * (rho_t / 3.0 + pr_t)
        profiles = {
            "rho": rho, "pr": pr, "pt": pt,
            "k_dot_matter_source": k_dot,
            "aa_dot_matter_source": aa_dot,
            "ktheta_dot_matter_source": kt_dot,
            "mass_dot_coordinate_from_ktheta_source": source_profile,
            "radial_pressure_work_composite": work_profile,
            "density_correction_composite": density_profile,
            "source_work_density_identity_residual_composite": identity_profile,
            "source_work_density_identity_kinematic_rhs_composite": identity_rhs_profile,
            "source_work_density_identity_algebra_residual_composite": (
                identity_profile - identity_rhs_profile
            ),
        }
        sector_rows[label] = {
            "present": bool(present),
            **{f"{name}_at_R_star": at(values) for name, values in profiles.items()},
            "mass_dot_coordinate_from_ktheta_source_factorized_at_R_star": float(factorized_source),
            "radial_pressure_work_factorized_at_R_star": float(factorized_work),
            "density_correction_factorized_at_R_star": float(factorized_density),
            "source_work_density_identity_residual_factorized_at_R_star": float(factorized_identity),
            "source_work_density_identity_kinematic_rhs_factorized_at_R_star": float(factorized_identity_rhs),
            "source_work_density_identity_algebra_residual_factorized_at_R_star": float(
                factorized_identity - factorized_identity_rhs
            ),
        }

    total_profiles = {
        name: _profile(f"total_{name}", getattr(total, name))
        for name in ("rho", "pr", "pt")
    }
    rho, pr, pt = total_profiles["rho"], total_profiles["pr"], total_profiles["pt"]
    total_k_dot = FOUR_PI * alpha_values * (rho + pr + 2.0 * pt)
    total_aa_dot = -(16.0 * math.pi / 3.0) * alpha_values * (pr - pt)
    total_kt_dot = total_k_dot / 3.0 - total_aa_dot / 2.0
    aggregate_profiles = {
        "k_dot_matter_source": total_k_dot,
        "aa_dot_matter_source": total_aa_dot,
        "ktheta_dot_matter_source": total_kt_dot,
        "mass_dot_coordinate_from_ktheta_source": radii**3 * Ktheta_values * total_kt_dot,
    }
    aggregate_at = {name: at(values) for name, values in aggregate_profiles.items()}
    sums, additive_residuals, additive_bounds = {}, {}, {}
    for name in aggregate_profiles:
        value = sum(row[f"{name}_at_R_star"] for row in sector_rows.values())
        residual = value - aggregate_at[name]
        scale = max(
            abs(value) + abs(aggregate_at[name]),
            sum(abs(row[f"{name}_at_R_star"]) for row in sector_rows.values()),
            float(np.finfo(np.float64).tiny),
        )
        bound = 64.0 * np.finfo(np.float64).eps * scale
        if abs(residual) > bound:
            raise RuntimeError(
                f"Ktheta source audit failed additive {name}: "
                f"residual={residual:.17g} (bound {bound:.17g})"
            )
        sums[name], additive_residuals[name], additive_bounds[name] = (
            float(value), float(residual), float(bound)
        )

    live_validation = {}
    if (live_primary_l3 is None) != (vacuum_primary_l3 is None):
        raise RuntimeError("Ktheta source audit requires both live and vacuum L3 RHS for source validation")
    if live_primary_l3 is not None:
        for field, expected in (("K", total_k_dot), ("Aa", total_aa_dot)):
            try:
                live = _profile(f"live_primary_l3_{field}", live_primary_l3[field])
                vacuum = _profile(f"vacuum_primary_l3_{field}", vacuum_primary_l3[field])
            except (KeyError, TypeError) as exc:
                raise RuntimeError(f"Ktheta source audit missing captured L3 field {field}") from exc
            observed = live - vacuum
            residual = observed - expected
            bound = 64.0 * np.finfo(np.float64).eps * np.maximum(
                np.abs(live) + np.abs(vacuum) + np.abs(expected),
                float(np.finfo(np.float64).tiny),
            )
            ratio_max = float(np.max(np.abs(residual) / bound))
            if not math.isfinite(ratio_max) or ratio_max > 1.0:
                raise RuntimeError(
                    f"live L3 {field} matter source mismatch: max bound ratio={ratio_max:.17g}"
                )
            live_validation[field] = {
                "live_minus_vacuum_at_R_star": float(at(observed)),
                "reconstructed_matter_source_at_R_star": float(at(expected)),
                "residual_at_R_star": float(at(residual)),
                "bound_at_R_star": float(at(bound)),
                "residual_absmax": float(np.max(np.abs(residual))),
                "bound_ratio_max": ratio_max,
            }

    total_work = -FOUR_PI * radii**2 * pr * Rdot_values
    total_density = (FOUR_PI / 3.0) * radii**2 * Rdot_values * rho
    total_identity = radii**3 * Ktheta_values * total_kt_dot - total_work + total_density
    total_identity_rhs = FOUR_PI * radii**2 * closure * (rho / 3.0 + pr)
    rho_t, pr_t, pt_t = at(rho), at(pr), at(pt)
    factorized_kt_dot = FOUR_PI * alpha_t * (rho_t / 3.0 + pr_t)
    factorized_source = target**3 * Ktheta_t * factorized_kt_dot
    factorized_work = -FOUR_PI * target**2 * pr_t * Rdot_t
    factorized_density = (FOUR_PI / 3.0) * target**2 * Rdot_t * rho_t
    factorized_identity = factorized_source - factorized_work + factorized_density
    factorized_closure = Rdot_t + alpha_t * target * Ktheta_t
    factorized_identity_rhs = FOUR_PI * target**2 * factorized_closure * (rho_t / 3.0 + pr_t)

    return {
        "sectors": sector_rows,
        "sum_at_R_star": sums,
        "aggregate_at_R_star": aggregate_at,
        "additive_residuals": additive_residuals,
        "additive_bounds": additive_bounds,
        "kinematic_closure_Rdot_plus_alpha_R_Ktheta_at_R_star": float(closure_t),
        "factorized_kinematic_closure_Rdot_plus_alpha_R_Ktheta_at_R_star": float(factorized_closure),
        "kinematic_closure_absmax": float(np.max(np.abs(closure))),
        "total_mass_source_composite_at_R_star": float(at(aggregate_profiles["mass_dot_coordinate_from_ktheta_source"])),
        "total_mass_source_factorized_at_R_star": float(factorized_source),
        "total_radial_pressure_work_composite_at_R_star": float(at(total_work)),
        "total_radial_pressure_work_factorized_at_R_star": float(factorized_work),
        "total_density_correction_composite_at_R_star": float(at(total_density)),
        "total_density_correction_factorized_at_R_star": float(factorized_density),
        "total_source_work_density_identity_residual_composite_at_R_star": float(at(total_identity)),
        "total_source_work_density_identity_kinematic_rhs_composite_at_R_star": float(at(total_identity_rhs)),
        "total_source_work_density_identity_algebra_residual_composite_at_R_star": float(at(total_identity - total_identity_rhs)),
        "total_source_work_density_identity_residual_factorized_at_R_star": float(factorized_identity),
        "total_source_work_density_identity_kinematic_rhs_factorized_at_R_star": float(factorized_identity_rhs),
        "total_source_work_density_identity_algebra_residual_factorized_at_R_star": float(
            factorized_identity - factorized_identity_rhs
        ),
        "live_l3_source_validation": live_validation,
    }

def _sector_fixed_R_rates(
    total,
    *,
    R,
    alpha,
    gamma_rr_inv,
    Rr,
    Rdot,
    target: float,
) -> dict:
    """Signed fixed-R stress-flux and pressure-work by TotalStressEnergy sector.

    Sectors are the three additive components already exposed by the canonical
    projection (scalar, fluid, electromagnetic). The sum of the three sector
    rates is required to recover the aggregate rate within a 64-epsilon
    scale-aware bound; any mismatch fails closed.
    """
    at = lambda profile: interpolate_at_areal_radius(R, profile, target)
    area = -FOUR_PI * target ** 2
    alpha_t = at(alpha)
    gamma_inv_t = at(gamma_rr_inv)
    Rr_t = at(Rr)
    Rdot_t = at(Rdot)

    def _safe_arr(arr):
        if arr is None:
            return np.zeros_like(R, dtype=float)
        out = np.asarray(arr, dtype=float)
        if out.shape != R.shape:
            raise RuntimeError("sector array shape does not match areal radius")
        if not np.all(np.isfinite(out)):
            raise RuntimeError("sector array contains non-finite values")
        return out

    def _one(j_arr, pr_arr, label: str) -> dict:
        present = j_arr is not None and pr_arr is not None
        j = _safe_arr(j_arr)
        pr = _safe_arr(pr_arr)
        j_t = at(j)
        pr_t = at(pr)
        flux = area * alpha_t * gamma_inv_t * j_t * Rr_t
        work = area * pr_t * Rdot_t
        if not (math.isfinite(flux) and math.isfinite(work)):
            raise RuntimeError(f"sector {label}: non-finite flux or work")
        return {
            "stress_flux": float(flux),
            "radial_pressure_work": float(work),
            "present": bool(present),
        }

    sectors = {
        "scalar": _one(total.scalar_j, total.scalar_pr, "scalar"),
        "fluid": _one(total.fluid_j, total.fluid_pr, "fluid"),
        "electromagnetic": _one(
            total.electromagnetic_j, total.electromagnetic_pr, "electromagnetic"
        ),
    }

    sum_flux = sum(s["stress_flux"] for s in sectors.values())
    sum_work = sum(s["radial_pressure_work"] for s in sectors.values())
    agg_j = at(_safe_arr(total.j))
    agg_pr = at(_safe_arr(total.pr))
    agg_flux = area * alpha_t * gamma_inv_t * agg_j * Rr_t
    agg_work = area * agg_pr * Rdot_t

    scale_flux = max(
        abs(agg_flux) + abs(sum_flux),
        abs(area) * abs(alpha_t) * abs(gamma_inv_t) * abs(Rr_t) * (
            abs(agg_j)
            + abs(at(_safe_arr(total.scalar_j)))
            + abs(at(_safe_arr(total.fluid_j)))
            + abs(at(_safe_arr(total.electromagnetic_j)))
        ),
        float(np.finfo(np.float64).tiny),
    )
    scale_work = max(
        abs(agg_work) + abs(sum_work),
        abs(area) * abs(Rdot_t) * (
            abs(agg_pr)
            + abs(at(_safe_arr(total.scalar_pr)))
            + abs(at(_safe_arr(total.fluid_pr)))
            + abs(at(_safe_arr(total.electromagnetic_pr)))
        ),
        float(np.finfo(np.float64).tiny),
    )
    bound_flux = 64.0 * np.finfo(np.float64).eps * scale_flux
    bound_work = 64.0 * np.finfo(np.float64).eps * scale_work
    residual_flux = sum_flux - agg_flux
    residual_work = sum_work - agg_work
    if abs(residual_flux) > bound_flux or abs(residual_work) > bound_work:
        raise RuntimeError(
            "sector fixed-R rates failed additive consistency: "
            f"flux residual={residual_flux:.17g} (bound {bound_flux:.17g}), "
            f"work residual={residual_work:.17g} (bound {bound_work:.17g})"
        )

    return {
        "sectors": sectors,
        "sum_stress_flux": float(sum_flux),
        "sum_radial_pressure_work": float(sum_work),
        "aggregate_stress_flux": float(agg_flux),
        "aggregate_radial_pressure_work": float(agg_work),
        "additive_flux_residual": float(residual_flux),
        "additive_work_residual": float(residual_work),
        "additive_flux_bound": float(bound_flux),
        "additive_work_bound": float(bound_work),
    }


def _explicit_metric_rate_profiles(grid, R, a, b, X, explicit_rhs) -> dict[str, np.ndarray]:
    """Reconstruct explicit metric-rate profiles with the production parity stencil.

    This is an observer-only path: it never advances or modifies geometry.
    """
    raw = {
        "R": R, "a": a, "b": b, "X": X,
        "a_dot": explicit_rhs["a"], "b_dot": explicit_rhs["b"], "X_dot": explicit_rhs["X"],
    }
    arrays = {name: np.asarray(value, dtype=float) for name, value in raw.items()}
    shape = arrays["R"].shape
    if arrays["R"].ndim != 1 or arrays["R"].size == 0:
        raise RuntimeError("explicit metric-rate R must be a nonempty 1-D profile")
    for name, value in arrays.items():
        if value.shape != shape:
            raise RuntimeError(f"explicit metric-rate {name} shape {value.shape} does not match {shape}")
        if not np.all(np.isfinite(value)):
            raise RuntimeError(f"explicit metric-rate {name} contains non-finite values")
    for name in ("a", "b", "X"):
        if np.any(arrays[name] == 0.0):
            raise RuntimeError(f"explicit metric-rate denominator {name} contains zero")
    Rdot = arrays["R"] * (
        0.5 * arrays["b_dot"] / arrays["b"] - arrays["X_dot"] / arrays["X"]
    )
    Rdot_r = np.asarray(grid.cell_derivative_fourth(Rdot, parity=-1), dtype=float)
    gamma_inv = arrays["X"]**2 / arrays["a"]
    gamma_inv_dot = gamma_inv * (
        2.0 * arrays["X_dot"] / arrays["X"] - arrays["a_dot"] / arrays["a"]
    )
    result = {
        "Rdot": Rdot, "Rdot_r": Rdot_r,
        "gamma_inv": gamma_inv, "gamma_inv_dot": gamma_inv_dot,
    }
    for name, value in result.items():
        if value.shape != shape or not np.all(np.isfinite(value)):
            raise RuntimeError(f"explicit metric-rate reconstruction {name} is invalid")
    return result


def _profile_roundoff_reconstruction(actual, reconstructed, *, label: str) -> dict[str, np.ndarray]:
    """Compare two paths with a term-scaled float64 roundoff bound (no absolute floor)."""
    left = np.asarray(actual, dtype=float)
    right = np.asarray(reconstructed, dtype=float)
    if left.ndim != 1 or right.shape != left.shape:
        raise RuntimeError(f"{label}: profiles have incompatible shapes")
    if not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
        raise RuntimeError(f"{label}: profiles contain non-finite values")
    residual = left - right
    scale = np.abs(left) + np.abs(right)
    bound = 64.0 * np.finfo(np.float64).eps * np.maximum(
        scale, np.finfo(np.float64).tiny
    )
    if np.any(np.abs(residual) > bound):
        raise RuntimeError(f"{label}: failed scale-aware roundoff reconstruction")
    return {"residual": residual, "roundoff_bound": bound}


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

    # Independently expose the five terms in the expanded derivative
    # Mdot = 1/2 Rdot(1-chi) - 1/2 R gamma_inv_dot Rr^2
    #        - R gamma_inv Rr Rdot_r + R^2 Rdot Ktheta^2
    #        + R^3 Ktheta Ktheta_dot.
    # These are algebraic witnesses only: the physical flux/work remainder is
    # deliberately not thresholded or altered by this reconstruction.
    term_Rdot = 0.5 * Rdot * (1.0 - chi)
    term_gamma_inv_dot = -0.5 * R * gamma_inv_dot * Rr**2
    term_Rdot_r = -R * gamma_inv * Rr * Rdot_r
    term_Ktheta_geometry = R**2 * Rdot * Ktheta**2
    term_Ktheta_dot = R**3 * Ktheta * Ktheta_dot
    term_sum = (
        term_Rdot
        + term_gamma_inv_dot
        + term_Rdot_r
        + term_Ktheta_geometry
        + term_Ktheta_dot
    )
    term_residual = mass_dot_coordinate - term_sum
    term_scale = (
        np.abs(term_Rdot)
        + np.abs(term_gamma_inv_dot)
        + np.abs(term_Rdot_r)
        + np.abs(term_Ktheta_geometry)
        + np.abs(term_Ktheta_dot)
        + np.abs(mass_dot_coordinate)
    )
    term_bound = 64.0 * np.finfo(np.float64).eps * np.maximum(
        term_scale, np.finfo(np.float64).tiny
    )
    outputs = {
        "chi": chi,
        "chi_dot": chi_dot,
        "mass": mass,
        "mass_dot_coordinate": mass_dot_coordinate,
        "mass_dot_coordinate_term_Rdot": term_Rdot,
        "mass_dot_coordinate_term_gamma_inv_dot": term_gamma_inv_dot,
        "mass_dot_coordinate_term_Rdot_r": term_Rdot_r,
        "mass_dot_coordinate_term_Ktheta_geometry": term_Ktheta_geometry,
        "mass_dot_coordinate_term_Ktheta_dot": term_Ktheta_dot,
        "mass_dot_coordinate_term_sum": term_sum,
        "mass_dot_coordinate_decomposition_residual": term_residual,
    }
    if not all(np.all(np.isfinite(value)) for value in outputs.values()):
        raise RuntimeError("Misner-Sharp chain-rule evaluation produced non-finite output")
    if np.any(np.abs(term_residual) > term_bound):
        raise RuntimeError(
            "Misner-Sharp expanded mass-rate terms failed roundoff reconstruction"
        )
    return outputs



def _ktheta_dot_component_profiles(
    primary_l2_K, primary_l2_Aa, vacuum_l3_K, vacuum_l3_Aa,
    einstein_source_l3_K, einstein_source_l3_Aa, *,
    expected_primary_l2, expected_vacuum_l3, expected_einstein_source_l3,
) -> dict[str, np.ndarray]:
    """Separate K/3 and -Aa/2 operands before recombination for each RHS lane."""
    raw = {
        "primary_l2_K": primary_l2_K,
        "primary_l2_Aa": primary_l2_Aa,
        "vacuum_l3_K": vacuum_l3_K,
        "vacuum_l3_Aa": vacuum_l3_Aa,
        "einstein_source_l3_K": einstein_source_l3_K,
        "einstein_source_l3_Aa": einstein_source_l3_Aa,
        "expected_primary_l2": expected_primary_l2,
        "expected_vacuum_l3": expected_vacuum_l3,
        "expected_einstein_source_l3": expected_einstein_source_l3,
    }
    arrays = {name: np.asarray(value, dtype=float) for name, value in raw.items()}
    first = arrays["primary_l2_K"]
    if first.ndim != 1 or first.size == 0:
        raise RuntimeError("Ktheta-dot component inputs must be nonempty 1-D profiles")
    for name, values in arrays.items():
        if values.shape != first.shape:
            raise RuntimeError(
                f"Ktheta-dot component {name} shape {values.shape} does not match {first.shape}"
            )
        if not np.all(np.isfinite(values)):
            raise RuntimeError(f"Ktheta-dot component {name} contains non-finite values")
    components = {
        "primary_l2_from_K": arrays["primary_l2_K"] / 3.0,
        "primary_l2_from_Aa": -arrays["primary_l2_Aa"] / 2.0,
        "vacuum_l3_from_K": arrays["vacuum_l3_K"] / 3.0,
        "vacuum_l3_from_Aa": -arrays["vacuum_l3_Aa"] / 2.0,
        "einstein_source_l3_from_K": arrays["einstein_source_l3_K"] / 3.0,
        "einstein_source_l3_from_Aa": -arrays["einstein_source_l3_Aa"] / 2.0,
    }
    groups = {
        "primary_l2": (
            "primary_l2_from_K", "primary_l2_from_Aa", arrays["expected_primary_l2"]
        ),
        "vacuum_l3": (
            "vacuum_l3_from_K", "vacuum_l3_from_Aa", arrays["expected_vacuum_l3"]
        ),
        "einstein_source_l3": (
            "einstein_source_l3_from_K", "einstein_source_l3_from_Aa",
            arrays["expected_einstein_source_l3"],
        ),
    }
    diagnostics = {}
    for group, (from_K, from_Aa, expected) in groups.items():
        reconstructed = components[from_K] + components[from_Aa]
        residual = reconstructed - expected
        scale = np.abs(components[from_K]) + np.abs(components[from_Aa]) + np.abs(expected)
        bound = 64.0 * np.finfo(np.float64).eps * np.maximum(
            scale, np.finfo(np.float64).tiny
        )
        if np.any(np.abs(residual) > bound):
            raise RuntimeError(f"Ktheta-dot {group} K/Aa split failed roundoff reconstruction")
        diagnostics[f"{group}_reconstruction_residual"] = residual
        diagnostics[f"{group}_reconstruction_bound"] = bound
    return {**components, **diagnostics}


def _misner_sharp_chain_rule_bracket_operands(
    R, Rdot, Rr, Rdot_r, gamma_inv, gamma_inv_dot, Ktheta, Ktheta_dot,
    Ktheta_dot_primary_l2, Ktheta_dot_vacuum_l3, Ktheta_dot_matter_source,
    *, target: float, ktheta_dot_component_profiles: dict[str, np.ndarray] | None = None,
    metric_rate_reconstruction_profiles: dict[str, np.ndarray] | None = None
) -> dict:
    """Retain signed stage operands and source/non-source Ktheta-dot terms at R*'s bracket cells."""
    raw = {
        "R": R,
        "Rdot": Rdot,
        "Rr": Rr,
        "Rdot_r": Rdot_r,
        "gamma_inv": gamma_inv,
        "gamma_inv_dot": gamma_inv_dot,
        "Ktheta": Ktheta,
        "Ktheta_dot": Ktheta_dot,
        "Ktheta_dot_primary_l2": Ktheta_dot_primary_l2,
        "Ktheta_dot_vacuum_l3": Ktheta_dot_vacuum_l3,
        "Ktheta_dot_matter_source": Ktheta_dot_matter_source,
    }
    arrays = {name: np.asarray(value, dtype=float) for name, value in raw.items()}
    radii = arrays["R"]
    if radii.ndim != 1 or radii.size < 2:
        raise RuntimeError("chain-rule bracket witness R must be 1-D with at least two cells")
    for name, values in arrays.items():
        if values.shape != radii.shape:
            raise RuntimeError(
                f"chain-rule bracket witness {name} shape {values.shape} does not match {radii.shape}"
            )
        if not np.all(np.isfinite(values)):
            raise RuntimeError(f"chain-rule bracket witness {name} contains non-finite values")
    if not np.all(np.diff(radii) > 0.0):
        raise RuntimeError("chain-rule bracket witness R must be strictly increasing")
    target = float(target)
    if not math.isfinite(target) or target < radii[0] or target > radii[-1]:
        raise RuntimeError("chain-rule bracket witness target is outside the sampled domain")

    right = int(np.searchsorted(radii, target, side="right"))
    right = min(max(right, 1), radii.size - 1)
    left = right - 1
    theta = float((target - radii[left]) / (radii[right] - radii[left]))

    normal_dR = radii * arrays["Ktheta"]
    chi = arrays["gamma_inv"] * arrays["Rr"]**2 - normal_dR**2
    term_Rdot = 0.5 * arrays["Rdot"] * (1.0 - chi)
    term_gamma_inv_dot = (
        -0.5 * radii * arrays["gamma_inv_dot"] * arrays["Rr"]**2
    )
    term_Rdot_r = (
        -radii * arrays["gamma_inv"] * arrays["Rr"] * arrays["Rdot_r"]
    )
    term_Ktheta_geometry = (
        radii**2 * arrays["Rdot"] * arrays["Ktheta"]**2
    )
    term_Ktheta_dot_primary_l2 = (
        radii**3 * arrays["Ktheta"] * arrays["Ktheta_dot_primary_l2"]
    )
    term_Ktheta_dot_vacuum_l3 = (
        radii**3 * arrays["Ktheta"] * arrays["Ktheta_dot_vacuum_l3"]
    )
    term_Ktheta_dot_matter_source = (
        radii**3 * arrays["Ktheta"] * arrays["Ktheta_dot_matter_source"]
    )
    ktheta_dot_from_captured_pieces = (
        arrays["Ktheta_dot_primary_l2"]
        + arrays["Ktheta_dot_vacuum_l3"]
        + arrays["Ktheta_dot_matter_source"]
    )
    ktheta_dot_reconstruction_residual = (
        arrays["Ktheta_dot"] - ktheta_dot_from_captured_pieces
    )
    ktheta_dot_reconstruction_scale = (
        np.abs(arrays["Ktheta_dot"])
        + np.abs(arrays["Ktheta_dot_primary_l2"])
        + np.abs(arrays["Ktheta_dot_vacuum_l3"])
        + np.abs(arrays["Ktheta_dot_matter_source"])
    )
    ktheta_dot_reconstruction_bound = (
        64.0 * np.finfo(np.float64).eps
        * np.maximum(ktheta_dot_reconstruction_scale, np.finfo(np.float64).tiny)
    )
    if np.any(np.abs(ktheta_dot_reconstruction_residual) > ktheta_dot_reconstruction_bound):
        raise RuntimeError(
            "captured primary-L2, vacuum-L3 and Einstein source Ktheta-dot pieces "
            "failed roundoff reconstruction"
        )
    term_Ktheta_dot_non_matter = (
        radii**3 * arrays["Ktheta"]
        * (arrays["Ktheta_dot"] - arrays["Ktheta_dot_matter_source"])
    )
    term_Ktheta_dot_total = radii**3 * arrays["Ktheta"] * arrays["Ktheta_dot"]
    component_profiles = None
    if ktheta_dot_component_profiles is not None:
        component_profiles = {}
        for name, value in ktheta_dot_component_profiles.items():
            array = np.asarray(value, dtype=float)
            if array.shape != radii.shape or not np.all(np.isfinite(array)):
                raise RuntimeError(
                    f"chain-rule bracket witness component {name} must be finite and match R"
                )
            component_profiles[name] = array
        required = {
            "primary_l2_from_K", "primary_l2_from_Aa",
            "vacuum_l3_from_K", "vacuum_l3_from_Aa",
            "einstein_source_l3_from_K", "einstein_source_l3_from_Aa",
            "primary_l2_reconstruction_residual", "primary_l2_reconstruction_bound",
            "vacuum_l3_reconstruction_residual", "vacuum_l3_reconstruction_bound",
            "einstein_source_l3_reconstruction_residual", "einstein_source_l3_reconstruction_bound",
        }
        if set(component_profiles) != required:
            raise RuntimeError("chain-rule bracket witness received an incomplete K/Aa split")
    metric_rate_profiles = None
    if metric_rate_reconstruction_profiles is not None:
        metric_rate_profiles = {}
        for name, value in metric_rate_reconstruction_profiles.items():
            array = np.asarray(value, dtype=float)
            if array.shape != radii.shape or not np.all(np.isfinite(array)):
                raise RuntimeError(
                    f"chain-rule bracket witness metric-rate {name} must be finite and match R"
                )
            metric_rate_profiles[name] = array
        required_metric = {
            "Rdot_residual", "Rdot_roundoff_bound",
            "gamma_inv_dot_residual", "gamma_inv_dot_roundoff_bound",
            "Rdot_r_residual", "Rdot_r_roundoff_bound",
        }
        if set(metric_rate_profiles) != required_metric:
            raise RuntimeError("chain-rule bracket witness received incomplete metric-rate witnesses")

    term_sum = (
        term_Rdot + term_gamma_inv_dot + term_Rdot_r
        + term_Ktheta_geometry + term_Ktheta_dot_total
    )
    chi_dot = (
        arrays["gamma_inv_dot"] * arrays["Rr"]**2
        + 2.0 * arrays["gamma_inv"] * arrays["Rr"] * arrays["Rdot_r"]
        - 2.0 * normal_dR * (
            arrays["Rdot"] * arrays["Ktheta"]
            + radii * arrays["Ktheta_dot"]
        )
    )
    mass_dot_coordinate = 0.5 * (
        arrays["Rdot"] * (1.0 - chi) - radii * chi_dot
    )
    residual = mass_dot_coordinate - term_sum
    scale = (
        np.abs(term_Rdot) + np.abs(term_gamma_inv_dot) + np.abs(term_Rdot_r)
        + np.abs(term_Ktheta_geometry) + np.abs(term_Ktheta_dot_matter_source)
        + np.abs(term_Ktheta_dot_non_matter) + np.abs(mass_dot_coordinate)
    )
    bound = 64.0 * np.finfo(np.float64).eps * np.maximum(
        scale, np.finfo(np.float64).tiny
    )
    if np.any(np.abs(residual) > bound):
        raise RuntimeError("chain-rule bracket witness failed roundoff decomposition")

    def cell(index):
        return {
            "index": int(index),
            "R": float(radii[index]),
            "Rdot_from_production_explicit_rhs": float(arrays["Rdot"][index]),
            "gamma_inv_dot_from_production_explicit_rhs": float(
                arrays["gamma_inv_dot"][index]
            ),
            "Rdot_r_from_cell_derivative_fourth_parity_minus_one": float(
                arrays["Rdot_r"][index]
            ),
            "Ktheta_dot_total_from_primary_L2_plus_L3": float(
                arrays["Ktheta_dot"][index]
            ),
            "Ktheta_dot_primary_l2_from_captured_production_L2": float(
                arrays["Ktheta_dot_primary_l2"][index]
            ),
            "Ktheta_dot_vacuum_l3_from_captured_vacuum_L3": float(
                arrays["Ktheta_dot_vacuum_l3"][index]
            ),
            "Ktheta_dot_matter_source_from_same_stage_live_minus_vacuum_L3": float(
                arrays["Ktheta_dot_matter_source"][index]
            ),
            "Ktheta_dot_non_matter_remainder": float(
                arrays["Ktheta_dot"][index]
                - arrays["Ktheta_dot_matter_source"][index]
            ),
            "Ktheta_dot_non_matter_from_L2_plus_vacuum_L3": float(
                arrays["Ktheta_dot_primary_l2"][index]
                + arrays["Ktheta_dot_vacuum_l3"][index]
            ),
            "Ktheta_dot_reconstruction_residual_total_minus_L2_vacuum_L3_and_source": float(
                ktheta_dot_reconstruction_residual[index]
            ),
            "Ktheta_dot_reconstruction_roundoff_bound": float(
                ktheta_dot_reconstruction_bound[index]
            ),
            "mass_dot_coordinate_term_Ktheta_dot_primary_l2": float(
                term_Ktheta_dot_primary_l2[index]
            ),
            "mass_dot_coordinate_term_Ktheta_dot_vacuum_l3": float(
                term_Ktheta_dot_vacuum_l3[index]
            ),
            "mass_dot_coordinate_term_Rdot": float(term_Rdot[index]),
            "mass_dot_coordinate_term_gamma_inv_dot": float(
                term_gamma_inv_dot[index]
            ),
            "mass_dot_coordinate_term_Rdot_r": float(term_Rdot_r[index]),
            "mass_dot_coordinate_term_Ktheta_geometry": float(
                term_Ktheta_geometry[index]
            ),
            "mass_dot_coordinate_term_Ktheta_dot_matter_source": float(
                term_Ktheta_dot_matter_source[index]
            ),
            "mass_dot_coordinate_term_Ktheta_dot_non_matter": float(
                term_Ktheta_dot_non_matter[index]
            ),
            "mass_dot_coordinate_term_Ktheta_dot_total": float(
                term_Ktheta_dot_total[index]
            ),
            "mass_dot_coordinate_from_chain_rule": float(
                mass_dot_coordinate[index]
            ),
            "mass_dot_coordinate_term_sum": float(term_sum[index]),
            "mass_dot_coordinate_decomposition_residual": float(residual[index]),
            "mass_dot_coordinate_decomposition_bound": float(bound[index]),
            "ktheta_dot_operand_split": (
                None if component_profiles is None else {
                    name: float(component_profiles[name][index])
                    for name in (
                        "primary_l2_from_K", "primary_l2_from_Aa",
                        "vacuum_l3_from_K", "vacuum_l3_from_Aa",
                        "einstein_source_l3_from_K", "einstein_source_l3_from_Aa",
                    )
                }
            ),
            "mass_dot_coordinate_term_ktheta_dot_operand_split": (
                None if component_profiles is None else {
                    name: float(
                        radii[index]**3 * arrays["Ktheta"][index]
                        * component_profiles[name][index]
                    )
                    for name in (
                        "primary_l2_from_K", "primary_l2_from_Aa",
                        "vacuum_l3_from_K", "vacuum_l3_from_Aa",
                        "einstein_source_l3_from_K", "einstein_source_l3_from_Aa",
                    )
                }
            ),
            "ktheta_dot_component_reconstruction": (
                None if component_profiles is None else {
                    group: {
                        "residual": float(component_profiles[f"{group}_reconstruction_residual"][index]),
                        "roundoff_bound": float(component_profiles[f"{group}_reconstruction_bound"][index]),
                    }
                    for group in ("primary_l2", "vacuum_l3", "einstein_source_l3")
                }
            ),
            "explicit_metric_rate_reconstruction": (
                None if metric_rate_profiles is None else {
                    "Rdot_residual": float(metric_rate_profiles["Rdot_residual"][index]),
                    "Rdot_roundoff_bound": float(metric_rate_profiles["Rdot_roundoff_bound"][index]),
                    "gamma_inv_dot_residual": float(metric_rate_profiles["gamma_inv_dot_residual"][index]),
                    "gamma_inv_dot_roundoff_bound": float(metric_rate_profiles["gamma_inv_dot_roundoff_bound"][index]),
                    "Rdot_r_residual": float(metric_rate_profiles["Rdot_r_residual"][index]),
                    "Rdot_r_roundoff_bound": float(metric_rate_profiles["Rdot_r_roundoff_bound"][index]),
                }
            ),
        }

    return {
        "target_R": target,
        "left_index": int(left),
        "right_index": int(right),
        "theta": theta,
        "left_cell": cell(left),
        "right_cell": cell(right),
    }

def _mass_gradient_identity_bracket_witness(
    R, mass_gradient, rho, current, Rr, Rdot, alpha,
    hamiltonian, momentum, *, target: float
) -> dict:
    """Expose the native mass-gradient identity at both cells bracketing R=target.

    The cell-local composite RHS is formed before interpolation:
      D_r M = 4*pi*R**2*(rho*D_r R + j*Rdot/alpha).
    The existing tube scalar instead interpolates each factor and then forms
    the RHS at the target radius. Both operators are recorded side-by-side;
    this helper measures their difference but never thresholds the physical
    identity residual or chooses an operator based on which closes better.
    """
    raw = {
        "R": R, "mass_gradient": mass_gradient, "rho": rho,
        "current": current, "Rr": Rr, "Rdot": Rdot, "alpha": alpha,
        "hamiltonian": hamiltonian, "momentum": momentum,
    }
    arrays = {name: np.asarray(value, dtype=float) for name, value in raw.items()}
    radii = arrays["R"]
    if radii.ndim != 1 or radii.size < 2:
        raise RuntimeError("mass-gradient witness R must be a 1-D array with at least two cells")
    for name, values in arrays.items():
        if values.shape != radii.shape:
            raise RuntimeError(
                f"mass-gradient witness {name} shape {values.shape} does not match {radii.shape}"
            )
        if not np.all(np.isfinite(values)):
            raise RuntimeError(f"mass-gradient witness {name} contains non-finite values")
    if not np.all(np.diff(radii) > 0.0):
        raise RuntimeError("mass-gradient witness R must be strictly increasing")
    if np.any(arrays["alpha"] == 0.0):
        raise RuntimeError("mass-gradient witness alpha contains zero")
    target = float(target)
    if not math.isfinite(target) or target < radii[0] or target > radii[-1]:
        raise RuntimeError("mass-gradient witness target is outside the sampled areal-radius domain")

    local_rhs = FOUR_PI * np.square(radii) * (
        arrays["rho"] * arrays["Rr"]
        + arrays["current"] * arrays["Rdot"] / arrays["alpha"]
    )
    residual = arrays["mass_gradient"] - local_rhs
    right = int(np.searchsorted(radii, target, side="right"))
    if right >= radii.size:
        right = radii.size - 1
    if right == 0:
        right = 1
    left = right - 1
    theta = float((target - radii[left]) / (radii[right] - radii[left]))

    def interp(values):
        return float((1.0 - theta) * values[left] + theta * values[right])

    rho_at = interp(arrays["rho"])
    current_at = interp(arrays["current"])
    Rr_at = interp(arrays["Rr"])
    Rdot_at = interp(arrays["Rdot"])
    alpha_at = interp(arrays["alpha"])
    factor_rhs = FOUR_PI * target**2 * (
        rho_at * Rr_at + current_at * Rdot_at / alpha_at
    )

    # Candidate continuum constraint identity, measured rather than enforced:
    # E_MS = R^2 R_r H/4 + R^2 U M_r/2, U = Rdot/alpha.
    # Keep the operator remainder visible; discrete product rules/operators
    # need not close pointwise and this physical residual is never thresholded.
    normal_velocity = arrays["Rdot"] / arrays["alpha"]
    hamiltonian_contribution = (
        np.square(radii) * arrays["Rr"] * arrays["hamiltonian"] / 4.0
    )
    momentum_contribution = (
        np.square(radii) * normal_velocity * arrays["momentum"] / 2.0
    )
    constraint_predicted_residual = hamiltonian_contribution + momentum_contribution
    constraint_identity_remainder = residual - constraint_predicted_residual

    def cell(index):
        return {
            "index": int(index),
            "normal_velocity_Rdot_over_alpha": float(normal_velocity[index]),
            "hamiltonian_constraint_contribution_to_mass_gradient_residual": float(
                hamiltonian_contribution[index]
            ),
            "momentum_constraint_contribution_to_mass_gradient_residual": float(
                momentum_contribution[index]
            ),
            "constraint_predicted_mass_gradient_residual": float(
                constraint_predicted_residual[index]
            ),
            "constraint_identity_operator_remainder": float(
                constraint_identity_remainder[index]
            ),
            "R": float(radii[index]),
            "mass_gradient_discrete": float(arrays["mass_gradient"][index]),
            "mass_gradient_rhs_local_composite": float(local_rhs[index]),
            "mass_gradient_residual_local_composite": float(residual[index]),
            "rho": float(arrays["rho"][index]),
            "j": float(arrays["current"][index]),
            "Rr": float(arrays["Rr"][index]),
            "Rdot": float(arrays["Rdot"][index]),
            "alpha": float(arrays["alpha"][index]),
            "hamiltonian_constraint": float(arrays["hamiltonian"][index]),
            "momentum_constraint": float(arrays["momentum"][index]),
        }

    result = {
        "target_R": target,
        "left_index": int(left),
        "right_index": int(right),
        "theta": theta,
        "left_cell": cell(left),
        "right_cell": cell(right),
        "mass_gradient_discrete_interpolated_at_target": interp(arrays["mass_gradient"]),
        "local_composite_rhs_interpolated_at_target": interp(local_rhs),
        "local_composite_residual_interpolated_at_target": interp(residual),
        "normal_velocity_Rdot_over_alpha_interpolated_at_target": interp(normal_velocity),
        "hamiltonian_constraint_contribution_interpolated_at_target": interp(
            hamiltonian_contribution
        ),
        "momentum_constraint_contribution_interpolated_at_target": interp(
            momentum_contribution
        ),
        "constraint_predicted_mass_gradient_residual_interpolated_at_target": interp(
            constraint_predicted_residual
        ),
        "constraint_identity_operator_remainder_interpolated_at_target": interp(
            constraint_identity_remainder
        ),
        "factor_interpolated_rhs_at_target": float(factor_rhs),
        "factor_interpolated_residual_at_target": float(
            interp(arrays["mass_gradient"]) - factor_rhs
        ),
        "factor_product_minus_composite_rhs": float(factor_rhs - interp(local_rhs)),
    }
    if not all(
        math.isfinite(value)
        for value in (
            result["mass_gradient_discrete_interpolated_at_target"],
            result["local_composite_rhs_interpolated_at_target"],
            result["local_composite_residual_interpolated_at_target"],
            result["normal_velocity_Rdot_over_alpha_interpolated_at_target"],
            result["hamiltonian_constraint_contribution_interpolated_at_target"],
            result["momentum_constraint_contribution_interpolated_at_target"],
            result["constraint_predicted_mass_gradient_residual_interpolated_at_target"],
            result["constraint_identity_operator_remainder_interpolated_at_target"],
            result["factor_interpolated_rhs_at_target"],
            result["factor_interpolated_residual_at_target"],
            result["factor_product_minus_composite_rhs"],
        )
    ):
        raise RuntimeError("mass-gradient witness produced a non-finite result")
    return result



def _hm_constraint_adjusted_fixed_R_rate(
    mass_gradient_witness: dict, *, tube_speed: float,
    composite_stress_flux: float, composite_pressure_work: float,
    measured_semidiscrete_rate: float, measured_interpolant_rate: float,
    target: float,
) -> dict:
    """Compare native and H/M-substituted tube sweeps at one exact stage.

    H_term and M_term belong to the spatial identity for D_r M_MS. They enter
    a fixed-R rate through the tube-coordinate speed, not as naked rate terms.
    Therefore P_HM uses F_composite + W_composite
    + (G_composite + H_term + M_term) * v_tube, while P_native uses
    F_composite + W_composite + (D_r M_MS)_native * v_tube.
    The unassigned operator remainder is retained explicitly; this helper
    neither thresholds it nor treats it as a new physical source.
    """
    if not isinstance(mass_gradient_witness, dict):
        raise RuntimeError("H/M fixed-R rate witness must be a mapping")
    required = (
        "target_R", "left_index", "right_index", "left_cell", "right_cell",
        "mass_gradient_discrete_interpolated_at_target",
        "local_composite_rhs_interpolated_at_target",
        "local_composite_residual_interpolated_at_target",
        "hamiltonian_constraint_contribution_interpolated_at_target",
        "momentum_constraint_contribution_interpolated_at_target",
        "constraint_identity_operator_remainder_interpolated_at_target",
    )
    missing = [key for key in required if key not in mass_gradient_witness]
    if missing:
        raise RuntimeError(f"H/M fixed-R rate witness missing fields: {missing}")

    def finite(value, label):
        try:
            number = float(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise RuntimeError(f"H/M fixed-R rate {label} must be finite numeric") from exc
        if not math.isfinite(number):
            raise RuntimeError(f"H/M fixed-R rate {label} must be finite")
        return number

    target_value = finite(target, "target")
    target_witness = finite(mass_gradient_witness["target_R"], "witness target")
    identity_scale = max(1.0, abs(target_value), abs(target_witness))
    if abs(target_value - target_witness) > 64.0 * np.finfo(float).eps * identity_scale:
        raise RuntimeError("H/M fixed-R rate target does not match mass-gradient witness")

    left = mass_gradient_witness["left_index"]
    right = mass_gradient_witness["right_index"]
    if (
        isinstance(left, bool) or isinstance(right, bool)
        or not isinstance(left, (int, np.integer))
        or not isinstance(right, (int, np.integer))
        or int(right) != int(left) + 1
    ):
        raise RuntimeError("H/M fixed-R rate witness has invalid bracket indices")
    for side, expected in (("left_cell", int(left)), ("right_cell", int(right))):
        cell = mass_gradient_witness[side]
        if not isinstance(cell, dict) or cell.get("index") != expected:
            raise RuntimeError(
                f"H/M fixed-R rate {side} identity does not match the bracket"
            )

    values = {
        key: finite(mass_gradient_witness[key], key)
        for key in (
            "mass_gradient_discrete_interpolated_at_target",
            "local_composite_rhs_interpolated_at_target",
            "local_composite_residual_interpolated_at_target",
            "hamiltonian_constraint_contribution_interpolated_at_target",
            "momentum_constraint_contribution_interpolated_at_target",
            "constraint_identity_operator_remainder_interpolated_at_target",
        )
    }
    speed = finite(tube_speed, "tube speed")
    flux = finite(composite_stress_flux, "composite stress flux")
    work = finite(composite_pressure_work, "composite radial pressure work")
    measured_chain = finite(measured_semidiscrete_rate, "semidiscrete measured rate")
    measured_interp = finite(measured_interpolant_rate, "interpolant measured rate")

    mass_gradient_native = values["mass_gradient_discrete_interpolated_at_target"]
    gradient_rhs = values["local_composite_rhs_interpolated_at_target"]
    gradient_residual = values["local_composite_residual_interpolated_at_target"]
    h_term = values["hamiltonian_constraint_contribution_interpolated_at_target"]
    m_term = values["momentum_constraint_contribution_interpolated_at_target"]
    eps_op = values["constraint_identity_operator_remainder_interpolated_at_target"]
    hm_sum = h_term + m_term
    gradient_from_hm = gradient_rhs + hm_sum
    gradient_from_full_identity = gradient_from_hm + eps_op
    identity_reconstruction_residual = mass_gradient_native - gradient_from_full_identity
    identity_scale = max(
        abs(mass_gradient_native), abs(gradient_rhs), abs(h_term), abs(m_term),
        abs(eps_op), float(np.finfo(float).tiny),
    )
    roundoff_bound = 256.0 * np.finfo(float).eps * identity_scale
    if abs(identity_reconstruction_residual) > roundoff_bound:
        raise RuntimeError("H/M fixed-R rate gradient identity reconstruction failed")

    base_flux_work = flux + work
    native_sweep = mass_gradient_native * speed
    h_sweep = h_term * speed
    m_sweep = m_term * speed
    eps_sweep = eps_op * speed
    constraint_sweep = gradient_from_hm * speed
    full_identity_sweep = gradient_from_full_identity * speed
    native_prediction = base_flux_work + native_sweep
    adjusted_prediction = base_flux_work + constraint_sweep
    full_identity_prediction = base_flux_work + full_identity_sweep
    adjusted_delta = adjusted_prediction - native_prediction
    expected_adjusted_delta = -eps_sweep
    rate_scale = max(
        abs(native_prediction), abs(adjusted_prediction), abs(expected_adjusted_delta),
        abs(base_flux_work), float(np.finfo(float).tiny),
    )
    rate_roundoff_bound = 256.0 * np.finfo(float).eps * rate_scale
    if abs(adjusted_delta - expected_adjusted_delta) > rate_roundoff_bound:
        raise RuntimeError("H/M fixed-R rate sweep substitution failed reconstruction")

    return {
        "target_R": target_value,
        "left_index": int(left),
        "right_index": int(right),
        "tube_coordinate_speed_for_fixed_R": speed,
        "base_composite_flux_plus_work_rate": base_flux_work,
        "composite_stress_flux_rate": flux,
        "composite_radial_pressure_work_rate": work,
        "native_mass_gradient_at_R_star": mass_gradient_native,
        "composite_mass_gradient_rhs_at_R_star": gradient_rhs,
        "hamiltonian_mass_gradient_term_at_R_star": h_term,
        "momentum_mass_gradient_term_at_R_star": m_term,
        "operator_remainder_mass_gradient_at_R_star": eps_op,
        "hamiltonian_term_contribution_to_sweep_rate": h_sweep,
        "momentum_term_contribution_to_sweep_rate": m_sweep,
        "operator_remainder_contribution_to_sweep_rate": eps_sweep,
        "constraint_substituted_mass_gradient_at_R_star": gradient_from_hm,
        "full_identity_reconstructed_mass_gradient_at_R_star": gradient_from_full_identity,
        "native_sweep_rate": native_sweep,
        "constraint_substituted_sweep_rate": constraint_sweep,
        "full_identity_reconstructed_sweep_rate": full_identity_sweep,
        "native_composite_fixed_R_prediction": native_prediction,
        "constraint_adjusted_fixed_R_prediction": adjusted_prediction,
        "full_identity_reconstructed_fixed_R_prediction": full_identity_prediction,
        "adjusted_minus_native_prediction": adjusted_delta,
        "expected_adjusted_minus_native_from_epsilon_sweep": expected_adjusted_delta,
        "measured_semidiscrete_fixed_R_rate": measured_chain,
        "measured_moving_interpolant_fixed_R_rate": measured_interp,
        "semidiscrete_measured_minus_native_prediction": measured_chain - native_prediction,
        "semidiscrete_measured_minus_adjusted_prediction": measured_chain - adjusted_prediction,
        "interpolant_measured_minus_native_prediction": measured_interp - native_prediction,
        "interpolant_measured_minus_adjusted_prediction": measured_interp - adjusted_prediction,
        "gradient_identity_reconstruction_residual": identity_reconstruction_residual,
        "gradient_identity_roundoff_bound": roundoff_bound,
        "rate_sweep_reconstruction_residual": adjusted_delta - expected_adjusted_delta,
        "rate_sweep_roundoff_bound": rate_roundoff_bound,
    }

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


def _fixed_R_interpolant_balance(interpolant_rate: dict, stress_flux: float, pressure_work: float):
    """Balance the fixed-R derivative against the *same* interpolant's sweep.

    The native-derivative tube sweep and the piecewise-linear M(R,t)
    interpolant sweep are distinct discrete operators. A direct-interpolant
    rate must be compared with flux + pressure work + its own sweep.
    """
    fixed_rate = float(interpolant_rate["mass_dot_fixed_R"])
    sweep = float(interpolant_rate["fixed_R_sweep"])
    stress_flux = float(stress_flux)
    pressure_work = float(pressure_work)
    if not all(math.isfinite(value) for value in (fixed_rate, sweep, stress_flux, pressure_work)):
        raise RuntimeError("fixed-R interpolant balance contains non-finite inputs")
    predicted = stress_flux + pressure_work + sweep
    return predicted, fixed_rate - predicted


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
    explicit_metric_rates = _explicit_metric_rate_profiles(
        grid, R, geom.a, geom.b, geom.X, explicit
    )
    Rdot_reconstruction = _profile_roundoff_reconstruction(
        Rdot, explicit_metric_rates["Rdot"], label="captured explicit Rdot"
    )

    total = _stress_projection(kernel, state)
    _, vacuum, _ = adapter.vendor_modules()
    raw_constraints = vacuum.constraints(grid, geom)
    hamiltonian_constraint_profile = (
        np.asarray(raw_constraints["hamiltonian"], dtype=float)
        - 16.0 * math.pi * np.asarray(total.rho, dtype=float)
    )
    momentum_constraint_profile = (
        np.asarray(raw_constraints["momentum"], dtype=float)
        - 8.0 * math.pi * np.asarray(total.j, dtype=float)
    )
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

    # Component-level attribution (scalar / fluid / electromagnetic).
    # Fail-closed additive consistency against the aggregate rates above.
    sector_rates = _sector_fixed_R_rates(
        total,
        R=R,
        alpha=geom.alpha,
        gamma_rr_inv=gamma_rr_inv,
        Rr=Rr,
        Rdot=Rdot,
        target=FIXED_AREAL_RADIUS,
    )

    # Independent named discriminators. The frozen-geometry expression uses
    # these same supplied snapshot profiles and does not re-evolve geometry.
    frozen_geometry_balance = frozen_geometry_stress_balance(
        R=R, alpha=geom.alpha, gamma_rr_inv=gamma_rr_inv, Rr=Rr,
        Rdot=Rdot, j=total.j, pr=total.pr, rho=total.rho,
        target_areal_radius=FIXED_AREAL_RADIUS,
    )
    sector_discriminator = pure_sector_rates(
        total, R=R, alpha=geom.alpha, gamma_rr_inv=gamma_rr_inv,
        Rr=Rr, Rdot=Rdot, target_areal_radius=FIXED_AREAL_RADIUS,
    )
    # Keep the legacy sector ledger and the extracted pure helper aligned.
    # This is an implementation-consistency check, not a physical residual test.
    for sector_name, legacy in sector_rates["sectors"].items():
        extracted = sector_discriminator["sectors"][sector_name]
        for term in ("stress_flux", "radial_pressure_work"):
            old_value = float(legacy[term])
            new_value = float(extracted[term])
            scale = max(abs(old_value) + abs(new_value),
                        float(np.finfo(np.float64).tiny))
            if (not math.isfinite(old_value) or not math.isfinite(new_value)
                    or abs(old_value - new_value)
                    > 64.0 * np.finfo(np.float64).eps * scale):
                raise RuntimeError(
                    "extracted pure-sector discriminator disagrees with legacy "
                    f"sector ledger: {sector_name}.{term}"
                )
    composite_flux_profile, composite_pressure_work_profile = (
        _fixed_R_composite_flux_work_profiles(
            R, geom.alpha, gamma_rr_inv, total.j, Rr, total.pr, Rdot
        )
    )
    composite_stress_flux = at(composite_flux_profile)
    composite_pressure_work = at(composite_pressure_work_profile)
    mass_gradient_rhs = FOUR_PI * FIXED_AREAL_RADIUS**2 * (
        rho_tube * Rr_tube + j_tube * Rdot_tube / alpha_tube
    )
    mass_gradient_residual = Mr_tube - mass_gradient_rhs
    mass_gradient_witness = _mass_gradient_identity_bracket_witness(
        R, Mr, total.rho, total.j, Rr, Rdot, geom.alpha,
        hamiltonian_constraint_profile, momentum_constraint_profile,
        target=FIXED_AREAL_RADIUS,
    )

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
        gamma_inv_dot_reconstruction = _profile_roundoff_reconstruction(
            gamma_inv_dot, explicit_metric_rates["gamma_inv_dot"],
            label="captured explicit gamma_inv_dot"
        )
        Ktheta = geom.K / 3.0 - geom.Aa / 2.0
        Ktheta_dot = K_dot / 3.0 - Aa_dot / 2.0
        vacuum_primary_l3 = split["vacuum_primary_l3"]
        ktheta_source_rates = _sector_ktheta_source_rates(
            total, R=R, alpha=geom.alpha, Ktheta=Ktheta, Rdot=Rdot,
            target=FIXED_AREAL_RADIUS,
            live_primary_l3=split["primary_l3"],
            vacuum_primary_l3=vacuum_primary_l3,
        )
        Rdot_r = np.asarray(
            grid.cell_derivative_fourth(Rdot, parity=-1), dtype=float
        )
        Rdot_r_reconstruction = _profile_roundoff_reconstruction(
            Rdot_r, explicit_metric_rates["Rdot_r"],
            label="captured explicit D_r(Rdot)"
        )
        metric_rate_reconstruction_profiles = {
            "Rdot_residual": Rdot_reconstruction["residual"],
            "Rdot_roundoff_bound": Rdot_reconstruction["roundoff_bound"],
            "gamma_inv_dot_residual": gamma_inv_dot_reconstruction["residual"],
            "gamma_inv_dot_roundoff_bound": gamma_inv_dot_reconstruction["roundoff_bound"],
            "Rdot_r_residual": Rdot_r_reconstruction["residual"],
            "Rdot_r_roundoff_bound": Rdot_r_reconstruction["roundoff_bound"],
        }
        chain_profiles = _misner_sharp_chain_rule_rate_profiles(
            R, Rdot, Rr, Rdot_r, gamma_rr_inv, gamma_inv_dot,
            Ktheta, Ktheta_dot,
        )
        chain_Mr = np.asarray(
            grid.cell_derivative_fourth(chain_profiles["mass"], parity=-1),
            dtype=float,
        )
        primary_l2_Ktheta_dot_profile = (
            np.asarray(split["primary_l2"]["K"], dtype=float) / 3.0
            - np.asarray(split["primary_l2"]["Aa"], dtype=float) / 2.0
        )
        vacuum_l3_Ktheta_dot_profile = (
            np.asarray(vacuum_primary_l3["K"], dtype=float) / 3.0
            - np.asarray(vacuum_primary_l3["Aa"], dtype=float) / 2.0
        )
        source_l3_K = (
            np.asarray(split["primary_l3"]["K"], dtype=float)
            - np.asarray(vacuum_primary_l3["K"], dtype=float)
        )
        source_l3_Aa = (
            np.asarray(split["primary_l3"]["Aa"], dtype=float)
            - np.asarray(vacuum_primary_l3["Aa"], dtype=float)
        )
        source_Ktheta_dot_profile = source_l3_K / 3.0 - 0.5 * source_l3_Aa
        ktheta_dot_component_profiles = _ktheta_dot_component_profiles(
            split["primary_l2"]["K"], split["primary_l2"]["Aa"],
            vacuum_primary_l3["K"], vacuum_primary_l3["Aa"],
            source_l3_K, source_l3_Aa,
            expected_primary_l2=primary_l2_Ktheta_dot_profile,
            expected_vacuum_l3=vacuum_l3_Ktheta_dot_profile,
            expected_einstein_source_l3=source_Ktheta_dot_profile,
        )
        chain_bracket_witness = _misner_sharp_chain_rule_bracket_operands(
            R, Rdot, Rr, Rdot_r, gamma_rr_inv, gamma_inv_dot,
            Ktheta, Ktheta_dot, primary_l2_Ktheta_dot_profile,
            vacuum_l3_Ktheta_dot_profile, source_Ktheta_dot_profile,
            target=FIXED_AREAL_RADIUS,
            ktheta_dot_component_profiles=ktheta_dot_component_profiles,
            metric_rate_reconstruction_profiles=metric_rate_reconstruction_profiles,
        )
        # Attach the actual-stage operands and local composite balance to both
        # cells, keeping their values adjacent to D_r M, H/M, and the RHS operators.
        for side in ("left_cell", "right_cell"):
            cell = mass_gradient_witness[side]
            index = int(cell["index"])
            cell["stage_chain_rule_operands"] = chain_bracket_witness[side]
            local_flux = (
                -FOUR_PI * R[index]**2 * geom.alpha[index]
                * gamma_rr_inv[index] * total.j[index] * Rr[index]
            )
            local_pressure_work = (
                -FOUR_PI * R[index]**2 * total.pr[index] * Rdot[index]
            )
            if Rr[index] == 0.0:
                raise RuntimeError("fixed-R bracket cell has zero D_r R")
            native_sweep = -Rdot[index] * Mr[index] / Rr[index]
            chain_sweep_cell = -Rdot[index] * chain_Mr[index] / Rr[index]
            coordinate_rate = chain_profiles["mass_dot_coordinate"][index]
            cell["stage_local_composite_balance"] = {
                "stress_flux": float(local_flux),
                "radial_pressure_work": float(local_pressure_work),
                "coordinate_rate_minus_flux_and_work": float(
                    coordinate_rate - local_flux - local_pressure_work
                ),
                "native_mass_gradient_tube_sweep": float(native_sweep),
                "chain_mass_gradient_tube_sweep": float(chain_sweep_cell),
                "native_minus_chain_mass_gradient": float(
                    Mr[index] - chain_Mr[index]
                ),
                "fixed_R_residual_using_native_mass_gradient": float(
                    coordinate_rate + native_sweep - local_flux - local_pressure_work
                ),
                "fixed_R_residual_using_chain_mass_gradient": float(
                    coordinate_rate + chain_sweep_cell - local_flux - local_pressure_work
                ),
            }
        # Attach the same-stage mass-rate tangent terms to the exact cells used
        # by the native D_r M identity witness. This keeps the local constraint
        # operands, mass-gradient residual, and mass-rate decomposition adjacent.
        chain_term_names = (
            "mass_dot_coordinate_term_Rdot",
            "mass_dot_coordinate_term_gamma_inv_dot",
            "mass_dot_coordinate_term_Rdot_r",
            "mass_dot_coordinate_term_Ktheta_geometry",
            "mass_dot_coordinate_term_Ktheta_dot",
        )
        for side in ("left_cell", "right_cell"):
            cell = mass_gradient_witness[side]
            cell_index = int(cell["index"])
            cell["mass_dot_coordinate_chain_rule"] = float(
                chain_profiles["mass_dot_coordinate"][cell_index]
            )
            cell["mass_dot_coordinate_terms"] = {
                name.removeprefix("mass_dot_coordinate_term_"): float(
                    chain_profiles[name][cell_index]
                )
                for name in chain_term_names
            }
            cell["mass_dot_coordinate_term_sum"] = float(
                chain_profiles["mass_dot_coordinate_term_sum"][cell_index]
            )
            cell["mass_dot_coordinate_decomposition_residual"] = float(
                chain_profiles["mass_dot_coordinate_decomposition_residual"][cell_index]
            )
        chain_Mdot_coord = at(chain_profiles["mass_dot_coordinate"])
        chain_Mr_tube = at(chain_Mr)
        chain_sweep = chain_Mr_tube * tube_speed
        chain_Mdot_fixed_R = chain_Mdot_coord + chain_sweep
        chain_coordinate_residual = chain_Mdot_coord - (stress_flux + pressure_work)
        chain_fixed_R_residual = chain_Mdot_fixed_R - predicted_fixed_R_rate
        ktheta_matter_source_mass_dot = float(
            ktheta_source_rates["aggregate_at_R_star"][
                "mass_dot_coordinate_from_ktheta_source"
            ]
        )

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
        interpolant_predicted_fixed_R_rate, interpolant_fixed_R_residual = (
            _fixed_R_interpolant_balance(interpolant_rate, stress_flux, pressure_work)
        )
        interpolant_composite_coordinate_residual = (
            interpolant_rate["interpolated_mass_dot_coordinate"]
            - (composite_stress_flux + composite_pressure_work)
        )
        interpolant_composite_predicted_fixed_R_rate, interpolant_composite_fixed_R_residual = (
            _fixed_R_interpolant_balance(
                interpolant_rate, composite_stress_flux, composite_pressure_work
            )
        )
        hm_constraint_adjusted_rate = _hm_constraint_adjusted_fixed_R_rate(
            mass_gradient_witness,
            tube_speed=tube_speed,
            composite_stress_flux=composite_stress_flux,
            composite_pressure_work=composite_pressure_work,
            measured_semidiscrete_rate=chain_Mdot_fixed_R,
            measured_interpolant_rate=interpolant_rate["mass_dot_fixed_R"],
            target=FIXED_AREAL_RADIUS,
        )
        chain_rule = {
            "tube_chain_rule_Rdot": at(Rdot),
            "tube_chain_rule_Rr": at(Rr),
            "tube_chain_rule_Rdot_r": at(Rdot_r),
            "tube_chain_rule_gamma_inv": at(gamma_rr_inv),
            "tube_chain_rule_gamma_inv_dot": at(gamma_inv_dot),
            "tube_chain_rule_Ktheta": at(Ktheta),
            "tube_chain_rule_Ktheta_dot": at(Ktheta_dot),
            "tube_ktheta_matter_source_rates": ktheta_source_rates,
            "tube_chain_rule_Aa_dot": at(Aa_dot),
            "tube_chain_rule_K_dot": at(K_dot),
            "tube_chain_rule_chi_dot": at(chain_profiles["chi_dot"]),
            "tube_chain_rule_chi_dot_absmax": float(np.max(np.abs(chain_profiles["chi_dot"]))),
            "tube_chain_rule_mass_dot_coordinate": chain_Mdot_coord,
            "tube_chain_rule_mass_dot_term_Rdot": at(
                chain_profiles["mass_dot_coordinate_term_Rdot"]
            ),
            "tube_chain_rule_mass_dot_term_gamma_inv_dot": at(
                chain_profiles["mass_dot_coordinate_term_gamma_inv_dot"]
            ),
            "tube_chain_rule_mass_dot_term_Rdot_r": at(
                chain_profiles["mass_dot_coordinate_term_Rdot_r"]
            ),
            "tube_chain_rule_mass_dot_term_Ktheta_geometry": at(
                chain_profiles["mass_dot_coordinate_term_Ktheta_geometry"]
            ),
            "tube_chain_rule_mass_dot_term_Ktheta_dot": at(
                chain_profiles["mass_dot_coordinate_term_Ktheta_dot"]
            ),
            "tube_chain_rule_mass_dot_coordinate_term_sum": at(
                chain_profiles["mass_dot_coordinate_term_sum"]
            ),
            "tube_chain_rule_mass_dot_coordinate_decomposition_residual": at(
                chain_profiles["mass_dot_coordinate_decomposition_residual"]
            ),
            "tube_chain_rule_mass_dot_term_Ktheta_matter_source": (
                ktheta_matter_source_mass_dot
            ),
            "tube_chain_rule_mass_dot_term_Ktheta_dot_minus_explicit_matter_source": (
                at(chain_profiles["mass_dot_coordinate_term_Ktheta_dot"])
                - ktheta_matter_source_mass_dot
            ),
            "tube_chain_rule_mass_gradient": chain_Mr_tube,
            "tube_chain_rule_fixed_R_sweep": chain_sweep,
            "tube_chain_rule_mass_dot_fixed_R": chain_Mdot_fixed_R,
            "tube_chain_rule_coordinate_rate_minus_flux_pressure": chain_coordinate_residual,
            "tube_chain_rule_fixed_R_residual": chain_fixed_R_residual,
            "tube_chain_rule_flux_work_operator_comparison": {
                "stress_flux_factor_interpolated": float(stress_flux),
                "stress_flux_local_composite_interpolated": float(composite_stress_flux),
                "radial_pressure_work_factor_interpolated": float(pressure_work),
                "radial_pressure_work_local_composite_interpolated": float(composite_pressure_work),
                "fixed_R_residual_factor_interpolated": float(chain_fixed_R_residual),
                "fixed_R_residual_local_composite_interpolated": float(
                    chain_Mdot_fixed_R
                    - composite_stress_flux - composite_pressure_work - chain_sweep
                ),
            },
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
            "tube_interpolant_predicted_fixed_R_rate": interpolant_predicted_fixed_R_rate,
            "tube_interpolant_composite_coordinate_rate_minus_flux_pressure": interpolant_composite_coordinate_residual,
            "tube_interpolant_composite_predicted_fixed_R_rate": interpolant_composite_predicted_fixed_R_rate,
            "tube_interpolant_composite_fixed_R_residual": interpolant_composite_fixed_R_residual,
            "tube_interpolant_mass_dot_fixed_R": interpolant_rate["mass_dot_fixed_R"],
            "tube_interpolant_coordinate_rate_minus_flux_pressure": interpolant_coordinate_residual,
            "tube_interpolant_fixed_R_residual": interpolant_fixed_R_residual,
            "tube_interpolant_witness_profiles": {
                "R": np.asarray(R, dtype=float).copy(),
                "Rdot": np.asarray(Rdot, dtype=float).copy(),
                "mass": np.asarray(chain_profiles["mass"], dtype=float).copy(),
                "mass_dot_coordinate": np.asarray(chain_profiles["mass_dot_coordinate"], dtype=float).copy(),
            },
            "tube_hm_constraint_adjusted_fixed_R_rate": hm_constraint_adjusted_rate,
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
        "tube_hamiltonian_constraint_at_R_star": at(hamiltonian_constraint_profile),
        "tube_momentum_constraint_at_R_star": at(momentum_constraint_profile),
        "tube_Rr": Rr_tube,
        "tube_Rdot_at_fixed_coordinate": Rdot_tube,
        "tube_coordinate_speed_for_fixed_R": tube_speed,
        "tube_mass_gradient": Mr_tube,
        "tube_stress_flux": stress_flux,
        "tube_radial_pressure_work": pressure_work,
        "tube_stress_flux_composite_interp": composite_stress_flux,
        "tube_radial_pressure_work_composite_interp": composite_pressure_work,
        "tube_coordinate_sweep": tube_sweep,
        "tube_predicted_fixed_R_mass_rate": predicted_fixed_R_rate,
        "tube_sector_rates": sector_rates,
        "frozen_geometry_balance": frozen_geometry_balance,
        "pure_sector_rates": sector_discriminator,
        "tube_native_mass_gradient_identity_rhs": mass_gradient_rhs,
        "tube_native_mass_gradient_identity_residual": mass_gradient_residual,
        "tube_mass_gradient_identity_witness": mass_gradient_witness,
        **(chain_rule or {}),
        "constraint_witnesses": constraint_witnesses,
    }


def _fixed_R_rate_snapshot(sample: dict) -> dict:
    """Keep signed fixed-R rate components for one actual PIRK RHS stage."""
    keys = (
        "tube_coordinate_r",
        "tube_mass_misner_sharp",
        "tube_hamiltonian_constraint_at_R_star",
        "tube_momentum_constraint_at_R_star",
        "tube_stress_flux",
        "tube_radial_pressure_work",
        "tube_stress_flux_composite_interp",
        "tube_radial_pressure_work_composite_interp",
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
        "tube_chain_rule_mass_dot_term_Rdot",
        "tube_chain_rule_mass_dot_term_gamma_inv_dot",
        "tube_chain_rule_mass_dot_term_Rdot_r",
        "tube_chain_rule_mass_dot_term_Ktheta_geometry",
        "tube_chain_rule_mass_dot_term_Ktheta_dot",
        "tube_chain_rule_mass_dot_coordinate_term_sum",
        "tube_chain_rule_mass_dot_coordinate_decomposition_residual",
        "tube_chain_rule_mass_dot_term_Ktheta_matter_source",
        "tube_chain_rule_mass_dot_term_Ktheta_dot_minus_explicit_matter_source",
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
        "tube_interpolant_predicted_fixed_R_rate",
        "tube_interpolant_composite_coordinate_rate_minus_flux_pressure",
        "tube_interpolant_composite_predicted_fixed_R_rate",
        "tube_interpolant_composite_fixed_R_residual",
        "tube_interpolant_mass_dot_fixed_R",
        "tube_interpolant_coordinate_rate_minus_flux_pressure",
        "tube_interpolant_fixed_R_residual",
    )
    snapshot = {"t": float(sample["t"]), **{key: float(sample[key]) for key in keys}}
    if "tube_interpolant_witness_profiles" in sample:
        snapshot["tube_interpolant_witness_profiles"] = sample["tube_interpolant_witness_profiles"]
    for key in (
        "tube_mass_gradient_identity_witness",
        "tube_chain_rule_flux_work_operator_comparison",
        "tube_hm_constraint_adjusted_fixed_R_rate",
    ):
        if key in sample:
            snapshot[key] = sample[key]
    if "tube_sector_rates" in sample:
        snapshot["tube_sector_rates"] = sample["tube_sector_rates"]
    for key in ("frozen_geometry_balance", "pure_sector_rates"):
        if key in sample:
            snapshot[key] = sample[key]
    if "tube_ktheta_matter_source_rates" in sample:
        snapshot["tube_ktheta_matter_source_rates"] = sample["tube_ktheta_matter_source_rates"]
    if "tube_mass_gradient_identity_witness" in sample:
        snapshot["tube_mass_gradient_identity_witness"] = sample["tube_mass_gradient_identity_witness"]
    return snapshot


def _pirk_accepted_update_rate_arrays(geometry_stages, final_primary_l2) -> dict:
    """Reconstruct the exact accepted-step coefficient weights for geometry variables.

    V55ProductionKernel.step advances explicit (a,b,X) with start/predictor
    trapezoid weights. Its accepted (Aa,K) update uses L2(start), L2(final
    mixed explicit / predictor implicit state), L3(start), and L3(primary
    predictor), each with weight 1/2. The predictor-only L2 evaluation is
    intentionally not substituted for the final mixed-state L2 call.
    """
    if len(geometry_stages) != 3:
        raise RuntimeError("accepted PIRK update needs the three captured geometry evaluations")
    start = geometry_stages[0]
    predictor = geometry_stages[2]

    def profile(value, *, name):
        array = np.asarray(value, dtype=float)
        if array.ndim != 1 or not np.all(np.isfinite(array)):
            raise RuntimeError(f"accepted PIRK update {name} must be a finite 1-D profile")
        return array

    explicit = {}
    for name in ("a", "b", "X"):
        try:
            left = profile(start["explicit"][name], name=f"start {name} RHS")
            right = profile(predictor["explicit"][name], name=f"predictor {name} RHS")
        except KeyError as exc:
            raise RuntimeError(f"accepted PIRK update is missing explicit profile: {exc}") from exc
        if left.shape != right.shape:
            raise RuntimeError(f"accepted PIRK update explicit {name} profiles have different shapes")
        explicit[name] = 0.5 * (left + right)

    primary = {}
    for name in ("Aa", "K"):
        try:
            terms = (
                profile(start["primary_l2"][name], name=f"start L2 {name}"),
                profile(final_primary_l2[name], name=f"final mixed-state L2 {name}"),
                profile(start["primary_l3"][name], name=f"start L3 {name}"),
                profile(predictor["primary_l3"][name], name=f"predictor L3 {name}"),
            )
        except KeyError as exc:
            raise RuntimeError(f"accepted PIRK update is missing primary profile: {exc}") from exc
        if any(term.shape != terms[0].shape for term in terms[1:]):
            raise RuntimeError(f"accepted PIRK update primary {name} profiles have different shapes")
        primary[name] = 0.5 * sum(terms)
    return {"explicit": explicit, "primary": primary}


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
            _, vacuum, _ = adapter.vendor_modules()
            vacuum_l3 = vacuum.primary_l3_rhs(grid, geometry)
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
                "vacuum_primary_l3": {
                    name: _array(vacuum_l3[name]) for name in ("Aa", "K")
                },
                "einstein_source_l3": {
                    name: _array(np.asarray(result["primary_l3"][name], dtype=float)
                                 - np.asarray(vacuum_l3[name], dtype=float))
                    for name in ("Aa", "K")
                },
                "einstein_source_roundoff_scale": {
                    name: _array(np.abs(np.asarray(result["primary_l3"][name], dtype=float))
                                 + np.abs(np.asarray(vacuum_l3[name], dtype=float)))
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
    _, vacuum, _ = adapter.vendor_modules()
    original_primary_l2 = vacuum.primary_l2_rhs

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

    def primary_l2_wrapper(grid, geometry, *args, **kwargs):
        result = original_primary_l2(grid, geometry, *args, **kwargs)
        active_step = kernel._active_step
        # The accepted step's fourth primary-L2 call occurs after all three
        # _geometry_stage_terms observations: its input mixes the final
        # explicit (a,b,X) update with the primary-predictor implicit fields.
        if (
            active_step is not None
            and len(active_step["geometry_stages"]) == 3
            and "final_primary_l2" not in active_step
        ):
            active_step["final_primary_l2"] = {
                "geometry": {
                    name: _array(getattr(geometry, name))
                    for name in ("a", "b", "X", "Aa", "K", "Lambda")
                },
                "Aa": _array(result["Aa"]),
                "K": _array(result["K"]),
            }
        return result

    matter_rhs_module._species_rhs_components = components_wrapper
    RadialMaxwellState.lorentz_source = lorentz_wrapper
    vacuum.primary_l2_rhs = primary_l2_wrapper
    try:
        yield
    finally:
        matter_rhs_module._species_rhs_components = original_components
        RadialMaxwellState.lorentz_source = original_lorentz
        vacuum.primary_l2_rhs = original_primary_l2


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


def _accepted_pirk_tangent_term_profiles(
    *,
    R, Rr, gamma_inv, Ktheta,
    Rdot, gamma_inv_dot, Ktheta_dot,
    Rdot_contributions: dict[str, np.ndarray],
    gamma_inv_dot_contributions: dict[str, np.ndarray],
    Ktheta_dot_contributions: dict[str, np.ndarray],
    grid,
) -> dict[str, np.ndarray]:
    """Expand accepted weighted PIRK tangent into signed RHS-piece contributions."""
    raw = {
        "R": R, "Rr": Rr, "gamma_inv": gamma_inv, "Ktheta": Ktheta,
        "Rdot": Rdot, "gamma_inv_dot": gamma_inv_dot, "Ktheta_dot": Ktheta_dot,
    }
    arrays = {name: np.asarray(value, dtype=float) for name, value in raw.items()}
    radius = arrays["R"]
    if radius.ndim != 1 or radius.size < 2:
        raise RuntimeError("accepted tangent term split R must be 1-D with at least two cells")
    for name, values in arrays.items():
        if values.shape != radius.shape or not np.all(np.isfinite(values)):
            raise RuntimeError(f"accepted tangent term split {name} has invalid shape or values")
    expected = {
        "Rdot": {"explicit_start", "explicit_predictor"},
        "gamma_inv_dot": {"explicit_start", "explicit_predictor"},
        "Ktheta_dot": {
            "primary_l2_start", "primary_l2_final_mixed_state",
            "vacuum_l3_start", "vacuum_l3_predictor",
            "einstein_source_l3_start", "einstein_source_l3_predictor",
        },
    }
    groups = {
        "Rdot": Rdot_contributions,
        "gamma_inv_dot": gamma_inv_dot_contributions,
        "Ktheta_dot": Ktheta_dot_contributions,
    }
    checked = {}
    for group, values in groups.items():
        if not isinstance(values, dict) or set(values) != expected[group]:
            raise RuntimeError(
                f"accepted tangent term split {group} requires exactly {sorted(expected[group])}"
            )
        checked[group] = {}
        for name, raw_profile in values.items():
            profile = np.asarray(raw_profile, dtype=float)
            if profile.shape != radius.shape or not np.all(np.isfinite(profile)):
                raise RuntimeError(
                    f"accepted tangent term split {group}.{name} has invalid shape or values"
                )
            checked[group][name] = profile
        summed = sum(checked[group].values())
        residual = arrays[group] - summed
        scale = np.abs(arrays[group]) + sum(np.abs(value) for value in checked[group].values())
        bound = 128.0 * np.finfo(np.float64).eps * np.maximum(scale, np.finfo(np.float64).tiny)
        if np.any(np.abs(residual) > bound):
            raise RuntimeError(f"accepted tangent {group} components fail roundoff reconstruction")
    chi = arrays["gamma_inv"] * arrays["Rr"]**2 - (radius * arrays["Ktheta"])**2
    terms: dict[str, np.ndarray] = {}
    for name, profile in checked["Rdot"].items():
        terms[f"mass_dot_coordinate_term_Rdot_{name}"] = 0.5 * profile * (1.0 - chi)
        terms[f"mass_dot_coordinate_term_Ktheta_geometry_{name}"] = (
            radius**2 * profile * arrays["Ktheta"]**2
        )
        derivative = np.asarray(grid.cell_derivative_fourth(profile, parity=-1), dtype=float)
        if derivative.shape != radius.shape or not np.all(np.isfinite(derivative)):
            raise RuntimeError(f"accepted tangent D_r(Rdot) component {name} is invalid")
        terms[f"mass_dot_coordinate_term_Rdot_r_{name}"] = (
            -radius * arrays["gamma_inv"] * arrays["Rr"] * derivative
        )
    for name, profile in checked["gamma_inv_dot"].items():
        terms[f"mass_dot_coordinate_term_gamma_inv_dot_{name}"] = (
            -0.5 * radius * profile * arrays["Rr"]**2
        )
    for name, profile in checked["Ktheta_dot"].items():
        terms[f"mass_dot_coordinate_term_Ktheta_dot_{name}"] = (
            radius**3 * arrays["Ktheta"] * profile
        )
    if not terms or not all(np.all(np.isfinite(value)) for value in terms.values()):
        raise RuntimeError("accepted tangent term split produced non-finite or empty terms")
    return terms


def _accepted_pirk_update_tangent(step: dict, stage_rates: dict) -> dict:
    """Reconcile the exact accepted-step RHS weights with the fixed-R mass observable."""
    geometry_stages = step["geometry_stages"]
    final_l2 = step.get("final_primary_l2")
    if not isinstance(final_l2, dict) or not {"Aa", "K", "geometry"} <= set(final_l2):
        raise RuntimeError("PIRK observer missed final mixed-state primary-L2 outputs")
    update_rates = _pirk_accepted_update_rate_arrays(geometry_stages, final_l2)
    update_einstein_source = _pirk_accepted_update_einstein_source_arrays(geometry_stages)
    if not step.get("rhs_stages"):
        raise RuntimeError("PIRK update tangent is missing the accepted left-endpoint state")

    state0 = step["rhs_stages"][0]["stage_state"]
    state1 = step["state_after"]
    grid = state0.grid
    geom = state0.geometry
    R = np.asarray(areal_radius(grid, geom), dtype=float)
    Rr = np.asarray(areal_radius_derivative(grid, geom), dtype=float)
    mass = np.asarray(misner_sharp(grid, geom), dtype=float)

    a_dot = update_rates["explicit"]["a"]
    b_dot = update_rates["explicit"]["b"]
    X_dot = update_rates["explicit"]["X"]
    Aa_dot = update_rates["primary"]["Aa"]
    K_dot = update_rates["primary"]["K"]
    Rdot = R * (0.5 * b_dot / geom.b - X_dot / geom.X)
    Rdot_r = np.asarray(grid.cell_derivative_fourth(Rdot, parity=-1), dtype=float)
    gamma_inv = geom.X**2 / geom.a
    gamma_inv_dot = gamma_inv * (2.0 * X_dot / geom.X - a_dot / geom.a)
    Ktheta = geom.K / 3.0 - geom.Aa / 2.0
    Ktheta_dot = K_dot / 3.0 - Aa_dot / 2.0
    source = update_einstein_source["source_l3"]
    source_ktheta_dot = source["K"] / 3.0 - source["Aa"] / 2.0
    source_mass_rate_profile = R**3 * Ktheta * source_ktheta_dot
    source_mass_rate_accepted_weights = interpolate_at_areal_radius(
        R, source_mass_rate_profile, FIXED_AREAL_RADIUS
    )
    chain_profiles = _misner_sharp_chain_rule_rate_profiles(
        R, Rdot, Rr, Rdot_r, gamma_inv, gamma_inv_dot, Ktheta, Ktheta_dot
    )
    direct = _fixed_target_interpolant_rate(
        R, Rdot, chain_profiles["mass"],
        chain_profiles["mass_dot_coordinate"], target=FIXED_AREAL_RADIUS,
    )

    start_geometry = geometry_stages[0]
    predictor_geometry = geometry_stages[2]
    explicit_contributions = {}
    for stage_name, stage in (
        ("explicit_start", start_geometry),
        ("explicit_predictor", predictor_geometry),
    ):
        explicit_contributions[stage_name] = {
            name: 0.5 * np.asarray(stage["explicit"][name], dtype=float)
            for name in ("a", "b", "X")
        }
    Rdot_contributions = {
        name: R * (0.5 * piece["b"] / geom.b - piece["X"] / geom.X)
        for name, piece in explicit_contributions.items()
    }
    gamma_inv_contributions = {
        name: gamma_inv * (2.0 * piece["X"] / geom.X - piece["a"] / geom.a)
        for name, piece in explicit_contributions.items()
    }

    def half_ktheta_dot(k_rhs, aa_rhs):
        return 0.5 * (
            np.asarray(k_rhs, dtype=float) / 3.0
            - np.asarray(aa_rhs, dtype=float) / 2.0
        )

    ktheta_dot_contributions = {
        "primary_l2_start": half_ktheta_dot(
            start_geometry["primary_l2"]["K"], start_geometry["primary_l2"]["Aa"]
        ),
        "primary_l2_final_mixed_state": half_ktheta_dot(
            final_l2["K"], final_l2["Aa"]
        ),
        "vacuum_l3_start": half_ktheta_dot(
            start_geometry["vacuum_primary_l3"]["K"],
            start_geometry["vacuum_primary_l3"]["Aa"],
        ),
        "vacuum_l3_predictor": half_ktheta_dot(
            predictor_geometry["vacuum_primary_l3"]["K"],
            predictor_geometry["vacuum_primary_l3"]["Aa"],
        ),
        "einstein_source_l3_start": half_ktheta_dot(
            start_geometry["einstein_source_l3"]["K"],
            start_geometry["einstein_source_l3"]["Aa"],
        ),
        "einstein_source_l3_predictor": half_ktheta_dot(
            predictor_geometry["einstein_source_l3"]["K"],
            predictor_geometry["einstein_source_l3"]["Aa"],
        ),
    }
    accepted_tangent_term_profiles = _accepted_pirk_tangent_term_profiles(
        R=R, Rr=Rr, gamma_inv=gamma_inv, Ktheta=Ktheta,
        Rdot=Rdot, gamma_inv_dot=gamma_inv_dot, Ktheta_dot=Ktheta_dot,
        Rdot_contributions=Rdot_contributions,
        gamma_inv_dot_contributions=gamma_inv_contributions,
        Ktheta_dot_contributions=ktheta_dot_contributions, grid=grid,
    )
    accepted_tangent_term_sum_profile = sum(accepted_tangent_term_profiles.values())
    accepted_tangent_term_reconstruction = (
        chain_profiles["mass_dot_coordinate"] - accepted_tangent_term_sum_profile
    )
    accepted_tangent_term_scale = (
        np.abs(chain_profiles["mass_dot_coordinate"])
        + sum(np.abs(value) for value in accepted_tangent_term_profiles.values())
    )
    accepted_tangent_term_bound = 128.0 * np.finfo(np.float64).eps * np.maximum(
        accepted_tangent_term_scale, np.finfo(np.float64).tiny
    )
    reconstruction_abs = np.abs(accepted_tangent_term_reconstruction)
    reconstruction_ratio = np.zeros_like(reconstruction_abs)
    nonzero_bound = accepted_tangent_term_bound > 0.0
    np.divide(
        reconstruction_abs, accepted_tangent_term_bound,
        out=reconstruction_ratio, where=nonzero_bound,
    )
    reconstruction_ratio[(~nonzero_bound) & (reconstruction_abs > 0.0)] = np.inf
    if np.any(reconstruction_abs > accepted_tangent_term_bound):
        worst_index = int(np.argmax(reconstruction_ratio))
        raise RuntimeError(
            "accepted PIRK tangent contributions do not reconstruct coordinate "
            f"mass-rate: max_abs={float(reconstruction_abs[worst_index]):.17e}, "
            f"roundoff_bound={float(accepted_tangent_term_bound[worst_index]):.17e}, "
            f"bound_ratio={float(reconstruction_ratio[worst_index]):.17e}, "
            f"cell_index={worst_index}, R={float(R[worst_index]):.17e}, "
            f"mass_dot_coordinate={float(chain_profiles['mass_dot_coordinate'][worst_index]):.17e}, "
            f"term_sum={float(accepted_tangent_term_sum_profile[worst_index]):.17e}"
        )

    at = lambda profile: interpolate_at_areal_radius(
        R, profile, FIXED_AREAL_RADIUS
    )
    mass_r = np.asarray(
        grid.cell_derivative_fourth(mass, parity=-1), dtype=float
    )
    Rr_tube = at(Rr)
    Rdot_tube = at(Rdot)
    mass_r_tube = at(mass_r)
    mass_dot_coordinate_tube = at(chain_profiles["mass_dot_coordinate"])
    tube_speed = -Rdot_tube / Rr_tube
    pointwise_sweep = mass_r_tube * tube_speed
    pointwise_rate = mass_dot_coordinate_tube + pointwise_sweep

    flux_average = 0.5 * (
        float(stage_rates["start"]["tube_stress_flux"])
        + float(stage_rates["primary_predictor"]["tube_stress_flux"])
    )
    pressure_average = 0.5 * (
        float(stage_rates["start"]["tube_radial_pressure_work"])
        + float(stage_rates["primary_predictor"]["tube_radial_pressure_work"])
    )
    predicted_with_pointwise_sweep = flux_average + pressure_average + pointwise_sweep
    pointwise_remainder = pointwise_rate - predicted_with_pointwise_sweep
    predicted_with_interpolant_sweep, interpolant_remainder = (
        _fixed_R_interpolant_balance(direct, flux_average, pressure_average)
    )
    composite_flux_average = 0.5 * (
        float(stage_rates["start"]["tube_stress_flux_composite_interp"])
        + float(stage_rates["primary_predictor"]["tube_stress_flux_composite_interp"])
    )
    composite_pressure_average = 0.5 * (
        float(stage_rates["start"]["tube_radial_pressure_work_composite_interp"])
        + float(stage_rates["primary_predictor"]["tube_radial_pressure_work_composite_interp"])
    )
    start_source = stage_rates["start"]["tube_ktheta_matter_source_rates"]
    predictor_source = stage_rates["primary_predictor"]["tube_ktheta_matter_source_rates"]
    source_stage_average = 0.5 * (
        float(start_source["total_mass_source_composite_at_R_star"])
        + float(predictor_source["total_mass_source_composite_at_R_star"])
    )
    pressure_stage_average = 0.5 * (
        float(start_source["total_radial_pressure_work_composite_at_R_star"])
        + float(predictor_source["total_radial_pressure_work_composite_at_R_star"])
    )
    density_stage_average = 0.5 * (
        float(start_source["total_density_correction_composite_at_R_star"])
        + float(predictor_source["total_density_correction_composite_at_R_star"])
    )
    identity_stage_average = 0.5 * (
        float(start_source["total_source_work_density_identity_residual_composite_at_R_star"])
        + float(predictor_source["total_source_work_density_identity_residual_composite_at_R_star"])
    )
    source_factor_stage_average = 0.5 * (
        float(start_source["total_mass_source_factorized_at_R_star"])
        + float(predictor_source["total_mass_source_factorized_at_R_star"])
    )
    pressure_factor_stage_average = 0.5 * (
        float(start_source["total_radial_pressure_work_factorized_at_R_star"])
        + float(predictor_source["total_radial_pressure_work_factorized_at_R_star"])
    )
    density_factor_stage_average = 0.5 * (
        float(start_source["total_density_correction_factorized_at_R_star"])
        + float(predictor_source["total_density_correction_factorized_at_R_star"])
    )
    identity_factor_stage_average = 0.5 * (
        float(start_source["total_source_work_density_identity_residual_factorized_at_R_star"])
        + float(predictor_source["total_source_work_density_identity_residual_factorized_at_R_star"])
    )
    predicted_with_composite_terms, interpolant_composite_remainder = (
        _fixed_R_interpolant_balance(
            direct, composite_flux_average, composite_pressure_average
        )
    )

    R1 = np.asarray(areal_radius(grid, state1.geometry), dtype=float)
    mass1 = np.asarray(misner_sharp(grid, state1.geometry), dtype=float)
    mass0_tube = interpolate_at_areal_radius(R, mass, FIXED_AREAL_RADIUS)
    mass1_tube = interpolate_at_areal_radius(R1, mass1, FIXED_AREAL_RADIUS)
    accepted_secant = (mass1_tube - mass0_tube) / float(step["dt"])

    tangent_terms_at_R_star = {
        name: at(profile) for name, profile in accepted_tangent_term_profiles.items()
    }
    tangent_coordinate_rate_from_terms = sum(tangent_terms_at_R_star.values())
    tangent_coordinate_rate_reference = float(direct["interpolated_mass_dot_coordinate"])
    tangent_coordinate_reconstruction_residual = (
        tangent_coordinate_rate_reference - tangent_coordinate_rate_from_terms
    )
    tangent_coordinate_scale = (
        abs(tangent_coordinate_rate_reference)
        + sum(abs(value) for value in tangent_terms_at_R_star.values())
    )
    tangent_coordinate_bound = 128.0 * np.finfo(np.float64).eps * max(
        tangent_coordinate_scale, float(np.finfo(np.float64).tiny)
    )
    if abs(tangent_coordinate_reconstruction_residual) > tangent_coordinate_bound:
        raise RuntimeError(
            "interpolated accepted tangent pieces do not reconstruct coordinate mass-rate"
        )

    tangent_fixed_R_prediction_composite = (
        composite_flux_average + composite_pressure_average + float(direct["fixed_R_sweep"])
    )
    tangent_fixed_R_residual_composite = (
        float(direct["mass_dot_fixed_R"]) - tangent_fixed_R_prediction_composite
    )
    tangent_fixed_R_prediction_factorized = (
        flux_average + pressure_average + float(direct["fixed_R_sweep"])
    )
    tangent_fixed_R_residual_factorized = (
        float(direct["mass_dot_fixed_R"]) - tangent_fixed_R_prediction_factorized
    )
    tangent_reconciliation = {
        "accepted_state_secant_mass_rate": float(accepted_secant),
        "accepted_update_interpolant_mass_rate": float(direct["mass_dot_fixed_R"]),
        "accepted_secant_minus_update_tangent": float(
            accepted_secant - float(direct["mass_dot_fixed_R"])
        ),
        "R_star": {
            "mass_dot_coordinate_terms": tangent_terms_at_R_star,
            "mass_dot_coordinate_from_terms": float(tangent_coordinate_rate_from_terms),
            "mass_dot_coordinate_from_direct_interpolant": tangent_coordinate_rate_reference,
            "coordinate_rate_reconstruction_residual": float(
                tangent_coordinate_reconstruction_residual
            ),
            "coordinate_rate_reconstruction_roundoff_bound": float(tangent_coordinate_bound),
            "fixed_R_sweep_from_same_interpolant": float(direct["fixed_R_sweep"]),
            "accepted_update_fixed_R_tangent": float(direct["mass_dot_fixed_R"]),
            "composite_stage_average_stress_flux": float(composite_flux_average),
            "composite_stage_average_pressure_work": float(composite_pressure_average),
            "factorized_stage_average_stress_flux": float(flux_average),
            "factorized_stage_average_pressure_work": float(pressure_average),
            "fixed_R_prediction_composite_terms": float(tangent_fixed_R_prediction_composite),
            "fixed_R_balance_remainder_composite_terms": float(tangent_fixed_R_residual_composite),
            "stored_fixed_R_residual_composite_interpolation": float(
                interpolant_composite_remainder
            ),
            "fixed_R_prediction_factorized_terms": float(tangent_fixed_R_prediction_factorized),
            "fixed_R_balance_remainder_factorized_terms": float(tangent_fixed_R_residual_factorized),
            "stored_fixed_R_residual_factor_interpolation": float(interpolant_remainder),
        },
        "bracket_cells": {},
    }
    for cell_name, cell_index in (
        ("left_cell", int(direct["left_index"])),
        ("right_cell", int(direct["right_index"])),
    ):
        cell_terms = {
            name: float(profile[cell_index])
            for name, profile in accepted_tangent_term_profiles.items()
        }
        cell_coordinate_rate = float(chain_profiles["mass_dot_coordinate"][cell_index])
        cell_term_sum = sum(cell_terms.values())
        tangent_reconciliation["bracket_cells"][cell_name] = {
            "index": cell_index,
            "R": float(R[cell_index]),
            "mass_dot_coordinate_terms": cell_terms,
            "mass_dot_coordinate_from_terms": float(cell_term_sum),
            "mass_dot_coordinate_from_chain_rule": cell_coordinate_rate,
            "coordinate_rate_reconstruction_residual": float(cell_coordinate_rate - cell_term_sum),
        }

    return {
        "t_start": float(step["t_before"]),
        "dt": float(step["dt"]),
        "accepted_state_tube_mass_before": float(mass0_tube),
        "accepted_state_tube_mass_after": float(mass1_tube),
        "accepted_state_secant_mass_rate": float(accepted_secant),
        "pirk_update_tangent_pointwise_mass_rate": float(pointwise_rate),
        "pirk_update_tangent_interpolant_mass_rate": float(direct["mass_dot_fixed_R"]),
        "accepted_secant_minus_pirk_tangent_pointwise_rate": float(accepted_secant - pointwise_rate),
        "accepted_secant_minus_pirk_tangent_interpolant_rate": float(
            accepted_secant - float(direct["mass_dot_fixed_R"])
        ),
        "stage_average_stress_flux": float(flux_average),
        "stage_average_radial_pressure_work": float(pressure_average),
        "pirk_update_tangent_pointwise_sweep": float(pointwise_sweep),
        "pirk_update_tangent_interpolant_sweep": float(direct["fixed_R_sweep"]),
        "pirk_update_tangent_predicted_rate_using_pointwise_sweep": float(
            predicted_with_pointwise_sweep
        ),
        "pirk_update_tangent_predicted_rate_using_interpolant_sweep": float(
            predicted_with_interpolant_sweep
        ),
        "pirk_update_tangent_interpolant_sweep_difference_from_pointwise": float(
            float(direct["fixed_R_sweep"]) - pointwise_sweep
        ),
        "pirk_update_tangent_pointwise_balance_remainder": float(pointwise_remainder),
        "pirk_update_tangent_interpolant_balance_remainder": float(interpolant_remainder),
        "pirk_update_tangent_composite_interpolated_flux_average": float(composite_flux_average),
        "pirk_update_tangent_composite_interpolated_pressure_work_average": float(composite_pressure_average),
        "pirk_update_tangent_predicted_rate_using_composite_interpolated_terms": float(
            predicted_with_composite_terms
        ),
        "pirk_update_tangent_interpolant_balance_remainder_composite_interpolation": float(
            interpolant_composite_remainder
        ),
        "accepted_update_tangent_termwise_reconciliation": tangent_reconciliation,
        "pirk_update_tangent_interpolant_coordinate_remainder": float(
            float(direct["interpolated_mass_dot_coordinate"]) - flux_average - pressure_average
        ),
        "pirk_update_tangent_einstein_source_mass_rate_accepted_weights": float(source_mass_rate_accepted_weights),
        "pirk_update_tangent_einstein_source_mass_rate_start_predictor_composite_average": float(source_stage_average),
        "pirk_update_tangent_einstein_pressure_work_start_predictor_composite_average": float(pressure_stage_average),
        "pirk_update_tangent_einstein_density_correction_start_predictor_composite_average": float(density_stage_average),
        "pirk_update_tangent_einstein_source_identity_start_predictor_composite_average": float(identity_stage_average),
        "pirk_update_tangent_einstein_source_mass_rate_start_predictor_factorized_average": float(source_factor_stage_average),
        "pirk_update_tangent_einstein_pressure_work_start_predictor_factorized_average": float(pressure_factor_stage_average),
        "pirk_update_tangent_einstein_density_correction_start_predictor_factorized_average": float(density_factor_stage_average),
        "pirk_update_tangent_einstein_source_identity_start_predictor_factorized_average": float(identity_factor_stage_average),
        "pirk_update_tangent_einstein_source_update_minus_stage_composite_average_diagnostic_only": float(source_mass_rate_accepted_weights - source_stage_average),
        "pirk_update_tangent_einstein_source_roundoff_scale_absmax": float(max(
            np.max(np.abs(update_einstein_source["roundoff_scale"]["Aa"])),
            np.max(np.abs(update_einstein_source["roundoff_scale"]["K"])),
        )),
        "interpolant_bracket_left_index": int(direct["left_index"]),
        "interpolant_bracket_right_index": int(direct["right_index"]),
        "interpolant_theta": float(direct["theta"]),
        "final_primary_l2_input_geometry": final_l2["geometry"],
        "final_primary_l2_Aa_profile": final_l2["Aa"],
        "final_primary_l2_K_profile": final_l2["K"],
        "rate_profiles": {
            "R": R,
            "Rr": Rr,
            "Rdot_from_accepted_pirk_coefficients": Rdot,
            "Rdot_r_from_accepted_pirk_coefficients": Rdot_r,
            "a_dot_from_accepted_pirk_coefficients": a_dot,
            "b_dot_from_accepted_pirk_coefficients": b_dot,
            "X_dot_from_accepted_pirk_coefficients": X_dot,
            "Aa_dot_from_accepted_pirk_coefficients": Aa_dot,
            "K_dot_from_accepted_pirk_coefficients": K_dot,
            "gamma_inv_dot_from_accepted_pirk_coefficients": gamma_inv_dot,
            "Ktheta_dot_from_accepted_pirk_coefficients": Ktheta_dot,
            "mass": chain_profiles["mass"],
            "mass_dot_coordinate_from_accepted_pirk_coefficients": chain_profiles[
                "mass_dot_coordinate"
            ],
        },
    }



def _pirk_accepted_update_einstein_source_arrays(geometry_stages) -> dict:
    """Accepted Aa/K Einstein source = 1/2 start L3 + 1/2 full predictor L3."""
    if len(geometry_stages) != 3:
        raise RuntimeError("accepted Einstein source update needs three geometry stages")
    out = {"source_l3": {}, "roundoff_scale": {}}
    shape = None
    for name in ("Aa", "K"):
        try:
            left = np.asarray(geometry_stages[0]["einstein_source_l3"][name], dtype=float)
            right = np.asarray(geometry_stages[2]["einstein_source_l3"][name], dtype=float)
            left_scale = np.asarray(geometry_stages[0]["einstein_source_roundoff_scale"][name], dtype=float)
            right_scale = np.asarray(geometry_stages[2]["einstein_source_roundoff_scale"][name], dtype=float)
        except (KeyError, TypeError) as exc:
            raise RuntimeError(f"accepted Einstein source update missing start/predictor {name} data") from exc
        if left.ndim != 1 or right.shape != left.shape or left_scale.shape != left.shape or right_scale.shape != left.shape:
            raise RuntimeError(f"accepted Einstein source update {name} arrays have incompatible shapes")
        if not all(np.all(np.isfinite(v)) for v in (left, right, left_scale, right_scale)):
            raise RuntimeError(f"accepted Einstein source update {name} arrays contain non-finite values")
        if np.any(left_scale < 0.0) or np.any(right_scale < 0.0):
            raise RuntimeError(f"accepted Einstein source update {name} has negative operand scale")
        if shape is None:
            shape = left.shape
        elif left.shape != shape:
            raise RuntimeError("accepted Einstein Aa/K source arrays have different shapes")
        out["source_l3"][name] = 0.5 * (left + right)
        out["roundoff_scale"][name] = 0.5 * (left_scale + right_scale)
    return out



def _serialize_step(kernel, step: dict):
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
    accepted_pirk_update_tangent = _accepted_pirk_update_tangent(
        step, fixed_R_stage_rates
    )
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
        "accepted_pirk_update_tangent": accepted_pirk_update_tangent,
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
            "tube_chain_rule_mass_dot_term_Rdot",
            "tube_chain_rule_mass_dot_term_gamma_inv_dot",
            "tube_chain_rule_mass_dot_term_Rdot_r",
            "tube_chain_rule_mass_dot_term_Ktheta_geometry",
            "tube_chain_rule_mass_dot_term_Ktheta_dot",
            "tube_chain_rule_mass_dot_coordinate_term_sum",
            "tube_chain_rule_mass_dot_coordinate_decomposition_residual",
            "tube_chain_rule_mass_dot_term_Ktheta_matter_source",
            "tube_chain_rule_mass_dot_term_Ktheta_dot_minus_explicit_matter_source",
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
        row = {
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
        }
        tangent = step.get("accepted_pirk_update_tangent")
        if isinstance(tangent, dict) and isinstance(
            tangent.get("accepted_update_tangent_termwise_reconciliation"), dict
        ):
            row["accepted_update_tangent_termwise_reconciliation"] = (
                tangent["accepted_update_tangent_termwise_reconciliation"]
            )
        rows.append(row)
    return rows


def _numeric_tree_difference(charged, neutral):
    """Subtract observables while preserving witness identity metadata."""
    if isinstance(charged, dict) and isinstance(neutral, dict):
        out = {}
        identity_keys = {"target_r", "left_index", "right_index", "index"}
        for key in sorted(set(charged) & set(neutral)):
            # Bounds and bracket/cell identifiers are per-arm metadata, not
            # signed observables. Keep identifiers only when both arms refer
            # to the same geometric target and the same native bracket cells.
            key_folded = key.casefold()
            if key_folded.endswith("_bound"):
                continue
            if key_folded in identity_keys:
                if charged[key] != neutral[key]:
                    raise RuntimeError(
                        "paired witness identity metadata mismatch at "
                        f"{key}: charged={charged[key]!r}, neutral={neutral[key]!r}"
                    )
                out[key] = charged[key]
                continue
            value = _numeric_tree_difference(charged[key], neutral[key])
            if value is not None:
                out[key] = value
        return out
    numeric = (int, float, np.number)
    if (
        isinstance(charged, numeric) and not isinstance(charged, (bool, np.bool_))
        and isinstance(neutral, numeric) and not isinstance(neutral, (bool, np.bool_))
    ):
        return float(charged) - float(neutral)
    return None


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
        "start_stage_tube_chain_rule_mass_dot_term_Rdot",
        "primary_predictor_stage_tube_chain_rule_mass_dot_term_Rdot",
        "start_stage_tube_chain_rule_mass_dot_term_gamma_inv_dot",
        "primary_predictor_stage_tube_chain_rule_mass_dot_term_gamma_inv_dot",
        "start_stage_tube_chain_rule_mass_dot_term_Rdot_r",
        "primary_predictor_stage_tube_chain_rule_mass_dot_term_Rdot_r",
        "start_stage_tube_chain_rule_mass_dot_term_Ktheta_geometry",
        "primary_predictor_stage_tube_chain_rule_mass_dot_term_Ktheta_geometry",
        "start_stage_tube_chain_rule_mass_dot_term_Ktheta_dot",
        "primary_predictor_stage_tube_chain_rule_mass_dot_term_Ktheta_dot",
        "start_stage_tube_chain_rule_mass_dot_coordinate_term_sum",
        "primary_predictor_stage_tube_chain_rule_mass_dot_coordinate_term_sum",
        "start_stage_tube_chain_rule_mass_dot_coordinate_decomposition_residual",
        "primary_predictor_stage_tube_chain_rule_mass_dot_coordinate_decomposition_residual",
        "start_stage_tube_chain_rule_mass_dot_term_Ktheta_matter_source",
        "primary_predictor_stage_tube_chain_rule_mass_dot_term_Ktheta_matter_source",
        "start_stage_tube_chain_rule_mass_dot_term_Ktheta_dot_minus_explicit_matter_source",
        "primary_predictor_stage_tube_chain_rule_mass_dot_term_Ktheta_dot_minus_explicit_matter_source",
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
        for stage_name in ("start", "primary_predictor"):
            charged_components = c[f"{stage_name}_stage_components"]
            neutral_components = n[f"{stage_name}_stage_components"]
            for witness_key in (
                "tube_mass_gradient_identity_witness",
                "tube_chain_rule_flux_work_operator_comparison",
                "tube_hm_constraint_adjusted_fixed_R_rate",
            ):
                charged_witness = charged_components.get(witness_key)
                neutral_witness = neutral_components.get(witness_key)
                if isinstance(charged_witness, dict) and isinstance(neutral_witness, dict):
                    row[
                        f"charged_minus_neutral_{stage_name}_stage_{witness_key}"
                    ] = _numeric_tree_difference(charged_witness, neutral_witness)
        charged_tangent = c.get("accepted_update_tangent_termwise_reconciliation")
        neutral_tangent = n.get("accepted_update_tangent_termwise_reconciliation")
        if isinstance(charged_tangent, dict) and isinstance(neutral_tangent, dict):
            row["charged_minus_neutral_accepted_update_tangent_termwise_reconciliation"] = (
                _numeric_tree_difference(charged_tangent, neutral_tangent)
            )
        elif charged_tangent is not None or neutral_tangent is not None:
            raise RuntimeError(
                f"{label}: accepted tangent termwise reconciliation missing in one charge arm"
            )
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
        "The live primary-L3 K and Aa matter-source increments are checked against independently reconstructed sector stress sources at a 64-epsilon operand-scaled bound.",
        "The Ktheta-source, radial-pressure-work, and density-correction identity is reported under both interpolate-after-product and multiply-interpolated-factors conventions. These discrete operators can differ; only like-for-like identity residuals and the Rdot + alpha*R*Ktheta closure are interpreted together.",
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
