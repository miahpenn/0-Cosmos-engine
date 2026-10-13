"""Bounded same-state A/B for archive-fluid stress normalization.

Runs one N=40 step from one frozen initializer-produced state twice. The only
difference is whether archived fluid stress is divided by KAPPA before it is
handed to Einstein/CMC geometry. This is a diagnostic, not a campaign or
physics-admission gate. It also measures the current production source-pair
identity in physical stress units; known scalar/DM source-law work remains
separate and is not patched here.
"""
from __future__ import annotations

import copy
import json
import math
import subprocess

import numpy as np

import engine.matter_system as matter_system
import engine.scalar_system as scalar_system
import engine.stress_energy as stress_energy
import engine.v55_matter as v55_matter
from engine.production_kernel import V55ProductionKernel, BETA_DM
from engine.matter_system import (
    Species,
    _species_rhs_components,
    _metric_arrays,
    geometry_metric_derivatives,
    project_species,
)
from engine.scalar_system import scalar_rhs_arrays, KAPPA
from engine.v55_matter import metric_slice_from_q, dm_density
from engine import v55_pirk_adapter as adapter

COMPONENTS = ("rho", "pr", "pt", "j")
SURFACE_AREAL_RADIUS = 10.0


def unscaled_archive_fluid_projection(metric, state):
    pieces = (
        project_species(metric, state.dark_matter, Species.DARK_MATTER),
        project_species(metric, state.baryons, Species.BARYON),
        project_species(metric, state.radiation, Species.RADIATION),
    )
    return {key: sum(piece[key] for piece in pieces) for key in COMPONENTS}


def assert_same_initial_state(left, right):
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        np.testing.assert_array_equal(
            getattr(left.geometry, name), getattr(right.geometry, name),
            err_msg=f"frozen initial geometry differs in {name}",
        )
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        np.testing.assert_array_equal(
            getattr(left.scalars, name), getattr(right.scalars, name),
            err_msg=f"frozen initial scalar differs in {name}",
        )
    for species_name in ("dark_matter", "baryons", "radiation"):
        for name in ("rest", "energy_t", "momentum_r"):
            np.testing.assert_array_equal(
                getattr(getattr(left.matter, species_name), name),
                getattr(getattr(right.matter, species_name), name),
                err_msg=f"frozen initial {species_name}.{name} differs",
            )


def radiation_margin(state):
    metric = metric_slice_from_q(state.grid, state.geometry)
    metrics = _metric_arrays(metric)
    sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics], dtype=float)
    gamma_rr_inv = np.asarray([m.gamma_rr_inv for m in metrics], dtype=float)
    energy = state.matter.radiation.energy_t / sqrt_gamma
    momentum = state.matter.radiation.momentum_r / sqrt_gamma
    abs_momentum = np.sqrt(np.maximum(gamma_rr_inv * momentum * momentum, 0.0))
    ratio = abs_momentum / np.maximum(np.abs(energy), 1.0e-300)
    return {
        "minimum_E_minus_absS": float(np.min(energy - abs_momentum)),
        "maximum_absS_over_E": float(np.max(ratio)),
        "outer_absS_over_E": float(ratio[-1]),
    }


def source_pair_in_physical_units(state):
    """Measure source increments from current-main equations, before flux terms."""
    grid, geom, fields = state.grid, state.geometry, state.scalars
    metric = metric_slice_from_q(grid, geom)
    rho_dm = dm_density(metric, state.matter)
    scalar_with = scalar_rhs_arrays(
        grid, geom, fields, beta_dm=BETA_DM, rho_dm=rho_dm
    )
    scalar_without = scalar_rhs_arrays(
        grid, geom, fields, beta_dm=0.0, rho_dm=rho_dm
    )
    delta_pi = scalar_with.Pi - scalar_without.Pi
    dphi_r = grid.cell_derivative_fourth(fields.phi, parity=1)
    derivatives = geometry_metric_derivatives(
        metric,
        lambda values, parity: grid.cell_derivative_fourth(values, parity=parity),
    )
    raw_with = _species_rhs_components(
        metric, derivatives, state.matter.dark_matter, Species.DARK_MATTER,
        dphi_t=scalar_with.phi, dphi_r=dphi_r, beta_dm=BETA_DM,
    )
    raw_without = _species_rhs_components(
        metric, derivatives, state.matter.dark_matter, Species.DARK_MATTER,
        dphi_t=scalar_with.phi, dphi_r=dphi_r, beta_dm=0.0,
    )
    q_energy = raw_with[3] - raw_without[3]
    q_radial = raw_with[4] - raw_without[4]
    sqrt_gamma = np.asarray([m.sqrt_gamma for m in _metric_arrays(metric)])
    # The source contribution to scalar rho is Pi * Delta(Pi); the source
    # contribution to scalar j=-Pi*d_r(phi) is -Delta(Pi)*d_r(phi).
    scalar_energy = sqrt_gamma * fields.Pi * delta_pi
    scalar_current = -sqrt_gamma * delta_pi * dphi_r
    net_energy = (scalar_energy + q_energy) / KAPPA
    net_current = (scalar_current + q_radial) / KAPPA

    def stats(values, scale):
        return {
            "linf": float(np.max(np.abs(values))),
            "relative_linf": float(
                np.max(np.abs(values))
                / max(float(np.max(np.abs(scale))), 1.0e-300)
            ),
        }

    return {
        "units": "physical Einstein-source units after common 1/KAPPA conversion",
        "scalar_energy_source_linf": float(np.max(np.abs(scalar_energy / KAPPA))),
        "dark_matter_energy_source_linf": float(np.max(np.abs(q_energy / KAPPA))),
        "net_energy_source": stats(net_energy, np.maximum(
            np.abs(scalar_energy / KAPPA), np.abs(q_energy / KAPPA)
        )),
        "scalar_radial_current_source_linf": float(np.max(np.abs(scalar_current / KAPPA))),
        "dark_matter_radial_current_source_linf": float(np.max(np.abs(q_radial / KAPPA))),
        "net_radial_current_source": stats(net_current, np.maximum(
            np.abs(scalar_current / KAPPA), np.abs(q_radial / KAPPA)
        )),
        "note": (
            "This records the current production source laws as written. "
            "A nonzero net is evidence of an unclosed source identity, not "
            "a reason to alter a coefficient in this normalization PR."
        ),
    }


def summarize(kernel, state, diag, target_areal_radius):
    R = np.asarray(diag["areal_radius"], dtype=float)
    M = np.asarray(diag["misner_sharp"], dtype=float)
    if np.any(np.diff(R) <= 0.0):
        raise AssertionError("areal radius is not strictly increasing on frozen A/B slice")
    mass_fixed_R = float(np.interp(target_areal_radius, R, M))
    idx = int(np.argmin(np.abs(R - target_areal_radius)))
    total = stress_energy.assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    _, _, moving = adapter.vendor_modules()
    explicit = moving.moving_puncture_explicit_rhs(state.grid, state.geometry)
    Rdot = R * (
        0.5 * explicit["b"] / state.geometry.b
        - explicit["X"] / state.geometry.X
    )
    coordinate_flux = 4.0 * math.pi * R[idx] ** 2 * (
        state.geometry.alpha[idx]
        * (state.geometry.X[idx] ** 2 / state.geometry.a[idx])
        * total.j[idx]
        - state.geometry.beta[idx] * total.rho[idx]
    )
    work = -4.0 * math.pi * R[idx] ** 2 * total.pr[idx] * Rdot[idx]
    return {
        "t": float(state.t),
        "H_full_absmax": float(diag["hamiltonian_full_absmax"]),
        "M_full_absmax": float(diag["momentum_full_absmax"]),
        "connection_full_absmax": float(diag["connection_full_absmax"]),
        "determinant_constraint_max": float(diag["determinant_constraint_max"]),
        "Lambda_linf": float(np.max(np.abs(state.geometry.Lambda))),
        "lapse_min": float(np.min(state.geometry.alpha)),
        "fixed_areal_R": float(target_areal_radius),
        "M_MS_at_fixed_areal_R": mass_fixed_R,
        "fixed_areal_R_current_flux": float(coordinate_flux),
        "fixed_areal_R_work_term": float(work),
        "surface_index_for_current": idx,
        "surface_areal_R_sample": float(R[idx]),
        "radiation": radiation_margin(state),
    }


def run_arm(kernel, state, projection, target_areal_radius, dt):
    v55_matter.total_fluid_projection = projection
    stress_energy.total_fluid_projection = projection
    before_diag = kernel.diagnostics(state, profiles=True)
    before = summarize(kernel, state, before_diag, target_areal_radius)
    advanced = kernel.step(state, dt)
    advanced.geometry.assert_finite_positive()
    after_diag = kernel.diagnostics(advanced, profiles=True)
    after = summarize(kernel, advanced, after_diag, target_areal_radius)
    if after["radiation"]["maximum_absS_over_E"] >= 1.0:
        raise AssertionError("radiation left the strict E>|S| cone in bounded A/B step")
    after["fixed_areal_R_M_MS_rate_one_step"] = (
        after["M_MS_at_fixed_areal_R"] - before["M_MS_at_fixed_areal_R"]
    ) / dt
    after["dt"] = float(dt)
    return {"initial": before, "after_one_step": after}, advanced


def state_difference(left, right):
    keys = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    return {
        f"geometry_{name}_max_abs": float(
            np.max(np.abs(getattr(left.geometry, name) - getattr(right.geometry, name)))
        )
        for name in keys
    }


def git_blob(path):
    try:
        return subprocess.check_output(
            ["git", "rev-parse", f"HEAD:{path}"], text=True
        ).strip()
    except Exception:
        return "unavailable"


def main():
    original_projection = v55_matter.total_fluid_projection
    initial = V55ProductionKernel().initialize(
        resolution=40, r_max=40.0, include_radiation=True
    )
    baseline_state = copy.deepcopy(initial)
    candidate_state = copy.deepcopy(initial)
    assert_same_initial_state(baseline_state, candidate_state)
    base_diag = V55ProductionKernel().diagnostics(initial, profiles=True)
    target_R = SURFACE_AREAL_RADIUS
    R_profile = np.asarray(base_diag["areal_radius"], dtype=float)
    if not (R_profile[0] <= target_R <= R_profile[-1]):
        target_R = float(R_profile[len(R_profile) // 4])
    dt = 0.03 * float(initial.grid.dr)
    kernel_base = V55ProductionKernel()
    kernel_candidate = V55ProductionKernel()

    try:
        base, base_final = run_arm(
            kernel_base, baseline_state, unscaled_archive_fluid_projection, target_R, dt
        )
        candidate, cand_final = run_arm(
            kernel_candidate, candidate_state, original_projection, target_R, dt
        )
    finally:
        v55_matter.total_fluid_projection = original_projection
        stress_energy.total_fluid_projection = original_projection

    report = {
        "status": "BOUNDED_FROZEN_STATE_AB_DIAGNOSTIC_ONLY",
        "comparison": {
            "baseline": "raw archive-normalized fluid projections summed without 1/KAPPA",
            "candidate": "all fluid rho, pr, pt, j projections divided by KAPPA once at Einstein boundary",
            "resolution": int(initial.grid.n),
            "r_max": float(initial.grid.r_max),
            "dt": float(dt),
            "common_initial_state": True,
            "initial_state_exact_match": True,
            "target_areal_radius": float(target_R),
            "no_long_evolution": True,
        },
        "source_identity_at_common_frozen_state": source_pair_in_physical_units(initial),
        "baseline": base,
        "candidate": candidate,
        "candidate_minus_baseline_after_one_step_geometry_max_abs": state_difference(
            cand_final, base_final
        ),
        "source_pins": {
            "HEAD": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
            "engine/v55_matter.py": git_blob("engine/v55_matter.py"),
            "engine/scalar_system.py": git_blob("engine/scalar_system.py"),
            "engine/matter_system.py": git_blob("engine/matter_system.py"),
            "engine/stress_energy.py": git_blob("engine/stress_energy.py"),
            "vendor_gitlink": git_blob("vendor/bb-palatini-unified-r0"),
        },
    }
    payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
    print("ARCHIVE_FLUID_STRESS_NORMALIZATION_AB_JSON_BEGIN", flush=True)
    print(payload, flush=True)
    print("ARCHIVE_FLUID_STRESS_NORMALIZATION_AB_JSON_END", flush=True)


if __name__ == "__main__":
    main()
