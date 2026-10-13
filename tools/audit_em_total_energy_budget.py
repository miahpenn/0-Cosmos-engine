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
from pathlib import Path
import subprocess
import tempfile

import numpy as np

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
DENSITY_PROFILE = ((0.0, 2.0e-5), (R_MAX, 2.0e-5))
CHARGE_TO_MASS_ARMS = (0.0, 0.25)
FIXED_AREAL_RADIUS = 10.0

# Blob IDs from main 5c51bf509504365df394e6cca96bf1d117f4614d. The audit is
# deliberately tied to these operators; update the pin only after a new source
# review and a new exact-head result, never silently on a moving main branch.
SOURCE_PINS = {
    "engine/campaign.py": "119942d9debf63b99997cda96c541754e7e8aded",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/matter_system.py": "05e23d2368120ea8a1a6aad3fa438625d5380750",
    "engine/production_kernel.py": "4e0bdccc2ac614c48c34521663a6d53b404b6747",
    "engine/radial_maxwell.py": "81269f1d1334c46c4c2f44ba816d792722c0a38e",
    "engine/run_production.py": "726f3fd09d46142862574d348aa0755f83bc0612",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/v55_matter.py": "acc8f8ceae28d82511c9a5ad9d2ca670e09fb0f7",
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
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    return {"candidate_head": head, "production_blob_pins": observed}


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


def _worldtube_sample(kernel, state) -> dict:
    """Evaluate native MS and EM-inclusive stress/flux/work on R=10."""
    grid, geom = state.grid, state.geometry
    R = np.asarray(areal_radius(grid, geom), dtype=float)
    if not np.all(np.diff(R) > 0.0):
        raise ValueError("native areal radius is not strictly increasing")
    Rr = np.asarray(areal_radius_derivative(grid, geom), dtype=float)
    chi = np.asarray(invariant_chi(grid, geom), dtype=float)
    mass = np.asarray(misner_sharp(grid, geom), dtype=float)
    Mr = np.asarray(grid.cell_derivative_fourth(mass, parity=-1), dtype=float)

    _, _, moving = adapter.vendor_modules()
    explicit = moving.moving_puncture_explicit_rhs(grid, geom)
    Rdot = R * (
        0.5 * np.asarray(explicit["b"], dtype=float) / geom.b
        - np.asarray(explicit["X"], dtype=float) / geom.X
    )

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
        "constraint_witnesses": constraint_witnesses,
    }


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
            expected_energy_rhs = (
                -np.diff(flux[:, 1]) / dr
                + components["source_energy"]
                + components["dm_exchange_energy"]
                + force_e
            )
            components["rhs_energy_profile"] = _array(rhs.energy_t)
            components["rhs_energy_from_named_terms_profile"] = _array(expected_energy_rhs)
            components["rhs_energy_component_residual_profile"] = (
                _array(rhs.energy_t) - _array(expected_energy_rhs)
            )
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
                    for name in ("a", "b", "X", "alpha", "beta")
                },
                "explicit": {
                    name: _array(result["explicit"][name])
                    for name in ("a", "b", "X")
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
        "matter_lorentz_power_cancels_in_matter_plus_em_sum": (
            0.5 * (sum(
                stages[k]["components"][name]["lorentz_energy_rate"]
                for name in stages[k]["components"]
            ) for k in (0, 1))
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
    }
    report = _json_safe(report)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return report


def write_step_summary(report: dict) -> None:
    summary_path = Path(str(__import__("os").environ.get("GITHUB_STEP_SUMMARY", "")))
    if not str(summary_path):
        return
    summary = report["summary"]
    lines = [
        "## Joined radial-EM / matter / Misner–Sharp ledger",
        "",
        "Source-pinned, research-only paired run; q/m=0 versus 0.25; N=24/48/96; t=0.03; Gaussian c=1 convention.",
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
    lines.extend([
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
        "candidate_head": report["provenance"]["candidate_head"],
        "source_pins_passed": all(
            value["match"]
            for value in report["provenance"]["production_blob_pins"].values()
        ),
        "paired_neutral_vs_charged": report["summary"]["paired_neutral_vs_charged"],
        "full_report_path": str(args.output.resolve()),
    }
    print("EM_TOTAL_ENERGY_LEDGER_REPORT=" + json.dumps(concise, sort_keys=True, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
