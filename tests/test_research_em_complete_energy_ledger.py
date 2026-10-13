"""Joined source ledger for the opt-in radial-Maxwell research path.

This observer records, without fitting a global conservation scalar:
- actual-stage species face fluxes and named RHS sources;
- the live radial-EM charge/Lorentz/metric/boundary work terms;
- the EM-inclusive Misner-Sharp mass-rate, stress-flux and pressure-work terms.

It is deliberately a one-step, zero-shift diagnostic. The report is not
scientific admission, a physical-normalization decision, or a claim that a
single global energy is conserved on a dynamical spacetime.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

import engine.matter_rhs as matter_rhs_module
from engine.matter_system import _metric_arrays, initialize_dust
from engine.production_kernel import V55ProductionKernel
from engine.radial_maxwell import RadialMaxwellState
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.worldtube import areal_radius_derivative


FOUR_PI = 4.0 * math.pi
ETA_ARMS = (0.0, 0.25)
RESOLUTIONS = (24, 48, 96)
R_MAX = 12.0
DENSITY = 2.0e-5
TARGET_AREAL_RADIUS = 10.0


def _max_abs(values) -> float:
    return float(np.max(np.abs(np.asarray(values, dtype=float))))


def _geometry_arrays(geometry) -> dict[str, np.ndarray]:
    return {
        name: np.asarray(getattr(geometry, name), dtype=float).copy()
        for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda")
    }


def _weight(grid, geometry) -> np.ndarray:
    r = np.asarray(grid.centers, dtype=float)
    return (
        np.asarray(geometry.X, dtype=float)
        * np.sqrt(np.asarray(geometry.a, dtype=float))
        / (np.square(r) * np.asarray(geometry.b, dtype=float))
    )


def _field_energy(grid, geometry, electromagnetic) -> float:
    q_faces = np.asarray(electromagnetic.enclosed_charge_faces, dtype=float)
    q_cell = 0.5 * (q_faces[:-1] + q_faces[1:])
    return 0.5 * grid.dr * float(np.sum(_weight(grid, geometry) * q_cell**2))


def _sample_at_areal_radius(radius: np.ndarray, values: np.ndarray, target: float):
    """Linear interpolation at target R when a positive ordered bracket exists."""
    R = np.asarray(radius, dtype=float)
    data = np.asarray(values, dtype=float)
    right = int(np.searchsorted(R, target, side="right"))
    if 0 < right < R.size and R[right] > R[right - 1] and R[right - 1] <= target <= R[right]:
        left = right - 1
        fraction = float((target - R[left]) / (R[right] - R[left]))
        return float((1.0 - fraction) * data[left] + fraction * data[right]), {
            "method": "linear_areal_radius_bracket",
            "left_index": left,
            "right_index": right,
            "left_R": float(R[left]),
            "right_R": float(R[right]),
            "fraction": fraction,
        }
    idx = int(np.argmin(np.abs(R - target)))
    return float(data[idx]), {
        "method": "nearest_cell_fallback",
        "index": idx,
        "sample_R": float(R[idx]),
    }


def _ms_rate_ledger(kernel, grid, geometry, scalars, matter, electromagnetic, terms):
    """Capture a signed quasi-local rate ledger from one real geometry stage.

    The finite-time fixed-R sweep and fixed-coordinate rate are both reported.
    Their residuals are observers only; no acceptance threshold is fitted to
    them and no claim of a timelike-Killing global energy is made.
    """
    explicit = terms["explicit"]
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    alpha = np.asarray(geometry.alpha, dtype=float)
    r = np.asarray(grid.centers, dtype=float)

    adot = np.asarray(explicit["a"], dtype=float)
    bdot = np.asarray(explicit["b"], dtype=float)
    Xdot = np.asarray(explicit["X"], dtype=float)
    Aadot = (
        np.asarray(terms["primary_l2"]["Aa"], dtype=float)
        + np.asarray(terms["primary_l3"]["Aa"], dtype=float)
    )
    Kdot = (
        np.asarray(terms["primary_l2"]["K"], dtype=float)
        + np.asarray(terms["primary_l3"]["K"], dtype=float)
    )

    R = r * np.sqrt(b) / X
    Rr = areal_radius_derivative(grid, geometry)
    Rdot = R * (0.5 * bdot / b - Xdot / X)
    Rr_dot = grid.cell_derivative_fourth(Rdot, parity=-1)

    Ktheta = np.asarray(geometry.K, dtype=float) / 3.0 - np.asarray(geometry.Aa, dtype=float) / 2.0
    Ktheta_dot = Kdot / 3.0 - Aadot / 2.0
    inv_rr = X**2 / a
    inv_rr_dot = 2.0 * X * Xdot / a - X**2 * adot / np.square(a)

    chi = inv_rr * np.square(Rr) - np.square(R * Ktheta)
    chi_dot = (
        inv_rr_dot * np.square(Rr)
        + 2.0 * inv_rr * Rr * Rr_dot
        - 2.0 * R * Ktheta * (Rdot * Ktheta + R * Ktheta_dot)
    )
    mass = 0.5 * R * (1.0 - chi)
    mass_dot = 0.5 * Rdot * (1.0 - chi) - 0.5 * R * chi_dot
    mass_r = grid.cell_derivative_fourth(mass, parity=-1)

    projection = kernel._electromagnetic_projection(grid, geometry, electromagnetic)
    total = assemble_total_stress_energy(
        grid, geometry, scalars, matter,
        electromagnetic_projection=projection,
    )
    flux_T = FOUR_PI * np.square(R) * (
        -alpha * inv_rr * np.asarray(total.j, dtype=float) * Rr
    )
    work_pR = -FOUR_PI * np.square(R) * np.asarray(total.pr, dtype=float) * Rdot
    fixed_coordinate_residual = mass_dot - (flux_T + work_pR)

    coordinate_speed_fixed_R = np.full_like(R, np.nan)
    valid = np.abs(Rr) > np.finfo(float).tiny
    coordinate_speed_fixed_R[valid] = -Rdot[valid] / Rr[valid]
    fixed_R_rate = mass_dot + mass_r * coordinate_speed_fixed_R
    fixed_R_residual = fixed_R_rate - (flux_T + work_pR)

    def sampled(name, values):
        value, bracket = _sample_at_areal_radius(R, values, TARGET_AREAL_RADIUS)
        return value, bracket

    fields = {
        "coordinate_r": r,
        "R": R,
        "Rr": Rr,
        "Rdot": Rdot,
        "mass": mass,
        "mass_dot_fixed_coordinate": mass_dot,
        "mass_radial_derivative": mass_r,
        "flux_T": flux_T,
        "work_pR": work_pR,
        "fixed_coordinate_residual": fixed_coordinate_residual,
        "coordinate_speed_fixed_R": coordinate_speed_fixed_R,
        "mass_rate_fixed_R": fixed_R_rate,
        "fixed_R_residual": fixed_R_residual,
        "total_rho_consumer_units": np.asarray(total.rho, dtype=float),
        "total_pr_consumer_units": np.asarray(total.pr, dtype=float),
        "total_j_consumer_units": np.asarray(total.j, dtype=float),
    }
    for name, values in fields.items():
        if not np.all(np.isfinite(np.asarray(values, dtype=float))):
            raise AssertionError(f"non-finite quasi-local field {name}")

    sampled_values = {}
    bracket = None
    for name in (
        "coordinate_r", "Rdot", "Rr", "mass", "mass_dot_fixed_coordinate",
        "mass_radial_derivative", "flux_T", "work_pR",
        "fixed_coordinate_residual", "coordinate_speed_fixed_R",
        "mass_rate_fixed_R", "fixed_R_residual",
    ):
        sampled_values[name], current_bracket = sampled(name, fields[name])
        if bracket is None:
            bracket = current_bracket

    return {
        "target_areal_radius": TARGET_AREAL_RADIUS,
        "sample_bracket": bracket,
        "sampled_terms": sampled_values,
        "interpretation": (
            "signed observer only; consumer-unit stress, physical EM normalization "
            "and the full donor+field+gravity conservation claim remain open"
        ),
    }


def _em_stage_powers(grid, stage, geometry_terms, mass_flux):
    electromagnetic = stage.electromagnetic
    if electromagnetic is None or stage.matter.charged_dust is None:
        raise AssertionError("feature-on stage lost its EM/carrier state")
    geometry = stage.geometry
    q_faces = np.asarray(electromagnetic.enclosed_charge_faces, dtype=float)
    q_cell = 0.5 * (q_faces[:-1] + q_faces[1:])
    weight = _weight(grid, geometry)
    current_faces = electromagnetic.charge_flux_faces(mass_flux)
    q_cell_dot = -0.5 * (current_faces[:-1] + current_faces[1:])
    charge_power = grid.dr * float(np.sum(weight * q_cell * q_cell_dot))

    explicit = geometry_terms["explicit"]
    log_weight_dot = (
        np.asarray(explicit["X"], dtype=float) / np.asarray(geometry.X, dtype=float)
        + 0.5 * np.asarray(explicit["a"], dtype=float) / np.asarray(geometry.a, dtype=float)
        - np.asarray(explicit["b"], dtype=float) / np.asarray(geometry.b, dtype=float)
    )
    metric_power = 0.5 * grid.dr * float(np.sum(weight * q_cell**2 * log_weight_dot))

    metric = metric_slice_from_q(grid, geometry)
    _, _, source = electromagnetic.lorentz_source(
        grid, geometry, metric, stage.matter.charged_dust
    )
    matter_lorentz_power = FOUR_PI * grid.dr * float(
        np.sum(np.asarray(source["source_energy"], dtype=float))
    )
    outer_current_power = (
        -0.5 * grid.dr * weight[-1] * q_cell[-1] * current_faces[-1]
    )

    values = (charge_power, metric_power, matter_lorentz_power, outer_current_power)
    if not all(math.isfinite(value) for value in values):
        raise AssertionError(f"non-finite EM stage powers: {values}")
    return {
        "charge_power": charge_power,
        "metric_power": metric_power,
        "matter_lorentz_power": matter_lorentz_power,
        "outer_current_power": outer_current_power,
        "outer_charge_current": float(current_faces[-1]),
    }


def _rhs_term_ledger(stage_index, stage_record, component_rows, grid):
    stage = stage_record["state"]
    dt_dr = 1.0 / float(grid.dr)
    species_rows = []
    for name in ("dark_matter", "baryons", "radiation", "charged_dust"):
        key = (stage_index, name)
        if key not in component_rows:
            if name == "charged_dust" and stage.matter.charged_dust is None:
                continue
            raise AssertionError(f"missing exact-stage component record for {key}")
        comp = component_rows[key]
        face_flux = comp["face_flux"]
        source_e = comp["source_e"]
        source_s = comp["source_s"]
        q_e = comp["q_e"]
        q_s = comp["q_s"]
        actual = stage_record["matter_rhs"][name]

        transport = -dt_dr * np.diff(face_flux, axis=0)
        force_e = np.zeros(grid.n, dtype=float)
        force_s = np.zeros(grid.n, dtype=float)
        if name == "charged_dust":
            _, _, force = stage.electromagnetic.lorentz_source(
                grid, stage.geometry, metric_slice_from_q(grid, stage.geometry),
                stage.matter.charged_dust,
            )
            force_e = np.asarray(force["source_energy"], dtype=float)
            force_s = np.asarray(force["source_momentum_r"], dtype=float)

        reconstructed = {
            "rest": transport[:, 0],
            "energy_t": transport[:, 1] + source_e + q_e + force_e,
            "momentum_r": transport[:, 2] + source_s + q_s + force_s,
        }
        actual_arrays = {
            "rest": np.asarray(actual.rest, dtype=float),
            "energy_t": np.asarray(actual.energy_t, dtype=float),
            "momentum_r": np.asarray(actual.momentum_r, dtype=float),
        }
        term_arrays = {
            "rest": (transport[:, 0],),
            "energy_t": (transport[:, 1], source_e, q_e, force_e),
            "momentum_r": (transport[:, 2], source_s, q_s, force_s),
        }
        residuals = {}
        for component, actual_values in actual_arrays.items():
            residual = actual_values - reconstructed[component]
            operands = [actual_values, reconstructed[component], *term_arrays[component]]
            scale = np.maximum.reduce([np.abs(x) for x in operands])
            roundoff_bound = (
                128.0 * np.finfo(float).eps
                * np.maximum(scale, np.finfo(float).tiny)
            )
            if np.any(np.abs(residual) > roundoff_bound):
                raise AssertionError(
                    f"{name} {component} RHS decomposition exceeds scale-aware "
                    f"roundoff at stage {stage_index}: max={_max_abs(residual):.6e}, "
                    f"bound={_max_abs(roundoff_bound):.6e}"
                )
            residuals[component] = _max_abs(residual)

        def full_sphere_rate(values):
            return FOUR_PI * grid.dr * float(np.sum(np.asarray(values, dtype=float)))

        species_rows.append({
            "species": name,
            "inner_face_flux_rest_energy_momentum": [float(x) for x in face_flux[0]],
            "outer_face_flux_rest_energy_momentum": [float(x) for x in face_flux[-1]],
            "outer_energy_flux_full_sphere": FOUR_PI * float(face_flux[-1, 1]),
            "transport_rate_full_sphere_rest_energy_momentum": [
                full_sphere_rate(transport[:, col]) for col in range(3)
            ],
            "geometric_energy_source_power_full_sphere": full_sphere_rate(source_e),
            "dm_scalar_exchange_energy_power_full_sphere": full_sphere_rate(q_e),
            "em_lorentz_energy_power_full_sphere": full_sphere_rate(force_e),
            "geometric_momentum_source_full_sphere": full_sphere_rate(source_s),
            "dm_scalar_exchange_momentum_source_full_sphere": full_sphere_rate(q_s),
            "em_lorentz_momentum_source_full_sphere": full_sphere_rate(force_s),
            "actual_rhs_rate_full_sphere_rest_energy_momentum": [
                full_sphere_rate(actual_arrays["rest"]),
                full_sphere_rate(actual_arrays["energy_t"]),
                full_sphere_rate(actual_arrays["momentum_r"]),
            ],
            "rhs_reconstruction_residual_max_abs": residuals,
        })

    return {
        "stage": "start" if stage_index == 0 else "primary_predictor",
        "scalar_rhs_absmax": _max_abs([
            getattr(stage_record["scalar_rhs"], name)
            for name in ("S", "PS", "D", "PD", "phi", "Pi")
        ]),
        "species": species_rows,
    }


def _run_probe(resolution, eta, dt):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution,
        r_max=R_MAX,
        include_radiation=True,
        radial_maxwell_charge_to_mass=eta,
        charged_dust_density=np.full(resolution, DENSITY, dtype=float),
    )

    # Match the explicit-carrier fixture used by the accepted local-work audit.
    r = np.asarray(state.grid.centers, dtype=float)
    velocity = 0.08 * (r / R_MAX) * (1.0 - r / R_MAX)
    carrier = initialize_dust(
        _metric_arrays(metric_slice_from_q(state.grid, state.geometry)),
        np.full(resolution, DENSITY, dtype=float),
        velocity,
    )
    state.matter = V55MatterState(
        state.matter.dark_matter,
        state.matter.baryons,
        state.matter.radiation,
        charged_dust=carrier,
    )
    state.electromagnetic = RadialMaxwellState.from_charged_dust(
        state.grid, carrier, eta
    )

    kernel.diagnostics(state)
    active = {"index": None, "state": None}
    rhs_records = []
    component_rows = {}
    mass_flux_records = []
    geometry_records = []

    original_components = matter_rhs_module._species_rhs_components
    original_rhs = kernel._rhs
    original_flux = kernel._charged_dust_mass_flux_faces
    original_geometry = kernel._geometry_stage_terms

    def record_components(metric, derivatives, species_state, species, *args, **kwargs):
        result = original_components(metric, derivatives, species_state, species, *args, **kwargs)
        if active["index"] is not None:
            stage = active["state"]
            label = next((
                name for name in ("dark_matter", "baryons", "radiation", "charged_dust")
                if getattr(stage.matter, name, None) is species_state
            ), None)
            if label is None:
                raise AssertionError("an active species RHS was not associated with its live stage")
            component_rows[(active["index"], label)] = {
                "face_flux": np.asarray(result[0], dtype=float).copy(),
                "source_e": np.asarray(result[1], dtype=float).copy(),
                "source_s": np.asarray(result[2], dtype=float).copy(),
                "q_e": np.asarray(result[3], dtype=float).copy(),
                "q_s": np.asarray(result[4], dtype=float).copy(),
            }
        return result

    def record_rhs(stage):
        index = len(rhs_records)
        active["index"] = index
        active["state"] = stage
        try:
            scalar_rhs, matter_rhs, metric_derivatives = original_rhs(stage)
        finally:
            active["index"] = None
            active["state"] = None
        rhs_records.append({
            "state": stage,
            "scalar_rhs": scalar_rhs,
            "matter_rhs": matter_rhs,
            "metric_derivatives": metric_derivatives,
        })
        return scalar_rhs, matter_rhs, metric_derivatives

    def record_flux(stage):
        flux = np.asarray(original_flux(stage), dtype=float).copy()
        mass_flux_records.append({"state": stage, "flux": flux})
        return flux

    def record_geometry(grid, geometry, scalars, matter, electromagnetic):
        terms = original_geometry(grid, geometry, scalars, matter, electromagnetic)
        geometry_records.append({
            "grid": grid,
            "geometry": geometry.copy(),
            "scalars": scalars,
            "matter": matter,
            "electromagnetic": None if electromagnetic is None else electromagnetic.copy(),
            "terms": terms,
        })
        return terms

    # These wrappers inspect the active production call sites; they do not
    # modify the numerical functions or their returned arrays.
    matter_rhs_module._species_rhs_components = record_components
    kernel._rhs = record_rhs
    kernel._charged_dust_mass_flux_faces = record_flux
    kernel._geometry_stage_terms = record_geometry
    try:
        accepted = kernel.step(state, dt)
    finally:
        matter_rhs_module._species_rhs_components = original_components
        kernel._rhs = original_rhs
        kernel._charged_dust_mass_flux_faces = original_flux
        kernel._geometry_stage_terms = original_geometry

    assert len(rhs_records) == 2
    assert len(mass_flux_records) == 2
    assert len(geometry_records) == 3
    assert all(np.all(np.isfinite(row["flux"])) for row in mass_flux_records)

    matter_ledger = [
        _rhs_term_ledger(index, rhs_records[index], component_rows, state.grid)
        for index in range(2)
    ]

    # Match PIRK's real stages: geometry call 0 is the accepted-start stage;
    # geometry call 2 is the primary predictor used for the trapezoidal update.
    em_stage_powers = []
    for stage_index, geometry_index in ((0, 0), (1, 2)):
        em_stage_powers.append(_em_stage_powers(
            state.grid,
            rhs_records[stage_index]["state"],
            geometry_records[geometry_index]["terms"],
            mass_flux_records[stage_index]["flux"],
        ))

    u_start = _field_energy(
        state.grid, rhs_records[0]["state"].geometry,
        rhs_records[0]["state"].electromagnetic,
    )
    u_end = _field_energy(state.grid, accepted.geometry, accepted.electromagnetic)
    field_energy_rate = (u_end - u_start) / dt
    average = lambda key: 0.5 * (
        em_stage_powers[0][key] + em_stage_powers[1][key]
    )
    em_residual = (
        field_energy_rate
        - average("metric_power")
        + average("matter_lorentz_power")
        - average("outer_current_power")
    )
    em_scale = max(
        abs(field_energy_rate),
        *(abs(average(key)) for key in (
            "charge_power", "metric_power", "matter_lorentz_power", "outer_current_power"
        )),
        np.finfo(float).tiny,
    )

    ms_ledger = []
    for stage_index, geometry_index in ((0, 0), (1, 2)):
        record = geometry_records[geometry_index]
        ms_ledger.append({
            "stage": "start" if stage_index == 0 else "primary_predictor",
            **_ms_rate_ledger(
                kernel, state.grid, record["geometry"], record["scalars"],
                record["matter"], record["electromagnetic"], record["terms"],
            ),
        })

    accepted_diag = accepted.history[-1]
    report_row = {
        "q_over_m": eta,
        "resolution": resolution,
        "r_max": R_MAX,
        "dt": dt,
        "t_start": float(state.t),
        "t_accepted": float(accepted.t),
        "accepted_steps": 1,
        "source_anchor": "5c51bf509504365df394e6cca96bf1d117f4614d",
        "status": "research_only_not_admitted",
        "matter_ledger": matter_ledger,
        "em_local_exchange": {
            "field_energy_start": u_start,
            "field_energy_accepted": u_end,
            "field_energy_rate": field_energy_rate,
            "start_stage_powers": em_stage_powers[0],
            "predictor_stage_powers": em_stage_powers[1],
            "trapezoid_averaged_powers": {
                key: average(key) for key in (
                    "charge_power", "metric_power", "matter_lorentz_power", "outer_current_power"
                )
            },
            "signed_residual_rate": em_residual,
            "relative_residual": abs(em_residual) / em_scale,
        },
        "misner_sharp_ledger": ms_ledger,
        "accepted_constraints": {
            name: float(accepted_diag[name])
            for name in (
                "hamiltonian_full_absmax", "momentum_full_absmax",
                "connection_full_absmax", "determinant_constraint_max",
            )
            if name in accepted_diag and accepted_diag[name] is not None
        },
        "scope_limits": [
            "one-step stage ledger only; not long-horizon stability or convergence",
            "zero-shift CMC branch only",
            "consumer stress arrays are reported in the current kernel's units; physical field-to-engine normalization remains open",
            "this is a component ledger, not a claim of one conserved global energy in a dynamical spacetime",
            "the current-main scalar-DM rest-current source contract remains a separate conditional gate (see PR #143)",
            "full donor+field+gravity closure, matter/geometry boundary bookkeeping, restart equivalence and scientific admission remain open",
        ],
    }
    serialized = json.dumps(report_row, sort_keys=True, allow_nan=False)
    assert json.loads(serialized)["status"] == "research_only_not_admitted"
    if eta == 0.0:
        em = report_row["em_local_exchange"]
        assert em["field_energy_start"] == 0.0
        assert em["field_energy_accepted"] == 0.0
        assert _max_abs([row["matter_lorentz_power"] for row in em_stage_powers]) == 0.0
    else:
        em = report_row["em_local_exchange"]
        assert em["field_energy_start"] > 0.0
        assert em["field_energy_accepted"] > 0.0
        assert max(
            abs(row["matter_lorentz_power"]) for row in em_stage_powers
        ) > 0.0
    assert math.isfinite(report_row["em_local_exchange"]["signed_residual_rate"])
    assert all(
        math.isfinite(value)
        for ms_row in ms_ledger
        for value in ms_row["sampled_terms"].values()
    )
    return report_row


def test_report_joined_opt_in_em_matter_and_misner_sharp_ledger(tmp_path, capsys):
    rows = []
    for resolution in RESOLUTIONS:
        dt = 0.002 * 24.0 / resolution
        for eta in ETA_ARMS:
            rows.append(_run_probe(resolution, eta, dt))

    report = {
        "classification": (
            "bounded one-step joined matter/EM/Misner-Sharp component ledger; "
            "diagnostic only, not scientific admission"
        ),
        "source_anchor": "5c51bf509504365df394e6cca96bf1d117f4614d",
        "settings": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "dt_rule": "0.002 * 24 / N",
            "explicit_separate_carrier_density": DENSITY,
            "charge_to_mass_arms": list(ETA_ARMS),
            "target_areal_radius": TARGET_AREAL_RADIUS,
        },
        "rows": rows,
    }
    payload = json.dumps(report, sort_keys=True, allow_nan=False)
    output_dir = Path(os.environ.get("EM_LEDGER_REPORT_DIR", str(tmp_path)))
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "joined_em_matter_misner_sharp_ledger.json"
    report_path.write_text(payload + "\n", encoding="utf-8")

    with capsys.disabled():
        print("EM_COMPLETE_LEDGER_REPORT_JSON=" + payload, flush=True)
    assert report_path.exists()
    assert json.loads(report_path.read_text(encoding="utf-8")) == report
