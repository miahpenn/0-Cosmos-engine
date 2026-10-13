#!/usr/bin/env python3
"""Bounded, diagnostic-only signed momentum-constraint tangent audit.

The audit evaluates the current production semi-discrete RHS on newly
initialized frozen slices. It does not evolve a campaign. It re-solves the CMC
lapse for both directional perturbations and splits the momentum residual rate
into its geometric constraint term and six independently projected current
sectors. A single small ordinary step and half-step are used only as a local
tangent consistency check.

A passing arithmetic closure is not physical/scientific admission.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.production_kernel import (
    BETA_DM,
    LAMBDA_M,
    ProductionState,
    V55ProductionKernel,
)
from engine.scalar_system import KAPPA, ScalarFields
from engine.matter_system import ConservedSpecies, Species, project_species
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter


SCHEMA = "current_source_signed_momentum_constraint_tangent_v1"
SOURCE_BASELINE = "29d828b4d72b017c12a20b6dc90549f2511f025c"
EXPECTED_VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_BLOBS = {
    "engine/production_kernel.py": "2954e7a79656849b5f088053ea0f4ac97d4bb443",
    "engine/matter_system.py": "4f99e3bb61736473d8b9220145aec9e4167a9cff",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/v55_matter.py": "0f71164494e534ca4226c56c10058845d9a14e1d",
    "engine/v55_pirk_adapter.py": "f671b838d9aefd8926702d9c2235209f52cd6815",
    "engine/cmc_gauge.py": "b1e6b3e95ff37702509e1b38ec6401d10e2482c3",
    "engine/v55_initial.py": "4cfbfe4be6ac3ce6ed492f4e435186997d94116a",
    "engine/stress_energy.py": "c31c647f74fab6893254b3dea1366eada15b3110",
    "engine/valencia.py": "9eba0855936e1974a15314eacdc2f6c378072596",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
}
SCALAR_NAMES = ("S", "PS", "D", "PD", "phi", "Pi")
GEOMETRY_RATE_NAMES = ("a", "b", "X", "Aa", "K", "Lambda")
SECTOR_NAMES = ("S", "D", "COSMOS_phi", "dark_matter", "baryons", "radiation")
GEOMETRIC_M_NAMES = (
    "D_r_Aa",
    "minus_2third_D_r_K",
    "6_Aa_chi_r",
    "anisotropy_2_over_r",
    "anisotropy_b_gradient",
)
MACHINE_EPS = float(np.finfo(np.float64).eps)


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def provenance() -> dict[str, Any]:
    actual_blobs = {
        path: _git("rev-parse", f"HEAD:{path}") for path in EXPECTED_BLOBS
    }
    wrong = {
        path: {"expected": expected, "actual": actual_blobs[path]}
        for path, expected in EXPECTED_BLOBS.items()
        if actual_blobs[path] != expected
    }
    if wrong:
        raise RuntimeError(
            "production source blob pin mismatch; refusing non-exact tangent "
            f"measurement: {json.dumps(wrong, sort_keys=True)}"
        )
    vendor = _git("-C", "vendor/bb-palatini-unified-r0", "rev-parse", "HEAD")
    if vendor != EXPECTED_VENDOR:
        raise RuntimeError(
            f"vendor pin mismatch: expected {EXPECTED_VENDOR}, got {vendor}"
        )
    return {
        "repository_commit": _git("rev-parse", "HEAD"),
        "production_source_baseline_anchor": SOURCE_BASELINE,
        "production_source_blobs": actual_blobs,
        "vendor_commit": vendor,
    }


def state_digest(state: ProductionState) -> str:
    h = hashlib.sha256()
    def add(name: str, values: Any) -> None:
        arr = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
        h.update(name.encode("utf-8"))
        h.update(str(arr.shape).encode("ascii"))
        h.update(arr.tobytes(order="C"))
    add("grid.centers", state.grid.centers)
    for name in ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B"):
        add(f"geometry.{name}", getattr(state.geometry, name))
    for name in SCALAR_NAMES:
        add(f"scalars.{name}", getattr(state.scalars, name))
    for species_name in ("dark_matter", "baryons", "radiation"):
        species = getattr(state.matter, species_name)
        for name in ("rest", "energy_t", "momentum_r"):
            add(f"matter.{species_name}.{name}", getattr(species, name))
    return h.hexdigest()


def finite(name: str, values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not bool(np.all(np.isfinite(arr))):
        raise FloatingPointError(f"non-finite values in {name}")
    return arr


def current_pieces(state: ProductionState) -> dict[str, np.ndarray]:
    """Rebuild each signed J contribution using production projection operators."""
    grid, geom, scalars, matter = (
        state.grid, state.geometry, state.scalars, state.matter
    )
    d1 = lambda f: grid.cell_derivative_fourth(f, parity=1)
    pieces = {
        "S": -np.asarray(scalars.PS) * d1(scalars.S),
        "D": -np.asarray(scalars.PD) * d1(scalars.D),
        "COSMOS_phi": (
            -np.asarray(scalars.Pi) * d1(scalars.phi) / KAPPA
        ),
    }
    metric = metric_slice_from_q(grid, geom)
    for key, state_name, species in (
        ("dark_matter", "dark_matter", Species.DARK_MATTER),
        ("baryons", "baryons", Species.BARYON),
        ("radiation", "radiation", Species.RADIATION),
    ):
        pieces[key] = np.asarray(
            project_species(metric, getattr(matter, state_name), species)["j"],
            dtype=np.float64,
        )
    return {key: finite(f"current {key}", value) for key, value in pieces.items()}


def geometric_momentum_pieces(state: ProductionState) -> dict[str, np.ndarray]:
    """Return signed terms from the exact pinned vendor momentum-constraint formula."""
    grid, geom = state.grid, state.geometry
    d1 = lambda f: grid.cell_derivative_fourth(f, parity=1)
    r = np.asarray(grid.centers, dtype=np.float64)
    Aa = np.asarray(geom.Aa, dtype=np.float64)
    Ab = -0.5 * Aa
    chip = -0.5 * d1(np.asarray(geom.X, dtype=np.float64)) / np.asarray(geom.X)
    bp = d1(np.asarray(geom.b, dtype=np.float64))
    return {
        "D_r_Aa": finite("M geometry D_r Aa", d1(Aa)),
        "minus_2third_D_r_K": finite(
            "M geometry -2/3 D_r K",
            -(2.0 / 3.0) * d1(np.asarray(geom.K, dtype=np.float64)),
        ),
        "6_Aa_chi_r": finite("M geometry 6 Aa chi_r", 6.0 * Aa * chip),
        "anisotropy_2_over_r": finite(
            "M geometry anisotropy 2/r", (Aa - Ab) * (2.0 / r),
        ),
        "anisotropy_b_gradient": finite(
            "M geometry anisotropy b gradient", (Aa - Ab) * bp / np.asarray(geom.b),
        ),
    }


def observe_constraint(state: ProductionState) -> dict[str, Any]:
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(state.grid, state.geometry)
    if "momentum" not in raw:
        raise KeyError("pinned vendor constraints output lacks 'momentum'")
    m_geom = finite("geometric momentum constraint", raw["momentum"])
    geom_pieces = geometric_momentum_pieces(state)
    geom_sum = np.zeros_like(m_geom)
    geom_operand_scale = np.abs(m_geom)
    for value in geom_pieces.values():
        geom_sum = geom_sum + value
        geom_operand_scale = geom_operand_scale + np.abs(value)
    geom_closure = m_geom - geom_sum
    geom_closure_bound = 64.0 * MACHINE_EPS * geom_operand_scale
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    j_total = finite("total projected current", total.j)
    pieces = current_pieces(state)
    j_sum = np.zeros_like(j_total)
    j_scale = np.abs(j_total)
    for value in pieces.values():
        j_sum = j_sum + value
        j_scale = j_scale + np.abs(value)
    j_residual = j_total - j_sum
    j_bound = 64.0 * MACHINE_EPS * j_scale
    m_total = m_geom - 8.0 * math.pi * j_total
    return {
        "M_geom": m_geom,
        "M_geom_pieces": geom_pieces,
        "M_geom_closure_abs": np.abs(geom_closure),
        "M_geom_closure_bound": geom_closure_bound,
        "M_geom_closure_ok": bool(np.all(np.abs(geom_closure) <= geom_closure_bound)),
        "J_total": j_total,
        "J_pieces": pieces,
        "M_total": finite("total momentum constraint", m_total),
        "J_closure_abs": np.abs(j_residual),
        "J_closure_bound": j_bound,
        "J_closure_ok": bool(np.all(np.abs(j_residual) <= j_bound)),
    }


def rhs_direction(
    kernel: V55ProductionKernel, state: ProductionState
) -> dict[str, Any]:
    srhs, mrhs, _ = kernel._rhs(state)
    terms = adapter.geometry_stage_terms(
        state.grid, state.geometry, state.scalars, state.matter,
        lambda_m=LAMBDA_M,
    )
    geom_rate: dict[str, np.ndarray] = {
        name: finite(f"explicit geometry rate {name}", terms["explicit"][name])
        for name in ("a", "b", "X")
    }
    for name in ("Aa", "K"):
        geom_rate[name] = finite(
            f"primary geometry rate {name}",
            np.asarray(terms["primary_l2"][name])
            + np.asarray(terms["primary_l3"][name]),
        )
    geom_rate["Lambda"] = finite(
        "Lambda geometry rate",
        np.asarray(terms["lambda_l2"]) + np.asarray(terms["lambda_l3"]),
    )
    scalar_rate = {
        name: finite(f"scalar rate {name}", getattr(srhs, name))
        for name in SCALAR_NAMES
    }
    matter_rate = {}
    for species_name in ("dark_matter", "baryons", "radiation"):
        r = mrhs[species_name]
        matter_rate[species_name] = {
            name: finite(
                f"matter rate {species_name}.{name}", getattr(r, name)
            )
            for name in ("rest", "energy_t", "momentum_r")
        }
    return {
        "geometry": geom_rate,
        "scalars": scalar_rate,
        "matter": matter_rate,
        "geometry_stage_terms": terms,
    }


def shifted_state(
    kernel: V55ProductionKernel,
    state: ProductionState,
    rates: dict[str, Any],
    delta: float,
) -> ProductionState:
    """Form a directional state +/- epsilon*RHS and re-solve the CMC lapse."""
    geom = state.geometry.copy()
    for name in GEOMETRY_RATE_NAMES:
        updated = (
            np.asarray(getattr(state.geometry, name), dtype=np.float64)
            + delta * rates["geometry"][name]
        )
        setattr(geom, name, finite(f"perturbed geometry {name}", updated))
    scalars = ScalarFields(
        *(
            finite(
                f"perturbed scalar {name}",
                getattr(state.scalars, name) + delta * rates["scalars"][name],
            )
            for name in SCALAR_NAMES
        )
    )
    species_out = {}
    for species_name in ("dark_matter", "baryons", "radiation"):
        base = getattr(state.matter, species_name)
        rate = rates["matter"][species_name]
        species_out[species_name] = ConservedSpecies(
            *(
                finite(
                    f"perturbed matter {species_name}.{name}",
                    getattr(base, name) + delta * rate[name],
                )
                for name in ("rest", "energy_t", "momentum_r")
            )
        )
    matter = V55MatterState(
        species_out["dark_matter"], species_out["baryons"],
        species_out["radiation"],
    )
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    geom.alpha = kernel._solve_lapse(
        state.grid, geom, scalars, matter
    )[0]
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    return ProductionState(
        grid=state.grid,
        geometry=geom,
        scalars=scalars,
        matter=matter,
        t=state.t,
        tau=state.tau,
        e_folds=state.e_folds,
    )


def _regions(r: np.ndarray) -> dict[str, np.ndarray]:
    i = np.arange(r.size)
    return {
        "cell_0": i == 0,
        "first_two_cells": i < 2,
        "first_three_cells": i < 3,
        "physical_annulus_2_le_r_lt_5": (r >= 2.0) & (r < 5.0),
        "physical_annulus_5_le_r_lt_10": (r >= 5.0) & (r < 10.0),
        "domain_r_ge_5": r >= 5.0,
    }


def summarize(values: Any, r: np.ndarray) -> dict[str, Any]:
    arr = finite("summary values", values)
    i = int(np.argmax(np.abs(arr)))
    return {
        "max_abs": float(abs(arr[i])),
        "signed_at_max": float(arr[i]),
        "index": i,
        "radius": float(r[i]),
        "rms": float(np.sqrt(np.mean(arr * arr))),
    }


def regional_summaries(values: Any, r: np.ndarray) -> dict[str, Any]:
    arr = finite("regional summary values", values)
    out = {}
    for label, mask in _regions(r).items():
        if bool(np.any(mask)):
            out[label] = summarize(arr[mask], r[mask])
        else:
            out[label] = None
    return out


def fixed_radius_samples(values: Any, r: np.ndarray) -> dict[str, float]:
    arr = finite("fixed-radius values", values)
    return {
        str(radius): float(np.interp(radius, r, arr))
        for radius in (2.0, 3.0, 5.0, 8.0)
        if float(r[0]) <= radius <= float(r[-1])
    }


def rhs_rate_witness(values: Any, r: np.ndarray) -> dict[str, Any]:
    """Compact cell-centred witnesses for an exact RHS rate profile."""
    arr = finite("RHS rate witness", values)
    return {
        "cell_0": float(arr[0]),
        "cell_1": float(arr[1]) if arr.size > 1 else None,
        "cell_2": float(arr[2]) if arr.size > 2 else None,
        "fixed_radius_samples": fixed_radius_samples(arr, r),
        "first_two_cells_max_abs": float(np.max(np.abs(arr[:2]))),
        "first_three_cells_max_abs": float(np.max(np.abs(arr[:3]))),
    }


def tangent_for_epsilon(
    kernel: V55ProductionKernel,
    state: ProductionState,
    rates: dict[str, Any],
    epsilon: float,
) -> dict[str, Any]:
    plus_state = shifted_state(kernel, state, rates, epsilon)
    minus_state = shifted_state(kernel, state, rates, -epsilon)
    plus = observe_constraint(plus_state)
    minus = observe_constraint(minus_state)
    denom = 2.0 * epsilon
    m_geom_dot = (plus["M_geom"] - minus["M_geom"]) / denom
    geom_term_dots = {
        key: (plus["M_geom_pieces"][key] - minus["M_geom_pieces"][key]) / denom
        for key in GEOMETRIC_M_NAMES
    }
    geom_term_sum = np.zeros_like(m_geom_dot)
    geom_term_sum_scale = np.zeros_like(m_geom_dot)
    geom_operand_scale = np.zeros_like(m_geom_dot)
    for key in GEOMETRIC_M_NAMES:
        geom_term_sum = geom_term_sum + geom_term_dots[key]
        geom_term_sum_scale = geom_term_sum_scale + np.abs(geom_term_dots[key])
        geom_operand_scale = geom_operand_scale + (
            np.abs(plus["M_geom_pieces"][key]) + np.abs(minus["M_geom_pieces"][key])
        )
    geom_tangent_closure = m_geom_dot - geom_term_sum
    geom_tangent_bound = (128.0 * MACHINE_EPS / denom) * geom_operand_scale
    geom_tangent_closure_ok = bool(
        np.all(np.abs(geom_tangent_closure) <= geom_tangent_bound)
    )
    j_dot = {
        key: (plus["J_pieces"][key] - minus["J_pieces"][key]) / denom
        for key in SECTOR_NAMES
    }
    sector_terms = {
        key: -8.0 * math.pi * value for key, value in j_dot.items()
    }
    m_dot = (plus["M_total"] - minus["M_total"]) / denom
    reconstructed = m_geom_dot.copy()
    for value in sector_terms.values():
        reconstructed = reconstructed + value
    closure = m_dot - reconstructed

    # Bound accounts for subtraction operands before differencing, not only
    # the usually much smaller derivative result.
    operand_scale = (
        np.abs(plus["M_total"]) + np.abs(minus["M_total"])
        + np.abs(plus["M_geom"]) + np.abs(minus["M_geom"])
    )
    for key in SECTOR_NAMES:
        operand_scale = operand_scale + 8.0 * math.pi * (
            np.abs(plus["J_pieces"][key]) + np.abs(minus["J_pieces"][key])
        )
    bound = (128.0 * MACHINE_EPS / denom) * operand_scale
    finite("momentum tangent", m_dot)
    finite("momentum tangent closure", closure)
    closure_ok = bool(np.all(np.abs(closure) <= bound))
    ratio = np.zeros_like(closure)
    positive = bound > 0.0
    ratio[positive] = np.abs(closure[positive]) / bound[positive]
    if bool(np.any((~positive) & (np.abs(closure) > 0.0))):
        max_ratio: float | None = None
    else:
        max_ratio = float(np.max(ratio))
    return {
        "epsilon": float(epsilon),
        "Mdot": m_dot,
        "Mdot_geom": m_geom_dot,
        "Mdot_geom_terms": geom_term_dots,
        "Mdot_geom_term_reconstructed": geom_term_sum,
        "Mdot_geom_tangent_closure_abs": np.abs(geom_tangent_closure),
        "Mdot_geom_tangent_closure_bound": geom_tangent_bound,
        "Mdot_geom_tangent_closure_ok": geom_tangent_closure_ok,
        "Mdot_geom_tangent_closure_max_abs": float(np.max(np.abs(geom_tangent_closure))),
        "Mdot_geom_tangent_closure_max_bound": float(np.max(geom_tangent_bound)),
        "Mdot_sector": sector_terms,
        "Mdot_reconstructed": reconstructed,
        "closure_abs": np.abs(closure),
        "closure_bound": bound,
        "closure_ok": closure_ok,
        "closure_max_abs": float(np.max(np.abs(closure))),
        "closure_max_bound": float(np.max(bound)),
        "closure_max_scaled_ratio": max_ratio,
        "plus_minus_current_projection_closure_ok": bool(
            plus["J_closure_ok"] and minus["J_closure_ok"]
        ),
        "plus_minus_current_projection_closure_max_abs": float(
            max(np.max(plus["J_closure_abs"]), np.max(minus["J_closure_abs"]))
        ),
    }


def compact_tangent_result(item: dict[str, Any], r: np.ndarray) -> dict[str, Any]:
    return {
        "epsilon": item["epsilon"],
        "Mdot": regional_summaries(item["Mdot"], r),
        "Mdot_fixed_radius_samples": fixed_radius_samples(item["Mdot"], r),
        "Mdot_geom": regional_summaries(item["Mdot_geom"], r),
        "Mdot_sector_contributions": {
            name: regional_summaries(value, r)
            for name, value in item["Mdot_sector"].items()
        },
        "Mdot_sector_fixed_radius_samples": {
            name: fixed_radius_samples(value, r)
            for name, value in item["Mdot_sector"].items()
        },
        "geometric_term_contributions": {
            name: regional_summaries(value, r)
            for name, value in item["Mdot_geom_terms"].items()
        },
        "geometric_term_fixed_radius_samples": {
            name: fixed_radius_samples(value, r)
            for name, value in item["Mdot_geom_terms"].items()
        },
        "geometric_term_tangent_closure": {
            "ok": item["Mdot_geom_tangent_closure_ok"],
            "max_abs": item["Mdot_geom_tangent_closure_max_abs"],
            "max_bound": item["Mdot_geom_tangent_closure_max_bound"],
        },
        "closure": {
            "ok": item["closure_ok"],
            "max_abs": item["closure_max_abs"],
            "max_bound": item["closure_max_bound"],
            "max_scaled_ratio": item["closure_max_scaled_ratio"],
        },
        "current_projection_closure": {
            "ok": item["plus_minus_current_projection_closure_ok"],
            "max_abs": item["plus_minus_current_projection_closure_max_abs"],
        },
    }


def one_step_check(
    kernel: V55ProductionKernel,
    state: ProductionState,
    m_initial: np.ndarray,
    tangent: np.ndarray,
    dt: float,
) -> dict[str, Any]:
    try:
        out = kernel.step(state, dt)
        m_after = observe_constraint(out)["M_total"]
        rate = (m_after - m_initial) / dt
        error = rate - tangent
        return {
            "completed": True,
            "dt": float(dt),
            "t_final": float(out.t),
            "M_change_over_dt": regional_summaries(rate, state.grid.centers),
            "difference_from_semi_discrete_tangent": regional_summaries(
                error, state.grid.centers
            ),
            "difference_fixed_radius_samples": fixed_radius_samples(
                error, state.grid.centers
            ),
        }
    except Exception as exc:
        return {
            "completed": False,
            "dt": float(dt),
            "failure_type": type(exc).__name__,
            "failure_message": str(exc),
        }


def audit_resolution(
    kernel: V55ProductionKernel,
    resolution: int,
    r_max: float,
    epsilons: tuple[float, float],
    dt: float,
) -> dict[str, Any]:
    state = kernel.initialize(
        resolution=resolution, r_max=r_max, include_radiation=True
    )
    original_digest = state_digest(state)
    rates = rhs_direction(kernel, state)
    w0 = observe_constraint(state)
    r = np.asarray(state.grid.centers, dtype=np.float64)
    tangent_results = [
        tangent_for_epsilon(kernel, state, rates, eps)
        for eps in epsilons
    ]
    e0, e1 = tangent_results
    eps_difference = e0["Mdot"] - e1["Mdot"]
    step_full = one_step_check(kernel, state, w0["M_total"], e0["Mdot"], dt)
    step_half = one_step_check(kernel, state, w0["M_total"], e0["Mdot"], 0.5 * dt)
    if state_digest(state) != original_digest:
        raise RuntimeError("frozen initial ProductionState was mutated by observer")

    arithmetic_ok = all(
        item["closure_ok"]
        and item["plus_minus_current_projection_closure_ok"]
        and item["Mdot_geom_tangent_closure_ok"]
        for item in tangent_results
    ) and bool(w0["J_closure_ok"] and w0["M_geom_closure_ok"])
    return {
        "resolution": int(resolution),
        "r_max": float(r_max),
        "dr": float(state.grid.dr),
        "initial_state_sha256": original_digest,
        "initial_rhs_rate_samples": {
            "geometry_full_rhs": {
                name: rhs_rate_witness(rates["geometry"][name], r)
                for name in GEOMETRY_RATE_NAMES
            },
            "Aa_primary_l2": rhs_rate_witness(
                rates["geometry_stage_terms"]["primary_l2"]["Aa"], r
            ),
            "Aa_primary_l3_including_matter": rhs_rate_witness(
                rates["geometry_stage_terms"]["primary_l3"]["Aa"], r
            ),
            "Aa_full_over_r_squared": rhs_rate_witness(
                rates["geometry"]["Aa"] / (r * r), r
            ),
            "K_primary_l2": rhs_rate_witness(
                rates["geometry_stage_terms"]["primary_l2"]["K"], r
            ),
            "K_primary_l3_including_matter": rhs_rate_witness(
                rates["geometry_stage_terms"]["primary_l3"]["K"], r
            ),
        },
        "initial_momentum_constraint": regional_summaries(w0["M_total"], r),
        "initial_geometric_momentum_term_closure": {
            "ok": w0["M_geom_closure_ok"],
            "max_abs": float(np.max(w0["M_geom_closure_abs"])),
            "max_bound": float(np.max(w0["M_geom_closure_bound"])),
        },
        "initial_current_projection_closure": {
            "ok": w0["J_closure_ok"],
            "max_abs": float(np.max(w0["J_closure_abs"])),
            "max_bound": float(np.max(w0["J_closure_bound"])),
        },
        "tangent_by_epsilon": [
            compact_tangent_result(item, r) for item in tangent_results
        ],
        "epsilon_sensitivity": {
            "max_abs_Mdot_difference": float(np.max(np.abs(eps_difference))),
            "regional_Mdot_difference": regional_summaries(eps_difference, r),
            "fixed_radius_samples": fixed_radius_samples(eps_difference, r),
        },
        "bounded_one_step_consistency": {
            "same_common_dt_all_resolutions": True,
            "dt": step_full,
            "dt_half": step_half,
            "interpretation": (
                "local tangent consistency only; a finite-step discrepancy "
                "is not by itself a continuum or physics verdict"
            ),
        },
        "arithmetic_gate_passed": bool(arithmetic_ok),
        "scientific_admission": False,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    prov = provenance()
    if len(args.epsilons) != 2 or any(e <= 0.0 for e in args.epsilons):
        raise ValueError("provide exactly two positive directional epsilons")
    if not args.resolutions or any(n < 4 for n in args.resolutions):
        raise ValueError("resolutions must each be at least four")
    if not math.isfinite(args.dt) or args.dt <= 0.0:
        raise ValueError("common dt must be finite and positive")
    kernel = V55ProductionKernel()
    per_resolution = [
        audit_resolution(
            kernel, n, args.r_max,
            (float(args.epsilons[0]), float(args.epsilons[1])),
            float(args.dt),
        )
        for n in args.resolutions
    ]
    arithmetic_ok = all(item["arithmetic_gate_passed"] for item in per_resolution)
    return {
        "schema": SCHEMA,
        "status": "completed_diagnostic_only" if arithmetic_ok else "arithmetic_gate_failed",
        "provenance": prov,
        "initialization": {
            "r_max": float(args.r_max),
            "amplitude": 0.01,
            "width": 7.0,
            "D_amplitude": 1.0e-10,
            "radiation": True,
            "beta_dm": BETA_DM,
            "lambda_m": LAMBDA_M,
            "common_dt": float(args.dt),
            "directional_epsilons": [float(x) for x in args.epsilons],
        },
        "resolutions": per_resolution,
        "arithmetic_gate_passed": bool(arithmetic_ok),
        "scientific_admission": False,
        "interpretation_limits": [
            "This is an initialized-state, semi-discrete tangent audit, not a long trajectory.",
            "The source-sector Mdot terms are total directional derivatives of each projected J sector along the full state and CMC-lapse direction; they are not frozen-metric partial derivatives.",
            "The signed decomposition diagnoses the implemented discrete RHS only. It does not identify a faulty physical equation by itself.",
            "The geometric subterm budget is a term-by-term derivative of the pinned vendor's exact momentum-constraint expression, with no alternate stencil."
            "The one-step and half-step checks are local finite-step comparisons, not production scientific admission.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[60, 80, 120])
    parser.add_argument("--r-max", type=float, default=80.0)
    parser.add_argument("--dt", type=float, default=0.015)
    parser.add_argument("--epsilons", nargs=2, type=float, default=[2.0e-5, 1.0e-5])
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/current-source-momentum-tangent/report.json"),
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        report = run(args)
        code = 0 if report["arithmetic_gate_passed"] else 2
    except Exception as exc:
        report = {
            "schema": SCHEMA,
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
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
