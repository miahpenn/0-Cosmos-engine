#!/usr/bin/env python3
"""Single-profile, no-evolution Issue #4 operator/parity discriminator.

Runs against one exact production-source commit and the pinned vendor. A smooth
regular manufactured state is evaluated with the production/vendor momentum
constraint and the first-party native diagnostic constraint. The tool records
full signed arrays, local-S current and matching Aa/K stress contributions,
instantaneous RHS tangents, central-cell witnesses, and fixed-radius annulus
summaries. It never calls kernel.step and never edits a production state.
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
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "0c40172f32b72a8149549adadf63d668a2f24f3e"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_ROOT = Path(
    os.environ.get(
        "COSMOS_PINNED_SOURCE_ROOT",
        str(ROOT / "source" if (ROOT / "source").is_dir() else ROOT),
    )
).resolve()
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from engine import v55_pirk_adapter as adapter  # noqa: E402
from engine.matter_system import ConservedSpecies, Species, project_species  # noqa: E402
from engine.production_kernel import ProductionState, V55ProductionKernel  # noqa: E402
from engine.scalar_system import KAPPA, ScalarFields  # noqa: E402
from engine.stress_energy import assemble_total_stress_energy  # noqa: E402
from engine.v55_matter import V55MatterState, metric_slice_from_q  # noqa: E402
from engine.native_geometry import (  # noqa: E402
    NativeSphericalCellGrid,
    native_bssn_constraints,
)

EPS = float(np.finfo(np.float64).eps)
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = ("dark_matter", "baryons", "radiation")
CONSERVED_FIELDS = ("rest", "energy_t", "momentum_r")
EPSILONS = (2.0e-5, 1.0e-5)
FIXED_RADII = (2.0, 3.0, 5.0, 8.0)


def git(path: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(path), *args], text=True).strip()


def digest_array_map(arrays: dict[str, np.ndarray]) -> str:
    h = hashlib.sha256()
    for name in sorted(arrays):
        value = np.ascontiguousarray(np.asarray(arrays[name]))
        h.update(name.encode("utf-8"))
        h.update(value.dtype.str.encode("ascii"))
        h.update(json.dumps(list(value.shape), separators=(",", ":")).encode("ascii"))
        h.update(value.tobytes(order="C"))
    return h.hexdigest()


def state_arrays(state: ProductionState) -> dict[str, np.ndarray]:
    if state.electromagnetic is not None or state.matter.charged_dust is not None:
        raise RuntimeError("this six-sector discriminator does not accept EM/charged-dust states")
    out: dict[str, np.ndarray] = {"r": np.asarray(state.grid.centers, dtype=np.float64).copy()}
    for name in GEOMETRY_FIELDS:
        out["geometry." + name] = np.asarray(getattr(state.geometry, name), dtype=np.float64).copy()
    for name in SCALAR_FIELDS:
        out["scalars." + name] = np.asarray(getattr(state.scalars, name), dtype=np.float64).copy()
    for sector in SPECIES:
        item = getattr(state.matter, sector)
        for name in CONSERVED_FIELDS:
            out[f"matter.{sector}.{name}"] = np.asarray(getattr(item, name), dtype=np.float64).copy()
    return out


def state_digest(state: ProductionState) -> str:
    return digest_array_map(state_arrays(state))


def d1(grid: Any, values: Any, parity: int) -> np.ndarray:
    out = np.asarray(grid.cell_derivative_fourth(values, parity=parity), dtype=np.float64)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError("non-finite first-derivative output")
    return out


def d2(grid: Any, values: Any, parity: int) -> np.ndarray:
    out = np.asarray(grid.cell_second_derivative_fourth(values, parity=parity), dtype=np.float64)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError("non-finite second-derivative output")
    return out


def roundoff_check(residual: np.ndarray, *operands: np.ndarray, multiplier: float = 512.0) -> dict:
    scale = np.zeros_like(np.asarray(residual, dtype=np.float64))
    for operand in operands:
        scale += np.abs(np.asarray(operand, dtype=np.float64))
    bound = multiplier * EPS * scale
    abs_residual = np.abs(np.asarray(residual, dtype=np.float64))
    okay = bool(np.all(abs_residual <= bound))
    ratio = np.divide(
        abs_residual, bound,
        out=np.zeros_like(abs_residual),
        where=bound > 0.0,
    )
    return {
        "passed": okay,
        "max_abs": float(np.max(abs_residual)),
        "max_bound": float(np.max(bound)),
        "max_scaled_ratio": float(np.max(ratio)),
    }


def witness(values: Any, r: np.ndarray) -> dict:
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape != r.shape or not np.all(np.isfinite(arr)):
        raise ValueError("witness must be a finite profile on the supplied grid")

    def peak(mask: np.ndarray) -> dict | None:
        ids = np.flatnonzero(mask)
        if ids.size == 0:
            return None
        i = int(ids[np.argmax(np.abs(arr[ids]))])
        return {
            "index": i,
            "radius": float(r[i]),
            "signed_value": float(arr[i]),
            "max_abs": float(abs(arr[i])),
        }

    fixed = {
        str(target): float(np.interp(target, r, arr))
        for target in FIXED_RADII if float(r[0]) <= target <= float(r[-1])
    }
    return {
        "signed_cells_0_1_2_6": {str(i): float(arr[i]) for i in (0, 1, 2, 6) if i < arr.size},
        "first_two_centre_cells": peak(np.arange(r.size) < 2),
        "r_le_10": peak(r <= 10.0 + 1.0e-12),
        "fixed_annulus_5_to_10": peak((r >= 5.0) & (r <= 10.0)),
        "fixed_radius_samples": fixed,
        "rms_full_grid": float(np.sqrt(np.mean(arr * arr))),
    }


def provenance() -> dict:
    source = git(SOURCE_ROOT, "rev-parse", "HEAD")
    source_vendor = git(SOURCE_ROOT / "vendor" / "bb-palatini-unified-r0", "rev-parse", "HEAD")
    if source != SOURCE_COMMIT:
        raise RuntimeError(f"production source commit mismatch: expected {SOURCE_COMMIT}, got {source}")
    if source_vendor != VENDOR_COMMIT:
        raise RuntimeError(f"production vendor mismatch: expected {VENDOR_COMMIT}, got {source_vendor}")
    # If a separate source checkout is used, the candidate checkout must carry
    # the same submodule commit. This protects the imported adapter from drift.
    candidate_vendor = git(ROOT / "vendor" / "bb-palatini-unified-r0", "rev-parse", "HEAD")
    if candidate_vendor != VENDOR_COMMIT:
        raise RuntimeError(f"candidate vendor mismatch: expected {VENDOR_COMMIT}, got {candidate_vendor}")
    adapter.vendor_modules()
    return {
        "production_source_commit": source,
        "production_source_vendor_commit": source_vendor,
        "diagnostic_commit": git(ROOT, "rev-parse", "HEAD"),
        "candidate_vendor_commit": candidate_vendor,
    }


def manufactured_state(kernel: V55ProductionKernel, resolution: int, r_max: float) -> ProductionState:
    state = kernel.initialize(resolution=resolution, r_max=r_max, include_radiation=True)
    grid = state.grid
    g = state.geometry.copy()
    r = np.asarray(grid.centers, dtype=np.float64)
    q = (r / 6.0) ** 2
    e = np.exp(-q)

    # Smooth regular spherical parity class: metric and scalar profiles are
    # even in r; Aa starts as O(r^2); Lambda is odd. This is a control profile,
    # not constraint-solved physical initial data.
    f = 2.0e-3 * q * e
    g.a = np.exp(f)
    g.b = np.exp(-0.5 * f)
    g.X = 1.0 + 1.5e-3 * e
    g.Aa = 2.0e-5 * r**2 * e
    g.K = -2.0e-2 + 2.0e-3 * e
    g.beta.fill(0.0)
    g.B.fill(0.0)

    g.Lambda = (
        d1(grid, g.a, +1) / (2.0 * g.a**2)
        - d1(grid, g.b, +1) / (g.a * g.b)
        + 2.0 / r * (1.0 / g.b - 1.0 / g.a)
    )
    state.scalars = ScalarFields(
        S=1.5e-2 * e,
        PS=5.0e-4 * e,
        D=1.0e-9 * e,
        PD=2.0e-9 * e,
        phi=0.179055 + 1.0e-3 * e,
        Pi=0.00341913 + 2.0e-4 * e,
    )
    g.alpha = kernel._solve_lapse_for_state(
        grid, g, state.scalars, state.matter, state.electromagnetic
    )[0]
    g.beta.fill(0.0)
    g.B.fill(0.0)
    g.assert_finite_positive()
    state.geometry = g
    state.t = 0.0
    return state


def momentum_pieces_production(grid: Any, g: Any) -> dict[str, np.ndarray]:
    r = np.asarray(grid.centers, dtype=np.float64)
    Aa = np.asarray(g.Aa, dtype=np.float64)
    Ab = -0.5 * Aa
    X = np.asarray(g.X, dtype=np.float64)
    b = np.asarray(g.b, dtype=np.float64)
    chi_r = -0.5 * d1(grid, X, +1) / X
    b_r = d1(grid, b, +1)
    return {
        "D_r_Aa": d1(grid, Aa, +1),
        "minus_2third_D_r_K": -(2.0 / 3.0) * d1(grid, g.K, +1),
        "6_Aa_chi_r": 6.0 * Aa * chi_r,
        "anisotropy_2_over_r": (Aa - Ab) * (2.0 / r),
        "anisotropy_b_gradient": (Aa - Ab) * b_r / b,
    }


def momentum_pieces_native(grid: Any, g: Any) -> dict[str, np.ndarray]:
    r = np.asarray(grid.centers, dtype=np.float64)
    Aa = np.asarray(g.Aa, dtype=np.float64)
    Ab = -0.5 * Aa
    areal = r * np.sqrt(np.asarray(g.b, dtype=np.float64)) / np.asarray(g.X, dtype=np.float64)
    return {
        "D_r_Aa": d1(grid, Aa, +1),
        "minus_2third_D_r_K": -(2.0 / 3.0) * d1(grid, g.K, +1),
        "areal_radius_chain": 2.0 * (Aa - Ab) * d1(grid, areal, -1) / areal,
    }


def current_pieces(state: ProductionState) -> dict[str, np.ndarray]:
    grid, g, s = state.grid, state.geometry, state.scalars
    metric = metric_slice_from_q(grid, g)
    return {
        "S": -np.asarray(s.PS) * d1(grid, s.S, +1),
        "D": -np.asarray(s.PD) * d1(grid, s.D, +1),
        "COSMOS_phi": -np.asarray(s.Pi) * d1(grid, s.phi, +1) / KAPPA,
        "dark_matter": np.asarray(project_species(metric, state.matter.dark_matter, Species.DARK_MATTER)["j"]),
        "baryons": np.asarray(project_species(metric, state.matter.baryons, Species.BARYON)["j"]),
        "radiation": np.asarray(project_species(metric, state.matter.radiation, Species.RADIATION)["j"]),
    }


def geometry_rates_and_source_terms(kernel: V55ProductionKernel, state: ProductionState) -> tuple[dict, dict]:
    g, grid = state.geometry, state.grid
    terms = kernel._geometry_stage_terms(
        grid, g, state.scalars, state.matter, state.electromagnetic
    )
    srhs, mrhs, _ = kernel._rhs(state)
    rates = {
        "geometry": {
            "a": np.asarray(terms["explicit"]["a"], dtype=np.float64),
            "b": np.asarray(terms["explicit"]["b"], dtype=np.float64),
            "X": np.asarray(terms["explicit"]["X"], dtype=np.float64),
            "Aa": np.asarray(terms["primary_l2"]["Aa"]) + np.asarray(terms["primary_l3"]["Aa"]),
            "K": np.asarray(terms["primary_l2"]["K"]) + np.asarray(terms["primary_l3"]["K"]),
            "Lambda": np.asarray(terms["lambda_l2"]) + np.asarray(terms["lambda_l3"]),
        },
        "scalars": {name: np.asarray(getattr(srhs, name), dtype=np.float64) for name in SCALAR_FIELDS},
        "matter": {
            sector: {
                name: np.asarray(getattr(mrhs[sector], name), dtype=np.float64)
                for name in CONSERVED_FIELDS
            } for sector in SPECIES
        },
    }
    return rates, {"stage_terms": terms, "scalar_rhs": srhs, "matter_rhs": mrhs}


def scalar_s_stress(state: ProductionState) -> dict[str, np.ndarray]:
    g, s, grid = state.geometry, state.scalars, state.grid
    grad = d1(grid, s.S, +1)
    grad_sq = (g.X**2 / g.a) * grad**2
    potential = 0.5 * s.S**2
    rho = 0.5 * s.PS**2 + 0.5 * grad_sq + potential
    pr = 0.5 * s.PS**2 + 0.5 * grad_sq - potential
    pt = 0.5 * s.PS**2 - 0.5 * grad_sq - potential
    return {"rho": rho, "pr": pr, "pt": pt, "j": -s.PS * grad}


def perturb_state(
    kernel: V55ProductionKernel, state: ProductionState, rates: dict, delta: float
) -> ProductionState:
    g = state.geometry.copy()
    for name in ("a", "b", "X", "Aa", "K", "Lambda"):
        value = np.asarray(getattr(state.geometry, name), dtype=np.float64) + delta * rates["geometry"][name]
        setattr(g, name, np.asarray(value, dtype=np.float64).copy())
    g.beta.fill(0.0)
    g.B.fill(0.0)
    scalars = ScalarFields(*(
        np.asarray(getattr(state.scalars, name), dtype=np.float64)
        + delta * rates["scalars"][name] for name in SCALAR_FIELDS
    ))
    items = {}
    for sector in SPECIES:
        base = getattr(state.matter, sector)
        items[sector] = ConservedSpecies(*(
            np.asarray(getattr(base, name), dtype=np.float64)
            + delta * rates["matter"][sector][name] for name in CONSERVED_FIELDS
        ))
    matter = V55MatterState(
        dark_matter=items["dark_matter"], baryons=items["baryons"],
        radiation=items["radiation"], charged_dust=None,
    )
    g.alpha = kernel._solve_lapse_for_state(
        state.grid, g, scalars, matter, None
    )[0]
    g.beta.fill(0.0)
    g.B.fill(0.0)
    return ProductionState(
        grid=state.grid, geometry=g, scalars=scalars, matter=matter,
        t=state.t, tau=state.tau, e_folds=state.e_folds,
        electromagnetic=None,
    )


def analytic_tangent_pieces(state: ProductionState, rates: dict) -> dict[str, np.ndarray]:
    grid, g, r = state.grid, state.geometry, np.asarray(state.grid.centers, dtype=np.float64)
    rd = rates["geometry"]
    Aa, X, b = np.asarray(g.Aa), np.asarray(g.X), np.asarray(g.b)
    Aadot, Xdot, bdot = np.asarray(rd["Aa"]), np.asarray(rd["X"]), np.asarray(rd["b"])
    chi = -0.5 * d1(grid, X, +1) / X
    chi_dot = -0.5 * (d1(grid, Xdot, +1) / X - d1(grid, X, +1) * Xdot / X**2)
    q = 2.0 / r + d1(grid, b, +1) / b
    q_dot = d1(grid, bdot, +1) / b - d1(grid, b, +1) * bdot / b**2
    prod_geom = (
        d1(grid, Aadot, +1) - (2.0 / 3.0) * d1(grid, rd["K"], +1)
        + 6.0 * Aadot * chi + 6.0 * Aa * chi_dot
        + 1.5 * Aadot * q + 1.5 * Aa * q_dot
    )

    areal = r * np.sqrt(b) / X
    areal_dot = areal * (0.5 * bdot / b - Xdot / X)
    native_geom = (
        d1(grid, Aadot, +1) - (2.0 / 3.0) * d1(grid, rd["K"], +1)
        + 3.0 * Aadot * d1(grid, areal, -1) / areal
        + 3.0 * Aa * (
            d1(grid, areal_dot, -1) / areal
            - d1(grid, areal, -1) * areal_dot / areal**2
        )
    )

    # Differentiate each sector's implemented current projection on this frozen
    # stage; fluid slots carry the sqrt(gamma) time-density term.
    s, sd = state.scalars, rates["scalars"]
    jdot = {
        "S": -(sd["PS"] * d1(grid, s.S, +1) + s.PS * d1(grid, sd["S"], +1)),
        "D": -(sd["PD"] * d1(grid, s.D, +1) + s.PD * d1(grid, sd["D"], +1)),
        "COSMOS_phi": -(
            sd["Pi"] * d1(grid, s.phi, +1) + s.Pi * d1(grid, sd["phi"], +1)
        ) / KAPPA,
    }
    volume_rate = 0.5 * rd["a"] / g.a + rd["b"] / g.b - 3.0 * rd["X"] / g.X
    sqrt_gamma = np.sqrt(g.a) * g.b * r**2 / g.X**3
    base_j = current_pieces(state)
    for sector in SPECIES:
        Udot = rates["matter"][sector]["momentum_r"]
        jdot[sector] = Udot / (KAPPA * sqrt_gamma) - base_j[sector] * volume_rate
    total_jdot = sum(jdot.values(), np.zeros_like(r))
    return {
        "production_M_geometry_dot": prod_geom,
        "native_M_geometry_dot": native_geom,
        **{"jdot_" + key: value for key, value in jdot.items()},
        "total_j_dot": total_jdot,
        "production_M_total_dot": prod_geom - KAPPA * total_jdot,
        "native_M_total_dot": native_geom - KAPPA * total_jdot,
    }


def observed_profiles(
    state: ProductionState, vendor: Any, native_grid: NativeSphericalCellGrid
) -> dict[str, Any]:
    raw = vendor.constraints(state.grid, state.geometry)
    production_geom_pieces = momentum_pieces_production(state.grid, state.geometry)
    native_geom_pieces = momentum_pieces_native(native_grid, state.geometry)
    prod_geom = np.asarray(raw["momentum"], dtype=np.float64)
    native_constraint = native_bssn_constraints(native_grid, state.geometry)
    native_geom = np.asarray(native_constraint["M_geometry"], dtype=np.float64)
    current = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    Jpieces = current_pieces(state)
    Jsum = sum(Jpieces.values(), np.zeros_like(prod_geom))
    return {
        "production_M_geometry": prod_geom,
        "native_M_geometry": native_geom,
        "production_M_geometry_pieces": production_geom_pieces,
        "native_M_geometry_pieces": native_geom_pieces,
        "J_total": np.asarray(current.j, dtype=np.float64),
        "J_pieces": Jpieces,
        "J_sum": Jsum,
        "production_M_total": prod_geom - KAPPA * np.asarray(current.j),
        "native_M_total": native_geom - KAPPA * np.asarray(current.j),
        "native_constraint_result": native_constraint,
    }


def run_audit(resolution: int = 64, r_max: float = 40.0) -> tuple[dict, dict[str, np.ndarray]]:
    pin = provenance()
    kernel = V55ProductionKernel()
    state = manufactured_state(kernel, resolution, r_max)
    before = state_digest(state)
    grid = state.grid
    native_grid = NativeSphericalCellGrid(resolution, r_max)
    if not np.allclose(native_grid.centers, grid.centers, rtol=0.0, atol=2.0e-14):
        raise RuntimeError("native and production grid centres differ")
    _, vendor, _ = adapter.vendor_modules()
    rates, rhs = geometry_rates_and_source_terms(kernel, state)
    after_rhs = state_digest(state)
    if before != after_rhs:
        raise RuntimeError("RHS extraction mutated the manufactured input state")

    # First establish whether shared field derivatives differ at all.
    r = np.asarray(grid.centers, dtype=np.float64)
    areal = r * np.sqrt(state.geometry.b) / state.geometry.X
    parity_fields = {
        "Aa_even": (state.geometry.Aa, +1),
        "K_even": (state.geometry.K, +1),
        "a_even": (state.geometry.a, +1),
        "b_even": (state.geometry.b, +1),
        "X_even": (state.geometry.X, +1),
        "S_even": (state.scalars.S, +1),
        "D_even": (state.scalars.D, +1),
        "phi_even": (state.scalars.phi, +1),
        "Lambda_odd": (state.geometry.Lambda, -1),
        "areal_radius_odd": (areal, -1),
    }
    operator_profiles: dict[str, np.ndarray] = {}
    operator_checks = {}
    for name, (values, parity) in parity_fields.items():
        prod = d1(grid, values, parity)
        native = d1(native_grid, values, parity)
        diff = prod - native
        operator_profiles["D1_production_" + name] = prod
        operator_profiles["D1_native_" + name] = native
        operator_profiles["D1_diff_" + name] = diff
        operator_checks[name] = roundoff_check(diff, prod, native, multiplier=128.0)
        operator_checks[name]["declared_input_parity"] = "even" if parity == 1 else "odd"
    for name, values, parity in (
        ("S_even", state.scalars.S, +1),
        ("alpha_even", state.geometry.alpha, +1),
        ("areal_radius_odd", areal, -1),
    ):
        prod = d2(grid, values, parity)
        native = d2(native_grid, values, parity)
        diff = prod - native
        operator_profiles["D2_production_" + name] = prod
        operator_profiles["D2_native_" + name] = native
        operator_profiles["D2_diff_" + name] = diff
        operator_checks["D2_" + name] = roundoff_check(diff, prod, native, multiplier=128.0)
        operator_checks["D2_" + name]["declared_input_parity"] = "even" if parity == 1 else "odd"
    operators_pass = all(item["passed"] for item in operator_checks.values())

    base = observed_profiles(state, vendor, native_grid)
    prod_piece_sum = sum(base["production_M_geometry_pieces"].values(), np.zeros_like(r))
    native_piece_sum = sum(base["native_M_geometry_pieces"].values(), np.zeros_like(r))
    prod_decomp = roundoff_check(
        base["production_M_geometry"] - prod_piece_sum,
        base["production_M_geometry"], *base["production_M_geometry_pieces"].values()
    )
    native_decomp = roundoff_check(
        base["native_M_geometry"] - native_piece_sum,
        base["native_M_geometry"], *base["native_M_geometry_pieces"].values()
    )
    current_projection_closure = roundoff_check(
        base["J_total"] - base["J_sum"], base["J_total"], *base["J_pieces"].values()
    )
    geometric_representation_delta = base["native_M_geometry"] - base["production_M_geometry"]
    source_closure_arrays = {}
    total_stress = assemble_total_stress_energy(
        grid, state.geometry, state.scalars, state.matter
    )
    total_terms = rhs["stage_terms"]
    vacuum_l3 = vendor.primary_l3_rhs(grid, state.geometry)
    aa_source_actual = np.asarray(total_terms["primary_l3"]["Aa"]) - np.asarray(vacuum_l3["Aa"])
    k_source_actual = np.asarray(total_terms["primary_l3"]["K"]) - np.asarray(vacuum_l3["K"])
    aa_source_expected = -(16.0 * math.pi / 3.0) * state.geometry.alpha * (
        total_stress.pr - total_stress.pt
    )
    k_source_expected = 4.0 * math.pi * state.geometry.alpha * (
        total_stress.rho + total_stress.pr + 2.0 * total_stress.pt
    )
    aa_source_closure = roundoff_check(
        aa_source_actual - aa_source_expected, aa_source_actual, aa_source_expected
    )
    k_source_closure = roundoff_check(
        k_source_actual - k_source_expected, k_source_actual, k_source_expected
    )

    s_stress = scalar_s_stress(state)
    aa_s_source = -(16.0 * math.pi / 3.0) * state.geometry.alpha * (
        s_stress["pr"] - s_stress["pt"]
    )
    k_s_source = 4.0 * math.pi * state.geometry.alpha * (
        s_stress["rho"] + s_stress["pr"] + 2.0 * s_stress["pt"]
    )
    g = state.geometry
    chi = -0.5 * d1(grid, g.X, +1) / g.X
    q = 2.0 / r + d1(grid, g.b, +1) / g.b
    local_s_production_geom_dot = (
        d1(grid, aa_s_source, +1) - (2.0 / 3.0) * d1(grid, k_s_source, +1)
        + aa_s_source * (6.0 * chi + 1.5 * q)
    )
    local_s_native_geom_dot = (
        d1(native_grid, aa_s_source, +1) - (2.0 / 3.0) * d1(native_grid, k_s_source, +1)
        + 3.0 * aa_s_source * d1(native_grid, areal, -1) / areal
    )
    local_s_jdot = -(
        np.asarray(rhs["scalar_rhs"].PS) * d1(grid, state.scalars.S, +1)
        + state.scalars.PS * d1(grid, rhs["scalar_rhs"].S, +1)
    )
    local_s_production_Mdot = local_s_production_geom_dot - KAPPA * local_s_jdot
    local_s_native_Mdot = local_s_native_geom_dot - KAPPA * local_s_jdot

    analytic = analytic_tangent_pieces(state, rates)
    measured = []
    arrays: dict[str, np.ndarray] = {
        "r": r.copy(),
        **operator_profiles,
        "production_M_geometry": base["production_M_geometry"],
        "native_M_geometry": base["native_M_geometry"],
        "production_minus_native_M_geometry": base["production_M_geometry"] - base["native_M_geometry"],
        "J_total": base["J_total"],
        "J_total_sum_from_sectors": base["J_sum"],
        "production_M_total": base["production_M_total"],
        "native_M_total": base["native_M_total"],
        **{"production_M_geometry_term_" + k: v for k, v in base["production_M_geometry_pieces"].items()},
        **{"native_M_geometry_term_" + k: v for k, v in base["native_M_geometry_pieces"].items()},
        **{"J_sector_" + k: v for k, v in base["J_pieces"].items()},
        **{"Jdot_analytic_" + k: v for k, v in analytic.items() if k.startswith("jdot_")},
        "Jdot_analytic_total": analytic["total_j_dot"],
        "Mdot_geometry_analytic_production": analytic["production_M_geometry_dot"],
        "Mdot_geometry_analytic_native": analytic["native_M_geometry_dot"],
        "Mdot_total_analytic_production": analytic["production_M_total_dot"],
        "Mdot_total_analytic_native": analytic["native_M_total_dot"],
        "local_S_rho": s_stress["rho"],
        "local_S_pr": s_stress["pr"],
        "local_S_pt": s_stress["pt"],
        "local_S_J": s_stress["j"],
        "local_S_Jdot_analytic": local_s_jdot,
        "local_S_Aa_rhs_source": aa_s_source,
        "local_S_K_rhs_source": k_s_source,
        "local_S_Mdot_geometry_production": local_s_production_geom_dot,
        "local_S_Mdot_geometry_native": local_s_native_geom_dot,
        "local_S_Mdot_total_production": local_s_production_Mdot,
        "local_S_Mdot_total_native": local_s_native_Mdot,
        "Aa_matter_source_actual": aa_source_actual,
        "Aa_matter_source_expected": aa_source_expected,
        "K_matter_source_actual": k_source_actual,
        "K_matter_source_expected": k_source_expected,
    }

    for epsilon in EPSILONS:
        plus = perturb_state(kernel, state, rates, epsilon)
        minus = perturb_state(kernel, state, rates, -epsilon)
        p = observed_profiles(plus, vendor, native_grid)
        m = observed_profiles(minus, vendor, native_grid)
        denom = 2.0 * epsilon
        fd = {
            "production_M_total": (p["production_M_total"] - m["production_M_total"]) / denom,
            "native_M_total": (p["native_M_total"] - m["native_M_total"]) / denom,
            "production_M_geometry": (p["production_M_geometry"] - m["production_M_geometry"]) / denom,
            "native_M_geometry": (p["native_M_geometry"] - m["native_M_geometry"]) / denom,
        }
        fd["J_total"] = (p["J_total"] - m["J_total"]) / denom
        for sector in SPECIES + ("S", "D", "COSMOS_phi"):
            fd["J_sector_" + sector] = (p["J_pieces"][sector] - m["J_pieces"][sector]) / denom
        fd["production_reassembled"] = fd["production_M_geometry"] - KAPPA * fd["J_total"]
        fd["native_reassembled"] = fd["native_M_geometry"] - KAPPA * fd["J_total"]
        fd["production_budget_closure"] = fd["production_M_total"] - fd["production_reassembled"]
        fd["native_budget_closure"] = fd["native_M_total"] - fd["native_reassembled"]
        for name, values in fd.items():
            arrays[f"fd_eps_{epsilon:.0e}_" + name] = np.asarray(values)
        measured.append({
            "epsilon": epsilon,
            "production_Mdot_measured": witness(fd["production_M_total"], r),
            "native_Mdot_measured": witness(fd["native_M_total"], r),
            "production_vs_analytic_tangent": witness(
                fd["production_M_total"] - analytic["production_M_total_dot"], r
            ),
            "native_vs_analytic_tangent": witness(
                fd["native_M_total"] - analytic["native_M_total_dot"], r
            ),
            "production_budget_closure": roundoff_check(
                fd["production_budget_closure"],
                fd["production_M_total"], fd["production_M_geometry"], KAPPA * fd["J_total"]
            ),
            "native_budget_closure": roundoff_check(
                fd["native_budget_closure"],
                fd["native_M_total"], fd["native_M_geometry"], KAPPA * fd["J_total"]
            ),
            "input_state_digest_plus": state_digest(plus),
            "input_state_digest_minus": state_digest(minus),
        })

    after_all = state_digest(state)
    if before != after_all:
        raise RuntimeError("observer mutated the frozen manufactured input state")
    if not operators_pass:
        raise RuntimeError("production/vendor and first-party shared D1/D2 operator parity checks failed")
    if not prod_decomp["passed"] or not native_decomp["passed"]:
        raise RuntimeError("signed momentum-constraint term decomposition failed roundoff closure")
    if not current_projection_closure["passed"]:
        raise RuntimeError("six-sector current projection does not reconstruct the canonical total current")
    if not aa_source_closure["passed"] or not k_source_closure["passed"]:
        raise RuntimeError("production Aa/K matter-source terms do not reconstruct their declared stress projections")

    report = {
        "schema": "issue4_manufactured_profile_rhs_stencil_parity_v1",
        "status": "diagnostic_only_not_scientific_admission",
        "provenance": pin,
        "configuration": {
            "resolution": resolution,
            "r_max": r_max,
            "dr": r_max / resolution,
            "profile": "smooth regular even/odd manufactured spherical profile; not constraint-solved",
            "directional_epsilons": list(EPSILONS),
            "production_steps_advanced": 0,
            "evolution_advanced": False,
            "production_equations_modified": False,
        },
        "input_state": {
            "digest_before": before,
            "digest_after_rhs": after_rhs,
            "digest_after_all": after_all,
            "immutable": before == after_rhs == after_all,
        },
        "decision": {
            "shared_operator_parity_check": {
                "passed": operators_pass,
                "checks": operator_checks,
                "interpretation": "Compares the actual production/vendor grid operators with the first-party native grid on identical arrays and declared input parity. A pass rules out a difference in these shared D1/D2 stencils as the explanation for any constraint-representation gap.",
            },
            "production_constraint_decomposition": prod_decomp,
            "native_constraint_decomposition": native_decomp,
            "current_projection_reconstruction": current_projection_closure,
            "primary_l3_stress_source_reconstruction": {
                "Aa": aa_source_closure,
                "K": k_source_closure,
            },
            "production_minus_native_geometry": witness(geometric_representation_delta, r),
            "production_minus_native_total_constraint": witness(
                base["production_M_total"] - base["native_M_total"], r
            ),
            "interpretation": (
                "A nonzero production/native difference with matched D1/D2 operators is a discrete formula/composition difference (not a stencil/parity mismatch). "
                "The signed local-S current and its paired Aa/K stress contributions are recorded separately. This control does not close the full matter-coupled propagation identity or authorize a production repair."
            ),
        },
        "base_state": {
            "production_M_geometry": witness(base["production_M_geometry"], r),
            "native_M_geometry": witness(base["native_M_geometry"], r),
            "production_M_total": witness(base["production_M_total"], r),
            "native_M_total": witness(base["native_M_total"], r),
            "sector_current_pieces": {name: witness(value, r) for name, value in base["J_pieces"].items()},
        },
        "local_S_source_budget": {
            "J_S": witness(s_stress["j"], r),
            "Jdot_S_analytic": witness(local_s_jdot, r),
            "Aa_rhs_stress_source": witness(aa_s_source, r),
            "K_rhs_stress_source": witness(k_s_source, r),
            "production_geometry_contribution_to_Mdot": witness(local_s_production_geom_dot, r),
            "native_geometry_contribution_to_Mdot": witness(local_s_native_geom_dot, r),
            "production_total_contribution_to_Mdot": witness(local_s_production_Mdot, r),
            "native_total_contribution_to_Mdot": witness(local_s_native_Mdot, r),
        },
        "source_projection": {
            "Aa_actual_minus_expected": aa_source_closure,
            "K_actual_minus_expected": k_source_closure,
        },
        "tangent": {
            "analytic_production_Mdot": witness(analytic["production_M_total_dot"], r),
            "analytic_native_Mdot": witness(analytic["native_M_total_dot"], r),
            "analytic_sector_Jdot": {
                name: witness(analytic["jdot_" + name], r)
                for name in SPECIES + ("S", "D", "COSMOS_phi")
            },
            "finite_difference_arms": measured,
        },
        "admission": {
            "full_matter_coupled_propagation_identity_closed": False,
            "total_stress_energy_conservation_closed": False,
            "scientific_admission": False,
        },
        "limitations": [
            "One manufactured profile at one resolution is a discriminator, not a convergence campaign.",
            "The native and production constraint formulas are different discrete representations; their difference is reported, not forced to zero.",
            "Arithmetic closure validates reconstruction of the implemented budget, not the physical validity of the source law.",
            "The saved accepted-state failures and source-convention obligations in PRs #420-#424 remain open and are not overridden by this control.",
        ],
    }
    return report, arrays


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=64)
    parser.add_argument("--r-max", type=float, default=40.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.resolution < 24:
        parser.error("resolution must be at least 24 for the fixed-radius and parity witnesses")
    if not math.isfinite(args.r_max) or args.r_max <= 10.0:
        parser.error("r-max must be finite and greater than 10")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report, arrays = run_audit(args.resolution, args.r_max)
    arrays_path = args.output_dir / "signed_arrays.npz"
    np.savez_compressed(arrays_path, **arrays)
    report["arrays"] = {
        "file": arrays_path.name,
        "sha256": hashlib.sha256(arrays_path.read_bytes()).hexdigest(),
        "array_count": len(arrays),
        "array_names": sorted(arrays),
    }
    report_path = args.output_dir / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "report.json": hashlib.sha256(report_path.read_bytes()).hexdigest(),
        "signed_arrays.npz": report["arrays"]["sha256"],
    }
    (args.output_dir / "SHA256SUMS.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": report["status"],
        "source": report["provenance"]["production_source_commit"],
        "vendor": report["provenance"]["production_source_vendor_commit"],
        "resolution": args.resolution,
        "shared_operator_parity_pass": report["decision"]["shared_operator_parity_check"]["passed"],
        "production_constraint_decomposition_pass": report["decision"]["production_constraint_decomposition"]["passed"],
        "native_constraint_decomposition_pass": report["decision"]["native_constraint_decomposition"]["passed"],
        "input_immutable": report["input_state"]["immutable"],
        "report": str(report_path),
        "arrays": str(arrays_path),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
