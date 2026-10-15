#!/usr/bin/env python3
"""Read-only manufactured-profile audit of the production discrete momentum tangent.

No production timestep is advanced. The audit checks:
  * parity-aware fourth-order derivatives on analytic regular even/odd profiles;
  * exact signed decomposition of the discrete momentum constraint;
  * local-S current and its correlated stress / geometric RHS contributions;
  * finite-difference of the complete production semi-discrete RHS direction,
    with the production CMC lapse re-solved at plus/minus perturbed states.

A manufactured profile is not asserted to satisfy the Einstein constraints.
A nonzero momentum tangent is therefore a diagnostic witness, not by itself
an error. The required identities are algebraic closures of the declared
discrete operators and projections.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

# Production API imports resolve against the separately checked-out exact source pin.
SOURCE_COMMIT = "d948c7ee692b72abf4cab759b76475e636838a1a"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EPSILONS = (2.0e-5, 1.0e-5)
RESOLUTION = 120
R_MAX = 40.0
R_CUT_INNER = 5.0
FIXED_RADII = (2.0, 3.0, 4.0, 5.0, 8.0)
SCALAR_NAMES = ("S", "PS", "D", "PD", "phi", "Pi")
GEOMETRY_RATE_NAMES = ("a", "b", "X", "Aa", "K", "Lambda")
SPECIES_PREFIX = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED_NAMES = ("rest", "energy_t", "momentum_r")

sys.path.insert(0, str(Path("source").resolve()))
from engine.production_kernel import V55ProductionKernel, ProductionState, LAMBDA_M
from engine.scalar_system import ScalarFields, scalar_projection
from engine.matter_system import ConservedSpecies, Species, project_species, validate_species_admissibility
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()


def finite(name: str, values: np.ndarray) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError(f"non-finite values in {name}")
    return out


def d1(grid, values: np.ndarray, parity: int = +1) -> np.ndarray:
    return finite("cell_derivative_fourth", grid.cell_derivative_fourth(values, parity=parity))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def profile_summary(values: np.ndarray, r: np.ndarray) -> dict:
    a = finite("profile_summary", values)
    i = int(np.argmax(np.abs(a)))
    inner = np.flatnonzero(r < R_CUT_INNER)
    ii = int(inner[np.argmax(np.abs(a[inner]))]) if inner.size else i
    return {
        "max_abs": float(np.max(np.abs(a))),
        "peak": {"i": i, "r": float(r[i]), "signed": float(a[i])},
        "inner_r_lt_5_peak": {"i": ii, "r": float(r[ii]), "signed": float(a[ii])},
        "rms": float(np.sqrt(np.mean(a * a))),
    }


def analytic_derivative_audit(grid) -> dict:
    """Analytic parity controls; smooth polynomials avoid fitting an order."""
    r = np.asarray(grid.centers, dtype=float)
    L = 40.0
    x = r / L
    even = 1.0 + 0.05 * x**2 + 0.01 * x**4
    even_exact = 0.1 * r / L**2 + 0.04 * r**3 / L**4
    odd = x + 0.05 * x**3
    odd_exact = 1.0 / L + 0.15 * r**2 / L**3
    records = {}
    for name, values, exact, parity in (
        ("even_quartic", even, even_exact, +1),
        ("odd_cubic", odd, odd_exact, -1),
    ):
        err = d1(grid, values, parity) - exact
        records[name] = {
            "parity": parity,
            "max_abs_derivative_error": float(np.max(np.abs(err))),
            "rms_derivative_error": float(np.sqrt(np.mean(err * err))),
            "first_two_signed_errors": [float(err[0]), float(err[1])],
            "last_two_signed_errors": [float(err[-2]), float(err[-1])],
            "max_error_location_r": float(r[int(np.argmax(np.abs(err)))]),
        }
    records["pass"] = bool(
        records["even_quartic"]["max_abs_derivative_error"] < 5.0e-11
        and records["odd_cubic"]["max_abs_derivative_error"] < 5.0e-11
    )
    if not records["pass"]:
        raise AssertionError(f"analytic parity derivative contract failed: {records}")
    return records


def metric_lambda(grid, geometry) -> np.ndarray:
    r = np.asarray(grid.centers, dtype=float)
    a = finite("metric a", geometry.a)
    b = finite("metric b", geometry.b)
    return d1(grid, a, +1) / (2.0 * a**2) - d1(grid, b, +1) / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)


def make_manufactured_state(kernel: V55ProductionKernel) -> ProductionState:
    """Smooth, regular parity-class data; deliberately not constraint-solved."""
    state = kernel.initialize(resolution=RESOLUTION, r_max=R_MAX, include_radiation=True)
    grid = state.grid
    r = np.asarray(grid.centers, dtype=float)
    q = (r / 6.0) ** 2
    e = np.exp(-q)
    g = state.geometry
    # A regular spherical metric has equal radial/angular coefficients at r=0.
    # Let f=O(r^2), then a=exp(f), b=exp(-f/2) also enforce a*b^2=1.
    f = 2.0e-3 * (r / 6.0) ** 2 * e
    g.a = np.exp(f)
    g.b = np.exp(-0.5 * f)
    g.X = 1.0 + 1.5e-3 * e
    g.Aa = 2.0e-6 * r**2 * e
    g.K = -2.0e-2 + 5.0e-4 * e
    g.beta.fill(0.0)
    g.B.fill(0.0)

    state.scalars = ScalarFields(
        S=1.5e-2 * e,
        PS=5.0e-4 * e,
        D=1.0e-9 * e,
        PD=2.0e-9 * e,
        phi=0.179055 + 1.0e-3 * e,
        Pi=0.00341913 + 2.0e-4 * e,
    )
    # Keep conservative matter fields from the valid initialized seed. The
    # small smooth metric deformation is checked before the CMC projection.
    validate_species_admissibility(
        metric_slice_from_q(grid, g), state.matter.radiation,
        Species.RADIATION, stage="manufactured-audit-base",
    )
    g.Lambda = metric_lambda(grid, g)
    g.alpha = np.asarray(kernel._solve_lapse(grid, g, state.scalars, state.matter)[0], dtype=float)
    g.beta.fill(0.0)
    g.B.fill(0.0)
    g.assert_finite_positive()
    validate_species_admissibility(
        metric_slice_from_q(grid, g), state.matter.radiation,
        Species.RADIATION, stage="manufactured-audit-cmc-base",
    )
    return state


def momentum_terms(grid, geometry, scalars, matter) -> tuple[dict[str, np.ndarray], dict]:
    """Signed terms assembled with the same vendor derivative and projections."""
    r = np.asarray(grid.centers, dtype=float)
    Aa = finite("Aa", geometry.Aa)
    Ab = -0.5 * Aa
    K = finite("K", geometry.K)
    X = finite("X", geometry.X)
    b = finite("b", geometry.b)
    chi_r = -0.5 * d1(grid, X, +1) / X

    geometric = {
        "geo_Dr_Aa": d1(grid, Aa, +1),
        "geo_minus_2third_Dr_K": -(2.0 / 3.0) * d1(grid, K, +1),
        "geo_6_Aa_chi_r": 6.0 * Aa * chi_r,
        "geo_Aa_radial_2_over_r": (Aa - Ab) * (2.0 / r),
        "geo_Aa_b_gradient": (Aa - Ab) * d1(grid, b, +1) / b,
    }
    geo_sum = sum(geometric.values(), np.zeros_like(r))

    # Expand the scalar current using the exact production projection's signs
    # and normalization. Phi is Einstein-normalized by 8*pi; S and D are not.
    j_s = -np.asarray(scalars.PS, dtype=float) * d1(grid, np.asarray(scalars.S, dtype=float), +1)
    j_d = -np.asarray(scalars.PD, dtype=float) * d1(grid, np.asarray(scalars.D, dtype=float), +1)
    j_phi = -np.asarray(scalars.Pi, dtype=float) * d1(grid, np.asarray(scalars.phi, dtype=float), +1) / (8.0 * math.pi)
    metric = metric_slice_from_q(grid, geometry)
    j_dm = np.asarray(project_species(metric, matter.dark_matter, Species.DARK_MATTER)["j"], dtype=float)
    j_baryon = np.asarray(project_species(metric, matter.baryons, Species.BARYON)["j"], dtype=float)
    j_radiation = np.asarray(project_species(metric, matter.radiation, Species.RADIATION)["j"], dtype=float)
    total = assemble_total_stress_energy(grid, geometry, scalars, matter)
    direct = np.asarray(vacuum.constraints(grid, geometry)["momentum"], dtype=float) - 8.0 * math.pi * np.asarray(total.j, dtype=float)

    sources = {
        "src_minus_8piJ_S": -8.0 * math.pi * j_s,
        "src_minus_8piJ_D": -8.0 * math.pi * j_d,
        "src_minus_8piJ_phi": -8.0 * math.pi * j_phi,
        "src_minus_8piJ_dark_matter": -8.0 * math.pi * j_dm,
        "src_minus_8piJ_baryons": -8.0 * math.pi * j_baryon,
        "src_minus_8piJ_radiation": -8.0 * math.pi * j_radiation,
    }
    source_sum = sum(sources.values(), np.zeros_like(r))
    direct_geom = np.asarray(vacuum.constraints(grid, geometry)["momentum"], dtype=float)
    arrays = {**geometric, **sources, "momentum_geometric_sum": geo_sum, "momentum_geometric_vendor": direct_geom,
              "momentum_source_sum": source_sum, "momentum_source_vendor": -8.0 * math.pi * np.asarray(total.j, dtype=float),
              "momentum_residual_direct": direct, "momentum_residual_rebuilt": geo_sum + source_sum}

    closure = {
        "geometry_decomposition_max_abs": float(np.max(np.abs(geo_sum - direct_geom))),
        "source_decomposition_max_abs": float(np.max(np.abs(source_sum + 8.0 * math.pi * np.asarray(total.j, dtype=float)))),
        "direct_vs_rebuilt_max_abs": float(np.max(np.abs(direct - (geo_sum + source_sum)))),
    }
    return {k: finite(k, v) for k, v in arrays.items()}, closure


def local_s_budget(grid, geometry, scalars) -> dict[str, np.ndarray]:
    """Local-S current and the matching source pieces in the A_a/K RHS."""
    invr = np.asarray(geometry.X, dtype=float)**2 / np.asarray(geometry.a, dtype=float)
    Sr = np.asarray(scalars.S, dtype=float)
    PS = np.asarray(scalars.PS, dtype=float)
    Sp = d1(grid, Sr, +1)
    gr2 = invr * Sp**2
    V = 0.5 * Sr**2
    rho_s = 0.5 * PS**2 + 0.5 * gr2 + V
    pr_s = 0.5 * PS**2 + 0.5 * gr2 - V
    pt_s = 0.5 * PS**2 - 0.5 * gr2 - V
    j_s = -PS * Sp
    alpha = np.asarray(geometry.alpha, dtype=float)
    return {
        "local_S_DrS": Sp,
        "local_S_current_J": j_s,
        "local_S_minus8piJ": -8.0 * math.pi * j_s,
        "local_S_rho": rho_s,
        "local_S_pr": pr_s,
        "local_S_pt": pt_s,
        "local_S_Aa_rhs_source": -(16.0 * math.pi / 3.0) * alpha * (pr_s - pt_s),
        "local_S_K_rhs_source": 4.0 * math.pi * alpha * (rho_s + pr_s + 2.0 * pt_s),
    }


def tangent_state(kernel, state, srhs, mrhs, geometry_rates, eps: float) -> ProductionState:
    g = state.geometry.copy()
    for name in GEOMETRY_RATE_NAMES:
        setattr(g, name, np.asarray(getattr(g, name), dtype=float) + eps * np.asarray(geometry_rates[name], dtype=float))
    g.beta.fill(0.0)
    g.B.fill(0.0)
    scalars = ScalarFields(*(
        np.asarray(getattr(state.scalars, name), dtype=float) + eps * np.asarray(getattr(srhs, name), dtype=float)
        for name in SCALAR_NAMES
    ))
    species = {}
    for name in SPECIES_PREFIX:
        base = getattr(state.matter, name)
        rhs = mrhs[name]
        species[name] = ConservedSpecies(*(
            np.asarray(getattr(base, field), dtype=float) + eps * np.asarray(getattr(rhs, field), dtype=float)
            for field in CONSERVED_NAMES
        ))
    matter = V55MatterState(
        dark_matter=species["dark_matter"],
        baryons=species["baryons"],
        radiation=species["radiation"],
    )
    validate_species_admissibility(
        metric_slice_from_q(state.grid, g), matter.radiation, Species.RADIATION,
        stage="manufactured-tangent-perturbed",
    )
    g.alpha = np.asarray(kernel._solve_lapse(state.grid, g, scalars, matter)[0], dtype=float)
    g.beta.fill(0.0)
    g.B.fill(0.0)
    g.assert_finite_positive()
    return ProductionState(
        grid=state.grid, geometry=g, scalars=scalars, matter=matter,
        t=state.t + eps, tau=state.tau, e_folds=state.e_folds,
        cycle=state.cycle, history=[], handoffs=[],
    )


def regional_peaks(arrays: dict[str, np.ndarray], r: np.ndarray) -> dict:
    masks = {
        "first_two_cells": np.arange(r.size) < 2,
        "center_adjacent_cells_0_to_4": np.arange(r.size) < min(5, r.size),
        "fixed_inner_annulus_2_le_r_lt_5": (r >= 2.0) & (r < 5.0),
        "full_domain": np.ones(r.size, dtype=bool),
    }
    out = {}
    for region, mask in masks.items():
        out[region] = {}
        for name, values in arrays.items():
            a = finite(name, values)
            ids = np.flatnonzero(mask)
            i = int(ids[np.argmax(np.abs(a[ids]))]) if ids.size else 0
            out[region][name] = {"i": i, "r": float(r[i]), "signed": float(a[i]), "abs": float(abs(a[i]))}
    return out


def cell_table(arrays: dict[str, np.ndarray], r: np.ndarray, indices: tuple[int, ...] = (0, 1)) -> dict:
    return {
        str(i): {"r": float(r[i]), **{name: float(values[i]) for name, values in arrays.items()}}
        for i in indices if 0 <= i < r.size
    }


def state_witness(kernel, state) -> dict:
    r = np.asarray(state.grid.centers, dtype=float)
    arrays, closure = momentum_terms(state.grid, state.geometry, state.scalars, state.matter)
    constraints = vacuum.constraints(state.grid, state.geometry)
    total = assemble_total_stress_energy(state.grid, state.geometry, state.scalars, state.matter)
    H = np.asarray(constraints["hamiltonian"], dtype=float) - 16.0 * math.pi * np.asarray(total.rho, dtype=float)
    M = arrays["momentum_residual_direct"]
    C = np.asarray(state.geometry.Lambda, dtype=float) - metric_lambda(state.grid, state.geometry)
    det = np.asarray(constraints["determinant"], dtype=float)
    return {
        "tangent_state_t": float(state.t),
        "constraint_residuals": {
            "H": profile_summary(H, r), "M": profile_summary(M, r),
            "C_lambda": profile_summary(C, r), "determinant": profile_summary(det, r),
        },
        "decomposition_closure": closure,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("results/manufactured-momentum-tangent"))
    args = parser.parse_args()
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    kernel = V55ProductionKernel()
    state = make_manufactured_state(kernel)
    grid = state.grid
    r = np.asarray(grid.centers, dtype=float)
    analytic = analytic_derivative_audit(grid)
    base_arrays, base_closure = momentum_terms(grid, state.geometry, state.scalars, state.matter)
    local_s = local_s_budget(grid, state.geometry, state.scalars)

    srhs, mrhs, _ = kernel._rhs(state)
    gterms = adapter.geometry_stage_terms(grid, state.geometry, state.scalars, state.matter, lambda_m=LAMBDA_M)
    geometry_rates = {
        "a": np.asarray(gterms["explicit"]["a"], dtype=float),
        "b": np.asarray(gterms["explicit"]["b"], dtype=float),
        "X": np.asarray(gterms["explicit"]["X"], dtype=float),
        "Aa": np.asarray(gterms["primary_l2"]["Aa"] + gterms["primary_l3"]["Aa"], dtype=float),
        "K": np.asarray(gterms["primary_l2"]["K"] + gterms["primary_l3"]["K"], dtype=float),
        "Lambda": np.asarray(gterms["lambda_l2"] + gterms["lambda_l3"], dtype=float),
    }
    tangent_by_eps = {}
    profile_arrays = {}
    for eps in EPSILONS:
        plus = tangent_state(kernel, state, srhs, mrhs, geometry_rates, +eps)
        minus = tangent_state(kernel, state, srhs, mrhs, geometry_rates, -eps)
        plus_arrays, plus_closure = momentum_terms(grid, plus.geometry, plus.scalars, plus.matter)
        minus_arrays, minus_closure = momentum_terms(grid, minus.geometry, minus.scalars, minus.matter)
        rates = {name: (plus_arrays[name] - minus_arrays[name]) / (2.0 * eps) for name in base_arrays}
        rates = {name: finite(f"tangent {name}", arr) for name, arr in rates.items()}
        key = format(eps, ".1e")
        tangent_by_eps[key] = rates
        for name, arr in rates.items():
            profile_arrays[f"eps_{key}_{name}"] = arr
        tangent_by_eps[key]["_closure"] = {
            "geometry_decomposition_rate_closure": float(np.max(np.abs(rates["momentum_geometric_vendor"] - sum((rates[f"geo_{k}"] for k in ("Dr_Aa", "minus_2third_Dr_K", "6_Aa_chi_r", "Aa_radial_2_over_r", "Aa_b_gradient")), np.zeros_like(r))))),
            "source_decomposition_rate_closure": float(np.max(np.abs(rates["momentum_source_vendor"] - sum((rates[k] for k in ("src_minus_8piJ_S", "src_minus_8piJ_D", "src_minus_8piJ_phi", "src_minus_8piJ_dark_matter", "src_minus_8piJ_baryons", "src_minus_8piJ_radiation")), np.zeros_like(r))))),
            "direct_M_tangent_vs_rebuilt_sum_max_abs": float(np.max(np.abs(rates["momentum_residual_direct"] - rates["momentum_residual_rebuilt"]))),
            "plus_state_decomposition": plus_closure,
            "minus_state_decomposition": minus_closure,
        }

    fine_key = format(EPSILONS[-1], ".1e")
    coarse_key = format(EPSILONS[0], ".1e")
    fine = {k: v for k, v in tangent_by_eps[fine_key].items() if not k.startswith("_")}
    coarse = {k: v for k, v in tangent_by_eps[coarse_key].items() if not k.startswith("_")}
    delta_epsilon = {name: finite(name, fine[name] - coarse[name]) for name in fine}
    for name, arr in delta_epsilon.items():
        profile_arrays[f"epsilon_delta_{name}"] = arr

    # Explicitly expose what the actual production geometry RHS does to the
    # local-S sector. This is attribution data, not an assertion that an S-only
    # term should independently preserve the total constraint.
    local_s["production_Aa_RHS_total"] = np.asarray(gterms["primary_l2"]["Aa"] + gterms["primary_l3"]["Aa"], dtype=float)
    local_s["production_K_RHS_total"] = np.asarray(gterms["primary_l2"]["K"] + gterms["primary_l3"]["K"], dtype=float)
    for name, arr in local_s.items():
        profile_arrays[name] = finite(name, arr)

    base_summary = state_witness(kernel, state)
    metric_det_residual = np.asarray(state.geometry.a, dtype=float) * np.asarray(state.geometry.b, dtype=float)**2 - 1.0
    metric_regularity = {
        "first_two_cell_abs_a_minus_b": [float(abs(state.geometry.a[i] - state.geometry.b[i])) for i in (0, 1)],
        "first_two_cell_abs_a_minus_b_over_r2": [float(abs(state.geometry.a[i] - state.geometry.b[i]) / r[i]**2) for i in (0, 1)],
        "max_abs_determinant_residual": float(np.max(np.abs(metric_det_residual))),
        "finite": bool(np.all(np.isfinite(metric_det_residual))),
    }
    if not metric_regularity["finite"] or metric_regularity["max_abs_determinant_residual"] > 1.0e-12:
        raise AssertionError(f"manufactured metric regularity/determinant contract failed: {metric_regularity}")
    tangent_summary = {
        "finite_difference_epsilons": list(EPSILONS),
        "finer_epsilon": fine_key,
        "regional_signed_peak_budgets": regional_peaks({
            name: fine[name] for name in (
                "geo_Dr_Aa", "geo_minus_2third_Dr_K", "geo_6_Aa_chi_r",
                "geo_Aa_radial_2_over_r", "geo_Aa_b_gradient",
                "src_minus_8piJ_S", "src_minus_8piJ_D", "src_minus_8piJ_phi",
                "src_minus_8piJ_dark_matter", "src_minus_8piJ_baryons",
                "src_minus_8piJ_radiation", "momentum_residual_direct",
            )
        }, r),
        "first_two_cells_signed_tangent_terms": cell_table({
            name: fine[name] for name in (
                "geo_Dr_Aa", "geo_minus_2third_Dr_K", "geo_6_Aa_chi_r",
                "geo_Aa_radial_2_over_r", "geo_Aa_b_gradient",
                "src_minus_8piJ_S", "src_minus_8piJ_D", "src_minus_8piJ_phi",
                "src_minus_8piJ_dark_matter", "src_minus_8piJ_baryons",
                "src_minus_8piJ_radiation", "momentum_residual_direct",
            )
        }, r, (0, 1)),
        "fixed_physical_radii_signed_samples": {},
        "termwise_epsilon_sensitivity_max_abs": {name: float(np.max(np.abs(delta_epsilon[name]))) for name in fine},
        "termwise_epsilon_sensitivity_rms": {name: float(np.sqrt(np.mean(delta_epsilon[name]**2))) for name in fine},
        "closure_by_epsilon": {key: value["_closure"] for key, value in tangent_by_eps.items()},
    }
    for radius in FIXED_RADII:
        idx = int(np.argmin(np.abs(r - radius)))
        tangent_summary["fixed_physical_radii_signed_samples"][str(radius)] = {
            "i": idx, "r_cell_center": float(r[idx]),
            **{name: float(fine[name][idx]) for name in (
                "geo_Dr_Aa", "geo_minus_2third_Dr_K", "geo_6_Aa_chi_r",
                "geo_Aa_radial_2_over_r", "geo_Aa_b_gradient",
                "src_minus_8piJ_S", "src_minus_8piJ_D", "src_minus_8piJ_phi",
                "src_minus_8piJ_dark_matter", "src_minus_8piJ_baryons",
                "src_minus_8piJ_radiation", "momentum_residual_direct",
            )},
        }

    report = {
        "schema": "v55_manufactured_momentum_tangent_parity_audit_v1",
        "status": "completed",
        "diagnostic_only": True,
        "production_physics_changed": False,
        "production_steps_run": 0,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "configuration": {"resolution": RESOLUTION, "r_max": R_MAX, "dr": float(grid.dr), "cell_centers": "r=(i+1/2)dr"},
        "manufactured_profile": {
            "constraint_solved": False,
            "description": "smooth even regular scalar/metric profiles; a=exp(f), b=exp(-f/2), f=2e-3*(r/6)^2*exp(-(r/6)^2), so a-b=O(r^2) and a*b^2=1; Aa proportional to r^2; Lambda recomputed from the metric connection; zero shift",
            "radiation": "conservative radiation retained from the initialized seed and validated on the deformed metric",
            "amplitude_definitions": {"f": "2e-3 (r/6)^2 exp(-(r/6)^2)", "a": "exp(f)", "b": "exp(-f/2)", "X_minus_1": "1.5e-3 exp(-(r/6)^2)", "Aa": "2e-6 r^2 exp(-(r/6)^2)", "K_plus_0.02": "5e-4 exp(-(r/6)^2)"},
        },
        "identities": {
            "vendor_momentum_constraint": "D_r Aa - (2/3)D_r K + 6 Aa chi_r + (Aa-Ab)(2/r + D_r b/b), chi_r=-D_r X/(2X), Ab=-Aa/2",
            "production_total_momentum_residual": "vendor momentum geometry - 8*pi*assembled total J",
            "source_parity": {"S": "even", "PS": "even", "D": "even", "PD": "even", "phi": "even", "Pi": "even", "a": "even", "b": "even", "X": "even", "K": "even", "Aa": "even O(r^2)", "Lambda": "odd from metric connection"},
            "derivative_operator": "Pinned vendor SphericalCellGrid.cell_derivative_fourth, same parity argument and ghost construction as production",
        },
        "analytic_parity_controls": analytic,
        "manufactured_metric_regularity": metric_regularity,
        "base_state": base_summary,
        "base_decomposition_closure": base_closure,
        "local_S_budget": {
            "signed_profiles": {name: profile_summary(value, r) for name, value in local_s.items()},
            "meaning": "The local-S current and its stress contributions are shown separately; the total geometric RHS couples to all stress sectors, so no one-sector compensation is assumed.",
        },
        "semi_discrete_tangent": tangent_summary,
        "interpretation_limits": [
            "The manufactured profile is not constraint-solved, so a nonzero M tangent is not by itself a conservation defect.",
            "The algebraic closures test the implementation's discrete split and source assembly, not the continuum Einstein equations.",
            "A source-law defect is not claimed unless the manufactured identity separates it from operator/parity mismatch.",
            "No production timestep or trajectory was advanced, and no production source was modified.",
        ],
    }
    report_path = outdir / "report.json"
    report_path.write_text(json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(outdir / "profiles.npz", r=r, **profile_arrays)
    print("MANUFACTURED_MOMENTUM_TANGENT_AUDIT " + json.dumps({
        "status": report["status"],
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "resolution": RESOLUTION,
        "analytic_parity_pass": analytic["pass"],
        "base_closure": base_closure,
        "closure_by_epsilon": tangent_summary["closure_by_epsilon"],
        "first_two_cells": tangent_summary["first_two_cells_signed_tangent_terms"],
        "production_steps_run": 0,
    }, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
