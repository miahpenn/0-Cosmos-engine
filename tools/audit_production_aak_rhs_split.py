#!/usr/bin/env python3
"""Frozen-state signed Aa/K RHS and matter-current tangent audit.

Inputs: pinned-production baseline and discrete-consistent initial-metric
checkpoints at t=0 and saved t=2.16, plus frozen t=0 states built with the
uniform-fourth-order radial quadrature candidate. No state is evolved.

The tool reports each signed contribution to the production Aa and K RHS,
the regular-centre quotients, and a directional finite-difference tangent
of geometric momentum terms and each scalar/fluid current using the actual
semi-discrete production RHS. The CMC lapse is re-solved on both perturbed
states. This is a numerical compatibility diagnostic, not an independent
proof of continuum stress-energy conservation.
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np

SOURCE = Path("source").resolve()
sys.path.insert(0, str(SOURCE))

from engine.production_kernel import V55ProductionKernel, ProductionState, LAMBDA_M
from engine.scalar_system import ScalarFields, cosmos_potential
from engine.matter_system import ConservedSpecies, Species, project_species
from engine.v55_matter import (
    V55MatterState, metric_slice_from_q,
    RHO_DM_PRESENT, RHO_B_PRESENT, RHO_R_PRESENT, initialize_from_archive,
)
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
AB_RUN_ID = 38284401548
AB_ROOT = Path("production-ab-artifact/runs/discrete-consistent-initial-production-test")
RESOLUTIONS = (60, 80, 120)
UNIFORM_T0_RESOLUTIONS = (40, 60, 80, 120, 160)
EPSILONS = (2.0e-5, 1.0e-5)
FIXED_RADII = (2.0, 3.0, 4.0, 5.0)
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES_PREFIX = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED_FIELDS = ("rest", "energy_t", "momentum_r")
ORIGIN_WEIGHTS = (401.0/720.0, -31.0/480.0, 11.0/1440.0)
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
R_MAX = 80.0


def d1(grid, values, parity=+1):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(values, dtype=float), parity=parity), dtype=float)


def summarize(arr, r, mask):
    vals = np.asarray(arr, dtype=float)
    ids = np.flatnonzero(np.asarray(mask, dtype=bool) & np.isfinite(vals))
    if not ids.size:
        return {"count": 0, "max_abs": None, "peak_signed": None, "peak_r": None, "rms": None}
    i = int(ids[np.argmax(np.abs(vals[ids]))])
    return {
        "count": int(ids.size),
        "max_abs": float(abs(vals[i])),
        "peak_signed": float(vals[i]),
        "peak_r": float(r[i]),
        "rms": float(np.sqrt(np.mean(vals[ids]**2))),
    }


def high_order_increments(src, dr):
    """Four-point cubic Lagrange integration over each centre interval."""
    n = len(src)
    if n < 4:
        raise ValueError("uniform radial quadrature requires >=4 cells")
    out = np.empty(n-1, dtype=float)
    # First interval: even-reflection ghost; cubic integral weights.
    out[0] = dr * (12.0*src[0] + 13.0*src[1] - src[2]) / 24.0
    for i in range(2, n-1):
        out[i-1] = dr * (-src[i-2] + 13.0*src[i-1] + 13.0*src[i] - src[i+1]) / 24.0
    # Final interval: one-sided cubic integral.
    out[-1] = dr * (src[-4] - 5.0*src[-3] + 19.0*src[-2] + 9.0*src[-1]) / 24.0
    return out


def uniform_fourth_order_initial_state(kernel, n):
    """Reconstruct the frozen uniform-fourth-order t=0 candidate exactly as audited."""
    grid_ops, _, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    r = np.asarray(grid.centers, dtype=float)
    dr = float(grid.dr)
    phi0, pi0 = 0.179055, 0.00341913
    dm0, b0, rho_r = RHO_DM_PRESENT, RHO_B_PRESENT, RHO_R_PRESENT
    rho0 = 0.5*pi0*pi0 + float(cosmos_potential(np.array([phi0]))[0]) + dm0 + b0 + rho_r
    H0 = math.sqrt(max(rho0/3.0, 0.0))
    S = AMPLITUDE*np.exp(-(r/WIDTH)**2)
    D = D_AMPLITUDE*np.exp(-(r/WIDTH)**2)
    PS = np.zeros_like(r)
    PD = D.copy()
    phi = np.full(n, phi0)
    Pi = np.full(n, pi0)
    Kt = np.full(n, -H0)
    Sp = d1(grid, S, +1)
    Dp = d1(grid, D, +1)
    Ktp = d1(grid, Kt, +1)
    j = -(PD*Dp)
    B = np.ones_like(r)
    for _ in range(6):
        rhoS = 0.5*(PS**2 + B*Sp**2) + 0.5*S**2
        rhoD = 0.5*(PD**2 + B*Dp**2) - 0.5*D**2
        rhoPhi = (0.5*Pi**2 + cosmos_potential(phi))/(8.0*math.pi)
        rhoMatter = (dm0+b0+rho_r)/(8.0*math.pi)
        rho = rhoS + rhoD + rhoPhi + rhoMatter
        Kr = Kt + r*Ktp + 4.0*math.pi*r*j
        src = 1.0 - 8.0*math.pi*r**2*rho + 2.0*r**2*Kr*Kt + r**2*Kt**2
        half_cell = dr*sum(ORIGIN_WEIGHTS[k]*src[k] for k in range(3))
        increments = high_order_increments(src, dr)
        integral = np.cumsum(np.r_[half_cell, increments])
        B = integral/r
        if not np.all(np.isfinite(B)) or np.any(B <= 0.0):
            raise FloatingPointError("uniform fourth-order integral produced invalid B")
    A = 1.0/np.sqrt(B)
    Kr = Kt + r*Ktp + 4.0*math.pi*r*j
    K = Kr + 2.0*Kt
    Aa = Kr - K/3.0
    a = A**(4.0/3.0)
    b = A**(-2.0/3.0)
    X = A**(-1.0/3.0)
    Lambda = d1(grid, a, +1)/(2.0*a*a) - d1(grid, b, +1)/(a*b) + 2.0/r*(1.0/b-1.0/a)
    geom = vacuum.flat_state(grid)
    geom.a[:], geom.b[:], geom.X[:] = a, b, X
    geom.alpha[:] = 1.0
    geom.beta[:] = 0.0
    geom.Aa[:], geom.K[:], geom.Lambda[:] = Aa, K, Lambda
    geom.B[:] = 0.75*Lambda
    scalars = ScalarFields(S, PS, D, PD, phi, Pi)
    matter = initialize_from_archive(grid, geom, include_radiation=True)
    geom.alpha = kernel._solve_lapse(grid, geom, scalars, matter)[0]
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    return ProductionState(grid=grid, geometry=geom, scalars=scalars, matter=matter, t=0., tau=0., e_folds=0.)


def restore_checkpoint(path, kernel):
    with np.load(path, allow_pickle=False) as z:
        data = {k: np.asarray(z[k]).copy() for k in z.files}
    n = int(data["r"].size)
    dr = float(data["r"][1]-data["r"][0])
    state = kernel.initialize(resolution=n, r_max=n*dr, amplitude=AMPLITUDE, width=WIDTH,
                              D_amplitude=D_AMPLITUDE, include_radiation=True)
    if not np.allclose(state.grid.centers, data["r"], rtol=0.0, atol=1.0e-12*max(1.0, n*dr)):
        raise ValueError(f"{path}: checkpoint grid mismatch")
    for name in GEOMETRY_FIELDS:
        setattr(state.geometry, name, np.asarray(data[name], dtype=float).copy())
    state.scalars = ScalarFields(*(np.asarray(data[name], dtype=float).copy() for name in SCALAR_FIELDS))
    species = {}
    for name, prefix in SPECIES_PREFIX.items():
        species[name] = ConservedSpecies(*(
            np.asarray(data[f"{prefix}_{field}"], dtype=float).copy() for field in CONSERVED_FIELDS
        ))
    state.matter = V55MatterState(
        dark_matter=species["dark_matter"], baryons=species["baryons"], radiation=species["radiation"]
    )
    state.t, state.tau, state.e_folds = float(data["t"]), float(data["tau"]), float(data["e_folds"])
    state.history, state.handoffs = [], []
    return state


def source_momentum_terms(grid, geom, scalars, matter):
    r = np.asarray(grid.centers, dtype=float)
    metric = metric_slice_from_q(grid, geom)
    Aa = np.asarray(geom.Aa, dtype=float)
    Ab = -0.5*Aa
    K = np.asarray(geom.K, dtype=float)
    X = np.asarray(geom.X, dtype=float)
    b = np.asarray(geom.b, dtype=float)
    total = assemble_total_stress_energy(grid, geom, scalars, matter)
    raw = vacuum.constraints(grid, geom)

    jS = -np.asarray(scalars.PS, dtype=float)*d1(grid, scalars.S, +1)
    jD = -np.asarray(scalars.PD, dtype=float)*d1(grid, scalars.D, +1)
    jPhi = -np.asarray(scalars.Pi, dtype=float)*d1(grid, scalars.phi, +1)/(8.0*math.pi)
    jDM = np.asarray(project_species(metric, matter.dark_matter, Species.DARK_MATTER)["j"], dtype=float)
    jB = np.asarray(project_species(metric, matter.baryons, Species.BARYON)["j"], dtype=float)
    jR = np.asarray(project_species(metric, matter.radiation, Species.RADIATION)["j"], dtype=float)

    tAa = d1(grid, Aa, +1)
    tK = -(2.0/3.0)*d1(grid, K, +1)
    tChi = 6.0*Aa*(-0.5*d1(grid, X, +1)/X)
    tB = (Aa-Ab)*(2.0/r+d1(grid, b, +1)/b)
    mGeom = np.asarray(raw["momentum"], dtype=float)
    residual = mGeom - 8.0*math.pi*np.asarray(total.j, dtype=float)
    signedSources = {
        "M_source_S": -8.0*math.pi*jS,
        "M_source_D": -8.0*math.pi*jD,
        "M_source_phi": -8.0*math.pi*jPhi,
        "M_source_dark_matter": -8.0*math.pi*jDM,
        "M_source_baryons": -8.0*math.pi*jB,
        "M_source_radiation": -8.0*math.pi*jR,
    }
    builtSources = sum(signedSources.values())
    return {
        "M_term_Dr_Aa": tAa,
        "M_term_minus_2third_Dr_K": tK,
        "M_term_6Aa_chi_r": tChi,
        "M_term_radial_b": tB,
        "M_geometry": mGeom,
        **signedSources,
        "M_source_total": -8.0*math.pi*np.asarray(total.j, dtype=float),
        "M_residual": residual,
        "M_geometry_closure": mGeom-(tAa+tK+tChi+tB),
        "M_source_closure": builtSources+8.0*math.pi*np.asarray(total.j, dtype=float),
    }


def production_rhs_terms(state):
    grid, geom = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=float)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    g = vacuum.geometry_terms(grid, geom)
    stages = adapter.geometry_stage_terms(grid, geom, state.scalars, state.matter, lambda_m=LAMBDA_M)
    pureL3 = vacuum.primary_l3_rhs(grid, geom)
    AaMatterExpected = -(16.0*math.pi/3.0)*np.asarray(geom.alpha)*(np.asarray(total.pr)-np.asarray(total.pt))
    AaAdvExtrinsic = np.asarray(geom.beta)*d1(grid, geom.Aa, +1) + np.asarray(geom.alpha)*np.asarray(geom.K)*np.asarray(geom.Aa)
    AaKO = np.asarray(pureL3["Aa"])-AaAdvExtrinsic
    AaHessTF = -(np.asarray(g["DrDr_alpha"])-np.asarray(g["lap_alpha"])/3.0)
    AaRicciTF = np.asarray(geom.alpha)*(np.asarray(g["Rr_r"])-np.asarray(g["R"])/3.0)
    AaMatterActual = np.asarray(stages["primary_l3"]["Aa"])-np.asarray(pureL3["Aa"])
    AaRhs = AaHessTF + AaRicciTF + AaAdvExtrinsic + AaKO + AaMatterActual
    AaRhsActual = np.asarray(stages["primary_l2"]["Aa"])+np.asarray(stages["primary_l3"]["Aa"])

    Ab = -0.5*np.asarray(geom.Aa)
    K_lapse = -np.asarray(g["lap_alpha"])
    KAdvExtrinsic = np.asarray(geom.beta)*d1(grid, geom.K, +1) + np.asarray(geom.alpha)*(
        np.asarray(geom.Aa)**2 + 2.0*Ab**2 + np.asarray(geom.K)**2/3.0
    )
    KKO = np.asarray(pureL3["K"])-KAdvExtrinsic
    KMatterExpected = 4.0*math.pi*np.asarray(geom.alpha)*(
        np.asarray(total.rho)+np.asarray(total.pr)+2.0*np.asarray(total.pt)
    )
    KMatterActual = np.asarray(stages["primary_l3"]["K"])-np.asarray(pureL3["K"])
    KRhs = K_lapse + KAdvExtrinsic + KKO + KMatterActual
    KRhsActual = np.asarray(stages["primary_l2"]["K"])+np.asarray(stages["primary_l3"]["K"])
    M = source_momentum_terms(grid, geom, state.scalars, state.matter)

    arrays = {
        "Aa_hessian_trace_free": AaHessTF,
        "Aa_ricci_trace_free": AaRicciTF,
        "Aa_advection_extrinsic": AaAdvExtrinsic,
        "Aa_KO_dissipation": AaKO,
        "Aa_matter_anisotropic_stress": AaMatterActual,
        "Aa_matter_expected_from_total_stress": AaMatterExpected,
        "Aa_rhs_reconstructed": AaRhs,
        "Aa_rhs_production": AaRhsActual,
        "Aa_rhs_closure": AaRhs-AaRhsActual,
        "K_lapse_laplacian": K_lapse,
        "K_advection_extrinsic": KAdvExtrinsic,
        "K_KO_dissipation": KKO,
        "K_matter_source": KMatterActual,
        "K_matter_expected_from_total_stress": KMatterExpected,
        "K_rhs_reconstructed": KRhs,
        "K_rhs_production": KRhsActual,
        "K_rhs_closure": KRhs-KRhsActual,
        **M,
    }
    # Quantify conservation-current derivatives per sector; matter source signs
    # in M_residual are already negative 8*pi times each current.
    if any(not np.all(np.isfinite(v)) for v in arrays.values()):
        bad = [name for name, v in arrays.items() if not np.all(np.isfinite(v))]
        raise FloatingPointError(f"nonfinite RHS/budget terms: {bad}")
    return arrays


def add_species(base, rhs, factor):
    return ConservedSpecies(
        np.asarray(base.rest, dtype=float)+factor*np.asarray(rhs.rest, dtype=float),
        np.asarray(base.energy_t, dtype=float)+factor*np.asarray(rhs.energy_t, dtype=float),
        np.asarray(base.momentum_r, dtype=float)+factor*np.asarray(rhs.momentum_r, dtype=float),
    )


def perturbed_state(kernel, state, scalar_rhs, matter_rhs, gt, factor):
    g = state.geometry.copy()
    explicit, l2, l3 = gt["explicit"], gt["primary_l2"], gt["primary_l3"]
    rates = {
        "a": explicit["a"], "b": explicit["b"], "X": explicit["X"],
        "Aa": l2["Aa"]+l3["Aa"], "K": l2["K"]+l3["K"],
        "Lambda": gt["lambda_l2"]+gt["lambda_l3"],
    }
    for name, rate in rates.items():
        setattr(g, name, np.asarray(getattr(g, name), dtype=float)+factor*np.asarray(rate, dtype=float))
    g.beta.fill(0.0)
    g.B.fill(0.0)
    scalars = ScalarFields(*(
        np.asarray(getattr(state.scalars, name), dtype=float)
        + factor*np.asarray(getattr(scalar_rhs, name), dtype=float) for name in SCALAR_FIELDS
    ))
    species = {name: add_species(getattr(state.matter, name), matter_rhs[name], factor) for name in SPECIES_PREFIX}
    matter = V55MatterState(dark_matter=species["dark_matter"], baryons=species["baryons"], radiation=species["radiation"])
    g.alpha = np.asarray(kernel._solve_lapse(state.grid, g, scalars, matter)[0], dtype=float)
    g.assert_finite_positive()
    return ProductionState(grid=state.grid, geometry=g, scalars=scalars, matter=matter,
                           t=state.t+factor, tau=state.tau, e_folds=state.e_folds,
                           cycle=state.cycle, history=[], handoffs=[])


def audit_state(kernel, state, label, arm, n, scheduled_t=None):
    grid = state.grid
    r = np.asarray(grid.centers, dtype=float)
    radius2 = r*r
    dt_rhs = {
        "Aa_rhs": None, "K_rhs": None
    }
    base = production_rhs_terms(state)
    scalar_rhs, matter_rhs, _ = kernel._rhs(state)
    gt = adapter.geometry_stage_terms(grid, state.geometry, state.scalars, state.matter, lambda_m=LAMBDA_M)
    tangent_eps = {}
    for eps in EPSILONS:
        plus = perturbed_state(kernel, state, scalar_rhs, matter_rhs, gt, +eps)
        minus = perturbed_state(kernel, state, scalar_rhs, matter_rhs, gt, -eps)
        p = production_rhs_terms(plus)
        m = production_rhs_terms(minus)
        tangent_eps[eps] = {key:(p[key]-m[key])/(2.0*eps) for key in base}
    fine, coarse = tangent_eps[EPSILONS[-1]], tangent_eps[EPSILONS[0]]

    masks = {
        "first_three_cells": np.arange(len(r)) < min(3, len(r)),
        "first_two_cells": np.arange(len(r)) < min(2, len(r)),
        "formal_stencil_i_ge_2": np.arange(len(r)) >= 2,
        "admitted_inner_i_ge_2_r_lt_5": (np.arange(len(r)) >= 2) & (r < 5.0),
        "fixed_annulus_2_le_r_lt_5": (r >= 2.0) & (r < 5.0),
        "exterior_r_ge_5": r >= 5.0,
        "full_domain": np.ones_like(r, dtype=bool),
    }
    focus = (
        "Aa_hessian_trace_free", "Aa_ricci_trace_free", "Aa_advection_extrinsic",
        "Aa_KO_dissipation", "Aa_matter_anisotropic_stress", "Aa_rhs_production",
        "Aa_rhs_closure", "K_lapse_laplacian", "K_advection_extrinsic",
        "K_KO_dissipation", "K_matter_source", "K_rhs_production", "K_rhs_closure",
        "M_term_Dr_Aa", "M_term_minus_2third_Dr_K", "M_term_6Aa_chi_r",
        "M_term_radial_b", "M_geometry", "M_source_S", "M_source_D", "M_source_phi",
        "M_source_dark_matter", "M_source_baryons", "M_source_radiation", "M_residual"
    )
    summaries = {}
    for name in focus:
        summaries[name] = {
            "base": {reg:summarize(base[name], r, mask) for reg, mask in masks.items()},
            "rhs_tangent": {reg:summarize(fine[name], r, mask) for reg, mask in masks.items()},
            "epsilon_difference_max_abs": float(np.max(np.abs(fine[name]-coarse[name]))),
        }
    cells = {}
    for i in range(min(6, len(r))):
        cells[str(i)] = {
            "r": float(r[i]),
            "base": {name: float(base[name][i]) for name in focus},
            "rhs_tangent": {name: float(fine[name][i]) for name in focus},
            "coarse_epsilon_tangent": {name: float(coarse[name][i]) for name in focus},
        }
        cells[str(i)]["series_witnesses"] = {
            "r2_Aa_rhs": float(radius2[i]*base["Aa_rhs_production"][i]),
            "Aa_rhs_over_r2": float(base["Aa_rhs_production"][i]/radius2[i]),
            "r2_K_rhs": float(radius2[i]*base["K_rhs_production"][i]),
            "Dr_K_rhs_over_r": float(d1(grid, base["K_rhs_production"], +1)[i]/r[i]),
        }
    fixed = {}
    for rr in FIXED_RADII:
        if r[0] <= rr <= r[-1]:
            fixed[str(rr)] = {
                group: {name:float(np.interp(rr, r, arr)) for name, arr in arrays.items()}
                for group, arrays in (("base", base), ("rhs_tangent", fine))
            }

    closures = {
        "Aa_rebuild_max_abs": float(np.max(np.abs(base["Aa_rhs_closure"]))),
        "K_rebuild_max_abs": float(np.max(np.abs(base["K_rhs_closure"]))),
        "M_geom_rebuild_max_abs": float(np.max(np.abs(base["M_geometry_closure"]))),
        "M_source_rebuild_max_abs": float(np.max(np.abs(base["M_source_closure"]))),
        "tangent_Aa_rebuild_max_abs": float(np.max(np.abs(fine["Aa_rhs_closure"]))),
        "tangent_K_rebuild_max_abs": float(np.max(np.abs(fine["K_rhs_closure"]))),
        "tangent_M_geom_rebuild_max_abs": float(np.max(np.abs(fine["M_geometry_closure"]))),
        "tangent_M_source_rebuild_max_abs": float(np.max(np.abs(fine["M_source_closure"]))),
        "max_epsilon_difference_by_term": {
            name: float(np.max(np.abs(fine[name]-coarse[name]))) for name in focus
        }
    }
    return {
        "label": label, "arm": arm, "resolution": int(n),
        "actual_t": float(state.t), "scheduled_t": scheduled_t,
        "tau": float(state.tau), "dr": float(grid.dr), "r_max": float(len(r)*grid.dr),
        "epsilons": list(EPSILONS), "terms": summaries, "first_cells": cells,
        "fixed_radius": fixed, "closures": closures,
        "interpretation_limit": (
            "RHS tangent uses the actual production semi-discrete vector field and a re-solved CMC lapse. "
            "Arithmetic closure is not an independent covariant conservation proof. No state is evolved."
        )
    }


def main():
    kernel = V55ProductionKernel()
    states = []
    # Baseline and discrete-consistent metric A/B at t=0 and saved t=2.16.
    for n in RESOLUTIONS:
        for arm in ("baseline", "candidate"):
            for target, filename in ((0.0, "checkpoint_t00000.000000.npz"), (2.16, "checkpoint_t00002.160000.npz")):
                path = AB_ROOT/arm/f"N{n}"/filename
                if not path.is_file():
                    raise FileNotFoundError(f"required paired checkpoint missing: {path}")
                states.append((restore_checkpoint(path, kernel), f"N{n}_{arm}_T{target:g}", arm, n, target))
    # The uniform fourth-order candidate is a frozen t=0 control only; it is not evolved.
    for n in UNIFORM_T0_RESOLUTIONS:
        states.append((uniform_fourth_order_initial_state(kernel, n), f"N{n}_uniform_fourth_order_T0",
                       "uniform_fourth_order_intervals", n, 0.0))

    records = [audit_state(kernel, state, label, arm, n, scheduled_t=target)
               for state, label, arm, n, target in states]
    report = {
        "schema": "production_Aa_K_rhs_and_matter_current_tangent_audit_v1",
        "status": "completed",
        "kind": "frozen_states_only_signed_rhs_decomposition_plus_directional_tangent_no_evolution",
        "production_physics_changed": False,
        "source_commit": SOURCE_COMMIT, "vendor_commit": VENDOR_COMMIT,
        "input_ab_run_id": AB_RUN_ID,
        "state_count": len(records),
        "states": records,
        "guardrails": [
            "The uniform-fourth-order initializer is measured only at t=0; no trajectory was continued from it.",
            "Baseline and discrete-consistent initial-metric arms use frozen t=0 and t=2.16 checkpoints from the paired A/B campaign.",
            "All signed terms are evaluated using the exact frozen production source and pinned vendor operators.",
            "No physics equation, coefficient, gauge, boundary, radiation admissibility rule, projection, damping, floor or clamp is changed.",
            "Term closure is an accounting check, not proof of a continuum identity or independent conservation law.",
        ]
    }
    out = Path("results/production-Aa-K-centre-source-tangent")
    out.mkdir(parents=True, exist_ok=True)
    (out/"report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+"\n", encoding="utf-8")
    print("PRODUCTION_AA_K_CENTRE_SOURCE_TANGENT "+json.dumps({
        "status": report["status"], "state_count": len(records),
        "cases": [{
            "label": c["label"], "N": c["resolution"], "arm": c["arm"], "t": c["actual_t"],
            "Aa_firstcell": c["first_cells"]["0"]["base"]["Aa_rhs_production"],
            "Aa_firstcell_over_r2": c["first_cells"]["0"]["series_witnesses"]["Aa_rhs_over_r2"],
            "Aa_hessian_firstcell": c["first_cells"]["0"]["base"]["Aa_hessian_trace_free"],
            "Aa_ricci_firstcell": c["first_cells"]["0"]["base"]["Aa_ricci_trace_free"],
            "Aa_matter_firstcell": c["first_cells"]["0"]["base"]["Aa_matter_anisotropic_stress"],
            "K_rhs_firstcell": c["first_cells"]["0"]["base"]["K_rhs_production"],
            "Dr_Kdot_over_r_firstcell": c["first_cells"]["0"]["series_witnesses"]["Dr_K_rhs_over_r"],
            "Mdot_i2": c["first_cells"].get("2",{}).get("rhs_tangent",{}).get("M_residual"),
            "closures": c["closures"]
        } for c in records],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
