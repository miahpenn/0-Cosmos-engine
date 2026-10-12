#!/usr/bin/env python3
"""Frozen-state multi-constraint tangent audit on matched initial-data A/B states.

Reads baseline/candidate t=0 checkpoints from the discrete-consistent initial
metric experiment. Evaluates the actual semi-discrete production RHS, perturbs
all evolved sectors +/- epsilon, re-solves the production CMC lapse, and
finite-differences signed H, M and Lambda-connection terms. One standard
production step is a local comparator; no trajectory is continued.
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path("source").resolve()))
from engine.production_kernel import V55ProductionKernel, ProductionState, LAMBDA_M
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()
SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
AB_RUN_ID = 38284401548
AB_ROOT = Path("production-ab-artifact/runs/discrete-consistent-initial-production-test")
RESOLUTIONS = (60, 80, 120)
ARMS = ("baseline", "candidate")
EPSILONS = (2.0e-5, 1.0e-5)
FIXED_RADII = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0)
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES_PREFIX = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED_FIELDS = ("rest", "energy_t", "momentum_r")


def d1(grid, x, parity=+1):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(x, dtype=float), parity=parity), dtype=float)


def metric_connection(grid, g):
    r = np.asarray(grid.centers, dtype=float)
    a, b = np.asarray(g.a, dtype=float), np.asarray(g.b, dtype=float)
    return d1(grid, a)/(2*a*a) - d1(grid, b)/(a*b) + 2/r*(1/b - 1/a)


def terms(state):
    grid, g, s, matter = state.grid, state.geometry, state.scalars, state.matter
    r = np.asarray(grid.centers, dtype=float)
    total = assemble_total_stress_energy(grid, g, s, matter)
    raw = vacuum.constraints(grid, g)

    Hgeom = np.asarray(raw["hamiltonian"], dtype=float)
    Hscalar = -16*np.pi*np.asarray(total.scalar_rho, dtype=float)
    Hfluid = -16*np.pi*np.asarray(total.fluid_rho, dtype=float)
    H = Hgeom - 16*np.pi*np.asarray(total.rho, dtype=float)

    Aa, Ab = np.asarray(g.Aa, dtype=float), -0.5*np.asarray(g.Aa, dtype=float)
    K, X, b = np.asarray(g.K, dtype=float), np.asarray(g.X, dtype=float), np.asarray(g.b, dtype=float)
    MAa = d1(grid, Aa, +1)
    MK = -(2/3)*d1(grid, K, +1)
    Mchi = 6*Aa*(-0.5*d1(grid, X, +1)/X)
    Mb = (Aa-Ab)*(2/r + d1(grid, b, +1)/b)
    Mgeom = np.asarray(raw["momentum"], dtype=float)
    Mscalar = -8*np.pi*np.asarray(total.scalar_j, dtype=float)
    Mfluid = -8*np.pi*np.asarray(total.fluid_j, dtype=float)
    M = Mgeom - 8*np.pi*np.asarray(total.j, dtype=float)

    lam = np.asarray(g.Lambda, dtype=float)
    lam_metric = metric_connection(grid, g)
    C = lam - lam_metric
    det = np.asarray(g.a, dtype=float)*np.asarray(g.b, dtype=float)**2 - 1.0
    out = {
        "H_geometry": Hgeom, "H_source_scalar": Hscalar, "H_source_fluid": Hfluid, "H_residual": H,
        "M_term_Dr_Aa": MAa, "M_term_minus_2third_Dr_K": MK, "M_term_6Aa_chi_r": Mchi,
        "M_term_radial_b": Mb, "M_geometry": Mgeom,
        "M_source_scalar": Mscalar, "M_source_fluid": Mfluid, "M_residual": M,
        "Lambda_evolved": lam, "Lambda_metric": lam_metric, "C_lambda": C, "determinant": det,
    }
    out["H_reconstruction_closure"] = Hgeom + Hscalar + Hfluid - H
    out["M_geometry_reconstruction_closure"] = MAa + MK + Mchi + Mb - Mgeom
    out["M_reconstruction_closure"] = Mgeom + Mscalar + Mfluid - M
    source_closure = {
        "energy_scalar_plus_fluid_minus_total_max_abs": float(np.max(np.abs(total.scalar_rho + total.fluid_rho - total.rho))),
        "current_scalar_plus_fluid_minus_total_max_abs": float(np.max(np.abs(total.scalar_j + total.fluid_j - total.j))),
    }
    return out, source_closure


def restore_checkpoint(path, kernel):
    with np.load(path, allow_pickle=False) as z:
        d = {k: np.asarray(z[k]).copy() for k in z.files}
    n, dr = int(d["r"].size), float(d["r"][1]-d["r"][0])
    state = kernel.initialize(resolution=n, r_max=n*dr, amplitude=.01, width=7.,
                              D_amplitude=1.e-10, include_radiation=True)
    if not np.allclose(state.grid.centers, d["r"], rtol=0., atol=1.e-12*max(1., n*dr)):
        raise ValueError(f"{path}: grid mismatch")
    for name in GEOMETRY_FIELDS:
        setattr(state.geometry, name, np.asarray(d[name], dtype=float).copy())
    state.scalars = ScalarFields(*(np.asarray(d[name], dtype=float).copy() for name in SCALAR_FIELDS))
    pieces = {}
    for name, prefix in SPECIES_PREFIX.items():
        pieces[name] = ConservedSpecies(*(np.asarray(d[f"{prefix}_{f}"], dtype=float).copy() for f in CONSERVED_FIELDS))
    state.matter = V55MatterState(dark_matter=pieces["dark_matter"], baryons=pieces["baryons"],
                                  radiation=pieces["radiation"])
    state.t, state.tau, state.e_folds = float(d["t"]), float(d["tau"]), float(d["e_folds"])
    state.history, state.handoffs = [], []
    return state


def add_species(base, rhs, factor):
    return ConservedSpecies(base.rest + factor*np.asarray(rhs.rest, dtype=float),
                            base.energy_t + factor*np.asarray(rhs.energy_t, dtype=float),
                            base.momentum_r + factor*np.asarray(rhs.momentum_r, dtype=float))


def perturbed_state(kernel, state, srhs, mrhs, gt, factor):
    g = state.geometry.copy()
    exp, l2, l3 = gt["explicit"], gt["primary_l2"], gt["primary_l3"]
    rates = {"a": exp["a"], "b": exp["b"], "X": exp["X"],
             "Aa": l2["Aa"]+l3["Aa"], "K": l2["K"]+l3["K"],
             "Lambda": gt["lambda_l2"]+gt["lambda_l3"]}
    for name, rate in rates.items():
        setattr(g, name, np.asarray(getattr(g, name), dtype=float) + factor*np.asarray(rate, dtype=float))
    g.beta.fill(0.0); g.B.fill(0.0)
    s = ScalarFields(*(np.asarray(getattr(state.scalars, name), dtype=float)
                       + factor*np.asarray(getattr(srhs, name), dtype=float) for name in SCALAR_FIELDS))
    parts = {name: add_species(getattr(state.matter, name), mrhs[name], factor) for name in SPECIES_PREFIX}
    matter = V55MatterState(dark_matter=parts["dark_matter"], baryons=parts["baryons"],
                            radiation=parts["radiation"])
    g.alpha = np.asarray(kernel._solve_lapse(state.grid, g, s, matter)[0], dtype=float)
    g.assert_finite_positive()
    return ProductionState(grid=state.grid, geometry=g, scalars=s, matter=matter,
                           t=state.t+factor, tau=state.tau, e_folds=state.e_folds,
                           cycle=state.cycle, history=[], handoffs=[])


def summary(arr, r, mask):
    ids = np.flatnonzero(np.asarray(mask, dtype=bool) & np.isfinite(arr))
    if not ids.size:
        return {"count": 0, "max_abs": None, "peak_signed": None, "peak_r": None, "rms": None}
    i = int(ids[np.argmax(np.abs(arr[ids]))])
    return {"count": int(ids.size), "max_abs": float(abs(arr[i])), "peak_signed": float(arr[i]),
            "peak_r": float(r[i]), "rms": float(np.sqrt(np.mean(np.asarray(arr[ids], dtype=float)**2)))}


def audit(kernel, state, arm, n):
    grid, r = state.grid, np.asarray(state.grid.centers, dtype=float)
    dt = .03*float(grid.dr)
    base, source_closure = terms(state)
    srhs, mrhs, _ = kernel._rhs(state)
    gt = adapter.geometry_stage_terms(grid, state.geometry, state.scalars, state.matter, lambda_m=LAMBDA_M)
    tangent_by_eps = {}
    for eps in EPSILONS:
        plus = perturbed_state(kernel, state, srhs, mrhs, gt, +eps)
        minus = perturbed_state(kernel, state, srhs, mrhs, gt, -eps)
        p, _ = terms(plus); m, _ = terms(minus)
        tangent_by_eps[eps] = {key: (p[key]-m[key])/(2*eps) for key in base}

    fine, coarse = tangent_by_eps[EPSILONS[-1]], tangent_by_eps[EPSILONS[0]]
    before = ProductionState(grid=state.grid, geometry=state.geometry.copy(), scalars=state.scalars,
                             matter=state.matter.copy(), t=state.t, tau=state.tau,
                             e_folds=state.e_folds, cycle=state.cycle, history=[], handoffs=[])
    after = kernel.step(before, dt)
    final, _ = terms(after)
    actual = {key: (final[key]-base[key])/dt for key in base}

    masks = {
        "first_two_cells": np.arange(len(r)) < 2,
        "formal_stencil_i_ge_2": np.arange(len(r)) >= 2,
        "admitted_inner_i_ge_2_r_lt_5": (np.arange(len(r)) >= 2) & (r < 5.),
        "fixed_annulus_2_le_r_lt_5": (r >= 2.) & (r < 5.),
        "exterior_r_ge_5": r >= 5.,
        "full_domain": np.ones_like(r, dtype=bool),
    }
    tracked = ("H_geometry", "H_source_scalar", "H_source_fluid", "H_residual",
               "M_term_Dr_Aa", "M_term_minus_2third_Dr_K", "M_term_6Aa_chi_r",
               "M_term_radial_b", "M_geometry", "M_source_scalar", "M_source_fluid",
               "M_residual", "Lambda_evolved", "Lambda_metric", "C_lambda", "determinant")
    base_summary = {key: {reg: summary(base[key], r, mask) for reg, mask in masks.items()} for key in tracked}
    tangent_summary = {
        key: {
            "rhs_tangent": {reg: summary(fine[key], r, mask) for reg, mask in masks.items()},
            "one_step_delta_over_dt": {reg: summary(actual[key], r, mask) for reg, mask in masks.items()},
            "epsilon_max_abs_difference": float(np.max(np.abs(fine[key]-coarse[key]))),
        } for key in tracked
    }
    cells = {}
    for i in range(min(6, len(r))):
        cells[str(i)] = {"r": float(r[i]),
            "base": {key: float(base[key][i]) for key in tracked},
            "rhs_tangent": {key: float(fine[key][i]) for key in tracked},
            "coarse_epsilon_tangent": {key: float(coarse[key][i]) for key in tracked},
            "one_step_delta_over_dt": {key: float(actual[key][i]) for key in tracked}}
    fixed = {}
    for target in FIXED_RADII:
        if r[0] <= target <= r[-1]:
            fixed[str(target)] = {
                group: {key: float(np.interp(target, r, arr)) for key, arr in dct.items()}
                for group, dct in (("base", base), ("rhs_tangent", fine), ("one_step_delta_over_dt", actual))
            }
    closures = {
        "base_H_reconstruction_max_abs": float(np.max(np.abs(base["H_reconstruction_closure"]))),
        "base_M_geometry_reconstruction_max_abs": float(np.max(np.abs(base["M_geometry_reconstruction_closure"]))),
        "base_M_reconstruction_max_abs": float(np.max(np.abs(base["M_reconstruction_closure"]))),
        "tangent_H_reconstruction_max_abs": float(np.max(np.abs(
            fine["H_geometry"]+fine["H_source_scalar"]+fine["H_source_fluid"]-fine["H_residual"]))),
        "tangent_M_geometry_reconstruction_max_abs": float(np.max(np.abs(
            fine["M_term_Dr_Aa"]+fine["M_term_minus_2third_Dr_K"]+
            fine["M_term_6Aa_chi_r"]+fine["M_term_radial_b"]-fine["M_geometry"]))),
        "tangent_M_reconstruction_max_abs": float(np.max(np.abs(
            fine["M_geometry"]+fine["M_source_scalar"]+fine["M_source_fluid"]-fine["M_residual"]))),
        "source_assembly": source_closure,
    }
    for collection in (base, fine, coarse, actual):
        for key, arr in collection.items():
            if not np.all(np.isfinite(arr)):
                raise FloatingPointError(f"N{n}/{arm}: nonfinite {key}")
    return {"arm": arm, "resolution": n, "actual_t": float(state.t), "tau": float(state.tau),
            "dr": float(grid.dr), "r_max": float(len(r)*grid.dr), "one_step_dt": dt,
            "starting_constraint_summaries": base_summary, "tangent_summaries": tangent_summary,
            "first_cells": cells, "fixed_radius": fixed, "closures": closures,
            "interpretation_limit": (
                "Frozen-state source-sector derivatives and assembly closures are not an independent "
                "proof of covariant stress-energy conservation. This diagnostic tests whether the "
                "near-constraint-satisfied initial slice develops a nonzero production-RHS tangent; "
                "it does not admit long-time constraint propagation."
            )}


def main():
    kernel = V55ProductionKernel()
    cases = []
    for n in RESOLUTIONS:
        for arm in ARMS:
            path = AB_ROOT / arm / f"N{n}" / "checkpoint_t00000.000000.npz"
            if not path.is_file():
                raise FileNotFoundError(f"missing paired initial checkpoint: {path}")
            cases.append(audit(kernel, restore_checkpoint(path, kernel), arm, n))

    paired = []
    for n in RESOLUTIONS:
        baseline = next(c for c in cases if c["resolution"] == n and c["arm"] == "baseline")
        candidate = next(c for c in cases if c["resolution"] == n and c["arm"] == "candidate")
        row = {"resolution": n, "fields": {}}
        for field in ("H_residual", "M_residual", "C_lambda"):
            row["fields"][field] = {}
            for region in ("first_two_cells", "formal_stencil_i_ge_2",
                           "admitted_inner_i_ge_2_r_lt_5", "fixed_annulus_2_le_r_lt_5", "exterior_r_ge_5"):
                b = baseline["tangent_summaries"][field]["rhs_tangent"][region]
                c = candidate["tangent_summaries"][field]["rhs_tangent"][region]
                row["fields"][field][region] = {
                    "baseline_tangent_max_abs": b["max_abs"], "candidate_tangent_max_abs": c["max_abs"],
                    "candidate_minus_baseline_tangent_max_abs": (
                        None if b["max_abs"] is None or c["max_abs"] is None else c["max_abs"]-b["max_abs"]),
                    "baseline_tangent_rms": b["rms"], "candidate_tangent_rms": c["rms"],
                    "candidate_minus_baseline_tangent_rms": (
                        None if b["rms"] is None or c["rms"] is None else c["rms"]-b["rms"]),
                }
        paired.append(row)

    report = {
        "schema": "v55_discrete_initial_multiconstraint_rhs_tangent_v1", "status": "completed",
        "diagnostic_only": True, "source_commit": SOURCE_COMMIT, "input_ab_run_id": AB_RUN_ID,
        "resolutions": list(RESOLUTIONS), "arms": list(ARMS), "epsilons": list(EPSILONS),
        "fixed_radii": list(FIXED_RADII), "state_count": len(cases), "one_step_comparators": len(cases),
        "cases": cases, "paired_tangent_comparison": paired,
        "guardrails": [
            "All inputs are frozen t=0 checkpoints from one matched initial-data A/B campaign.",
            "The candidate initial metric closes the discrete H residual by its Newton solve; that is reported, not treated as independent admission.",
            "The production source commit and pinned vendor are read-only and exact.",
            "No long trajectory, production equation edit, gauge change, coefficient change, or center projection is performed.",
            "Source-sector tangent and assembly closures are not an independent covariant stress-energy conservation proof.",
        ],
    }
    out = Path("results/discrete-initial-constraint-tangent")
    out.mkdir(parents=True, exist_ok=True)
    (out/"report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+"\n", encoding="utf-8")
    print("DISCRETE_INITIAL_CONSTRAINT_TANGENT "+json.dumps({
        "status": report["status"], "state_count": len(cases),
        "paired_tangent_comparison": paired,
        "closure_maxima": {
            key: max(c["closures"][key] for c in cases if isinstance(c["closures"].get(key), float))
            for key in ("base_H_reconstruction_max_abs", "base_M_geometry_reconstruction_max_abs",
                        "base_M_reconstruction_max_abs", "tangent_H_reconstruction_max_abs",
                        "tangent_M_geometry_reconstruction_max_abs", "tangent_M_reconstruction_max_abs")
        },
    }, sort_keys=True))


if __name__ == "__main__":
    main()
