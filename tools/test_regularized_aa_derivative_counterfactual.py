#!/usr/bin/env python3
"""One-step counterfactual for tensor-regularized D_r Aa in BSSN constraints.

Uses unmodified production source at a pinned commit and four frozen inputs:
N=60/80 at t=0 and saved t~10. The only diagnostic-process substitution is to
evaluate D_r Aa via the vendor's existing regularized_aa_prime identity in the
momentum constraint and Lambda L2 RHS. It does not project fields or modify
any production file. Both direct and regularized momentum residuals are reported
on each output state so a change of diagnostic definition cannot masquerade as
an evolution improvement.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import sys

import numpy as np

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
CAMPAIGN_RUN_ID = 38208253687
R_CUT = 5.0
LAMBDA_M = 2.0
FIXED_RADII = (1.0, 2.0, 2.5, 3.0, 4.0, 5.0)
SOURCE = Path("source").resolve()
sys.path.insert(0, str(SOURCE))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()

GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED = ("rest", "energy_t", "momentum_r")

ORIGINAL_MOMENTUM = vacuum.momentum_constraint
ORIGINAL_LAMBDA_L2 = vacuum.lambda_l2_rhs
REGULARIZED_DERIVATIVE = vacuum.regularized_aa_prime


def d1(grid, values, parity):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(values, dtype=float), parity=parity), dtype=float)


def d2(grid, values, parity):
    return np.asarray(grid.cell_second_derivative_fourth(np.asarray(values, dtype=float), parity=parity), dtype=float)


def momentum_constraint_regularized(grid, state):
    """Pinned momentum expression; only D_r Aa uses the tensor identity."""
    r = np.asarray(grid.centers, dtype=float)
    Aa = np.asarray(state.Aa, dtype=float)
    Ab = -0.5 * Aa
    chip = -0.5 * d1(grid, state.X, +1) / np.asarray(state.X, dtype=float)
    return (
        REGULARIZED_DERIVATIVE(grid, Aa)
        - (2.0 / 3.0) * d1(grid, state.K, +1)
        + 6.0 * Aa * chip
        + (Aa - Ab) * (2.0 / r + d1(grid, state.b, +1) / np.asarray(state.b, dtype=float))
    )


def lambda_l2_rhs_regularized(grid, state, *, lambda_m=2.0):
    """Exact pinned Lambda L2 expression with regularized D_r Aa only."""
    r = np.asarray(grid.centers, dtype=float)
    div_beta = vacuum.shift_divergence(grid, state)
    Aa = np.asarray(state.Aa, dtype=float)
    Ab = -0.5 * Aa
    alphap = d1(grid, state.alpha, +1)
    Aap = REGULARIZED_DERIVATIVE(grid, Aa)
    beta_over_r = np.asarray(state.beta, dtype=float) / r
    return (
        d2(grid, state.beta, -1) / np.asarray(state.a, dtype=float)
        + 2.0 * d1(grid, beta_over_r, +1) / np.asarray(state.b, dtype=float)
        + d1(grid, div_beta, +1) / (3.0 * np.asarray(state.a, dtype=float))
        - 2.0 * (Aa * alphap + np.asarray(state.alpha, dtype=float) * Aap) / np.asarray(state.a, dtype=float)
        - 4.0 * np.asarray(state.alpha, dtype=float) * (Aa - Ab) / (r * np.asarray(state.b, dtype=float))
        + lambda_m * np.asarray(state.alpha, dtype=float) * momentum_constraint_regularized(grid, state) / np.asarray(state.a, dtype=float)
    )


def restore_checkpoint(path: Path, kernel):
    with np.load(path, allow_pickle=False) as z:
        data = {key: z[key].copy() for key in z.files}
    n = int(data["r"].size)
    dr = float(data["r"][1] - data["r"][0])
    rmax = float(n * dr)
    state = kernel.initialize(resolution=n, r_max=rmax, include_radiation=True)
    if not np.allclose(state.grid.centers, data["r"], rtol=0.0, atol=1e-12 * max(1.0, rmax)):
        raise ValueError(f"{path}: reconstructed grid mismatch")
    for name in GEOMETRY_FIELDS:
        setattr(state.geometry, name, np.asarray(data[name], dtype=float).copy())
    state.scalars = ScalarFields(*(np.asarray(data[name], dtype=float).copy() for name in SCALAR_FIELDS))
    parts = {}
    for name, prefix in SPECIES.items():
        parts[name] = ConservedSpecies(*(np.asarray(data[f"{prefix}_{k}"], dtype=float).copy() for k in CONSERVED))
    state.matter = V55MatterState(
        dark_matter=parts["dark_matter"], baryons=parts["baryons"], radiation=parts["radiation"]
    )
    state.t = float(data["t"])
    state.tau = float(data["tau"])
    state.e_folds = float(data["e_folds"])
    state.history = []
    state.handoffs = []
    return state


def summarise(values, radius):
    values = np.asarray(values, dtype=float)
    masks = {
        "full_domain": np.ones(radius.shape, dtype=bool),
        "stencil_i_ge_2": np.arange(radius.size) >= 2,
        "r_lt_5": radius < R_CUT,
        "stencil_i_ge_2_and_r_lt_5": (np.arange(radius.size) >= 2) & (radius < R_CUT),
        "r_ge_5": radius >= R_CUT,
    }
    out = {}
    for label, mask in masks.items():
        ids = np.flatnonzero(mask & np.isfinite(values))
        if not ids.size:
            out[label] = {"count": 0}
            continue
        idx = int(ids[np.argmax(np.abs(values[ids]))])
        out[label] = {
            "count": int(ids.size),
            "max_abs": float(abs(values[idx])),
            "peak_signed": float(values[idx]),
            "peak_r": float(radius[idx]),
            "peak_i": idx,
            "rms": float(np.sqrt(np.mean(values[ids] ** 2))),
        }
    out["fixed_radius_linear_samples"] = {
        str(target): float(np.interp(target, radius, values))
        for target in FIXED_RADII if radius[0] <= target <= radius[-1]
    }
    out["first_five_cells"] = {
        str(i): {"r": float(radius[i]), "value": float(values[i])}
        for i in range(min(5, radius.size))
    }
    return out


def state_witness(state, label):
    grid, geom = state.grid, state.geometry
    radius = np.asarray(grid.centers, dtype=float)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    h_raw = np.asarray(vacuum.constraints(grid, geom)["hamiltonian"], dtype=float)
    source_H = 16.0 * math.pi * np.asarray(total.rho, dtype=float)
    H = h_raw - source_H
    J8 = 8.0 * math.pi * np.asarray(total.j, dtype=float)
    # Always evaluate both operators explicitly for apples-to-apples checking.
    M_geo_direct = np.asarray(ORIGINAL_MOMENTUM(grid, geom), dtype=float)
    M_geo_regularized = np.asarray(momentum_constraint_regularized(grid, geom), dtype=float)
    M_direct = M_geo_direct - J8
    M_regularized = M_geo_regularized - J8
    metric_lambda = (
        d1(grid, geom.a, +1) / (2.0 * np.asarray(geom.a, dtype=float) ** 2)
        - d1(grid, geom.b, +1) / (np.asarray(geom.a, dtype=float) * np.asarray(geom.b, dtype=float))
        + 2.0 / radius * (1.0 / np.asarray(geom.b, dtype=float) - 1.0 / np.asarray(geom.a, dtype=float))
    )
    C = np.asarray(geom.Lambda, dtype=float) - metric_lambda
    det = np.asarray(geom.a, dtype=float) * np.asarray(geom.b, dtype=float) ** 2 - 1.0
    return {
        "label": label, "t": float(state.t),
        "H": summarise(H, radius),
        "M_direct": summarise(M_direct, radius),
        "M_regularized_operator": summarise(M_regularized, radius),
        "difference_M_regularized_minus_direct": summarise(M_regularized - M_direct, radius),
        "C_lambda": summarise(C, radius),
        "determinant": summarise(det, radius),
        "raw_derivative_Aa": summarise(d1(grid, geom.Aa, +1), radius),
        "regularized_derivative_Aa": summarise(REGULARIZED_DERIVATIVE(grid, geom.Aa), radius),
        "derivative_difference_regularized_minus_raw": summarise(REGULARIZED_DERIVATIVE(grid, geom.Aa) - d1(grid, geom.Aa, +1), radius),
    }


def compare_case(state, label, kernel):
    dt = 0.03 * float(state.grid.dr)
    baseline_in = copy.deepcopy(state)
    probe_in = copy.deepcopy(state)
    pre = state_witness(state, f"{label}_INPUT")

    baseline_out = kernel.step(baseline_in, dt)
    baseline_witness = state_witness(baseline_out, f"{label}_BASELINE_OUTPUT")

    try:
        vacuum.momentum_constraint = momentum_constraint_regularized
        vacuum.lambda_l2_rhs = lambda_l2_rhs_regularized
        probe_out = kernel.step(probe_in, dt)
        probe_witness = state_witness(probe_out, f"{label}_REGULARIZED_OUTPUT")
    finally:
        vacuum.momentum_constraint = ORIGINAL_MOMENTUM
        vacuum.lambda_l2_rhs = ORIGINAL_LAMBDA_L2

    # Compare changed state under both definitions, not only its native diagnostic.
    delta = {}
    for name in ("H", "M_direct", "M_regularized_operator", "C_lambda", "determinant"):
        b = baseline_witness[name]
        p = probe_witness[name]
        for region in ("stencil_i_ge_2", "stencil_i_ge_2_and_r_lt_5", "r_ge_5"):
            bv, pv = b[region], p[region]
            delta.setdefault(name, {})[region] = {
                "baseline_max_abs": bv.get("max_abs"),
                "probe_max_abs": pv.get("max_abs"),
                "delta_max_abs": None if bv.get("max_abs") is None or pv.get("max_abs") is None else pv["max_abs"] - bv["max_abs"],
                "baseline_rms": bv.get("rms"),
                "probe_rms": pv.get("rms"),
                "delta_rms": None if bv.get("rms") is None or pv.get("rms") is None else pv["rms"] - bv["rms"],
            }
    return {
        "case": label, "resolution": int(len(state.grid.centers)), "r_max": float(len(state.grid.centers) * state.grid.dr),
        "base_t": float(state.t), "dt": dt, "production_steps_per_path": 1,
        "production_files_modified": False, "pre_step": pre,
        "baseline_output": baseline_witness, "regularized_operator_output": probe_witness,
        "probe_state_compared_using_production_operator": {
            "M_direct": probe_witness["M_direct"],
            "M_regularized_operator": probe_witness["M_regularized_operator"],
        },
        "probe_minus_baseline": delta,
    }


def main():
    kernel = V55ProductionKernel()
    cases = []
    for n in (60, 80):
        initial = kernel.initialize(resolution=n, r_max=80.0, include_radiation=True)
        cases.append(compare_case(initial, f"N{n}_INITIAL_T0", kernel))
        path = Path(f"runs/t12/N{n}/checkpoint_t00010.000000.npz")
        if not path.exists():
            raise FileNotFoundError(f"required checkpoint missing: {path}")
        saved = restore_checkpoint(path, kernel)
        cases.append(compare_case(saved, f"N{n}_SAVED_T10", kernel))

    report = {
        "schema": "regularized_Aa_derivative_one_step_counterfactual_v1",
        "status": "completed",
        "kind": "diagnostic_only_no_projection_no_production_equation_change",
        "source_commit": SOURCE_COMMIT,
        "campaign_run_id": CAMPAIGN_RUN_ID,
        "settings": {
            "resolutions": [60, 80], "r_max": 80.0, "CFL": 0.03,
            "inputs": ["production initializer t=0", "saved checkpoint t~10"],
            "counterfactual": "use pinned vendor regularized_aa_prime for D_r Aa in momentum constraint and Lambda L2 only",
            "full_trajectory": False,
        },
        "cases": cases,
        "interpretation_guardrail": (
            "The alternative is a frozen one-step operator counterfactual, not an approved repair. "
            "Compare alternative-state constraints under both direct and regularized definitions. "
            "Do not treat a better residual under the alternative definition alone as evidence of a better state."
        ),
    }
    out = Path("results/regularized-aa-derivative-counterfactual")
    out.mkdir(parents=True, exist_ok=True)
    (out / "regularized_aa_derivative_counterfactual.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print("REGULARIZED_AA_DERIVATIVE_COUNTERFACTUAL " + json.dumps({
        "status": report["status"],
        "cases": [{
            "case": c["case"], "dt": c["dt"],
            "input_M_direct_i2": c["pre_step"]["M_direct"]["first_five_cells"]["2"]["value"],
            "input_M_regularized_i2": c["pre_step"]["M_regularized_operator"]["first_five_cells"]["2"]["value"],
            "baseline_final_M_direct_i2": c["baseline_output"]["M_direct"]["first_five_cells"]["2"]["value"],
            "probe_final_M_direct_i2": c["regularized_operator_output"]["M_direct"]["first_five_cells"]["2"]["value"],
            "probe_final_M_regularized_i2": c["regularized_operator_output"]["M_regularized_operator"]["first_five_cells"]["2"]["value"],
            "baseline_final_H_i2": c["baseline_output"]["H"]["first_five_cells"]["2"]["value"],
            "probe_final_H_i2": c["regularized_operator_output"]["H"]["first_five_cells"]["2"]["value"],
            "baseline_final_C_i2": c["baseline_output"]["C_lambda"]["first_five_cells"]["2"]["value"],
            "probe_final_C_i2": c["regularized_operator_output"]["C_lambda"]["first_five_cells"]["2"]["value"],
            "H_formal_stencil_delta_max": c["probe_minus_baseline"]["H"]["stencil_i_ge_2"]["delta_max_abs"],
            "M_direct_formal_stencil_delta_max": c["probe_minus_baseline"]["M_direct"]["stencil_i_ge_2"]["delta_max_abs"],
            "C_formal_stencil_delta_max": c["probe_minus_baseline"]["C_lambda"]["stencil_i_ge_2"]["delta_max_abs"],
        } for c in cases]
    }, sort_keys=True))


if __name__ == "__main__":
    main()
