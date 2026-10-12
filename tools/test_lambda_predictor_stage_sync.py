#!/usr/bin/env python3
"""Compare baseline PIRK predictor ordering with a Lambda-RHS stage-sync probe.

All source inputs are frozen t=0/t~10 states from the admitted production
campaign. Each baseline and probe path takes exactly one step. The only
counterfactual difference is: for the second geometry_stage_terms call, recompute
lambda_l2 after applying the predictor-stage Aa/K update implied by the already
computed PIRK primary blocks. No production file or physics equation is edited.
"""
from __future__ import annotations

import copy
import json
import math
from contextlib import contextmanager
from pathlib import Path
import sys

import numpy as np

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
CAMPAIGN_RUN_ID = 38208253687
LAMBDA_M = 2.0
R_CUT = 5.0
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
CONSERVED_FIELDS = ("rest", "energy_t", "momentum_r")


def d1(grid, values, parity=1):
    return np.asarray(grid.cell_derivative_fourth(
        np.asarray(values, dtype=float), parity=parity
    ), dtype=float)


def metric_connection(grid, geometry):
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    return d1(grid, a, +1) / (2.0*a*a) - d1(grid, b, +1)/(a*b) + 2.0/r*(1.0/b-1.0/a)


def state_metrics(kernel, state):
    grid = state.grid
    geometry = state.geometry
    total = assemble_total_stress_energy(grid, geometry, state.scalars, state.matter)
    raw = vacuum.constraints(grid, geometry)
    radius = np.asarray(grid.centers, dtype=float)
    H = np.asarray(raw["hamiltonian"], dtype=float) - 16.0*math.pi*np.asarray(total.rho, dtype=float)
    M = np.asarray(raw["momentum"], dtype=float) - 8.0*math.pi*np.asarray(total.j, dtype=float)
    C = np.asarray(geometry.Lambda, dtype=float) - metric_connection(grid, geometry)
    det = np.asarray(geometry.a, dtype=float)*np.asarray(geometry.b, dtype=float)**2 - 1.0
    fields = {"H": H, "M": M, "C": C, "determinant": det}
    masks = {
        "full_domain": np.ones(radius.shape, dtype=bool),
        "formal_stencil_i_ge_2": np.arange(radius.size) >= 2,
        "r_lt_5": radius < R_CUT,
        "r_ge_5": radius >= R_CUT,
    }
    summary = {}
    for name, values in fields.items():
        summary[name] = {}
        for region, mask in masks.items():
            inds = np.flatnonzero(mask & np.isfinite(values))
            if not inds.size:
                summary[name][region] = {"count": 0}
                continue
            k = int(inds[np.argmax(np.abs(values[inds]))])
            summary[name][region] = {
                "count": int(inds.size),
                "max_abs": float(abs(values[k])),
                "peak_signed": float(values[k]),
                "peak_r": float(radius[k]),
                "peak_i": k,
                "rms": float(np.sqrt(np.mean(values[inds]**2))),
            }
    first = {
        str(i): {
            "r": float(radius[i]),
            **{name: float(values[i]) for name, values in fields.items()},
        }
        for i in range(min(5, radius.size))
    }
    return {"t": float(state.t), "summary": summary, "first_five_cells": first}


def restore_checkpoint(path: Path, kernel):
    with np.load(path, allow_pickle=False) as archive:
        data = {key: archive[key].copy() for key in archive.files}
    n = int(data["r"].size)
    dr = float(data["r"][1]-data["r"][0])
    r_max = float(n*dr)
    state = kernel.initialize(resolution=n, r_max=r_max, include_radiation=True)
    if not np.allclose(state.grid.centers, data["r"], rtol=0.0,
                       atol=1.0e-12*max(1.0, r_max)):
        raise ValueError(f"{path}: checkpoint grid differs from reconstructed grid")
    for name in GEOMETRY_FIELDS:
        setattr(state.geometry, name, np.asarray(data[name], dtype=float).copy())
    state.scalars = ScalarFields(*(
        np.asarray(data[name], dtype=float).copy() for name in SCALAR_FIELDS
    ))
    species = {}
    for name, prefix in SPECIES.items():
        species[name] = ConservedSpecies(*(
            np.asarray(data[f"{prefix}_{field}"], dtype=float).copy()
            for field in CONSERVED_FIELDS
        ))
    state.matter = V55MatterState(
        dark_matter=species["dark_matter"],
        baryons=species["baryons"],
        radiation=species["radiation"],
    )
    state.t = float(data["t"])
    state.tau = float(data["tau"])
    state.e_folds = float(data["e_folds"])
    state.history = []
    state.handoffs = []
    return state


@contextmanager
def synchronized_predictor_lambda_l2(dt):
    """Patch only the diagnostic process's 2nd predictor RHS result; always restore."""
    original = adapter.geometry_stage_terms
    calls = {"count": 0, "first_terms": None}
    record = {}

    def wrapper(grid, geometry, fields, matter, *args, **kwargs):
        calls["count"] += 1
        terms = original(grid, geometry, fields, matter, *args, **kwargs)
        if calls["count"] == 1:
            calls["first_terms"] = terms
            return terms
        if calls["count"] != 2:
            return terms

        first = calls["first_terms"]
        if first is None:
            raise RuntimeError("predictor trace called before start-stage RHS")

        # This matches the production g1.Aa/K PIRK primary-stage formulas
        # immediately after this call returns, using the already-computed terms.
        sync_geometry = geometry.copy()
        sync_geometry.Aa = (
            np.asarray(geometry.Aa, dtype=float)
            + dt*(
                0.5*np.asarray(first["primary_l2"]["Aa"], dtype=float)
                + 0.5*np.asarray(terms["primary_l2"]["Aa"], dtype=float)
                + np.asarray(first["primary_l3"]["Aa"], dtype=float)
            )
        )
        sync_geometry.K = (
            np.asarray(geometry.K, dtype=float)
            + dt*(
                0.5*np.asarray(first["primary_l2"]["K"], dtype=float)
                + 0.5*np.asarray(terms["primary_l2"]["K"], dtype=float)
                + np.asarray(first["primary_l3"]["K"], dtype=float)
            )
        )
        lambda_m = kwargs.get("lambda_m", LAMBDA_M)
        corrected_lambda_l2 = vacuum.lambda_l2_rhs(
            grid, sync_geometry, lambda_m=lambda_m
        )
        record.update({
            "stage": "explicit_predictor_after_predicted_Aa_K",
            "lambda_m": float(lambda_m),
            "max_abs_Aa_predictor_change": float(np.max(np.abs(sync_geometry.Aa - geometry.Aa))),
            "max_abs_K_predictor_change": float(np.max(np.abs(sync_geometry.K - geometry.K))),
            "lambda_l2_stale_first_five": np.asarray(terms["lambda_l2"], dtype=float)[:5].tolist(),
            "lambda_l2_synced_first_five": np.asarray(corrected_lambda_l2, dtype=float)[:5].tolist(),
            "lambda_l2_difference_first_five": (
                np.asarray(corrected_lambda_l2, dtype=float)
                - np.asarray(terms["lambda_l2"], dtype=float)
            )[:5].tolist(),
            "lambda_l2_max_abs_difference": float(np.max(np.abs(
                np.asarray(corrected_lambda_l2, dtype=float)
                - np.asarray(terms["lambda_l2"], dtype=float)
            ))),
        })
        replacement = dict(terms)
        replacement["lambda_l2"] = corrected_lambda_l2
        return replacement

    adapter.geometry_stage_terms = wrapper
    try:
        yield record
    finally:
        adapter.geometry_stage_terms = original


def compare_case(state, label, kernel):
    dt = 0.03*float(state.grid.dr)
    baseline_input = copy.deepcopy(state)
    probe_input = copy.deepcopy(state)

    base_initial = state_metrics(kernel, baseline_input)
    base_output = kernel.step(baseline_input, dt)
    base_final = state_metrics(kernel, base_output)

    with synchronized_predictor_lambda_l2(dt) as stage_record:
        probe_output = kernel.step(probe_input, dt)
    probe_final = state_metrics(kernel, probe_output)

    delta = {}
    for field in ("H", "M", "C", "determinant"):
        delta[field] = {}
        for region in ("formal_stencil_i_ge_2", "r_ge_5", "r_lt_5"):
            a = base_final["summary"][field][region]
            b = probe_final["summary"][field][region]
            delta[field][region] = {
                "baseline_max_abs": a["max_abs"],
                "synced_max_abs": b["max_abs"],
                "change_in_max_abs": b["max_abs"] - a["max_abs"],
                "baseline_rms": a["rms"],
                "synced_rms": b["rms"],
                "change_in_rms": b["rms"] - a["rms"],
            }

    first_cell_delta = {}
    for i in range(min(5, len(state.grid.centers))):
        key = str(i)
        first_cell_delta[key] = {
            field: {
                "baseline": base_final["first_five_cells"][key][field],
                "synced": probe_final["first_five_cells"][key][field],
                "delta": (
                    probe_final["first_five_cells"][key][field]
                    - base_final["first_five_cells"][key][field]
                ),
            }
            for field in ("H", "M", "C", "determinant")
        }

    return {
        "case": label,
        "resolution": int(len(state.grid.centers)),
        "r_max": float(len(state.grid.centers)*state.grid.dr),
        "base_t": float(state.t),
        "dt": dt,
        "steps_per_path": 1,
        "source_equations_changed": False,
        "baseline_initial": base_initial,
        "baseline_final": base_final,
        "synced_final": probe_final,
        "synced_minus_baseline_final": delta,
        "first_five_cell_baseline_vs_sync": first_cell_delta,
        "stage_sync_probe": stage_record,
    }


def main():
    kernel = V55ProductionKernel()
    cases = []
    for n in (60, 80):
        initial = kernel.initialize(resolution=n, r_max=80.0, include_radiation=True)
        cases.append(compare_case(initial, f"N{n}_INITIAL_T0", kernel))
        checkpoint = Path(f"runs/t12/N{n}/checkpoint_t00010.000000.npz")
        if not checkpoint.exists():
            raise FileNotFoundError(f"required checkpoint missing: {checkpoint}")
        saved = restore_checkpoint(checkpoint, kernel)
        cases.append(compare_case(saved, f"N{n}_SAVED_T10", kernel))

    report = {
        "schema": "lambda_pirk_stage_sync_counterfactual_v1",
        "status": "completed",
        "kind": "diagnostic_only_single_step_counterfactual",
        "source_commit": SOURCE_COMMIT,
        "campaign_run_id": CAMPAIGN_RUN_ID,
        "settings": {
            "resolutions": [60, 80],
            "r_max": 80.0,
            "CFL": 0.03,
            "inputs": ["production initializer t=0", "saved t~10 checkpoints"],
            "baseline": "unchanged production stage ordering",
            "counterfactual": "recompute predictor lambda_l2 after predicted Aa/K update only",
            "full_evolution": False,
        },
        "cases": cases,
        "interpretation_guardrail": (
            "A one-step counterfactual change does not confirm root cause or authorize a "
            "production patch. The reference implementation has a related stage ordering, "
            "but production CMC and outer-boundary operations differ. Compare H, M, C, "
            "determinant and fixed-physical-radius witnesses before deciding."
        ),
    }
    out = Path("results/lambda-pirk-stage-sync")
    out.mkdir(parents=True, exist_ok=True)
    (out/"lambda_pirk_stage_sync.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+"\n"
    )
    print("LAMBDA_PIRK_STAGE_SYNC " + json.dumps({
        "status": report["status"],
        "cases": [{
            "case": c["case"],
            "dt": c["dt"],
            "i2_H_baseline": c["baseline_final"]["first_five_cells"]["2"]["H"],
            "i2_H_synced": c["synced_final"]["first_five_cells"]["2"]["H"],
            "i2_M_baseline": c["baseline_final"]["first_five_cells"]["2"]["M"],
            "i2_M_synced": c["synced_final"]["first_five_cells"]["2"]["M"],
            "i2_C_baseline": c["baseline_final"]["first_five_cells"]["2"]["C"],
            "i2_C_synced": c["synced_final"]["first_five_cells"]["2"]["C"],
            "H_formal_stencil": c["synced_minus_baseline_final"]["H"]["formal_stencil_i_ge_2"],
            "M_formal_stencil": c["synced_minus_baseline_final"]["M"]["formal_stencil_i_ge_2"],
            "C_formal_stencil": c["synced_minus_baseline_final"]["C"]["formal_stencil_i_ge_2"],
            "predictor_lambda_l2_max_delta": c["stage_sync_probe"]["lambda_l2_max_abs_difference"],
        } for c in cases]
    }, sort_keys=True))


if __name__ == "__main__":
    main()
