"""Trace the production initial-data construction before the first evolution step.

Diagnostic-only: build the existing strong-D initial data, then record the
same Hamiltonian terms before centre regularity projection, after projection,
after the CMC lapse solve, and from the public kernel.initialize result. This
tests whether the O(1e-3) centre Hamiltonian residual is introduced by the
projection/lapse stages or is already present in the unprojected initializer.
No production equation, parameter, gauge, projection, or default is changed.
"""
import json
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import cosmos_potential
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_initial import build_initial_data
from diagnostics.curvature_consistency.independent_metric_ricci import (
    finite_array,
    metric_connection,
    metric_derived_ricci,
    metric_summary,
)

grid_ops, vacuum, _ = adapter.vendor_modules()

N = 160
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-4
INCLUDE_RADIATION = True
EPS = np.finfo(float).eps


def eval_stage(grid, geom, scalars, matter, label):
    radius = np.asarray(grid.centers, dtype=float)
    gt = vacuum.geometry_terms(grid, geom)
    raw_constraints = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, scalars, matter)

    R_vendor = finite_array("vendor_R", np.asarray(gt["R"], dtype=float))
    R_metric = finite_array(
        "metric_R", metric_derived_ricci(grid, geom)
    )
    aa = finite_array("Aa", np.asarray(geom.Aa, dtype=float))
    ab = -0.5 * aa
    term_A = -(aa**2 + 2.0 * ab**2)
    term_K = (2.0 / 3.0) * np.asarray(geom.K, dtype=float)**2
    rho_s = finite_array("scalar_rho", np.asarray(total.scalar_rho, dtype=float))
    rho_f = finite_array("fluid_rho", np.asarray(total.fluid_rho, dtype=float))
    rho = rho_s + rho_f
    term_s = -16.0 * math.pi * rho_s
    term_f = -16.0 * math.pi * rho_f
    term_matter = term_s + term_f
    H_raw = finite_array(
        "vendor H", np.asarray(raw_constraints["hamiltonian"], dtype=float) + term_matter
    )
    H_reconstructed = R_vendor + term_A + term_K + term_matter
    H_metric = R_metric + term_A + term_K + term_matter

    lambda_metric = finite_array(
        "metric Lambda", metric_connection(grid, geom)
    )
    C_lambda = np.asarray(geom.Lambda, dtype=float) - lambda_metric
    b_over_X2_minus_one = np.asarray(geom.b, dtype=float) / np.asarray(geom.X, dtype=float)**2 - 1.0

    gate_values = {
        "vendor_H_term_closure": float(np.max(np.abs(H_raw - H_reconstructed))),
        "curvature_replacement_closure": float(np.max(np.abs((H_raw - H_metric) - (R_vendor - R_metric)))),
    }
    scale = max(1.0, float(np.max(np.abs(R_vendor))), float(np.max(np.abs(R_metric))),
                float(np.max(np.abs(term_matter))), float(np.max(np.abs(term_K))))
    tolerance = 2048.0 * EPS * scale
    failures = [{"stage": label, "gate": k, "value": v, "tolerance": tolerance}
                for k,v in gate_values.items() if v > tolerance]

    fields = {
        "R_vendor": R_vendor,
        "R_metric": R_metric,
        "term_A": term_A,
        "term_K": term_K,
        "term_scalar": term_s,
        "term_fluid": term_f,
        "term_matter": term_matter,
        "H_vendor": H_raw,
        "H_metric": H_metric,
        "C_lambda": C_lambda,
        "a": np.asarray(geom.a, dtype=float),
        "b": np.asarray(geom.b, dtype=float),
        "X": np.asarray(geom.X, dtype=float),
        "Lambda_evolved": np.asarray(geom.Lambda, dtype=float),
        "Lambda_metric": lambda_metric,
        "Aa": np.asarray(geom.Aa, dtype=float),
        "K": np.asarray(geom.K, dtype=float),
        "alpha": np.asarray(geom.alpha, dtype=float),
        "B": np.asarray(geom.B, dtype=float),
    }
    summary = {}
    for key, values in fields.items():
        if key in ("a", "b", "X", "Lambda_evolved", "Lambda_metric", "Aa", "K", "alpha", "B"):
            v = np.asarray(values, dtype=float)
            summary[key] = {
                "cell0": float(v[0]),
                "first_five": [float(z) for z in v[:5]],
                "max_abs_all": float(np.max(np.abs(v))),
            }
        else:
            summary[key] = metric_summary(values, radius)

    return {
        "label": label,
        "summary": summary,
        "central_terms": {
            "R_vendor": float(R_vendor[0]),
            "R_metric": float(R_metric[0]),
            "term_A": float(term_A[0]),
            "term_K": float(term_K[0]),
            "term_scalar": float(term_s[0]),
            "term_fluid": float(term_f[0]),
            "term_matter": float(term_matter[0]),
            "rho_scalar": float(rho_s[0]),
            "rho_fluid": float(rho_f[0]),
            "rho_total": float(rho[0]),
            "H_vendor": float(H_raw[0]),
            "H_metric": float(H_metric[0]),
            "curvature_gap": float(R_vendor[0]-R_metric[0]),
            "C_lambda": float(C_lambda[0]),
            "b_over_X2_minus_one": float(b_over_X2_minus_one[0]),
        },
        "first_five_H_vendor": [float(x) for x in H_raw[:5]],
        "first_five_H_metric": [float(x) for x in H_metric[:5]],
        "gates": {**gate_values, "tolerance": tolerance, "failures": failures},
        "_arrays": {
            "R_vendor": R_vendor,
            "term_A": term_A,
            "term_K": term_K,
            "term_scalar": term_s,
            "term_fluid": term_f,
            "term_matter": term_matter,
            "H_vendor": H_raw,
        },
    }


def serializable_stage(stage):
    return {k:v for k,v in stage.items() if k != "_arrays"}


def main():
    kernel = V55ProductionKernel()
    grid = grid_ops.SphericalCellGrid(N, R_MAX)
    init = build_initial_data(
        grid, vacuum.flat_state, amplitude=AMPLITUDE, width=WIDTH,
        D_amplitude=D_AMPLITUDE, include_radiation=INCLUDE_RADIATION,
    )
    scalars = init.scalars
    matter = init.matter

    stages = []
    failures = []
    # Stage A: exact geometry returned by build_initial_data, before production projection.
    g_raw = init.geometry.copy()
    s_raw = eval_stage(grid, g_raw, scalars, matter, "A_build_initial_data_raw")
    stages.append(s_raw)
    failures.extend(s_raw["gates"]["failures"])

    # Stage B: apply the exact existing production centre projection.
    g_projected = kernel._enforce_center_regularity(grid, g_raw.copy())
    s_proj = eval_stage(grid, g_projected, scalars, matter, "B_after_centre_projection")
    stages.append(s_proj)
    failures.extend(s_proj["gates"]["failures"])

    # Stage C: same lapse and shift preparation used by kernel.initialize.
    g_gauge = g_projected.copy()
    g_gauge.alpha = kernel._solve_lapse(grid, g_gauge, scalars, matter)[0]
    g_gauge.beta.fill(0.0)
    g_gauge.B.fill(0.0)
    s_gauge = eval_stage(grid, g_gauge, scalars, matter, "C_after_CMC_lapse_and_gauge_reset")
    stages.append(s_gauge)
    failures.extend(s_gauge["gates"]["failures"])

    # Stage D: public production initialization must reproduce Stage C.
    prod_state = kernel.initialize(
        resolution=N, r_max=R_MAX, amplitude=AMPLITUDE, width=WIDTH,
        D_amplitude=D_AMPLITUDE, include_radiation=INCLUDE_RADIATION,
    )
    s_public = eval_stage(
        prod_state.grid, prod_state.geometry, prod_state.scalars, prod_state.matter,
        "D_public_kernel_initialize",
    )
    stages.append(s_public)
    failures.extend(s_public["gates"]["failures"])

    # Attribute raw->projected, projected->gauge and raw->public H delta exactly
    # to curvature, extrinsic-A/K and scalar/fluid matter terms on the same grid.
    comparisons = []
    pairs = (
        ("raw_to_projection", s_raw, s_proj),
        ("projection_to_gauge", s_proj, s_gauge),
        ("raw_to_public_initialize", s_raw, s_public),
    )
    for name, before, after in pairs:
        parts = {}
        for field in ("R_vendor", "term_A", "term_K", "term_scalar", "term_fluid"):
            delta = after["_arrays"][field] - before["_arrays"][field]
            parts[field] = {
                "central_delta": float(delta[0]),
                "max_abs_delta_all": float(np.max(np.abs(delta))),
            }
        delta_h = after["_arrays"]["H_vendor"] - before["_arrays"]["H_vendor"]
        part_sum = sum(
            after["_arrays"][f] - before["_arrays"][f]
            for f in ("R_vendor","term_A","term_K","term_scalar","term_fluid")
        )
        comparisons.append({
            "name": name,
            "H_vendor_central_before": float(before["_arrays"]["H_vendor"][0]),
            "H_vendor_central_after": float(after["_arrays"]["H_vendor"][0]),
            "delta_H_vendor_central": float(delta_h[0]),
            "delta_H_vendor_max_abs_all": float(np.max(np.abs(delta_h))),
            "term_deltas": parts,
            "component_delta_sum_closure_max_abs": float(np.max(np.abs(delta_h-part_sum))),
        })

    # Confirm the initializer's archived budget values and identify initial S/D/phi densities.
    r = np.asarray(grid.centers, dtype=float)
    Sp = grid.cell_derivative_fourth(scalars.S, parity=1)
    Dp = grid.cell_derivative_fourth(scalars.D, parity=1)
    B_archive = np.asarray(init.geometry.X, dtype=float)**6
    rho_S0 = 0.5*(np.asarray(scalars.PS)**2 + B_archive*Sp**2) + 0.5*np.asarray(scalars.S)**2
    rho_D0 = 0.5*(np.asarray(scalars.PD)**2 + B_archive*Dp**2) - 0.5*np.asarray(scalars.D)**2
    rho_phi0 = (0.5*np.asarray(scalars.Pi)**2 + cosmos_potential(scalars.phi))/(8.0*math.pi)

    report = {
        "schema": "strong_D_production_initial_stage_trace_v1",
        "status": "completed" if not failures else "completed_with_diagnostic_gate_failures",
        "diagnostic_only": True,
        "configuration": {
            "resolution": N, "r_max": R_MAX, "dr": float(grid.dr),
            "amplitude": AMPLITUDE, "width": WIDTH, "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "initializer": "build_initial_data / V55ProductionKernel.initialize, unchanged",
            "projection": "existing vacuum.enforce_algebraic_regularity",
            "production_physics_changed": False,
        },
        "initial_budget": {
            "H0": float(init.H0),
            "initial_scalar_S_rho_center": float(rho_S0[0]),
            "initial_scalar_D_rho_center": float(rho_D0[0]),
            "initial_cosmos_phi_rho_center": float(rho_phi0[0]),
            "initial_fluid_rho_center_postprojection_metric": float(s_public["central_terms"]["rho_fluid"]),
        },
        "stages": [serializable_stage(s) for s in stages],
        "stage_comparisons": comparisons,
        "diagnostic_gate_failures": failures,
        "interpretation_guardrail": (
            "This tests the origin of the central constraint residual within the "
            "initialization pipeline only. If the residual precedes projection, that "
            "rules out projection as its origin but does not by itself prove whether "
            "the radial integral, central stencil, stress-energy normalization, or "
            "other initializer convention is responsible. No correction or production "
            "admission is authorized by this report."
        ),
    }
    out = ROOT / "runs" / "strong-D-production-initial-stage-trace"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "report.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("[STRONG_D_INITIAL_STAGE_TRACE] report=" + str(path))
    for stage in report["stages"]:
        print("[STRONG_D_INITIAL_STAGE_TRACE] STAGE " +
              json.dumps({"label":stage["label"],"central_terms":stage["central_terms"],
                          "H_vendor":stage["summary"]["H_vendor"],
                          "gates":stage["gates"]}, sort_keys=True))
    print("[STRONG_D_INITIAL_STAGE_TRACE] COMPARISONS " + json.dumps(comparisons, sort_keys=True))
    print("[STRONG_D_INITIAL_STAGE_TRACE] gate_failures=" + json.dumps(failures, sort_keys=True))
    print("[STRONG_D_INITIAL_STAGE_TRACE] status=" + report["status"])
    print("[STRONG_D_INITIAL_STAGE_TRACE] report_json=" + json.dumps(report, sort_keys=True, allow_nan=False))
    if failures:
        raise RuntimeError("Initial stage diagnostic gate failed; report preserved.")


if __name__ == "__main__":
    main()
