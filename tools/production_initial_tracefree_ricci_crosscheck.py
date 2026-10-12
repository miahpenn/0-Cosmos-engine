#!/usr/bin/env python3
"""Same-state radial trace-free Ricci cross-check on actual production initial data.

For each resolution, compare the pinned vendor's R^r_r-R/3 with an independent
metric-derived form computed from areal radius and proper radial distance on the
same initialized slice. Split the difference into the evolved-Lambda contribution
and the remaining route/discretization contribution. Also record the actual
Aa RHS split to establish whether the Ricci term is the leading first-cell
tendency. No evolution or production edit is made.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path("source").resolve()))
from engine.production_kernel import V55ProductionKernel, LAMBDA_M
from engine import v55_pirk_adapter as adapter
from engine.stress_energy import assemble_total_stress_energy

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
RESOLUTIONS = (40, 60, 80, 120)
R_MAX = 80.0
ANNULUS = (2.0, 5.0)
RADII = (2.0, 2.5, 3.0, 4.0, 5.0)
CENTRE_CELLS = (0, 1, 2, 3, 4)


def d1(grid, values, parity):
    return np.asarray(grid.cell_derivative_fourth(values, parity=parity), dtype=float)


def metric_connection(grid, g):
    r = np.asarray(grid.centers, dtype=float)
    a, b = np.asarray(g.a), np.asarray(g.b)
    ap, bp = d1(grid, a, +1), d1(grid, b, +1)
    return ap / (2.0 * a**2) - bp / (a*b) + 2.0/r * (1.0/b - 1.0/a)


def independent_metric_ricci(grid, g):
    """Mixed Ricci components from ds²=dℓ²+R(ℓ)²dΩ²."""
    r = np.asarray(grid.centers, dtype=float)
    areal = r * np.sqrt(np.asarray(g.b)) / np.asarray(g.X)
    d_ell_dr = np.asarray(g.X) / np.sqrt(np.asarray(g.a))
    areal_ell = d_ell_dr * d1(grid, areal, -1)
    areal_ell_ell = d_ell_dr * d1(grid, areal_ell, +1)
    rrr = -2.0 * areal_ell_ell / areal
    rtheta = (1.0 - areal_ell**2 - areal * areal_ell_ell) / areal**2
    scalar = rrr + 2.0*rtheta
    return {
        "Rr_r": rrr, "Rtheta_theta": rtheta, "R": scalar,
        "R_tf": rrr - scalar/3.0,
        "areal_radius": areal,
        "areal_ell": areal_ell,
        "areal_ell_ell": areal_ell_ell,
    }


def peak_rms(v, mask):
    a = np.asarray(v, dtype=float)
    ids = np.flatnonzero(mask)
    if not ids.size:
        return {"count": 0, "rms": None, "max_abs": None}
    return {
        "count": int(ids.size),
        "rms": float(np.sqrt(np.mean(a[ids]**2))),
        "max_abs": float(np.max(np.abs(a[ids]))),
    }


def audit(n):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=n, r_max=R_MAX, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True,
    )
    grid, g = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=float)
    _, vacuum, _ = adapter.vendor_modules()

    connection_metric = metric_connection(grid, g)
    connection_residual = np.asarray(g.Lambda) - connection_metric
    metric_geometry = g.copy()
    metric_geometry.Lambda = connection_metric.copy()

    # Keep all other state variables and derivatives the same; only replace
    # Lambda in a diagnostic copy, exactly as the previously admitted scalar-R audit.
    terms_evolved = vacuum.geometry_terms(grid, g)
    terms_metric_lambda = vacuum.geometry_terms(grid, metric_geometry)
    metric_ricci = independent_metric_ricci(grid, g)

    tf_evolved = np.asarray(terms_evolved["Rr_r"]) - np.asarray(terms_evolved["R"])/3.0
    tf_metric_lambda = np.asarray(terms_metric_lambda["Rr_r"]) - np.asarray(terms_metric_lambda["R"])/3.0
    tf_metric = np.asarray(metric_ricci["R_tf"])
    tf_connection_piece = tf_evolved - tf_metric_lambda
    tf_route_piece = tf_metric_lambda - tf_metric
    tf_total_piece = tf_evolved - tf_metric
    split_closure = tf_total_piece - (tf_connection_piece + tf_route_piece)

    # Exact source-decomposed Aa RHS at t=0.
    total = assemble_total_stress_energy(grid, g, state.scalars, state.matter)
    stages = adapter.geometry_stage_terms(grid, g, state.scalars, state.matter, lambda_m=LAMBDA_M)
    vendor_l3 = vacuum.primary_l3_rhs(grid, g)
    aa_lapse = -(np.asarray(terms_evolved["DrDr_alpha"]) - np.asarray(terms_evolved["lap_alpha"])/3.0)
    aa_tf_evolved = np.asarray(g.alpha)*tf_evolved
    aa_tf_metric = np.asarray(g.alpha)*tf_metric
    aa_base_l3 = np.asarray(vendor_l3["Aa"])
    aa_matter = np.asarray(stages["primary_l3"]["Aa"]) - aa_base_l3
    aa_rhs_evolved = np.asarray(stages["primary_l2"]["Aa"]) + np.asarray(stages["primary_l3"]["Aa"])
    aa_rhs_metric_counterfactual = aa_lapse + aa_tf_metric + aa_base_l3 + aa_matter
    aa_rhs_rebuild = aa_lapse + aa_tf_evolved + aa_base_l3 + aa_matter

    ann = (r >= ANNULUS[0]) & (r < ANNULUS[1])
    profiles = {
        "Rr_r_evolved_Lambda": np.asarray(terms_evolved["Rr_r"]),
        "Rr_r_metric_Lambda": np.asarray(terms_metric_lambda["Rr_r"]),
        "Rr_r_metric": np.asarray(metric_ricci["Rr_r"]),
        "R_evolved_Lambda": np.asarray(terms_evolved["R"]),
        "R_metric_Lambda": np.asarray(terms_metric_lambda["R"]),
        "R_metric": np.asarray(metric_ricci["R"]),
        "R_tf_evolved_Lambda": tf_evolved,
        "R_tf_metric_Lambda": tf_metric_lambda,
        "R_tf_metric": tf_metric,
        "R_tf_connection_piece": tf_connection_piece,
        "R_tf_route_piece": tf_route_piece,
        "R_tf_total_difference": tf_total_piece,
        "connection_residual": connection_residual,
        "Aa_rhs_lapse_hessian": aa_lapse,
        "Aa_rhs_tf_Ricci_evolved": aa_tf_evolved,
        "Aa_rhs_tf_Ricci_metric": aa_tf_metric,
        "Aa_rhs_vacuum_L3": aa_base_l3,
        "Aa_rhs_matter_L3": aa_matter,
        "Aa_rhs_total_evolved": aa_rhs_evolved,
        "Aa_rhs_total_metric_Ricci_counterfactual": aa_rhs_metric_counterfactual,
        "Aa_rhs_reconstruction": aa_rhs_rebuild,
    }
    for key, values in profiles.items():
        if not np.all(np.isfinite(values)):
            raise FloatingPointError(f"N={n}: non-finite {key}")
    if float(np.max(np.abs(split_closure))) > 2.0e-12:
        raise AssertionError(f"N={n}: Ricci split failed to close")
    if float(np.max(np.abs(aa_rhs_rebuild-aa_rhs_evolved))) > 2.0e-12:
        raise AssertionError(f"N={n}: Aa RHS term split failed to close")

    first = []
    for i in CENTRE_CELLS:
        first.append({
            "i": i, "r": float(r[i]),
            **{k: float(v[i]) for k,v in profiles.items()},
            "R_tf_evolved_over_r2": float(tf_evolved[i]/r[i]**2),
            "R_tf_metric_Lambda_over_r2": float(tf_metric_lambda[i]/r[i]**2),
            "R_tf_metric_over_r2": float(tf_metric[i]/r[i]**2),
            "a_over_b_minus1_over_r2": float((g.a[i]/g.b[i]-1.0)/r[i]**2),
            "Lambda_evolved_over_r": float(g.Lambda[i]/r[i]),
            "Aa_over_r2": float(g.Aa[i]/r[i]**2),
        })
    fixed = {}
    for rr in RADII:
        fixed[str(rr)] = {k: float(np.interp(rr, r, v)) for k,v in profiles.items()}
    ann_summary = {
        k: peak_rms(v,ann) for k,v in profiles.items()
        if k.startswith("R_tf") or k.startswith("Aa_rhs")
    }
    return {
        "resolution": n, "dr": float(grid.dr), "r_max": R_MAX, "time": 0.0,
        "source_commit": SOURCE_COMMIT, "initialization": {
            "amplitude": 0.01, "width": 7.0, "D_amplitude": 1.0e-10, "radiation": True
        },
        "first_cells": first, "fixed_radius": fixed, "fixed_annulus": ann_summary,
        "full_domain_max_abs_connection_residual": float(np.max(np.abs(connection_residual))),
        "ricci_split_closure_max_abs": float(np.max(np.abs(split_closure))),
        "aa_rhs_split_closure_max_abs": float(np.max(np.abs(aa_rhs_rebuild-aa_rhs_evolved))),
        "interpretation_limit": "Same-state initial-slice geometry cross-check. Metric-derived curvature uses the pinned derivative operators but a separate metric formula. The counterfactual swaps only the radial trace-free Ricci term in a diagnostic RHS sum; no RHS is evolved.",
    }


def main():
    reports=[audit(n) for n in RESOLUTIONS]
    out=Path("results/production-initial-tracefree-ricci")
    out.mkdir(parents=True,exist_ok=True)
    result={
        "schema":"v55_production_initial_tracefree_ricci_crosscheck_v1",
        "status":"completed","diagnostic_only":True,"source_commit":SOURCE_COMMIT,
        "resolutions":list(RESOLUTIONS),"r_max":R_MAX,"reports":reports,
        "guardrails":[
            "No production evolution step was taken.",
            "The actual production initializer and CMC lapse solve were called unchanged.",
            "Vendor evolved-Lambda, vendor metric-Lambda and independent metric-derived Ricci are compared on the same initialized slice.",
            "The Aa RHS is decomposed without applying the metric-Ricci counterfactual to any trajectory.",
        ],
    }
    (out/"report.json").write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    summary=[]
    for x in reports:
        c=x["first_cells"][0]
        summary.append({
            "N":x["resolution"],"dr":x["dr"],
            "Rtf_vendor_i0":c["R_tf_evolved_Lambda"],
            "Rtf_metric_i0":c["R_tf_metric"],
            "Rtf_vendor_over_r2_i0":c["R_tf_evolved_over_r2"],
            "Rtf_metric_over_r2_i0":c["R_tf_metric_over_r2"],
            "connection_piece_i0":c["R_tf_connection_piece"],
            "route_piece_i0":c["R_tf_route_piece"],
            "Aa_rhs_total_i0":c["Aa_rhs_total_evolved"],
            "Aa_rhs_metric_counterfactual_i0":c["Aa_rhs_total_metric_Ricci_counterfactual"],
            "Aa_tf_evolved_i0":c["Aa_rhs_tf_Ricci_evolved"],
            "Aa_tf_metric_i0":c["Aa_rhs_tf_Ricci_metric"],
            "annulus_Rtf_route_error":x["fixed_annulus"]["R_tf_route_piece"]["rms"],
            "annulus_Rtf_connection_piece":x["fixed_annulus"]["R_tf_connection_piece"]["rms"],
        })
    print("PRODUCTION_INITIAL_TRACEFREE_RICCI_CROSSCHECK "+json.dumps(summary,sort_keys=True))


if __name__=="__main__":
    main()
