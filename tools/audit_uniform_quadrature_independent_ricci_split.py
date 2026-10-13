#!/usr/bin/env python3
"""Independent sampled-metric Ricci/Hessian cross-check on frozen uniform-quadrature initial states.

This diagnostic reconstructs, but never evolves or edits, the existing uniform
fourth-order radial-quadrature initializer. It compares the pinned vendor BSSN
TF Ricci/lapse-Hessian pieces with a separately assembled areal-radius/proper-
radial metric formula on the same sampled geometry.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from tools.audit_uniform_radial_quadrature_initial_profile import (
    candidate_uniform_initial_state,
)

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (40, 60, 80, 120, 160)
R_MAX = 80.0
FIXED_RADII = (2.0, 3.0, 4.0, 5.0, 8.0)
OUT = Path("results/uniform-quadrature-independent-ricci-split")
grid_ops, vacuum, _ = adapter.vendor_modules()


def _max_abs(values, mask):
    v = np.asarray(values, dtype=float)
    m = np.asarray(mask, dtype=bool) & np.isfinite(v)
    if not np.any(m):
        return None
    return float(np.max(np.abs(v[m])))


def _rms(values, mask):
    v = np.asarray(values, dtype=float)
    m = np.asarray(mask, dtype=bool) & np.isfinite(v)
    if not np.any(m):
        return None
    return float(np.sqrt(np.mean(v[m] ** 2)))


def _interp(values, r, target):
    return float(np.interp(float(target), r, np.asarray(values, dtype=float)))


def _state_digest(state):
    h = hashlib.sha256()
    for group, names in (
        (state.geometry, ("a", "b", "X", "alpha", "Aa", "K", "Lambda")),
        (state.scalars, ("S", "PS", "D", "PD", "phi", "Pi")),
    ):
        for name in names:
            arr = np.ascontiguousarray(np.asarray(getattr(group, name), dtype="<f8"))
            h.update(name.encode("ascii") + b"\0")
            h.update(arr.tobytes())
    for name in ("dark_matter", "baryons", "radiation"):
        species = getattr(state.matter, name)
        for component in ("rest", "energy_t", "momentum_r"):
            arr = np.ascontiguousarray(
                np.asarray(getattr(species, component), dtype="<f8")
            )
            h.update((name + "." + component).encode("ascii") + b"\0")
            h.update(arr.tobytes())
    return h.hexdigest()


def audit_one(n, kernel):
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    state = candidate_uniform_initial_state(kernel, grid)
    if state.t != 0.0:
        raise RuntimeError(f"N={n}: expected frozen initial slice t=0, got {state.t!r}")
    digest_before = _state_digest(state)
    r = np.asarray(grid.centers, dtype=float)
    g = state.geometry
    metric_factor = np.asarray(g.X, dtype=float) / np.sqrt(np.asarray(g.a, dtype=float))

    vendor_geometry = vacuum.geometry_terms(grid, g)
    vendor_l2 = vacuum.primary_l2_rhs(grid, g)
    vendor_l3 = vacuum.primary_l3_rhs(grid, g)
    vendor_ricci_tf = np.asarray(
        vendor_geometry["Rr_r"] - vendor_geometry["R"] / 3.0, dtype=float
    )
    vendor_hessian_tf = np.asarray(
        vendor_geometry["DrDr_alpha"] - vendor_geometry["lap_alpha"] / 3.0,
        dtype=float,
    )
    vendor_ricci_contribution = np.asarray(g.alpha, dtype=float) * vendor_ricci_tf
    vendor_hessian_contribution = -vendor_hessian_tf
    vendor_l2_rhs = np.asarray(vendor_l2["Aa"], dtype=float)
    vendor_l3_rhs = np.asarray(vendor_l3["Aa"], dtype=float)
    vendor_full_rhs = vendor_l2_rhs + vendor_l3_rhs

    # Independent geometric construction for dl^2 + R(l)^2 dOmega^2.
    # It uses sampled a,b,X,alpha plus the grid derivative stencil, not
    # vacuum.geometry_terms or its conformal-connection input.
    areal = r * np.sqrt(np.asarray(g.b, dtype=float)) / np.asarray(g.X, dtype=float)
    areal_r = np.asarray(grid.cell_derivative_fourth(areal, parity=-1), dtype=float)
    areal_l = metric_factor * areal_r
    areal_l_r = np.asarray(grid.cell_derivative_fourth(areal_l, parity=+1), dtype=float)
    areal_ll = metric_factor * areal_l_r
    independent_ricci_tf = -(2.0 / 3.0) * (
        areal_ll / areal + (1.0 - areal_l**2) / areal**2
    )

    alpha_r = np.asarray(
        grid.cell_derivative_fourth(np.asarray(g.alpha, dtype=float), parity=+1),
        dtype=float,
    )
    alpha_l = metric_factor * alpha_r
    alpha_l_r = np.asarray(grid.cell_derivative_fourth(alpha_l, parity=-1), dtype=float)
    alpha_ll = metric_factor * alpha_l_r
    independent_hessian_tf = (2.0 / 3.0) * (
        alpha_ll - (areal_l / areal) * alpha_l
    )

    independent_ricci_contribution = np.asarray(g.alpha, dtype=float) * independent_ricci_tf
    independent_hessian_contribution = -independent_hessian_tf
    independent_l2_rhs = independent_hessian_contribution + independent_ricci_contribution

    # Verify the reported vendor L2 identity separately from the cross-metric test.
    vendor_identity_rhs = vendor_hessian_contribution + vendor_ricci_contribution
    vendor_identity_closure = vendor_l2_rhs - vendor_identity_rhs
    ricci_delta = vendor_ricci_tf - independent_ricci_tf
    hessian_delta = vendor_hessian_tf - independent_hessian_tf
    l2_delta = vendor_l2_rhs - independent_l2_rhs
    radii = {
        str(target): {
            "vendor_ricci_tf": _interp(vendor_ricci_tf, r, target),
            "independent_metric_ricci_tf": _interp(independent_ricci_tf, r, target),
            "ricci_delta_vendor_minus_metric": _interp(ricci_delta, r, target),
            "vendor_hessian_tf": _interp(vendor_hessian_tf, r, target),
            "independent_metric_hessian_tf": _interp(independent_hessian_tf, r, target),
            "hessian_delta_vendor_minus_metric": _interp(hessian_delta, r, target),
            "vendor_L2_Aa_rhs": _interp(vendor_l2_rhs, r, target),
            "independent_metric_L2_Aa_rhs": _interp(independent_l2_rhs, r, target),
            "L2_delta_vendor_minus_metric": _interp(l2_delta, r, target),
            "full_vendor_Aa_rhs_including_L3": _interp(vendor_full_rhs, r, target),
            "full_vendor_Aa_rhs_over_r2": _interp(vendor_full_rhs / r**2, r, target),
        }
        for target in FIXED_RADII if r[0] <= target <= r[-1]
    }
    first_cells = {}
    for i in range(min(8, n)):
        first_cells[str(i)] = {
            "r": float(r[i]),
            "vendor_ricci_tf": float(vendor_ricci_tf[i]),
            "independent_metric_ricci_tf": float(independent_ricci_tf[i]),
            "ricci_delta_vendor_minus_metric": float(ricci_delta[i]),
            "vendor_ricci_contribution_alpha_times_Rtf": float(vendor_ricci_contribution[i]),
            "independent_metric_ricci_contribution": float(independent_ricci_contribution[i]),
            "vendor_hessian_tf": float(vendor_hessian_tf[i]),
            "independent_metric_hessian_tf": float(independent_hessian_tf[i]),
            "hessian_delta_vendor_minus_metric": float(hessian_delta[i]),
            "vendor_hessian_contribution_minus_hessian": float(vendor_hessian_contribution[i]),
            "independent_metric_hessian_contribution": float(independent_hessian_contribution[i]),
            "vendor_L2_Aa_rhs": float(vendor_l2_rhs[i]),
            "independent_metric_L2_Aa_rhs": float(independent_l2_rhs[i]),
            "L2_delta_vendor_minus_metric": float(l2_delta[i]),
            "vendor_L3_Aa_rhs": float(vendor_l3_rhs[i]),
            "full_vendor_Aa_rhs_including_L3": float(vendor_full_rhs[i]),
            "full_vendor_Aa_rhs_over_r2": float(vendor_full_rhs[i] / r[i]**2),
            "independent_metric_L2_Aa_rhs_over_r2": float(independent_l2_rhs[i] / r[i]**2),
            "vendor_identity_closure": float(vendor_identity_closure[i]),
            "connection_constraint": float(
                vacuum.constraints(grid, g)["connection"][i]
            ),
            "determinant_constraint": float(
                vacuum.constraints(grid, g)["determinant"][i]
            ),
        }

    region_masks = {
        "first_four_cells": np.arange(n) < 4,
        "first_eight_cells": np.arange(n) < min(8, n),
        "interior_r_lt_5": r < 5.0,
        "fixed_annulus_3_to_6": (r >= 3.0) & (r <= 6.0),
        "exterior_r_ge_5": r >= 5.0,
        "full_domain": np.ones(n, dtype=bool),
    }
    summary = {}
    for region, mask in region_masks.items():
        summary[region] = {
            "vendor_ricci_tf_max_abs": _max_abs(vendor_ricci_tf, mask),
            "independent_metric_ricci_tf_max_abs": _max_abs(independent_ricci_tf, mask),
            "ricci_delta_max_abs": _max_abs(ricci_delta, mask),
            "ricci_delta_rms": _rms(ricci_delta, mask),
            "vendor_hessian_tf_max_abs": _max_abs(vendor_hessian_tf, mask),
            "independent_metric_hessian_tf_max_abs": _max_abs(independent_hessian_tf, mask),
            "hessian_delta_max_abs": _max_abs(hessian_delta, mask),
            "hessian_delta_rms": _rms(hessian_delta, mask),
            "vendor_L2_Aa_rhs_max_abs": _max_abs(vendor_l2_rhs, mask),
            "independent_metric_L2_Aa_rhs_max_abs": _max_abs(independent_l2_rhs, mask),
            "L2_delta_max_abs": _max_abs(l2_delta, mask),
            "L2_delta_rms": _rms(l2_delta, mask),
            "vendor_identity_closure_max_abs": _max_abs(vendor_identity_closure, mask),
        }

    arrays = (
        vendor_ricci_tf, independent_ricci_tf, ricci_delta, vendor_hessian_tf,
        independent_hessian_tf, hessian_delta, vendor_l2_rhs, independent_l2_rhs,
        l2_delta, vendor_full_rhs,
    )
    finite = all(bool(np.all(np.isfinite(x))) for x in arrays)
    if not finite:
        raise FloatingPointError(f"N={n}: non-finite decomposition field")
    if _state_digest(state) != digest_before or state.t != 0.0:
        raise RuntimeError(f"N={n}: frozen initializer state was mutated")
    if maxabs := _max_abs(vendor_identity_closure, np.ones(n, dtype=bool)):
        if maxabs > 1.0e-12:
            raise AssertionError(f"N={n}: vendor RHS signed identity closure is {maxabs:g}")

    return {
        "N": n,
        "dr": float(grid.dr),
        "state_sha256": digest_before,
        "finite": finite,
        "vendor_identity_closure_full_domain_max_abs": _max_abs(
            vendor_identity_closure, np.ones(n, dtype=bool)
        ),
        "summary": summary,
        "first_cells": first_cells,
        "fixed_radius": radii,
        "interpretation_limit": (
            "Frozen t=0 initializer only. Independent metric formula shares the "
            "sampled radial derivative stencil, but does not use vendor Ricci "
            "or the evolved conformal connection. No evolution or production edit."
        ),
    }


def main():
    kernel = V55ProductionKernel()
    records = [audit_one(n, kernel) for n in RESOLUTIONS]
    report = {
        "schema": "uniform_quadrature_independent_metric_ricci_hessian_v1",
        "status": "completed_finite_closure_checked",
        "production_source_commit": SOURCE_COMMIT,
        "pinned_vendor_commit": VENDOR_COMMIT,
        "diagnostic_branch": "diag/uniform-quadrature-independent-ricci-split-2026-10-13",
        "kind": "reconstructed_uniform_quadrature_initial_data_no_evolution",
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "independent_formula": (
            "dl=sqrt(a)/X dr; R_areal=r*sqrt(b)/X; "
            "TF_Ricci^r_r=-(2/3)*(R_ll/R+(1-R_l^2)/R^2); "
            "TF_Hessian^r_r=(2/3)*(alpha_ll-(R_l/R)*alpha_l)"
        ),
        "records": records,
        "guardrails": [
            "Uses the archived uniform-quadrature initializer exactly; no production source edits.",
            "Does not modify or evolve the reconstructed states.",
            "No fit, projection, damping, clipping, source change, or physics promotion.",
            "Workflow PASS means finite decomposition and algebraic closure only, not scientific admission.",
            "The independent metric cross-check shares the radial derivative stencil, so it is not a formal continuum proof.",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    compact = [
        {
            "N": rec["N"],
            "state_sha256": rec["state_sha256"],
            "first_cells": {
                k: {
                    "r": v["r"],
                    "ricci_vendor": v["vendor_ricci_tf"],
                    "ricci_metric": v["independent_metric_ricci_tf"],
                    "ricci_delta": v["ricci_delta_vendor_minus_metric"],
                    "hessian_vendor": v["vendor_hessian_tf"],
                    "hessian_metric": v["independent_metric_hessian_tf"],
                    "hessian_delta": v["hessian_delta_vendor_minus_metric"],
                    "l2_rhs_delta": v["L2_delta_vendor_minus_metric"],
                    "vendor_rhs_over_r2": v["full_vendor_Aa_rhs_over_r2"],
                }
                for k, v in list(rec["first_cells"].items())[:4]
            },
            "interior": rec["summary"]["interior_r_lt_5"],
            "exterior": rec["summary"]["exterior_r_ge_5"],
            "vendor_identity_closure": rec["vendor_identity_closure_full_domain_max_abs"],
        }
        for rec in records
    ]
    print("UNIFORM_QUADRATURE_INDEPENDENT_RICCI_SPLIT " + json.dumps(compact, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
