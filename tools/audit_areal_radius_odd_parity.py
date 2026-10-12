#!/usr/bin/env python3
"""Frozen-checkpoint audit of areal-radius parity in invariant diagnostics.

The areal radius R=r*sqrt(b)/X is odd under regular radial reflection when b
and X are even. Compare the old even-parity derivative against the correct
odd-parity derivative on saved production checkpoints. No kernel steps, source
edits, or field modifications occur.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys
import numpy as np

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
CAMPAIGN_RUN_ID = 38208253687
RESOLUTIONS = (40, 60, 80)
TIMES = (2, 6, 12)
SOURCE = Path("source").resolve()
sys.path.insert(0, str(SOURCE))

from engine import v55_pirk_adapter as adapter

grid_ops, _, _ = adapter.vendor_modules()


def summarize(values, radius, mask):
    values = np.asarray(values, dtype=float)
    selected = np.asarray(mask, dtype=bool) & np.isfinite(values)
    ids = np.flatnonzero(selected)
    if not ids.size:
        return {"count": 0}
    i = int(ids[np.argmax(np.abs(values[ids]))])
    return {
        "count": int(ids.size),
        "max_abs": float(abs(values[i])),
        "peak_signed": float(values[i]),
        "peak_r": float(radius[i]),
        "peak_i": i,
        "rms": float(np.sqrt(np.mean(values[ids] ** 2))),
    }


def roots(radius, chi):
    out = []
    for i in np.flatnonzero(np.asarray(chi[:-1]) * np.asarray(chi[1:]) <= 0.0):
        if chi[i] != chi[i + 1]:
            out.append(float(radius[i] - chi[i] * (radius[i + 1] - radius[i])/(chi[i + 1]-chi[i])))
    return out


def audit(path, resolution, time):
    with np.load(path, allow_pickle=False) as z:
        data = {k: np.asarray(z[k]).copy() for k in z.files}
    r = np.asarray(data["r"], dtype=float)
    n = int(r.size)
    dr = float(r[1]-r[0])
    grid = grid_ops.SphericalCellGrid(n, n*dr)
    if not np.allclose(grid.centers, r, rtol=0.0, atol=1e-12*max(1.0,n*dr)):
        raise ValueError(f"{path}: checkpoint grid mismatch")

    a = np.asarray(data["a"], dtype=float)
    b = np.asarray(data["b"], dtype=float)
    X = np.asarray(data["X"], dtype=float)
    K = np.asarray(data["K"], dtype=float)
    Aa = np.asarray(data["Aa"], dtype=float)
    R = r*np.sqrt(b)/X

    # Previous implementation treated R as even; regular signed-radius R is odd.
    dR_old = np.asarray(grid.cell_derivative_fourth(R, parity=+1), dtype=float)
    dR_new = np.asarray(grid.cell_derivative_fourth(R, parity=-1), dtype=float)
    Ktheta = K/3.0 - Aa/2.0
    normal_dR = -R*Ktheta
    gamma_rr_inv = X*X/a
    chi_old = gamma_rr_inv*dR_old*dR_old-normal_dR*normal_dR
    chi_new = gamma_rr_inv*dR_new*dR_new-normal_dR*normal_dR
    mass_old = 0.5*R*(1.0-chi_old)
    mass_new = 0.5*R*(1.0-chi_new)

    masks = {
        "first_two_cells": np.arange(n)<2,
        "first_five_cells": np.arange(n)<5,
        "formal_stencil_i_ge_2": np.arange(n)>=2,
        "r_lt_5": r<5.0,
        "r_ge_5": r>=5.0,
        "full_domain": np.ones(n,dtype=bool),
    }
    summary = {}
    for name, values_old, values_new in (
        ("dR_dr",dR_old,dR_new),("chi",chi_old,chi_new),("misner_sharp_mass",mass_old,mass_new)
    ):
        summary[name] = {
            "old_even_parity": {region:summarize(values_old,r,mask) for region,mask in masks.items()},
            "new_odd_parity": {region:summarize(values_new,r,mask) for region,mask in masks.items()},
            "new_minus_old": {region:summarize(values_new-values_old,r,mask) for region,mask in masks.items()},
        }

    cells = {}
    for i in range(min(5,n)):
        cells[str(i)] = {
            "r":float(r[i]),"R":float(R[i]),
            "dR_dr_old_even":float(dR_old[i]),"dR_dr_new_odd":float(dR_new[i]),
            "chi_old_even":float(chi_old[i]),"chi_new_odd":float(chi_new[i]),
            "mass_old_even":float(mass_old[i]),"mass_new_odd":float(mass_new[i]),
            "chi_delta":float(chi_new[i]-chi_old[i]),
        }
    # Cell values i>=2 use the same stencil for both parities because no ghost
    # points enter there; their pointwise changes must be bitwise zero.
    interior_delta = float(np.max(np.abs(chi_new[2:]-chi_old[2:])))
    if interior_delta != 0.0:
        raise AssertionError(f"{path}: parity correction unexpectedly changed i>=2 chi by {interior_delta}")

    return {
        "resolution":resolution,"time":time,"checkpoint":str(path),
        "grid":{"n":n,"dr":dr,"r_max":n*dr},
        "chi_min_old":float(np.min(chi_old)),"chi_min_old_r":float(r[int(np.argmin(chi_old))]),
        "chi_min_new":float(np.min(chi_new)),"chi_min_new_r":float(r[int(np.argmin(chi_new))]),
        "marginal_roots_old":roots(r,chi_old),"marginal_roots_new":roots(r,chi_new),
        "chi_max_abs_delta_i_ge_2":interior_delta,
        "summary":summary,"first_five_cells":cells,
        "interpretation_guardrail":"This changes only how the odd radial areal radius is differentiated in invariant diagnostics. It does not alter the evolved state, H/M/C constraints, or source terms."
    }


def main():
    cases=[]
    for n in RESOLUTIONS:
        for t in TIMES:
            p=Path(f"baseline/N{n}/checkpoint_t{t:012.6f}.npz")
            if not p.is_file(): raise FileNotFoundError(p)
            cases.append(audit(p,n,t))
    out=Path("results/areal-radius-odd-parity-audit")
    out.mkdir(parents=True,exist_ok=True)
    report={
        "schema":"areal_radius_odd_parity_invariant_audit_v1",
        "status":"completed","kind":"saved_checkpoint_diagnostic_only",
        "source_commit":SOURCE_COMMIT,"campaign_run_id":CAMPAIGN_RUN_ID,
        "cases":cases,
        "summary_guardrail":"The previous even-parity derivative is not admissible at a regular center. Use odd parity for R; retain the corrected H/M constraint findings unchanged."
    }
    (out/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("AREAL_RADIUS_ODD_PARITY_AUDIT "+json.dumps({
        "status":report["status"],
        "cases":[{
            "N":c["resolution"],"t":c["time"],
            "chi_min_old":c["chi_min_old"],"chi_min_new":c["chi_min_new"],
            "chi_min_old_r":c["chi_min_old_r"],"chi_min_new_r":c["chi_min_new_r"],
            "roots_old":c["marginal_roots_old"],"roots_new":c["marginal_roots_new"],
            "dR_i0_old":c["first_five_cells"]["0"]["dR_dr_old_even"],
            "dR_i0_new":c["first_five_cells"]["0"]["dR_dr_new_odd"],
            "chi_i0_old":c["first_five_cells"]["0"]["chi_old_even"],
            "chi_i0_new":c["first_five_cells"]["0"]["chi_new_odd"],
            "chi_delta_i_ge_2":c["chi_max_abs_delta_i_ge_2"],
            "chi_delta_first_two_max":c["summary"]["chi"]["new_minus_old"]["first_two_cells"]["max_abs"]
        } for c in cases]},sort_keys=True))


if __name__=="__main__":
    main()
