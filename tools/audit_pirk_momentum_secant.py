#!/usr/bin/env python3
"""Post-process a saved PIRK stage capture into its exact discrete momentum secant.

This diagnostic imports no production kernel, calls no RHS, and advances no
state. It consumes report.json and stage_arrays.npz from the prior stage-capture
artifact and verifies the accepted split using archived operands only.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from typing import Any
import numpy as np

KAPPA = 8.0 * np.pi


def d1(values: np.ndarray, dr: float, parity: int = 1) -> np.ndarray:
    """Pinned SphericalCellGrid fourth-order first derivative."""
    if parity not in (-1, 1):
        raise ValueError("parity must be +1 or -1")
    u = np.asarray(values, dtype=float)
    n = u.size
    if u.ndim != 1 or n < 5:
        raise ValueError("values must be a 1-D array with at least five cells")
    result = np.empty_like(u)
    extended = np.empty(n + 4)
    extended[2:-2] = u
    extended[1] = parity * u[0]
    extended[0] = parity * u[1]
    for i in range(n - 2):
        q = extended[i : i + 5]
        result[i] = (q[0] - 8.0*q[1] + 8.0*q[3] - q[4]) / (12.0*dr)
    for i in (n - 2, n - 1):
        indices = np.arange(n - 5, n)
        offsets = indices.astype(float) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5)
        target[1] = 1.0
        result[i] = float(np.linalg.solve(matrix, target) @ u[indices]) / dr
    return result


def d2(values: np.ndarray, dr: float, parity: int = 1) -> np.ndarray:
    """Pinned SphericalCellGrid fourth-order second derivative."""
    if parity not in (-1, 1):
        raise ValueError("parity must be +1 or -1")
    u = np.asarray(values, dtype=float)
    n = u.size
    if u.ndim != 1 or n < 5:
        raise ValueError("values must be a 1-D array with at least five cells")
    result = np.empty_like(u)
    extended = np.empty(n + 4)
    extended[2:-2] = u
    extended[1] = parity * u[0]
    extended[0] = parity * u[1]
    for i in range(n - 2):
        q = extended[i : i + 5]
        result[i] = (-q[0] + 16.0*q[1] - 30.0*q[2] + 16.0*q[3] - q[4]) / (12.0*dr**2)
    for i in (n - 2, n - 1):
        indices = np.arange(n - 5, n)
        offsets = indices.astype(float) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5)
        target[2] = 2.0
        result[i] = float(np.linalg.solve(matrix, target) @ u[indices]) / dr**2
    return result


def momentum_constraint_geom(
    r: np.ndarray, b: np.ndarray, X: np.ndarray,
    Aa: np.ndarray, K: np.ndarray, dr: float,
) -> np.ndarray:
    """Accepted geometric M_r with A_b=-A_a/2 and chi_r=-D_1 X/(2X)."""
    chip = -0.5 * d1(X, dr, 1) / X
    return (
        d1(Aa, dr, 1) - (2.0/3.0)*d1(K, dr, 1)
        + 6.0*Aa*chip + 1.5*Aa*(2.0/r + d1(b, dr, 1)/b)
    )


def total_current(arrays: Any, prefix: str, r: np.ndarray, dr: float) -> np.ndarray:
    """Rebuild normalized six-sector j from the archived profiles."""
    g = lambda name: arrays[f"{prefix}__geometry__{name}"]
    s = lambda name: arrays[f"{prefix}__scalars__{name}"]
    sqrt_gamma = np.sqrt(g("a")) * g("b") * r**2 / g("X")**3
    j = (
        -s("PS")*d1(s("S"), dr, 1)
        -s("PD")*d1(s("D"), dr, 1)
        -s("Pi")*d1(s("phi"), dr, 1)/KAPPA
    )
    for species in ("dark_matter", "baryons", "radiation"):
        j += arrays[f"{prefix}__matter__{species}__momentum_r"] / sqrt_gamma / KAPPA
    return j


def momentum_constraint_total(
    arrays: Any, prefix: str, r: np.ndarray, dr: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    g = lambda name: arrays[f"{prefix}__geometry__{name}"]
    geom = momentum_constraint_geom(r, g("b"), g("X"), g("Aa"), g("K"), dr)
    j = total_current(arrays, prefix, r, dr)
    return geom - KAPPA*j, geom, j


def constraint_secant_components(
    r: np.ndarray, dr: float,
    b0: np.ndarray, b1: np.ndarray, X0: np.ndarray, X1: np.ndarray,
    A0: np.ndarray, A1: np.ndarray, K0: np.ndarray, K1: np.ndarray,
    delta_j: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return linear and exact nonlinear-cross terms in Delta M_r.

    Delta M = D1(Delta A)-2/3 D1(Delta K)
      + 6(Delta A chi0 + A0 Delta chi + Delta A Delta chi)
      + 3/2(Delta A q0 + A0 Delta q + Delta A Delta q)
      - kappa Delta j, with chi=-D1(X)/(2X), q=2/r+D1(b)/b.
    """
    dA, dK = A1-A0, K1-K0
    chi0, chi1 = -0.5*d1(X0, dr, 1)/X0, -0.5*d1(X1, dr, 1)/X1
    q0, q1 = 2.0/r + d1(b0, dr, 1)/b0, 2.0/r + d1(b1, dr, 1)/b1
    dchi, dq = chi1-chi0, q1-q0
    linear = (
        d1(dA, dr, 1) - (2.0/3.0)*d1(dK, dr, 1)
        + 6.0*(dA*chi0 + A0*dchi)
        + 1.5*(dA*q0 + A0*dq) - KAPPA*delta_j
    )
    nonlinear = 6.0*dA*dchi + 1.5*dA*dq
    return linear, nonlinear


def primary_l2_from_geometry(
    arrays: Any, prefix: str, r: np.ndarray, dr: float,
) -> dict[str, np.ndarray]:
    """Rebuild vendor primary L2 RHS at the archived final-primary geometry."""
    g = lambda name: arrays[f"{prefix}__geometry__{name}"]
    a,b,X,alpha,Lambda = (g(k) for k in ("a","b","X","alpha","Lambda"))
    ap,bp = d1(a,dr,1),d1(b,dr,1)
    app,bpp = d2(a,dr,1),d2(b,dr,1)
    chip = -0.5*d1(X,dr,1)/X
    chipp = -0.5*(d2(X,dr,1)/X - (d1(X,dr,1)/X)**2)
    Lambda_p = d1(Lambda,dr,-1)
    metric_lambda = (1.0-a/b)/r**2
    invrad = X**2/a
    radial = (
        0.5*app/a - a*Lambda_p - 0.75*(ap/a)**2 + 0.5*(bp/b)**2
        - 0.5*Lambda*ap + ap/(r*b) + 2.0*metric_lambda*(1.0+r*bp/b)
        + 4.0*chipp - 2.0*chip*(ap/a-bp/b-2.0/r)
    )
    scalar = (
        0.5*app/a + bpp/b - a*Lambda_p - (ap/a)**2 + 0.5*(bp/b)**2
        + 2.0*(3.0-a/b)*bp/(r*b) + 4.0*metric_lambda
        + 8.0*(chipp+chip**2) - 8.0*chip*(0.5*ap/a-bp/b-2.0/r)
    )
    Rrr, R = -invrad*radial, -invrad*scalar
    alphap, alphapp = d1(alpha,dr,1), d2(alpha,dr,1)
    DrDr_alpha = invrad*(alphapp-alphap*(0.5*ap/a+2.0*chip))
    lap_alpha = invrad*(alphapp-alphap*(0.5*ap/a-bp/b-2.0*chip-2.0/r))
    return {"Aa": -(DrDr_alpha-lap_alpha/3.0)+alpha*(Rrr-R/3.0), "K": -lap_alpha}


def _max_abs(values: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(values,dtype=float))))


def _region(values: np.ndarray, r: np.ndarray) -> dict[str, Any]:
    values = np.asarray(values,dtype=float)
    mask = r <= 10.0+1e-12
    loc = np.where(mask)[0][int(np.argmax(np.abs(values[mask])))]
    return {
        "max_abs_r_le_10": _max_abs(values[mask]),
        "max_index_r_le_10": int(loc),
        "radius_at_max_r_le_10": float(r[loc]),
        "signed_at_max_r_le_10": float(values[loc]),
        "cells_0_1_2_6": {str(i): float(values[i]) for i in (0,1,2,6) if i < values.size},
    }


def _verify_bundle(bundle: Path) -> dict[str,str]:
    sums_path = bundle/"SHA256SUMS.txt"
    if not sums_path.exists():
        raise FileNotFoundError(f"missing artifact checksum manifest: {sums_path}")
    verified: dict[str,str] = {}
    for line in sums_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, filename = line.split(maxsplit=1)
        filename = filename.strip().lstrip("*")
        actual = hashlib.sha256((bundle/filename).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"SHA256 mismatch for {filename}: {actual} != {expected}")
        verified[filename] = actual
    for name in ("report.json","stage_arrays.npz"):
        if name not in verified:
            raise ValueError(f"checksum manifest does not include {name}")
    return verified


def audit_bundle(bundle_dir: str | Path) -> dict[str,Any]:
    """Calculate endpoint M_r secant and reconstruct saved PIRK field increments."""
    bundle = Path(bundle_dir)
    verified = _verify_bundle(bundle)
    report = json.loads((bundle/"report.json").read_text(encoding="utf-8"))
    step = report["step"]
    dt = float(step["dt"])
    arrays = np.load(bundle/"stage_arrays.npz",allow_pickle=False)
    r = np.asarray(arrays["r"],dtype=float)
    if r.ndim != 1 or r.size < 5:
        raise ValueError("invalid radial profile in captured stage arrays")
    dr = float(r[1]-r[0])
    if not np.allclose(np.diff(r),dr,rtol=0.0,atol=64*np.finfo(float).eps*max(1.0,abs(dr))):
        raise ValueError("captured grid is not uniform")
    if not (dt>0.0 and np.isfinite(dt)):
        raise ValueError("invalid captured step dt")

    M0,geom0,j0 = momentum_constraint_total(arrays,"input",r,dr)
    M1,geom1,j1 = momentum_constraint_total(arrays,"copied_output",r,dr)
    j_stage0 = arrays["source_start_lapse_refreshed__j_total"]
    geom_stage0 = arrays["source_start_lapse_refreshed__M_r_geometry"]
    crosscheck = {
        "j_start_max_abs_difference": _max_abs(j0-j_stage0),
        "M_geom_start_max_abs_difference": _max_abs(geom0-geom_stage0),
    }
    deltaM, delta_geom, delta_j = M1-M0, geom1-geom0, j1-j0
    secant, geom_secant, j_secant = deltaM/dt, delta_geom/dt, delta_j/dt
    exact_sum_residual = deltaM-(delta_geom-KAPPA*delta_j)
    A0,A1 = arrays["input__geometry__Aa"],arrays["copied_output__geometry__Aa"]
    K0,K1 = arrays["input__geometry__K"],arrays["copied_output__geometry__K"]
    b0,b1 = arrays["input__geometry__b"],arrays["copied_output__geometry__b"]
    X0,X1 = arrays["input__geometry__X"],arrays["copied_output__geometry__X"]
    linear_delta, nonlinear_cross = constraint_secant_components(
        r,dr,b0,b1,X0,X1,A0,A1,K0,K1,delta_j
    )
    nonlinear_identity_residual = deltaM-linear_delta-nonlinear_cross

    # Reconstruct the archived stage assignment without invoking a kernel/RHS.
    start,predictor,final_l2 = (
        "stage_start_lapse_refreshed","stage_primary_predictor","lambda2_final_primary_L2"
    )
    l2_final = primary_l2_from_geometry(arrays,final_l2,r,dr)
    increment_residuals: dict[str,float] = {}
    mask = r <= 10.0+1e-12
    for name in ("a","b","X"):
        predicted = 0.5*dt*(
            arrays[f"{start}__term__explicit__{name}"]+
            arrays[f"{predictor}__term__explicit__{name}"]
        )
        actual = arrays[f"copied_output__geometry__{name}"]-arrays[f"input__geometry__{name}"]
        increment_residuals[name] = _max_abs((actual-predicted)[mask])
    for name in ("Aa","K"):
        predicted = 0.5*dt*(
            arrays[f"{start}__term__primary_l2__{name}"]+l2_final[name]+
            arrays[f"{start}__term__primary_l3__{name}"]+
            arrays[f"{predictor}__term__primary_l3__{name}"]
        )
        actual = arrays[f"copied_output__geometry__{name}"]-arrays[f"input__geometry__{name}"]
        increment_residuals[name] = _max_abs((actual-predicted)[mask])
    lambda_predicted = 0.5*dt*(
        arrays[f"{start}__term__lambda_l2"]+
        arrays[f"{final_l2}__term__lambda_l2"]+
        arrays[f"{start}__term__lambda_l3"]+
        arrays[f"{predictor}__term__lambda_l3"]
    )
    lambda_actual = arrays["copied_output__geometry__Lambda"]-arrays["input__geometry__Lambda"]
    increment_residuals["Lambda"] = _max_abs((lambda_actual-lambda_predicted)[mask])

    # Instantaneous geometric tangent at the refreshed start, not stale input alpha.
    g = lambda name: arrays[f"{start}__geometry__{name}"]
    rate_A = arrays[f"{start}__term__primary_l2__Aa"]+arrays[f"{start}__term__primary_l3__Aa"]
    rate_K = arrays[f"{start}__term__primary_l2__K"]+arrays[f"{start}__term__primary_l3__K"]
    rate_b,rate_X = arrays[f"{start}__term__explicit__b"],arrays[f"{start}__term__explicit__X"]
    bstart,Xstart,Astart = g("b"),g("X"),g("Aa")
    chi = -0.5*d1(Xstart,dr,1)/Xstart
    chi_rate = -0.5*(d1(rate_X,dr,1)/Xstart-d1(Xstart,dr,1)*rate_X/Xstart**2)
    q = 2.0/r+d1(bstart,dr,1)/bstart
    q_rate = d1(rate_b,dr,1)/bstart-d1(bstart,dr,1)*rate_b/bstart**2
    tangent = (
        d1(rate_A,dr,1)-(2.0/3.0)*d1(rate_K,dr,1)
        +6.0*rate_A*chi+6.0*Astart*chi_rate
        +1.5*rate_A*q+1.5*Astart*q_rate
    )
    alpha_change = g("alpha")-arrays["input__geometry__alpha"]
    return {
        "diagnostic_only": True,
        "no_rhs_evaluation_or_state_advance": True,
        "resolution": int(r.size),
        "dt": dt,
        "dr": dr,
        "input_time": float(step.get("input_time",float(step["copied_output_time"])-dt)),
        "copied_output_time": float(step["copied_output_time"]),
        "bundle_member_sha256": verified,
        "source_start_crosscheck": crosscheck,
        "max_abs_alpha_refresh_r_le_10": _max_abs(alpha_change[mask]),
        "stage_increment_residuals_max_abs_r_le_10": increment_residuals,
        "constraint_endpoint_max_abs_r_le_10": {
            "input": _max_abs(M0[mask]),"copied_output": _max_abs(M1[mask]),
        },
        "M_r_secant_per_time": _region(secant,r),
        "geometric_M_r_secant_per_time": _region(geom_secant,r),
        "total_j_secant_per_time": _region(j_secant,r),
        "M_r_geometric_tangent_at_refreshed_start": _region(tangent,r),
        "exact_deltaM_minus_deltaMgeom_plus_kappa_deltaj_max_abs_r_le_10": _max_abs(exact_sum_residual[mask]),
        "nonlinear_secant_cross_terms_per_time": _region(nonlinear_cross/dt,r),
        "nonlinear_secant_identity_max_abs_r_le_10": _max_abs(nonlinear_identity_residual[mask]/dt),
        "interpretation": (
            "The endpoint secant is the accepted finite-step operator. Brown's continuum identity is "
            "instantaneous and must be compared using refreshed-start/final-primary/predictor stage "
            "ownership; the old saved-state tangent is not the same starting RHS. This does not close "
            "the full Brown/source reconciliation."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle_dir",type=Path,help="unzipped SHA-verified archived stage bundle")
    parser.add_argument("--output",type=Path,help="optional JSON output path")
    args = parser.parse_args()
    result = audit_bundle(args.bundle_dir)
    rendered = json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(rendered,encoding="utf-8")
    print(rendered,end="")


if __name__ == "__main__":
    main()
