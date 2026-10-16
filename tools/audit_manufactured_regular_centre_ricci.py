#!/usr/bin/env python3
"""Analytic regular-centre Ricci oracle for Issue #4 (diagnostic only).

The target-state artifacts are not read. A fixed smooth regular metric is
prescribed analytically, so its continuum scalar spatial Ricci is known
without differentiating sampled data. It is compared with (a) the pinned BSSN
Ricci representation and (b) a metric-only areal-radius/proper-radius formula.
The copied D1/D2 stencil is the repository's pinned cell-centred fourth-order
stencil with parity ghosts. This tests its centre error class; it does not
verify the matter-coupled Bianchi identity or scientifically admit the model.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "27da39e1e4d3a37ce0f7ef70b228009323183048"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
DEFAULT_RESOLUTIONS = (80, 120, 160, 240, 320)
R_MAX = 40.0
AMPLITUDE = 0.1
LENGTH_SCALE = 8.0
ANNULUS = (5.0, 10.0)


def _finite(name: str, values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise FloatingPointError(f"non-finite values in {name}")
    return arr


def cell_derivative_fourth(values: np.ndarray, dr: float, parity: int) -> np.ndarray:
    """Pinned SphericalCellGrid D1, copied explicitly for audit transparency."""
    u = _finite("D1 input", values)
    if u.ndim != 1 or u.size < 5 or dr <= 0.0 or parity not in (-1, 1):
        raise ValueError("invalid fourth-order D1 inputs")
    n = u.size
    extended = np.empty(n + 4, dtype=np.float64)
    extended[2:-2] = u
    extended[1] = parity * u[0]
    extended[0] = parity * u[1]
    out = np.empty_like(u)
    for i in range(n - 2):
        s = extended[i:i + 5]
        out[i] = (s[0] - 8.0 * s[1] + 8.0 * s[3] - s[4]) / (12.0 * dr)
    indices = np.arange(n - 5, n)
    for i in (n - 2, n - 1):
        offsets = indices.astype(np.float64) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5, dtype=np.float64)
        target[1] = 1.0
        out[i] = float(np.linalg.solve(matrix, target) @ u[indices]) / dr
    return _finite("D1 output", out)


def cell_second_derivative_fourth(values: np.ndarray, dr: float, parity: int) -> np.ndarray:
    """Pinned SphericalCellGrid D2, copied explicitly for audit transparency."""
    u = _finite("D2 input", values)
    if u.ndim != 1 or u.size < 5 or dr <= 0.0 or parity not in (-1, 1):
        raise ValueError("invalid fourth-order D2 inputs")
    n = u.size
    extended = np.empty(n + 4, dtype=np.float64)
    extended[2:-2] = u
    extended[1] = parity * u[0]
    extended[0] = parity * u[1]
    out = np.empty_like(u)
    for i in range(n - 2):
        s = extended[i:i + 5]
        out[i] = (-s[0] + 16.0 * s[1] - 30.0 * s[2] + 16.0 * s[3] - s[4]) / (12.0 * dr**2)
    indices = np.arange(n - 5, n)
    for i in (n - 2, n - 1):
        offsets = indices.astype(np.float64) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5, dtype=np.float64)
        target[2] = 2.0
        out[i] = float(np.linalg.solve(matrix, target) @ u[indices]) / dr**2
    return _finite("D2 output", out)


def exact_profile(r: np.ndarray, amplitude: float, length_scale: float) -> dict[str, np.ndarray]:
    """Analytic regular profile with gamma_rr=1 and areal radius R=r*q(r)."""
    e = np.exp(-((r / length_scale) ** 2))
    q = 1.0 + amplitude * (1.0 - e)
    q_r = 2.0 * amplitude * r * e / length_scale**2
    q_rr = 2.0 * amplitude * e / length_scale**2 * (1.0 - 2.0 * r**2 / length_scale**2)
    areal = r * q
    areal_r = q + r * q_r
    areal_rr = 2.0 * q_r + r * q_rr
    # For dl^2 + R(l)^2 dOmega^2, R^(3)=-4 R_ll/R+2(1-R_l^2)/R^2.
    ricci_exact = -4.0 * areal_rr / areal + 2.0 * (1.0 - areal_r**2) / areal**2
    a = q ** (-4.0 / 3.0)
    b = q ** (2.0 / 3.0)
    X = q ** (-2.0 / 3.0)
    return {
        "q": _finite("q", q),
        "a": _finite("a", a),
        "b": _finite("b", b),
        "X": _finite("X", X),
        "areal_radius_exact": _finite("exact areal radius", areal),
        "proper_factor_exact": np.ones_like(r),
        "ricci_exact": _finite("analytic scalar Ricci", ricci_exact),
    }


def bssn_scalar_ricci(a: np.ndarray, b: np.ndarray, X: np.ndarray,
                      Lambda: np.ndarray, r: np.ndarray, dr: float) -> np.ndarray:
    """Pinned vendor scalar-Ricci expression, kept separate from metric Ricci."""
    d1 = lambda v, p=1: cell_derivative_fourth(v, dr, p)
    d2 = lambda v, p=1: cell_second_derivative_fourth(v, dr, p)
    ap, bp = d1(a, +1), d1(b, +1)
    app, bpp = d2(a, +1), d2(b, +1)
    chi_r = -0.5 * d1(X, +1) / X
    chi_rr = -0.5 * (d2(X, +1) / X - (d1(X, +1) / X) ** 2)
    lambda_r = d1(Lambda, -1)
    metric_lambda = (1.0 - a / b) / r**2
    inverse_physical_radial = X**2 / a
    scalar_bracket = (
        0.5 * app / a + bpp / b - a * lambda_r - (ap / a) ** 2
        + 0.5 * (bp / b) ** 2
        + 2.0 * (3.0 - a / b) * bp / (r * b)
        + 4.0 * metric_lambda + 8.0 * (chi_rr + chi_r**2)
        - 8.0 * chi_r * (0.5 * ap / a - bp / b - 2.0 / r)
    )
    return _finite("BSSN scalar Ricci", -inverse_physical_radial * scalar_bracket)


def observed_orders(rows: list[dict[str, Any]], key: str) -> list[dict[str, float]]:
    out = []
    for coarse, fine in zip(rows, rows[1:]):
        e0, e1 = abs(float(coarse[key])), abs(float(fine[key]))
        ratio_n = float(fine["N"]) / float(coarse["N"])
        p = math.log(e0 / e1) / math.log(ratio_n) if e0 > 0.0 and e1 > 0.0 else math.inf
        out.append({"from_N": int(coarse["N"]), "to_N": int(fine["N"]), "order": float(p)})
    return out


def analyze_resolution(n: int, *, r_max: float = R_MAX, amplitude: float = AMPLITUDE,
                       length_scale: float = LENGTH_SCALE,
                       annulus: tuple[float, float] = ANNULUS) -> dict[str, Any]:
    if n < 16 or r_max <= 0.0 or amplitude <= 0.0 or length_scale <= 0.0:
        raise ValueError("invalid analytic-oracle grid/profile parameters")
    if not 0.0 < annulus[0] < annulus[1] < r_max:
        raise ValueError("fixed physical annulus must lie inside the domain")
    dr = r_max / n
    r = (np.arange(n, dtype=np.float64) + 0.5) * dr
    p = exact_profile(r, amplitude, length_scale)
    a, b, X, exact = p["a"], p["b"], p["X"], p["ricci_exact"]

    # Reconstruct the connection target with the same pinned D1/parity contract
    # used by the initializer. This remains an evolved-connection representation witness.
    Lambda_metric = (
        cell_derivative_fourth(a, dr, +1) / (2.0 * a**2)
        - cell_derivative_fourth(b, dr, +1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    ricci_bssn = bssn_scalar_ricci(a, b, X, Lambda_metric, r, dr)

    # Independent geometric representation: direct metric quantities are
    # reconstructed from areal radius R=r*sqrt(b)/X and dl/dr=sqrt(a)/X.
    areal = r * np.sqrt(b) / X
    proper = np.sqrt(a) / X
    areal_r = cell_derivative_fourth(areal, dr, -1)
    areal_rr = cell_second_derivative_fourth(areal, dr, -1)
    proper_r = cell_derivative_fourth(proper, dr, +1)
    R_l = areal_r / proper
    R_ll = (areal_rr - areal_r * proper_r / proper) / proper**2
    ricci_metric = -4.0 * R_ll / areal + 2.0 * (1.0 - R_l**2) / areal**2

    # Deliberately invalid parity: areal radius is odd, not even. This must fail loudly.
    wrong_R_r = cell_derivative_fourth(areal, dr, +1)
    wrong_R_rr = cell_second_derivative_fourth(areal, dr, +1)
    wrong_R_l = wrong_R_r / proper
    wrong_R_ll = (wrong_R_rr - wrong_R_r * proper_r / proper) / proper**2
    wrong_ricci = -4.0 * wrong_R_ll / areal + 2.0 * (1.0 - wrong_R_l**2) / areal**2

    err_bssn = _finite("BSSN Ricci error", ricci_bssn - exact)
    err_metric = _finite("metric Ricci error", ricci_metric - exact)
    gap = _finite("BSSN-minus-metric Ricci gap", ricci_bssn - ricci_metric)
    wrong_err = _finite("wrong parity Ricci error", wrong_ricci - exact)
    ann_mask = (r >= annulus[0]) & (r <= annulus[1])
    if int(np.count_nonzero(ann_mask)) < 4:
        raise ValueError("annulus contains fewer than four cell centres")
    determinant_error = _finite("determinant error", a * b**2 - 1.0)
    if float(np.max(np.abs(determinant_error))) > 2.0e-14:
        raise RuntimeError("manufactured determinant-one profile failed")

    return {
        "N": int(n), "dr": float(dr), "r_cell0": float(r[0]),
        "analytic_ricci_cell0": float(exact[0]),
        "bssn_ricci_cell0": float(ricci_bssn[0]),
        "metric_ricci_cell0": float(ricci_metric[0]),
        "bssn_error_cell0": float(err_bssn[0]),
        "metric_error_cell0": float(err_metric[0]),
        "route_gap_cell0": float(gap[0]),
        "bssn_error_cell1": float(err_bssn[1]),
        "metric_error_cell1": float(err_metric[1]),
        "route_gap_cell1": float(gap[1]),
        "bssn_first2_max_abs_error": float(np.max(np.abs(err_bssn[:2]))),
        "metric_first2_max_abs_error": float(np.max(np.abs(err_metric[:2]))),
        "wrong_parity_cell0_error_abs": float(abs(wrong_err[0])),
        "bssn_fixed_annulus_max_abs_error": float(np.max(np.abs(err_bssn[ann_mask]))),
        "metric_fixed_annulus_max_abs_error": float(np.max(np.abs(err_metric[ann_mask]))),
        "route_gap_fixed_annulus_max_abs": float(np.max(np.abs(gap[ann_mask]))),
        "determinant_max_abs": float(np.max(np.abs(determinant_error))),
        "all_finite": bool(
            np.all(np.isfinite(ricci_bssn)) and np.all(np.isfinite(ricci_metric))
            and np.all(np.isfinite(wrong_ricci))
        ),
    }


def analyze(resolutions: tuple[int, ...] = DEFAULT_RESOLUTIONS, *, r_max: float = R_MAX,
            amplitude: float = AMPLITUDE, length_scale: float = LENGTH_SCALE,
            annulus: tuple[float, float] = ANNULUS) -> dict[str, Any]:
    if len(resolutions) < 3 or tuple(sorted(set(resolutions))) != tuple(resolutions):
        raise ValueError("provide at least three unique, strictly increasing resolutions")
    rows = [analyze_resolution(n, r_max=r_max, amplitude=amplitude,
                               length_scale=length_scale, annulus=annulus) for n in resolutions]
    orders = {
        "bssn_cell0_error": observed_orders(rows, "bssn_error_cell0"),
        "metric_cell0_error": observed_orders(rows, "metric_error_cell0"),
        "route_gap_cell0": observed_orders(rows, "route_gap_cell0"),
        "bssn_annulus_error": observed_orders(rows, "bssn_fixed_annulus_max_abs_error"),
        "metric_annulus_error": observed_orders(rows, "metric_fixed_annulus_max_abs_error"),
        "route_gap_annulus": observed_orders(rows, "route_gap_fixed_annulus_max_abs"),
    }
    # Bounds calibrate this fixed MMS family; they are not target-state acceptance tolerances.
    centre_bssn = [x["order"] for x in orders["bssn_cell0_error"]]
    centre_metric = [x["order"] for x in orders["metric_cell0_error"]]
    centre_gap = [x["order"] for x in orders["route_gap_cell0"]]
    annulus_bssn = [x["order"] for x in orders["bssn_annulus_error"]]
    annulus_metric = [x["order"] for x in orders["metric_annulus_error"]]
    all_finite = all(row["all_finite"] for row in rows)
    negative_control = all(
        row["wrong_parity_cell0_error_abs"] > 1.0e3 * row["metric_first2_max_abs_error"]
        for row in rows
    )
    contract = bool(
        all_finite and negative_control
        and all(3.5 <= x <= 4.5 for x in centre_bssn)
        and all(1.7 <= x <= 2.3 for x in centre_metric)
        and all(1.7 <= x <= 2.3 for x in centre_gap)
        and all(3.5 <= x <= 4.5 for x in annulus_bssn)
        and all(3.5 <= x <= 4.5 for x in annulus_metric)
    )
    if not contract:
        raise RuntimeError(f"manufactured regular-centre calibration contract failed: {orders}")
    return {
        "schema": "issue4_manufactured_regular_centre_ricci_oracle_v1",
        "status": "completed_diagnostic_only",
        "diagnostic_contract_passed": contract,
        "scientific_admission": False,
        "production_evolution_steps": 0,
        "target_state_artifacts_read": False,
        "production_or_vendor_files_modified": False,
        "source_anchor": SOURCE_COMMIT,
        "vendor_anchor": VENDOR_COMMIT,
        "profile": {
            "q": "1 + A*(1-exp(-(r/L)^2))",
            "A": float(amplitude), "L": float(length_scale), "r_max": float(r_max),
            "resolutions": [int(x) for x in resolutions],
            "annulus": [float(annulus[0]), float(annulus[1])],
            "metric": "a=q^(-4/3), b=q^(2/3), X=q^(-2/3); a*b^2=1, gamma_rr=1, R_areal=r*q",
            "continuum_Ricci": "-4*R''/R + 2*(1-(R')^2)/R^2 for dl=dr and areal radius R(r)",
        },
        "rows": rows,
        "descriptive_orders_on_manufactured_profile": orders,
        "interpretation": (
            "For this fixed analytic regular metric, the BSSN scalar-Ricci error at cell 0 is fourth-order, "
            "while the metric-only areal/proper-radius reconstruction and their route gap are second-order at "
            "the first cell; both routes are approximately fourth-order on the fixed physical annulus. "
            "The incorrect even parity for the odd areal radius is strongly rejected. This is a calibrated "
            "manufactured-profile discretization result only: it does not identify the origin of PR #283's "
            "target-profile gap, prove the matter-coupled propagation identity, validate total stress-energy "
            "conservation, or admit scientific correctness."
        ),
    }


def _git(path: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(path), *args], text=True).strip()


def verify_pins() -> dict[str, str]:
    source_root = Path(os.environ.get("COSMOS_PINNED_SOURCE_ROOT", ROOT / "source")).resolve()
    vendor_rel = Path("vendor/bb-palatini-unified-r0")
    if _git(source_root, "rev-parse", "HEAD") != SOURCE_COMMIT:
        raise RuntimeError("pinned production source checkout mismatch")
    if _git(source_root / vendor_rel, "rev-parse", "HEAD") != VENDOR_COMMIT:
        raise RuntimeError("pinned source vendor revision mismatch")
    if _git(ROOT / vendor_rel, "rev-parse", "HEAD") != VENDOR_COMMIT:
        raise RuntimeError("diagnostic checkout vendor revision mismatch")
    production_paths = (
        "engine/production_kernel.py", "engine/matter_system.py", "engine/matter_rhs.py",
        "engine/scalar_system.py", "engine/v55_matter.py", "engine/v55_pirk_adapter.py",
        "engine/cmc_gauge.py", "engine/v55_initial.py", "engine/stress_energy.py",
        "engine/valencia.py", "engine/worldtube.py",
    )
    for path in production_paths:
        if _git(ROOT, "rev-parse", f"HEAD:{path}") != _git(source_root, "rev-parse", f"HEAD:{path}"):
            raise RuntimeError(f"diagnostic checkout production source mismatch: {path}")
    return {
        "production_source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "production_blob_scope_checked": ",".join(production_paths),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=list(DEFAULT_RESOLUTIONS))
    parser.add_argument("--r-max", type=float, default=R_MAX)
    parser.add_argument("--amplitude", type=float, default=AMPLITUDE)
    parser.add_argument("--length-scale", type=float, default=LENGTH_SCALE)
    parser.add_argument("--annulus", nargs=2, type=float, default=list(ANNULUS))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    pins = verify_pins()
    report = analyze(tuple(args.resolutions), r_max=args.r_max, amplitude=args.amplitude,
                     length_scale=args.length_scale, annulus=tuple(args.annulus))
    report["provenance"] = pins
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "schema": report["schema"], "status": report["status"],
        "diagnostic_contract_passed": report["diagnostic_contract_passed"],
        "scientific_admission": report["scientific_admission"],
        "source_anchor": pins["production_source_commit"],
        "vendor_anchor": pins["vendor_commit"],
        "centre_gap_orders": report["descriptive_orders_on_manufactured_profile"]["route_gap_cell0"],
        "annulus_gap_orders": report["descriptive_orders_on_manufactured_profile"]["route_gap_annulus"],
        "output": str(args.output),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
