#!/usr/bin/env python3
"""Diagnostic-only source-pinned test of the analytic origin half-cell s2 coefficient.

This script compares the existing fitted source-moment controls against a new
fixed-s0/fixed-s2 arm on copied initializer recurrences. It performs no
production evolution and edits no engine, vendor, dependency, or campaign file.
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_pirk_adapter import vendor_modules
from engine.scalar_system import cosmos_potential
from engine.v55_matter import RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT
from tools.audit_current_main_centre_proper_radial_ricci import (
    _state_digest,
    proper_radial_ricci,
)

SOURCE_COMMIT = "27da39e1e4d3a37ce0f7ef70b228009323183048"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_BLOBS = {
    "engine/production_kernel.py": "9685e76bb5a6f26d6bf383e240de24f37523f265",
    "engine/matter_system.py": "e52a117690253ddad7b24d116145a65474a351d2",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/v55_pirk_adapter.py": "eeee0492d80889d70e52b0946563eee5054c0f11",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
}
RESOLUTIONS = (160, 240, 320)
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
FIT_RADIUS = 3.5
TOLERANCE = 2.0e-14
MAX_ITERATIONS = 100
OUT = Path("results/centre-ricci-analytic-s2-halfcell-2026-10-16")
VARIANTS = (
    ("fit_degree3_r2p5", "fitted", 3, 2.5),
    ("fit_degree4_r3", "fitted", 4, 3.0),
    ("fit_degree5_r3p5", "fitted", 5, 3.5),
    ("fixed_s0_s2_fit_s4_s6_s8", "fixed_s2", 3, FIT_RADIUS),
    ("parity_composite_simpson", "parity_simpson", 0, FIT_RADIUS),
)


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def _max_abs(values) -> float:
    a = np.asarray(values, dtype=np.float64)
    return float(np.max(np.abs(a))) if a.size else 0.0


def _finite(values) -> bool:
    return bool(np.all(np.isfinite(np.asarray(values, dtype=np.float64))))


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _verify_source_pins() -> dict:
    source_head = _git("rev-parse", "HEAD")
    source_vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if source_vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor pin mismatch: expected {VENDOR_COMMIT}, found {source_vendor}")
    pinned = {}
    for path, expected in SOURCE_BLOBS.items():
        base_blob = _git("rev-parse", f"{SOURCE_COMMIT}:{path}")
        head_blob = _git("rev-parse", f"HEAD:{path}")
        match = base_blob == expected and head_blob == expected
        pinned[path] = {
            "expected_blob": expected,
            "source_commit_blob": base_blob,
            "diagnostic_head_blob": head_blob,
            "matches": match,
        }
        if not match:
            raise RuntimeError(f"source blob mismatch for {path}: {pinned[path]}")
    changed = _git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--",
        "engine", "vendor", ".gitmodules", "requirements.txt",
    )
    if changed:
        raise RuntimeError("production/vendor/dependency paths changed:\n" + changed)
    return {
        "diagnostic_head": source_head,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": source_vendor,
        "source_blobs": pinned,
        "engine_vendor_dependency_diff_empty": True,
    }


def _fixed_s2_fit(source: np.ndarray, r: np.ndarray, r_fit: float) -> dict:
    """Fit only s4/s6/s8 after imposing analytically known s0=1, s2=-4*pi*A^2."""
    source = np.asarray(source, dtype=np.float64)
    r = np.asarray(r, dtype=np.float64)
    s0 = 1.0
    s2 = -4.0 * np.pi * AMPLITUDE**2
    mask = r <= r_fit + 16.0 * np.finfo(float).eps * r_fit
    rr = r[mask]
    if rr.size < 5:
        raise ValueError(f"fixed-s2 fit requires at least five cells; got {rr.size}")
    x = rr**2
    z = x / (r_fit * r_fit)
    y = (source[mask] - s0 - s2 * x) / (x * x)
    design = np.column_stack((np.ones_like(z), z, z * z))
    scaled, *_ = np.linalg.lstsq(design, y, rcond=None)
    s4 = float(scaled[0])
    s6 = float(scaled[1] / (r_fit * r_fit))
    s8 = float(scaled[2] / (r_fit**4))
    predicted = s0 + s2 * x + s4 * x**2 + s6 * x**3 + s8 * x**4
    residual = source[mask] - predicted
    loo = []
    for omitted in range(rr.size):
        keep = np.arange(rr.size) != omitted
        coeff, *_ = np.linalg.lstsq(design[keep], y[keep], rcond=None)
        loo.append(float(y[omitted] - design[omitted] @ coeff))
    return {
        "method": "fix s0=1 and s2=-4*pi*A^2; fit only s4,s6,s8",
        "fit_r_max": float(r_fit),
        "fit_cells": int(rr.size),
        "s0_fixed": s0,
        "s2_fixed": float(s2),
        "s4_fitted": s4,
        "s6_fitted": s6,
        "s8_fitted": s8,
        "fit_max_abs_source_residual": _max_abs(residual),
        "fit_rms_source_residual": float(np.sqrt(np.mean(residual * residual))),
        "loo_scaled_ratio_r4_rms": float(np.sqrt(np.mean(np.square(loo)))),
        "loo_scaled_ratio_r4_max_abs": _max_abs(loo),
        "coefficients_finite": bool(np.all(np.isfinite([s4, s6, s8]))),
    }


def _fitted_even_polynomial(source: np.ndarray, r: np.ndarray, r0: float,
                            r_fit: float, degree: int) -> dict:
    """Fit an unconstrained even source polynomial for retained comparison arms."""
    mask = r <= r_fit + 16.0 * np.finfo(float).eps * r_fit
    rr = r[mask]
    if rr.size < degree + 2:
        raise ValueError(
            f"degree-{degree} comparison fit at r<={r_fit} needs degree+2 cells; got {rr.size}"
        )
    z = (rr / r_fit) ** 2
    design = np.column_stack([z**power for power in range(degree + 1)])
    scaled, *_ = np.linalg.lstsq(design, source[mask], rcond=None)
    coeff = np.asarray(
        [scaled[power] / (r_fit ** (2 * power)) for power in range(degree + 1)],
        dtype=np.float64,
    )
    residual = source[mask] - design @ scaled
    loo = []
    for omitted in range(rr.size):
        keep = np.arange(rr.size) != omitted
        loo_coeff, *_ = np.linalg.lstsq(design[keep], source[mask][keep], rcond=None)
        loo.append(float(source[mask][omitted] - design[omitted] @ loo_coeff))
    coefficient_names = [f"s{2 * power}" for power in range(degree + 1)]
    return {
        "method": "unconstrained even-polynomial fit",
        "fit_r_max": float(r_fit),
        "fit_cells": int(rr.size),
        "degree_in_r_squared": int(degree),
        "coefficient_names": coefficient_names,
        "coefficients_even_powers": [float(v) for v in coeff],
        "fit_max_abs_source_residual": _max_abs(residual),
        "fit_rms_source_residual": float(np.sqrt(np.mean(residual * residual))),
        "loo_source_rms_residual": float(np.sqrt(np.mean(np.square(loo)))),
        "loo_source_max_abs_residual": _max_abs(loo),
    }


def _moment_from_fit(fit: dict, r0: float, mode: str) -> tuple[float, dict]:
    if mode == "fixed_s2":
        coeff = [
            fit["s0_fixed"], fit["s2_fixed"], fit["s4_fitted"],
            fit["s6_fitted"], fit["s8_fitted"],
        ]
    else:
        coeff = fit["coefficients_even_powers"]
    coefficient_names = (
        ["s0", "s2", "s4", "s6", "s8"]
        if mode == "fixed_s2" else fit["coefficient_names"]
    )
    moment = float(sum(
        value * r0 ** (2 * power + 1) / (2 * power + 1)
        for power, value in enumerate(coeff)
    ))
    return moment, {
        "integrated_coefficient_names": coefficient_names,
        "integrated_coefficients_even_powers": [float(v) for v in coeff],
        "half_cell_integral": moment,
    }


def _make_fit(source: np.ndarray, r: np.ndarray, r0: float,
              mode: str, degree: int, fit_r_max: float) -> tuple[float, dict]:
    if mode == "fixed_s2":
        fit = _fixed_s2_fit(source, r, fit_r_max)
    else:
        fit = _fitted_even_polynomial(source, r, r0, fit_r_max, degree)
    moment, meta = _moment_from_fit(fit, r0, mode)
    point_sample_moment = float(r0 * source[0])
    meta.update(fit)
    meta.update({
        "point_sample_half_cell_integral": point_sample_moment,
        "delta_integral_vs_point_sample": float(moment - point_sample_moment),
        "analytic_s0_s2_leading_delta": float(
            r0 + (-4.0 * np.pi * AMPLITUDE**2) * r0**3 / 3.0
            - point_sample_moment
        ),
        "fit_residual_and_coefficients_finite": _finite(
            [moment, point_sample_moment, *meta["integrated_coefficients_even_powers"]]
        ),
    })
    return moment, meta


def _source_for_B(state, B: np.ndarray) -> np.ndarray:
    grid, fields = state.grid, state.scalars
    r = np.asarray(grid.centers, dtype=np.float64)
    d1 = lambda values: grid.cell_derivative_fourth(values, parity=1)
    S, PS = np.asarray(fields.S), np.asarray(fields.PS)
    D, PD = np.asarray(fields.D), np.asarray(fields.PD)
    phi, Pi = np.asarray(fields.phi), np.asarray(fields.Pi)
    Sp, Dp = d1(S), d1(D)
    current = -(PD * Dp)
    rhoS = 0.5 * (PS * PS + B * Sp * Sp) + 0.5 * S * S
    rhoD = 0.5 * (PD * PD + B * Dp * Dp) - 0.5 * D * D
    rho_phi = (0.5 * Pi * Pi + cosmos_potential(phi)) / (8.0 * np.pi)
    rho_matter = (
        float(RHO_DM_PRESENT) + float(RHO_B_PRESENT) + float(RHO_R_PRESENT)
    ) / (8.0 * np.pi)
    rho = rhoS + rhoD + rho_phi + rho_matter
    phi0, pi0 = 0.179055, 0.00341913
    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += float(RHO_DM_PRESENT) + float(RHO_B_PRESENT) + float(RHO_R_PRESENT)
    H0 = np.sqrt(max(rho0 / 3.0, 0.0))
    Kt = np.full_like(r, -H0)
    Ktp = d1(Kt)
    Kr = Kt + r * Ktp + 4.0 * np.pi * r * current
    source = 1.0 - 8.0 * np.pi * r * r * rho + 2.0 * r * r * Kr * Kt + r * r * Kt * Kt
    if not _finite(source):
        raise FloatingPointError("nonfinite source in copied initializer recurrence")
    return source


def _integrate_point_sample(source: np.ndarray, dr: float, r: np.ndarray) -> np.ndarray:
    r0 = float(r[0])
    return np.cumsum(np.r_[r0 * source[0], 0.5 * dr * (source[1:] + source[:-1])]) / r


def _parity_composite_simpson_integral(
    source: np.ndarray, r: np.ndarray, dr: float, origin_value: float = 1.0
) -> np.ndarray:
    """Integrate centre samples through a parity-even Simpson reconstruction.

    For a target centre r_i=(i+1/2)dr, construct the even extension on
    [-r_i,r_i], integrate with Simpson 1/3 panels plus a final Simpson 3/8
    panel, and halve the symmetric integral. The first half-cell uses the
    known origin value and Simpson's two-panel formula. This is diagnostic
    only; the production initializer remains unchanged.
    """
    samples = np.asarray(source, dtype=np.float64)
    integral = np.zeros_like(samples)
    for i in range(samples.size):
        if i == 0:
            r0 = float(r[0])
            integral[i] = r0*(samples[i] + 2.0*origin_value)/3.0
            continue
        signed_indices = np.arange(-i-1, i+1, dtype=int)
        sample_indices = np.where(signed_indices < 0, -signed_indices-1, signed_indices)
        values = samples[sample_indices]
        intervals = 2*i + 1
        simpson13_intervals = intervals - 3
        symmetric_integral = 0.0
        if simpson13_intervals > 0:
            first = values[:simpson13_intervals+1]
            symmetric_integral += (dr/3.0)*(
                first[0] + first[-1]
                + 4.0*np.sum(first[1:-1:2])
                + 2.0*np.sum(first[2:-1:2])
            )
        start = simpson13_intervals
        symmetric_integral += (3.0*dr/8.0)*(
            values[start] + 3.0*values[start+1]
            + 3.0*values[start+2] + values[start+3]
        )
        integral[i] = 0.5*symmetric_integral
    return integral


def _integrate_variant(source: np.ndarray, r: np.ndarray, dr: float,
                       mode: str, degree: int, fit_r_max: float) -> tuple[np.ndarray, dict]:
    r0 = float(r[0])
    point_sample_moment = float(r0*source[0])
    if mode == "parity_simpson":
        integral = _parity_composite_simpson_integral(source, r, dr, origin_value=1.0)
        Bnext = integral/r
        moment = float(integral[0])
        fit_meta = {
            "method": "parity-even-extension composite Simpson 1/3 + 3/8",
            "half_cell_integral": moment,
            "point_sample_half_cell_integral": point_sample_moment,
            "delta_integral_vs_point_sample": float(moment-point_sample_moment),
            "analytic_s0_s2_leading_delta": float(
                r0 + (-4.0*np.pi*AMPLITUDE**2)*r0**3/3.0 - point_sample_moment
            ),
            "fit_applicable": False,
            "fit_max_abs_source_residual": None,
            "fit_rms_source_residual": None,
            "integrated_coefficient_names": [],
            "integrated_coefficients_even_powers": [],
            "coefficients_finite": bool(_finite([moment, point_sample_moment])),
        }
    else:
        moment, fit_meta = _make_fit(source, r, r0, mode, degree, fit_r_max)
        Bnext = np.cumsum(np.r_[moment, 0.5*dr*(source[1:]+source[:-1])])/r
    fit_meta["iteration_B_positive_finite"] = bool(_finite(Bnext) and np.all(Bnext > 0.0))
    return Bnext, fit_meta


def _fit_rows_for_report(source: np.ndarray, r: np.ndarray) -> dict:
    rows = {}
    for tag, mode, degree, fit_r_max in VARIANTS:
        if mode == "parity_simpson":
            rows[tag] = {
                "method": "no source polynomial fit; parity-aware composite Simpson quadrature",
                "fit_applicable": False,
                "fit_max_abs_source_residual": None,
                "fit_rms_source_residual": None,
                "coefficients_finite": True,
            }
        else:
            _, fit = _make_fit(source, r, float(r[0]), mode, degree, fit_r_max)
            rows[tag] = fit
    return rows


def _matter_subtracted_constraint_profiles(
    raw_constraints: dict,
    total_rho: np.ndarray,
    total_j: np.ndarray,
) -> dict[str, np.ndarray]:
    """Separate geometric constraint sides from matter-subtracted residuals.

    The pinned vendor's vacuum.constraints() returns the geometry-only sides.
    In matter-filled data these are not expected to vanish by themselves:
    the total ADM constraints subtract 16*pi*rho and 8*pi*j using the canonical
    production stress-energy projection. Record both operands and residuals;
    this is not an independent stress-energy-conservation proof.
    """
    h_metric = np.asarray(raw_constraints["hamiltonian"], dtype=np.float64)
    m_metric = np.asarray(raw_constraints["momentum"], dtype=np.float64)
    rho = np.asarray(total_rho, dtype=np.float64)
    j = np.asarray(total_j, dtype=np.float64)
    if h_metric.shape != rho.shape or m_metric.shape != j.shape:
        raise ValueError("metric constraints and projected matter sources must share the grid shape")
    if not all(_finite(v) for v in (h_metric, m_metric, rho, j)):
        raise FloatingPointError("non-finite metric constraint or projected matter source")
    h_source = 16.0 * np.pi * rho
    m_source = 8.0 * np.pi * j
    return {
        "matter_hamiltonian_source_16pi_rho": h_source,
        "hamiltonian_total_residual": h_metric - h_source,
        "matter_momentum_source_8pi_j": m_source,
        "momentum_total_residual": m_metric - m_source,
    }


def _metric_witness(grid, state, vacuum, B: np.ndarray) -> tuple[dict, dict[str, np.ndarray]]:
    original_geometry = state.geometry
    r = np.asarray(grid.centers, dtype=np.float64)
    A = 1.0 / np.sqrt(B)
    a = A ** (4.0 / 3.0)
    b = A ** (-2.0 / 3.0)
    X = A ** (-1.0 / 3.0)
    d1 = lambda values: grid.cell_derivative_fourth(values, parity=1)
    Lambda = d1(a) / (2.0 * a * a) - d1(b) / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)
    copied = copy.deepcopy(original_geometry)
    copied.a, copied.b, copied.X, copied.Lambda = a.copy(), b.copy(), X.copy(), Lambda.copy()
    copied.B = 0.75 * Lambda
    terms = vacuum.geometry_terms(grid, copied)
    vendor_tf = np.asarray(terms["Rr_r"] - terms["R"] / 3.0, dtype=np.float64)
    direct = proper_radial_ricci(grid, r * np.sqrt(b) / X, np.sqrt(a) / X)
    direct_tf = np.asarray(direct["tf_ricci"], dtype=np.float64)
    gap = vendor_tf - direct_tf
    raw_constraints = vacuum.constraints(grid, copied)
    constraints = {
        name: np.asarray(raw_constraints[name], dtype=np.float64)
        for name in ("hamiltonian", "momentum", "connection")
        if name in raw_constraints
    }
    total = assemble_total_stress_energy(
        grid, copied, state.scalars, state.matter
    )
    constraints.update(
        _matter_subtracted_constraint_profiles(raw_constraints, total.rho, total.j)
    )
    determinant_constraint = np.asarray(raw_constraints["determinant"], dtype=np.float64)
    metric_only_names = ("hamiltonian", "momentum", "connection")
    total_residual_names = ("hamiltonian_total_residual", "momentum_total_residual")
    matter_source_names = ("matter_hamiltonian_source_16pi_rho", "matter_momentum_source_8pi_j")
    annulus = (r >= 5.0) & (r <= 10.0)
    centre = r <= 2.0
    invariants = {
        "a_times_b_squared_minus_one_max_abs": _max_abs(a * b * b - 1.0),
        "b_minus_X_squared_max_abs": _max_abs(b - X * X),
        "a_minus_X_minus_four_max_abs": _max_abs(a - X ** (-4.0)),
        "X6_minus_B_max_abs": _max_abs(X**6 - B),
        "areal_radius_minus_r_max_abs": _max_abs(r * np.sqrt(b) / X - r),
        "all_metric_profiles_finite": all(_finite(v) for v in (B, a, b, X, Lambda, vendor_tf, direct_tf, gap)),
    }
    witness = {
        "cell0": float(gap[0]),
        "cell1": float(gap[1]),
        "first_two_max_abs": _max_abs(gap[:2]),
        "centre_r_le_2_max_abs": _max_abs(gap[centre]) if np.any(centre) else None,
        "fixed_annulus_r5_to10_max_abs": _max_abs(gap[annulus]) if np.any(annulus) else None,
        "fixed_annulus_r5_to10_signed_at_max": float(gap[annulus][np.argmax(np.abs(gap[annulus]))]) if np.any(annulus) else None,
        "constraint_max_abs": {k: _max_abs(constraints[k]) for k in metric_only_names},
        "constraint_first_two": {k: [float(x) for x in constraints[k][:2]] for k in metric_only_names},
        "matter_constraint_source_max_abs": {k: _max_abs(constraints[k]) for k in matter_source_names},
        "matter_constraint_source_first_two": {k: [float(x) for x in constraints[k][:2]] for k in matter_source_names},
        "total_constraint_residual_max_abs": {k: _max_abs(constraints[k]) for k in total_residual_names},
        "total_constraint_residual_first_two": {k: [float(x) for x in constraints[k][:2]] for k in total_residual_names},
        "determinant_constraint_max_abs": _max_abs(determinant_constraint),
        "determinant_constraint_first_two": [float(x) for x in determinant_constraint[:2]],
        "constraint_interpretation": {
            "hamiltonian": "metric-only vacuum Hamiltonian side; matter source 16*pi*rho is reported separately",
            "momentum": "metric-only vacuum momentum side; matter source 8*pi*j is reported separately",
            "total_residuals": "metric-only side minus canonical projected matter source; not a stress-energy conservation test",
        },
        "invariants": invariants,
    }
    arrays = {
        "B": B.copy(), "a": a, "b": b, "X": X, "Lambda": Lambda,
        "vendor_tf_ricci": vendor_tf, "direct_tf_ricci": direct_tf,
        "ricci_gap": gap, "ricci_gap_over_r2": gap / (r*r),
        **{f"constraint_{k}": v for k, v in constraints.items()},
    }
    return witness, arrays


def _run_variant(state, vacuum, tag: str, mode: str, degree: int, fit_r_max: float) -> tuple[dict, dict[str, np.ndarray]]:
    grid, r = state.grid, np.asarray(state.grid.centers, dtype=np.float64)
    dr = float(grid.dr)
    state_before = _state_digest(state)
    B_baseline = np.ones_like(r)
    baseline_updates = []
    for _ in range(6):
        Bnext = _integrate_point_sample(_source_for_B(state, B_baseline), dr, r)
        baseline_updates.append(_max_abs(Bnext - B_baseline))
        B_baseline = Bnext
    B_prod = np.asarray(state.geometry.X, dtype=np.float64) ** 6
    baseline_error = _max_abs(B_baseline - B_prod)
    baseline_scale = max(1.0, _max_abs(B_prod))
    if baseline_error > 2e-12 * baseline_scale:
        raise RuntimeError(f"{tag}: copied point-sampled recurrence misses production X^6 by {baseline_error:.6e}")

    B = np.ones_like(r)
    history = []
    B_six = None
    converged = False
    last_update = float("inf")
    last_fit = None
    for iteration in range(1, MAX_ITERATIONS + 1):
        source = _source_for_B(state, B)
        Bnext, meta = _integrate_variant(source, r, dr, mode, degree, fit_r_max)
        last_update = _max_abs(Bnext - B)
        history.append({
            "iteration": iteration,
            "max_abs_B_update": float(last_update),
            "half_cell_delta_integral": meta["delta_integral_vs_point_sample"],
            "fit_max_abs_source_residual": meta["fit_max_abs_source_residual"],
            "fit_rms_source_residual": meta["fit_rms_source_residual"],
        })
        B = Bnext
        last_fit = meta
        if iteration == 6:
            B_six = B.copy()
        if iteration >= 6 and last_update <= TOLERANCE * max(1.0, _max_abs(B)):
            converged = True
            break
    if B_six is None:
        raise RuntimeError(f"{tag}: copied recurrence did not reach its six-step witness")
    source_fixed = _source_for_B(state, B)
    B_check, check_meta = _integrate_variant(source_fixed, r, dr, mode, degree, fit_r_max)
    fixed_point_defect = _max_abs(B_check - B)
    state_after = _state_digest(state)
    input_unchanged = state_after == state_before
    if not converged:
        raise RuntimeError(f"{tag}: fixed-point iteration did not converge in {MAX_ITERATIONS} iterations")
    fit_scan = _fit_rows_for_report(source_fixed, r)
    fixed_witness, fixed_arrays = _metric_witness(grid, state, vacuum, B)
    six_witness, six_arrays = _metric_witness(grid, state, vacuum, B_six)
    fixed_witness.update({
        "convergence": {
            "converged": converged,
            "iterations_used": len(history),
            "last_B_update_max_abs": float(last_update),
            "fixed_point_defect_max_abs": float(fixed_point_defect),
            "baseline_point_sample_six_step_reconstruction_max_abs": float(baseline_error),
            "baseline_six_step_max_abs_updates": [float(x) for x in baseline_updates],
        },
        "half_cell": {
            "applied_method": last_fit["method"],
            "applied_delta_integral": last_fit["delta_integral_vs_point_sample"],
            "analytic_s0_s2_leading_delta": last_fit["analytic_s0_s2_leading_delta"],
            "applied_integral": last_fit["half_cell_integral"],
            "point_sample_integral": last_fit["point_sample_half_cell_integral"],
            "selected_fit": last_fit,
            "source_fit_family_controls_at_fixed_point": fit_scan,
            "fixed_point_check_fit": check_meta,
            "iteration_history": history,
        },
        "six_step_witness": six_witness,
        "input_state_unchanged": input_unchanged,
        "status": "completed_diagnostic_only",
    })
    arrays_out = {
        **{f"fixed_{k}": v for k, v in fixed_arrays.items()},
        **{f"six_{k}": v for k, v in six_arrays.items()},
        "production_B": B_prod.copy(),
        "source_fixed_point": source_fixed.copy(),
    }
    return fixed_witness, arrays_out


def main() -> int:
    provenance = _verify_source_pins()
    _, vacuum, _ = vendor_modules()
    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "schema": "issue4_centre_ricci_analytic_s2_halfcell_v1",
        "status": "completed_diagnostic_only",
        "scientific_admission": False,
        "source_and_vendor_pins": provenance,
        "configuration": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "amplitude": AMPLITUDE,
            "width": WIDTH,
            "D_amplitude": D_AMPLITUDE,
            "fit_radius_fixed_s2": FIT_RADIUS,
            "fixed_analytic_coefficients": {
                "s0": 1.0,
                "s2": float(-4.0 * np.pi * AMPLITUDE**2),
                "derivation_scope": "standard uncharged initializer; not generalized to optional charged dust or radial EM",
            },
            "evolution_steps": 0,
            "production_state_modified": False,
        },
        "resolutions": {},
        "checks": {},
        "interpretation_limits": [
            "The exact s2 coefficient is derived for the declared standard uncharged initializer and its current archive normalization.",
            "The metric/Ricci calculation uses the same pinned cell-centred derivative family; it is an independent formula route, not an independent discretization.",
            "A small fit residual or reduced centre gap does not select the continuum solution or authorize a production patch.",
            "The raw vacuum.constraints Hamiltonian and momentum values are metric-only sides; matter-subtracted total residuals are separately reported using the canonical local projection.",
            "These same-path matter-subtracted residuals are initial-data witnesses, not independent total stress-energy conservation or source-normalization proof.",
            "The fitted s4/s6/s8 coefficients are model-dependent and are reported with leave-one-out diagnostics.",
            "No production/vendor/dependency file or acceptance tolerance is changed; scientific admission remains false.",
        ],
    }
    all_converged = True
    all_finite = True
    all_input_unchanged = True
    all_invariants_ok = True
    all_baselines_reproduced = True
    arrays_for_npz = {}
    rows_for_csv = []
    for n in RESOLUTIONS:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=n, r_max=R_MAX, amplitude=AMPLITUDE, width=WIDTH,
            D_amplitude=D_AMPLITUDE, include_radiation=True,
        )
        digest_before = _state_digest(state)
        r = np.asarray(state.grid.centers, dtype=np.float64)
        baseline_initial_source = _source_for_B(state, np.ones_like(r))
        variants = {}
        for tag, mode, degree, fit_rmax in VARIANTS:
            witness, arrays = _run_variant(state, vacuum, tag, mode, degree, fit_rmax)
            variants[tag] = {
                "variant": {
                    "tag": tag, "mode": mode, "degree": degree,
                    "fit_r_max": fit_rmax,
                },
                **witness,
            }
            for key, values in arrays.items():
                arrays_for_npz[f"N{n}_{tag}_{key}"] = np.asarray(values, dtype=np.float64)
            all_converged &= bool(witness["convergence"]["converged"])
            all_input_unchanged &= bool(witness["input_state_unchanged"])
            all_baselines_reproduced &= witness["convergence"]["baseline_point_sample_six_step_reconstruction_max_abs"] <= 2.0e-12
            all_finite &= bool(witness["invariants"]["all_metric_profiles_finite"])
            all_invariants_ok &= all(
                witness["invariants"][k] <= 1.0e-12
                for k in (
                    "a_times_b_squared_minus_one_max_abs",
                    "b_minus_X_squared_max_abs",
                    "a_minus_X_minus_four_max_abs",
                    "X6_minus_B_max_abs",
                    "areal_radius_minus_r_max_abs",
                )
            )
            rows_for_csv.append({
                "N": n, "variant": tag, "cell0_gap": witness["cell0"],
                "cell1_gap": witness["cell1"], "first_two_gap_max_abs": witness["first_two_max_abs"],
                "centre_r_le_2_gap_max_abs": witness["centre_r_le_2_max_abs"],
                "annulus_r5_to10_gap_max_abs": witness["fixed_annulus_r5_to10_max_abs"],
                "hamiltonian_metric_only_cell0": witness["constraint_first_two"]["hamiltonian"][0],
                "hamiltonian_metric_only_cell1": witness["constraint_first_two"]["hamiltonian"][1],
                "hamiltonian_total_residual_cell0": witness["total_constraint_residual_first_two"]["hamiltonian_total_residual"][0],
                "hamiltonian_total_residual_cell1": witness["total_constraint_residual_first_two"]["hamiltonian_total_residual"][1],
                "hamiltonian_total_residual_max_abs": witness["total_constraint_residual_max_abs"]["hamiltonian_total_residual"],
                "momentum_total_residual_cell0": witness["total_constraint_residual_first_two"]["momentum_total_residual"][0],
                "momentum_total_residual_cell1": witness["total_constraint_residual_first_two"]["momentum_total_residual"][1],
                "momentum_total_residual_max_abs": witness["total_constraint_residual_max_abs"]["momentum_total_residual"],
                "iterations": witness["convergence"]["iterations_used"],
                "delta_integral": witness["half_cell"]["applied_delta_integral"],
                "s2_fixed": witness["half_cell"]["selected_fit"].get("s2_fixed"),
            })
        digest_after = _state_digest(state)
        if digest_after != digest_before:
            raise RuntimeError(f"N={n}: original initialized state mutated during copied diagnostic")
        report["resolutions"][str(n)] = {
            "resolution": n,
            "dr": float(state.grid.dr),
            "first_cell_radius": float(r[0]),
            "baseline_initial_source_first_two": [float(v) for v in baseline_initial_source[:2]],
            "variants": variants,
            "original_state_sha256_before": digest_before,
            "original_state_sha256_after": digest_after,
            "original_state_unchanged": digest_after == digest_before,
        }
    report["checks"] = {
        "source_pins_verified": True,
        "production_vendor_dependency_diff_empty": True,
        "all_fixed_points_converged": bool(all_converged),
        "all_profiles_finite": bool(all_finite),
        "all_input_states_unchanged": bool(all_input_unchanged),
        "all_metric_invariants_within_1e-12": bool(all_invariants_ok),
        "all_production_point_sample_reconstructions_within_2e-12": bool(all_baselines_reproduced),
    }
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    npz_path = OUT / "profiles.npz"
    np.savez_compressed(npz_path, **arrays_for_npz)
    csv_path = OUT / "summary.csv"
    with csv_path.open("w", encoding="utf-8") as stream:
        stream.write("N,variant,cell0_gap,cell1_gap,first_two_gap_max_abs,centre_r_le_2_gap_max_abs,annulus_r5_to10_gap_max_abs,iterations,delta_integral,s2_fixed\n")
        for row in rows_for_csv:
            values = [row[k] for k in ("N","variant","cell0_gap","cell1_gap","first_two_gap_max_abs","centre_r_le_2_gap_max_abs","annulus_r5_to10_gap_max_abs","iterations","delta_integral","s2_fixed")]
            stream.write(",".join(str(v) for v in values) + "\n")
    manifest_files = (report_path, npz_path, csv_path)
    (OUT / "SHA256SUMS.txt").write_text(
        "\n".join(f"{_sha256(p)}  {p.name}" for p in manifest_files) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "schema": report["schema"],
        "status": report["status"],
        "scientific_admission": report["scientific_admission"],
        "source_and_vendor_pins": {
            "diagnostic_head": provenance["diagnostic_head"],
            "source_commit": provenance["source_commit"],
            "vendor_commit": provenance["vendor_commit"],
        },
        "checks": report["checks"],
        "output_dir": str(OUT),
    }, indent=2, allow_nan=False), flush=True)
    if not all_converged or not all_finite or not all_input_unchanged or not all_invariants_ok or not all_baselines_reproduced:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
