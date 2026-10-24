#!/usr/bin/env python3
"""Source-pinned production evolution audit of regular-centre A_a/K behaviour.

Research-only observer. Runs the unchanged production kernel at N=60/120/240
to t=0.6, captures the actual CMC/PIRK geometry stages, decomposes the signed
A_a L2 source, and measures even-series origin extrapolants on fixed physical
support. It neither modifies production equations nor rejects a state because
of the measured regularity residual.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

SOURCE_COMMIT = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_BLOBS = {
    "engine/production_kernel.py": "1720615628a98a8d7cfc9fcadcc984fe28c88be1",
    "engine/v55_pirk_adapter.py": "be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3",
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/matter_system.py": "0969f4d996e4f8e0fc5b70a7dff147f024907302",
    "engine/matter_rhs.py": "0153c232ee15b5a9d3dbac37ea858e4d2d730011",
    "engine/scalar_system.py": "e2c2dae9a016b4d6454cf74a9eae8b0167d5733d",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
    "engine/native_geometry/__init__.py": "796d338d7cf2ed87bcb6ae63808dcf08110aeff8",
    "engine/native_geometry/evolution.py": "7b47d700680b12ab3c82cb6a556d40580b103ce9",
    "engine/native_geometry/operators.py": "2748be304b58328e2620f0a7fde49b01c0794ae2",
    "engine/native_geometry/geometry.py": "1bf12c1139572ff4597eaec0570aa43f2563fd91",
    "engine/native_geometry/constraints.py": "124d1f237e54d082b4d76213f9ebe144e02b02ba",
    "engine/handoff.py": "1e0b8a8a126c67dafaa5d370babdd519cbedf974",
    "engine/production_contract.py": "0a3b6f0a4223f93a5e33ab30706dc803f69643ca",
    "engine/invariant_diagnostics.py": "4bbcc6174c978b07c4f2175d1dcfeccb2d0c12b1",
    "engine/cosmology_observables.py": "4d6f39308685bae5e5bad2e17baa9c4043b67b4e",
    "engine/atomic_io.py": "484020bf0e2b8b1a37fe85a0f54b27d4cec5d2dd",
    "engine/radial_maxwell.py": "81269f1d1334c46c4c2f44ba816d792722c0a38e",
    "requirements.txt": "1e2a5bcbfd897f6465081a8794e14966f6e5db09",
    ".gitmodules": "502d49645859dc22c848efee0f654ed2170e18b9",
}

RESOLUTION = int(os.environ.get("REGULARITY_RESOLUTION", "120"))
R_MAX = float(os.environ.get("REGULARITY_R_MAX", "80.0"))
TARGET_TIME = float(os.environ.get("REGULARITY_TARGET_TIME", "0.6"))
CFL = float(os.environ.get("REGULARITY_CFL", "0.03"))
FIT_R_MAX = float(os.environ.get("REGULARITY_FIT_R_MAX", "8.0"))
FIT_DEGREE = int(os.environ.get("REGULARITY_FIT_DEGREE", "3"))
OUT = Path(os.environ.get(
    "REGULARITY_OUT",
    f"results/evolved-origin-aa-regularity-N{RESOLUTION}-t{TARGET_TIME:g}-2026-10-24",
))
GEOMETRY_STAGE_NAMES = (
    "start_lapse_refreshed",
    "explicit_predictor",
    "primary_predictor",
)


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def verify_source_pins() -> dict:
    base_vendor = _git("rev-parse", f"{SOURCE_COMMIT}:vendor/bb-palatini-unified-r0")
    head_vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if base_vendor != VENDOR_COMMIT or head_vendor != VENDOR_COMMIT:
        raise RuntimeError(
            f"vendor gitlink mismatch: base={base_vendor}, head={head_vendor}, expected={VENDOR_COMMIT}"
        )
    pin_rows = {}
    for path, expected in SOURCE_BLOBS.items():
        base_blob = _git("rev-parse", f"{SOURCE_COMMIT}:{path}")
        head_blob = _git("rev-parse", f"HEAD:{path}")
        matched = base_blob == expected and head_blob == expected
        pin_rows[path] = {
            "expected_blob": expected,
            "base_blob": base_blob,
            "head_blob": head_blob,
            "matches": matched,
        }
        if not matched:
            raise RuntimeError(f"source pin mismatch for {path}: {pin_rows[path]}")
    changed = _git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--",
        "engine", "vendor", ".gitmodules", "requirements.txt",
    )
    if changed:
        raise RuntimeError(
            "Production/source/vendor/dependency paths changed from the baseline: " + changed
        )
    return {
        "baseline_commit": SOURCE_COMMIT,
        "audit_head": _git("rev-parse", "HEAD"),
        "vendor_commit": head_vendor,
        "production_vendor_dependency_diff_empty": True,
        "blobs": pin_rows,
    }


def _arr(values) -> np.ndarray:
    return np.asarray(values, dtype=float).copy()


def _max_abs(values) -> float:
    values = _arr(values)
    return float(np.max(np.abs(values))) if values.size else 0.0


def _state_digest(state) -> str:
    arrays = {}
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        arrays["geometry." + name] = getattr(state.geometry, name)
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        arrays["scalars." + name] = getattr(state.scalars, name)
    for sector in ("dark_matter", "baryons", "radiation", "charged_dust"):
        item = getattr(state.matter, sector, None)
        if item is None:
            continue
        for name in ("rest", "energy_t", "momentum_r"):
            arrays["matter." + sector + "." + name] = getattr(item, name)
    digest = hashlib.sha256()
    for name in sorted(arrays):
        values = np.ascontiguousarray(np.asarray(arrays[name]))
        digest.update(name.encode("utf-8"))
        digest.update(str(values.dtype).encode("ascii"))
        digest.update(repr(values.shape).encode("ascii"))
        digest.update(values.tobytes())
    if state.electromagnetic is not None:
        for name in ("enclosed_charge_faces", "electric_field_faces"):
            if hasattr(state.electromagnetic, name):
                values = np.ascontiguousarray(np.asarray(getattr(state.electromagnetic, name)))
                digest.update(("electromagnetic." + name).encode("utf-8"))
                digest.update(str(values.dtype).encode("ascii"))
                digest.update(repr(values.shape).encode("ascii"))
                digest.update(values.tobytes())
    return digest.hexdigest()


def _even_origin_fit(r: np.ndarray, values: np.ndarray) -> dict:
    """Fit c0+c2*r^2+... on fixed physical support, excluding cell zero."""
    r = _arr(r)
    y = _arr(values)
    if r.shape != y.shape or r.ndim != 1:
        raise ValueError("even-origin fit requires matching one-dimensional arrays")
    mask = (r > 0.0) & (r <= FIT_R_MAX + 1.0e-12)
    mask[0] = False
    count = int(np.count_nonzero(mask))
    degree = min(FIT_DEGREE, count - 1)
    if degree < 1:
        raise ValueError(f"too few cells for origin fit: Nfit={count}")
    x = (r[mask] / FIT_R_MAX) ** 2
    design = np.column_stack([x**order for order in range(degree + 1)])
    coeff_scaled, *_ = np.linalg.lstsq(design, y[mask], rcond=None)
    fitted = design @ coeff_scaled
    residual = y[mask] - fitted
    physical_coeff = [
        float(coeff_scaled[order] / (FIT_R_MAX ** (2 * order)))
        for order in range(degree + 1)
    ]
    return {
        "origin_intercept": physical_coeff[0],
        "coefficients_c0_c2_etc": physical_coeff,
        "degree": degree,
        "fit_cells": count,
        "fit_r_max": FIT_R_MAX,
        "fit_max_abs_residual": _max_abs(residual),
        "fit_rms_residual": float(np.sqrt(np.mean(residual * residual))),
        "profile_max_abs_on_fit_support": _max_abs(y[mask]),
        "cell0": float(y[0]),
        "cell0_minus_fit_intercept": float(y[0] - physical_coeff[0]),
        "cell0_cell1_cell2_cell3": [float(v) for v in y[:4]],
    }


def _local_metric_regularization_counterfactual(
    grid, geometry, vacuum, actual_geometry_terms: dict[str, np.ndarray]
) -> dict:
    """Initial-slice-only regular continuation test; never impose it on evolved fields.

    q=(1-a/b)/r^2 is fitted through cells 1,2,3 in r^2. At initialization only,
    the fit defines the same determinant-consistent ratio a/b=1-r^2*q and the
    initializer's corresponding X and metric-derived Lambda. The unchanged
    pinned discrete Ricci operator is reevaluated on a copied state. After
    evolution X and Lambda are independent evolved fields, so this relation
    would not be an admissible same-state physical counterfactual. This helper
    is diagnostic only; it never mutates the accepted state or production equations.
    """
    r = _arr(grid.centers)
    if len(r) < 5:
        raise ValueError("local origin counterfactual requires at least five cells")
    r0 = float(r[0])
    if r0 <= 0.0:
        raise ValueError("cell-centred grid must have positive first radius")
    q = (1.0 - _arr(geometry.a) / _arr(geometry.b)) / r**2
    fit_x = r[1:4] ** 2
    coeff = np.polynomial.polynomial.polyfit(fit_x, q[1:4], deg=2)
    x0 = r0**2
    q0_fit = float(np.polynomial.polynomial.polyval(x0, coeff))
    qx0_fit = float(coeff[1] + 2.0 * coeff[2] * x0)
    ratio = 1.0 - x0 * q0_fit
    if not np.isfinite(ratio) or ratio <= 0.0:
        raise FloatingPointError("regularized centre metric ratio is not positive")

    a_reg = float(ratio ** (2.0 / 3.0))
    b_reg = float(ratio ** (-1.0 / 3.0))
    X_reg = float(ratio ** (-1.0 / 6.0))
    s_x = -(q0_fit + x0 * qx0_fit)
    a_x = (2.0 / 3.0) * ratio ** (-1.0 / 3.0) * s_x
    b_x = (-1.0 / 3.0) * ratio ** (-4.0 / 3.0) * s_x
    lambda_reg = float(
        r0 * a_x / a_reg**2
        - 2.0 * r0 * b_x / (a_reg * b_reg)
        + 2.0 / r0 * (1.0 / b_reg - 1.0 / a_reg)
    )

    alpha0 = float(geometry.alpha[0])
    actual_rtf = alpha0 * (
        float(actual_geometry_terms["Rr_r"][0])
        - float(actual_geometry_terms["R"][0]) / 3.0
    )
    actual_lapse_hessian_tf = -(
        float(actual_geometry_terms["DrDr_alpha"][0])
        - float(actual_geometry_terms["lap_alpha"][0]) / 3.0
    )
    actual_l2 = actual_lapse_hessian_tf + actual_rtf

    counterfactual = geometry.copy()
    counterfactual.a[0] = a_reg
    counterfactual.b[0] = b_reg
    counterfactual.X[0] = X_reg
    counterfactual.Lambda[0] = lambda_reg
    counter_terms = vacuum.geometry_terms(grid, counterfactual)
    counter_rtf = alpha0 * (
        float(counter_terms["Rr_r"][0]) - float(counter_terms["R"][0]) / 3.0
    )
    counter_l2 = -(
        float(counter_terms["DrDr_alpha"][0])
        - float(counter_terms["lap_alpha"][0]) / 3.0
    ) + counter_rtf

    return {
        "role": "offline_same_discrete_operator_counterfactual_not_a_production_fix",
        "metric_quotient": "q=(1-a/b)/r^2",
        "fit_support_cells": [1, 2, 3],
        "fit_variable": "q as a quadratic in r^2",
        "q_cell0_actual": float(q[0]),
        "q_cell0_fit_from_cells_1_to_3": q0_fit,
        "q_cell0_actual_minus_fit": float(q[0] - q0_fit),
        "a_cell0_actual": float(geometry.a[0]),
        "a_cell0_regularized": a_reg,
        "b_cell0_actual": float(geometry.b[0]),
        "b_cell0_regularized": b_reg,
        "X_cell0_actual": float(geometry.X[0]),
        "X_cell0_regularized": X_reg,
        "Lambda_cell0_actual": float(geometry.Lambda[0]),
        "Lambda_cell0_from_fitted_metric": lambda_reg,
        "tracefree_Ricci_cell0_actual": actual_rtf,
        "tracefree_Ricci_cell0_counterfactual": counter_rtf,
        "tracefree_Ricci_counterfactual_minus_actual": float(counter_rtf - actual_rtf),
        "Aa_L2_cell0_actual": actual_l2,
        "Aa_L2_cell0_counterfactual": float(counter_l2),
        "counterfactual_metric_identity_ab2": float(a_reg * b_reg**2),
        "accepted_state_mutated": False,
    }


def _capture_arrays(
    *,
    label: str,
    step_index: int,
    t_start: float,
    grid,
    geometry,
    vacuum,
    terms: dict | None = None,
    l2_only: dict | None = None,
) -> tuple[dict, dict[str, np.ndarray]]:
    r = _arr(grid.centers)
    g = geometry
    field_names = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda")
    arrays: dict[str, np.ndarray] = {
        "r": r.copy(),
        **{name: _arr(getattr(g, name)) for name in field_names},
    }
    geom_terms = vacuum.geometry_terms(grid, g)
    alpha = _arr(g.alpha)
    aa_lapse = -(
        _arr(geom_terms["DrDr_alpha"]) - _arr(geom_terms["lap_alpha"]) / 3.0
    )
    aa_ricci = alpha * (
        _arr(geom_terms["Rr_r"]) - _arr(geom_terms["R"]) / 3.0
    )
    arrays["Aa_lapse_hessian_tf"] = aa_lapse
    arrays["Aa_tracefree_ricci"] = aa_ricci
    l2_aa = None
    l2_k = None
    l3_aa = None
    l3_k = None
    if terms is not None:
        l2_aa = _arr(terms["primary_l2"]["Aa"])
        l2_k = _arr(terms["primary_l2"]["K"])
        l3_aa = _arr(terms["primary_l3"]["Aa"])
        l3_k = _arr(terms["primary_l3"]["K"])
        arrays["Aa_l3_total"] = l3_aa
        arrays["K_l3_total"] = l3_k
        arrays["a_rhs"] = _arr(terms["explicit"]["a"])
        arrays["b_rhs"] = _arr(terms["explicit"]["b"])
        arrays["X_rhs"] = _arr(terms["explicit"]["X"])
        vacuum_l3 = vacuum.primary_l3_rhs(grid, g)
        arrays["Aa_l3_vacuum"] = _arr(vacuum_l3["Aa"])
        arrays["Aa_l3_matter_increment"] = l3_aa - arrays["Aa_l3_vacuum"]
    elif l2_only is not None:
        l2_aa = _arr(l2_only["Aa"])
        l2_k = _arr(l2_only["K"])
    if l2_aa is not None:
        arrays["Aa_l2_total"] = l2_aa
    if l2_k is not None:
        arrays["K_l2_total"] = l2_k
    if l2_aa is not None:
        arrays["Aa_l2_closure"] = l2_aa - aa_lapse - aa_ricci
    if l2_aa is not None and l3_aa is not None:
        arrays["Aa_rhs_l2_plus_l3"] = l2_aa + l3_aa
        arrays["K_rhs_l2_plus_l3"] = l2_k + l3_k

    def fit_if_present(name: str) -> dict | None:
        return _even_origin_fit(r, arrays[name]) if name in arrays else None

    even_fields = {
        name: _even_origin_fit(r, arrays[name])
        for name in ("a", "b", "X", "alpha", "Aa", "K")
    }
    source_fits = {
        name: fit_if_present(name)
        for name in (
            "Aa_lapse_hessian_tf", "Aa_tracefree_ricci",
            "Aa_l2_total", "Aa_l3_total", "Aa_l3_vacuum",
            "Aa_l3_matter_increment", "Aa_rhs_l2_plus_l3", "K_rhs_l2_plus_l3",
        )
    }
    source_fits = {key: value for key, value in source_fits.items() if value is not None}

    identities = {}
    if terms is not None:
        a, b, X = _arr(g.a), _arr(g.b), _arr(g.X)
        Aa, K = _arr(g.Aa), _arr(g.K)
        identities = {
            "max_abs_a_rhs_plus_2_alpha_a_Aa": _max_abs(
                arrays["a_rhs"] + 2.0 * alpha * a * Aa
            ),
            "max_abs_b_rhs_minus_alpha_b_Aa": _max_abs(
                arrays["b_rhs"] - alpha * b * Aa
            ),
            "max_abs_X_rhs_minus_alpha_X_K_over3": _max_abs(
                arrays["X_rhs"] - alpha * X * K / 3.0
            ),
            "beta_max_abs": _max_abs(g.beta),
        }
    closure = _max_abs(arrays["Aa_l2_closure"]) if "Aa_l2_closure" in arrays else None
    record = {
        "label": label,
        "local_metric_regularization_counterfactual": (
            _local_metric_regularization_counterfactual(
                grid, g, vacuum, geom_terms
            )
            if label == "accepted_000" and r.size >= 5 else None
        ),
        "step_index": int(step_index),
        "accepted_time_at_step_start": float(t_start),
        "regularity_fits_even_fields": even_fields,
        "regularity_fits_Aa_source_channels": source_fits,
        "metric_rhs_identities": identities,
        "Aa_L2_signed_decomposition": {
            "formula": "Aa_L2 = -(DrDr_alpha - lap_alpha/3) + alpha*(Rr_r - R/3)",
            "cell0_lapse_hessian": float(aa_lapse[0]),
            "cell0_tracefree_ricci": float(aa_ricci[0]),
            "cell0_primary_l2": float(l2_aa[0]) if l2_aa is not None else None,
            "cell0_signed_closure": (
                float(arrays["Aa_l2_closure"][0]) if "Aa_l2_closure" in arrays else None
            ),
            "max_abs_closure_full_grid": closure,
        },
        "available_profiles": sorted(k for k in arrays if k != "r"),
    }
    return record, arrays


BASELINE_REFERENCE_RUN_ID = 40127716287
BASELINE_ARTIFACT_PREFIX = "evolved-origin-aa-regularity-N"
ORIGIN_HALF_CELL_WEIGHTS = np.asarray(
    [401.0 / 720.0, -31.0 / 480.0, 11.0 / 1440.0], dtype=float
)


def _integrate_origin_half_cell(src: np.ndarray, dr: float) -> float:
    """Exact origin half-cell integral for an even source through r^4.

    This is the preserved PR #84/#106 candidate rule. All inter-centre
    intervals remain the production trapezoid; no other source or coefficient
    is changed.
    """
    values = _arr(src)
    if values.ndim != 1 or values.size < 3 or not np.all(np.isfinite(values[:3])):
        raise ValueError("origin half-cell rule requires three finite source samples")
    if not np.isfinite(dr) or dr <= 0.0:
        raise ValueError("origin half-cell rule requires a positive finite dr")
    return float(dr * np.dot(ORIGIN_HALF_CELL_WEIGHTS, values[:3]))


def _rebuild_origin_half_cell_candidate(kernel, state) -> dict:
    """Apply the archived origin-half-cell initializer candidate to this copy.

    Reconstructs the six-pass initializer's Hamiltonian-source recurrence from
    the accepted baseline initialized data. The baseline point-sampled
    recurrence must independently reproduce X^6 before the candidate rule is
    applied. Only the initial metric recurrence changes; the physical scalar
    profiles and archive matter densities are reinitialized consistently on
    the candidate metric. Production code and the supplied state object are
    not changed outside this diagnostic arm.
    """
    from engine.scalar_system import cosmos_potential
    from engine.v55_matter import (
        RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT, initialize_from_archive,
    )

    grid = state.grid
    r = _arr(grid.centers)
    dr = float(grid.dr)
    g = state.geometry
    fields = state.scalars
    S, PS = _arr(fields.S), _arr(fields.PS)
    D, PD = _arr(fields.D), _arr(fields.PD)
    phi, Pi = _arr(fields.phi), _arr(fields.Pi)
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    Sp, Dp = d1(S, 1), d1(D, 1)
    current = -(PD * Dp)

    phi0, pi0 = 0.179055, 0.00341913
    rho_r = float(RHO_R_PRESENT)
    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.asarray([phi0]))[0])
    rho0 += float(RHO_DM_PRESENT) + float(RHO_B_PRESENT) + rho_r
    H0 = float(np.sqrt(max(rho0 / 3.0, 0.0)))
    Kt, Ktp = np.full_like(r, -H0), d1(np.full_like(r, -H0), 1)
    rho_matter = (
        float(RHO_DM_PRESENT) + float(RHO_B_PRESENT) + rho_r
    ) / (8.0 * np.pi)
    rho_phi = (0.5 * Pi**2 + cosmos_potential(phi)) / (8.0 * np.pi)

    def source_for_B(B: np.ndarray) -> np.ndarray:
        rhoS = 0.5 * (PS**2 + B * Sp**2) + 0.5 * S**2
        rhoD = 0.5 * (PD**2 + B * Dp**2) - 0.5 * D**2
        rho = rhoS + rhoD + rho_phi + rho_matter
        Kr = Kt + r * Ktp + 4.0 * np.pi * r * current
        source = (
            1.0 - 8.0 * np.pi * r**2 * rho
            + 2.0 * r**2 * Kr * Kt
            + r**2 * Kt**2
        )
        if not np.all(np.isfinite(source)):
            raise FloatingPointError("non-finite source in copied origin initializer")
        return source

    def integrate(src: np.ndarray, origin_rule: str) -> np.ndarray:
        if origin_rule == "point_sampled":
            origin_integral = 0.5 * dr * float(src[0])
        elif origin_rule == "origin_half_cell_candidate":
            origin_integral = _integrate_origin_half_cell(src, dr)
        else:
            raise ValueError(f"unknown origin integral rule: {origin_rule}")
        integral = np.cumsum(np.r_[
            origin_integral, 0.5 * dr * (src[1:] + src[:-1])
        ])
        return integral / r

    B_baseline = np.ones_like(r)
    B_candidate = np.ones_like(r)
    baseline_updates = []
    candidate_updates = []
    for _ in range(6):
        src_baseline = source_for_B(B_baseline)
        Bnext_baseline = integrate(src_baseline, "point_sampled")
        baseline_updates.append(_max_abs(Bnext_baseline - B_baseline))
        B_baseline = Bnext_baseline

        src_candidate = source_for_B(B_candidate)
        Bnext_candidate = integrate(src_candidate, "origin_half_cell_candidate")
        candidate_updates.append(_max_abs(Bnext_candidate - B_candidate))
        B_candidate = Bnext_candidate

    baseline_B_from_state = _arr(g.X)**6
    baseline_reconstruction_error = _max_abs(B_baseline - baseline_B_from_state)
    baseline_scale = max(1.0, _max_abs(baseline_B_from_state))
    if baseline_reconstruction_error > 2.0e-12 * baseline_scale:
        raise RuntimeError(
            "baseline point-sampled recurrence failed X^6 reconstruction: "
            f"{baseline_reconstruction_error:.6e}"
        )
    if not np.all(np.isfinite(B_candidate)) or np.any(B_candidate <= 0.0):
        raise FloatingPointError("origin-half-cell candidate produced invalid B")

    # Reconstruct the same initial-data map used by the archived candidate.
    A = 1.0 / np.sqrt(B_candidate)
    a_candidate = A ** (4.0 / 3.0)
    b_candidate = A ** (-2.0 / 3.0)
    X_candidate = A ** (-1.0 / 3.0)
    Lambda_candidate = (
        d1(a_candidate, 1) / (2.0 * a_candidate**2)
        - d1(b_candidate, 1) / (a_candidate * b_candidate)
        + 2.0 / r * (1.0 / b_candidate - 1.0 / a_candidate)
    )
    Kr = Kt + r * Ktp + 4.0 * np.pi * r * current
    K_candidate = Kr + 2.0 * Kt
    Aa_candidate = Kr - K_candidate / 3.0

    original_digest = _state_digest(state)
    candidate_state = state
    gc = candidate_state.geometry
    gc.a[:] = a_candidate
    gc.b[:] = b_candidate
    gc.X[:] = X_candidate
    gc.Lambda[:] = Lambda_candidate
    gc.K[:] = K_candidate
    gc.Aa[:] = Aa_candidate
    gc.beta.fill(0.0)
    gc.B.fill(0.0)
    candidate_state.matter = initialize_from_archive(
        grid, gc, include_radiation=True
    )
    # The baseline initialize() has already run the CMC solve once. Re-solve
    # on the candidate metric and its reinitialized conservative matter state.
    gc.alpha[:] = kernel._solve_lapse_for_state(
        grid, gc, candidate_state.scalars, candidate_state.matter,
        candidate_state.electromagnetic,
    )[0]
    gc.beta.fill(0.0)
    gc.B.fill(0.0)

    after_digest = _state_digest(candidate_state)
    return {
        "method": "PR #84/#106 origin-half-cell-only initializer candidate",
        "source_rule": "dr * dot([401/720, -31/480, 11/1440], source_cells_0_to_2)",
        "source_rule_polynomial_exactness": "even source polynomial in r through r^4",
        "production_intercentre_rule": "unchanged trapezoidal quadrature",
        "initializer_iterations": 6,
        "baseline_point_sampled_reconstruction_max_abs_X6": baseline_reconstruction_error,
        "baseline_reconstruction_matches_X6": True,
        "baseline_point_sampled_max_abs_iteration_updates": baseline_updates,
        "candidate_max_abs_iteration_updates": candidate_updates,
        "candidate_B_cell0_minus_baseline_B_cell0": float(B_candidate[0] - B_baseline[0]),
        "candidate_B_first_three": [float(v) for v in B_candidate[:3]],
        "baseline_B_first_three": [float(v) for v in B_baseline[:3]],
        "q=(1-a/b)/r^2_cell0_baseline": float(
            (1.0 - state.geometry.a[0] / state.geometry.b[0]) / r[0]**2
        ) if False else None,
        "candidate_initial_metric_algebraic_max_abs_ab2_minus_one": _max_abs(
            a_candidate * b_candidate**2 - 1.0
        ),
        "baseline_state_digest_before_candidate_build": original_digest,
        "candidate_initial_state_digest_after_cmc_refresh": after_digest,
        "candidate_initial_state_reconstruction_only": True,
        "production_or_vendor_files_modified": False,
        "scientific_admission": False,
    }


def _load_baseline_reference() -> tuple[dict, dict[str, np.ndarray], dict]:
    """Verify and load the already-accepted baseline artifact from run #40127716287."""
    baseline_dir_value = os.environ.get("REGULARITY_BASELINE_REFERENCE_DIR")
    if not baseline_dir_value:
        raise RuntimeError("REGULARITY_BASELINE_REFERENCE_DIR is required for candidate A/B")
    directory = Path(baseline_dir_value)
    required = ("report.json", "stage_profiles.npz", "final_state.npz", "SHA256SUMS.json")
    missing = [name for name in required if not (directory / name).is_file()]
    if missing:
        raise FileNotFoundError(f"baseline reference artifact missing {missing} in {directory}")

    sums = json.loads((directory / "SHA256SUMS.json").read_text(encoding="utf-8"))
    for name in ("report.json", "stage_profiles.npz", "final_state.npz"):
        actual = hashlib.sha256((directory / name).read_bytes()).hexdigest()
        expected = sums.get(name)
        if not expected or actual != expected:
            raise RuntimeError(
                f"baseline reference checksum mismatch for {name}: expected={expected}, actual={actual}"
            )
    baseline = json.loads((directory / "report.json").read_text(encoding="utf-8"))
    if baseline.get("schema") != "source_pinned_evolved_origin_Aa_regularity_v1":
        raise RuntimeError(f"unexpected baseline report schema: {baseline.get('schema')}")
    if baseline.get("source_pins", {}).get("baseline_commit") != SOURCE_COMMIT:
        raise RuntimeError("baseline reference uses a different production source commit")
    if baseline.get("source_pins", {}).get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError("baseline reference uses a different vendor pin")
    pins = baseline.get("source_pins", {}).get("blobs", {})
    for path, expected in SOURCE_BLOBS.items():
        row = pins.get(path)
        if not row or row.get("expected_blob") != expected:
            raise RuntimeError(f"baseline reference source pin missing/mismatched: {path}")
        if row.get("base_blob") != expected or row.get("head_blob") != expected or not row.get("matches"):
            raise RuntimeError(f"baseline reference does not verify production source blob: {path}")
    config = baseline.get("configuration", {})
    if (
        int(config.get("resolution", -1)) != RESOLUTION
        or not np.isclose(float(config.get("r_max", np.nan)), R_MAX, rtol=0.0, atol=1e-12)
        or not np.isclose(float(config.get("target_time", np.nan)), TARGET_TIME, rtol=0.0, atol=1e-12)
        or not np.isclose(float(config.get("cfl", np.nan)), CFL, rtol=0.0, atol=1e-12)
        or not np.isclose(float(config.get("fit", {}).get("fit_r_max", FIT_R_MAX), FIT_R_MAX), FIT_R_MAX)
        or int(config.get("fit", {}).get("degree", FIT_DEGREE)) != FIT_DEGREE
    ):
        raise RuntimeError(f"baseline reference configuration mismatch: {config}")
    with np.load(directory / "stage_profiles.npz", allow_pickle=False) as archive:
        profiles = {name: np.asarray(archive[name], dtype=float).copy() for name in archive.files}
    reference = {
        "workflow_run_id": BASELINE_REFERENCE_RUN_ID,
        "artifact_name": f"{BASELINE_ARTIFACT_PREFIX}{RESOLUTION}-t0.6-2026-10-24",
        "artifact_report_audit_head": baseline.get("source_pins", {}).get("audit_head"),
        "report_sha256": sums["report.json"],
        "stage_profiles_sha256": sums["stage_profiles.npz"],
        "final_state_sha256": sums["final_state.npz"],
        "reference_directory": str(directory),
        "checksums_verified": True,
        "production_source_pins_verified": True,
    }
    return baseline, profiles, reference


def _record_metrics(record: dict) -> dict:
    fields = record["regularity_fits_even_fields"]
    src = record["Aa_L2_signed_decomposition"]
    source_fits = record.get("regularity_fits_Aa_source_channels", {})
    return {
        "time": float(record["accepted_time_at_step_start"]),
        "Aa_cell0": float(fields["Aa"]["cell0"]),
        "Aa_even_fit_origin": float(fields["Aa"]["origin_intercept"]),
        "Aa_cell0_minus_even_fit": float(fields["Aa"]["cell0_minus_fit_intercept"]),
        "K_cell0": float(fields["K"]["cell0"]),
        "K_cell0_minus_even_fit": float(fields["K"]["cell0_minus_fit_intercept"]),
        "Aa_L2_cell0": float(src["cell0_primary_l2"]) if src["cell0_primary_l2"] is not None else None,
        "Aa_lapse_hessian_cell0": float(src["cell0_lapse_hessian"]),
        "Aa_tracefree_Ricci_cell0": float(src["cell0_tracefree_ricci"]),
        "Aa_L2_even_fit_origin": source_fits.get("Aa_l2_total", {}).get("origin_intercept"),
        "Aa_tracefree_Ricci_even_fit_origin": source_fits.get("Aa_tracefree_ricci", {}).get("origin_intercept"),
    }


def _compare_candidate_to_baseline(
    baseline_report: dict, baseline_profiles: dict[str, np.ndarray],
    candidate_report: dict, candidate_profiles: dict[str, np.ndarray],
) -> dict:
    baseline_records = [r for r in baseline_report["records"] if r["label"].startswith("accepted_")]
    candidate_records = [r for r in candidate_report["records"] if r["label"].startswith("accepted_")]
    if len(baseline_records) != len(candidate_records):
        raise RuntimeError(
            f"accepted-state history length mismatch: baseline={len(baseline_records)}, candidate={len(candidate_records)}"
        )
    r = _arr(candidate_profiles["r"])
    support = r <= 10.0
    if not np.any(support):
        raise RuntimeError("comparison physical support r<=10 is empty")
    history = []
    for index, (base_rec, cand_rec) in enumerate(zip(baseline_records, candidate_records)):
        if base_rec["label"] != cand_rec["label"]:
            raise RuntimeError(
                f"accepted-state label mismatch at {index}: {base_rec['label']} vs {cand_rec['label']}"
            )
        base_metrics = _record_metrics(base_rec)
        cand_metrics = _record_metrics(cand_rec)
        if abs(base_metrics["time"] - cand_metrics["time"]) > 1.0e-12:
            raise RuntimeError(
                f"accepted-state time mismatch at {index}: {base_metrics['time']} vs {cand_metrics['time']}"
            )
        history.append({
            "accepted_index": index,
            "time": cand_metrics["time"],
            "baseline": base_metrics,
            "origin_half_cell_candidate": cand_metrics,
            "candidate_minus_baseline": {
                key: (
                    float(cand_metrics[key] - base_metrics[key])
                    if cand_metrics.get(key) is not None and base_metrics.get(key) is not None
                    else None
                )
                for key in base_metrics if key != "time"
            },
        })

    def endpoint_profile_comparison(index: int) -> dict:
        base_rec, cand_rec = baseline_records[index], candidate_records[index]
        base_prefix, cand_prefix = base_rec["profile_prefix"], cand_rec["profile_prefix"]
        items = {}
        for field in (
            "a", "b", "X", "alpha", "Aa", "K",
            "Aa_lapse_hessian_tf", "Aa_tracefree_ricci", "Aa_l2_total", "Aa_l3_total",
        ):
            base_key, cand_key = base_prefix + "__" + field, cand_prefix + "__" + field
            if base_key not in baseline_profiles or cand_key not in candidate_profiles:
                continue
            bvalues = _arr(baseline_profiles[base_key])
            cvalues = _arr(candidate_profiles[cand_key])
            delta = cvalues - bvalues
            items[field] = {
                "cell0_baseline": float(bvalues[0]),
                "cell0_candidate": float(cvalues[0]),
                "cell0_candidate_minus_baseline": float(delta[0]),
                "max_abs_candidate_minus_baseline_r_le_10": _max_abs(delta[support]),
                "rms_candidate_minus_baseline_r_le_10": float(np.sqrt(np.mean(delta[support]**2))),
                "max_abs_candidate_minus_baseline_full_grid": _max_abs(delta),
            }
        return {
            "time": candidate_records[index]["accepted_time_at_step_start"],
            "baseline_profile_prefix": base_prefix,
            "candidate_profile_prefix": cand_prefix,
            "fields": items,
        }

    initial = history[0]
    final = history[-1]
    initial_base_l2 = initial["baseline"]["Aa_L2_cell0"]
    initial_cand_l2 = initial["origin_half_cell_candidate"]["Aa_L2_cell0"]
    final_base_l2 = final["baseline"]["Aa_L2_cell0"]
    final_cand_l2 = final["origin_half_cell_candidate"]["Aa_L2_cell0"]
    initial_source_abs_reduction = (
        1.0 - abs(initial_cand_l2) / abs(initial_base_l2)
        if initial_base_l2 not in (None, 0.0) and initial_cand_l2 is not None else None
    )
    final_source_abs_change_ratio = (
        abs(final_cand_l2) / abs(final_base_l2)
        if final_base_l2 not in (None, 0.0) and final_cand_l2 is not None else None
    )
    return {
        "schema": "origin_half_cell_candidate_vs_existing_baseline_v1",
        "baseline_reference_run_id": BASELINE_REFERENCE_RUN_ID,
        "initialization_method": "historical PR #84/#106 exact three-sample origin half-cell quadrature candidate; production trapezoids retained elsewhere",
        "comparison_contract": {
            "same_production_source_commit": SOURCE_COMMIT,
            "same_vendor_commit": VENDOR_COMMIT,
            "same_resolution_and_domain": True,
            "same_cfl_and_target_time": True,
            "same_accepted_state_times": True,
            "baseline_reused_from_existing_artifact_without_re_evolving": True,
            "candidate_evolved_by_unchanged_pinned_production_kernel": True,
            "physical_profile_comparison_support": "r<=10",
            "fits_are_diagnostic_only": True,
            "scientific_admission": False,
        },
        "baseline_state_digest_initial": baseline_report.get("initial_state_sha256"),
        "baseline_state_digest_final": baseline_report.get("final_state_sha256"),
        "candidate_state_digest_initial": candidate_report.get("initial_state_sha256"),
        "candidate_state_digest_final": candidate_report.get("final_state_sha256"),
        "baseline_reference_integrity": baseline_report.get("accepted_state_observer_mutation_check_passed"),
        "key_outcomes": {
            "baseline_initial": initial["baseline"],
            "candidate_initial": initial["origin_half_cell_candidate"],
            "initial_Aa_L2_abs_reduction_fraction": initial_source_abs_reduction,
            "baseline_final": final["baseline"],
            "candidate_final": final["origin_half_cell_candidate"],
            "final_Aa_L2_abs_candidate_to_baseline_ratio": final_source_abs_change_ratio,
            "final_Aa_cell0_candidate_minus_baseline": final["candidate_minus_baseline"]["Aa_cell0"],
            "final_Aa_cell0_minus_even_fit_candidate_minus_baseline": final["candidate_minus_baseline"]["Aa_cell0_minus_even_fit"],
            "final_K_cell0_minus_even_fit_candidate_minus_baseline": final["candidate_minus_baseline"]["K_cell0_minus_even_fit"],
            "initial_profile_delta": endpoint_profile_comparison(0),
            "final_profile_delta": endpoint_profile_comparison(len(history)-1),
        },
        "accepted_state_centre_source_history": history,
        "interpretation_limits": [
            "The baseline is the already completed exact-source run artifact, verified internally and against every production source/vendor blob; baseline was not re-evolved.",
            "The candidate differs only in origin half-cell quadrature inside a copied six-pass initializer reconstruction; inter-centre trapezoids, evolution equations, coefficients, gauge, boundary rules and production/vendor files are unchanged.",
            "The archived candidate's quadrature is exact only for even source polynomials through r^4; real-source truncation and matter/metric coupling are not eliminated by the quadrature formula.",
            "A lower initial A_a source does not establish a correct continuum centre, full constraint closure, or long-time scientific admission.",
            "The final fit intercepts use the same fixed support/degree as the existing baseline report and are shape indicators, not pass thresholds.",
        ],
        "scientific_admission": False,
    }


def main() -> int:
    if RESOLUTION < 16 or R_MAX <= 0.0 or TARGET_TIME <= 0.0 or CFL <= 0.0:
        raise ValueError("invalid resolution / domain / time / CFL")
    if FIT_R_MAX <= 0.0:
        raise ValueError("fit radius must be positive")
    provenance = verify_source_pins()

    from engine.production_kernel import V55ProductionKernel
    from engine.v55_pirk_adapter import vendor_modules

    _, vacuum, _ = vendor_modules()
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    dt = CFL * R_MAX / RESOLUTION
    steps = int(round(TARGET_TIME / dt))
    if steps < 1 or not np.isclose(steps * dt, TARGET_TIME, rtol=0.0, atol=1.0e-12):
        raise ValueError(
            f"target time {TARGET_TIME} is not an integer number of steps at dt={dt}"
        )

    out_profiles: dict[str, np.ndarray] = {"r": _arr(state.grid.centers)}
    records: list[dict] = []
    step_budgets: list[dict] = []
    active = {
        "enabled": False,
        "inside_geometry_stage": False,
        "step_index": None,
        "t_start": None,
        "geometry_calls": 0,
        "stage_records": [],
        "final_l2_record": None,
    }
    original_geometry_stage = kernel._geometry_stage_terms
    original_primary_l2 = vacuum.primary_l2_rhs

    def capture(label, geometry, terms=None, l2_only=None):
        rec, arrays = _capture_arrays(
            label=label,
            step_index=active["step_index"] if active["step_index"] is not None else -1,
            t_start=active["t_start"] if active["t_start"] is not None else float(state.t),
            grid=state.grid,
            geometry=geometry,
            vacuum=vacuum,
            terms=terms,
            l2_only=l2_only,
        )
        prefix = f"{len(records):04d}_{label}"
        for name, values in arrays.items():
            if name == "r":
                continue
            out_profiles[prefix + "__" + name] = _arr(values)
        rec["profile_prefix"] = prefix
        rec["capture_index"] = len(records)
        records.append(rec)
        return rec, arrays

    def wrapped_geometry_stage(grid, geometry, scalars, matter, electromagnetic):
        if not active["enabled"]:
            return original_geometry_stage(grid, geometry, scalars, matter, electromagnetic)
        active["inside_geometry_stage"] = True
        try:
            terms = original_geometry_stage(grid, geometry, scalars, matter, electromagnetic)
        finally:
            active["inside_geometry_stage"] = False
        call_index = active["geometry_calls"]
        if call_index >= len(GEOMETRY_STAGE_NAMES):
            raise RuntimeError(f"unexpected geometry stage call {call_index}")
        label = GEOMETRY_STAGE_NAMES[call_index]
        active["geometry_calls"] += 1
        rec, arrays = capture("step%03d_%s" % (active["step_index"], label), geometry, terms=terms)
        active["stage_records"].append((label, rec, arrays))
        return terms

    def wrapped_primary_l2(grid, geometry):
        result = original_primary_l2(grid, geometry)
        # CMC lapse projection also calls primary_l2_rhs outside the geometry
        # stage. Only the exact production step() call is the final implicit
        # A_a/K operand needed by the accepted-state update reconstruction.
        caller = sys._getframe(1)
        is_production_step_operand = (
            caller.f_code.co_name == "step"
            and caller.f_code.co_filename.endswith("/engine/production_kernel.py")
        )
        if (
            active["enabled"]
            and not active["inside_geometry_stage"]
            and is_production_step_operand
        ):
            if active["final_l2_record"] is not None:
                raise RuntimeError("duplicate production-step final primary_l2 operand")
            rec, arrays = capture(
                "step%03d_final_u_old_v" % active["step_index"],
                geometry,
                l2_only=result,
            )
            active["final_l2_record"] = (rec, arrays)
        return result

    kernel._geometry_stage_terms = wrapped_geometry_stage
    vacuum.primary_l2_rhs = wrapped_primary_l2
    try:
        # Accepted-state observer: same grid and operators, no evolution, then
        # verify that observation leaves every saved state array byte-identical.
        for step_index in range(steps + 1):
            active.update({
                "enabled": False,
                "step_index": step_index,
                "t_start": float(state.t),
                "inside_geometry_stage": False,
            })
            before = _state_digest(state)
            terms = original_geometry_stage(
                state.grid, state.geometry, state.scalars, state.matter, state.electromagnetic
            )
            rec, arrays = capture("accepted_%03d" % step_index, state.geometry, terms=terms)
            after = _state_digest(state)
            if before != after:
                raise RuntimeError(f"accepted-state observer mutated state at index {step_index}")
            rec["state_sha256_before_observation"] = before
            rec["state_sha256_after_observation"] = after
            rec["observer_left_state_unchanged"] = True
            if step_index == 0:
                initial_r = _arr(state.grid.centers)
                initial_q0 = float(
                    (1.0 - state.geometry.a[0] / state.geometry.b[0])
                    / initial_r[0]**2
                )
                # Initial rho_S=0.5*S^2 and src contains -8*pi*r^2*rho_S.
                # Since B[0]=src[0], this source alone predicts q0 ~= -4*pi*S[0]^2.
                pulse_q0 = float(-4.0 * np.pi * np.asarray(state.scalars.S)[0]**2)
                rec["initial_data_origin_quadrature_audit"] = {
                    "initializer_rule": "B[0]=src[0] because integ[0]=0.5*dr*src[0] and r[0]=0.5*dr",
                    "observed_q0": initial_q0,
                    "local_source_prediction_minus_4pi_S0_squared": pulse_q0,
                    "observed_minus_pulse_prediction": float(initial_q0 - pulse_q0),
                    "source_term": "src contains -8*pi*r^2*rho_S with rho_S=0.5*S^2",
                    "interpretation": "diagnostic correlation; compare with the correct even-source integral average before any production change",
                }

            if step_index == steps:
                break

            old_state = state
            active.update({
                "enabled": True,
                "inside_geometry_stage": False,
                "step_index": step_index,
                "t_start": float(old_state.t),
                "geometry_calls": 0,
                "stage_records": [],
                "final_l2_record": None,
            })
            state = kernel.step(old_state, dt)
            active["enabled"] = False
            if active["geometry_calls"] != 3 or active["final_l2_record"] is None:
                raise RuntimeError(
                    "production step capture incomplete: "
                    f"geometry calls={active['geometry_calls']}, final L2={active['final_l2_record'] is not None}"
                )
            stages = {name: arrays for name, _, arrays in active["stage_records"]}
            stage_metadata = {name: rec for name, rec, _ in active["stage_records"]}
            start = stages["start_lapse_refreshed"]
            primary = stages["primary_predictor"]
            final_l2 = active["final_l2_record"][1]
            r = _arr(state.grid.centers)
            interior = r <= 10.0
            actual_Aa_delta = _arr(state.geometry.Aa) - _arr(old_state.geometry.Aa)
            actual_K_delta = _arr(state.geometry.K) - _arr(old_state.geometry.K)
            predicted_Aa_delta = 0.5 * dt * (
                _arr(start["Aa_l2_total"]) + _arr(final_l2["Aa_l2_total"])
                + _arr(start["Aa_l3_total"]) + _arr(primary["Aa_l3_total"])
            )
            predicted_K_delta = 0.5 * dt * (
                _arr(start["K_l2_total"]) + _arr(final_l2["K_l2_total"])
                + _arr(start["K_l3_total"]) + _arr(primary["K_l3_total"])
            )
            closure_Aa = actual_Aa_delta - predicted_Aa_delta
            closure_K = actual_K_delta - predicted_K_delta
            budget = {
                "step_index": step_index,
                "time_start": float(old_state.t),
                "time_end": float(state.t),
                "dt": dt,
                "Aa_cell0_start": float(old_state.geometry.Aa[0]),
                "Aa_cell0_end": float(state.geometry.Aa[0]),
                "Aa_cell0_actual_delta": float(actual_Aa_delta[0]),
                "Aa_cell0_predicted_delta": float(predicted_Aa_delta[0]),
                "Aa_cell0_signed_update_closure": float(closure_Aa[0]),
                "Aa_update_closure_max_abs_r_le_10": _max_abs(closure_Aa[interior]),
                "Aa_update_closure_max_abs_full_grid": _max_abs(closure_Aa),
                "K_cell0_start": float(old_state.geometry.K[0]),
                "K_cell0_end": float(state.geometry.K[0]),
                "K_cell0_actual_delta": float(actual_K_delta[0]),
                "K_cell0_predicted_delta": float(predicted_K_delta[0]),
                "K_cell0_signed_update_closure": float(closure_K[0]),
                "K_update_closure_max_abs_r_le_10": _max_abs(closure_K[interior]),
                "K_update_closure_max_abs_full_grid": _max_abs(closure_K),
                "Aa_regular_origin_fit_start": stage_metadata["start_lapse_refreshed"]["regularity_fits_even_fields"]["Aa"],
                "Aa_regular_origin_fit_primary_predictor": stage_metadata["primary_predictor"]["regularity_fits_even_fields"]["Aa"],
            }
            # Preserve the signed pointwise residuals in the compact report;
            # full profiles below carry every cell's closure and source channels.
            for name, values in (
                ("Aa_actual_delta", actual_Aa_delta),
                ("Aa_predicted_delta", predicted_Aa_delta),
                ("Aa_update_closure", closure_Aa),
                ("K_actual_delta", actual_K_delta),
                ("K_predicted_delta", predicted_K_delta),
                ("K_update_closure", closure_K),
            ):
                prefix = f"step{step_index:03d}_budget__{name}"
                out_profiles[prefix] = _arr(values)
            step_budgets.append(budget)
            if budget["Aa_update_closure_max_abs_r_le_10"] > 2.0e-12:
                raise RuntimeError(
                    "production A_a update did not close at r<=10 after exact captured-stage reconstruction: "
                    f"step {step_index}, residual={budget['Aa_update_closure_max_abs_r_le_10']:.6e}"
                )
            if budget["K_update_closure_max_abs_r_le_10"] > 2.0e-12:
                raise RuntimeError(
                    "production K update did not close at r<=10 after exact captured-stage reconstruction: "
                    f"step {step_index}, residual={budget['K_update_closure_max_abs_r_le_10']:.6e}"
                )
    finally:
        kernel._geometry_stage_terms = original_geometry_stage
        vacuum.primary_l2_rhs = original_primary_l2

    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "schema": "source_pinned_evolved_origin_Aa_regularity_v1",
        "status": "completed_bounded_production_origin_regularity_audit",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "source_pins": provenance,
        "configuration": {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "dr": float(R_MAX / RESOLUTION),
            "first_cell_center": float(R_MAX / RESOLUTION / 2.0),
            "cfl": CFL,
            "dt": dt,
            "target_time": TARGET_TIME,
            "accepted_steps": steps,
            "fit": {
                "parity_assumption": "even Taylor series in r: c0+c2*r^2+c4*r^4+...",
                "physical_support": f"0 < r <= {FIT_R_MAX}",
                "cell_zero_excluded_from_fit": True,
                "degree": FIT_DEGREE,
                "purpose": "extrapolation diagnostic only; no regularity acceptance threshold",
            },
        },
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "initial_state_sha256": None,
        "final_state_sha256": _state_digest(state),
        "accepted_state_observer_mutation_check_passed": all(
            row.get("observer_left_state_unchanged", False)
            for row in records if row["label"].startswith("accepted_")
        ),
        "records": records,
        "accepted_step_Aa_K_update_budgets": step_budgets,
        "interpretation_limits": [
            "The fit window and polynomial degree are explicit diagnostics, not a physical regularity bound.",
            "A nonzero extrapolated Aa intercept flags incompatibility with the expected even, trace-free centre series for this discrete profile; it does not identify the continuum equation responsible.",
            "The A_a L2 split is the pinned production source decomposition; its algebraic closure alone is not a continuum regularity proof.",
            "The cell-0 metric/connection counterfactual is restricted to accepted_000: after evolution X and Lambda are independent evolved fields, so imposing initializer relations on later states would be invalid.",
            "The initial q0 versus -4*pi*S0^2 comparison isolates the scalar-pulse contribution in the initializer's first source sample; it does not substitute for a fully reconstructed regular radial integral.",
            "Stage histories distinguish the start-lapse-refreshed slice, explicit predictor, primary predictor, final implicit L2 operand, and accepted state. This does not supply internal sub-stages outside the documented production PIRK sequence.",
            "The accepted A_a/K update closure is checked only on r<=10 to avoid conflating the interior source budget with the finite outer-boundary reconstruction.",
            "No production physics, equation, boundary, coefficient, threshold, or scientific-admission flag is changed by this diagnostic.",
        ],
    }
    # Record the first digest separately after all source-preserving runs.
    init_record = next(row for row in records if row["label"] == "accepted_000")
    report["initial_state_sha256"] = init_record["state_sha256_before_observation"]
    report_path = OUT / "report.json"
    OUT.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    npz_path = OUT / "stage_profiles.npz"
    np.savez_compressed(npz_path, **out_profiles)
    checkpoint_path = OUT / "final_state.npz"
    kernel.checkpoint(state, checkpoint_path)
    sums = {}
    for path in (report_path, npz_path, checkpoint_path):
        sums[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    (OUT / "SHA256SUMS.json").write_text(json.dumps(sums, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    final_record = next(row for row in records if row["label"] == f"accepted_{steps:03d}")
    print(json.dumps({
        "status": report["status"],
        "resolution": RESOLUTION,
        "target_time": TARGET_TIME,
        "dt": dt,
        "accepted_steps": steps,
        "initial_Aa_cell0": init_record["regularity_fits_even_fields"]["Aa"]["cell0"],
        "initial_Aa_even_fit_origin": init_record["regularity_fits_even_fields"]["Aa"]["origin_intercept"],
        "final_Aa_cell0": final_record["regularity_fits_even_fields"]["Aa"]["cell0"],
        "final_Aa_even_fit_origin": final_record["regularity_fits_even_fields"]["Aa"]["origin_intercept"],
        "final_Aa_cell0_minus_fit": final_record["regularity_fits_even_fields"]["Aa"]["cell0_minus_fit_intercept"],
        "final_K_cell0_minus_fit": final_record["regularity_fits_even_fields"]["K"]["cell0_minus_fit_intercept"],
        "max_Aa_update_closure_r_le_10": max(x["Aa_update_closure_max_abs_r_le_10"] for x in step_budgets),
        "max_K_update_closure_r_le_10": max(x["K_update_closure_max_abs_r_le_10"] for x in step_budgets),
        "artifact_checksums": sums,
        "scientific_admission": False,
    }, indent=2, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
