#!/usr/bin/env python3
"""No-evolution audit of regular-centre geometry on production initial slices.

This observer compares the pinned BSSN Ricci source with an independent
spherical warped-product reconstruction from areal radius and proper radial
distance. It intentionally reports discrepancies without changing the engine
or imposing a scientific-admission threshold.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import math
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

SOURCE_COMMIT = "285a2ee0c1a3b97cddeb4519d8e749109b115aab"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_BLOBS = {
    "engine/production_kernel.py": "9685e76bb5a6f26d6bf383e240de24f37523f265",
    "engine/matter_system.py": "05e23d2368120ea8a1a6aad3fa438625d5380750",
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
CENTRE_CELLS = 6
OUT = Path("results/current-main-centre-proper-radial-ricci-highres-2026-10-15")


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def _arr(values) -> np.ndarray:
    return np.asarray(values, dtype=float).copy()


def _max_abs(values) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.max(np.abs(values))) if values.size else 0.0


def _regional_summary(r, fields, mask):
    """Return signed, scale-preserving summaries for a fixed radial region."""
    r = np.asarray(r, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    indices = np.flatnonzero(mask)
    if indices.size == 0:
        return {"cells": 0, "profiles": {}}
    output = {
        "cells": int(indices.size),
        "r_min": float(r[indices[0]]),
        "r_max": float(r[indices[-1]]),
        "profiles": {},
    }
    for name, values in fields.items():
        values = np.asarray(values, dtype=float)
        selected = values[indices]
        if not np.all(np.isfinite(selected)):
            output["profiles"][name] = {"finite": False}
            continue
        local = int(np.argmax(np.abs(selected)))
        source_index = int(indices[local])
        output["profiles"][name] = {
            "max_abs": float(abs(selected[local])),
            "signed_at_max_abs": float(selected[local]),
            "r_at_max_abs": float(r[source_index]),
            "min": float(np.min(selected)),
            "max": float(np.max(selected)),
            "rms": float(np.sqrt(np.mean(selected * selected))),
        }
    return output


def _finite(values) -> bool:
    return bool(np.all(np.isfinite(np.asarray(values, dtype=float))))


def _state_digest(state) -> str:
    digest = hashlib.sha256()
    arrays = {}
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        arrays["geometry." + name] = getattr(state.geometry, name)
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        arrays["scalars." + name] = getattr(state.scalars, name)
    for sector in ("dark_matter", "baryons", "radiation"):
        item = getattr(state.matter, sector)
        for name in ("rest", "energy_t", "momentum_r"):
            arrays["matter." + sector + "." + name] = getattr(item, name)
    for name in sorted(arrays):
        values = np.ascontiguousarray(np.asarray(arrays[name]))
        digest.update(name.encode("utf-8"))
        digest.update(str(values.dtype).encode("ascii"))
        digest.update(repr(values.shape).encode("ascii"))
        digest.update(values.tobytes())
    return digest.hexdigest()


def proper_radial_ricci(grid, areal_radius, proper_factor):
    """Compute mixed Ricci eigenvalues from dl=f(r)dr and R_areal(r).

    For ds^2 = dl^2 + R(l)^2 dOmega^2:
      Ric^l_l = -2 R_ll/R
      Ric^theta_theta = (1 - R_l^2 - R R_ll)/R^2.
    The construction does not use the evolved conformal connection Lambda or
    the vendor's BSSN Ricci brackets. Derivatives use the pinned grid operator
    with odd parity for R and even parity for dl/dr.
    """
    R = _arr(areal_radius)
    f = _arr(proper_factor)
    if R.shape != np.asarray(grid.centers).shape or f.shape != R.shape:
        raise ValueError("areal radius and proper factor must match the grid")
    if np.any(R <= 0.0) or np.any(f <= 0.0):
        raise ValueError("areal radius and proper factor must be positive")
    R_r = grid.cell_derivative_fourth(R, parity=-1)
    R_rr = grid.cell_second_derivative_fourth(R, parity=-1)
    f_r = grid.cell_derivative_fourth(f, parity=1)
    R_l = R_r / f
    R_ll = (R_rr - R_r * f_r / f) / (f * f)
    Ric_l_l = -2.0 * R_ll / R
    Ric_theta_theta = (1.0 - R_l * R_l - R * R_ll) / (R * R)
    tf = (2.0 / 3.0) * (Ric_l_l - Ric_theta_theta)
    return {
        "R_r": R_r,
        "R_rr": R_rr,
        "f_r": f_r,
        "R_l": R_l,
        "R_ll": R_ll,
        "Ric_l_l": Ric_l_l,
        "Ric_theta_theta": Ric_theta_theta,
        "tf_ricci": tf,
    }


def even_taylor_fit(r, values, count=CENTRE_CELLS):
    """Fit c0 + c2*r^2 + c4*r^4 on the first cell-centred samples."""
    r = _arr(r)
    values = _arr(values)
    if r.shape != values.shape or count < 3 or r.size < count:
        raise ValueError("invalid centre profile for even Taylor fit")
    rr = r[:count]
    yy = values[:count]
    scale = float(rr[-1])
    x = (rr / scale) ** 2
    design = np.column_stack((np.ones(count), x, x * x))
    coeff_scaled, *_ = np.linalg.lstsq(design, yy, rcond=None)
    fitted = design @ coeff_scaled
    residual = yy - fitted
    physical = [
        float(coeff_scaled[0]),
        float(coeff_scaled[1] / (scale * scale)),
        float(coeff_scaled[2] / (scale**4)),
    ]
    return {
        "basis": "c0 + c2*r^2 + c4*r^4",
        "fit_cells": int(count),
        "c0": physical[0],
        "c2": physical[1],
        "c4": physical[2],
        "max_abs_residual": _max_abs(residual),
        "max_abs_profile": _max_abs(yy),
        "signed_residuals": [float(x) for x in residual],
        "fitted_values": [float(x) for x in fitted],
    }



def _even_half_cell_integral(coefficients: np.ndarray, r0: float) -> float:
    """Integrate an even polynomial in r² analytically from 0 to r0."""
    coeff = np.asarray(coefficients, dtype=float).reshape(-1)
    if coeff.size not in (3, 4) or np.any(~np.isfinite(coeff)):
        raise ValueError("even half-cell integral expects finite s0/s2/s4[/s6] coefficients")
    return float(sum(value * r0 ** (2 * order + 1) / (2 * order + 1)
                     for order, value in enumerate(coeff)))


def _origin_half_cell_frozen_source(grid, state, vacuum) -> tuple[dict[str, np.ndarray], dict]:
    """Reconstruct the pinned initializer source and apply only an origin-cell moment correction.

    This is a *frozen-source counterfactual*: its source is taken from the final
    iteration of the production initializer recurrence, and the correction to
    the half-cell integral is propagated as delta_I/r to each cell. It does not
    iterate the source after changing B, evolve the machine, or edit production
    code. All inter-centre trapezoid intervals are preserved.
    """
    from engine.scalar_system import cosmos_potential
    from engine.v55_matter import RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT

    r = _arr(grid.centers)
    dr = float(grid.dr)
    geom = state.geometry
    fields = state.scalars
    S, PS, D, PD = (_arr(fields.S), _arr(fields.PS), _arr(fields.D), _arr(fields.PD))
    phi, Pi = _arr(fields.phi), _arr(fields.Pi)
    d1 = lambda v, parity=1: grid.cell_derivative_fourth(v, parity=parity)
    Sp, Dp = d1(S, 1), d1(D, 1)
    j = -(PD * Dp)

    # Reproduce the initializer's declared homogeneous Kt exactly rather
    # than recovering it by subtracting K/Aa (which would amplify roundoff in
    # Ktp on the large-r grid).
    phi0, pi0 = 0.179055, 0.00341913
    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += float(RHO_DM_PRESENT) + float(RHO_B_PRESENT) + float(RHO_R_PRESENT)
    H0 = math.sqrt(max(rho0 / 3.0, 0.0))
    Kt = np.full_like(r, -H0)
    Ktp = d1(Kt, 1)
    rho_matter = (
        float(RHO_DM_PRESENT) + float(RHO_B_PRESENT) + float(RHO_R_PRESENT)
    ) / (8.0 * math.pi)
    rho_phi = (0.5 * Pi * Pi + cosmos_potential(phi)) / (8.0 * math.pi)
    B_iter = np.ones_like(r)
    source_last = np.zeros_like(r)
    for _ in range(6):
        rhoS = 0.5 * (PS * PS + B_iter * Sp * Sp) + 0.5 * S * S
        rhoD = 0.5 * (PD * PD + B_iter * Dp * Dp) - 0.5 * D * D
        rho = rhoS + rhoD + rho_phi + rho_matter
        Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
        source_last = (
            1.0 - 8.0 * math.pi * r * r * rho
            + 2.0 * r * r * Kr * Kt
            + r * r * Kt * Kt
        )
        # Reproduce the production recurrence: 0.5*dr*src[0] == r0*src[0].
        integ = np.cumsum(
            np.r_[0.5 * dr * source_last[0],
                  0.5 * dr * (source_last[1:] + source_last[:-1])]
        )
        B_iter = integ / r

    B_from_geometry = _arr(geom.X) ** 6
    reconstruction_error = _max_abs(B_iter - B_from_geometry)
    r0 = float(r[0])

    def fit_even_source(fit_r_max: float) -> dict:
        # Fix support in physical radius, not cell count, so N=160/240/320
        # extrapolate the same source interval to the origin.
        mask = r <= float(fit_r_max)
        x = r[mask] ** 2
        ncells = int(np.count_nonzero(mask))
        if ncells < 4:
            raise ValueError("even-source fit needs at least four cells inside its fixed-radius support")
        design = np.column_stack((np.ones(ncells), x, x * x, x * x * x))
        coeff, *_ = np.linalg.lstsq(design, source_last[mask], rcond=None)
        residual = source_last[mask] - design @ coeff
        integral = _even_half_cell_integral(coeff, r0)
        delta = float(integral - r0 * source_last[0])
        return {
            "fit_r_max": float(fit_r_max),
            "fit_cells": ncells,
            "coefficients_s0_s2_s4_s6": [float(x) for x in coeff],
            "fit_max_abs_residual": _max_abs(residual),
            "fit_rms_residual": float(np.sqrt(np.mean(residual * residual))),
            "even_series_half_cell_integral": float(integral),
            "production_half_cell_integral": float(r0 * source_last[0]),
            "delta_integral": delta,
        }

    fit25 = fit_even_source(2.5)
    fit30 = fit_even_source(3.0)
    delta_integral = float(fit25["delta_integral"])
    B_cf = B_from_geometry + delta_integral / r
    if np.any(~np.isfinite(B_cf)) or np.any(B_cf <= 0.0):
        raise FloatingPointError("origin half-cell frozen-source counterfactual produced invalid B")
    A_cf = 1.0 / np.sqrt(B_cf)
    a_cf = A_cf ** (4.0 / 3.0)
    b_cf = A_cf ** (-2.0 / 3.0)
    X_cf = A_cf ** (-1.0 / 3.0)
    Lambda_cf = (
        d1(a_cf, 1) / (2.0 * a_cf * a_cf)
        - d1(b_cf, 1) / (a_cf * b_cf)
        + 2.0 / r * (1.0 / b_cf - 1.0 / a_cf)
    )

    # Recompute curvature on a deep copy. The production state and its matter,
    # K/Aa, and every inter-centre integral interval stay untouched.
    geom_cf = copy.deepcopy(geom)
    geom_cf.a = a_cf.copy()
    geom_cf.b = b_cf.copy()
    geom_cf.X = X_cf.copy()
    geom_cf.Lambda = Lambda_cf.copy()
    geom_cf.B = 0.75 * Lambda_cf
    vendor_terms_cf = vacuum.geometry_terms(grid, geom_cf)
    vendor_tf_cf = _arr(vendor_terms_cf["Rr_r"] - vendor_terms_cf["R"] / 3.0)
    areal_cf = r * np.sqrt(b_cf) / X_cf
    proper_factor_cf = np.sqrt(a_cf) / X_cf
    direct_cf = proper_radial_ricci(grid, areal_cf, proper_factor_cf)
    direct_tf_cf = _arr(direct_cf["tf_ricci"])
    gap_cf = vendor_tf_cf - direct_tf_cf

    vendor_terms_base = vacuum.geometry_terms(grid, geom)
    vendor_tf_base = _arr(vendor_terms_base["Rr_r"] - vendor_terms_base["R"] / 3.0)
    areal_base = r * np.sqrt(_arr(geom.b)) / _arr(geom.X)
    proper_factor_base = np.sqrt(_arr(geom.a)) / _arr(geom.X)
    direct_base = proper_radial_ricci(grid, areal_base, proper_factor_base)
    direct_tf_base = _arr(direct_base["tf_ricci"])
    gap_base = vendor_tf_base - direct_tf_base
    annulus = (r >= 5.0) & (r <= 10.0)

    first_two = []
    for i in range(min(2, len(r))):
        first_two.append({
            "i": i,
            "r": float(r[i]),
            "B_production": float(B_from_geometry[i]),
            "B_frozen_source_counterfactual": float(B_cf[i]),
            "X_production": float(geom.X[i]),
            "X_counterfactual": float(X_cf[i]),
            "vendor_tf_ricci_production": float(vendor_tf_base[i]),
            "direct_tf_ricci_production": float(direct_tf_base[i]),
            "signed_gap_production": float(gap_base[i]),
            "gap_over_r2_production": float(gap_base[i] / (r[i] * r[i])),
            "vendor_tf_ricci_counterfactual": float(vendor_tf_cf[i]),
            "direct_tf_ricci_counterfactual": float(direct_tf_cf[i]),
            "signed_gap_counterfactual": float(gap_cf[i]),
            "gap_over_r2_counterfactual": float(gap_cf[i] / (r[i] * r[i])),
        })
    arrays = {
        "initializer_source_final_iteration": source_last,
        "initializer_B_reconstructed": B_iter,
        "initializer_B_from_X6": B_from_geometry,
        "origin_cf_B": B_cf,
        "origin_cf_a": a_cf,
        "origin_cf_b": b_cf,
        "origin_cf_X": X_cf,
        "origin_cf_Lambda": Lambda_cf,
        "origin_cf_areal_radius": areal_cf,
        "origin_cf_proper_factor": proper_factor_cf,
        "origin_cf_vendor_tf_ricci": vendor_tf_cf,
        "origin_cf_direct_tf_ricci": direct_tf_cf,
        "origin_cf_tf_ricci_difference": gap_cf,
        "origin_cf_vendor_tf_ricci_over_r2": vendor_tf_cf / (r * r),
        "origin_cf_direct_tf_ricci_over_r2": direct_tf_cf / (r * r),
        "origin_cf_tf_ricci_difference_over_r2": gap_cf / (r * r),
    }
    summary = {
        "method": "frozen-source analytic even half-cell integral; final pinned initializer source; delta_I/r propagated to every cell; inter-centre trapezoids unchanged",
        "evolution_steps": 0,
        "source_reconstruction_B_vs_X6_max_abs": float(reconstruction_error),
        "source_reconstruction_matches_geometry": bool(reconstruction_error <= 2.0e-12 * max(1.0, _max_abs(B_from_geometry))),
        "source_final_iteration_first_two": [float(x) for x in source_last[:2]],
        "even_fit_r_le_2_5": fit25,
        "even_fit_r_le_3_0": fit30,
        "fit_sensitivity_delta_integral_abs": float(abs(fit30["delta_integral"] - fit25["delta_integral"])),
        "production_half_cell_integral": float(r0 * source_last[0]),
        "corrected_half_cell_integral": float(fit25["even_series_half_cell_integral"]),
        "delta_integral_applied": delta_integral,
        "first_two_cells": first_two,
        "fixed_annulus_r5_to10": {
            "cells": int(np.count_nonzero(annulus)),
            "production_gap_max_abs": _max_abs(gap_base[annulus]) if np.any(annulus) else None,
            "counterfactual_gap_max_abs": _max_abs(gap_cf[annulus]) if np.any(annulus) else None,
        },
        "first_two_cell_gap_reduction_factor": float(
            _max_abs(gap_base[:2]) / max(_max_abs(gap_cf[:2]), np.finfo(float).tiny)
        ),
        "counterfactual_profiles_finite": all(_finite(v) for v in arrays.values()),
        "interpretation_limits": [
            "This is a frozen-source counterfactual, not a self-consistent re-solve of the six-step initializer fixed point.",
            "Only the origin half-cell integral is changed; the same delta integral is propagated as delta_I/r, while all inter-centre trapezoid intervals are unchanged.",
            "Metric, Lambda, and Ricci are recomputed on a copied geometry; no production-state field is mutated and no evolution step is taken.",
            "A reduced Ricci gap would support the initialization quadrature hypothesis, not authorize promotion of the quadrature rule or prove evolved convergence."
        ],
    }
    if not summary["source_reconstruction_matches_geometry"]:
        raise RuntimeError(
            "reconstructed initializer B does not match the actual X^6 profile: "
            f"{reconstruction_error:.6g}"
        )
    return arrays, summary



def _lambda_chain_rule_counterfactual(
    grid,
    r: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    X: np.ndarray,
    Lambda: np.ndarray,
    vendor_terms: dict,
    direct_tf_ricci: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict]:
    """Offline-only Lambda derivative chain-rule comparison at fixed metric samples.

    The pinned vendor BSSN Ricci uses D4(Lambda) because Lambda is an independently
    evolved variable. This helper computes a counterfactual using the derivative
    of the displayed metric-to-Lambda relation, with the same discrete D1/D2
    operands for a and b. It does not alter state.Lambda or production outputs.
    """
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    d2 = lambda values, parity=1: grid.cell_second_derivative_fourth(values, parity=parity)
    ap = d1(a, 1)
    bp = d1(b, 1)
    app = d2(a, 1)
    bpp = d2(b, 1)
    lambda_p_discrete = d1(Lambda, -1)
    lambda_p_chain = (
        app / (2.0 * a * a)
        - ap * ap / (a * a * a)
        - bpp / (a * b)
        + bp * (ap * b + a * bp) / (a * a * b * b)
        - 2.0 * (1.0 / b - 1.0 / a) / (r * r)
        + 2.0 * (-bp / (b * b) + ap / (a * a)) / r
    )
    lambda_p_delta = lambda_p_chain - lambda_p_discrete
    ricci_scalar_shift = X * X * lambda_p_delta
    vendor_r_cf = _arr(vendor_terms["Rr_r"]) + ricci_scalar_shift
    vendor_scalar_cf = _arr(vendor_terms["R"]) + ricci_scalar_shift
    vendor_tf_cf = vendor_r_cf - vendor_scalar_cf / 3.0
    direct_tf = _arr(direct_tf_ricci)
    gap_base = _arr(vendor_terms["Rr_r"] - vendor_terms["R"] / 3.0) - direct_tf
    gap_cf = vendor_tf_cf - direct_tf
    fixed = (r >= 5.0) & (r <= 10.0)
    first_two = []
    for i in range(min(2, len(r))):
        first_two.append({
            "i": int(i),
            "r": float(r[i]),
            "lambda_prime_discrete_D4": float(lambda_p_discrete[i]),
            "lambda_prime_chain_rule": float(lambda_p_chain[i]),
            "lambda_prime_chain_minus_D4": float(lambda_p_delta[i]),
            "induced_trace_free_Ricci_shift": float((2.0 / 3.0) * ricci_scalar_shift[i]),
            "baseline_vendor_minus_direct_tf_ricci": float(gap_base[i]),
            "chain_rule_vendor_minus_direct_tf_ricci": float(gap_cf[i]),
            "baseline_gap_over_r2": float(gap_base[i] / (r[i] * r[i])),
            "chain_rule_gap_over_r2": float(gap_cf[i] / (r[i] * r[i])),
        })
    arrays = {
        "lambda_prime_discrete_D4": lambda_p_discrete,
        "lambda_prime_chain_rule": lambda_p_chain,
        "lambda_prime_chain_minus_D4": lambda_p_delta,
        "lambda_chain_rule_induced_tf_ricci_shift": (2.0 / 3.0) * ricci_scalar_shift,
        "lambda_chain_rule_vendor_tf_ricci": vendor_tf_cf,
        "lambda_chain_rule_vendor_minus_direct_tf_ricci": gap_cf,
        "lambda_chain_rule_vendor_minus_direct_tf_ricci_over_r2": gap_cf / (r * r),
    }
    summary = {
        "method": "offline chain-rule derivative of the displayed metric Lambda relation using the same D1/D2 a,b operands; only the -a*Lambda' Ricci contribution is counterfactually changed",
        "production_or_vendor_code_changed": False,
        "evolution_steps": 0,
        "lambda_prime_D4_max_abs": _max_abs(lambda_p_discrete),
        "lambda_prime_chain_rule_max_abs": _max_abs(lambda_p_chain),
        "lambda_prime_chain_minus_D4_max_abs": _max_abs(lambda_p_delta),
        "induced_tf_ricci_shift_max_abs": _max_abs((2.0 / 3.0) * ricci_scalar_shift),
        "baseline_gap_max_abs_full": _max_abs(gap_base),
        "chain_rule_gap_max_abs_full": _max_abs(gap_cf),
        "baseline_gap_max_abs_first_two": _max_abs(gap_base[:min(2, len(r))]),
        "chain_rule_gap_max_abs_first_two": _max_abs(gap_cf[:min(2, len(r))]),
        "baseline_gap_max_abs_fixed_annulus_r5_to10": _max_abs(gap_base[fixed]) if np.any(fixed) else None,
        "chain_rule_gap_max_abs_fixed_annulus_r5_to10": _max_abs(gap_cf[fixed]) if np.any(fixed) else None,
        "first_two_cells": first_two,
        "counterfactual_profiles_finite": all(_finite(v) for v in arrays.values()),
        "interpretation_limits": [
            "This is an offline representation counterfactual, not the production BSSN operator.",
            "Lambda is independently evolved in the vendor BSSN formulation; chain-rule substitution may not preserve its off-constraint principal structure.",
            "Gap collapse supports discrete product/chain-rule mismatch on this initialized slice, not a production fix or scientific admission.",
        ],
    }
    return arrays, summary


def _origin_half_cell_self_consistent_counterfactual(
    grid, state, vacuum, *, frozen_delta_integral: float,
    max_iterations: int = 100, tolerance: float = 2.0e-14,
    source_fit_degree: int = 3, source_fit_r_max: float = 2.5,
    include_fit_scan: bool = True,
) -> tuple[dict[str, np.ndarray], dict]:
    """Iterate corrected origin moment and its metric/source map on a copied state.

    Every iterate recomputes the Hamiltonian source from its current B metric,
    fits an even source series on fixed physical support, corrects only the
    origin half-cell integral, and preserves inter-centre trapezoids. Reports
    both a six-step counterpart to the pinned initializer and a converged
    fixed point. The initialized engine state is never modified or evolved.
    """
    if max_iterations < 6 or tolerance <= 0.0:
        raise ValueError("invalid self-consistent origin iteration controls")
    if source_fit_degree not in (3, 4, 5) or source_fit_r_max not in (2.5, 3.0, 3.5):
        raise ValueError("origin fit-family comparison expects degree 3/4/5 and radius 2.5/3.0/3.5")
    from engine.scalar_system import cosmos_potential
    from engine.v55_matter import RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT

    state_before = _state_digest(state)
    r = _arr(grid.centers)
    dr, r0 = float(grid.dr), float(grid.centers[0])
    geom, fields = state.geometry, state.scalars
    S, PS, D, PD = (_arr(fields.S), _arr(fields.PS), _arr(fields.D), _arr(fields.PD))
    phi, Pi = _arr(fields.phi), _arr(fields.Pi)
    d1 = lambda v, parity=1: grid.cell_derivative_fourth(v, parity=parity)
    Sp, Dp = d1(S, 1), d1(D, 1)
    current = -(PD * Dp)

    phi0, pi0 = 0.179055, 0.00341913
    rho0 = 0.5*pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += float(RHO_DM_PRESENT) + float(RHO_B_PRESENT) + float(RHO_R_PRESENT)
    H0 = math.sqrt(max(rho0/3.0, 0.0))
    Kt, Ktp = np.full_like(r, -H0), d1(np.full_like(r, -H0), 1)
    rho_matter = (float(RHO_DM_PRESENT) + float(RHO_B_PRESENT) + float(RHO_R_PRESENT))/(8.0*math.pi)
    rho_phi = (0.5*Pi*Pi + cosmos_potential(phi))/(8.0*math.pi)

    def source_for_B(B: np.ndarray) -> np.ndarray:
        rhoS = 0.5*(PS*PS + B*Sp*Sp) + 0.5*S*S
        rhoD = 0.5*(PD*PD + B*Dp*Dp) - 0.5*D*D
        rho = rhoS + rhoD + rho_phi + rho_matter
        Kr = Kt + r*Ktp + 4.0*math.pi*r*current
        source = 1.0 - 8.0*math.pi*r*r*rho + 2.0*r*r*Kr*Kt + r*r*Kt*Kt
        if not _finite(source):
            raise FloatingPointError("non-finite source in copied initializer recurrence")
        return source

    def fit_even(source: np.ndarray, fit_r_max: float, degree: int = source_fit_degree) -> dict:
        mask = r <= float(fit_r_max)
        x, nfit = r[mask]**2, int(np.count_nonzero(mask))
        if nfit < degree + 1:
            raise ValueError(
                f"even-source degree-{degree} fit at r<={fit_r_max} has only {nfit} samples"
            )
        design = np.column_stack([x**power for power in range(degree + 1)])
        coeff, *_ = np.linalg.lstsq(design, source[mask], rcond=None)
        residual = source[mask] - design@coeff
        integral = float(sum(
            value * r0 ** (2 * power + 1) / (2 * power + 1)
            for power, value in enumerate(coeff)
        ))
        result = {
            "fit_r_max": float(fit_r_max), "fit_cells": nfit,
            "degree_in_r_squared": int(degree),
            "overdetermined": bool(nfit >= degree + 2),
            "coefficients_s_even": [float(v) for v in coeff],
            "fit_max_abs_residual": _max_abs(residual),
            "fit_rms_residual": float(np.sqrt(np.mean(residual*residual))),
            "half_cell_integral": float(integral),
            "delta_integral_vs_point_sample": float(integral-r0*source[0]),
        }
        if degree == 3:
            result["coefficients_s0_s2_s4_s6"] = [float(v) for v in coeff]
        return result

    def integrate_original(source: np.ndarray) -> np.ndarray:
        return np.cumsum(np.r_[r0*source[0], 0.5*dr*(source[1:]+source[:-1])])/r

    def integrate_corrected(source: np.ndarray) -> tuple[np.ndarray, dict]:
        fit25 = fit_even(source, 2.5) if np.count_nonzero(r <= 2.5) >= source_fit_degree + 1 else None
        fit30 = fit_even(source, 3.0) if np.count_nonzero(r <= 3.0) >= source_fit_degree + 1 else None
        selected = fit_even(source, source_fit_r_max, source_fit_degree)
        alternative = fit30 if source_fit_r_max == 2.5 else fit25
        sensitivity = (
            abs(fit30["half_cell_integral"] - fit25["half_cell_integral"])
            if fit25 is not None and fit30 is not None else None
        )
        Bnext = np.cumsum(np.r_[
            selected["half_cell_integral"], 0.5*dr*(source[1:]+source[:-1])
        ])/r
        return Bnext, {
            "selected_fit": selected,
            "fit_r_le_2_5": fit25,
            "fit_r_le_3_0": fit30,
            "comparison_fit": alternative,
            "fit_window_integral_sensitivity_abs": float(sensitivity) if sensitivity is not None else None,
        }

    # Internal control: reconstruct the production point-sampled recurrence.
    B_baseline = np.ones_like(r)
    baseline_updates = []
    for _ in range(6):
        Bnext = integrate_original(source_for_B(B_baseline))
        baseline_updates.append(_max_abs(Bnext-B_baseline))
        B_baseline = Bnext
    B_production = _arr(geom.X)**6
    baseline_error = _max_abs(B_baseline-B_production)
    baseline_scale = max(1.0, _max_abs(B_production))
    if baseline_error > 2.0e-12*baseline_scale:
        raise RuntimeError(f"copied point-sampled recurrence failed production X^6 closure: {baseline_error:.6e}")

    # The corrected iteration starts at B=1 and recomputes the source every pass.
    B = np.ones_like(r)
    history: list[dict] = []
    B_six, fit_meta_six = None, None
    converged, last_update = False, float("inf")
    for iteration in range(1, max_iterations+1):
        src = source_for_B(B)
        Bnext, metadata = integrate_corrected(src)
        last_update = _max_abs(Bnext-B)
        scale = max(1.0, _max_abs(Bnext))
        history.append({
            "iteration": iteration, "max_abs_B_update": float(last_update),
            "max_relative_B_update": float(last_update/scale),
            "source_fit_degree_in_r_squared": int(source_fit_degree),
            "source_fit_r_max": float(source_fit_r_max),
            "delta_integral_applied": float(metadata["selected_fit"]["delta_integral_vs_point_sample"]),
            "selected_fit_max_abs_residual": float(metadata["selected_fit"]["fit_max_abs_residual"]),
            "fit_r_le_2_5_max_abs_residual": (
                float(metadata["fit_r_le_2_5"]["fit_max_abs_residual"])
                if metadata["fit_r_le_2_5"] is not None else None
            ),
            "fit_window_integral_sensitivity_abs": metadata["fit_window_integral_sensitivity_abs"],
        })
        B = Bnext
        if iteration == 6:
            B_six, fit_meta_six = B.copy(), metadata
        if iteration >= 6 and last_update <= tolerance*scale:
            converged = True
            break
    if B_six is None or fit_meta_six is None:
        raise RuntimeError("corrected-source iteration missed the six-step checkpoint")

    # Fresh source evaluations expose the fixed-point defects after iteration.
    src_six, src_fixed = source_for_B(B_six), source_for_B(B)
    B_six_map, six_meta = integrate_corrected(src_six)
    B_fixed_map, fixed_meta = integrate_corrected(src_fixed)
    six_defect, fixed_defect = _max_abs(B_six_map-B_six), _max_abs(B_fixed_map-B)

    def fit_sensitivity_scan(source: np.ndarray) -> dict:
        """Show how source-moment delta-I depends on polynomial degree/window."""
        rows = []
        for fit_r_max in (1.5, 2.0, 2.5, 3.0, 3.5):
            mask = r <= fit_r_max
            x = r[mask]**2
            nfit = int(np.count_nonzero(mask))
            for degree in (1, 2, 3, 4, 5):
                if nfit < degree + 1:
                    continue
                design = np.column_stack([x**power for power in range(degree + 1)])
                coeff, *_ = np.linalg.lstsq(design, source[mask], rcond=None)
                residual = source[mask] - design@coeff
                integral = float(sum(
                    value * r0 ** (2 * power + 1) / (2 * power + 1)
                    for power, value in enumerate(coeff)
                ))
                delta_i = float(integral - r0*source[0])
                loo = []
                if nfit >= degree + 2:
                    for omitted in range(nfit):
                        keep = np.arange(nfit) != omitted
                        loo_coeff, *_ = np.linalg.lstsq(design[keep], source[mask][keep], rcond=None)
                        loo.append(float(source[mask][omitted] - design[omitted]@loo_coeff))
                rows.append({
                    "fit_r_max": float(fit_r_max), "degree_in_r_squared": int(degree),
                    "fit_cells": nfit, "overdetermined": bool(nfit >= degree + 2),
                    "delta_integral": delta_i,
                    "fit_max_abs_residual": _max_abs(residual),
                    "fit_rms_residual": float(np.sqrt(np.mean(residual*residual))),
                    "leave_one_cell_out_rms_residual": float(np.sqrt(np.mean(np.square(loo)))) if loo else None,
                    "leave_one_cell_out_max_abs_residual": float(np.max(np.abs(loo))) if loo else None,
                })
        reference = next(
            row["delta_integral"] for row in rows
            if row["degree_in_r_squared"] == 3 and row["fit_r_max"] == 2.5
        )
        for row in rows:
            row["relative_delta_integral_percent_vs_degree3_r2_5"] = float(
                100.0*(row["delta_integral"]/reference - 1.0)
            )
        stable_degree3 = [
            row for row in rows
            if row["degree_in_r_squared"] == 3 and 1.5 <= row["fit_r_max"] <= 3.0
        ]
        all_degree3 = [row for row in rows if row["degree_in_r_squared"] == 3]
        overdetermined_higher = [
            row for row in rows
            if row["degree_in_r_squared"] in (3, 4, 5)
            and row["overdetermined"] and row["fit_r_max"] <= 3.0
        ]
        higher_spread = (
            100.0 * (max(row["delta_integral"] for row in overdetermined_higher)
                     - min(row["delta_integral"] for row in overdetermined_higher))
            / abs(reference)
        ) if overdetermined_higher else None
        return {
            "reference_fit": {"degree_in_r_squared": 3, "fit_r_max": 2.5,
                              "delta_integral": float(reference)},
            "fits": rows,
            "degree3_max_abs_delta_percent_over_windows_r1p5_to_r3": float(max(
                abs(row["relative_delta_integral_percent_vs_degree3_r2_5"])
                for row in stable_degree3
            )),
            "degree3_max_abs_delta_percent_over_windows_r1p5_to_r3p5": float(max(
                abs(row["relative_delta_integral_percent_vs_degree3_r2_5"])
                for row in all_degree3
            )),
            "degree3_to5_overdetermined_delta_spread_percent_windows_to_r3": higher_spread,
            "loo_qualified_fit_count": int(sum(row["leave_one_cell_out_rms_residual"] is not None for row in rows)),
        }

    fit_scan = fit_sensitivity_scan(src_fixed) if include_fit_scan else None

    def make_metric(Bprofile: np.ndarray, tag: str, source_profile: np.ndarray):
        if np.any(~np.isfinite(Bprofile)) or np.any(Bprofile <= 0.0):
            raise FloatingPointError(f"{tag}: invalid B profile")
        A = 1.0/np.sqrt(Bprofile)
        a, b, X = A**(4.0/3.0), A**(-2.0/3.0), A**(-1.0/3.0)
        Lambda = d1(a, 1)/(2.0*a*a) - d1(b, 1)/(a*b) + 2.0/r*(1.0/b-1.0/a)
        copied = copy.deepcopy(geom)
        copied.a, copied.b, copied.X, copied.Lambda = a.copy(), b.copy(), X.copy(), Lambda.copy()
        copied.B = 0.75*Lambda
        vt = vacuum.geometry_terms(grid, copied)
        vendor_tf = _arr(vt["Rr_r"]-vt["R"]/3.0)
        areal, proper_factor = r*np.sqrt(b)/X, np.sqrt(a)/X
        direct_tf = _arr(proper_radial_ricci(grid, areal, proper_factor)["tf_ricci"])
        gap = vendor_tf-direct_tf
        arrays = {
            f"{tag}_B": Bprofile.copy(), f"{tag}_a": a, f"{tag}_b": b, f"{tag}_X": X,
            f"{tag}_Lambda": Lambda, f"{tag}_areal_radius": areal,
            f"{tag}_proper_factor": proper_factor, f"{tag}_vendor_tf_ricci": vendor_tf,
            f"{tag}_direct_tf_ricci": direct_tf, f"{tag}_tf_ricci_difference": gap,
            f"{tag}_tf_ricci_difference_over_r2": gap/(r*r),
            f"{tag}_delta_over_r2": (1.0/b-1.0/a)/(r*r),
            f"{tag}_Lambda_over_r": Lambda/r,
            f"{tag}_source_recomputed_at_profile": source_profile,
        }
        invariants = {
            "a_times_b_squared_minus_one_max_abs": _max_abs(a*b*b-1.0),
            "b_minus_X_squared_max_abs": _max_abs(b-X*X),
            "a_minus_X_minus_four_max_abs": _max_abs(a-X**(-4.0)),
            "X6_minus_B_max_abs": _max_abs(X**6-Bprofile),
            "areal_radius_minus_r_max_abs": _max_abs(areal-r),
            "source_finite": _finite(source_profile),
            "all_profiles_finite": all(_finite(v) for v in arrays.values()),
        }
        cells = []
        for i in range(min(2, len(r))):
            cells.append({
                "i": i, "r": float(r[i]), "B": float(Bprofile[i]),
                "delta_over_r2": float(arrays[f"{tag}_delta_over_r2"][i]),
                "Lambda_over_r": float(arrays[f"{tag}_Lambda_over_r"][i]),
                "vendor_tf_ricci": float(vendor_tf[i]), "direct_tf_ricci": float(direct_tf[i]),
                "signed_gap": float(gap[i]), "signed_gap_over_r2": float(gap[i]/(r[i]*r[i])),
                "source_recomputed_at_profile": float(source_profile[i]),
            })
        return arrays, invariants, cells

    six_arrays, six_invariants, six_cells = make_metric(B_six, "origin_sc_six", src_six)
    fixed_arrays, fixed_invariants, fixed_cells = make_metric(B, "origin_sc_fixed", src_fixed)
    arrays = {**six_arrays, **fixed_arrays}
    arrays["origin_sc_baseline_production_B"] = B_production
    arrays["origin_sc_frozen_source_cf_B"] = B_production + float(frozen_delta_integral)/r
    arrays["origin_sc_initializer_source_baseline_final"] = source_for_B(B_production)

    baseline_vt = vacuum.geometry_terms(grid, geom)
    baseline_vendor = _arr(baseline_vt["Rr_r"]-baseline_vt["R"]/3.0)
    baseline_direct = _arr(proper_radial_ricci(
        grid, r*np.sqrt(_arr(geom.b))/_arr(geom.X), np.sqrt(_arr(geom.a))/_arr(geom.X)
    )["tf_ricci"])
    baseline_gap = baseline_vendor-baseline_direct
    annulus = (r >= 5.0) & (r <= 10.0)
    state_after = _state_digest(state)
    summary = {
        "method": "copied initializer recurrence: re-evaluate source at every B iterate; integrate its fitted even half-cell series; preserve inter-centre trapezoids",
        "production_or_vendor_code_changed": False, "evolution_steps": 0,
        "input_state_sha256_before": state_before, "input_state_sha256_after": state_after,
        "input_state_unchanged": state_before == state_after,
        "resolutions": {"N": int(len(r)), "r_max": float(r[-1]+dr/2.0), "dr": dr},
        "source_pins": {"engine_source_anchor": SOURCE_COMMIT, "vendor_pin": VENDOR_COMMIT},
        "baseline_six_iteration_reconstruction_max_abs": float(baseline_error),
        "baseline_six_iteration_reconstruction_matches_production": bool(baseline_error <= 2.0e-12*baseline_scale),
        "baseline_six_iteration_max_abs_updates": [float(v) for v in baseline_updates],
        "convergence": {
            "maximum_iterations": max_iterations, "tolerance": tolerance, "converged": converged,
            "iterations_used": len(history), "last_B_update_max_abs": float(last_update),
            "recomputed_six_step_fixed_point_defect_max_abs": float(six_defect),
            "recomputed_fixed_point_defect_max_abs": float(fixed_defect),
            "six_step_history": history[:6], "tail_history": history[-min(5, len(history)):],
        },
        "corrected_source": {
            "source_fit_degree_in_r_squared": int(source_fit_degree),
            "source_fit_r_max": float(source_fit_r_max),
            "six_step_delta_integral": float(fit_meta_six["selected_fit"]["delta_integral_vs_point_sample"]),
            "fixed_point_delta_integral": float(fixed_meta["selected_fit"]["delta_integral_vs_point_sample"]),
            "fixed_point_fit_window_sensitivity_abs": fixed_meta["fit_window_integral_sensitivity_abs"],
            "fixed_point_fit_selected": fixed_meta["selected_fit"],
            "fixed_point_fit_r_le_2_5": fixed_meta["fit_r_le_2_5"],
            "fixed_point_fit_r_le_3_0": fixed_meta["fit_r_le_3_0"],
            "fit_degree_window_sensitivity_scan": fit_scan,
        },
        "metric_invariants_six_step": six_invariants,
        "metric_invariants_fixed_point": fixed_invariants,
        "production_first_two_cells": [
            {"i": i, "r": float(r[i]), "B": float(B_production[i]),
             "delta_over_r2": float((1.0/_arr(geom.b)[i]-1.0/_arr(geom.a)[i])/(r[i]*r[i])),
             "Lambda_over_r": float(_arr(geom.Lambda)[i]/r[i]), "signed_gap": float(baseline_gap[i]),
             "signed_gap_over_r2": float(baseline_gap[i]/(r[i]*r[i]))}
            for i in range(min(2, len(r)))
        ],
        "six_step_corrected_first_two_cells": six_cells,
        "fixed_point_corrected_first_two_cells": fixed_cells,
        "fixed_annulus_r5_to10": {
            "production_gap_max_abs": _max_abs(baseline_gap[annulus]) if np.any(annulus) else None,
            "six_step_corrected_gap_max_abs": _max_abs(six_arrays["origin_sc_six_tf_ricci_difference"][annulus]) if np.any(annulus) else None,
            "fixed_point_corrected_gap_max_abs": _max_abs(fixed_arrays["origin_sc_fixed_tf_ricci_difference"][annulus]) if np.any(annulus) else None,
        },
        "first_two_cell_gap_reduction": {
            "six_step_vs_production": float(_max_abs(six_arrays["origin_sc_six_tf_ricci_difference"][:2])/max(_max_abs(baseline_gap[:2]), np.finfo(float).tiny)),
            "fixed_point_vs_production": float(_max_abs(fixed_arrays["origin_sc_fixed_tf_ricci_difference"][:2])/max(_max_abs(baseline_gap[:2]), np.finfo(float).tiny)),
        },
        "all_counterfactual_profiles_finite": all(_finite(v) for v in arrays.values()),
        "interpretation_limits": [
            "This is a copied-state diagnostic recurrence, not an edited production initializer.",
            "The baseline recurrence is reconstructed using the original point-sampled half-cell as a consistency control.",
            "The corrected source is refitted and recomputed from each candidate B; the converged fixed point is separate from the six-step initializer analogue.",
            "Scalar/matter fields stay at initialized values while their energy/current are reevaluated with the candidate metric. No field is evolved.",
            "A reduced Ricci mismatch supports an origin-quadrature hypothesis only; it does not authorize a production change or scientific admission.",
        ],
    }
    if not summary["input_state_unchanged"]:
        raise RuntimeError("self-consistent origin counterfactual mutated the input state")
    if not summary["all_counterfactual_profiles_finite"]:
        raise FloatingPointError("self-consistent origin counterfactual produced non-finite arrays")
    return arrays, summary


def _combined_origin_and_lambda_counterfactual(
    grid,
    r: np.ndarray,
    origin_cf_arrays: dict[str, np.ndarray],
) -> tuple[dict[str, np.ndarray], dict]:
    """Compose the two offline discriminators on the copied half-cell metric.

    First use the frozen-source origin half-cell metric from
    _origin_half_cell_frozen_source; then, on that copied metric only, compare
    D4(Lambda) with the chain-rule derivative of the metric-to-Lambda relation.
    No production state, source file, or vendor operator is changed.
    """
    r = _arr(r)
    a = _arr(origin_cf_arrays["origin_cf_a"])
    b = _arr(origin_cf_arrays["origin_cf_b"])
    X = _arr(origin_cf_arrays["origin_cf_X"])
    Lambda = _arr(origin_cf_arrays["origin_cf_Lambda"])
    vendor_tf_origin = _arr(origin_cf_arrays["origin_cf_vendor_tf_ricci"])
    direct_tf_origin = _arr(origin_cf_arrays["origin_cf_direct_tf_ricci"])
    gap_origin = _arr(origin_cf_arrays["origin_cf_tf_ricci_difference"])

    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    d2 = lambda values, parity=1: grid.cell_second_derivative_fourth(values, parity=parity)
    ap = d1(a, 1)
    bp = d1(b, 1)
    app = d2(a, 1)
    bpp = d2(b, 1)
    lambda_p_discrete = d1(Lambda, -1)
    lambda_p_chain = (
        app / (2.0 * a * a)
        - ap * ap / (a * a * a)
        - bpp / (a * b)
        + bp * (ap * b + a * bp) / (a * a * b * b)
        - 2.0 * (1.0 / b - 1.0 / a) / (r * r)
        + 2.0 * (-bp / (b * b) + ap / (a * a)) / r
    )
    lambda_p_delta = lambda_p_chain - lambda_p_discrete
    induced_tf_shift = (2.0 / 3.0) * X * X * lambda_p_delta
    combined_vendor_tf = vendor_tf_origin + induced_tf_shift
    combined_gap = combined_vendor_tf - direct_tf_origin
    fixed = (r >= 5.0) & (r <= 10.0)
    first_two = []
    for i in range(min(2, len(r))):
        first_two.append({
            "i": int(i),
            "r": float(r[i]),
            "half_cell_only_gap": float(gap_origin[i]),
            "chain_rule_induced_tf_ricci_shift_on_half_cell_metric": float(induced_tf_shift[i]),
            "combined_gap": float(combined_gap[i]),
            "half_cell_only_gap_over_r2": float(gap_origin[i] / (r[i] * r[i])),
            "combined_gap_over_r2": float(combined_gap[i] / (r[i] * r[i])),
        })
    arrays = {
        "combined_origin_lambda_cf_lambda_prime_discrete_D4": lambda_p_discrete,
        "combined_origin_lambda_cf_lambda_prime_chain_rule": lambda_p_chain,
        "combined_origin_lambda_cf_lambda_prime_chain_minus_D4": lambda_p_delta,
        "combined_origin_lambda_cf_induced_tf_ricci_shift": induced_tf_shift,
        "combined_origin_lambda_cf_vendor_tf_ricci": combined_vendor_tf,
        "combined_origin_lambda_cf_vendor_minus_direct_tf_ricci": combined_gap,
        "combined_origin_lambda_cf_vendor_minus_direct_tf_ricci_over_r2": combined_gap / (r * r),
    }
    summary = {
        "method": "compose the frozen-source half-cell metric counterfactual with an offline Lambda chain-rule derivative comparison; all sampled geometry values and operators remain diagnostic-only",
        "production_or_vendor_code_changed": False,
        "evolution_steps": 0,
        "half_cell_only_gap_max_abs_full": _max_abs(gap_origin),
        "combined_gap_max_abs_full": _max_abs(combined_gap),
        "half_cell_only_gap_max_abs_first_two": _max_abs(gap_origin[:min(2, len(r))]),
        "combined_gap_max_abs_first_two": _max_abs(combined_gap[:min(2, len(r))]),
        "half_cell_only_gap_max_abs_fixed_annulus_r5_to10": _max_abs(gap_origin[fixed]) if np.any(fixed) else None,
        "combined_gap_max_abs_fixed_annulus_r5_to10": _max_abs(combined_gap[fixed]) if np.any(fixed) else None,
        "first_two_cells": first_two,
        "counterfactual_profiles_finite": all(_finite(v) for v in arrays.values()) and _finite(gap_origin),
        "interpretation_limits": [
            "This composes two frozen-slice counterfactuals; it is not a self-consistent initializer re-solve.",
            "The Lambda chain-rule substitution is not the production BSSN operator and may alter the off-constraint principal structure if promoted to evolution.",
            "A small combined gap localizes representation compatibility on this initialized slice; it does not prove a unique cause, evolved stability, or scientific admission.",
        ],
    }
    return arrays, summary



def _manufactured_lambda_parity_discriminator(grid, vacuum) -> tuple[dict[str, np.ndarray], dict]:
    """Compare continuum, discrete metric-derived, and off-constraint Lambda arms.

    Two smooth-even manufactured metric profiles are sampled on the existing
    cell-centred grid. Their radial contracted connection Lambda is odd. For
    each profile, the same fixed metric is passed through the pinned vendor
    Ricci route with (A) the analytic continuum metric-derived Lambda sampled
    on the grid, (B) the initializer-style discrete metric-derived Lambda,
    (C) an odd off-constraint perturbation of B, and (D) a separate even
    parity-breaking negative control. This is a no-evolution diagnostic only.
    """
    r = _arr(grid.centers)
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    profiles = (
        {"tag": "manufactured_q1", "c3": 1.0e-5, "c5": 2.0e-8, "epsilon": 2.0e-5, "width": 4.0},
        {"tag": "manufactured_q2", "c3": -8.0e-6, "c5": 1.0e-8, "epsilon": -1.5e-5, "width": 3.5},
    )
    annulus = (r >= 5.0) & (r <= 10.0)
    arrays: dict[str, np.ndarray] = {}
    profile_reports = []

    for spec in profiles:
        tag = spec["tag"]
        c3 = float(spec["c3"])
        c5 = float(spec["c5"])
        epsilon = float(spec["epsilon"])
        width = float(spec["width"])
        q = 1.0 + c3 * r * r + c5 * r**4
        if np.any(~np.isfinite(q)) or np.any(q <= 0.0):
            raise FloatingPointError(f"{tag}: manufactured q is not finite and positive")
        qp = 2.0 * c3 * r + 4.0 * c5 * r**3
        # Preserve the production metric identities on the manufactured slice:
        # b = X^2, a = X^-4, and a*b^2 = 1.
        a = q ** (-4.0 / 3.0)
        b = q ** (2.0 / 3.0)
        X = q ** (1.0 / 3.0)
        ap_cont = (-4.0 / 3.0) * a * qp / q
        bp_cont = (2.0 / 3.0) * b * qp / q
        lambda_cont = (
            ap_cont / (2.0 * a * a)
            - bp_cont / (a * b)
            + 2.0 / r * (1.0 / b - 1.0 / a)
        )
        lambda_discrete = (
            d1(a, 1) / (2.0 * a * a)
            - d1(b, 1) / (a * b)
            + 2.0 / r * (1.0 / b - 1.0 / a)
        )
        odd_perturbation = epsilon * r * np.exp(-((r / width) ** 2))
        arms = {
            "continuum": lambda_cont,
            "discrete": lambda_discrete,
            "off_constraint_odd": lambda_discrete + odd_perturbation,
            "negative_control_even": lambda_discrete + epsilon,
        }
        areal = r * np.sqrt(b) / X
        proper_factor = np.sqrt(a) / X
        expected_proper_factor = 1.0 / q
        direct = proper_radial_ricci(grid, areal, proper_factor)
        direct_tf = _arr(direct["tf_ricci"])
        ap_discrete = _arr(d1(a, 1))
        # Reference response for the initializer-style discrete Lambda arm.
        reference_geo = vacuum.flat_state(grid)
        reference_geo.a[:] = a
        reference_geo.b[:] = b
        reference_geo.X[:] = X
        reference_geo.alpha[:] = 1.0
        reference_geo.beta[:] = 0.0
        reference_geo.Lambda[:] = lambda_discrete
        reference_terms = vacuum.geometry_terms(grid, reference_geo)
        reference_tf = _arr(reference_terms["Rr_r"] - reference_terms["R"] / 3.0)
        reference_gap = reference_tf - direct_tf
        arm_reports = []
        for arm_name, lambda_values in arms.items():
            geo = vacuum.flat_state(grid)
            geo.a[:] = a
            geo.b[:] = b
            geo.X[:] = X
            geo.alpha[:] = 1.0
            geo.beta[:] = 0.0
            geo.Lambda[:] = lambda_values
            vendor_terms = vacuum.geometry_terms(grid, geo)
            vendor_tf = _arr(vendor_terms["Rr_r"] - vendor_terms["R"] / 3.0)
            gap = vendor_tf - direct_tf
            chain_arrays, chain_summary = _lambda_chain_rule_counterfactual(
                grid, r, a, b, X, _arr(lambda_values), vendor_terms, direct_tf
            )

            # Exact fixed-metric response identity of the pinned vendor TF-Ricci:
            # delta(Rr_r - R/3) = (2/3) X^2 D4(delta Lambda)
            #   + (1/2) X^2 [D1(a)/a] delta Lambda.
            # Keep the derivative and algebraic Lambda*a' channels separate.
            lambda_delta = _arr(lambda_values) - _arr(lambda_discrete)
            d4_lambda_delta = _arr(d1(lambda_delta, -1))
            observed_gap_delta = gap - reference_gap
            derivative_response = (2.0 / 3.0) * X * X * d4_lambda_delta
            algebraic_response = 0.5 * X * X * (ap_discrete / a) * lambda_delta
            predicted_response = derivative_response + algebraic_response
            response_closure_error = observed_gap_delta - predicted_response

            prefix = f"{tag}_{arm_name}"
            arrays[prefix + "_a"] = _arr(a).copy()
            arrays[prefix + "_b"] = _arr(b).copy()
            arrays[prefix + "_X"] = _arr(X).copy()
            arrays[prefix + "_areal_radius"] = _arr(areal).copy()
            arrays[prefix + "_proper_factor"] = _arr(proper_factor).copy()
            arrays[prefix + "_Lambda"] = _arr(lambda_values)
            arrays[prefix + "_Lambda_over_r"] = _arr(lambda_values) / r
            arrays[prefix + "_vendor_Rr_r"] = _arr(vendor_terms["Rr_r"])
            arrays[prefix + "_vendor_R_scalar"] = _arr(vendor_terms["R"])
            arrays[prefix + "_vendor_tf_ricci"] = vendor_tf
            arrays[prefix + "_vendor_tf_ricci_over_r2"] = vendor_tf / (r * r)
            arrays[prefix + "_signed_gap"] = gap
            arrays[prefix + "_signed_gap_over_r2"] = gap / (r * r)
            arrays[prefix + "_chain_rule_gap"] = _arr(
                chain_arrays["lambda_chain_rule_vendor_minus_direct_tf_ricci"]
            )
            arrays[prefix + "_chain_rule_induced_tf_shift"] = _arr(
                chain_arrays["lambda_chain_rule_induced_tf_ricci_shift"]
            )
            arrays[prefix + "_delta_Lambda_from_discrete"] = lambda_delta.copy()
            arrays[prefix + "_D4_delta_Lambda_from_discrete"] = d4_lambda_delta.copy()
            arrays[prefix + "_observed_gap_change_from_discrete"] = observed_gap_delta.copy()
            arrays[prefix + "_predicted_gap_change_derivative_only"] = derivative_response.copy()
            arrays[prefix + "_predicted_gap_change_algebraic_Lambda_ap"] = algebraic_response.copy()
            arrays[prefix + "_predicted_gap_change_total"] = predicted_response.copy()
            arrays[prefix + "_response_closure_error"] = response_closure_error.copy()
            annulus_signed = {}
            for stat_name, stat_values in {
                "vendor_Rr_r": arrays[prefix + "_vendor_Rr_r"],
                "vendor_R_scalar": arrays[prefix + "_vendor_R_scalar"],
                "vendor_tf_ricci": vendor_tf,
                "direct_tf_ricci": direct_tf,
                "signed_gap": gap,
                "signed_gap_over_r2": gap / (r * r),
                "observed_gap_change_from_discrete": observed_gap_delta,
                "predicted_gap_change_derivative_only": derivative_response,
                "predicted_gap_change_algebraic_Lambda_ap": algebraic_response,
                "predicted_gap_change_total": predicted_response,
                "response_closure_error": response_closure_error,
            }.items():
                selected = _arr(stat_values)[annulus]
                annulus_signed[stat_name] = {
                    "signed_min": float(np.min(selected)) if selected.size else None,
                    "signed_max": float(np.max(selected)) if selected.size else None,
                    "max_abs": _max_abs(selected) if selected.size else None,
                    "cells": int(selected.size),
                }
            arm_reports.append({
                "arm": arm_name,
                "lambda_parity_witness_max_abs_first6_lambda_over_r": _max_abs(
                    (_arr(lambda_values) / r)[:CENTRE_CELLS]
                ),
                "lambda_first_two": [float(v) for v in _arr(lambda_values)[:2]],
                "vendor_Rr_r_first_two": [float(v) for v in _arr(vendor_terms["Rr_r"])[:2]],
                "vendor_R_scalar_first_two": [float(v) for v in _arr(vendor_terms["R"])[:2]],
                "vendor_tf_ricci_first_two": [float(v) for v in vendor_tf[:2]],
                "direct_tf_ricci_first_two": [float(v) for v in direct_tf[:2]],
                "signed_gap_first_two": [float(v) for v in gap[:2]],
                "signed_gap_over_r2_first_two": [
                    float(v) for v in (gap / (r * r))[:2]
                ],
                "observed_gap_change_from_discrete_first_two": [
                    float(v) for v in observed_gap_delta[:2]
                ],
                "predicted_gap_change_derivative_only_first_two": [
                    float(v) for v in derivative_response[:2]
                ],
                "predicted_gap_change_algebraic_Lambda_ap_first_two": [
                    float(v) for v in algebraic_response[:2]
                ],
                "predicted_gap_change_total_first_two": [
                    float(v) for v in predicted_response[:2]
                ],
                "response_closure_error_first_two": [
                    float(v) for v in response_closure_error[:2]
                ],
                "gap_max_abs_first_two": _max_abs(gap[:2]),
                "gap_max_abs_fixed_annulus_r5_to10": (
                    _max_abs(gap[annulus]) if np.any(annulus) else None
                ),
                "fixed_annulus_signed_budget": annulus_signed,
                "response_closure_error_max_abs_full": _max_abs(response_closure_error),
                "response_closure_error_max_abs_first_two": _max_abs(response_closure_error[:2]),
                "response_closure_error_max_abs_fixed_annulus": (
                    _max_abs(response_closure_error[annulus]) if np.any(annulus) else None
                ),
                "chain_rule_gap_max_abs_full": chain_summary["chain_rule_gap_max_abs_full"],
                "chain_rule_gap_max_abs_fixed_annulus_r5_to10": (
                    chain_summary["chain_rule_gap_max_abs_fixed_annulus_r5_to10"]
                ),
                "profiles_finite": all(_finite(v) for v in (
                    lambda_values, vendor_tf, direct_tf, gap,
                    _arr(vendor_terms["Rr_r"]), _arr(vendor_terms["R"]),
                    observed_gap_delta, derivative_response, algebraic_response,
                    predicted_response, response_closure_error,
                    chain_arrays["lambda_chain_rule_vendor_minus_direct_tf_ricci"],
                )),
            })
        profile_reports.append({
            "tag": tag,
            "c3": c3,
            "c5": c5,
            "determinant_one_max_abs": _max_abs(a * b * b - 1.0),
            "b_minus_X_squared_max_abs": _max_abs(b - X * X),
            "a_minus_X_minus_four_max_abs": _max_abs(a - X ** (-4.0)),
            "proper_factor_minus_q_inverse_max_abs": _max_abs(proper_factor - expected_proper_factor),
            "areal_radius_minus_r_max_abs": _max_abs(areal - r),
            "off_constraint_odd_amplitude": epsilon,
            "off_constraint_width": width,
            "arms": arm_reports,
        })

    all_finite = all(_finite(values) for values in arrays.values()) and all(
        arm["profiles_finite"]
        for profile in profile_reports
        for arm in profile["arms"]
    )
    summary = {
        "method": "two smooth-even manufactured determinant-one metric profiles; continuum metric Lambda, discrete metric Lambda, odd off-constraint Lambda perturbation, and separate even-parity negative control; no evolution",
        "evolution_steps": 0,
        "production_or_vendor_code_changed": False,
        "profiles_finite": all_finite,
        "fixed_annulus_coordinate_r5_to10": {"cells": int(np.count_nonzero(annulus))},
        "profiles": profile_reports,
        "interpretation_limits": [
            "Lambda has odd parity for smooth-even spherical metric components; the off-constraint perturbation is smooth-odd and the parity-breaking control is deliberately separate.",
            "The metric and proper-radial profile are identical across Lambda arms; only the supplied connection field changes.",
            "The metric-chain-rule derivative is a diagnostic comparison, not the evolved BSSN Lambda equation and not a production substitution recommendation.",
            "The Lambda response closure includes both D4(delta Lambda) and the algebraic Lambda*a' term in the pinned trace-free Ricci expression.",
            "Manufactured profiles test discrete operator sensitivity and parity, not the cause of the initialized production mismatch or scientific admission.",
        ],
    }
    if not all_finite:
        raise FloatingPointError("manufactured Lambda parity discriminator produced non-finite values")
    return arrays, summary


def _verify_source_pins() -> dict:
    source_head = _git("rev-parse", "HEAD")
    source_vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if source_vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor pin mismatch: expected {VENDOR_COMMIT}, found {source_vendor}")
    pins = {}
    for path, expected in SOURCE_BLOBS.items():
        baseline_blob = _git("rev-parse", f"{SOURCE_COMMIT}:{path}")
        head_blob = _git("rev-parse", f"HEAD:{path}")
        pins[path] = {
            "expected_blob": expected,
            "baseline_blob": baseline_blob,
            "audit_head_blob": head_blob,
            "matches": baseline_blob == expected and head_blob == expected,
        }
        if not pins[path]["matches"]:
            raise RuntimeError(f"source blob mismatch for {path}: {pins[path]}")
    changed = _git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--",
        "engine", "vendor", ".gitmodules", "requirements.txt",
    )
    if changed:
        raise RuntimeError("engine/vendor/dependency paths changed from source pin:\n" + changed)
    return {
        "audit_head": source_head,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": source_vendor,
        "source_blobs": pins,
        "engine_vendor_dependency_diff_empty": True,
    }


def main() -> int:
    provenance = _verify_source_pins()

    from engine import v55_pirk_adapter as adapter
    from engine.production_kernel import V55ProductionKernel

    _, vacuum, _ = adapter.vendor_modules()
    OUT.mkdir(parents=True, exist_ok=True)
    reports = []
    profiles: dict[str, np.ndarray] = {}
    problems = []

    for resolution in RESOLUTIONS:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=resolution,
            r_max=R_MAX,
            amplitude=0.01,
            width=7.0,
            D_amplitude=1.0e-10,
            include_radiation=True,
        )
        grid = state.grid
        g = state.geometry
        r = _arr(grid.centers)
        before = _state_digest(state)

        a, b, X = _arr(g.a), _arr(g.b), _arr(g.X)
        alpha, Lambda = _arr(g.alpha), _arr(g.Lambda)
        areal = r * np.sqrt(b) / X
        proper_factor = np.sqrt(a) / X
        direct = proper_radial_ricci(grid, areal, proper_factor)
        vendor_terms = vacuum.geometry_terms(grid, g)
        vendor_tf = _arr(vendor_terms["Rr_r"] - vendor_terms["R"] / 3.0)
        vendor_l2_tf = alpha * vendor_tf
        direct_l2_tf = alpha * direct["tf_ricci"]
        tf_delta = vendor_tf - direct["tf_ricci"]
        l2_delta = vendor_l2_tf - direct_l2_tf
        lambda_chain_arrays, lambda_chain_summary = _lambda_chain_rule_counterfactual(
            grid, r, a, b, X, Lambda, vendor_terms, direct["tf_ricci"]
        )
        manufactured_lambda_arrays, manufactured_lambda_summary = (
            _manufactured_lambda_parity_discriminator(grid, vacuum)
        )
        if not manufactured_lambda_summary["profiles_finite"]:
            problems.append(
                f"N={resolution}: manufactured Lambda parity discriminator is non-finite"
            )

        # Signed source budget on the initialized profile; keep each driver
        # separate so a local pole is not hidden by summing unlike terms.
        aa_lapse_hessian = -(
            _arr(vendor_terms["DrDr_alpha"])
            - _arr(vendor_terms["lap_alpha"]) / 3.0
        )
        aa_ricci_tf = alpha * vendor_tf
        l3_vacuum = vacuum.primary_l3_rhs(grid, g)
        l3_with_matter = adapter.primary_l3_with_matter(
            grid, g, state.scalars, state.matter
        )
        aa_l3_vacuum = _arr(l3_vacuum["Aa"])
        aa_l3_matter = _arr(l3_with_matter["Aa"]) - aa_l3_vacuum
        aa_l2 = _arr(aa_lapse_hessian + aa_ricci_tf)
        aa_rhs = aa_l2 + _arr(l3_with_matter["Aa"])
        aa_l2_contract = _arr(vacuum.primary_l2_rhs(grid, g)["Aa"])
        aa_l2_reconstruction_error = aa_l2_contract - aa_l2
        aa_rhs_over_r2 = aa_rhs / (r * r)
        vendor_tf_over_r2 = vendor_tf / (r * r)
        direct_tf_over_r2 = _arr(direct["tf_ricci"]) / (r * r)
        tf_delta_over_r2 = tf_delta / (r * r)

        # Reconstruct the origin half-cell moment with a frozen final source.
        # This counterfactual keeps the actual initialized state read-only.
        origin_cf_arrays, origin_cf_summary = _origin_half_cell_frozen_source(
            grid, state, vacuum
        )
        origin_sc_arrays, origin_sc_summary = _origin_half_cell_self_consistent_counterfactual(
            grid, state, vacuum,
            frozen_delta_integral=float(origin_cf_summary["delta_integral_applied"]),
        )
        # Bound the origin-moment effect on several defensible, overdetermined
        # source-fit families. Each variant independently re-solves the copied
        # source/metric recurrence and recomputes its geometry/Ricci witnesses.
        origin_sc_fit_family_arrays = {}
        origin_sc_fit_family_summaries = {}
        fit_family_specs = (
            ("degree3_r3", 3, 3.0),
            ("degree4_r3", 4, 3.0),
            ("degree5_r3p5", 5, 3.5),
        )
        for family_tag, family_degree, family_rmax in fit_family_specs:
            family_arrays, family_summary = _origin_half_cell_self_consistent_counterfactual(
                grid, state, vacuum,
                frozen_delta_integral=float(origin_cf_summary["delta_integral_applied"]),
                source_fit_degree=family_degree,
                source_fit_r_max=family_rmax,
                include_fit_scan=False,
            )
            if not family_summary["convergence"]["converged"]:
                problems.append(
                    f"N={resolution}/{family_tag}: corrected-source fit-family iteration did not converge"
                )
            if not family_summary["input_state_unchanged"]:
                problems.append(
                    f"N={resolution}/{family_tag}: fit-family counterfactual mutated input state"
                )
            if not family_summary["all_counterfactual_profiles_finite"]:
                problems.append(
                    f"N={resolution}/{family_tag}: fit-family counterfactual produced non-finite values"
                )
            if not family_summary["baseline_six_iteration_reconstruction_matches_production"]:
                problems.append(
                    f"N={resolution}/{family_tag}: baseline recurrence failed to reproduce production X^6"
                )
            invariant_names = (
                "a_times_b_squared_minus_one_max_abs", "b_minus_X_squared_max_abs",
                "a_minus_X_minus_four_max_abs", "X6_minus_B_max_abs",
                "areal_radius_minus_r_max_abs",
            )
            for invariant_stage in ("metric_invariants_six_step", "metric_invariants_fixed_point"):
                invariants = family_summary[invariant_stage]
                if not invariants["all_profiles_finite"] or not invariants["source_finite"]:
                    problems.append(f"N={resolution}/{family_tag}: {invariant_stage} has non-finite values")
                for invariant_name in invariant_names:
                    if invariants[invariant_name] > 1.0e-12:
                        problems.append(
                            f"N={resolution}/{family_tag}: {invariant_stage}/{invariant_name} "
                            f"exceeds copied-metric identity tolerance"
                        )
            for array_name, values in family_arrays.items():
                if array_name.startswith("origin_sc_six_"):
                    output_name = (
                        f"origin_sc_{family_tag}_six_" +
                        array_name[len("origin_sc_six_"):]
                    )
                elif array_name.startswith("origin_sc_fixed_"):
                    output_name = (
                        f"origin_sc_{family_tag}_fixed_" +
                        array_name[len("origin_sc_fixed_"):]
                    )
                else:
                    output_name = (
                        f"origin_sc_{family_tag}_" +
                        array_name[len("origin_sc_"):]
                    )
                origin_sc_fit_family_arrays[output_name] = _arr(values)
            origin_sc_fit_family_summaries[family_tag] = {
                "source_fit_degree_in_r_squared": family_degree,
                "source_fit_r_max": family_rmax,
                "baseline_six_iteration_reconstruction_max_abs":
                    family_summary["baseline_six_iteration_reconstruction_max_abs"],
                "convergence": family_summary["convergence"],
                "corrected_source": family_summary["corrected_source"],
                "metric_invariants_six_step": family_summary["metric_invariants_six_step"],
                "metric_invariants_fixed_point": family_summary["metric_invariants_fixed_point"],
                "production_first_two_cells": family_summary["production_first_two_cells"],
                "six_step_corrected_first_two_cells": family_summary["six_step_corrected_first_two_cells"],
                "fixed_point_corrected_first_two_cells": family_summary["fixed_point_corrected_first_two_cells"],
                "fixed_annulus_r5_to10": family_summary["fixed_annulus_r5_to10"],
                "first_two_cell_gap_reduction": family_summary["first_two_cell_gap_reduction"],
                "input_state_unchanged": family_summary["input_state_unchanged"],
                "all_counterfactual_profiles_finite": family_summary["all_counterfactual_profiles_finite"],
            }
        # Re-test the independent Lambda derivative counterfactual on the
        # self-consistent origin profile too, rather than inferring it only
        # from the frozen-source metric.
        self_consistent_metric_as_origin_cf = {
            "origin_cf_a": origin_sc_arrays["origin_sc_fixed_a"],
            "origin_cf_b": origin_sc_arrays["origin_sc_fixed_b"],
            "origin_cf_X": origin_sc_arrays["origin_sc_fixed_X"],
            "origin_cf_Lambda": origin_sc_arrays["origin_sc_fixed_Lambda"],
            "origin_cf_vendor_tf_ricci": origin_sc_arrays["origin_sc_fixed_vendor_tf_ricci"],
            "origin_cf_direct_tf_ricci": origin_sc_arrays["origin_sc_fixed_direct_tf_ricci"],
            "origin_cf_tf_ricci_difference": origin_sc_arrays["origin_sc_fixed_tf_ricci_difference"],
        }
        combined_sc_raw_arrays, combined_sc_summary = _combined_origin_and_lambda_counterfactual(
            grid, r, self_consistent_metric_as_origin_cf
        )
        combined_sc_arrays = {
            name.replace("combined_origin_lambda_cf_", "combined_origin_sc_lambda_cf_"): values
            for name, values in combined_sc_raw_arrays.items()
        }
        if not origin_sc_summary["convergence"]["converged"]:
            problems.append(f"N={resolution}: corrected-source origin iteration did not converge")
        if not origin_sc_summary["all_counterfactual_profiles_finite"]:
            problems.append(f"N={resolution}: corrected-source origin counterfactual has non-finite profiles")
        combined_cf_arrays, combined_cf_summary = _combined_origin_and_lambda_counterfactual(
            grid, r, origin_cf_arrays
        )
        if not combined_cf_summary["counterfactual_profiles_finite"]:
            problems.append(f"N={resolution}: combined counterfactual produced non-finite profiles")

        aa_budget_scale = max(
            _max_abs(aa_l2),
            _max_abs(aa_l2_contract),
            _max_abs(aa_lapse_hessian),
            _max_abs(aa_ricci_tf),
            np.finfo(float).tiny,
        )
        aa_l2_roundoff_bound = 512.0 * np.finfo(np.float64).eps * aa_budget_scale
        if _max_abs(aa_l2_reconstruction_error) > aa_l2_roundoff_bound:
            problems.append(
                f"N={resolution}: signed L2 Aa source pieces do not close "
                "at the operand-scaled roundoff bound"
            )

        d1 = lambda v, parity=1: grid.cell_derivative_fourth(v, parity=parity)
        lambda_from_metric = (
            d1(a, 1) / (2.0 * a * a)
            - d1(b, 1) / (a * b)
            + 2.0 / r * (1.0 / b - 1.0 / a)
        )
        a_minus_b_over_r2 = (a - b) / (r * r)
        lambda_over_r = Lambda / r

        metric_arrays = {
            "a": a,
            "b": b,
            "X": X,
            "alpha": alpha,
            "Lambda": Lambda,
            "a_minus_b_over_r2": a_minus_b_over_r2,
            "Lambda_over_r": lambda_over_r,
            "areal_radius": areal,
            "proper_factor": proper_factor,
            "vendor_tf_ricci": vendor_tf,
            "direct_tf_ricci": direct["tf_ricci"],
            "tf_ricci_difference": tf_delta,
            "vendor_tf_ricci_over_r2": vendor_tf_over_r2,
            "direct_tf_ricci_over_r2": direct_tf_over_r2,
            "tf_ricci_difference_over_r2": tf_delta_over_r2,
            "vendor_l2_tf_ricci": vendor_l2_tf,
            "direct_l2_tf_ricci": direct_l2_tf,
            "l2_tf_ricci_difference": l2_delta,
            "aa_lapse_hessian_tf": aa_lapse_hessian,
            "aa_ricci_tf": aa_ricci_tf,
            "aa_l3_vacuum": aa_l3_vacuum,
            "aa_l3_matter": aa_l3_matter,
            "aa_l2": aa_l2,
            "aa_l2_reconstruction_error": aa_l2_reconstruction_error,
            "aa_rhs": aa_rhs,
            "aa_rhs_over_r2": aa_rhs_over_r2,
            "R_r": direct["R_r"],
            "R_rr": direct["R_rr"],
            "f_r": direct["f_r"],
            "R_l": direct["R_l"],
            "R_ll": direct["R_ll"],
            "Ric_l_l": direct["Ric_l_l"],
            "Ric_theta_theta": direct["Ric_theta_theta"],
            "a_prime_over_r": d1(a, 1) / r,
            "b_prime_over_r": d1(b, 1) / r,
            "X_prime_over_r": d1(X, 1) / r,
            "alpha_prime_over_r": d1(alpha, 1) / r,
            "lambda_metric_identity_residual": Lambda - lambda_from_metric,
            "areal_radius_minus_coordinate_radius": areal - r,
        }
        metric_arrays.update({name: _arr(values) for name, values in lambda_chain_arrays.items()})
        metric_arrays.update({name: _arr(values) for name, values in origin_cf_arrays.items()})
        metric_arrays.update({name: _arr(values) for name, values in combined_cf_arrays.items()})
        metric_arrays.update({name: _arr(values) for name, values in origin_sc_arrays.items()})
        metric_arrays.update({name: _arr(values) for name, values in origin_sc_fit_family_arrays.items()})
        metric_arrays.update({name: _arr(values) for name, values in combined_sc_arrays.items()})
        metric_arrays.update({
            name: _arr(values) for name, values in manufactured_lambda_arrays.items()
        })

        taylor = {}
        for name in ("a", "b", "X", "alpha", "a_minus_b_over_r2", "Lambda_over_r"):
            taylor[name] = even_taylor_fit(r, metric_arrays[name])

        identities = {
            "a_times_b_squared_minus_one_max_abs": _max_abs(a * b * b - 1.0),
            "b_minus_X_squared_max_abs": _max_abs(b - X * X),
            "a_minus_X_minus_four_max_abs": _max_abs(a - X ** (-4.0)),
            "Lambda_metric_reconstruction_max_abs": _max_abs(Lambda - lambda_from_metric),
            "areal_radius_minus_r_max_abs": _max_abs(areal - r),
        }

        finite_names = list(metric_arrays)
        all_finite = all(_finite(metric_arrays[name]) for name in finite_names)
        after = _state_digest(state)
        unchanged = before == after
        if not all_finite:
            problems.append(f"N={resolution}: non-finite profile in initialized-state audit")
        if not unchanged:
            problems.append(f"N={resolution}: input state mutated by observer")

        first = {
            "i": 0,
            "r": float(r[0]),
            "vendor_tf_ricci": float(vendor_tf[0]),
            "direct_tf_ricci": float(direct["tf_ricci"][0]),
            "signed_difference": float(tf_delta[0]),
            "vendor_l2_tf_ricci": float(vendor_l2_tf[0]),
            "direct_l2_tf_ricci": float(direct_l2_tf[0]),
            "l2_signed_difference": float(l2_delta[0]),
            "vendor_tf_ricci_over_r2": float(vendor_tf_over_r2[0]),
            "direct_tf_ricci_over_r2": float(direct_tf_over_r2[0]),
            "tf_ricci_difference_over_r2": float(tf_delta_over_r2[0]),
            "aa_lapse_hessian_tf": float(aa_lapse_hessian[0]),
            "aa_ricci_tf": float(aa_ricci_tf[0]),
            "aa_l3_vacuum": float(aa_l3_vacuum[0]),
            "aa_l3_matter": float(aa_l3_matter[0]),
            "aa_rhs": float(aa_rhs[0]),
            "aa_rhs_over_r2": float(aa_rhs_over_r2[0]),
            "Lambda_over_r": float(lambda_over_r[0]),
            "a_minus_b_over_r2": float(a_minus_b_over_r2[0]),
        }
        centre = []
        for i in range(min(CENTRE_CELLS, len(r))):
            centre.append({
                "i": i,
                "r": float(r[i]),
                "a_minus_b_over_r2": float(a_minus_b_over_r2[i]),
                "Lambda_over_r": float(lambda_over_r[i]),
                "a_prime_over_r": float(metric_arrays["a_prime_over_r"][i]),
                "b_prime_over_r": float(metric_arrays["b_prime_over_r"][i]),
                "X_prime_over_r": float(metric_arrays["X_prime_over_r"][i]),
                "alpha_prime_over_r": float(metric_arrays["alpha_prime_over_r"][i]),
                "vendor_tf_ricci": float(vendor_tf[i]),
                "direct_tf_ricci": float(direct["tf_ricci"][i]),
                "signed_difference": float(tf_delta[i]),
                "vendor_l2_tf_ricci": float(vendor_l2_tf[i]),
                "direct_l2_tf_ricci": float(direct_l2_tf[i]),
                "l2_signed_difference": float(l2_delta[i]),
                "vendor_tf_ricci_over_r2": float(vendor_tf_over_r2[i]),
                "direct_tf_ricci_over_r2": float(direct_tf_over_r2[i]),
                "tf_ricci_difference_over_r2": float(tf_delta_over_r2[i]),
                "aa_lapse_hessian_tf": float(aa_lapse_hessian[i]),
                "aa_ricci_tf": float(aa_ricci_tf[i]),
                "aa_l3_vacuum": float(aa_l3_vacuum[i]),
                "aa_l3_matter": float(aa_l3_matter[i]),
                "aa_rhs": float(aa_rhs[i]),
                "aa_rhs_over_r2": float(aa_rhs_over_r2[i]),
                "Ric_l_l": float(direct["Ric_l_l"][i]),
                "Ric_theta_theta": float(direct["Ric_theta_theta"][i]),
            })

        denom = max(
            _max_abs(vendor_tf[:CENTRE_CELLS]),
            _max_abs(direct["tf_ricci"][:CENTRE_CELLS]),
            np.finfo(float).tiny,
        )
        first_two_mask = np.zeros(r.size, dtype=bool)
        first_two_mask[:min(2, r.size)] = True
        near_centre_annulus = (r >= 2.0) & (r < 5.0)
        fixed_annulus = (r >= 5.0) & (r <= 10.0)
        region_fields = {
            "vendor_tf_ricci": vendor_tf,
            "direct_tf_ricci": _arr(direct["tf_ricci"]),
            "tf_ricci_difference": tf_delta,
            "vendor_tf_ricci_over_r2": vendor_tf_over_r2,
            "direct_tf_ricci_over_r2": direct_tf_over_r2,
            "tf_ricci_difference_over_r2": tf_delta_over_r2,
            "aa_lapse_hessian_tf": aa_lapse_hessian,
            "aa_ricci_tf": aa_ricci_tf,
            "aa_l3_vacuum": aa_l3_vacuum,
            "aa_l3_matter": aa_l3_matter,
            "aa_rhs": aa_rhs,
            "aa_rhs_over_r2": aa_rhs_over_r2,
        }
        report = {
            "resolution": resolution,
            "r_max": R_MAX,
            "dr": float(grid.dr),
            "evolution_steps": 0,
            "coordinate_time": float(state.t),
            "proper_time": float(state.tau),
            "radiation_enabled": True,
            "input_state_sha256_before": before,
            "input_state_sha256_after": after,
            "input_state_unchanged": unchanged,
            "all_profiles_finite": all_finite,
            "algebraic_and_connection_identities": identities,
            "first_cell": first,
            "first_two_cells": centre[:2],
            "first_six_cells": centre,
            "centre_vs_fixed_annulus": {
                "first_two_centre_cells": _regional_summary(
                    r, region_fields, first_two_mask
                ),
                "near_centre_annulus_r2_to5": _regional_summary(
                    r, region_fields, near_centre_annulus
                ),
                "fixed_annulus_r5_to10": _regional_summary(
                    r, region_fields, fixed_annulus
                ),
            },
            "signed_Aa_source_budget_closure": {
                "l2_operator_minus_lapse_plus_ricci_max_abs": _max_abs(
                    aa_l2_reconstruction_error
                ),
                "roundoff_bound": float(aa_l2_roundoff_bound),
                "bound_ratio": float(
                    _max_abs(aa_l2_reconstruction_error)
                    / max(aa_l2_roundoff_bound, np.finfo(float).tiny)
                ),
                "lapse_hessian_max_abs_first_six": _max_abs(
                    aa_lapse_hessian[:CENTRE_CELLS]
                ),
                "trace_free_ricci_max_abs_first_six": _max_abs(
                    aa_ricci_tf[:CENTRE_CELLS]
                ),
                "vacuum_L3_max_abs_first_six": _max_abs(
                    aa_l3_vacuum[:CENTRE_CELLS]
                ),
                "matter_L3_delta_max_abs_first_six": _max_abs(
                    aa_l3_matter[:CENTRE_CELLS]
                ),
                "full_Aa_rhs_over_r2_first_six": [
                    float(value) for value in aa_rhs_over_r2[:CENTRE_CELLS]
                ],
            },
            "ricci_comparison_first_six": {
                "vendor_tf_max_abs": _max_abs(vendor_tf[:CENTRE_CELLS]),
                "direct_tf_max_abs": _max_abs(direct["tf_ricci"][:CENTRE_CELLS]),
                "signed_difference_max_abs": _max_abs(tf_delta[:CENTRE_CELLS]),
                "relative_difference_to_larger_signal": float(
                    _max_abs(tf_delta[:CENTRE_CELLS]) / denom
                ),
                "vendor_l2_tf_max_abs": _max_abs(vendor_l2_tf[:CENTRE_CELLS]),
                "direct_l2_tf_max_abs": _max_abs(direct_l2_tf[:CENTRE_CELLS]),
                "l2_signed_difference_max_abs": _max_abs(l2_delta[:CENTRE_CELLS]),
            },
            "even_taylor_fits_first_six_cells": taylor,
            "frozen_source_origin_half_cell_counterfactual": origin_cf_summary,
            "self_consistent_origin_half_cell_counterfactual": origin_sc_summary,
            "self_consistent_origin_fit_family_counterfactuals": origin_sc_fit_family_summaries,
            "self_consistent_origin_and_lambda_chain_rule_counterfactual": combined_sc_summary,
            "lambda_chain_rule_counterfactual": lambda_chain_summary,
            "combined_origin_and_lambda_chain_rule_counterfactual": combined_cf_summary,
            "manufactured_lambda_parity_discriminator": manufactured_lambda_summary,
            "interpretation_limits": [
                "A successful observer workflow is not scientific admission.",
                "The direct warped-product reconstruction is independent of evolved Lambda and the vendor Ricci brackets, but shares the pinned radial grid derivative operator.",
                "Taylor fits summarize the first six cell-centred values; fit residuals are evidence, not a post-hoc pass/fail tolerance.",
                "A discrepancy localizes a representation/operator compatibility question but does not alone prove which equation is wrong.",
                "The source budget is an initialized-slice semi-discrete RHS decomposition, not an evolved-stage or trajectory result.",
                "The first-two-cell versus r=5..10 comparison is a coordinate-radial localization diagnostic, not a proper-volume or continuum-convergence proof.",
            ],
        }
        reports.append(report)

        prefix = f"N{resolution}_"
        for name, values in metric_arrays.items():
            profiles[prefix + name] = _arr(values)
        for name, result in taylor.items():
            profiles[prefix + name + "_taylor_residual_first6"] = _arr(result["signed_residuals"])

    status = "completed_diagnostic_only" if not problems else "observer_integrity_failure"
    result = {
        "schema": "v55_current_main_centre_proper_radial_ricci_v1",
        "status": status,
        "source_pins": provenance,
        "resolutions": reports,
        "errors": problems,
        "no_evolution": True,
        "production_or_vendor_files_modified": False,
        "scientific_admission": False,
    }

    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    npz_path = OUT / "profiles.npz"
    np.savez_compressed(npz_path, **profiles)
    csv_path = OUT / "profiles.csv"
    # Export the full signed profiles in a deterministic, flat table.
    all_names = list(reports and [
        "a", "b", "X", "alpha", "Lambda", "a_minus_b_over_r2", "Lambda_over_r",
        "areal_radius", "proper_factor", "vendor_tf_ricci", "direct_tf_ricci",
        "tf_ricci_difference", "vendor_tf_ricci_over_r2", "direct_tf_ricci_over_r2",
        "tf_ricci_difference_over_r2", "vendor_l2_tf_ricci", "direct_l2_tf_ricci",
        "l2_tf_ricci_difference", "aa_lapse_hessian_tf", "aa_ricci_tf",
        "aa_l3_vacuum", "aa_l3_matter", "aa_l2", "aa_l2_reconstruction_error",
        "aa_rhs", "aa_rhs_over_r2", "R_r", "R_rr", "f_r", "R_l", "R_ll",
        "Ric_l_l", "Ric_theta_theta", "a_prime_over_r", "b_prime_over_r",
        "X_prime_over_r", "alpha_prime_over_r", "lambda_metric_identity_residual",
        "areal_radius_minus_coordinate_radius", "initializer_source_final_iteration",
        "initializer_B_reconstructed", "initializer_B_from_X6", "origin_cf_B",
        "origin_cf_a", "origin_cf_b", "origin_cf_X", "origin_cf_Lambda",
        "origin_cf_areal_radius", "origin_cf_proper_factor",
        "origin_cf_vendor_tf_ricci", "origin_cf_direct_tf_ricci",
        "origin_cf_tf_ricci_difference", "origin_cf_vendor_tf_ricci_over_r2",
        "origin_cf_direct_tf_ricci_over_r2", "origin_cf_tf_ricci_difference_over_r2",
        "origin_sc_six_B", "origin_sc_six_a", "origin_sc_six_b", "origin_sc_six_X",
        "origin_sc_six_Lambda", "origin_sc_six_areal_radius", "origin_sc_six_proper_factor",
        "origin_sc_six_vendor_tf_ricci", "origin_sc_six_direct_tf_ricci",
        "origin_sc_six_tf_ricci_difference", "origin_sc_six_tf_ricci_difference_over_r2",
        "origin_sc_six_delta_over_r2", "origin_sc_six_Lambda_over_r",
        "origin_sc_six_source_recomputed_at_profile", "origin_sc_fixed_B", "origin_sc_fixed_a",
        "origin_sc_fixed_b", "origin_sc_fixed_X", "origin_sc_fixed_Lambda",
        "origin_sc_fixed_areal_radius", "origin_sc_fixed_proper_factor",
        "origin_sc_fixed_vendor_tf_ricci", "origin_sc_fixed_direct_tf_ricci",
        "origin_sc_fixed_tf_ricci_difference", "origin_sc_fixed_tf_ricci_difference_over_r2",
        "origin_sc_fixed_delta_over_r2", "origin_sc_fixed_Lambda_over_r",
        "origin_sc_fixed_source_recomputed_at_profile", "origin_sc_baseline_production_B",
        "origin_sc_frozen_source_cf_B", "origin_sc_initializer_source_baseline_final",
        "combined_origin_sc_lambda_cf_lambda_prime_discrete_D4",
        "combined_origin_sc_lambda_cf_lambda_prime_chain_rule",
        "combined_origin_sc_lambda_cf_lambda_prime_chain_minus_D4",
        "combined_origin_sc_lambda_cf_induced_tf_ricci_shift",
        "combined_origin_sc_lambda_cf_vendor_tf_ricci",
        "combined_origin_sc_lambda_cf_vendor_minus_direct_tf_ricci",
        "combined_origin_sc_lambda_cf_vendor_minus_direct_tf_ricci_over_r2",
        "lambda_prime_discrete_D4", "lambda_prime_chain_rule", "lambda_prime_chain_minus_D4",
        "lambda_chain_rule_induced_tf_ricci_shift", "lambda_chain_rule_vendor_tf_ricci",
        "lambda_chain_rule_vendor_minus_direct_tf_ricci",
        "lambda_chain_rule_vendor_minus_direct_tf_ricci_over_r2",
        "combined_origin_lambda_cf_lambda_prime_discrete_D4",
        "combined_origin_lambda_cf_lambda_prime_chain_rule",
        "combined_origin_lambda_cf_lambda_prime_chain_minus_D4",
        "combined_origin_lambda_cf_induced_tf_ricci_shift",
        "combined_origin_lambda_cf_vendor_tf_ricci",
        "combined_origin_lambda_cf_vendor_minus_direct_tf_ricci",
        "combined_origin_lambda_cf_vendor_minus_direct_tf_ricci_over_r2",
    ])
    for manufactured_tag in ("manufactured_q1", "manufactured_q2"):
        for manufactured_arm in (
            "continuum", "discrete", "off_constraint_odd", "negative_control_even"
        ):
            manufactured_prefix = f"{manufactured_tag}_{manufactured_arm}"
            all_names.extend([
                manufactured_prefix + "_a",
                manufactured_prefix + "_b",
                manufactured_prefix + "_X",
                manufactured_prefix + "_areal_radius",
                manufactured_prefix + "_proper_factor",
                manufactured_prefix + "_Lambda",
                manufactured_prefix + "_Lambda_over_r",
                manufactured_prefix + "_vendor_tf_ricci",
                manufactured_prefix + "_vendor_tf_ricci_over_r2",
                manufactured_prefix + "_signed_gap",
                manufactured_prefix + "_signed_gap_over_r2",
                manufactured_prefix + "_chain_rule_gap",
                manufactured_prefix + "_chain_rule_induced_tf_shift",
                manufactured_prefix + "_vendor_Rr_r",
                manufactured_prefix + "_vendor_R_scalar",
                manufactured_prefix + "_delta_Lambda_from_discrete",
                manufactured_prefix + "_D4_delta_Lambda_from_discrete",
                manufactured_prefix + "_observed_gap_change_from_discrete",
                manufactured_prefix + "_predicted_gap_change_derivative_only",
                manufactured_prefix + "_predicted_gap_change_algebraic_Lambda_ap",
                manufactured_prefix + "_predicted_gap_change_total",
                manufactured_prefix + "_response_closure_error",
            ])
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["resolution", "i", "r"] + all_names)
        for n in RESOLUTIONS:
            source = {key[len(f"N{n}_"):]: value for key, value in profiles.items()
                      if key.startswith(f"N{n}_") and not key.endswith("_taylor_residual_first6")}
            for i in range(source["a"].size):
                writer.writerow([n, i, float(i * (R_MAX / n) + 0.5 * (R_MAX / n))] +
                                [float(source[name][i]) for name in all_names])

    checksum_paths = (report_path, npz_path, csv_path)
    (OUT / "SHA256SUMS.txt").write_text(
        "\n".join(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}"
                  for path in checksum_paths) + "\n",
        encoding="utf-8",
    )
    (OUT / "SOURCE_COMMIT.txt").write_text(SOURCE_COMMIT + "\n", encoding="utf-8")
    (OUT / "VENDOR_COMMIT.txt").write_text(VENDOR_COMMIT + "\n", encoding="utf-8")

    print(json.dumps({
        "status": status,
        "head": provenance["audit_head"],
        "resolutions": [
            {
                "N": item["resolution"],
                "centre_vendor_tf": item["first_cell"]["vendor_tf_ricci"],
                "centre_direct_tf": item["first_cell"]["direct_tf_ricci"],
                "centre_difference": item["first_cell"]["signed_difference"],
                "first_six_relative_difference": item["ricci_comparison_first_six"]["relative_difference_to_larger_signal"],
                "annulus_r5_to10_difference": (
                    item["centre_vs_fixed_annulus"]["fixed_annulus_r5_to10"]["profiles"]["tf_ricci_difference"]["max_abs"]
                    if item["centre_vs_fixed_annulus"]["fixed_annulus_r5_to10"]["cells"] else None
                ),
                "first_cell_aa_rhs_over_r2": item["first_cell"]["aa_rhs_over_r2"],
                "state_unchanged": item["input_state_unchanged"],
            }
            for item in reports
        ],
        "errors": problems,
    }, indent=2, allow_nan=False), flush=True)
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())