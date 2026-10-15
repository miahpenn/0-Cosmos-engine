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
        metric_arrays.update({name: _arr(values) for name, values in origin_cf_arrays.items()})

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
