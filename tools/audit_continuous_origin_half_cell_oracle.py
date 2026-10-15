"""Research-only source-specific continuous origin half-cell integral oracle."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.scalar_system import cosmos_potential
from engine.v55_matter import RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT
from tools.audit_current_main_centre_proper_radial_ricci import (
    _finite, _max_abs, _state_digest, _verify_source_pins,
)

PI = float(np.pi)
PHI0 = 0.179055
PI0 = 0.00341913
SOURCE_COMMIT = "285a2ee0c1a3b97cddeb4519d8e749109b115aab"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SUBSTEPS = (256, 512, 1024, 2048, 4096, 8192)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def continuous_coefficients(r, amplitude=0.01, width=7.0, D_amplitude=1.0e-10):
    """Return A,q for the initializer source s=A-qB under its analytic fields."""
    r = np.asarray(r, dtype=float)
    if r.ndim != 1 or not np.all(np.isfinite(r)) or np.any(r < 0):
        raise ValueError("r must be a finite nonnegative one-dimensional array")
    if width <= 0 or not np.isfinite(width):
        raise ValueError("width must be finite and positive")
    if not np.isfinite(amplitude) or not np.isfinite(D_amplitude):
        raise ValueError("field amplitudes must be finite")
    rho_r = RHO_R_PRESENT
    rho_phi = (0.5 * PI0**2 + float(cosmos_potential(np.array([PHI0]))[0])) / (8.0 * PI)
    rho_matter = (RHO_DM_PRESENT + RHO_B_PRESENT + rho_r) / (8.0 * PI)
    rho0 = 0.5 * PI0**2 + float(cosmos_potential(np.array([PHI0]))[0])
    H0 = float(np.sqrt(max((rho0 + RHO_DM_PRESENT + RHO_B_PRESENT + rho_r) / 3.0, 0.0)))
    Kt = -H0
    S = amplitude * np.exp(-((r / width) ** 2))
    D = D_amplitude * np.exp(-((r / width) ** 2))
    Sp = -2.0 * r * S / width**2
    Dp = -2.0 * r * D / width**2
    j = -(D * Dp)
    Kr = Kt + 4.0 * PI * r * j
    A = (
        1.0 - 8.0 * PI * r**2 * (0.5 * S**2 + rho_phi + rho_matter)
        + 2.0 * r**2 * Kr * Kt + r**2 * Kt**2
    )
    q = 4.0 * PI * r**2 * (Sp**2 + Dp**2)
    if not (_finite(A) and _finite(q)):
        raise FloatingPointError("non-finite continuous-source coefficients")
    return A, q




def _background_parameters():
    """Return Kt and the homogeneous coefficient for the declared uncharged setup."""
    rho_r = RHO_R_PRESENT
    potential = float(cosmos_potential(np.array([PHI0]))[0])
    rho_phi = (0.5 * PI0**2 + potential) / (8.0 * PI)
    rho_matter = (RHO_DM_PRESENT + RHO_B_PRESENT + rho_r) / (8.0 * PI)
    rho0 = 0.5 * PI0**2 + potential + RHO_DM_PRESENT + RHO_B_PRESENT + rho_r
    H0 = math.sqrt(max(rho0 / 3.0, 0.0))
    Kt = -H0
    c_bg = 3.0 * Kt**2 - 8.0 * PI * (rho_phi + rho_matter)
    return Kt, c_bg


def _gaussian_even_moment_integral(r: float, power: int, a: float) -> float:
    """Evaluate integral_0^r t**power exp(-a*t**2) dt stably for power 2 or 4."""
    if power not in (2, 4) or not np.isfinite(r) or r < 0.0 or a <= 0.0:
        raise ValueError("expected power 2 or 4, r >= 0, and a > 0")
    if r == 0.0:
        return 0.0
    z = a * r * r
    if z <= 0.5:
        term = r ** (power + 1) / (power + 1)
        total = term
        for n in range(1, 100):
            term *= (-z / n) * ((power + 2 * n - 1) / (power + 2 * n + 1))
            total += term
            if abs(term) <= np.finfo(float).eps * max(abs(total), np.finfo(float).tiny):
                break
        return float(total)
    x = math.sqrt(a) * r
    erf_term = math.sqrt(math.pi) * math.erf(x)
    if power == 2:
        return float(erf_term / (4.0 * a**1.5) - r * math.exp(-z) / (2.0 * a))
    return float(
        3.0 * erf_term / (8.0 * a**2.5)
        - math.exp(-z) * (r**3 / (2.0 * a) + 3.0 * r / (4.0 * a**2))
    )


def analytic_A_integral(r, amplitude=0.01, width=7.0, D_amplitude=1.0e-10):
    """Closed-form primitive integral_0^r A(t) dt for the declared continuous source."""
    r = float(r)
    if not np.isfinite(r) or r < 0.0:
        raise ValueError("r must be finite and nonnegative")
    if not np.isfinite(width) or width <= 0.0:
        raise ValueError("width must be finite and positive")
    if not np.isfinite(amplitude) or not np.isfinite(D_amplitude):
        raise ValueError("field amplitudes must be finite")
    Kt, c_bg = _background_parameters()
    a = 2.0 / width**2
    F2 = _gaussian_even_moment_integral(r, 2, a)
    F4 = _gaussian_even_moment_integral(r, 4, a)
    return float(
        r + c_bg * r**3 / 3.0
        - 4.0 * PI * amplitude**2 * F2
        + 16.0 * PI * Kt * D_amplitude**2 / width**2 * F4
    )


def analytic_half_cell_envelope(
    r_half, amplitude=0.01, width=7.0, D_amplitude=1.0e-10
):
    """Analytic integral/error enclosure for the declared uncharged Gaussian source.

    With s=A-qB, I=rB, and I'=A-qB, let I_A=integral(A) and
    Q=integral(qB). Positivity and the bounds B<=A_max and q<=C*r^4
    imply 0<=Q<=A_max*C*r_half^5/5, hence
    I_A-Q_max<=I<=I_A. The midpoint defect has a corresponding enclosure.
    """
    r_half = float(r_half)
    if not np.isfinite(r_half) or r_half <= 0.0:
        raise ValueError("r_half must be finite and positive")
    if not np.isfinite(width) or width <= 0.0:
        raise ValueError("width must be finite and positive")
    if not np.isfinite(amplitude) or not np.isfinite(D_amplitude):
        raise ValueError("field amplitudes must be finite")
    Kt, c_bg = _background_parameters()
    if Kt >= 0.0:
        raise ValueError("the analytic enclosure assumes the declared Kt < 0 branch")
    a = 2.0 / width**2
    A_upper = 1.0 + max(c_bg, 0.0) * r_half**2
    A_lower = (
        1.0 + min(c_bg, 0.0) * r_half**2
        - 4.0 * PI * amplitude**2 * r_half**2
        - 16.0 * PI * abs(Kt) * D_amplitude**2 * r_half**4 / width**2
    )
    if not (np.isfinite(A_lower) and np.isfinite(A_upper) and A_lower > 0.0):
        raise ValueError("source positivity could not be established on this half-cell")
    q_coefficient = 16.0 * PI * (amplitude**2 + D_amplitude**2) / width**4
    q_absorption_bound = A_upper * q_coefficient * r_half**5 / 5.0
    I_A = analytic_A_integral(
        r_half, amplitude=amplitude, width=width, D_amplitude=D_amplitude
    )
    A_end, q_end = continuous_coefficients(
        np.array([r_half]), amplitude=amplitude, width=width, D_amplitude=D_amplitude
    )
    A_end = float(A_end[0])
    q_end = float(q_end[0])
    integral_lower = I_A - q_absorption_bound
    integral_upper = I_A
    delta_A = I_A - r_half * A_end
    defect_lower = delta_A - q_absorption_bound
    defect_upper = delta_A + r_half * q_end * A_upper
    return {
        "A_integral_closed_form": float(I_A),
        "A_lower_bound_on_half_cell": float(A_lower),
        "A_upper_bound_on_half_cell": float(A_upper),
        "q_coefficient_C": float(q_coefficient),
        "Q_absorption_upper_bound": float(q_absorption_bound),
        "continuous_integral_lower_bound": float(integral_lower),
        "continuous_integral_upper_bound": float(integral_upper),
        "midpoint_defect_A_only": float(delta_A),
        "continuous_midpoint_defect_lower_bound": float(defect_lower),
        "continuous_midpoint_defect_upper_bound": float(defect_upper),
        "assumptions": [
            "Default uncharged source with the declared analytic Gaussian S and D fields.",
            "Constant phi, Pi, and Kt; radiation enabled; no charged dust or Maxwell source.",
            "Kt < 0 and A(r) > 0 on [0, r_half].",
            "Analytic inequalities are evaluated in float64; no directed-rounding interval arithmetic is used.",
        ],
    }


def integrate_half_cell(r_half, n, coeff_fn):
    """Implicit trapezoid solve of I'=A-q I/r with regular I(0)=0, B(0)=1."""
    if not np.isfinite(r_half) or r_half <= 0:
        raise ValueError("r_half must be finite and positive")
    if isinstance(n, bool) or int(n) != n or n < 2:
        raise ValueError("n must be an integer of at least two")
    r = np.linspace(0.0, float(r_half), int(n) + 1)
    h = float(r_half) / int(n)
    A, q = coeff_fn(r)
    A = np.asarray(A, dtype=float)
    q = np.asarray(q, dtype=float)
    if A.shape != r.shape or q.shape != r.shape:
        raise ValueError("coefficient arrays must match the integration grid")
    if not (_finite(A) and _finite(q)) or abs(A[0] - 1.0) > 1.0e-12 or abs(q[0]) > 1.0e-12:
        raise ValueError("non-finite coefficients or invalid regular-origin limits")
    if np.any(q < -1.0e-15):
        raise ValueError("q must be nonnegative")
    I = np.zeros_like(r)
    B = np.empty_like(r)
    s = np.empty_like(r)
    B[0] = 1.0
    s[0] = 1.0
    for i in range(1, len(r)):
        I[i] = (I[i - 1] + 0.5 * h * (s[i - 1] + A[i])) / (1.0 + 0.5 * h * q[i] / r[i])
        B[i] = I[i] / r[i]
        s[i] = A[i] - q[i] * B[i]
    if not (_finite(I) and _finite(B) and _finite(s)):
        raise FloatingPointError("non-finite continuous integration output")
    return {
        "moment": float(I[-1]),
        "B_at_half_cell": float(B[-1]),
        "source_at_half_cell": float(s[-1]),
        "substeps": int(n),
        "step_size": h,
    }


def continuous_reference(r_half, n_levels=SUBSTEPS):
    coeff_fn = continuous_coefficients
    raw = {n: integrate_half_cell(r_half, n, coeff_fn) for n in n_levels}
    rich = []
    for coarse, fine in zip(n_levels, n_levels[1:]):
        if fine != 2 * coarse:
            continue
        c, f = raw[coarse]["moment"], raw[fine]["moment"]
        rich.append({
            "coarse_substeps": coarse, "fine_substeps": fine,
            "coarse_moment": c, "fine_moment": f,
            "richardson_moment": f + (f - c) / 3.0,
            "raw_difference": f - c,
        })
    if len(rich) < 2:
        raise ValueError("need at least two consecutive doubling pairs")
    ref = rich[-1]["richardson_moment"]
    change = abs(ref - rich[-2]["richardson_moment"])
    A_end, q_end = coeff_fn(np.array([r_half]))
    B_ref = ref / r_half
    source_end = float(A_end[0] - q_end[0] * B_ref)
    diffs = [abs(x["raw_difference"]) for x in rich]
    orders = [
        float(np.log(a / b) / np.log(2.0)) if a > 0 and b > 0 else None
        for a, b in zip(diffs, diffs[1:])
    ]
    return {
        "reference_moment": float(ref),
        "reference_B": float(B_ref),
        "reference_endpoint_source": source_end,
        "estimated_numeric_uncertainty": float(max(change / 15.0, 8.0 * np.finfo(float).eps * abs(ref))),
        "estimated_roundoff_floor": float(8.0 * np.finfo(float).eps * abs(ref)),
        "change_between_final_extrapolants": float(change),
        "raw_levels": [raw[n] for n in n_levels],
        "richardson_levels": rich,
        "observed_raw_orders": orders,
        "method": "implicit trapezoid plus second-order Richardson extrapolation",
        "limit": "convergence estimate only, not a rigorous analytic remainder bound",
    }


def sampled_initializer_source(state):
    """Replay the six point-sampled source iterations used by build_initial_data."""
    g = state.grid
    r = np.asarray(g.centers, dtype=float)
    dr = float(g.dr)
    f = state.scalars
    S, D = np.asarray(f.S), np.asarray(f.D)
    PS, PD = np.asarray(f.PS), np.asarray(f.PD)
    phi, Pi = np.asarray(f.phi), np.asarray(f.Pi)
    Sp = g.cell_derivative_fourth(S, parity=1)
    Dp = g.cell_derivative_fourth(D, parity=1)
    Kt = np.asarray(state.geometry.K) / 3.0 - np.asarray(state.geometry.Aa) / 2.0
    Ktp = g.cell_derivative_fourth(Kt, parity=1)
    rho_matter = (RHO_DM_PRESENT + RHO_B_PRESENT + RHO_R_PRESENT) / (8.0 * PI)
    rho_phi = (0.5 * Pi**2 + cosmos_potential(phi)) / (8.0 * PI)
    j = -(PD * Dp)
    B = np.ones_like(r)
    src_last = None
    for _ in range(6):
        rhoS = 0.5 * (PS**2 + B * Sp**2) + 0.5 * S**2
        rhoD = 0.5 * (PD**2 + B * Dp**2) - 0.5 * D**2
        rho = rhoS + rhoD + rho_phi + rho_matter
        Kr = Kt + r * Ktp + 4.0 * PI * r * j
        src_last = 1.0 - 8.0 * PI * r**2 * rho + 2.0 * r**2 * Kr * Kt + r**2 * Kt**2
        integ = np.cumsum(np.r_[
            0.5 * dr * src_last[0],
            0.5 * dr * (src_last[1:] + src_last[:-1]),
        ])
        B = integ / r
    B_actual = np.asarray(state.geometry.X) ** 6
    return {
        "r": r,
        "source_last_iteration": np.asarray(src_last),
        "B_replayed": np.asarray(B),
        "B_actual": B_actual,
        "B_replay_max_abs_error": _max_abs(B - B_actual),
        "production_point_moment": float(r[0] * src_last[0]),
    }


def fit_even_polynomial_integral(r, source, r_half, window, degree, anchor=False):
    """Fit a polynomial in z=(r/window)^2 and integrate it analytically."""
    r = np.asarray(r, dtype=float)
    source = np.asarray(source, dtype=float)
    if r.ndim != 1 or r.shape != source.shape:
        raise ValueError("r and source must be matching one-dimensional arrays")
    if degree < 1 or window <= 0:
        raise ValueError("degree and window must be positive")
    mask = r <= window + 16.0 * np.finfo(float).eps * max(1.0, window)
    rs, ys = r[mask], source[mask]
    if rs.size < degree + 2:
        raise ValueError("fit must be overdetermined by at least one sample")
    z = (rs / window) ** 2
    if anchor:
        design = np.column_stack([z**k for k in range(1, degree + 1)])
        shifted_y = ys - 1.0
        tail, _, rank, _ = np.linalg.lstsq(design, shifted_y, rcond=None)
        c = np.r_[1.0, tail]
        loo_intercept = 1.0
    else:
        design = np.column_stack([z**k for k in range(degree + 1)])
        shifted_y = ys
        c, _, rank, _ = np.linalg.lstsq(design, shifted_y, rcond=None)
        loo_intercept = 0.0
    if rank < design.shape[1]:
        raise ValueError("fit design is rank deficient")
    pred = sum(c[k] * z**k for k in range(degree + 1))
    moment = sum(
        c[k] * r_half**(2 * k + 1) / ((2 * k + 1) * window**(2 * k))
        for k in range(degree + 1)
    )
    residual = ys - pred

    # Leave-one-out prediction: refit on all but one source sample and predict
    # the omitted sample. Report the worst fold condition number as well as error.
    loo_residuals = []
    loo_condition_numbers = []
    loo_training_residual_dof = []
    for holdout in range(rs.size):
        keep = np.ones(rs.size, dtype=bool)
        keep[holdout] = False
        train_design = design[keep]
        train_y = shifted_y[keep]
        if train_design.shape[0] < design.shape[1]:
            raise ValueError("leave-one-out training set is underdetermined")
        loo_coeff, _, loo_rank, _ = np.linalg.lstsq(train_design, train_y, rcond=None)
        if loo_rank < design.shape[1]:
            raise ValueError("leave-one-out fit design is rank deficient")
        held_prediction = loo_intercept + float(design[holdout] @ loo_coeff)
        loo_residuals.append(float(ys[holdout] - held_prediction))
        loo_condition_numbers.append(float(np.linalg.cond(train_design)))
        loo_training_residual_dof.append(int(train_design.shape[0] - design.shape[1]))
    loo_residuals = np.asarray(loo_residuals, dtype=float)
    if not (_finite(c) and _finite(residual) and _finite(loo_residuals) and np.isfinite(moment)):
        raise FloatingPointError("non-finite fit or leave-one-out output")
    return {
        "degree": int(degree), "window": float(window), "anchor_origin": bool(anchor),
        "sample_count": int(rs.size), "coefficient_count": int(degree + 1),
        "condition_number": float(np.linalg.cond(design)),
        "fit_coefficients_z": [float(v) for v in c],
        "max_abs_sample_residual": _max_abs(residual),
        "rms_sample_residual": float(np.sqrt(np.mean(residual**2))),
        "loo_residuals": [float(v) for v in loo_residuals],
        "loo_max_abs_prediction_error": _max_abs(loo_residuals),
        "loo_rms_prediction_error": float(np.sqrt(np.mean(loo_residuals**2))),
        "loo_max_condition_number": max(loo_condition_numbers),
        "loo_min_training_residual_dof": min(loo_training_residual_dof),
        "fitted_moment": float(moment),
    }


def run_resolution(n, r_max=80.0):
    from engine.production_kernel import V55ProductionKernel
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=n, r_max=r_max, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True,
    )
    digest_before = _state_digest(state)
    samples = sampled_initializer_source(state)
    if samples["B_replay_max_abs_error"] > 5.0e-12:
        raise RuntimeError(f"production source replay mismatch at N={n}: {samples['B_replay_max_abs_error']}")
    r_half = 0.5 * float(state.grid.dr)
    ref = continuous_reference(r_half)
    envelope = analytic_half_cell_envelope(r_half)
    I_cont = ref["reference_moment"]
    enclosure_pad = 64.0 * np.finfo(float).eps * max(1.0, abs(envelope["A_integral_closed_form"]))
    if (
        I_cont < envelope["continuous_integral_lower_bound"] - enclosure_pad
        or I_cont > envelope["continuous_integral_upper_bound"] + enclosure_pad
    ):
        raise RuntimeError(f"Richardson estimate lies outside analytic integral enclosure at N={n}")
    I_mid_cont = r_half * ref["reference_endpoint_source"]
    I_prod = samples["production_point_moment"]
    q_delta = I_cont - I_mid_cont
    source_delta = I_mid_cont - I_prod
    total_delta = I_cont - I_prod
    closure = q_delta + source_delta - total_delta
    if abs(closure) > 64.0 * np.finfo(float).eps * max(abs(I_cont), 1.0e-300):
        raise RuntimeError(f"moment decomposition failed closure at N={n}: {closure}")

    specs = [
        (3, 2.5, False, "degree3_window2.5_unanchored"),
        (3, 2.5, True, "degree3_window2.5_origin_anchored"),
        (3, 3.0, False, "degree3_window3.0_unanchored"),
        (4, 3.0, False, "degree4_window3.0_unanchored"),
        (5, 3.5, False, "degree5_window3.5_unanchored"),
    ]
    fits = []
    for degree, window, anchor, name in specs:
        try:
            fit = fit_even_polynomial_integral(
                samples["r"], samples["source_last_iteration"], r_half,
                window, degree, anchor,
            )
            fit["name"] = name
            fit["fit_minus_continuous_reference"] = fit["fitted_moment"] - I_cont
            fit["fit_error_abs"] = abs(fit["fit_minus_continuous_reference"])
            fit["fit_error_rel_to_integral"] = fit["fit_error_abs"] / abs(I_cont)
            fit["fit_minus_production_point_moment"] = fit["fitted_moment"] - I_prod
            err_lo = fit["fitted_moment"] - envelope["A_integral_closed_form"]
            err_hi = err_lo + envelope["Q_absorption_upper_bound"]
            fit["analytic_signed_error_lower_bound"] = float(err_lo)
            fit["analytic_signed_error_upper_bound"] = float(err_hi)
            fit["analytic_absolute_error_bound"] = float(max(abs(err_lo), abs(err_hi)))
            estimate_error = fit["fitted_moment"] - I_cont
            fit["richardson_error_inside_analytic_enclosure"] = bool(
                estimate_error >= err_lo - enclosure_pad
                and estimate_error <= err_hi + enclosure_pad
            )
            if not fit["richardson_error_inside_analytic_enclosure"]:
                raise RuntimeError(f"fit error lies outside analytic enclosure for {name} at N={n}")
            fits.append(fit)
        except ValueError as exc:
            fits.append({
                "name": name, "degree": degree, "window": window,
                "anchor_origin": anchor, "status": "not_eligible", "reason": str(exc),
            })
    digest_after = _state_digest(state)
    if digest_after != digest_before:
        raise RuntimeError("observer mutated initialized state")
    return {
        "resolution": int(n), "r_max": float(r_max), "dr": float(state.grid.dr),
        "r_half": r_half, "continuous_reference": ref,
        "analytic_enclosure": {
            **envelope,
            "richardson_reference_inside_enclosure": True,
            "floating_point_enclosure_padding": float(enclosure_pad),
        },
        "production_discrete_source": {
            "source_last_iteration_cell0": float(samples["source_last_iteration"][0]),
            "production_point_moment": I_prod,
            "B_replay_max_abs_error": float(samples["B_replay_max_abs_error"]),
            "production_B0": float(samples["B_actual"][0]),
            "replayed_B0": float(samples["B_replayed"][0]),
        },
        "moment_decomposition": {
            "continuous_half_cell_integral_reference": I_cont,
            "continuous_midpoint_moment": I_mid_cont,
            "production_point_sampled_moment": I_prod,
            "continuous_quadrature_delta": q_delta,
            "continuum_vs_discrete_sampled_source_delta": source_delta,
            "total_continuum_minus_production_delta": total_delta,
            "decomposition_closure_error": closure,
            "continuous_B0": ref["reference_B"],
            "production_B0": float(samples["B_actual"][0]),
            "estimated_reference_uncertainty": ref["estimated_numeric_uncertainty"],
        },
        "fit_family": fits,
        "state_provenance": {
            "initialized_state_sha256_before": digest_before,
            "initialized_state_sha256_after": digest_after,
            "observer_left_state_unchanged": digest_before == digest_after,
            "evolution_steps": 0,
        },
    }


def main():
    provenance = _verify_source_pins()
    if provenance.get("source_commit") != SOURCE_COMMIT or provenance.get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError("source/vendor commit pins differ from declared audit pins")
    if not provenance.get("engine_vendor_dependency_diff_empty"):
        raise RuntimeError("engine/vendor/dependency diff is not empty")
    if not all(item.get("matches") for item in provenance.get("source_blobs", {}).values()):
        raise RuntimeError("one or more source blob pins failed")
    resolutions = tuple(int(x.strip()) for x in os.environ.get("ORACLE_RESOLUTIONS", "160,240,320").split(",") if x.strip())
    if len(resolutions) < 2 or any(n < 32 for n in resolutions):
        raise ValueError("at least two resolutions of 32 or greater are required")
    r_max = float(os.environ.get("ORACLE_R_MAX", "80.0"))
    if not np.isfinite(r_max) or r_max <= 10.0:
        raise ValueError("r_max must be finite and greater than 10")
    output = Path(os.environ.get("ORACLE_OUT", "results/continuous-origin-half-cell-oracle-2026-10-15"))
    rows = [run_resolution(n, r_max=r_max) for n in resolutions]
    report = {
        "schema": "v55_continuous_origin_half_cell_oracle_v2",
        "status": "completed_research_only_continuous_integral_enclosure",
        "scientific_admission": False, "production_or_vendor_files_modified": False,
        "source_pins": provenance,
        "declared_pins": {"engine_source": SOURCE_COMMIT, "vendor": VENDOR_COMMIT},
        "configuration": {
            "resolutions": list(resolutions), "r_max": r_max, "amplitude": 0.01,
            "width": 7.0, "D_amplitude": 1.0e-10, "radiation": True,
            "charged_dust": False, "radial_maxwell": False, "evolution_steps": 0,
        },
        "resolutions": rows,
        "all_observers_left_state_unchanged": all(x["state_provenance"]["observer_left_state_unchanged"] for x in rows),
        "all_B_replays_match_production": all(x["production_discrete_source"]["B_replay_max_abs_error"] <= 5.0e-12 for x in rows),
        "all_profiles_finite": all(
            _finite([
                x["continuous_reference"]["reference_moment"],
                x["continuous_reference"]["reference_B"],
                x["continuous_reference"]["reference_endpoint_source"],
                x["moment_decomposition"]["continuous_quadrature_delta"],
                x["moment_decomposition"]["continuum_vs_discrete_sampled_source_delta"],
                x["moment_decomposition"]["total_continuum_minus_production_delta"],
                x["moment_decomposition"]["decomposition_closure_error"],
                x["moment_decomposition"]["estimated_reference_uncertainty"],
                x["analytic_enclosure"]["A_integral_closed_form"],
                x["analytic_enclosure"]["A_lower_bound_on_half_cell"],
                x["analytic_enclosure"]["A_upper_bound_on_half_cell"],
                x["analytic_enclosure"]["q_coefficient_C"],
                x["analytic_enclosure"]["Q_absorption_upper_bound"],
                x["analytic_enclosure"]["continuous_integral_lower_bound"],
                x["analytic_enclosure"]["continuous_integral_upper_bound"],
                x["analytic_enclosure"]["midpoint_defect_A_only"],
                x["analytic_enclosure"]["continuous_midpoint_defect_lower_bound"],
                x["analytic_enclosure"]["continuous_midpoint_defect_upper_bound"],
            ])
            and all(
                fit.get("status") == "not_eligible"
                or _finite([
                    fit["fitted_moment"], fit["fit_minus_continuous_reference"],
                    fit["fit_error_abs"], fit["fit_error_rel_to_integral"],
                    fit["condition_number"], fit["max_abs_sample_residual"],
                    fit["rms_sample_residual"], fit["loo_max_abs_prediction_error"],
                    fit["loo_rms_prediction_error"], fit["loo_max_condition_number"],
                    fit["analytic_signed_error_lower_bound"],
                    fit["analytic_signed_error_upper_bound"],
                    fit["analytic_absolute_error_bound"],
                ])
                and fit["richardson_error_inside_analytic_enclosure"]
                for fit in x["fit_family"]
            )
            for x in rows
        ),
        "interpretation_limits": [
            "The continuous source assumes the declared analytic Gaussian S/D fields and constant phi, Pi, K_t.",
            "A closed-form integral of A and a derived analytic inequality enclose the full source integral and any fit estimator error under the declared uncharged model.",
            "The analytic inequality is evaluated in float64, not directed-rounding interval arithmetic; it is conditional on the declared continuous initialization.",
            "This is the default uncharged initialization and does not include optional charged-dust or Maxwell source terms.",
            "Richardson extrapolation provides a convergence estimate, not a rigorous analytic remainder bound.",
            "Fit-to-continuum errors are conditional on this continuous extension and do not by themselves establish a production correction.",
            "No production or vendor code changed; the audit takes zero evolution steps.",
        ],
    }
    if not report["all_observers_left_state_unchanged"]:
        raise RuntimeError("one or more observers changed their initialized state")
    if not report["all_B_replays_match_production"]:
        raise RuntimeError("one or more source replays failed against production X^6")
    if not report["all_profiles_finite"]:
        raise FloatingPointError("one or more report profiles or fit candidates are non-finite")
    output.mkdir(parents=True, exist_ok=True)
    paths = {}
    paths["report.json"] = output / "report.json"
    paths["report.json"].write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    paths["moments.csv"] = output / "moments.csv"
    with paths["moments.csv"].open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "N", "dr", "r_half", "I_cont_Richardson", "I_prod", "I_A_closed_form",
            "I_true_lower", "I_true_upper", "Q_absorption_bound",
            "midpoint_defect_lower", "midpoint_defect_upper",
            "quadrature_delta_Richardson", "sampled_source_delta", "total_delta",
            "oracle_uncertainty",
        ])
        for x in rows:
            m = x["moment_decomposition"]
            e = x["analytic_enclosure"]
            w.writerow([
                x["resolution"], x["dr"], x["r_half"],
                m["continuous_half_cell_integral_reference"],
                m["production_point_sampled_moment"],
                e["A_integral_closed_form"],
                e["continuous_integral_lower_bound"],
                e["continuous_integral_upper_bound"],
                e["Q_absorption_upper_bound"],
                e["continuous_midpoint_defect_lower_bound"],
                e["continuous_midpoint_defect_upper_bound"],
                m["continuous_quadrature_delta"],
                m["continuum_vs_discrete_sampled_source_delta"],
                m["total_continuum_minus_production_delta"],
                m["estimated_reference_uncertainty"],
            ])
    paths["fit_family.csv"] = output / "fit_family.csv"
    with paths["fit_family.csv"].open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "N", "name", "status", "degree", "window", "anchor_origin",
            "sample_count", "fitted_moment", "fit_minus_Richardson_reference",
            "fit_error_abs_vs_Richardson", "analytic_signed_error_lower",
            "analytic_signed_error_upper", "analytic_absolute_error_bound",
            "max_abs_sample_residual", "condition_number", "LOO_max_abs_prediction_error",
            "LOO_rms_prediction_error", "LOO_max_condition_number",
            "LOO_min_training_residual_dof",
        ])
        for x in rows:
            for fit in x["fit_family"]:
                w.writerow([
                    x["resolution"], fit["name"], fit.get("status", "ok"),
                    fit["degree"], fit["window"], fit["anchor_origin"],
                    fit.get("sample_count"), fit.get("fitted_moment"),
                    fit.get("fit_minus_continuous_reference"), fit.get("fit_error_abs"),
                    fit.get("analytic_signed_error_lower_bound"),
                    fit.get("analytic_signed_error_upper_bound"),
                    fit.get("analytic_absolute_error_bound"),
                    fit.get("max_abs_sample_residual"), fit.get("condition_number"),
                    fit.get("loo_max_abs_prediction_error"),
                    fit.get("loo_rms_prediction_error"),
                    fit.get("loo_max_condition_number"),
                    fit.get("loo_min_training_residual_dof"),
                ])
    (output / "SHA256SUMS.txt").write_text("\n".join(f"{_sha256(p)}  {p.name}" for p in paths.values()) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"], "audit_head": provenance.get("audit_head"),
        "results": [{
            "N": x["resolution"],
            "I_cont": x["moment_decomposition"]["continuous_half_cell_integral_reference"],
            "I_prod": x["moment_decomposition"]["production_point_sampled_moment"],
            "delta_total": x["moment_decomposition"]["total_continuum_minus_production_delta"],
            "delta_quadrature": x["moment_decomposition"]["continuous_quadrature_delta"],
            "delta_sampled_source": x["moment_decomposition"]["continuum_vs_discrete_sampled_source_delta"],
            "B_replay_error": x["production_discrete_source"]["B_replay_max_abs_error"],
            "uncertainty_estimate": x["moment_decomposition"]["estimated_reference_uncertainty"],
            "analytic_I_A": x["analytic_enclosure"]["A_integral_closed_form"],
            "analytic_integral_lower": x["analytic_enclosure"]["continuous_integral_lower_bound"],
            "analytic_integral_upper": x["analytic_enclosure"]["continuous_integral_upper_bound"],
            "analytic_fit_errors_enclosed": all(
                fit.get("status") == "not_eligible" or fit["richardson_error_inside_analytic_enclosure"]
                for fit in x["fit_family"]
            ),
        } for x in rows],
        "state_unchanged": report["all_observers_left_state_unchanged"],
        "finite": report["all_profiles_finite"], "errors": [],
    }, indent=2, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
