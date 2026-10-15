"""Research-only source-specific continuous origin half-cell integral oracle."""
from __future__ import annotations

import csv
import hashlib
import json
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
        "estimated_numeric_uncertainty": float(change / 15.0),
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
        tail, _, _, _ = np.linalg.lstsq(design, ys - 1.0, rcond=None)
        c = np.r_[1.0, tail]
    else:
        design = np.column_stack([z**k for k in range(degree + 1)])
        c, _, _, _ = np.linalg.lstsq(design, ys, rcond=None)
    pred = sum(c[k] * z**k for k in range(degree + 1))
    moment = sum(
        c[k] * r_half**(2 * k + 1) / ((2 * k + 1) * window**(2 * k))
        for k in range(degree + 1)
    )
    residual = ys - pred
    if not (_finite(c) and _finite(residual) and np.isfinite(moment)):
        raise FloatingPointError("non-finite fit output")
    return {
        "degree": int(degree), "window": float(window), "anchor_origin": bool(anchor),
        "sample_count": int(rs.size), "coefficient_count": int(degree + 1),
        "condition_number": float(np.linalg.cond(design)),
        "fit_coefficients_z": [float(v) for v in c],
        "max_abs_sample_residual": _max_abs(residual),
        "rms_sample_residual": float(np.sqrt(np.mean(residual**2))),
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
    I_cont = ref["reference_moment"]
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
    resolutions = tuple(int(x.strip()) for x in os.environ.get("ORACLE_RESOLUTIONS", "160,240,320").split(",") if x.strip())
    if len(resolutions) < 2 or any(n < 32 for n in resolutions):
        raise ValueError("at least two resolutions of 32 or greater are required")
    r_max = float(os.environ.get("ORACLE_R_MAX", "80.0"))
    if not np.isfinite(r_max) or r_max <= 10.0:
        raise ValueError("r_max must be finite and greater than 10")
    output = Path(os.environ.get("ORACLE_OUT", "results/continuous-origin-half-cell-oracle-2026-10-15"))
    rows = [run_resolution(n, r_max=r_max) for n in resolutions]
    report = {
        "schema": "v55_continuous_origin_half_cell_oracle_v1",
        "status": "completed_research_only_continuous_integral_oracle",
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
        "all_profiles_finite": all(np.isfinite(x["continuous_reference"]["reference_moment"]) for x in rows),
        "interpretation_limits": [
            "The continuous source assumes the declared analytic Gaussian S/D fields and constant phi, Pi, K_t.",
            "This is the default uncharged initialization and does not include optional charged-dust or Maxwell source terms.",
            "Richardson extrapolation provides a convergence estimate, not a rigorous analytic remainder bound.",
            "Fit-to-continuum errors are conditional on this continuous extension and do not by themselves establish a production correction.",
            "No production or vendor code changed; the audit takes zero evolution steps.",
        ],
    }
    output.mkdir(parents=True, exist_ok=True)
    paths = {}
    paths["report.json"] = output / "report.json"
    paths["report.json"].write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    paths["moments.csv"] = output / "moments.csv"
    with paths["moments.csv"].open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["N", "dr", "r_half", "I_cont", "I_mid_cont", "I_prod", "quadrature_delta", "sampled_source_delta", "total_delta", "oracle_uncertainty"])
        for x in rows:
            m = x["moment_decomposition"]
            w.writerow([x["resolution"], x["dr"], x["r_half"], m["continuous_half_cell_integral_reference"], m["continuous_midpoint_moment"], m["production_point_sampled_moment"], m["continuous_quadrature_delta"], m["continuum_vs_discrete_sampled_source_delta"], m["total_continuum_minus_production_delta"], m["estimated_reference_uncertainty"]])
    paths["fit_family.csv"] = output / "fit_family.csv"
    with paths["fit_family.csv"].open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["N", "name", "status", "degree", "window", "anchor_origin", "sample_count", "fitted_moment", "fit_minus_continuous_reference", "fit_error_abs", "max_abs_sample_residual", "condition_number"])
        for x in rows:
            for fit in x["fit_family"]:
                w.writerow([x["resolution"], fit["name"], fit.get("status", "ok"), fit["degree"], fit["window"], fit["anchor_origin"], fit.get("sample_count"), fit.get("fitted_moment"), fit.get("fit_minus_continuous_reference"), fit.get("fit_error_abs"), fit.get("max_abs_sample_residual"), fit.get("condition_number")])
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
        } for x in rows],
        "state_unchanged": report["all_observers_left_state_unchanged"],
        "finite": report["all_profiles_finite"], "errors": [],
    }, indent=2, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
