#!/usr/bin/env python3
"""Research-only semidiscrete fixed-R Misner-Sharp mass-tangent budget.

Unlike the endpoint finite-difference probes, this evaluates the instantaneous
mass tangent from the exact pinned metric and K/Aa RHS at accepted/stage states.
For chi = gamma^rr (D_r R)^2 - (R K_theta)^2, its coordinate-time tangent is
split into gamma^rr_dot, D_r(R_t), R_t K_theta^2, K_theta_dot(K_dot,Aa_dot),
then transported to fixed areal radius:

  m_dot|R = -R/2 * [chi_dot|r - (R_t/R_r) D_r(chi)].

Both direct D_r(chi) and the chain-rule reconstruction of chi_r are reported.
The same fixed-R samples are compared with F=4*pi*R^2*T_t^R. This separates
the semidiscrete equation-of-motion identity from finite-difference time
truncation and exposes the discrete product-rule remainder. Diagnostic only:
no production equations or physical-source convention are changed.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))

from tools import audit_bssn_areal_misner_sharp_current as areal
from tools import audit_bssn_areal_stage_current_balance as stage_balance
from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
N = int(os.environ.get("AREAL_MASS_TANGENT_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_MASS_TANGENT_OUT", f"results/bssn-areal-semidiscrete-mass-tangent-N{N}-2026-10-25"))
TARGET_TIME = 0.60
CFL = 0.03
FOUR_PI = 4.0 * np.pi
FIXED_RADII = (0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0)
STAGE_LABELS = (
    "accepted_start",
    "kernel_start_stage",
    "forward_euler_predictor",
    "kernel_corrected_stage",
    "accepted_output",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(name: str, values: Any) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(out)):
        raise RuntimeError(f"{name} contains non-finite values")
    return out


def scalar(name: str, value: Any) -> float:
    out = float(value)
    if not np.isfinite(out):
        raise RuntimeError(f"{name} is non-finite")
    return out


def sample(profile: dict[str, Any], values: np.ndarray, radii: np.ndarray, order: str, name: str) -> np.ndarray:
    coordinates = finite("areal-coordinate sample locations", profile["R"])
    values = finite(f"{name} native profile", values)
    if order == "linear":
        out = np.interp(radii, coordinates, values)
    elif order == "cubic":
        out = np.asarray([areal.interp4(float(R), coordinates, values) for R in radii], dtype=np.float64)
    else:
        raise ValueError(f"unsupported interpolation order {order}")
    return finite(f"{order} fixed-R samples of {name}", out)


def error_summary(values: np.ndarray, radii: np.ndarray) -> dict[str, float]:
    values = finite("summary values", values)
    inner = radii <= 3.0
    annulus = (radii >= 2.0) & (radii <= 10.0)
    if not inner.any() or not annulus.any():
        raise RuntimeError("common fixed-R samples do not cover the inner and annular regions")
    return {
        "max_abs_R_le_3": float(np.max(np.abs(values[inner]))),
        "rms_R_2_to_10": float(np.sqrt(np.mean(values[annulus] ** 2))),
        "max_abs_R_2_to_10": float(np.max(np.abs(values[annulus]))),
    }


def semidiscrete_tangent(kernel, state) -> dict[str, Any]:
    # Reuse the pinned observer's gauge refresh and native geometric definition,
    # then independently form the analytic semidiscrete tangent from actual RHSs.
    profile = areal.geometry_rate_and_mass(kernel, state)
    q = profile["state"]
    grid = q.grid
    g = q.geometry
    r = finite("coordinate centers", grid.centers)
    R = finite("areal radius", profile["R"])
    Rr = finite("areal-radius radial derivative", profile["Rr"])
    Rt = finite("areal-radius time derivative", profile["Rt"])
    alpha = finite("lapse", g.alpha)
    a = finite("a", g.a)
    b = finite("b", g.b)
    X = finite("X", g.X)
    K = finite("K", g.K)
    Aa = finite("Aa", g.Aa)
    beta = finite("beta", g.beta)
    if np.max(np.abs(beta)) > 1e-14:
        raise RuntimeError("semidiscrete fixed-R tangent requires the pinned zero-shift lane")

    terms = kernel._geometry_stage_terms(
        grid, g, q.scalars, q.matter, q.electromagnetic
    )
    explicit = terms["explicit"]
    l2 = terms["primary_l2"]
    l3 = terms["primary_l3"]
    adot = finite("a RHS", explicit["a"])
    bdot = finite("b RHS", explicit["b"])
    Xdot = finite("X RHS", explicit["X"])
    Kdot = finite("K RHS L2+L3", finite("K L2", l2["K"]) + finite("K L3", l3["K"]))
    Aadot = finite("Aa RHS L2+L3", finite("Aa L2", l2["Aa"]) + finite("Aa L3", l3["Aa"]))

    gamma_inv = finite("gamma^rr", X * X / a)
    gamma_inv_dot = finite(
        "gamma^rr tangent",
        gamma_inv * (2.0 * Xdot / X - adot / a),
    )
    gamma_inv_r = finite("D_r gamma^rr", grid.cell_derivative_fourth(gamma_inv, parity=1))
    Rrr = finite("D_r^2 R", grid.cell_second_derivative_fourth(R, parity=-1))
    Rt_r = finite("D_r R_t", grid.cell_derivative_fourth(Rt, parity=-1))
    Ktheta = finite("Ktheta", K / 3.0 - Aa / 2.0)
    Ktheta_r = finite("D_r Ktheta", grid.cell_derivative_fourth(Ktheta, parity=1))
    Ktheta_dot = finite("Ktheta tangent", Kdot / 3.0 - Aadot / 2.0)

    chi_native = finite("chi reconstructed from Rr and Ktheta", gamma_inv * Rr**2 - (R * Ktheta)**2)
    chi_gap = finite("native chi reconstruction gap", chi_native - finite("observer chi", profile["chi"]))
    if float(np.max(np.abs(chi_gap))) > 1e-10:
        raise RuntimeError(f"chi kinematic reassembly failed: max gap={np.max(np.abs(chi_gap))}")

    # Exact RHS derivative of the continuum expression, evaluated with the
    # pinned discrete metric RHS and derivative operators.
    chi_dot_terms = {
        "gamma_inv_dot": finite("chi_dot gamma_inv term", gamma_inv_dot * Rr**2),
        "D_r_Rt": finite("chi_dot D_r(Rt) term", 2.0 * gamma_inv * Rr * Rt_r),
        "Rt_Ktheta_geometry": finite("chi_dot Rt Ktheta geometry term", -2.0 * R * Rt * Ktheta**2),
        "Ktheta_dot": finite("chi_dot Ktheta_dot term", -2.0 * R**2 * Ktheta * Ktheta_dot),
    }
    chi_dot = finite("semidiscrete chi tangent", sum(chi_dot_terms.values()))

    # Two independent spatial derivatives of chi: D4 applied to the composed
    # chi field, versus a chain-rule reconstruction from its primitive factors.
    chi_r_direct = finite("D_r chi direct", grid.cell_derivative_fourth(chi_native, parity=1))
    chi_r_chain_terms = {
        "gamma_inv_r": finite("chi_r gamma_inv term", gamma_inv_r * Rr**2),
        "Rrr": finite("chi_r Rrr term", 2.0 * gamma_inv * Rr * Rrr),
        "R_Ktheta_geometry": finite("chi_r R Ktheta geometry term", -2.0 * R * Rr * Ktheta**2),
        "Ktheta_r": finite("chi_r Ktheta_r term", -2.0 * R**2 * Ktheta * Ktheta_r),
    }
    chi_r_chain = finite("chain-rule chi_r", sum(chi_r_chain_terms.values()))
    chi_r_product_rule_gap = finite("D_r chi minus chain-rule chi_r", chi_r_direct - chi_r_chain)

    # Fixed-R derivative: dt|R = dt|r - (R_t/R_r) dr|t.
    coordinate_velocity_at_fixed_R = finite("coordinate velocity at fixed R", -Rt / Rr)
    mass_rate_terms = {
        "gamma_inv_dot": finite("mass tangent gamma_inv term", -0.5 * R * chi_dot_terms["gamma_inv_dot"]),
        "D_r_Rt": finite("mass tangent D_r(Rt) term", -0.5 * R * chi_dot_terms["D_r_Rt"]),
        "Rt_Ktheta_geometry": finite("mass tangent Rt Ktheta term", -0.5 * R * chi_dot_terms["Rt_Ktheta_geometry"]),
        "Kdot": finite("mass tangent Kdot term", (R**3) * Ktheta * Kdot / 3.0),
        "Aadot": finite("mass tangent Aadot term", -0.5 * (R**3) * Ktheta * Aadot),
        "fixed_R_advection_direct_D4_chi": finite(
            "mass tangent fixed-R advection term", 0.5 * R * (Rt / Rr) * chi_r_direct
        ),
        "fixed_R_advection_chain_rule_chi_r": finite(
            "mass tangent chain fixed-R advection term", 0.5 * R * (Rt / Rr) * chi_r_chain
        ),
    }
    mass_rate_direct = finite(
        "semidiscrete fixed-R mass rate direct D4 chi",
        -0.5 * R * (chi_dot - (Rt / Rr) * chi_r_direct),
    )
    mass_rate_chain = finite(
        "semidiscrete fixed-R mass rate chain chi_r",
        -0.5 * R * (chi_dot - (Rt / Rr) * chi_r_chain),
    )
    # The term ledger's direct lane must reconstruct the mass-rate profile.
    direct_sum = finite("mass-rate component sum", sum(
        mass_rate_terms[k] for k in (
            "gamma_inv_dot", "D_r_Rt", "Rt_Ktheta_geometry", "Kdot", "Aadot",
            "fixed_R_advection_direct_D4_chi",
        )
    ))
    ledger_gap = finite("mass-rate component closure", direct_sum - mass_rate_direct)
    ledger_scale = max(float(np.max(np.abs(direct_sum))), float(np.max(np.abs(mass_rate_direct))), 1e-300)
    if float(np.max(np.abs(ledger_gap))) > 2e-10 * ledger_scale:
        raise RuntimeError(f"semidiscrete mass-rate component ledger does not close: {np.max(np.abs(ledger_gap))}")

    current = finite("transformed current F", profile["rhs"])
    rhs_component_sum = finite("sum of transformed current components", sum(profile["rhs_components"].values()))
    current_component_gap = finite("current component closure", current - rhs_component_sum)

    result = {
        "R": R.copy(),
        "chi": chi_native,
        "chi_dot": chi_dot,
        "chi_dot_terms": chi_dot_terms,
        "chi_r_direct": chi_r_direct,
        "chi_r_chain": chi_r_chain,
        "chi_r_chain_terms": chi_r_chain_terms,
        "chi_r_product_rule_gap": chi_r_product_rule_gap,
        "Rr": Rr,
        "Rt": Rt,
        "Rt_r": Rt_r,
        "Ktheta": Ktheta,
        "Ktheta_dot": Ktheta_dot,
        "Kdot": Kdot,
        "Aadot": Aadot,
        "gamma_inv_dot": gamma_inv_dot,
        "mass_rate_terms": mass_rate_terms,
        "mass_rate_direct": mass_rate_direct,
        "mass_rate_chain": mass_rate_chain,
        "mass_rate_direct_minus_chain": finite("mass rate direct minus chain", mass_rate_direct - mass_rate_chain),
        "mass_rate_component_closure_gap": ledger_gap,
        "current": current,
        "current_component_gap": current_component_gap,
        "profile": profile,
    }
    return result


def compare_arm(kernel, accepted, dt: float) -> dict[str, Any]:
    accepted_digest = saved.tangent.state_digest(accepted)
    candidate, captured, captured_digests = stage_balance.capture_stage_states(kernel, accepted, dt)
    raw_states = [accepted, captured[0], captured[1], captured[2], candidate]
    profiles = {label: semidiscrete_tangent(kernel, state) for label, state in zip(STAGE_LABELS, raw_states)}
    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("accepted state mutated during semidiscrete tangent audit")
    for label, state, digest in zip(STAGE_LABELS[1:4], captured, captured_digests):
        if saved.tangent.state_digest(state) != digest:
            raise RuntimeError(f"{label} mutated during semidiscrete tangent observation")

    lo = max(float(item["R"][0]) for item in profiles.values())
    hi = min(min(float(item["R"][-1]) for item in profiles.values()), 10.0)
    if not hi > lo:
        raise RuntimeError("no common fixed-areal-radius range across captured stages")
    radii = np.linspace(lo, hi, max(300, int(accepted.grid.n) * 2))

    report: dict[str, Any] = {
        "dt": float(dt),
        "accepted_start_time": float(accepted.t),
        "accepted_output_time": float(candidate.t),
        "captured_stage_labels": list(STAGE_LABELS[1:4]),
        "accepted_input_immutable": True,
        "captured_stages_immutable": True,
        "common_areal_radius_interval": [float(lo), float(hi)],
        "common_samples": int(radii.size),
        "by_stage": {},
        "fixed_R": {},
        "endpoint_vs_tangent_integral": {},
    }

    for label, item in profiles.items():
        one = {"linear": {}, "cubic": {}}
        for order in ("linear", "cubic"):
            mass_direct = sample(item, item["mass_rate_direct"], radii, order, f"{label}.mass_rate_direct")
            mass_chain = sample(item, item["mass_rate_chain"], radii, order, f"{label}.mass_rate_chain")
            current = sample(item, item["current"], radii, order, f"{label}.current")
            product_gap = sample(item, item["mass_rate_direct_minus_chain"], radii, order, f"{label}.direct_minus_chain")
            one[order] = {
                "direct_tangent_vs_current": error_summary(mass_direct - current, radii),
                "chain_tangent_vs_current": error_summary(mass_chain - current, radii),
                "direct_minus_chain_mass_rate": error_summary(product_gap, radii),
                "Aadot_mass_rate_contribution": error_summary(
                    sample(item, item["mass_rate_terms"]["Aadot"], radii, order, f"{label}.Aadot_term"), radii
                ),
                "Kdot_mass_rate_contribution": error_summary(
                    sample(item, item["mass_rate_terms"]["Kdot"], radii, order, f"{label}.Kdot_term"), radii
                ),
                "D_r_Rt_mass_rate_contribution": error_summary(
                    sample(item, item["mass_rate_terms"]["D_r_Rt"], radii, order, f"{label}.DrRt_term"), radii
                ),
                "fixed_R_advection_mass_rate_contribution": error_summary(
                    sample(item, item["mass_rate_terms"]["fixed_R_advection_direct_D4_chi"], radii, order, f"{label}.advection_term"), radii
                ),
                "chi_r_direct_vs_chain": error_summary(
                    sample(item, item["chi_r_product_rule_gap"], radii, order, f"{label}.chi_r_product_rule_gap"), radii
                ),
                "current_component_closure": error_summary(
                    sample(item, item["current_component_gap"], radii, order, f"{label}.current_component_gap"), radii
                ),
            }
        report["by_stage"][label] = one

    # Integrate the instantaneous tangent using the production step's start and
    # corrected-stage pairing, compare with accepted endpoint mass change, and
    # separately compare that prediction with the actual stage-current average.
    start = profiles["accepted_start"]
    end = profiles["accepted_output"]
    stage0 = profiles["kernel_start_stage"]
    stage1 = profiles["kernel_corrected_stage"]
    for order in ("linear", "cubic"):
        m0 = sample(start, start["chi"], radii, order, f"{order}.accepted_start.chi")
        m1 = sample(end, end["chi"], radii, order, f"{order}.accepted_output.chi")
        R_fixed = radii
        mass_rate_fd = finite(f"{order}.endpoint mass rate", -0.5 * R_fixed * (m1 - m0) / dt)
        tangent0 = sample(start, start["mass_rate_direct"], radii, order, f"{order}.start tangent")
        tangent1 = sample(stage1, stage1["mass_rate_direct"], radii, order, f"{order}.corrected-stage tangent")
        tangent_trap = finite(f"{order}.stage-integrated tangent", 0.5 * (tangent0 + tangent1))
        F0 = sample(stage0, stage0["current"], radii, order, f"{order}.start-stage F")
        F1 = sample(stage1, stage1["current"], radii, order, f"{order}.corrected-stage F")
        current_trap = finite(f"{order}.stage-integrated current", 0.5 * (F0 + F1))
        report["endpoint_vs_tangent_integral"][order] = {
            "accepted_endpoint_mass_rate": error_summary(mass_rate_fd, radii),
            "start_corrected_tangent_trapezoid_vs_endpoint_rate": error_summary(tangent_trap - mass_rate_fd, radii),
            "actual_stage_current_trapezoid_vs_endpoint_rate": error_summary(current_trap - mass_rate_fd, radii),
            "tangent_trapezoid_vs_stage_current_trapezoid": error_summary(tangent_trap - current_trap, radii),
        }

    for radius in FIXED_RADII:
        if not lo <= radius <= hi:
            continue
        report["fixed_R"][str(radius)] = {}
        for label, item in profiles.items():
            row = {}
            for order in ("linear", "cubic"):
                def at(field: str, field_name: str) -> float:
                    arr = sample(item, item[field], np.asarray([radius]), order, f"{label}.{field_name}")
                    return scalar(f"{label}.{field_name}@R={radius}", arr[0])
                row[order] = {
                    "mass_rate_direct": at("mass_rate_direct", "mass_rate_direct"),
                    "mass_rate_chain": at("mass_rate_chain", "mass_rate_chain"),
                    "current": at("current", "current"),
                    "residual_direct": at("mass_rate_direct", "mass_rate_direct") - at("current", "current"),
                    "residual_chain": at("mass_rate_chain", "mass_rate_chain") - at("current", "current"),
                    "Aadot_mass_rate_contribution": at_from_array(
                        item, "mass_rate_terms", "Aadot", radius, order, f"{label}.Aadot"
                    ),
                    "Kdot_mass_rate_contribution": at_from_array(
                        item, "mass_rate_terms", "Kdot", radius, order, f"{label}.Kdot"
                    ),
                    "D_r_Rt_mass_rate_contribution": at_from_array(
                        item, "mass_rate_terms", "D_r_Rt", radius, order, f"{label}.D_r_Rt"
                    ),
                    "fixed_R_advection_mass_rate_contribution": at_from_array(
                        item, "mass_rate_terms", "fixed_R_advection_direct_D4_chi", radius, order, f"{label}.advection"
                    ),
                }
            report["fixed_R"][str(radius)][label] = row

    return report


def at_from_array(item: dict[str, Any], parent: str, key: str, radius: float, order: str, name: str) -> float:
    values = item[parent][key]
    result = sample(item, values, np.asarray([radius]), order, name)
    return scalar(f"{name}@R={radius}", result[0])


def main() -> int:
    if N not in (60, 120, 240):
        raise ValueError("resolution must be 60, 120 or 240")
    pins = origin_shadow.verify_source_pins()
    if pins.get("source_commit") != SOURCE or pins.get("vendor_commit") != VENDOR:
        raise RuntimeError("source/vendor pin mismatch")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("protected production/vendor/dependency diff is not empty")
    manifest = saved.check_artifact(FULL_DIR, N, full_profile=True)
    stored = saved.load_profile_arrays(FULL_DIR)
    prefix = f"N{N}_full_continuous_metric_profile_final_"
    kernel, accepted = saved.reconstruct_state(N, stored, prefix)
    origin_shadow._install_regular_row(kernel)
    if abs(float(accepted.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("accepted input is not at t=0.60")
    digest = saved.tangent.state_digest(accepted)
    dt = CFL * float(accepted.grid.dr)
    arms = [{"arm": label, **compare_arm(kernel, accepted, h)} for label, h in (("dt", dt), ("dt_half", dt / 2.0))]
    if saved.tangent.state_digest(accepted) != digest:
        raise RuntimeError("accepted state changed across semidiscrete tangent arms")

    report = {
        "schema": "bssn-areal-semidiscrete-mass-tangent/v1",
        "status": "completed_diagnostic_only",
        "source_pins": {
            "source_commit": SOURCE,
            "vendor_commit": VENDOR,
            "live_pin_check": pins,
            "accepted_artifact_manifest": manifest,
            "accepted_state_digest_before_after_match": True,
        },
        "runtime": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform()},
        "configuration": {
            "resolution": N,
            "accepted_time": TARGET_TIME,
            "dt_rule": "dt=0.03*dr; dt and dt/2",
            "R_t": "R*(0.5*b_dot/b-X_dot/X), checked against -alpha*R*(K/3-Aa/2)",
            "chi": "gamma^rr*(D_r R)^2-(R*Ktheta)^2, Ktheta=K/3-Aa/2",
            "chi_t": "gamma^rr_dot*Rr^2+2*gamma^rr*Rr*D_r(R_t)-2*R*R_t*Ktheta^2-2*R^2*Ktheta*Ktheta_dot",
            "fixed_R_transport": "m_dot|R=-R/2*(chi_dot|r-(R_t/R_r)*D_r(chi))",
            "Ktheta_dot": "K_dot/3-Aa_dot/2 from the same pinned primary-L2+primary-L3 stage RHS",
            "spatial_derivative_controls": "direct D_r(chi) versus chain-rule reconstruction of chi_r",
            "current": "same pinned total-matter transform F=4*pi*R^2*T_t^R",
        },
        "arms": arms,
        "interpretation_guard": [
            "This is an instantaneous semidiscrete identity test; it is not a finite-difference-in-time residual.",
            "The K/Aa RHS is represented by the captured pinned L2+L3 semidiscrete stage sum. The PIRK accepted-step update has a split temporal quadrature, separately audited in endpoint_vs_tangent_integral.",
            "The direct and chain-rule chi_r variants are both discrete operators; neither is declared authoritative by its residual.",
            "No production source, initializer, vendor, dependency, source normalization, or scientific-admission gate was changed.",
        ],
        "acceptance": {
            "source_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_input_immutable": True,
            "captured_stages_immutable": True,
            "semidiscrete_mass_tangent_reconstructed": True,
            "Aadot_and_Kdot_contributions_explicit": True,
            "direct_vs_chain_spatial_gradient_compared": True,
            "integrated_current_identity_closed": False,
            "full_matter_coupled_propagation_identity_closed": False,
            "total_stress_energy_conservation_closed": False,
            "scientific_admission": False,
        },
        "scientific_admission": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "report.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "SHA256SUMS.txt").write_text(f"{sha256(path)}  report.json\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "resolution": N,
        "scientific_admission": False,
        "arms": [{
            "arm": arm["arm"],
            "by_stage": {
                stage: {order: {
                    "direct_tangent_vs_current": result["direct_tangent_vs_current"],
                    "chain_tangent_vs_current": result["chain_tangent_vs_current"],
                    "direct_minus_chain": result["direct_minus_chain_mass_rate"],
                    "Aadot_contribution": result["Aadot_mass_rate_contribution"],
                    "D_r_Rt_contribution": result["D_r_Rt_mass_rate_contribution"],
                    "chi_r_product_rule_gap": result["chi_r_direct_vs_chain"],
                } for order, result in methods.items()}
                for stage, methods in arm["by_stage"].items()
            },
            "endpoint_vs_tangent_integral": arm["endpoint_vs_tangent_integral"],
            "fixed_R": arm["fixed_R"],
        } for arm in arms],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
