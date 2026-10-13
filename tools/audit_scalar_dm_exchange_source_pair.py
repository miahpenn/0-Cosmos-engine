#!/usr/bin/env python3
"""Frozen-state scalar/DM exchange-source identity audit; no evolution step.

The audit exposes the coupling terms separately at the existing production
call boundaries. It contrasts the current implementation with lapse factors
implied by the documented covariant source
  partial_mu(sqrt(-g) T^mu_nu) = geometric_nu + sqrt(-g) Q_nu,
  Q_nu = beta_DM rho_DM partial_nu(phi).
No counterfactual is passed back into the production kernel.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.production_kernel import BETA_DM, V55ProductionKernel
from engine.scalar_system import KAPPA, scalar_rhs_arrays
from engine.matter_system import (
    Species,
    _metric_arrays,
    _species_rhs_components,
    primitives,
)
from engine.valencia import dark_matter_covector_source
from engine.v55_matter import dm_density, metric_slice_from_q

SOURCE_BASELINE = "c01f2467347e8aebd57d220374fd7bc4e93045bc"
EXPECTED_BLOBS = {
    "engine/production_kernel.py": "2954e7a79656849b5f088053ea0f4ac97d4bb443",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/matter_system.py": "4f99e3bb61736473d8b9220145aec9e4167a9cff",
    "engine/valencia.py": "9eba0855936e1974a15314eacdc2f6c378072596",
    "engine/v55_matter.py": "0f71164494e534ca4226c56c10058845d9a14e1d",
    "engine/stress_energy.py": "c31c647f74fab6893254b3dea1366eada15b3110",
}
EXPECTED_VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (60, 80, 120)
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
EPS = np.finfo(np.float64).eps
OUT = Path("results/scalar-dm-exchange-source-pair")
FIXED_RADII = (2.0, 3.0, 5.0, 8.0)


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def arr(value: Any) -> np.ndarray:
    return np.asarray(value, dtype=np.float64)


def max_abs(value: Any) -> float:
    x = arr(value)
    return float(np.max(np.abs(x))) if x.size else 0.0


def finite(name: str, value: Any) -> np.ndarray:
    x = arr(value)
    if not bool(np.all(np.isfinite(x))):
        raise FloatingPointError(f"non-finite values in {name}")
    return x


def digest_state(state) -> str:
    h = hashlib.sha256()
    named: dict[str, Any] = {"grid.centers": state.grid.centers}
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        named["geometry." + name] = getattr(state.geometry, name)
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        named["scalars." + name] = getattr(state.scalars, name)
    for sec in ("dark_matter", "baryons", "radiation"):
        bundle = getattr(state.matter, sec)
        for name in ("rest", "energy_t", "momentum_r"):
            named["matter." + sec + "." + name] = getattr(bundle, name)
    for name, value in sorted(named.items()):
        x = np.ascontiguousarray(arr(value))
        h.update(name.encode("utf-8"))
        h.update(repr(x.shape).encode("ascii"))
        h.update(x.tobytes(order="C"))
    return h.hexdigest()


def provenance() -> dict[str, Any]:
    blobs = {}
    for path, expected in EXPECTED_BLOBS.items():
        source_blob = git("rev-parse", f"{SOURCE_BASELINE}:{path}")
        branch_blob = git("rev-parse", f"HEAD:{path}")
        if source_blob != expected or branch_blob != expected:
            raise RuntimeError(
                f"source blob mismatch {path}: expected={expected}, "
                f"source={source_blob}, branch={branch_blob}"
            )
        blobs[path] = branch_blob
    vendor = git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if vendor != EXPECTED_VENDOR:
        raise RuntimeError(f"vendor mismatch: expected={EXPECTED_VENDOR}, actual={vendor}")
    forbidden_diff = git(
        "diff", "--name-only", SOURCE_BASELINE, "HEAD", "--",
        "engine", "vendor", ".gitmodules", "requirements.txt"
    )
    if forbidden_diff:
        raise RuntimeError("production/vendor/dependency paths changed:\n" + forbidden_diff)
    return {
        "source_baseline": SOURCE_BASELINE,
        "diagnostic_head": git("rev-parse", "HEAD"),
        "production_blobs": blobs,
        "vendor": vendor,
        "forbidden_diff_empty": True,
    }


def profile_summary(values: Any, radii: np.ndarray) -> dict[str, Any]:
    x = finite("summary profile", values)
    i = int(np.argmax(np.abs(x)))
    regions = {
        "full": np.ones(x.size, dtype=bool),
        "centre_i0_i2": np.arange(x.size) < min(3, x.size),
        "r_ge_5": radii >= 5.0,
    }
    out: dict[str, Any] = {
        "max_abs": float(abs(x[i])),
        "signed_at_max": float(x[i]),
        "peak_index": i,
        "peak_r": float(radii[i]),
        "rms": float(np.sqrt(np.mean(x*x))),
        "cell_0": float(x[0]),
        "cell_1": float(x[1]) if x.size > 1 else None,
        "cell_2": float(x[2]) if x.size > 2 else None,
        "fixed_coordinate_radius_samples": {
            str(rad): float(np.interp(rad, radii, x))
            for rad in FIXED_RADII if radii[0] <= rad <= radii[-1]
        },
        "regions": {},
    }
    for name, mask in regions.items():
        if not np.any(mask):
            out["regions"][name] = None
            continue
        idx = np.flatnonzero(mask)
        j = int(idx[np.argmax(np.abs(x[idx]))])
        out["regions"][name] = {
            "max_abs": float(abs(x[j])),
            "signed_at_max": float(x[j]),
            "index": j,
            "r": float(radii[j]),
            "rms": float(np.sqrt(np.mean(x[mask]*x[mask]))),
        }
    return out


def audit_resolution(resolution: int) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution, r_max=R_MAX, amplitude=AMPLITUDE,
        width=WIDTH, D_amplitude=D_AMPLITUDE, include_radiation=True,
    )
    before = digest_state(state)
    grid, geom, scalars, matter = state.grid, state.geometry, state.scalars, state.matter
    r = arr(grid.centers)
    alpha = finite("alpha", geom.alpha)
    beta_shift = finite("shift", geom.beta)
    pi = finite("COSMOS Pi", scalars.Pi)
    phi = finite("COSMOS phi", scalars.phi)
    phi_r = finite("D_r phi", grid.cell_derivative_fourth(phi, parity=1))
    metric = metric_slice_from_q(grid, geom)
    metric_cells = _metric_arrays(metric)
    sqrt_gamma = finite("sqrt_gamma", np.asarray([m.sqrt_gamma for m in metric_cells]))

    # Match the actual scalar call in V55ProductionKernel._rhs.
    rho_dm = finite("rho_DM input", dm_density(metric, matter))
    scalar_with_coupling, matter_rhs, metric_derivatives = kernel._rhs(state)
    scalar_without_coupling = scalar_rhs_arrays(
        grid, geom, scalars, beta_dm=BETA_DM, rho_dm=None
    )
    pi_source_input = finite("Pi-source input beta*rho_DM", BETA_DM * rho_dm)
    pi_source_observed = finite(
        "Pi RHS coupled-minus-uncoupled",
        scalar_with_coupling.Pi - scalar_without_coupling.Pi,
    )
    pi_source_observed_residual = pi_source_observed - pi_source_input
    pi_source_observed_bound = 64.0 * EPS * (
        np.abs(scalar_with_coupling.Pi) + np.abs(scalar_without_coupling.Pi)
        + np.abs(pi_source_input)
    )

    # Match the production matter RHS call and separately retain Q_E / Q_S.
    dm_state = matter.dark_matter
    primitives_dm = primitives(metric_cells, dm_state, Species.DARK_MATTER)
    rho_dm_primitive = finite("DM primitive rho", np.asarray([q.rho for q in primitives_dm]))
    dphi_t = finite("production dphi_t", scalar_with_coupling.phi)
    q_t, q_r = [], []
    for i, primitive in enumerate(primitives_dm):
        qt, qr = dark_matter_covector_source(
            BETA_DM, float(primitive.rho), float(dphi_t[i]), float(phi_r[i])
        )
        q_t.append(qt)
        q_r.append(qr)
    q_t = finite("DM Q_t", q_t)
    q_r = finite("DM Q_r", q_r)

    face_flux, source_e_geom, source_s_geom, q_e_actual, q_s_actual = _species_rhs_components(
        metric, metric_derivatives, dm_state, Species.DARK_MATTER,
        dphi_t=dphi_t, dphi_r=phi_r, beta_dm=BETA_DM,
    )
    q_e_actual = finite("actual DM exchange energy source", q_e_actual)
    q_s_actual = finite("actual DM exchange momentum source", q_s_actual)

    q_e_expected = -sqrt_gamma * (q_t - beta_shift*q_r)
    # The covariant identity has sqrt(-g) Q_r = alpha*sqrt(gamma)*Q_r.
    q_s_expected_covariant = alpha * sqrt_gamma * q_r
    q_s_expected_current_code = sqrt_gamma * q_r
    q_e_reconstruction_residual = q_e_actual - q_e_expected
    q_s_code_residual = q_s_actual - q_s_expected_current_code
    q_s_covariant_residual = q_s_actual - q_s_expected_covariant

    q_pi_actual = pi_source_input
    q_pi_expected_covariant = alpha * BETA_DM * rho_dm
    kappa = float(KAPPA)

    # Source contribution to the normalized scalar Einstein projections.
    scalar_energy_actual = pi * q_pi_actual / kappa
    scalar_current_actual = -q_pi_actual * phi_r / kappa
    scalar_energy_covariant = pi * q_pi_expected_covariant / kappa
    scalar_current_covariant = -q_pi_expected_covariant * phi_r / kappa

    # Actual covariant source increments in the DM Eulerian energy/current channels.
    dm_energy_channel_actual = q_e_actual / (sqrt_gamma * kappa)
    dm_current_channel_actual = q_s_actual / (sqrt_gamma * kappa)
    dm_current_channel_covariant = q_s_expected_covariant / (sqrt_gamma * kappa)

    energy_residual_actual = scalar_energy_actual + dm_energy_channel_actual
    current_residual_actual = scalar_current_actual + dm_current_channel_actual
    energy_residual_scalar_lapse_only = scalar_energy_covariant + dm_energy_channel_actual
    current_residual_scalar_lapse_only = scalar_current_covariant + dm_current_channel_actual
    energy_residual_matter_lapse_only = scalar_energy_actual + dm_energy_channel_actual
    current_residual_matter_lapse_only = scalar_current_actual + dm_current_channel_covariant
    energy_residual_both_lapse = scalar_energy_covariant + dm_energy_channel_actual
    current_residual_both_lapse = scalar_current_covariant + dm_current_channel_covariant

    pi_source_identity_ok = bool(np.all(np.abs(pi_source_observed_residual) <= pi_source_observed_bound))
    q_e_identity_ok = bool(np.all(np.abs(q_e_reconstruction_residual) <=
        64.0*EPS*(np.abs(q_e_actual)+np.abs(q_e_expected))))
    q_s_code_identity_ok = bool(np.all(np.abs(q_s_code_residual) <=
        64.0*EPS*(np.abs(q_s_actual)+np.abs(q_s_expected_current_code))))
    q_s_covariant_identity_ok = bool(np.all(np.abs(q_s_actual-q_s_expected_covariant) <=
        64.0*EPS*(np.abs(q_s_actual)+np.abs(q_s_expected_covariant))))

    final = digest_state(state)
    if before != final:
        raise RuntimeError("frozen state was mutated during the source audit")

    profiles = {
        "r": r, "alpha": alpha, "sqrt_gamma": sqrt_gamma,
        "rho_dm_projection": rho_dm, "rho_dm_primitive": rho_dm_primitive,
        "Pi": pi, "phi_r": phi_r, "dphi_t": dphi_t,
        "Q_t": q_t, "Q_r": q_r,
        "Pi_source_actual": q_pi_actual,
        "Pi_source_covariant_comparator": q_pi_expected_covariant,
        "Pi_source_observed_rhs_difference": pi_source_observed,
        "q_E_actual": q_e_actual, "q_E_covariant_formula": q_e_expected,
        "q_S_actual": q_s_actual, "q_S_current_code_formula": q_s_expected_current_code,
        "q_S_covariant_formula": q_s_expected_covariant,
        "scalar_energy_source_actual": scalar_energy_actual,
        "dm_energy_channel_actual": dm_energy_channel_actual,
        "energy_net_actual": energy_residual_actual,
        "scalar_current_source_actual": scalar_current_actual,
        "dm_current_channel_actual": dm_current_channel_actual,
        "current_net_actual": current_residual_actual,
        "energy_net_scalar_lapse_only": energy_residual_scalar_lapse_only,
        "current_net_scalar_lapse_only": current_residual_scalar_lapse_only,
        "energy_net_matter_lapse_only": energy_residual_matter_lapse_only,
        "current_net_matter_lapse_only": current_residual_matter_lapse_only,
        "energy_net_both_lapse": energy_residual_both_lapse,
        "current_net_both_lapse": current_residual_both_lapse,
        "Pi_source_observed_residual": pi_source_observed_residual,
        "q_E_formula_residual": q_e_reconstruction_residual,
        "q_S_covariant_residual": q_s_covariant_residual,
        "q_S_code_formula_residual": q_s_code_residual,
    }

    report = {
        "resolution": int(resolution), "r_max": R_MAX, "dr": float(grid.dr),
        "t": float(state.t), "steps_taken": 0,
        "initial_state_sha256_before": before, "initial_state_sha256_after": final,
        "state_unchanged": before == final,
        "lapse": {
            "min": float(np.min(alpha)), "max": float(np.max(alpha)),
            "cell_0": float(alpha[0]), "cell_1": float(alpha[1]),
        },
        "dm_density_consistency": {
            "max_abs_dm_density_api_vs_primitive": max_abs(rho_dm-rho_dm_primitive),
            "source_input_density": "dm_density(metric,state.matter), matching the production scalar_rhs_arrays input",
        },
        "actual_rhs_contract": {
            "Pi_source": "scalar_rhs_arrays appends beta_DM*rho_DM to Pi_t, outside the alpha-weighted wave/potential terms",
            "DM_Q_t": "beta_DM*rho_DM*dphi_t; production dphi_t=alpha*Pi at zero shift",
            "DM_Q_r": "beta_DM*rho_DM*dphi_r",
            "q_E": "-sqrt_gamma*(Q_t-beta_shift*Q_r)",
            "q_S_actual": "sqrt_gamma*Q_r",
            "covariant_momentum_source_comparator": "sqrt(-g)*Q_r=alpha*sqrt_gamma*Q_r",
        },
        "identity_checks": {
            "Pi source observed vs input": {
                "max_abs_residual": max_abs(pi_source_observed_residual),
                "max_roundoff_bound": float(np.max(pi_source_observed_bound)),
                "within_64eps_operand_bound": pi_source_identity_ok,
            },
            "DM energy source formula": {
                "max_abs_residual": max_abs(q_e_reconstruction_residual),
                "within_64eps_operand_bound": q_e_identity_ok,
            },
            "DM current source vs current code formula": {
                "max_abs_residual": max_abs(q_s_code_residual),
                "within_64eps_operand_bound": q_s_code_identity_ok,
            },
            "DM current source vs covariant sqrt(-g)Q_r": {
                "max_abs_residual": max_abs(q_s_covariant_residual),
                "within_64eps_operand_bound": q_s_covariant_identity_ok,
                "note": "A failing check is a measured difference from the covariant RHS convention, not a runtime arithmetic failure.",
            },
        },
        "source_pair_closure": {
            "actual_code_energy": profile_summary(energy_residual_actual, r),
            "actual_code_current": profile_summary(current_residual_actual, r),
            "scalar_lapse_only_energy": profile_summary(energy_residual_scalar_lapse_only, r),
            "scalar_lapse_only_current": profile_summary(current_residual_scalar_lapse_only, r),
            "matter_lapse_only_energy": profile_summary(energy_residual_matter_lapse_only, r),
            "matter_lapse_only_current": profile_summary(current_residual_matter_lapse_only, r),
            "both_lapse_weighted_energy": profile_summary(energy_residual_both_lapse, r),
            "both_lapse_weighted_current": profile_summary(current_residual_both_lapse, r),
            "interpretation": (
                "Source-only algebraic comparators. Lapse-weighted entries are counterfactual identities only; "
                "they are not executed by production and do not authorize a source change."
            ),
        },
        "scientific_admission": False,
    }
    return report, profiles


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    try:
        prov = provenance()
        results = []
        profiles_all: dict[str, np.ndarray] = {}
        for n in RESOLUTIONS:
            result, profiles = audit_resolution(n)
            results.append(result)
            for name, value in profiles.items():
                profiles_all[f"N{n}_{name}"] = arr(value)
        report = {
            "schema": "current_main_scalar_dm_separated_source_pair_v1",
            "status": "completed_diagnostic_only",
            "provenance": prov,
            "initialization": {
                "resolutions": list(RESOLUTIONS), "r_max": R_MAX,
                "amplitude": AMPLITUDE, "width": WIDTH, "D_amplitude": D_AMPLITUDE,
                "radiation": True, "beta_DM": BETA_DM, "KAPPA": KAPPA,
            },
            "resolution_reports": results,
            "scientific_admission": False,
            "interpretation_limits": [
                "No kernel.step call is made; all checks are on the initialized accepted state at t=0.",
                "The report separates actual source injections from lapse-weighted covariant comparators.",
                "A source-only algebraic closure is not a full stress-energy conservation or trajectory test.",
                "Counterfactual factors are evaluated only in report arrays; production code is never patched or monkeypatched.",
                "No production-repair conclusion is authorized from this frozen-state identity audit alone.",
            ],
        }
        report["artifact_manifest"] = {
            "source_commit": SOURCE_BASELINE,
            "vendor_commit": EXPECTED_VENDOR,
            "state_hashes": {f"N{x['resolution']}": x["initial_state_sha256_before"] for x in results},
        }
        report_path = OUT/"report.json"
        report_path.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
        profile_paths = []
        for n in RESOLUTIONS:
            fields = sorted(k[len(f"N{n}_"):] for k in profiles_all if k.startswith(f"N{n}_") and k != f"N{n}_r")
            radii = profiles_all[f"N{n}_r"]
            ppath = OUT/f"profiles_N{n}.csv"
            with ppath.open("w",newline="",encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(["i","r"]+fields)
                for i,rv in enumerate(radii):
                    writer.writerow([i,float(rv)]+[float(profiles_all[f"N{n}_{field}"][i]) for field in fields])
            profile_paths.append(ppath)
        index_path=OUT/"profiles_index.csv"
        with index_path.open("w",newline="",encoding="utf-8") as stream:
            writer=csv.writer(stream)
            writer.writerow(["resolution","profile_file","cells","profile_columns"])
            for n,ppath in zip(RESOLUTIONS,profile_paths):
                writer.writerow([n,ppath.name,len(profiles_all[f"N{n}_r"]),len([k for k in profiles_all if k.startswith(f"N{n}_")])])
        files=[report_path,index_path,*profile_paths]
        (OUT/"SHA256SUMS.txt").write_text(
            "".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in files),
            encoding="utf-8",
        )
        print("SCALAR_DM_EXCHANGE_SOURCE_PAIR_PASS")
        print("source=" + SOURCE_BASELINE + " vendor=" + EXPECTED_VENDOR + " head=" + prov["diagnostic_head"])
        print("summary=" + json.dumps([
            {
                "N":x["resolution"],"alpha0":x["lapse"]["cell_0"],
                "actual_energy_net_max":x["source_pair_closure"]["actual_code_energy"]["max_abs"],
                "actual_current_net_max":x["source_pair_closure"]["actual_code_current"]["max_abs"],
                "qS_covariant_mismatch":x["identity_checks"]["DM current source vs covariant sqrt(-g)Q_r"]["max_abs_residual"],
                "joint_lapse_energy_net_max":x["source_pair_closure"]["both_lapse_weighted_energy"]["max_abs"],
                "joint_lapse_current_net_max":x["source_pair_closure"]["both_lapse_weighted_current"]["max_abs"],
            } for x in results
        ],sort_keys=True))
        print("report_sha256="+hashlib.sha256(report_path.read_bytes()).hexdigest())
        return 0
    except Exception as exc:
        failure={"schema":"current_main_scalar_dm_separated_source_pair_v1","status":"failed_closed",
                 "error_type":type(exc).__name__,"error":str(exc),"scientific_admission":False}
        OUT.mkdir(parents=True,exist_ok=True)
        (OUT/"failure.json").write_text(json.dumps(failure,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
        print("SCALAR_DM_EXCHANGE_SOURCE_PAIR_FAIL_CLOSED "+json.dumps(failure,sort_keys=True))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
