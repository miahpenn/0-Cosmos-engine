#!/usr/bin/env python3
"""Capture centre-series and primary-Aa RHS witnesses on initialized production slices.

This is an initialization-only audit: it takes no time step, runs no campaign,
and never modifies engine or vendor source. It records the actual initialized
geometry/gauge fields and exact current production RHS components.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine.production_kernel import V55ProductionKernel
from engine import v55_pirk_adapter as adapter
from engine.stress_energy import assemble_total_stress_energy

SOURCE_BASELINE = "29d828b4d72b017c12a20b6dc90549f2511f025c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_BLOBS = {
    "engine/cmc_gauge.py": "b1e6b3e95ff37702509e1b38ec6401d10e2482c3",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/matter_system.py": "4f99e3bb61736473d8b9220145aec9e4167a9cff",
    "engine/production_kernel.py": "2954e7a79656849b5f088053ea0f4ac97d4bb443",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/stress_energy.py": "c31c647f74fab6893254b3dea1366eada15b3110",
    "engine/v55_initial.py": "4cfbfe4be6ac3ce6ed492f4e435186997d94116a",
    "engine/v55_matter.py": "0f71164494e534ca4226c56c10058845d9a14e1d",
    "engine/v55_pirk_adapter.py": "f671b838d9aefd8926702d9c2235209f52cd6815",
    "engine/valencia.py": "9eba0855936e1974a15314eacdc2f6c378072596",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
}
RESOLUTIONS = (60, 80, 120)
R_MAX = 80.0
SAMPLE_RADII = (0.5, 1.0, 2.0, 3.0, 5.0, 8.0)
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "Aa", "K", "Lambda", "beta", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES_NAMES = ("dark_matter", "baryons", "radiation")
SPECIES_FIELDS = ("rest", "energy_t", "momentum_r")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def digest_state(state) -> str:
    digest = hashlib.sha256()
    arrays = [("r", state.grid.centers)]
    arrays.extend((f"geometry.{key}", getattr(state.geometry, key)) for key in GEOMETRY_FIELDS)
    arrays.extend((f"scalar.{key}", getattr(state.scalars, key)) for key in SCALAR_FIELDS)
    for species_name in SPECIES_NAMES:
        species = getattr(state.matter, species_name)
        arrays.extend((f"matter.{species_name}.{key}", getattr(species, key)) for key in SPECIES_FIELDS)
    for name, value in arrays:
        a = np.ascontiguousarray(np.asarray(value, dtype=np.float64))
        digest.update(name.encode() + b"\0")
        digest.update(str(a.shape).encode() + b"\0")
        digest.update(a.tobytes())
    return digest.hexdigest()


def stats(values: np.ndarray, r: np.ndarray, count: int = 6) -> dict:
    v = np.asarray(values, dtype=float)
    if v.shape != r.shape or not np.all(np.isfinite(v)):
        raise ValueError("profile must be finite and match the radius array")
    n = min(count, len(v))
    samples = [
        {"i": int(i), "r": float(r[i]), "value": float(v[i])}
        for i in range(n)
    ]
    return {
        "cell_0": float(v[0]),
        "first_cells": samples,
        "fixed_radius_samples": {
            f"{target:g}": {
                "i": int(i),
                "r": float(r[i]),
                "value": float(v[i]),
            }
            for target in SAMPLE_RADII
            for i in [int(np.argmin(np.abs(r - target)))]
        },
        "max_abs": float(np.max(np.abs(v))),
    }


def require_finite_profiles(named_arrays: dict[str, np.ndarray]) -> None:
    for name, values in named_arrays.items():
        a = np.asarray(values, dtype=float)
        if a.ndim != 1 or not np.all(np.isfinite(a)):
            raise ValueError(f"non-finite or non-vector diagnostic profile: {name}")


def audit_resolution(n: int, outdir: Path) -> dict:
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=n,
        r_max=R_MAX,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    grid, g = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=float)
    if state.t != 0.0:
        raise RuntimeError("initializer unexpectedly advanced coordinate time")
    if len(state.history) != 0:
        raise RuntimeError("initialization unexpectedly appended evolution history")
    before = digest_state(state)

    _, vacuum, _ = adapter.vendor_modules()
    geom = vacuum.geometry_terms(grid, g)
    l2 = vacuum.primary_l2_rhs(grid, g)
    l3_vac = vacuum.primary_l3_rhs(grid, g)
    l3_total = adapter.primary_l3_with_matter(grid, g, state.scalars, state.matter)
    matter_l3_aa = np.asarray(l3_total["Aa"]) - np.asarray(l3_vac["Aa"])

    lapse_hessian_aa = -(
        np.asarray(geom["DrDr_alpha"]) - np.asarray(geom["lap_alpha"]) / 3.0
    )
    ricci_tracefree_aa = np.asarray(g.alpha) * (
        np.asarray(geom["Rr_r"]) - np.asarray(geom["R"]) / 3.0
    )
    reconstructed_l2 = lapse_hessian_aa + ricci_tracefree_aa
    primary_l2_aa = np.asarray(l2["Aa"], dtype=float)
    primary_l3_aa = np.asarray(l3_total["Aa"], dtype=float)
    full_aa = primary_l2_aa + primary_l3_aa

    stress = assemble_total_stress_energy(grid, g, state.scalars, state.matter)
    constraint = vacuum.constraints(grid, g)
    H = np.asarray(constraint["hamiltonian"]) - 16.0 * math.pi * np.asarray(stress.rho)
    M = np.asarray(constraint["momentum"]) - 8.0 * math.pi * np.asarray(stress.j)
    C = np.asarray(constraint["connection"])
    det = np.asarray(constraint["determinant"])

    d1 = lambda values, parity: grid.cell_derivative_fourth(values, parity=parity)
    even_quotients = {
        f"{name}_radial_derivative_over_r": d1(np.asarray(getattr(g, name)), 1) / r
        for name in ("a", "b", "X", "alpha", "K")
    }
    even_quotients.update({
        f"scalar_{name}_radial_derivative_over_r": d1(np.asarray(getattr(state.scalars, name)), 1) / r
        for name in ("S", "D", "phi")
    })
    radius2 = r * r
    regularity = {
        "metric_ratio_even_quotient_(1-a/b)/r2": (1.0 - g.a / g.b) / radius2,
        "Aa_even_quotient_Aa/r2": np.asarray(g.Aa) / radius2,
        "Lambda_odd_quotient_Lambda/r": np.asarray(g.Lambda) / r,
        **even_quotients,
        "tracefree_Ricci_source_over_r2": ricci_tracefree_aa / radius2,
        "lapse_Hessian_Aa_source_over_r2": lapse_hessian_aa / radius2,
        "primary_L2_Aa_over_r2": primary_l2_aa / radius2,
        "primary_L3_with_matter_Aa_over_r2": primary_l3_aa / radius2,
        "full_Aa_RHS_over_r2": full_aa / radius2,
    }
    profiles = {
        "r": r,
        **{f"geometry__{name}": np.asarray(getattr(g, name), dtype=float) for name in GEOMETRY_FIELDS},
        **{f"scalar__{name}": np.asarray(getattr(state.scalars, name), dtype=float) for name in SCALAR_FIELDS},
        **{
            f"matter__{species_name}__{field}": np.asarray(getattr(getattr(state.matter, species_name), field), dtype=float)
            for species_name in SPECIES_NAMES for field in SPECIES_FIELDS
        },
        "R": np.asarray(geom["R"]),
        "Rr_r": np.asarray(geom["Rr_r"]),
        "lapse_hessian_Aa": lapse_hessian_aa,
        "ricci_tracefree_Aa": ricci_tracefree_aa,
        "primary_L2_Aa": primary_l2_aa,
        "primary_L3_vacuum_Aa": np.asarray(l3_vac["Aa"]),
        "primary_L3_with_matter_Aa": primary_l3_aa,
        "matter_L3_Aa_increment": matter_l3_aa,
        "full_Aa_RHS": full_aa,
        "constraint_H": H,
        "constraint_M": M,
        "constraint_C_Lambda": C,
        "constraint_determinant": det,
        **{f"regularity__{name}": values for name, values in regularity.items()},
    }
    require_finite_profiles(profiles)

    l2_closure = primary_l2_aa - reconstructed_l2
    source_closure = full_aa - (lapse_hessian_aa + ricci_tracefree_aa + primary_l3_aa)
    source_scale = np.abs(primary_l2_aa) + np.abs(lapse_hessian_aa) + np.abs(ricci_tracefree_aa)
    l2_bound = 64.0 * np.finfo(float).eps * np.maximum(source_scale, np.finfo(float).tiny)
    if np.any(np.abs(l2_closure) > l2_bound):
        raise RuntimeError("L2 Aa decomposition does not close within operand-scaled roundoff")
    if float(np.max(np.abs(source_closure))) > 1e-14:
        raise RuntimeError("full Aa source decomposition closure exceeded 1e-14")
    if digest_state(state) != before:
        raise RuntimeError("an analysis call mutated the initialized state")

    indices = range(min(6, n))
    center_rows = []
    for i in indices:
        center_rows.append({
            "i": int(i),
            "r": float(r[i]),
            "a_over_b_minus_one_over_r2": float((g.a[i] / g.b[i] - 1.0) / (r[i] * r[i])),
            "Aa_over_r2": float(g.Aa[i] / (r[i] * r[i])),
            "Lambda_over_r": float(g.Lambda[i] / r[i]),
            "X_prime_over_r": float(even_quotients["X_radial_derivative_over_r"][i]),
            "alpha_prime_over_r": float(even_quotients["alpha_radial_derivative_over_r"][i]),
            "K_prime_over_r": float(even_quotients["K_radial_derivative_over_r"][i]),
            "Rr_r_minus_R_over_3": float(ricci_tracefree_aa[i] / g.alpha[i]),
            "ricci_tracefree_Aa": float(ricci_tracefree_aa[i]),
            "lapse_hessian_Aa": float(lapse_hessian_aa[i]),
            "L2_Aa": float(primary_l2_aa[i]),
            "L3_Aa_total": float(primary_l3_aa[i]),
            "full_Aa_RHS": float(full_aa[i]),
            "full_Aa_RHS_over_r2": float(full_aa[i] / (r[i] * r[i])),
        })
    annuli = {}
    for lo, hi in ((0.0, 2.0), (2.0, 5.0), (5.0, 10.0), (10.0, 20.0)):
        mask = (r >= lo) & (r < hi)
        annuli[f"{lo:g}-{hi:g}"] = {
            "count": int(np.count_nonzero(mask)),
            "r_min": float(np.min(r[mask])) if np.any(mask) else None,
            "r_max": float(np.max(r[mask])) if np.any(mask) else None,
            "full_Aa_RHS_max_abs": float(np.max(np.abs(full_aa[mask]))) if np.any(mask) else None,
            "full_Aa_RHS_over_r2_max_abs": float(np.max(np.abs((full_aa / radius2)[mask]))) if np.any(mask) else None,
            "ricci_tracefree_Aa_max_abs": float(np.max(np.abs(ricci_tracefree_aa[mask]))) if np.any(mask) else None,
        }

    report = {
        "resolution": int(n),
        "r_max": R_MAX,
        "dr": float(grid.dr),
        "t": float(state.t),
        "tau": float(state.tau),
        "history_rows": int(len(state.history)),
        "state_digest_before": before,
        "state_digest_after": digest_state(state),
        "state_unchanged_by_observer": digest_state(state) == before,
        "initial_constraint_witnesses": {
            "H_cell_0": float(H[0]),
            "H_full_absmax": float(np.max(np.abs(H))),
            "M_cell_0": float(M[0]),
            "M_full_absmax": float(np.max(np.abs(M))),
            "C_Lambda_cell_0": float(C[0]),
            "C_Lambda_full_absmax": float(np.max(np.abs(C))),
            "determinant_cell_0": float(det[0]),
            "determinant_full_absmax": float(np.max(np.abs(det))),
        },
        "regularity_quotients": {
            name: stats(values, r)
            for name, values in regularity.items()
        },
        "source_decomposition": {
            "max_abs_L2_closure": float(np.max(np.abs(l2_closure))),
            "max_abs_full_Aa_closure": float(np.max(np.abs(source_closure))),
            "max_operand_scaled_L2_closure": float(np.max(np.abs(l2_closure) / np.maximum(l2_bound, np.finfo(float).tiny))),
            "initial_L2_Aa": stats(primary_l2_aa, r),
            "initial_L3_vacuum_Aa": stats(np.asarray(l3_vac["Aa"]), r),
            "initial_L3_matter_increment_Aa": stats(matter_l3_aa, r),
            "initial_full_Aa_RHS": stats(full_aa, r),
            "initial_full_Aa_RHS_over_r2": stats(full_aa / radius2, r),
            "initial_tracefree_Ricci_Aa": stats(ricci_tracefree_aa, r),
            "initial_lapse_Hessian_Aa": stats(lapse_hessian_aa, r),
        },
        "first_six_cells": center_rows,
        "fixed_physical_annuli": annuli,
    }
    return report, profiles


def main() -> None:
    outdir = Path("artifacts/initial-centre-series")
    outdir.mkdir(parents=True, exist_ok=True)
    kernel = V55ProductionKernel()
    reports, all_profiles = [], {}
    for n in RESOLUTIONS:
        report, profiles = audit_resolution(n, outdir)
        reports.append(report)
        all_profiles[f"N{n}"] = profiles
        for name, values in profiles.items():
            all_profiles[f"N{n}__{name}"] = values
    archive = outdir / "initial_centre_series_profiles.npz"
    np.savez_compressed(archive, **{
        key: values for key, values in all_profiles.items()
    })
    report = {
        "schema": "v55_initial_centre_series_capture_v1",
        "status": "completed_initialization_only",
        "diagnostic_only": True,
        "scientific_admission": False,
        "evolution_steps": 0,
        "source_baseline": SOURCE_BASELINE,
        "source_blobs": SOURCE_BLOBS,
        "vendor_commit": VENDOR_COMMIT,
        "initialization_parameters": {
            "resolutions": list(RESOLUTIONS), "r_max": R_MAX,
            "amplitude": 0.01, "width": 7.0, "D_amplitude": 1.0e-10,
            "include_radiation": True, "beta_dm": -0.04, "lambda_m": 2.0,
        },
        "resolutions": reports,
        "interpretation_limits": [
            "No time step, evolution, retry, or campaign was run; this captures only the initialized CMC-solved slices.",
            "Regularity quotients are diagnostics, not a pass/fail rule or a fitted physical coefficient.",
            "The first-cell Aa RHS is compared to r^2 as a centre-series witness; finite-resolution values alone do not identify a wrong continuum equation.",
            "A smooth manufactured conformally-flat control separately shows the pinned Ricci operator can reproduce O(r^2) behavior; this actual-state capture tests whether the initialized fields share the required centre-compatible profile.",
            "No production equation, vendor code, physical source, boundary, gauge, initializer, or campaign parameter was changed.",
        ],
    }
    json_path = outdir / "report.json"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    hashes = [f"{sha256_file(path)}  {path.name}\n" for path in (json_path, archive)]
    (outdir / "SHA256SUMS.txt").write_text("".join(hashes), encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "scientific_admission": report["scientific_admission"],
        "evolution_steps": report["evolution_steps"],
        "resolutions": [{
            "N": item["resolution"],
            "r0": item["dr"] / 2.0,
            "H0": item["initial_constraint_witnesses"]["H_cell_0"],
            "M0": item["initial_constraint_witnesses"]["M_cell_0"],
            "Rtf_Aa0": item["first_six_cells"][0]["ricci_tracefree_Aa"],
            "AaDot0": item["source_decomposition"]["initial_full_Aa_RHS"]["cell_0"],
            "AaDotOverR20": item["source_decomposition"]["initial_full_Aa_RHS_over_r2"]["cell_0"],
            "LambdaOverR0": item["regularity_quotients"]["Lambda_odd_quotient_Lambda/r"]["cell_0"],
            "AaOverR20": item["regularity_quotients"]["Aa_even_quotient_Aa/r2"]["cell_0"],
            "state_unchanged": item["state_unchanged_by_observer"],
        } for item in reports],
        "output_hashes": {path.name: sha256_file(path) for path in (json_path, archive)},
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
