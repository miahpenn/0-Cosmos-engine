#!/usr/bin/env python3
"""Zero-evolution audit of the native-cell Misner-Sharp radial identity.

The observer initializes the exact pinned V5.5 production machine, advances no
time steps, and compares the production derivative of the odd Misner-Sharp mass
with the matter/geometry identity

    D_r M_MS = 4 pi R^2 (rho_total D_r R + j_total Rdot / alpha).

It records signed native-cell values around areal R=10 and at the centre, plus
the same-state H/M/connection constraints. The independent centered fourth-
order stencil is only a derivative-representation control; it is not a
replacement evolution law or a physical conservation prescription.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine import production_kernel as kernel_module
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.run_production import assert_finite_state
from engine.stress_energy import assemble_total_stress_energy
from engine.worldtube import areal_radius_derivative

SOURCE_COMMIT = "6df02e1f91dd4aa414bd760b7ee0d7d86f1605e3"
SOURCE_BLOBS = {
    "engine/production_kernel.py": "2954e7a79656849b5f088053ea0f4ac97d4bb443",
    "engine/v55_initial.py": "4cfbfe4be6ac3ce6ed492f4e435186997d94116a",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
    "engine/stress_energy.py": "c31c647f74fab6893254b3dea1366eada15b3110",
    "engine/v55_pirk_adapter.py": "f671b838d9aefd8926702d9c2235209f52cd6815",
}
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (80, 120, 160, 240, 320)
R_MAX = 80.0
TARGET_AREAL_R = 10.0
OUT = Path("results/native-ms-gradient-identity")


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def verify_source_pins() -> dict:
    observed = {}
    for path, expected in SOURCE_BLOBS.items():
        source_blob = _git("rev-parse", f"{SOURCE_COMMIT}:{path}")
        head_blob = _git("rev-parse", f"HEAD:{path}")
        if source_blob != expected or head_blob != expected:
            raise RuntimeError(
                f"source blob mismatch for {path}: source={source_blob}, "
                f"head={head_blob}, expected={expected}"
            )
        observed[path] = {"source_blob": source_blob, "head_blob": head_blob}

    vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor pin mismatch: {vendor}")

    engine_diff = _git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--", "engine"
    )
    if engine_diff:
        raise RuntimeError("engine files changed from source baseline:\n" + engine_diff)
    forbidden_diff = _git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--",
        "vendor", ".gitmodules", "requirements.txt", "pyproject.toml",
    )
    if forbidden_diff:
        raise RuntimeError(
            "vendor/dependency/config drift is forbidden:\n" + forbidden_diff
        )

    return {
        "source_commit": SOURCE_COMMIT,
        "source_blobs": observed,
        "vendor_commit": vendor,
        "engine_diff_empty": True,
        "vendor_dependency_diff_empty": True,
        "diagnostic_head": _git("rev-parse", "HEAD"),
    }


def state_digest(state) -> str:
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
        values = np.ascontiguousarray(np.asarray(arrays[name], dtype=np.float64))
        digest.update(name.encode("utf-8") + b"\0")
        digest.update(repr(values.shape).encode("ascii") + b"\0")
        digest.update(values.tobytes())
    return digest.hexdigest()


def centered_fourth_derivative(values, dr: float) -> np.ndarray:
    """Independent interior centered fourth-order derivative on cell centres.

    The first and last two entries are NaN because this helper deliberately
    avoids reproducing boundary/ghost-parity rules; comparisons use only
    interior entries where the standard five-point stencil applies.
    """
    x = np.asarray(values, dtype=float)
    step = float(dr)
    if x.ndim != 1 or x.size < 5:
        raise ValueError("fourth-order derivative requires at least five samples")
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("dr must be finite and positive")
    out = np.full_like(x, np.nan)
    out[2:-2] = (
        x[:-4] - 8.0 * x[1:-3] + 8.0 * x[3:-1] - x[4:]
    ) / (12.0 * step)
    return out


def _finite_array(values, label: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(arr)):
        bad = int(np.flatnonzero(~np.isfinite(arr))[0])
        raise FloatingPointError(f"{label} has non-finite value at index {bad}")
    return arr


def _native_cell_row(
    i: int,
    r: np.ndarray,
    R: np.ndarray,
    Rr: np.ndarray,
    Rdot: np.ndarray,
    mass: np.ndarray,
    mass_r_prod: np.ndarray,
    mass_r_ind: np.ndarray,
    mass_r_stress: np.ndarray,
) -> dict:
    return {
        "i": int(i),
        "r_coordinate": float(r[i]),
        "R_areal": float(R[i]),
        "R_r": float(Rr[i]),
        "Rdot": float(Rdot[i]),
        "tube_coordinate_speed_at_fixed_R": float(-Rdot[i] / Rr[i]),
        "M_MS": float(mass[i]),
        "D_r_M_MS_production_odd_parity": float(mass_r_prod[i]),
        "D_r_M_MS_independent_centered4": (
            float(mass_r_ind[i]) if math.isfinite(float(mass_r_ind[i])) else None
        ),
        "D_r_M_MS_stress_identity": float(mass_r_stress[i]),
        "production_signed_identity_error": float(mass_r_prod[i] - mass_r_stress[i]),
        "independent_signed_identity_error": (
            float(mass_r_ind[i] - mass_r_stress[i])
            if math.isfinite(float(mass_r_ind[i])) else None
        ),
        "production_minus_independent_stencil": (
            float(mass_r_prod[i] - mass_r_ind[i])
            if math.isfinite(float(mass_r_ind[i])) else None
        ),
    }


def audit_resolution(n: int, kernel: V55ProductionKernel) -> dict:
    state = kernel.initialize(
        resolution=n,
        r_max=R_MAX,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    assert_finite_state(state)
    state.geometry.assert_finite_positive()
    digest_before = state_digest(state)

    grid, geom = state.grid, state.geometry
    r = _finite_array(grid.centers, f"N{n}.r")
    diagnostic = kernel.diagnostics(state, profiles=True)
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)

    H = _finite_array(
        np.asarray(raw["hamiltonian"], dtype=float)
        - 16.0 * math.pi * np.asarray(total.rho, dtype=float),
        f"N{n}.H",
    )
    M = _finite_array(
        np.asarray(raw["momentum"], dtype=float)
        - 8.0 * math.pi * np.asarray(total.j, dtype=float),
        f"N{n}.M",
    )
    C = _finite_array(raw["connection"], f"N{n}.C_lambda")
    determinant = _finite_array(raw["determinant"], f"N{n}.determinant")

    R = _finite_array(diagnostic["areal_radius"], f"N{n}.R")
    if not np.all(np.diff(R) > 0.0):
        raise RuntimeError(f"N={n}: areal radius is not strictly increasing")
    Rr = _finite_array(areal_radius_derivative(grid, geom), f"N{n}.Rr")
    mass = _finite_array(diagnostic["misner_sharp"], f"N{n}.M_MS")

    stage = adapter.geometry_stage_terms(
        grid,
        geom,
        state.scalars,
        state.matter,
        lambda_m=kernel_module.LAMBDA_M,
    )
    explicit = stage["explicit"]
    Rdot = _finite_array(
        R * (
            0.5 * np.asarray(explicit["b"], dtype=float) / np.asarray(geom.b, dtype=float)
            - np.asarray(explicit["X"], dtype=float) / np.asarray(geom.X, dtype=float)
        ),
        f"N{n}.Rdot",
    )
    rho = _finite_array(diagnostic["total_rho"], f"N{n}.rho_total")
    j = _finite_array(diagnostic["total_j"], f"N{n}.j_total")
    mass_r_stress = _finite_array(
        4.0 * math.pi * R**2 * (
            rho * Rr + j * Rdot / np.asarray(geom.alpha, dtype=float)
        ),
        f"N{n}.M_r_stress",
    )
    # Misner-Sharp mass is odd under radial reflection: M=R(1-chi)/2,
    # with odd R and even chi. Use parity=-1 for the production derivative.
    mass_r_prod = _finite_array(
        grid.cell_derivative_fourth(mass, parity=-1),
        f"N{n}.D_r_M_MS",
    )
    mass_r_ind = centered_fourth_derivative(mass, float(grid.dr))

    interior = slice(2, n - 2)
    stencil_delta = mass_r_prod[interior] - mass_r_ind[interior]
    stencil_scale = max(
        float(np.max(np.abs(mass_r_prod[interior]))),
        float(np.max(np.abs(mass_r_ind[interior]))),
        float(np.max(np.abs(mass[interior]))) / float(grid.dr),
        np.finfo(float).tiny,
    )
    stencil_delta_max = float(np.max(np.abs(stencil_delta)))
    stencil_roundoff_bound = 256.0 * np.finfo(float).eps * stencil_scale

    target_i = int(np.searchsorted(R, TARGET_AREAL_R, side="left"))
    if target_i <= 0 or target_i >= n:
        raise RuntimeError(f"N={n}: R=10 not bracketed by native areal-radius cells")
    lower = target_i - 1
    native_indices = sorted(
        set(i for i in (lower - 1, lower, target_i, target_i + 1) if 0 <= i < n)
    )
    annulus = (R >= 8.0) & (R <= 12.0)
    annulus_error = mass_r_prod - mass_r_stress
    l2 = vacuum.primary_l2_rhs(grid, geom)
    l3 = adapter.primary_l3_with_matter(grid, geom, state.scalars, state.matter)
    aa_rhs = _finite_array(
        np.asarray(l2["Aa"], dtype=float) + np.asarray(l3["Aa"], dtype=float),
        f"N{n}.Aa_rhs",
    )

    digest_after = state_digest(state)
    if digest_after != digest_before:
        raise RuntimeError(f"N={n}: observer mutated the initialized state")

    rows = {
        "native_cells_bracketing_R10": [
            _native_cell_row(
                i, r, R, Rr, Rdot, mass, mass_r_prod, mass_r_ind, mass_r_stress
            )
            for i in native_indices
        ],
        "first_five_center_cells": [
            _native_cell_row(
                i, r, R, Rr, Rdot, mass, mass_r_prod, mass_r_ind, mass_r_stress
            )
            for i in range(min(5, n))
        ],
    }
    return {
        "N": int(n),
        "dr": float(grid.dr),
        "coordinate_dt_if_evolved": float(0.03 * grid.dr),
        "coordinate_time": float(state.t),
        "time_steps_taken": 0,
        "radiation_enabled": True,
        "state_digest_before_diagnostics": digest_before,
        "state_digest_after_diagnostics": digest_after,
        "state_unchanged_by_diagnostics": True,
        "central_constraints": {
            "H_cell_0": float(H[0]),
            "M_cell_0": float(M[0]),
            "C_lambda_cell_0": float(C[0]),
            "determinant_cell_0": float(determinant[0]),
            "H_full_absmax": float(np.max(np.abs(H))),
            "M_full_absmax": float(np.max(np.abs(M))),
            "C_lambda_full_absmax": float(np.max(np.abs(C))),
        },
        "central_Aa_rhs_over_r0_squared": float(aa_rhs[0] / r[0]**2),
        "native_mass_gradient_identity": {
            "target_R": TARGET_AREAL_R,
            "bracketing_indices": [int(lower), int(target_i)],
            "bracketing_R": [float(R[lower]), float(R[target_i])],
            "annulus_R_8_to_12_signed_error_max_abs": float(
                np.max(np.abs(annulus_error[annulus]))
            ),
            "annulus_R_8_to_12_signed_error_rms": float(
                np.sqrt(np.mean(annulus_error[annulus]**2))
            ),
            # Report both native cells bracketing R=10. Do not use argmin:
            # when R=10 is exactly halfway between cells, last-bit rounding can
            # select different sides at different resolutions.
            "lower_bracket_coordinate_r": float(r[lower]),
            "upper_bracket_coordinate_r": float(r[target_i]),
            "signed_identity_error_lower_R10_bracket_cell": float(annulus_error[lower]),
            "signed_identity_error_upper_R10_bracket_cell": float(annulus_error[target_i]),
            "production_vs_independent_centered4": {
                "interior_cells": int(n - 4),
                "max_abs_difference": stencil_delta_max,
                "operand_scaled_roundoff_bound": stencil_roundoff_bound,
                "within_operand_scaled_roundoff": bool(
                    stencil_delta_max <= stencil_roundoff_bound
                ),
            },
            "native_cells": rows["native_cells_bracketing_R10"],
        },
        "center_mass_gradient_samples": rows["first_five_center_cells"],
    }


def main() -> int:
    provenance = verify_source_pins()
    kernel = V55ProductionKernel()
    records = [audit_resolution(n, kernel) for n in RESOLUTIONS]
    report = {
        "schema": "native_ms_gradient_identity_zero_evolution_v1",
        "status": "completed_zero_evolution_audit",
        "diagnostic_only": True,
        "scientific_admission": False,
        "source_pins": provenance,
        "method": {
            "evolution": "No time steps; each resolution is initialized from the exact current-main production initializer.",
            "mass_parity": "Odd. The production D_r M_MS operator is called with parity=-1.",
            "stress_identity": "4*pi*R^2*(rho_total*R_r + j_total*Rdot/alpha), from the current production diagnostic arrays and geometry RHS.",
            "independent_stencil": "Explicit centred fourth-order derivative on interior cell-centred grid points, compared with the production derivative to check operator representation.",
            "native_sampling": "Native cells immediately around areal R=10 and first five centre cells; no interpolation is used for the reported identity errors.",
            "unchanged": "No engine/vendor/source, physical equation, coefficient, initializer, boundary or campaign setting changes.",
        },
        "configuration": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "amplitude": 0.01,
            "width": 7.0,
            "D_amplitude": 1.0e-10,
            "include_radiation": True,
            "target_areal_radius": TARGET_AREAL_R,
        },
        "records": records,
        "interpretation_limits": [
            "This measures the initialized semi-discrete state only; it does not establish trajectory conservation or long-horizon stability.",
            "A matching centred-stencil implementation is a derivative-representation control, not an independent Einstein-equation derivation.",
            "The radial identity error and centre constraints are separate signed witnesses; correlation does not uniquely identify a faulty equation.",
            "Scientific admission remains false.",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    digest = hashlib.sha256(report_path.read_bytes()).hexdigest()
    (OUT / "SHA256SUMS.txt").write_text(
        f"{digest}  report.json\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": report["status"],
        "head": provenance["diagnostic_head"],
        "source_commit": SOURCE_COMMIT,
        "report_sha256": digest,
        "rows": [{
            "N": row["N"],
            "time_steps_taken": row["time_steps_taken"],
            "H0": row["central_constraints"]["H_cell_0"],
            "M0": row["central_constraints"]["M_cell_0"],
            "C0": row["central_constraints"]["C_lambda_cell_0"],
            "Aa_rhs_over_r0_squared": row["central_Aa_rhs_over_r0_squared"],
            "mass_gradient_error_R8_to_12_max_abs": row["native_mass_gradient_identity"]["annulus_R_8_to_12_signed_error_max_abs"],
            "mass_gradient_error_lower_R10_bracket": row["native_mass_gradient_identity"]["signed_identity_error_lower_R10_bracket_cell"],
            "mass_gradient_error_upper_R10_bracket": row["native_mass_gradient_identity"]["signed_identity_error_upper_R10_bracket_cell"],
            "stencil_delta_max_abs": row["native_mass_gradient_identity"]["production_vs_independent_centered4"]["max_abs_difference"],
        } for row in records],
    }, indent=2, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
