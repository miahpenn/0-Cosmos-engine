#!/usr/bin/env python3
"""Current-source, single-configuration, no-evolution origin-moment audit.

The production initializer and all vendor/runtime code are read-only. The sole
counterfactual changes the first half-cell moment on the final production
source profile while preserving all centre-to-centre trapezoids. A separate
continuous integral is admitted only when this initialization satisfies its
declared uncharged Gaussian source-family premises.
"""
from __future__ import annotations

import csv
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SOURCE_COMMIT = "ea2dabc17426b714eb445871c3c6b880638eee38"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTION = 80
R_MAX = 40.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
PHI0 = 0.179055
PI0 = 0.00341913
OUTPUT_DIR = Path("results/current-source-origin-initializer-2026-10-17")

SOURCE_BLOBS = {
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/production_kernel.py": "9685e76bb5a6f26d6bf383e240de24f37523f265",
    "engine/matter_system.py": "e52a117690253ddad7b24d116145a65474a351d2",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/v55_pirk_adapter.py": "eeee0492d80889d70e52b0946563eee5054c0f11",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
    "engine/native_geometry/__init__.py": "796d338d7cf2ed87bcb6ae63808dcf08110aeff8",
    "engine/native_geometry/operators.py": "2748be304b58328e2620f0a7fde49b01c0794ae2",
    "engine/native_geometry/geometry.py": "9e916b562a137ada0a8ba0268a5d2140162654d3",
    "engine/native_geometry/constraints.py": "124d1f237e54d082b4d76213f9ebe144e02b02ba",
    "engine/native_geometry/evolution.py": "7b47d700680b12ab3c82cb6a556d40580b103ce9",
    "tools/audit_current_main_centre_proper_radial_ricci.py": "da485d220664269cd52563728acfc9441678cda0",
    "tools/audit_continuous_origin_half_cell_oracle.py": "6e6f70cc4e57db67815f3cab41558d14a7c63797",
}
ROOT_BLOBS = {
    ".gitmodules": "502d49645859dc22c848efee0f654ed2170e18b9",
    "requirements.txt": "1e2a5bcbfd897f6465081a8794e14966f6e5db09",
}
VENDOR_BLOBS = {
    "r0_operators.py": "7f92fd9857981561dc2dff568addae87790c989f",
    "vacuum_bssn.py": "099ce52ba9e6dcd38d50c2915248c76c4236f7f4",
}


def _git(*args: str, cwd: Path | None = None) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=str(cwd or ROOT), text=True
    ).strip()


def verify_source_pins() -> dict[str, Any]:
    """Fail closed unless the current source and the untouched vendor are pinned."""
    head = _git("rev-parse", "HEAD")
    gitlink = _git("rev-parse", f"{head}:vendor/bb-palatini-unified-r0")
    if gitlink != VENDOR_COMMIT:
        raise RuntimeError(
            f"vendor gitlink mismatch: expected {VENDOR_COMMIT}, found {gitlink}"
        )

    source_rows: dict[str, Any] = {}
    for path, expected in {**SOURCE_BLOBS, **ROOT_BLOBS}.items():
        source_blob = _git("rev-parse", f"{SOURCE_COMMIT}:{path}")
        head_blob = _git("rev-parse", f"{head}:{path}")
        row = {
            "expected_blob": expected,
            "source_commit_blob": source_blob,
            "audit_head_blob": head_blob,
            "matches": source_blob == expected and head_blob == expected,
        }
        source_rows[path] = row
        if not row["matches"]:
            raise RuntimeError(f"source blob pin mismatch for {path}: {row}")

    vendor_root = ROOT / "vendor" / "bb-palatini-unified-r0"
    vendor_head = _git("rev-parse", "HEAD", cwd=vendor_root)
    if vendor_head != VENDOR_COMMIT:
        raise RuntimeError(
            f"checked-out vendor mismatch: expected {VENDOR_COMMIT}, found {vendor_head}"
        )
    vendor_rows: dict[str, Any] = {}
    for path, expected in VENDOR_BLOBS.items():
        actual = _git("rev-parse", f"HEAD:{path}", cwd=vendor_root)
        vendor_rows[path] = {"expected_blob": expected, "actual_blob": actual,
                             "matches": actual == expected}
        if actual != expected:
            raise RuntimeError(f"vendor file pin mismatch for {path}: {actual}")

    changed = _git(
        "diff", "--name-only", SOURCE_COMMIT, head, "--",
        "engine", "vendor", ".gitmodules", "requirements.txt",
    )
    if changed:
        raise RuntimeError(
            "unexpected engine/vendor/dependency changes relative to pinned source:\n"
            + changed
        )
    return {
        "source_commit": SOURCE_COMMIT,
        "audit_head": head,
        "vendor_commit": gitlink,
        "source_blobs": source_rows,
        "vendor_file_blobs": vendor_rows,
        "engine_vendor_dependency_diff_empty": True,
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _max_abs(values: Any) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.max(np.abs(array))) if array.size else 0.0


def _state_digest(state: Any) -> str:
    """Hash all initialized fields used here, plus optional charge sectors."""
    digest = hashlib.sha256()

    def add(label: str, values: Any) -> None:
        array = np.ascontiguousarray(np.asarray(values))
        digest.update(label.encode("utf-8") + b"\0")
        digest.update(array.dtype.str.encode("ascii") + b"\0")
        digest.update(repr(array.shape).encode("ascii") + b"\0")
        digest.update(array.tobytes(order="C"))

    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        add("geometry." + name, getattr(state.geometry, name))
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        add("scalars." + name, getattr(state.scalars, name))
    for species_name in ("dark_matter", "baryons", "radiation", "charged_dust"):
        species = getattr(state.matter, species_name, None)
        if species is None:
            digest.update(("matter." + species_name + ".NONE").encode("ascii"))
            continue
        for name in ("rest", "energy_t", "momentum_r"):
            add("matter." + species_name + "." + name, getattr(species, name))
    electromagnetic = getattr(state, "electromagnetic", None)
    if electromagnetic is None:
        digest.update(b"electromagnetic.NONE")
    else:
        digest.update(repr(type(electromagnetic).__qualname__).encode("utf-8"))
        for name, value in sorted(vars(electromagnetic).items()):
            add("electromagnetic." + name, value)
    add("grid.centers", state.grid.centers)
    digest.update(repr((state.grid.n, state.grid.r_max, state.grid.dr,
                        state.t, state.tau, state.e_folds)).encode("ascii"))
    return digest.hexdigest()


def _source_family_checks(state, r: np.ndarray) -> dict[str, Any]:
    """Check field-family premises; do not equate discrete derivatives to exact ones."""
    from engine.scalar_system import cosmos_potential
    from engine.v55_matter import RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT

    del cosmos_potential  # documented source-family constant; used by the initializer
    geometry, fields, matter = state.geometry, state.scalars, state.matter
    expected_s = AMPLITUDE * np.exp(-((r / WIDTH) ** 2))
    expected_d = D_AMPLITUDE * np.exp(-((r / WIDTH) ** 2))
    rho0 = 0.5 * PI0**2
    # The initializer's phi is constant; the potential is checked by the same
    # source function used by build_initial_data rather than inferred from output.
    from engine.scalar_system import cosmos_potential as potential
    rho0 += float(potential(np.array([PHI0]))[0])
    rho0 += float(RHO_DM_PRESENT + RHO_B_PRESENT + RHO_R_PRESENT)
    h0 = float(np.sqrt(max(rho0 / 3.0, 0.0)))
    kt = np.asarray(geometry.K, dtype=float) / 3.0 - np.asarray(geometry.Aa, dtype=float) / 2.0

    def profile_check(name: str, actual: Any, expected: Any, atol: float = 2.0e-14) -> dict[str, Any]:
        actual_array = np.asarray(actual, dtype=float)
        expected_array = np.asarray(expected, dtype=float)
        error = _max_abs(actual_array - expected_array)
        return {"passed": bool(actual_array.shape == expected_array.shape and error <= atol),
                "max_abs_difference": error, "absolute_check_tolerance": atol}

    checks = {
        "S_is_declared_gaussian": profile_check("S", fields.S, expected_s),
        "D_is_declared_gaussian": profile_check("D", fields.D, expected_d),
        "PS_is_zero": profile_check("PS", fields.PS, np.zeros_like(r)),
        "PD_equals_D": profile_check("PD", fields.PD, fields.D),
        "phi_is_declared_constant": profile_check("phi", fields.phi, np.full_like(r, PHI0)),
        "Pi_is_declared_constant": profile_check("Pi", fields.Pi, np.full_like(r, PI0)),
        "Kt_is_declared_constant": profile_check("K/3-Aa/2", kt, np.full_like(r, -h0)),
        "radiation_enabled": {"passed": matter.radiation is not None},
        "charged_dust_absent": {"passed": getattr(matter, "charged_dust", None) is None},
        "electromagnetic_sector_absent": {"passed": getattr(state, "electromagnetic", None) is None},
        "metric_is_default_areal_branch": profile_check(
            "R_areal-r", r * np.sqrt(np.asarray(geometry.b)) / np.asarray(geometry.X), r
        ),
    }
    return {
        "premises_satisfied": all(bool(check["passed"]) for check in checks.values()),
        "checks": checks,
        "H0_from_declared_background": h0,
        "interpretation": (
            "These checks establish the declared field family and absence of the charged/Maxwell sectors. "
            "The continuous formula is still a continuum comparator: production uses finite-grid "
            "derivatives for S', D' and Kt', so pointwise equality with the discrete source is not presumed."
        ),
    }


def _source_components(
    B: np.ndarray, *, r: np.ndarray, grid: Any,
    S: np.ndarray, PS: np.ndarray, D: np.ndarray, PD: np.ndarray,
    phi: np.ndarray, Pi: np.ndarray, H0: float,
) -> dict[str, np.ndarray]:
    """Replay the actual uncharged source at a supplied copied B profile."""
    from engine.scalar_system import cosmos_potential
    from engine.v55_matter import RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT

    Sp = np.asarray(grid.cell_derivative_fourth(S, parity=1), dtype=float)
    Dp = np.asarray(grid.cell_derivative_fourth(D, parity=1), dtype=float)
    Kt = np.full_like(r, -H0)
    Ktp = np.asarray(grid.cell_derivative_fourth(Kt, parity=1), dtype=float)
    current = -(PD * Dp)
    rho_s = 0.5 * (PS**2 + B * Sp**2) + 0.5 * S**2
    rho_d = 0.5 * (PD**2 + B * Dp**2) - 0.5 * D**2
    rho_phi = (0.5 * Pi**2 + cosmos_potential(phi)) / (8.0 * np.pi)
    rho_matter = np.full_like(
        r, (RHO_DM_PRESENT + RHO_B_PRESENT + RHO_R_PRESENT) / (8.0 * np.pi)
    )
    rho_total = rho_s + rho_d + rho_phi + rho_matter
    Kr = Kt + r * Ktp + 4.0 * np.pi * r * current
    source = (
        1.0 - 8.0 * np.pi * r**2 * rho_total
        + 2.0 * r**2 * Kr * Kt + r**2 * Kt**2
    )
    if not all(np.all(np.isfinite(x)) for x in (source, rho_total, Kr, Sp, Dp, Ktp)):
        raise FloatingPointError("source replay produced a non-finite profile")
    return {
        "source": source, "rho_total": rho_total, "rho_s": rho_s,
        "rho_d": rho_d, "rho_phi": rho_phi, "rho_matter": rho_matter,
        "Kr": Kr, "Sp_discrete": Sp, "Dp_discrete": Dp, "Kt_prime_discrete": Ktp,
    }


def _metric_from_B(geometry: Any, B: np.ndarray, grid: Any) -> Any:
    """Make an isolated geometry copy from a candidate B; preserve K and Aa."""
    from copy import deepcopy

    copied = deepcopy(geometry)
    A = 1.0 / np.sqrt(B)
    a, b, X = A ** (4.0 / 3.0), A ** (-2.0 / 3.0), A ** (-1.0 / 3.0)
    r = np.asarray(grid.centers, dtype=float)
    ap = np.asarray(grid.cell_derivative_fourth(a, parity=1), dtype=float)
    bp = np.asarray(grid.cell_derivative_fourth(b, parity=1), dtype=float)
    lam = ap / (2.0 * a**2) - bp / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)
    copied.a = a.copy()
    copied.b = b.copy()
    copied.X = X.copy()
    copied.Lambda = lam.copy()
    copied.B = 0.75 * lam
    return copied


def _metric_witnesses(
    geometry: Any, *, grid: Any, native_grid: Any, vacuum: Any,
    rho_total: np.ndarray, tag: str,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    from engine.native_geometry import native_bssn_constraints, native_bssn_geometry_terms
    from tools.audit_current_main_centre_proper_radial_ricci import proper_radial_ricci

    r = np.asarray(grid.centers, dtype=float)
    vendor_terms = vacuum.geometry_terms(grid, geometry)
    vendor_tf = np.asarray(vendor_terms["Rr_r"] - vendor_terms["R"] / 3.0, dtype=float)
    native_terms = native_bssn_geometry_terms(native_grid, geometry)
    native_tf = np.asarray(native_terms["Rr_r"] - native_terms["R"] / 3.0, dtype=float)
    areal = r * np.sqrt(np.asarray(geometry.b)) / np.asarray(geometry.X)
    proper_factor = np.sqrt(np.asarray(geometry.a)) / np.asarray(geometry.X)
    independent = proper_radial_ricci(grid, areal, proper_factor)
    proper_tf = np.asarray(independent["tf_ricci"], dtype=float)

    vendor_constraints = vacuum.constraints(grid, geometry)
    vendor_h = np.asarray(vendor_constraints["hamiltonian"], dtype=float) - 16.0 * np.pi * rho_total
    native_constraints = native_bssn_constraints(native_grid, geometry, rho=rho_total)
    native_h = np.asarray(native_constraints["H"], dtype=float)

    profiles = {
        f"{tag}_vendor_tf_ricci": vendor_tf,
        f"{tag}_native_tf_ricci": native_tf,
        f"{tag}_proper_radial_tf_ricci": proper_tf,
        f"{tag}_vendor_hamiltonian_full": vendor_h,
        f"{tag}_native_hamiltonian_full": native_h,
        f"{tag}_rho_total_recomputed": np.asarray(rho_total, dtype=float),
        f"{tag}_areal_radius": areal,
        f"{tag}_proper_factor": proper_factor,
        f"{tag}_connection_constraint_native": np.asarray(native_terms["connection_constraint"]),
        f"{tag}_determinant_constraint_native": np.asarray(native_constraints["determinant"]),
    }
    for name, values in profiles.items():
        if not np.all(np.isfinite(values)):
            raise FloatingPointError(f"non-finite witness profile {name}")

    cells = []
    for i in (0, 1):
        cells.append({
            "index": i, "r": float(r[i]),
            "vendor_tf_ricci": float(vendor_tf[i]),
            "native_tf_ricci": float(native_tf[i]),
            "proper_radial_tf_ricci": float(proper_tf[i]),
            "vendor_minus_proper_radial_tf_ricci": float(vendor_tf[i] - proper_tf[i]),
            "native_minus_proper_radial_tf_ricci": float(native_tf[i] - proper_tf[i]),
            "vendor_hamiltonian_full": float(vendor_h[i]),
            "native_hamiltonian_full": float(native_h[i]),
            "native_minus_vendor_hamiltonian": float(native_h[i] - vendor_h[i]),
            "rho_total_recomputed": float(rho_total[i]),
            "connection_constraint_native": float(native_terms["connection_constraint"][i]),
            "determinant_constraint_native": float(native_constraints["determinant"][i]),
        })
    return profiles, {
        "cells_0_and_1": cells,
        "interpretation": (
            "Signed same-configuration witnesses. The vendor and native H values both include "
            "the explicitly replayed 16*pi*rho source. They are diagnostics, not acceptance thresholds."
        ),
    }


def build_report(provenance: dict[str, Any] | None = None) -> tuple[dict[str, Any], dict[str, np.ndarray], list[dict[str, Any]]]:
    from engine.native_geometry import NativeSphericalCellGrid
    from engine.production_kernel import V55ProductionKernel
    from engine.v55_pirk_adapter import vendor_modules
    from tools.audit_continuous_origin_half_cell_oracle import (
        analytic_half_cell_envelope, continuous_reference,
    )

    grid_ops, vacuum, _ = vendor_modules()
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=RESOLUTION, r_max=R_MAX, amplitude=AMPLITUDE, width=WIDTH,
        D_amplitude=D_AMPLITUDE, include_radiation=True,
    )
    before = _state_digest(state)
    grid = state.grid
    native_grid = NativeSphericalCellGrid(RESOLUTION, R_MAX)
    r = np.asarray(grid.centers, dtype=float)
    dr = float(grid.dr)
    r0 = float(r[0])
    if not np.allclose(native_grid.centers, r, rtol=0.0, atol=1.0e-14):
        raise RuntimeError("native/reference grid centres differ for the fixed configuration")

    premise = _source_family_checks(state, r)
    if not premise["premises_satisfied"]:
        after = _state_digest(state)
        return ({
            "schema": "current_source_origin_initializer_audit_v1",
            "status": "UNKNOWN_CONTINUOUS_SOURCE_PREMISES_NOT_SATISFIED",
            "source_law_premises": premise,
            "configuration": {
                "resolution": RESOLUTION, "r_max": R_MAX, "dr": dr, "r0": r0,
                "amplitude": AMPLITUDE, "width": WIDTH, "D_amplitude": D_AMPLITUDE,
                "include_radiation": True,
            },
            "provenance": provenance,
            "input_state_sha256_before": before,
            "input_state_sha256_after": after,
            "input_state_unchanged": before == after,
            "evolution_steps": 0,
            "production_or_vendor_code_changed": False,
            "scientific_admission": False,
            "continuous_moment": None,
            "counterfactual": None,
            "interpretation": "Stopped fail-closed; no continuous moment was transplanted.",
        }, {}, [])

    fields = state.scalars
    S, PS, D, PD = (np.asarray(fields.S, dtype=float), np.asarray(fields.PS, dtype=float),
                    np.asarray(fields.D, dtype=float), np.asarray(fields.PD, dtype=float))
    phi, Pi = np.asarray(fields.phi, dtype=float), np.asarray(fields.Pi, dtype=float)
    H0 = float(premise["H0_from_declared_background"])

    # Exact replay of the six source-to-metric updates in build_initial_data.
    B_replay = np.ones_like(r)
    iteration_rows: list[dict[str, Any]] = []
    source_for_last_update: dict[str, np.ndarray] | None = None
    for iteration in range(1, 7):
        source_data = _source_components(
            B_replay, r=r, grid=grid, S=S, PS=PS, D=D, PD=PD, phi=phi, Pi=Pi, H0=H0
        )
        src = source_data["source"]
        moment = float(0.5 * dr * src[0])
        trapezoids = 0.5 * dr * (src[1:] + src[:-1])
        B_next = np.cumsum(np.r_[moment, trapezoids]) / r
        iteration_rows.append({
            "iteration": iteration, "applied_first_half_cell_moment": moment,
            "source_at_first_cell": float(src[0]),
            "max_abs_B_update": _max_abs(B_next - B_replay),
        })
        B_replay = B_next
        source_for_last_update = source_data

    if source_for_last_update is None:
        raise RuntimeError("the six-step initializer replay did not run")
    B_production = np.asarray(state.geometry.X, dtype=float) ** 6
    replay_error = _max_abs(B_replay - B_production)
    if replay_error > 5.0e-12:
        raise RuntimeError(
            f"current six-step source replay does not close to initialized X^6: {replay_error:.6e}"
        )

    # The integral is evaluated only after the field-family gate above passes.
    envelope = analytic_half_cell_envelope(r0, amplitude=AMPLITUDE, width=WIDTH,
                                           D_amplitude=D_AMPLITUDE)
    reference = continuous_reference(r0)
    I_cont = float(reference["reference_moment"])
    I_point = float(iteration_rows[-1]["applied_first_half_cell_moment"])
    bound_pad = 64.0 * np.finfo(float).eps * max(
        1.0, abs(float(envelope["A_integral_closed_form"]))
    )
    enclosed = (
        I_cont >= float(envelope["continuous_integral_lower_bound"]) - bound_pad
        and I_cont <= float(envelope["continuous_integral_upper_bound"]) + bound_pad
    )
    if not enclosed:
        raise RuntimeError("continuous integral estimate fell outside its model-conditional analytic enclosure")

    # One frozen-source control changes only the first moment. Every centre-to-centre
    # trapezoid remains byte-for-byte the value from the final production source profile.
    source_last = source_for_last_update["source"]
    B_counterfactual = np.cumsum(np.r_[
        I_cont, 0.5 * dr * (source_last[1:] + source_last[:-1])
    ]) / r

    geometry_baseline = copy.deepcopy(state.geometry)
    geometry_counterfactual = _metric_from_B(state.geometry, B_counterfactual, grid)
    source_baseline_now = _source_components(
        B_production, r=r, grid=grid, S=S, PS=PS, D=D, PD=PD, phi=phi, Pi=Pi, H0=H0
    )
    source_counterfactual_now = _source_components(
        B_counterfactual, r=r, grid=grid, S=S, PS=PS, D=D, PD=PD, phi=phi, Pi=Pi, H0=H0
    )
    baseline_profiles, baseline_witness = _metric_witnesses(
        geometry_baseline, grid=grid, native_grid=native_grid, vacuum=vacuum,
        rho_total=source_baseline_now["rho_total"], tag="production",
    )
    counter_profiles, counter_witness = _metric_witnesses(
        geometry_counterfactual, grid=grid, native_grid=native_grid, vacuum=vacuum,
        rho_total=source_counterfactual_now["rho_total"], tag="origin_moment_counterfactual",
    )

    def profile_summary(geometry: Any, B_metric: np.ndarray,
                        source_data: dict[str, np.ndarray], tag: str) -> tuple[dict[str, np.ndarray], list[dict[str, Any]]]:
        a, b, X = (np.asarray(geometry.a, dtype=float),
                   np.asarray(geometry.b, dtype=float), np.asarray(geometry.X, dtype=float))
        lam = np.asarray(geometry.Lambda, dtype=float)
        profiles = {
            f"{tag}_a": a, f"{tag}_b": b, f"{tag}_X": X, f"{tag}_B_metric": B_metric,
            f"{tag}_delta_over_r2": (1.0 / b - 1.0 / a) / (r * r),
            f"{tag}_Lambda_over_r": lam / r,
            f"{tag}_source_at_metric": source_data["source"],
        }
        rows = []
        for i in range(min(6, len(r))):
            rows.append({
                "index": i, "r": float(r[i]),
                "a": float(a[i]), "b": float(b[i]), "X": float(X[i]),
                "B_metric_X6": float(X[i]**6),
                "B_profile_used": float(B_metric[i]),
                "delta_over_r2": float((1.0 / b[i] - 1.0 / a[i]) / (r[i] * r[i])),
                "Lambda_over_r": float(lam[i] / r[i]),
                "source_recomputed_at_metric": float(source_data["source"][i]),
                "rho_recomputed_at_metric": float(source_data["rho_total"][i]),
            })
        return profiles, rows

    base_metric_profiles, base_rows = profile_summary(
        geometry_baseline, B_production, source_baseline_now, "production"
    )
    cf_metric_profiles, cf_rows = profile_summary(
        geometry_counterfactual, B_counterfactual, source_counterfactual_now,
        "origin_moment_counterfactual",
    )
    cells = []
    for base_row, cf_row in zip(base_rows, cf_rows):
        cells.append({
            "index": base_row["index"], "r": base_row["r"],
            **{f"production_{k}": v for k, v in base_row.items() if k not in ("index", "r")},
            **{f"counterfactual_{k}": v for k, v in cf_row.items() if k not in ("index", "r")},
        })

    arrays: dict[str, np.ndarray] = {"r": r.copy(), "B_production_replayed": B_production.copy(),
                                     "B_origin_moment_counterfactual": B_counterfactual.copy(),
                                     "source_last_production_update": source_last.copy()}
    arrays.update(base_metric_profiles)
    arrays.update(cf_metric_profiles)
    arrays.update(baseline_profiles)
    arrays.update(counter_profiles)
    for name, values in arrays.items():
        if not np.all(np.isfinite(np.asarray(values, dtype=float))):
            raise FloatingPointError(f"output array {name} contains a non-finite value")

    after = _state_digest(state)
    report = {
        "schema": "current_source_origin_initializer_audit_v1",
        "status": "COMPLETED_SOURCE_CURRENT_NO_EVOLUTION_AUDIT",
        "configuration": {
            "resolution": RESOLUTION, "r_max": R_MAX, "dr": dr, "r0": r0,
            "amplitude": AMPLITUDE, "width": WIDTH, "D_amplitude": D_AMPLITUDE,
            "include_radiation": True, "charged_dust_density": None,
            "radial_maxwell_charge_to_mass": None,
        },
        "provenance": provenance,
        "source_law_premises": premise,
        "production_replay": {
            "method": "exact six-iteration replay of current engine/v55_initial.py",
            "first_half_cell_rule": "I0 = 0.5*dr*source[0] = r0*source(r0)",
            "subsequent_interval_rule": "centre-to-centre trapezoids, 0.5*dr*(source[i]+source[i+1])",
            "iteration_history": iteration_rows,
            "replayed_B_vs_initialized_X6_max_abs": replay_error,
            "replay_closure_limit_from_existing_audit_contract": 5.0e-12,
            "replay_closed": bool(replay_error <= 5.0e-12),
            "applied_final_iteration_point_sample_moment": I_point,
            "source_at_first_cell_in_final_applied_update": float(source_last[0]),
            "final_production_source_first_six": [float(v) for v in source_last[:6]],
            "discrete_derivative_vs_continuum_field_derivative": {
                "S_prime_max_abs_difference": _max_abs(
                    source_last["Sp_discrete"] if isinstance(source_last, dict) else
                    grid.cell_derivative_fourth(S, parity=1) - (-2.0 * r * S / WIDTH**2)
                ),
            } if False else {
                "S_prime_max_abs_difference": _max_abs(
                    grid.cell_derivative_fourth(S, parity=1) - (-2.0 * r * S / WIDTH**2)
                ),
                "D_prime_max_abs_difference": _max_abs(
                    grid.cell_derivative_fourth(D, parity=1) - (-2.0 * r * D / WIDTH**2)
                ),
                "interpretation": "descriptive discretization differences; no tolerance or source correction inferred",
            },
        },
        "continuous_moment": {
            "method": "source-family-specific implicit-trapezoid solve with Richardson convergence estimate",
            "nominal_moment": I_cont,
            "point_sample_moment": I_point,
            "nominal_minus_point_sample": float(I_cont - I_point),
            "analytic_enclosure": {
                "lower": float(envelope["continuous_integral_lower_bound"]),
                "upper": float(envelope["continuous_integral_upper_bound"]),
                "width": float(envelope["continuous_integral_upper_bound"] - envelope["continuous_integral_lower_bound"]),
                "richardson_estimate_inside_enclosure": bool(enclosed),
                "floating_point_comparison_pad": bound_pad,
            },
            "reference": reference,
            "premise_assumptions": envelope["assumptions"],
            "interpretation_limit": (
                "The analytic enclosure is model-conditional and evaluated in float64, not directed-rounding interval arithmetic. "
                "The Richardson value is a convergence estimate. It is not assumed to equal the discrete initializer's true error."
            ),
        },
        "counterfactual": {
            "method": "frozen final production source; replace only I0 with the continuous nominal moment; preserve all centre-to-centre trapezoids",
            "origin_moment_used": I_cont,
            "changed_first_half_cell_moment_by": float(I_cont - I_point),
            "production_vs_counterfactual_B_first_six": [
                {"index": int(i), "B_production": float(B_production[i]),
                 "B_counterfactual": float(B_counterfactual[i]),
                 "signed_delta_B": float(B_counterfactual[i] - B_production[i])}
                for i in range(min(6, len(r))
            ],
            "signed_witnesses": {
                "production": baseline_witness,
                "origin_moment_counterfactual": counter_witness,
            },
            "interpretation_limit": (
                "This is a copied, frozen-source, one-moment intervention. Density is recomputed on each reported metric for H, "
                "but the counterfactual metric is not evolved or used to authorize a production integration change."
            ),
        },
        "input_state_sha256_before": before,
        "input_state_sha256_after": after,
        "input_state_unchanged": bool(before == after),
        "evolution_steps": 0,
        "production_or_vendor_code_changed": False,
        "scientific_admission": False,
        "interpretation_limits": [
            "One current-source configuration only: N=80, r_max=40, dr=0.5; no resolution sweep or evolution.",
            "The current production source replay is required to close against initialized X^6 before any counterfactual is reported.",
            "The continuous source model is a separate comparator, not an exact pointwise reconstruction of finite-grid derivative operations.",
            "Changing the first half-cell moment alone does not correct the unchanged centre-to-centre trapezoid error.",
            "Signed Ricci and Hamiltonian witnesses are descriptive diagnostics; no new acceptance threshold is introduced.",
        ],
    }
    if before != after:
        raise RuntimeError("the observer mutated the initialized state")
    return report, arrays, cells


def main() -> int:
    provenance = verify_source_pins()
    report, arrays, cells = build_report(provenance)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    report_path = OUTPUT_DIR / "audit.json"
    arrays_path = OUTPUT_DIR / "profiles.npz"
    csv_path = OUTPUT_DIR / "first_six.csv"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    np.savez_compressed(arrays_path, **arrays)
    headers = sorted({key for row in cells for key in row})
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=headers)
        writer.writeheader()
        writer.writerows(cells)
    manifest = {
        "schema": "current_source_origin_initializer_manifest_v1",
        "sha256": {
            "audit.json": _sha256(report_path),
            "profiles.npz": _sha256(arrays_path),
            "first_six.csv": _sha256(csv_path),
        },
        "source_commit": SOURCE_COMMIT,
        "audit_head": provenance["audit_head"],
        "vendor_commit": VENDOR_COMMIT,
        "evolution_steps": 0,
        "scientific_admission": False,
    }
    (OUTPUT_DIR / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": report["status"], "output_dir": str(OUTPUT_DIR),
        "input_state_unchanged": report["input_state_unchanged"],
        "evolution_steps": report["evolution_steps"],
        "replay_error": report.get("production_replay", {}).get("replayed_B_vs_initialized_X6_max_abs"),
        "point_sample_moment": report.get("production_replay", {}).get("applied_final_iteration_point_sample_moment"),
        "continuous_moment": report.get("continuous_moment", {}).get("nominal_moment")
        if report.get("continuous_moment") else None,
        "scientific_admission": False,
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
