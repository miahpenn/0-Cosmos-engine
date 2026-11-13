#!/usr/bin/env python3
"""Bounded, diagnostic-only signed momentum-constraint tangent audit.

The audit evaluates the current production semi-discrete RHS on newly
initialized frozen slices. It does not evolve a campaign. It re-solves the CMC
lapse for both directional perturbations and splits the momentum residual rate
into its geometric constraint term and six independently projected current
sectors. A single small ordinary step and half-step are used only as a local
tangent consistency check.

The analytic geometric chain rule cross-checks the finite-difference tangent.
A passing arithmetic closure is not physical/scientific admission.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
# In CI, import only from the separately checked-out historical source snapshot.
# Outside that workflow, default to the script's own checkout for importability.
_default_source_root = REPOSITORY_ROOT / "source"
PRODUCTION_SOURCE_ROOT = Path(
    os.environ.get(
        "COSMOS_PINNED_SOURCE_ROOT",
        str(_default_source_root if _default_source_root.is_dir() else REPOSITORY_ROOT),
    )
).resolve()
# Import the candidate implementation from this checkout. The separately
# pinned PRODUCTION_SOURCE_ROOT is retained for provenance checks only; putting
# it first would silently replace the candidate modules with historical code.
sys.path.insert(0, str(REPOSITORY_ROOT))

from engine.production_kernel import (
    BETA_DM,
    LAMBDA_M,
    ProductionState,
    V55ProductionKernel,
)
from engine.scalar_system import KAPPA, ScalarFields
from engine.matter_system import ConservedSpecies, Species, project_species
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter
from tools.momentum_constraint_chain_rule import analytic_geometric_momentum_tangent
from engine.native_geometry import native_bssn_constraints


SCHEMA = "current_source_signed_momentum_constraint_tangent_v1"
SOURCE_BASELINE = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
EXPECTED_VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_BLOBS = {
    "engine/production_kernel.py": "9685e76bb5a6f26d6bf383e240de24f37523f265",
    "engine/matter_system.py": "e52a117690253ddad7b24d116145a65474a351d2",
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
SCALAR_NAMES = ("S", "PS", "D", "PD", "phi", "Pi")
GEOMETRY_RATE_NAMES = ("a", "b", "X", "Aa", "K", "Lambda")
SECTOR_NAMES = ("S", "D", "COSMOS_phi", "dark_matter", "baryons", "radiation")
GEOMETRIC_M_NAMES = (
    "D_r_Aa",
    "minus_2third_D_r_K",
    "6_Aa_chi_r",
    "anisotropy_2_over_r",
    "anisotropy_b_gradient",
)
MACHINE_EPS = float(np.finfo(np.float64).eps)


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def provenance() -> dict[str, Any]:
    actual_blobs = {
        path: _git(
            "-C", str(PRODUCTION_SOURCE_ROOT), "rev-parse", f"HEAD:{path}"
        )
        for path in EXPECTED_BLOBS
    }
    wrong = {
        path: {"expected": expected, "actual": actual_blobs[path]}
        for path, expected in EXPECTED_BLOBS.items()
        if actual_blobs[path] != expected
    }
    if wrong:
        raise RuntimeError(
            "production source blob pin mismatch; refusing non-exact tangent "
            f"measurement: {json.dumps(wrong, sort_keys=True)}"
        )
    vendor = _git(
        "-C", str(PRODUCTION_SOURCE_ROOT / "vendor" / "bb-palatini-unified-r0"),
        "rev-parse", "HEAD",
    )
    if vendor != EXPECTED_VENDOR:
        raise RuntimeError(
            f"vendor pin mismatch: expected {EXPECTED_VENDOR}, got {vendor}"
        )
    return {
        "diagnostic_commit": _git(
            "-C", str(REPOSITORY_ROOT), "rev-parse", "HEAD"
        ),
        "production_source_commit": _git(
            "-C", str(PRODUCTION_SOURCE_ROOT), "rev-parse", "HEAD"
        ),
        "production_source_baseline_anchor": SOURCE_BASELINE,
        "production_source_root": str(PRODUCTION_SOURCE_ROOT),
        "production_source_blobs": actual_blobs,
        "vendor_commit": vendor,
    }


def assert_supported_baseline_source_scope(state: ProductionState) -> None:
    """Fail closed when the six-sector baseline audit would omit opt-in sources.

    The explicit J-piece inventory and state digest cover S, D, COSMOS phi,
    dark matter, baryons, and radiation only. Charged dust and radial Maxwell
    must be audited with their own current/stress witnesses; silently passing
    either state here would make the total-versus-sector closure incomplete.
    """
    if getattr(state.matter, "charged_dust", None) is not None:
        raise RuntimeError(
            "baseline six-sector tangent audit does not support charged_dust; "
            "use an audit that includes its projected J/rho/pr/pt and tangent"
        )
    if getattr(state, "electromagnetic", None) is not None:
        raise RuntimeError(
            "baseline six-sector tangent audit does not support radial EM; "
            "pass and audit the normalized EM stress/current projection explicitly"
        )


def state_digest(state: ProductionState) -> str:
    assert_supported_baseline_source_scope(state)
    h = hashlib.sha256()
    def add(name: str, values: Any) -> None:
        arr = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
        h.update(name.encode("utf-8"))
        h.update(str(arr.shape).encode("ascii"))
        h.update(arr.tobytes(order="C"))
    add("grid.centers", state.grid.centers)
    for name in ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B"):
        add(f"geometry.{name}", getattr(state.geometry, name))
    for name in SCALAR_NAMES:
        add(f"scalars.{name}", getattr(state.scalars, name))
    for species_name in ("dark_matter", "baryons", "radiation"):
        species = getattr(state.matter, species_name)
        for name in ("rest", "energy_t", "momentum_r"):
            add(f"matter.{species_name}.{name}", getattr(species, name))
    return h.hexdigest()


def finite(name: str, values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not bool(np.all(np.isfinite(arr))):
        raise FloatingPointError(f"non-finite values in {name}")
    return arr


def current_pieces(state: ProductionState) -> dict[str, np.ndarray]:
    """Rebuild each signed J contribution using production projection operators."""
    assert_supported_baseline_source_scope(state)
    assert_supported_baseline_source_scope(state)
    grid, geom, scalars, matter = (
        state.grid, state.geometry, state.scalars, state.matter
    )
    d1 = lambda f: grid.cell_derivative_fourth(f, parity=1)
    pieces = {
        "S": -np.asarray(scalars.PS) * d1(scalars.S),
        "D": -np.asarray(scalars.PD) * d1(scalars.D),
        "COSMOS_phi": (
            -np.asarray(scalars.Pi) * d1(scalars.phi) / KAPPA
        ),
    }
    metric = metric_slice_from_q(grid, geom)
    for key, state_name, species in (
        ("dark_matter", "dark_matter", Species.DARK_MATTER),
        ("baryons", "baryons", Species.BARYON),
        ("radiation", "radiation", Species.RADIATION),
    ):
        pieces[key] = np.asarray(
            project_species(metric, getattr(matter, state_name), species)["j"],
            dtype=np.float64,
        )
    return {key: finite(f"current {key}", value) for key, value in pieces.items()}


def _d1_stencil_operand_scale(grid, values: Any, parity: int = 1) -> np.ndarray:
    """Pointwise sum of absolute fourth-order D1 stencil operands.

    The derivative output can be tiny by cancellation even when the values
    entering its stencil are not. Roundoff budgets therefore use this operand
    scale, not abs(D1(values)) alone.
    """
    u = np.asarray(values, dtype=np.float64)
    if u.ndim != 1 or u.size < 5 or parity not in (-1, 1):
        raise ValueError("invalid D1 stencil-scale inputs")
    dr = float(grid.dr)
    if not math.isfinite(dr) or dr <= 0.0:
        raise FloatingPointError("D1 stencil-scale grid spacing must be positive")
    au = np.abs(u)
    n = u.size
    extended = np.empty(n + 4, dtype=np.float64)
    extended[2:-2] = au
    extended[1] = au[0]
    extended[0] = au[1]
    scale = np.empty_like(u)
    for i in range(n - 2):
        stencil = extended[i:i + 5]
        scale[i] = (
            stencil[0] + 8.0 * stencil[1]
            + 8.0 * stencil[3] + stencil[4]
        ) / (12.0 * dr)
    indices = np.arange(n - 5, n)
    for i in (n - 2, n - 1):
        offsets = indices.astype(float) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5, dtype=np.float64)
        target[1] = 1.0
        weights = np.linalg.solve(matrix, target)
        scale[i] = float(np.sum(np.abs(weights) * au[indices])) / dr
    return finite("D1 stencil operand scale", scale)


def _geometric_constraint_stencil_operand_scales_by_field(
    state: ProductionState,
) -> dict[str, np.ndarray]:
    """Expose raw D1 input scales separately for local failure diagnostics."""
    grid, geom = state.grid, state.geometry
    return {
        "Aa": _d1_stencil_operand_scale(grid, geom.Aa),
        "K": _d1_stencil_operand_scale(grid, geom.K),
        "b": _d1_stencil_operand_scale(grid, geom.b),
        "X": _d1_stencil_operand_scale(grid, geom.X),
    }


def _geometric_constraint_operand_scale(state: ProductionState) -> np.ndarray:
    """Estimate roundoff operands for evaluating M_geom at every cell."""
    grid, geom = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=np.float64)
    aa = np.asarray(geom.Aa, dtype=np.float64)
    kk = np.asarray(geom.K, dtype=np.float64)
    bb = np.asarray(geom.b, dtype=np.float64)
    xx = np.asarray(geom.X, dtype=np.float64)
    d_aa = grid.cell_derivative_fourth(aa, parity=1)
    d_kk = grid.cell_derivative_fourth(kk, parity=1)
    d_bb = grid.cell_derivative_fourth(bb, parity=1)
    d_xx = grid.cell_derivative_fourth(xx, parity=1)
    chi_r = -0.5 * d_xx / xx
    scale = (
        _d1_stencil_operand_scale(grid, aa)
        + (2.0 / 3.0) * _d1_stencil_operand_scale(grid, kk)
        + 3.0 * np.abs(aa / xx) * _d1_stencil_operand_scale(grid, xx)
        + 1.5 * np.abs(aa / bb) * _d1_stencil_operand_scale(grid, bb)
    )
    # Include post-stencil products/divisions and the sums in M_geom itself.
    scale = scale + (
        np.abs(d_aa)
        + (2.0 / 3.0) * np.abs(d_kk)
        + 6.0 * np.abs(aa * chi_r)
        + 1.5 * np.abs(aa * (2.0 / r + d_bb / bb))
    )
    return finite("geometric constraint roundoff operand scale", scale)


def geometric_momentum_pieces(state: ProductionState) -> dict[str, np.ndarray]:
    """Return signed terms from the exact pinned vendor momentum-constraint formula."""
    grid, geom = state.grid, state.geometry
    d1 = lambda f: grid.cell_derivative_fourth(f, parity=1)
    r = np.asarray(grid.centers, dtype=np.float64)
    Aa = np.asarray(geom.Aa, dtype=np.float64)
    Ab = -0.5 * Aa
    chip = -0.5 * d1(np.asarray(geom.X, dtype=np.float64)) / np.asarray(geom.X)
    bp = d1(np.asarray(geom.b, dtype=np.float64))
    return {
        "D_r_Aa": finite("M geometry D_r Aa", d1(Aa)),
        "minus_2third_D_r_K": finite(
            "M geometry -2/3 D_r K",
            -(2.0 / 3.0) * d1(np.asarray(geom.K, dtype=np.float64)),
        ),
        "6_Aa_chi_r": finite("M geometry 6 Aa chi_r", 6.0 * Aa * chip),
        "anisotropy_2_over_r": finite(
            "M geometry anisotropy 2/r", (Aa - Ab) * (2.0 / r),
        ),
        "anisotropy_b_gradient": finite(
            "M geometry anisotropy b gradient", (Aa - Ab) * bp / np.asarray(geom.b),
        ),
    }


def observe_constraint(state: ProductionState) -> dict[str, Any]:
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(state.grid, state.geometry)
    if "momentum" not in raw:
        raise KeyError("pinned vendor constraints output lacks 'momentum'")
    m_geom = finite("geometric momentum constraint", raw["momentum"])
    geom_pieces = geometric_momentum_pieces(state)
    geom_sum = np.zeros_like(m_geom)
    geom_operand_scale = np.abs(m_geom)
    for value in geom_pieces.values():
        geom_sum = geom_sum + value
        geom_operand_scale = geom_operand_scale + np.abs(value)
    geom_closure = m_geom - geom_sum
    geom_closure_bound = 64.0 * MACHINE_EPS * geom_operand_scale
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    j_total = finite("total projected current", total.j)
    pieces = current_pieces(state)
    j_sum = np.zeros_like(j_total)
    j_scale = np.abs(j_total)
    for value in pieces.values():
        j_sum = j_sum + value
        j_scale = j_scale + np.abs(value)
    j_residual = j_total - j_sum
    j_bound = 64.0 * MACHINE_EPS * j_scale
    m_total = m_geom - 8.0 * math.pi * j_total
    return {
        "M_geom": m_geom,
        "M_geom_pieces": geom_pieces,
        "M_geom_closure_abs": np.abs(geom_closure),
        "M_geom_closure_bound": geom_closure_bound,
        "M_geom_closure_ok": bool(np.all(np.abs(geom_closure) <= geom_closure_bound)),
        "J_total": j_total,
        "J_pieces": pieces,
        "M_total": finite("total momentum constraint", m_total),
        "J_closure_abs": np.abs(j_residual),
        "J_closure_bound": j_bound,
        "J_closure_ok": bool(np.all(np.abs(j_residual) <= j_bound)),
    }


def rhs_direction(
    kernel: V55ProductionKernel,
    state: ProductionState,
    *,
    geometry_terms_override: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate a signed source direction; optionally pin terms to the kernel's exact stage interface.

    The default preserves the historical initialized-slice audit. Accepted-state
    consumers should pass kernel._geometry_stage_terms(...) so native explicit
    a/b/X RHS overrides are not silently bypassed.
    """
    srhs, mrhs, _ = kernel._rhs(state)
    terms = (
        geometry_terms_override
        if geometry_terms_override is not None
        else adapter.geometry_stage_terms(
            state.grid, state.geometry, state.scalars, state.matter,
            lambda_m=LAMBDA_M,
        )
    )

    # Same-state decomposition of the vendor's primary L2 Aa rate. The
    # production L2 formula is the sum of the trace-free lapse Hessian and
    # lapse-weighted trace-free Ricci terms; the split is an accounting check,
    # not an independent alternative evolution operator.
    _, vacuum, _ = adapter.vendor_modules()
    vendor_geometry = vacuum.geometry_terms(state.grid, state.geometry)
    alpha = finite("Aa L2 split lapse", state.geometry.alpha)
    aa_lapse_hessian = -(
        finite("Aa L2 DrDr(alpha)", vendor_geometry["DrDr_alpha"])
        - finite("Aa L2 lap(alpha)", vendor_geometry["lap_alpha"]) / 3.0
    )
    aa_tracefree_ricci = alpha * (
        finite("Aa L2 Rr_r", vendor_geometry["Rr_r"])
        - finite("Aa L2 scalar R", vendor_geometry["R"]) / 3.0
    )
    aa_l2_operator = finite("Aa primary L2 operator", terms["primary_l2"]["Aa"])
    aa_l2_reconstructed = aa_lapse_hessian + aa_tracefree_ricci
    aa_l2_split_residual = aa_l2_operator - aa_l2_reconstructed
    aa_l2_split_scale = (
        np.abs(aa_l2_operator)
        + np.abs(aa_lapse_hessian)
        + np.abs(aa_tracefree_ricci)
    )
    aa_l2_split_bound = 512.0 * MACHINE_EPS * aa_l2_split_scale
    aa_l2_split = {
        "lapse_hessian": finite("Aa L2 lapse Hessian contribution", aa_lapse_hessian),
        "tracefree_ricci": finite("Aa L2 trace-free Ricci contribution", aa_tracefree_ricci),
        "operator": aa_l2_operator,
        "reconstructed": finite("Aa L2 reconstructed contribution", aa_l2_reconstructed),
        "residual": finite("Aa L2 split residual", aa_l2_split_residual),
        "closure": {
            "ok": bool(np.all(np.abs(aa_l2_split_residual) <= aa_l2_split_bound)),
            "max_abs": float(np.max(np.abs(aa_l2_split_residual))),
            "max_bound": float(np.max(aa_l2_split_bound)),
            "max_scaled_ratio": float(np.max(
                np.divide(
                    np.abs(aa_l2_split_residual),
                    aa_l2_split_bound,
                    out=np.zeros_like(aa_l2_split_residual),
                    where=aa_l2_split_bound > 0.0,
                )
            )),
        },
    }

    geom_rate: dict[str, np.ndarray] = {
        name: finite(f"explicit geometry rate {name}", terms["explicit"][name])
        for name in ("a", "b", "X")
    }
    for name in ("Aa", "K"):
        geom_rate[name] = finite(
            f"primary geometry rate {name}",
            np.asarray(terms["primary_l2"][name])
            + np.asarray(terms["primary_l3"][name]),
        )
    geom_rate["Lambda"] = finite(
        "Lambda geometry rate",
        np.asarray(terms["lambda_l2"]) + np.asarray(terms["lambda_l3"]),
    )
    scalar_rate = {
        name: finite(f"scalar rate {name}", getattr(srhs, name))
        for name in SCALAR_NAMES
    }
    matter_rate = {}
    for species_name in ("dark_matter", "baryons", "radiation"):
        r = mrhs[species_name]
        matter_rate[species_name] = {
            name: finite(
                f"matter rate {species_name}.{name}", getattr(r, name)
            )
            for name in ("rest", "energy_t", "momentum_r")
        }
    return {
        "geometry": geom_rate,
        "scalars": scalar_rate,
        "matter": matter_rate,
        "geometry_stage_terms": terms,
        "Aa_l2_split": aa_l2_split,
    }


def shifted_state(
    kernel: V55ProductionKernel,
    state: ProductionState,
    rates: dict[str, Any],
    delta: float,
) -> ProductionState:
    """Form a directional state +/- epsilon*RHS and re-solve the CMC lapse."""
    assert_supported_baseline_source_scope(state)
    assert_supported_baseline_source_scope(state)
    geom = state.geometry.copy()
    for name in GEOMETRY_RATE_NAMES:
        updated = (
            np.asarray(getattr(state.geometry, name), dtype=np.float64)
            + delta * rates["geometry"][name]
        )
        setattr(geom, name, finite(f"perturbed geometry {name}", updated))
    scalars = ScalarFields(
        *(
            finite(
                f"perturbed scalar {name}",
                getattr(state.scalars, name) + delta * rates["scalars"][name],
            )
            for name in SCALAR_NAMES
        )
    )
    species_out = {}
    for species_name in ("dark_matter", "baryons", "radiation"):
        base = getattr(state.matter, species_name)
        rate = rates["matter"][species_name]
        species_out[species_name] = ConservedSpecies(
            *(
                finite(
                    f"perturbed matter {species_name}.{name}",
                    getattr(base, name) + delta * rate[name],
                )
                for name in ("rest", "energy_t", "momentum_r")
            )
        )
    matter = V55MatterState(
        species_out["dark_matter"], species_out["baryons"],
        species_out["radiation"],
    )
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    geom.alpha = kernel._solve_lapse(
        state.grid, geom, scalars, matter
    )[0]
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    return ProductionState(
        grid=state.grid,
        geometry=geom,
        scalars=scalars,
        matter=matter,
        t=state.t,
        tau=state.tau,
        e_folds=state.e_folds,
    )


def _regions(r: np.ndarray) -> dict[str, np.ndarray]:
    i = np.arange(r.size)
    return {
        "cell_0": i == 0,
        "first_two_cells": i < 2,
        "first_three_cells": i < 3,
        "physical_annulus_2_le_r_lt_5": (r >= 2.0) & (r < 5.0),
        "physical_annulus_5_le_r_lt_10": (r >= 5.0) & (r < 10.0),
        "domain_r_ge_5": r >= 5.0,
    }


def summarize(values: Any, r: np.ndarray) -> dict[str, Any]:
    arr = finite("summary values", values)
    i = int(np.argmax(np.abs(arr)))
    return {
        "max_abs": float(abs(arr[i])),
        "signed_at_max": float(arr[i]),
        "index": i,
        "radius": float(r[i]),
        "rms": float(np.sqrt(np.mean(arr * arr))),
    }


def regional_summaries(values: Any, r: np.ndarray) -> dict[str, Any]:
    arr = finite("regional summary values", values)
    out = {}
    for label, mask in _regions(r).items():
        if bool(np.any(mask)):
            out[label] = summarize(arr[mask], r[mask])
        else:
            out[label] = None
    return out


def fixed_radius_samples(values: Any, r: np.ndarray) -> dict[str, float]:
    arr = finite("fixed-radius values", values)
    return {
        str(radius): float(np.interp(radius, r, arr))
        for radius in (2.0, 3.0, 5.0, 8.0)
        if float(r[0]) <= radius <= float(r[-1])
    }


def rhs_rate_witness(values: Any, r: np.ndarray) -> dict[str, Any]:
    """Compact cell-centred witnesses for an exact RHS rate profile."""
    arr = finite("RHS rate witness", values)
    return {
        "cell_0": float(arr[0]),
        "cell_1": float(arr[1]) if arr.size > 1 else None,
        "cell_2": float(arr[2]) if arr.size > 2 else None,
        "fixed_radius_samples": fixed_radius_samples(arr, r),
        "first_two_cells_max_abs": float(np.max(np.abs(arr[:2]))),
        "first_three_cells_max_abs": float(np.max(np.abs(arr[:3]))),
    }


def tangent_for_epsilon(
    kernel: V55ProductionKernel,
    state: ProductionState,
    rates: dict[str, Any],
    epsilon: float,
) -> dict[str, Any]:
    plus_state = shifted_state(kernel, state, rates, epsilon)
    minus_state = shifted_state(kernel, state, rates, -epsilon)
    plus = observe_constraint(plus_state)
    minus = observe_constraint(minus_state)
    denom = 2.0 * epsilon
    m_geom_dot = (plus["M_geom"] - minus["M_geom"]) / denom
    geom_term_dots = {
        key: (plus["M_geom_pieces"][key] - minus["M_geom_pieces"][key]) / denom
        for key in GEOMETRIC_M_NAMES
    }
    geom_term_sum = np.zeros_like(m_geom_dot)
    geom_term_sum_scale = np.zeros_like(m_geom_dot)
    geom_operand_scale = np.zeros_like(m_geom_dot)
    for key in GEOMETRIC_M_NAMES:
        geom_term_sum = geom_term_sum + geom_term_dots[key]
        geom_term_sum_scale = geom_term_sum_scale + np.abs(geom_term_dots[key])
        geom_operand_scale = geom_operand_scale + (
            np.abs(plus["M_geom_pieces"][key]) + np.abs(minus["M_geom_pieces"][key])
        )
    geom_tangent_closure = m_geom_dot - geom_term_sum
    geom_tangent_bound = (128.0 * MACHINE_EPS / denom) * geom_operand_scale
    geom_tangent_closure_ok = bool(
        np.all(np.abs(geom_tangent_closure) <= geom_tangent_bound)
    )
    j_dot = {
        key: (plus["J_pieces"][key] - minus["J_pieces"][key]) / denom
        for key in SECTOR_NAMES
    }
    sector_terms = {
        key: -8.0 * math.pi * value for key, value in j_dot.items()
    }
    m_dot = (plus["M_total"] - minus["M_total"]) / denom
    reconstructed = m_geom_dot.copy()
    for value in sector_terms.values():
        reconstructed = reconstructed + value
    closure = m_dot - reconstructed

    # Bound accounts for subtraction operands before differencing, not only
    # the usually much smaller derivative result.
    operand_scale = (
        np.abs(plus["M_total"]) + np.abs(minus["M_total"])
        + np.abs(plus["M_geom"]) + np.abs(minus["M_geom"])
    )
    for key in SECTOR_NAMES:
        operand_scale = operand_scale + 8.0 * math.pi * (
            np.abs(plus["J_pieces"][key]) + np.abs(minus["J_pieces"][key])
        )
    bound = (128.0 * MACHINE_EPS / denom) * operand_scale
    finite("momentum tangent", m_dot)
    finite("momentum tangent closure", closure)
    closure_ok = bool(np.all(np.abs(closure) <= bound))
    ratio = np.zeros_like(closure)
    positive = bound > 0.0
    ratio[positive] = np.abs(closure[positive]) / bound[positive]
    if bool(np.any((~positive) & (np.abs(closure) > 0.0))):
        max_ratio: float | None = None
    else:
        max_ratio = float(np.max(ratio))
    fd_stencil_scales_plus = _geometric_constraint_stencil_operand_scales_by_field(
        plus_state
    )
    fd_stencil_scales_minus = _geometric_constraint_stencil_operand_scales_by_field(
        minus_state
    )
    fd_stencil_scales_by_field = {
        name: finite(
            f"finite-difference {name} D1 stencil operand scale",
            fd_stencil_scales_plus[name] + fd_stencil_scales_minus[name],
        )
        for name in ("Aa", "K", "b", "X")
    }
    return {
        "epsilon": float(epsilon),
        "Mdot": m_dot,
        "Mdot_geom": m_geom_dot,
        "M_geom_fd_roundoff_operand_scale": finite(
            "finite-difference M_geom operand scale",
            _geometric_constraint_operand_scale(plus_state)
            + _geometric_constraint_operand_scale(minus_state),
        ),
        "M_geom_fd_stencil_operand_scales_by_field": fd_stencil_scales_by_field,
        "Mdot_geom_terms": geom_term_dots,
        "Mdot_geom_term_reconstructed": geom_term_sum,
        "Mdot_geom_tangent_closure_abs": np.abs(geom_tangent_closure),
        "Mdot_geom_tangent_closure_bound": geom_tangent_bound,
        "Mdot_geom_tangent_closure_ok": geom_tangent_closure_ok,
        "Mdot_geom_tangent_closure_max_abs": float(np.max(np.abs(geom_tangent_closure))),
        "Mdot_geom_tangent_closure_max_bound": float(np.max(geom_tangent_bound)),
        "Mdot_sector": sector_terms,
        "Mdot_reconstructed": reconstructed,
        "closure_abs": np.abs(closure),
        "closure_bound": bound,
        "closure_ok": closure_ok,
        "closure_max_abs": float(np.max(np.abs(closure))),
        "closure_max_bound": float(np.max(bound)),
        "closure_max_scaled_ratio": max_ratio,
        "plus_minus_current_projection_closure_ok": bool(
            plus["J_closure_ok"] and minus["J_closure_ok"]
        ),
        "plus_minus_current_projection_closure_max_abs": float(
            max(np.max(plus["J_closure_abs"]), np.max(minus["J_closure_abs"]))
        ),
    }


def compact_tangent_result(item: dict[str, Any], r: np.ndarray) -> dict[str, Any]:
    return {
        "epsilon": item["epsilon"],
        "Mdot": regional_summaries(item["Mdot"], r),
        "Mdot_fixed_radius_samples": fixed_radius_samples(item["Mdot"], r),
        "Mdot_geom": regional_summaries(item["Mdot_geom"], r),
        "Mdot_sector_contributions": {
            name: regional_summaries(value, r)
            for name, value in item["Mdot_sector"].items()
        },
        "Mdot_sector_fixed_radius_samples": {
            name: fixed_radius_samples(value, r)
            for name, value in item["Mdot_sector"].items()
        },
        "geometric_term_contributions": {
            name: regional_summaries(value, r)
            for name, value in item["Mdot_geom_terms"].items()
        },
        "geometric_term_fixed_radius_samples": {
            name: fixed_radius_samples(value, r)
            for name, value in item["Mdot_geom_terms"].items()
        },
        "geometric_term_tangent_closure": {
            "ok": item["Mdot_geom_tangent_closure_ok"],
            "max_abs": item["Mdot_geom_tangent_closure_max_abs"],
            "max_bound": item["Mdot_geom_tangent_closure_max_bound"],
        },
        "closure": {
            "ok": item["closure_ok"],
            "max_abs": item["closure_max_abs"],
            "max_bound": item["closure_max_bound"],
            "max_scaled_ratio": item["closure_max_scaled_ratio"],
        },
        "current_projection_closure": {
            "ok": item["plus_minus_current_projection_closure_ok"],
            "max_abs": item["plus_minus_current_projection_closure_max_abs"],
        },
    }


def one_step_check(
    kernel: V55ProductionKernel,
    state: ProductionState,
    m_initial: np.ndarray,
    tangent: np.ndarray,
    dt: float,
) -> dict[str, Any]:
    try:
        out = kernel.step(state, dt)
        m_after = observe_constraint(out)["M_total"]
        rate = (m_after - m_initial) / dt
        error = rate - tangent
        return {
            "completed": True,
            "dt": float(dt),
            "t_final": float(out.t),
            "M_change_over_dt": regional_summaries(rate, state.grid.centers),
            "difference_from_semi_discrete_tangent": regional_summaries(
                error, state.grid.centers
            ),
            "difference_fixed_radius_samples": fixed_radius_samples(
                error, state.grid.centers
            ),
        }
    except Exception as exc:
        return {
            "completed": False,
            "dt": float(dt),
            "failure_type": type(exc).__name__,
            "failure_message": str(exc),
        }


def audit_resolution(
    kernel: V55ProductionKernel,
    resolution: int,
    r_max: float,
    epsilons: tuple[float, float],
    dt: float,
    include_accepted_state: bool = False,
) -> dict[str, Any]:
    state = kernel.initialize(
        resolution=resolution, r_max=r_max, include_radiation=True
    )
    original_digest = state_digest(state)
    rates = rhs_direction(kernel, state)
    w0 = observe_constraint(state)
    r = np.asarray(state.grid.centers, dtype=np.float64)
    tangent_results = [
        tangent_for_epsilon(kernel, state, rates, eps)
        for eps in epsilons
    ]

    # Independent analytic chain rule for the exact pinned discrete M_geom.
    # This is deliberately separate from the continuum/Brown propagation oracle.
    chain_rule = analytic_geometric_momentum_tangent(
        lambda values: state.grid.cell_derivative_fourth(values, parity=1),
        r,
        state.geometry.Aa,
        state.geometry.K,
        state.geometry.b,
        state.geometry.X,
        rates["geometry"]["Aa"],
        rates["geometry"]["K"],
        rates["geometry"]["b"],
        rates["geometry"]["X"],
    )
    geom_fd_epsilon_delta = np.abs(
        tangent_results[0]["Mdot_geom"] - tangent_results[1]["Mdot_geom"]
    )
    analytic_stencil_operand_scales_by_field = {
        "Aa": _d1_stencil_operand_scale(state.grid, state.geometry.Aa),
        "K": _d1_stencil_operand_scale(state.grid, state.geometry.K),
        "b": _d1_stencil_operand_scale(state.grid, state.geometry.b),
        "X": _d1_stencil_operand_scale(state.grid, state.geometry.X),
        "dot_Aa": _d1_stencil_operand_scale(state.grid, rates["geometry"]["Aa"]),
        "dot_K": _d1_stencil_operand_scale(state.grid, rates["geometry"]["K"]),
        "dot_b": _d1_stencil_operand_scale(state.grid, rates["geometry"]["b"]),
        "dot_X": _d1_stencil_operand_scale(state.grid, rates["geometry"]["X"]),
    }
    chain_rule_checks = []
    for item in tangent_results:
        difference = finite(
            "analytic-vs-directional geometric tangent difference",
            item["Mdot_geom"] - chain_rule["total"],
        )
        # Bound D1 roundoff from the raw stencil operands, including cases
        # where D1(field) is near zero due to cancellation. Then add a separately
        # reported analytic-evaluation roundoff term and an epsilon-halving
        # truncation estimate. This envelope is for tangent accounting only.
        fd_roundoff_bound = (
            64.0 * MACHINE_EPS
            / (2.0 * float(item["epsilon"]))
            * item["M_geom_fd_roundoff_operand_scale"]
        )
        analytic_operand_scale = (
            _d1_stencil_operand_scale(state.grid, rates["geometry"]["Aa"])
            + (2.0 / 3.0) * _d1_stencil_operand_scale(state.grid, rates["geometry"]["K"])
            + 3.0 * np.abs(state.geometry.Aa / state.geometry.X)
                * _d1_stencil_operand_scale(state.grid, rates["geometry"]["X"])
            + 1.5 * np.abs(state.geometry.Aa / state.geometry.b)
                * _d1_stencil_operand_scale(state.grid, rates["geometry"]["b"])
        )
        # Extra sensitivity terms for D1(X) and D1(b) inside the product rule.
        dot_aa = np.asarray(rates["geometry"]["Aa"], dtype=np.float64)
        dot_bb = np.asarray(rates["geometry"]["b"], dtype=np.float64)
        dot_xx = np.asarray(rates["geometry"]["X"], dtype=np.float64)
        aa = np.asarray(state.geometry.Aa, dtype=np.float64)
        bb = np.asarray(state.geometry.b, dtype=np.float64)
        xx = np.asarray(state.geometry.X, dtype=np.float64)
        analytic_operand_scale = analytic_operand_scale + (
            3.0 * np.abs(dot_aa / xx) * _d1_stencil_operand_scale(state.grid, state.geometry.X)
            + 3.0 * np.abs(aa * dot_xx / (xx * xx))
                * _d1_stencil_operand_scale(state.grid, state.geometry.X)
            + 1.5 * np.abs(dot_aa / bb) * _d1_stencil_operand_scale(state.grid, state.geometry.b)
            + 1.5 * np.abs(aa * dot_bb / (bb * bb))
                * _d1_stencil_operand_scale(state.grid, state.geometry.b)
            + sum(np.abs(value) for value in chain_rule["terms"].values())
        )
        analytic_roundoff_bound = 64.0 * MACHINE_EPS * analytic_operand_scale
        truncation_bound = 4.0 * geom_fd_epsilon_delta
        fd_stencil_scales_by_field = item[
            "M_geom_fd_stencil_operand_scales_by_field"
        ]
        bound = finite(
            "analytic tangent cross-check bound",
            fd_roundoff_bound + analytic_roundoff_bound + truncation_bound,
        )
        scaled = np.zeros_like(difference)
        positive_bound = bound > 0.0
        scaled[positive_bound] = (
            np.abs(difference[positive_bound]) / bound[positive_bound]
        )
        zero_bound_mismatch = (~positive_bound) & (np.abs(difference) > 0.0)
        scaled[zero_bound_mismatch] = np.inf
        worst_index = int(np.argmax(scaled))
        chain_rule_checks.append({
            "epsilon": float(item["epsilon"]),
            "max_abs_difference": float(np.max(np.abs(difference))),
            "max_bound": float(np.max(bound)),
            "max_scaled_ratio": float(np.max(scaled)),
            "worst_scaled_ratio": float(scaled[worst_index]),
            "worst_index": worst_index,
            "worst_radius": float(r[worst_index]),
            "worst_difference": float(difference[worst_index]),
            "worst_bound": float(bound[worst_index]),
            "worst_fd_roundoff_bound": float(fd_roundoff_bound[worst_index]),
            "worst_analytic_roundoff_bound": float(analytic_roundoff_bound[worst_index]),
            "worst_epsilon_sensitivity_bound": float(truncation_bound[worst_index]),
            "worst_fd_total_operand_scale": float(
                item["M_geom_fd_roundoff_operand_scale"][worst_index]
            ),
            "worst_fd_stencil_operand_scales_by_field": {
                name: float(values[worst_index])
                for name, values in fd_stencil_scales_by_field.items()
            },
            "worst_analytic_stencil_operand_scales_by_field": {
                name: float(values[worst_index])
                for name, values in analytic_stencil_operand_scales_by_field.items()
            },
            "worst_analytic_total_operand_scale": float(analytic_operand_scale[worst_index]),
            "ok": bool(np.all(np.abs(difference) <= bound)),
            "_difference": difference,
            "_bound": bound,
        })

    # Negative control: deliberately flip the sign of the largest analytic
    # signed term. Both measured epsilon tangents must reject that mutation.
    mutated_term = max(
        chain_rule["terms"],
        key=lambda name: float(np.max(np.abs(chain_rule["terms"][name]))),
    )
    mutated_terms = dict(chain_rule["terms"])
    mutated_terms[mutated_term] = -mutated_terms[mutated_term]
    mutated_prediction = np.zeros_like(chain_rule["total"])
    for values in mutated_terms.values():
        mutated_prediction = mutated_prediction + values
    negative_control_records = []
    for item, check in zip(tangent_results, chain_rule_checks):
        mutation_difference = finite(
            "signed-term negative-control difference",
            item["Mdot_geom"] - mutated_prediction,
        )
        detected = bool(np.any(np.abs(mutation_difference) > check["_bound"]))
        negative_control_records.append({
            "epsilon": float(item["epsilon"]),
            "mutated_term": mutated_term,
            "max_abs_difference": float(np.max(np.abs(mutation_difference))),
            "max_bound": float(np.max(check["_bound"])),
            "detected": detected,
        })
    chain_rule_ok = all(check["ok"] for check in chain_rule_checks)
    negative_control_ok = all(record["detected"] for record in negative_control_records)
    chain_rule_report = {
        "basis": (
            "analytic derivative of the pinned discrete M_geom definition, using "
            "the same fixed parity-aware D1 for fields and RHS rates"
        ),
        "Mdot_geom": regional_summaries(chain_rule["total"], r),
        "terms": {
            name: rhs_rate_witness(values, r)
            for name, values in chain_rule["terms"].items()
        },
        "chi_r": rhs_rate_witness(chain_rule["chi_r"], r),
        "dot_chi_r": rhs_rate_witness(chain_rule["dot_chi_r"], r),
        "finite_difference_crosscheck": [
            {key: value for key, value in record.items() if not key.startswith("_")}
            for record in chain_rule_checks
        ],
        "crosscheck_envelope": {
            "finite_difference_roundoff_bound": "64 machine eps times plus/minus M_geom raw D1 stencil operand scales and post-stencil operation scales, divided by the two-sided epsilon",
            "analytic_roundoff_bound": "64 machine eps times the analytic tangent term and D1 operand scale",
            "finite_difference_truncation_bound": "4x pointwise difference between the two directional-epsilon Mdot_geom estimates",
            "scope_limit": "this envelope checks the algebraic tangent only; it is not a continuum propagation error bound"
        },
        "signed_term_negative_control": {
            "mutated_term": mutated_term,
            "mutation": "flip the selected analytic signed term",
            "per_epsilon": negative_control_records,
            "detected_at_both_epsilons": bool(negative_control_ok),
        },
        "arithmetic_gate_passed": bool(chain_rule_ok and negative_control_ok),
    }
    # Do not leak arrays through the report structure.
    for record in chain_rule_checks:
        record.pop("_difference", None)
        record.pop("_bound", None)

    e0, e1 = tangent_results
    eps_difference = e0["Mdot"] - e1["Mdot"]
    step_full = one_step_check(kernel, state, w0["M_total"], e0["Mdot"], dt)
    step_half = one_step_check(kernel, state, w0["M_total"], e0["Mdot"], 0.5 * dt)
    if state_digest(state) != original_digest:
        raise RuntimeError("frozen initial ProductionState was mutated by observer")

    # One accepted production step only, on the smallest requested resolution.
    # Reuse the existing signed tangent machinery at that exact accepted state;
    # this is not a campaign or a new physical-law/closure assumption.
    accepted_state_report = None
    if include_accepted_state:
        accepted_state = kernel.step(state, dt)
        accepted_state_digest = state_digest(accepted_state)
        accepted_rates = rhs_direction(kernel, accepted_state)
        accepted_witness = observe_constraint(accepted_state)
        accepted_tangents = [
            tangent_for_epsilon(kernel, accepted_state, accepted_rates, eps)
            for eps in epsilons
        ]
        accepted_r = np.asarray(accepted_state.grid.centers, dtype=np.float64)
        accepted_geom = accepted_state.geometry
        native_constraint = native_bssn_constraints(
            accepted_state.grid, accepted_geom
        )
        native_m_geom = finite(
            "accepted-state native geometric momentum",
            native_constraint["M_geometry"],
        )
        vendor_m_geom = finite(
            "accepted-state vendor geometric momentum",
            accepted_witness["M_geom"],
        )
        native_vendor_delta = finite(
            "accepted-state native minus vendor M geometry",
            native_m_geom - vendor_m_geom,
        )
        accepted_aa = finite("accepted-state Aa", accepted_geom.Aa)
        accepted_b = finite("accepted-state b", accepted_geom.b)
        accepted_x = finite("accepted-state X", accepted_geom.X)
        accepted_areal_radius = finite(
            "accepted-state areal radius",
            accepted_r * np.sqrt(accepted_b) / accepted_x,
        )
        accepted_d1 = lambda values, parity=1: finite(
            "accepted-state D1",
            accepted_state.grid.cell_derivative_fourth(values, parity=parity),
        )
        native_log_derivative = (
            accepted_d1(accepted_areal_radius, -1) / accepted_areal_radius
        )
        vendor_log_derivative = (
            1.0 / accepted_r
            + 0.5 * accepted_d1(accepted_b) / accepted_b
            - accepted_d1(accepted_x) / accepted_x
        )
        predicted_remainder = finite(
            "accepted-state signed product quotient remainder",
            3.0 * accepted_aa * (native_log_derivative - vendor_log_derivative),
        )
        remainder_closure = finite(
            "accepted-state remainder closure",
            native_vendor_delta - predicted_remainder,
        )
        vendor_term_operand_scale = np.zeros_like(vendor_m_geom)
        for values in accepted_witness["M_geom_pieces"].values():
            vendor_term_operand_scale = vendor_term_operand_scale + np.abs(values)
        native_term_operand_scale = (
            np.abs(accepted_d1(accepted_geom.Aa))
            + (2.0 / 3.0) * np.abs(accepted_d1(accepted_geom.K))
            + 2.0 * np.abs(accepted_geom.Aa - (-0.5 * accepted_geom.Aa))
                * np.abs(native_log_derivative)
        )
        # Bound floating-point cancellation from the actual formula operands,
        # not only from the small post-cancellation residuals.
        closure_scale = (
            np.abs(native_m_geom) + np.abs(vendor_m_geom)
            + np.abs(native_vendor_delta) + np.abs(predicted_remainder)
            + vendor_term_operand_scale + native_term_operand_scale
        )
        closure_bound = 128.0 * MACHINE_EPS * closure_scale
        remainder_closure_ok = bool(
            np.all(np.abs(remainder_closure) <= closure_bound)
        )
        accepted_digest_after = state_digest(accepted_state)
        if accepted_state_digest != accepted_digest_after:
            raise RuntimeError("accepted ProductionState was mutated by observer")
        accepted_arithmetic_ok = all(
            item["closure_ok"]
            and item["plus_minus_current_projection_closure_ok"]
            and item["Mdot_geom_tangent_closure_ok"]
            for item in accepted_tangents
        ) and bool(
            accepted_witness["J_closure_ok"]
            and accepted_witness["M_geom_closure_ok"]
            and accepted_rates["Aa_l2_split"]["closure"]["ok"]
            and remainder_closure_ok
        )
        accepted_state_report = {
            "scope": "one ordinary accepted production step from the initialized state; only the smallest requested resolution",
            "accepted_time": float(accepted_state.t),
            "accepted_state_sha256_before_observation": accepted_state_digest,
            "accepted_state_sha256_after_observation": accepted_digest_after,
            "accepted_state_unchanged": accepted_state_digest == accepted_digest_after,
            "momentum_constraint_total": regional_summaries(
                accepted_witness["M_total"], accepted_r
            ),
            "native_vendor_geometry_crosswalk": {
                "native_M_geometry": regional_summaries(native_m_geom, accepted_r),
                "vendor_M_geometry": regional_summaries(vendor_m_geom, accepted_r),
                "native_minus_vendor": regional_summaries(native_vendor_delta, accepted_r),
                "native_minus_vendor_fixed_radius_samples": fixed_radius_samples(
                    native_vendor_delta, accepted_r
                ),
                "predicted_signed_product_quotient_remainder": regional_summaries(
                    predicted_remainder, accepted_r
                ),
                "remainder_closure": {
                    "ok": remainder_closure_ok,
                    "max_abs": float(np.max(np.abs(remainder_closure))),
                    "max_bound": float(np.max(closure_bound)),
                    "rms": float(np.sqrt(np.mean(remainder_closure**2))),
                },
                "parity": {
                    "native_D_areal_radius": -1,
                    "vendor_D_b_X": 1,
                },
                "source_terms_excluded_from_crosswalk": True,
            },
            "tangent_by_epsilon": [
                compact_tangent_result(item, accepted_r)
                for item in accepted_tangents
            ],
            "geometric_term_signed_tangent_closure": [
                {
                    "epsilon": float(item["epsilon"]),
                    "ok": bool(item["Mdot_geom_tangent_closure_ok"]),
                    "max_abs": float(item["Mdot_geom_tangent_closure_max_abs"]),
                    "max_bound": float(item["Mdot_geom_tangent_closure_max_bound"]),
                }
                for item in accepted_tangents
            ],
            "current_projection_closure": {
                "ok": bool(accepted_witness["J_closure_ok"]),
                "max_abs": float(np.max(accepted_witness["J_closure_abs"])),
                "max_bound": float(np.max(accepted_witness["J_closure_bound"])),
            },
            "arithmetic_gate_passed": accepted_arithmetic_ok,
        }

    arithmetic_ok = all(
        item["closure_ok"]
        and item["plus_minus_current_projection_closure_ok"]
        and item["Mdot_geom_tangent_closure_ok"]
        for item in tangent_results
    ) and bool(
        w0["J_closure_ok"]
        and w0["M_geom_closure_ok"]
        and rates["Aa_l2_split"]["closure"]["ok"]
        and chain_rule_ok
        and negative_control_ok
        and (
            accepted_state_report is None
            or accepted_state_report["arithmetic_gate_passed"]
        )
    )
    return {
        "resolution": int(resolution),
        "r_max": float(r_max),
        "dr": float(state.grid.dr),
        "initial_state_sha256": original_digest,
        "initial_rhs_rate_samples": {
            "geometry_full_rhs": {
                name: rhs_rate_witness(rates["geometry"][name], r)
                for name in GEOMETRY_RATE_NAMES
            },
            "Aa_primary_l2": rhs_rate_witness(
                rates["geometry_stage_terms"]["primary_l2"]["Aa"], r
            ),
            "Aa_primary_l2_split": {
                "lapse_hessian": rhs_rate_witness(rates["Aa_l2_split"]["lapse_hessian"], r),
                "tracefree_ricci": rhs_rate_witness(rates["Aa_l2_split"]["tracefree_ricci"], r),
                "operator": rhs_rate_witness(rates["Aa_l2_split"]["operator"], r),
                "reconstructed": rhs_rate_witness(rates["Aa_l2_split"]["reconstructed"], r),
                "residual": rhs_rate_witness(rates["Aa_l2_split"]["residual"], r),
                "closure": rates["Aa_l2_split"]["closure"],
            },
            "Aa_primary_l3_including_matter": rhs_rate_witness(
                rates["geometry_stage_terms"]["primary_l3"]["Aa"], r
            ),
            "Aa_full_over_r_squared": rhs_rate_witness(
                rates["geometry"]["Aa"] / (r * r), r
            ),
            "K_primary_l2": rhs_rate_witness(
                rates["geometry_stage_terms"]["primary_l2"]["K"], r
            ),
            "K_primary_l3_including_matter": rhs_rate_witness(
                rates["geometry_stage_terms"]["primary_l3"]["K"], r
            ),
        },
        "initial_momentum_constraint": regional_summaries(w0["M_total"], r),
        "initial_geometric_momentum_term_closure": {
            "ok": w0["M_geom_closure_ok"],
            "max_abs": float(np.max(w0["M_geom_closure_abs"])),
            "max_bound": float(np.max(w0["M_geom_closure_bound"])),
        },
        "initial_current_projection_closure": {
            "ok": w0["J_closure_ok"],
            "max_abs": float(np.max(w0["J_closure_abs"])),
            "max_bound": float(np.max(w0["J_closure_bound"])),
        },
        "semi_discrete_chain_rule": chain_rule_report,
        "tangent_by_epsilon": [
            compact_tangent_result(item, r) for item in tangent_results
        ],
        "epsilon_sensitivity": {
            "max_abs_Mdot_difference": float(np.max(np.abs(eps_difference))),
            "regional_Mdot_difference": regional_summaries(eps_difference, r),
            "fixed_radius_samples": fixed_radius_samples(eps_difference, r),
        },
        "bounded_one_step_consistency": {
            "same_common_dt_all_resolutions": True,
            "dt": step_full,
            "dt_half": step_half,
            "interpretation": (
                "local tangent consistency only; a finite-step discrepancy "
                "is not by itself a continuum or physics verdict"
            ),
        },
        "accepted_state_one_step_observer": accepted_state_report,
        "arithmetic_gate_passed": bool(arithmetic_ok),
        "scientific_admission": False,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    prov = provenance()
    if len(args.epsilons) != 2 or any(e <= 0.0 for e in args.epsilons):
        raise ValueError("provide exactly two positive directional epsilons")
    if not args.resolutions or any(n < 4 for n in args.resolutions):
        raise ValueError("resolutions must each be at least four")
    if not math.isfinite(args.dt) or args.dt <= 0.0:
        raise ValueError("common dt must be finite and positive")
    kernel = V55ProductionKernel()
    per_resolution = [
        audit_resolution(
            kernel, n, args.r_max,
            (float(args.epsilons[0]), float(args.epsilons[1])),
            float(args.dt),
            include_accepted_state=(n == min(args.resolutions)),
        )
        for n in args.resolutions
    ]
    arithmetic_ok = all(item["arithmetic_gate_passed"] for item in per_resolution)
    return {
        "schema": SCHEMA,
        "status": "completed_diagnostic_only" if arithmetic_ok else "arithmetic_gate_failed",
        "provenance": prov,
        "initialization": {
            "r_max": float(args.r_max),
            "amplitude": 0.01,
            "width": 7.0,
            "D_amplitude": 1.0e-10,
            "radiation": True,
            "beta_dm": BETA_DM,
            "lambda_m": LAMBDA_M,
            "common_dt": float(args.dt),
            "directional_epsilons": [float(x) for x in args.epsilons],
        },
        "resolutions": per_resolution,
        "arithmetic_gate_passed": bool(arithmetic_ok),
        "scientific_admission": False,
        "interpretation_limits": [
            "This is an initialized-state, semi-discrete tangent audit, not a long trajectory.",
            "The source-sector Mdot terms are total directional derivatives of each projected J sector along the full state and CMC-lapse direction; they are not frozen-metric partial derivatives.",
            "The signed decomposition diagnoses the implemented discrete RHS only. It does not identify a faulty physical equation by itself.",
            "The geometric subterm budget is a term-by-term derivative of the pinned vendor's exact momentum-constraint expression, with no alternate stencil. "
            "The one-step and half-step checks are local finite-step comparisons, not production scientific admission.",
            "The accepted-state observer, when enabled, evaluates the signed tangent and native/vendor geometry-only remainder after one ordinary step at the smallest requested resolution; it is not a multi-step trajectory or source-law admission test.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[60, 80, 120])
    parser.add_argument("--r-max", type=float, default=80.0)
    parser.add_argument("--dt", type=float, default=0.015)
    parser.add_argument("--epsilons", nargs=2, type=float, default=[2.0e-5, 1.0e-5])
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/current-source-momentum-tangent/report.json"),
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        report = run(args)
        code = 0 if report["arithmetic_gate_passed"] else 2
    except Exception as exc:
        report = {
            "schema": SCHEMA,
            "status": "failed_closed",
            "scientific_admission": False,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        code = 1
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
