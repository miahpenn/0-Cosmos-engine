#!/usr/bin/env python3
"""RHS-only same-stage matter/Brown propagation audit from a saved PIRK capture.

Reads the checksum-verified stage-0 arrays of the existing copied-step artifact,
evaluates the accepted instantaneous scalar/fluid RHS at that frozen refreshed
state, and compares the derived total M_r tangent with Brown's source-aware
candidate including the pointwise Bianchi source. It never calls kernel.step,
re-solves a lapse, or advances a state.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.matter_system import ConservedSpecies
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.cmc_gauge import solve_cmc_lapse
from engine.v55_matter import total_matter_projection
from scipy.linalg import solve_banded
from engine.scalar_system import ScalarFields
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import V55MatterState
from tools import audit_brown_vacuum_momentum_propagation as brown

KAPPA = 8.0 * math.pi
R_MAX = 80.0
SOURCE_COMMIT = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
STAGE = "stage_start_lapse_refreshed"
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = ("dark_matter", "baryons", "radiation")
SPECIES_FIELDS = ("rest", "energy_t", "momentum_r")


def finite(name: str, values: Any) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError(f"non-finite values in {name}")
    return out


def d1(grid: Any, values: Any, parity: int = 1) -> np.ndarray:
    return finite("D1", grid.cell_derivative_fourth(finite("D1 input", values), parity=parity))


def state_digest(state: ProductionState) -> str:
    """Stable digest of every base state array used by this same-stage audit."""
    arrays: dict[str, np.ndarray] = {}
    for name in GEOMETRY_FIELDS:
        arrays["geometry." + name] = getattr(state.geometry, name)
    for name in SCALAR_FIELDS:
        arrays["scalars." + name] = getattr(state.scalars, name)
    for sector in (*SPECIES, "charged_dust"):
        item = getattr(state.matter, sector, None)
        if item is None:
            continue
        for name in SPECIES_FIELDS:
            arrays["matter." + sector + "." + name] = getattr(item, name)
    digest = hashlib.sha256()
    for name in sorted(arrays):
        values = np.ascontiguousarray(np.asarray(arrays[name]))
        digest.update(name.encode("utf-8"))
        digest.update(str(values.dtype).encode("ascii"))
        digest.update(repr(values.shape).encode("ascii"))
        digest.update(values.tobytes())
    return digest.hexdigest()


def write_signed_artifacts(out_dir: Path, report: dict[str, Any], arrays: dict[str, np.ndarray]) -> dict[str, str]:
    """Write the complete signed radial arrays and a SHA-256 manifest."""
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "report.json"
    arrays_path = out_dir / "signed_arrays.npz"
    manifest_path = out_dir / "SHA256SUMS.txt"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    np.savez_compressed(arrays_path, **{
        key: finite("saved signed array " + key, value)
        for key, value in sorted(arrays.items())
    })
    sums = {
        "report.json": hashlib.sha256(report_path.read_bytes()).hexdigest(),
        "signed_arrays.npz": hashlib.sha256(arrays_path.read_bytes()).hexdigest(),
    }
    manifest_path.write_text(
        "".join(f"{digest}  {name}\n" for name, digest in sorted(sums.items())),
        encoding="utf-8",
    )
    return sums


def solve_cmc_lapse_parity_center_pde_row(
    grid: Any, c: np.ndarray, m: np.ndarray, rhs: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Solve the same second-order CMC interior rows with a parity-ghost centre PDE row.

    For a cell-centred even lapse, alpha[-1]=alpha[0] is the reflected ghost.
    Thus the first-cell centred differences are D1(alpha)_0=(alpha[1]-alpha[0])/(2h)
    and D2(alpha)_0=(alpha[1]-alpha[0])/h^2. This diagnostic-only alternative
    replaces the pinned equality row alpha[0]-alpha[1]=0.
    """
    h = float(grid.dr)
    n = int(grid.n)
    c = finite("counterfactual CMC c coefficient", c)
    m = finite("counterfactual CMC m coefficient", m)
    rhs = finite("counterfactual CMC right hand side", rhs)
    if c.shape != (n,) or m.shape != (n,) or rhs.shape != (n,) or n < 5:
        raise ValueError("counterfactual CMC arrays do not match the cell-centred grid")
    if not math.isfinite(h) or h <= 0.0:
        raise ValueError("counterfactual CMC grid spacing must be positive")
    band = np.zeros((3, n), dtype=np.float64)
    bvec = rhs.copy()
    # Centre PDE row using the even-reflection ghost alpha[-1]=alpha[0].
    band[1, 0] = -1.0 / h**2 - c[0] / (2.0 * h) - m[0]
    band[0, 1] = 1.0 / h**2 + c[0] / (2.0 * h)
    for i in range(1, n - 1):
        ci = c[i] / (2.0 * h)
        band[2, i - 1] = 1.0 / h**2 - ci
        band[1, i] = -2.0 / h**2 - m[i]
        band[0, i + 1] = 1.0 / h**2 + ci
    # Preserve the production outer normalization.
    band[1, n - 1] = 1.0
    bvec[n - 1] = 1.0
    alpha = finite(
        "counterfactual CMC lapse",
        solve_banded((1, 1), band, bvec, check_finite=False),
    )
    if float(np.min(alpha)) <= 0.0:
        raise ValueError("counterfactual CMC lapse is not positive")
    centre_residual = (
        (alpha[1] - alpha[0]) / h**2
        + c[0] * (alpha[1] - alpha[0]) / (2.0 * h)
        - m[0] * alpha[0] - rhs[0]
    )
    return alpha, float(centre_residual)


def primary_l2_decomposition(
    alpha: np.ndarray,
    vendor_geometry: dict[str, Any],
    primary_l2_aa: np.ndarray,
    primary_l2_k: np.ndarray,
) -> dict[str, np.ndarray]:
    """Recover the pinned vendor's exact Aa/K L2 split on this frozen stage."""
    alpha = finite("primary L2 lapse", alpha)
    drdr_alpha = finite("vendor DrDr_alpha", vendor_geometry["DrDr_alpha"])
    lap_alpha = finite("vendor lap_alpha", vendor_geometry["lap_alpha"])
    rr = finite("vendor Rr_r", vendor_geometry["Rr_r"])
    scalar_r = finite("vendor R", vendor_geometry["R"])
    aa_operator = finite("saved primary L2 Aa", primary_l2_aa)
    k_operator = finite("saved primary L2 K", primary_l2_k)
    aa_lapse_hessian = -(drdr_alpha - lap_alpha / 3.0)
    aa_tracefree_ricci = alpha * (rr - scalar_r / 3.0)
    aa_reconstructed = aa_lapse_hessian + aa_tracefree_ricci
    k_reconstructed = -lap_alpha
    return {
        "DrDr_alpha": drdr_alpha,
        "lap_alpha": lap_alpha,
        "Rr_r": rr,
        "R": scalar_r,
        "Aa_lapse_hessian": finite("Aa L2 lapse Hessian", aa_lapse_hessian),
        "Aa_tracefree_ricci": finite("Aa L2 trace-free Ricci", aa_tracefree_ricci),
        "Aa_operator": aa_operator,
        "Aa_reconstructed": finite("Aa L2 reconstructed", aa_reconstructed),
        "Aa_split_residual": finite("Aa L2 split residual", aa_operator - aa_reconstructed),
        "K_operator": k_operator,
        "K_reconstructed_minus_lap_alpha": finite("K L2 reconstructed", k_reconstructed),
        "K_split_residual": finite("K L2 split residual", k_operator - k_reconstructed),
    }


def summary(values: Any, r: np.ndarray) -> dict[str, Any]:
    arr = finite("summary input", values)
    mask = r <= 10.0 + 1.0e-12
    idxs = np.flatnonzero(mask)
    idx = int(idxs[np.argmax(np.abs(arr[mask]))])
    return {
        "signed_cells_0_1_2_6": {str(i): float(arr[i]) for i in (0, 1, 2, 6) if i < arr.size},
        "max_abs_full_grid": float(np.max(np.abs(arr))),
        "max_abs_r_le_10": float(np.max(np.abs(arr[mask]))),
        "signed_at_max_r_le_10": float(arr[idx]),
        "index_at_max_r_le_10": idx,
        "radius_at_max_r_le_10": float(r[idx]),
        "fixed_radius_linear_samples": {
            str(x): float(np.interp(x, r, arr)) for x in (2.0, 3.0, 5.0, 8.0)
            if r[0] <= x <= r[-1]
        },
    }


def restore_stage_state(
    arrays: Any, report: dict[str, Any], resolution: int,
) -> ProductionState:
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(resolution, R_MAX)
    if not np.allclose(arrays["r"], grid.centers, rtol=0.0, atol=1.0e-12):
        raise ValueError("stage artifact radii disagree with the pinned grid")

    geometry = vacuum.VacuumState(**{
        name: finite(f"stage geometry {name}", arrays[f"{STAGE}__geometry__{name}"]).copy()
        for name in GEOMETRY_FIELDS
    })
    scalars = ScalarFields(*(
        finite(f"stage scalar {name}", arrays[f"{STAGE}__scalars__{name}"]).copy()
        for name in SCALAR_FIELDS
    ))
    items = {}
    for species in SPECIES:
        items[species] = ConservedSpecies(*(
            finite(
                f"stage matter {species}.{name}",
                arrays[f"{STAGE}__matter__{species}__{name}"],
            ).copy()
            for name in SPECIES_FIELDS
        ))
    matter = V55MatterState(
        dark_matter=items["dark_matter"],
        baryons=items["baryons"],
        radiation=items["radiation"],
        charged_dust=None,
    )
    return ProductionState(
        grid=grid,
        geometry=geometry,
        scalars=scalars,
        matter=matter,
        t=float(report["input_artifact"]["input_time"]),
        tau=0.0,
        e_folds=0.0,
        electromagnetic=None,
    )


def accepted_geometry_rates(arrays: Any, grid: Any) -> dict[str, np.ndarray]:
    """Recover the exact stage-0 coordinate-time rates from saved operands."""
    def term(block: str, field: str) -> np.ndarray:
        return finite(
            f"saved stage term {block}.{field}",
            arrays[f"{STAGE}__term__{block}__{field}"],
        )
    rates = {name: term("explicit", name) for name in ("a", "b", "X")}
    for name in ("Aa", "K"):
        rates[name] = term("primary_l2", name) + term("primary_l3", name)
    rates["Lambda"] = arrays[f"{STAGE}__term__lambda_l2"] + arrays[f"{STAGE}__term__lambda_l3"]
    return {name: finite(f"accepted rate {name}", value) for name, value in rates.items()}


def total_current_tangent(
    state: ProductionState,
    rates: dict[str, np.ndarray],
    scalar_rates: ScalarFields,
    matter_rates: dict[str, ConservedSpecies],
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """Differentiate each of the six accepted j_s projections on one frozen state."""
    g, s, grid = state.geometry, state.scalars, state.grid
    d = lambda x: d1(grid, x, 1)
    pieces: dict[str, np.ndarray] = {}

    for name, field, momentum, field_rate, momentum_rate, norm in (
        ("S", s.S, s.PS, scalar_rates.S, scalar_rates.PS, 1.0),
        ("D", s.D, s.PD, scalar_rates.D, scalar_rates.PD, 1.0),
        ("COSMOS_phi", s.phi, s.Pi, scalar_rates.phi, scalar_rates.Pi, KAPPA),
    ):
        pieces[name] = -(
            momentum_rate * d(field) + momentum * d(field_rate)
        ) / norm

    r = finite("r", grid.centers)
    sqrt_gamma = np.sqrt(g.a) * g.b * r**2 / g.X**3
    volume_rate = (
        0.5 * rates["a"] / g.a + rates["b"] / g.b - 3.0 * rates["X"] / g.X
    )
    for name in SPECIES:
        U_r = finite(f"{name} conservative covector", getattr(state.matter, name).momentum_r)
        U_r_dot = finite(f"{name} conservative covector rate", matter_rates[name].momentum_r)
        j = U_r / (sqrt_gamma * KAPPA)
        pieces[name] = U_r_dot / (sqrt_gamma * KAPPA) - j * volume_rate

    total = np.zeros_like(r)
    for value in pieces.values():
        total += value
    return pieces, finite("total current tangent", total), finite("spatial volume rate", volume_rate)


def geometric_constraint_tangent_terms(
    state: ProductionState, rates: dict[str, np.ndarray],
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Signed term-by-term chain rule for the implemented discrete geometric M_r."""
    g, grid = state.geometry, state.grid
    d = lambda x: d1(grid, x, 1)
    r = finite("r", grid.centers)
    aa, b, X = g.Aa, g.b, g.X
    aadot, kdot, bdot, xdot = rates["Aa"], rates["K"], rates["b"], rates["X"]
    Db, DX = d(b), d(X)
    terms = {
        "D1_Aa_dot": d(aadot),
        "minus_two_thirds_D1_K_dot": -(2.0 / 3.0) * d(kdot),
        "six_Aa_dot_chi": 6.0 * aadot * (-0.5 * DX / X),
        "minus_three_Aa_D1_X_dot_over_X": -3.0 * aa * d(xdot) / X,
        "plus_three_Aa_D1_X_X_dot_over_X_squared": 3.0 * aa * DX * xdot / X**2,
        "three_Aa_dot_over_r": 3.0 * aadot / r,
        "one_point_five_Aa_dot_D1_b_over_b": 1.5 * aadot * Db / b,
        "one_point_five_Aa_D1_b_dot_over_b": 1.5 * aa * d(bdot) / b,
        "minus_one_point_five_Aa_D1_b_b_dot_over_b_squared": -1.5 * aa * Db * bdot / b**2,
    }
    total = np.zeros_like(r)
    for value in terms.values():
        total += value
    return (
        {name: finite("Mdot geometry term " + name, value) for name, value in terms.items()},
        finite("geometric M_r tangent", total),
    )


def geometric_constraint_tangent(
    state: ProductionState, rates: dict[str, np.ndarray],
) -> np.ndarray:
    """Compatibility wrapper for the full signed chain-rule budget."""
    return geometric_constraint_tangent_terms(state, rates)[1]


def pointwise_total_divergence_r(
    state: ProductionState,
    stress: Any,
    j_dot: np.ndarray,
    volume_rate: np.ndarray,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Observer-matched pointwise nabla_mu T^mu_r on the zero-shift stage."""
    g, grid = state.geometry, state.grid
    r, alpha = finite("r", grid.centers), finite("lapse", g.alpha)
    sqrt_gamma = np.sqrt(g.a) * g.b * r**2 / g.X**3
    gamma_rr = g.a / g.X**2
    gamma_thth = g.b * r**2 / g.X**2

    # At beta=0, sqrt(-g) T^t_r = sqrt(gamma) j_r and
    # sqrt(-g) T^r_r = alpha sqrt(gamma) p_r.
    dt_density_r = sqrt_gamma * (j_dot + stress.j * volume_rate)
    radial_flux_r = d1(grid, alpha * sqrt_gamma * stress.pr, 1)
    Ttt = stress.rho / alpha**2
    Trr = stress.pr / gamma_rr
    Tthth = stress.pt / gamma_thth
    source_r = -0.5 * (
        Ttt * d1(grid, -(alpha**2), 1)
        + Trr * d1(grid, gamma_rr, 1)
        + 2.0 * Tthth * d1(grid, gamma_thth, 1)
    )
    div = (dt_density_r + radial_flux_r) / (alpha * sqrt_gamma) + source_r
    return finite("pointwise total stress divergence radial", div), {
        "time_density": finite("time density divergence component", dt_density_r / (alpha * sqrt_gamma)),
        "radial_flux": finite("radial flux divergence component", radial_flux_r / (alpha * sqrt_gamma)),
        "metric_source": finite("radial metric source component", source_r),
    }


def evaluate_center_row_counterfactual(
    state: ProductionState,
    arrays: Any,
    kernel: V55ProductionKernel,
    baseline_residual: np.ndarray,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Test a parity-centre CMC row on the same saved pre-refresh profile only."""
    grid = state.grid
    h = float(grid.dr)
    r = finite("counterfactual radii", grid.centers)
    stage_g = state.geometry
    pre_g = stage_g.copy()
    pre_g.alpha = finite("saved pre-refresh alpha", arrays["input__geometry__alpha"]).copy()
    pre_state = ProductionState(
        grid=grid, geometry=pre_g, scalars=state.scalars, matter=state.matter,
        t=state.t, tau=state.tau, e_folds=state.e_folds,
        electromagnetic=state.electromagnetic,
    )

    # Confirm lapse refresh changed alpha only. Do not splice state families.
    unchanged_differences: dict[str, float] = {}
    for name in GEOMETRY_FIELDS:
        if name == "alpha":
            continue
        unchanged_differences["geometry." + name] = float(np.max(np.abs(
            np.asarray(arrays["input__geometry__" + name], dtype=float)
            - np.asarray(getattr(stage_g, name), dtype=float)
        )))
    for name in SCALAR_FIELDS:
        unchanged_differences["scalars." + name] = float(np.max(np.abs(
            np.asarray(arrays["input__scalars__" + name], dtype=float)
            - np.asarray(getattr(state.scalars, name), dtype=float)
        )))
    for sector in SPECIES:
        for name in SPECIES_FIELDS:
            unchanged_differences[f"matter.{sector}.{name}"] = float(np.max(np.abs(
                np.asarray(arrays[f"input__matter__{sector}__{name}"], dtype=float)
                - np.asarray(getattr(getattr(state.matter, sector), name), dtype=float)
            )))
    other_profile_max_difference = max(unchanged_differences.values(), default=0.0)
    if other_profile_max_difference > 1.0e-13:
        raise ValueError("saved lapse-refresh stage differs from input profile outside alpha")

    # Reconstruct the pinned solve from the exact pre-refresh input. Its target
    # kdot is then held fixed while only the first row is changed.
    baseline_alpha_rebuilt, kdot = solve_cmc_lapse(
        grid, pre_g.copy(), state.scalars, state.matter
    )
    baseline_alpha_reconstruction = float(np.max(np.abs(baseline_alpha_rebuilt - stage_g.alpha)))

    a = finite("counterfactual lapse metric a", pre_g.a)
    b = finite("counterfactual lapse metric b", pre_g.b)
    X = finite("counterfactual lapse metric X", pre_g.X)
    D1 = lambda values: finite(
        "counterfactual CMC D1", grid.cell_derivative_fourth(values, parity=1)
    )
    c = -0.5 * D1(a) / a + D1(b) / b - D1(X) / X + 2.0 / r
    inv = X * X / a
    total = total_matter_projection(grid, pre_g, state.scalars, state.matter)
    Q = (
        1.5 * np.asarray(pre_g.Aa, dtype=np.float64) ** 2
        + np.asarray(pre_g.K, dtype=np.float64) ** 2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total["rho"], dtype=np.float64)
            + np.asarray(total["pr"], dtype=np.float64)
            + 2.0 * np.asarray(total["pt"], dtype=np.float64)
        )
    )
    m = Q / inv
    rhs = -float(kdot) / inv
    alpha_cf, cf_centre_residual = solve_cmc_lapse_parity_center_pde_row(grid, c, m, rhs)

    # Residual of the same centre PDE row evaluated on the saved baseline lapse.
    alpha_base = finite("saved baseline refreshed alpha", stage_g.alpha)
    base_centre_residual = (
        (alpha_base[1] - alpha_base[0]) / h**2
        + c[0] * (alpha_base[1] - alpha_base[0]) / (2.0 * h)
        - m[0] * alpha_base[0] - rhs[0]
    )
    centre_profile = pre_g.copy()
    centre_profile.alpha = alpha_cf.copy()
    cf_state = ProductionState(
        grid=grid, geometry=centre_profile, scalars=state.scalars, matter=state.matter,
        t=state.t, tau=state.tau, e_folds=state.e_folds,
        electromagnetic=state.electromagnetic,
    )

    # Instantaneous RHS only: no accepted step, no trajectory advancement.
    scalar_rhs_cf, matter_rhs_cf, _ = kernel._rhs(cf_state)
    terms_cf = kernel._geometry_stage_terms(
        grid, cf_state.geometry, cf_state.scalars, cf_state.matter, cf_state.electromagnetic
    )
    rates_cf = {
        name: finite("counterfactual explicit rate " + name, terms_cf["explicit"][name])
        for name in ("a", "b", "X")
    }
    for name in ("Aa", "K"):
        rates_cf[name] = finite(
            "counterfactual primary rate " + name,
            terms_cf["primary_l2"][name] + terms_cf["primary_l3"][name],
        )
    rates_cf["Lambda"] = finite(
        "counterfactual Lambda rate", terms_cf["lambda_l2"] + terms_cf["lambda_l3"]
    )
    _, jdot_cf, volume_rate_cf = total_current_tangent(
        cf_state, rates_cf, scalar_rhs_cf, matter_rhs_cf
    )
    stress_cf = assemble_total_stress_energy(
        grid, cf_state.geometry, cf_state.scalars, cf_state.matter
    )
    _, vacuum, _ = adapter.vendor_modules()
    constraints_cf = vacuum.constraints(grid, cf_state.geometry)
    Mgeom_cf = finite("counterfactual geometric M", constraints_cf["momentum"])
    Mtotal_cf = Mgeom_cf - KAPPA * finite("counterfactual total j", stress_cf.j)
    geom_terms_cf, Mdot_geom_cf = geometric_constraint_tangent_terms(cf_state, rates_cf)
    Mdot_total_cf = Mdot_geom_cf - KAPPA * jdot_cf

    brown_constraints_cf = {
        "momentum": Mtotal_cf,
        "connection": finite("counterfactual connection", constraints_cf["connection"]),
    }
    _, brown_base_cf = brown.brown_lower_momentum_rhs(
        grid, cf_state.geometry, brown_constraints_cf
    )
    Hmatter_cf = -2.0 * KAPPA * finite("counterfactual total rho", stress_cf.rho)
    Hmatter_gradient_cf = (
        -(D1(cf_state.geometry.alpha) / 3.0) * Hmatter_cf
        + (cf_state.geometry.alpha / 6.0) * D1(Hmatter_cf)
    )
    divergence_cf, _ = pointwise_total_divergence_r(
        cf_state, stress_cf, jdot_cf, volume_rate_cf
    )
    bianchi_cf = -KAPPA * cf_state.geometry.alpha * divergence_cf
    predicted_cf = brown_base_cf + Hmatter_gradient_cf + bianchi_cf
    residual_cf = Mdot_total_cf - predicted_cf
    vendor_geometry_cf = vacuum.geometry_terms(grid, cf_state.geometry)
    l2_cf = primary_l2_decomposition(
        cf_state.geometry.alpha, vendor_geometry_cf,
        terms_cf["primary_l2"]["Aa"], terms_cf["primary_l2"]["K"],
    )

    mask = r <= 10.0 + 1.0e-12
    base_max = float(np.max(np.abs(baseline_residual[mask])))
    cf_max = float(np.max(np.abs(residual_cf[mask])))
    report_cf = {
        "status": "completed_diagnostic_only",
        "input_state": "saved pre-refresh profile; geometry/scalar/matter unchanged",
        "baseline_row": "alpha[0] - alpha[1] = 0",
        "counterfactual_row": "D2(alpha)[0] + c[0]*D1(alpha)[0] - m[0]*alpha[0] = rhs[0], with even ghost alpha[-1]=alpha[0]",
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "local_pre_refresh_state_digest": state_digest(pre_state),
        "pinned_baseline_lapse_reconstruction_max_abs": baseline_alpha_reconstruction,
        "frozen_cmc_target_kdot": float(kdot),
        "max_abs_non_lapse_profile_difference": other_profile_max_difference,
        "non_lapse_profile_equal": other_profile_max_difference == 0.0,
        "baseline_alpha_0_minus_1": float(alpha_base[0] - alpha_base[1]),
        "counterfactual_alpha_0_minus_1": float(alpha_cf[0] - alpha_cf[1]),
        "baseline_centre_pde_row_residual": float(base_centre_residual),
        "counterfactual_centre_pde_row_residual": float(cf_centre_residual),
        "alpha_baseline": summary(alpha_base, r),
        "alpha_counterfactual": summary(alpha_cf, r),
        "alpha_counterfactual_minus_baseline": summary(alpha_cf - alpha_base, r),
        "baseline_max_abs_source_aware_residual_r_le_10": base_max,
        "counterfactual_max_abs_source_aware_residual_r_le_10": cf_max,
        "residual_max_ratio_baseline_over_counterfactual": base_max / cf_max if cf_max > 0.0 else None,
        "baseline_residual": summary(baseline_residual, r),
        "counterfactual_residual": summary(residual_cf, r),
        "counterfactual_Aa_lapse_hessian": summary(l2_cf["Aa_lapse_hessian"], r),
        "counterfactual_Aa_tracefree_ricci": summary(l2_cf["Aa_tracefree_ricci"], r),
        "counterfactual_K_L2": summary(l2_cf["K_operator"], r),
        "no_step_or_trajectory": True,
        "scientific_admission": False,
    }
    arrays_cf = {
        "counterfactual__alpha": alpha_cf.copy(),
        "counterfactual__alpha_minus_baseline": (alpha_cf - alpha_base).copy(),
        "counterfactual__DrDr_alpha": l2_cf["DrDr_alpha"].copy(),
        "counterfactual__lap_alpha": l2_cf["lap_alpha"].copy(),
        "counterfactual__Aa_lapse_hessian": l2_cf["Aa_lapse_hessian"].copy(),
        "counterfactual__Aa_tracefree_ricci": l2_cf["Aa_tracefree_ricci"].copy(),
        "counterfactual__Aa_L2": l2_cf["Aa_operator"].copy(),
        "counterfactual__K_L2": l2_cf["K_operator"].copy(),
        "counterfactual__Mdot_geometry": Mdot_geom_cf.copy(),
        "counterfactual__Mdot_total": Mdot_total_cf.copy(),
        "counterfactual__source_aware_prediction": predicted_cf.copy(),
        "counterfactual__residual": residual_cf.copy(),
        "counterfactual__Mdot_Aa_rate": rates_cf["Aa"].copy(),
        "counterfactual__Mdot_K_rate": rates_cf["K"].copy(),
    }
    for name, values in geom_terms_cf.items():
        arrays_cf["counterfactual__Mdot_geometry_term__" + name] = np.asarray(values).copy()
    return report_cf, arrays_cf


def verify_stage_bundle(bundle_dir: Path, resolution: int) -> tuple[dict[str, Any], Any, dict[str, Any]]:
    required = ("report.json", "stage_arrays.npz", "SHA256SUMS.txt")
    for filename in required:
        if not (bundle_dir / filename).is_file():
            raise FileNotFoundError(f"stage artifact missing {filename}: {bundle_dir}")
    expected: dict[str, str] = {}
    for line in (bundle_dir / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, filename = line.split(maxsplit=1)
        expected[filename.strip().lstrip("*")] = digest
    verified = {}
    for filename in ("report.json", "stage_arrays.npz"):
        if filename not in expected:
            raise ValueError(f"checksum manifest omits {filename}")
        actual = hashlib.sha256((bundle_dir / filename).read_bytes()).hexdigest()
        if actual != expected[filename]:
            raise ValueError(f"SHA-256 mismatch for {filename}: {actual} != {expected[filename]}")
        verified[filename] = actual
    report = json.loads((bundle_dir / "report.json").read_text(encoding="utf-8"))
    if report.get("schema") != "issue4-pirk-stage-resolved-momentum-source/v1":
        raise ValueError("unexpected saved-stage artifact schema")
    if report.get("status") != "completed_diagnostic_only":
        raise ValueError("saved-stage artifact is not a completed diagnostic-only capture")
    if report.get("admission_flags", {}).get("scientific_admission") is not False:
        raise ValueError("saved-stage artifact admission state is not fail-closed")
    pins = report.get("pins", {})
    if pins.get("source_pin_commit") != SOURCE_COMMIT or pins.get("source_commit") != SOURCE_COMMIT:
        raise ValueError("saved-stage artifact has the wrong accepted source anchor")
    if pins.get("vendor_pin") != VENDOR_COMMIT or pins.get("vendor_commit") != VENDOR_COMMIT:
        raise ValueError("saved-stage artifact has the wrong vendor pin")
    if pins.get("protected_blobs_match") is not True or int(pins.get("protected_source_blobs_checked", 0)) != 11:
        raise ValueError("saved-stage artifact does not attest all protected source blobs")
    input_artifact = report.get("input_artifact", {})
    if int(input_artifact.get("profile_resolution", -1)) != resolution:
        raise ValueError("saved-stage artifact resolution mismatch")
    if not math.isclose(float(input_artifact.get("input_time", math.nan)), 0.60, rel_tol=0.0, abs_tol=1.0e-10):
        raise ValueError("saved-stage artifact does not match the accepted t=0.60 state")
    with np.load(bundle_dir / "stage_arrays.npz", allow_pickle=False) as npz:
        arrays = {key: np.asarray(npz[key]).copy() for key in npz.files}
    if "r" not in arrays:
        raise ValueError("saved stage arrays have no radial coordinate")
    return report, arrays, {"sha256": verified, "report_sha256": verified["report.json"]}


def audit_bundle(bundle_dir: str | Path, resolution: int) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    report, arrays, bundle = verify_stage_bundle(Path(bundle_dir), resolution)
    state = restore_stage_state(arrays, report, resolution)
    stage_state_digest = state_digest(state)
    grid, geom = state.grid, state.geometry
    r, dr = finite("r", grid.centers), float(grid.dr)
    kernel = V55ProductionKernel()

    # This is a single frozen RHS evaluation at the captured refreshed stage,
    # not a step, copied evolution, or lapse re-solve.
    scalar_rhs, matter_rhs, _ = kernel._rhs(state)
    rates = accepted_geometry_rates(arrays, grid)
    current_parts_dot, j_dot, volume_rate = total_current_tangent(
        state, rates, scalar_rhs, matter_rhs
    )

    stress = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    _, vacuum, _ = adapter.vendor_modules()
    vendor_geometry = vacuum.geometry_terms(grid, geom)
    constraints = vacuum.constraints(grid, geom)
    l2_arrays = primary_l2_decomposition(
        geom.alpha,
        vendor_geometry,
        arrays[f"{STAGE}__term__primary_l2__Aa"],
        arrays[f"{STAGE}__term__primary_l2__K"],
    )
    M_geom = finite("geometric momentum constraint", constraints["momentum"])
    M_total = M_geom - KAPPA * finite("total stress current", stress.j)
    geometric_Mdot_terms, tangent_geom = geometric_constraint_tangent_terms(state, rates)
    geometric_Mdot_terms_sum = np.zeros_like(r)
    for values in geometric_Mdot_terms.values():
        geometric_Mdot_terms_sum += values
    geometric_Mdot_terms_residual = tangent_geom - geometric_Mdot_terms_sum
    Mdot_total = tangent_geom - KAPPA * j_dot

    # Brown's published structural operator uses the metric-only Hamiltonian.
    # Matter translation follows H_total=H_geom-2*kappa*rho and M_total=M_geom-kappa*j.
    # Its Hamiltonian-gradient matter difference is derived algebraically, not fitted.
    brown_constraints = {
        "momentum": M_total,
        "connection": finite("connection constraint", constraints["connection"]),
    }
    brown_terms, brown_Hgeom_Mtotal = brown.brown_lower_momentum_rhs(
        grid, geom, brown_constraints
    )
    Hmatter = -2.0 * KAPPA * finite("total Eulerian energy density", stress.rho)
    alpha_r = d1(grid, geom.alpha, 1)
    Hmatter_gradient = -(alpha_r / 3.0) * Hmatter + (geom.alpha / 6.0) * d1(grid, Hmatter, 1)

    divergence_r, divergence_components = pointwise_total_divergence_r(
        state, stress, j_dot, volume_rate
    )
    bianchi_source = -KAPPA * geom.alpha * divergence_r
    brown_source_aware = brown_Hgeom_Mtotal + Hmatter_gradient + bianchi_source
    residual = Mdot_total - brown_source_aware

    # Independent projection cross-check: conservative fluid momentum / sqrt(gamma)/kappa
    # must reproduce the canonical per-sector j under this zero-shift representation.
    j_rebuilt = (
        -state.scalars.PS * d1(grid, state.scalars.S, 1)
        -state.scalars.PD * d1(grid, state.scalars.D, 1)
        -state.scalars.Pi * d1(grid, state.scalars.phi, 1) / KAPPA
    )
    for name in SPECIES:
        j_rebuilt += getattr(state.matter, name).momentum_r / (np.sqrt(geom.a) * geom.b * r**2 / geom.X**3 * KAPPA)

    stage_state_digest_after = state_digest(state)
    if stage_state_digest_after != stage_state_digest:
        raise RuntimeError("frozen input stage state mutated during RHS-only audit")

    # Preserve all signed cell profiles used to form the same-stage budget.
    # The upstream stage capture remains a separately checksummed input artifact.
    signed_arrays: dict[str, np.ndarray] = {"r": r.copy()}
    for name in GEOMETRY_FIELDS:
        signed_arrays["stage_state__geometry__" + name] = np.asarray(getattr(geom, name)).copy()
    for name in SCALAR_FIELDS:
        signed_arrays["stage_state__scalars__" + name] = np.asarray(getattr(state.scalars, name)).copy()
    for sector in SPECIES:
        item = getattr(state.matter, sector)
        for name in SPECIES_FIELDS:
            signed_arrays[f"stage_state__matter__{sector}__{name}"] = np.asarray(getattr(item, name)).copy()
    for name, value in rates.items():
        signed_arrays["accepted_geometry_rate__" + name] = np.asarray(value).copy()
    for name in SCALAR_FIELDS:
        signed_arrays["accepted_scalar_rhs__" + name] = np.asarray(getattr(scalar_rhs, name)).copy()
    for sector in SPECIES:
        for name in SPECIES_FIELDS:
            signed_arrays[f"accepted_matter_rhs__{sector}__{name}"] = np.asarray(
                getattr(matter_rhs[sector], name)
            ).copy()
    for name, value in current_parts_dot.items():
        signed_arrays["Jdot_sector__" + name] = np.asarray(value).copy()
    signed_arrays.update({
        "spatial_volume_rate": np.asarray(volume_rate).copy(),
        "stress__rho": np.asarray(stress.rho).copy(),
        "stress__j": np.asarray(stress.j).copy(),
        "stress__pr": np.asarray(stress.pr).copy(),
        "stress__pt": np.asarray(stress.pt).copy(),
        "M_geometry": np.asarray(M_geom).copy(),
        "M_total": np.asarray(M_total).copy(),
        "Mdot_geometry": np.asarray(tangent_geom).copy(),
        "Mdot_geometry_term_reconstruction_residual": np.asarray(geometric_Mdot_terms_residual).copy(),
        "Mdot_total_accepted_RHS": np.asarray(Mdot_total).copy(),
        "primary_l2_DrDr_alpha": l2_arrays["DrDr_alpha"].copy(),
        "primary_l2_lap_alpha": l2_arrays["lap_alpha"].copy(),
        "primary_l2_Rr_r": l2_arrays["Rr_r"].copy(),
        "primary_l2_R": l2_arrays["R"].copy(),
        "primary_l2_Aa_lapse_hessian": l2_arrays["Aa_lapse_hessian"].copy(),
        "primary_l2_Aa_tracefree_ricci": l2_arrays["Aa_tracefree_ricci"].copy(),
        "primary_l2_Aa_operator": l2_arrays["Aa_operator"].copy(),
        "primary_l2_Aa_reconstructed": l2_arrays["Aa_reconstructed"].copy(),
        "primary_l2_Aa_split_residual": l2_arrays["Aa_split_residual"].copy(),
        "primary_l2_K_operator": l2_arrays["K_operator"].copy(),
        "primary_l2_K_reconstructed_minus_lap_alpha": l2_arrays["K_reconstructed_minus_lap_alpha"].copy(),
        "primary_l2_K_split_residual": l2_arrays["K_split_residual"].copy(),
        "connection_Cr": np.asarray(constraints["connection"]).copy(),
        "connection_D1_Cr": d1(grid, constraints["connection"], 1).copy(),
        "connection_D2_Cr": finite(
            "connection D2", grid.cell_second_derivative_fourth(constraints["connection"], parity=1)
        ).copy(),
        "H_matter": np.asarray(Hmatter).copy(),
        "H_matter_gradient": np.asarray(Hmatter_gradient).copy(),
        "div_T_r_pointwise": np.asarray(divergence_r).copy(),
        "Bianchi_source_minus_kappa_alpha_D": np.asarray(bianchi_source).copy(),
        "Brown_homogeneous_source_aware_base": np.asarray(brown_Hgeom_Mtotal).copy(),
        "source_aware_Brown_prediction": np.asarray(brown_source_aware).copy(),
        "accepted_minus_Brown_source_aware_residual": np.asarray(residual).copy(),
        "j_dot_total": np.asarray(j_dot).copy(),
        "j_projection_rebuilt": np.asarray(j_rebuilt).copy(),
    })
    for name, value in geometric_Mdot_terms.items():
        signed_arrays["Mdot_geometry_term__" + name] = np.asarray(value).copy()
    for name, value in divergence_components.items():
        signed_arrays["div_T_component__" + name] = np.asarray(value).copy()
    for name, value in brown_terms.items():
        if isinstance(value, np.ndarray):
            signed_arrays["Brown_term__" + name] = np.asarray(value).copy()

    member_hashes = bundle.get("sha256", {})
    result = {
        "schema": "issue4-same-stage-matter-brown-propagation/v1",
        "status": "completed_rhs_only_diagnostic",
        "scientific_admission": False,
        "resolution": resolution,
        "stage": STAGE,
        "time": float(state.t),
        "dt_from_copied_capture": float(report["step"]["dt"]),
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "source_profile_state_digest": report.get("input_artifact", {}).get("input_state_digest_recomputed"),
        "stage_state_digest": stage_state_digest,
        "stage_state_digest_after_rhs": stage_state_digest_after,
        "stage_state_immutable": stage_state_digest_after == stage_state_digest,
        "input_artifact_report_sha256": member_hashes.get("report.json"),
        "input_stage_arrays_sha256": member_hashes.get("stage_arrays.npz"),
        "artifact_member_sha256": member_hashes,
        "signed_arrays_schema": "issue4-same-stage-signed-arrays/v1",
        "signed_arrays_fields": sorted(signed_arrays),
        "no_step_or_trajectory": True,
        "no_lapse_resolve": True,
        "source_projection_reconstruction_max_abs": float(np.max(np.abs(j_rebuilt - stress.j))),
        "source_j_stage_crosscheck_max_abs": float(np.max(np.abs(stress.j - arrays["source_start_lapse_refreshed__j_total"]))),
        "geometric_M_stage_crosscheck_max_abs": float(np.max(np.abs(M_geom - arrays["source_start_lapse_refreshed__M_r_geometry"]))),
        "current_tangent_by_sector": {k: summary(v, r) for k, v in current_parts_dot.items()},
        "j_dot_total": summary(j_dot, r),
        "geometric_M_dot": summary(tangent_geom, r),
        "geometric_M_dot_term_budget": {k: summary(v, r) for k, v in geometric_Mdot_terms.items()},
        "geometric_M_dot_term_reconstruction_residual": summary(geometric_Mdot_terms_residual, r),
        "primary_L2_same_stage_decomposition": {
            name: summary(values, r) for name, values in l2_arrays.items()
        },
        "total_M_dot_accepted_RHS": summary(Mdot_total, r),
        "total_stress_divergence_pointwise": summary(divergence_r, r),
        "pointwise_divergence_components": {k: summary(v, r) for k, v in divergence_components.items()},
        "Bianchi_source_minus_kappa_alpha_D": summary(bianchi_source, r),
        "Brown_homogeneous_source_aware_base": summary(brown_Hgeom_Mtotal, r),
        "Brown_Hamiltonian_matter_gradient": summary(Hmatter_gradient, r),
        "Brown_terms_metric_H": {k: summary(v, r) for k, v in brown_terms.items()},
        "source_aware_Brown_prediction": summary(brown_source_aware, r),
        "accepted_minus_Brown_source_aware_residual": summary(residual, r),
        "interpretation": (
            "Same-state semi-discrete comparison only. The accepted side uses exact captured "
            "stage-0 geometry rates plus the accepted scalar/fluid RHS. Brown's side uses the "
            "metric-only Hamiltonian translated by -2*kappa*rho, total M_r, and the declared "
            "pointwise Bianchi source. Remaining residual is not a fit target and does not itself "
            "select the off-constraint connection convention."
        ),
    }
    return result, signed_arrays


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--resolution", type=int, choices=(60, 120, 240), required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    result, signed_arrays = audit_bundle(args.input_dir, args.resolution)
    sums = write_signed_artifacts(args.out_dir, result, signed_arrays)
    print(json.dumps({
        "resolution": args.resolution,
        "status": result["status"],
        "j_dot_total_max_abs_r_le_10": result["j_dot_total"]["max_abs_r_le_10"],
        "bianchi_source_max_abs_r_le_10": result["Bianchi_source_minus_kappa_alpha_D"]["max_abs_r_le_10"],
        "residual_max_abs_r_le_10": result["accepted_minus_Brown_source_aware_residual"]["max_abs_r_le_10"],
        "stage_state_digest": result["stage_state_digest"],
        "signed_arrays_sha256": sums["signed_arrays.npz"],
        "scientific_admission": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
