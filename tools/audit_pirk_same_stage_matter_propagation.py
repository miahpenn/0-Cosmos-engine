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


def geometric_constraint_tangent(
    state: ProductionState, rates: dict[str, np.ndarray],
) -> np.ndarray:
    """Exact chain-rule derivative of the implemented discrete geometric M_r."""
    g, grid = state.geometry, state.grid
    d = lambda x: d1(grid, x, 1)
    r = finite("r", grid.centers)
    aa, b, X = g.Aa, g.b, g.X
    aadot, kdot, bdot, xdot = rates["Aa"], rates["K"], rates["b"], rates["X"]
    chi = -0.5 * d(X) / X
    chi_dot = -0.5 * (d(xdot) / X - d(X) * xdot / X**2)
    q = 2.0 / r + d(b) / b
    q_dot = d(bdot) / b - d(b) * bdot / b**2
    return finite(
        "geometric M_r tangent",
        d(aadot) - (2.0 / 3.0) * d(kdot)
        + 6.0 * aadot * chi + 6.0 * aa * chi_dot
        + 1.5 * aadot * q + 1.5 * aa * q_dot,
    )


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


def audit_bundle(bundle_dir: str | Path, resolution: int) -> dict[str, Any]:
    report, arrays, bundle = verify_stage_bundle(Path(bundle_dir), resolution)
    state = restore_stage_state(arrays, report, resolution)
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
    constraints = vacuum.constraints(grid, geom)
    M_geom = finite("geometric momentum constraint", constraints["momentum"])
    M_total = M_geom - KAPPA * finite("total stress current", stress.j)
    tangent_geom = geometric_constraint_tangent(state, rates)
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

    member_hashes = bundle.get("sha256", {})
    return {
        "schema": "issue4-same-stage-matter-brown-propagation/v1",
        "status": "completed_rhs_only_diagnostic",
        "scientific_admission": False,
        "resolution": resolution,
        "stage": STAGE,
        "time": float(state.t),
        "dt_from_copied_capture": float(report["step"]["dt"]),
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "artifact_member_sha256": member_hashes,
        "no_step_or_trajectory": True,
        "no_lapse_resolve": True,
        "source_projection_reconstruction_max_abs": float(np.max(np.abs(j_rebuilt - stress.j))),
        "source_j_stage_crosscheck_max_abs": float(np.max(np.abs(stress.j - arrays["source_start_lapse_refreshed__j_total"]))),
        "geometric_M_stage_crosscheck_max_abs": float(np.max(np.abs(M_geom - arrays["source_start_lapse_refreshed__M_r_geometry"]))),
        "current_tangent_by_sector": {k: summary(v, r) for k, v in current_parts_dot.items()},
        "j_dot_total": summary(j_dot, r),
        "geometric_M_dot": summary(tangent_geom, r),
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--resolution", type=int, choices=(60, 120, 240), required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    result = audit_bundle(args.input_dir, args.resolution)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.out_dir / "report.json"
    report_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "resolution": args.resolution,
        "status": result["status"],
        "j_dot_total_max_abs_r_le_10": result["j_dot_total"]["max_abs_r_le_10"],
        "bianchi_source_max_abs_r_le_10": result["Bianchi_source_minus_kappa_alpha_D"]["max_abs_r_le_10"],
        "residual_max_abs_r_le_10": result["accepted_minus_Brown_source_aware_residual"]["max_abs_r_le_10"],
        "scientific_admission": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
