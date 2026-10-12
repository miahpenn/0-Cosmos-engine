"""Strong-D snapshot CMC map/Jacobian audit (diagnostic only).

Loads the admitted boundary-matched R80 full-state snapshot at t=12, preserves
its recorded outer lapse normalization, and measures:
  F(alpha)=solve_cmc_lapse(alpha-dependent target)-alpha,
the incoming-state discrete CMC residual, and a finite-difference Jacobian on
four independent radial lapse modes. No evolution or production edits occur.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import engine.cmc_gauge as cmc
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q, total_matter_projection
from engine import v55_pirk_adapter as adapter

EPSILONS = (1.0e-4, 1.0e-5, 1.0e-6)
TARGET_TIME = 12.0
EXPECTED_N = 160
EXPECTED_RMAX = 80.0
EXPECTED_D_AMPLITUDE = 1.0e-4
OUT = ROOT / "runs" / "cmc_strongD_snapshot_jacobian"


def rms(x):
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def max_abs(x):
    return float(np.max(np.abs(np.asarray(x, dtype=float))))


def profile_summary(values):
    values = np.asarray(values, dtype=float)
    idx = int(np.argmax(np.abs(values)))
    return {
        "rms": rms(values),
        "max_abs": max_abs(values),
        "center": float(values[0]),
        "outer": float(values[-1]),
        "max_abs_index": idx,
        "max_abs_radius": float((idx + 0.5) * EXPECTED_RMAX / EXPECTED_N),
        "max_abs_signed_value": float(values[idx]),
    }


def cosine_similarity(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    denom = float(np.linalg.norm(left) * np.linalg.norm(right))
    return None if denom == 0.0 else float(np.dot(left, right) / denom)


def read_snapshot(path: Path):
    selected = None
    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            if abs(float(row.get("t", -1.0)) - TARGET_TIME) <= 1.0e-8:
                if row.get("state_snapshot") is not None:
                    selected = row
                    break
    if selected is None:
        raise ValueError(f"no full state_snapshot at t={TARGET_TIME} in {path}")
    domain = selected.get("domain", {})
    if int(domain.get("N", -1)) != EXPECTED_N:
        raise ValueError(f"snapshot N mismatch: {domain}")
    if abs(float(domain.get("r_max", -1.0)) - EXPECTED_RMAX) > 1.0e-12:
        raise ValueError(f"snapshot r_max mismatch: {domain}")
    if not selected.get("diagnostic_only", False):
        raise ValueError("expected the archived diagnostic-only snapshot")
    return selected


def load_state(row):
    snap = row["state_snapshot"]
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=EXPECTED_N,
        r_max=EXPECTED_RMAX,
        amplitude=1.0e-2,
        width=7.0,
        D_amplitude=EXPECTED_D_AMPLITUDE,
        include_radiation=True,
    )
    for name, values in snap["geometry"].items():
        arr = np.asarray(values, dtype=float)
        if arr.shape != (state.grid.n,):
            raise ValueError(f"geometry.{name} shape mismatch: {arr.shape}")
        getattr(state.geometry, name)[:] = arr
    for name, values in snap["scalars"].items():
        arr = np.asarray(values, dtype=float)
        if arr.shape != (state.grid.n,):
            raise ValueError(f"scalars.{name} shape mismatch: {arr.shape}")
        getattr(state.scalars, name)[:] = arr
    for species, fields in snap["matter"].items():
        target = getattr(state.matter, species)
        for name, values in fields.items():
            arr = np.asarray(values, dtype=float)
            if arr.shape != (state.grid.n,):
                raise ValueError(f"matter.{species}.{name} shape mismatch: {arr.shape}")
            getattr(target, name)[:] = arr

    radius = np.asarray(row["radius"], dtype=float)
    if radius.shape != (state.grid.n,) or not np.allclose(
        radius, np.asarray(state.grid.centers, dtype=float), rtol=0.0, atol=1.0e-12
    ):
        raise ValueError("snapshot radius does not match native grid centers")
    state.t = float(row["t"])
    state.tau = float(row["tau"])
    return state


def solve_map(state, frozen_geometry, alpha, outer_alpha):
    geom = frozen_geometry.copy()
    geom.alpha = np.asarray(alpha, dtype=float).copy()
    out, target = cmc.solve_cmc_lapse(
        state.grid, geom, state.scalars, state.matter,
        outer_alpha=float(outer_alpha),
    )
    return np.asarray(out, dtype=float), float(target)


def build_modes(grid):
    r = np.asarray(grid.centers, dtype=float)
    x = r / float(r[-1])
    raw = [
        ("central_gaussian", np.exp(-(r / 3.0) ** 2)),
        ("broad_gaussian", np.exp(-(r / 15.0) ** 2)),
        ("off_center_shell", np.exp(-0.5 * ((r - 25.0) / 7.0) ** 2)),
        ("signed_radial_mode", np.cos(2.0 * np.pi * x) * (1.0 - x) ** 2),
    ]
    modes = []
    for name, values in raw:
        v = np.asarray(values, dtype=float)
        v[-1] = 0.0  # keep the snapshot's outer lapse normalization fixed
        scale = rms(v)
        if not np.isfinite(scale) or scale <= 0.0:
            raise ValueError(f"degenerate lapse mode {name}")
        modes.append((name, v / scale))
    return modes


def c_coefficients(grid, geometry):
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    r = np.asarray(grid.centers, dtype=float)
    ap = grid.cell_derivative_fourth(a, parity=1)
    bp = grid.cell_derivative_fourth(b, parity=1)
    Xp = grid.cell_derivative_fourth(X, parity=1)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r
    return c


def residual_for_target(state, coefficient_geometry, alpha_eval, kdot, outer_alpha):
    """Evaluate the exact native fourth-order discrete CMC operator residual."""
    grid = state.grid
    total = total_matter_projection(
        grid, coefficient_geometry, state.scalars, state.matter
    )
    a = np.asarray(coefficient_geometry.a, dtype=float)
    X = np.asarray(coefficient_geometry.X, dtype=float)
    inv = X * X / a
    Q = (
        1.5 * coefficient_geometry.Aa * coefficient_geometry.Aa
        + coefficient_geometry.K * coefficient_geometry.K / 3.0
        + 4.0 * np.pi * (total["rho"] + total["pr"] + 2.0 * total["pt"])
    )
    m = Q / inv
    c = c_coefficients(grid, coefficient_geometry)
    band = cmc._cmc_operator_banded(grid, c, m)
    n = grid.n
    u = 4
    dense = np.zeros((n, n), dtype=float)
    for j in range(n):
        for i in range(max(0, j - u), min(n, j + u + 1)):
            dense[i, j] = band[u + i - j, j]
    rhs = np.zeros(n, dtype=float)
    rhs[:-1] = -float(kdot) / inv[:-1]
    rhs[-1] = float(outer_alpha)
    residual = dense @ np.asarray(alpha_eval, dtype=float) - rhs
    scale = max(
        max_abs(dense @ np.asarray(alpha_eval, dtype=float)),
        max_abs(rhs),
        max_abs(np.asarray(alpha_eval, dtype=float)),
        1.0,
    )
    return residual, float(max_abs(residual) / scale), scale


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: python diagnostics/cmc_strongD_snapshot_jacobian.py <profiles.jsonl>")
    profile_path = Path(sys.argv[1])
    row = read_snapshot(profile_path)
    state = load_state(row)
    grid = state.grid
    geometry0 = state.geometry.copy()
    alpha0 = np.asarray(geometry0.alpha, dtype=float).copy()
    outer_alpha = float(alpha0[-1])

    # Reference provenance from the archived matched-boundary trajectory.
    if not all(np.all(np.isfinite(np.asarray(a, dtype=float))) for a in
               [alpha0, geometry0.a, geometry0.b, geometry0.X, geometry0.Aa, geometry0.K]):
        raise ValueError("snapshot has non-finite geometry values")
    if min(float(np.min(geometry0.a)), float(np.min(geometry0.b)),
           float(np.min(geometry0.X)), float(np.min(alpha0))) <= 0.0:
        raise ValueError("snapshot has non-positive metric or lapse")

    alpha_base, kdot_used = solve_map(state, geometry0, alpha0, outer_alpha)
    base_shift = alpha_base - alpha0
    input_target = float(cmc.target_kdot(grid, geometry0, state.scalars, state.matter))
    geometry_base = geometry0.copy()
    geometry_base.alpha = alpha_base.copy()
    output_target = float(cmc.target_kdot(grid, geometry_base, state.scalars, state.matter))

    r_input, rel_input, scale_input = residual_for_target(
        state, geometry0, alpha0, input_target, outer_alpha
    )
    r_same, rel_same, scale_same = residual_for_target(
        state, geometry0, alpha_base, input_target, outer_alpha
    )
    r_recomputed, rel_recomputed, scale_recomputed = residual_for_target(
        state, geometry_base, alpha_base, output_target, outer_alpha
    )

    modes = build_modes(grid)
    directional = []
    finite = True
    for name, mode in modes:
        eps_results = []
        for eps in EPSILONS:
            plus, target_plus = solve_map(state, geometry0, alpha0 + eps * mode, outer_alpha)
            minus, target_minus = solve_map(state, geometry0, alpha0 - eps * mode, outer_alpha)
            jv = (plus - minus) / (2.0 * eps)
            item = {
                "epsilon": eps,
                "J_G_v": profile_summary(jv),
                "response_rms_gain": rms(jv) / max(rms(mode), 1.0e-300),
                "identity_response_error_rms_gain": rms(jv - mode) / max(rms(mode), 1.0e-300),
                "rayleigh_quotient": float(np.dot(mode, jv) / np.dot(mode, mode)),
                "response_cosine_with_mode": cosine_similarity(mode, jv),
                "dkdot_d_direction": float((target_plus - target_minus) / (2.0 * eps)),
                "finite": bool(np.all(np.isfinite(plus)) and np.all(np.isfinite(minus)) and np.all(np.isfinite(jv))),
            }
            eps_results.append(item)
            finite = finite and item["finite"]
        directional.append({"mode": name, "input_mode": profile_summary(mode), "by_epsilon": eps_results})

    raw = np.column_stack([v for _, v in modes])
    q, _ = np.linalg.qr(raw / np.sqrt(grid.n))
    basis = q * np.sqrt(grid.n)
    epsilon = EPSILONS[1]
    Jcols = []
    off_subspace = []
    for j in range(basis.shape[1]):
        v = basis[:, j]
        gp, _ = solve_map(state, geometry0, alpha0 + epsilon * v, outer_alpha)
        gm, _ = solve_map(state, geometry0, alpha0 - epsilon * v, outer_alpha)
        jv = (gp - gm) / (2.0 * epsilon)
        coeff = (basis.T @ jv) / grid.n
        Jcols.append(coeff)
        off_subspace.append(rms(jv - basis @ coeff))
        finite = finite and bool(np.all(np.isfinite(jv)) and np.all(np.isfinite(coeff)))

    J = np.column_stack(Jcols)
    JF = J - np.eye(J.shape[0])
    eig_G = np.linalg.eigvals(J)
    eig_F = np.linalg.eigvals(JF)
    svals_F = np.linalg.svd(JF, compute_uv=False)
    target_digest = __import__("hashlib").sha256(profile_path.read_bytes()).hexdigest()

    payload = {
        "test": "strong-D snapshot CMC map Jacobian and discrete residual",
        "provenance": {
            "source_run": "https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38162318255",
            "source_artifact_id": 11710858125,
            "source_commit": "35072830f9e4d7029bb10243a0520b9a273a17dd",
            "profile_file": profile_path.name,
            "profile_file_sha256": target_digest,
            "snapshot_label": row.get("label"),
            "schema": row.get("schema"),
            "diagnostic_only_source_snapshot": bool(row.get("diagnostic_only")),
        },
        "state": {
            "t": float(state.t),
            "tau": float(state.tau),
            "N": int(grid.n),
            "r_max": float(grid.centers[-1] + 0.5 * grid.dr),
            "dr": float(grid.dr),
            "D_center": float(state.scalars.D[0]),
            "PD_center": float(state.scalars.PD[0]),
            "min_alpha": float(np.min(alpha0)),
            "max_alpha": float(np.max(alpha0)),
            "outer_alpha_fixed": outer_alpha,
            "min_a": float(np.min(geometry0.a)),
            "min_b": float(np.min(geometry0.b)),
            "min_X": float(np.min(geometry0.X)),
        },
        "production_equations_changed": False,
        "fixed_point_map": "G(alpha)=solve_cmc_lapse(alpha-dependent proper-volume target), F(alpha)=G(alpha)-alpha; outer alpha held at snapshot alpha(R)",
        "base_state": {
            "input_target_kdot": input_target,
            "target_used_by_map": float(kdot_used),
            "output_recomputed_target_kdot": output_target,
            "target_output_minus_input": float(output_target - input_target),
            "map_offset_F": profile_summary(base_shift),
            "input_slice_cmc_residual": profile_summary(r_input),
            "input_slice_cmc_residual_relative_max": rel_input,
            "mapped_lapse_same_target_residual": profile_summary(r_same),
            "mapped_lapse_same_target_residual_relative_max": rel_same,
            "mapped_lapse_recomputed_target_residual": profile_summary(r_recomputed),
            "mapped_lapse_recomputed_target_residual_relative_max": rel_recomputed,
            "input_residual_max_scale": scale_input,
            "mapped_same_target_residual_max_scale": scale_same,
            "mapped_recomputed_residual_max_scale": scale_recomputed,
            "map_shift_residual_cosine_descriptive_only": cosine_similarity(base_shift, r_input),
        },
        "finite_difference_epsilons": list(EPSILONS),
        "directional_modes": directional,
        "reduced_jacobian": {
            "basis": "RMS-orthonormal span of four smooth center-regular, outer-zero radial profiles",
            "epsilon": epsilon,
            "J_G_matrix": J.tolist(),
            "J_G_eigenvalues": [{"real": float(np.real(z)), "imag": float(np.imag(z))} for z in eig_G],
            "J_F_matrix": JF.tolist(),
            "J_F_eigenvalues": [{"real": float(np.real(z)), "imag": float(np.imag(z))} for z in eig_F],
            "J_F_singular_values": [float(v) for v in svals_F],
            "min_J_F_singular_value": float(np.min(svals_F)),
            "off_subspace_response_rms_by_basis_direction": [float(v) for v in off_subspace],
            "max_off_subspace_response_rms": float(max(off_subspace)),
            "interpretation_limit": "local finite differences on one snapshot and four radial directions only",
        },
        "proper_time_context": {
            "archived_R80_minus_R160_tau_at_t12": 0.0106210937742,
            "archived_baseline_R80_minus_R160_tau_at_t12": 0.0117426496868,
            "warning": "A same-slice CMC defect is a gauge consistency diagnostic; it does not alone explain an accumulated trajectory proper-time difference.",
        },
        "checks": {
            "all_finite": bool(finite and np.all(np.isfinite(r_input)) and np.all(np.isfinite(r_same)) and np.all(np.isfinite(r_recomputed))),
            "mapped_lapse_closes_same_target_equation": bool(rel_same < 1.0e-10),
            "positive_snapshot_and_mapped_lapse": bool(np.min(alpha0) > 0.0 and np.min(alpha_base) > 0.0),
        },
    }
    payload["pass"] = all(payload["checks"].values())
    OUT.mkdir(parents=True, exist_ok=True)
    out_path = OUT / "cmc_strongD_snapshot_jacobian.json"
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    print(f"RESULT_JSON={out_path}", flush=True)
    if not payload["pass"]:
        raise SystemExit("FAIL: non-finite, non-positive, or same-target equation closure failed")
    print("PASS: strong-D snapshot Jacobian/residual diagnostic complete; no production physics changed")


if __name__ == "__main__":
    main()
