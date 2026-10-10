"""Opt-in diagnostic: discrete-consistent initial metric function.

The production initializer builds B by a radial cumulative integral. This
module solves for B so that the discrete Hamiltonian used by evolution
(vendor constraint minus 16 pi rho) vanishes cell by cell.

Opt-in only. The default kernel initializer is unchanged. This changes initial
data and is not admitted for production use; it must pass the admission gates.
"""
import math

import numpy as np

from .stress_energy import assemble_total_stress_energy
from . import v55_pirk_adapter as adapter


class ConvergenceError(RuntimeError):
    """The discrete-consistent initial-data solve did not safely converge."""

    def __init__(self, message, diagnostics=None):
        super().__init__(message)
        self.diagnostics = list(diagnostics or [])


def _residual(kernel, state, B, K_base, Aa_base):
    grid, G = state.grid, state.geometry
    r = np.asarray(grid.centers)
    _, vacuum, _ = adapter.vendor_modules()
    A = 1.0 / np.sqrt(B)
    G.a[:] = A ** (4.0 / 3.0)
    G.b[:] = A ** (-2.0 / 3.0)
    G.X[:] = A ** (-1.0 / 3.0)
    G.Lambda[:] = (
        np.asarray(grid.cell_derivative_fourth(G.a, parity=1)) / (2.0 * G.a * G.a)
        - np.asarray(grid.cell_derivative_fourth(G.b, parity=1)) / (G.a * G.b)
        + 2.0 / r * (1.0 / G.b - 1.0 / G.a)
    )
    G.K[:] = K_base
    G.Aa[:] = Aa_base
    raw = vacuum.constraints(grid, G)
    tot = assemble_total_stress_energy(grid, G, state.scalars, state.matter)
    return (
        np.asarray(raw["hamiltonian"])
        - 16.0 * math.pi * np.asarray(tot.rho)
    )


def discrete_consistent_state(
    kernel,
    resolution,
    r_max,
    amplitude=0.01,
    width=7.0,
    D_amplitude=1.0e-10,
    include_radiation=True,
    max_iter=12,
    tol=1.0e-13,
    max_cond=1.0e12,
    return_info=False,
):
    """Return (state, B, residual_history), optionally with solver diagnostics.

    Newton updates are accepted only when the candidate metric B remains
    finite and strictly positive. The solve fails closed on invalid inputs,
    non-finite residuals, singular/ill-conditioned Jacobians, or non-convergence.
    """
    if isinstance(max_iter, (bool, np.bool_)) or not isinstance(
        max_iter, (int, np.integer)
    ) or max_iter < 0:
        raise ValueError("max_iter must be a non-negative integer")
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive")
    if not np.isfinite(max_cond) or max_cond <= 1.0:
        raise ValueError("max_cond must be finite and greater than one")

    state = kernel.initialize(
        resolution=resolution,
        r_max=r_max,
        amplitude=amplitude,
        width=width,
        D_amplitude=D_amplitude,
        include_radiation=include_radiation,
    )
    K_base = np.asarray(state.geometry.K).copy()
    Aa_base = np.asarray(state.geometry.Aa).copy()
    B = np.asarray(state.geometry.X, dtype=float) ** 6

    if B.ndim != 1 or B.size == 0 or not np.all(np.isfinite(B)) or np.any(B <= 0.0):
        raise ConvergenceError("initial metric B must be finite and strictly positive")

    def evaluate(candidate):
        value = np.asarray(
            _residual(kernel, state, candidate, K_base, Aa_base), dtype=float
        )
        if value.shape != B.shape:
            raise ConvergenceError("Hamiltonian residual has an unexpected shape")
        if not np.all(np.isfinite(value)):
            raise ConvergenceError("Hamiltonian residual contains non-finite values")
        return value

    H = evaluate(B)
    hist = [float(np.max(np.abs(H)))]
    max_condition_seen = None
    iterations = 0
    diagnostics = []

    for iteration in range(max_iter):
        if hist[-1] <= tol:
            break

        record = {
            "iteration": iteration + 1,
            "residual_before": float(hist[-1]),
            "jacobian_relative_eps": 1.0e-7,
            "jacobian_condition": None,
            "update_linf": None,
            "update_l2": None,
            "candidate_B_min": None,
            "residual_after": None,
            "status": "building_jacobian",
        }
        diagnostics.append(record)
        J = np.zeros((B.size, B.size), dtype=float)
        eps = record["jacobian_relative_eps"]
        for j in range(B.size):
            Bp = B.copy()
            Bp[j] *= 1.0 + eps
            denominator = Bp[j] - B[j]
            if not np.isfinite(denominator) or denominator == 0.0:
                raise ConvergenceError("finite-difference Jacobian step is invalid")
            Hp = evaluate(Bp)
            J[:, j] = (Hp - H) / denominator

        if not np.all(np.isfinite(J)):
            raise ConvergenceError("Jacobian contains non-finite values")
        try:
            condition = float(np.linalg.cond(J))
        except np.linalg.LinAlgError as exc:
            raise ConvergenceError("could not estimate Jacobian condition") from exc
        if max_condition_seen is None or (
            np.isfinite(condition) and condition > max_condition_seen
        ):
            max_condition_seen = condition
        record["jacobian_condition"] = condition
        record["status"] = "jacobian_conditioned"
        if not np.isfinite(condition) or condition > max_cond:
            record["status"] = "jacobian_condition_rejected"
            raise ConvergenceError(
                f"Jacobian condition {condition!r} exceeds limit {max_cond:g}",
                diagnostics=diagnostics,
            )

        try:
            delta = np.linalg.solve(J, -H)
        except np.linalg.LinAlgError as exc:
            raise ConvergenceError("Jacobian solve is singular") from exc
        record["update_linf"] = float(np.max(np.abs(delta)))
        record["update_l2"] = float(np.linalg.norm(delta))
        if not np.all(np.isfinite(delta)):
            record["status"] = "nonfinite_update"
            raise ConvergenceError(
                "Newton update contains non-finite values", diagnostics=diagnostics
            )

        candidate = B + delta
        record["candidate_B_min"] = float(np.min(candidate))
        if not np.all(np.isfinite(candidate)) or np.any(candidate <= 0.0):
            record["status"] = "candidate_B_invalid"
            raise ConvergenceError(
                "Newton update would make metric B non-finite or non-positive",
                diagnostics=diagnostics,
            )

        B = candidate
        H = evaluate(B)
        hist.append(float(np.max(np.abs(H))))
        record["residual_after"] = hist[-1]
        record["status"] = "updated"
        iterations += 1

    if hist[-1] > tol:
        raise ConvergenceError(
            f"did not converge in {max_iter} iterations; "
            f"final max |H|={hist[-1]:.17g}, tolerance={tol:.17g}",
            diagnostics=diagnostics,
        )

    # Leave the returned state's geometry at the converged B, not at the last
    # finite-difference probe evaluated while assembling the Jacobian.
    evaluate(B)
    info = {
        "iterations": iterations,
        "max_condition_number": max_condition_seen,
        "final_max_residual": hist[-1],
        "tolerance": float(tol),
    }
    if return_info:
        return state, B, hist, info
    return state, B, hist
