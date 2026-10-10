"""Opt-in diagnostic: discrete-consistent initial metric function.

The production initializer builds B by a radial cumulative integral. That
satisfies the continuum radial equation to O(dr^2) at the first cell, and the
Hamiltonian divides by r0^2, so the first-cell residual stays O(1) under
refinement. This module instead solves for B so that the discrete Hamiltonian
the evolution uses (vendor constraints minus 16 pi rho) vanishes at every cell.

Opt-in only. The default kernel initializer is unchanged. Not admitted for
production use: it changes initial data and must pass the admission gates.
"""
import math

import numpy as np

from .stress_energy import assemble_total_stress_energy
from . import v55_pirk_adapter as adapter


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
    return np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(tot.rho)


def discrete_consistent_state(kernel, resolution, r_max, amplitude=0.01, width=7.0,
                              D_amplitude=1.0e-10, include_radiation=True,
                              max_iter=12, tol=1.0e-13):
    """Return (state, B, residual_history). Newton solve in B at every cell."""
    state = kernel.initialize(resolution=resolution, r_max=r_max, amplitude=amplitude,
                              width=width, D_amplitude=D_amplitude,
                              include_radiation=include_radiation)
    K_base = np.asarray(state.geometry.K).copy()
    Aa_base = np.asarray(state.geometry.Aa).copy()
    B = np.asarray(state.geometry.X) ** 6
    H = _residual(kernel, state, B, K_base, Aa_base)
    hist = [float(np.max(np.abs(H)))]
    for _ in range(max_iter):
        if hist[-1] < tol:
            break
        J = np.zeros((len(B), len(B)))
        eps = 1.0e-7
        for j in range(len(B)):
            Bp = B.copy()
            Bp[j] *= 1.0 + eps
            J[:, j] = (_residual(kernel, state, Bp, K_base, Aa_base) - H) / (Bp[j] - B[j])
        B = B + np.linalg.solve(J, -H)
        H = _residual(kernel, state, B, K_base, Aa_base)
        hist.append(float(np.max(np.abs(H))))
    _residual(kernel, state, B, K_base, Aa_base)  # leave geometry at the solved B
    return state, B, hist
