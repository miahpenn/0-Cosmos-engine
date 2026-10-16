"""Project-owned spherical-grid operators for the native-kernel migration.

This module is deliberately independent of the optional/pinned vendor package.
The formulas here follow the cell-centred spherical finite-volume geometry and
polynomial moment conditions documented in the method comments. The vendor is
used only by separate tests as a compatibility witness during migration.

Phase 1 contract: this module is not wired into the production adapter.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral

import numpy as np


def _real_profile(values: np.ndarray, shape: tuple[int, ...], name: str) -> np.ndarray:
    """Return a real numeric array with the exact requested shape."""
    raw = np.asarray(values)
    if raw.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    if raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numeric values")
    return np.asarray(raw, dtype=float)


def _parity_sign(parity: int) -> int:
    if isinstance(parity, (bool, np.bool_)) or parity not in (-1, 1):
        raise ValueError("parity must be +1 (even) or -1 (odd)")
    return int(parity)


@dataclass(frozen=True)
class NativeSphericalCellGrid:
    """Uniform cell-centred radial grid with exact spherical shell weights.

    Cell centres are r_i=(i+1/2)h and faces are r_(i+1/2)=i h, where
    h=r_max/n. Shell volumes are the exact integral of r^2 dr over each cell,
    per unit solid angle: (r_outer^3-r_inner^3)/3. A common factor 4*pi may be
    applied by callers when integrating over a full sphere.
    """

    n: int
    r_max: float

    def __post_init__(self) -> None:
        if isinstance(self.n, (bool, np.bool_)) or not isinstance(self.n, Integral):
            raise TypeError("n must be an integer cell count")
        if self.n < 4:
            raise ValueError("n must be at least four")
        try:
            radius = float(self.r_max)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("r_max must be finite and positive") from exc
        if not math.isfinite(radius) or radius <= 0.0:
            raise ValueError("r_max must be finite and positive")
        object.__setattr__(self, "n", int(self.n))
        object.__setattr__(self, "r_max", radius)

    @property
    def dr(self) -> float:
        return self.r_max / self.n

    @property
    def faces(self) -> np.ndarray:
        return np.arange(self.n + 1, dtype=float) * self.dr

    @property
    def centers(self) -> np.ndarray:
        return (np.arange(self.n, dtype=float) + 0.5) * self.dr

    @property
    def face_areas(self) -> np.ndarray:
        return self.faces**2

    @property
    def volumes(self) -> np.ndarray:
        faces = self.faces
        return (faces[1:]**3 - faces[:-1]**3) / 3.0

    def gradient_faces(
        self,
        values: np.ndarray,
        *,
        inner_flux_gradient: float = 0.0,
        outer_flux_gradient: float = 0.0,
    ) -> np.ndarray:
        """Second-order difference from cell values to interior faces.

        Boundary gradients are explicit inputs rather than hidden extrapolation.
        """
        u = _real_profile(values, (self.n,), "values")
        result = np.empty(self.n + 1, dtype=float)
        result[0] = inner_flux_gradient
        result[-1] = outer_flux_gradient
        result[1:-1] = (u[1:] - u[:-1]) / self.dr
        return result

    def gradient_faces_fourth(
        self,
        values: np.ndarray,
        *,
        parity: int,
        inner_flux_gradient: float = 0.0,
        outer_flux_gradient: float = 0.0,
    ) -> np.ndarray:
        """Fourth-order interior cell-to-face gradient with declared parity.

        The first internal face uses a reflected cell value. The last internal
        face retains the documented first-difference closure, so convergence
        claims must mask/exclude that boundary cell where appropriate.
        """
        sign = _parity_sign(parity)
        u = _real_profile(values, (self.n,), "values")
        result = np.empty(self.n + 1, dtype=float)
        result[0] = inner_flux_gradient
        result[-1] = outer_flux_gradient
        result[1] = (sign * u[0] - 27.0 * u[0] + 27.0 * u[1] - u[2]) / (24.0 * self.dr)
        result[2:-2] = (u[:-3] - 27.0 * u[1:-2] + 27.0 * u[2:-1] - u[3:]) / (24.0 * self.dr)
        result[-2] = (u[-1] - u[-2]) / self.dr
        return result

    def interpolate_faces_fourth(self, values: np.ndarray, *, parity: int) -> np.ndarray:
        """Fourth-order interior cell-to-face interpolation.

        Boundary values use the declared parity at the origin and the existing
        constant outer-face closure. That outer closure is not fourth order.
        """
        sign = _parity_sign(parity)
        u = _real_profile(values, (self.n,), "values")
        result = np.empty(self.n + 1, dtype=float)
        result[0] = 0.5 * (u[0] + sign * u[0])
        result[1] = (-sign * u[0] + 9.0 * u[0] + 9.0 * u[1] - u[2]) / 16.0
        result[2:-2] = (-u[:-3] + 9.0 * u[1:-2] + 9.0 * u[2:-1] - u[3:]) / 16.0
        result[-2] = 0.5 * (u[-2] + u[-1])
        result[-1] = u[-1]
        return result

    def divergence_faces(self, flux: np.ndarray) -> np.ndarray:
        """Conservative spherical divergence from face fluxes to cell values."""
        f = _real_profile(flux, (self.n + 1,), "flux")
        area_weighted = self.face_areas * f
        return (area_weighted[1:] - area_weighted[:-1]) / self.volumes

    def diffusion(
        self,
        values: np.ndarray,
        coefficient_faces: np.ndarray | float = 1.0,
        *,
        inner_flux_gradient: float = 0.0,
        outer_flux_gradient: float = 0.0,
    ) -> np.ndarray:
        """Conservative divergence of a face-centred coefficient times gradient."""
        gradient = self.gradient_faces(
            values,
            inner_flux_gradient=inner_flux_gradient,
            outer_flux_gradient=outer_flux_gradient,
        )
        coefficient_raw = np.asarray(coefficient_faces)
        if coefficient_raw.dtype.kind not in "iuf":
            raise ValueError("coefficient_faces must be real numeric values")
        coefficient = np.asarray(coefficient_raw, dtype=float)
        if coefficient.ndim == 0:
            coefficient = np.full(self.n + 1, float(coefficient))
        if coefficient.shape != (self.n + 1,):
            raise ValueError("coefficient_faces must be scalar or shape (n+1,)")
        return self.divergence_faces(coefficient * gradient)

    def integration_by_parts_residual(self, values: np.ndarray, flux: np.ndarray) -> float:
        """Return the discrete spherical summation-by-parts residual.

        With shell weights V_i and face areas A_f, the volume divergence and
        cell-jump terms reduce exactly to the two boundary flux terms.
        """
        u = _real_profile(values, (self.n,), "values")
        f = _real_profile(flux, (self.n + 1,), "flux")
        volume_term = np.sum(self.volumes * u * self.divergence_faces(f))
        jump_term = np.sum(self.face_areas[1:-1] * f[1:-1] * (u[1:] - u[:-1]))
        boundary_term = self.face_areas[-1] * f[-1] * u[-1] - self.face_areas[0] * f[0] * u[0]
        return float(volume_term + jump_term - boundary_term)

    def cell_derivative(self, values: np.ndarray, *, parity: int) -> np.ndarray:
        """Second-order cell derivative, reflecting the first cell through r=0."""
        sign = _parity_sign(parity)
        u = _real_profile(values, (self.n,), "values")
        result = np.empty_like(u)
        result[0] = (u[1] - sign * u[0]) / (2.0 * self.dr)
        result[1:-1] = (u[2:] - u[:-2]) / (2.0 * self.dr)
        result[-1] = (3.0 * u[-1] - 4.0 * u[-2] + u[-3]) / (2.0 * self.dr)
        return result

    def cell_second_derivative(self, values: np.ndarray, *, parity: int) -> np.ndarray:
        """Second-order second derivative with a parity-reflected centre ghost."""
        sign = _parity_sign(parity)
        u = _real_profile(values, (self.n,), "values")
        result = np.empty_like(u)
        result[0] = (u[1] - 2.0 * u[0] + sign * u[0]) / self.dr**2
        result[1:-1] = (u[2:] - 2.0 * u[1:-1] + u[:-2]) / self.dr**2
        result[-1] = (2.0 * u[-1] - 5.0 * u[-2] + 4.0 * u[-3] - u[-4]) / self.dr**2
        return result

    @staticmethod
    def _polynomial_derivative_weights(offsets: np.ndarray, order: int) -> np.ndarray:
        """Solve polynomial moment equations for a one-sided derivative stencil.

        For nodes x_j, weights w_j satisfy sum_j w_j x_j^p = p! delta_(p,order).
        These moment conditions make the stencil exact for polynomials through
        degree len(nodes)-1 without relying on a vendor-specific stencil table.
        """
        x = np.asarray(offsets, dtype=float)
        if x.ndim != 1 or x.size < 2 or order < 0 or order >= x.size:
            raise ValueError("invalid finite-difference stencil")
        moment_matrix = np.vstack([x**power for power in range(x.size)])
        target = np.zeros(x.size, dtype=float)
        target[order] = float(math.factorial(order))
        return np.linalg.solve(moment_matrix, target)

    def _cell_derivative_fourth(self, values: np.ndarray, *, parity: int, order: int) -> np.ndarray:
        """Fourth-order cell derivative from a centred stencil and polynomial moments."""
        sign = _parity_sign(parity)
        u = _real_profile(values, (self.n,), "values")
        if self.n < 5:
            raise ValueError("fourth-order derivative requires at least five cells")
        if order not in (1, 2):
            raise ValueError("order must be one or two")

        result = np.empty_like(u)
        extended = np.empty(self.n + 4, dtype=float)
        extended[2:-2] = u
        extended[1] = sign * u[0]
        extended[0] = sign * u[1]

        # The centred five-point formula follows directly from polynomial
        # moment conditions on offsets -2,-1,0,1,2.
        for i in range(self.n - 2):
            stencil = extended[i:i + 5]
            if order == 1:
                result[i] = (stencil[0] - 8.0 * stencil[1] + 8.0 * stencil[3] - stencil[4]) / (12.0 * self.dr)
            else:
                result[i] = (-stencil[0] + 16.0 * stencil[1] - 30.0 * stencil[2] + 16.0 * stencil[3] - stencil[4]) / (12.0 * self.dr**2)

        # The outer two cells use a five-node one-sided polynomial stencil.
        # Construct weights from coordinates relative to the target cell, so
        # the equations—not a copied coefficient table—define the closure.
        indices = np.arange(self.n - 5, self.n, dtype=int)
        for i in (self.n - 2, self.n - 1):
            weights = self._polynomial_derivative_weights(indices.astype(float) - float(i), order)
            result[i] = float(weights @ u[indices]) / self.dr**order
        return result

    def cell_derivative_fourth(self, values: np.ndarray, *, parity: int) -> np.ndarray:
        """Fourth-order first derivative with parity-aware regular-origin stencil."""
        return self._cell_derivative_fourth(values, parity=parity, order=1)

    def cell_second_derivative_fourth(self, values: np.ndarray, *, parity: int) -> np.ndarray:
        """Fourth-order second derivative with parity-aware regular-origin stencil."""
        return self._cell_derivative_fourth(values, parity=parity, order=2)
