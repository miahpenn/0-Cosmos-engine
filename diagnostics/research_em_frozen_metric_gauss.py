"""Research-only finite-volume Gauss solve on frozen production geometry.

The source charge is a manufactured test profile, not a charge assignment to
any GEAR field/species and not a claim that the engine has electromagnetism.
The input is a read-only, initialized V5.5 ProductionState; this module does
not step it, change its arrays, or enter the production RHS.

For the active spherical BSSN variables:
    gamma_rr = a / X**2
    gamma_thth = b * r**2 / X**2
    sqrt_gamma (without sin(theta)) = sqrt(a) * b * r**2 / X**3
    R_areal = r * sqrt(b) / X

The discrete Gauss law is d_r(sqrt_gamma E^r) = 4*pi*rho_e*sqrt_gamma.
The finite-volume field below is reconstructed from shell source charges,
so exact shell closure is an algebraic constraint-solve check, not evidence
of Maxwell time evolution or matter/Einstein source coupling.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


FOUR_PI = 4.0 * math.pi


def _finite_real(value: object, name: str) -> float:
    numeric_types = (int, float, np.integer, np.floating)
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numeric_types):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


def _positive_array(values: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    try:
        result = np.asarray(values, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite positive array") from exc
    if result.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {result.shape}")
    if not np.all(np.isfinite(result)) or np.any(result <= 0.0):
        raise ValueError(f"{name} must be finite and strictly positive")
    return result


def _interpolate_even_cell_centred(
    values: np.ndarray, centers: np.ndarray, points: np.ndarray
) -> np.ndarray:
    """Linear interior interpolation, even reflection at r=0, linear outer extrapolation.

    The production arrays are cell-centred and regular spherical metric fields
    are even at the origin. Setting the reflected ghost value equal to the
    first centre value imposes zero first derivative at the origin; because
    the data are sampled at dr/2, this gives a second-order face value there.
    """
    flat = np.asarray(points, dtype=float).reshape(-1)
    result = np.interp(flat, centers, values)
    left = flat < centers[0]
    result[left] = values[0]
    right = flat > centers[-1]
    if np.any(right):
        slope = (values[-1] - values[-2]) / (centers[-1] - centers[-2])
        result[right] = values[-1] + slope * (flat[right] - centers[-1])
    return result.reshape(np.shape(points))


def _metric_at(
    radii: np.ndarray,
    centers: np.ndarray,
    a_cells: np.ndarray,
    b_cells: np.ndarray,
    X_cells: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    a = _interpolate_even_cell_centred(a_cells, centers, radii)
    b = _interpolate_even_cell_centred(b_cells, centers, radii)
    X = _interpolate_even_cell_centred(X_cells, centers, radii)
    for name, values in (("a", a), ("b", b), ("X", X)):
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError(
                f"interpolated metric {name} must be finite and strictly positive"
            )

    sqrt_gamma = np.sqrt(a) * b * radii**2 / X**3
    gamma_rr = a / X**2
    areal_radius = radii * np.sqrt(b) / X
    if (
        not np.all(np.isfinite(sqrt_gamma))
        or not np.all(np.isfinite(gamma_rr))
        or not np.all(np.isfinite(areal_radius))
        or np.any(sqrt_gamma < 0.0)
        or np.any(gamma_rr <= 0.0)
        or np.any(areal_radius < 0.0)
    ):
        raise ValueError("interpolated metric produced invalid volume/area factors")
    return sqrt_gamma, gamma_rr, areal_radius, a, b


@dataclass(frozen=True)
class FrozenProductionGaussResult:
    """Read-only-by-convention output for a frozen-slice Gauss constraint solve."""

    faces: np.ndarray
    cell_charges: np.ndarray
    face_charges: np.ndarray
    face_sqrt_gamma: np.ndarray
    face_gamma_rr: np.ndarray
    face_areal_radius: np.ndarray
    face_electric_field_r: np.ndarray
    face_electric_field_hat: np.ndarray
    cell_residuals: np.ndarray
    total_charge: float
    midpoint_total_charge: float
    max_abs_cell_residual: float


def solve_frozen_production_metric_gauss(
    state: object,
    *,
    rho0: float = 0.75,
    source_radius: float = 2.0,
    quadrature_order: int = 8,
) -> FrozenProductionGaussResult:
    """Solve the static discrete Gauss constraint on a frozen production slice.

    A manufactured compact profile rho_e = rho0*(1-(r/R0)**2)**2 for r<R0
    is integrated against the interpolated proper-volume factor of the actual
    initialized production geometry. R0 must coincide with a cell face so no
    cell straddles the compact-support edge. An independent centre-sampled
    shell-volume estimate is returned to expose refinement behavior.
    """
    rho0_value = _finite_real(rho0, "rho0")
    radius_value = _finite_real(source_radius, "source_radius")
    if radius_value <= 0.0:
        raise ValueError("source_radius must be positive")
    if (
        isinstance(quadrature_order, bool)
        or not isinstance(quadrature_order, int)
        or not 2 <= quadrature_order <= 32
    ):
        raise ValueError("quadrature_order must be an integer from 2 through 32")

    try:
        grid = state.grid
        geometry = state.geometry
        faces = np.asarray(grid.faces, dtype=float)
        centers = np.asarray(grid.centers, dtype=float)
        n = int(grid.n)
        a_cells = _positive_array(geometry.a, (n,), "a")
        b_cells = _positive_array(geometry.b, (n,), "b")
        X_cells = _positive_array(geometry.X, (n,), "X")
        grid_volumes = np.asarray(grid.volumes, dtype=float)
    except AttributeError as exc:
        raise ValueError("state must expose production grid and geometry arrays") from exc

    if faces.shape != (n + 1,) or centers.shape != (n,):
        raise ValueError("grid faces/centres do not match grid.n")
    if (
        not np.all(np.isfinite(faces))
        or not np.all(np.isfinite(centers))
        or faces[0] != 0.0
        or np.any(np.diff(faces) <= 0.0)
        or np.any(centers <= faces[:-1])
        or np.any(centers >= faces[1:])
    ):
        raise ValueError("grid faces/centres must be finite and strictly ordered")
    expected_centers = 0.5 * (faces[:-1] + faces[1:])
    if not np.allclose(
        centers,
        expected_centers,
        rtol=0.0,
        atol=64.0 * np.finfo(float).eps * max(1.0, abs(float(faces[-1]))),
    ):
        raise ValueError("grid must be uniform cell-centred spherical geometry")
    if (
        grid_volumes.shape != (n,)
        or not np.all(np.isfinite(grid_volumes))
        or np.any(grid_volumes <= 0.0)
    ):
        raise ValueError("grid spherical shell volumes must be finite and positive")

    step = float(faces[1] - faces[0])
    source_face = int(round(radius_value / step))
    tolerance = (
        64.0
        * np.finfo(float).eps
        * max(1.0, abs(radius_value), abs(float(faces[-1])))
    )
    if (
        source_face < 1
        or source_face >= len(faces)
        or abs(float(faces[source_face]) - radius_value) > tolerance
    ):
        raise ValueError("source_radius must align with an interior cell face")

    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    half_widths = 0.5 * np.diff(faces)
    midpoints = 0.5 * (faces[:-1] + faces[1:])
    quadrature_r = midpoints[:, None] + half_widths[:, None] * nodes[None, :]
    quadrature_sqrt_gamma, _, _, _, _ = _metric_at(
        quadrature_r, centers, a_cells, b_cells, X_cells
    )
    if np.any(quadrature_r < radius_value):
        source_active = quadrature_r < radius_value
        rho = np.zeros_like(quadrature_r)
        scaled_r = quadrature_r[source_active] / radius_value
        rho[source_active] = rho0_value * (1.0 - scaled_r**2) ** 2
    else:
        rho = np.zeros_like(quadrature_r)
    charge_integrand = FOUR_PI * rho * quadrature_sqrt_gamma
    cell_charges = half_widths * np.sum(
        weights[None, :] * charge_integrand, axis=1
    )
    if not np.all(np.isfinite(cell_charges)):
        raise ValueError("proper-volume cell charges are non-finite")

    face_sqrt_gamma, face_gamma_rr, face_areal_radius, _, _ = _metric_at(
        faces, centers, a_cells, b_cells, X_cells
    )
    face_charges = np.concatenate(([0.0], np.cumsum(cell_charges)))
    face_electric_field_r = np.zeros_like(faces)
    np.divide(
        face_charges,
        face_sqrt_gamma,
        out=face_electric_field_r,
        where=faces > 0.0,
    )
    face_electric_field_hat = np.sqrt(face_gamma_rr) * face_electric_field_r
    if (
        not np.all(np.isfinite(face_charges))
        or not np.all(np.isfinite(face_electric_field_r))
        or not np.all(np.isfinite(face_electric_field_hat))
    ):
        raise ValueError("finite-volume charge/field reconstruction is non-finite")

    face_flux = face_sqrt_gamma * face_electric_field_r
    cell_residuals = np.diff(face_flux) - cell_charges
    total_charge = float(np.sum(cell_charges))
    midpoint_factor = np.sqrt(a_cells) * b_cells / X_cells**3
    rho_centers = np.zeros_like(centers)
    active_centers = centers < radius_value
    scaled_centers = centers[active_centers] / radius_value
    rho_centers[active_centers] = rho0_value * (1.0 - scaled_centers**2) ** 2
    midpoint_total_charge = float(
        FOUR_PI
        * np.sum(rho_centers * midpoint_factor * grid_volumes)
    )
    if (
        not math.isfinite(total_charge)
        or not math.isfinite(midpoint_total_charge)
        or not np.all(np.isfinite(cell_residuals))
    ):
        raise ValueError("Gauss diagnostic produced non-finite charge/residual")

    return FrozenProductionGaussResult(
        faces=faces.copy(),
        cell_charges=cell_charges,
        face_charges=face_charges,
        face_sqrt_gamma=face_sqrt_gamma,
        face_gamma_rr=face_gamma_rr,
        face_areal_radius=face_areal_radius,
        face_electric_field_r=face_electric_field_r,
        face_electric_field_hat=face_electric_field_hat,
        cell_residuals=cell_residuals,
        total_charge=total_charge,
        midpoint_total_charge=midpoint_total_charge,
        max_abs_cell_residual=float(np.max(np.abs(cell_residuals))),
    )
