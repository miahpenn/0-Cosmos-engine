"""Finite-volume Gauss diagnostic on a frozen production-initialized metric.

Research-only: this snapshots the actual V5.5 initializer's spherical metric,
then solves a manufactured compact charge profile on that fixed background.
It does not add charge to production matter, evolve Maxwell fields, or modify
the production state/source graph.

The metric samples are cell-centred. We interpolate a,b,X piecewise-linearly
between their frozen samples, use even reflection at the regular centre, and
linear extrapolation only over the outer half-cell. Any nonfinite/nonpositive
interpolated metric value fails closed. These are diagnostic interpolation
rules, not production geometry/evolution changes.
"""
from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
import math

import numpy as np

from tools.research_em_radial_reference import RadialCoulombReference


_GAUSS8_X, _GAUSS8_W = np.polynomial.legendre.leggauss(8)
_GAUSS16_X, _GAUSS16_W = np.polynomial.legendre.leggauss(16)


def _finite_real(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


@dataclass(frozen=True)
class FrozenProductionMetric:
    """Immutable metric sample captured from a real initialized V5.5 state."""

    centers: tuple[float, ...]
    a: tuple[float, ...]
    b: tuple[float, ...]
    X: tuple[float, ...]
    dr: float
    r_max: float

    def __post_init__(self) -> None:
        centers = tuple(_finite_real(v, "center") for v in self.centers)
        a = tuple(_finite_real(v, "a") for v in self.a)
        b = tuple(_finite_real(v, "b") for v in self.b)
        X = tuple(_finite_real(v, "X") for v in self.X)
        dr = _finite_real(self.dr, "dr")
        r_max = _finite_real(self.r_max, "r_max")
        n = len(centers)
        if n < 3 or len(a) != n or len(b) != n or len(X) != n:
            raise ValueError("metric samples must have matching length >= 3")
        if dr <= 0.0 or r_max <= 0.0:
            raise ValueError("dr and r_max must be positive")
        if any(v <= 0.0 for v in a + b + X):
            raise ValueError("a, b, and X must be finite and positive")
        # Scale tolerance to the actual coordinate domain, not to unity.
        # A unit-sized absolute tolerance would accept malformed tiny grids.
        scale = max(abs(r_max), abs(dr * n), abs(dr))
        tol = 64.0 * math.ulp(scale)
        if abs(r_max - dr * n) > tol:
            raise ValueError("r_max must equal resolution * dr")
        for i, center in enumerate(centers):
            expected = (i + 0.5) * dr
            if abs(center - expected) > tol:
                raise ValueError("production grid must be uniform and cell-centred")
            if i and not center > centers[i - 1]:
                raise ValueError("metric centers must be strictly increasing")
        object.__setattr__(self, "centers", centers)
        object.__setattr__(self, "a", a)
        object.__setattr__(self, "b", b)
        object.__setattr__(self, "X", X)
        object.__setattr__(self, "dr", dr)
        object.__setattr__(self, "r_max", r_max)

    @classmethod
    def from_production_state(cls, state) -> "FrozenProductionMetric":
        """Copy the live initializer's a,b,X arrays; later state changes cannot affect it."""
        try:
            centers = np.asarray(state.grid.centers, dtype=float)
            dr = float(state.grid.dr)
            a = np.asarray(state.geometry.a, dtype=float)
            b = np.asarray(state.geometry.b, dtype=float)
            X = np.asarray(state.geometry.X, dtype=float)
        except (AttributeError, TypeError, ValueError) as exc:
            raise ValueError("state must expose grid centers/dr and geometry a,b,X") from exc
        if any(values.ndim != 1 for values in (centers, a, b, X)):
            raise ValueError("grid and metric arrays must be one-dimensional")
        if len(centers) == 0:
            raise ValueError("production metric has no cells")
        return cls(
            centers=tuple(centers),
            a=tuple(a),
            b=tuple(b),
            X=tuple(X),
            dr=dr,
            r_max=float(len(centers) * dr),
        )

    def _radius(self, radius: float) -> float:
        value = _finite_real(radius, "radius")
        if value < 0.0 or value > self.r_max:
            raise ValueError("radius must lie inside the frozen metric domain")
        return value

    def _interpolate(self, radius: float, values: tuple[float, ...]) -> float:
        r = self._radius(radius)
        # Metric fields are regular even fields. Reflection maps the first
        # half-cell to equal samples at +/- centers[0], giving a parity-respecting
        # constant interpolation there without inventing a centre value.
        if r <= self.centers[0]:
            value = values[0]
        elif r >= self.centers[-1]:
            # Only the last half-cell is extrapolated; reject if that linear
            # continuation yields an invalid metric value.
            slope = (values[-1] - values[-2]) / (
                self.centers[-1] - self.centers[-2]
            )
            value = values[-1] + slope * (r - self.centers[-1])
        else:
            i = bisect_right(self.centers, r) - 1
            fraction = (r - self.centers[i]) / (
                self.centers[i + 1] - self.centers[i]
            )
            value = values[i] + fraction * (values[i + 1] - values[i])
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(
                f"interpolated metric value must be finite and positive at r={r:.17e}"
            )
        return value

    def metric_values(self, radius: float) -> tuple[float, float, float]:
        r = self._radius(radius)
        return (
            self._interpolate(r, self.a),
            self._interpolate(r, self.b),
            self._interpolate(r, self.X),
        )

    def sqrt_gamma(self, radius: float) -> float:
        """Return sqrt(gamma)=sqrt(a)*b*r^2/X^3, without angular sin(theta)."""
        r = self._radius(radius)
        if r == 0.0:
            return 0.0
        a, b, X = self.metric_values(r)
        value = math.sqrt(a) * b * r * r / (X**3)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"invalid proper-volume factor at r={r:.17e}")
        return value

    def areal_radius(self, radius: float) -> float:
        """Return the production spherical areal radius R=r*sqrt(b)/X."""
        r = self._radius(radius)
        if r == 0.0:
            return 0.0
        _, b, X = self.metric_values(r)
        value = r * math.sqrt(b) / X
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"invalid areal radius at r={r:.17e}")
        return value

    def _breakpoints(self, lower: float, upper: float) -> tuple[float, ...]:
        lo = self._radius(lower)
        hi = self._radius(upper)
        if hi < lo:
            raise ValueError("upper integration bound must not precede lower bound")
        if hi == lo:
            return (lo, hi)
        first = bisect_right(self.centers, lo)
        last = bisect_right(self.centers, hi)
        interior = tuple(c for c in self.centers[first:last] if lo < c < hi)
        return (lo, *interior, hi)

    @staticmethod
    def _gauss_integral(function, lower: float, upper: float, nodes, weights) -> float:
        midpoint = 0.5 * (lower + upper)
        half_width = 0.5 * (upper - lower)
        result = half_width * math.fsum(
            float(weight) * function(midpoint + half_width * float(node))
            for node, weight in zip(nodes, weights)
        )
        if not math.isfinite(result):
            raise ValueError("nonfinite quadrature result")
        return result

    def _integrate_piecewise(self, function, lower: float, upper: float, *, order: int) -> float:
        breaks = self._breakpoints(lower, upper)
        nodes, weights = (_GAUSS8_X, _GAUSS8_W) if order == 8 else (_GAUSS16_X, _GAUSS16_W)
        pieces = (
            self._gauss_integral(function, left, right, nodes, weights)
            for left, right in zip(breaks, breaks[1:])
        )
        result = math.fsum(pieces)
        if not math.isfinite(result):
            raise ValueError("nonfinite piecewise integral")
        return result

    def proper_shell_volume(self, lower: float, upper: float) -> float:
        """Independently integrate 4*pi*sqrt(gamma) dr over one target shell."""
        lo = self._radius(lower)
        hi = self._radius(upper)
        if hi < lo:
            raise ValueError("upper shell face must not precede lower face")
        if hi == lo:
            return 0.0
        volume = 4.0 * math.pi * self._integrate_piecewise(
            self.sqrt_gamma, lo, hi, order=8
        )
        if not math.isfinite(volume) or volume <= 0.0:
            raise ValueError("proper shell volume must be finite and positive")
        return volume

    def proper_shell_charge(
        self,
        source: RadialCoulombReference,
        lower: float,
        upper: float,
    ) -> float:
        """Integrate 4*pi*rho_e*sqrt(gamma) over one target cell.

        The manufactured source is integrated together with the metric's
        proper-volume factor instead of sampling rho_e only at the cell centre.
        Piecewise integration splits at metric sample centres, where the linear
        interpolant's slope may change. Source support is clipped at R0; the
        solver separately requires R0 to coincide with a target face.
        """
        if not isinstance(source, RadialCoulombReference):
            raise ValueError("source must be a RadialCoulombReference")
        lo = self._radius(lower)
        hi = self._radius(upper)
        if hi < lo:
            raise ValueError("upper shell face must not precede lower face")
        support_hi = min(hi, source.source_radius)
        if support_hi <= lo:
            return 0.0

        def integrand(radius: float) -> float:
            return FOUR_PI * source.charge_density(radius) * self.sqrt_gamma(radius)

        charge = self._integrate_piecewise(integrand, lo, support_hi, order=8)
        if not math.isfinite(charge):
            raise ValueError("proper-volume shell charge must be finite")
        return charge

    def reference_enclosed_charge(
        self, source: RadialCoulombReference, radius: float
    ) -> float:
        """High-order proper-volume reference integral on this same frozen metric.

        This intentionally does not call source.enclosed_charge(): the latter
        is a flat-space oracle and is not valid on this production geometry.
        """
        if not isinstance(source, RadialCoulombReference):
            raise ValueError("source must be a RadialCoulombReference")
        r = self._radius(radius)
        upper = min(r, source.source_radius)
        if upper == 0.0:
            return 0.0

        def integrand(s: float) -> float:
            return 4.0 * math.pi * source.charge_density(s) * self.sqrt_gamma(s)

        return self._integrate_piecewise(integrand, 0.0, upper, order=16)


@dataclass(frozen=True)
class ProductionMetricGaussSolution:
    """Static cell-source / face-field solution over the frozen production metric."""

    faces: tuple[float, ...]
    cell_centers: tuple[float, ...]
    cell_charges: tuple[float, ...]
    face_charges: tuple[float, ...]
    face_electric_field_hat: tuple[float, ...]
    face_electric_field_r: tuple[float, ...]
    cell_residuals: tuple[float, ...]

    @property
    def total_charge(self) -> float:
        return self.face_charges[-1]

    @property
    def max_abs_cell_residual(self) -> float:
        return max((abs(value) for value in self.cell_residuals), default=0.0)


def solve_discrete_production_metric_gauss(
    metric: FrozenProductionMetric,
    source: RadialCoulombReference,
    *,
    resolution: int,
) -> ProductionMetricGaussSolution:
    """Solve the static radial Gauss constraint on a frozen real production metric.

    The source density is sampled at target cell centres; proper shell volumes
    are separately integrated over the frozen production metric. Face charge is
    accumulated from those source cells, never from the flat or curved analytic
    enclosed-charge functions. The coordinate field follows
    E^r=Q_h/sqrt(gamma), while E_hat_r=Q_h/R_areal^2. This remains a static
    diagnostic, not field time evolution or Einstein-source coupling.
    """
    if not isinstance(metric, FrozenProductionMetric):
        raise ValueError("metric must be a FrozenProductionMetric")
    if not isinstance(source, RadialCoulombReference):
        raise ValueError("source must be a RadialCoulombReference")
    if isinstance(resolution, bool) or not isinstance(resolution, int) or resolution < 2:
        raise ValueError("resolution must be an integer >= 2")
    dr = metric.r_max / resolution
    if not math.isfinite(dr) or dr <= 0.0:
        raise ValueError("invalid target grid spacing")
    if source.source_radius >= metric.r_max:
        raise ValueError("source_radius must lie inside the frozen metric domain")
    face_coordinate = source.source_radius / dr
    if not math.isfinite(face_coordinate):
        raise ValueError("source_radius must align with an interior target face")
    source_face = round(face_coordinate)
    nearest_face_radius = source_face * dr
    align_error = abs(nearest_face_radius - source.source_radius)
    # Only roundoff at the actual coordinate scale counts as face alignment.
    align_tol = 64.0 * max(
        math.ulp(source.source_radius),
        math.ulp(nearest_face_radius),
    )
    if not (
        1 <= source_face < resolution
        and math.isfinite(nearest_face_radius)
        and math.isfinite(align_tol)
        and align_error <= align_tol
    ):
        raise ValueError(
            "source_radius must align with an interior target face; "
            f"source_radius={source.source_radius:.17e}, dr={dr:.17e}"
        )

    faces = tuple(i * dr for i in range(resolution + 1))
    centers = tuple(0.5 * (faces[i] + faces[i + 1]) for i in range(resolution))
    cell_charges = tuple(
        metric.proper_shell_charge(source, left, right)
        for left, right in zip(faces, faces[1:])
    )
    if not all(math.isfinite(value) for value in cell_charges):
        raise ValueError("finite-volume cell charges must remain finite")

    face_charges = [0.0]
    running = 0.0
    compensation = 0.0
    for charge in cell_charges:
        corrected = charge - compensation
        updated = running + corrected
        if not (math.isfinite(corrected) and math.isfinite(updated)):
            raise ValueError("cumulative finite-volume charge must remain finite")
        next_compensation = (updated - running) - corrected
        if not math.isfinite(next_compensation):
            raise ValueError("finite-volume charge compensation must remain finite")
        compensation = next_compensation
        running = updated
        face_charges.append(running)

    fields_hat = []
    fields_r = []
    for radius, charge in zip(faces, face_charges):
        if radius == 0.0:
            fields_hat.append(0.0)
            fields_r.append(0.0)
            continue
        areal = metric.areal_radius(radius)
        sqrt_gamma = metric.sqrt_gamma(radius)
        areal_squared = areal * areal
        if not (
            math.isfinite(areal_squared)
            and areal_squared > 0.0
            and math.isfinite(sqrt_gamma)
            and sqrt_gamma > 0.0
        ):
            raise ValueError("production-metric face geometry must remain finite and positive")
        field_hat = charge / areal_squared
        field_r = charge / sqrt_gamma
        if not (math.isfinite(field_hat) and math.isfinite(field_r)):
            raise ValueError("reconstructed production-metric electric field must remain finite")
        fields_hat.append(field_hat)
        fields_r.append(field_r)

    fluxes = tuple(
        metric.sqrt_gamma(radius) * field
        for radius, field in zip(faces, fields_r)
    )
    if not all(math.isfinite(value) for value in fluxes):
        raise ValueError("production-metric Gauss fluxes must remain finite")
    residuals = tuple(
        fluxes[i + 1] - fluxes[i] - cell_charges[i]
        for i in range(resolution)
    )
    if any(not math.isfinite(value) for value in residuals):
        raise ValueError("production-metric Gauss residuals must remain finite")
    return ProductionMetricGaussSolution(
        faces=faces,
        cell_centers=centers,
        cell_charges=cell_charges,
        face_charges=tuple(face_charges),
        face_electric_field_hat=tuple(fields_hat),
        face_electric_field_r=tuple(fields_r),
        cell_residuals=residuals,
    )
