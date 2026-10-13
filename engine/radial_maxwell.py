"""Feature-gated spherical Einstein--Maxwell state for a separate charged-dust sector.

Convention: Gaussian units with c=G=1, signature (-,+,+,+), and a radial
electric field only. Exact spherical symmetry has no propagating transverse
Maxwell degree of freedom: the field is reconstructed from face-enclosed
charge by Gauss' law. The evolved electromagnetic degree of freedom is Q on
radial faces, advanced from the *same numerical Valencia rest-mass face flux*
as a separately specified dust carrier with constant caller-supplied q/m.

This module is not enabled by default and does not assign charge to any
existing GEAR species. The caller must create the additional charged-dust
state and explicitly supply q/m. A successful numerical run of this extension
does not prove that the original GEAR archive derives this carrier or action.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    Species,
    _metric_arrays,
    primitives,
)


FOUR_PI = 4.0 * math.pi


def _finite_real_array(values, name: str, *, length: int | None = None) -> np.ndarray:
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real one-dimensional array") from exc
    if raw.ndim != 1 or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a finite real one-dimensional array")
    out = np.array(raw, dtype=float, copy=True)
    if out.size < 2 or (length is not None and out.size != length):
        expected = "at least two" if length is None else f"exactly {length}"
        raise ValueError(f"{name} must contain {expected} values")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _finite_real_scalar(value, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise ValueError(f"{name} must be a finite real scalar")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be a finite real scalar")
    return out


@dataclass(frozen=True)
class RadialMaxwellState:
    """Face-enclosed charge for an opt-in radial electric sector."""

    charge_to_rest_mass: float
    enclosed_charge_faces: np.ndarray

    def __post_init__(self) -> None:
        eta = _finite_real_scalar(self.charge_to_rest_mass, "charge_to_rest_mass")
        q = _finite_real_array(self.enclosed_charge_faces, "enclosed_charge_faces")
        if q[0] != 0.0:
            raise ValueError("enclosed charge at the regular centre must be zero")
        q.setflags(write=False)
        object.__setattr__(self, "charge_to_rest_mass", eta)
        object.__setattr__(self, "enclosed_charge_faces", q)

    @classmethod
    def from_charged_dust(
        cls,
        grid,
        charged_dust: ConservedSpecies,
        charge_to_rest_mass: float,
    ) -> "RadialMaxwellState":
        """Initialize Gauss-consistent face charge from one dust rest state.

        The per-solid-angle rest variable is sqrt(gamma)*D. For constant q/m,
        each shell's charge is 4*pi*(q/m)*rest*dr, matching the numerical
        Valencia rest-mass continuity variable and the declared Gaussian units.
        """
        eta = _finite_real_scalar(charge_to_rest_mass, "charge_to_rest_mass")
        rest = _finite_real_array(charged_dust.rest, "charged_dust.rest")
        if rest.size != int(grid.n):
            raise ValueError("charged dust rest state must match grid cell count")
        if np.any(rest < 0.0):
            raise ValueError("charged dust rest state must be nonnegative")
        dr = _finite_real_scalar(grid.dr, "grid.dr")
        if dr <= 0.0:
            raise ValueError("grid.dr must be positive")
        q_faces = np.empty(rest.size + 1, dtype=float)
        q_faces[0] = 0.0
        q_faces[1:] = np.cumsum(FOUR_PI * eta * dr * rest)
        if not np.all(np.isfinite(q_faces)):
            raise ValueError("initial enclosed charge is non-finite")
        return cls(eta, q_faces)

    def copy(self) -> "RadialMaxwellState":
        return RadialMaxwellState(
            self.charge_to_rest_mass, self.enclosed_charge_faces.copy()
        )

    def charge_flux_faces(self, rest_mass_flux_faces) -> np.ndarray:
        """Use the exact same-stage Valencia face flux for charge continuity."""
        flux = _finite_real_array(
            rest_mass_flux_faces, "rest_mass_flux_faces",
            length=self.enclosed_charge_faces.size,
        )
        if flux[0] != 0.0:
            raise ValueError("regular-centre Valencia rest-mass flux must be zero")
        result = FOUR_PI * self.charge_to_rest_mass * flux
        if not np.all(np.isfinite(result)):
            raise ValueError("charge-current face flux is non-finite")
        return result

    def advance(self, dt: float, rest_mass_flux_faces) -> "RadialMaxwellState":
        if isinstance(dt, (bool, np.bool_)):
            raise ValueError("dt must be finite and positive")
        try:
            step = float(dt)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("dt must be finite and positive") from exc
        if not math.isfinite(step) or step <= 0.0:
            raise ValueError("dt must be finite and positive")
        current = self.charge_flux_faces(rest_mass_flux_faces)
        updated = self.enclosed_charge_faces - step * current
        if not np.all(np.isfinite(updated)):
            raise ValueError("radial Maxwell charge update produced non-finite values")
        updated[0] = 0.0
        return RadialMaxwellState(self.charge_to_rest_mass, updated)

    def areal_radius_faces(self, grid, geometry) -> np.ndarray:
        """Interpolate the live areal radius R=r*sqrt(b)/X to radial faces."""
        r_faces = _finite_real_array(grid.faces, "grid.faces", length=self.enclosed_charge_faces.size)
        b = _finite_real_array(geometry.b, "geometry.b", length=grid.n)
        X = _finite_real_array(geometry.X, "geometry.X", length=grid.n)
        if np.any(b <= 0.0) or np.any(X <= 0.0):
            raise ValueError("metric b and X must be strictly positive")
        bf = np.empty(grid.n + 1, dtype=float)
        Xf = np.empty(grid.n + 1, dtype=float)
        bf[0], Xf[0] = b[0], X[0]
        if grid.n > 1:
            bf[1:-1] = 0.5 * (b[:-1] + b[1:])
            Xf[1:-1] = 0.5 * (X[:-1] + X[1:])
        bf[-1], Xf[-1] = b[-1], X[-1]
        result = r_faces * np.sqrt(bf) / Xf
        if result[0] != 0.0 or not np.all(np.isfinite(result)):
            raise ValueError("areal-radius faces violate centre/finiteness invariant")
        if np.any(np.diff(result) <= 0.0):
            raise ValueError("areal radius must be strictly increasing on the face grid")
        return result

    def electric_field_hat(self, grid, geometry) -> np.ndarray:
        """Return cell-centred orthonormal Ehat_r from the face Gauss state."""
        radius = _finite_real_array(grid.centers, "grid.centers", length=grid.n)
        q_centres = 0.5 * (
            self.enclosed_charge_faces[:-1] + self.enclosed_charge_faces[1:]
        )
        areal_centres = radius * np.sqrt(
            _finite_real_array(geometry.b, "geometry.b", length=grid.n)
        ) / _finite_real_array(geometry.X, "geometry.X", length=grid.n)
        if np.any(areal_centres <= 0.0):
            raise ValueError("cell-centred areal radius must be positive")
        field = q_centres / np.square(areal_centres)
        if not np.all(np.isfinite(field)):
            raise ValueError("Gauss reconstruction produced non-finite electric field")
        return field

    def stress_projection(self, grid, geometry) -> dict[str, np.ndarray]:
        """Radial pure-electric stress in Gaussian units: (u,-u,+u,0)."""
        Ehat = self.electric_field_hat(grid, geometry)
        energy = np.square(Ehat) / (8.0 * math.pi)
        return {
            "rho": energy,
            "pr": -energy,
            "pt": energy,
            "j": np.zeros_like(energy),
        }

    def lorentz_source(
        self,
        grid,
        geometry,
        metric: BSSNMetricSlice,
        charged_dust: ConservedSpecies,
    ) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
        """Return Eulerian Valencia energy/momentum sources for the charged dust.

        E_r = sqrt(gamma_rr)*Ehat_r, rho_e=(q/m)*rho0*W and j_e^r=rho_e*v^r.
        The matter receives (+alpha*sqrt(gamma)*E_r*j^r,
        +alpha*sqrt(gamma)*rho_e*E_r). The field charge state is advanced from
        the same dust's *numerical* rest-mass face flux in the PIRK stage.
        """
        if self.enclosed_charge_faces.size != grid.n + 1:
            raise ValueError("face-charge state must have grid.n + 1 entries")
        metrics = _metric_arrays(metric)
        primitives_by_cell = primitives(metrics, charged_dust, Species.BARYON)
        Ehat = self.electric_field_hat(grid, geometry)
        root_gamma_rr = np.sqrt(np.asarray([m.gamma_rr for m in metrics]))
        E_r = root_gamma_rr * Ehat
        rho0 = np.asarray([q.rho for q in primitives_by_cell], dtype=float)
        velocity = np.asarray([q.v_r for q in primitives_by_cell], dtype=float)
        W = np.asarray([q.lorentz() for q in primitives_by_cell], dtype=float)
        rho_e = self.charge_to_rest_mass * rho0 * W
        current_r = rho_e * velocity
        alpha = np.asarray([m.alpha for m in metrics], dtype=float)
        sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics], dtype=float)
        source_energy = alpha * sqrt_gamma * E_r * current_r
        source_momentum = alpha * sqrt_gamma * E_r * rho_e
        values = {
            "Ehat_r": Ehat,
            "E_r": E_r,
            "rho_e": rho_e,
            "j_e_r": current_r,
            "source_energy": source_energy,
            "source_momentum_r": source_momentum,
        }
        if not all(np.all(np.isfinite(v)) for v in values.values()):
            raise ValueError("Lorentz source produced non-finite stage data")
        return source_energy, source_momentum, values

    def gauss_shell_residual(
        self, grid, charged_dust: ConservedSpecies
    ) -> np.ndarray:
        """Discrete Gauss/charge-to-mass invariant, shell by shell."""
        rest = _finite_real_array(charged_dust.rest, "charged_dust.rest", length=grid.n)
        dr = _finite_real_scalar(grid.dr, "grid.dr")
        return np.diff(self.enclosed_charge_faces) - (
            FOUR_PI * self.charge_to_rest_mass * dr * rest
        )
