"""Source-agnostic radial Maxwell state/RHS primitives (candidate only).

This module is intentionally not imported by ProductionState or the active
V5.5 stepper. It implements the standard spherical Gauss/charge-continuity
identities and conditional Eulerian Valencia Lorentz-source projection, but
does not choose a charge carrier, current-closure law, or field normalization
for the 0-Cosmos/GEAR model.

Conventions: Gaussian units, c=1, signature (-,+,+,+), Eulerian charge density
rho_e, Eulerian contravariant radial current j_e^r, areal radius R, and spatial
volume density sqrt(gamma) including the spherical r**2 factor but excluding
sin(theta). See docs/EM_LITERATURE_AND_MODEL_GATES_2026-10-13.md.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

FOUR_PI = 4.0 * math.pi
EIGHT_PI = 8.0 * math.pi


def _array(values, name: str, *, ndim: int = 1) -> np.ndarray:
    """Copy a finite real numeric array; reject bools, complex and object data."""
    raw = np.asarray(values)
    if raw.dtype.kind not in "iuf" or raw.ndim != ndim:
        raise ValueError(f"{name} must be a {ndim}-D real numeric array")
    result = np.array(raw, dtype=float, copy=True)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _same_shape(named_arrays: dict[str, np.ndarray]) -> None:
    shapes = {value.shape for value in named_arrays.values()}
    if len(shapes) != 1:
        details = ", ".join(f"{name}={value.shape}" for name, value in named_arrays.items())
        raise ValueError(f"array shapes must agree: {details}")


@dataclass(frozen=True)
class RadialMaxwellState:
    """Finite-volume enclosed-charge state on N radial cells.

    cell_charge[i] is the total electric charge in shell i, including the 4*pi
    angular integral. It is the conservative variable advanced by the discrete
    charge-continuity equation; it is not a pointwise charge density.
    """

    cell_charge: np.ndarray

    def __post_init__(self) -> None:
        charge = _array(self.cell_charge, "cell_charge")
        if charge.size < 1:
            raise ValueError("cell_charge must contain at least one radial cell")
        charge.setflags(write=False)
        object.__setattr__(self, "cell_charge", charge)

    @property
    def enclosed_charge_faces(self) -> np.ndarray:
        """Return Q at N+1 faces, with Q(0)=0 and Q[k]=sum(cell_charge[:k])."""
        result = np.empty(self.cell_charge.size + 1, dtype=float)
        result[0] = 0.0
        total = 0.0
        compensation = 0.0
        for i, charge in enumerate(self.cell_charge, start=1):
            corrected = float(charge) - compensation
            updated = total + corrected
            if not math.isfinite(updated):
                raise ValueError("cumulative enclosed charge must remain finite")
            compensation = (updated - total) - corrected
            if not math.isfinite(compensation):
                raise ValueError("cumulative charge compensation must remain finite")
            total = updated
            result[i] = total
        return result

    @classmethod
    def from_cell_density(
        cls, charge_density: np.ndarray, proper_shell_volume: np.ndarray
    ) -> "RadialMaxwellState":
        """Build shell charge from cell-averaged rho_e and proper shell volume.

        proper_shell_volume is the complete spherical cell volume, including
        the 4*pi angular integral. The caller is responsible for computing it
        from the same stage's metric.
        """
        rho = _array(charge_density, "charge_density")
        volume = _array(proper_shell_volume, "proper_shell_volume")
        _same_shape({"charge_density": rho, "proper_shell_volume": volume})
        if np.any(volume <= 0.0):
            raise ValueError("proper_shell_volume must be positive in every cell")
        charge = rho * volume
        if not np.all(np.isfinite(charge)):
            raise ValueError("cell charge must remain finite")
        return cls(charge)


@dataclass(frozen=True)
class RadialElectricStress:
    """Eulerian orthonormal stress projections for a purely radial electric field."""

    rho: np.ndarray
    pr: np.ndarray
    pt: np.ndarray
    j: np.ndarray


@dataclass(frozen=True)
class ValenciaLorentzSources:
    """Matter-side increments for active Eulerian energy and covariant momentum."""

    q_energy: np.ndarray
    q_momentum_r: np.ndarray


def radial_electric_field_faces(
    enclosed_charge_faces: np.ndarray,
    areal_radius_faces: np.ndarray,
) -> np.ndarray:
    """Reconstruct Ehat_r=Q/R**2 from charge and areal radius at radial faces.

    The first face must be the regular centre (R=Q=0); noncentral areal radii
    must be strictly positive and strictly increasing. A point charge at the
    centre is deliberately not represented by this regular-source state.
    """
    charge = _array(enclosed_charge_faces, "enclosed_charge_faces")
    radius = _array(areal_radius_faces, "areal_radius_faces")
    _same_shape({"enclosed_charge_faces": charge, "areal_radius_faces": radius})
    if charge.size < 2:
        raise ValueError("at least two radial faces are required")
    if radius[0] != 0.0 or charge[0] != 0.0:
        raise ValueError("the regular centre requires first face R=0 and Q=0")
    if np.any(radius[1:] <= 0.0) or np.any(np.diff(radius) <= 0.0):
        raise ValueError("noncentral areal radii must be positive and increasing")
    field = np.zeros_like(charge)
    field[1:] = charge[1:] / np.square(radius[1:])
    if not np.all(np.isfinite(field)):
        raise ValueError("radial electric field must remain finite")
    return field


def charge_current_flux_faces(
    sqrt_gamma_faces: np.ndarray,
    alpha_faces: np.ndarray,
    beta_r_faces: np.ndarray,
    rho_e_faces: np.ndarray,
    j_e_r_faces: np.ndarray,
) -> np.ndarray:
    """Return angular-integrated radial charge flux at N+1 faces.

    With J^mu=rho_e*n^mu+j^mu and n^mu=(1/alpha,-beta^i/alpha), charge
    continuity is
      d_t(sqrt(gamma)*rho_e)
        + d_r[sqrt(gamma)*(alpha*j_e^r-beta^r*rho_e)] = 0.
    The returned flux includes 4*pi because cell_charge is the full shell
    charge. Inputs must all belong to the same integrator stage and face grid.
    """
    values = {
        "sqrt_gamma_faces": _array(sqrt_gamma_faces, "sqrt_gamma_faces"),
        "alpha_faces": _array(alpha_faces, "alpha_faces"),
        "beta_r_faces": _array(beta_r_faces, "beta_r_faces"),
        "rho_e_faces": _array(rho_e_faces, "rho_e_faces"),
        "j_e_r_faces": _array(j_e_r_faces, "j_e_r_faces"),
    }
    _same_shape(values)
    rootg = values["sqrt_gamma_faces"]
    alpha = values["alpha_faces"]
    beta = values["beta_r_faces"]
    rho = values["rho_e_faces"]
    current = values["j_e_r_faces"]
    if np.any(rootg < 0.0):
        raise ValueError("sqrt_gamma_faces must be nonnegative")
    if np.any(alpha <= 0.0):
        raise ValueError("alpha_faces must be positive")
    flux = FOUR_PI * rootg * (alpha * current - beta * rho)
    if not np.all(np.isfinite(flux)):
        raise ValueError("radial charge flux must remain finite")
    return flux


def charge_continuity_rhs(charge_flux_faces: np.ndarray) -> np.ndarray:
    """Return conservative d(cell_charge)/dt = -(F[i+1/2]-F[i-1/2])."""
    flux = _array(charge_flux_faces, "charge_flux_faces")
    if flux.size < 2:
        raise ValueError("charge_flux_faces must have at least two faces")
    rhs = -(flux[1:] - flux[:-1])
    if not np.all(np.isfinite(rhs)):
        raise ValueError("charge-continuity RHS must remain finite")
    return rhs


def radial_electric_stress_projection(
    electric_field_hat_cells: np.ndarray,
) -> RadialElectricStress:
    """Project cell-sampled pure radial E to (rho, pr, pt, j) in Gaussian units.

    This function assumes zero magnetic field and zero Poynting flux. It must
    not be used for a general electromagnetic configuration.
    """
    field = _array(electric_field_hat_cells, "electric_field_hat_cells")
    rho = np.square(field) / EIGHT_PI
    if not np.all(np.isfinite(rho)):
        raise ValueError("electromagnetic energy density must remain finite")
    return RadialElectricStress(
        rho=rho,
        pr=-rho.copy(),
        pt=rho.copy(),
        j=np.zeros_like(rho),
    )


def valencia_lorentz_sources(
    alpha: np.ndarray,
    beta_r: np.ndarray,
    sqrt_gamma: np.ndarray,
    electric_covector_r: np.ndarray,
    rho_e: np.ndarray,
    j_e_r: np.ndarray,
) -> ValenciaLorentzSources:
    """Map the conventional radial Lorentz force into active Valencia sources.

    Conditional on the convention documented in the EM literature note, the
    covector is Q_nu=F_nu_lambda J^lambda, with
      Q_r=rho_e*E_r,
      Q_t=-alpha*E_r*j_e^r+beta^r*rho_e*E_r.
    For active Eulerian energy U_E=sqrt(gamma)*E and radial momentum
    U_Sr=sqrt(gamma)*S_r, the charged-matter source is
      q_E=-sqrt(gamma)*(Q_t-beta^r*Q_r),
      q_Sr=alpha*sqrt(gamma)*Q_r.
    This is an algebraic source map only: the caller must have an approved
    current law and inject equal-and-opposite EM/matter exchange at the same
    real integrator stage. It does not choose which species carries charge.
    """
    values = {
        "alpha": _array(alpha, "alpha"),
        "beta_r": _array(beta_r, "beta_r"),
        "sqrt_gamma": _array(sqrt_gamma, "sqrt_gamma"),
        "electric_covector_r": _array(electric_covector_r, "electric_covector_r"),
        "rho_e": _array(rho_e, "rho_e"),
        "j_e_r": _array(j_e_r, "j_e_r"),
    }
    _same_shape(values)
    if np.any(values["alpha"] <= 0.0):
        raise ValueError("alpha must be positive")
    if np.any(values["sqrt_gamma"] < 0.0):
        raise ValueError("sqrt_gamma must be nonnegative")
    alpha_a = values["alpha"]
    beta = values["beta_r"]
    rootg = values["sqrt_gamma"]
    e_r = values["electric_covector_r"]
    rho = values["rho_e"]
    current = values["j_e_r"]
    q_r = rho * e_r
    q_t = -alpha_a * e_r * current + beta * rho * e_r
    q_energy = -rootg * (q_t - beta * q_r)
    q_momentum = alpha_a * rootg * q_r
    if not (
        np.all(np.isfinite(q_energy))
        and np.all(np.isfinite(q_momentum))
    ):
        raise ValueError("Valencia Lorentz source must remain finite")
    return ValenciaLorentzSources(q_energy=q_energy, q_momentum_r=q_momentum)
