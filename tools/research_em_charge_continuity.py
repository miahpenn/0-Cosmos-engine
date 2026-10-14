"""Research-only conservative charge-continuity kernel for spherical Maxwell.

This module advances face-enclosed charge from a supplied, same-stage,
sphere-integrated outward current flux. It is not imported by the production
kernel and does not prescribe a carrier, current constitutive law, or
matter-field interaction.

Under the declared Gaussian convention:
    I_face = 4*pi*sqrt(gamma)*(alpha*j_e^r - beta^r*rho_e)
    dQ_face/dt = -I_face
    E_hat_r(face) = Q_face / R_areal(face)**2

Here sqrt(gamma) excludes sin(theta) and includes the spherical r**2 factor.
The current density is Eulerian spatial current with its radial index raised,
j_e^r. All arrays must describe the same integrator stage.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


FOUR_PI = 4.0 * math.pi


def _finite_array(values, name: str, *, length: int | None = None) -> np.ndarray:
    """Copy a finite real 1-D array, rejecting booleans and complex inputs."""
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real 1-D array") from exc
    if raw.ndim != 1 or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a finite real 1-D array")
    result = np.array(raw, dtype=float, copy=True)
    if result.size < 2:
        raise ValueError(f"{name} must contain at least two radial faces")
    if length is not None and result.size != length:
        raise ValueError(f"{name} must contain exactly {length} face values")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _positive_dt(dt: float) -> float:
    if isinstance(dt, (bool, np.bool_)):
        raise ValueError("dt must be finite and positive")
    try:
        value = float(dt)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("dt must be finite and positive") from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("dt must be finite and positive")
    return value


@dataclass(frozen=True)
class RadialChargeState:
    """Cumulative charge Q on the radial faces, with Q(0)=0.

    enclosed_charge[i] is the charge enclosed by face i. Thus the charge in
    shell i is enclosed_charge[i+1] - enclosed_charge[i]. The array is
    defensively copied and made read-only. Charge units and current convention
    must be chosen by the caller.
    """

    enclosed_charge: np.ndarray

    def __post_init__(self) -> None:
        charge = _finite_array(self.enclosed_charge, "enclosed_charge")
        if charge[0] != 0.0:
            raise ValueError("enclosed charge at the regular centre must be zero")
        charge.setflags(write=False)
        object.__setattr__(self, "enclosed_charge", charge)

    @property
    def cell_charge(self) -> np.ndarray:
        """Return shell-integrated charge, without assuming a volume rule."""
        return np.diff(self.enclosed_charge)

    def advance(
        self,
        dt: float,
        integrated_outward_current_faces,
    ) -> "RadialChargeState":
        """Advance enclosed charge using one supplied stage-local face flux.

        The flux is the outward charge-current integrated over the sphere,
        already including 4*pi; positive flux reduces enclosed charge.
        Regularity requires exactly zero integrated flux through the centre.
        The outer-face flux may be nonzero and sets the global charge change,
        so charge is not artificially conserved across an open outer boundary.
        """
        step = _positive_dt(dt)
        flux = _finite_array(
            integrated_outward_current_faces,
            "integrated_outward_current_faces",
            length=self.enclosed_charge.size,
        )
        if flux[0] != 0.0:
            raise ValueError("integrated charge-current flux at the centre must be zero")
        with np.errstate(over="ignore", invalid="ignore"):
            updated = self.enclosed_charge - step * flux
        if not np.all(np.isfinite(updated)):
            raise ValueError("charge update produced a non-finite face value")
        # Both operands are exactly zero at the centre; keep that invariant
        # explicit rather than allowing future arithmetic refactors to drift.
        updated[0] = 0.0
        return RadialChargeState(updated)


def integrated_outward_current_flux_faces(
    alpha_faces,
    beta_faces,
    sqrt_gamma_faces,
    charge_density_faces,
    spatial_current_r_faces,
) -> np.ndarray:
    """Compute the Gaussian-convention charge flux integrated over each sphere.

    Inputs are face samples from one stage:
    - alpha: positive lapse;
    - beta: contravariant radial shift beta^r;
    - sqrt_gamma: sqrt(det(gamma_ij)) without angular sin(theta), including r^2;
    - charge_density: Eulerian charge density rho_e = -n_mu J^mu;
    - spatial_current_r: contravariant Eulerian spatial current j_e^r.

    Returns I = 4*pi*sqrt(gamma)*(alpha*j_e^r - beta^r*rho_e).
    At a regular spherical origin sqrt(gamma)=0, so I(0)=0.
    """
    alpha = _finite_array(alpha_faces, "alpha_faces")
    count = alpha.size
    beta = _finite_array(beta_faces, "beta_faces", length=count)
    sqrt_gamma = _finite_array(sqrt_gamma_faces, "sqrt_gamma_faces", length=count)
    rho_e = _finite_array(charge_density_faces, "charge_density_faces", length=count)
    current_r = _finite_array(
        spatial_current_r_faces, "spatial_current_r_faces", length=count
    )
    if np.any(alpha <= 0.0):
        raise ValueError("lapse must be strictly positive at every face")
    if np.any(sqrt_gamma < 0.0):
        raise ValueError("sqrt_gamma must be nonnegative")
    if sqrt_gamma[0] != 0.0:
        raise ValueError("spherical sqrt_gamma must vanish at the regular centre")
    with np.errstate(over="ignore", invalid="ignore"):
        flux = FOUR_PI * sqrt_gamma * (alpha * current_r - beta * rho_e)
    if not np.all(np.isfinite(flux)):
        raise ValueError("integrated charge-current flux is non-finite")
    flux[0] = 0.0
    return flux


def radial_electric_field_faces(
    state: RadialChargeState,
    areal_radius_faces,
) -> np.ndarray:
    """Recover E_hat_r=Q/R_areal^2 under the declared Gaussian convention.

    The centre value is set to zero by regularity. Noncentral areal radii must
    be positive; they are not required to be monotone, because this helper
    must not impose an extra geometric restriction on a strong-field slice.
    This is a Gauss-law reconstruction, not an independent wave evolution.
    """
    radius = _finite_array(
        areal_radius_faces,
        "areal_radius_faces",
        length=state.enclosed_charge.size,
    )
    if radius[0] != 0.0:
        raise ValueError("areal radius at the regular centre must be zero")
    if np.any(radius[1:] <= 0.0):
        raise ValueError("noncentral areal radii must be strictly positive")
    field = np.zeros_like(radius)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        field[1:] = state.enclosed_charge[1:] / (radius[1:] * radius[1:])
    if not np.all(np.isfinite(field)):
        raise ValueError("radial electric field reconstruction is non-finite")
    return field
