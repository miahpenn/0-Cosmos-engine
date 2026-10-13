"""Analytic flat-space oracle for the research-only spherical Maxwell/Coulomb gate.

This module intentionally uses only the Python standard library and is not
connected to the production kernel. It defines a manufactured compact charge
profile and exact reference values for testing a future radial Maxwell solver.
It is not a field-evolution algorithm and makes no claim that GEAR generates
charge or electromagnetism.

Convention: Gaussian units with c=1, so
    (1/r^2) d(r^2 E_r)/dr = 4*pi*rho_e,
    E_r(r) = Q_enclosed(r)/r^2,
    Q_enclosed(r) = 4*pi integral_0^r rho_e(s) s^2 ds.
"""

from __future__ import annotations

from math import isfinite, pi
from typing import NamedTuple


class RadialElectricStress(NamedTuple):
    """Eulerian orthonormal projections for a pure radial electric field."""

    energy_density: float
    radial_pressure: float
    tangential_pressure: float
    momentum_density: float


def _profile_inputs(radius: float, rho0: float, support_radius: float) -> tuple[float, float, float]:
    """Validate scalar inputs and return normalized floats."""
    try:
        r = float(radius)
        rho = float(rho0)
        r0 = float(support_radius)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("radius, rho0, and support_radius must be finite real numbers") from exc
    if not (isfinite(r) and isfinite(rho) and isfinite(r0)):
        raise ValueError("radius, rho0, and support_radius must be finite")
    if r < 0.0:
        raise ValueError("radius must be non-negative")
    if r0 <= 0.0:
        raise ValueError("support_radius must be positive")
    return r, rho, r0


def charge_density(radius: float, rho0: float, support_radius: float) -> float:
    r"""Return rho_e(r)=rho0(1-r^2/R0^2)^2 inside R0, zero outside."""
    r, rho, r0 = _profile_inputs(radius, rho0, support_radius)
    if r >= r0:
        return 0.0
    x = r / r0
    return rho * (1.0 - x * x) ** 2


def total_charge(rho0: float, support_radius: float) -> float:
    r"""Exact total charge, Q_total = 32*pi*rho0*R0^3/105."""
    _, rho, r0 = _profile_inputs(0.0, rho0, support_radius)
    return (32.0 * pi / 105.0) * rho * r0**3


def enclosed_charge(radius: float, rho0: float, support_radius: float) -> float:
    r"""Exact enclosed charge Q(r)=4*pi*integral_0^r rho_e(s)s^2 ds."""
    r, rho, r0 = _profile_inputs(radius, rho0, support_radius)
    if r >= r0:
        return total_charge(rho, r0)
    x = r / r0
    # Factor r^3 to avoid subtracting high powers of a small dimensional r.
    shape = (1.0 / 3.0) - (2.0 * x * x / 5.0) + (x**4 / 7.0)
    return 4.0 * pi * rho * r**3 * shape


def radial_electric_field(radius: float, rho0: float, support_radius: float) -> float:
    """Exact physical E_hat_r for the manufactured charge profile."""
    r, rho, r0 = _profile_inputs(radius, rho0, support_radius)
    if r == 0.0:
        return 0.0
    if r >= r0:
        return total_charge(rho, r0) / (r * r)
    x = r / r0
    # Algebraically Q(r)/r^2, evaluated in a center-regular form.
    shape = (1.0 / 3.0) - (2.0 * x * x / 5.0) + (x**4 / 7.0)
    return 4.0 * pi * rho * r * shape


def radial_electric_stress(
    radius: float, rho0: float, support_radius: float
) -> RadialElectricStress:
    """Return pure-radial-E stress projections in an orthonormal frame.

    For zero magnetic field/Poynting flux, rho_EM=E^2/(8*pi), p_r=-rho_EM,
    p_t=+rho_EM, and radial momentum density is zero.
    """
    field = radial_electric_field(radius, rho0, support_radius)
    energy = field * field / (8.0 * pi)
    return RadialElectricStress(
        energy_density=energy,
        radial_pressure=-energy,
        tangential_pressure=energy,
        momentum_density=0.0,
    )


def shell_gauss_residual(
    inner_radius: float,
    outer_radius: float,
    inner_field: float,
    outer_field: float,
    rho0: float,
    support_radius: float,
) -> float:
    r"""Return the integrated flat-space shell Gauss residual.

    The signed residual is
        [r^2 E_r]_(inner)^outer - [Q(outer)-Q(inner)].
    It should be zero for the exact field, because d(r^2 E_r)/dr=4*pi*rho_e*r^2.
    Supplying fields independently makes this diagnostic useful as a negative
    control; this function does not silently replace them by the analytic E.
    """
    ri, rho, r0 = _profile_inputs(inner_radius, rho0, support_radius)
    ro, _, _ = _profile_inputs(outer_radius, rho, r0)
    try:
        ei = float(inner_field)
        eo = float(outer_field)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("field values must be finite real numbers") from exc
    if not (isfinite(ei) and isfinite(eo)):
        raise ValueError("field values must be finite")
    if ro <= ri:
        raise ValueError("outer_radius must be greater than inner_radius")
    flux_difference = ro * ro * eo - ri * ri * ei
    charge_difference = enclosed_charge(ro, rho, r0) - enclosed_charge(ri, rho, r0)
    return flux_difference - charge_difference
