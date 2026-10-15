"""Research-only algebraic diagnostics for selected COSMOS/GEAR identities.

These helpers intentionally do not import or modify the production engine.
They evaluate declared formulae only; they do not certify a physical model,
solve its equations of motion, or decide that a curvature wall is active.
"""
from __future__ import annotations

import math
from collections.abc import Sequence

Triple = tuple[float, float, float]


def _finite(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _triple(name: str, values: Sequence[float]) -> Triple:
    if len(values) != 3:
        raise ValueError(f"{name} must contain exactly three values")
    return tuple(_finite(f"{name}[{i}]", values[i]) for i in range(3))  # type: ignore[return-value]


def _positive(name: str, value: float) -> float:
    value = _finite(name, value)
    if value <= 0.0:
        raise ValueError(f"{name} must be strictly positive")
    return value


def _nonzero(name: str, value: float) -> float:
    value = _finite(name, value)
    if value == 0.0:
        raise ValueError(f"{name} must be nonzero; the normalized identity is undefined here")
    return value


def _require_traceless(sigma_i: Triple, *, rtol: float) -> None:
    rtol = _finite("rtol", rtol)
    if rtol < 0.0:
        raise ValueError("rtol must be nonnegative")
    scale = max(abs(x) for x in sigma_i)
    if scale > 0.0 and abs(sum(sigma_i)) > rtol * scale:
        raise ValueError("directional shear eigenvalues must sum to zero")


def omega_sigma_from_directional_shear(
    sigma_i: Sequence[float], H: float, *, trace_rtol: float = 1e-10
) -> float:
    """Return Omega_sigma = sum(sigma_i**2)/(6 H**2), requiring H != 0.

    This is the archive convention sum(sigma_i**2) = 6 Omega_sigma H**2.
    It is a normalized observable and must not be sampled at a turnaround H=0.
    A result outside [0, 1) is returned raw rather than clamped; it can signal
    that the assumed positive non-shear-density/Friedmann constraint is absent.
    """
    shear = _triple("sigma_i", sigma_i)
    _require_traceless(shear, rtol=trace_rtol)
    hubble = _nonzero("H", H)
    return sum(component * component for component in shear) / (6.0 * hubble * hubble)


def directional_kasner_exponents(
    sigma_i: Sequence[float], H: float, *, trace_rtol: float = 1e-10
) -> Triple:
    """Return q_i = 1/3 + sigma_i/(3H), with a traceless shear triple and H != 0.

    The dangerous extremal orientation changes when the sign of H changes.
    This function deliberately makes no orientation label or bounce judgement.
    """
    shear = _triple("sigma_i", sigma_i)
    _require_traceless(shear, rtol=trace_rtol)
    hubble = _nonzero("H", H)
    return tuple(1.0 / 3.0 + component / (3.0 * hubble) for component in shear)  # type: ignore[return-value]


def preserve_achieve_slope(w: float, omega_sigma: float) -> float:
    """Return d Omega_sigma / d ln(a) under the unsourced Bianchi-I assumptions.

    Formula: 3 (w - 1) Omega_sigma (1 - Omega_sigma). Requires the equation
    of state w=p/rho for a separately conserved, positive non-shear sector,
    isotropic pressure/no shear source, and 0 <= Omega_sigma < 1. It is invalid
    at H=0 because ln(a) is not a good time parameter there.
    """
    w = _finite("w", w)
    omega = _finite("omega_sigma", omega_sigma)
    if not (0.0 <= omega < 1.0):
        raise ValueError("omega_sigma must satisfy 0 <= omega_sigma < 1")
    return 3.0 * (w - 1.0) * omega * (1.0 - omega)


def raychaudhuri_bounce_margin(rho_plus_p: float, sigma2: float) -> float:
    """Return -(rho+p)-2 sigma^2 for Hdot = -(rho+p)/2 - sigma^2.

    A strictly positive result is necessary and sufficient for Hdot>0 (a
    local classical minimum of a under this stated convention). Zero is the
    marginal case; negative means the local bounce inequality is not met.
    This does not cover effective quantum corrections or other conventions.
    """
    rho_plus_p = _finite("rho_plus_p", rho_plus_p)
    sigma2 = _finite("sigma2", sigma2)
    if sigma2 < 0.0:
        raise ValueError("sigma2 must be nonnegative")
    return -rho_plus_p - 2.0 * sigma2


def heron_discriminant(a: float, b: float, c: float) -> float:
    """Return the Heron discriminant T for positive Bianchi-IX scale factors.

    T = 2 a^2 b^2 + 2 b^2 c^2 + 2 c^2 a^2 - a^4 - b^4 - c^4.
    Equivalently T=(a+b+c)(-a+b+c)(a-b+c)(a+b-c). Its sign classifies the
    triangle inequality only; T<0 is NOT by itself a wall-activity detector.
    """
    a = _positive("a", a)
    b = _positive("b", b)
    c = _positive("c", c)
    return (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c)


def bianchi_ix_trace_term(
    volume: float, rho: float, pressures: Sequence[float]
) -> float:
    """Return V**2 * (3 rho - sum_i P_i), the matter term in the GEAR-47 S identity."""
    volume = _positive("volume", volume)
    rho = _finite("rho", rho)
    p_i = _triple("pressures", pressures)
    return volume * volume * (3.0 * rho - sum(p_i))


def bianchi_ix_S_tau_derivative(
    scale_factors: Sequence[float], rho: float, pressures: Sequence[float]
) -> float:
    """Return S'=1/2[-T+V**2(3rho-sum P_i)] for stated diagonal Bianchi IX equations.

    Scale factors must be positive. The result is the derivative in the archive's
    tau-time gauge dt=V d tau. It is not a COSMOS production observable.
    """
    a, b, c = _triple("scale_factors", scale_factors)
    a = _positive("a", a)
    b = _positive("b", b)
    c = _positive("c", c)
    volume = a * b * c
    T = heron_discriminant(a, b, c)
    trace_term = bianchi_ix_trace_term(volume, rho, pressures)
    return 0.5 * (-T + trace_term)


def scalar_kasner_exponent_square_sum(phi_prime: float, S: float) -> float:
    """Return sum(p_i**2)=1-phi_prime**2/S**2 in the declared scalar Bianchi-IX regime.

    Formula-only helper: caller must independently establish that the scalar
    reduced system and its conserved phi_prime assumptions apply.
    """
    phi_prime = _finite("phi_prime", phi_prime)
    S = _positive("S", S)
    return 1.0 - (phi_prime / S) ** 2
