"""Covariant spherical fluid state and flux layer for V5.5.

The conservative variables are the densitized rest-mass current and the
mixed stress-energy components

    U_D = sqrt(-g) J^t
    U_t = sqrt(-g) T^t_t
    U_r = sqrt(-g) T^t_r.

Their radial fluxes and metric source are direct consequences of

    ∂_μ(sqrt(-g) T^μ_ν)
      = 1/2 sqrt(-g) T^{ab} ∂_ν g_ab + sqrt(-g) Q_ν.

For the archived dark-sector exchange,

    Q_DM,ν = + beta rho_DM ∂_ν phi,

while baryons and radiation have Q_ν = 0.

No Riemann solver, reconstruction, numerical dissipation, or artificial
interface source is selected here; those remain production-integrator
choices.
"""
from dataclasses import dataclass
import math
import numpy as np

KAPPA = 8.0 * math.pi


@dataclass(frozen=True)
class SphericalMetric:
    alpha: float
    beta: float
    gamma_rr: float
    gamma_rr_inv: float
    sqrt_gamma: float

    @property
    def shift_cov(self) -> float:
        return self.gamma_rr * self.beta

    @property
    def sqrt_minus_g(self) -> float:
        return self.alpha * self.sqrt_gamma


@dataclass(frozen=True)
class FluidPrimitive:
    rho: float
    pressure: float
    v_r: float

    @property
    def v2(self) -> float:
        return self.gamma_inv * self.v_r * self.v_r

    gamma_inv: float = 1.0

    def lorentz(self) -> float:
        v2 = self.v2
        if v2 < 0.0 or v2 >= 1.0:
            raise ValueError("|v| must be < 1 in the Eulerian frame")
        return 1.0 / math.sqrt(1.0 - v2)


@dataclass(frozen=True)
class FluidConserved:
    rest: float
    energy_t: float
    momentum_r: float


def primitive_to_conserved(metric: SphericalMetric, q: FluidPrimitive) -> FluidConserved:
    W = q.lorentz()
    h = q.rho + q.pressure
    S_r = h * W * W * metric.gamma_rr * q.v_r
    E = h * W * W - q.pressure
    Tt_t = -metric.alpha * E + metric.beta * S_r
    Tt_r = S_r
    Jt = q.rho * W / metric.alpha
    rootg = metric.sqrt_minus_g
    return FluidConserved(
        rest=rootg * Jt,
        energy_t=rootg * Tt_t,
        momentum_r=rootg * Tt_r,
    )


def recover_dust(metric: SphericalMetric, U: FluidConserved) -> FluidPrimitive:
    if U.rest <= 0.0:
        raise ValueError("dust rest-current must be positive")
    D = U.rest / metric.sqrt_gamma
    S_r = U.momentum_r / metric.alpha / metric.sqrt_gamma
    if D <= 0.0:
        raise ValueError("dust conserved density must be positive")
    s2 = metric.gamma_rr_inv * S_r * S_r
    v2 = s2 / (D * D + s2)
    v_r = math.copysign(
        math.sqrt(v2 / metric.gamma_rr_inv) if v2 else 0.0, S_r
    )
    W = 1.0 / math.sqrt(1.0 - v2)
    rho = D / W
    return FluidPrimitive(rho=rho, pressure=0.0, v_r=v_r,
                          gamma_inv=metric.gamma_rr_inv)


def conserved_energy(metric: SphericalMetric, U: FluidConserved) -> tuple[float, float]:
    D = U.rest / metric.sqrt_gamma
    S_r = U.momentum_r / (metric.alpha * metric.sqrt_gamma)
    E = (metric.alpha * metric.beta * S_r
         - U.energy_t / metric.sqrt_gamma) / (metric.alpha * metric.alpha)
    return D, E


def recover_radiation(metric: SphericalMetric, U: FluidConserved,
                       tol: float = 1.0e-12,
                       max_iter: int = 80) -> FluidPrimitive:
    """Recover the p=rho/3 radiation primitive from conservative variables."""
    D, E = conserved_energy(metric, U)
    if D < 0.0 or E < 0.0:
        raise ValueError("radiation conservative state is nonphysical")
    S_r = U.momentum_r / (metric.alpha * metric.sqrt_gamma)
    S2 = metric.gamma_rr_inv * S_r * S_r
    if S2 == 0.0:
        rho = E
        return FluidPrimitive(rho=rho, pressure=rho / 3.0, v_r=0.0,
                              gamma_inv=metric.gamma_rr_inv)

    def residual(p: float) -> tuple[float, float]:
        # Radiation: h = 4 rho / 3 and E = h W^2 - p.
        rho = 3.0 * p
        denom = E + p
        v2 = S2 / (denom * denom)
        if v2 >= 1.0:
            return math.inf, v2
        # h W^2 = E+p, so rho = 3p and the identity must match E.
        W2 = 1.0 / (1.0 - v2)
        implied_E = 4.0 * rho * W2 / 3.0 - p
        return implied_E - E, v2

    lo = max(0.0, 0.5 * (E - math.sqrt(E * E + S2)))
    hi = max(E, math.sqrt(E * E + S2)) + 1.0
    # A safeguarded bisection is used because the EOS inversion is part of the
    # conservative state definition, not an additional physical closure.
    f_lo, _ = residual(lo)
    f_hi, _ = residual(hi)
    if not math.isfinite(f_lo) or not math.isfinite(f_hi) or f_lo * f_hi > 0.0:
        lo = 1.0e-16
        hi = max(E + math.sqrt(S2), 1.0e-12)
        for _ in range(20):
            f_hi, _ = residual(hi)
            if math.isfinite(f_hi) and f_lo * f_hi <= 0.0:
                break
            hi *= 2.0
        else:
            raise ValueError("radiation primitive inversion has no bracket")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid, v2 = residual(mid)
        if abs(f_mid) <= tol * max(E, 1.0):
            rho = 3.0 * mid
            v = math.sqrt(v2 / metric.gamma_rr_inv)
            v = math.copysign(v, S_r)
            return FluidPrimitive(rho=rho, pressure=rho / 3.0,
                                  v_r=v, gamma_inv=metric.gamma_rr_inv)
        if f_lo * f_mid <= 0.0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    mid = 0.5 * (lo + hi)
    f_mid, v2 = residual(mid)
    if not math.isfinite(f_mid) or v2 >= 1.0:
        raise ValueError("radiation primitive inversion did not converge")
    rho = 3.0 * mid
    v = math.sqrt(v2 / metric.gamma_rr_inv)
    v = math.copysign(v, S_r)
    return FluidPrimitive(rho=rho, pressure=rho / 3.0,
                          v_r=v, gamma_inv=metric.gamma_rr_inv)


def four_velocity(metric: SphericalMetric, q: FluidPrimitive) -> tuple[float, float]:
    W = q.lorentz()
    u_t = W * (-metric.alpha + metric.shift_cov * q.v_r)
    u_r = W * q.v_r
    return u_t, u_r


def stress_tensor_contravariant(
    metric: SphericalMetric, q: FluidPrimitive
) -> dict[str, float]:
    W = q.lorentz()
    h = q.rho + q.pressure
    v = q.v_r
    u_t, _ = four_velocity(metric, q)
    u0 = W / metric.alpha
    ur = W * (v - metric.beta / metric.alpha)
    gtt = -1.0 / (metric.alpha * metric.alpha)
    gtr = metric.beta / (metric.alpha * metric.alpha)
    grr = metric.gamma_rr_inv - (metric.beta / metric.alpha) ** 2
    Ttt = h * u0 * u0 + q.pressure * gtt
    Ttr = h * u0 * ur + q.pressure * gtr
    Trr = h * ur * ur + q.pressure * grr
    return {"tt": Ttt, "tr": Ttr, "rr": Trr}


def mixed_flux(metric: SphericalMetric, q: FluidPrimitive) -> FluidConserved:
    """Return sqrt(-g) times the radial fluxes (J^r, T^r_t, T^r_r)."""
    W = q.lorentz()
    h = q.rho + q.pressure
    u0 = W / metric.alpha
    ur = W * (q.v_r - metric.beta / metric.alpha)
    ut, _ = four_velocity(metric, q)

    Jr = q.rho * ur
    Tr_t = h * ur * ut
    Tr_r = h * ur * W * metric.gamma_rr * q.v_r + q.pressure
    rootg = metric.sqrt_minus_g
    return FluidConserved(
        rest=rootg * Jr,
        energy_t=rootg * Tr_t,
        momentum_r=rootg * Tr_r,
    )


def metric_mixed_components(metric: SphericalMetric) -> dict[str, float]:
    gtt = -metric.alpha * metric.alpha + metric.gamma_rr * metric.beta**2
    gtr = metric.gamma_rr * metric.beta
    grr = metric.gamma_rr
    return {"tt": gtt, "tr": gtr, "rr": grr}


def geometric_source(
    metric: SphericalMetric,
    q: FluidPrimitive,
    dg_tt: float,
    dg_tr: float,
    dg_rr: float,
) -> tuple[float, float]:
    """Return (source_t, source_r) for the mixed stress equations."""
    T = stress_tensor_contravariant(metric, q)
    pref = 0.5 * metric.sqrt_minus_g
    source_t = pref * T["tt"] * dg_tt + pref * 2.0 * T["tr"] * dg_tr
    source_r = pref * T["tt"] * 0.0 + pref * T["tr"] * 0.0 + pref * T["rr"] * dg_rr
    return source_t, source_r


def dark_matter_covector_source(beta_dm: float, rho_dm: float,
                                dphi_t: float, dphi_r: float) -> tuple[float, float]:
    """Locked V5.5 DM exchange covector Q_nu."""
    return beta_dm * rho_dm * dphi_t, beta_dm * rho_dm * dphi_r


def spherical_metric_from_bssn(
    r: float, a: float, b: float, X: float, alpha: float, beta: float
) -> SphericalMetric:
    gamma_rr = a / (X * X)
    gamma_rr_inv = (X * X) / a
    # sqrt(det gamma) = sqrt(a) b r^2 / X^3 times sin(theta).
    # The radial cell measure is supplied separately by the cell-volume grid;
    # this scalar is the angular-free radial metric factor.
    sqrt_gamma = math.sqrt(max(a, 0.0)) * b * r * r / (X**3)
    return SphericalMetric(alpha, beta, gamma_rr, gamma_rr_inv, sqrt_gamma)
