"""Conservative spherical matter transport for the V5.5 production adapter.

This module evolves the covariant mixed-tensor variables defined in the
Valencia layer. The discretization is conservative finite-volume transport
with piecewise-linear reconstruction and HLL numerical fluxes.

The physical source identity remains

    partial_mu(sqrt(-g) T^mu_nu)
      = 1/2 sqrt(-g) T^(ab) partial_nu g_ab + sqrt(-g) Q_nu.

For dark matter:
    Q_nu = beta rho_DM partial_nu phi.

For baryons and radiation:
    Q_nu = 0.

No interface reservoir, feedback closure, or fitted physical coefficient is
introduced.
"""
from dataclasses import dataclass
from enum import Enum
import math
import numpy as np

from .valencia import (
    FluidConserved,
    FluidPrimitive,
    SphericalMetric,
    dark_matter_covector_source,
    recover_dust,
    recover_radiation,
    spherical_metric_from_bssn,
    mixed_flux,
    primitive_to_conserved,
)

KAPPA = 8.0 * math.pi


class Species(str, Enum):
    DARK_MATTER = "dark_matter"
    BARYON = "baryon"
    RADIATION = "radiation"


@dataclass
class ConservedSpecies:
    rest: np.ndarray
    energy_t: np.ndarray
    momentum_r: np.ndarray

    def copy(self) -> "ConservedSpecies":
        return ConservedSpecies(
            self.rest.copy(), self.energy_t.copy(), self.momentum_r.copy()
        )


@dataclass(frozen=True)
class BSSNMetricSlice:
    r: np.ndarray
    a: np.ndarray
    b: np.ndarray
    X: np.ndarray
    alpha: np.ndarray
    beta: np.ndarray
    Aa: np.ndarray
    K: np.ndarray
    Lambda: np.ndarray
    B: np.ndarray


@dataclass(frozen=True)
class MetricDerivativeSet:
    time: dict[str, np.ndarray]
    radial: dict[str, np.ndarray]


def _metric_arrays(g: BSSNMetricSlice) -> list[SphericalMetric]:
    return [
        spherical_metric_from_bssn(
            float(r), float(a), float(b), float(X), float(alpha), float(beta)
        )
        for r, a, b, X, alpha, beta in zip(
            g.r, g.a, g.b, g.X, g.alpha, g.beta
        )
    ]


def geometry_metric_derivatives(
    g: BSSNMetricSlice,
    d1,
) -> MetricDerivativeSet:
    """Construct partial_t g_ab and partial_r g_ab from ADM kinematics."""
    grr = g.a / (g.X * g.X)
    gthth = g.b * g.r * g.r / (g.X * g.X)

    dr_a = d1(g.a, 1)
    dr_b = d1(g.b, 1)
    dr_X = d1(g.X, 1)
    dr_alpha = d1(g.alpha, 1)
    dr_beta = d1(g.beta, -1)

    dr_grr = grr * (dr_a / g.a - 2.0 * dr_X / g.X)
    dr_gthth = gthth * (
        dr_b / g.b + 2.0 / g.r - 2.0 * dr_X / g.X
    )
    dr_gtt = (
        -2.0 * g.alpha * dr_alpha
        + dr_grr * g.beta**2
        + 2.0 * grr * g.beta * dr_beta
    )
    dr_gtr = dr_grr * g.beta + grr * dr_beta

    Krr = grr * (g.K / 3.0 + g.Aa)
    Kthth = gthth * (g.K / 3.0 - g.Aa / 2.0)

    dt_grr = -2.0 * g.alpha * Krr + g.beta * dr_grr + 2.0 * grr * dr_beta
    dt_gthth = -2.0 * g.alpha * Kthth + g.beta * dr_gthth

    dt_alpha = -2.0 * g.alpha * g.K + g.beta * dr_alpha
    dt_beta = g.B
    dt_gtt = (
        -2.0 * g.alpha * dt_alpha
        + dt_grr * g.beta**2
        + 2.0 * grr * g.beta * dt_beta
    )
    dt_gtr = dt_grr * g.beta + grr * dt_beta

    return MetricDerivativeSet(
        time={"tt": dt_gtt, "tr": dt_gtr, "rr": dt_grr, "thth": dt_gthth},
        radial={"tt": dr_gtt, "tr": dr_gtr, "rr": dr_grr, "thth": dr_gthth,
                 "alpha": dr_alpha, "beta": dr_beta},
    )


def _state_at(U: ConservedSpecies, i: int) -> FluidConserved:
    return FluidConserved(
        rest=float(U.rest[i]),
        energy_t=float(U.energy_t[i]),
        momentum_r=float(U.momentum_r[i]),
    )


def primitives(
    metrics: list[SphericalMetric],
    U: ConservedSpecies,
    species: Species,
) -> list[FluidPrimitive]:
    out: list[FluidPrimitive] = []
    for i, metric in enumerate(metrics):
        state = _state_at(U, i)
        if species in (Species.DARK_MATTER, Species.BARYON):
            out.append(recover_dust(metric, state))
        else:
            try:
                out.append(recover_radiation(metric, state))
            except ValueError as exc:
                raise ValueError(
                    f"radiation inversion failure at cell i={i}: "
                    f"alpha={metric.alpha:.17e}, beta={metric.beta:.17e}, "
                    f"gamma_rr={metric.gamma_rr:.17e}; {exc}"
                ) from exc
    return out


def _minmod(a: float, b: float) -> float:
    if a * b <= 0.0:
        return 0.0
    return math.copysign(min(abs(a), abs(b)), a)


def _reconstruct(values: np.ndarray, i: int) -> tuple[float, float]:
    n = len(values)
    if i <= 0:
        return float(values[0]), float(values[0])
    if i >= n - 1:
        return float(values[-1]), float(values[-1])
    slope = _minmod(
        float(values[i] - values[i - 1]),
        float(values[i + 1] - values[i]),
    )
    return float(values[i] - 0.5 * slope), float(values[i] + 0.5 * slope)


def _reconstructed_primitive(
    prim: list[FluidPrimitive],
    i: int,
    side: str,
    gamma_rr: float,
) -> FluidPrimitive:
    if side not in ("left", "right"):
        raise ValueError("side must be left or right")
    attrs = {}
    for name in ("rho", "pressure", "v_r"):
        vals = np.asarray([getattr(q, name) for q in prim])
        lo, hi = _reconstruct(vals, i)
        attrs[name] = lo if side == "left" else hi
    rho = max(0.0, attrs["rho"])
    pressure = max(0.0, attrs["pressure"])
    v = attrs["v_r"]
    vmax = (1.0 - 1.0e-12) / math.sqrt(gamma_rr)
    v = max(-vmax, min(vmax, v))
    return FluidPrimitive(
        rho=rho, pressure=pressure, v_r=v, gamma_rr=gamma_rr
    )


def _cons_vector(U: FluidConserved) -> np.ndarray:
    return np.asarray([U.rest, U.energy_t, U.momentum_r], dtype=float)


def _valencia_conserved(metric: SphericalMetric, q: FluidPrimitive) -> np.ndarray:
    W = q.lorentz()
    h = q.rho + q.pressure
    D = q.rho * W
    E = h * W * W - q.pressure
    S_r = h * W * W * metric.gamma_rr * q.v_r
    return np.asarray([
        metric.sqrt_gamma * D,
        metric.sqrt_gamma * E,
        metric.sqrt_gamma * S_r,
    ], dtype=float)


def _valencia_flux(metric: SphericalMetric, q: FluidPrimitive) -> np.ndarray:
    W = q.lorentz()
    h = q.rho + q.pressure
    D = q.rho * W
    E = h * W * W - q.pressure
    S_r = h * W * W * metric.gamma_rr * q.v_r
    S_up = metric.gamma_rr_inv * S_r
    return np.asarray([
        metric.sqrt_gamma * D * (metric.alpha * q.v_r - metric.beta),
        metric.sqrt_gamma * (metric.alpha * S_up - metric.beta * E),
        metric.sqrt_gamma * (
            metric.alpha * S_r * q.v_r
            + metric.alpha * q.pressure
            - metric.beta * S_r
        ),
    ], dtype=float)


def _valencia_source(
    metric: SphericalMetric,
    q: FluidPrimitive,
    K: float,
    Aa: float,
    dr_alpha: float,
    dr_beta: float,
    dr_grr: float,
    dr_gthth: float,
) -> tuple[float, float]:
    W = q.lorentz()
    h = q.rho + q.pressure
    E = h * W * W - q.pressure
    S_r = h * W * W * metric.gamma_rr * q.v_r
    S_up = metric.gamma_rr_inv * S_r
    Srr = h * W * W * q.v_r * q.v_r + q.pressure * metric.gamma_rr_inv
    Sthth = q.pressure * metric.gamma_thth_inv
    Krr = metric.gamma_rr * (K / 3.0 + Aa)
    Kthth = metric.gamma_thth * (K / 3.0 - Aa / 2.0)
    energy = metric.sqrt_gamma * (
        metric.alpha * (Krr * Srr + 2.0 * Kthth * Sthth)
        - S_up * dr_alpha
    )
    momentum = metric.sqrt_gamma * (
        -E * dr_alpha
        + S_r * dr_beta
        + 0.5 * metric.alpha * (Srr * dr_grr + 2.0 * Sthth * dr_gthth)
    )
    return energy, momentum


def _hll_flux(
    metric: SphericalMetric,
    left: FluidPrimitive,
    right: FluidPrimitive,
    species: Species,
) -> np.ndarray:
    fl = _valencia_flux(metric, left)
    fr = _valencia_flux(metric, right)
    ul = _valencia_conserved(metric, left)
    ur = _valencia_conserved(metric, right)

    c_s = 0.0 if species != Species.RADIATION else 1.0 / math.sqrt(3.0)

    def speed(q: FluidPrimitive, sign: float) -> float:
        vhat = math.sqrt(metric.gamma_rr) * q.v_r
        denom = 1.0 + sign * vhat * c_s
        num = q.v_r + sign * c_s / math.sqrt(metric.gamma_rr)
        return metric.alpha * num / denom - metric.beta

    s_minus = min(0.0, speed(left, -1.0), speed(right, -1.0))
    s_plus = max(0.0, speed(left, 1.0), speed(right, 1.0))
    if s_plus <= 0.0:
        return fr
    if s_minus >= 0.0:
        return fl
    return (
        s_plus * fl
        - s_minus * fr
        + s_minus * s_plus * (ur - ul)
    ) / (s_plus - s_minus)


def species_projection(
    metrics: list[SphericalMetric],
    prim: list[FluidPrimitive],
) -> dict[str, np.ndarray]:
    """Project a perfect fluid onto the Eulerian stress-energy components."""
    n = len(metrics)
    rho = np.empty(n)
    pr = np.empty(n)
    pt = np.empty(n)
    j = np.empty(n)
    for i, (m, q) in enumerate(zip(metrics, prim)):
        W = q.lorentz()
        h = q.rho + q.pressure
        vhat2 = m.gamma_rr * q.v_r * q.v_r
        rho[i] = (h * W * W - q.pressure) / KAPPA
        pr[i] = (h * W * W * vhat2 + q.pressure) / KAPPA
        pt[i] = q.pressure / KAPPA
        j[i] = h * W * W * m.gamma_rr * q.v_r / KAPPA
    return {"rho": rho, "pr": pr, "pt": pt, "j": j}


def initialize_dust(
    metrics: list[SphericalMetric],
    rho: np.ndarray,
    v_r: np.ndarray | None = None,
) -> ConservedSpecies:
    v = np.zeros_like(rho) if v_r is None else np.asarray(v_r, dtype=float)
    if rho.shape != v.shape or rho.shape != (len(metrics),):
        raise ValueError("rho and v_r must match the metric grid")
    rest = np.empty_like(rho)
    energy = np.empty_like(rho)
    momentum = np.empty_like(rho)
    for i, m in enumerate(metrics):
        q = FluidPrimitive(
            rho=float(rho[i]), pressure=0.0,
            v_r=float(v[i]), gamma_rr=m.gamma_rr,
        )
        W = q.lorentz()
        E = q.rho * W * W
        D = q.rho * W
        S_r = q.rho * W * W * m.gamma_rr * q.v_r
        rest[i] = m.sqrt_gamma * D
        energy[i] = m.sqrt_gamma * E
        momentum[i] = m.sqrt_gamma * S_r
    return ConservedSpecies(rest, energy, momentum)


def initialize_radiation(
    metrics: list[SphericalMetric],
    rho: np.ndarray,
    v_r: np.ndarray | None = None,
) -> ConservedSpecies:
    v = np.zeros_like(rho) if v_r is None else np.asarray(v_r, dtype=float)
    if rho.shape != v.shape or rho.shape != (len(metrics),):
        raise ValueError("rho and v_r must match the metric grid")
    rest = np.zeros_like(rho)
    energy = np.empty_like(rho)
    momentum = np.empty_like(rho)
    for i, m in enumerate(metrics):
        q = FluidPrimitive(
            rho=float(rho[i]),
            pressure=float(rho[i]) / 3.0,
            v_r=float(v[i]),
            gamma_rr=m.gamma_rr,
        )
        W = q.lorentz()
        h = q.rho + q.pressure
        E = h * W * W - q.pressure
        S_r = h * W * W * m.gamma_rr * q.v_r
        rest[i] = 0.0
        energy[i] = m.sqrt_gamma * E
        momentum[i] = m.sqrt_gamma * S_r
    return ConservedSpecies(rest, energy, momentum)


def _radiation_conserved(metric: SphericalMetric, q: FluidPrimitive) -> np.ndarray:
    W = q.lorentz()
    h = q.rho + q.pressure
    E = h * W * W - q.pressure
    S_r = h * W * W * metric.gamma_rr * q.v_r
    return np.asarray([0.0, metric.sqrt_gamma * E, metric.sqrt_gamma * S_r], dtype=float)


def _radiation_flux(metric: SphericalMetric, q: FluidPrimitive) -> np.ndarray:
    W = q.lorentz()
    h = q.rho + q.pressure
    E = h * W * W - q.pressure
    S_r = h * W * W * metric.gamma_rr * q.v_r
    S_up = metric.gamma_rr_inv * S_r
    return np.asarray([
        0.0,
        metric.sqrt_gamma * (metric.alpha * S_up - metric.beta * E),
        metric.sqrt_gamma * (metric.alpha * S_r * q.v_r
                              + metric.alpha * q.pressure
                              - metric.beta * S_r),
    ], dtype=float)


def _radiation_source(metric, q, K, Aa, dr_alpha, dr_beta, dr_grr, dr_gthth):
    W = q.lorentz()
    h = q.rho + q.pressure
    E = h * W * W - q.pressure
    S_r = h * W * W * metric.gamma_rr * q.v_r
    S_up = metric.gamma_rr_inv * S_r
    Srr = h * W * W * q.v_r * q.v_r + q.pressure * metric.gamma_rr_inv
    Sthth = q.pressure * metric.gamma_thth_inv
    Krr = metric.gamma_rr * (K / 3.0 + Aa)
    Kthth = metric.gamma_thth * (K / 3.0 - Aa / 2.0)
    energy = metric.sqrt_gamma * (
        metric.alpha * (Krr * Srr + 2.0 * Kthth * Sthth) - S_up * dr_alpha
    )
    momentum = metric.sqrt_gamma * (
        -E * dr_alpha + S_r * dr_beta
        + 0.5 * metric.alpha * (Srr * dr_grr + 2.0 * Sthth * dr_gthth)
    )
    return energy, momentum


def evolve_species(
    metric: BSSNMetricSlice,
    metric_derivatives: MetricDerivativeSet,
    state: ConservedSpecies,
    species: Species,
    dt: float,
    dphi_t: np.ndarray | None = None,
    dphi_r: np.ndarray | None = None,
    beta_dm: float = -0.04,
) -> ConservedSpecies:
    """Advance one conservative finite-volume step."""
    n = len(metric.r)
    if any(np.shape(x) != (n,) for x in (
        state.rest, state.energy_t, state.momentum_r
    )):
        raise ValueError("conserved state shape mismatch")

    metrics = _metric_arrays(metric)
    prim = primitives(metrics, state, species)

    face_flux = np.zeros((n + 1, 3), dtype=float)

    # Regular center: the integrated spherical flux through r=0 vanishes.
    face_flux[0] = 0.0
    for i in range(n - 1):
        m0, m1 = metrics[i], metrics[i + 1]
        mf = SphericalMetric(
            alpha=0.5 * (m0.alpha + m1.alpha),
            beta=0.5 * (m0.beta + m1.beta),
            gamma_rr=0.5 * (m0.gamma_rr + m1.gamma_rr),
            gamma_rr_inv=0.5 * (m0.gamma_rr_inv + m1.gamma_rr_inv),
            gamma_thth=0.5 * (m0.gamma_thth + m1.gamma_thth),
            gamma_thth_inv=0.5 * (m0.gamma_thth_inv + m1.gamma_thth_inv),
            sqrt_gamma=0.5 * (m0.sqrt_gamma + m1.sqrt_gamma),
        )
        ql = _reconstructed_primitive(prim, i, "right", mf.gamma_rr)
        qr = _reconstructed_primitive(prim, i + 1, "left", mf.gamma_rr)
        face_flux[i + 1] = _hll_flux(mf, ql, qr, species)

    # Causal/outflow outer closure: continue the last physical state.
    face_flux[-1] = _valencia_flux(metrics[-1], prim[-1])

    out = state.copy()
    inv_dr = 1.0 / (metric.r[1] - metric.r[0])
    for i, m in enumerate(metrics):
        source_e, source_s = _valencia_source(
            m, prim[i], float(metric.K[i]), float(metric.Aa[i]),
            float(metric_derivatives.radial["alpha"][i]),
            float(metric_derivatives.radial["beta"][i]),
            float(metric_derivatives.radial["rr"][i]),
            float(metric_derivatives.radial["thth"][i]),
        )

        q_e = q_s = 0.0
        if species == Species.DARK_MATTER:
            if dphi_t is None or dphi_r is None:
                raise ValueError("DM evolution requires dphi_t and dphi_r")
            q_t, q_r = dark_matter_covector_source(
                beta_dm, prim[i].rho, float(dphi_t[i]), float(dphi_r[i])
            )
            q_e = -m.sqrt_gamma * (q_t - m.beta * q_r)
            q_s = m.sqrt_gamma * q_r

        out.rest[i] -= dt * inv_dr * (
            face_flux[i + 1, 0] - face_flux[i, 0]
        )
        out.energy_t[i] -= dt * inv_dr * (
            face_flux[i + 1, 1] - face_flux[i, 1]
        )
        out.momentum_r[i] -= dt * inv_dr * (
            face_flux[i + 1, 2] - face_flux[i, 2]
        )
        out.energy_t[i] += dt * (source_e + q_e)
        out.momentum_r[i] += dt * (source_s + q_s)

    return out


def project_species(
    metric: BSSNMetricSlice,
    state: ConservedSpecies,
    species: Species,
) -> dict[str, np.ndarray]:
    """Recover primitives and return Einstein-source projections."""
    metrics = _metric_arrays(metric)
    return species_projection(metrics, primitives(metrics, state, species))
