"""Conditional GEAR-63 non-minimal-coupling algebra oracle.

Research-only helper. This module encodes one explicit action and is not
imported by the production engine. Here kappa means 8*pi*G, and mass_sq is
the positive parameter alpha**2 in V(d) = -mass_sq*d**2/2.
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class ExplicitXiState:
    """Algebraic witnesses for a nonsingular homogeneous state."""

    ricci_scalar: float
    curvature_force: float
    net_scalar_force: float
    effective_gravity_prefactor: float
    trace_denominator: float


def _finite_inputs(**values: float) -> None:
    for name, value in values.items():
        if not math.isfinite(value):
            raise ValueError(f"{name} must be finite")


def equation_denominators(kappa: float, xi: float, d: float) -> tuple[float, float]:
    """Return (effective-gravity prefactor, trace-elimination denominator)."""
    _finite_inputs(kappa=kappa, xi=xi, d=d)
    if kappa <= 0.0:
        raise ValueError("kappa must be positive")
    d2 = d * d
    effective = 1.0 - kappa * xi * d2
    trace = 1.0 - kappa * xi * (1.0 - 6.0 * xi) * d2
    _finite_inputs(effective_gravity_prefactor=effective, trace_denominator=trace)
    return effective, trace


def evaluate_state(
    kappa: float,
    xi: float,
    mass_sq: float,
    d: float,
    velocity: float,
    *,
    pole_tolerance: float = 1.0e-12,
) -> ExplicitXiState:
    """Evaluate exact R and forces, rejecting either singular prefactor.

    Assumptions: canonical scalar, signature (-+++), standard Ricci convention,
    V(d) = -mass_sq*d**2/2, no other matter and no constant-potential offset.
    """
    _finite_inputs(
        kappa=kappa, xi=xi, mass_sq=mass_sq, d=d, velocity=velocity,
        pole_tolerance=pole_tolerance,
    )
    if mass_sq < 0.0:
        raise ValueError("mass_sq must be nonnegative")
    if pole_tolerance <= 0.0:
        raise ValueError("pole_tolerance must be positive")

    effective, trace = equation_denominators(kappa, xi, d)
    if abs(effective) <= pole_tolerance:
        raise ValueError("effective-gravity prefactor is singular")
    if abs(trace) <= pole_tolerance:
        raise ValueError("trace-elimination denominator is singular")

    d2 = d * d
    v2 = velocity * velocity
    numerator = (1.0 - 6.0 * xi) * v2 + (2.0 - 6.0 * xi) * mass_sq * d2
    ricci = -kappa * numerator / trace
    curvature_force = -xi * ricci * d
    net_force = mass_sq * d + curvature_force
    _finite_inputs(ricci_scalar=ricci, curvature_force=curvature_force, net_scalar_force=net_force)
    return ExplicitXiState(
        ricci_scalar=ricci,
        curvature_force=curvature_force,
        net_scalar_force=net_force,
        effective_gravity_prefactor=effective,
        trace_denominator=trace,
    )


def closed_form_net_scalar_force(
    kappa: float, xi: float, mass_sq: float, d: float, velocity: float
) -> float:
    """Independently written eliminated form of mass_sq*d - xi*R*d."""
    state = evaluate_state(kappa, xi, mass_sq, d, velocity)
    d2 = d * d
    v2 = velocity * velocity
    return (
        mass_sq * d * (1.0 + kappa * xi * d2)
        + kappa * xi * (1.0 - 6.0 * xi) * v2 * d
    ) / state.trace_denominator


def trace_equation_residual(
    kappa: float, xi: float, mass_sq: float, d: float, velocity: float
) -> float:
    """Residual of the original trace equation after using the scalar equation.

    For this potential, T = v**2 + 2*mass_sq*d**2 and
    box(d**2) = -2*mass_sq*d**2 + 2*xi*R*d**2 - 2*v**2.
    The residual is (1-kappa*xi*d**2)*R + kappa*T
    + 3*kappa*xi*box(d**2).
    """
    state = evaluate_state(kappa, xi, mass_sq, d, velocity)
    d2 = d * d
    v2 = velocity * velocity
    trace_T = v2 + 2.0 * mass_sq * d2
    box_d2 = -2.0 * mass_sq * d2 + 2.0 * xi * state.ricci_scalar * d2 - 2.0 * v2
    return (
        state.effective_gravity_prefactor * state.ricci_scalar
        + kappa * trace_T
        + 3.0 * kappa * xi * box_d2
    )


def scalar_equation_residual(
    acceleration: float,
    hubble: float,
    velocity: float,
    d: float,
    xi: float,
    mass_sq: float,
    ricci_scalar: float,
) -> float:
    """Homogeneous scalar-EOM residual for the stated potential."""
    _finite_inputs(
        acceleration=acceleration, hubble=hubble, velocity=velocity, d=d,
        xi=xi, mass_sq=mass_sq, ricci_scalar=ricci_scalar,
    )
    return acceleration + 3.0 * hubble * velocity - mass_sq * d + xi * ricci_scalar * d


def friedmann_constraint_residual(
    kappa: float,
    xi: float,
    mass_sq: float,
    d: float,
    velocity: float,
    hubble: float,
) -> float:
    """Flat-FLRW constraint residual, including the 6*xi*H*d*velocity term."""
    _finite_inputs(hubble=hubble)
    state = evaluate_state(kappa, xi, mass_sq, d, velocity)
    potential = -0.5 * mass_sq * d * d
    return (
        3.0 * hubble * hub * state.effective_gravity_prefactor
        - kappa * (0.5 * velocity * velocity + potential + 6.0 * xi * hubble * d * velocity)
    )


def friedmann_hubble_roots(
    kappa: float,
    xi: float,
    mass_sq: float,
    d: float,
    velocity: float,
    *,
    pole_tolerance: float = 1.0e-12,
) -> tuple[float, ...]:
    """Real H roots of the full flat-FLRW quadratic constraint, if any."""
    state = evaluate_state(
        kappa, xi, mass_sq, d, velocity, pole_tolerance=pole_tolerance
    )
    potential = -0.5 * mass_sq * d * d
    a = 3.0 * state.effective_gravity_prefactor
    b = -6.0 * kappa * xi * d * velocity
    c = -kappa * (0.5 * velocity * velocity + potential)
    discriminant = b * b - 4.0 * a * c
    scale = max(1.0, b * b, abs(4.0 * a * c))
    if discriminant < -pole_tolerance * scale:
        return ()
    discriminant = max(0.0, discriminant)
    root = math.sqrt(discriminant)
    roots = ((-b - root) / (2.0 * a), (-b + root) / (2.0 * a))
    if math.isclose(roots[0], roots[1], rel_tol=1.0e-14, abs_tol=1.0e-14):
        return (roots[0],)
    return tuple(sorted(roots))


def nonzero_rest_root_d_squared(kappa: float, xi: float) -> float:
    """Return d_*^2 = -1/(kappa*xi), defined for the xi<0 reduced-force root."""
    _finite_inputs(kappa=kappa, xi=xi)
    if kappa <= 0.0:
        raise ValueError("kappa must be positive")
    if xi >= 0.0:
        raise ValueError("the nonzero real rest-force root requires xi < 0")
    return -1.0 / (kappa * xi)
