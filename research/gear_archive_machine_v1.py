"""Research-only top-level runner for the archived two-lane GEAR machine.

The default local lane follows the reduced local-mode equations in the preserved
GEAR archive DOCX; an explicit comparison mode calls the current engine.local RHS.
The homogeneous COSMOS lane uses the corrected, normalized engine.cosmos equations.
No cross-lane source is introduced. Local time tau and cosmic time t stay separate.
This is a bounded research harness, not a production or admission claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import atan2, isclose, isfinite, sqrt
from typing import Callable, Sequence

from engine.cosmos import (
    CosmosParams,
    CosmosState,
    construct_present_day,
    cosmos_rhs,
    cosmos_rho_p,
    normalized_budgets,
)
from engine.local import LocalParams, local_rhs


ARCHIVE_OPEN_GATES = {
    "local_proper_time_to_cosmic_time": "OPEN — no derived clock map",
    "D_mode_to_localized_matter": "OPEN — no conversion/source law",
    "local_source_to_global_cosmos_Q": "OPEN — no reciprocal Q law",
    "full_GR_oscillaton_endpoint": "OPEN — not represented by this reduced runner",
    "first_principles_bounce": "OPEN — no forced or assumed bounce",
    "global_anisotropic_three_curvature": "OPEN — not represented by exact spherical lane",
}

DEFAULT_LOCAL_INITIAL_STATE = (0.02, 0.0, 1.0e-8, 0.0)


def archive_docx_local_rhs(y: Sequence[float], p: LocalParams) -> tuple[float, ...]:
    """Reduced local RHS transcribed from the preserved GEAR archive DOCX.

    This is the archive's mode-space equation, not the full spherical PDE:
    Sdd = -(kappa+g) S
    Ddd = (g-kappa) D - eta_G D**2 Ddot
          - (2*pi/3) G (g-kappa) D**3
    """
    S, Sd, D, Dd = _finite_vector(y, name="archive local state", size=4)
    omega2 = p.kappa + p.g
    alpha2 = p.g - p.kappa
    barrier = (2.0 * 3.141592653589793 / 3.0) * p.G * alpha2 * D**3
    damping = p.eta_G * D**2 * Dd
    return (Sd, -omega2 * S, Dd, alpha2 * D - damping - barrier)


@dataclass(frozen=True)
class LaneTrace:
    """Immutable sampled trace with its own lane-specific time coordinate."""

    time_name: str
    time: tuple[float, ...]
    state: tuple[tuple[float, ...], ...]


def _finite_vector(values: Sequence[float], *, name: str, size: int) -> tuple[float, ...]:
    out = tuple(float(v) for v in values)
    if len(out) != size:
        raise ValueError(f"{name} must contain exactly {size} values")
    if not all(isfinite(v) for v in out):
        raise ValueError(f"{name} values must all be finite")
    return out


def _validate_time_controls(final_time: float, dt: float, sample_every: int, lane: str) -> None:
    if not isfinite(final_time) or final_time < 0.0:
        raise ValueError(f"{lane} final time must be finite and nonnegative")
    if not isfinite(dt) or dt <= 0.0:
        raise ValueError(f"{lane} dt must be finite and positive")
    if isinstance(sample_every, bool) or not isinstance(sample_every, int) or sample_every < 1:
        raise ValueError("sample_every must be an integer >= 1")


def _add_scaled(y: tuple[float, ...], k: tuple[float, ...], scale: float) -> tuple[float, ...]:
    return tuple(yi + scale * ki for yi, ki in zip(y, k, strict=True))


def rk4_step(
    rhs: Callable[[tuple[float, ...]], Sequence[float]],
    state: Sequence[float],
    dt: float,
) -> tuple[float, ...]:
    """One checked classical RK4 step; this helper adds no model terms."""
    y = _finite_vector(state, name="state", size=len(state))
    if not isfinite(dt) or dt <= 0.0:
        raise ValueError("dt must be finite and positive")

    def evaluate(x: tuple[float, ...]) -> tuple[float, ...]:
        return _finite_vector(rhs(x), name="RHS", size=len(y))

    k1 = evaluate(y)
    k2 = evaluate(_add_scaled(y, k1, 0.5 * dt))
    k3 = evaluate(_add_scaled(y, k2, 0.5 * dt))
    k4 = evaluate(_add_scaled(y, k3, dt))
    result = tuple(
        yi + (dt / 6.0) * (a + 2.0 * b + 2.0 * c + d)
        for yi, a, b, c, d in zip(y, k1, k2, k3, k4, strict=True)
    )
    return _finite_vector(result, name="RK4 result", size=len(y))


def _integrate_lane(
    initial: Sequence[float],
    rhs: Callable[[tuple[float, ...]], Sequence[float]],
    *,
    final_time: float,
    dt: float,
    sample_every: int,
    time_name: str,
) -> LaneTrace:
    y = _finite_vector(initial, name=f"{time_name} initial state", size=len(initial))
    times = [0.0]
    states = [y]
    t = 0.0
    nstep = 0

    while t < final_time:
        remaining = final_time - t
        h = min(dt, remaining)
        next_y = rk4_step(rhs, y, h)
        next_t = final_time if h >= remaining else t + h
        if not next_t > t:
            raise FloatingPointError(f"{time_name} integration step did not advance time")
        y = next_y
        t = next_t
        nstep += 1
        if nstep % sample_every == 0 or t >= final_time:
            times.append(t)
            states.append(y)

    return LaneTrace(time_name=time_name, time=tuple(times), state=tuple(states))


def _local_energy_ledger(y: Sequence[float], p: LocalParams) -> dict[str, float]:
    """Archive mode-energy split; positive excursion energy is only a proxy."""
    S, Sd, D, Dd = y
    omega2 = p.kappa + p.g
    alpha2 = p.g - p.kappa
    e_s = 0.5 * Sd * Sd + 0.5 * omega2 * S * S
    e_d_raw = 0.5 * Dd * Dd - 0.5 * alpha2 * D * D
    e_d_kinetic = 0.5 * Dd * Dd
    e_d_excursion = 0.5 * alpha2 * D * D
    return {
        "E_S": e_s,
        "E_D_raw": e_d_raw,
        "E_D_kinetic": e_d_kinetic,
        "E_D_excursion": e_d_excursion,
        "E_total_raw": e_s + e_d_raw,
        "generated_activity_proxy": e_d_kinetic + e_d_excursion,
        "transport_activity_proxy": abs(D * Dd),
    }


def _local_cycle_phase(trace: LaneTrace, p: LocalParams) -> tuple[float, ...]:
    omega2 = p.kappa + p.g
    if omega2 <= 0.0:
        raise ValueError("GEAR-03 local phase requires kappa + g > 0")
    omega = sqrt(omega2)
    phases = [atan2(-y[1] / omega, y[0]) for y in trace.state]
    out = [phases[0]]
    for phase in phases[1:]:
        delta = phase - out[-1]
        while delta > 3.141592653589793:
            phase -= 2.0 * 3.141592653589793
            delta = phase - out[-1]
        while delta < -3.141592653589793:
            phase += 2.0 * 3.141592653589793
            delta = phase - out[-1]
        out.append(phase)
    return tuple(out)


def _cosmos_state_tuple(state: CosmosState) -> tuple[float, ...]:
    return (
        state.a,
        state.H,
        state.phi,
        state.pi_phi,
        state.rho_dm,
        state.rho_b,
        state.rho_r,
        state.rho_shear,
    )


def _cosmos_state_dict(y: Sequence[float]) -> dict[str, float]:
    keys = ("a", "H", "phi", "pi_phi", "rho_dm", "rho_b", "rho_r", "rho_shear")
    return dict(zip(keys, (float(v) for v in y), strict=True))


def run_machine(
    *,
    local_tau_final: float = 0.1,
    local_dt: float = 1.0e-3,
    cosmos_t_final: float = 0.1,
    cosmos_dt: float = 1.0e-3,
    cosmic_branch: int = 1,
    sample_every: int = 1,
    local_params: LocalParams | None = None,
    cosmos_params: CosmosParams | None = None,
    local_initial_state: Sequence[float] = DEFAULT_LOCAL_INITIAL_STATE,
    local_equation_set: str = "archive_docx",
) -> dict[str, object]:
    """Run the archived local and homogeneous lanes independently.

    cosmic_branch selects only the initial sign of H (+1 expansion, -1
    contraction); thereafter H is integrated from cosmos_rhs. No clock
    synchronization, exchange_Q, D-to-matter conversion, or bounce is added.
    """
    _validate_time_controls(local_tau_final, local_dt, sample_every, "local")
    _validate_time_controls(cosmos_t_final, cosmos_dt, sample_every, "COSMOS")
    if isinstance(cosmic_branch, bool) or not isinstance(cosmic_branch, int) or cosmic_branch not in (-1, 1):
        raise ValueError("cosmic_branch must be +1 (expansion) or -1 (contraction)")
    if local_equation_set not in ("archive_docx", "current_source"):
        raise ValueError("local_equation_set must be 'archive_docx' or 'current_source'")

    lp = local_params or LocalParams()
    cp = cosmos_params or CosmosParams()
    y0_local = _finite_vector(local_initial_state, name="local_initial_state", size=4)

    omega2 = lp.kappa + lp.g
    alpha2 = lp.g - lp.kappa
    if omega2 <= 0.0 or alpha2 <= 0.0:
        raise ValueError("archive GEAR-03 lock requires kappa+g > 0 and g-kappa > 0")
    if not isclose(lp.alpha_D * lp.alpha_D, alpha2, rel_tol=1.0e-12, abs_tol=1.0e-14):
        raise ValueError("archive GEAR-03 lock requires alpha_D**2 == g-kappa")

    if local_equation_set == "archive_docx":
        selected_local_rhs = lambda y: archive_docx_local_rhs(y, lp)
        local_equation_source = "preserved GEAR archive DOCX reduced local_rhs"
    else:
        selected_local_rhs = lambda y: local_rhs(y, lp, H=0.0)
        local_equation_source = "engine.local.local_rhs(..., H=0.0); eta_G unused by source"

    local_trace = _integrate_lane(
        y0_local,
        selected_local_rhs,
        final_time=local_tau_final,
        dt=local_dt,
        sample_every=sample_every,
        time_name="tau_local",
    )

    cosmic0_obj = construct_present_day(cp)
    cosmic0 = _cosmos_state_tuple(cosmic0_obj)
    cosmic0 = (cosmic0[0], cosmic_branch * abs(cosmic0[1]), *cosmic0[2:])
    cosmos_trace = _integrate_lane(
        cosmic0,
        lambda y: cosmos_rhs(y, cp),
        final_time=cosmos_t_final,
        dt=cosmos_dt,
        sample_every=sample_every,
        time_name="t_cosmos",
    )

    local_states = [
        {"S": y[0], "Sd": y[1], "D": y[2], "Dd": y[3]}
        for y in local_trace.state
    ]
    local_ledgers = [_local_energy_ledger(y, lp) for y in local_trace.state]

    cosmic_states = [_cosmos_state_dict(y) for y in cosmos_trace.state]
    cosmic_ledgers = []
    for y in cosmos_trace.state:
        rho, pressure = cosmos_rho_p(y, cp)
        state = CosmosState(*y)
        cosmic_ledgers.append({
            "rho_total": float(rho),
            "pressure_total": float(pressure),
            "friedmann_residual_H2_minus_rho_over_3": float(y[1] * y[1] - rho / 3.0),
            "normalized_budgets": normalized_budgets(state, cp),
        })

    return {
        "schema_version": 1,
        "status": "research_only_independent_lanes",
        "scientific_admission": "not_claimed",
        "lane_contract": {
            "local_time": "tau_local; evolved independently",
            "cosmos_time": "t_cosmos; evolved independently",
            "cross_lane_source_terms": False,
            "cosmic_branch_selection": "sets only initial sign of H; H then evolves",
        },
        "open_gates": dict(ARCHIVE_OPEN_GATES),
        "local": {
            "equation_set": local_equation_set,
            "equation_source": local_equation_source,
            "time_name": local_trace.time_name,
            "time": list(local_trace.time),
            "state": local_states,
            "energy_and_activity_ledgers": local_ledgers,
            "oscillator_phase_theta": list(_local_cycle_phase(local_trace, lp)),
        },
        "cosmos": {
            "time_name": cosmos_trace.time_name,
            "time": list(cosmos_trace.time),
            "state": cosmic_states,
            "component_and_constraint_ledgers": cosmic_ledgers,
            "equation_source": "engine.cosmos.cosmos_rhs",
            "initial_branch_sign": cosmic_branch,
            "normalization_source": "engine.cosmos.construct_present_day",
        },
        "provenance": {
            "local_equation_module": "engine/local.py",
            "cosmos_equation_module": "engine/cosmos.py",
            "machine_module": "research/gear_archive_machine_v1.py",
            "note": "The archived DOCX baseline is preserved; this runner composes the current source-mapped sector RHS modules and does not alter them.",
        },
    }


__all__ = [
    "ARCHIVE_OPEN_GATES",
    "DEFAULT_LOCAL_INITIAL_STATE",
    "LaneTrace",
    "archive_docx_local_rhs",
    "rk4_step",
    "run_machine",
]
