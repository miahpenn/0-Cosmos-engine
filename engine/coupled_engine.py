"""One-machine V5.5 coupled evolution driver.

The local scalar, COSMOS scalar, dark matter, and baryons share one metric
and one total stress-energy tensor.  The interface ledger observes that
solution; it does not inject a reciprocal phenomenological source.

The optional step limit is only a computational safety bound.  It is not a
physical stop condition and does not alter the equations.
"""
from dataclasses import dataclass
import numpy as np

from . import reference_pirk_unified as q
from .true_cmc_reference import cmc_lapse, step_true_cmc
from .bidirectional_production import ledger_observables
from .handoff import CycleLedger, handoff_from_ledger


@dataclass
class CoupledRun:
    grid: object
    state: object
    fields: tuple
    t: float
    tau: float
    cycle: CycleLedger
    handoffs: list
    history: list


def initialize(n=80, rmax=40.0, amplitude=.01, width=7.0):
    g = q.Grid(n, rmax)
    s, fields = q.make_initial(g, amplitude, width)
    s.alpha = cmc_lapse(g, s, fields)[0]
    s.beta.fill(0.0)
    s.B.fill(0.0)
    return CoupledRun(g, s, tuple(fields), 0.0, 0.0,
                      CycleLedger(), [], [])


def advance(run, dt):
    previous_H = -float(np.mean(run.state.K)) / 3.0
    run.state, run.fields = step_true_cmc(
        run.grid, run.state, run.fields, dt
    )
    run.t += dt
    run.tau += dt * float(run.state.alpha[0])

    obs = ledger_observables(run.grid, run.state, run.fields)
    event = run.cycle.observe(
        len(run.history), run.t, run.tau, previous_H, obs["H_eff"]
    )

    # Handoff is a measurement of the actual solved state.
    S, Sd, D, Dd = run.fields[:4]
    hp = handoff_from_ledger(
        run.t, run.tau, obs,
        float(S[0]), float(Sd[0]), float(D[0]), float(Dd[0])
    )
    run.handoffs.append(hp)
    run.history.append(obs)
    return event


def evolve(run, dt, max_steps=None):
    """Advance continuously; only non-finite state terminates evolution.

    max_steps is a host-computation safety bound, never a physics event.
    """
    steps = 0
    while max_steps is None or steps < max_steps:
        event = advance(run, dt)
        steps += 1
        if not np.all(np.isfinite(run.state.alpha)):
            return {"status": "numerical_failure", "steps": steps,
                    "event": event}
        if not np.all(np.isfinite(run.state.a)):
            return {"status": "numerical_failure", "steps": steps,
                    "event": event}
    return {"status": "host_limit", "steps": steps,
            "event": run.cycle.events[-1] if run.cycle.events else None}


def trajectory(run):
    if not run.history:
        return {}
    keys = run.history[0].keys()
    return {k: np.asarray([row[k] for row in run.history]) for k in keys}
