"""One-machine V5.5 coupled evolution driver.

The local scalar, COSMOS scalar, dark matter, and baryons share one metric
and one total stress-energy tensor.  The interface ledger observes that
solution; it does not inject a reciprocal phenomenological source.

Evolution, proper time, the cycle ledger and handoffs are those of
``V55ProductionKernel``; this module is a thin driver over it and owns no
equations. It formerly stepped its own copy of the V5.5 system
(``reference_pirk_unified`` + ``true_cmc_reference``), kept for provenance only;
see docs/PROVENANCE_OLDER_V55_PATHS.md.

The optional step limit is only a computational safety bound.  It is not a
physical stop condition and does not alter the equations.
"""
from dataclasses import dataclass
import numpy as np

from .production_kernel import V55ProductionKernel
from .bidirectional_production import ledger_observables


@dataclass
class CoupledRun:
    kernel: V55ProductionKernel
    state: object

    @property
    def grid(self):
        return self.state.grid

    @property
    def t(self):
        return self.state.t

    @property
    def tau(self):
        return self.state.tau

    @property
    def cycle(self):
        return self.state.cycle

    @property
    def handoffs(self):
        return self.state.handoffs

    @property
    def history(self):
        return self.state.history


def initialize(n=80, rmax=40.0, amplitude=.01, width=7.0, **initial_data):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=n, r_max=rmax, amplitude=amplitude, width=width, **initial_data
    )
    return CoupledRun(kernel, state)


def advance(run, dt):
    """Take one production step; return the cycle event it observed, if any."""
    events_before = len(run.state.cycle.events)
    run.state = run.kernel.step(run.state, dt)
    events = run.state.cycle.events
    return events[-1] if len(events) > events_before else None


def evolve(run, dt, max_steps=None):
    """Advance continuously; only non-finite state terminates evolution.

    max_steps is a host-computation safety bound, never a physics event.
    """
    steps = 0
    event = None
    while max_steps is None or steps < max_steps:
        event = advance(run, dt)
        steps += 1
        geom = run.state.geometry
        if not np.all(np.isfinite(geom.alpha)) or not np.all(np.isfinite(geom.a)):
            return {"status": "numerical_failure", "steps": steps,
                    "event": event}
    return {"status": "host_limit", "steps": steps,
            "event": run.cycle.events[-1] if run.cycle.events else None}


def ledger(run):
    """Interface-ledger row of the run's current state."""
    return ledger_observables(run.kernel, run.state)


def trajectory(run):
    if not run.history:
        return {}
    keys = [k for k in run.history[0].keys()
            if all(np.isscalar(row.get(k)) for row in run.history)]
    return {k: np.asarray([row[k] for row in run.history]) for k in keys}
