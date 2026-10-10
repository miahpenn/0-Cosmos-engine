"""Capture failure-path diagnostics for the opt-in N=320 initial-data solve.

This is an initial-data-only diagnostic. It does not evolve the state and does
not relax the solver tolerance or any admission criterion. A completed report
means diagnostics were captured, not that the initial data were admitted.
"""
import json
import sys

from engine.discrete_consistent_initial import (
    ConvergenceError,
    discrete_consistent_state,
)
from engine.production_kernel import V55ProductionKernel


def main():
    settings = {
        "resolution": 320,
        "r_max": 40.0,
        "amplitude": 0.01,
        "width": 7.0,
        "D_amplitude": 1.0e-10,
        "include_radiation": True,
        "max_iter": 20,
        "tol": 1.0e-13,
        "max_cond": 1.0e12,
    }
    print("[INITIAL_SOLVER_DIAGNOSTIC] mode=initial_data_only")
    print("[INITIAL_SOLVER_DIAGNOSTIC] admission_tolerance_unchanged=1e-13")
    print("[INITIAL_SOLVER_DIAGNOSTIC] evolution_started=false")
    print("[INITIAL_SOLVER_DIAGNOSTIC] settings=" + json.dumps(settings, sort_keys=True))

    try:
        _, _, residual_history, info = discrete_consistent_state(
            V55ProductionKernel(), return_info=True, **settings
        )
    except ConvergenceError as exc:
        print("[INITIAL_SOLVER_DIAGNOSTIC] outcome=NOT_ADMITTED")
        print("[INITIAL_SOLVER_DIAGNOSTIC] exception=" + str(exc))
        print("[INITIAL_SOLVER_DIAGNOSTIC] iteration_records=")
        print(json.dumps(exc.diagnostics, indent=2, sort_keys=True))
        print("[INITIAL_SOLVER_DIAGNOSTIC] interpretation=solver failure captured; no physical conclusion")
        return 0

    print("[INITIAL_SOLVER_DIAGNOSTIC] outcome=CONVERGED")
    print("[INITIAL_SOLVER_DIAGNOSTIC] info=" + json.dumps(info, sort_keys=True))
    print("[INITIAL_SOLVER_DIAGNOSTIC] residual_history=" + json.dumps(residual_history))
    print("[INITIAL_SOLVER_DIAGNOSTIC] admission=solver tolerance passed; downstream gates still required")
    return 0


if __name__ == "__main__":
    sys.exit(main())
