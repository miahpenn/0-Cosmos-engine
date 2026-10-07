# V5.5 strong-field gate — 2026-10-07

## Objective

Resolve the numerical strong-field representation failure without changing the V5.5 matter model, stress-energy definitions, interface semantics, or physical source terms.

## Captured defect

The pre-repair moving-puncture PIRK campaigns were captured through the first numerical failure.

At CFL=0.015, N=40 failed at t≈27.18 when the lapse became nonpositive. At CFL=0.015, N=60 (dr=2/3, so dt=0.01) failed at t≈27.36.

Immediately before failure:

| N | t | min alpha | location | max H | outer H L2 | det constraint |
|---:|---:|---:|---:|---:|---:|---:|
| 40 | 27.18 | 1.76e-21 | r=0.5 | 7.34 | 2.85 | ~3.3e-16 |
| 60 | 27.36 | 7.24e-22 | r=0.333 | 6.15 | 1.97 | ~3.3e-16 |

The dominant constraint growth was in the outer shell while the conformal determinant identity stayed exact. Radiation remained physically admissible in the captured ledgers. The failure is therefore treated as a moving-gauge/finite-radius representation problem, not a radiation-model failure.

## Repairs

1. The legacy finite-radius characteristic reconstruction was removed from the moving-gauge PIRK step. That reconstruction belongs to the archived nonadvective-lapse/algebraic-B characteristic system and was not allowed to be mixed with advective 1+log plus independently evolved B.

2. The moving-gauge centre projection now preserves the independently evolved B field and imposes the required odd spherical regularity on B directly, instead of overwriting B with the algebraic B=3/4 Lambda relation.

3. Real RK-stage matter admissibility is now checked after physical stage updates. No clipping or matter repair is performed.

4. A separate stage-aware CMC kernel was added. It solves the archive-anchored elliptic lapse at the beginning, predictor, and completed stages while retaining the current V5.5 conservative matter system. The first gauge gate holds beta=0 to isolate the foliation layer.

## Gates now running

### Repaired moving gauge

N=40/60, Rmax=40, CFL=0.03, final_time=30.

Workflow:
https://github.com/miahpenn/0-Cosmos-engine/actions/workflows/zero_star_moving_pirk_repaired_gate.yml

### Stage-aware true CMC

N=40/60/80, Rmax=40, CFL=0.06, final_time=24.

Workflow:
https://github.com/miahpenn/0-Cosmos-engine/actions/workflows/zero_star_true_cmc_pirk_gate.yml

### Next gate, prepared but not launched

Same-CFL N=40/60/80 CMC evolution to t=40, requiring observed H_eff turnaround without manual sign reversal.

Workflow:
https://github.com/miahpenn/0-Cosmos-engine/actions/workflows/zero_star_true_cmc_turnaround.yml

## Acceptance rules

No fitted coefficient, artificial source, lapse floor, shell, bounce, reset, stop command, branch sign flip, or silent matter reinterpretation is permitted.

A numerical continuation is not by itself evidence of a physical bounce, completed cycle, black hole, or singularity. Those require independent invariant and resolution evidence.

## Remaining sequence

1. Complete the repaired moving-gauge and true-CMC t=24 gates.
2. Select the representation that is numerically healthy by controlled resolution evidence, not by visual continuation time.
3. Run same-CFL CMC turnaround convergence.
4. Extend through the post-turnaround regime with Hamiltonian/momentum/connection/determinant and Misner-Sharp checks.
5. Reconnect the full bidirectional COSMOS handoff only after the strong-field representation passes.
6. Run the unrestricted cycle search with no imposed bounce or stop condition.
