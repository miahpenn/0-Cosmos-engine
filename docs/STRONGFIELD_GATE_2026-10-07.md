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


## 2026-10-06 CMC operator repair

The first repository true-CMC gate reached t=24 at N=40/60/80, but its Hamiltonian maximum was
center-dominated and increased with resolution. Audit of the actual numerical operators identified
the gauge-layer defect: the CMC elliptic solve used a mixed second-order three-point interior
operator while the evolved K equation used the pinned fourth-order cell-centered derivative,
with a separate first-cell alpha equality row. The two discrete operators therefore did not
enforce the same CMC condition.

The repair replaces that mixed operator with the native fourth-order/parity derivative rows,
including the regular first cell and the one-sided penultimate cell, while retaining the outer
alpha(R)=1 normalization. No physical source, equation of state, coefficient, lapse floor,
boundary prescription, bounce rule, or interface term was added.

The isolated true-CMC kernel now calls the same shared production CMC solver, so the production
and diagnostic branches cannot silently diverge in their gauge operator or CMC target.

The next controlled gate is the same N=40/60/80, Rmax=40, CFL=0.06, t=24 campaign. Promotion
requires actual resolution improvement in the invariant/constraint witnesses, not merely a
successful run.

## 2026-10-07: first same-CFL turnaround result

The first same-CFL True-CMC turnaround campaign used Rmax=40, CFL=0.06 and N=40/60/80. All three resolutions crossed H_eff=0 before the later finite-radius breakdown:

| N | turnaround t | turnaround tau | later numerical boundary |
|---:|---:|---:|---|
| 40 | 27.6000 | 21.21158 | t≈34.44 |
| 60 | 27.5600 | 21.21044 | t≈34.28 |
| 80 | 27.5700 | 21.21253 | t≈34.23 |

The turnaround times agree to approximately 0.15% across the three resolutions. This is a strong same-CFL dynamical convergence witness for the H_eff=0 crossing itself.

The later failures are qualitatively different from the earlier lapse-collapse failure. The lapse remains positive, while the outer retained shell develops very large conformal-geometry derivatives and Hamiltonian/momentum/connection residuals. The failure is therefore classified as finite-radius outer-domain contamination, not as a failed turnaround or a demonstrated physical singularity.

At the recorded failure states, the radiation characteristic state remains admissible at the outermost cell. The divergent constraint terms are dominated by the outer geometry variables. This does not justify a physical interpretation of the late trapped roots.

The next gate enlarges Rmax from 40 to 80 while preserving the same physical spacing: N=80/120/160 gives dr=1, 2/3, 1/2, matching the earlier N=40/60/80 spacings. This directly tests whether the late failure follows the finite outer boundary.

## 2026-10-07: event-ledger repair

The production kernel previously appended a HandoffPoint at every timestep. That was a bookkeeping defect: a handoff is defined by an actual cycle crossing, while the continuous trajectory belongs in history. The kernel now records handoffs only when CycleLedger.observe() emits a real turnaround or re-expansion event.