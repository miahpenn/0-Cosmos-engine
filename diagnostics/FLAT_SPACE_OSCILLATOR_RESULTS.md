# Flat-space scalar oscillator smoke test

## Status

- Diagnostic branch: `diag/bond-sweep-self-consistent`
- Initial workflow: [GitHub Actions run 38201499478](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38201499478)
- Initial artifact: [flat-oscillator-smoke](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38201499478/artifacts/11725521660)
- Initial run result: **PASS** for the scalar RHS under a fixed flat background, using an independent RK4 driver.
- Follow-up code correction: JSON output newline was changed from a literal backslash-n suffix to a real newline. A clean rerun should be recorded before treating the artifact as the canonical machine-readable result.

## Physical setup and analytic answer

Set the geometry to flat and static (`alpha=1`, `beta=0`, `K=0`, `a=b=X=1`). Set the S field spatially uniform, with `S(0)=A=10^-3` and `PS(0)=0`; set the other scalar fields to zero.

For the implemented S subsystem, the spatial derivatives and curvature terms vanish, leaving

- `dS/dt = PS`
- `dPS/dt = -S`

The analytic solution is `S(t)=A cos(t)`, `PS(t)=-A sin(t)`. After one period, `T=2 pi`, the expected state returns to `S=A`, `PS=0`.

## Initial observed results

| Steps per period | Maximum absolute error | Error reduction vs. next finer run |
|---:|---:|---:|
| 400 | 3.18743e-12 | 15.9986 |
| 800 | 1.99231e-13 | 15.9971 |
| 1600 | 1.24542e-14 | — |

Halving the timestep reduced the error by about 16 each time, consistent with fourth-order convergence for the RK4 integration. The numerical state returned close to `S=0.001`, `PS=0`.

## What this establishes

It is a useful, inexpensive unit-level check that the scalar RHS gives the expected harmonic-oscillator dynamics in this homogeneous, fixed-background limit and that the chosen test driver converges at fourth order.

## What it does not establish

This test **does not validate the production time integrator, evolving geometry, Hamiltonian/momentum constraint propagation, inhomogeneous scalar gradients, matter coupling, the S-D bond, or cosmological evolution**. The RK4 driver is independent of the engine's full production evolution path; that makes the test a focused check of the RHS, not a full-machine certification.

## Next cheap known-physics tests

1. Run a clean repeat after the JSON newline fix and preserve its artifact.
2. Test exact static flat vacuum through the production geometry update path, with all matter and couplings off; record metric/lapse drift and both constraints.
3. If that passes, test a small-amplitude inhomogeneous scalar mode on fixed flat geometry against a known linear-wave solution, then test a simple coupled geometry case.

Do not infer anything about the 0* hypothesis from this smoke test. Keep all changes on the diagnostic branch; production `main` remains untouched.
