# Canonical production proper-time admission request

This file triggers a controlled numerical run; it does not modify equations.

- Source branch: diag/production-proper-time-admission
- Configuration: N=120, r_max=120.0, dr=1.0, CFL=0.03
- Requested central proper-time target: tau=22.1
- Coordinate-time safety cap: t=100
- Production entry point: engine.run_production.run_campaign
- Radiation admissibility retry: unchanged from current canonical runner; timestep bisection only, no clipping.
- Required evidence: native per-step ledger, checkpoints, explicit result classification for target reached versus coordinate cap reached first, and SHA-verified artifact.
- Purpose: determine whether the full production runner reproduces the same outer-boundary radiation failure seen in the proper-time-matched trace. This is not a physics change or production promotion.
