# Proper-time matched radiation-domain comparison request

This file is an execution trigger, not a physics or numerical-source change.

Frozen run:
- Branch: diag/radiation-proper-time-domain-control
- Resolution: N=120
- Domain: r_max=120.0, dr=1.0
- CFL: 0.03
- Central proper-time target: tau=22.1
- Coordinate-time safety cap: t=100
- Radiation admissibility retry: enabled, timestep bisection only
- Radiation recovery metric: original stage/predictor metric
- Expected result classification:
  - target_proper_time_completed_without_failure
  - radiation_admissibility_step_underflow / expected failure
  - coordinate_cap_reached_before_proper_time_target
- No clipping, floors, artificial damping, source edits, or physics changes.
- Purpose: compare the expanded domain at the same central proper time as the prior r_max=80 failure near tau=22.02484.
