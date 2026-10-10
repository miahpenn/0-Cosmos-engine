# Extended-B post-fix diagnostic rerun

- Purpose: validate the patched diagnostic path and stage/cell failure attribution.
- Parent code revision: `9ecef141a1cbacba31986f5a65d9f2b41379a9a3`.
- Branch: `diag/extended-B-recovery-metric`.
- The Extended-B recovery selector remains diagnostic-only and defaults OFF in production code. The workflow explicitly enables it for this trace.
- Required gates: full repository `pytest -q`, then N=80, r_max=80, dr=1, CFL=0.03, radiation ON, D amplitude=1e-10, target t=50.
- Goal: verify exception stage, caller, grid index, and radius are recorded when radiation primitive recovery rejects the predictor state; retain RHS-budget closure checks.
- No production branch changes, equation changes, clipping, floors, damping, timestep adjustments, or resolution campaign authorization.
- Interpretation rule: a green workflow means the expected failure was captured after tests passed; it does not mean the physical evolution reached t=50.
