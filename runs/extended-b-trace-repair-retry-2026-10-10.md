# Extended-B post-fix test and trace retry

- Parent code revision: `7750b3980db8df34ac3e75099ae677684141008c`.
- Branch: `diag/extended-B-recovery-metric`.
- Reason for retry: the previous trigger commit was created before the recovery-stash regression fixture was corrected to use a mutable geometry object implementing `copy()`. The full test gate must be green before the numerical trace step can execute.
- Current trace instrumentation also records the throwing RHS stage, caller, cell index, and radius; the focused test guards that contract.
- Run setup remains N=80, r_max=80, dr=1, CFL=0.03, radiation ON, D amplitude=1e-10, target t=50, Extended-B explicitly ON only inside this diagnostic workflow.
- No production branch change, equation change, clipping, floor, damping, timestep adjustment, or resolution-study authorization.
- A successful workflow means the expected radiation failure was captured after pytest passed; it does not mean the evolution reached t=50.
