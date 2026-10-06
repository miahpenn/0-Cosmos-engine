# 0* D-mode synchronization — housekeeping status
## 2026-10-06

### Repository state
- Experimental branch: `0star-central-clock`
- Pull request: #1, intentionally **draft** and not merged.
- Production branch remains unchanged by the central-clock experiments.
- No new physical source, fitted coefficient, damping, floor, boundary condition, bounce rule, or reset mechanism has been introduced.

### Established result
The D-mode source-order and phase-alignment campaigns support a phase-coordinate interpretation:
- different initial inverted-D amplitudes reach nearly the same local geometric/dynamical state at different coordinate/proper times;
- proper-time phase offsets were nearly constant across multiple D_active levels;
- the common-state geometry collapsed much more tightly than the coordinate lapse;
- this does not yet establish a complete feedback loop, turnaround/bounce, or 0* cycle.

Reference:
- `docs/0STAR_D_MODE_TIME_SHAPE_RESULT_2026-10-06.md`
- `docs/0STAR_D_MODE_PHASE_ALIGNMENT_2026-10-06.md`

### Current tests
Two follow-up campaigns are the active discriminators:
1. **Early-phase long synchronization:** Dhalf/Dbase/Ddouble, synchronize at D_active=0.002, then continue to t=40.
2. **Early-phase resolution check:** Dbase at N=80/120/160, synchronize at D_active=0.002, then continue to t=30.

The synchronization is a diagnostic rephasing of already evolved trajectories. It does not modify the state or dynamics.

### Workflow housekeeping
The experimental synchronization workflows are now **manual-dispatch only**. They no longer trigger automatically on every branch commit. This prevents duplicate expensive campaigns while preserving the ability to launch the intended test explicitly.

The current running jobs were triggered before this housekeeping change and are left untouched.

### Interpretation gate
The next decision depends on the completed early-sync runs:
- persistent post-synchronization collapse -> strengthen the phase-trajectory interpretation and test its domain;
- measurable divergence -> identify the first state variable carrying memory;
- resolution dependence -> stop physical interpretation and isolate the numerical layer.

No 0* gate is being paid or declared by these synchronization tests. 0* remains the independent diagnostic/scoring layer.
