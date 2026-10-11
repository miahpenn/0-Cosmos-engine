# D-mode time-shape source-order campaign — preregistration

Status: preregistered before the campaign trigger commit.

## Question
Does the existing D-mode source organize the interior clock/lapse response more tightly when trajectories are aligned by instantaneous D-active state than by coordinate time, and what is the measured ordering among D growth, existing sector source contributions, the matter/geometry K-source split, CMC response, lapse gradient, and central proper-time rate?

## Fixed configuration
- Production branch: `0star-central-clock`
- Existing production kernel and CMC gauge; no evolution-equation, boundary, or physical-parameter changes
- Resolution: N=160
- Domain: Rmax=160
- Coordinate CFL: 0.0075
- Final coordinate time: 22.5
- Diagnostic checkpoint interval: 0.25
- Radiation: enabled
- Sequential initial-D cases: D0=0, Dhalf=5e-11, Dbase=1e-10, Ddouble=2e-10
- D0 remains the mandatory control.

## Measurements
Record existing state/source diagnostics at the center and the instantaneous maximum-|alpha_r| probe, including D, PD, D-active, existing S/D/COSMOS scalar-sector stress projections, dark-matter/baryon/radiation projections, vacuum and matter contributions to the existing geometry K source, CMC Kdot, H/H_eff, lapse and lapse gradient, and central proper-time rate.

## Analysis fixed before execution
1. Compare histories against coordinate time.
2. Compare the same histories against instantaneous D-active state; do not fit a new coefficient or impose a source law.
3. Determine onset ordering after the run from recorded histories. Any onset thresholds are analysis markers only and never feed back into evolution.
4. Retain all four cases and raw ledgers. Report missing/nonfinite data or failed runs rather than silently excluding them.
5. If the D-state alignment and ordering are ambiguous, report that ambiguity; do not infer one-way causality from correlation alone.

## Admission and interpretation limits
The workflow's repository tests must pass and all four trajectories must complete for a complete campaign. Completion alone does not admit a physical claim. This is a diagnostic source decomposition at one resolution and one domain; it cannot establish convergence, causal propagation, a bounce, a turnaround, or a trapped surface. No physics change follows automatically from its outcome.
