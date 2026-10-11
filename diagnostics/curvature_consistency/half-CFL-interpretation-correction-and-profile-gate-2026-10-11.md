# Half-CFL interpretation correction and domain-profile decision gate

**Date:** 2026-10-11  
**Branch:** `diag/independent-metric-ricci-audit`  
**Status:** interpretation correction; diagnostic only.

## Review from Lumen

Lumen independently checked the half-CFL report against the repository. The reported matched-time magnitudes and small timestep dependence are consistent with the source data, but two interpretive limits need to be explicit.

### 1. Late-time agreement does not erase initializer dependence

At (t=4), off-centre max-|H| is about (1.24\times10^{-4}) for production and (4.8\times10^{-6}) for the regular-F shadow (roughly 25-fold smaller in the shadow). At (t=8), production is about (1.01\times10^{-4}), while the shadow is about (5.94\times10^{-5}) (roughly six-fold smaller). At (t=12), both are about (2.829\times10^{-4}), agreeing to about 0.04%.

Permitted conclusion: timestep dependence is small at the sampled targets; the off-centre maxima happen to agree closely at (t=12). Do **not** conclude that the initializer does not matter. Possible late convergence to a common state, lapse collapse, or coincidental equality remain open.

### 2. A peak-radius jump does not establish a moving front

The half-CFL report gives the shadow's full-domain maximum at (r=6.25) at (t=8), while production's maximum at the same time is at the inner peak (r=1.25). Production's maximum later appears at (r=7.25) near (t=9). These values are compatible with an argmax switch between an inner and outer local peak; they do not prove a single feature propagated continuously.

The half-CFL artifact contains sparse-time summaries, not full radial residual profiles. Therefore it cannot settle the peak-identity question.

## Domain-size workflow status

Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38108719482  
Commit: `9087dbea4c8507894088cc2603c6608b6cfc8623`

At the latest status check, the workflow's validation, compile, fixed-spacing checks and dependency installation had passed. The four-trajectory numerical evolution step was still in progress, and no report artifact was yet published. Do not describe the domain run as complete until GitHub reports a terminal result and its artifact is inspected.

## Required profile analysis when artifact is ready

For each initializer and each domain ((N=160,r_{max}=80); (N=320,r_{max}=160), both (Delta r=0.5)):

1. Read full (H(r,t)) and lapse profiles every 0.25 time units, not just a global argmax.
2. Identify inner and outer local peaks separately from at least (t=4) onward; record their radii, amplitudes, prominence and evolution.
3. Track (alpha(0,t)) and (alpha_{min}(t)) beside both peaks, since the elliptic lapse solution may itself depend on the outer domain.
4. Compare the two domains at matched actual times. Distinguish genuine domain dependence from a maximum switching between peaks.
5. Interpret the initializer comparison independently of the timestep comparison.

Diagnostic signatures:
- Similar outer-peak emergence in both initializers supports a shared feature, but does not prove initializer independence in the full trajectory.
- Different onset times/radii support sensitivity to the initial state.
- A decaying inner peak plus a separately growing outer peak supports an argmax-switch explanation rather than a single travelling maximum.
- Domain-dependent (alpha(0,t)) indicates that the elliptic lapse solution changed with the domain and complicates direct causal attribution.

## Decision

No physics patch. No bounce/turnaround or travelling-front claim. The narrow finding is: half-CFL changes the sampled off-centre magnitudes by less than about 0.2%, but initializer dependence remains substantial at (t=4) and (t=8), and the peak trajectory remains unresolved until full profiles are examined.
