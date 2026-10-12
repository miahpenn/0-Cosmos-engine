# Radiation centre-state capture probe (diagnostic only)

Status: design / not yet executed  
Target branch: `diag/radiation-four-force-counterfactual-2026-10-12`  
Parent experiment: run 38315678068, commit `1042e75396d8aca889d28e0a75608b041354ae9b`

## Purpose

Localize the first-RHS radiation admissibility failure at the central cell for the already-passing N=80 cases and failing N=160 cases. This is a diagnostic capture only: no changes to the radiation equations, initial data, reconstruction, boundary treatment, recovery algorithm, guard, CFL, or coefficients.

## Fixed case matrix

Use the same default initial data, `r_max=40`, `CFL=0.03`, radiation enabled, and one ordinary `kernel.step` per case:

- N=80, beta_rad=0.0
- N=80, beta_rad=-0.04
- N=160, beta_rad=0.0
- N=160, beta_rad=-0.04

Retain the existing focused regression and complete test suite as preconditions. Preserve per-case errors and write/hash/upload the JSON even if one or more cases fail.

## Capture points

### A. Immediately before the first RHS call

At the start of `kernel.step`, before invoking `_rhs`, capture for the radiation species:

- conserved variables at cells i=0, 1, 2: `rest`, `energy_t`, `momentum_r`;
- metric values at i=0, 1, 2: `alpha`, `beta`, `a`, `b`, `X`, `K`, `Aa`, and derived `gamma_rr`, `gamma_rr_inv`, `sqrt_gamma`;
- recovered primitives at those cells: `rho`, `pressure`, `v_r`, `W`, `h=rho+pressure`, `vhat=sqrt(gamma_rr)*v_r`, and `|S|/E`;
- the actual grid centres and `dr`.

Do not infer primitive values from a separate approximate formula when the production recovery routine can be called directly on a copy of the unchanged pre-step state.

### B. During the first radiation RHS evaluation

Instrument the existing call path without changing its return values. Record:

- conserved state passed to the radiation RHS and recovered primitives at i=0, 1, 2;
- each geometric source component at i=0, including the energy and radial-momentum sources already used by the update;
- left/right reconstructed primitive states at the first interior face (between cells 0 and 1), the face metric, characteristic speeds, and the resulting HLL flux components;
- central flux boundary value (the existing regular-centre zero-flux treatment), first interior face flux, and the resulting finite-volume flux differences/tendencies at i=0;
- source and flux contributions separately, plus the post-update conserved state and admissibility ratio immediately before the existing guard raises;
- any exception raised by primitive recovery, preserving its complete traceback.

Instrumentation must be observational: wrappers log arguments/results and return the original values unchanged. No added projection, clamp, retry, alternate recovery, or altered centre parity.

### C. Comparison and interpretation

Compare N=80 and N=160 by physical cell index and normalized location, preserving full-precision values. Separate the following hypotheses rather than choosing one prematurely:

1. Initial conserved state or primitive recovery already near/over the bound.
2. Reconstruction / HLL face state or flux imbalance at the centre.
3. Geometric source discretization at i=0.
4. Primitive recovery inconsistency or conditioning.
5. Centre parity/regularity handling or grid-centre convention.

The existing exception's flux tendency and source are useful, but they are not a substitute for recording the actual face flux and pre-update conserved state.

## Admission criteria

- All four cases complete their one-step attempt and are represented in the artifact.
- Pre-step captures are explicitly distinguished from RHS-entry and rejection-point captures.
- Captures include all requested fields or an explicit capture error; missing fields must not silently be represented as zero.
- JSON SHA-256 is emitted and the artifact is uploaded regardless of case failure.
- No physics interpretation or candidate promotion follows from this probe alone.
- PR #9 remains draft; no guard relaxation or phenomenological fix.

## Known parent result

Run 38315678068 on `1042e75396d8aca889d28e0a75608b041354ae9b`:
- focused tests: 3 passed;
- existing tests: 25 passed;
- N=80 null and candidate: passed;
- N=160 null: failed at i=0 with `|S|/E ~= 1.3585922116`;
- N=160 candidate: failed at i=0 with `|S|/E ~= 1.3584029299`;
- JSON SHA-256: `f5b9cdf9e0e267249821adc8c87d04c860d719ebdd6aee58bfff7a74c2b6ad62`.

This parent result identifies a resolution-dependent first-RHS admissibility failure but does not yet establish which sublayer is responsible.
