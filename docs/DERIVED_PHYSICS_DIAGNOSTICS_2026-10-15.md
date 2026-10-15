# Derived-physics diagnostic helpers — research-only

**Date:** 15 October 2026  
**Status:** Isolated algebraic utilities and unit tests; no production engine integration.  
**Main/source anchor when claimed:** \`89b577e2cd98d755d653c0f2d32104a15f0f571c\`.  
**Related audit records:** [theorem audit #260](https://github.com/miahpenn/0-Cosmos-engine/issues/260), [A6 derivation PR #261](https://github.com/miahpenn/0-Cosmos-engine/pull/261), [A6 numerical benchmark PR #263](https://github.com/miahpenn/0-Cosmos-engine/pull/263), and [A10/A11 scope PR #262](https://github.com/miahpenn/0-Cosmos-engine/pull/262).

## Purpose

\`tools/derived_physics_diagnostics.py\` provides pure, standard-library-only calculations for selected equations independently recovered in the theorem audit. It is intended as a small algebraic cross-checking instrument, not as a solver, model selector, or statement that the full cyclic model is physically admitted.

The A6 preserve/achieve slope and numerical benchmark are owned by PR #263. This helper set deliberately does **not** implement that slope formula again. It focuses on the distinct directional-shear/orientation calculation, bounce inequality, and Bianchi IX identities.

There are no imports from the production engine, no changed physical parameters, no source-law assumptions hidden behind defaults, and no integration/trajectory routine.

## Available calculations

| Function | Formula and intended scope |
|---|---|
| \`omega_sigma_from_directional_shear\` | \(\\Omega_\\sigma=\\sum_i\\sigma_i^2/(6H^2)\). Requires traceless shear and \(H\\ne0\). The raw value is returned without clamping. |
| \`directional_kasner_exponents\` | \(q_i=1/3+\\sigma_i/(3H)\). Makes no orientation label; changing the sign of H changes which shear pattern is dangerous. |
| \`raychaudhuri_bounce_margin\` | \(- (\\rho+p)-2\\sigma^2\). Under \(\\dot H=-(\\rho+p)/2-\\sigma^2\), positive margin means \(\\dot H>0\). The helper does not cover an effective-quantum bounce equation. |
| \`heron_discriminant\` | \(T=(a+b+c)(-a+b+c)(a-b+c)(a+b-c)\) for positive scale factors. Returns a discriminant, not a boolean wall-activity decision. |
| \`bianchi_ix_trace_term\` | \(V^2(3\\rho-\\sum_iP_i)\), the matter-trace term in the declared GEAR-47 identity. |
| \`bianchi_ix_S_tau_derivative\` | \(S'=\\tfrac12[-T+V^2(3\\rho-\\sum_iP_i)]\) under the stated diagonal Bianchi IX equations and \(dt=V\\,d\\tau\). |
| \`scalar_kasner_exponent_square_sum\` | \(\\sum_i p_i^2=1-\\phi'^2/S^2\), only as a formula helper for the declared scalar Bianchi IX reduction. |

## Guard behavior

- Non-finite numeric inputs are rejected.
- The shear fraction and directional exponents reject \(H=0\); the H-normalized formulas are undefined at a turnaround and are not regularized.
- The shear triple must be traceless to a declared numerical relative tolerance.
- Scale factors and volume must be strictly positive; \(\\sigma^2\) must be nonnegative.
- Values are not clamped into expected physical ranges.
- The Heron discriminant sign classifies the triangle inequality only. The code deliberately does not translate \(T<0\) into “wall active”: \(|T|\) and the dynamics matter separately.
- The Bianchi IX helper retains the matter-trace correction. It does not silently assume \(\\sum_iP_i=3\\rho\); that is a model condition, not a universal identity.

## Tests

The companion test module checks the exact shear normalization, the expansion/contraction orientation inversion, the signed Raychaudhuri margin, positive/negative/degenerate Heron cases, the trace-term cancellation when the total pressure equals three times the density, plus a nonzero anisotropic-pressure-sum control, the full \(S'\) identity, and the scalar/vacuum limit.

Local isolated test command:

\`\`\`bash
python -m unittest discover -s tests -p 'test_derived_physics_diagnostics.py' -v
\`\`\`

Result on the authored local copy: **8 tests passed**. This is a unit test of the declared algebra and guard behavior only. It is not a full repository run, a numerical-relativity validation, proof of finite BKL bounce count, validation of a global cycle, or scientific admission.

## Follow-up acceptance gate

The diagnostic set must remain research-only until independent review checks the GEAR normalization, Bianchi IX time gauge, and pressure signs. If the helpers are later used on simulation data, the caller must preserve provenance (input source SHA, sign/unit convention, realized shear fraction, branch sign, and excluded \(H\\approx0\) samples) and report constraint/validity failures rather than filtering them away.

**No production equations or GEAR archive files are modified by this contribution. Scientific admission remains NOT ADMITTED.**

## Coordination note

PR #263 owns the bounded A6 preserve/achieve numerical benchmark, including the free and source-aware normalized-shear slopes and the H=0 exclusion. This helper set intentionally does not duplicate that slope implementation. It retains the independent directional-shear normalization/orientation calculation, the conditional Raychaudhuri margin, and the Bianchi IX/Heron and scalar Kasner identities.
