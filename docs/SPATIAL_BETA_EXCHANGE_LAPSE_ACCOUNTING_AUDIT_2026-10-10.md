# Spatial beta-coupled scalar/DM source accounting audit

**Status:** suspected source-consistency defect; equation-level verification required; no physics repair approved.
**Date:** 2026-10-10. **Branch:** research/gear-kernel-integration-gates.
**Scope:** read-only source audit. No source edits or simulation.

## Finding

The live spatial source path appears to apply the lapse to the dark-matter energy-transfer term but not to the paired scalar momentum-equation source. This mismatch vanishes in the homogeneous alpha=1 limit. Treat this as a candidate defect pending derivation from the governing action, the stress-energy divergence convention, and the ADM definition/sign of Pi. Do not patch from this report alone.

## Source evidence

- [Scalar RHS](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-kernel-integration-gates/engine/scalar_system.py)
- [Conservative fluid RHS](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-kernel-integration-gates/engine/matter_system.py)
- [Covector source and coordinate/normal derivative](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-kernel-integration-gates/engine/valencia.py)
- [Live source assembly](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-kernel-integration-gates/engine/production_kernel.py)

In scalar_rhs_arrays(), the code defines the field derivative as

    partial_t(phi) = alpha * Pi + beta^r * partial_r(phi)

and first constructs the Pi RHS with lapse-weighted wave, curvature, and potential terms, including -alpha * Vp. It then adds the DM coupling separately:

    if rho_dm is not None:
        Pit = Pit + beta_dm * rho_dm

This is an additive coordinate-time source +beta_dm*rho_dm without an explicit lapse.

The live production RHS passes that scalar time derivative and the same rho_dm slice into the conservative DM species RHS. It forms

    Q_t = beta_dm * rho_dm * partial_t(phi)
    Q_r = beta_dm * rho_dm * partial_r(phi)

and adds to the conservative DM energy equation

    q_E = -sqrt(gamma) * (Q_t - beta^r * Q_r).

Substituting the implemented partial_t(phi) expression cancels the shift terms, leaving

    q_E = -sqrt(gamma) * alpha * beta_dm * rho_dm * Pi.

## Conditional source-balance implication

For the conventional normal-derivative Pi and the paired covariant interaction suggested by these sources, the scalar energy-transfer contribution should be proportional to +sqrt(gamma)*alpha*beta_dm*rho_dm*Pi and cancel the DM contribution. The current scalar source appears instead to contribute +sqrt(gamma)*beta_dm*rho_dm*Pi. Under those source/sign conventions, the pair's coordinate-time residual would be

    Delta_exchange = sqrt(gamma) * (1 - alpha) * beta_dm * rho_dm * Pi.

This is an algebraic implication of the source expressions, not a measured production residual. Before accepting it, an independent derivation must confirm the action, field equation, metric signature, Pi convention, source index/sign, and density normalization. A nonzero remainder might only be called a bug once that derivation fixes the expected identity.

The homogeneous lane has alpha=1 by construction, so the homogeneous scalar/DM continuity test cannot probe this lapse-dependent mismatch.

## Test-coverage limitation

tests/test_exchange.py tests the separate v55_matter.exchange_pair() helper. The live spatial path computes its scalar and DM source contributions separately in scalar_system.py and matter_system.py. The helper's opposite-sign arrays therefore do not prove cancellation in the actual production RHS. Existing center-regularity and outer-constraint tests passed in the normal suite, but they do not test exchange-source conservation.

## Required next verification

1. Derive the scalar and DM sources from the governing action / archived covariant equations, fixing metric signature, Pi, units, index placement, and Q_nu sign.
2. Derive the coordinate-time scalar energy balance, explicitly tracking alpha, sqrt(gamma), shift, and the same DM density primitive.
3. Build a deterministic, non-evolving source-pair check at alpha != 1 and nonzero shift, using the actual RHS source assembly. It must show cancellation, or explain from the derived equations why a remainder is physical. Do not test only the disconnected helper.
4. If derivation confirms the lapse factor is missing, make the smallest correction on a separate physics branch, keep prior commits/artifacts intact, and run deterministic gates before any trajectory.
5. Do not launch more expensive runs or attempt cross-lane coupling until the source accounting is resolved.

## Explicitly not done

- No alpha factor, sign, coefficient, or source term changed.
- No damping, fit, floor, clamp, or radiation boundary change.
- No simulation or long campaign launched.
- No archive or REV16 criterion changed; no merge to main.

**Provisional conclusion:** sharing the same beta parameter does not prove a conserved pair. Resolve the lapse weighting first, before admitting spatial trajectories as physically valid.
