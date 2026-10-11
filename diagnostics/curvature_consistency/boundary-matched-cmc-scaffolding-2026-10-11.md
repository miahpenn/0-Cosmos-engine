# Boundary-matched CMC experiment — scaffolding — 2026-10-11

**Branch:** diag/independent-metric-ricci-audit

## Scope

This diagnostic scaffold adds an optional outer Dirichlet value to the CMC lapse solve. When unset, the default remains alpha(R)=1.0; the operator coefficients are unchanged. The wrapper forwards the optional argument, but existing production call sites do not activate the override.

The reference series contains 49 R=160 samples at index 159 (r=79.75), the cell centre that aligns with the R=80 outer cell centre. It is not the R=160 outermost cell, which is pinned to 1.0. Source: frozen-state CMC crop audit run 38117001920. The uploaded source JSON SHA-256 was 0641d5387b054a2aa97b7dfdf77a6ffed9b18ceed823bd70af9824405b9ea5ba; this repository copy preserves the sample values but is JSON-reserialized.

## Confound

geom.alpha[-1] is read by radiation characteristic speeds and outer-face flux/metric assembly and by dust outer flux/source diagnostics. Changing the CMC Dirichlet can therefore change matter-boundary inputs. A joint intervention allows those reads to follow the changed alpha. A separated CMC-only arm must freeze baseline R=80 matter-boundary alpha/metric values and assert they remain unchanged. That arm is not yet implemented.

## Admission gates

1. Diagnostic-off initialization matches the archived R=80 capture exactly across geometry, scalar, and matter arrays.
2. Default / explicit None solves are identical; default outer alpha is exactly 1.0.
3. A finite override imposes that outer value and demonstrably changes the interior lapse. Non-finite inputs fail closed.
4. The interpolation helper is validated for schema, sample ordering, finiteness, endpoints, and midpoint interpolation.
5. No full trajectory until stage-time values are threaded through all five CMC solves and the separated matter-boundary freeze assertions are implemented and passing.

Linear interpolation of the 49-slice reference (~0.25 spacing) is not a stage-level error bound for dt=0.00375. The held-out 0.5-spacing max error was ~1.5e-4 (RMS ~5.6e-5). Do not interpret the residual ridge as a physical front from this experiment alone.
