# Issue #4 — radiation recovery round-trip sweep findings

Date: 2026-10-19
Status: diagnostic-only; current-identity cause not yet isolated.

## Provenance

- PR: [#408](https://github.com/miahpenn/0-Cosmos-engine/pull/408), remains draft.
- Workflow: [run #39190152670](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/39190152670), completed successfully.
- Artifact: `accepted-state-raw-momentum-witness`.
- Historical production source: `7f47669a1c0edfafef9fc5f27b2db80840d3f035`.
- Vendor pin: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Initial digest: `a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0`.
- Accepted digest: `c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a`.
- Configuration: resolution 60, `r_max=80`, `dt=0.015`, one production step; recovery tolerances `1e-8`, `1e-10`, `1e-12`, `1e-13`, `1e-14`.

## Results

Maximum relative round-trip errors across 60 cells:

| State | Tolerance | Energy error | Momentum error | Exceptions |
|---|---:|---:|---:|---:|
| Accepted | 1e-8 | 7.4506e-9 | 1.2668e-13 | 0 |
| Accepted | 1e-10 | 9.9186e-11 | 1.2903e-15 | 0 |
| Accepted | 1e-12 | 9.9055e-13 | 1.3432e-17 | 0 |
| Accepted | 1e-13 | 9.9633e-14 | 1.6000e-18 | 0 |
| Accepted | 1e-14 | 9.5155e-15 | 1.2811e-19 | 0 |
| Manufactured | 1e-8 | 9.6942e-9 | 1.1657e-10 | 0 |
| Manufactured | 1e-10 | 9.9654e-11 | 1.1977e-12 | 0 |
| Manufactured | 1e-12 | 9.9635e-13 | 1.1976e-14 | 0 |
| Manufactured | 1e-13 | 9.4694e-14 | 1.1372e-15 | 0 |
| Manufactured | 1e-14 | 9.7173e-15 | 1.1953e-16 | 0 |

Every tested output was finite, and there were no recovery exceptions. The recovery round-trip errors decrease substantially as tolerance tightens.

## What this establishes — and does not establish

The data establish a tolerance-dependent recovery round-trip error trend. They do not establish that recovery alone causes the radiation projected-current identity mismatch. The sweep does not rerun the projected-versus-direct identity at each tolerance and does not record per-cell pressure residual or solver exit class. Therefore no recovery-aware identity bound is justified yet, and the strict roundoff-only identity gate remains unchanged.

## Next bounded test

For each of the five tolerances, rerun the current identity on both accepted and manufactured states. Preserve per-cell signed `E`, `S_r`, `S^2`, recovered pressure `p`, `v^2`, residual `f(p)`, iteration count/exit class, measured and predicted signed `j_projected - j_direct`, and arithmetic remainder. Include explicit zero-velocity handling and a wrong-sign negative control. Fail closed on nonfinite values, failed recovery, missing cells/fields, digest/pin mismatch, or a failed negative control. Define any residual-to-current error enclosure before evaluating the identity results; do not loosen the original roundoff gate post hoc.

Only after the recovery/current discrepancy is independently understood should work proceed to the separate independently assembled matter + geometry + scalar momentum-propagation closure.

## Guardrails

PR #408 remains draft. PR #406 remains untouched. Production equations, source/vendor pins, expected state digests, and scientific-admission rules are unchanged. Keep `independent_full_propagation_identity_closed=false`, `scientific_admission=false`, and `acceptance_oracle_satisfied=false`.
