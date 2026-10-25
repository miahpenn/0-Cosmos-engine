# Issue #4 — accepted-state momentum witness inspection

Date: 2026-10-19  
Status: diagnostic evidence verified; independent propagation closure remains open.

## Provenance and replay

- PR: [#408](https://github.com/miahpenn/0-Cosmos-engine/pull/408), draft and diagnostic-only.
- Replay commit: `5d982a2b2e9991fb3534c64941745fc6ee3d8b80`.
- Diagnostic base/main recorded by PR: `4aec88d0e3ff028e0a65a8ac36b30948625df792`.
- Historical production source pin: `7f47669a1c0edfafef9fc5f27b2db80840d3f035`.
- Vendor pin: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Configuration: resolution 60, (r_{max}=80), one production step, (Delta t=0.015).
- Expected initial-state SHA-256: `a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0`.
- Expected accepted-state SHA-256: `c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a`.

The successful witness workflow records exact state identity, checks that observation leaves the accepted state unchanged, emits raw arrays and a SHA-256 manifest. The report explicitly marks `independent_full_propagation_identity_closed=false`, `scientific_admission=false`, and `acceptance_oracle_satisfied=false`.

## Signed centre witness (reported from artifact inspection)

| Quantity | Cell 0 (`r=2/3`) | Cell 1 (`r=2`) |
|---|---:|---:|
| Total momentum constraint (M_r) | +5.9050e-5 | -1.1199e-5 |
| Geometric contribution (M_{geom}) | +5.7908e-5 | -1.4145e-5 |
| Signed tangent (dot M_r) | +3.9337e-3 | -7.4620e-4 |
| Geometric tangent contribution | +3.8576e-3 | -9.4251e-4 |
| S-sector tangent contribution | +7.6060e-5 | +1.9632e-4 |

The S-sector term offsets part of the geometric tangent at both centre cells. This is a signed decomposition only; it does not establish causation or identify a faulty sector. All current-sector terms must be retained in the independent identity.

The two finite-difference directions, (epsilon=2e-5) and (1e-5), agree closely at the centre (reported differences about (2.22e-14) at cell 0 and (1.16e-13) at cell 1). The reported internal tangent-decomposition residual is at most (3.47e-17) for the first perturbation. These are consistency checks of the existing tangent path, not independent propagation closure.

For the sampled fixed-coordinate (r=5)–(10) annulus, the reported maximum absolute total tangent is approximately (1.696e-6) at (r=7.333); the geometric tangent reaches approximately (1.521e-4) at (r=6), with substantial signed cancellation. This annulus is coordinate-fixed, not a fixed-areal-radius worldtube.

## Interpretation and limits

1. State identity, pin checks, finite-array checks, and artifact hashes passed for the recorded run.
2. The accepted state was not mutated by observation.
3. The finite-difference tangent is internally consistent, but is not an independent continuum/semi-discrete propagation oracle.
4. Do not infer causation from term magnitudes or tune coefficients to the observed residual.
5. No production equations, source/vendor pins, or admission rules are changed by this record.

## Next gate: independent signed closure

Build the propagation RHS independently on this exact accepted state and its recorded rates. Preserve index placement, signs, parity-aware discrete derivatives, lapse/shift conventions, total stress-energy divergence, and every matter/scalar exchange term. Compare predicted and measured signed (dot M_r) at cells 0 and 1 and separately over fixed-coordinate annuli. Include manufactured controls and preregistered numerical tolerances; fail closed on missing terms, non-finite arrays, provenance mismatch, or unsuccessful signed closure.

Only after this independent closure succeeds should the team decide whether further residual attribution is warranted. PR #406 remains untouched. Both scientific-admission gates remain false.
