# High-resolution centre Ricci replay — 15 October 2026

## Question

Do the existing first-cell vendor-minus-warped-product trace-free Ricci differences at N=60/80/120 enter a stable refinement regime when the same production initializer is sampled more finely?

This is a bounded no-evolution diagnostic. It is not a production physics correction or a scientific-admission claim.

## Pins and configuration

- PR base / live main snapshot at branch creation: `3971d9a7e959e4128c030319ee7305d390df105b`.
- Production-source comparison anchor: `285a2ee0c1a3b97cddeb4519d8e749109b115aab`.
- Vendor submodule: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Resolutions: N=160, 240, 320; outer radius r_max=80.
- The audit validates the 11 expected engine blob hashes, vendor pin, and empty engine/vendor/dependency diff before it evaluates any profiles.
- The observer reports state digests and zero evolution steps. Its result bundle contains signed CSV/NPZ profiles, strict JSON, source-pin metadata, and checksums.

## Why this is distinct

The merged #218 artifact already records N=60/80/120 on the same declared engine-source anchor and vendor. PR #242 tests a synthetic c5 profile and an exactly predicted direct-metric stencil truncation term; it does not sample the initialized production geometry. This run extends only the initialized production-profile resolutions and writes to a fresh output path, preserving prior artifacts.

## Reported comparisons

For each N, inspect the signed vendor and direct warped-product TF-Ricci values for cells i=0 and i=1; their difference and difference divided by r²; the lapse-weighted Ricci contribution to the A_a RHS; and the maximum signed/absolute gap in the fixed coordinate annulus 5≤r≤10. Keep these center and annulus results separate and do not assume a convergence order ahead of the measurements.

## Acceptance and limits

A useful result requires all source/vendor pins to match, finite profiles, unchanged state digests, zero evolution steps, and a reproducible checksum manifest. Passing those integrity gates means the diagnostic completed, not that either formula is correct. No production equation, vendor code, coefficient, gauge, boundary, threshold, dependency, or campaign setting is changed. Scientific admission remains **NOT ADMITTED**.
