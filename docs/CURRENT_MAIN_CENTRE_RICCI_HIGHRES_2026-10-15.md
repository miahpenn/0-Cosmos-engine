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


## Primary discriminator: offline chain-rule closure for the BSSN connection `Lambda`

A newer independent offline decomposition of the original high-resolution artifact (team handoff: [Issue #5 comment 6185063664](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6185063664)) suggests the strongest current lead is the discrete derivative of the composite metric-to-`Lambda` relation, not the origin half-cell quadrature alone.

The pinned initializer sets
[
\Lambda = \frac{D_1(a)}{2a^2}-\frac{D_1(b)}{ab}+\frac{2}{r}\left(\frac1b-\frac1a\right).
]
The vendor BSSN Ricci calculation uses a separate (D_4(\Lambda)) term because (\Lambda) is an independently evolved connection variable. For a controlled **offline-only counterfactual**, this observer computes the ordinary chain-rule derivative of the displayed metric definition using the same sampled (a',b',a'',b'') stencil operands:
[
\Lambda'_{\rm chain}=\frac{a''}{2a^2}-\frac{(a')^2}{a^3}-\frac{b''}{ab}
+\frac{b'(a'b+ab')}{a^2b^2}
-\frac{2(1/b-1/a)}{r^2}
+\frac{2(-b'/b^2+a'/a^2)}{r}.
]
The sole Ricci-term substitution in this counterfactual is (D_4(\Lambda)\mapsto\Lambda'_{\rm chain}); the induced radial trace-free Ricci shift is (	frac23 X^2(\Lambda'_{\rm chain}-D_4(\Lambda))). The new report records both derivatives, their signed difference, the induced shift, and baseline-versus-counterfactual gaps for cells 0/1, the full grid and fixed (5\le r\le10) annulus, with all original source/vendor/state/checksum gates retained.

The team's read-only calculation on the preserved artifact reported a full-grid gap reduction from approximately (1.59\times10^{-3}) at cell 0 to (2.45\times10^{-10}), (5.06\times10^{-11}), and (1.63\times10^{-11}) at N=160/240/320. The refreshed observer must reproduce the decomposition at its exact head before relying on those prior values.

**Critical limit:** this is strong evidence of a discrete product/chain-rule closure mismatch between two representations on the initialized slice. It does not authorize replacing the evolved BSSN connection derivative in production. The vendor explicitly relies on independent `Lambda` as part of the BSSN principal structure, and direct metric Ricci need not equal the evolved-connection expression away from the connection-constraint surface. No vendor or production change is part of this PR.

## Added frozen-source origin half-cell counterfactual

The pinned production initializer computes its first cumulative integral entry as `0.5*dr*source[0]` at the first cell centre (r_0=\Delta r/2). For a smooth even source
[
s(r)=s_0+s_2r^2+s_4r^4+O(r^6),
]
the regular half-cell integral is
[
I_0=s_0r_0+\frac{s_2r_0^3}{3}+\frac{s_4r_0^5}{5}+O(r_0^7),
]
but the point-sampled implementation yields (r_0s(r_0)=s_0r_0+s_2r_0^3+s_4r_0^5+O(r_0^7)). The quadratic moment is therefore overweighted by a factor of three. Because (X=B^{1/6}), the leading first-cell coefficient in (X) inherits the same factor-three tendency.

The expanded observer reconstructs the final source of the pinned six-iteration initializer and verifies the reconstructed radial integral (B) against the actual geometry value (X^6). It fits an even (s_0+s_2r^2+s_4r^4+s_6r^6) model over fixed physical radii r≤2.5 and r≤3.0, with the difference reported as a fit-width sensitivity check. It then corrects only the origin half-cell integral and propagates its (delta I/r) effect to a **copied metric**. All inter-centre trapezoids remain unchanged, the source is frozen rather than iterated again, and no production state or source file is edited. The report records baseline and counterfactual signed Ricci gaps for the first two cells and fixed annulus.

The origin-half-cell mechanism remains a secondary hypothesis. This is deliberately a local mechanism discriminator, not a self-consistent initializer candidate. Earlier origin-half-cell experiments improved initial constraints but did not restore evolved centre convergence; their source blobs differ from this current source anchor. A reduced frozen-source gap therefore supports the quadrature hypothesis but does not authorize a production patch, prove which curvature representation is wrong, or establish scientific admission.

## Acceptance and limits

A useful result requires all source/vendor pins to match, finite profiles, unchanged state digests, zero evolution steps, and a reproducible checksum manifest. Passing those integrity gates means the diagnostic completed, not that either formula is correct. No production equation, vendor code, coefficient, gauge, boundary, threshold, dependency, or campaign setting is changed. Scientific admission remains **NOT ADMITTED**.
