# Test of the bound-pocket virial in an expanding background (research only)

Derivation: chat, Drive "DERIVATION-bound-pocket-virial-in-an-expanding-background … Nov15-2026". Stored states only
(x_c 2.0 at r_max 16, 24, 36; x_c 1.5). Script `research/core_virial_v1.py`. Outputs: `core_virial.txt` (form 1),
`core_virial_form2.txt` (form 2, superseded), `core_virial_form3.txt` (form 3, the test), `core_flux.txt` (superseded:
its flux formula assumed a diluting core).

## Summary

- **The identity holds once the boundary flux it neglects is small (Read).** In the central observer's frame,
  ΣT^ii = ½Ï − (Ḣ+H²)I − W_x + S_b closes to residual/T^ii ≈ 10–25%, i.e. |resid|/(3E) = 0.002–0.009, for core
  radii r_b = 8–10. At those radii the energy flux through the boundary is ≤ 0.02 HE. At r_b = 4, energy flows into the
  core at 0.7–0.9 HE and the identity misses by 0.05–0.16. At r_b = 6 (flux 0.12 HE) it misses by 0.010–0.015. **The
  residual tracks the neglected flux, as the derivation says it should.**
- **It does not yet predict the core's w.** The prediction is the small difference of large terms. ½Ï and −(Ḣ+H²)I
  are each 0.1–0.2 per 3E and cancel to 0.01–0.03. The residual (0.002–0.009) is as large as the w being predicted
  (+0.004 to +0.009). **Open.**
- **Read from the terms:** the steadiest piece is the self-gravity support, −W_x/(3E) = +0.012 to +0.019 (x_c 2.0),
  +0.008 to +0.013 (x_c 1.5): larger than the core's w itself. The boundary term S_b is small and **positive**
  (+0.0002 to +0.004 at r_b 6–10), so it does not look like a held, w = −1 piece at the core's edge. **Jeremiah's
  guess that the leftover held step is the sea-pressure term is not supported here** (Read, at this accuracy).

## Two refinements to the derivation (for chat)

1. **The surface term is the radial stress, not the isotropic pressure.** For a spherical boundary,
   ∮ x_i T^ij n_j dS = ℓ_b T^rr(ℓ_b) × area. For a field, T^rr = K + G − V (gradient counted fully) differs from
   p = K − G/3 − V. Used here.
2. **Frames.** The cell's slice observers move with the Hubble flow, so the stored stresses are slice-frame. The
   identity is in the central observer's frame. To first order in the boost:
   T^ii_F = T^ii_slice + 2H Σ j ℓ dV + H² Σ (e+p) ℓ² dV, with j the field's radial momentum density. The sign of j is
   fixed by the core's energy balance (correlation 0.97–1.00 at r_b 4–6). Without these terms (forms 1–2) the test
   fails by 0.02–0.2.

**On the expected scaling:** the derivation expects errors to shrink as (R_b H)² from r_max 16 to 36. They can't: the
core's size doesn't depend on r_max. The three x_c 2.0 cells give the same terms to ~0.001. Here ℓ_b H ≈ 0.4–1.2,
which is not small. The first-order boost is the weak point at the larger radii (boost terms 0.02–0.2 at r_b = 10).

## Results (each term / 3E, t = 13–54 mean; full tables in `core_virial_form3.txt`)

| cell | r_b | flux in / HE | slice w | T^ii_F | ½Ï | −(Ḣ+H²)I | −W_x | S_b | residual |
|---|---|---|---|---|---|---|---|---|---|
| x_c 2.0, r_max 24 | 6 | 0.12 | +0.008 | +0.017 | −0.059 | +0.069 | +0.016 | +0.004 | −0.012 |
| | 8 | 0.02 | +0.007 | +0.031 | −0.065 | +0.087 | +0.017 | +0.000 | −0.009 |
| | 10 | 0.00 | +0.007 | +0.040 | −0.070 | +0.099 | +0.017 | +0.000 | −0.007 |
| x_c 1.5, r_max 24 | 6 | 0.13 | +0.004 | +0.036 | −0.059 | +0.092 | +0.011 | +0.004 | −0.012 |
| | 8 | 0.02 | +0.004 | +0.053 | −0.069 | +0.116 | +0.012 | +0.000 | −0.005 |
| | 10 | 0.00 | +0.004 | +0.066 | −0.077 | +0.132 | +0.012 | +0.000 | −0.002 |

(r_b = 10 rows with j's sign fixed at +1 from r_b 4–6; there the energy balance can't fix it, because the flux is ~0.)

## What this means for the contract

- The core's pressure is no longer unexplained in kind: it is the support against its own gravity plus the tidal and
  relaxation terms, with the expansion handled. It is **not yet predicted in size**. The contract keeps the w(a) fit
  for the S–D species (Amendment 2), now with a named physical content.
- To get a sharp prediction: either (a) keep the boundary flux terms exactly (the Leibniz terms of a comoving region:
  energy and momentum through ℓ_b), so that r_b = 6 can be used where the boost is small; or (b) evaluate Ï and the
  tidal term together analytically for a bound core (I nearly constant in proper terms), so that the large
  cancellation is done on paper. (b) is a chat derivation; (a) is a short script change.

Unpaid inputs: smoothing (two boxcars of π); the core's mean lapse as the central clock; the first-order boost.
