# D_r's return paths on an expanding background (research only)

Queued item from HANDOFF Nov 15 (1), and GEAR-105's open "same analysis on an expanding background".
Script: `research/dr_return_toy_v1.py` (homogeneous, linear, scaled units, a few seconds). Output: `toy.txt`.
H is in units of mu (S's mass scale). Production today: H/mu = 0.0057.

## Path 1: S ringing in the sea makes D_r (the sea's answer), now with redshift

- **Derived.** A D wave redshifts through the S^2 resonance once. In the wave's number (its adiabatic invariant)
  the crossing multiplies it by exp(pi mu^2 / H), mu = (3/32) x_A at the crossing.
- **Measured.** The toy (constant H, S diluting as a^-3, seed phase averaged) matches this to 3–30%, tightening as H falls:
  at x_A 0.33, H = 0.001: toy ln gain 2.74, predicted 2.97; at x_A 0.13, H = 0.001: 0.449 vs 0.439.
- **Thresholds (Derived).** An order-one gain needs x_A > 6.0 sqrt(H/mu). Draining most of S's ringing from a seed of
  relative size sigma needs pi mu^2/H > ln(1/sigma^2), i.e. x_A > 6.0 sqrt(H/mu) sqrt(ln 1/sigma^2).

  | H/mu | x_A for gain e | x_A to drain from sigma = 1e-10 |
  |---|---|---|
  | 0.1 | 1.9 | 12.9 |
  | 0.0057 (production today) | 0.45 | 3.07 |
  | 1e-10 | 6e-5 | 4e-4 |

- **Read.** At production's H/mu, the whole hold band (x_A 1.27–2.67) lies below the drain line for a 1e-10 seed: S's
  ringing is pumped into D_r by up to e^35 in number, but whether that finishes the job depends on the seed. The
  result is **seed-limited**, and the seed (numerical noise, initial-data fluctuations, or anything quantum) is an
  **unpaid input**. At the very small H/mu of a physical universe the drain is essentially unconditional.

## Path 2: D_r on its own (GEAR-105's long-wave breakup)

- A uniform D ring of relative amplitude eps breaks into long D waves at ~1.1 eps^2 (GEAR-105), against dilution
  eps ~ a^-3/2. Accumulated exponent 1.1 eps0^2/(3H). At production H it needs eps0 ~ 0.84 to reach 46: **D_r made
  by S does not break up; it dilutes like matter** (mass sqrt2, w = 0). (Derived, using GEAR-105's rate.)
- GEAR-105 already showed D_r cannot decay back into S outside matter. So once made, D_r's only returns are
  expansion (dilution) and matter: being reabsorbed where S still rings.

## Path 3: back through matter (distance)

- **Derived.** D_r made by S moves at k/omega = sqrt2/2 = 0.71. It is pumped back inside a pocket of radius R only
  if mu R / 0.71 > 1, i.e. x_A R > 7.5 (R in 1/mu). Toy pockets (R ~ 8) need x_A > 0.95 — the hold band sits above
  this, so inside held pockets the exchange runs both ways; in the thin sea it is one way, outward.

## What this means for the lane's S–D species (Read; Open)

- D_r from the sea's answer is a **cold, matter-like** component (w = 0, mass sqrt2) that does not return to S in the sea.
  The lane's S–D species therefore has a third piece besides held and released: drained D_r, counted as matter.
  Its size is seed-limited at production H/mu and needs a run with a stated seed to fix. No lane change proposed.
- Nothing here changes the two-component lane's fit (its residual is already matter-like); it changes the label of
  part of "released".

## Correction to SEA_ANSWER (see AMENDMENT-1 below)

Unpaid inputs: the seed sigma; constant H in the toy (a slowly varying H changes the crossing rate locally, not the form);
linear waves (no backreaction on S); homogeneous S (pockets only through path 3's distance estimate).
