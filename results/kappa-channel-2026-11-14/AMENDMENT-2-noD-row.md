# Amendment 2 to the channel note (14 Nov 2026): the no-D row; decisions adopted

## Correction

The B3 table's row "no D, r_max 24" is wrong. `research/cosmos_channel_v1.py` left the restored hand floor
(D_POT = 0.10 in the no-D control) out of φ's potential. The missing constant went into the fluid remainder and
inflated F for that row only. Found by B6's fully independent build on a fresh slice. Script corrected;
`channel_frozen_states.txt` regenerated. No other row changes (they have no hand floor).

| no D, r_max 24, t = 15 … 63 | was | corrected |
|---|---|---|
| C/F | +0.003 … +0.0002 | +0.035 … +0.032 |
| S sector | +0.003 … +0.0002 | +0.033 … +0.029 |
| shear | 1e-5 | 2e-4 … 1.3e-3 |
| spread of Π | 2e-5 … 5e-6 | 4e-4 … 9e-4 |
| Jensen | −8e-6 | −1e-4 … −6e-5 |
| lapse on cosmos terms | 9e-5 … 2e-6 | 1e-3 … 3e-4 |

Consequence: B6's third control ("no-D pocket: C falls to 0.003 F and below") was a wrong prediction. A pocket of S
without D carries about 40% of the full channel. The note's statement that the leading part is "the S–D pocket's own
ρ + 3p" stands; S alone already supplies much of it.

## Decisions adopted (Jeremiah, Nov 14)

1. The cosmos lane's t is the proper time of its own fluids, rest-mass weighted (derived from its conservation laws).
2. The pocket's inner redshift is its own labelled quantity, sized only on a slicing-independent clock.
3. B6 runs on the dust clock: `results/b6-dust-clock-2026-11-14/README.md`.

Flag: the dust is imported, so t rests on an imported species until S–D matter enters the cosmos lane.
