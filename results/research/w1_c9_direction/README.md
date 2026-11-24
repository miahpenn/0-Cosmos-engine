# W1 input tails: direction and reversed-tail fault test (diagnostic only)

Stored native.npz (sha256 16151c44...a98ea) only; no simulation; no gate changed. #577's input-check result stands
on the gates as declared in AMENDMENT 2. A direction check 'e is a future proposal (Lumen², Amendment 3).
Script: research/w1_c9_direction/run.py. Full record: report.json.

**Direction.** Signed flux F(r) = -4 pi r^2 f_t f_r, with f_t = P - H0 r f_r (the chart term removed); positive means
outward. Window means (D r in [20, 60-5dr], S r in [35, 60-5dr]):

| | D mean (share > 0) | S mean (share > 0) |
|---|---|---|
| dr 0.025 | +4.275e-10 (99.8%) | +5.432e-13 (99.8%) |
| dr 0.05 | +4.221e-10 (99.7%) | +5.271e-13 (99.8%) |

Both tails are outward on both grids. The D mean sits within 1.6% / 0.3% of the recorded P_D 4.207e-10.

**Reversed-tail fault** (rest momentum sign flipped for r beyond ~12, producer's chart term kept):

| | caught by 'a | 'b | 'c | 'd |
|---|---|---|---|---|
| S reversed, both grids | no | no | no | no |
| D reversed, both grids | no | no | no | **yes** (last-5 PD at 1.62x / 1.77x its bound) |

'a-'c square the chart-corrected momentum, so they cannot see direction. 'd catches a reversed D only because it
bounds |PD| with the chart term in, where the Doppler factor 1 - H0 r k / Omega becomes 1 + H0 r k / Omega. 'd reads
no S field, so a reversed S tail passes all four (Aster's finding, confirmed).
