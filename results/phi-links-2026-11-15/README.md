# φ's links to S and D, from the machine's own coupling, against the unpaid items (research only)

Chat (Nov 15), reading 0*: everything the S–D sector can't supply alone (floor, seed, second reset, κ*/λ*, dust)
needs a moving player outside S and D, and φ is the half-wired candidate. Its links to S and D are missing in production
(interaction matrix, Nov 13). Task: build them from the machine's own coupling, without forcing, and see whether the
loop sea → S → D → sea closes.

The machine's own coupling is GEAR-91 (research switch, PR #465):
V_SD = e^(−2β_S φ)S²/2 − e^(−2β_D φ)D²/2 + (λ/8)(S⁴ + 6S²D² + D⁴), and φ's equation gains 8π ∂V_SD/∂φ. Here
β_S = 0 (Cassini) and β_D = −0.04 (locked with dark matter, borrowed). No runs; algebra checked numerically on the
switch's own `potential()`. Script `research/phi_links_v1.py` → `phi_links.txt`.

## Scorecard (Derived unless marked)

| unpaid item | does φ pay it through GEAR-91? | why |
|---|---|---|
| **floor** | **yes, this is the floor mechanism** | In the sea, D's vacuum depends on φ: V_vac(φ) = −e^(−4β_D φ)/(2λ). φ feels it through the 8π ∂V_SD/∂φ term, so D's step *is* φ's floor, with no hand floor. It matches the hand floor at φ = 0 for λ = 4π/(D_POT V0) = 1.52e6. |
| φ parks at its node? | **no** | V_eff = V0 e^(−φ) − (4π/λ) e^(0.16φ) falls everywhere (both slopes negative). It has a node, no minimum. The node moves from ln 10 = 2.303 (hand floor) to 2.303/1.16 = **1.985** (D's floor). The floor's slope is 16% of the exponential's at the node. |
| **seed** | **no** | V_SD is even in D. ∂V/∂D = 0 at D = 0 for every S and φ (checked: exactly 0). Nothing φ does can start D off the saddle. The seed has to be in D already. |
| **second reset** | **no** | D's mass² at the saddle is −e^(−2β_D φ) < 0 for every φ. φ can never make the saddle stable for D. Only S's energy can, via (3λ/2)S² > e^(−2β_D φ). With β_S = 0 φ has no channel into S. With β_S ≠ 0 it changes S's mass only at its own rate (~H, far below S's mass), so adiabatically, and makes no S. |
| **κ*/λ*** | **partly: φ sets κ, not λ** | GEAR-03's split κ = (K_S − K_D)/2 = (1 − e^(−2β_D φ))/2 ≈ β_D φ: −0.02 at φ = 0.5, −0.04 at 1, −0.086 at the node. This is C6's "hierarchy relocated to φ near its node", now from the coupling. λ is still unpaid. |
| **dust** | **no** | φ → dark matter exists (β_dm), but it moves imported dust; it doesn't make it. |

**So the loop does not close through φ in GEAR-91's form.** sea → φ closes (the floor). φ → D is adiabatic (D_s and
the step deepen as φ rolls: +38% in the step by φ = 2). D → sea is the floor again. But **φ → S, and φ → back to the
saddle, are forbidden by the coupling's form**: it is even in D, and the exponential never changes the sign of D's
mass². Paying the seed or a second reset with φ would need a coupling that isn't in the gears. That would be forcing.

## What this means

- **φ pays the floor.** That is the floor question we had queued next, and it now arrives with a derived
  consequence: φ's node moves to 1.985 and the step grows as φ rolls. The test is a short cell or lane run with the
  GEAR-91 switch on, the hand floor off, λ at the matching value (or the toy's 5000 with V0 rescaled), checking that
  geometry counts the floor once (the double count of AMENDMENT-REV20-2) and that φ rolls past the moving node.
- **The seed and the second reset need another player, or are features rather than gaps.** Within the machine's own
  couplings, φ can't supply them. The D_s sea elsewhere (chat's other candidate) is the same field. Its D is not at
  the saddle, so it can only seed through gradients at a boundary, and the cell already shows that: the release starts
  from D's tail out of the surrounding sea (handoff note). That is a candidate for the seed, not for a second reset.

## Two stones (scope)

The cell with its reflecting wall is already a lattice of identical stones: each cell is one Wigner–Seitz cell of a
population. So the population question (STONE_SPLIT) is partly answered by the cell as it stands. Two *unequal* stones
(different sizes or times) need two or three dimensions, and the radial code can't do that. A centre-plus-shell
configuration is possible in spherical symmetry, but it is a shell, not a stone. **I'd hold the widening: it doesn't help
any item above, and the φ floor run can be done in the existing cell.**

Unpaid inputs: β_D = −0.04 (borrowed), β_S = 0, λ, V0.
