# The production lane with the coupled floor: how far today is from φ's node (research only)

Follow-up to the floor run (README.md here), approved Nov 15. The production cosmos lane (flat, H² = ρ/3), forward from
the archive's present-day state. Script `research/production_floor_lane_v1.py` → `production_lane.txt`. The research
right-hand side reproduces `engine.cosmos.cosmos_rhs` to 1e-13 when the coupling is off
(`production_lane_validation.txt`). **Chart-dependent:** every time below is counted from the imported φ_today = 0.179.

| variant | φ's node | φ reaches it | turnaround (H = 0) | crunch | φ at turnaround |
|---|---|---|---|---|---|
| A hand floor (production, engine) | 2.303 | 3.75 Hubble times | 5.42 (1.99 e-folds from today) | 7.59 | 3.10 |
| B coupled floor, matched at φ = 0 (λ = 1.52e6) | 1.985 | 2.90 | 4.19 (1.72 e-folds) | 5.82 | 2.76 |
| B′ coupled floor, matched today (same dark energy today as A; λ = 1.57e6) | 2.010 | 2.95 | 4.26 (1.74 e-folds) | 5.91 | 2.79 |

Times are from today in units of today's Hubble time, 1/H_today = 176.4 lane units. At a Hubble time of ~14 Gyr
(an outside conversion, not derived) that is a turnaround in ~76 Gyr with the hand floor and ~60 Gyr with the coupled
floor.

## Reading

- **The turnaround is not new with the coupling: production already has it.** With the hand floor, φ reaches its node
  in 3.75 Hubble times and the lane turns around at 5.42. That follows from φ's potential having no minimum and a
  negative floor (Derived; chat's first verdict). The coupled floor brings everything forward by ~1.2 Hubble times
  (~22%), because the floor deepens as φ rolls.
- **Matching the floor at φ = 0 or today changes the times by under 2%.** The result doesn't hinge on that choice.
- **The crunch follows ~1.6–2.2 Hubble times after the turnaround.** The lane is flat and homogeneous, so it has no
  bounce: H reaches 0 only where ρ = 0, and after that ρ + p > 0 drives H downward without limit. That agrees with
  chat's derivation. The lane is integrated to a = 10⁻³ a_today, short of the singularity itself.
- In the research cell (c²-rescaled, one-budget), the same sequence ran in t ≈ 42 → 58.8. Here it is the production
  chart's version of it.

## For the next derivation (chat)

Does S reach the band before the singularity? Two notes from the lane:
- In a homogeneous contraction, any adiabatic x grows as a^-3 and so passes 4/3 at some finite a > 0. Reaching the
  band before a = 0 is therefore guaranteed *if S is still oscillating* (|H| ≪ S's mass at that moment). The real
  condition is adiabaticity at the crossing. Near the crunch |H| grows without bound, and once |H| ≳ S's mass, S freezes
  and x stops growing as a^-3.
- In the one-stone cosmos the sea's own S is tiny, and the relevant x is the stones'. The cell's core went from
  0.03 to 0.39 while |H| rose to 0.4 (S's mass is 2 there), with the band at 1.27. A factor of 3.3 more needs a to
  fall by 1.5 more, which in that cell happens within the last ~1 time unit before the singularity, where |H| ~ mass.
  **Marginal in the cell; a derivation should state the adiabaticity condition explicitly.**
