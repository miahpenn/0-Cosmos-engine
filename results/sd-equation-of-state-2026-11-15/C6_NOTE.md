# C6 candidate: the reset sets κ*/λ*. Derivation note with 0* verdicts

HANDOFF Nov 15 2026 (1), Task 3. RESEARCH ONLY. Algebra only; no runs. `research/c6_reset_v1.py` → `c6_reset.txt`.
Sources read: AMENDMENT-CORE-SHEET-C6 (Nov 15), GEAR-104, AMENDMENT-BRIEF-FOR-FABLE (A2: the clock-stop band).

Marks: **Derived**, **Read**, **Picture** (Jeremiah's, labelled), **Hypothesis** (the archive's clock-stop reset),
**Open**. Scaled units: μ = S's mass at the saddle = 1, energies × 1/λ.

## Summary

- **The closure count holds.** C6 plus the one-budget rule are two relations for two unknowns (κ/λ and H/μ). C6 alone
  fixes κ*/λ*, without circularity.
- **κ\*/λ\* ≈ 0.02–0.5, λ\* ≈ 50–1 200** (clock-stop variants). The literal balance H = ω_S gives λ* ≈ 8–12.
  Ruler-dependent through how "balance" is defined, by up to a factor 100. Order fixed: at the reset H is within an
  order of magnitude of matter's clock, by construction of the criterion.
- **It does not pay for the hierarchy.** Under the one-budget rule this λ* makes φ's potential scale V0 ≈ 0.1–2.3 μ².
  Today's H ≪ μ then needs φ's potential and D's step to nearly cancel: φ within δ ~ ε² of the potential's node. The
  unpaid number moves from a ratio of scales to a near-cancellation: the cosmological-constant problem again.
- **The region's size is the whole question (3d).** The uniform version has no size because the region is a Hubble
  patch; that is why no hierarchy appears. A finite pocket would need its size in matter's units, which nothing in the
  machine selects yet.

## 3a. The energy at the saddle (Derived)

A uniform region of S–D matter at the saddle, D_s held out, S ringing with turning point x_A:

    ρ̃ = B + U + G + ½          (bond, each element's quartic, gradients, the held vacuum step ΔV = ½)
      = x_A/2 + x_A²/8 + ½       (uniform: G = 0; S's ring gives x_A/2 + x_A²/8)

- Every term is λ-free in scaled units.
- **Gravity's share** is the only place κ*/λ* enters: H² = (κ/λ) s ρ̃ / 3, s = S–D matter's share of the total density
  at the reset.
- "The uncoupled stiffness, felt as expansion": in the identity of the EOS note the quartic's ringing averages to
  radiation, and radiation's pressure is what makes H fall faster than dust's. That much is derived. That it *is* the
  breathing is the picture.

## 3b. The balance, and κ*/λ* (Derived, given the criterion)

"Breathing H against matter's clock at the saddle" needs a rule for when they balance. Three variants:

- **A.** S–D matter dominates at the reset (s = 1). The archive's clock-stop rule applied with the region's own w (the
  held w, −0.32 to −0.07 across the band), not w = 1.
- **B.** A stiff component dominates (the archive's w = 1 contraction). The archive's band 0.139–0.373 unchanged.
  S–D matter's share s is unknown, so the result is κ*/λ* · s.
- **C.** The literal balance H = ω_S (S's clock at the saddle, 1.21–1.41).

The clock-stop rule generalised (**Derived**): in a contraction with H = −p/τ, p = 2/(3(1 + w)), S rings as
τ^(−ν)[c₁J_ν(τ) + c₂Y_ν(τ)], ν = (3p − 1)/2. The archive's band is p/j_ν1 to p/y_ν1. At w = 1 that reproduces
0.139–0.373. At w = 0: 0.21–0.42. At w = −0.16: 0.23–0.44.

| x_A | ρ̃ | A: h | A: κ*/λ* | A: λ* | B: κ*/λ*·s | B: λ*/s | C: κ*/λ* | C: λ* |
|---|---|---|---|---|---|---|---|---|
| 1.28 (band edge) | 1.35 | 0.26–0.45 | 0.15–0.46 | 55–170 | 0.043–0.31 | 80–590 | 3.3 | 7.6 |
| 2.00 | 2.00 | 0.23–0.44 | 0.08–0.29 | 87–310 | 0.029–0.21 | 120–870 | 2.6 | 9.6 |
| 2.68 (band edge) | 2.74 | 0.22–0.43 | 0.05–0.20 | 124–470 | 0.021–0.15 | 165–1 190 | 2.2 | 11.6 |

λ* = 8π/(κ*/λ*). The research runs (λ = 5000) lie above every value here.

**One unit check.** The archive's band is in S's bare mass μ = 1, which is S's clock at the saddle. In the sea S's
clock is 2. If the reset clock is S in the sea, h doubles and κ*/λ* goes up fourfold.

## 3c. Against the clock-stop band, and GEAR-104

- The band is an input to variants A and B, so they agree with it by construction. What is new is the factor the
  region's own w makes (band 0.22–0.45 instead of 0.14–0.37), and the conversion to κ*/λ*.
- **GEAR-104's 1/[x(1 + 9x/32)] is not the coupled-to-uncoupled ratio.** It is the held vacuum step over S's energy,
  0.55 → 0.22 across the band (falling). The bond-to-quartic ratio U/B at the saddle is 0.24 → 0.50 across the band
  (rising). They overlap in range by coincidence and run in opposite directions.
- **Neither can set H.** Both are λ-free and scale-free, so they can fix *where in x* ejection happens (for example
  the band's upper edge, x_A = 2.68, where D_s can no longer be held out), but not κ*/λ*. Gravity enters only through
  H at that x.

## 3d. 0*

| question | answer |
|---|---|
| circular through the one-budget floor? | **No**, for κ*/λ*: it uses the S–D sector and the clock-stop rule only; no f, no V0, no "today". The step to today's H (below) does use f (FORCED) and Ω_de (BORROWED, enters as an O(1) factor). |
| derived or borrowed? | Derived **given** the C6 criterion (Picture) and the clock-stop rule (Hypothesis for G_reset). |
| region's size at ejection? | **Not in the uniform version**, because the region is a Hubble patch: its compactness ~ H²R² ~ 1 by construction. That is exactly why no hierarchy appears. A finite pocket of radius R and compactness Φ would give κ*/λ* = 3Φ/(ρ̃R²): a small κ*/λ* is a large R in matter's units. The hold band does not select R (shorter D waves grow at any x, PR #470). **Open.** |
| ruler dependence | **RULER-DEPENDENT** through the balance's definition (A, B, C: up to ×100 in κ*/λ*), S–D's share s (B), and which of S's clocks is meant (×4). Invariant: the order, λ* ~ 10–1 000. |

**The consequence through the one-budget rule (Derived; uses f FORCED):**

- The rule sets D's step equal to f·V0, so V0 = (κ*/λ*)/(2f) = **0.11–2.3 μ²**. The lane's own H scale,
  √(V0/3), is 0.19–0.87 μ. φ's potential lives at the reset scale.
- Today's dark energy is V0(e^(−φ) − f) + π²/2. For H_today ≪ μ the potential and D's step must nearly cancel:
  δ = e^(−φ_today) − f must be small.

| ε_today = H_today/μ | δ/f required (potential alone) |
|---|---|
| 0.099 (research runs) | 0.09–1.9: no cancellation needed |
| 0.00567 (production) | 3e-4 to 6e-3 |
| 1e-40 (S as a particle, order of magnitude) | 1e-79 to 2e-78 |

- So C6 plus the one-budget rule moves the unpaid number from H/μ into **how close φ sits to the node** of
  V0(e^(−φ) − f), φ = ln 10. With ε ~ 1e-40 that is an 80-decade cancellation: the cosmological-constant problem,
  restated, not solved.
- **Contrast with production (Read):** production's φ_today = 0.179 is far from the node (δ/f = 7.4); its small ε
  comes from a small V0 (8.24e-5 μ², i.e. λ ≈ 1.5 million). C6 plus the rule describe a different universe: λ not
  large, ε small because φ is near the node. Stored data cannot tell them apart.
- **One possible way out (Open):** if the lane's own dynamics parks φ at the node, the cancellation would be
  dynamical rather than tuned. In the present lane φ rolls through the node on its way to turnaround (the potential
  turns negative beyond it), so near the node the dark energy is small but only transiently. Whether anything holds φ
  there is a question for chat.

**Verdict.** The C6 candidate is a non-circular way to set κ*/λ* (λ* ~ 10–1 000, ruler-dependent). Together with the
one-budget rule it closes the system, but the hierarchy reappears as φ's closeness to the node. **Not a derivation of
the scale ratio.** The two things that could change that are both open: a finite region size selected by the dynamics,
or a mechanism that parks φ at the node.

## Not done

A finite pocket at ejection (needs R). The reset with S in the sea (h doubles). φ's kinetic share at "today" (left out of
the δ table; it only tightens the requirement). GEAR-91's mass factor e^(−2βφ): between the reset and the node φ moves
by at most about 2.3, so the factor is at most 1.2; harmless here.
