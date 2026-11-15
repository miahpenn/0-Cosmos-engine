# 2c, derived: the pocket's source spectrum against the sea's thresholds

Owed by amendment 1 to the leftover note (PR #489). RESEARCH ONLY. Stored states (x_c 2.0 and 1.5, r_max 24); no
runs. `research/sea_answer_v1.py` → `sea_answer.txt`.

Marks: **Derived**, **Read**, **Open**. Scaled units; frequencies in local proper time.

## Summary

- **Correction to amendment 1 (kept):** the line that matters for pumping D is not S's own frequency but S², which
  rings at twice S's clock. In the sea that is **4.0, above D's parametric threshold 2√2 = 2.83 everywhere S rings**,
  not only at the knot.
- **The sea is open by frequency, closed by rate, in these runs.** The growth of D_r from S's pumping is 4–10 times
  smaller than D_r's Hubble damping after t = 12. Nothing grows parametrically; the D_r seen in the sea is radiated
  from the pocket.
- **The balance is a ruler effect.** Growth is λ-free (it depends on x); damping is ∝ H/μ. At production's ratio
  the channel opens for x_A above about 0.1. **The toy universe hides this channel.**
- **The leftover verdict stands:** whatever the sea answers with is paid by S's ringing, not by D's step.

## 1. The thresholds (Derived)

About the sea (D_s at rest, d = √2, s = 0): S waves have ω² = 4 + k², D waves ω² = 2 + k².

| coupling (from V) | what it pumps | forced wave needs | parametric needs |
|---|---|---|---|
| (3/2) s² d | D, at the frequencies of s² | pump ≥ √2 | pump ≥ 2√2 = 2.83 |
| (3/2) s d² | S, at the frequencies of d² | pump ≥ 2 | pump ≥ 4 |
| d³/2 | D's own harmonics | — | — |

## 2. The source lines (Read, t = 18–54)

| place | S | S² (pump on D) | D's ring | D² (pump on S) |
|---|---|---|---|---|
| pocket r = 3 | 2.01 | **4.03** (96%) | 1.36 | 1.38 |
| pocket edge r = 8 | 1.99 | **4.01** | 1.41 | 1.41 |
| the sea r = 15 | 2.03 | **4.07** | 1.44 | 1.44 |
| knot (x_c 2.0) | 1.83, sidebands 0.63, 3.07 | 1.18, 2.49, **3.67**, **4.90** | 1.21, 2.46 | 1.23, 2.44, 3.67 |
| centre, x_c 1.5 (no knot) | 2.01 | **4.04** | 1.43 | 1.43 |

- **S² is above D's parametric threshold everywhere** (bold). This is the Lamé conversion already on file: S rings
  D_s into D_r (PR #470; core-sheet amendment §3).
- D² lines (1.2–1.44, and 3.67 at the knot) stay below S's parametric threshold 4. **S is not pumped by D**, except
  marginally at the knot.
- S itself (2.0) is just above S's propagating gap: S radiates weakly. D's ring at 1.41–1.44 sits at D's gap: D_r
  near the pocket is barely propagating.

## 3. Rate against damping (Derived; Read)

**The rate (Derived).** With s = A cos 2t in the sea, D's perturbation obeys a Mathieu equation:

    δ'' + [2 + k² + (3/4) x_A (1 + cos 4t)] δ = 0

The first resonance is at 2 + k² + (3/4)x_A = 4, i.e. D waves at ω = 2, k ≈ √2. The amplitude grows at

    μ = (3/32) x_A per unit time.

D_r's amplitude in the expanding cell is damped at (3/2)H.

**Read** (pocket core r < 6, x_c 2.0; x_c 1.5 is the same within 30%):

| t | x_A | μ | (3/2)H | net |
|---|---|---|---|---|
| 12 | 0.18 | 0.017 | 0.086 | −0.069 |
| 24 | 0.044 | 0.004 | 0.054 | −0.050 |
| 40 | 0.014 | 0.001 | 0.026 | −0.024 |
| 50 | 0.010 | 0.001 | 0.009 | −0.008 |

- **Closed by rate throughout** the settled runs. S's concentration dilutes faster than H falls.
- **Opening condition (Derived): x_A > 16 H/μ.** At the research ratio (H ≈ 0.1 early) that needs x_A > 1.6, i.e.
  inside the hold band, where D_s is held out and there is no sea to ring. At production's H_today/μ = 0.0057 it
  needs only x_A > 0.09.
- **So the answer is RULER-DEPENDENT through H/μ, and in the opposite direction to most toy effects.** The research
  universe expands too fast for the sea to answer. At small H/μ a settled pocket with x_A ≳ 0.1 would ring D_s into
  D_r. This belongs in class D of PR #486's classification ("cannot say from scaling") and should move to: "opens as
  H/μ falls".

**D_r in the sea (Read, crude).** Comoving Σ ring² over r = 10–20 rises 0.036 → 0.116 (t = 12 → 24), then falls to
0.055 by t = 48 (x_c 2.0). With growth closed, this is D_r radiated from the pocket and then diluting, not grown in
place.

## 4. Verdicts

| question | verdict |
|---|---|
| does the pocket's source reach the sea's thresholds? | **yes**: S² at 4 > 2√2 everywhere S rings in the sea (correcting amendment 1) |
| does the sea answer in these runs? | **no** (Read): growth 4–10× below Hubble damping |
| in general | **RULER-DEPENDENT**: opens when x_A > 16 H/μ; closed in the toy, open at small H/μ |
| does the answer remove energy from D's step? | **no** (Derived): D_r is paid by S's ringing; D_s at rest is the zero of energy |
| route to the leftover | **not supplied**, unchanged |

**Open:**

- Whether, at small H/μ, the sea's answer drains S's ringing (matter turning into D_r) fast enough to matter for the
  cosmos lane's S–D species. That is "D_r's return paths", parked.
- Whether this connects to Jeremiah's "where dark energy is rung".
