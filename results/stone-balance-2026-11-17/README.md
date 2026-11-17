# Stone balance data: what a stone rings out against what it takes back (research only)

HANDOFF Nov 17 (1), task 4. Target question (AMENDMENT-CORE-SHEET-2 §3): is there a stable stone amplitude x_A* where
the D_r a stone rings out into the sea equals what it takes back into S inside the pocket?
Script: `research/stone_balance_v1.py`. Output: `balance.txt` (all cells, both rulers, every sample time). Stored cells
only; no new runs; no fitting; nothing in `engine/` or `vendor/` changed.

Marks: **Derived**, **Read** (measured from the stored cells), **Open**. Scaled units: S's bare mass μ = 1, so H/μ = H.
μ_M = (3/32) x_A is the Mathieu rate (PR #491, SEA_ANSWER), not S's mass.

## Headline

- **No crossing with a stable slope is found.** Across x_A = 0.02–0.35 at H/μ = 0.008–0.065, with the
  line-resolving window, gain − loss is **positive** in every sample at x_c 1.5 down to x_A ≈ 0.026 and at x_c 2.0 down
  to x_A ≈ 0.07. It falls *with* x_A (positive slope in 9–11 of 12–13 neighbour pairs per cell and core radius): the **unstable**
  sign by the handoff's criterion. It is a time track, though, along which x_A and H fall together (below).
- **The one zero-crossing is not a function of x_A.** At x_c 2.0, gain − loss changes sign at t = 27–30, between
  x_A = 0.075 and 0.066 (R_c 6) or 0.068 and 0.061 (R_c 8). The slope is +0.008 to +0.013 (unstable sign) and all
  three x_c 2.0 cells agree. But x_c 1.5 passes the same x_A (0.066, t = 21) with gain − loss = +4e-4. At x_c 2.0 the
  sign change comes as the central knot grows (x(0) 0.5 → 0.8 while the core's weighted x_A falls). **Read: it follows
  the knot, not x_A.**
- **The stored cells never reach the regime the question is about.** After the release, which the windows can
  reach only from t ≈ 6–9, the core's x_A is ≤ 0.35. That is far below the hold band (1.27–2.67) and below PR #491's
  two-way condition x_A R_c > 7.5 (x_A > 1.25 at R_c 6, > 0.94 at R_c 8). μ_M R_c/v here is 0.015–0.31. The 2c channel
  is also **closed by rate** throughout: 2μ_M = 0.01–0.065 < 3H = 0.12–0.19 (as SEA_ANSWER found). **What the stored
  data measure is the balance of a released, sub-threshold remnant, not of a held stone.**

## What was measured (definitions: Derived; values: Read)

Per cell, core radius R_c, averaging window W, and time (rates per unit of the core's S energy E_S, per proper time):

| quantity | definition |
|---|---|
| x_A | x_loc = 2⟨s²⟩ over one S period (`handoff_contract_v1.x_local`), weighted by S's energy over r < R_c. x(0) = the centre's. |
| E_S | S's ring with its in-sea mass: ½ṡ² + ½(∂_l s)² + ½s²(1 + (3/2)D_slow²) + s⁴/8, summed over r < R_c |
| **loss l** | outward flux of D's fast part d_f = D − D_slow (`reset_handoff_v1.d_slow`, box2) through the sphere r = R_c: −ḋ_f ∂_l d_f R², averaged over W. Split by line (complex demodulation as `handoff_contract_v1.line_energy`): **l_2c**, the line at ω = 2 (S² at 4 rings D at 2, moving at 0.71: PR #491's D_r), and l_gap (the gap line √2). |
| **gain g** | the coupling's power into S from D's fluctuations inside the core: −(3/2) s ṡ d_f² (the coupling term of V is (3/4)s²d²), summed over r < R_c and averaged over W. The two-way exchange of PR #491. |
| also shown | g_x (cross term −3 s ṡ D_slow d_f, the forced ring), g_h (work of the moving held D on S's mass, (3/4)s² ∂_τ D_slow²), l_S (S's own flux out through R_c), and the closure residual dE_S/dτ − (g + g_x + g_h − l_S) |

ṡ = PS and ḋ = PD in proper time (checked: dS/dt = α PS and dD/dt = α PD, correlation 0.99998).

Cells (Read from the files): f_x20, f_x15 (x_c 2.0 / 1.5, N 240, r_max 24, closure ii, cold, output every 0.1,
t ≤ 30); u_x20_cold, u_x15_cold (same at output every 0.5, t ≤ 70); u36_x20_cold (r_max 36). At overlapping times
the 0.5-cadence cells reproduce the 0.1-cadence values to 7–25% in g and to ≤ 2% in l (e.g. t = 9, R_c 6:
g 7.4e-4 vs 6.9e-4, l −1.97e-3 vs −1.94e-3). r_max 36 agrees with r_max 24 to ≤ 5% wherever |g| > 1e-4.

## Table (W = 10.7; rates per unit E_S per unit time; Read, except the #491 columns, which are Derived forms)

| cell | R_c | t | x_A | x(0) | H/μ | gain g | loss l | g − l | g_h | l_S | 2μ_M | 3H | μ_M R_c/v | l_2c/E_2c | 3v/R_c | g/E_2c |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| x_c 2.0 | 6 | 9 | 0.342 | 0.26 | 0.065 | +6.9e-4 | −1.9e-3 | +2.6e-3 | +2.4e-2 | −5.2e-2 | 0.064 | 0.194 | 0.27 | −0.013 | 0.35 | +0.023 |
| | 6 | 15 | 0.173 | 0.18 | 0.051 | +5.0e-4 | −2.0e-4 | +7.0e-4 | +8.2e-3 | −1.4e-2 | 0.032 | 0.153 | 0.14 | −0.002 | 0.35 | +0.016 |
| | 6 | 21 | 0.105 | 0.39 | 0.041 | +1.6e-4 | −5.1e-5 | +2.2e-4 | +3.6e-3 | −5.3e-3 | 0.020 | 0.122 | 0.08 | −0.0005 | 0.35 | +0.005 |
| | 6 | 30 | 0.066 | 0.79 | 0.028 | −2.5e-5 | −2.2e-5 | −2.8e-6 | +1.5e-3 | −2.2e-3 | 0.013 | 0.085 | 0.05 | −0.0002 | 0.35 | −0.0007 |
| | 6 | 39 | 0.050 | 0.64 | 0.018 | −1.9e-4 | −1.3e-5 | −1.8e-4 | +6.7e-4 | −1.2e-3 | 0.009 | 0.054 | 0.04 | −0.0001 | 0.35 | −0.005 |
| | 8 | 9 | 0.287 | 0.26 | 0.065 | +5.5e-4 | −1.7e-4 | +7.2e-4 | +2.2e-2 | −5.6e-3 | 0.054 | 0.194 | 0.30 | −0.001 | 0.27 | +0.023 |
| | 8 | 15 | 0.151 | 0.18 | 0.051 | +4.2e-4 | −2.6e-5 | +4.5e-4 | +7.2e-3 | −1.6e-3 | 0.028 | 0.153 | 0.16 | −0.0002 | 0.27 | +0.015 |
| | 8 | 21 | 0.094 | 0.39 | 0.041 | +1.4e-4 | −1.1e-5 | +1.5e-4 | +3.2e-3 | −7.5e-4 | 0.018 | 0.122 | 0.10 | −0.0001 | 0.27 | +0.005 |
| | 8 | 30 | 0.061 | 0.79 | 0.028 | −2.2e-5 | −4.5e-6 | −1.7e-5 | +1.4e-3 | −3.5e-4 | 0.011 | 0.085 | 0.06 | −0.0000 | 0.27 | −0.0007 |
| x_c 1.5 | 6 | 9 | 0.246 | 0.17 | 0.065 | +1.2e-3 | −2.0e-3 | +3.2e-3 | +1.3e-2 | −4.3e-2 | 0.046 | 0.194 | 0.20 | −0.018 | 0.35 | +0.083 |
| | 6 | 15 | 0.118 | 0.07 | 0.051 | +5.6e-4 | −3.8e-4 | +9.4e-4 | +4.9e-3 | −1.4e-2 | 0.022 | 0.153 | 0.09 | −0.003 | 0.35 | +0.017 |
| | 6 | 21 | 0.066 | 0.05 | 0.041 | +3.2e-4 | −9.1e-5 | +4.1e-4 | +2.2e-3 | −6.0e-3 | 0.012 | 0.122 | 0.05 | −0.0006 | 0.35 | +0.007 |
| | 8 | 9 | 0.205 | 0.17 | 0.065 | +8.8e-4 | −1.5e-4 | +1.0e-3 | +1.2e-2 | −5.2e-3 | 0.038 | 0.194 | 0.22 | −0.002 | 0.27 | +0.071 |
| | 8 | 15 | 0.100 | 0.07 | 0.051 | +4.6e-4 | −2.9e-5 | +4.9e-4 | +4.3e-3 | −1.7e-3 | 0.019 | 0.153 | 0.11 | −0.0002 | 0.27 | +0.016 |
| | 8 | 21 | 0.057 | 0.05 | 0.041 | +2.7e-4 | −1.2e-5 | +2.8e-4 | +1.9e-3 | −8.4e-4 | 0.011 | 0.122 | 0.06 | −0.0001 | 0.27 | +0.007 |

Rows t ≤ 21 are the 0.1-cadence cells; t = 30 and 39 are u_x20_cold (0.5 cadence). E_2c is the 2c line's energy inside the core. It is
1–4% of E_S, but the window leaks the √2 line into it (≈ 25% in energy, handoff Amendment 1), so the per-E_2c columns
are upper-bounded in size, not exact.

## Readings

1. **Loss is not outward (Read).** The D fluctuations' flux through R_c is **inward** in every sample, ruler and cell:
   −2e-3 to −1e-5 of E_S per unit time. The 2c line carries ≤ 30% of it, also inward. Per unit of its own energy the
   2c line moves at |l_2c/E_2c| ≤ 0.018, against 3v/R_c = 0.27–0.35 for a freely escaping 2c wave field (this note's
   geometric step: uniform fill, speed v = 0.71). **Nothing in the core is ringing out at v.** With the channel closed
   by rate (2μ_M < 3H), the core makes no fresh D_r, and the inward flux is D carried in with the infall. STONE_SPLIT
   found the same lag (ν ≈ 0.7) and infall absorption for S.
2. **Gain is small and runs D → S (Read).** g = +1e-4 to +1.2e-3 of E_S per unit time while x_A ≳ 0.07 (both
   cells). That is the opposite sign to Mathieu pumping (S → D at 2μ_M per unit D_r energy). Per unit 2c energy it is
   +0.005 to +0.09, the same order as 2μ_M but reversed. At x_c 2.0 after t ≈ 27 it turns S → D (−2e-5 to −2e-4) as the
   knot grows. **Read:** inside these remnants the net exchange is D handing energy back to S, except at the knot.
3. **Both are small against the stone's other channels (Read).** g and l are 10–100× smaller than g_h (the held D's
   work on S's mass as D_slow relaxes: +1e-3 to +2.4e-2) and l_S (S's own energy flowing **in** through R_c: −2e-4 to
   −5e-2, the infall). The stone's S energy is set by infall and the held step's relaxation, not by the D_r exchange.
   The closure residual of S's budget is 4e-4 to 1.2e-2 of E_S per unit time, the same order as g_h, so S's budget
   does not close at the level of g. The residual likely includes the dilution work of S's quartic (radiation-like)
   part, which is not in the bookkeeping (**Open**).
4. **Window ruler (Read).** With W = 2π, g's sign follows a beat of period ≈ 5–6 (the 2 and √2 lines are not
   separated) and gain − loss changes sign 6–16 times per record with alternating slopes. W = 2π does not resolve the
   exchange: only W = 10.7 rows are used for the headline. **Core-radius ruler:** R_c 6 and 8 agree in every sign
   and in the crossing time. Loss falls 5–12× from R_c 6 to 8 (the infall front lies between them); gain changes
   ≤ 30%.
5. **x_A ruler (unpaid).** Core-weighted x_A and the centre's x(0) differ by up to ×12 late at x_c 2.0 (knot). Any
   "x_A*" depends on which one is meant.

## PR #491's forms beside the measured rates (Derived forms evaluated at the measured x_A and H)

| form | value range here | measured counterpart |
|---|---|---|
| 2μ_M = (3/16) x_A, the 2c line's in-band energy growth | 0.003–0.065 | g/E_2c = −0.005 to +0.09 (sign reversed while x_A ≳ 0.07) |
| 3H, D_r's dilution | 0.025–0.19 | above 2μ_M everywhere: closed by rate |
| π μ_M²/H, the crossing exponent | 0.0006–0.05 | ≪ 1: no crossing gain |
| μ_M R_c / v, the return condition (> 1 ⇔ x_A R_c > 7.5) | 0.015–0.31 | < 1 everywhere: #491 says one way (outward), but no outward flow is measured |
| 3v/R_c, escape of a filled sphere (this note's step) | 0.27, 0.35 | \|l_2c/E_2c\| ≤ 0.018, inward |

## What is settled and what is not

- **Settled (Read):** in the stored cells' post-release remnants (x_A ≤ 0.35), D_r is not rung out through the core
  boundary; the exchange runs weakly D → S; gain − loss > 0 down to x_A ≈ 0.026 (x_c 1.5) and 0.07 (x_c 2.0) with the unstable slope sign; the
  only sign change follows the knot, not x_A. Ruler-safe in R_c; needs W ≥ one beat (10.7).
- **Not settled (Open):** the balance in the regime the question names: x_A in the hold band or above #491's return
  threshold, and 2μ_M > 3H. The stored cells don't have it at any windowable time. x_A and H/μ also fall together
  along each record, so their effects can't be separated here.
- **The minimal run that would settle it** (not made; proposed): a core held at x_A ≳ 1 (R_c 8) for at least two
  beat windows (≈ 21 time units) after D_slow settles, output every 0.1, r_max 24. At research H/μ ≈ 0.06 that is
  also what opening the 2c channel needs (x_A > 16 H/μ ≈ 1). No stored initial data does this: x_c 1.5 and 2.0 both
  release by t ≈ 6–8. Whether any initial data holds a stone that long is itself **Open**, and choosing one is for
  Jeremiah / chat. A second H/μ, at fixed x_A, would separate the two dependences.

## Unpaid inputs

The core radii (6, 8); the windows (2π, 10.7); D_slow's separation rule (box2, two boxcars TD and T4); the split of
the coupling energy (D_slow's share in S's mass; D_fast's as exchange); the core weighting of x_A (S's energy); the
demodulation's √2 leak into the 2c line (≈ 25%); the 3v/R_c escape form (uniform fill, this note); H read as
−⟨K⟩/3; samples with H ≤ 0.005 (turnaround near t = 55) excluded from the crossing count.
