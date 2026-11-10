# Teeth test on a slicing-independent clock (HANDOFF-to-projects Nov 10, Task 4; backlog item 6)

`research/teeth_invariant_clock_v1.py`. Research only. Production engine as is: λ = 0, V0 = 0 (φ off), CMC, β = 0.
G_N = 1; G_reset not involved.

**Runs.** AMENDMENT-GEAR-68's no-wall set rebuilt here: S = A e^(−(r/7)²), D seed 1e-10 (same width), time-symmetric,
N = 160, r_max = 80, T = 45, A = 0.005, 0.01, 0.02. All three ran to t = 45.

**Clock.** In spherical symmetry the slicing-independent clock rate is the Kodama norm N_K = √(1 − 2m/R)
(|K|² = −(∇R)², `engine.worldtube.invariant_chi`; m the Misner–Sharp mass, R the areal radius): the redshift relative
to Kodama time, which is the Killing time wherever the region is static. Points: D-dominated (|D| > 1e-3, |D| > 5|S|,
r < 60) and untrapped (1 − 2m/R > 0, where K is timelike; ≈ 4% of D-dominated points are trapped and excluded).

## Result (`teeth_output.txt`; columns A = 0.01 / 0.02 / 0.005)

| | CMC lapse ln α | Kodama clock ln N_K |
|---|---|---|
| spread (std) | 1.41 / 1.40 / 1.43 | **0.164 / 0.148 / 0.152** |
| η² on local D² | 0.853 / 0.831 / 0.856 | **0.368 / 0.385 / 0.422** |
| η² on act_D = PD² + D² | 0.852 / 0.829 / 0.856 | 0.368 / 0.385 / 0.423 |
| η² on enclosed active source | 0.249 / 0.298 / 0.205 | 0.031 / 0.036 / 0.041 |
| η² on 1 − 2m/R (N_K's own definition; binning limit) | 0.194 / 0.117 / 0.129 | 0.737 / 0.860 / 0.737 |
| cross-run error, curve from A = 0.01 predicting A = 0.02 / 0.005 (their spread) | 0.629 (1.40) / 0.546 (1.43) | 0.119 (0.148) / 0.109 (0.152) |

1. **Reproduction.** The CMC-lapse numbers reproduce AMENDMENT-2-GEAR-68 (η² 0.850 / 0.831 / 0.856; cross-run
   0.62 / 0.54) to the third digit, now on the untrapped subset.
2. **Most of the lapse's range is slicing.** The CMC lapse spans ln α from −4.9 to +0.5 (spread 1.4); the invariant
   clock in the same points spans −2.0 to +0.37 with spread 0.15–0.16, about a ninth.
3. **On the invariant clock D does not write the clock one-to-one.** Local D² explains 37–42% of ln N_K (against 83–86%
   of ln α), and the curve from one run predicts the others with error 0.11–0.12 against a spread of 0.15 (a 20–30%
   reduction). D's local value is still the best of the local variables, and the enclosed active source and position
   explain almost nothing.
4. **By construction the invariant clock is a function of the enclosed mass ratio 2m/R** (Misner–Sharp m is an integral
   of the local energy). So in spherical symmetry no slicing-independent clock can be a purely local function of D;
   the strong local correlation of AMENDMENT-2-GEAR-68 belongs to the CMC slicing.

## Limits

One invariant clock (Kodama redshift); proper time along matter worldlines would be another, needing tracers. Trapped
points excluded (no timelike Kodama vector there). Flat-space density formulas as in AMENDMENT-2-GEAR-68; S kinetic
term now included in the enclosed source. One resolution, one width, one seed.
