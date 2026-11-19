# Floquet scan of the k=0 D mode above the band (research only)
Script research/floquet_above_band_v1.py; raw output summary.txt (1.7 s). Toy S s''=-(s+0.5s^3), D linear d''=(1-1.5 s^2)d, real S orbit, no friction in the Floquet exponent. Design and amendment (exact adiabatic descent, dust and radiation bounds) fixed in the coordination log before the run.
- Method check: mu = 0 at x=1.5, 2.0; onset at the band edge (trace -2.004 at x=2.67, mu 0.029), i.e. the k=0 instability starts at x=8/3 from the Hill equation itself.
- (a) POSITIVE-ALL above the band (min 0.27 at x=3.0). mu: 0.35 at 3.27, 0.55 at 5, 0.74 at 16, 0.81 at 157; saturates near 0.81.
- (b) MATCH: mu(3.27)=0.346 vs GEAR-100's 0.35 (amplitude reading). The 0.35 is the toy S, not a double count.
- (c) Descent 157 -> 8/3, H=H_b*sqrt(x/x_b): SATURATES ABOVE for all three H_b and all three laws. N_net (exact law): 82.9 (H_b 0.0057), 46.0 (0.01), 21.6 (0.02); dust-law 59.9/33.3/15.6, radiation-law 89.9/49.9/23.4; descent time exact/dust = 1.36, radiation/dust = 1.50. Never damped (min mu 0.09 vs 1.5H <= 0.03).
- (d) mu ~ constant (spread 2.3% over x=48..160), mu/sqrt(x) falls (spread 79%): the paper estimate (marginal quartic limit) holds; "growth rises with sqrt(x)" does not.
Caveats: small-amplitude exponent only; no back-reaction (saturated D is paid for by S); k=0 only; x0=157 is an inherited initial value, no population shown to exist there; seed needed to saturate = e^-N_net (about 4e-10 at H_b 0.02, below 1e-16 rounding only for H_b <= ~0.01).
