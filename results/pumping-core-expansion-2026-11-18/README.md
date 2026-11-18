# Pumping under expansion, homogeneous 2-DOF core (Nov 18 2026)

Script `research/pumping_core_expansion_v1.py` (`python research/pumping_core_expansion_v1.py run`), tables in `summary.txt`. Design and tolerances: coordination log, run request v2 (Nov 18). Research only; nothing in engine/ touched.

## Result (Measured)
1. **The k = 0 core does not pump at all.** At H = 0, tiny D, sigma is within +-0.002 of zero for x0 = 1.4-2.67 (D stays 1e-8; max change 15%). Gamma0 = 0 here, not the 0.16-0.28 of the Floquet analysis.
2. **Reason (Derived, checked by item 4):** in the band D's k = 0 frequency is sqrt(1.5 n - 1) with n = x/2 (0.7 at x = 2), far from S's 2 and the ring's 4, so the first Mathieu bands (D frequency 2 or 4) are not reached. Growth needs a mode with omega_k = omega_S, i.e. k ~ 2.
3. **Main grid (S and D friction, k = 0):** total D gain G = 0.5-2.3 everywhere, no trend with H, x0 or d0 (small or large D: the same to ~10%). Q1 as posed has no H_G10: G never reaches 10, even at H = 0. The variation is the envelope reading's own jitter on a non-growing oscillation. Q2 is moot (nothing to shift).
4. **Method check (C3-style):** frozen S, sigma at H = 0.0057 is -0.0103 vs the prediction Gamma0 - 1.5H = -0.0086 (within the 0.01 absolute tolerance); the H = 0 value is -0.002. The envelope reading works, but the growth it measures is 0.
5. **Extension (not in the approved design; exploratory; H = 0, tiny D, a finite-k mode, k^2 scanned near resonance 5 - 0.375x):** max sigma = 0.140 (x0 1.6), 0.172 (2.0), 0.224 (2.67) at k ~ 2.0-2.1. This reproduces the Floquet window (0.16-0.28) and fixes where it lives: a finite-k, omega_k = omega_S phenomenon.
6. t_exit at x0 = 2.0 (S friction, x < 4/3): 50.1, 28.1, 16.4, 9.2, 7.7, 4.9, 3.5 for H = 0.003, 0.0057, 0.01, 0.02, 0.03, 0.05, 0.11 (S's x falls as exp(-3Ht) as estimated; at H = 0.11 the first S period already leaves the band).

## What this does and does not settle
- The pump the SD closure cares about is not in the k = 0 core. Whether it develops under expansion needs the finite-k mode with S friction, where there is one more effect the k = 0 case lacks: as x falls, the resonant k shifts (k^2 = 5 - 0.375x), so the mode drifts off resonance as S decays. Not run (needs a go).
- Q1 and Q2 stay open; the core cannot answer them. The paper estimate H ~ Gamma/17..Gamma/11 for G = 10 is untested.
- Caveats from the design stand: homogeneous (the finite-k row is a single linear mode with no back-reaction on S), constant H.
