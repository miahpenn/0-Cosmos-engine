# Frozen-S amplitude scan, x0 2.67, H 0.0057 (research only)
Script research/frozen_s_amp_scan_v1.py (frozen_s_trapped_v1 run() + rms output). Raw: summary.txt.
d0: ratio G/G(1e-8) | final D rms: 1e-7 1.000|0.0014; 1e-6 1.001|0.014; 3e-6 1.011|0.043; 1e-5 1.119|0.129; 3e-5 1.092|0.64.
Verdict by pre-fixed rule: INTERMEDIATE. No row with rms<0.1 exceeds 1.05 (not an early excess), but the >5% excess appears at rms 0.13, below the 0.3 the pass clause required. So: a smooth nonlinear-amplitude effect that switches on at final rms ~0.1, 12% peak, 9% at 3e-5 (k* shifts 1.325->1.276 there). Correction: earlier note said rms ~1.3 at d0 1e-5; the measured rms is 0.13.
Does not touch the k_res deviation, which is read at linear amplitude.
