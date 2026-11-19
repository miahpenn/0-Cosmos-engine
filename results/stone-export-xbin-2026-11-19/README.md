# Offline x-binning of the stone-export windows (research only)

Tests chat's hypothesis: near-gap and mid are one pair channel read at different x(t) (omega_S = sqrt(4 - 1.875x), x < 1.07); broad is the
single channel 2 omega_S. Same 12 runs (deterministic rerun; probe series and x(t) kept in memory, not saved: runs take ~15 s).
Each window's band-passed outward flux is binned by x at emission (t - r/v_g, v_g from the window's mean omega) at the centre (r < 1) and core (r < 5).

**Result: not supported.**
- Flux-weighted instantaneous omega in each window is flat in x (R 24 f 0, centre: near-gap 1.46-1.47, mid 1.65-1.68, broad 2.24-2.35 over x 0.3-5). The mid range over all 12 runs is 1.64-1.79 with no trend; the pair line would sweep 1.85 -> 1.41 over x 0.3 -> 1.07.
- Most flux is emitted at x > 1.5 (centre: ~60-80% of each window), where neither formula applies. About 15-25% of a window's flux is emitted at x < 1.07.
- Caveats: Hilbert omega is confined by the rectangular mask (near-gap [1.35, 1.55] is narrow, so flatness there is partly by construction; mid and broad are wide enough to show a sweep). The centre's x (r < 1 mean) is not one clock for a stone with no single omega_S. Emission-time uses the dispersion v_g.
- **Total w_stress** of the whole field (1.35-3.6): 0.055-0.155 over runs; 3R: 0.104-0.155 (mostly 0.105-0.13); 1.5R: 0.055-0.126.
