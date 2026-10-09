# Preregistration: domain-comparability run for the D time-shape campaign (Lumen, 2026-10-09)

Written BEFORE any run. Purpose: test whether the amplitude-dependent final-time trend and the
D-driven timing structure from campaign 37873115509 (r_max = 160) survive at r_max = 80.

## Fixed configuration
- Cases, in this order, one at a time: D0 (amplitude 0), Dbase (amplitude 1e-10).
- N = 160, r_max = 80.0, CFL = 0.0075, radiation ON, final time 22.5, checkpoint interval 0.25.
- Production kernel and CMC gauge at the pinned commit of the launch; no equation, boundary, or
  parameter change. Same time-shape recorder and export as the r_max = 160 campaign.
- Onset rules exactly as in D_TIME_SHAPE_ANALYSIS_2026-10-09.md: deviation-from-initial marker,
  primary f = 0.25, sensitivity f = 0.10 and 0.50; lags relative to D_active_at_source_probe.

## Already-admitted context (not new runs)
- D0 at r_max = 80, CFL 0.0075, N = 160 (admitted four-way file): tau_central = 22.349398.
- Campaign D0 at r_max = 160: final tau = 22.343354.
- Control domain effect at fixed D = 0: 0.006 in final tau. Recorded as context, not a new result.

## Required artifacts per case
- Ledger with 3,000 rows, 90 checkpoints, required fields finite
  (D_center, PD_center, D_at_source_probe, PD_at_source_probe, D_active_at_source_probe).
- Pinned commit, ZIP SHA-256, summary SHA-256, per-case ledger SHA-256.
- Same acceptance gate as run 37873115509.

## Questions and decision rules (fixed now)
1. Trend survival. Final-time gap G = tau(D0) - tau(Dbase) at r_max = 80. Record G and compare with
   the r_max = 160 gap of 1.503. Report the value; no pass/fail threshold on the magnitude.
   Decision: if G > 0 (the control is above Dbase), the amplitude-dependent suppression survives the
   domain change in sign. If G <= 0, the r_max = 160 trend is a domain artifact and is NOT carried forward.
2. Timing structure. For Dbase at r_max = 80: are the D-driven onsets (D_center, D_K_source_center,
   D_growth_rate_center, alpha_r_max, tau_rate) within about one time unit of the D-active onset,
   as at r_max = 160? Report the lags and the spread. Decision: the structure is carried forward only if
   the lags stay within about one time unit of D-active onset at the primary marker.
3. Marker sensitivity. Report f = 0.10, 0.25, 0.50 in full. Do not select the marker after seeing this run.

## What this run does NOT decide
- Whether the amplitude trend has a physical mechanism.
- Whether the local-to-cosmic time map or any cycle exists.
- Resolution dependence (N = 320 is a separate, later, conditional step).

## Launch and gate
- Launch by a new workflow with its own tag (not [0star-d]), added in a separate commit that changes
  only the workflow and the case list. No physics code change. The non-interference gate must pass first.
- One case at a time: D0 first, then Dbase after the D0 artifact is verified.
- Owner (Jeremiah) controls launches; no other campaign runs during these two.
