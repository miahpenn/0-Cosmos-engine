# Gear-to-code discrepancy map (read-only, 11 November 2026)

HANDOFF v2 to projects, task T3; REV20 open item 6. **No code was changed to produce this map.** Repo main at 30cd67f.
Each row: where the engine departs from the gears (or the gears from the engine), file and line, the gear, the
current status, and what it affects. Ranked by how many current results depend on it. "Production" means the default
code path (`engine/`); research settings that already exist are named.

Unpaid inputs behind the rows: λ, V0, φ_today, β_D, β_S, D's seed, Ω_σ,today, η (REV20 §7). G_N = 1 throughout;
G_reset (the 0\* reset) has no code at all (row 9).

| # | Departure | Code | Gear | Status | Affects |
|---|---|---|---|---|---|
| 1 | **Two floors.** The COSMOS lane's negative floor is written by hand (`D = 0.10`, so F = 0.10·V0) and the lane carries no D, while the gears put the floor in D's ejected vacuum −ΔV(φ). With λ > 0 and φ on, both would count. | `engine/cosmos.py:27`, `:66`–`:71` (potential and its derivative); `engine/scalar_system.py:12` (V0) | GEAR-91 §2, §4(d); GEAR-99; HANDOFF v2 | One-budget lane built as research (`research/one_budget_cosmos_v1.py`, PR #469): reproduces production to 10 digits at β_D = 0; double count moves turnaround 1.988 → 1.37 e-folds. Production unchanged. Adoption: Jeremiah, after β_D. | every run "at production V0"; the all-on pocket |
| 2 | **GEAR-67 quartic off.** Production λ = 0 (bit-identical to the pre-#456 engine). | `engine/scalar_system.py:24`; used in `:55`–`:76`, `engine/v55_initial.py:145`–`:147` | GEAR-67 3(e) | Explicit since #456; research runs use λ = 5000; λ unpaid (tied to the floor under the identity: 1.52458e6 at φ_ref = 0). | all S–D pocket, hold, emission and bag results |
| 3 | **No GEAR-91 coupling in production.** S/D masses do not depend on φ; φ's equation has no 8π ∂V_SD/∂φ. | `engine/scalar_system.py:79` (RHS), `:122` (projection) | GEAR-91 and amendments 1–4 | Research switch `research/gear91_coupling_v1.py` (PR #465; exchange balanced to 3e-9). β_D = −0.04 borrowed, β_S = 0 by choice. | the β_D-dependent turnaround (1.739 vs 1.988); dark-matter pocket |
| 4 | **COSMOS lane is one-way.** `exchange_Q` is accepted and ignored; nothing from the local lane reaches the cosmos. | `engine/cosmos.py:128`–`:150` | Oct 7/9 objective; REV19 §4c; backlog 3 | Open; the largest structural job. | every "self-sustaining loop" statement |
| 5 | **Two arrest closures.** The reduced lane carries the GEAR-63/65 cubic barrier −(2π/3)Gα²D³ (conditional on D-to-matter conversion); production carries the GEAR-67 quartic slot. #457 stops them silently disagreeing at λ ≠ 0 but they remain two theories. | `engine/local.py:1`–`:13`, `:38` (guard) | GEAR-63/65, GEAR-67, AMENDMENT-GEAR-73 | **Decided Nov 11** (DECISION-map-row-5): quartic is the one closure; cubic barrier and damping retired to provenance (owner-reviewed PR). | reduced-lane results; Term B |
| 6 | **Term A damping (η) absent.** GEAR-61's η_G D² Ḋ is imposed in the research archive mode only; production has no dissipation for D. The released ΔV therefore floods (GEAR-95; bag controls, PR #470). | `engine/local.py:25` (declared, unused); none in production | GEAR-61, GEAR-73, REV19 §3 | η WALLED; candidate origin (S radiated by an oscillating D) not derived. | whether D can rest; every hot-start pocket run |
| 7 | **Outer boundary not specified by the gears.** Production applies the vendor's incoming-light constraint at a finite radius. | `engine/production_kernel.py:273`–`:283`, `:343`–`:360` | none (gears silent) | Open (backlog 7): 6% mass loss when D reached r_max; Oct 7 outer-cell radiation failure. | any run long enough for D to reach the edge |
| 8 | **Slicing.** Production uses CMC with β = 0; the gears quote lapse-based clocks. | `engine/production_kernel.py:194`–`:199`; `engine/cmc_gauge.py:156` | GEAR-68 (teeth), GEAR-82 | Partly answered: on the Kodama clock D² explains 37–42% (vs 83–86% on the lapse) (PR #466). Lapse statements are slicing statements. | clock claims; stall depths |
| 9 | **No G_reset, no pole swap, no bounce in code.** The cycle module only observes turnarounds ("no imposed reset rule"). | `engine/cycle.py:1`–`:14` | REV18/19 §1 (G_reset), GEAR 50, GEAR-94; Oct 2 pole swap | By design (always derive). Clock-stop candidate (|H| = 0.14–0.37 μ) is analytic only. What the pole swap does to φ: open (chat). | cycle closure; the D–φ ledger |
| 10 | **"Today" is imported.** The present-day state is hard-coded from Ω_de = 0.69. | `engine/cosmos.py:83`–`:105` | GEAR-92 §3, GEAR-93 §3 | φ_today and the 1.988 e-folds to turnaround are not derived. | everything quoted "today" |
| 11 | **Contraction endpoint.** `n_contract = 47.203` is unused; production's cosmos halts ~9.3 e-folds below turnaround; the density endpoint lives in research. | `engine/cosmos.py:33` | GEAR 57; REV19 §2; AMENDMENT-GEAR-93 | Research tools: `research/cosmos_efold_v1.py`, `research/cosmos_endpoint_depths_v1.py`. Planck convention open; moot if the endpoint is the clock stop. | depth statements |
| 12 | **Baryon/dark-matter split.** `f_b = 0.157` is not read; the live split is 0.157 baryons/DM at φ_ref = 0; Planck's 0.157 is baryons/all matter. | `engine/cosmos.py:32`; densities `:88`–`:89`; `engine/v55_matter.py` (RHO_*_PRESENT) | AMENDMENT-GEAR-19, 2, 3; GEAR-75 | Deferred by Jeremiah (one PR with a before/after table later). | COSMOS budget (~0.9% in 1 + z†) |
| 13 | **Shear input.** Ω_σ,today = 1.745e-8 is an assumed input; the shear lane is background only. | `engine/cosmos.py:31`, `:92` | GEAR-80; REV19 §2 | Origin not found; condition s untestable in the machine. | condition s |
| 14 | **Scale assertion.** Sharing one length unit between lanes with S at μ = 1 makes production assert H_today = 5.67e-3 μ (μ ≈ 2.55e-31 eV). | `engine/scalar_system.py:12` with `engine/v55_initial.py` (local lane initialised from the cosmos numbers) | GEAR-92 §2 | V0/μ² not fixed by the equations (HANDOFF-REV20-ready R1). | any reading of S as tangible matter |
| 15 | **φ normalisation.** φ is at 8πG = 1 (stress divided by KAPPA), S/D at G = 1; κ placement in the source action is an owner question. | `engine/scalar_system.py:11`; exchange `:118` | REV19 §4a (Gate A item 4) | Consistent in code (GEAR-91's 8π carried in research); owner question open. | exchange terms |
| 16 | **Saddle ratio g/κ.** Production has κ = 0 (Ω_S = Γ_D = 1). | `engine/scalar_system.py` (mass terms ±½) | GEAR-89, GEAR-90 | Bridge: matter supplies κ = 9x/16 (averaged, single-harmonic). Exact edge 1.270 (PR #470); fails for finite-wavelength modes. | the "main clock" readings |

## Not discrepancies (recorded so they are not re-opened)

- **D as content that curves spacetime.** Correct bookkeeping for a matter-specific field; D is not gravity's scalar (GEAR-97, AMENDMENT-2-GEAR-97; PR #468 checks D, G).
- **Dark matter not coupled to D through A(D).** Matches Jeremiah's recorded position (D matter-specific); A(D) exists only as a chat sandbox script (GEAR-98).
- **Misner–Sharp ledger shift term.** Closed by #454.
- **Older V5.5 drivers.** Re-pointed by #443; provenance in `docs/PROVENANCE_OLDER_V55_PATHS.md`.

## Gear statements contradicted by machine runs (the gears, not the code, need amending)

- "Matter holds D out where it is concentrated": only the uniform mode, only for 1.270 < x < 2.670; no upper band; finite-wavelength modes grow inside the band (PR #470, exact Floquet).
- Bag pockets in the ejected vacuum: with the released energy removed, oscillating S still fills the pocket with D and, with back-reaction, loses 50–95% of its energy (PR #470). Engine bags (x_c = 6) collapse with or without D.
