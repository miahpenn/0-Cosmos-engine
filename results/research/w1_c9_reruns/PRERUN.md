# C9 W1 reruns: pre-run declarations (written before any run)

Go: Jeremiah, "Go" (projects timeline, 2026-11-24 01:51 UTC), relayed by the coordinator. Engine main 6b5728d (#575 merged), vendor d6052d6.
Gates: GATES-W1-reruns-C9-Nov24-2026 (Drive 1J7rzIr0Rn8rlyKOlFrRmYcuiUmZWktP3uUtQVMuvoR0) with its Amendment 1 (11qsWXifjhhLziO4F8Kvy3AA0nYSfd6l6wOQI32L8dRU); G3.1's control and G3.2 replaced by DESIGN-W1-admission-control-Nov24-2026 (1TaWODevd-4CIgcMjCcqhOEYHQJqOG-9RLJH1dyJTL70) with its Amendment 1 (1aSXD5T_P3i21o8AB9d-Li06CS32tIeQ8OIgj-1c0jok). G3.4's 0.970 band lapses (design section 8); K5 replaces it.
Input: native.npz, sha256 16151c445c9dcfcc34cf61e28e268e991ab2bfea629beaa3aa52fadc4c5a98ea (Drive 1OUObtOwD7B51pDMiTacKJEjbavVBqv6D). H0_chart = the initializer's H0.

Run B (eight_node_areal on) is the admission run. Run A (engine default integral) is run for the readings C9 defines only on Run A (G3.3, G4.2, G4.3).

Implementation choices, fixed now (the gate documents leave them open):
1. G0.1 label fields: u from members.json drhat/dr squared; theta0 0 and snapshot phase 0 (native.py calls snapshot(x, 0.0)); matching rule and toy radius from Aster's README/members.json.
2. G0.2: the slice's r is the grid's centres (exact, by construction); the producer's nodes H.r/sqrt(u) were asserted within 2e-14 of them by Aster's native.py (Read). Recorded beside the PASS.
3. G1.5 (Kr frozen): Kr is not returned by the initializer. Checked as: the momentum projection that forms j is called exactly once per build (Kr = Kt + r Ktp + 4 pi r j has no other loop dependence).
4. G2.1: E_ring = integral of 4 pi r^2 (rho_ring - rho_none) over r < 25 and over the box, flat weight and areal weight (1/sqrt(B)), rho from the engine's total stress-energy, Run B. f_ring is NOT computed: its definition (coupled share, derivation section 8.1) is not in the repo or the gate documents. The 4.207e-10 P_D value is a radiated power per unit time; max |P_D| on the slice is a field amplitude, so that row is reported with this mismatch stated, not read as a pass or fail of the toy power.
5. G3 arms (design section 6): none (F), K sech^2 k 10, K sech^2 k 100, K Gaussian k 10; injection wraps the engine's own integral function for that run so f becomes f - 8 pi r^2 drho_c in all six iterations; readout = research b7_w2_origin_check_v1.evaluate post-lapse H. K3 tolerance at dr 0.05 computed now by the amendment's rule: 8.268e-13/2.2306e-9 + 0.009 = 0.00937 (report-only).
6. G4.1 gate V needs #568's manufactured test slice, which is built from the v7 member file (not available here). It is built by #568's same fixed rule from the aligned dr 0.05 member instead, and labelled so.
7. G4.4 reading: PASS if max over cells 0-40 of |B_engine - B_ref| <= max over cells 0-40 of the reference's uncertainty u_B = |B5 - B3| + |B5 - B5b| (#568's definition). The cellwise count of cells where |dB| > u_B is reported beside it.
