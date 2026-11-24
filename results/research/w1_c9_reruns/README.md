# C9 W1 reruns: results

Pre-run declarations: `PRERUN.md` (commit 2e99a1b, pushed before the run). Engine main 6b5728d, vendor d6052d6. Run B = eight-node on (admission run); Run A = engine default integral. Full numbers: `report.json`. Script: `research/w1_c9_reruns/run.py`.

## Headline

- **G0.5 FAILS on both grids, and G0 is a precondition ("refuse the run if any fails").** Over the last five cells PS and PD sit at 5.1e-4 to 2.9e-3 of their own grid maxima (dr 0.025: PS 5.08e-4, PD 2.62e-3; dr 0.05: PS 7.40e-4, PD 2.86e-3); D - D_rest also fails (5.5e-5, 6.2e-5); only S passes (6.7e-7, 6.4e-7). The tail is the radiating D field plus the Hubble chart term H0 r f_r at r = 60. Under C9 as written this input is not admitted, whatever the readings below say.
- Sequencing (my error): G0 was evaluated in the same pass as G1-G4 instead of before them. The readings below were therefore produced; they are on record, but none of them admits W1 while G0.5 stands.
- **G4c with the known-perturbation control (Run B, dr 0.025): rule PASS, all K checks pass.** R = 9.542e-12 (cell 0; second largest cell 2399), C = 1.410e-11, F = 2.906e-12: (R + F)/(10 C) = 0.088, R/C = 0.677 (Gaussian 0.480), R/F = 3.28. Design prediction 0.68.

## Gates

| Gate | dr 0.025 | dr 0.05 |
|---|---|---|
| G0.1 label, G0.2 alignment, G0.3 H0_chart = H0 (0.0056689681106181795), G0.4 background | PASS | PASS |
| **G0.5 edge** (<= 1e-6 of grid max) | **FAIL**: S 6.7e-7, D 5.5e-5, PS 5.1e-4, PD 2.6e-3 | **FAIL**: S 6.4e-7, D 6.2e-5, PS 7.4e-4, PD 2.9e-3 |
| G1.1 momentum residual | PASS 1.32e-16 (<= 3.89e-14) | PASS 2.05e-16 (<= 3.90e-14) |
| G1.2 P_S detectability | PASS 2.9e8 | PASS 1.9e8 |
| G1.3 j_S vs #574 | PASS, bit-identical | PASS, bit-identical |
| G1.4 loop convergence | PASS 0.0 | PASS 0.0 |
| G1.5 Kr frozen (projection called once) | PASS | PASS |
| G2.1 E_ring vs 1.1525e-4 (1%) | PASS 1.153e-4 (r < 25, flat and areal) | PASS 1.152e-4 |
| G2.1 f_ring | NOT COMPUTED (definition not available) | same |
| G2.1 max abs P_D | 3.26e-5 field amplitude; the 4.207e-10 is a radiated power, not comparable | 3.24e-5 |
| K0 floor | PASS 2.906e-12 | PASS 8.269e-13 |
| K1 response | PASS 1.00e-7 (<= 2e-7), argmax cell 0 | PASS 1.59e-6 (<= 3.2e-6) |
| K2 separation | PASS 48.5 F | PASS 2698 F |
| K3 linearity | PASS 0.0020 (<= 0.03) | 0.0056 (<= 0.0094, report-only) |
| K4 order | PASS 15.82 (in [12, 20]) | |
| K5 reproduction | PASS 1.4102e-10 | PASS 2.2306e-9 |
| **G4c reading (Run B)** | **PASS** (R/C 0.677) | PASS, recorded only (R/C 0.291; ring max cell 38) |
| G3.3 Run A reading | UNDETERMINED: K0, K1, K5 fail on the old integral (rule alone PASS, R/C 0.427 against its own control) | UNDETERMINED, R/C 0.452 |
| G3.5 old control (ring / rescaled pocket), recorded | 0.970 | 62.50 |
| G3.6 direct areal ring/pocket, recorded | 7-node 0.379, 9-node 0.365 (Aster: 0.671, 0.797; not traced) | 1.60, 1.92 (Aster: 1.91, 2.41) |
| G4.1 gate V (manufactured slice per #568's rule, built from the aligned dr 0.05 member) | PASS both runs | PASS both runs |
| G4.2 Run A asymptotic | cells 1-5 share 0.982, U 0.033 | share 0.897, U 0.016 |
|  | **ASYMPTOTIC**: change 0.086 <= 0.1, U <= 0.1 on both grids | |
| G4.3 Run A defect / ring imprint (cells 1-40), recorded | 3.8e-6 | 1.5e-5 |
| G4.4 Run B centre defect <= reference uncertainty (cells 0-40) | PASS: 1.9e-14 <= 1.8e-12; 0 cells over | PASS: 1.5e-12 <= 2.4e-11; 0 cells over |

## Observations (post hoc, not gates)

- The new control cannot read the old integral at all: its own K arms fail K0, K1 and K5, so Run A is UNDETERMINED rather than the 31.4x FAIL in design Amendment 1 (that figure used the old pocket control). Against a control built with the same old integral, the old ring sits at 0.43 C, because the control inherits the integral's error.
- G4.4's max-versus-max reading also passes Run A (1.6e-13 <= 1.8e-12), so on its own it does not separate the integrals; the cellwise count does (Run A 16 cells over, Run B 0).
- Run B gate U fails (u_h / dEF 3.6) because the eight-node B is already within the reference's uncertainty; gate U is a Run A gate.
- Research reconstruct and engine supplied_slice give the same B to 5.6e-16 on both grids.
