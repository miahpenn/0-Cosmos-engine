# W1 corrected slice: 'f under Amendments 5 and 6, with S3 separated (declared before any build or reading)

**Go:** Jeremiah, "Go", Nov 24 03:47 UTC in the Lumen² channel (verified by the coordinator session). He was answering
"Go to add the S3-only input to the gate script and then take the reading?".

**Governing documents:**
- W1 derivation AMENDMENT 6, Drive 1fnORgBkyOnowBr9gqP350xniDPZMe6icQYAYCJxO5tM.
- AMENDMENT 5, Drive 1U59VB-y99ivjRYx-zO0PaRjTaKt1_YHzIX8-Jfw0gj4.
- Both are post hoc and apply to future 'f readings only. #580's readings stay on record as read.

**Admission rule:** dr 0.025 is the admission grid. dr 0.05 is record only (A3 section 5).

## Committed before anything is built

1. **Sum-back tolerance:** the S3-only and S1-only arrays are built from the stored harmonic vectors through the
   producer's unchanged `mapped()`. Per grid, they must satisfy

       max|S3 + S1 - S| and max|PS3 + PS1 - PS| <= 1e-6 x edge_tolerance

   That is 1.943e-13 at dr 0.025 and 1.949e-13 at dr 0.05.
   - **Why this size:** the split differs from the slice only by floating-point reassociation, which is about 1e-18
     on |S| <= 1.3e-3. Anything above a millionth of the edge scale is a real build difference.
   - **If it fails:** STOP; no 'f reading is taken.
2. **k_d per grid** (A6 2.1): k_d = arccos(1 - (k dr)^2/2) / dr. The solver's (k-hat h-hat) equals k dr in engine
   units.
   - D2: 3.544330 at dr 0.025, 3.547821 at dr 0.05.
   - S3: 5.367633 at dr 0.025, 5.379812 at dr 0.05.
   - These match the fitted phase slope of the stored harmonic vectors.
3. **Operator** (A6 2.2): the gate differentiates with the engine's 4th-order centred stencil (vendor
   `cell_derivative_fourth`). On a k_d wave it reads k_d (1 - (k_d h)^4 / 30), which is -1.1e-5 / -1.7e-4 for S3.
   So k_d stands, and the windows are unchanged.
4. **Windows** (A5 2.2 rule with k_d; H0 = 0.0056689681106181795, the G0.3 value):

   | field | dr 0.025 | dr 0.05 |
   |---|---|---|
   | D2 | [20, 59.850], 17 wavelengths, 1594 cells | [20, 57.050], 16 wavelengths, 741 cells |
   | S3 | [35, 59.600], 15 wavelengths, 984 cells | [35, 59.500], 15 wavelengths, 490 cells |

   The rule is recomputed in the script from the engine H0 and must reproduce these edges.
5. **Estimator:**
   - F_st = sum(kappa_i g_i) / sum(kappa_i), with kappa = k_d - Omega H0 r and g = 4 pi r^2 (-P_st G_st).
   - P_st = P - H0 r f_r and G_st = f_r - H0 r P_st, using engine derivatives.
   - The ratio is F_st / P_tail, with P_tail = 2 pi Omega k_d A. A is the constant term of the a^2 fit over the
     'b window, fitted at 2 k_d.
   - **For S,** every one of these (a^2, P_st, G_st, the flux, P_tail) uses the S3-only arrays.
   - **For D,** the total field is used (D4 is about 2e-4 of D2).
6. **Bands** (A5 3.2, unchanged by A6): [0.999, 1.001] at dr 0.025 and [0.9985, 1.0015] at dr 0.05, for D and S.
7. **Other checks unchanged:** G0.1-G0.4, the engine entry, and 'a-'e are as predeclared (A3 5, with the 'd PD bound
   1.1331e-7). They are evaluated first, on the total fields. The A3 version of 'f is recorded beside the new one.
8. **Faults rerun alongside:**
   - uncorrected (admitted) slice;
   - v^2 piece missing (watch D: 0.992 under A5);
   - wrong-sign phase.

   Each fault's S split uses that fault's own construction, so that 'f, not the sum-back check, is what is tested.
   Each must fail 'f on the admission grid.

Paper expectations (A6, not gates): worst |'f - 1| of 1.9e-4 (D) and 1.1e-4 (S) at dr 0.025, and 1.2e-4 (D) and
4.1e-4 (S) at dr 0.05, plus beats of about 3e-4. S1 contributes nothing once S3 is separated.
