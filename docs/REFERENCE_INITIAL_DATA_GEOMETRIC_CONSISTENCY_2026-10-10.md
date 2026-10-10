# Reference initial-data geometric consistency audit — 2026-10-10

**Branch:** `diag/reference-initial-data-geometric-consistency`  
**Scope:** deterministic initial-slice diagnostics only; no trajectory integration  
**Reference:** `engine/reference_pirk_unified.py`, blob `a1636a2a26670a2c4dc9ba5f23985cd53c6b737d`  
**Status:** OPEN numerical inconsistency; no physical or production patch proposed

## Question

Does the legacy radial metric-source reconstruction agree with its discrete radial derivative, and does the reference BSSN Ricci scalar agree with an independent polar-areal expression? The comparison separates the first five center cells from the rest of the grid and runs both fixed-domain refinement and fixed-spacing/domain growth.

For the metric used by this reference initializer,
`ds^2 = dr^2/B + r^2 dOmega^2`, the independent polar-coordinate identity is
`R = 2(1-B)/r^2 - 2 B'/r`. The source-balance witness is
`B + r D(B) - source(B)`. Both use the reference's stated discrete derivative for `D(B)`; the polar Ricci expression is independent of the reference BSSN Ricci function.

## Preliminary numerical replay

The 16-case control matrix uses N = 40, 80, 160, 320 with (a) fixed `r_max=40` and (b) fixed `dr=1` with `r_max=N`; each is run with S amplitude 0 and 0.01.

### Control: S amplitude = 0

The uniform-background control closes: source balance, independent polar Ricci, and both Hamiltonian evaluations are at floating-point round-off (maximum residual about `8.5e-15` in the N=80, `r_max=40` replay). This supports the algebraic H0 cancellation in the reference background but does not repair the perturbed case.

### Perturbed slice: S amplitude = 0.01

At N=40, `r_max=40`, `dr=1`:
- Source balance maximum over first five cells: `1.6599733e-4`; away from center: `1.3014123e-4`.
- BSSN-minus-polar Ricci maximum over first five cells: `4.1451664e-4`; away from center: `1.6395210e-6`.
- First-cell BSSN Hamiltonian residual: `+1.0247393e-3`; independent-polar Hamiltonian residual: `+6.1022268e-4`.
- BSSN Hamiltonian residual near r=2.5: `-3.6980660e-5`.

At N=320, fixed `r_max=40`, `dr=0.125`:
- Source balance maximum over first five cells: `4.1281510e-6`; away from center: `2.9162622e-6`.
- BSSN-minus-polar Ricci maximum over first five cells: `3.7485921e-4`; away from center: `5.8119208e-8`.
- First-cell BSSN Hamiltonian residual: `+1.0648262e-3`; independent-polar Hamiltonian residual: `+6.8996697e-4`.
- BSSN Hamiltonian residual near r=2.5: `-4.8433701e-8`.

At fixed `dr=1`, domain growth from `r_max=40` to `r_max=320` leaves the local value near r=2.5 at `-3.6980660e-5` to numerical precision. This is consistent with a local discrete residual rather than a simple outer-radius effect.

## Interpretation and guardrails

1. The H0 algebraic cancellation passes the uniform, S=0 control.
2. The source-balance and Ricci mismatch are exposed by the perturbed initial slice.
3. Off-center residuals improve substantially under fixed-domain refinement, while the center Ricci mismatch and center Hamiltonian residual persist.
4. The exact source coefficient/stencil/center handling responsible remains undiagnosed. This audit intentionally does not select or patch one.
5. The independent polar formula is an identity for the stated areal metric, not a new evolution equation.
6. The reference initializer retains its historical variable-mapping mismatch and radiation omission. These are not silently changed by this audit.
7. No expensive trajectory is involved. No claim about a bounce, turnaround, completed cycle, or physical validity is made.

## Reproduction

On the diagnostic branch:
- `python -m pytest -q tests/test_reference_initial_data_consistency.py`
- `python -m engine.reference_initial_data_consistency`

The workflow uploads `runs/reference-initial-data-geometric-consistency/report.json`. Its source commit and the artifact should be retained with the audit; the preliminary values above are baseline witnesses, not a replacement for the generated artifact.
