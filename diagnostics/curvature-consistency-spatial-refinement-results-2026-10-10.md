# Independent metric-Ricci spatial refinement — results

Date: 2026-10-10  
Branch: `diag/independent-metric-ricci-audit`  
Source workflow: [run 38096180929](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38096180929)  
Source commit: `8b528199e84f05e6bcf2876b5e0f5878f4be82a2`  
Artifact: [independent-metric-ricci-fixed-dt-spatial](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38096180929/artifacts/11685444977)  
Uploaded ZIP SHA-256: `863eb4d679661e52fcf1df10b8cc465d434c7f607de6d34ad6201b878ef02bfc`

## Question

Does the discrepancy between the vendor spatial Ricci scalar and an independent metric-derived reconstruction shrink when spatial resolution is refined while coordinate timestep is held fixed?

## Configuration

- Resolutions: (N=40,80,160)
- Fixed domain: (r_{\max}=40)
- Fixed coordinate timestep: (\Delta t=0.0025)
- Effective CFL by resolution: (0.0025, 0.005, 0.01)
- Initial matter amplitude (0.01), width (7), D amplitude (10^{-10})
- Radiation on
- Sample times: (t=0,1,2,3)
- Initial state constructed with the existing discrete-consistent initializer at each resolution
- Vendor submodule pinned via the diagnostic branch checkout

Independent curvature identity:
`R_3=-4 R_{,\ell\ell}/R+2(1-R_{,\ell}^2)/R^2`, with (R=r\sqrt{b}/X) and (d/d\ell=(X/\sqrt{a})d/dr). It does not use evolved `Lambda` and does not assume coordinate radius equals areal radius.

## Results at (t=3)

| (N) | (\Delta r) | max (|R_{vendor}-R_{metric}|) | max connection contribution | max metric-route difference |
|---:|---:|---:|---:|---:|
| 40 | 1.00 | (1.90560\times10^{-5}) | (2.52810\times10^{-5}) | (6.74370\times10^{-6}) |
| 80 | 0.50 | (3.24264\times10^{-6}) | (4.21468\times10^{-6}) | (9.72039\times10^{-7}) |
| 160 | 0.25 | (5.77301\times10^{-7}) | (7.11605\times10^{-7}) | (1.34304\times10^{-7}) |

The maximum total curvature gap is in cell 0 at all three resolutions at (t=3), at radii (0.5,0.25,0.125), respectively. Its magnitude falls by factors (5.88) (N40→N80) and (5.62) (N80→N160). The connection contribution falls by factors (6.00) and (5.92). The remaining metric-route difference falls by factors (6.94) and (7.24).

At (t=0), the maximum total gap is instead around (r\approx3.3\) and falls from (4.25226\times10^{-6}) (N40) to (2.77650\times10^{-7}) (N80) to (1.75085\times10^{-8}) (N160), factors (15.32) and (15.86).

Algebraic split closure is exactly zero in the reported discrete arithmetic. The connection-identity closure is at roundoff scale (roughly (10^{-18})); vendor connection-constraint reconstruction mismatch is zero. No diagnostic gates failed.

## Interpretation guardrails

- The discrepancy shrinks strongly under refinement; these data do not support a fixed non-vanishing curvature offset in this experiment.
- This is not a formal convergence-order result. The effective CFL varies as (N) changes to preserve fixed (\Delta t), and the metric-derived route uses the repository's finite-difference operator family.
- Curvature agreement alone does not explain a Hamiltonian-constraint drift: the curvature gap must be compared directly with the full Hamiltonian residual and its other terms on the same states.
- This diagnostic does not admit the candidate initializer or change physical interpretation.
- No production equation, gauge, source, projection, or default was changed.

## Next diagnostic

Run same-state Hamiltonian residual comparison using the vendor Ricci scalar versus the independent metric-derived Ricci scalar, keeping the extrinsic-curvature and matter terms identical. Record the center and whole-grid residuals at matching sample times and resolutions. This isolates how much of the measured Hamiltonian residual could be attributed to the curvature-route gap without applying a correction.
