# Accepted-state source-mapping ratio witness, cells 0–5 — 2026-10-25

**Request:** Lumen's Slack handoff (message timestamp `1792942493.056289`, thread rooted at `1792887830.446859`) requested the full-profile per-cell arrays for (dot M^{matter}-8pidot J) and the provisional (-8pialpha Q_r), N=60/120/240 at accepted t=0.60. The arrays and this script/path were shared in the Slack thread.

**Source run:** [40159027760](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40159027760), observer head `ddc560964067c0e2d1569e401ff6cc92fe969adf). Production anchor `4aec88d0e3ff028e0a65a8ac36b30948625df792`; vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. The source/vendor/protected-file guards passed and accepted-state artifacts are manifest-hash verified.

## Definition and extraction path

In `tools/audit_hybrid_origin_integral_patch_newton.py`:
- lines 520–539 form the matter contribution to the geometric momentum-constraint tangent from the explicit (A_a) and (K) matter source rates;
- lines 544–547 define `matter_source_total_composite = Mgeom_dot_matter_injection - 8*pi*j_cov_dot`;
- lines 548–555 define `stress_divergence_source = -8*pi*alpha*Q_r`.

The values below are extracted from `profiles.npz` fields `N{N}_full_continuous_profile_eps1e-05_propagation_actual_matter_source_composite`, `..._propagation_stress_divergence_source`, and `..._propagation_Qr_total_stress_divergence`. The run's directional (epsilon=2\times10^{-5}) arrays agree closely with the listed (epsilon=10^{-5}) values. Cells are the native coordinate-grid cell centres, not matched fixed areal radii.

## Per-cell arrays

### N=60

```
Mdot_matter - 8*pi*Jdot:
[+2.242612792715e-05, +5.830097790691e-05, +7.295880885762e-05,
 +6.644673984452e-05, +4.815360817734e-05, +2.874894243519e-05]

-8*pi*alpha*Qr:
[+9.539374709231e-08, +2.109797260556e-07, +1.608133412568e-07,
 +4.251358282756e-10, -1.231790576072e-07, -1.316586195663e-07]

Qr:
[-3.800343798666e-09, -8.404193616703e-09, -6.404699016653e-09,
 -1.692847738438e-11, +4.904005928753e-09, +5.240894635133e-09]

composite / (-8*pi*alpha*Qr):
[+2.350901e+02, +2.763345e+02, +4.536863e+02, +1.562953e+05,
 -3.909237e+02, -2.183597e+02]
```

### N=120

```
Mdot_matter - 8*pi*Jdot:
[+1.141877555821e-05, +3.304158139479e-05, +5.123825900548e-05,
 +6.438487115194e-05, +7.167584348096e-05, +7.317723347867e-05]

-8*pi*alpha*Qr:
[+3.308476028663e-09, +9.254477804813e-09, +1.329189605637e-08,
 +1.453433935053e-08, +1.271487741127e-08, +8.347190467794e-09]

Qr:
[-1.318021739449e-10, -3.686670794894e-10, -5.294759438079e-10,
 -5.789253028367e-10, -5.064086663093e-10, -3.324204025758e-10]

composite / (-8*pi*alpha*Qr):
[+3.451370e+03, +3.570335e+03, +3.854850e+03, +4.429845e+03,
 +5.637164e+03, +8.766690e+03]
```

### N=240

```
Mdot_matter - 8*pi*Jdot:
[+5.731697718258e-06, +1.703699385888e-05, +2.788518897619e-05,
 +3.799813646536e-05, +4.712437822211e-05, +5.505842624889e-05]

-8*pi*alpha*Qr:
[+1.060486659553e-10, +3.132020627462e-10, +5.039111102551e-10,
 +6.685747077679e-10, +7.983783342843e-10, +8.858049958655e-10]

Qr:
[-4.224718447976e-12, -1.247711352937e-11, -2.007416192033e-11,
 -2.663327509719e-11, -3.180327832890e-11, -3.528478057203e-11]

composite / (-8*pi*alpha*Qr):
[+5.404781e+04, +5.439617e+04, +5.533752e+04, +5.683454e+04,
 +5.902512e+04, +6.215637e+04]
```

## Interpretation guard

- The ratio is not a resolution-independent factor such as (4pi) or (8pi); however, a ratio by itself is ill-conditioned where (Q_r) is near zero. In N=60, (Q_r) crosses sign between cells 3 and 4, producing large/sign-changing ratios there. Compare the signed **absolute arrays** first.
- The matter composite is dominated by the S sector; the per-sector decomposition is already stored in the same profile NPZ.
- The subsequent Brown source-bridge report establishes that (S_{act}-S_Q) is almost exactly balanced by the matter-dependent difference (B_{tot}-B_{geom}), with signed reassembly closure (O(10^{-20})). Thus this ratio is a diagnostic of the provisional split, not independent proof that a physical source term is missing. See [Brown bridge report](brown_matter_source_bridge_absolute_residual_2026-10-25.md).
- No source law or production equation is changed and `scientific_admission=false`.

## Artifact hashes

| N | Artifact ID | SHA-256 of `profiles.npz` |
|---:|---:|---|
| 60 | 12584544612 | `41c9a4c4c1c26eec6c51ad662f371152f1e46656c27b90df2121e41313a455c5` |
| 120 | 12585144988 | `623931dc52cfd0923b1b96c737f48a7ec9059808fd0e587450cc49d6bc2108ab` |
| 240 | 12585018461 | `48ee3b91a9ffb428bd2d253f65bb1632281664967b9a31a3a4a5b0efb46f08ea` |
