# Build Status — 2026-10-06

## Current state

The canonical repository contains a unified production-facing state graph:

1. corrected archive COSMOS background equations;
2. local GEAR-03 S/D sector with GEAR-63/65 barrier;
3. scalar S/D/COSMOS evolution and stress projections;
4. conservative spherical dark-matter, baryon, and radiation matter sectors;
5. locked beta DM/scalar four-force exchange;
6. one total Einstein stress-energy assembly;
7. pinned reference-metric BSSN/PIRK numerical geometry;
8. synchronized geometry/scalar/matter evolution state;
9. dynamically derived proper time, H_eff, e-folds, trapping, and Misner-Sharp observables;
10. cycle/handoff ledgers without a synthetic lambda_cycle;
11. checkpoint and long-horizon multi-resolution campaign infrastructure;
12. an isolated experimental central-proper-time CMC branch for time-interface testing.

## Explicit exclusions

- Bianchi-I shear remains outside the exact spherical local source.
- No D-to-matter identification law.
- No shell EOS/action.
- No phenomenological interface source.
- No fitted coefficient.
- No invented clock conversion.
- No bounce/reset/branch flip/lapse floor/ejection threshold/physical stop.

## 0* status

0* is being used as a measurement and scoring layer, not as a source of model terms or gauge prescriptions.

Authoritative criteria:

- 0*_eff iff t and b and (w >= 1)
- 0*_smooth iff t and b and (w > 1)
- 0*_global iff 0*_eff and G

Turnaround t is H=0 with dot H<0; bounce b is H=0 with dot H>0. Contraction is subsequently parameterized in N=ln(a). The diagnostic must use the prescribed invariant shear test rather than manufacturing rho_tot.

## Time-interface investigation

The production CMC solve retains its archive-supported alpha(R)=1 normalization. An experimental branch tests the archive relation d tau = alpha(0) dt as a full coordinate transformation.

The first lapse-only implementation failed before turnaround at an inner radiation cell. The result is being treated as an interface-consistency witness, not as evidence against the central proper-time coordinate itself.

The current branch therefore tests whether all evolution components that interpret dt are transformed consistently with the new time parameter. No physical equation, source, coefficient, boundary condition, or bounce rule is being added.

## Validation state

Substantial numerical diagnostics have already isolated the earlier outer-worldtube failure as a boundary-layer/gauge interaction and established that the effect changes with worldtube radius. The present engineering gate is the time-interface consistency test.

Do not merge the experimental gauge branch into production until:

- the complete repository test suite passes;
- the central-clock transformation is internally consistent;
- the full numerical campaign is rerun;
- turnaround is located in the correct physical/cosmic time variable;
- contraction is analyzed in N=ln(a);
- shear/anti-BKL and bounce conditions are scored through 0*;
- resolution and worldtube dependence are understood.

## Standing rule

Run the complete machine, expose the actual defect, repair only that layer, and rerun. No hand tuning, arbitrary damping/floors, fitted feedback, or invented boundary physics.

No completed multi-cycle or physical-bounce claim is made until invariant, resolution-independent evidence supports it.

## 2026-10-07 center/clock diagnostic gate

The widened strong-D constraint scan (N=160/240/320/480, r_max=N, Δr≈1) shows stable proper-time suppression and slowly decreasing central Hamiltonian residual, but it is a domain-size/fixed-Δr scan rather than formal fixed-domain convergence.

The apparent trapping root appears immediately and at a nearly fixed physical radius as the outer domain is enlarged; it is therefore treated as a domain/pre-existing-geometry diagnostic candidate until explicitly tested.

Next diagnostic sequence:
1. Center-regularity repair A/B at N=160, r_max=80 for D=0 and D=1e-4, repair ON/OFF.
2. On those same runs, report radial τ(r) and α(r) profiles at common times.
3. Only then spend N=320 fixed-r_max=80 runs.

The A/B campaign is diagnostic-only and changes no physical evolution equation, source term, fitted coefficient, or boundary physics.



## 2026-10-10 engineering update — isolated source/CMC branch

This is a dated addendum to the 2026-10-06 record above. It does not replace that historical status or merge the current branch into production.

- Active branch: `physics/spatial-beta-covariant-source-closure`.
- Current tested head at this update: `1ad4b6219bfd6899d5f064cd90f19f6a2dee32c3`.
- `main` remains `48f06d09a5b88b831ac86f6bdd3d9ecc94f092f6`; no merge was made.
- [CI #572](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068034703) passes **124 tests in 60.02 s**. This commit adds fail-fast CMC spatial-metric validation in both the projection and full lapse-solver entry point. [The preceding #571 failure](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38067950792) is preserved and documents a fixed stale variable reference.
- Also verified: [#568, multi-configuration CMC residual](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38067318842); [#569, actual predictor/final CMC residual](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38067488971); [#570, invalid-metric projection gate](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38067741630).
- The current CI runs tested code, not long physical trajectories; optional campaign jobs were skipped. The radiation outer-boundary failure still needs reproduction against a pinned current-head run.
- Open questions remain: full-domain constraint/convergence behavior, radiation outer-boundary realizability, reciprocal local-COSMOS coupling, D-to-matter identification, and global shear closure. No turnaround, bounce, completed cycle, or validated cosmology is claimed.

Standing rule unchanged: run the full machine, expose the actual defect, repair only the affected layer, and rerun. Keep code diagnostics, physics changes, and long-run campaigns independently identifiable.

## 2026-10-10 radiation predictor trace — captured, not repaired

- Source baseline: `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`; trigger commit: `c8a9967126ba07cd12e4da9c86b7066a703d3419`; branch remains isolated. No production physics source changed in the instrumentation commits.
- [Trace run #7](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938182) completed its diagnostic capture; [CI #576](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938171) passed 127 tests. Artifact [11676486591](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938182/artifacts/11676486591).
- Captured-state metric-cone regression added in `tests/test_radiation_predictor_metric_cone_attribution.py`; [CI #579](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38070237509) passed 128 tests in 35.48 s. The diagnostic trace job was correctly skipped on this non-launch commit.
- N=80, r_max=80, dr=1, dt=CFL=0.03, radiation ON, D=1e-10. Last accepted t=46.98, τ=22.0161154; attempted predictor t=47.01 fails radiation `E>=|S|` at cell i=79 (r=79.5), ratio 1.0000569941.
- Ordered cone budget at i=79: accepted `C=+4.6971e-7`; flux-updated on accepted metric `+1.4812e-6`; source-updated on accepted metric `+2.2171e-6`; predictor metric contribution `−2.2340e-6`, giving `C=−1.6900e-8`. The predictor metric changes `gamma_rr_inv` upward by 1.5239% at this cell; flux/source terms in this ordered budget improve its margin.
- Detailed data and artifact hashes: [`docs/RADIATION_PREDICTOR_TRACE_2026-10-10.md`](RADIATION_PREDICTOR_TRACE_2026-10-10.md).
- Code-order refinement: the failing first CMC lapse solve sees the explicit Euler geometry predictor after `enforce_algebraic_regularity`; it fails before the later `_apply_outer_light_boundary` call. Earlier boundary effects on the accepted state/RHS remain possible, but that later boundary call is not the direct failing operation.
- Next: record raw explicit RHS (`a`, `b`, `X`), unprojected Euler predictor, and post-projection values at cells 75–79; reconstruct `Δgamma_rr_inv`. Do not add clipping, damping, floors, or a speculative boundary law. The upstream defect remains open.