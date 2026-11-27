# V1-I1 — Engine-owned homogeneous state contract (review candidate)

Baseline: `5623c11c5f477de2502d8bdf98548cc829f2df2e` (PR #614, all five prior CI checks passed).
Frozen checkpoint SHA256 `d83c71a7a4d9b50a65b4682af9f906145370ec35decb4f335afefe991eba446c`.
Static package transfer `b49cb4976911747e72283ff11f876c999587f86e`; unchanged.

## Scope and representation

This is not a V55 ProductionKernel implementation. Do not advertise BSSN/PIRK2/moving gauge capabilities or change default dispatch. V1RunSpec accepts only N=8/16, stride=1/2, radius=1, tick spacing .000125, final tick 8. An accepted-state tick must be a multiple of stride. Time is derived from the integer tick, never repeatedly accumulated. Constructor checks representation, not physical admissibility or authorization. The eventual driver alone must decide acceptance.

| Channel | Shape |
|---|---|
| p, k, ell | 3N+4 |
| J, RV | 3N+5 |
| L0, K, kh, eh, w | scalar |
| legacy geometry incoming | N |
| dark_matter, baryons, radiation | each three raw N-vectors |
| S, P_S, D, P_D, phi, P_phi | six N-vectors |
| actual augmented packed incoming gauge z | 2N+4, last entry exactly 1 |

The legacy geometry incoming field is retained unchanged in the ordered physical snapshot for provenance; rate() replaces it on its internal copy from z and bR. It is **not** the actual gauge. Do not remove it or feed it as z. The state owns immutable float64 backing bytes; the exported kernel input is a detached ordered copy with original scalar values, never a live view. No implicit conversion from float32 or invalid shape, no floors, no coefficient edits. Grid assembly is not accepted state; it must be constructed deterministically by the future adapter from the run spec and checked for shape and resource identity.

## Evidence (local, isolated transferred source; not a full checkout)

Python 3.13.5 / NumPy 2.3.5 / SciPy 1.17.0, single-thread OPENBLAS and OMP; pytest 9.0.2. Local targeted suite: 11 passed with V1_FROZEN_REPLAY=1; N8/N16 direct vs round-tripped initial records and frozen JSON exactly equal. Mutation isolation and invalid scope/tick/shape/nonfinite/gauge controls pass. Separately read 52 valid frozen formal records (4 initial, 24 start, 24 endpoint) from the verified checkpoint: all physical input digests and packed gauge bytes round-trip exactly, using incoming z for initial/start and outgoing z for accepted endpoints. No midpoint records treated as accepted state. No new evolution or reference integration executed. The tests normally compare direct vs round-trip in the same environment; the dedicated pinned workflow also requests exact archived comparison. Hosted results for this candidate pending.

Local blob pins:
- engine/native_v1_contract.py: Git `10aefd6d18cef35fb1fb06238065829351c87534`, SHA256 `215a9c56d04b12858a2eaf75319ac66c43acbc8f30daac0aeee0997f2316e27f`.
- tests/test_native_v1_contract.py: Git `a6530f16200192977a31f03dfbd0047c562fd12a`, SHA256 `ce266d1917fa2337bbfed63822ba69f9581019a914b74849563f96ef80a1d98d`.

## I2/I3 seam requirements (not implemented)

- Build a research-only explicit engine adapter, not a V55 ProductionKernel. Preserve the frozen package and its resource files A8_certificate.json, basis.json, moments.json. Do not expose inherited diagnostic main/run helpers as engine APIs.
- Initial evaluation uses declared compatible z; its outgoing gauge is diagnostic, not a replacement for initial incoming z. Use original accepted z for both start and midpoint. With F0 and Fm evaluated on the proper physical and gauge slices: xm=x+hF0/2, xp=x+hFm; endpoint incoming gauge=2zm-z, endpoint rate yields zp. Only after endpoint guard and finite-step work ledger pass may the driver commit xp,zp. Midpoint is a transient candidate, not an accepted state.
- Explicit boundary inputs only: geometry boundary g=0 and scalar flux alpha(Y)*Y*J(Y)*sqrt(s2(Y))*(P_bg(t)-P_interior(Y)), on the same pre-solved lapse and original incoming z. Reject implicit scalar_boundary=None fallback. Ensure same-incoming pre-solve/rate outgoing gauge equality and immutability. Do not duplicate a competing guard implementation; preserve original attempt record schema and stop ordering. Reference DOP853 dense solution validation, both original tolerance pairs and 17 scheduled samples, mandatory before a separately authorized run; bind reference identity to run metadata.
- The frozen guard thresholds H/M <=1e-8, physical homogeneous departure <=1e-7, positive lapse, matter branch/cone, finite rates/ledgers remain unchanged. No retries or resets on failure. Record attempted and failed stages, keep accepted state unmodified, do not silently restart or promote a failed candidate.
- I4 restart/checkpoint contract must include spec, accepted tick, all physical bytes, actual accepted z, reference identity, source/environment pins and trace/ledger identity, with atomic persistence and fail-closed readback. No restart claim in I1.
- No general CFL or coupled stability claim, no inhomogeneous input or larger domains, no production or admission changes.

The existing historical Ricci source guard must fail closed on this new engine path until separately reviewed and explicitly authorized by exact blob. Do not extend the prior 18-file authorization or use a wildcard. Draft PR only. A8 unchanged, Stage B held, scientific admission false.
