# V5.5 radiation-sector ablation — 2026-10-06

## Purpose

Separate the radiation matter-sector admissibility boundary from the underlying
strong-field geometry/finite-radius continuation problem without changing any
production equation, parameter, lapse floor, clipping rule, or boundary source.

Both runs use the same stage-aware True-CMC kernel, N=80, Rmax=80, CFL=0.06,
with radiation either included or excluded at initialization.

## Result

| case | turnaround | first failure | proper time at failure | trapped roots |
|---|---:|---|---:|---:|
| radiation ON | 1 | radiation inversion at outer cell | 22.02348 | 2 |
| radiation OFF | 1 | nonpositive conformal metric ratio | 22.71259 | 4 |

### Radiation ON

The outermost radiation state reaches the characteristic limit:

- E = 4.3817515e-8
- |S| = 4.3854674e-8
- |S|/E = 1.000848

The failure is confined to the outer transport layer at the same strong-field
continuation stage already seen in the multi-resolution runs.

### Radiation OFF

Removing radiation does **not** restore healthy continuation. The evolution
continues farther in proper time, then fails with:

**FloatingPointError: nonpositive conformal metric ratio**

At the final valid slice:

- H_eff = -0.71184
- Hamiltonian outer L2 = 40.8331
- momentum outer L2 = 4.49724
- connection max = 83.9888
- Misner-Sharp mass at the sampled worldtube = 140.712
- trapping indicator is deeply negative
- two outer/inner trapped roots are reported.

## Adjudication

Radiation is therefore **not the fundamental strong-field continuation blocker**.

Radiation changes where the finite-radius computation first becomes invalid:
the radiation-on branch reaches its matter-cone boundary first. The radiation-off
branch exposes the underlying geometric deterioration later.

This supports the archived disposition that the primary unresolved layer is the
strong-field spacetime representation / finite-radius continuation, not an
unmodeled radiation feedback term.

## Next gate

Do not add a radiation clip, lapse floor, fitted transport coefficient, or
phenomenological boundary source.

Proceed with the archived reference-metric / native-PIRK strong-field route:

1. use the pinned reference BSSN/PIRK operators rather than another hand-reduced
   tensor implementation;
2. retain the frozen V5.5 stress-energy and scalar source semantics;
3. reproduce the t~22.576 trapped event;
4. continue through the trapped regime;
5. compare invariant witnesses and resolution convergence before reconnecting the
   full long bidirectional campaign.
