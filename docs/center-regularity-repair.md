# Center-regularity repair and validation plan

## Finding that triggered the repair

The completed pre-repair resolution campaign showed a constraint maximum that moved toward the spherical center as resolution increased:

- N=80: Hamiltonian max near r=1.0, about 1.58e-2 at t=30.
- N=120: Hamiltonian max near r=2/3, about 3.41e-2.
- N=160: Hamiltonian max near r=1/2, about 5.97e-2.

At the same time, the CMC elliptic residual remained at roughly 1e-9 and the lapse minimum became strongly resolution-sensitive. The location pattern is tied to the first retained cell, not to a fixed physical radius.

Inspection of the pinned spherical BSSN kernel exposed an unused numerical regularity projection, `enforce_algebraic_regularity`, which explicitly enforces the exact center identities:

- a*b^2 = 1 everywhere;
- a/b = 1 + O(r^2) at the first cell;
- Aa = O(r^2) at the center;
- odd connection/shift variables = O(r).

The production branch was evolving these fields without reapplying that projection.

## Repair

The production kernel now reapplies the pinned regularity projection at initialization and at every geometry stage needed by the PIRK/CMC sequence. The projection is algebraic/numerical only; it adds no source, coefficient, damping, lapse floor, boundary condition, bounce rule, or matter-identification law.

## Validation gates

1. Re-run the N=80/120/160 Dbase resolution campaign to see whether the inward-moving Hamiltonian spike and lapse sensitivity collapse.
2. Re-run the Dhalf/Dbase/Ddouble synchronized long campaign at N=160 to test whether the D-phase collapse survives the repair.
3. Compare CFL=0.0075 against 0.00375 at N=160 and compare same-dr Rmax=80 against Rmax=160.
4. Only after these gates pass should a longer turnaround/bounce search be treated as physically interpretable.

The key rule remains: repair the numerical layer first; do not compensate for a constraint defect with physical tuning.
