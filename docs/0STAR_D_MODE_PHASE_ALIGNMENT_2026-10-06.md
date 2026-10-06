# D-mode phase alignment — 2026-10-06

## Experiment

Rephase the existing N=160, Rmax=160 source-order trajectories by the instantaneous central

[
D_{active}=2P_D^2+D^2
]

rather than by coordinate time. No engine physics, gauge prescription, boundary condition, coefficient, damping, or timestep rule was changed.

The common-state crossings were evaluated at D_active = 0.001, 0.002, 0.003, and 0.004 for Dhalf = 5e-11, Dbase = 1e-10, and Ddouble = 2e-10.

## Main result

The three trajectories are almost a single trajectory with different phase offsets.

At every tested D_active level, the inferred coordinate-time offsets are nearly constant:

- Dhalf − Dbase: mean 0.72675, standard deviation 0.00346
- Dbase − Ddouble: mean 0.74035, standard deviation 0.00458
- Dhalf − Ddouble: mean 1.46710, standard deviation 0.00722

The proper-time offsets are even tighter:

- Dhalf − Dbase: mean 0.71803, standard deviation 0.00106
- Dbase − Ddouble: mean 0.71802, standard deviation 0.00074
- Dhalf − Ddouble: mean 1.43605, standard deviation 0.00179

This is the key signature: changing the initial D amplitude primarily changes how quickly the system reaches the same local state.

## Common D_active = 0.004 slice

Crossing times:

| case | t | tau | D | P_D |
|---|---:|---:|---:|---:|
| Dhalf | 22.48388 | 21.42907 | 0.0371858 | 0.0361666 |
| Dbase | 21.75194 | 20.71228 | 0.0372025 | 0.0361652 |
| Ddouble | 21.00467 | 19.99490 | 0.0372204 | 0.0361547 |

The local spatial geometry is very tightly aligned at this common D state. Over r <= 20:

- a profile RMS relative spread: 0.19%
- b profile RMS relative spread: 0.095%
- X profile RMS relative spread: 0.75%
- K profile RMS relative spread: 0.55%

The largest local K spread is about 1.2%.

The lapse itself is less tightly collapsed: RMS relative spread about 3.2%, maximum about 4.0%. That distinction matters. The spatial geometry/extrinsic geometry is much closer to a common state than the coordinate lapse.

## Interpretation

This is stronger than simply observing that the three runs eventually pass through similar values.

The data support a phase-coordinate picture:

[
D_{active} longmapsto 	ext{local dynamical state}
]

while the initial amplitude determines approximately the phase shift needed to reach that state.

The near-constant time offsets across four different D_active values are especially important. They are difficult to explain as a single accidental coincidence at one checkpoint.

The result does **not** yet prove a one-way causal chain. The source-order campaign already showed that D-source growth precedes the strongest lapse-gradient response, but the matter K-source peak is nearly coincident with the D-source peak in some cases. The system remains a coupled dynamical system.

## Next clean test

The next useful analysis is to replot/recompare the three histories using D_active as the independent variable and test whether the state-vector components collapse continuously, rather than only at four crossing points.

The most informative state vector is:

[
(D,P_D,K,alpha,alpha_r,dot K_{m CMC},
	ext{sector K-sources})
]

with the spatial geometry ((a,b,X)) kept as the invariant-style geometry check.

No production change should be made from this result alone.
