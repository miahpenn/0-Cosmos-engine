# W1 Stage 0 origin trace

This diagnostic replays the corrected ring (A) and background-only control (B) to the common coordinate time \(t=0.015\). The grids are \(dr=0.05\) and \(0.025\), with the open edge, eight-node areal initializer and KO off. Each fine run also uses half its nominal timestep. The trace samples the six innermost and three outermost cells through the first PIRK2 step, then at each accepted slice. Cell 0 has coordinate radius \(dr/2\).

## First step

The initial lapse and gauge refresh leaves the sampled \(H\) and \(M_r\) unchanged. On the fine ring run, the explicit predictor changes the cell-0 scalar momentum source \(8\pi j_r^{\mathrm{scalar}}\) by \(-4.1480\times10^{-10}\), with \(A_a\) and \(K\) still at their initial values. The primary predictor changes \(D_r A_a+3A_a/r\) by \(-4.3752\times10^{-10}\), bringing the geometry and source back into partial cancellation. The corrector reduces the remaining difference again. The values below show the signed cell-0 momentum on that run.

| Slice | \(M_r\), cell 0 |
|---|---:|
| After initial lapse refresh | \(-2.0720\times10^{-16}\) |
| Explicit predictor | \(+4.1479\times10^{-10}\) |
| Primary predictor, before outer write | \(-2.2623\times10^{-11}\) |
| Corrected, before outer write | \(-9.9839\times10^{-13}\) |
| First accepted state | \(-9.9839\times10^{-13}\) |

Neither outer write changes the sampled centre values. The fine background-only run has a separate outer-cell signal. At the corrector boundary write its outermost momentum residual changes from \(-2.3118\times10^{-14}\) to \(-9.1029\times10^{-14}\). Its cell-0 reading stays at \(-7.5000\times10^{-15}\) through that write. At \(t=0.015\), its largest momentum residual is at cell 0, so the initial outer maximum alone does not explain the later background signal.

## Matched-time fine-grid comparison

| Ring, cell 0 | Nominal \(\Delta t=7.5\times10^{-4}\) | Half \(\Delta t=3.75\times10^{-4}\) |
|---|---:|---:|
| \(H\) at \(t=0.0075\) | \(+5.7594\times10^{-11}\) | \(+4.5659\times10^{-12}\) |
| \(H\) at \(t=0.015\) | \(+1.4096\times10^{-10}\) | \(+2.7064\times10^{-11}\) |
| \(M_r\) at \(t=0.00075\) | \(-9.9839\times10^{-13}\) | \(+1.3268\times10^{-13}\) |
| \(M_r\) at \(t=0.015\) | \(+1.8413\times10^{-12}\) | \(+1.4796\times10^{-12}\) |

At \(t=0.015\), the fine ring's Hamiltonian residual reaches its maximum at cell 0 with either timestep. Halving the timestep reduces that maximum by \(5.21\times\). The coarse nominal run has \(H=+1.4312\times10^{-10}\) at cell 0, close to the fine nominal result. The same-grid budget shows where the timestep difference appears at this endpoint. The change in the cell-0 Ricci scalar from initialization is \(+4.3920\times10^{-8}\) at nominal step and \(+4.3806\times10^{-8}\) at half step. Changes in \(16\pi\rho\) and the kinematic terms agree to four displayed significant digits between the two timesteps. This identifies the Ricci contribution on two evolved states. An independent curvature evaluation of the same states is needed to distinguish a change in the evolved metric from sensitivity of the Ricci evaluator.

The ring's \(A_a/r^2\) is nearly flat across fine cells 0–2 on the first accepted slice and at the endpoint. At \(t=0.015\), these values are \((-1.8440,-1.8447,-1.8412)\times10^{-7}\) at the nominal timestep and \((-1.8448,-1.8442,-1.8410)\times10^{-7}\) at half step. This sample does not establish origin regularity as \(r\to0\). The late-time cell-0 momentum changes weakly with halving, even as the fine ring's grid-wide maximum falls from \(1.1853\times10^{-11}\) to \(3.5574\times10^{-12}\), both at cell 1. These mixed ratios do not establish a timestep convergence order.

The next diagnostic can use copies of the same accepted states to compare the production and direct areal curvature residuals, with their source and kinematic terms held fixed. The Stage 0 result remains STOP, including the independent runtime limit.

## Verification

All six runs reached the matched endpoint, completing 140 accepted steps in total without retries. The nominal replays reproduced all 128 stored Hamiltonian and momentum max-norm values through \(t=0.015\) exactly. The full pytest suite finished with 1,202 passed, two skipped and two expected failures.

## Reproduction

Engine and smoke inputs come from PR #584 head `05c2a09aae5fc1224f110c055c0c5b4b3acf8b21`. The vendor pin is `d6052d605673ce9d82cdc99b0df79855fa2ea215` and the corrected NPZ SHA-256 is `339652a2a7e1379b19812893e43bafcf9ab0211422f2b8685c951107b863a295`. The observer checks all three and records its own SHA-256. The runtime used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and OpenBLAS 0.3.30 (SkylakeX, one thread).

From a clean worktree at the smoke head with the pinned vendor checked out, copy `research/w1_fring_smoke/origin_trace.py` from the research revision. Run the following from the root of that worktree.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for arm in A B; do
    for dr in .05 .025; do
        python -m research.w1_fring_smoke.origin_trace "$dr" "$arm"
    done
    python -m research.w1_fring_smoke.origin_trace .025 "$arm" --half
done
```

`M_Aa` means the even-parity radial derivative of \(A_a\) plus (3A_a/r). `M_K` is the \(-2/3\) multiple of the even-parity radial derivative of \(K\). `M_other_geometry` is the remainder of the geometry operator after subtracting those two terms. The signed source is subtracted in both constraints. Negative sample indices count from the outer edge. No production or acceptance configuration files are part of this changeset.
