# W1 diagnostic spatial refinement at matched timesteps

Compare the same corrected W1 ring at \(dr=0.025\) and \(dr=0.0125\), separately at \(\Delta t=0.000375\) and \(0.0001875\). Those are the half and quarter timesteps from the earlier fine-grid work. The four trajectories reach \(t=0.015\) in 40 and 80 steps per grid. They use the open edge, eight-node areal initializer and KO off. The run is a research comparison. It does not set a new entry gate or change A6.

## Fine-grid input

The harmonic-balance equations were solved on a fresh 9,532-node toy grid aligned to the 4,800-cell engine grid. The prior \(dr=0.025\) member was used to construct a Newton seed. The evolved fine input is the newly solved member after the same static-time phase map, rather than a resampling of the archived fields. The pinned W1 producer reproduces the archived corrected \(dr=0.025\) fields bit for bit. The fine solve converges in two iterations, with maximum residual \(6.8541\times10^{-12}\) and maximum node discrepancy \(7.1054\times10^{-15}\). The diagnostic uses the archived \(1.943\times10^{-7}\) edge threshold for both supplied inputs. Both pass engine entry.

## Readbacks

Production \(H\) and \(M\) include their matter sources. The seven- and nine-point direct metrics use the same kinematic term and matter density as production. The tracefree witness uses the integrated stage-sampled \(\alpha A_a\). The write witness subtracts the accumulated signed binary64 storage errors from the frozen \(a,b,X\) profiles, holding computed stage rates fixed. Every corrector write is reconstructed bit for bit from the actual native stage increments. Changing the storage representation in an evolved trajectory could also change those rates.

Cell 0 is shown separately at its own radius on each grid. A second comparison samples 17 fixed physical areal radii \(R=0.05\) to \(0.25\) using cubic interpolation for the readbacks. The JSON also records the difference between linear and cubic interpolation as a sensitivity check. The evolution itself uses each grid's native fields. Signed residuals are archived at the initial slice and every \(0.00375\) through \(0.015\). Both the raw endpoints and the changes from each grid's own initial residual matter because the initial fine-grid \(H\) already differs from the reference.

The initial production \(H\) maxima are \(7.1373\times10^{-12}\) at \(dr=0.025\) and \(1.6469\times10^{-11}\) at \(dr=0.0125\). The signed cell-0 initial values are \(+7.1373\times10^{-12}\) and \(-1.3666\times10^{-11}\). The initial global \(M\) maxima are \(4.00\times10^{-16}\) and \(8.89\times10^{-16}\), respectively.

## Endpoint comparison

At \(t=0.015\), all table entries are in engine units. Band maxima are over the same 17 physical radii. Grid maxima use native cells.

| Matched \(\Delta t\) | \(dr\) | Signed cell-0 \(H\) | Peak \(\lvert H\rvert\) | Band \(\max\lvert M\rvert\) | Peak \(\lvert M\rvert\) |
|---|---|---|---|---|---|
| 0.000375 | 0.025 | \(+2.7064\times10^{-11}\) | \(2.7064\times10^{-11}\) | \(2.3333\times10^{-12}\) | \(3.5574\times10^{-12}\) |
| 0.000375 | 0.0125 | \(-1.6789\times10^{-11}\) | \(8.1575\times10^{-11}\) | \(4.6901\times10^{-12}\) | \(2.8602\times10^{-11}\) |
| 0.0001875 | 0.025 | \(-2.1073\times10^{-11}\) | \(2.5262\times10^{-11}\) | \(1.3087\times10^{-12}\) | \(6.8067\times10^{-12}\) |
| 0.0001875 | 0.0125 | \(-8.4649\times10^{-11}\) | \(1.4908\times10^{-10}\) | \(4.6738\times10^{-12}\) | \(8.2393\times10^{-11}\) |

Both refinements raise the global \(H\) and \(M\) maxima. On the finer grid the worst \(H\) is at the outer edge and the worst \(M\) at cell 0. The endpoint first-cell areal radii are \(0.01250106\) and \(0.00625053\). The momentum increase persists at fixed radii, where the band maximum rises by \(2.01\times\) and \(3.57\times\) at the two matched timesteps. The interpolation sensitivity is below the increase in both pairs. The centre \(H\) changes from its own initial value by \(+1.993\times10^{-11}\) and \(-2.821\times10^{-11}\) on \(dr=0.025\), and by \(-3.124\times10^{-12}\) and \(-7.098\times10^{-11}\) on \(dr=0.0125\).

The maximum cross-grid differences in \(S,D\) on the physical band stay near their initial magnitudes. They are \(3.29\times10^{-8}\) and \(2.62\times10^{-8}\) at initialization and \(3.28\times10^{-8}\) and \(2.67\times10^{-8}\) at the endpoint for either matched timestep. These are comparisons of newly solved same-branch initial data and their evolutions. They do not measure spatial order with identical initial arrays.

## Curvature write and tracefree witnesses

Here each difference is larger timestep minus smaller timestep on the indicated grid, using nine-point direct curvature at cell 0.

| \(dr\) | Direct \(H\) difference | Tracefree curvature witness | Accumulated write witness |
|---|---|---|---|
| 0.025 | \(+3.6945\times10^{-11}\) | \(+1.8598\times10^{-11}\) | \(+1.8513\times10^{-11}\) |
| 0.0125 | \(+1.8719\times10^{-10}\) | \(+8.4670\times10^{-11}\) | \(+1.0395\times10^{-10}\) |

The centre write witness grows about \(5.6\times\) and exceeds the tracefree witness on the finer grid. At fixed physical radii the respective maximum magnitudes of the timestep differences in the tracefree and write witnesses are \(2.18\times10^{-12}\) and \(1.22\times10^{-11}\) on \(dr=0.025\), versus \(1.39\times10^{-11}\) and \(1.73\times10^{-11}\) on \(dr=0.0125\). These are finite changes through absolute field substitutions and are not a complete additive Hamiltonian budget. The seven-point witness gives the same qualitative result. The matched comparisons show no spatial improvement in the two constraint maxima for this window. A convergence order and the cause of the increased momentum residual remain open.

## Verification and reproduction

All 240 steps reached their prescribed endpoints without retries. Both reference endpoint production \(H\) values match the archived ledger exactly. All 37 targeted tests passed. The engine and A6 settings were preserved.

Use the pinned W1 ringkit v10 archive (archive SHA-256 b498ef42e68e778d536b1be0b559df9d65eb758348c3006b9b3a15dbcd3ced91) and the aligned producer directory with the hashes in input.json. The producer members.json SHA-256 is 697056dc28d6871c11ab95858dd301bb0bc124b5aa96739b5d8a5e83870ab4f0. The fine NPZ SHA-256 is 0cc47fb17b149e36c6ac6c45abb58ac240ed632912549f538d9400a48d46c39b. Runs used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and one OpenBLAS 0.3.30 SkylakeX thread.

    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    python -m research.w1_fring_smoke.prepare_spatial "$RINGKIT_V10" "$PRODUCER" results/research/w1_fring_spatial
    output=$(mktemp -d)
    for dr in 0.025 0.0125; do
        for dt in 0.000375 0.0001875; do
            python -m research.w1_fring_smoke.spatial_replay "$dr" "$dt" "$output/$dr-$dt.npz"
        done
    done
    python -m research.w1_fring_smoke.spatial_compare "$output"
