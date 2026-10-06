# 0* D-Mode Time-Shape Source-Order Result — 2026-10-06

Campaign: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37524557434
Branch: `0star-central-clock`
Resolution: N=160
Worldtube: Rmax=160
CFL: 0.0075
Cases: D0=0, Dhalf=5e-11, Dbase=1e-10, Ddouble=2e-10
All four cases completed to t=22.5 without the earlier interior radiation failure.

## Main result

The decisive diagnostic is state alignment. When the three nonzero-D trajectories are re-parameterized by the instantaneous central

D_active = 2 PD^2 + D^2,

their local geometric/time response nearly collapses onto one trajectory, even though the same D_active state is reached at substantially different coordinate times.

At D_active=0.002:
- Dhalf t=21.574, max|alpha_r|=0.026421, H_center=-0.004463, tau_rate=0.53207
- Dbase t=20.852, max|alpha_r|=0.026756, H_center=-0.004560, tau_rate=0.53192
- Ddouble t=20.114, max|alpha_r|=0.026992, H_center=-0.004354, tau_rate=0.53560

At D_active=0.004:
- Dhalf t=22.483, max|alpha_r|=0.037476, H_center=-0.009918, tau_rate=0.27054
- Dbase t=21.752, max|alpha_r|=0.037511, H_center=-0.010038, tau_rate=0.27663
- Ddouble t=21.005, max|alpha_r|=0.037671, H_center=-0.009993, tau_rate=0.28182

At D_active=0.004 the relative spread of max|alpha_r| is about 0.5%, H_center about 1.2%, and CMC kdot about 0.9%.

## Source ordering

For the D scalar, the existing stress projection gives exactly

K_D = 4 pi alpha (rho_D + p_r,D + 2 p_t,D)
    = 4 pi alpha (2 PD^2 + D^2)
    = 4 pi alpha D_active.

The direct D geometric source therefore rises with D_active but is simultaneously modulated by the falling central lapse.

For the three nonzero-D cases, the strongest positive growth in the D K-source occurs before the strongest growth in the lapse-gradient signal. Approximate derivative-peak ordering:

Dhalf:
D-source 20.963 -> matter K-source 21.038 -> alpha_r 21.248 -> central lapse-rate decline 21.413

Dbase:
D-source 20.265 -> matter K-source 20.280 -> alpha_r 20.505 -> central lapse-rate decline 20.708

Ddouble:
D-source 19.515 -> matter K-source 19.395 -> alpha_r 19.703 -> central lapse-rate decline 19.913

The matter-source peak is essentially coincident with the D-source peak, not robustly one-way after it. This favors a coupled interior response rather than a clean one-directional chain.

## Control

The D=0 control reaches t=22.5 with:
- H_center=+9.34e-4
- H_eff=5.245e-3
- lapse_min=0.9900
- max|alpha_r|=0.00127
- D_active=0

The nonzero-D cases instead develop negative H_center and strong local lapse deformation.

## Global-versus-local clock

At t=22.5 the global H_eff remains approximately +0.00524 in every case, while H_center is:
- D0: +0.000934
- Dhalf: -0.010014
- Dbase: -0.014008
- Ddouble: -0.017177

Thus the global proper-volume Hubble diagnostic can remain in expansion while the central region is already contracting. The earlier worldtube/time-shape discrepancy is therefore plausibly a consequence of a localized interior dynamical phase being viewed through a global CMC observable.

## Interpretation status

Strongly supported:
1. The interior response is D-dependent.
2. The response is organized much more tightly by instantaneous D_active than by coordinate time.
3. The same local state produces nearly the same lapse-gradient and central-H response across different initial D amplitudes.
4. The direct D source grows before the strongest lapse-gradient response.
5. The D=0 control does not reproduce the effect.

Not yet established:
- strict one-way causality;
- whether the D/matter sector closes a genuine feedback loop;
- turnaround/bounce behavior or any completed 0* cycle;
- worldtube-independent long-span evolution.

Working picture:
An unstable D-mode drives a localized stress/curvature state. The geometry/lapse responds to that state, and the changing lapse then feeds back into the D evolution through its existing time derivatives. Initial D amplitude primarily changes how quickly the system reaches the same local dynamical state, rather than selecting a different local state trajectory. This is the current best-supported meaning of the emerging "time-shape" hypothesis.

No production physics, gauge prescription, boundary condition, damping, floor, or fitted coefficient was changed by this campaign.
