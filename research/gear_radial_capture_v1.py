"""Research-only capture runner for the archived face-consistent GEAR solver.

Source equation/parameter anchor: Library file
GEAR_radial_validated_face_consistent_0125_2026-10-03.py, SHA-256
0f2fc685f90b901275b33670c24aee3d2f1dbe8712d6e8de95d426ad2e1f5f4b.
This is an instrumented copy, not a byte-identical canonical source. The
additional work records profiles, RK4-consistent central proper time, and a
separate geometric mass/face-flux ledger. It does not add evolution terms.

Example: python -m research.gear_radial_capture_v1 --ngrid 300 --final-time 10
Scientific admission is explicitly not claimed.
"""

import numpy as np

from research.gear_radial_derived_observables_v1 import (
    CORRECTED_ARCHIVE_SOURCE_SHA256,
    LEGACY_HELPER_SOURCE_SHA256,
    balance_residual,
    face_flux_from_centres,
    misner_sharp_face_mass,
    profile_peak_observables,
    rk4_proper_time_increment,
)

# 0*/GEAR two-mode Einstein-scalar engine
# Locked inputs from the archive:
#   g = 1, kappa = 0  ->  omega_S^2 = kappa+g = 1
#                         mu_D^2     = g-kappa = 1
# Stable/coupled mode S: V_S = +1/2 S^2
# Unstable/ejected mode D: V_D = -1/2 D^2
# No phenomenological matter-production law, damping coefficient,
# lifetime, or arrest threshold is inserted.

def run_engine(G=1.0, R=50.0, ngrid=300, A=0.01, eps=1e-10,
               width=7.0, T=25.0, cfl=0.12, sample_every=25):
    ms = md = 1.0
    h = R/ngrid
    r = (np.arange(ngrid)+0.5)*h
    rf = (np.arange(ngrid)+1.0)*h
    dt = cfl*h
    nsteps = int(T/dt)

    S = A*np.exp(-(r/width)**2)
    D = eps*np.exp(-(r/width)**2)
    PiS = np.zeros(ngrid)
    PiD = md*D                 # growing saddle branch, D_dot = +D

    def V(S,D):
        return 0.5*ms**2*S**2 - 0.5*md**2*D**2

    def dVS(S,D): return ms**2*S
    def dVD(S,D): return -md**2*D

    def metric(Y):
        S,D,PiS,PiD = Y
        PhiS = (np.r_[S[1:],-S[-1]] - np.r_[S[0],S[:-1]])/(2*h)
        PhiD = (np.r_[D[1:],-D[-1]] - np.r_[D[0],D[:-1]])/(2*h)

        a = np.ones(ngrid)
        for _ in range(80):
            rho = (PiS**2+PhiS**2+PiD**2+PhiD**2)/(2*a**2) + V(S,D)
            m = np.cumsum(4*np.pi*G*r*r*rho)*h - .5*h*4*np.pi*G*r*r*rho
            arg = 1-2*m/r
            if np.any(arg <= 1e-5): return None
            anew = 1/np.sqrt(arg)
            if np.max(np.abs(anew-a)) < 1e-9:
                a = anew
                break
            a = .125*anew + .875*a

        rho = (PiS**2+PhiS**2+PiD**2+PhiD**2)/(2*a**2) + V(S,D)
        mR = np.sum(4*np.pi*G*r*r*rho)*h
        if 1-2*mR/R <= 1e-5: return None
        aR = 1/np.sqrt(1-2*mR/R)

        Sg = ((a*a-1)/(2*r)
              + 2*np.pi*G*r*(PiS**2+PiD**2+PhiS**2+PhiD**2)
              - 4*np.pi*G*r*a*a*V(S,D))
        tail = np.cumsum((Sg*h)[::-1])[::-1] - .5*Sg*h
        lapse = np.exp(np.log(1/aR)-tail)
        return a,lapse,PhiS,PhiD,mR,aR

    def rhs(Y):
        M = metric(Y)
        if M is None: return None
        a,lapse,PhiS,PhiD,mR,aR = M
        last_center_lapse[0] = float(lapse[0])
        S,D,PiS,PiD = Y
        if balance_face < 0 or balance_face >= ngrid - 1:
            raise ValueError("fixed R=10 balance face must be an interior numerical face")
        kk = balance_face
        Rf_balance = rf[kk]
        mass_balance = misner_sharp_face_mass(Rf_balance, a[kk], a[kk+1])
        flux_balance = face_flux_from_centres(
            Rf_balance, h, lapse[kk], lapse[kk+1], a[kk], a[kk+1],
            S[kk], S[kk+1], PiS[kk], PiS[kk+1],
            D[kk], D[kk+1], PiD[kk], PiD[kk+1], G=G,
        )
        last_balance[0] = (float(mass_balance), float(flux_balance))
        St = lapse*PiS/a
        Dt = lapse*PiD/a

        dSf = np.r_[np.diff(S)/h, -S[-1]/(.5*h)]
        dDf = np.r_[np.diff(D)/h, -D[-1]/(.5*h)]
        af = np.r_[.5*(a[1:]+a[:-1]), aR]
        lf = np.r_[.5*(lapse[1:]+lapse[:-1]), 1/aR]

        fS = rf**2*lf*dSf/af
        fD = rf**2*lf*dDf/af
        return (St, Dt,
                np.r_[fS[0],np.diff(fS)]/(h*r**2) - lapse*a*dVS(S,D),
                np.r_[fD[0],np.diff(fD)]/(h*r**2) - lapse*a*dVD(S,D))

    Y = (S,D,PiS,PiD)
    history = []
    proper_profiles = []
    t = 0.0
    tau = 0.0
    last_center_lapse = [None]
    last_balance = [None]
    balance_history = []
    peak_history = []
    balance_face = int(round(10.0/h)) - 1

    for step in range(nsteps+1):
        if step % sample_every == 0:
            M = metric(Y)
            if M is None:
                return {"failed": True, "t": t, "history": np.asarray(history), "state": Y, "grid": (r,rf)}

            a,lapse,PhiS,PhiD,mR,aR = M
            S,D,PiS,PiD = Y

            rhoS = (PiS**2+PhiS**2)/(2*a*a) + .5*ms**2*S**2
            rhoD = (PiD**2+PhiD**2)/(2*a*a) - .5*md**2*D**2
            ES = 4*np.pi*np.sum(r*r*rhoS)*h
            ED = 4*np.pi*np.sum(r*r*rhoD)*h

            # The pre-existing scalar diagnostic samples are preserved.
            # These flux entries remain as saved-history context only.
            flux = []
            for rr in (10,20,30):
                if rr >= R:
                    flux.append(np.nan)
                    continue
                jf = int(round(rr/h)) - 1
                jf = max(0, min(jf, ngrid-2))
                Nf = 0.5*(lapse[jf] + lapse[jf+1])
                af = 0.5*(a[jf] + a[jf+1])
                PfS = 0.5*(PhiS[jf] + PhiS[jf+1])
                PfD = 0.5*(PhiD[jf] + PhiD[jf+1])
                PnS = 0.5*(PiS[jf] + PiS[jf+1])
                PnD = 0.5*(PiD[jf] + PiD[jf+1])
                Ftot = -4*np.pi*G*rr**2*Nf/af**3*(PfS*PnS + PfD*PnD)
                flux.append(Ftot)

            history.append((t,S[0],D[0],lapse[0]*aR,mR,
                            ES,ED,*flux))
            ell=(np.cumsum(a)-0.5*a)*h
            proper_profiles.append((t, tau, ell.copy(), a.copy(), lapse.copy(), S.copy(), D.copy(), PiS.copy(), PiD.copy(), PhiS.copy(), PhiD.copy()))
            peak = profile_peak_observables(r, a, lapse, PhiS, PiS, PhiD, PiD, G=G)
            peak_history.append((float(t), float(tau), peak["d_peak_areal_radius"],
                                 peak["d_peak_proper_radius"], peak["d_peak_signed_flux"]))

        # Always include the terminal state in the balance ledger, even when
        # sample_every does not divide the number of accepted steps.
        if step == nsteps and (not balance_history or balance_history[-1][0] < t):
            M_end = metric(Y)
            if M_end is None:
                return {"failed": True, "t": t, "history": np.asarray(history), "state": Y, "grid": (r,rf)}
            a_end, lapse_end, phi_s_end, phi_d_end, _, _ = M_end
            s_end, d_end, pi_s_end, pi_d_end = Y
            kk = balance_face
            Rf_balance = rf[kk]
            mass_balance = misner_sharp_face_mass(Rf_balance, a_end[kk], a_end[kk+1])
            flux_balance = face_flux_from_centres(
                Rf_balance, h, lapse_end[kk], lapse_end[kk+1], a_end[kk], a_end[kk+1],
                s_end[kk], s_end[kk+1], pi_s_end[kk], pi_s_end[kk+1],
                d_end[kk], d_end[kk+1], pi_d_end[kk], pi_d_end[kk+1], G=G,
            )
            balance_history.append((float(t), float(mass_balance), float(flux_balance)))

        if step == nsteps: break

        k1 = rhs(Y)
        if k1 is None: return {"failed": True, "t": t, "history": np.asarray(history), "state": Y, "grid": (r,rf)}
        l1 = last_center_lapse[0]
        balance_history.append((float(t), *last_balance[0]))
        Y2 = tuple(y+.5*dt*k for y,k in zip(Y,k1))
        k2 = rhs(Y2)
        l2 = last_center_lapse[0]
        if k2 is None: return {"failed": True, "t": t+.5*dt, "history": np.asarray(history)}
        Y3 = tuple(y+.5*dt*k for y,k in zip(Y,k2))
        k3 = rhs(Y3)
        l3 = last_center_lapse[0]
        if k3 is None: return {"failed": True, "t": t+.5*dt, "history": np.asarray(history)}
        Y4 = tuple(y+dt*k for y,k in zip(Y,k3))
        k4 = rhs(Y4)
        l4 = last_center_lapse[0]
        if k4 is None: return {"failed": True, "t": t+dt, "history": np.asarray(history)}

        Y = tuple(y + dt*(a1+2*a2+2*a3+a4)/6
                  for y,a1,a2,a3,a4 in zip(Y,k1,k2,k3,k4))
        t += dt
        tau += rk4_proper_time_increment(dt, (l1, l2, l3, l4))

    balance_array = np.asarray(balance_history, dtype=np.float64)
    if balance_array.shape[0] >= 2:
        balance_C = balance_residual(balance_array[:, 0], balance_array[:, 1], balance_array[:, 2])
    else:
        balance_C = float("nan")
    return {
        "failed": False, "t": t, "history": np.asarray(history),
        "proper_profiles": proper_profiles, "balance_history": balance_history,
        "balance_residual_R10": balance_C, "peak_history": np.asarray(peak_history),
        "provenance": {
            "corrected_archive_source_sha256": CORRECTED_ARCHIVE_SOURCE_SHA256,
            "legacy_helper_source_sha256": LEGACY_HELPER_SOURCE_SHA256,
            "diagnostic_mode": "observational_only_no_RHS_or_metric_equation_change",
        },
        "state": Y, "grid": (r,rf),
    }


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ngrid", type=int, default=300)
    parser.add_argument("--final-time", type=float, default=10.0)
    parser.add_argument("--sample-every", type=int, default=10)
    parser.add_argument("--eps", type=float, default=1e-10)
    args = parser.parse_args()
    result = run_engine(ngrid=args.ngrid, T=args.final_time,
                        sample_every=args.sample_every, eps=args.eps)
    final_peak = result.get("peak_history")
    summary = {
        "failed": result["failed"],
        "t": result["t"],
        "ngrid": args.ngrid,
        "final_time_requested": args.final_time,
        "balance_residual_R10": result.get("balance_residual_R10"),
        "final_peak": None if final_peak is None or len(final_peak) == 0 else {
            "t": float(final_peak[-1, 0]),
            "tau_center": float(final_peak[-1, 1]),
            "r_peak": float(final_peak[-1, 2]),
            "ell_peak": float(final_peak[-1, 3]),
            "signed_D_flux_peak": float(final_peak[-1, 4]),
        },
        "provenance": result.get("provenance"),
        "scientific_admission": "research_only_not_claimed",
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
