"""RESEARCH ONLY. The candidate local-to-cosmos channel of HANDOFF Nov 14 2026 (5), Task B3, evaluated on stored
frozen states of the lattice-cell runs (no new runs). Derivation: results/kappa-channel-2026-11-14/CHANNEL_NOTE.md.

On a CMC slice of the cell (zero lapse gradient at the wall):  dK/dt = <alpha Q>,  < > = proper-volume average,
    Q = (3/2) Aa^2 + K^2/3 + 4 pi (rho + S),   4 pi (rho + S) = [Pi^2 - V_c(phi)]            (phi, Einstein-normalised)
                                                              + 8 pi [PS^2 + PD^2 - V_SD(S, D)] (S and D)
                                                              + fluids (dust rho_E W^2 (1 + v^2)/2, radiation)
The cosmos lane evolves averages f_bar and so uses Q(f_bar) = K^2/3 + Pi_bar^2 - V_c(phi_bar) + fluids(f_bar)
(+ the constant of D at rest in its vacuum, one-budget runs). The channel is C = <alpha Q> - Q(f_bar):
    C = shear + S sector + D sector (above its vacuum) + phi variance (Var Pi - Jensen gap of V_c) + fluid part
        + lapse covariance <(alpha - 1) Q>.
Printed in units of Q(f_bar) - K^2/3 = (3/2)(rho_E + p_E)(f_bar) - ... i.e. as a fraction of the cosmos lane's own
source term. Second measure (B5): the same with S's action density as the weight in f_bar.
Usage: python research/cosmos_channel_v1.py RUN.npz [RUN.npz ...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

LAM = 5000.0
C2 = 304.92
V0_PROD = 8.242522415500654e-5
KAP = 8 * math.pi


def vsd(S, D):
    return 0.5 * S**2 - 0.5 * D**2 + LAM / 8 * (S**4 + 6 * S**2 * D**2 + D**4)


def vc(phi):
    return C2 * V0_PROD * np.exp(-phi)                       # one-budget rule: hand floor off (D_POT = 0)


def channel(d, i, weight=None):
    """Components of C at snapshot i. weight: None = proper volume; else an array used for the averaged fields."""
    dV = d["dV"][i]; al = d["alpha"][i]
    av = lambda f: float(np.sum(f * dV) / np.sum(dV))                                   # noqa: E731
    w = dV if weight is None else weight * dV
    bar = lambda f: float(np.sum(f * w) / np.sum(w))                                    # noqa: E731
    K = float(d["K"][i].mean())
    S, PS, D, PD, phi, Pi = (d[k][i] for k in ("S", "PS", "D", "PD", "phi", "Pi"))
    tot = 4 * math.pi * (d["rho"][i] + d["pr"][i] + 2 * d["pt"][i])
    floor = 0.0 if np.any(D != 0) else 0.10 * C2 * V0_PROD   # the no-D control restores the hand floor (amendment 2)
    q_phi = Pi**2 - (vc(phi) - floor)
    q_S = KAP * (PS**2 - (0.5 * S**2 + LAM / 8 * S**4))
    q_D = KAP * (PD**2 - (-0.5 * D**2 + LAM / 8 * (6 * S**2 * D**2 + D**4)))
    q_fl = tot - q_phi - q_S - q_D
    shear = 1.5 * d["Aa"][i]**2
    Q = shear + K * K / 3 + tot
    Dstar = math.sqrt(2 / LAM)
    q_vac = KAP * (-vsd(0.0, Dstar)) if np.any(D != 0) else 0.0                         # D at rest in its vacuum
    Qbar_src = bar(Pi)**2 - (vc(bar(phi)) - floor) + bar(q_fl) + q_vac                            # cosmos lane's source at f_bar
    parts = dict(
        shear=av(al * shear),
        S=av(al * q_S),
        D_above_vacuum=av(al * q_D) - q_vac * av(al),
        phi_var_Pi=av(Pi**2) - bar(Pi)**2,
        phi_Jensen_V=-(av(vc(phi)) - vc(bar(phi))),
        lapse_cov=av((al - 1) * (K * K / 3 + q_phi + q_fl)) + q_vac * (av(al) - 1),
    )
    total = av(al * Q) - (K * K / 3 + Qbar_src)
    parts["SD"] = parts["S"] + parts["D_above_vacuum"]
    parts["fluid_measure"] = av(q_fl) - bar(q_fl)                                      # zero for proper volume
    return dict(fl=av(q_fl), K=K, kdot=float(d["kdot"][i]), aQ=av(al * Q), meanQ=av(Q), src=Qbar_src, total=total, parts=parts,
                alpha_mean=av(al), alpha_wall=float(al[-1]))


if __name__ == "__main__":
    for path in sys.argv[1:]:
        d = np.load(path, allow_pickle=True); t = d["t"]; r = d["r"]
        print(f"\n=== {path.split('/')[-1]}: x_c = {float(d['x_c']):g}, {d['start']} start, r_max = {r[-1] + 0.5*(r[1]-r[0]):g}, closure {d['closure']}")
        print("  identity dK/dt = <alpha Q>: largest |kdot - <alpha Q>|/kdot over the run: "
              f"{max(abs(channel(d, i)['kdot'] - channel(d, i)['aQ']) / abs(channel(d, i)['kdot']) for i in range(0, len(t), 8)):.1e}")
        print("  channel C = <alpha Q> - Q(averaged fields), in units of the fluids' source term F = <(rho_dm + rho_b)/2 + rho_r> (always positive);")
        print("  each row is the mean over the snapshots in a window of 6 centred on t (S's and D's parts change sign within a cycle)")
        print("      t     alpha_wall  <alpha>   cosmos source/F   C/F     |  shear      S          D-vac      S+D        VarPi      JensenV    lapse     ||  C/F with S's action as the averaging weight")
        keys = ("shear", "S", "D_above_vacuum", "SD", "phi_var_Pi", "phi_Jensen_V", "lapse_cov")
        for tq in (3, 9, 15, 27, 39, 51, 63):
            sel = np.where((t >= tq - 3) & (t < tq + 3))[0]
            jw = np.sqrt((d["PS"][sel]**2).mean(0) * (d["S"][sel]**2).mean(0)) + 1e-300
            cs = [channel(d, i) for i in sel]; cj = [channel(d, i, weight=jw) for i in sel]
            m = lambda f: float(np.mean([f(c) for c in cs]))                            # noqa: E731
            print(f"    {tq:5.1f}   {m(lambda c: c['alpha_wall']):.4f}    {m(lambda c: c['alpha_mean']):.4f}   {m(lambda c: c['src']/c['fl']):+.4f}        {m(lambda c: c['total']/c['fl']):+.5f} | "
                  + " ".join(f"{m(lambda c, k=k: c['parts'][k]/c['fl']):+.3e}" for k in keys)
                  + f"  ||  {float(np.mean([c['total']/c['fl'] for c in cj])):+.4f}")
