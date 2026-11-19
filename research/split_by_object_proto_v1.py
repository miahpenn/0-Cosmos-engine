"""RESEARCH ONLY. Stage 1 of the split-by-object design (coordination log, Nov 19 2026; claude/split-by-object-design.md).
Reads engine/ without changing it. Two objects, stones and free D_r, each with its own coupling and its own transport, on a FROZEN kernel state
(metric and phi frozen, exchange OFF, D_r pressure not fed to the metric: FORCED). The kernel's own dust routines (species_rhs, evolve_species) are used for every object.
Parts: A lane-level identities (r(phi), w_r); B frozen-state tests (i)-(v) on the kernel's state with an IMPOSED phi bump (FORCED, so that grad phi != 0)."""
import math, sys, dataclasses
sys.path.insert(0, ".")
import numpy as np
from engine import model_constants as mc, scalar_system as ss
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import dm_density, metric_slice_from_q
from engine.matter_system import geometry_metric_derivatives, ConservedSpecies, Species, evolve_species
from engine.matter_rhs import species_rhs

def f_st(phi): u = np.exp(-2 * mc.BETA_D * phi); return 3 * u / (1 + 3 * u)
def f_D(phi, r=mc.SD_DR_SHARE_TODAY, pt=mc.PHI_TODAY): return mc.sd_matter_coupled_fraction(np.asarray(phi, float), r_today=r, phi_today=pt)
def r_of(phi, r=mc.SD_DR_SHARE_TODAY, pt=mc.PHI_TODAY):
    fs = f_st(phi); return (f_D(phi, r, pt) - fs) / (1 - fs)
def scale(U, c): return ConservedSpecies(U.rest * c, U.energy_t * c, U.momentum_r * c)
def add(A, B): return ConservedSpecies(A.rest + B.rest, A.energy_t + B.energy_t, A.momentum_r + B.momentum_r)
def comps(U): return np.array([U.rest, U.energy_t, U.momentum_r])

def frozen_state(N, amp, width, r0=20.0):
    k = V55ProductionKernel(); st = k.initialize(resolution=N, r_max=40.0, outer_edge="open_edge")
    r = np.asarray(st.grid.centers); phi = np.asarray(st.scalars.phi) + amp * np.exp(-((r - r0) / width) ** 2)
    sc = dataclasses.replace(st.scalars, phi=phi)
    st = dataclasses.replace(st, scalars=sc)
    metric = metric_slice_from_q(st.grid, st.geometry)
    md = geometry_metric_derivatives(metric, lambda v, p: k._d1(st.grid, v, p))
    return k, st, metric, md

def objects(k, st, metric, md, r_today=mc.SD_DR_SHARE_TODAY, phi_today=mc.PHI_TODAY):
    phi = np.asarray(st.scalars.phi); rho = dm_density(metric, st.matter); U = st.matter.dark_matter
    rc = r_of(phi, r_today, phi_today); fs = f_st(phi); fd = f_D(phi, r_today, phi_today)
    b_eff, b_st, b_dr = mc.BETA_D * fd, mc.BETA_D * fs, mc.BETA_D * np.ones_like(phi)
    dphi_r = k._d1(st.grid, phi, 1)
    def srhs(beta, rho_i): return ss.scalar_rhs_arrays(st.grid, st.geometry, st.scalars, beta_dm=beta, rho_dm=rho_i)
    dphi_t = srhs(b_eff, rho).phi
    return dict(phi=phi, rho=rho, U=U, rc=rc, fs=fs, fd=fd, b_eff=b_eff, b_st=b_st, b_dr=b_dr, dphi_r=dphi_r, dphi_t=dphi_t, srhs=srhs)

def rhs(o, metric, md, U, beta):
    return species_rhs(metric, md, U, Species.DARK_MATTER, dphi_t=o["dphi_t"], dphi_r=o["dphi_r"], beta_dm=beta)

def part_B(N, amp, width, r_today=mc.SD_DR_SHARE_TODAY, phi_today=mc.PHI_TODAY, verbose=True, steps=(1, 5, 20)):
    k, st, metric, md = frozen_state(N, amp, width); o = objects(k, st, metric, md, r_today, phi_today)
    rc = o["rc"]; U = o["U"]
    # shares by cell: the same state scaled (same velocity); rest, energy and momentum are linear in density at fixed velocity
    Ust = ConservedSpecies(U.rest * (1 - rc), U.energy_t * (1 - rc), U.momentum_r * (1 - rc))
    Udr = ConservedSpecies(U.rest * rc, U.energy_t * rc, U.momentum_r * rc)
    R_L, R_S, R_D = rhs(o, metric, md, U, o["b_eff"]), rhs(o, metric, md, Ust, o["b_st"]), rhs(o, metric, md, Udr, o["b_dr"])
    R_sum = add(R_S, R_D); a, b = comps(R_L), comps(R_sum)
    i_max = float(np.max(np.abs(a - b)) / np.max(np.abs(a)))
    # (ii) pull per unit density: source part = RHS(beta) - RHS(0)
    z = np.zeros_like(o["b_eff"])
    P_L = comps(R_L) - comps(rhs(o, metric, md, U, z)); P_S = comps(R_S) - comps(rhs(o, metric, md, Ust, z)); P_D = comps(R_D) - comps(rhs(o, metric, md, Udr, z))
    mom = 2; m = np.abs(P_L[mom]) > 1e-3 * np.max(np.abs(P_L[mom]))
    lump_on_stone = (P_L[mom][m] / (1 - rc[m])) / (P_S[mom][m] / (1 - rc[m])) * (1 - rc[m]) - 1   # lumped pull per unit density / stones' own pull per unit density, -1
    # per unit density: lumped (beta_eff), stones (beta_st), D_r (beta_dr); same velocity so pull/rho ~ beta
    ps = P_S[mom][m] / (1 - rc[m]); pd = P_D[mom][m] / rc[m]; pl = P_L[mom][m]
    err_st = float(np.max(np.abs(pl / ps - 1))); err_dr = float(np.max(np.abs(pl / pd - 1)))
    split_check_st = float(np.max(np.abs(ps / (o["b_st"][m] * (pl / o["b_eff"][m])) - 1))); split_check_dr = float(np.max(np.abs(pd / (o["b_dr"][m] * (pl / o["b_eff"][m])) - 1)))
    # (iii) phi's source
    base = o["srhs"](z, z * 0 + 1e-300); L = o["srhs"](o["b_eff"], o["rho"]); S1 = o["srhs"](o["b_st"], (1 - rc) * o["rho"]); S2 = o["srhs"](o["b_dr"], rc * o["rho"])
    dL = L.Pi - base.Pi; dS = (S1.Pi - base.Pi) + (S2.Pi - base.Pi); iii = float(np.max(np.abs(dL - dS)) / max(np.max(np.abs(dL)), 1e-300))
    # (iv) negative control: beta_dr := beta_st
    R_Dn = rhs(o, metric, md, Udr, o["b_st"]); Pn = comps(R_Dn) - comps(rhs(o, metric, md, Udr, z)); pdn = Pn[mom][m] / rc[m]
    neg_dr_def = float(np.median(1 - pdn / pd)); lump_dr_def = float(np.median(1 - pl / pd)); lump_st_err = float(np.median(pl / ps - 1))
    # few-step frozen evolution of the three: velocity of each object at the cell of strongest pull
    dt = 0.06 * float(st.grid.dr); evo = {}
    cL, cS, cD = U, Ust, Udr
    j = int(np.argmax(np.abs(pl))); out = []
    for n in range(1, max(steps) + 1):
        cL = evolve_species(metric, md, cL, Species.DARK_MATTER, dt, dphi_t=o["dphi_t"], dphi_r=o["dphi_r"], beta_dm=o["b_eff"], validate_admissibility=False)
        cS = evolve_species(metric, md, cS, Species.DARK_MATTER, dt, dphi_t=o["dphi_t"], dphi_r=o["dphi_r"], beta_dm=o["b_st"], validate_admissibility=False)
        cD = evolve_species(metric, md, cD, Species.DARK_MATTER, dt, dphi_t=o["dphi_t"], dphi_r=o["dphi_r"], beta_dm=o["b_dr"], validate_admissibility=False)
        if n in steps:
            vL = cL.momentum_r[j] / cL.energy_t[j]; vS = cS.momentum_r[j] / cS.energy_t[j]; vD = cD.momentum_r[j] / cD.energy_t[j]
            tot = (cS.momentum_r + cD.momentum_r); out.append((n, float(vL / vS - 1), float(vL / vD - 1), float(vD / vS - 1), float(np.max(np.abs(tot - cL.momentum_r)) / max(np.max(np.abs(cL.momentum_r)), 1e-300))))
    res = dict(N=N, amp=amp, width=width, r_today=r_today, phi_today=phi_today, i_rhs_sum=i_max, err_st=err_st, err_dr=err_dr, split_check_st=split_check_st, split_check_dr=split_check_dr,
               iii=iii, neg_dr_def=neg_dr_def, lump_dr_def=lump_dr_def, lump_st_err=lump_st_err, evo=out, phi_range=(float(o["phi"].min()), float(o["phi"].max())))
    return res

def part_A():
    print("A. LANE-LEVEL IDENTITIES")
    A_ = (1 - mc.SD_DR_SHARE_TODAY) / math.sqrt(1 + 3 * math.exp(-2 * mc.BETA_D * mc.PHI_TODAY)); B_ = mc.SD_DR_SHARE_TODAY / math.sqrt(math.exp(-2 * mc.BETA_D * mc.PHI_TODAY))
    worst = 0.0
    for phi in np.linspace(0, 3, 13):
        u = math.exp(-2 * mc.BETA_D * phi); ms, md_ = A_ * math.sqrt(1 + 3 * u), B_ * math.sqrt(u); rs = md_ / (ms + md_)
        worst = max(worst, abs(rs - float(r_of(phi))))
    print(f"(vi) r(phi) from two number-conserving objects (rest energy A sqrt(1+3u), B sqrt(u)) vs the lump's r(phi): max |difference| {worst:.2e}")
    print("     r(phi) at phi = 0, today, 2, 3:", [round(float(r_of(p)), 5) for p in (0, mc.PHI_TODAY, 2, 3)], "(READ today 0.07, bounds 0.06-0.08)")
    print("(vii) D_r pressure: w = k^2/(3 omega^2), m = sqrt2, k ~ 1/a: energy per quantum falls from omega0 to m, share r falls by q = m/omega0 then holds")
    for w0 in (2.0, 4.0):
        q = math.sqrt(2) / w0; w_birth = (1 - q * q) / 3
        row = [f"omega0 {w0:g}: q {q:.3f}, w at birth {w_birth:.3f}"]
        for ra in (0.06, 0.07, 0.08): row.append(f"r_after {ra} <- r_birth {ra / (q + ra * (1 - q)):.3f}")
        print("     " + "; ".join(row))
    print("     AMEND-3-GEAR-100 travelling D w 0.13-0.27: omega0 2 gives w = (1-q^2)/3 = %.3f (inside), omega0 4 gives %.3f (outside)" % ((1-0.5)/3, (1-0.125)/3))

if __name__ == "__main__":
    part_A()
    print("B. FROZEN KERNEL STATE (phi bump FORCED; exchange OFF; D_r pressure not fed to metric)")
    rows = []
    for N in (20, 40):
        for amp, width in ((0.05, 5.0), (0.05, 2.5), (0.15, 5.0)):
            r = part_B(N, amp, width); rows.append(r)
            print(f"N {N} bump {amp}/{width}: phi {r['phi_range'][0]:.3f}-{r['phi_range'][1]:.3f} | (i) RHS sum vs lump {r['i_rhs_sum']:.2e} | (ii) lumped-vs-true pull: stones {r['err_st']:.3%} D_r {r['err_dr']:.3%} (max over cells); split-vs-own-beta {r['split_check_st']:.1e}/{r['split_check_dr']:.1e} | (iii) phi source sum vs lump {r['iii']:.2e} | (iv) lump stone err {r['lump_st_err']:.3%}, lump D_r deficit {r['lump_dr_def']:.3%}, negative-control D_r deficit {r['neg_dr_def']:.3%}", flush=True)
            for n, a, b, c, d in r["evo"]: print(f"     steps {n}: v_lump/v_stones-1 {a:+.3%}, v_lump/v_D_r-1 {b:+.3%}, v_D_r/v_stones-1 {c:+.3%}, |sum of objects - lump| momentum {d:.2e}")
    print("(v) RULERS: r_today and phi_today (N 40, bump 0.05/5):")
    for rt in (0.06, 0.07, 0.08):
        r = part_B(40, 0.05, 5.0, r_today=rt); print(f"   r_today {rt}: lump stone err {r['lump_st_err']:.3%}, lump D_r deficit {r['lump_dr_def']:.3%}, negative control D_r deficit {r['neg_dr_def']:.3%}")
    for pt in (mc.PHI_TODAY - 0.1, mc.PHI_TODAY, mc.PHI_TODAY + 0.1):
        r = part_B(40, 0.05, 5.0, phi_today=pt); print(f"   phi_today {pt:.3f}: lump stone err {r['lump_st_err']:.3%}, lump D_r deficit {r['lump_dr_def']:.3%}")
