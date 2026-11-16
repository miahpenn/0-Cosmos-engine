"""RESEARCH ONLY. Probe of the production kernel drift (results/floor-swap-2026-11-15/KERNEL_DRIFT.md): K-dot projected
with the stale lapse vs Q, at the floor-swap initial state. Usage: python research/kernel_drift_probe_v1.py"""
import sys, math, dataclasses
sys.path.insert(0, ".")
import numpy as np
import research.dm_one_budget_v1 as dm
from engine import cmc_gauge as cg
from engine.v55_matter import total_matter_projection
from engine.cosmos import HAND_FLOOR_REFERENCE_PARAMS, CosmosParams, construct_present_day, cosmos_rho_p
dm.LAM = 1.5246e6; dm.C2 = 1.0
k, st, r, H0 = dm.initial("noD", 0.0, 0.0, 20, 40.0)
g = st.geometry
w, W = cg._proper_volume_weights(st.grid, g)
tot = total_matter_projection(st.grid, g, st.scalars, st.matter)
Q = 1.5 * g.Aa**2 + g.K**2 / 3 + 4 * math.pi * (tot["rho"] + tot["pr"] + 2 * tot["pt"])
kd = cg._project_kdot_from_weights(st.grid, g, st.scalars, st.matter, w, W)
print("H0", H0, "K", g.K[0], "-> H", -g.K[0] / 3)
print("Q (uniform?)", Q[[0, 10, 19]], " kdot", kd, " kdot/Q", kd / Q[0])
# lane: K = -3H, dK/dt = -3 dH/dt = 1.5 (rho+p) with rho_E = 8 pi rho
p = HAND_FLOOR_REFERENCE_PARAMS; s0 = construct_present_day(p)  # p: hand-floor lane (labelled reference since step 2, Nov 16 2026)
y = (s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear)
rho, pr = cosmos_rho_p(y, p)
print("lane: H", s0.H, " dK/dt = 1.5(rho+p) =", 1.5 * (rho + pr), "  K^2/3 + (rho+3p)/2 =", (3*s0.H)**2/3 + (rho + 3*pr)/2)
print("cell rho_E = 8 pi rho:", 8*math.pi*tot["rho"][0], " lane rho:", rho, "; cell p_E:", 8*math.pi*(tot["pr"][0]), " lane p:", pr)
from engine import v55_pirk_adapter as adapter
import engine.scalar_system as ss
_, vac, _ = adapter.vendor_modules()
l2 = vac.primary_l2_rhs(st.grid, g)["K"]; l3 = adapter.primary_l3_with_matter(st.grid, g, st.scalars, st.matter)["K"]
print("L2 K-RHS", l2[[0,10,19]], " L3 K-RHS", l3[[0,10,19]], " sum", (l2+l3)[[0,10,19]])
print("K^2/3", g.K[0]**2/3, "  4pi(rho+pr+2pt)", 4*math.pi*(tot['rho'][0]+tot['pr'][0]+2*tot['pt'][0]), " 4pi(rho+S) alt  4pi(rho + 3pr)", 4*math.pi*(tot['rho'][0]+3*tot['pr'][0]))
print("alpha stored", g.alpha[[0,10,19]])
# which source? zero phi's potential and kinetic
for lab, mod in (("Pi=0", dict(Pi=np.zeros_like(st.scalars.Pi))),):
    sc = dataclasses.replace(st.scalars, **mod)
    l3b = adapter.primary_l3_with_matter(st.grid, g, sc, st.matter)["K"]; tb = total_matter_projection(st.grid, g, sc, st.matter)
    print(lab, "L3", l3b[0], " 4pi(rho+3p)", 4*math.pi*(tb['rho'][0]+tb['pr'][0]+2*tb['pt'][0]))
V0save, Dsave = ss.V0, ss.D_POT
ss.V0 = 0.0
l3c = adapter.primary_l3_with_matter(st.grid, g, st.scalars, st.matter)["K"]; tc = total_matter_projection(st.grid, g, st.scalars, st.matter)
print("V0=0:", "L3", l3c[0], " 4pi(rho+3p)", 4*math.pi*(tc['rho'][0]+tc['pr'][0]+2*tc['pt'][0]))
ss.V0 = V0save
