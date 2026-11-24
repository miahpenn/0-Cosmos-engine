"""W1 entry checks C1-C5, C7 (Lumen^2 design 16eY2t9Mdz4NxUpY_SHcucH7lCyQOxa9bMVtkOhKcars, section 6).

Gates were posted in the projects thread before this ran. Each check also runs a discrimination probe: a mutant
or perturbed input the check must FAIL, so a PASS means something.
Usage (repo root): python3 research/w1_supplied_slice_entry/checks.py BEFORE.npz MAIN.npz native.npz OUTDIR
  BEFORE = dump.py on 1a49290 (pre-#572), MAIN = dump.py on 3213444 (main), native.npz = Aster's W1 members.
"""
import json
import math
import os
import sys
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.getcwd())
import numpy as np  # noqa: E402

from engine import scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402
from engine.stress_energy import assemble_total_stress_energy  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import dump as D  # noqa: E402

PHI0, PI0 = 0.179055, 0.00341913
GRIDS = ((0.05, 1200), (0.025, 2400))
go, vac, _ = ad.vendor_modules()
out = {}


def fresh(mut=None):
    """Entry-branch dump in memory; mut patches the initializer for a discrimination probe."""
    data = {}
    if mut is None:
        D.main_into(data)
    else:
        with mut:
            D.main_into(data)
    return data


def rel(a, b):
    """max|a - b| / max|b| per array (several arrays cross zero, so a per-cell ratio is not defined there)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if np.array_equal(a, b):
        return 0.0
    return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), np.finfo(float).tiny))


def compare(new, old, prefix):
    keys = [k for k in old if k.split("/")[2] == prefix] if prefix else list(old)
    differ = sorted({k for k in keys if not np.array_equal(new[k], old[k])})
    worst = max((rel(new[k], old[k]) for k in keys), default=0.0)
    return differ, worst


before, mainref = np.load(sys.argv[1]), np.load(sys.argv[2])
z = np.load(sys.argv[3])
outdir = Path(sys.argv[4])
outdir.mkdir(parents=True, exist_ok=True)
entry = fresh()

# ---- C1 baseline identity
d_main, _ = compare(entry, mainref, None)
d_bef, w_bef = compare(entry, before, None)
bad_scal_matter = [k for k in d_bef if k.split("/")[2] == "scal"]   # matter follows from the geometry
c1 = dict(vs_1a49290_differ_keys=d_bef, vs_main_bitwise_differ=d_main, vs_1a49290_differ=len(d_bef), vs_1a49290_max_rel=w_bef,
          vs_1a49290_scalars_differ=bad_scal_matter,
          PASS=(w_bef <= 1e-13 and not bad_scal_matter))
# probe: the uniform-phi guard (5b) moves Kr at the binary64 floor; bitwise C1 vs main must catch it
real_proj = ss.scalar_momentum_projection
guard = mock.patch.object(vi._scalar_system, "scalar_momentum_projection",
                          lambda f, d: real_proj(f, (d[0], d[1], np.zeros_like(d[2]))))
d_probe, w_probe = compare(fresh(guard), mainref, None)
# probe for the 1e-13 branch: H0 * (1 + 1e-12) must exceed it
real_sqrt = math.sqrt
h0bump = mock.patch.object(vi.math, "sqrt", lambda x: real_sqrt(x) * (1 + 1e-12))
_, w_h0 = compare(fresh(h0bump), before, "geom")
c1["probe_uniform_phi_guard_flagged"] = bool(d_probe)
c1["probe_H0_1e-12_max_rel"] = w_h0
c1["probe_H0_flagged"] = w_h0 > 1e-13
out["C1"] = c1

# ---- C2 diagnostics identity (scalar_projection on the same slices)
d2m, _ = compare(entry, mainref, "proj")
d2b, w2b = compare(entry, before, "proj")
# probe: a one-ulp change in S in one cell must show
grid = go.SphericalCellGrid(1200, 60.0)
init = vi.build_initial_data(grid, vac.flat_state, initial_structure="rescaled_pocket", pocket_x=0.5)
p0 = ss.scalar_projection(grid, init.geometry, init.scalars)
S1 = init.scalars.S.copy()
S1[10] *= 1 + 1e-10
p1 = ss.scalar_projection(grid, init.geometry, ss.ScalarFields(S1, *[getattr(init.scalars, k)
                                                                     for k in ("PS", "D", "PD", "phi", "Pi")]))
out["C2"] = dict(vs_main_differ=d2m, vs_1a49290_differ=d2b, vs_1a49290_max_rel=w2b,
                 probe_cell_1e10_flagged=not all(np.array_equal(a, b) for a, b in zip(p0, p1)),
                 PASS=(not d2m and not d2b))


# ---- C3 round trip
def as_slice(sc, r, scale=1.0):
    return vi.SuppliedSlice(ss.ScalarFields(sc.S * scale, sc.PS, sc.D, sc.PD, sc.phi, sc.Pi), np.array(r), 0.0,
                            1e-3, "round trip of rescaled_pocket x=0.5")


c3 = {}
for dr, n in GRIDS:
    grid = go.SphericalCellGrid(n, 60.0)
    ref = vi.build_initial_data(grid, vac.flat_state, initial_structure="rescaled_pocket", pocket_x=0.5)
    got = vi.build_initial_data(grid, vac.flat_state, initial_structure="supplied_slice",
                                supplied_slice=as_slice(ref.scalars, grid.centers))
    mut = vi.build_initial_data(grid, vac.flat_state, initial_structure="supplied_slice",
                                supplied_slice=as_slice(ref.scalars, grid.centers, 1 + 1e-12))
    same = all(np.array_equal(getattr(got.geometry, k), getattr(ref.geometry, k)) for k in D.GEOM) and \
        all(np.array_equal(getattr(got.scalars, k), getattr(ref.scalars, k)) for k in vi.SUPPLIED_SLICE_FIELDS) and \
        got.H0 == ref.H0
    probe = all(np.array_equal(getattr(mut.geometry, k), getattr(ref.geometry, k)) for k in D.GEOM)
    c3[f"dr{dr}"] = dict(bit_identical=same, probe_rescale_1e12_flagged=not probe)
c3["PASS"] = all(v["bit_identical"] for v in c3.values())
out["C3"] = c3


# ---- W1 slices
def w1(dr, n):
    t = f"dr{dr}"
    grid = go.SphericalCellGrid(n, 60.0)
    f = ss.ScalarFields(z[t + "_S"], z[t + "_PS"], z[t + "_D"], z[t + "_PD"], np.full(n, PHI0), np.full(n, PI0))
    label = ("Aster B7-W1-grid-aligned-comparison.zip (Drive 1OUObtOwD7B51pDMiTacKJEjbavVBqv6D) native.npz sha256 "
             "16151c44...; grid-aligned harmonic-balance member, engine centres, r_max 60; chart-term flag not "
             "recorded in the archive (H0_chart supplied as 0.0)")
    return grid, vi.SuppliedSlice(f, np.array(grid.centers), 0.0, 1e-7, label)


kernel = V55ProductionKernel()
c4, c7 = {}, {}
for dr, n in GRIDS:
    grid, sl = w1(dr, n)
    st = kernel.initialize(resolution=n, r_max=60.0, initial_structure="supplied_slice", supplied_slice=sl)
    assert st.grid.dr == dr
    mom = vac.constraints(st.grid, st.geometry)["momentum"]
    total = assemble_total_stress_energy(st.grid, st.geometry, st.scalars, st.matter)
    res = float(np.max(np.abs(mom - 8 * np.pi * total.j)))
    Sp, Dp = vi._d1(st.grid, st.scalars.S, 1), vi._d1(st.grid, st.scalars.D, 1)
    jS, jD = -st.scalars.PS * Sp, -st.scalars.PD * Dp
    i = int(np.argmax(np.abs(jS)))
    res_noS = float(np.max(np.abs(mom - 8 * np.pi * (total.j - jS))))   # probe: drop the P_S S_r term
    c4[f"dr{dr}"] = dict(residual=res, PASS_residual=res <= 1e-12, cell=i, r=float(grid.centers[i]),
                         jS=float(jS[i]), jD=float(jD[i]), opposite_sign=bool(jS[i] * jD[i] < 0),
                         max_abs_jS=float(np.max(np.abs(jS))), probe_residual_without_jS=res_noS,
                         probe_flagged=res_noS > 1e-12)
    # C7: record every areal step
    steps, real = [], vi._cumulative_origin_integral

    def rec(g, r, h, f):
        F = real(g, r, h, f)
        steps.append(F / r)
        return F
    with mock.patch.object(vi, "_cumulative_origin_integral", rec):
        ini = vi.build_initial_data(grid, vac.flat_state, initial_structure="supplied_slice", supplied_slice=sl)
    seq = [float(np.max(np.abs(steps[0] - 1.0)))] + [float(np.max(np.abs(b - a))) for a, b in zip(steps, steps[1:])]
    c7[f"dr{dr}"] = dict(last=ini.areal_last_step, recorded_matches=ini.areal_last_step == seq[-1],
                         steps=seq, PASS=ini.areal_last_step <= 1e-13, probe_second_step_flagged=seq[1] > 1e-13)
c4["PASS"] = all(v["PASS_residual"] and v["opposite_sign"] for k, v in c4.items() if k != "PASS")
c7["PASS"] = all(v["PASS"] for k, v in c7.items() if k != "PASS")
out["C4"], out["C7"] = c4, c7

# ---- C5 refusals
grid, sl = w1(0.05, 1200)
f = sl.scalars
fields = {k: getattr(f, k) for k in vi.SUPPLIED_SLICE_FIELDS}


def variant(**kw):
    ff = dict(fields)
    ff.update({k: v for k, v in kw.items() if k in ff})
    return vi.SuppliedSlice(ss.ScalarFields(*[ff[k] for k in vi.SUPPLIED_SLICE_FIELDS]), kw.get("r", sl.r),
                            kw.get("H0_chart", 0.0), 1e-7, "C5 variant")


bump = fields["phi"].copy()
bump[3] += 1e-12
nan = fields["D"].copy()
nan[7] = np.nan
H0 = vi.build_initial_data(grid, vac.flat_state).H0
cases = {
    "misaligned_r": lambda: vi.build_initial_data(grid, vac.flat_state, initial_structure="supplied_slice",
                                                  supplied_slice=variant(r=sl.r + 1e-12)),
    "wrong_shape": lambda: vi.build_initial_data(grid, vac.flat_state, initial_structure="supplied_slice",
                                                 supplied_slice=variant(PS=fields["PS"][:-1])),
    "non_finite": lambda: vi.build_initial_data(grid, vac.flat_state, initial_structure="supplied_slice",
                                                supplied_slice=variant(D=nan)),
    "non_uniform_phi": lambda: vi.build_initial_data(grid, vac.flat_state, initial_structure="supplied_slice",
                                                     supplied_slice=variant(phi=bump)),
    "non_uniform_Pi": lambda: vi.build_initial_data(grid, vac.flat_state, initial_structure="supplied_slice",
                                                    supplied_slice=variant(Pi=fields["Pi"] * (1 + 1e-12))),
    "H0_chart_mismatch": lambda: vi.build_initial_data(grid, vac.flat_state, initial_structure="supplied_slice",
                                                       supplied_slice=variant(H0_chart=H0 * (1 + 1e-9))),
    "pocket_args": lambda: vi.build_initial_data(grid, vac.flat_state, initial_structure="supplied_slice",
                                                 supplied_slice=sl, amplitude=0.01),
    "lattice_wall": lambda: kernel.initialize(resolution=1200, r_max=60.0, initial_structure="supplied_slice",
                                              supplied_slice=sl, outer_edge="lattice_wall"),
}
c5 = {}
for name, fn in cases.items():
    try:
        fn()
        c5[name] = "no error"
    except ValueError as e:
        c5[name] = "ValueError: " + str(e)[:90]
    except Exception as e:  # noqa: BLE001
        c5[name] = f"{type(e).__name__}: {str(e)[:90]}"
# probe: the valid slice and an H0_chart equal to H0 both go through
vi.build_initial_data(grid, vac.flat_state, initial_structure="supplied_slice", supplied_slice=variant(H0_chart=H0))
c5["probe_valid_and_matching_H0_chart_accepted"] = True
c5["PASS"] = all(str(v).startswith("ValueError") for k, v in c5.items() if k in cases)
out["C5"] = c5

out["ALL_PASS"] = all(out[c]["PASS"] for c in ("C1", "C2", "C3", "C4", "C5", "C7"))
(outdir / "checks.json").write_text(json.dumps(out, indent=1))
print(json.dumps({c: out[c]["PASS"] for c in ("C1", "C2", "C3", "C4", "C5", "C7")}), "ALL", out["ALL_PASS"])
