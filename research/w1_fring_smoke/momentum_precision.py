"""Zero-step manufactured and saved-state momentum precision readback."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import io
import json
from pathlib import Path
import subprocess

import numpy as np

from engine.scalar_system import ScalarFields
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_pirk_adapter import vendor_modules
from research.w1_fring_smoke.frozen_centre_derivative import BAND, TEAM, peak, shell_rms
from research.w1_fring_smoke.origin_trace import VENDOR

ARMS = Path("results/research/w1_fring_a6/arms")
TERMS = ("Aa_prime", "three_Aa_over_r", "Aa_b_prime_over_b",
         "six_Aa_chi_prime", "K_part", "j_part")
COEFFICIENTS = tuple(Decimal(v) for v in ("2e-7", "-3e-7", "4e-7"))


def wide(value):
    return Decimal.from_float(float(value))


def derivative(values, i, h):
    # Interior rows and the same even ghost mapping used by the vendor.
    def at(j):
        return values[j if j >= 0 else -j-1]
    return (at(i-2)-8*at(i-1)+8*at(i+1)-at(i+2))/(12*h)


def aa_derivative(aa, r, quotient, i, h, regularized):
    if regularized:
        return 2*aa[i]/r[i] + r[i]**2*derivative(quotient,i,h)
    return derivative(aa,i,h)


def stats(r, values):
    band = (r >= .05) & (r < .25)
    return {"cell0_signed":float(values[0]),"band_max_abs":peak(values[band])}


def manufactured(go, vacuum):
    A,B,C = COEFFICIENTS
    output = {}
    for dr in (.025,.0125):
        grid = go.SphericalCellGrid(round(1/dr),1.)
        r = grid.centers
        indices = range(np.searchsorted(r,.25))
        h = Decimal(str(dr))
        ideal_r = [(Decimal(i)+Decimal("0.5"))*h for i in range(grid.n)]
        ideal_aa = [A*x**2+B*x**4+C*x**6 for x in ideal_r]
        ideal_j = [5*A*x+7*B*x**3+9*C*x**5 for x in ideal_r]  # 8 pi j
        state = vacuum.flat_state(grid)
        state.K.fill(-.017)
        state.Aa[:] = np.array(ideal_aa,dtype=float)
        j = np.array(ideal_j,dtype=float)
        float_current = vacuum.momentum_constraint(grid,state)-j
        stored_r = list(map(wide,r))
        stored_aa = list(map(wide,state.Aa))
        stored_j = list(map(wide,j))
        arms = {}
        for regularized in (False,True):
            quotient = [a/x**2 for a,x in zip(ideal_aa,ideal_r)]
            stored_q = [a/x**2 for a,x in zip(stored_aa,stored_r)]
            binary = (float_current - vacuum._d(grid,state.Aa,1)
                      + vacuum.regularized_aa_prime(grid,state.Aa)) if regularized else float_current
            ideal, stored = [], []
            for i in indices:
                x = ideal_r[i]
                actual = aa_derivative(ideal_aa,ideal_r,quotient,i,h,regularized) + 3*ideal_aa[i]/x - ideal_j[i]
                expected = Decimal(0) if regularized else -24*C*h**4*x
                if abs(actual-expected) > Decimal("1e-58"):
                    raise AssertionError("manufactured derivative fails its analytic residual")
                ideal.append(actual)
                x = stored_r[i]
                stored.append(aa_derivative(stored_aa,stored_r,stored_q,i,wide(dr),regularized)
                              + 3*stored_aa[i]/x-stored_j[i])
            ideal_values = np.array(ideal,dtype=float)
            stored_values = np.array(stored,dtype=float)
            band_r = r[:len(indices)]
            arms["regularized" if regularized else "current"] = {
                "binary64":stats(band_r,binary[:len(indices)]),
                "wide_stored_binary64_inputs":stats(band_r,stored_values),
                "wide_ideal_inputs":stats(band_r,ideal_values),
                "max_abs_readback_rounding":peak(binary[:len(indices)]-stored_values),
                "max_abs_input_representation":peak(stored_values-ideal_values),
                "ideal_residual_over_r_at_cell0":float(ideal[0]/ideal_r[0]),
            }
        output[dr] = arms
    return output


def read_git(path):
    return subprocess.check_output(["git","show",f"{TEAM}:{path}"])


def team_fine(go, vacuum):
    result = {}
    for dr in (.025,.0125):
        for dt in (.000375,.0001875):
            name = f"A-dr{dr:g}-dt{dt:g}"
            content = read_git(BAND/f"{name}-split.npz")
            with np.load(io.BytesIO(content),allow_pickle=False) as data:
                r = data["r"]
                grid = go.SphericalCellGrid(len(r),60.)
                np.testing.assert_array_equal(r,grid.centers)
                # This archive saves 3 Aa/r and the already computed momentum terms.
                # Invert the stored radial term at wide precision for the derivative-only readback.
                radii = list(map(wide,r))
                aa = [wide(t)*x/3 for t,x in zip(data["three_Aa_over_r"],radii)]
                band = (r >= .05) & (r < .25)
                frozen_sum = np.zeros_like(r)
                frozen_stencil = np.zeros_like(r)
                archive_terms = [data[key] for key in TERMS]
                recovered_aa = np.array(aa,dtype=float)
                np.testing.assert_allclose(vacuum._d(grid,recovered_aa,1)[band],data["Aa_prime"][band],
                                           atol=48*np.finfo(float).eps*peak(recovered_aa)/dr,rtol=0)
                for i in np.flatnonzero(band):
                    frozen_sum[i] = float(sum(wide(term[i]) for term in archive_terms))
                    fresh = derivative(aa,i,wide(dr))
                    frozen_stencil[i] = float(fresh + sum(wide(term[i]) for term in archive_terms[1:]))
                baseline = data["M"]
                result[name] = {
                    "archive_sha256":hashlib.sha256(content).hexdigest(),
                    "archived":shell_rms(r,baseline),
                    "wide_sum":shell_rms(r,frozen_sum),
                    "wide_sum_and_Aa_stencil":shell_rms(r,frozen_stencil),
                    "max_abs_difference_wide_sum":peak((frozen_sum-baseline)[band]),
                    "max_abs_difference_wide_stencil":peak((frozen_stencil-baseline)[band]),
                }
    return result


def full_coarse(go, vacuum):
    grid = go.SphericalCellGrid(1200,60.)
    r = grid.centers
    radii = list(map(wide,r))
    h = wide(grid.dr)
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=1200,r_max=60.,outer_edge="open_edge",
                              eight_node_areal=True,initial_structure="none")
    result = {}
    for arm in ("A1","A2"):
        rec = json.loads(read_git(ARMS/f"{arm}.json"))
        content = read_git(ARMS/f"{arm}-t0.015.npz")
        with np.load(io.BytesIO(content),allow_pickle=False) as data:
            for key in vacuum.ALL_FIELDS:
                np.copyto(getattr(state.geometry,key),data[f"geometry_{key}"])
            state.scalars = ScalarFields(*(data[f"scalars_{key}"] for key in ("S","PS","D","PD","phi","Pi")))
            for species in ("dark_matter","baryons","radiation"):
                field = getattr(state.matter,species)
                for key in ("rest","energy_t","momentum_r"):
                    np.copyto(getattr(field,key),data[f"matter_{species}_{key}"])
        g = state.geometry
        j = 8*np.pi*assemble_total_stress_energy(grid,g,state.scalars,state.matter).j
        current = vacuum.momentum_constraint(grid,g)-j
        np.testing.assert_array_equal(current[:6],rec["saves"]["0.015"]["row"]["M"][:6])
        aa,K,b,X = (list(map(wide,getattr(g,key))) for key in ("Aa","K","b","X"))
        band = (r >= .05) & (r < .25)
        selected = np.flatnonzero(band | (np.arange(grid.n)==0))
        values = np.zeros_like(r)
        for i in selected:
            values[i] = float(
                derivative(aa,i,h) - 2*derivative(K,i,h)/3
                -3*aa[i]*derivative(X,i,h)/X[i]
                + 3*aa[i]/2*(2/radii[i]+derivative(b,i,h)/b[i]) - wide(j[i])
            )
        result[arm] = {
            "archive_sha256":hashlib.sha256(content).hexdigest(),
            "binary64":stats(r,current),
            "wide_geometry_same_binary64_j":stats(r,values),
            "max_abs_readback_difference":peak((values-current)[selected]),
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    root = "vendor/bb-palatini-unified-r0"
    pin = subprocess.check_output(["git","-C",root,"rev-parse","HEAD"],text=True).strip()
    if pin != VENDOR:
        raise ValueError("unexpected vendor revision")
    go,vacuum,_ = vendor_modules()
    with localcontext(prec=60):
        result = {"vendor":pin,"archive_revision":TEAM,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  "decimal_digits":60,"coefficients":[str(v) for v in COEFFICIENTS],
                  "manufactured":manufactured(go,vacuum), "fine_band":team_fine(go,vacuum),
                  "complete_coarse":full_coarse(go,vacuum)}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
