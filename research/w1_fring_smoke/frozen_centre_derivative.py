"""Evaluate a frozen Aa derivative swap on archived W1 momentum traces."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess

import numpy as np

from engine.v55_pirk_adapter import vendor_modules
from research.w1_fring_smoke.origin_trace import VENDOR
from research.w1_fring_smoke.prepare_spatial import sha

CENTRE = Path("results/research/w1_fring_centre_projection")
TEAM = "b46c545a617dd52018c4ba61c8d5e43de11be604"
BAND = Path("results/research/w1_fring_spatial_rerun/runs")
WEIGHTS = np.array([4, -14/3, 2, 0])


def peak(values):
    return np.max(abs(values)).item()


def split(grid, vacuum, aa):
    """Split the change into the cell-0 departure and the remaining profile."""
    r = grid.centers
    smooth = aa.copy()
    smooth[0] = r[0]**2 * vacuum._quadratic_extrapolate_first(grid, aa/r**2)
    delta = aa[0] - smooth[0]
    difference = vacuum.regularized_aa_prime(grid, aa)-vacuum._d(grid, aa, 1)
    regular = vacuum.regularized_aa_prime(grid, smooth)-vacuum._d(grid, smooth, 1)
    spike = delta * WEIGHTS / grid.dr
    return difference, regular, spike, delta


def check_split(difference, regular, spike, atol):
    np.testing.assert_allclose(difference[:4], regular[:4] + spike, atol=atol, rtol=0)


def centre_traces(go, vacuum):
    source = json.loads((CENTRE/"summary.json").read_text())
    result = {}
    for case, group in source.items():
        for arm, expected in group.items():
            path = CENTRE/f"{case}-{arm}.npz"
            if sha(path) != expected["sha256"]:
                raise ValueError(f"archived centre trace checksum mismatch: {path}")
            with np.load(path, allow_pickle=False) as data:
                dt = json.loads(data["metadata"].item())
                dr = dt["dr"]
                grid = go.SphericalCellGrid(8,8*dr)
                r = grid.centers
                np.testing.assert_array_equal(r[:6], data["r"])
                # The first four stencil rows read only the six available native cells.
                aa = np.pad(data["Aa"], ((0,0),(0,2)), constant_values=np.nan)
                m = data["M"][:,:4]
                changes, regular, spikes, departures = [], [], [], []
                for i, profile in enumerate(aa):
                    old = vacuum._d(grid, profile, 1) + 3*profile/r
                    np.testing.assert_array_equal(old[:4], data["M_Aa"][i,:4])
                    dm, reg, spike, delta = split(grid,vacuum,profile)
                    # Bound evaluation error from the scale of the original derivative stencils.
                    scale = peak(profile[:6])/dr
                    check_split(dm, reg, spike, 32*np.finfo(float).eps*scale)
                    changes.append(dm[:4]); regular.append(reg[:4]); spikes.append(spike)
                    departures.append(delta)
                changes = np.stack(changes)
                adjusted = m + changes
                result[f"{case}-{arm}"] = {
                    "dr":dr, "dt":dt["dt"], "sha256":expected["sha256"], "t":data["t"].tolist(),
                    "M_current":m.tolist(), "M_regularized":adjusted.tolist(),
                    "delta_derivative":changes.tolist(), "delta_regular":np.stack(regular).tolist(),
                    "delta_spike":np.stack(spikes).tolist(), "delta_Aa0":departures,
                }
    return result


def shell_rms(r, values):
    # Use the coordinate-radius nested shells proposed in the band amendment.
    coarse = 0.025
    edges = np.arange(2, 11)*coarse
    counts = []
    average = []
    for a,b in zip(edges[:-1],edges[1:]):
        use = (r >= a) & (r < b)
        weight = r[use]**2
        counts.append(int(np.count_nonzero(use)))
        average.append(np.average(values[use],weights=weight))
    if counts != [round(coarse/(r[2]-r[1]))]*8:
        raise ValueError("unexpected nested shell occupancy")
    volumes = np.diff(edges**3)/3
    radius = (edges[:-1]+edges[1:])/2
    average = np.asarray(average)
    design = np.column_stack((np.ones(8), radius, radius**2))
    root = np.sqrt(volumes)
    smooth = design@np.linalg.lstsq(root[:,None]*design,average*root,rcond=None)[0]
    def rms(v):
        return float(np.sqrt(np.average(v**2,weights=volumes)))
    return {"rms":rms(average),"smooth_rms":rms(smooth),"rough_rms":rms(average-smooth),
            "shell_values":average.tolist()}


def team_band(go, vacuum):
    result = {}
    for mode in ("A","B"):
        for dr in (.025,.0125):
            for dt in (.000375,.0001875):
                name = f"{mode}-dr{dr:g}-dt{dt:g}"
                p = BAND/f"{name}-split.npz"
                content = subprocess.check_output(["git","show",f"{TEAM}:{p}"])
                with np.load(io.BytesIO(content),allow_pickle=False) as data:
                    r = data["r"]
                    grid = go.SphericalCellGrid(len(r),60.)
                    np.testing.assert_array_equal(grid.centers,r)
                    aa = r*data["three_Aa_over_r"]/3
                    dm, reg, spike, _ = split(grid,vacuum,aa)
                    scale = peak(aa)/dr
                    atol = 48*np.finfo(float).eps*scale
                    check_split(dm, reg, spike, atol)
                    np.testing.assert_allclose(dm[3:],reg[3:],atol=atol,rtol=0)
                    stored_derivative = (vacuum._d(grid,aa,1) if mode == "A"
                                         else vacuum.regularized_aa_prime(grid,aa))
                    np.testing.assert_allclose(stored_derivative,data["Aa_prime"],atol=atol,rtol=0)
                    # B already used the regularized derivative in its evolution and readback.
                    current = data["M"] if mode == "A" else data["M"]-dm
                    regularized = data["M"] if mode == "B" else data["M"]+dm
                    band = (r >= .05) & (r <= .25)
                    result[name] = {
                        "npz_sha256":hashlib.sha256(content).hexdigest(),
                        "archived_operator":"current" if mode=="A" else "regularized",
                        "current_native_band_max":peak(current[band]),
                        "regularized_native_band_max":peak(regularized[band]),
                        "derivative_change_native_band_max":peak(dm[band]),
                        "current":shell_rms(r,current), "regularized":shell_rms(r,regularized),
                        "derivative_change":shell_rms(r,dm),
                    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    vendor = subprocess.check_output(["git","-C","vendor/bb-palatini-unified-r0","rev-parse","HEAD"],text=True).strip()
    if vendor != VENDOR:
        raise ValueError("incorrect vendor pin")
    go,vacuum,_ = vendor_modules()
    report = {"vendor":vendor,"team_revision":TEAM,"script_sha256":sha(Path(__file__)),
              "centre_traces_revision":"fed7459769cada709466fc0338af60cfae6ac3c7",
              "centre_traces":centre_traces(go,vacuum),"team_band":team_band(go,vacuum)}
    args.output.write_text(json.dumps(report,allow_nan=False,separators=(",",":"))+"\n")


if __name__ == "__main__":
    main()
