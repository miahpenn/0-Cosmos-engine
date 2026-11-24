"""Dump the four archived initial-data structures (open edge, dr 0.05 and 0.025, r_max 60) for checks C1/C2.

Run from a checkout's root: python3 research/w1_supplied_slice_entry/dump.py <out.npz>. Works on any commit
since #503 (it uses only build_initial_data's existing arguments), so it is run on 1a49290, 3213444 and the
entry branch. Writes every state array, the scalars, the matter state and scalar_projection's outputs.
"""
import contextlib
import os
import sys

sys.path.insert(0, os.getcwd())
import numpy as np  # noqa: E402

from engine import scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad  # noqa: E402

CASES = (("none", {}), ("rescaled_pocket", {"pocket_x": 0.5}), ("inherited_pocket", {}),
         ("hand_floor_reference", {}))
GEOM = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")


def dump(init, grid, tag, out):
    for k in GEOM:
        out[f"{tag}/geom/{k}"] = np.array(getattr(init.geometry, k))
    for k in ("S", "PS", "D", "PD", "phi", "Pi"):
        out[f"{tag}/scal/{k}"] = np.array(getattr(init.scalars, k))
    for sp in ("dark_matter", "baryons", "radiation"):
        for f in ("rest", "energy_t", "momentum_r"):
            out[f"{tag}/matter/{sp}.{f}"] = np.array(getattr(getattr(init.matter, sp), f))
    for k, v in zip(("e", "pr", "pt", "j"), ss.scalar_projection(grid, init.geometry, init.scalars)):
        out[f"{tag}/proj/{k}"] = np.array(v)
    out[f"{tag}/H0"] = np.array(init.H0)


def main_into(out):
    go, vac, _ = ad.vendor_modules()
    for dr in (0.05, 0.025):
        grid = go.SphericalCellGrid(round(60 / dr), 60.0)
        for name, kw in CASES:
            ctx = ss.hand_floor_reference() if name == "hand_floor_reference" else contextlib.nullcontext()
            with ctx:
                init = vi.build_initial_data(grid, vac.flat_state, initial_structure=name, **kw)
                dump(init, grid, f"dr{dr}/{name}", out)
    return out


def main(path):
    np.savez(path, **main_into({}))


if __name__ == "__main__":
    main(sys.argv[1])
