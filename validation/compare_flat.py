"""
Build the flat with triplespec_reducer and compare it with IDL's
(validation/idl/ref_flat.pro, written to validation/idl_out/).
"""

import os
import sys
import time

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from triplespec_reducer import flat as F            # noqa: E402
from triplespec_reducer.frames import rotate        # noqa: E402

RAW = os.path.join(os.path.dirname(HERE), "data", "raw", "2025-11-02")
REF = os.path.join(HERE, "idl_out")
OUT = os.path.join(HERE, "py_out")
os.makedirs(OUT, exist_ok=True)

on = [os.path.join(RAW, f"SPEC_DFLAT_{i:04d}.fits") for i in range(2, 52)]
off = [os.path.join(RAW, f"SPEC_DFLAT_{i:04d}.fits") for i in range(52, 104)]

t0 = time.time()
res = F.make_flat(on, off, progress=lambda s: print("  ", s))
print(f"python flat: {time.time() - t0:.1f} s")
F.write_flat(os.path.join(OUT, "flat.fits"), res)


def stat(name, a, b, mask=None):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b) if mask is None else mask & np.isfinite(a) & np.isfinite(b)
    d = np.abs(a - b)[m]
    rel = d / np.maximum(np.abs(b[m]), 1e-30)
    print(f"  {name:22s} max|d| {d.max():.3e}  median|d| {np.median(d):.3e}  "
          f"max rel {rel.max():.3e}   n={m.sum()}")


# IDL WRITEFITS of an [x, y] array is already numpy's rotated [row, col].
ref_mean = fits.getdata(os.path.join(REF, "flat_mean.fits")).astype(float)
stat("scaled flat mean", res["combined"], ref_mean)

ref = F.read_flat(os.path.join(REF, "flat.fits"))
print("  edge coefficients, python - IDL, evaluated (pixels):")
x = np.arange(2048.)
for i, o in enumerate(ref["orders"]):
    for j, nm in ((0, "bottom"), (1, "top")):
        dp = F.I.poly(x, res["edgecoeffs"][i, j]) - F.I.poly(x, ref["edgecoeffs"][i, j])
        print(f"    order {o} {nm:6s} max |dy| {np.abs(dp).max():.2e} px")
print("  xranges python", res["xranges"].tolist(), " IDL", ref["xranges"].tolist())
inorder = ref["omask"] > 0
stat("normalized flat", res["flat"], ref["flat"], inorder)
stat("flat variance", res["var"], ref["var"], inorder)
print("  rms python", np.round(res["rms"], 7), "\n  rms IDL   ", np.round(ref["rms"], 7))
model_ref = fits.getdata(os.path.join(REF, "flat_model.fits")).astype(float)
stat("normalization model", res["model"], model_ref, inorder)
