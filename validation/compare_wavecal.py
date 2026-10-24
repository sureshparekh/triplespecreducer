"""
Wavelength-calibrate NGC134 with triplespec_reducer and compare every
intermediate with IDL's (validation/idl/ref_wavecal.pro). Uses IDL's flat,
so only the wavecal step is being tested.
"""

import os
import sys
import time

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from triplespec_reducer import flat as F, wavecal as W   # noqa: E402

RAW = os.path.join(os.path.dirname(HERE), "data", "raw", "2025-11-02")
REF = os.path.join(HERE, "idl_out")
OUT = os.path.join(HERE, "py_out")
os.makedirs(OUT, exist_ok=True)

files = [os.path.join(RAW, f"SPEC_SO2025D-006-{i:04d}.fits") for i in range(120, 132)]
flat = F.read_flat(os.path.join(REF, "flat.fits"))
t0 = time.time()
r = W.make_wavecal(files, flat, progress=lambda s: print("  ", s))
print(f"python wavecal: {time.time() - t0:.1f} s")
W.write_wavecal(os.path.join(OUT, "wavecal_NGC134.fits"), r)


def stat(name, a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    d = np.abs(a - b)[m]
    nan_same = np.array_equal(np.isfinite(a), np.isfinite(b))
    print(f"  {name:26s} max|d| {d.max() if d.size else 0:.3e}  median {np.median(d) if d.size else 0:.3e}"
          f"  n={m.sum()}  finite-mask {'same' if nan_same else 'DIFFERS'}")


for i, o in enumerate(r["orders"]):
    ref = fits.getdata(os.path.join(REF, f"wc_spec1d_{o}.fits")).astype(float)  # [3, nx]
    stat(f"1D spectrum order {o}", r["spec1d"][i][:2], ref[:2])
print(f"  offset python {r['offset']:.6f}   IDL -0.656294")
print("  nlines python", r["sol1d"]["nlines"], "  IDL [31, 64, 63, 76, 32]")
print(f"  rms1d python {r['sol1d']['rms']:.6f}   IDL 0.34187767")
idl_c = [1.3863339, 0.00023600114, -1.9454762e-08, 6.3993170e-12, 0.013314428,
         -1.3572513e-05, 1.2849736e-08, -4.2275137e-12, -0.0022233863, 3.0858513e-06,
         -2.7446726e-09, 8.5297285e-13, 0.00014450317, -2.1819946e-07, 1.9004792e-10,
         -5.6745503e-14]
print("  1D coeff max rel diff", np.max(np.abs((r["sol1d"]["coeffs"] - idl_c) / np.array(idl_c))))
for i, o in enumerate(r["orders"]):
    ref = fits.getdata(os.path.join(REF, f"wc_lines2d_{o}.fits")).astype(float)
    ours = r["lines2d"][i]
    n_ref, n_py = ref.shape[1], ours["sx"].size
    print(f"  order {o}: traced points python {n_py}  IDL {n_ref}", end="")
    if n_ref == n_py:
        d = [np.max(np.abs(ours[k] - ref[j])) for j, k in enumerate(("sx", "sy", "ss", "sl"))]
        print("   max |d| sx %.2e sy %.2e ss %.2e sl %.2e" % tuple(d))
    else:
        print()
ref = W.read_wavecal(os.path.join(REF, "wavecal_NGC134.fits"))
ours = W.read_wavecal(os.path.join(OUT, "wavecal_NGC134.fits"))
stat("wavelength image (um)", ours["wavecal"], ref["wavecal"])
stat("spatial image (arcsec)", ours["spatcal"], ref["spatcal"])
print("  dispersion python", np.round(ours["disp"], 9), "\n  dispersion IDL   ", np.round(ref["disp"], 9))
