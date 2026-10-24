"""
Point-source (HD225200 116-117) and extended-source (NGC134 120-121)
extraction, compared with IDL (validation/idl/ref_extract.pro). Uses IDL's
flat and wavecal, so only extraction is being tested.
"""

import os
import sys
import time

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from triplespec_reducer import flat as F, wavecal as W, extract as E   # noqa: E402
from triplespec_reducer.instrument import read_atran                    # noqa: E402

RAW = os.path.join(os.path.dirname(HERE), "data", "raw", "2025-11-02")
REF = os.path.join(HERE, "idl_out")
fl = F.read_flat(os.path.join(REF, "flat.fits"))
wc = W.read_wavecal(os.path.join(REF, "wavecal_NGC134.fits"), rotation=fl["rotation"])
aw, at = read_atran(fl["rp"])
orders = fl["orders"]


def stat(name, a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.shape != b.shape:
        print(f"  {name:28s} SHAPE {a.shape} vs IDL {b.shape}")
        return
    m = np.isfinite(a) & np.isfinite(b)
    d = np.abs(a - b)[m]
    rel = d / np.maximum(np.abs(b[m]), 1e-30)
    print(f"  {name:28s} max|d| {d.max():.3e}  max rel {np.median(rel):.1e} (median) "
          f"{np.percentile(rel, 99):.1e} (99%)  finite {'same' if np.array_equal(np.isfinite(a), np.isfinite(b)) else 'DIFFERS'}")


# ---------------------------------------------------------------- point source
f = [os.path.join(RAW, "SPEC_SO2025D-006-0116.fits"), os.path.join(RAW, "SPEC_SO2025D-006-0117.fits")]
work, var, hdrs, itot = E.load_image(f, fl)
stat("work image", work, fits.getdata(os.path.join(REF, "ps_work.fits")))
t0 = time.time()
prof = E.mkspatprof2d(work, fl["omask"], wc["wavecal"], wc["spatcal"],
                      fl["slith_arc"] / fl["slith_pix"], 1, aw, at, 0.4)
print(f"  profiles {time.time() - t0:.1f} s")
for o in orders:
    ref = fits.getdata(os.path.join(REF, f"ps_prof_{o}.fits"))
    stat(f"profile order {o}", np.array(prof[o]), ref)
pos, sign = E.findpeaks(prof, 2, mode="auto")
print("  positions python\n", np.round(pos, 5), "\n  IDL\n", np.round(fits.getdata(os.path.join(REF, "ps_positions.fits")), 5))
t0 = time.time()
tc, pts = E.tracespec2d(work, wc["wavecal"], wc["spatcal"], fl["edgecoeffs"], fl["xranges"],
                        pos, 5, 5, 3.0, 0.1, 1)
print(f"  trace {time.time() - t0:.1f} s")
stat("trace coefficients", tc, fits.getdata(os.path.join(REF, "ps_trace.fits")))
apsign = np.zeros(2, int)
for i, o in enumerate(orders):
    x, p = prof[o]
    med = E.I.median(p)
    for j in range(2):
        apsign[j] = 1 if p[E.I.itrunc(E.I.tabinv(x, pos[i, j]))] > med else -1
print("  apsign", apsign)
t0 = time.time()
spec = E.extspec2d(work, var, fl["omask"], orders, wc["wavecal"], wc["spatcal"], fl["ps"],
                   apsign, tc, np.ones((len(orders), 2)), psbginfo=[2.5, 5.0], bgdeg=1)
print(f"  extraction {time.time() - t0:.1f} s")
k = 0
for o in orders:
    for a in (1, 2):
        ref = fits.getdata(os.path.join(REF, f"ps_spec_{k}.fits"))
        stat(f"PS order {o} ap {a} (w,f,e)", spec[(o, a)], ref)
        k += 1

# ---------------------------------------------------------------- extended source
f = [os.path.join(RAW, "SPEC_SO2025D-006-0120.fits"), os.path.join(RAW, "SPEC_SO2025D-006-0121.fits")]
work, var, hdrs, itot = E.load_image(f, fl)
tcx = np.array([[20.1, 0.0]] * len(orders))
spec = E.extspec2d(work, var, fl["omask"], orders, wc["wavecal"], wc["spatcal"], fl["ps"],
                   np.array([1]), tcx, np.ones((len(orders), 1)), xsbginfo=None, bgdeg=None)
for k, o in enumerate(orders):
    ref = fits.getdata(os.path.join(REF, f"xs_spec_{k}.fits"))
    stat(f"XS order {o} (w,f,e)", spec[(o, 1)], ref)
