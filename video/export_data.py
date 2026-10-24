"""
Film 1's data: real arrays and images from the tutorial night, made with the
package's own routines, for the film's animated plots.

    python video/export_data.py          (after docs/make_shots.py has made reduced/tutorial)

Writes video/film/public/explainer/data.json and its images.
"""

import json
import os
import sys

import numpy as np
from astropy.io import fits
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from triplespec_reducer import flat as F, wavecal as W, extract as E, telluric as T   # noqa: E402
from triplespec_reducer.frames import read_ts4                                         # noqa: E402
from triplespec_reducer.instrument import read_atran                                   # noqa: E402

RED = os.path.join(ROOT, "reduced", "tutorial")
RAW = os.path.join(ROOT, "data", "raw", "2025-11-02")
OUT = os.path.join(ROOT, "video", "film", "public", "explainer")
os.makedirs(OUT, exist_ok=True)
SCI, STD = "NGC134", "HD225200_1"


def raw(n):
    return os.path.join(RAW, f"SPEC_SO2025D-006-{n:04d}.fits")


def thin(x, y, n=1400):
    x, y = np.asarray(x, float), np.asarray(y, float)
    g = np.isfinite(x) & np.isfinite(y)
    x, y = x[g], y[g]
    if x.size > n:
        k = np.linspace(0, x.size - 1, n).astype(int)
        x, y = x[k], y[k]
    return [[round(float(a), 6), float(f"{b:.6g}")] for a, b in zip(x, y)]


def spec(name):
    with fits.open(os.path.join(RED, "proc", name)) as h:
        d = h[0].data
        return d.astype(d.dtype.newbyteorder("=")), h[0].header


def save_img(a, name, lo, hi, cmap="gray", w=None):
    import matplotlib
    m = matplotlib.colormaps[cmap]
    z = np.clip((np.nan_to_num(a, nan=lo) - lo) / (hi - lo), 0, 1)
    rgb = (m(z)[..., :3] * 255).astype(np.uint8)[::-1]          # origin at the bottom
    im = Image.fromarray(rgb)
    if w:
        im = im.resize((w, int(im.height * w / im.width)), Image.LANCZOS)
    im.save(os.path.join(OUT, name))
    return [im.width, im.height]


D = {}
fl = F.read_flat(os.path.join(RED, "cal", "flat.fits"))
wc = W.read_wavecal(os.path.join(RED, "cal", f"{SCI}_wavecal.fits"), rotation=fl["rotation"])
aw, at = read_atran(fl["rp"])
fr, _, _, _ = read_ts4([raw(116), raw(117)], rotation=fl["rotation"])
diff = (fr[0] - fr[1]) / fl["flat"]

# the raw frame with its order edges
lo, hi = np.nanpercentile(fr[0], [5, 99.5])
D["raw"] = dict(size=save_img(fr[0], "raw.png", lo, hi, w=1400), ny=fr[0].shape[0], nx=fr[0].shape[1],
                edges={str(o): [[float(np.polyval(fl["edgecoeffs"][i, j][::-1], x)) for x in range(0, 2048, 32)]
                                for j in (0, 1)] for i, o in enumerate(fl["orders"])})
# A, B and A-B, detail
sl = (slice(370, 560), slice(900, 1300))
lo, hi = np.nanpercentile(fr[0][sl], [2, 99.5])
save_img(fr[0][sl], "A.png", lo, hi)
save_img(fr[1][sl], "B.png", lo, hi)
v = np.nanpercentile(np.abs(diff[sl]), 99)
save_img(diff[sl], "AB.png", -v, v, cmap="RdBu_r")
# the flat: a lamp frame and the normalised flat
lamp, _, _, _ = read_ts4([os.path.join(RAW, "SPEC_DFLAT_0002.fits")], rotation=fl["rotation"])
sl2 = (slice(150, 560), slice(800, 1250))
lo, hi = np.nanpercentile(lamp[0][sl2], [2, 99.5])
save_img(lamp[0][sl2], "lamp.png", lo, hi)
save_img(fl["flat"][sl2], "flat.png", 0.9, 1.1)
# the wavelength image
D["wave_img"] = dict(size=save_img(wc["wavecal"], "wavecal.png", 0.92, 2.5, cmap="viridis", w=1400))
# the atmosphere and the orders
D["atran"] = thin(aw, at, 2400)
comb_sci, h_sci = spec(f"{SCI}_comb.fits")
orders = [int(t) for t in str(h_sci["ORDERS"]).split(",")]
D["orders"] = {str(o): [float(np.nanmin(comb_sci[k, 0])), float(np.nanmax(comb_sci[k, 0]))]
               for k, o in enumerate(orders)}
# the sky: OH lines in order 4
sky, _, _, _ = read_ts4([raw(121)], rotation=fl["rotation"])
sky = sky[0] / fl["flat"]
m = (fl["omask"] == 4) & (wc["spatcal"] > 14) & (wc["spatcal"] < 24)
wl, fx = [], []
for c in np.unique(np.nonzero(m)[1]):
    r = m[:, c]
    wl.append(np.nanmedian(wc["wavecal"][r, c]))
    fx.append(np.nanmedian(sky[r, c]))
D["oh"] = thin(wl, fx, 2000)
# the spatial profile of the standard, order 5, with apertures and background
prof = E.mkspatprof2d(diff, fl["omask"], wc["wavecal"], wc["spatcal"], fl["slith_arc"] / fl["slith_pix"],
                      1, aw, at, 0.4)
pos, signs = E.findpeaks(prof, 2)
y, p = prof[5]
D["profile"] = dict(curve=thin(y, p, 400), peaks=[float(x) for x in pos[fl["orders"].index(5)]],
                    signs=[int(s) for s in signs], radius=1.0, bg=[2.5, 7.5])
# combining: the four spectra of the standard and their combination, a stretch of order 5
comb_std, h_std = spec(f"{STD}_comb.fits")
k5 = orders.index(5)
D["combine"] = dict(spectra=[], combined=None)
for f in sorted(x for x in os.listdir(os.path.join(RED, "proc")) if x.startswith(f"{STD}_spec")):
    d, _ = spec(f)
    z = (d[k5, 0] > 1.20) & (d[k5, 0] < 1.27)
    D["combine"]["spectra"].append(thin(d[k5, 0][z], d[k5, 1][z], 500))
z = (comb_std[k5, 0] > 1.20) & (comb_std[k5, 0] < 1.27)
D["combine"]["combined"] = thin(comb_std[k5, 0][z], comb_std[k5, 1][z], 500)
# the kernel on Pa delta
st = T.load(os.path.join(RED, "proc", f"{STD}_comb.fits"), os.path.join(RED, "proc", f"{SCI}_comb.fits"), 6.40, 6.39)
T.conkernel(st, order=6, normreg=[(0.990, 0.997), (1.013, 1.020)], normdegree=3, linereg=(0.997, 1.013))
kp = st["kernel_plot"]
D["kernel"] = dict(star=thin(kp["w"], kp["data"]), vega=thin(kp["w"], kp["vega"]), conv=thin(kp["w"], kp["conv"]),
                   resid=thin(kp["w"], 100 * kp["resid"]), maxdev=float(st["maxdev"]), rmsdev=float(st["rmsdev"]),
                   vshift=float(st["vshift"]))
# the telluric correction, order 4
tell, _ = spec(f"{SCI}_tellspec.fits")
tc, _ = spec(f"{SCI}_tc.fits")
k4 = orders.index(4)
z = (comb_std[k4, 0] > 1.45) & (comb_std[k4, 0] < 1.80)


def sl_(a):
    return a[z] if a.shape == z.shape else a


D["tell"] = dict(std=thin(comb_std[k4, 0][z], comb_std[k4, 1][z]),
                 tell=thin(tell[k4, 0][z], tell[k4, 1][z]),
                 before=thin(comb_sci[k4, 0][z], comb_sci[k4, 1][z] / np.nanmedian(comb_sci[k4, 1][z])),
                 after=thin(tc[k4, 0][z], tc[k4, 1][z] / np.nanmedian(tc[k4, 1][z])))
# merging
mg, _ = spec(f"{SCI}_merged.fits")
mg = mg if mg.ndim == 2 else mg[0]
# in units of the merged spectrum's median, so the plot reads in plain numbers
mm = np.nanmedian(mg[1])
# the opaque bands between J, H and K left out: there the corrected flux is only noise
def clear(w, f):
    w, f = np.asarray(w, float), np.array(f, float)
    f[((w > 1.35) & (w < 1.44)) | ((w > 1.81) & (w < 1.94))] = np.nan
    return w, f
D["merge"] = dict(orders={str(o): thin(*clear(tc[k, 0], tc[k, 1] / mm), 700) for k, o in enumerate(orders)},
                  merged=thin(*clear(mg[0], mg[1] / mm), 2000))
json.dump(D, open(os.path.join(OUT, "data.json"), "w"))
print("written", OUT, f"{os.path.getsize(os.path.join(OUT, 'data.json')) / 1e6:.1f} MB")
