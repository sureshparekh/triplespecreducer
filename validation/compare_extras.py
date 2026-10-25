"""
The options beyond the default path, compared with IDL
(validation/idl/ref_extras.pro): combined images, A-mode extraction of a
combined image, xtellcor units, edited H-line scales, combine apertures.
"""

import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from triplespec_reducer import (flat as F, wavecal as W, extract as E, combineimgs as CI,   # noqa: E402
                                telluric as T, combine as C, idl as I)
from triplespec_reducer.instrument import read_atran                                         # noqa: E402

RAW = os.path.join(os.path.dirname(HERE), "data", "raw", "2025-11-02")
REF = os.path.join(HERE, "idl_out")
OUT = os.path.join(HERE, "py_out")


def stat(name, a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.shape != b.shape:
        print(f"  {name:32s} SHAPE {a.shape} vs IDL {b.shape}")
        return
    m = np.isfinite(a) & np.isfinite(b) & (np.abs(b) > 0)
    rel = np.abs(a - b)[m] / np.abs(b[m])
    print(f"  {name:32s} rel. diff median {np.median(rel):.1e}  99% {np.percentile(rel, 99):.1e}  "
          f"max {rel.max():.1e}  finite-mask {'same' if np.array_equal(np.isfinite(a), np.isfinite(b)) else 'DIFFERS'}")


fl = F.read_flat(os.path.join(REF, "flat.fits"))
wc = W.read_wavecal(os.path.join(REF, "wavecal_NGC134.fits"), rotation=fl["rotation"])
files = [os.path.join(RAW, f"SPEC_SO2025D-006-{n:04d}.fits") for n in (116, 117, 118, 119)]
p = CI.combine_images(files, os.path.join(OUT, "HD225200_combined.fits"), mode="A-B",
                      abba2abab=True, statistic="Robust Mean (Std Error)", thresh=8.0)
with fits.open(p) as h:
    stat("combined image", h[1].data, fits.getdata(os.path.join(REF, "ci_mean.fits")))
    stat("combined variance", h[2].data, fits.getdata(os.path.join(REF, "ci_var.fits")))
    print("   ITOT", h[0].header["ITOT"], "(IDL 10.0)")

pb = CI.combine_images(files, os.path.join(OUT, "HD225200_combined_bg.fits"), mode="A-B",
                       abba2abab=True, statistic="Robust Mean (Std Error)", thresh=8.0,
                       subtract_sky=True, flatinfo=fl)
from triplespec_reducer.frames import rotate as _rot
with fits.open(pb) as h:
    stat("combined image, sky subtracted", _rot(h[1].data, fl["rotation"]),
         fits.getdata(os.path.join(REF, "ci_bg_mean.fits")))
    stat("its variance", _rot(h[2].data, fl["rotation"]), fits.getdata(os.path.join(REF, "ci_bg_var.fits")))

r = E.reduce_pair([], fl, wc, [os.path.join(OUT, "HD225200_ci_spec.fits")], source="point",
                  atran=read_atran(fl["rp"]), mode="A", combined=p,
                  optimal=False)              # ref_extras.pro: sum with background
print("   positions\n", np.round(r["positions"], 5), "\n   signs", r["signs"],
      "(IDL CI_POS 9.87706 28.0891 ... signs 1 -1)")
k = 0
for o in fl["orders"]:
    for a in (1, 2):
        stat(f"A-mode spectrum order {o} ap {a}", r["spectra"][(o, a)],
             fits.getdata(os.path.join(REF, f"ci_spec_{k}.fits")))
        k += 1

t = fits.getdata(os.path.join(REF, "tc_tellspec.fits")).astype(np.float32)   # READFITS keeps FLOAT
jy, ejy = T.chfunits(t[1, 0], t[1, 1], 0, 1, 4, error=t[1, 2])
stat("telluric in Jy (f, e)", np.array([jy, ejy]), fits.getdata(os.path.join(REF, "un_jy.fits")))
wm, _ = T.chfunits(t[1, 0], t[1, 1], 0, 1, 0)
stat("telluric in W m-2 um-1", wm, fits.getdata(os.path.join(REF, "un_wm.fits")))

st = T.load(os.path.join(REF, "comb_HD225200.fits"), os.path.join(REF, "comb_NGC134.fits"), 6.33, 6.31)
st["scale"], st["vshift"] = np.float32(0.74261242), np.float32(33.111367)   # FLOAT literals in ref_extras.pro
T.initial_points(st)
cp = st["cpoints"][4]
k11 = np.argmin(np.abs(cp[:, 0] - 1.6811))
cp[k11, 1] = np.float32(0.8) * st["scale"]
stat("H-line control points", cp, fits.getdata(os.path.join(REF, "hl_cpoints.fits")))
sc = T.scale_spectrum(st, 4)
stat("tension-spline scales", sc, fits.getdata(os.path.join(REF, "hl_scales.fits")))
kern = fits.getdata(os.path.join(REF, "tc_kern_4.fits")).astype(np.float32)
s = st["std"]["spec"][1]
tc, te, _ = T.mktellspec(s[0], s[1], s[2], st["vmag"], st["bmag"] - st["vmag"], kern, sc, st["wvega"], st["fvega"],
                         st["cfvega"], st["cf2vega"], st["vshift"])
stat("telluric with edited scales", np.array([tc, te]), fits.getdata(os.path.join(REF, "hl_tell4.fits")))

# xmc_scalelines 'estimate scale' around Br 11, from the untouched control points
st["kernels"] = {o: kern for o in st["stdorders"]}
T.initial_points(st)
est = T.estimate_scale(st, 4, (1.6770, 1.6850))
idl_est = [l.split()[1:] for l in open(os.path.join(REF, "ref_extras.log")) if l.startswith("ESTSCALE")][0]
print("   estimate scale", repr(float(est[0])), "(IDL", idl_est[0] + ")",
      "identical" if float(est[0]) == float(idl_est[0]) else "DIFFERS")
stat("estimate-scale Lorentzian", st["estimate_fit"], fits.getdata(os.path.join(REF, "hl_est_gcoeff.fits")))

# xmergeorders with auto scale, trims and a gap (see ref_extras.pro)
from triplespec_reducer import merge as M
cs = fits.getdata(os.path.join(REF, "tc_NGC134.fits"))
cs = cs.astype(cs.dtype.newbyteorder("="))
m, used, done = M.merge_orders(cs, [3, 4, 5, 6, 7], anchor=4, sequence=[5, 3, 6],
                               scales={5: "auto"}, regions={5: (1.42, 1.47)},
                               trims={5: ("Green", "Right", 1.47), 6: ("Blue", "Right", 1.30)})
idl_ms = [l.split()[1] for l in open(os.path.join(REF, "ref_extras.log")) if l.startswith("MSCALE")][0]
print("   auto scale", repr(float(used[5])), "(IDL", idl_ms + ")",
      "identical" if np.float32(used[5]) == np.float32(float(idl_ms)) else "DIFFERS")
stat("merged with trims and scale", m, fits.getdata(os.path.join(REF, "m_trims.fits")))

# xcleanspec: fix a region, then remove one (cursor order reversed)
from triplespec_reducer import clean as CLN
mspec = fits.getdata(os.path.join(REF, "m_NGC134.fits"))
mspec = mspec.astype(mspec.dtype.newbyteorder("="))[None]
ed = CLN.fix_region(mspec, 0, 1.27, 1.29)
ed = CLN.remove_region(ed, 0, 1.82, 1.80)
stat("xcleanspec fix + remove", ed[0], fits.getdata(os.path.join(REF, "cl_edit.fits")))

# xcombspec statistics 1-8, order 5, unscaled
L5 = C.load([os.path.join(REF, f"spectra{n:04d}.fits") for n in (116, 117, 118, 119)])
rs = fits.getdata(os.path.join(REF, "cs_stats.fits"))          # [9, 2, npix]
for k, name in enumerate(C.STATISTICS[1:], start=1):
    comb = C.combine(L5, statistic=name)
    stat(f"combine: {name}", comb[2, 1:3], rs[k])

# mc_findpeaks with GUESS and FIXED positions
wk = fits.getdata(os.path.join(REF, "ps_work.fits")).astype(np.float32)
aw_, at_ = read_atran(fl["rp"])
pr = E.mkspatprof2d(wk, fl["omask"], wc["wavecal"], wc["spatcal"], fl["slith_arc"] / fl["slith_pix"], 1, aw_, at_, 0.4)
gpos = [[9.5, 28.3]] * len(fl["orders"])
pg, sg = E.findpeaks(pr, 2, mode="guess", guess=gpos)
pf, sf_ = E.findpeaks(pr, 2, mode="fixed", fixed=gpos)
stat("findpeaks GUESS positions", pg, fits.getdata(os.path.join(REF, "fp_guess.fits")))
stat("findpeaks FIXED positions", pf, fits.getdata(os.path.join(REF, "fp_fixed.fits")))
print("   signs guess/fixed", sg, sf_, "(IDL", [l.split()[1:] for l in open(os.path.join(REF, "ref_extras.log")) if l.startswith("FPSIGN")][0], ")")

# A-Sky: raw 120 minus the combined image as a sky
wa, va, _, _ = E.load_image([os.path.join(RAW, "SPEC_SO2025D-006-0120.fits")], fl, mode="A-Sky",
                            sky=os.path.join(OUT, "HD225200_combined.fits"))
stat("A-Sky work image", wa, fits.getdata(os.path.join(REF, "as_work.fits")))
sa = E.extspec2d(wa, va, fl["omask"], fl["orders"], wc["wavecal"], wc["spatcal"], fl["ps"],
                 np.array([1]), np.array([[20.1, 0.0]] * 5), np.ones((5, 1)))
for k, o in enumerate(fl["orders"]):
    stat(f"A-Sky extended spectrum {o}", sa[(o, 1)], fits.getdata(os.path.join(REF, f"as_spec_{k}.fits")))

L = C.load([os.path.join(OUT, "HD225200_ci_spec.fits")] * 2, combine_apertures=True)
print("   combine apertures: one aperture out of", L["spc"][0][0].shape[0], "spectra, ITOT", L["itot"])
