"""
Run the same inputs through IDL 7.1 and triplespec_reducer.idl, and report
the largest difference for each primitive.

    python validation/primitives_parity.py
"""

import os
import subprocess
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from triplespec_reducer import idl as I  # noqa: E402

WORK = os.path.join(HERE, "idl_out", "prims")
os.makedirs(WORK, exist_ok=True)
rng = np.random.default_rng(42)


def w(name, arr):
    fits.writeto(os.path.join(WORK, name + ".fits"), np.asarray(arr), overwrite=True)


def r(name):
    return fits.getdata(os.path.join(WORK, name + ".fits")).astype(float)


# ---- inputs (numpy order; IDL reads FITS as [x, y] = transposed view)
v = rng.normal(size=40).astype(np.float32)
v[[3, 17]] = np.nan
w("v", v)
img = rng.normal(size=(30, 45)).astype(np.float32)     # [row, col]
img[5, 7] = np.nan
w("img", img)
cube = rng.normal(5, 1, size=(7, 20, 25)).astype(np.float32)   # [k, row, col]
cube[2, 3, 4] = 40.0
w("cube", cube)
w("cubevar", (np.abs(cube) * 0.1 + 0.05).astype(np.float32))
ker_odd = np.array([0.1, 0.25, 0.4, 0.2, 0.05], np.float32)
ker_even = np.array([0.1, 0.3, 0.4, 0.2], np.float32)
sig = np.cumsum(rng.normal(size=200)).astype(np.float32)
w("sig", sig); w("kodd", ker_odd); w("keven", ker_even)
xi = (rng.uniform(-1.5, 45.5, 300)).astype(np.float32)
yi = (rng.uniform(-1.5, 30.5, 300)).astype(np.float32)
w("xi", xi); w("yi", yi)
xt = np.sort(rng.uniform(0, 10, 25)).astype(np.float32)
yt = np.sin(xt).astype(np.float32)
xq = np.linspace(-1, 11, 77).astype(np.float32)
w("xt", xt); w("yt", yt); w("xq", xq)
gx = np.arange(31, dtype=np.float32)
gy = (50 * np.exp(-0.5 * ((gx - 14.3) / 2.1) ** 2) + 3 + rng.normal(0, 0.5, 31)).astype(np.float32)
w("gx", gx); w("gy", gy)
px = rng.uniform(0, 100, 60); py = rng.uniform(0, 3, 60)
pz = 1 + 0.01 * px + 0.3 * py + 0.001 * px * py + rng.normal(0, 0.01, 60)
pz[5] += 1.0
w("px", px); w("py", py); w("pz", pz)
sx = rng.uniform(0, 20, 400); sy = rng.uniform(0, 10, 400)
sv = np.sin(sx / 3) + sy ** 1.3
w("sx", sx); w("sy", sy); w("sv", sv)
polyx = np.array([2.3, 5.7, 6.1, 2.9]); polyy = np.array([1.2, 1.6, 4.4, 3.9])
w("polyx", polyx); w("polyy", polyy)

script = f"""
!quiet=1
sp = os.environ.get('SPEXTOOL_DIR', 'Spextool4.1.0')
!path = expand_path('+'+sp) + ':' + !path
d = '{WORK}/'
v = readfits(d+'v.fits',/SILENT) & img = readfits(d+'img.fits',/SILENT)
cube = readfits(d+'cube.fits',/SILENT) & cvar = readfits(d+'cubevar.fits',/SILENT)
writefits, d+'o_med.fits', [median(v), median(v,/EVEN), median(v[0:9]), median(v[0:9],/EVEN)]
writefits, d+'o_medfilt.fits', median(img,5)
writefits, d+'o_meddim.fits', median(img,DIMENSION=2)
writefits, d+'o_meddimE.fits', median(img,DIMENSION=1,/EVEN)
sig = readfits(d+'sig.fits',/SILENT) & ko=readfits(d+'kodd.fits',/SILENT) & ke=readfits(d+'keven.fits',/SILENT)
writefits, d+'o_convodd.fits', convol(sig,ko)
writefits, d+'o_conveven.fits', convol(sig,ke)
writefits, d+'o_convscale.fits', convol(sig,ko,total(ko))
xi=readfits(d+'xi.fits',/SILENT) & yi=readfits(d+'yi.fits',/SILENT)
im2 = img & im2[where(finite(im2) eq 0)] = 0.0
writefits, d+'o_cubic.fits', interpolate(im2,xi,yi,CUBIC=-0.5)
writefits, d+'o_cubicmiss.fits', interpolate(im2,xi,yi,CUBIC=-0.5,MISSING=!values.f_nan)
writefits, d+'o_lin2.fits', interpolate(im2,xi,yi)
xt=readfits(d+'xt.fits',/SILENT) & yt=readfits(d+'yt.fits',/SILENT) & xq=readfits(d+'xq.fits',/SILENT)
writefits, d+'o_interpol.fits', interpol(yt,xt,xq)
linterp, xt, yt, xq, yl & writefits, d+'o_linterp.fits', yl
linterp, xt, yt, xq, yl2, MISSING=-99 & writefits, d+'o_linterpm.fits', yl2
tabinv, xt, xq, ie & writefits, d+'o_tabinv.fits', ie
writefits, d+'o_spline.fits', spline(xt,yt,xq[5:70],10.)
writefits, d+'o_ccor.fits', c_correlate(yt,shift(yt,3),indgen(11)-5)
gx=readfits(d+'gx.fits',/SILENT) & gy=readfits(d+'gy.fits',/SILENT)
junk = gaussfit(gx,gy,a,NTERMS=4) & writefits, d+'o_gfit.fits', a
junk = mpfitpeak(double(gx),double(gy),b,NTERMS=4) & writefits, d+'o_mpfp.fits', b
junk = mpfitpeak(double(gx),double(gy),b2,NTERMS=4,ESTIMATES=[40d,13d,1.5d,0d]) & writefits, d+'o_mpfpe.fits', b2
px=readfits(d+'px.fits',/SILENT) & py=readfits(d+'py.fits',/SILENT) & pz=readfits(d+'pz.fits',/SILENT)
c1 = poly_fit1d(px,pz,2,/SILENT) & writefits, d+'o_pf1.fits', c1
c2 = robustpoly1d(px,pz,2,3,0.01,/SILENT,OGOODBAD=g1) & writefits, d+'o_rp1.fits', [c2,g1]
c3 = robustpoly2d(px,py,pz,2,1,3,0.01,/SILENT,OGOODBAD=g2) & writefits, d+'o_rp2.fits', [c3,g2]
mc_meancomb, cube, m1, v1 & writefits, d+'o_mc1.fits', [[[m1]],[[v1]]]
mc_meancomb, cube, m2, v2, DATAVAR=cvar & writefits, d+'o_mc2.fits', [[[m2]],[[v2]]]
mc_meancomb, cube, m3, v3, ROBUST=3 & writefits, d+'o_mc3.fits', [[[m3]],[[v3]]]
mc_moments, pz, mn, vr, sd, ROBUST=5, /SILENT & writefits, d+'o_mom.fits', [mn,vr,sd]
sx=readfits(d+'sx.fits',/SILENT) & sy=readfits(d+'sy.fits',/SILENT) & sv=readfits(d+'sv.fits',/SILENT)
triangulate, sx, sy, tri
g = trigrid(sx,sy,sv,tri,[0.37,0.21],[min(sx),min(sy),max(sx),max(sy)],MISSING=!values.f_nan,XGRID=xg,YGRID=yg)
writefits, d+'o_trigrid.fits', g & writefits, d+'o_trixg.fits', xg & writefits, d+'o_triyg.fits', yg
polyx=readfits(d+'polyx.fits',/SILENT) & polyy=readfits(d+'polyy.fits',/SILENT)
inds = polyfillaa(polyx,polyy,10,10,AREAS=areas)
writefits, d+'o_pfa.fits', [[double(inds)],[areas]]
writefits, d+'o_round.fits', round([0.5,1.5,2.5,-0.5,-2.5])
exit
"""
with open(os.path.join(WORK, "prims.pro"), "w") as fh:
    fh.write(script)
res = subprocess.run(["/usr/local/bin/idl", "-e", "@" + os.path.join(WORK, "prims.pro")],
                     capture_output=True, text=True, cwd=WORK)
errs = [l for l in (res.stdout + res.stderr).splitlines() if l.startswith("%")]
if errs:
    print("\n".join(errs[:20]))


def report(name, ours, theirs, tol=1e-5):
    ours = np.asarray(ours, float)
    theirs = np.asarray(theirs, float)
    if ours.shape != theirs.shape:
        print(f"  {name:14s} SHAPE {ours.shape} vs IDL {theirs.shape}")
        return
    nan_ok = np.array_equal(np.isnan(ours), np.isnan(theirs))
    d = np.nanmax(np.abs(ours - theirs)) if np.isfinite(theirs).any() else 0.0
    scale = max(np.nanmax(np.abs(theirs)), 1e-30)
    ok = nan_ok and d <= tol * max(scale, 1.0)
    print(f"  {name:14s} max|diff| {d:10.3e}  (scale {scale:9.3e})  "
          f"NaNs {'same' if nan_ok else 'DIFFER'}  {'OK' if ok else 'CHECK'}")


print("IDL parity of primitives")
report("median", [I.median(v), I.median(v, True), I.median(v[:10]), I.median(v[:10], True)], r("o_med"))
report("median_filt", I.median_filter(img, 5), r("o_medfilt"))
report("median dim2", I.median(img, axis=0), r("o_meddim"))
report("median dim1E", I.median(img, even=True, axis=1), r("o_meddimE"))
report("convol odd", I.convol(sig, ker_odd), r("o_convodd"))
report("convol even", I.convol(sig, ker_even), r("o_conveven"))
report("convol scale", I.convol(sig, ker_odd, ker_odd.sum()), r("o_convscale"))
im2 = np.nan_to_num(img)
report("cubic", I.interpolate2d(im2, xi, yi, cubic=-0.5), r("o_cubic"))
report("cubic miss", I.interpolate2d(im2, xi, yi, cubic=-0.5, missing=np.nan), r("o_cubicmiss"))
report("bilinear", I.interpolate2d(im2, xi, yi), r("o_lin2"))
report("interpol", I.interpol(yt, xt, xq), r("o_interpol"))
report("linterp", I.linterp(xt, yt, xq), r("o_linterp"))
report("linterp miss", I.linterp(xt, yt, xq, missing=-99), r("o_linterpm"))
report("tabinv", I.tabinv(xt, xq), r("o_tabinv"))
report("spline", I.spline(xt, yt, xq[5:71], 10.0)[(xq[5:71]>xt[0])&(xq[5:71]<xt[-1])], r("o_spline")[(xq[5:71]>xt[0])&(xq[5:71]<xt[-1])])
report("c_correlate", I.c_correlate(yt, np.roll(yt, 3), np.arange(11) - 5), r("o_ccor"))
report("gaussfit", I.gaussfit(gx, gy, nterms=4)[1], r("o_gfit"), tol=1e-4)
report("mpfitpeak", I.mpfitpeak(gx.astype(float), gy.astype(float), nterms=4)[1], r("o_mpfp"), tol=1e-6)
report("mpfitpeak est", I.mpfitpeak(gx.astype(float), gy.astype(float), nterms=4, estimates=[40., 13., 1.5, 0.])[1], r("o_mpfpe"), tol=1e-6)
report("poly_fit1d", I.poly_fit1d(px, pz, 2), r("o_pf1"), tol=1e-8)
c, g = I.robustpoly1d(px, pz, 2, 3, 0.01)
report("robustpoly1d", np.concatenate([c, g]), r("o_rp1"), tol=1e-8)
c, g = I.robustpoly2d(px, py, pz, 2, 1, 3, 0.01)
report("robustpoly2d", np.concatenate([c, g]), r("o_rp2"), tol=1e-8)
m, vv, _ = I.meancomb(cube)
report("meancomb", np.stack([m, vv]), r("o_mc1"))
m, vv, _ = I.meancomb(cube, datavar=fits.getdata(os.path.join(WORK, "cubevar.fits")))
report("meancomb wt", np.stack([m, vv]), r("o_mc2"))
m, vv, _ = I.meancomb(cube, robust=3)
report("meancomb rob", np.stack([m, vv]), r("o_mc3"))
mn, vr, sd, _ = I.moments(pz, robust=5)
report("moments", [mn, vr, sd], r("o_mom"), tol=1e-8)
xg, yg = r("o_trixg"), r("o_triyg")
ourx, oury = I.trigrid_axes(sx.min(), sx.max(), I.FL(0.37), sy.min(), sy.max(), I.FL(0.21))   # IDL GS=[0.37,0.21] is FLOAT
report("trigrid xgrid", ourx, xg)
report("trigrid ygrid", oury, yg)
(gg,), _ = I.triangulate_grid(sx, sy, [sv], xg, yg)
report("trigrid", gg, r("o_trigrid"), tol=1e-4)
ci, ri, ar = I.polyfillaa(polyx, polyy, 10, 10)
ours = np.column_stack([ri * 10 + ci, ar]).T
report("polyfillaa", ours, r("o_pfa"), tol=1e-6)
report("round", I.iround([0.5, 1.5, 2.5, -0.5, -2.5]), r("o_round"))
