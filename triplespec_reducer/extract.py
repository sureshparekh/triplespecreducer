"""
Spectral extraction: what xspextool does once the calibrations exist.

    load the pair            XSPEXTOOL_LOADIMAGES   (A-B)/flat, var/flat^2
    spatial profiles         MC_MKSPATPROF2D        median along lambda of
                                                    the rectified order
    aperture positions       MC_FINDPEAKS           /AUTO, GUESS or FIXED
    trace                    MC_TRACESPEC2D         (point sources), or the
                                                    positions held fixed
    aperture signs           XSPEXTOOL_PSDEFINEAP
    extraction               MC_EXTSPEC2D           exact-area resampling of
                                                    every (lambda, slit) cell
    the file                 MC_WRITESPEC

Point sources are extracted optimally by default, as the TS4_SOAR panel
does: ts4_soar.dat greys out the 'Optimal' box (OPTIMALEXTRACTION=0 0),
but the startup 'PS BG Subtraction' event (PSBGSUB=1) ticks it, so
XSPEXTOOL_PSEXTRACTSPEC hands MC_EXTSPEC2D the PSF radius and the spatial
map. Unticking BG subtraction gives a plain sum with no background.

Extended sources: in the IDL, MC_EXTSPEC2D asks MC_MKAPMASK1D for the
background regions as XSBGINFO=xsbgino (a misspelling of xsbginfo), and
XSPEXTOOL_XSEXTRACTSPEC passes BGDEG=bgorder where the degree was read into
bgfit. Either alone switches the background fit off, so in 2D mode the
regions a user types for an extended source are drawn but never used. That
is reproduced by default. `xs_background=True` applies them exactly as the
point-source path does -- the subtraction Spextool was written to make.
"""

import math
import os

import numpy as np
from astropy.io import fits

from triplespec_reducer import idl as I
from triplespec_reducer.frames import read_ts4
from triplespec_reducer.instrument import fsextract_index

try:
    from numba import njit
except ImportError:                    # the loops still run, only slower
    def njit(*a, **k):
        def wrap(f):
            return f
        return wrap if not (a and callable(a[0])) else a[0]


# ------------------------------------------------------------- loading

def load_image(files, flatinfo, mode="A-B", sky=None, flatfield=True, errorprop=True,
               combined=None):
    """
    XSPEXTOOL_LOADIMAGES: the work image and its variance.

    A-B differences the two frames (A/t - B/t). A takes one raw frame, or,
    with `combined`, an image XSPEXTOOL_COMBINEIMGS wrote. A-Sky subtracts
    a sky image: `sky` is the path of a combined image in cal/. The image is
    divided by the flat and the variance by its square.
    """
    rot = flatinfo["rotation"]
    from triplespec_reducer.combineimgs import read_combined
    if combined:
        data, var, flags, h, itot0 = read_combined(combined, rot)
        hdrs, itot = [h], [itot0]
    elif mode == "A-B":
        data, var, hdrs, itot = read_ts4(files, rotation=rot, pair=True)
        data, var = data[0], var[0]
    else:
        data, var, hdrs, itot = read_ts4(files[:1], rotation=rot)
        data, var = data[0], var[0]
    if mode == "A-Sky":
        sd, sv, _, sh, _ = read_combined(sky, rot)
        data = data - sd
        var = var + sv
        hdrs = list(hdrs) + [sh]
    if flatfield:
        data = data / flatinfo["flat"]
        var = var / I.sq(flatinfo["flat"])
    if not errorprop:
        var[:] = 1.0
    return data.astype(np.float32), var.astype(np.float32), hdrs, itot


def load_bitmask(files, flatinfo, mode="A-B", sky=None, combined=None, lincormax=80000):
    """
    The linearity bitmask XSPEXTOOL_LOADIMAGES keeps (KNA_READTS4FITS'
    BITINFO {lincormax, bit 0}): bit 0 set where a raw frame exceeds
    LINCORMAX (the LONG 80000), OR-ed over the pair; from a combined image
    its 3rd extension; A-Sky adds the sky's. Only the 1D extraction uses it.
    """
    from triplespec_reducer.combineimgs import read_combined
    from triplespec_reducer.frames import _read_one, rotate
    rot = flatinfo["rotation"]
    lim = np.float64(lincormax)                       # DOUBLE(bitinfo.lincormax)
    if combined:
        return read_combined(combined, rot)[2]
    n = 2 if mode == "A-B" else 1
    m = np.zeros(1, dtype=np.uint8)
    for f in files[:n]:
        img = _read_one(f)[0]
        m = m | (img > lim).astype(np.uint8)
    m = rotate(m, rot)
    if mode == "A-Sky":
        m = m | read_combined(sky, rot)[2]
    return m


def load_bdpxmask(rotation):
    """xspextool's bad pixel mask: BADPIXMASK (fakebmask.fits) as BYTE, rotated."""
    from triplespec_reducer._paths import INSTRUMENT_DIR
    from triplespec_reducer.frames import rotate
    m = fits.getdata(os.path.join(INSTRUMENT_DIR, "fakebmask.fits"))
    return rotate(m.astype(np.uint8), rotation)


def _strtrim(v):
    """STRTRIM(v,2) of a number, in IDL's default format: DOUBLE %#.8g, FLOAT %#.6g."""
    return f"{float(v):#.8g}" if not I.isf32(np.asarray(v)) else f"{float(v):#.6g}"


def write_trace(path, orders, naps, tracecoeffs, obsmode="XD"):
    """
    The trace file XSPEXTOOL_PSTRACESPEC writes to cal/: the selected
    orders, and each coefficient as STRTRIM prints it (a traced, DOUBLE
    fit to 8 significant figures; the FLTARR of Use Ap Positions to 6).
    """
    tracecoeffs = np.asarray(tracecoeffs)
    with open(path, "w") as fh:
        fh.write(f"OBSMODE={obsmode}\n")
        fh.write("ORDERS=" + ",".join(str(o) for o in orders) + "\n")
        fh.write(f"NAPS={naps}\n")
        k = 0
        for o in orders:
            for j in range(naps):
                fh.write(f"ORDER{o:02d}_AP{j + 1:02d}= " +
                         "  ".join(_strtrim(c) for c in tracecoeffs[k]) + "\n")
                k += 1


def read_trace(path):
    """
    MC_READTRACE: the orders the file lists, the number of apertures, and
    the coefficients as FLOAT() reads them, one row per order and aperture
    in the file's order.
    """
    lines = [ln for ln in open(path).read().splitlines() if "=" in ln]
    orders = fsextract_index(lines[1].split("=", 1)[1])
    naps = int(lines[2].split("=", 1)[1].strip())
    tc = [[np.float32(t) for t in ln.split("=", 1)[1].split()] for ln in lines[3:]]
    return orders, naps, np.array(tc, dtype=np.float32)


# ------------------------------------------------------------- profiles

_GRIDS = {}


def _order_grid(omask, wavecal, spatcal, order, ds):
    """
    The (lambda, slit) grid TRIGRID puts an order on, with pixel x, y there.

    It depends only on the calibration images and the step, never on the
    frame being extracted, so it is made once per wavelength solution and
    reused for every pair: the same arrays, so the same numbers, without
    triangulating 150 000 pixels again for each frame.
    """
    key = (id(omask), id(wavecal), id(spatcal), int(order), float(ds))
    hit = _GRIDS.get(key)
    if hit is not None and hit[0] is omask and hit[1] is wavecal and hit[2] is spatcal:
        return hit[3]
    nrows, ncols = omask.shape
    yy, xx = np.mgrid[0:nrows, 0:ncols]
    z = omask == order
    w, s = wavecal[z], spatcal[z]
    wmin, wmax = np.nanmin(w), np.nanmax(w)
    smin, smax = np.nanmin(s), np.nanmax(s)
    xmin, xmax = xx[z].min(), xx[z].max()
    dw = (wmax - wmin) / (xmax - xmin)
    xgrid, ygrid = I.trigrid_axes(wmin, wmax, dw, smin, smax, ds)
    (ix, iy), _ = I.triangulate_grid(w, s, [xx[z].astype(float), yy[z].astype(float)],
                                     xgrid, ygrid)
    out = (xgrid, ygrid, ix, iy, dw)     # ix, iy are [len(ygrid), len(xgrid)]
    if len(_GRIDS) > 64:
        _GRIDS.clear()
    _GRIDS[key] = (omask, wavecal, spatcal, out)
    return out


def mkspatprof2d(img, omask, wavecal, spatcal, ds, ybuffer, wtrans=None,
                 ttrans=None, atmosthresh=None, orders=None, coeffs=False):
    """
    MC_MKSPATPROF2D: the mean spatial profile of each order.

    The order is rectified onto (lambda, arcsec); each wavelength has its
    median along the slit removed; wavelengths where the atmosphere
    transmits less than `atmosthresh` are blanked; the profile is the
    median along wavelength, normalized by the sum of its absolute values.
    Returns {order: (ygrid, profile)}.

    With `coeffs`, also SPATCOEFFS and MEDSPATCOEFFS, the spatial maps
    optimal extraction rebuilds the profile from at every wavelength:
    returns (profiles, spatcoeffs, medspatcoeffs), the last two
    {order: (ygrid, c)} with c an FLTARR [ny, 5] (IDL [5, ny]).
    """
    orders = orders or sorted(int(o) for o in np.unique(omask) if o != 0)
    out, spat, med = {}, {}, {}
    for o in orders:
        xgrid, ygrid, ix, iy, _ = _order_grid(omask, wavecal, spatcal, o, ds)
        rect = I.interpolate2d(img, ix, iy, cubic=-0.5, missing=np.nan)
        ny = ygrid.size
        ygrid = ygrid[ybuffer:ny - ybuffer]
        if I.isf32(spatcal) and I.isf32(wavecal):
            # MC_SIMWAVECAL2D's FLOAT images (no wavecal file): TRIGRID's LIMITS
            # and so its YGRID are FLOAT (a DOUBLE wavecal makes them DOUBLE)
            ygrid = ygrid.astype(np.float32)
        rect = rect[ybuffer:ny - ybuffer, :]
        bg = I.median(rect, axis=0)
        rect = rect - bg[None, :]
        if wtrans is not None:
            mask = np.ones(xgrid.size, dtype=np.float32)            # FLTARR
            rtrans = I.linterp(wtrans, ttrans, xgrid, missing=1)
            # ATMOSTHRESH comes from its widget as a FLOAT
            mask[rtrans < np.float32(atmosthresh)] = np.nan
            rect = rect * mask[None, :]
        prof = I.median(rect, axis=1)
        a = np.abs(prof)
        prof = (prof / I.ftotal(a[np.isfinite(a)])).astype(np.float32)   # TOTAL(/NAN), float
        # IDL keeps [[ygrid],[prof]]: DOUBLE, so the slit grid keeps its precision.
        out[o] = (ygrid, prof)
        if coeffs:
            spat[o], med[o] = _spatcoeffs(xgrid, ygrid, rect, prof)
    return (out, spat, med) if coeffs else out


def _spatcoeffs(xgrid, ygrid, rect, prof, ndeg=4):
    """
    The SPATCOEFFS / MEDSPATCOEFFS part of MC_MKSPATPROF2D.

    The rectified order (atmosphere mask applied) is divided by the mean
    profile, each wavelength scaled by its median (/EVEN) along the slit,
    and every slit row then fitted along wavelength with a robust degree-4
    polynomial (5 sigma, GAUSSJ). A row the fit refuses gets NaNs, as
    IDL's CANCEL branch gives it. MEDSPATCOEFFS keep only the constant
    term: the mean profile itself.
    """
    ny = ygrid.size
    with np.errstate(divide="ignore", invalid="ignore"):
        norms = I.median(rect / prof[:, None], even=True, axis=0)     # / meanmap, DIMENSION=2
        rect = rect / norms[None, :]
    c = np.zeros((ny, ndeg + 1), dtype=np.float32)                   # FLTARR
    for j in range(ny):
        try:
            cj, _ = I.robustpoly1d(xgrid, rect[j], ndeg, 5, 0.1)
        except (ValueError, np.linalg.LinAlgError):
            cj = np.full(ndeg + 1, np.nan)
        c[j] = cj
    a = np.zeros((ny, ndeg + 1), dtype=np.float32)
    a[:, 0] = prof
    return (ygrid, c), (ygrid, a)


def _sinc_weights(xa, x, dampfac=3.25, ksize=21):
    """
    The part of MC_SINCINTERP that depends only on the abscissae: for every
    x, its TABINV index, and the damped-sinc kernel and the (edge-clamped)
    table indices it multiplies. In double, as the IDL computes it (FINDGEN
    minus a DOUBLE fraction; !PI is the FLOAT pi; EXP and SIN the libm
    double functions). Optimal extraction interpolates every wavelength of
    an order on the same slit grid, so these are made once per order.
    """
    import math
    idx = np.atleast_1d(I.tabinv(xa, x)).astype(np.float64)
    ix = np.trunc(idx).astype(np.int64)                              # FIX
    fx = idx - ix
    nf = np.asarray(xa).size
    pi32 = float(np.float32(np.pi))
    base = np.arange(ksize, dtype=np.float32) - np.float32(10)       # FINDGEN(ksize)-10
    w = np.zeros((idx.size, ksize))
    for p in np.where(fx != 0)[0]:
        for k in range(ksize):
            xk = float(base[k]) - float(fx[p])
            u1 = xk / dampfac
            u2 = pi32 * xk
            w[p, k] = math.exp(-(u1 * u1)) * math.sin(u2) / u2
    lobe = np.clip(np.arange(ksize)[None, :] - 10 + ix[:, None], 0, nf - 1)
    return ix, fx, w, lobe


def _sinc_apply(weights, ya):
    """MC_SINCINTERP's sums for one table of values: TOTAL(sinc*vals), in order."""
    ix, fx, w, lobe = weights
    f = np.asarray(ya)
    vals = f[lobe].astype(np.float32)                                # FLTARR(ksize)
    r = np.zeros(ix.size, dtype=w.dtype)                             # xx*0: TABINV's type
    for k in range(w.shape[1]):
        r = (r + w[:, k] * vals[:, k]).astype(w.dtype)
    z = fx == 0
    r[z] = f[ix[z]]
    return r


def sincinterp(xa, ya, x):
    """MC_SINCINTERP: damped-sinc interpolation of the table (xa, ya) at x."""
    return _sinc_apply(_sinc_weights(xa, x), ya)


def findpeaks(profiles, npeaks, doorders=None, mode="auto", guess=None, fixed=None):
    """
    MC_FINDPEAKS. `mode` is 'auto', 'guess' or 'fixed'; guess/fixed are
    positions in arcsec, one list of `npeaks` per order. Returns
    (positions, an array norders x npeaks; apsign, one sign per peak).
    """
    orders = list(profiles)
    no = len(orders)
    doorders = [1] * no if doorders is None else doorders
    positions = np.full((no, npeaks), np.nan, dtype=np.float32)    # REPLICATE(!VALUES.F_NAN)
    apsign = np.zeros((no, npeaks), dtype=int)
    for i, o in enumerate(orders):
        if not doorders[i]:
            continue
        # The IDL profile is [[ygrid],[prof]]: DOUBLE, unless the slit grid
        # is FLOAT (no wavecal file), and the fits are made in that type.
        P = np.float32 if I.isf32(np.asarray(profiles[o][0])) else np.float64
        x, p = np.asarray(profiles[o][0]).astype(P), np.asarray(profiles[o][1]).astype(P)
        if mode == "guess":
            y = p - I.median(p, even=True)
            for j in range(npeaks):
                xpix = I.tabinv(x, guess[i][j])
                _, fit = I.mpfitpeak(x, y, nterms=4, estimates=np.array(
                    [y[I.itrunc(xpix)], guess[i][j], 0.5, 0], dtype=P))
                fit = fit.astype(P)                   # MPFITPEAK's A: the data's type
                positions[i, j] = np.float32(fit[1])
                xpix = I.tabinv(x, fit[1])
                apsign[i, j] = 1 if y[I.itrunc(xpix)] > 0 else -1
                y = y - (fit[3] + I.gaussian(x, fit[0:3]))
        elif mode == "auto":
            y = np.abs(p - I.median(p, even=True))
            z = np.isfinite(y)
            xx, y = x[z], y[z]
            for j in range(npeaks):
                k = int(np.argmax(y))
                _, fit = I.mpfitpeak(xx, y, nterms=4,
                                     estimates=np.array([y[k], xx[k], 1, 0], dtype=P))
                fit = fit.astype(P)
                positions[i, j] = np.float32(fit[1])
                # IDL looks the sign up in the unfiltered profile with an
                # index found on the NaN-filtered abscissa; kept as written.
                xpix = I.tabinv(xx, fit[1])
                apsign[i, j] = 1 if p[I.itrunc(xpix)] > 0 else -1
                y = y - (fit[3] + I.gaussian(xx, fit[0:3]))
        else:
            for j in range(npeaks):
                positions[i, j] = np.float32(fixed[i][j])
                xpix = I.tabinv(x, fixed[i][j])
                apsign[i, j] = 1 if p[I.itrunc(xpix)] > 0 else -1
    for i in range(no):
        if doorders[i]:
            s = np.argsort(positions[i], kind="stable")
            positions[i] = positions[i][s]
            apsign[i] = apsign[i][s]
    cnt = int(np.sum(doorders))
    sign = np.sum(apsign, axis=0) / float(cnt) if no > 1 else apsign[0].astype(float)
    sign = np.where(sign >= 0, 1, -1)
    return positions, sign


def tracespec2d(image, wavecal, spatcal, edgecoeffs, xranges, appos, step,
                sumap, winthresh, sigthresh, fitorder):
    """
    MC_TRACESPEC2D: follow each aperture along its order.

    Every `step` columns a Gaussian is fitted across the slit (to a median
    of `sumap` columns) near the aperture; the centroids, converted to
    (lambda, arcsec), are fitted with a robust polynomial of slit position
    in wavelength. Returns coefficients [norders*naps, fitorder+1] and the
    traced points for display.
    """
    nrows, ncols = image.shape
    naps = appos.shape[1]
    half = I.itrunc(sumap / 2.)
    y = np.arange(nrows, dtype=np.float32)
    coeffs = np.zeros((len(xranges) * naps, fitorder + 1))
    points = []
    for i in range(len(xranges)):
        starts = xranges[i][0] + step - 1
        stops = xranges[i][1] - step + 1
        numstep = I.itrunc((stops - starts) / step) + 1
        column = np.arange(numstep, dtype=np.float32) * step + starts
        bot = I.poly(column, edgecoeffs[i, 0])
        top = I.poly(column, edgecoeffs[i, 1])
        peaks = np.full((numstep, naps), np.nan)
        for j in range(numstep):
            if bot[j] < -0.5 or top[j] > nrows - 0.5:
                continue
            yb = max(0.0, bot[j])
            yt = min(nrows - 1.0, top[j])
            c = int(column[j])
            r0, r1 = I.iround(yb), I.iround(yt)
            zimg = I.median(image[r0:r1 + 1, max(c - half, 0):min(c + half, ncols - 1) + 1],
                            even=True, axis=1).astype(np.float32)
            yarc = spatcal[r0:r1 + 1, c]
            yidx = y[r0:r1 + 1]
            z = np.isfinite(yarc)
            zimg, yarc, yidx = zimg[z], yarc[z], yidx[z]
            for k in range(naps):
                if yarc.size - 3 <= 0:
                    continue
                guessz = I.linterp(yarc, zimg, appos[i, k])
                guessy = I.linterp(yarc, yidx, appos[i, k])
                if guessy >= nrows - 1:
                    continue
                _, fit = I.gaussfit(yidx, zimg, nterms=4,
                                    estimates=np.array([guessz, guessy, 1, 0], dtype=np.float32))
                if fit[2] > sigthresh and abs(fit[1] - guessy) <= winthresh:
                    peaks[j, k] = fit[1]
        for k in range(naps):
            wave = I.interpolate2d(wavecal, column, peaks[:, k], cubic=-1.0)
            spat = I.interpolate2d(spatcal, column, peaks[:, k], cubic=-1.0)
            c, gb = I.robustpoly1d(wave, spat, fitorder, 3, 0.01)
            coeffs[i * naps + k] = c
            points.append(dict(order_index=i, aperture=k, column=column,
                               row=peaks[:, k], goodbad=gb))
    return coeffs, points


def tracetoap2d(omask, wavecal, spatcal, tracecoeffs, naps, orders, doorders=None):
    """
    MC_TRACETOAP2D: an aperture's path back on the detector, for display,
    and its mean slit position (MEAN of the DOUBLE trace, kept in a FLOAT
    array), which Import Trace takes as the aperture position.
    """
    nrows, ncols = omask.shape
    yy, xx = np.mgrid[0:nrows, 0:ncols]
    doorders = [1] * len(orders) if doorders is None else doorders
    out, l = [], 0
    for i, o in enumerate(orders):
        if not doorders[i]:
            continue
        z = omask == o
        wmin, wmax = np.nanmin(wavecal[z]), np.nanmax(wavecal[z])
        xmin, xmax = xx[z].min(), xx[z].max()
        dw = (wmax - wmin) / (xmax - xmin)
        wave = np.arange(xmax - xmin - 1) * dw + wmin
        for j in range(naps):
            tr = I.poly(wave, tracecoeffs[l])
            (ix, iy), _ = I.trigrid_points(wavecal[z], spatcal[z],
                                           [xx[z].astype(float), yy[z].astype(float)], wave, tr)
            out.append(dict(order=o, aperture=j + 1, x=ix, y=iy, appos=np.float32(I.itotal(tr) / tr.size)))
            l += 1
    return out


def mkapmask1d(slit, positions, apradius, psbginfo=None, xsbginfo=None):
    """
    MC_MKAPMASK1D: aperture k is marked k+1 (fractional at its ends),
    background -1, the rest 0. Raises if apertures overlap, as IDL refuses.

    The positions (FLTARR trace), FLOAT(apradius) and the background limits
    are single precision; TABINV on the double slit grid is double. All the
    limits go through one TABINV call (it works element by element), and the
    mask is then written in the IDL's order.
    """
    slit = np.asarray(slit)
    slit = slit.astype(np.float32 if I.isf32(slit) else np.float64)   # 1D slits are FLOAT
    npix = slit.size
    naps = len(positions)
    f32 = np.float32
    positions = np.asarray(positions).astype(f32)
    apradius = np.asarray(apradius).astype(f32)
    pts = [positions - apradius, positions + apradius]                 # apertures
    if psbginfo is not None:
        bgstart, bgwidth = f32(psbginfo[0]), f32(psbginfo[1])
        pts += [positions - bgstart - bgwidth, positions - bgstart,     # bottom bg
                positions + bgstart, positions + bgstart + bgwidth]     # top bg
    xs = [] if xsbginfo is None else [np.array([r[0], r[1]], f32) for r in xsbginfo]
    flat = np.concatenate([np.asarray(p_, f32).ravel() for p_ in pts] + xs) if (pts or xs) else np.zeros(0, f32)
    t = np.atleast_1d(I.tabinv(slit, flat.astype(f32)))
    lo, hi = t[0:naps], t[naps:2 * naps]
    tr = I.itrunc
    mask = np.zeros(npix, dtype=np.float32)
    for i in range(naps):
        mask[tr(lo[i]):tr(hi[i]) + 1] += 1
    if (mask > 1).any():
        raise ValueError("The extraction apertures overlap. Please lower the aperture radii.")
    mask[:] = 0
    k = 2 * naps
    if psbginfo is not None:
        b0, b1 = t[k:k + naps], t[k + naps:k + 2 * naps]
        t0, t1 = t[k + 2 * naps:k + 3 * naps], t[k + 3 * naps:k + 4 * naps]
        for i in range(naps):
            mask[tr(b0[i]):tr(b1[i]) + 1] = -1
            mask[tr(t0[i]):tr(t1[i]) + 1] = -1
        for i in range(naps):
            # the cleared zone [pos - bgstart, pos + bgstart]: the same limits
            mask[tr(b1[i]):tr(t0[i]) + 1] = 0
        k += 4 * naps
    for j in range(len(xs)):
        mask[tr(t[k + 2 * j]):tr(t[k + 2 * j + 1]) + 1] = -1
    for i in range(naps):
        a = (lo[i], hi[i])
        a0, a1 = tr(a[0]), tr(a[1])
        mask[a0:a1 + 1] = i + 1
        if a[0] - np.floor(a[0]) >= 0.5:
            mask[a0] = 0
            mask[tr(a[0] + 1)] = (0.5 + I.iround(a[0]) - a[0]) + i
        else:
            mask[a0] = (0.5 - (a[0] - np.floor(a[0]))) + i
        if a[1] - np.floor(a[1]) >= 0.5:
            mask[tr(a[1] + 1)] = (0.5 - (I.iround(a[1]) - a[1])) + i
        else:
            mask[a1] = (0.5 + (a[1] - I.iround(a[1]))) + i
    return mask


# ------------------------------------------------------------- exact-area cells

@njit(cache=True)
def _clip_area(px, py, i, j):
    """
    Area of the polygon (px, py) inside pixel [i, i+1] x [j, j+1], as the
    compiled clipper POLYFILLAA calls (Spextool's other/misc/polyclip.c)
    computes it: vertices in single precision, Sutherland-Hodgman against
    the left, right, top and bottom edges in turn, each intersection
    S + (P - S) slope * (edge - S) in float, the shoelace sum in double,
    returned in float.
    """
    f32 = np.float32
    xs = np.empty(24, dtype=np.float32)
    ys = np.empty(24, dtype=np.float32)
    n = px.shape[0]
    for k in range(n):
        xs[k] = px[k]
        ys[k] = py[k]
    fi = f32(i)
    fi1 = f32(i + 1)
    fj = f32(j)
    fj1 = f32(j + 1)
    for side in range(4):
        ox = np.empty(24, dtype=np.float32)
        oy = np.empty(24, dtype=np.float32)
        m = 0
        sx_ = f32(0.)
        sy_ = f32(0.)
        last_in = False
        for k in range(n + 1):
            closing = k == n
            if closing:
                if n == 0:
                    break
                cx, cy = xs[0], ys[0]
            else:
                cx, cy = xs[k], ys[k]
            if side == 0:
                cin = cx >= fi
            elif side == 1:
                cin = cx <= fi1
            elif side == 2:
                cin = cy <= fj1
            else:
                cin = cy >= fj
            if k > 0 and (last_in != cin):
                # every operation rounded to float, as the C code does
                # (numba would otherwise divide float32s in double)
                if side < 2:
                    edge = fi if side == 0 else fi1
                    q = f32(f32(cy - sy_) / f32(cx - sx_))
                    ox[m] = edge
                    oy[m] = f32(sy_ + f32(q * f32(edge - sx_)))
                else:
                    edge = fj1 if side == 2 else fj
                    q = f32(f32(cx - sx_) / f32(cy - sy_))
                    ox[m] = f32(sx_ + f32(q * f32(edge - sy_)))
                    oy[m] = edge
                m += 1
            if closing:
                break
            if cin:
                ox[m] = cx
                oy[m] = cy
                m += 1
            sx_ = cx
            sy_ = cy
            last_in = cin
        n = m
        if n == 0:
            return f32(0.)
        for k in range(n):
            xs[k] = ox[k]
            ys[k] = oy[k]
    a = 0.0
    for k in range(n):
        kk = 0 if k == n - 1 else k + 1
        # np.float64: numba's float() of a float32 multiplies in single precision
        a += np.float64(xs[k]) * np.float64(ys[kk]) - np.float64(ys[k]) * np.float64(xs[kk])
    if a < 0.:
        a = -a
    return f32(a / 2.0)


@njit(cache=True)
def _cell_sums(px, py, img, var):
    """
    POLYFILLAA over one column of cells, summed as MC_EXTSPEC2D sums:
    total(img[pixels]*area) in single precision, columns outermost.
    """
    ncell = px.shape[0]
    sy, sx = img.shape
    flx = np.empty(ncell, dtype=np.float32)
    vr = np.empty(ncell, dtype=np.float32)
    for c in range(ncell):
        qx = px[c]
        qy = py[c]
        ok = True
        for k in range(4):
            if not (np.isfinite(qx[k]) and np.isfinite(qy[k])):
                ok = False
        if not ok:
            flx[c] = np.nan
            vr[c] = np.nan
            continue
        # The bounding box from the double vertices, the clipping on their
        # FLOAT copies ("No double please"), as POLYFILLAA does.
        i0 = max(int(np.floor(qx.min())), 0)
        i1 = min(int(np.floor(qx.max())), sx - 1)
        j0 = max(int(np.floor(qy.min())), 0)
        j1 = min(int(np.floor(qy.max())), sy - 1)
        qx32 = qx.astype(np.float32)
        qy32 = qy.astype(np.float32)
        f = np.float32(0.0)
        v = np.float32(0.0)
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                a = _clip_area(qx32, qy32, i, j)
                if a != 0.0:
                    f = np.float32(f + img[j, i] * a)
                    v = np.float32(v + var[j, i] * a)
        flx[c] = f
        vr[c] = v
    return flx, vr


def _idlmax(a, b):
    """IDL's a > b (maximum): b where a is NaN."""
    return np.where(a > b, a, b)


def _idlmin(a, b):
    """IDL's a < b (minimum): a where a is NaN."""
    return np.where(a > b, b, a)


def extspec2d(img, var, omask, orders, wavecal, spatcal, ds, apsign, tracecoeffs,
              apradii, psbginfo=None, xsbginfo=None, bgdeg=None, psfradius=None,
              smapinfo=None):
    """
    MC_EXTSPEC2D: sum or optimal extraction.

    Each order is cut into cells of one resampled wavelength step by `ds`
    arcsec; every cell's flux is the exact-area sum of the pixels under it.
    The background fit, when asked for, is a robust polynomial of degree
    `bgdeg` across the slit at each wavelength, its uncertainty added to
    the variance. Returns {(order, aperture): [wave, flux, error]}.

    Without `psfradius` the flux is the sum over the aperture (fractional
    weights at its ends). With `psfradius` and `smapinfo` (one (ygrid,
    coeffs) pair per extracted order, from mkspatprof2d's SPATCOEFFS or
    MEDSPATCOEFFS) it is the optimal estimate xspextool makes for TS4 point
    sources: the profile is rebuilt at each wavelength from the spatial
    map, sinc-interpolated onto the slit, scaled to the slit with a robust
    linear fit whose 6-sigma rejects are left out, forced to the aperture's
    sign, normalized within the PSF radius, and the flux is the
    variance-weighted mean of data/profile over the aperture.
    """
    optext = bool(psfradius)                         # KEYWORD_SET(PSFRADIUS)
    if smapinfo is not None and not optext:
        # 'Fix Bad Pixels' (FIXBADPIXELS=0 0): insensitive for TS4
        raise ValueError("A spatial map without a PSF radius (Fix Bad Pixels) is not "
                         "available for TS4_SOAR.")
    if optext and (smapinfo is None or bgdeg is None):
        # xspextool switches optimal extraction on only with PS BG Subtraction,
        # and always hands it the spatial map.
        raise ValueError("Optimal extraction needs the spatial map and the background fit.")
    nrows, ncols = img.shape
    naps = apradii.shape[1]
    f32 = np.float32
    out = {}
    for i, o in enumerate(orders):
        xgrid, ygrid, ix, iy, dw = _order_grid(omask, wavecal, spatcal, o, ds)
        # [lambda, slit] in IDL; numpy [slit, lambda] here, so transpose.
        ix, iy = ix.T, iy.T
        nx, ny = xgrid.size, ygrid.size
        yb = 3
        ix, iy = ix[:, yb:ny - yb - 1 + 1], iy[:, yb:ny - yb - 1 + 1]
        ygrid = ygrid[yb:ny - yb]
        with np.errstate(invalid="ignore"):
            ixnan = np.sum(ix, axis=1)
            zx = np.where(np.isfinite(ixnan))[0]
            xgood = [0, ixnan.size - 1] if zx.size == 0 else [zx[0], zx[-1]]
            iynan = np.sum(ix, axis=0)
            zy = np.where(np.isfinite(iynan))[0]
            ygood = [0, iynan.size - 1] if zy.size == 0 else [zy[0], zy[-1]]
        xgrid = xgrid[xgood[0]:xgood[1] + 1]
        ygrid = ygrid[ygood[0]:ygood[1] + 1]
        ix = ix[xgood[0]:xgood[1] + 1, ygood[0]:ygood[1] + 1]
        iy = iy[xgood[0]:xgood[1] + 1, ygood[0]:ygood[1] + 1]
        xmask = np.sum((ix <= 2).astype(int) + (ix >= ncols - 2).astype(int), axis=1)
        zx = np.where(xmask == 0)[0]
        ymask = np.sum((iy <= 2).astype(int) + (iy >= nrows - 2).astype(int), axis=0)
        zy = np.where(ymask == 0)[0]
        xgrid = xgrid[zx[0]:zx[-1] + 1]
        ygrid = ygrid[zy[0]:zy[-1] + 1]
        ix = ix[zx[0]:zx[-1] + 1, zy[0]:zy[-1] + 1]
        iy = iy[zx[0]:zx[-1] + 1, zy[0]:zy[-1] + 1]
        nw = xgrid.size - 1
        ns = ygrid.size - 1
        slit_arc = (ygrid + ds / 2.)[0:ns]
        owave = (xgrid + dw / 2.)[0:nw]
        ofspec = np.full((nw, naps), np.nan, dtype=np.float32)
        ovspec = np.full((nw, naps), np.nan, dtype=np.float32)
        trace_arc = np.array([[I.poly(owave[j], tracecoeffs[i * naps + k]) for k in range(naps)]
                              for j in range(nw)], dtype=np.float32)       # FLTARR
        imgf = np.ascontiguousarray(img, dtype=np.float32)
        varf = np.ascontiguousarray(var, dtype=np.float32)
        if optext:
            # The spatial map: every slit row's polynomial evaluated on the
            # output wavelengths (DOUBLE POLY stored in an FLTARR).
            yarc, cmap = smapinfo[i]
            cmap = np.asarray(cmap)
            spatmap = np.zeros((nw, cmap.shape[0]), dtype=f32)
            for r in range(cmap.shape[0]):
                spatmap[:, r] = I.poly(owave, cmap[r])
            sincw = _sinc_weights(yarc, slit_arc)
            psfr = np.full(naps, psfradius, dtype=f32)              # REPLICATE(psfradius,naps)
        for j in range(nw):
            px = np.stack([ix[j, 0:ns], ix[j, 1:ns + 1], ix[j + 1, 1:ns + 1], ix[j + 1, 0:ns]], axis=1) + 0.5
            py = np.stack([iy[j, 0:ns], iy[j, 1:ns + 1], iy[j + 1, 1:ns + 1], iy[j + 1, 0:ns]], axis=1) + 0.5
            slit_flx, slit_var = _cell_sums(px.astype(np.float64), py.astype(np.float64), imgf, varf)
            slitmask = mkapmask1d(slit_arc, trace_arc[j], apradii[i], psbginfo=psbginfo,
                                  xsbginfo=xsbginfo)
            if optext:
                psfmask = mkapmask1d(slit_arc, trace_arc[j], psfr)
            if bgdeg is not None:
                z = np.where((slitmask == -1.0) & (slit_var[:ns] != 0))[0]
                if z.size >= 2:
                    c, _, cvar = I.robustpoly1d(slit_arc[z], slit_flx[z], bgdeg, 4, 0.1,
                                                yerr=np.sqrt(slit_var[z]), want_covar=True)
                    bg, yvar = I.poly1d(slit_arc, c, cvar)
                    # FLOAT minus the DOUBLE fit: the slit stays double, as in IDL
                    slit_flx = slit_flx[:ns] - bg
                    slit_var = slit_var[:ns] + yvar
            if optext:
                # Without a background fit at this wavelength the IDL slit
                # stays ny long and ROBUSTPOLY1D pairs it with an ns-long
                # profile; the ns values that exist are used here.
                sflx, svar = slit_flx[:ns], slit_var[:ns]
                sprofile = _sinc_apply(sincw, spatmap[j])
                _, ogoodbad = I.robustpoly1d(sprofile, sflx, 1, 6, 0.1)
                for k in range(naps):
                    zpsf = (psfmask > f32(k)) & (psfmask <= f32(k + 1))
                    ps_ = sprofile[zpsf]
                    # enforce the aperture's sign (IDL's > and < operators)
                    ps_ = _idlmax(ps_, 0.0) if apsign[k] == 1 else _idlmin(ps_, 0.0)
                    ps_ = apsign[k] * np.abs(ps_ / I.itotal(ps_[np.isfinite(ps_)]))
                    sprofile[zpsf] = ps_
                    zs = ((slitmask > f32(k)) & (slitmask <= f32(k + 1)) &
                          (sprofile != 0.0) & (ogoodbad == 1))
                    if zs.any():
                        with np.errstate(divide="ignore", invalid="ignore"):
                            flux = sflx[zs] / sprofile[zs]
                            fvar = svar[zs] / I.sq(sprofile[zs])
                            m, v, _ = I.meancomb(flux, datavar=fvar)
                        ofspec[j, k], ovspec[j, k] = m, v
                    else:
                        ofspec[j, k] = ovspec[j, k] = np.nan
                continue
            for k in range(naps):
                z = (slitmask > np.float32(k)) & (slitmask <= np.float32(k + 1))
                w = (slitmask[z] - np.float32(k)).astype(np.float32)
                tot = I.ftotal if I.isf32(slit_flx) else I.itotal
                ofspec[j, k] = apsign[k] * tot(slit_flx[:ns][z] * w) if z.any() else 0.0
                tot = I.ftotal if I.isf32(slit_var) else I.itotal
                ovspec[j, k] = tot(slit_var[:ns][z] * I.sq(w)) if z.any() else 0.0
        for k in range(naps):
            out[(o, k + 1)] = np.array([owave, ofspec[:, k], np.sqrt(ovspec[:, k])])
    return out


# ------------------------------------------------------------- 1D (no wavecal file)

@njit(cache=True)
def _sinc_kernel(fx, dampfac, pi32):
    """MC_SINCINTERP's damped sinc for each fractional offset (DOUBLE, libm exp/sin)."""
    n = fx.shape[0]
    w = np.zeros((n, 21))
    for p in range(n):
        if fx[p] == 0.0:
            continue
        for k in range(21):
            xk = np.float64(np.float32(k) - np.float32(10)) - fx[p]
            u1 = xk / dampfac
            u2 = pi32 * xk
            w[p, k] = math.exp(-(u1 * u1)) * math.sin(u2) / u2
    return w


@njit(cache=True)
def _sinc_kernel32(fx, dampfac, pi32):
    """The same in FLOAT (both abscissae FLOAT): expf, sinf, every step rounded."""
    n = fx.shape[0]
    w = np.zeros((n, 21), dtype=np.float32)
    for p in range(n):
        if fx[p] == np.float32(0):
            continue
        for k in range(21):
            xk = np.float32(np.float32(np.float32(k) - np.float32(10)) - fx[p])
            u1 = np.float32(xk / dampfac)
            u2 = np.float32(pi32 * xk)
            e = np.float32(math.exp(np.float32(-np.float32(u1 * u1))))
            sn = np.float32(math.sin(u2))
            w[p, k] = np.float32(np.float32(e * sn) / u2)
    return w


def _sinc_weights_fast(xa, x, dampfac=3.25, ksize=21):
    """
    _sinc_weights with the kernel loop compiled: the 1D extraction needs one
    per column. There the table and the slit are both FLOAT, so TABINV and
    the whole kernel are single precision.
    """
    single = I.isf32(np.asarray(xa)) and I.isf32(np.asarray(x))
    idx = np.atleast_1d(I.tabinv(xa, x))
    ix = np.trunc(idx).astype(np.int64)
    if single:
        fx = (idx.astype(np.float32) - ix.astype(np.float32)).astype(np.float32)
        w = _sinc_kernel32(fx, np.float32(dampfac), np.float32(np.pi))
    else:
        fx = idx.astype(np.float64) - ix
        w = _sinc_kernel(fx, float(dampfac), float(np.float32(np.pi)))
    nf = np.asarray(xa).size
    lobe = np.clip(np.arange(ksize)[None, :] - 10 + ix[:, None], 0, nf - 1)
    return ix, fx, w, lobe


def _moments_goodbad(data, igoodbad, robust):
    """
    MC_MOMENTS' OGOODBAD with ROBUST (the moments themselves unused here).
    As written, the flags of the points kept are set in a count_initgood
    array at their positions in the full data (GOOD_BAD[Z[GOODDATAS]]),
    out-of-range subscripts clipped to the last element, and that array is
    then copied back onto the initially good points.
    """
    d = np.asarray(data)
    gb = np.array(igoodbad, dtype=np.int64).copy()
    gb[~np.isfinite(d)] = 2
    z = np.where(gb == 1)[0]
    if z.size <= 1:
        raise ValueError("mc_moments: fewer than two good points")
    dz = d[z]
    med = I.median(dz, even=True)
    P = np.float32 if I.isf32(d) else np.float64
    mad = P(np.float32(1.482) * I.median(np.abs(dz - med), even=True))
    with np.errstate(divide="ignore", invalid="ignore"):
        good = np.where(np.abs((dz - med) / mad) <= robust)[0]
    if good.size == 0:
        good = np.arange(dz.size)
    good_bad = np.zeros(z.size, dtype=np.int64)
    good_bad[np.minimum(z[good], z.size - 1)] = 1
    gb[z] = good_bad
    return gb


def _nantrim2(a):
    """MC_NANTRIM(array, 2): drop the leading and trailing NaNs (and zeros)."""
    a = np.asarray(a, dtype=float)
    c = np.cumsum(np.where(np.isnan(a), 0.0, a))
    if not np.isnan(a).any():
        return np.arange(a.size)
    lo = int(np.where(c != 0)[0][0])
    cr = np.cumsum(np.where(np.isnan(a[::-1]), 0.0, a[::-1]))
    hi = a.size - 1 - int(np.where(cr != 0)[0][0])
    return np.arange(lo, hi + 1)


def _sqrt_store(v):
    """SQRT(TOTAL(...)) stored in an FLTARR: a FLOAT root, or a DOUBLE root rounded."""
    return np.sqrt(v) if I.isf32(np.asarray(v)) else np.float32(np.sqrt(np.float64(v)))


def _spatmap1d(owave, spatcoeffs_i):
    yarc, c = spatcoeffs_i
    c = np.asarray(c)
    smap = np.zeros((owave.size, c.shape[0]), dtype=np.float32)
    for r in range(c.shape[0]):
        smap[:, r] = I.poly(owave, c[r])
    return yarc, smap


def extpsspec(image, var, wavecal, edgecoeffs, tracecoeffs, naps, xranges, slith_arc,
              apradii, apsign, psfwidth=None, spatcoeffs=None, bgorder=None, psbginfo=None,
              bdpxmask=None, bitmask=None, bdpixthresh=7.0):
    """
    MC_EXTPSSPEC: the SpeX-type point-source extraction xspextool uses when
    there is no wavecal file. Every detector column of an order is one
    output pixel (the 'wavelength' is the column number); the slit is the
    column's pixels between the order edges, in arcsec by linear scaling.
    With the background, a robust fit (5 sigma MC_MOMENTS clip, then
    ROBUSTPOLY1D) over the background pixels is subtracted. With
    `psfwidth` (and the spatial map) the flux is the optimal estimate, bad
    pixels found by a `bdpixthresh`-sigma fit of the profile; otherwise the
    sum. A 4th row carries the flags: 1 linearity, 8 optimal extraction
    failed. Returns a list, one {aperture: [4, n]} per order.
    """
    nrows, ncols = image.shape
    f32 = np.float32
    optextract = psfwidth is not None
    fixbdpx = spatcoeffs is not None and not optextract
    if fixbdpx:
        raise ValueError("Fix Bad Pixels (FIXBADPIXELS=0 0) is not available for TS4_SOAR.")
    if bdpxmask is None:
        bdpxmask = np.ones(image.shape, dtype=np.int16)
    linmask = (np.asarray(bitmask) & 1) if bitmask is not None else np.zeros(image.shape, np.uint8)
    y = np.arange(nrows, dtype=np.float32)
    apradii = np.asarray(apradii, dtype=f32)
    out = []
    for i in range(len(xranges)):
        start, stop = int(xranges[i][0]), int(xranges[i][1])
        x = np.arange(stop - start + 1, dtype=f32) + f32(start)
        bot = I.poly(x, edgecoeffs[i][0])
        top = I.poly(x, edgecoeffs[i][1])
        owave = wavecal[I.itrunc((bot + top) / 2.), I.itrunc(x)]
        n = owave.size
        ofspec = np.full((n, naps), np.nan, dtype=f32)
        oespec = np.full((n, naps), np.nan, dtype=f32)
        obspec = np.zeros((n, naps), dtype=f32)
        trace = np.zeros((n, naps), dtype=f32)
        for j in range(naps):
            trace[:, j] = I.poly(owave, tracecoeffs[i * naps + j])
        pixtoarc = np.zeros((n, 2), dtype=f32)
        pixtoarc[:, 1] = np.float64(f32(slith_arc)) / (top - bot)
        pixtoarc[:, 0] = -1. * (pixtoarc[:, 1] * bot)
        if optextract:
            yarc, spatmap = _spatmap1d(owave, spatcoeffs[i])
            psfr = np.full(naps, psfwidth, dtype=f32)
        for j in range(n):
            if top[j] > nrows - 0.5 or bot[j] < -0.5:
                continue
            b, t, c = I.itrunc(bot[j]), I.itrunc(top[j]), I.itrunc(x[j])
            slity_arc = I.poly(y[b:t + 1], pixtoarc[j])
            slitz = image[b:t + 1, c]
            vslitz = var[b:t + 1, c]
            bdpxmkz = np.asarray(bdpxmask[b:t + 1, c]).astype(np.int64)
            linmkz = linmask[b:t + 1, c]
            mask = mkapmask1d(slity_arc, trace[j], apradii, psbginfo=psbginfo)
            if optextract:
                psfmask = mkapmask1d(slity_arc, trace[j], psfr)
            if bgorder is not None:
                z = np.where((mask == -1) & (vslitz != 0.0))[0]
                if z.size >= 2:
                    og = _moments_goodbad(slitz[z], bdpxmkz[z], 5)
                    coeff, _, cvar = I.robustpoly1d(slity_arc[z], slitz[z], bgorder, 4, 0.1,
                                                    yerr=np.sqrt(vslitz[z]), igoodbad=og,
                                                    want_covar=True)
                    bg, yvar = I.poly1d(slity_arc, coeff, cvar)
                    slitz = slitz - bg
                    vslitz = vslitz + yvar
            if optextract:
                sprofile = _sinc_apply(_sinc_weights_fast(yarc, slity_arc), spatmap[j])
                bdpx = bdpxmkz.copy()
                _, ogoodbad = I.robustpoly1d(sprofile, slitz, 1, bdpixthresh, 0.1, igoodbad=bdpx)
                # ROBUSTPOLY1D marks the NaNs 2 in the IGOODBAD array it is given
                bdpx[~(np.isfinite(sprofile) & np.isfinite(slitz))] = 2
            for k in range(naps):
                if optextract:
                    zpsf = (psfmask > f32(k)) & (psfmask <= f32(k + 1))
                    ps_ = sprofile[zpsf]
                    ps_ = _idlmax(ps_, np.float32(0.0)) if apsign[k] == 1 else _idlmin(ps_, np.float32(0.0))
                    tot = I.ftotal if I.isf32(ps_) else I.itotal
                    ps_ = (apsign[k] * np.abs(ps_ / tot(ps_[np.isfinite(ps_)]))).astype(ps_.dtype)
                    sprofile[zpsf] = ps_
                    # IGOODBAD=bdpxmkz is the caller's array: its NaNs became 2 in it
                    goodbad = ogoodbad * bdpx
                    zs = ((mask > f32(k)) & (mask <= f32(k + 1)) & (sprofile != 0.0) &
                          (goodbad == 1))
                    with np.errstate(divide="ignore", invalid="ignore"):
                        if zs.any():
                            be = slitz[zs] / sprofile[zs]
                            bv = vslitz[zs] / I.sq(sprofile[zs])
                            if (linmkz[zs] == 1).any():
                                obspec[j, k] += 1
                        else:
                            be = slitz * 0.0
                            bv = slitz / slitz
                            obspec[j, k] += 8
                        m, v, _ = I.meancomb(be, datavar=bv)
                    ofspec[j, k] = m
                    oespec[j, k] = np.sqrt(f32(v))
                else:
                    z = (mask > f32(k)) & (mask <= f32(k + 1))
                    w = (mask[z] - f32(k)).astype(f32)
                    tot = I.ftotal if I.isf32(slitz) else I.itotal
                    ofspec[j, k] = f32(apsign[k]) * tot(slitz[z] * w)
                    tot = I.ftotal if I.isf32(vslitz) else I.itotal
                    oespec[j, k] = _sqrt_store(tot(vslitz[z] * I.sq(w)))
                    if (linmkz[z] == 1).any():
                        obspec[j, k] += 1
        keep = _nantrim2(owave)
        out.append({k + 1: np.array([owave[keep], ofspec[keep, k], oespec[keep, k],
                                     obspec[keep, k]], dtype=f32) for k in range(naps)})
    return out


def extxsspec(image, var, wavecal, edgecoeffs, tracecoeffs, naps, xranges, slith_arc,
              apradii, bgr=None, bgorder=None, spatcoeffs=None, bdpxmask=None, bitmask=None,
              bdpixthresh=7.0, errorprop=False):
    """
    MC_EXTXSSPEC: the SpeX-type extended-source extraction (no wavecal
    file). The background regions `bgr` are fitted and subtracted (in 1D
    the background does work); XSPEXTOOL_XSEXTRACTSPEC does not pass
    ERRORPROP, so after a background fit the variance is set to 1. With
    the spatial map (xspextool hands it over whenever the point-source
    Optimal flag is on, the TS4 default) pixels off a `bdpixthresh`-sigma
    fit of the profile are replaced by the scaled profile and their
    variance set to NaN. Flags: 1 linearity, 2 a bad pixel in the
    aperture. Returns a list, one {aperture: [4, n]} per order.
    """
    nrows, ncols = image.shape
    f32 = np.float32
    fixbdpx = spatcoeffs is not None
    if bdpxmask is None:
        bdpxmask = np.ones(image.shape, dtype=np.int16)
    linmask = (np.asarray(bitmask) & 1) if bitmask is not None else np.zeros(image.shape, np.uint8)
    y = np.arange(nrows, dtype=np.float32)
    apradii = np.asarray(apradii, dtype=f32)
    out = []
    for i in range(len(xranges)):
        start, stop = int(xranges[i][0]), int(xranges[i][1])
        x = np.arange(stop - start + 1, dtype=f32) + f32(start)
        bot = I.poly(x, edgecoeffs[i][0])
        top = I.poly(x, edgecoeffs[i][1])
        owave = wavecal[I.itrunc((bot + top) / 2.), I.itrunc(x)]
        n = owave.size
        ofspec = np.full((n, naps), np.nan, dtype=f32)
        oespec = np.full((n, naps), np.nan, dtype=f32)
        obspec = np.zeros((n, naps), dtype=f32)
        trace = np.zeros((n, naps), dtype=f32)
        for j in range(naps):
            trace[:, j] = I.poly(owave, tracecoeffs[i * naps + j])
        pixtoarc = np.zeros((n, 2), dtype=f32)
        pixtoarc[:, 1] = np.float64(f32(slith_arc)) / (top - bot)
        pixtoarc[:, 0] = -1. * (pixtoarc[:, 1] * bot)
        if fixbdpx:
            yarc, spatmap = _spatmap1d(owave, spatcoeffs[i])
        for j in range(n):
            if top[j] > nrows - 0.5 or bot[j] < -0.5:
                continue
            lo = max(0.0, float(bot[j]))
            hi = min(float(nrows - 1), float(top[j]))
            b, t, c = I.itrunc(lo), I.itrunc(hi), I.itrunc(x[j])
            slity_arc = I.poly(y[b:t + 1], pixtoarc[j])
            slitz = image[b:t + 1, c]
            vslitz = var[b:t + 1, c].copy()
            bb, tt = I.itrunc(bot[j]), I.itrunc(top[j])
            bdpxmkz = np.asarray(bdpxmask[bb:tt + 1, c]).astype(np.int64)
            linmkz = linmask[bb:tt + 1, c]
            mask = mkapmask1d(slity_arc, trace[j], apradii, xsbginfo=bgr)
            if bgorder is not None:
                z = np.where((mask == -1) & (vslitz != 0.0))[0]
                og = _moments_goodbad(slitz[z], bdpxmkz[z], 5)
                coeff, _, cvar = I.robustpoly1d(slity_arc[z], slitz[z], bgorder, 4, 0.1,
                                                igoodbad=og, yerr=np.sqrt(vslitz[z]),
                                                want_covar=True)
                if errorprop:
                    bg, yvar = I.poly1d(slity_arc, coeff, cvar)
                    slitz = slitz - bg
                    vslitz = vslitz + yvar
                else:
                    bg, _ = I.poly1d(slity_arc, coeff)
                    slitz = slitz - bg
                    vslitz[:] = 1.0
            if fixbdpx:
                sprofile = I.linterp(yarc, spatmap[j], slity_arc)
                try:
                    cf, ogoodbad = I.robustpoly1d(sprofile, slitz, 1, bdpixthresh, 0.1)
                except (ValueError, np.linalg.LinAlgError):
                    continue                          # CANCEL: GOTO, cont2
                scs = I.poly(sprofile, cf)
                bad = (ogoodbad == 0) | (bdpxmkz == 0)
                if bad.any():
                    slitz = np.array(slitz)
                    slitz[bad] = scs[bad]
                    vslitz[bad] = np.nan
            for k in range(naps):
                z = (mask > f32(k)) & (mask <= f32(k + 1))
                w = (mask[z] - f32(k)).astype(f32)
                tot = I.ftotal if I.isf32(slitz) else I.itotal
                ofspec[j, k] = tot(slitz[z] * w)
                tot = I.ftotal if I.isf32(vslitz) else I.itotal
                oespec[j, k] = _sqrt_store(tot(vslitz[z] * I.sq(w)))
                if (linmkz[z] == 1).any():
                    obspec[j, k] += 1
                if fixbdpx:
                    goodbad = ogoodbad * bdpxmkz
                    if ((goodbad == 0) & z).any():
                        obspec[j, k] += 2
        keep = _nantrim2(owave)
        out.append({k + 1: np.array([owave[keep], ofspec[keep, k], oespec[keep, k],
                                     obspec[keep, k]], dtype=f32) for k in range(naps)})
    return out


# ------------------------------------------------------------- writing

_WRITESPEC_HISTORY = (
    "The output FITS files contain a 3-D array of data consisting of sets of triplet arrays "
    "of data for each aperture and each order, where each triplet is composed of an array "
    "for the wavelength, an array for the flux, and an array for the error. The triplets for "
    "each aperture are stacked behind one another, followed by the triplets for the next "
    "order, etc. If no orders have been skipped or deslected in the extraction process, the "
    "contents of aperture Y in order X can be found as follows:lambda = array[*,0,( X - "
    "(smallest order number))*naps + (Y-1)], flux   = array[*,1,( X - (smallest order "
    "number))*naps + (Y-1)], error  = array[*,2,( X - (smallest order number))*naps + "
    "(Y-1)]  For example, for an SXD file with two apertures, the wavelength array for "
    "aperture 1 in order 3 is located in array [*,0,0], the flux is in array [*,1,0] and the "
    "error is in array [*,2,0]. For aperture 2, the wavelength is in array [*,0,1], the flux "
    "is in array [*,1,1] and the error is in array [*,2,1]. For order 4, the flux for "
    "aperture 1 is in array [*,1,2], while the flux for aperture 2 is in array [*,1,3]. For "
    "order 5, the fluxes for the two apertures are in arrays [*,1,4] and [*,1,5], etc.")


def write_spectrum(path, spectra, orders, naps, hdr0, info):
    """
    MC_WRITESPEC: a [norders*naps, 3, nwave] cube (IDL [nwave, 3, n]) of
    wavelength, flux and error, NaN-padded to the longest order, with the
    keywords every later Spextool step reads. The 1D extraction (no wavecal
    file, info['pixels']) adds a 4th row of flags and is in pixels, with no
    WAVECAL, WAVETYPE or dispersions.
    """
    keys = [(o, a) for o in orders for a in range(1, naps + 1)]
    npix = max(spectra[k].shape[1] for k in keys)
    nrow = max(spectra[k].shape[0] for k in keys)
    pixels = bool(info.get("pixels"))
    cube = np.full((len(keys), nrow, npix), np.nan)
    for n, k in enumerate(keys):
        s = spectra[k]
        cube[n, :, :s.shape[1]] = s
    hdr = fits.Header()
    f32 = np.float32
    # The instrument file's keywords, as MC_GETHDRINFO read them (FXPAR), then
    # MC_WRITESPEC's own; every number written in its IDL type by FXADDPAR.
    from triplespec_reducer.instrument import read_instrfile
    for key in read_instrfile().keywords:
        if key in hdr0:
            I.fxaddpar(hdr, key, I.fxpar(hdr0, key), I.fxcomment(hdr0, key))
    I.fxaddpar(hdr, "FILENAME", os.path.basename(path),
               I.fxcomment(hdr0, "FILENAME") if "FILENAME" in hdr0 else "")
    I.fxaddpar(hdr, "INSTR", "TS4_SOAR", " Instrument")
    I.fxaddpar(hdr, "MODENAME", info.get("modename", "XD"), " Spectroscopy mode")
    I.fxaddpar(hdr, "AIMAGE", info.get("aimage", ""), " A image")
    I.fxaddpar(hdr, "ITOT", f32(info.get("itot", 0.0)), " Total integration time (sec)")   # FLTARR
    I.fxaddpar(hdr, "LINCRMAX", 80000, " Maximum of linearity correction (DN)")   # LONG, instrument file
    I.fxaddpar(hdr, "SKY", info.get("sky", "None"), " Sky image")
    I.fxaddpar(hdr, "FLAT", info.get("flat", ""), " Flat field image")
    if not pixels:
        I.fxaddpar(hdr, "WAVECAL", info.get("wavecal", ""), " Wavecal file used for wavelength calibration")
        I.fxaddpar(hdr, "WAVETYPE", "Vacuum", " Wavelength type (air or vacuum)")
    if info.get("tracefile"):
        # XSPEXTOOL_WRITESPEC names the file only in an 'Import Coeffs' mode
        # the TS4 point-source panel does not have, so IDL never writes it.
        I.fxaddpar(hdr, "TRACEFIL", info["tracefile"], " Trace file name")
    I.fxaddpar(hdr, "NAPS", int(naps), " Number of apertures")
    I.fxaddpar(hdr, "NORDERS", len(orders), " Number of orders")
    I.fxaddpar(hdr, "ORDERS", ",".join(str(o) for o in orders), " Order numbers")
    I.fxaddpar(hdr, "SLTH_PIX", info["slith_pix"], " Slit height in pixels", fmt="G0.5")
    I.fxaddpar(hdr, "SLTH_ARC", info["slith_arc"], " Slit height in arcseconds", fmt="G0.5")
    I.fxaddpar(hdr, "SLTW_PIX", info["slitw_pix"], " Slit width in pixels", fmt="G0.5")
    I.fxaddpar(hdr, "SLTW_ARC", info["slitw_arc"], " Slit width in arcseconds", fmt="G0.5")
    I.fxaddpar(hdr, "RP", int(info["rp"]), " Resolving power")
    for i, o in enumerate(orders):
        I.fxaddpar(hdr, f"APPOSO{o:02d}", ",".join(sigfig(p, 3) for p in info["appos"][i]),
                   f" Aperture positions (arcsec) for order {o:02d}")
    if info.get("psfradius") is not None:
        I.fxaddpar(hdr, "PSFRAD", f32(info["psfradius"]),
                   " Optimal extraction PSF radius in arcseconds", fmt="f4.2")
    for i, r in enumerate(info["apradius"]):
        I.fxaddpar(hdr, f"AP{i + 1:02d}RAD", f32(r), " Aperture radius in arcseconds", fmt="f4.2")
    if info.get("psbginfo") is not None:
        b = np.asarray(info["psbginfo"], dtype=f32)
        I.fxaddpar(hdr, "BGSTART", b[0], " Background start radius in arcseconds", fmt="f4.2")
        I.fxaddpar(hdr, "BGWIDTH", b[1], " Background width in arcseconds", fmt="G0.5")
        I.fxaddpar(hdr, "BGORDER", b[2], " Background polynomial fit degree", fmt="I1.1")
    if info.get("xsbginfo") is not None:
        I.fxaddpar(hdr, "BGR", info["xsbginfo"][0], " Background regions in arcseconds")
        I.fxaddpar(hdr, "BGORDER", int(info["xsbginfo"][1]), " Background polynomial fit degree")
    if info.get("xs_background"):             # the opt-in, not an IDL keyword
        I.fxaddpar(hdr, "XSBGSUB", True, " Extended-source background actually subtracted")
    I.fxaddpar(hdr, "XUNITS", "pixels" if pixels else "um", " Units of the X axis")
    I.fxaddpar(hdr, "YUNITS", "DN / s", " Units of the Y axis")
    I.fxaddpar(hdr, "XTITLE", "!7k!5 (pixels)" if pixels else "!7k!5 (!7l!5m)", " IDL X title")
    I.fxaddpar(hdr, "YTITLE", "f (!5DN s!u-1!N)", " IDL Y title")
    I.fxaddpar(hdr, "VERSION", info.get("version", "TSR-1.0"), " Spextool version")
    for i, o in enumerate(orders if not pixels else []):
        # SXADDPAR of the DOUBLE dispersions MC_RDWAVECAL2D read
        I.fxaddpar(hdr, f"DISPO{o:02d}", np.float64(info["disp"][i]),
                   f" Dispersion (um pixel-1) for order {o:02d}")
    I.addhistlabel(hdr, "Xspextool History")
    I.sxaddhist(hdr, _WRITESPEC_HISTORY)
    fits.PrimaryHDU(cube.astype(np.float64), header=hdr).writeto(path, overwrite=True)


def caaspecflags(spec):
    """
    MC_CAASPECFLAGS: a spectrum with three columns gets a fourth (flags) by
    being copied into an FLTARR, so it becomes single precision; one with
    four is returned as it is.
    """
    spec = np.asarray(spec)
    if spec.shape[1] == 4:
        return spec
    new = np.zeros((spec.shape[0], 4) + spec.shape[2:], dtype=np.float32)
    new[:, 0:3] = spec
    return new


def read_spectrum(path):
    """
    MC_READSPEC: the cube [norders*naps, ncol, nwave] and its header, the
    numbers typed by FXPAR from their text as the IDL reads them.
    """
    with fits.open(path) as hdul:
        d = hdul[0].data
        spc = d.astype(d.dtype.newbyteorder("="))     # READFITS keeps the stored type
        hdr = hdul[0].header
    orders = [int(t) for t in str(hdr["ORDERS"]).split(",")]
    get = lambda k, default: I.fxpar(hdr, k) if k in hdr else default
    airmass = get("AIRMASS", 0)
    if airmass == 0:
        airmass = get("TCS_AM", 0)
    vals = [I.fxpar(hdr, k) for k in hdr if k.startswith("DISP")]        # MC_FXPAR 'DISP*'
    return dict(spec=spc, header=hdr, orders=orders, naps=I.fxpar(hdr, "NAPS"),
                norders=I.fxpar(hdr, "NORDERS"), airmass=airmass, rp=get("RP", 0),
                slith_pix=get("SLTH_PIX", 0), slith_arc=get("SLTH_ARC", 0),
                slitw_pix=get("SLTW_PIX", 0), slitw_arc=get("SLTW_ARC", 0),
                itot=get("ITOT", 0),
                disp=np.array(vals, dtype=np.result_type(*vals)) if vals else None)


def sigfig(x, n):
    """MC_SIGFIG's string: `n` significant figures, rounded from IDL's (e) text."""
    t = f"{abs(float(x)):e}"
    mant, ex = t.split("e")
    osign = -1 if ex[0] == "-" else 1
    order = int(ex[1:3])
    nodec = mant[0] + mant[2:]
    move = min(n - 1, len(mant) - 1)
    num0 = nodec[:1 + move] + "." + nodec[1 + move:]
    num1 = str(int(np.floor(abs(float(num0)) + 0.5)))       # ROUND, half away from zero
    order_inc = int(len(num1) > n)
    num = num1[0] + "." + num1[1:n]
    order = order + osign * order_inc
    new = float(num) * 10.0 ** (osign * order)
    if osign * order + 1 > n:
        s = str(int(np.round(new)))
    else:
        d = max(n - osign * order - 1, 0)
        s = f"{new:.{d}f}" + ("." if d == 0 else "")        # F40.0 keeps the point
    return ("-" if x < 0 else "") + s


# ------------------------------------------------------------- one pair, start to finish

def psdefine_signs(profiles, positions, doorders=None):
    """
    XSPEXTOOL_PSDEFINEAP's signs: an aperture is positive where the profile
    there lies above the profile's median. The loop runs over the orders
    and each overwrites the last, so the signs are the last order's.
    """
    orders = list(profiles)
    doorders = [1] * len(orders) if doorders is None else doorders
    naps = positions.shape[1]
    sign = np.zeros(naps, dtype=int)
    for i, o in enumerate(orders):
        if not doorders[i]:
            continue
        x, p = profiles[o]
        med = I.median(p)
        for j in range(naps):
            sign[j] = 1 if p[I.itrunc(I.tabinv(x, positions[i, j]))] > med else -1
    return sign


def check_ps_inputs(slith_arc, apradius, psbgsub=True, psbg=(2.5, 5.0, 1), optimal=True,
                    psfradius=2.2, naps=2, positions=None, tracestep=None, sumap=None,
                    tracedeg=None):
    """
    The range checks xspextool makes on the point-source inputs (MC_CRANGE
    in XSPEXTOOL_PSAPFIND, _PSTRACESPEC, _PSDEFINEAP and _PSEXTRACTSPEC),
    with its messages. Raises ValueError on the first one that fails.
    """
    f32 = np.float32
    if not 1 <= naps <= 4:
        raise ValueError("Cannot have more than four apertures per order.")
    if positions is not None:
        p = np.asarray(positions, dtype=f32)
        if (p < 0).any() or (p > slith_arc).any():
            raise ValueError(f"Pos must be between 0 and {slith_arc:g}.")
    if tracestep is not None and not 1 <= tracestep <= 20:
        raise ValueError("Step Size must be between 1 and 20.")
    if sumap is not None and (sumap < 1 or sumap % 2 == 0):
        raise ValueError("Colums Add must be odd.")
    if tracedeg is not None and not 0 <= tracedeg <= 4:
        raise ValueError("Trace Polynomial Degree must be between 0 and 4.")
    if optimal and not 0 <= f32(psfradius) <= slith_arc:
        raise ValueError(f"PSF Radius must be between 0 and {slith_arc:g}.")
    rng = f32(psfradius) if optimal else slith_arc
    if f32(apradius) > rng:
        raise ValueError(f"Ap Radius must be <= {rng:g}.")
    if psbgsub:
        if not f32(psbg[0]) > f32(apradius):
            raise ValueError(f"Background Start must be > {f32(apradius):g}.")
        if not 0 <= f32(psbg[1]) <= slith_arc:
            raise ValueError(f"Background Width must be between 0 and {slith_arc:g}.")
        if not 0 <= psbg[2] <= 7:
            raise ValueError("Background Fit Order must be between 0 and 7.")


def reduce_pair(files, flatinfo, wcinfo, opaths, source="point", atran=None,
                npeaks=2, findmode="auto", guess=None, apradius=1.0,
                psbg=(2.5, 5.0, 1), psbgsub=True, optimal=None, psfradius=2.2,
                medprofile=False, ps_orders=None, xs_orders=None,
                xs_positions=None, xs_radii=None,
                xsbg=((2., 6.), (32., 36.)), xsbgdeg=0, xs_background=False,
                trace=True, mode="A-B", tracestep=5, sumap=5, winthresh=3.0,
                sigthresh=0.0, tracedeg=1, ybuffer=1, atmosthresh=0.4,
                flatname="", wavecalname="", progress=None, signs=None,
                combined=None, sky=None, trace_out=None, trace_in=None,
                flatfield=True, errorprop=True, xsbgsub=True, bdpixthresh=7.0):
    """
    What xspextool does for one pair of frames with "Do All Steps", for a
    point or an extended source. Writes opaths[0] (and opaths[1] for the
    negative apertures of an A-B point source) and returns everything
    computed on the way, for display.

    Point sources, with the TS4_SOAR defaults of the IDL panel:

    psbgsub      'PS BG Subtraction'. In the TS4 panel this box also sets
                 the (greyed) 'Optimal' box, so on it means optimal
                 extraction with the background fit, off a plain sum with
                 no background.
    optimal      None follows psbgsub, as the panel does. False with
                 psbgsub on gives the background-subtracted sum
                 MC_EXTSPEC2D makes without a PSF radius: NOT something the
                 TS4 panel can produce (its Optimal box is insensitive).
    psfradius    PSPSFRAD, arcsec; the aperture radius may not exceed it.
    medprofile   'Median Profile': the optimal profile is the mean spatial
                 profile at every wavelength (MEDSPATCOEFFS) instead of the
                 fitted spatial map (SPATCOEFFS).
    findmode     'auto' (npeaks peaks), 'guess' or 'fixed'. For the last two
                 `guess` is the one list of positions typed in the panel:
                 it is used for every order and its length is the number
                 of apertures. 'fixed' also ticks Use Ap Positions, so the
                 apertures are not traced.
    trace        False is 'Use Ap Positions': a constant trace at the
                 found positions.
    ps_orders    'PS Orders': the order numbers to extract (all if None);
                 xs_orders is the same for an extended source.
    trace_in     'Import Trace': the apertures, and the orders, come from
                 a stored trace file; the signs are then set from the
                 profile, as Do All Steps does.
    npeaks       PSNAPS (1-4), in every reduction mode.

    Extended sources: xsbgsub is 'XS BG Subtraction' (XSBGSUB=1). In 2D it
    only decides whether A/A-Sky files name the regions (the IDL fit is
    dead code; xs_background=True is the opt-in that applies it).

    With no wavecal file (wcinfo from wavecal.no_wavecal, wctype '1D') the
    extraction is MC_EXTPSSPEC / MC_EXTXSSPEC on detector columns, with
    the bad pixel mask, the linearity bitmask and `bdpixthresh`
    (BADPIXELTHRESH 7); there the XS background is subtracted when
    xsbgsub is on, and, as xspextool hands MC_EXTXSSPEC the spatial map
    whenever the point-source Optimal flag is on, bad pixels of an
    extended source are fixed by default.
    """
    say = progress or (lambda *_: None)
    fl, wc = flatinfo, wcinfo
    f32 = np.float32
    allorders = list(fl["orders"])
    no = len(allorders)
    if optimal is None:
        optimal = bool(psbgsub)
    if source == "point" and optimal and not psbgsub:
        raise ValueError("Optimal extraction needs PS BG Subtraction (the TS4 panel "
                         "switches them on and off together).")
    say("Loading the images")
    work, var, hdrs, itot = load_image(files, fl, mode=mode, sky=sky, combined=combined,
                                       flatfield=flatfield, errorprop=errorprop)
    one_d = wc.get("wctype") == "1D"
    # no ATRAN without a wavecal file (XSPEXTOOL_MKPROFILES)
    aw, at = atran if (atran is not None and not one_d) else (None, None)
    say("Making the spatial profiles")
    # the spatial maps are needed by optimal extraction only (the fits take a while)
    want = optimal and (source == "point" or one_d)
    prof = mkspatprof2d(work, fl["omask"], wc["wavecal"], wc["spatcal"],
                        fl["slith_arc"] / fl["slith_pix"], ybuffer, aw, at, atmosthresh,
                        coeffs=want)
    prof, spc, mspc = prof if want else (prof, None, None)
    info = dict(modename=fl["modename"], slith_pix=fl["slith_pix"], slith_arc=fl["slith_arc"],
                slitw_pix=fl["slitw_pix"], slitw_arc=fl["slitw_arc"], rp=fl["rp"],
                flat=flatname, wavecal=wavecalname)
    doorders = np.ones(no, dtype=int)
    if source == "point":
        apradius, psfradius = f32(apradius), f32(psfradius)       # MC_CFLD(...,4)
        psbg = (f32(psbg[0]), f32(psbg[1]), psbg[2])
        if ps_orders is not None:
            doorders = np.array([1 if o in ps_orders else 0 for o in allorders])
            if doorders.sum() < 1:
                raise ValueError("Must select at least one Order.")
        pts = []
        if trace_in:
            # 'Import Trace' (XSPEXTOOL_PSAPFIND): orders and coefficients
            # from the file, positions their mean, signs from the profile.
            torders, npeaks, tc = read_trace(trace_in)
            doorders = np.array([1 if o in torders else 0 for o in allorders])
            paths = tracetoap2d(fl["omask"], wc["wavecal"], wc["spatcal"], tc, npeaks,
                                allorders, doorders)
            pos = np.full((no, npeaks), np.nan, dtype=f32)
            l = 0
            for i in range(no):
                if doorders[i]:
                    for j in range(npeaks):
                        pos[i, j] = paths[l]["appos"]
                        l += 1
            check_ps_inputs(fl["slith_arc"], apradius, psbgsub, psbg, optimal, psfradius, npeaks)
        else:
            say("Finding the apertures")
            if findmode == "auto":
                check_ps_inputs(fl["slith_arc"], apradius, psbgsub, psbg, optimal, psfradius,
                                npeaks, tracestep=tracestep if trace else None,
                                sumap=sumap if trace else None, tracedeg=tracedeg if trace else None)
                pos, _ = findpeaks(prof, npeaks, doorders, mode="auto")
            else:
                g = np.asarray(guess, dtype=f32)                  # FLOAT(STRSPLIT(...))
                if g.ndim == 1:                                   # one list for every order
                    g = np.repeat(g[None, :], no, axis=0)
                npeaks = g.shape[1]
                if findmode == "fixed":
                    trace = False                                 # ticks Use Ap Positions
                check_ps_inputs(fl["slith_arc"], apradius, psbgsub, psbg, optimal, psfradius,
                                npeaks, positions=g, tracestep=tracestep if trace else None,
                                sumap=sumap if trace else None, tracedeg=tracedeg if trace else None)
                if findmode == "guess":
                    pos, _ = findpeaks(prof, npeaks, doorders, mode="guess", guess=g)
                else:
                    pos, _ = findpeaks(prof, npeaks, doorders, mode="fixed", fixed=g)
            sel = np.where(doorders == 1)[0]
            if trace:
                say("Tracing the apertures")
                tc, pts = tracespec2d(work, wc["wavecal"], wc["spatcal"], fl["edgecoeffs"][sel],
                                      [fl["xranges"][i] for i in sel], pos[sel], tracestep,
                                      sumap, f32(winthresh), f32(sigthresh), tracedeg)
            else:
                # Use Ap Positions: an FLTARR [pos, 0] per aperture, taking
                # appos[*, i] with i counting the selected orders (so with
                # an order deselected below a selected one it reads the
                # wrong order's positions, as the IDL does).
                tc = np.array([[pos[i, j], 0.0] for i in range(len(sel)) for j in range(npeaks)],
                              dtype=f32)
                if np.isnan(tc).any():
                    raise ValueError("Use Ap Positions with these PS Orders reads the positions "
                                     "of a deselected order (NaN); IDL's extraction never "
                                     "returns from that. Select the orders from the lowest up.")
            if trace_out:
                write_trace(trace_out, [allorders[i] for i in sel], npeaks, tc, fl["modename"])
        # XSPEXTOOL_PSDEFINEAP: apertures inside the slit, signs from the profile
        for i in range(no):
            if doorders[i]:
                t = np.concatenate([pos[i] + apradius, pos[i] - apradius])
                if (t < 0).any() or (t > fl["slith_arc"]).any():
                    raise ValueError("Some of the extraction apertures extend beyond the slit. "
                                     "Please lower the aperture radius.")
        sign = psdefine_signs(prof, pos, doorders)
        if signs is not None:                 # XSPEXTOOL_OVERRIDESIGNS
            sign = np.asarray(signs, dtype=int)
        orders = [allorders[i] for i in range(no) if doorders[i]]
        apradii = np.full((len(orders), npeaks), apradius, dtype=f32)
        smap = None
        if optimal:
            c = mspc if medprofile else spc
            smap = [c[o] for o in orders]
        say("Extracting")
        if one_d:
            sel = [i for i in range(no) if doorders[i]]
            r1 = extpsspec(work, var, wc["wavecal"], [fl["edgecoeffs"][i] for i in sel], tc, npeaks,
                           [fl["xranges"][i] for i in sel], fl["slith_arc"],
                           np.full(npeaks, apradius, dtype=f32), sign,
                           psfwidth=psfradius if optimal else None, spatcoeffs=smap,
                           bgorder=psbg[2] if psbgsub else None,
                           psbginfo=[psbg[0], psbg[1]] if psbgsub else None,
                           bdpxmask=load_bdpxmask(fl["rotation"]),
                           bitmask=load_bitmask(files, fl, mode, sky, combined),
                           bdpixthresh=f32(bdpixthresh))
            spec = {(o, a): r1[n][a] for n, o in enumerate(orders) for a in r1[n]}
        else:
            spec = extspec2d(work, var, fl["omask"], orders, wc["wavecal"], wc["spatcal"], fl["ps"],
                             sign, tc, apradii,
                             psbginfo=[psbg[0], psbg[1]] if psbgsub else None,
                             bgdeg=psbg[2] if psbgsub else None,
                             psfradius=psfradius if optimal else None, smapinfo=smap)
        info.update(apradius=[apradius], psbginfo=list(psbg) if psbgsub else None,
                    psfradius=psfradius if optimal else None)
    else:
        if xs_orders is not None:
            doorders = np.array([1 if o in xs_orders else 0 for o in allorders])
            if doorders.sum() < 1:
                raise ValueError("Must select at least one Order.")
        orders = [allorders[i] for i in range(no) if doorders[i]]
        pos = np.array([list(xs_positions)] * no, dtype=float)
        naps = pos.shape[1]
        radii = list(xs_radii) if xs_radii is not None else [1.0] * naps
        tc = np.array([[pos[i, j], 0.0] for i in range(len(orders)) for j in range(naps)])
        pts = []
        sign = np.ones(naps, dtype=int)
        apradii = np.array([radii] * len(orders), dtype=float)
        say("Extracting")
        if one_d:
            sel = [i for i in range(no) if doorders[i]]
            smap = None
            if optimal:
                # built from the PS Orders (XSPEXTOOL_XSEXTRACTSPEC), used by position
                c = mspc if medprofile else spc
                psel = allorders if ps_orders is None else [o for o in allorders if o in ps_orders]
                smap = [c[o] for o in psel]
            r1 = extxsspec(work, var, wc["wavecal"], [fl["edgecoeffs"][i] for i in sel],
                           tc.astype(f32), naps, [fl["xranges"][i] for i in sel], fl["slith_arc"],
                           np.asarray(radii, dtype=f32),
                           bgr=[np.asarray(r_, dtype=f32) for r_ in xsbg] if xsbgsub else None,
                           bgorder=int(xsbgdeg) if xsbgsub else None, spatcoeffs=smap,
                           bdpxmask=load_bdpxmask(fl["rotation"]),
                           bitmask=load_bitmask(files, fl, mode, sky, combined),
                           bdpixthresh=f32(bdpixthresh))
            spec = {(o, a): r1[n][a] for n, o in enumerate(orders) for a in r1[n]}
        else:
            spec = extspec2d(work, var, fl["omask"], orders, wc["wavecal"], wc["spatcal"], fl["ps"],
                             sign, tc, apradii,
                             xsbginfo=[list(r) for r in xsbg] if xs_background else None,
                             bgdeg=xsbgdeg if xs_background else None)
        bgtxt = ",".join(f"{a:g}-{b:g}" for a, b in xsbg)
        # XSPEXTOOL_WRITESPEC passes the background regions in A and A-Sky only,
        # and only with XS BG Subtraction on
        info.update(apradius=radii,
                    xsbginfo=[bgtxt, str(xsbgdeg)] if (mode != "A-B" and xsbgsub) else None,
                    xs_background=xs_background)
    good = [i for i in range(no) if doorders[i]]
    info["disp"] = [wc["disp"][i] for i in good] if not one_d else None
    info["pixels"] = one_d
    written = []
    names = [os.path.basename(f) for f in files] if files else [os.path.basename(combined or "")]
    if mode != "A-B":
        # A and A-Sky write every aperture, sign applied, to one file.
        sub = {(o, j + 1): spec[(o, j + 1)] for o in orders for j in range(len(sign))}
        d = dict(info)
        d.update(aimage=names[0], sky=os.path.basename(sky) if (mode == "A-Sky" and sky) else "None",
                 itot=itot[0], appos=[list(pos[i]) for i in good])
        write_spectrum(opaths[0], sub, orders, len(sign), hdrs[0], d)
        written.append(opaths[0])
    for k, (want, hidx) in enumerate(((1, 0), (-1, 1)) if mode == "A-B" else ()):
        aps = [j for j in range(len(sign)) if sign[j] == want]
        if not aps:
            continue
        sub = {}
        for o in orders:
            for n, j in enumerate(aps):
                sub[(o, n + 1)] = spec[(o, j + 1)]
        d = dict(info)
        d.update(aimage=names[hidx], sky=names[1 - hidx],
                 itot=itot[hidx], appos=[[pos[i, j] for j in aps] for i in good])
        write_spectrum(opaths[k], sub, orders, len(aps), hdrs[hidx], d)
        written.append(opaths[k])
    return dict(written=written, profiles=prof, positions=pos, signs=sign, trace=tc,
                trace_points=pts, work=work, var=var, spectra=spec, headers=hdrs,
                orders=orders, spatcoeffs=spc, medspatcoeffs=mspc)
