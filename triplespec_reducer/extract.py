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

Extended sources: in the IDL, MC_EXTSPEC2D asks MC_MKAPMASK1D for the
background regions as XSBGINFO=xsbgino (a misspelling of xsbginfo), and
XSPEXTOOL_XSEXTRACTSPEC passes BGDEG=bgorder where the degree was read into
bgfit. Either alone switches the background fit off, so in 2D mode the
regions a user types for an extended source are drawn but never used. That
is reproduced by default. `xs_background=True` applies them exactly as the
point-source path does -- the subtraction Spextool was written to make.
"""

import os

import numpy as np
from astropy.io import fits

from triplespec_reducer import idl as I
from triplespec_reducer.frames import read_ts4

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


def write_trace(path, orders, naps, tracecoeffs, obsmode="XD"):
    """The trace file XSPEXTOOL_PSTRACESPEC writes to cal/."""
    with open(path, "w") as fh:
        fh.write(f"OBSMODE={obsmode}\n")
        fh.write("ORDERS=" + ",".join(str(o) for o in orders) + "\n")
        fh.write(f"NAPS={naps}\n")
        k = 0
        for o in orders:
            for j in range(naps):
                fh.write(f"ORDER{o:02d}_AP{j + 1:02d}= " + "  ".join(repr(float(c)) for c in tracecoeffs[k]) + "\n")
                k += 1


def read_trace(path):
    """MC_READTRACE."""
    d = {}
    with open(path) as fh:
        for line in fh:
            if "=" in line:
                k, v = line.split("=", 1)
                d[k.strip()] = v.strip()
    orders = [int(t) for t in d["ORDERS"].split(",")]
    naps = int(d["NAPS"])
    tc = [[float(t) for t in d[f"ORDER{o:02d}_AP{j + 1:02d}"].split()] for o in orders for j in range(naps)]
    return orders, naps, np.array(tc)


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
                 ttrans=None, atmosthresh=None, orders=None):
    """
    MC_MKSPATPROF2D: the mean spatial profile of each order.

    The order is rectified onto (lambda, arcsec); each wavelength has its
    median along the slit removed; wavelengths where the atmosphere
    transmits less than `atmosthresh` are blanked; the profile is the
    median along wavelength, normalized by the sum of its absolute values.
    Returns {order: (ygrid, profile)}.
    """
    orders = orders or sorted(int(o) for o in np.unique(omask) if o != 0)
    out = {}
    for o in orders:
        xgrid, ygrid, ix, iy, _ = _order_grid(omask, wavecal, spatcal, o, ds)
        rect = I.interpolate2d(img, ix, iy, cubic=-0.5, missing=np.nan)
        ny = ygrid.size
        ygrid = ygrid[ybuffer:ny - ybuffer]
        rect = rect[ybuffer:ny - ybuffer, :]
        bg = I.median(rect, axis=0)
        rect = rect - bg[None, :]
        if wtrans is not None:
            mask = np.ones(xgrid.size, dtype=np.float32)            # FLTARR
            rtrans = I.linterp(wtrans, ttrans, xgrid, missing=1)
            mask[rtrans < atmosthresh] = np.nan
            rect = rect * mask[None, :]
        prof = I.median(rect, axis=1)
        a = np.abs(prof)
        prof = (prof / I.ftotal(a[np.isfinite(a)])).astype(np.float32)   # TOTAL(/NAN), float
        # IDL keeps [[ygrid],[prof]]: DOUBLE, so the slit grid keeps its precision.
        out[o] = (ygrid, prof)
    return out


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
        # The IDL profile is [[ygrid],[prof]]: DOUBLE.
        x, p = np.asarray(profiles[o][0], float), np.asarray(profiles[o][1], float)
        if mode == "guess":
            y = p - I.median(p, even=True)
            for j in range(npeaks):
                xpix = I.tabinv(x, guess[i][j])
                _, fit = I.mpfitpeak(x, y, nterms=4, estimates=np.array(
                    [y[I.itrunc(xpix)], guess[i][j], 0.5, 0], dtype=float))
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
                                     estimates=np.array([y[k], xx[k], 1, 0], dtype=float))
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
    """MC_TRACETOAP2D: an aperture's path back on the detector, for display."""
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
            out.append(dict(order=o, aperture=j + 1, x=ix, y=iy, appos=float(np.mean(tr))))
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
    slit = np.asarray(slit, dtype=float)
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


def extspec2d(img, var, omask, orders, wavecal, spatcal, ds, apsign, tracecoeffs,
              apradii, psbginfo=None, xsbginfo=None, bgdeg=None):
    """
    MC_EXTSPEC2D (sum extraction; optimal extraction is off for TS4).

    Each order is cut into cells of one resampled wavelength step by `ds`
    arcsec; every cell's flux is the exact-area sum of the pixels under it.
    The background fit, when asked for, is a robust polynomial of degree
    `bgdeg` across the slit at each wavelength, its uncertainty added to
    the variance. Returns {(order, aperture): [wave, flux, error]}.
    """
    nrows, ncols = img.shape
    naps = apradii.shape[1]
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
        for j in range(nw):
            px = np.stack([ix[j, 0:ns], ix[j, 1:ns + 1], ix[j + 1, 1:ns + 1], ix[j + 1, 0:ns]], axis=1) + 0.5
            py = np.stack([iy[j, 0:ns], iy[j, 1:ns + 1], iy[j + 1, 1:ns + 1], iy[j + 1, 0:ns]], axis=1) + 0.5
            slit_flx, slit_var = _cell_sums(px.astype(np.float64), py.astype(np.float64), imgf, varf)
            slitmask = mkapmask1d(slit_arc, trace_arc[j], apradii[i], psbginfo=psbginfo,
                                  xsbginfo=xsbginfo)
            if bgdeg is not None:
                z = np.where((slitmask == -1.0) & (slit_var[:ns] != 0))[0]
                if z.size >= 2:
                    c, _, cvar = I.robustpoly1d(slit_arc[z], slit_flx[z], bgdeg, 4, 0.1,
                                                yerr=np.sqrt(slit_var[z]), want_covar=True)
                    bg, yvar = I.poly1d(slit_arc, c, cvar)
                    # FLOAT minus the DOUBLE fit: the slit stays double, as in IDL
                    slit_flx = slit_flx[:ns] - bg
                    slit_var = slit_var[:ns] + yvar
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
    keywords every later Spextool step reads.
    """
    keys = [(o, a) for o in orders for a in range(1, naps + 1)]
    npix = max(spectra[k].shape[1] for k in keys)
    cube = np.full((len(keys), 3, npix), np.nan)
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
    I.fxaddpar(hdr, "WAVECAL", info.get("wavecal", ""), " Wavecal file used for wavelength calibration")
    I.fxaddpar(hdr, "WAVETYPE", "Vacuum", " Wavelength type (air or vacuum)")
    if info.get("tracefile"):
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
        if info.get("xs_background"):         # the opt-in, not an IDL keyword
            I.fxaddpar(hdr, "XSBGSUB", True, " Extended-source background actually subtracted")
    I.fxaddpar(hdr, "XUNITS", "um", " Units of the X axis")
    I.fxaddpar(hdr, "YUNITS", "DN / s", " Units of the Y axis")
    I.fxaddpar(hdr, "XTITLE", "!7k!5 (!7l!5m)", " IDL X title")
    I.fxaddpar(hdr, "YTITLE", "f (!5DN s!u-1!N)", " IDL Y title")
    I.fxaddpar(hdr, "VERSION", info.get("version", "TSR-1.0"), " Spextool version")
    for i, o in enumerate(orders):
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


def reduce_pair(files, flatinfo, wcinfo, opaths, source="point", atran=None,
                npeaks=2, findmode="auto", guess=None, apradius=1.0,
                psbg=(2.5, 5.0, 1), xs_positions=None, xs_radii=None,
                xsbg=((2., 6.), (32., 36.)), xsbgdeg=0, xs_background=False,
                trace=True, mode="A-B", tracestep=5, sumap=5, winthresh=3.0,
                sigthresh=0.1, tracedeg=1, ybuffer=1, atmosthresh=0.4,
                flatname="", wavecalname="", progress=None, signs=None,
                combined=None, sky=None, trace_out=None, trace_in=None,
                flatfield=True, errorprop=True):
    """
    What xspextool does for one pair of frames with "Do All Steps", for a
    point or an extended source. Writes opaths[0] (and opaths[1] for the
    negative apertures of an A-B point source) and returns everything
    computed on the way, for display.
    """
    say = progress or (lambda *_: None)
    fl, wc = flatinfo, wcinfo
    orders = fl["orders"]
    say("Loading the images")
    work, var, hdrs, itot = load_image(files, fl, mode=mode, sky=sky, combined=combined,
                                       flatfield=flatfield, errorprop=errorprop)
    aw, at = atran if atran is not None else (None, None)
    say("Making the spatial profiles")
    prof = mkspatprof2d(work, fl["omask"], wc["wavecal"], wc["spatcal"],
                        fl["slith_arc"] / fl["slith_pix"], ybuffer, aw, at, atmosthresh)
    info = dict(modename=fl["modename"], slith_pix=fl["slith_pix"], slith_arc=fl["slith_arc"],
                slitw_pix=fl["slitw_pix"], slitw_arc=fl["slitw_arc"], rp=fl["rp"],
                disp=wc["disp"], flat=flatname, wavecal=wavecalname)
    if source == "point" and trace_in:
        # 'Import Trace': the apertures follow a stored trace; their signs
        # come from the profile at the trace's mean position.
        _, npeaks, tc = read_trace(trace_in)
        paths = tracetoap2d(fl["omask"], wc["wavecal"], wc["spatcal"], tc, npeaks, orders)
        pos = np.array([[paths[i * npeaks + j]["appos"] for j in range(npeaks)]
                        for i in range(len(orders))])
        _, sgn = findpeaks(prof, npeaks, mode="fixed", fixed=pos)
        sign = np.asarray(signs if signs is not None else sgn, dtype=int)
        pts = []
        apradii = np.full((len(orders), npeaks), float(apradius))
        say("Extracting")
        spec = extspec2d(work, var, fl["omask"], orders, wc["wavecal"], wc["spatcal"], fl["ps"],
                         sign, tc, apradii, psbginfo=list(psbg[:2]), bgdeg=psbg[2])
        info.update(apradius=[float(apradius)], psbginfo=list(psbg), tracefile=os.path.basename(trace_in))
    elif source == "point":
        say("Finding the apertures")
        if findmode == "auto":
            pos, _ = findpeaks(prof, npeaks, mode="auto")
        elif findmode == "guess":
            pos, _ = findpeaks(prof, npeaks, mode="guess", guess=guess)
        else:
            pos, _ = findpeaks(prof, npeaks, mode="fixed", fixed=guess)
        if trace:
            say("Tracing the apertures")
            tc, pts = tracespec2d(work, wc["wavecal"], wc["spatcal"], fl["edgecoeffs"],
                                  fl["xranges"], pos, tracestep, sumap, winthresh,
                                  sigthresh, tracedeg)
        else:
            tc = np.array([[pos[i, j], 0.0] for i in range(len(orders)) for j in range(npeaks)])
            pts = []
        sign = psdefine_signs(prof, pos)
        if signs is not None:                 # XSPEXTOOL_OVERRIDESIGNS
            sign = np.asarray(signs, dtype=int)
        if trace_out:
            write_trace(trace_out, orders, npeaks, tc, fl["modename"])
        apradii = np.full((len(orders), npeaks), float(apradius))
        say("Extracting")
        spec = extspec2d(work, var, fl["omask"], orders, wc["wavecal"], wc["spatcal"], fl["ps"],
                         sign, tc, apradii, psbginfo=list(psbg[:2]), bgdeg=psbg[2])
        info.update(apradius=[float(apradius)], psbginfo=list(psbg))
    else:
        pos = np.array([list(xs_positions)] * len(orders), dtype=float)
        naps = pos.shape[1]
        radii = list(xs_radii) if xs_radii is not None else [1.0] * naps
        tc = np.array([[pos[i, j], 0.0] for i in range(len(orders)) for j in range(naps)])
        pts = []
        sign = np.ones(naps, dtype=int)
        apradii = np.array([radii] * len(orders), dtype=float)
        say("Extracting")
        spec = extspec2d(work, var, fl["omask"], orders, wc["wavecal"], wc["spatcal"], fl["ps"],
                         sign, tc, apradii,
                         xsbginfo=[list(r) for r in xsbg] if xs_background else None,
                         bgdeg=xsbgdeg if xs_background else None)
        bgtxt = ",".join(f"{a:g}-{b:g}" for a, b in xsbg)
        info.update(apradius=radii, xsbginfo=[bgtxt, str(xsbgdeg)], xs_background=xs_background)
    written = []
    names = [os.path.basename(f) for f in files] if files else [os.path.basename(combined or "")]
    if mode != "A-B":
        # A and A-Sky write every aperture, sign applied, to one file.
        sub = {(o, j + 1): spec[(o, j + 1)] for o in orders for j in range(len(sign))}
        d = dict(info)
        d.update(aimage=names[0], sky=os.path.basename(sky) if (mode == "A-Sky" and sky) else "None",
                 itot=itot[0], appos=[list(pos[i]) for i in range(len(orders))])
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
                 itot=itot[hidx], appos=[[pos[i, j] for j in aps] for i in range(len(orders))])
        write_spectrum(opaths[k], sub, orders, len(aps), hdrs[hidx], d)
        written.append(opaths[k])
    return dict(written=written, profiles=prof, positions=pos, signs=sign, trace=tc,
                trace_points=pts, work=work, var=var, spectra=spec, headers=hdrs)
