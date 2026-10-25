"""
The normalized flat field: TS41CALS_MKFLAT and the routines it calls.

    lamps-on and lamps-off frames  ->  KNA_READTS4FITS
    lamps-off combined, subtracted ->  MC_COMBIMGS (mean; no clipping, as
                                       the IDL passes an undefined THRESH)
    orders found on the mean       ->  TS4_FINDORDERS (the outrigger edges)
    frames scaled, combined        ->  KA_SCALEFLATS, MC_COMBIMGS
    each order normalized          ->  MC_NORMSPECFLAT (MC_RECTORDER,
                                       MC_FITERPOLATE)
    written                        ->  MC_WRITEFLAT

Images here are numpy [row, col] in the rotated orientation: columns are
the dispersion axis (IDL x) and rows the spatial axis (IDL y).
"""

import numpy as np
from astropy.io import fits

from triplespec_reducer import idl as I
from triplespec_reducer.frames import read_ts4, unrotate, rotate
from triplespec_reducer.instrument import (read_modefile, SLITW_ARC, SLITW_PIX,
                                           RESOLVINGPOWER)


# ------------------------------------------------------------- order edges

def ts4_findorders(image, guesspos, sranges, step, slith_pix, slith_pix_range,
                   degree, bufpix, comwidth, qa=None):
    """
    TS4_FINDORDERS: trace the two outrigger peaks of every order.

    The flat of each TripleSpec order carries a bright narrow peak just
    outside each edge. Every `step` columns a Gaussian is fitted to each
    peak, starting from the order's guess position and walking right then
    left; a pair is kept when the peaks are slith_pix_range apart. The
    peaks, robustly fitted, are then moved inward by 13% of the slit height
    to give the edges actually used. Returns edgecoeffs [norders, 2, deg+1]
    (bottom, top) and xranges [norders, 2].

    With `qa` a list, one dict per order is appended holding what the IDL
    draws on the combined flat while it works (WID=): guess, the guess
    position (red dot); top / bot, the fitted outrigger centres it plots
    (green '+', NaN where a walk had stopped); and for the bottom and top
    edge (0, 1) the fitted curve over the order's columns (red) and the
    points the robust fit rejected (blue '+').
    """
    nrows, ncols = image.shape
    norders = len(guesspos)
    edgecoeffs = np.zeros((norders, 2, degree + 1))
    xranges = np.zeros((norders, 2), dtype=int)
    for i in range(norders):
        start, stop = int(sranges[i][0]), int(sranges[i][1])
        starts = start + step - 1
        stops = stop - step + 1
        scols = (np.arange(I.itrunc((stops - starts) / step) + 1, dtype=np.float32)
                 * np.float32(step) + np.float32(starts))
        nscols = scols.size
        # FLTARR and REPLICATE(!VALUES.F_NAN) in the IDL: every traced value
        # is stored in single precision.
        oedges = np.full((nscols, 2), np.nan, dtype=np.float32)
        cen = np.full(nscols, np.nan, dtype=np.float32)
        idx = I.tabinv(scols, float(guesspos[i][0]))
        gidx = I.iround(idx)
        cen[gidx - degree:gidx + degree + 1] = guesspos[i][1]
        oedges[gidx - degree:gidx + degree + 1, 0] = guesspos[i][1] - slith_pix / 2.0
        oedges[gidx - degree:gidx + degree + 1, 1] = guesspos[i][1] + slith_pix / 2.0

        pts_top, pts_bot = [], []

        def walk(ks, pdeg):
            dotop = dobot = True
            for k in ks:
                fcol = image[:, int(scols[k])]
                coeff = I.poly_fit1d(scols, cen, pdeg)
                y_guess = min(max(bufpix, I.poly(scols[k], coeff)), nrows - bufpix - 1)
                ctop = I.poly_fit1d(scols, oedges[:, 1], pdeg)
                cbot = I.poly_fit1d(scols, oedges[:, 0], pdeg)
                slith_guess = I.poly(scols[k], ctop) - I.poly(scols[k], cbot)
                com_top = com_bot = np.nan
                if dotop:
                    pix = (np.arange(comwidth) - I.iround((comwidth - 1) / 2.)
                           + I.iround(y_guess + slith_guess / 2.0))
                    pix = pix[(pix > bufpix) & (pix < nrows - bufpix - 1)]
                    com_top = I.gaussfit(pix, fcol[pix], nterms=4)[1][1]
                if dobot:
                    pix = (np.arange(comwidth) - I.iround((comwidth - 1) / 2.)
                           + I.iround(y_guess - slith_guess / 2.0))
                    pix = pix[(pix > bufpix) & (pix < nrows - bufpix - 1)]
                    com_bot = I.gaussfit(pix, fcol[pix], nterms=4)[1][1]
                if dotop and dobot and np.isfinite(com_bot) and np.isfinite(com_top):
                    sep = abs(com_bot - com_top)
                    if slith_pix_range[0] < sep < slith_pix_range[1]:
                        oedges[k] = [com_bot, com_top]
                        cen[k] = (com_bot + com_top) / 2.
                    else:
                        continue                       # IDL: goto, cont
                else:
                    com_bot = com_top = np.nan
                    oedges[k] = [np.nan, np.nan]
                    cen[k] = y_guess
                if dotop:
                    pts_top.append((float(scols[k]), float(com_top)))
                if dobot:
                    pts_bot.append((float(scols[k]), float(com_bot)))
                if com_top <= bufpix or com_top >= (nrows - 1 - bufpix):
                    dotop = False
                if com_bot <= bufpix or com_bot > (nrows - 1 - bufpix):
                    dobot = False
                if not dotop and not dobot:
                    break

        # Rightward the guess polynomials are degree-1 (>= 1), leftward
        # degree (>= 1), exactly as the two loops in the IDL differ.
        walk(range(gidx + 1, nscols), max(1, degree - 1))
        walk([gidx - j for j in range(gidx + 1)], max(1, degree))

        fedges = np.zeros_like(oedges)
        for j in range(2):
            c, _ = I.robustpoly1d(scols, oedges[:, j], degree, 3, 0.01)
            fedges[:, j] = I.poly(scols, c)
        outrigger_offset = np.float32(0.13) * (fedges[:, 1] - fedges[:, 0])
        edges = np.empty_like(oedges)
        edges[:, 1] = oedges[:, 1] - outrigger_offset
        edges[:, 0] = oedges[:, 0] + outrigger_offset
        x = np.arange(stop - start + 1, dtype=float) + start
        curves, rejected = [], []
        for j in range(2):
            c, gb = I.robustpoly1d(scols, edges[:, j], degree, 3, 0.01)
            edgecoeffs[i, j] = c
            curves.append(I.poly(x, c))
            z = gb == 0
            rejected.append((scols[z].astype(float), edges[z, j].astype(float)))
        bot = I.poly(x, edgecoeffs[i, 0])
        top = I.poly(x, edgecoeffs[i, 1])
        z = (top > 0.0) & (top < nrows - 1) & (bot > 0) & (bot < nrows - 1)
        xranges[i] = [int(x[z].min()), int(x[z].max())]
        if qa is not None:
            qa.append(dict(guess=(float(guesspos[i][0]), float(guesspos[i][1])),
                           top=np.array(pts_top).reshape(-1, 2),
                           bot=np.array(pts_bot).reshape(-1, 2),
                           x=x, edges=curves, rejected=rejected))
    return edgecoeffs, xranges


def order_mask(ncols, nrows, orders, edgecoeffs, xranges):
    """The OMASK MC_READFLAT builds: each pixel holds its order number."""
    omask = np.zeros((nrows, ncols), dtype=int)
    for i, o in enumerate(orders):
        x = np.arange(xranges[i][0], xranges[i][1] + 1)
        bot = I.poly1d(x, edgecoeffs[i, 0])[0]
        top = I.poly1d(x, edgecoeffs[i, 1])[0]
        for j, xx in enumerate(x):
            if top[j] > nrows - 0.5 or bot[j] < -0.5:
                continue
            omask[I.itrunc(bot[j]):I.itrunc(top[j]) + 1, xx] = o
    return omask


def ka_scaleflats(data, var, edgecoeffs, xranges):
    """KA_SCALEFLATS: scale every flat to the median in-order level."""
    n, nrows, ncols = data.shape
    mask = np.zeros((nrows, ncols), dtype=bool)
    for i in range(len(xranges)):
        x = np.arange(xranges[i][0], xranges[i][1] + 1)
        bot = I.poly1d(x, edgecoeffs[i, 0])[0]
        top = I.poly1d(x, edgecoeffs[i, 1])[0]
        for j, xx in enumerate(x):
            if top[j] > nrows - 0.5 or bot[j] < -0.5:
                continue
            mask[I.itrunc(bot[j]):I.itrunc(top[j]) + 1, xx] = True
    medflux = np.array([I.median(d[mask], even=True) for d in data], dtype=np.float32)
    scale = (np.float32(I.median(medflux, even=True)) / medflux).astype(np.float32)
    data = data * scale[:, None, None]
    var = var * (I.sq(scale))[:, None, None]
    return data, var


# ------------------------------------------------------------- normalization

def rectorder(image, edgecoeffs, xrange, slith_arc, step, oversamp, ybuffer):
    """
    MC_RECTORDER: resample one order onto (column, slit arcsec).

    Returns map [nslit, ncolumns] and yarc, the slit coordinate of each row.
    """
    start, stop = int(xrange[0]), int(xrange[1])
    x = np.arange(stop - start + 1, dtype=float) + start
    starts = start + step - 1
    stops = stop - step + 1
    numstep = I.itrunc((stops - starts) / step) + 1
    scolumn = np.arange(numstep, dtype=float) * step + starts
    botedge = I.poly1d(scolumn, edgecoeffs[0])[0]
    topedge = I.poly1d(scolumn, edgecoeffs[1])[0]
    dif = topedge - botedge
    pixtoarc = dif / slith_arc
    avedif = I.iround(np.nanmin(dif))
    nslit = int(oversamp * avedif)
    # FINDGEN times the FLOAT slit height: single precision, as in the IDL.
    yarc_map = (np.arange(nslit, dtype=np.float32) * np.float32(slith_arc)
                / np.float32(nslit - 1)).astype(np.float32)
    tmp = yarc_map.copy()
    tmp[0:ybuffer] = tmp[ybuffer]
    tmp[(nslit - 1 - ybuffer - 1):nslit] = tmp[nslit - 1 - ybuffer - 2]
    ypix = tmp[:, None] * pixtoarc[None, :] + botedge[None, :]
    xpix = np.broadcast_to(x[None, :], ypix.shape)
    m = I.interpolate2d(image, xpix, ypix, cubic=-0.5)
    return m.astype(np.float32), yarc_map


_BICU = np.array([
    [1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
    [0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
    [-3., 0., 0., 3., 0., 0., 0., 0., -2., 0., 0., -1., 0., 0., 0., 0.],
    [2., 0., 0., -2., 0., 0., 0., 0., 1., 0., 0., 1., 0., 0., 0., 0.],
    [0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.],
    [0., 0., 0., 0., -3., 0., 0., 3., 0., 0., 0., 0., -2., 0., 0., -1.],
    [0., 0., 0., 0., 2., 0., 0., -2., 0., 0., 0., 0., 1., 0., 0., 1.],
    [-3., 3., 0., 0., -2., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
    [0., 0., 0., 0., 0., 0., 0., 0., -3., 3., 0., 0., -2., -1., 0., 0.],
    [9., -9., 9., -9., 6., 3., -3., -6., 6., -6., -3., 3., 4., 2., 1., 2.],
    [-6., 6., -6., 6., -4., -2., 2., 4., -3., 3., 3., -3., -2., -1., -1., -2.],
    [2., -2., 0., 0., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
    [0., 0., 0., 0., 0., 0., 0., 0., 2., -2., 0., 0., 1., 1., 0., 0.],
    [-6., 6., -6., 6., -3., -3., 3., 3., -4., 4., 2., -2., -2., -2., -1., -1.],
    [4., -4., 4., -4., 2., 2., -2., -2., 2., -2., -2., 2., 1., 1., 1., 1.]])


def bicuval(z, dz1, dz2, dz12, xl, xu, yl, yu, x, y):
    """
    MC_BICUVAL with MC_BICUCOEFFS, in IDL's precision: single for float
    nodes and pixel grids (FINDGEN), as MC_FITERPOLATE passes them.
    """
    P = np.float64 if I.isf64(x) else np.float32
    nd1, nd2 = (xu - xl), (yu - yl)
    v = np.concatenate([np.asarray(z, P), (np.asarray(dz1, P) * P(nd1)).astype(P),
                        (np.asarray(dz2, P) * P(nd2)).astype(P),
                        (np.asarray(dz12, P) * P(nd1 * nd2)).astype(P)])
    w = _BICU.astype(P)
    c = np.zeros(16, dtype=P)
    for m in range(16):
        c = (c + w[:, m] * v[m]).astype(P)
    c = c.reshape(4, 4)
    t = ((np.asarray(x, P) - P(xl)) / P(nd1)).astype(P)
    u = ((np.asarray(y, P) - P(yl)) / P(nd2)).astype(P)
    nz = np.zeros_like(t)
    for i in range(3, -1, -1):
        nz = t * nz + c[i, 0] + u * (c[i, 1] + u * (c[i, 2] + u * c[i, 3]))
    return nz


_QEXP = [(0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1)]


def quadfit(img):
    """
    MC_QUADFIT /JUSTFIT: c0 + c1 x + c2 y + c3 x^2 + c4 y^2 + c5 x y on an [y, x] patch.

    The design matrix is built in float, as the IDL builds it from FINDGEN, and
    the normal equations are solved by MC_GAUSSJ; numpy's solver rounds them
    differently.
    """
    ysize, xsize = img.shape
    x, y = np.meshgrid(np.arange(xsize, dtype=np.float32), np.arange(ysize, dtype=np.float32))
    x, y, z = x.ravel(), y.ravel(), img.ravel()
    g = np.isfinite(z)
    x, y, z = x[g], y[g], z[g].astype(float)
    at = [(x ** ex * y ** ey).astype(np.float32).astype(float) for ex, ey in _QEXP]
    alpha = np.array([[I.itotal(at[i] * at[j]) for j in range(6)] for i in range(6)])
    beta = np.array([I.itotal(a * z) for a in at])
    return I._gaussj(alpha, beta)[0]


def fiterpolate(image, ncg, nrg):
    """
    MC_FITERPOLATE: a smooth surface through an ncg x nrg grid of cells.

    At every grid node a quadratic surface fitted to the surrounding block
    gives the value, both slopes and the cross derivative; bicubic patches
    join the nodes. `image` is [rows, cols]; IDL's x is the column.
    """
    nrows, ncols = image.shape
    fitimg = np.zeros((nrows, ncols), dtype=np.float32)
    nxgrid, nygrid = ncg + 1, nrg + 1
    gx = I.iround(np.arange(nxgrid) * (ncols - 1) / float(nxgrid - 1))
    gy = I.iround(np.arange(nygrid) * (nrows - 1) / float(nygrid - 1))
    vals = np.zeros((nxgrid * nygrid, 4), dtype=np.float32)
    idx = 0
    for i in range(nxgrid):
        x1 = (gx[max(i - 1, 0)] + gx[i]) // 2
        x2 = (gx[i] + gx[min(nxgrid - 1, i + 1)]) // 2
        for j in range(nygrid):
            y1 = (gy[max(j - 1, 0)] + gy[j]) // 2
            y2 = (gy[j] + gy[min(nygrid - 1, j + 1)]) // 2
            c = quadfit(image[y1:y2 + 1, x1:x2 + 1])
            x = gx[i] - x1
            y = gy[j] - y1
            vals[idx] = [c[0] + c[1] * x + c[2] * y + c[3] * I.sq(x) + c[4] * I.sq(y) + c[5] * x * y,
                         c[1] + 2. * c[3] * x + c[5] * y,
                         c[2] + 2. * c[4] * y + c[5] * x,
                         c[5]]
            idx += 1
    for i in range(ncg):
        for j in range(nrg):
            ids = [i * nygrid + j, (i + 1) * nygrid + j, (i + 1) * nygrid + j + 1,
                   (1 + i * nygrid) + j]
            xr = (gx[i], gx[i + 1])
            yr = (gy[j], gy[j + 1])
            xsize = xr[1] - xr[0] + 1
            ysize = yr[1] - yr[0] + 1
            yy, xx = np.meshgrid(np.arange(ysize, dtype=np.float32),
                                 np.arange(xsize, dtype=np.float32), indexing="ij")
            fitimg[yr[0]:yr[1] + 1, xr[0]:xr[1] + 1] = bicuval(
                vals[ids, 0], vals[ids, 1], vals[ids, 2], vals[ids, 3],
                0.0, xsize, 0.0, ysize, xx, yy)
    return fitimg


def normspecflat(flat, edgecoeffs, xranges, slith_arc, nxgrid, nygrid,
                 ybuffer, ivar):
    """
    MC_NORMSPECFLAT: divide each order by a smooth model of itself.

    The model is MC_FITERPOLATE of a 5-pixel median of the order rectified
    at twice the slit sampling, interpolated back onto each column.
    Returns (normalized flat, normalized variance, rms per order, model).
    As in the IDL, the variance is divided by the model, not its square.
    """
    nrows, ncols = flat.shape
    y = np.arange(nrows, dtype=np.float32)
    nflat = np.ones((nrows, ncols), dtype=np.float32)
    ovar = np.full((nrows, ncols), np.nan, dtype=np.float32)
    model = np.ones((nrows, ncols), dtype=np.float32)
    rms = np.zeros(len(xranges), dtype=np.float32)
    for i in range(len(xranges)):
        m, yarc_map = rectorder(flat, edgecoeffs[i], xranges[i], slith_arc, 1, 2, ybuffer)
        m = fiterpolate(I.median_filter(m, 5), nxgrid, nygrid)
        start, stop = int(xranges[i][0]), int(xranges[i][1])
        x = np.arange(stop - start + 1, dtype=float) + start
        botedge = I.poly1d(x, edgecoeffs[i, 0])[0]
        topedge = I.poly1d(x, edgecoeffs[i, 1])[0]
        dif = topedge - botedge
        mask = np.zeros((nrows, ncols), dtype=bool)
        for j in range(stop - start + 1):
            if topedge[j] > nrows - 0.5 or botedge[j] < -0.5:
                continue
            # pixtoarc is FLTARR(2) in the IDL
            p1 = np.float32(np.float32(slith_arc) / np.float32(dif[j]))
            p0 = np.float32(np.float32(-1. * p1) * botedge[j])
            r0 = I.iround(botedge[j] + ybuffer)
            r1 = I.iround(topedge[j] - ybuffer)
            yy = y[r0:r1 + 1]
            nslit = I.linterp(yarc_map, m[:, j], I.poly(yy, np.array([p0, p1])))
            col = j + start
            nflat[r0:r1 + 1, col] = flat[r0:r1 + 1, col] / nslit
            ovar[r0:r1 + 1, col] = ivar[r0:r1 + 1, col] / nslit
            mask[r0:r1 + 1, col] = True
            model[r0:r1 + 1, col] = nslit
        rms[i] = I.moments(nflat[mask], robust=5)[2]
    nflat[nflat <= 0.0] = 1.0
    return nflat, ovar, rms, model


# ------------------------------------------------------------- the whole step

def make_flat(on_files, off_files=None, modeinfo=None, badpixmask="instrument",
              progress=None):
    """
    TS41CALS_MKFLAT. Returns a dict with the normalized flat and everything
    MC_WRITEFLAT records.

    on_files    the lamps-on frames (one or more; one frame is used as it
                is, neither scaled nor combined, as in the IDL)
    off_files   the lamps-off frames, or none (the IDL's blank field)
    badpixmask  "instrument": the mask ts4_soar.dat's BADPIXMASK names, as
                xspextool loads it (none at all for BADPIXMASK=None); an
                array; or None, no mask. Pixels where the mask is 0 are set
                to 1 in the normalized flat.

    result['findorders_qa'] holds the order-finding overlay for the combined
    flat (see ts4_findorders), result['combined_for_orders'] that image.
    """
    if not on_files:
        raise ValueError("no lamps-on flat frames given")
    if isinstance(badpixmask, str) and badpixmask == "instrument":
        from triplespec_reducer.instrument import read_badpixmask
        badpixmask = read_badpixmask()
    mi = modeinfo or read_modefile()
    say = progress or (lambda *_: None)
    say("Loading the flats")
    data, var, hdrs, _ = read_ts4(on_files, rotation=mi.rotation)
    n = data.shape[0]
    if off_files:
        say("Loading the darks")
        ddata, dvar, _, _ = read_ts4(off_files, rotation=mi.rotation)
        if ddata.shape[0] != 1:
            dark, darkvar, _ = I.meancomb(ddata)
        else:
            dark, darkvar = ddata[0], dvar[0]
        data = data - dark[None]
        var = var + darkvar[None]
    say("Finding the orders")
    med = I.meancomb(data)[0] if n != 1 else data[0]
    foqa = []
    edgecoeffs, xranges = ts4_findorders(
        med, mi.guesspos, mi.findrange, mi.step, mi.slith_pix,
        mi.slith_pix_range, mi.edgedeg, mi.norm_ybuffer, mi.comwin, qa=foqa)
    if n != 1:
        say("Scaling the flats")
        data, var = ka_scaleflats(data, var, edgecoeffs, xranges)
        say("Combining the flats")
        mean, var, _ = I.meancomb(data)
    else:
        mean, var = data[0], var[0]
    say("Normalizing the flat")
    flat, ovar, rms, model = normspecflat(mean, edgecoeffs, xranges, mi.slith_arc,
                                          mi.norm_nxg, mi.norm_nyg,
                                          mi.norm_ybuffer, var)
    if badpixmask is not None:
        flat[rotate(np.asarray(badpixmask), mi.rotation) == 0] = 1.0
    return dict(flat=flat, var=ovar, rms=rms, model=model, combined=mean,
                combined_for_orders=med, findorders_qa=foqa,
                edgecoeffs=edgecoeffs, xranges=xranges, orders=list(mi.orders),
                modeinfo=mi, header=hdrs[0], headers=hdrs,
                files=[str(f) for f in on_files],
                darks=[str(f) for f in (off_files or [])])


def flat_output_path(calpath, flatoname):
    """The file the IDL writes for Output Flat Name `flatoname`: cal/<name>.fits."""
    import os
    name = flatoname[:-5] if flatoname.endswith(".fits") else flatoname
    return os.path.join(calpath, name + ".fits")


def _hdrinfo(header, keywords):
    """
    MC_GETHDRINFO: each keyword's value as FXPAR types it, and its comment;
    a missing one becomes "Doesn't Exist" with no comment. Tag names lose
    their '-' (CREATE_STRUCT).
    """
    out = []
    for k in keywords:
        tag = "".join(k.split("-")).upper()
        if k in header:
            out.append((tag, I.fxpar(header, k), I.fxcomment(header, k)))
        else:
            out.append((tag, "Doesn't Exist", ""))
    return out


def write_flat(path, result, version=None, instrinfo=None):
    """
    MC_WRITEFLAT as TS41CALS_MKFLAT calls it: ext 0 the header, ext 1-3
    flat, variance, flags (unrotated).

    The header starts with the instrument file's KEYWORD list taken from
    the first lamps-on frame, its IRAFNAME (FILENAME) keyword set to the
    output name; VERSION is the Spextool version (data/Version.dat, '4.1')
    unless `version` is given; LINCRMAX is the instrument file's LINCORMAX.
    `instrinfo` is read_instrfile()'s structure (default: ts4_soar.dat).
    """
    import os
    from triplespec_reducer.instrument import read_instrfile, spextool_version
    info = instrinfo or read_instrfile()
    version = spextool_version() if version is None else version
    mi = result["modeinfo"]
    orders = result["orders"]
    ec = result["edgecoeffs"]
    hdr = fits.Header()
    src = result["header"]
    f32 = np.float32
    # MC_WRITEFLAT: the raw keywords as FXPAR typed them, then its own, every
    # number written by FXADDPAR/SXADDPAR in its IDL type (the FLOATs below).
    irafname = "".join(info.irafname.split("-")).upper()
    for tag, val, com in _hdrinfo(src, info.keywords):
        if tag == irafname:
            val = os.path.basename(path)
        I.fxaddpar(hdr, tag, val, com)
    I.fxaddpar(hdr, "MODENAME", mi.modename, " Spectroscopy mode")
    I.fxaddpar(hdr, "EDGEDEG", ec.shape[2] - 1, " Degree of the polynomial fit to order edges")
    I.fxaddpar(hdr, "NORDERS", len(orders), " Number of orders identified")
    I.fxaddpar(hdr, "ORDERS", ",".join(str(o) for o in orders), " Orders identified")
    I.fxaddpar(hdr, "PLTSCALE", f32(mi.ps), " Plate scale (arcseconds per pixel)")
    I.fxaddpar(hdr, "SLTH_PIX", mi.slith_pix, " Slit height in pixels")
    I.fxaddpar(hdr, "SLTH_ARC", f32(mi.slith_arc), " Slit height in arcseconds")
    I.fxaddpar(hdr, "SLTW_PIX", f32(SLITW_PIX), " Slit width in pixels")
    I.fxaddpar(hdr, "SLTW_ARC", f32(SLITW_ARC), " Slit width in arcseconds")
    I.fxaddpar(hdr, "RP", int(RESOLVINGPOWER), " Resolving power")
    I.fxaddpar(hdr, "ROTATION", mi.rotation, " IDL ROTATE value")
    I.fxaddpar(hdr, "VERSION", version, " Spextool version")
    I.fxaddpar(hdr, "LINCRMAX", int(info.lincormax), " Linearity maximum (DN)")    # LONG
    I.fxaddpar(hdr, "DARK", "None", " Dark frame")
    for i, o in enumerate(orders):
        I.fxaddpar(hdr, f"ODR{o:02d}RMS", f32(result["rms"][i]),
                   f" RMS deviation of the surface fit to order {o:02d}")
    for i, o in enumerate(orders):
        I.fxaddpar(hdr, f"ODR{o:02d}_XR", f"{result['xranges'][i][0]},{result['xranges'][i][1]}",
                   f" Extraction range for order {o:02d}")
    for i, o in enumerate(orders):
        for k in range(ec.shape[2]):
            I.fxaddpar(hdr, f"ODR{o:02d}_B{k + 1}", np.float64(ec[i, 0, k]),
                       f" a{k} edge coefficient for bottom of order {o:02d}")
        for k in range(ec.shape[2]):
            I.fxaddpar(hdr, f"ODR{o:02d}_T{k + 1}", np.float64(ec[i, 1, k]),
                       f" a{k} edge coefficient for top of order {o:02d}")
    I.addhistlabel(hdr, "Xspextool History")
    names = [str(I.fxpar(h, info.irafname)) if info.irafname in h else "Doesn't Exist"
             for h in result.get("headers", [result["header"]])]
    I.sxaddhist(hdr, "This flat was created by 1. scaling " + ", ".join(names)
                + " using the median and then 2. Robust Mean combining the resulting images.")
    rot = mi.rotation
    hdul = fits.HDUList([fits.PrimaryHDU(header=hdr),
                         fits.ImageHDU(unrotate(result["flat"], rot).astype(np.float32)),
                         fits.ImageHDU(unrotate(result["var"], rot).astype(np.float32)),
                         fits.ImageHDU(np.zeros_like(unrotate(result["flat"], rot), dtype=np.uint8))])
    hdul.writeto(path, overwrite=True)


def read_flat(path):
    """MC_READFLAT: the flat, rotated back, and its order geometry."""
    with fits.open(path) as hdul:
        hdr = hdul[0].header
        rot = int(hdr["ROTATION"])
        image = rotate(hdul[1].data.astype(np.float32), rot)
        variance = rotate(hdul[2].data.astype(np.float32), rot)
    # MC_READFLAT reads every number back with FXPAR, which types it from its text
    norders = I.fxpar(hdr, "NORDERS")
    edgedeg = I.fxpar(hdr, "EDGEDEG")
    orders = [int(t) for t in str(hdr["ORDERS"]).split(",")]
    ec = np.zeros((norders, 2, edgedeg + 1))
    xr = np.zeros((norders, 2), dtype=int)
    rms = np.zeros(norders, dtype=np.float32)
    for i, o in enumerate(orders):
        ec[i, 1] = [I.fxpar(hdr, f"ODR{o:02d}_T{k + 1}") for k in range(edgedeg + 1)]
        ec[i, 0] = [I.fxpar(hdr, f"ODR{o:02d}_B{k + 1}") for k in range(edgedeg + 1)]
        xr[i] = [int(t) for t in str(hdr[f"ODR{o:02d}_XR"]).split(",")]
        rms[i] = I.fxpar(hdr, f"ODR{o:02d}RMS") if f"ODR{o:02d}RMS" in hdr else np.nan
    nrows, ncols = image.shape
    return dict(flat=image, var=variance, header=hdr, orders=orders,
                edgecoeffs=ec, xranges=xr, rms=rms, rotation=rot,
                ps=I.fxpar(hdr, "PLTSCALE"), slith_arc=I.fxpar(hdr, "SLTH_ARC"),
                slith_pix=I.fxpar(hdr, "SLTH_PIX"), slitw_arc=I.fxpar(hdr, "SLTW_ARC"),
                slitw_pix=I.fxpar(hdr, "SLTW_PIX"), rp=I.fxpar(hdr, "RP"),
                modename=str(hdr["MODENAME"]).strip(),
                omask=order_mask(ncols, nrows, orders, ec, xr))
