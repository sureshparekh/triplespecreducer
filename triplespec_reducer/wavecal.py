"""
Wavelength calibration on sky (OH) lines or the CuHeAr lamp: TS41CALS_MKWAVECAL.

    frames read                     KNA_READTS4FITS
    median sky / flat               (the 2D image the lines are traced on)
    a pixel-only "wavecal"          MC_SIMWAVECAL2D
    a 1D spectrum of each order     KNA_MEDEXTSPEC1D  (+-5" about slit centre)
    offset from the template        C_CORRELATE + GAUSSFIT on the home order
    1D solution over all orders     MC_XDWAVECAL1D   lambda*m/4 = P(x, m)
    each line traced along the slit MC_FINDLINES2D_KNAHACK
    2D distortion and the images    MC_WAVECAL2D_KNAHACK

Two things are reproduced from the IDL as it is, not as it may have been
meant. The 1D spectra are taken from the first frame only, unflattened:
KNA_MEDEXTSPEC1D is handed the whole image cube and indexes it with two
subscripts, which in IDL addresses the first plane. And the aperture mask
there is an integer array, so its fractional edge weights truncate to 0.
"""

import os

import numpy as np
from astropy.io import fits

from triplespec_reducer import idl as I
from triplespec_reducer.frames import read_ts4, rotate, unrotate
from triplespec_reducer.instrument import read_modefile, INSTRUMENT_DIR


# ------------------------------------------------------------- inputs

def read_linecal(path):
    """MC_READLINECALFILE."""
    with fits.open(path) as hdul:
        spec = hdul[0].data.astype(np.float32)     # [norders, 3, npix], FLOAT as stored
        hdr = hdul[0].header
    from triplespec_reducer.instrument import fsextract_index
    orders = [int(t) for t in str(hdr["ORDERS"]).split(",")]
    w2pdeg = [int(t) for t in str(hdr["W2PDEG"]).split(",")]
    w2p = np.zeros((len(orders), max(w2pdeg) + 1))
    for i, o in enumerate(orders):
        w2p[i, :w2pdeg[i] + 1] = [hdr[f"W2P{o:02d}_A{k}"] for k in range(w2pdeg[i] + 1)]
    xcororder = int(hdr["XCORORDR"])
    z = orders.index(xcororder)
    return dict(xcorspec=spec[z], xcororder=xcororder,
                ordertype=[t.strip() for t in str(hdr["ORDRTYPE"]).split(",")],
                w2pcoeffs=w2p, homeorder=int(hdr["HOMEORDR"]),
                dispdeg=int(hdr["DISPDEG"]), ordrdeg=int(hdr["ORDRDEG"]),
                orders=orders, spectra=spec)


def read_linelist(path):
    """READCOL,file,lines,FORMAT='D': the first number of every line that has one."""
    out = []
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if not t:
                continue
            try:
                out.append(float(t[0]))
            except ValueError:
                pass
    return np.array(out)


def rdohll(path, res=None, thresh=1.0):
    """MC_RDOHLL: drop lines closer than thresh*lambda/R/2.534 to another."""
    lines = read_linelist(path)
    if res is None:
        return lines
    n = lines.size
    gb = np.ones(n, dtype=int)
    for i in range(1, n - 1):
        delta = lines[i] / float(res)
        z = (lines >= lines[i] - thresh * delta / 2.534) & \
            (lines <= lines[i] + thresh * delta / 2.534)
        if z.sum() > 1:
            gb[z] = 0
    return lines[gb == 1]


# ------------------------------------------------------------- 1D

def simwavecal2d(ncols, nrows, edgecoeffs, xranges, slith_arc):
    """MC_SIMWAVECAL2D: column number as 'wavelength', arcsec along the slit."""
    wavecal = np.full((nrows, ncols), np.nan, dtype=np.float32)
    spatcal = np.full((nrows, ncols), np.nan, dtype=np.float32)
    y = np.arange(nrows, dtype=float)
    for i in range(len(xranges)):
        start, stop = int(xranges[i][0]), int(xranges[i][1])
        x = np.arange(stop - start + 1, dtype=float) + start
        bot = I.poly1d(x, edgecoeffs[i, 0])[0]
        top = I.poly1d(x, edgecoeffs[i, 1])[0]
        p1 = slith_arc / (top - bot)
        p0 = -1. * (p1 * bot)
        for j in range(stop - start + 1):
            r0, r1 = I.itrunc(bot[j]), I.itrunc(top[j])
            wavecal[r0:r1 + 1, int(x[j])] = x[j]
            spatcal[r0:r1 + 1, int(x[j])] = p0[j] + p1[j] * y[r0:r1 + 1]
    return wavecal, spatcal


def medextspec1d(img, var, omask, orders, wavecal, spatcal, apsign, apranges):
    """
    KNA_MEDEXTSPEC1D with APRANGES: the median of each column's aperture.

    Returns one [3, nx] array (wave, flux, variance) per order.
    """
    nrows, ncols = img.shape
    out = []
    for i, o in enumerate(orders):
        cols = np.where((omask == o).any(axis=0))[0]
        x0, x1 = cols.min(), cols.max()
        nx = x1 - x0 + 1
        ow = np.full(nx, np.nan, dtype=np.float32)
        of = np.full(nx, np.nan, dtype=np.float32)
        ov = np.full(nx, np.nan, dtype=np.float32)
        for k in range(x0, x1 + 1):
            z = np.where(omask[:, k] == o)[0]
            cnt = z.size
            slitimg = img[z, k]
            slitspat = spatcal[z, k].astype(float)
            slitwave = wavecal[z, k]
            apmask = np.zeros(cnt, dtype=int)
            aprange = [apranges[i][0], apranges[i][1]]
            apidx = I.tabinv(slitspat, aprange)
            a0, a1 = I.itrunc(apidx[0]), I.itrunc(apidx[1])
            apmask[a0:a1 + 1] = 1
            # Integer mask: every fractional weight below truncates to 0.
            if apidx[0] - np.floor(apidx[0]) >= 0.5:
                apmask[a0] = 0
                apmask[I.itrunc(apidx[0] + 1)] = I.itrunc(0.5 + I.iround(apidx[0]) - apidx[0])
            else:
                apmask[a0] = I.itrunc(0.5 - (apidx[0] - np.floor(apidx[0])))
            if apidx[1] - np.floor(apidx[1]) >= 0.5:
                apmask[I.itrunc(apidx[1] + 1)] = I.itrunc(0.5 - (I.iround(apidx[1]) - apidx[1]))
            else:
                apmask[a1] = I.itrunc(0.5 + (apidx[1] - I.iround(apidx[1])))
            zz = apmask > 0
            ow[k - x0] = slitwave[cnt // 2]
            of[k - x0] = I.median(slitimg[zz])
            sel = slitimg[zz].astype(float)
            ov[k - x0] = I.sq((np.std(sel, ddof=1) / float(sel.size))) if sel.size > 1 else np.nan
        out.append(np.array([ow, float(apsign[0]) * of, ov]))
    return out


def xcorr_offset(calinfo, spec_home, ncols):
    """The pixel offset of the home order from the template, as autoxcorr does."""
    whome = calinfo["xcorspec"][0]
    fhome = calinfo["xcorspec"][1]
    idx = nantrim(whome, 2)
    whome, fhome = whome[idx], fhome[idx]
    wobs, fobs = spec_home[0], spec_home[1]
    ifobs = I.linterp(wobs, fobs, whome, missing=0)
    lag = np.arange(ncols, dtype=np.float32) - I.itrunc(ncols / 2.)
    corr = I.c_correlate(fhome, ifobs, lag).astype(np.float32)
    k = int(np.nanargmax(corr))
    _, a = I.gaussfit(lag[k - 20:k + 21], corr[k - 20:k + 21], nterms=4)
    return float(a[1]), (lag, corr)


def nantrim(array, flag):
    """MC_NANTRIM: indices that drop leading (1), trailing (0), both (2) NaNs, or all (3)."""
    a = np.asarray(array, dtype=float)
    n = a.size
    good = np.isfinite(a)
    if good.all():
        return np.arange(n)
    if flag == 3:
        return np.where(good)[0]
    nz = np.where(np.cumsum(np.nan_to_num(a)) != 0)[0]
    rz = np.where(np.cumsum(np.nan_to_num(a[::-1])) != 0)[0]
    lo = nz[0] if flag in (1, 2) else 0
    hi = (n - 1 - rz[0]) if flag in (0, 2) else n - 1
    return np.arange(lo, hi + 1)


def findlines1d(x, y, xguess, sigma):
    """
    MC_FINDLINES1D: a Gaussian on a constant at every guessed line.

    The start is FLOAT([...]) and the spectrum is single precision, so
    MPFITEXPR fits in float; mpfit.fit_gauss_const does the same.
    """
    from triplespec_reducer.mpfit import fit_gauss_const
    x = np.asarray(x, dtype=np.float32)
    y = np.asarray(y, dtype=np.float32)
    nx = x.size
    hwin = I.itrunc(np.float32(sigma) * np.float32(5.))
    goodbad = np.zeros(len(xguess), dtype=int)
    lines = np.zeros(len(xguess), dtype=np.float32)          # FLTARR
    for i, xg in enumerate(xguess):
        idx = I.tabinv(x, xg)
        lo = max(0, I.iround(idx) - hwin)
        hi = min(nx - 1, I.iround(idx) + hwin)
        xwin = x[lo:hi + 1]
        ywin = y[lo:hi + 1]
        if not np.all(np.isfinite(ywin)):
            continue
        idx = I.tabinv(xwin, xg)
        start = np.array([I.median(ywin), ywin[I.iround(idx)], xg, sigma], dtype=np.float32)
        fit, _ = fit_gauss_const(xwin, ywin, start)
        lines[i] = fit[2]
        if fit[2] <= xg + 1 and fit[2] >= xg - 1 and fit[3] > 0.0 and fit[1] > 0:
            goodbad[i] = 1
    return lines, goodbad


def xdwavecal1d(spectra, orders, w2pcoeffs, offset, lines, ordertypes,
                linetypes, slitw_pix, homeorder, dispdeg, ordrdeg, robustthresh):
    """
    MC_XDWAVECAL1D: one polynomial in column and order for every order.

    lambda * m / homeorder = P(x, m). Returns the coefficients, the ranges,
    lines found per order, the rms in Angstrom, and the fitted lines.
    """
    nlines = []
    rows = []
    for i, o in enumerate(orders):
        wspec = spectra[i][0]
        fspec = spectra[i][1]
        idx = nantrim(wspec, 2)
        wspec, fspec = wspec[idx], fspec[idx]
        goodl = np.array([t == ordertypes[i].strip() for t in linetypes])
        plines = lines[goodl]
        xguess = I.poly(plines, w2pcoeffs[i]) + offset
        z = (xguess > wspec.min()) & (xguess < wspec.max())
        xguess, plines = xguess[z], plines[z]
        slines = plines * o / homeorder
        xpos, gb = findlines1d(wspec, fspec, xguess, slitw_pix / I.FL(2.354))
        g = gb == 1
        nlines.append(int(g.sum()))
        for xp, sl in zip(xpos[g], slines[g]):
            rows.append((xp, o, sl))
    # SARRAY is DOUBLE (the scaled wavelengths are), holding FLOAT centroids.
    a = np.array(rows, dtype=float)
    s = np.argsort(a[:, 2], kind="stable")
    sx, so, sl = a[s, 0], a[s, 1], a[s, 2]
    coeffs, ogb, rms = I.robustpoly2d(sx, so, sl, dispdeg, ordrdeg, robustthresh, 0.01,
                                      want_rms=True)
    rms = rms * 10000.
    wranges = np.zeros((len(orders), 2), dtype=np.float32)      # FLTARR
    for i, o in enumerate(orders):
        w = spectra[i][0].astype(float)
        xmin, xmax = np.min(w), np.max(w)
        wranges[i] = I.poly2d(np.array([xmin, xmax]), np.array([o, o], float),
                              dispdeg, ordrdeg, coeffs) / o * homeorder
    return dict(coeffs=coeffs, wranges=wranges, nlines=nlines, rms=rms,
                fitlines=dict(x=sx, order=so, slines=sl, goodbad=ogb))


# ------------------------------------------------------------- 2D

_FITLOG = None        # debugging: set to a list to record every line fit


def findlines2d(img, edgecoeffs, xranges, orders, wlines, w2pcoeffs, offset,
                wranges, slitw_pix, slith_arc, linedeg, linewin, ystep, ysum):
    """
    MC_FINDLINES2D_KNAHACK: follow every line up and down the slit.

    For each line, Gaussians are fitted to 5-row medians every `ystep` rows
    from mid-slit towards the top, then the bottom; the centroids give the
    line's column at each row. Each point is labelled with the line's
    column at mid-slit (sl) and its arc-length position along the line,
    as a fraction of the slit (ss). Returns one dict per order.
    """
    nrows, ncols = img.shape
    x = np.arange(ncols, dtype=np.float32)          # FINDGEN: xtop, xbot, xmid come out FLOAT
    out = []
    for i, o in enumerate(orders):
        z = (wlines >= wranges[i][0]) & (wlines <= wranges[i][1])
        wline = wlines[z]
        xgline = I.poly(wline, w2pcoeffs[i]) + offset
        top = I.poly(x, edgecoeffs[i, 1])
        bot = I.poly(x, edgecoeffs[i, 0])
        mid = bot + (top - bot) / 2.
        idx = I.tabinv(x, xgline)
        acc = dict(sx=[], sy=[], ss=[], sl=[], sw=[])
        for j in range(np.size(idx)):
            guessx = x[I.itrunc(idx[j])]
            if guessx < 4 or guessx > ncols - 5:
                continue
            midslit = mid[I.itrunc(idx[j])]
            tx, ty = [], []
            skip = False

            def fit_at(gx, row):
                lo = I.itrunc(max(gx - linewin, 0))
                hi = I.itrunc(min(gx + linewin, ncols - 1))
                xsub = x[lo:hi + 1]
                r0 = I.itrunc(row - ysum)
                r1 = I.itrunc(row + ysum)
                ysub = I.median(img[r0:r1 + 1, lo:hi + 1], axis=0)
                zidx = I.tabinv(xsub, gx)
                guessz = ysub[I.itrunc(zidx)]
                est = [guessz, gx, slitw_pix / I.FL(2.354), np.nanmin(ysub)]
                _, a = I.mpfitpeak(xsub.astype(float), ysub.astype(float), nterms=4,
                                   estimates=est)
                if _FITLOG is not None:
                    _FITLOG.append((o, gx, row, xsub.copy(), ysub.copy(), np.array(est, float), a.copy()))
                return a

            for k in range(nrows):
                topidx = I.tabinv(x, guessx)
                if midslit + k * ystep + 2 >= top[I.itrunc(topidx)]:
                    break
                if k == 5 and len(tx) + 1 < 3:
                    skip = True
                    break
                a = fit_at(guessx, midslit + k * ystep)
                if a[0] < 0.0 or abs(a[1] - guessx) > 2:
                    continue
                tx.append(a[1])
                ty.append(midslit + k * ystep)
                if len(tx) + 1 > 3:
                    c, _ = I.robustpoly1d(np.array(ty), np.array(tx),
                                          min(len(tx), linedeg), 4, 0.1)
                    guessx = I.poly(midslit + k * ystep, c)
                else:
                    guessx = a[1]
            if skip:
                continue
            if len(tx) < 2:
                # IDL's POLY_FIT stops the whole calibration here; a line
                # with fewer than two points is left out instead.
                continue
            for k in range(1, nrows):
                c = I.poly_fit(np.array(ty), np.array(tx), linedeg)
                guessx = I.poly(midslit - k * ystep, c)
                if guessx < 4 or guessx > ncols - 5:
                    continue
                botidx = I.tabinv(x, guessx)
                if midslit - k * ystep - 2 <= bot[I.itrunc(botidx)]:
                    break
                a = fit_at(guessx, midslit - k * ystep)
                if a[0] < 0.0 or abs(a[1] - guessx) > 2:
                    continue
                tx.append(a[1])
                ty.append(midslit - k * ystep)
            tx = np.array(tx)
            ty = np.array(ty)
            s = np.argsort(tx, kind="stable")
            tx, ty = tx[s], ty[s]
            lc, ogb = I.robustpoly1d(ty, tx, linedeg, 3.5, 0.01)
            xtop = I.linterp(x - I.poly(top, lc), x, 0.0)
            ytop = I.poly(xtop, edgecoeffs[i, 1])
            xbot = I.linterp(x - I.poly(bot, lc), x, 0.0)
            ybot = I.poly(xbot, edgecoeffs[i, 0])
            xmid = I.linterp(x - I.poly(mid, lc), x, 0.0)
            deriv = np.array([k * lc[k] for k in range(1, lc.size)])
            func = lambda yy: np.sqrt(1.0 + I.sq(I.poly(yy, deriv)))
            linelength = I.qsimp(func, ybot, ytop)
            g = np.where(ogb == 1)[0]
            tss = np.array([I.qsimp(func, ybot, ty[kk]) / linelength * slith_arc for kk in g])
            ss_order = np.argsort(tss, kind="stable")
            acc["sx"] += list(tx[g][ss_order])
            acc["sy"] += list(ty[g][ss_order])
            acc["sl"] += [xmid] * g.size
            acc["ss"] += list(tss[ss_order])
            # The IDL labels these with wlines[j], indexing the full list
            # with a position in this order's sublist; the label is only
            # used in a printed report, and is kept as the IDL writes it.
            acc["sw"] += [wlines[j]] * g.size
        out.append({k: np.array(v, dtype=float) for k, v in acc.items()})
    return out


def wavecal2d(wcinfo1d, wcinfo2d, edgecoeffs, omask, orders, xranges, wranges,
              wdx, wdy, sdx, sdy):
    """MC_WAVECAL2D_KNAHACK: the wavelength and slit-position images."""
    nrows, ncols = omask.shape
    warr = np.full((nrows, ncols), np.nan)
    sarr = np.full((nrows, ncols), np.nan, dtype=np.float32)
    ximg, yimg = np.meshgrid(np.arange(ncols, dtype=float), np.arange(nrows, dtype=float))
    wcoeffs, scoeffs, disp, stats = [], [], [], []
    for i, o in enumerate(orders):
        L = wcinfo2d[i]
        oyimg = yimg - I.poly(np.arange(ncols, dtype=float), edgecoeffs[i, 0])[None, :]
        sx = L["sx"]
        sy = L["sy"] - I.poly(sx, edgecoeffs[i, 0])
        sc, sgb = I.robustpoly2d(sx, sy, L["ss"], sdx, sdy, 4, 0.01)
        wc, wgb = I.robustpoly2d(sx, sy, L["sl"], wdx, wdy, 4, 0.01)
        z = omask == o
        xcen = I.poly2d(ximg[z], oyimg[z], wdx, wdy, wc)
        warr[z] = I.poly2d(xcen, np.full(xcen.size, float(o)), wcinfo1d["dispdeg"],
                           wcinfo1d["ordrdeg"], wcinfo1d["coeffs"]) / float(o) * wcinfo1d["homeorder"]
        sarr[z] = I.poly2d(ximg[z], oyimg[z], sdx, sdy, sc)
        wcoeffs.append(wc)
        scoeffs.append(sc)
        # FLOAT / LONG: a single-precision division
        disp.append(np.float32(wranges[i][1] - wranges[i][0])
                    / np.float32(xranges[i][1] - xranges[i][0]))
        sres = L["ss"] - I.poly2d(sx, sy, sdx, sdy, sc)
        lres = L["sl"] - I.poly2d(sx, sy, wdx, wdy, wc)
        stats.append(dict(order=o, nlines=int(np.unique(L["sl"]).size),
                          spatial_rms=float(np.std(sres[sgb == 1], ddof=1)),
                          wave_rms_pix=float(np.std(lres[wgb == 1], ddof=1)),
                          sx=sx, sy=L["sy"], sres=sres, lres=lres, sgb=sgb, wgb=wgb))
    return warr, sarr, wcoeffs, scoeffs, disp, stats


# ------------------------------------------------------------- the whole step

def make_wavecal(files, flatinfo, arctype="Sky", modeinfo=None, progress=None, offset=None):
    """
    TS41CALS_MKWAVECAL. `flatinfo` is read_flat()'s dict. `offset` replaces
    the automatic cross-correlation (xspextool's 'Auto X-Correlate' off,
    where XGETOFFSET has the user mark one line in both spectra).
    """
    mi = modeinfo or read_modefile()
    say = progress or (lambda *_: None)
    rot = flatinfo["rotation"]
    orders = flatinfo["orders"]
    ec, xr = flatinfo["edgecoeffs"], flatinfo["xranges"]
    nrows, ncols = flatinfo["flat"].shape
    say("Loading the sky frames" if arctype == "Sky" else "Loading the arc frames")
    data, var, hdrs, _ = read_ts4(files, rotation=rot)
    arc = (I.median(data, even=True, axis=0) if data.shape[0] != 1 else data[0])
    arc = (arc / flatinfo["flat"]).astype(np.float32)
    say("Wavelength calibrating")
    wfake, sfake = simwavecal2d(ncols, nrows, ec, xr, flatinfo["slith_arc"])
    sa = flatinfo["slith_arc"]
    spec = medextspec1d(data[0], var[0], flatinfo["omask"], orders, wfake, sfake,
                        [1], [[sa / 2. - 5, sa / 2 + 5]] * len(orders))
    calname = "CuHeArXD_LineCal.fits" if arctype != "Sky" else "OHXD_LineCal.fits"
    calinfo = read_linecal(os.path.join(INSTRUMENT_DIR, calname))
    zh = orders.index(calinfo["xcororder"])
    if offset is None:
        offset, xc = xcorr_offset(calinfo, spec[zh], ncols)
    else:
        offset, xc = float(offset), None
    llname = "Arll.dat" if arctype != "Sky" else "OHll.dat"
    lines = read_linelist(os.path.join(INSTRUMENT_DIR, llname))
    lt = "Ar" if arctype != "Sky" else "OH"
    w1 = xdwavecal1d(spec, orders, calinfo["w2pcoeffs"], offset, lines,
                     calinfo["ordertype"], [lt] * lines.size, flatinfo["slitw_pix"],
                     calinfo["homeorder"], mi.dispdeg, mi.ordrdeg, 3.5)
    say("Tracing the lines along the slit")
    plines = rdohll(os.path.join(INSTRUMENT_DIR, llname), res=flatinfo["rp"], thresh=3)
    w2 = findlines2d(arc, ec, xr, orders, plines, calinfo["w2pcoeffs"], offset,
                     w1["wranges"], flatinfo["slitw_pix"], sa, mi.linedeg,
                     mi.linereg / 2., mi.ystep, mi.ysum)
    say("Fitting the 2D solution")
    wcinfo1d = dict(coeffs=w1["coeffs"], dispdeg=mi.dispdeg, ordrdeg=mi.ordrdeg,
                    homeorder=mi.homeorder)
    warr, sarr, wc, sc, disp, stats = wavecal2d(
        wcinfo1d, w2, ec, flatinfo["omask"], orders, xr, w1["wranges"],
        mi.wxdeg, mi.wydeg, mi.sxdeg, mi.sydeg)
    return dict(wavecal=warr, spatcal=sarr, wcinfo1d=wcinfo1d, wcoeffs=wc,
                scoeffs=sc, disp=disp, orders=orders, rotation=rot, offset=offset,
                xcorr=xc, spec1d=spec, sol1d=w1, lines2d=w2, stats=stats, arc=arc,
                modeinfo=mi, files=[str(f) for f in files], arctype=arctype)


def write_wavecal(path, r):
    """The file MC_WAVECAL2D_KNAHACK writes, read back by MC_RDWAVECAL2D."""
    mi = r["modeinfo"]
    w1 = r["wcinfo1d"]
    hdr = fits.Header()
    # numbers as FXADDPAR/SXADDPAR write them; in MC_WAVECAL2D_KNAHACK's order
    I.fxaddpar(hdr, "WCTYPE", "2D", " Wavelength calibration type")
    I.fxaddpar(hdr, "WDEG", int(w1["dispdeg"]), " X Polynomial degree for 1D wavelength calibration")
    I.fxaddpar(hdr, "ODEG", int(w1["ordrdeg"]), " Y Polynomial degree for 1D wavelength calibration")
    I.fxaddpar(hdr, "WXDEG", int(mi.wxdeg), ' X Polynomial degree for 2D "wavelength" distortion calibration')
    I.fxaddpar(hdr, "WYDEG", int(mi.wydeg), ' Y Polynomial degree for 2D "wavelength" distortion calibration')
    I.fxaddpar(hdr, "SXDEG", int(mi.sxdeg), ' X Polynomial degree for 2D "spatial" distortion calibration')
    I.fxaddpar(hdr, "SYDEG", int(mi.sydeg), ' Y Polynomial degree for 2D "spatial" distortion calibration')
    I.fxaddpar(hdr, "NORDERS", len(r["orders"]), " Number of orders")
    I.fxaddpar(hdr, "ORDERS", ",".join(str(o) for o in r["orders"]), " Order numbers")
    I.fxaddpar(hdr, "HOMEORDR", f"{int(w1['homeorder']):02d}", " Home order for 1D wavelength calibration")
    for j, c in enumerate(w1["coeffs"]):
        I.fxaddpar(hdr, f"1DWC{j + 1:02d}", np.float64(c), f" a{j + 1:02d} 1D wavecal coefficient")
    for i, o in enumerate(r["orders"]):
        for j, c in enumerate(r["wcoeffs"][i]):
            I.fxaddpar(hdr, f"OR{o:02d}WC{j + 1:02d}", np.float64(c),
                       f' a{j + 1:02d} 2D "wavecal" distortion coefficient for order {o:02d}')
        for j, c in enumerate(r["scoeffs"][i]):
            I.fxaddpar(hdr, f"OR{o:02d}SC{j + 1:02d}", np.float64(c),
                       f' a{j + 1:02d} 2D "spatial" distortion coefficient for order {o:02d}')
        # FLOAT: the FLOAT wavelength range over the LONG column range
        I.fxaddpar(hdr, f"DISPO{o:02d}", np.float32(r["disp"][i]),
                   f" Dispersion (microns pixel-1) for order {o:02d}")
    hdr["TSRTYPE"] = (r["arctype"], " Sky (OH) or CuHeAr arc")
    hdr["TSROFFS"] = (float(r["offset"]), " Cross-correlation offset (pixels)")
    hdr["TSRRMS1D"] = (float(r["sol1d"]["rms"]), " RMS of the 1D solution (Angstrom)")
    for f in r["files"]:
        hdr.add_history("wavecal image: " + os.path.basename(f))
    rot = r["rotation"]
    fits.HDUList([fits.PrimaryHDU(header=hdr),
                  fits.ImageHDU(unrotate(r["wavecal"], rot)),
                  fits.ImageHDU(unrotate(r["spatcal"], rot).astype(np.float32))]
                 ).writeto(path, overwrite=True)


def read_wavecal(path, rotation=3):
    """MC_RDWAVECAL2D."""
    with fits.open(path) as hdul:
        hdr = hdul[0].header
        # MRDFITS keeps the stored types: DOUBLE wavelengths, FLOAT slit positions.
        d1, d2 = hdul[1].data, hdul[2].data
        wavecal = rotate(d1.astype(d1.dtype.newbyteorder("=")), rotation)
        spatcal = rotate(d2.astype(d2.dtype.newbyteorder("=")), rotation)
    orders = [int(t) for t in str(hdr["ORDERS"]).split(",")]
    # MC_FXPAR(hdr, 'DISP*'): each typed by FXPAR, concatenated
    vals = [I.fxpar(hdr, k) for k in hdr if k.startswith("DISP")]
    disp = np.array(vals, dtype=np.result_type(*vals))
    return dict(wavecal=wavecal, spatcal=spatcal, header=hdr, wctype=str(hdr["WCTYPE"]).strip(),
                orders=orders, disp=disp)
