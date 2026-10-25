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
    """
    MC_RDOHLL: drop lines closer than thresh*lambda/R/2.534 to another.

    2.534 is a FLOAT literal in a DOUBLE expression, so it enters as the
    single-precision value promoted (I.FL).
    """
    lines = read_linelist(path)
    if res is None:
        return lines
    n = lines.size
    gb = np.ones(n, dtype=int)
    for i in range(1, n - 1):
        delta = lines[i] / float(res)
        z = (lines >= lines[i] - thresh * delta / I.FL(2.534)) & \
            (lines <= lines[i] + thresh * delta / I.FL(2.534))
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
        # PIXTOARC is FLTARR: both coefficients are stored in single precision
        p1 = (np.float64(np.float32(slith_arc)) / (top - bot)).astype(np.float32)
        p0 = (-1. * (p1.astype(float) * bot)).astype(np.float32).astype(float)
        p1 = p1.astype(float)
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
    # The <oname>_Wavecal.eps panels: residuals in Angstrom, their MOMENT
    # over the lines kept, and the lines the fit rejected.
    resid = (sl - I.poly2d(sx, so, dispdeg, ordrdeg, coeffs)) * 10000.
    m = I.moments(resid[ogb == 1])
    qa = dict(order=so, column=sx, resid=resid, goodbad=ogb, mean=float(m[0]),
              var=float(m[1]), title="RMS = %5.3f Angstroms" % np.sqrt(m[1]))
    return dict(coeffs=coeffs, wranges=wranges, nlines=nlines, rms=rms,
                fitlines=dict(x=sx, order=so, slines=sl, goodbad=ogb), qa=qa)


# ------------------------------------------------------------- 2D

_FITLOG = None        # debugging: set to a list to record every line fit


class ShortLineError(RuntimeError):
    """
    A line was traced at too few rows for the fit that starts the downward
    walk. MC_FINDLINES2D_KNAHACK's POLY_FIT stops there (a subscript error
    with no point, 'Singular matrix detected' with fewer than LINEDEG+1), so
    the IDL wavelength calibration ends without writing anything.
    """


def findlines2d(img, edgecoeffs, xranges, orders, wlines, w2pcoeffs, offset,
                wranges, slitw_pix, slith_arc, linedeg, linewin, ystep, ysum,
                short_lines="abort"):
    """
    MC_FINDLINES2D_KNAHACK: follow every line up and down the slit.

    For each line, Gaussians are fitted to 5-row medians every `ystep` rows
    from mid-slit towards the top, then the bottom; the centroids give the
    line's column at each row. Each point is labelled with the line's
    column at mid-slit (sl) and its arc-length position along the line,
    as a fraction of the slit (ss). Returns one dict per order.

    A line that reaches the top of the slit with fewer than LINEDEG+1
    points stops the IDL calibration (POLY_FIT fails). `short_lines` is
    "abort" (the IDL behaviour: ShortLineError) or "skip" (not IDL: the
    line is left out and the calibration goes on).

    Besides sx, sy, ss, sl, sw, each dict holds what the IDL draws on the
    arc while tracing (WID=): wid_x, wid_y, every accepted centroid,
    including those of lines dropped later; and, not drawn by the IDL,
    rej_x, rej_y, the centroids the line fit (ROBUSTPOLY1D) rejected, and
    skipped, the lines left out with short_lines="skip".
    """
    if short_lines not in ("abort", "skip"):
        raise ValueError("short_lines must be 'abort' (IDL) or 'skip'")
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
        wid_x, wid_y, rej_x, rej_y, skipped = [], [], [], [], []
        for j in range(np.size(idx)):
            guessx = x[I.itrunc(idx[j])]
            if guessx < 4 or guessx > ncols - 5:
                continue
            guessx0 = guessx
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
                wid_x.append(a[1])
                wid_y.append(midslit + k * ystep)
                if len(tx) + 1 > 3:
                    c, _ = I.robustpoly1d(np.array(ty), np.array(tx),
                                          min(len(tx), linedeg), 4, 0.1)
                    guessx = I.poly(midslit + k * ystep, c)
                else:
                    guessx = a[1]
            if skip:
                continue
            if len(tx) < linedeg + 1:
                if short_lines == "skip":
                    skipped.append(dict(wline=float(wline[j]), column=float(guessx0),
                                        npoints=len(tx)))
                    continue
                raise ShortLineError(
                    f"Order {o}: the line at {wline[j]:.6f} um (column {guessx0:.1f}) was "
                    f"found at only {len(tx)} row(s) before the top of the slit; "
                    f"{linedeg + 1} are needed to start tracing it downwards. "
                    "Spextool's MC_FINDLINES2D_KNAHACK stops the whole wavelength "
                    "calibration here (POLY_FIT), writing nothing. To leave such lines "
                    "out instead (not what Spextool does), use short_lines='skip'.")
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
                wid_x.append(a[1])
                wid_y.append(midslit - k * ystep)
            tx = np.array(tx)
            ty = np.array(ty)
            s = np.argsort(tx, kind="stable")
            tx, ty = tx[s], ty[s]
            lc, ogb = I.robustpoly1d(ty, tx, linedeg, 3.5, 0.01)
            rej_x += list(tx[ogb == 0])
            rej_y += list(ty[ogb == 0])
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
        d = {k: np.array(v, dtype=float) for k, v in acc.items()}
        d.update(wid_x=np.array(wid_x, float), wid_y=np.array(wid_y, float),
                 rej_x=np.array(rej_x, float), rej_y=np.array(rej_y, float),
                 skipped=skipped)
        out.append(d)
    return out


def _idl_uniq(a):
    """UNIQ on an unsorted array: the last index of every run of equal values."""
    a = np.asarray(a)
    if a.size == 0:
        return np.zeros(0, dtype=int)
    return np.where(np.r_[a[1:] != a[:-1], True])[0]


def _strtrim_double(v):
    """STRTRIM(double, 2): IDL's default G16.8 text."""
    return f"{float(v):#.8g}"


def wavecal2d(wcinfo1d, wcinfo2d, edgecoeffs, omask, orders, xranges, wranges,
              wdx, wdy, sdx, sdy, ps=None):
    """
    MC_WAVECAL2D_KNAHACK: the wavelength and slit-position images.

    Each stats entry also carries what the routine plots in
    <oname>_Distortion.ps (the residuals, their MOMENTs, the RMS titles as
    printed, with the plate scale `ps`) and the line labels of the points
    each fit rejected. The IDL's printed reports and overlay points are
    collected by distortion_report().
    """
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
        sm = I.moments(sres[sgb == 1])
        lm = I.moments(lres[wgb == 1])
        sw = L["sw"]
        sb, lb = np.where(sgb == 0)[0], np.where(wgb == 0)[0]
        st = dict(order=o, nlines=int(np.unique(L["sl"]).size),
                  spatial_rms=float(np.std(sres[sgb == 1], ddof=1)),
                  wave_rms_pix=float(np.std(lres[wgb == 1], ddof=1)),
                  sx=sx, sy=L["sy"], sres=sres, lres=lres, sgb=sgb, wgb=wgb,
                  sy_rel=sy, xrange=[int(xranges[i][0]), int(xranges[i][1])],
                  smoment=(float(sm[0]), float(sm[1])), lmoment=(float(lm[0]), float(lm[1])),
                  bad_spatial_lines=[_strtrim_double(v) for v in sw[sb][_idl_uniq(sw[sb])]],
                  bad_wave_lines=[_strtrim_double(v) for v in sw[lb][_idl_uniq(sw[lb])]],
                  wave_title="RMS = %6.3f pixels" % np.sqrt(lm[1]))
        if ps is not None:
            # sqrt(variance / ps), as the IDL title has it (not sqrt(variance) / ps)
            st["ps"] = float(I.FL(ps))
            st["spatial_title"] = "RMS = %6.4f pixels" % np.sqrt(sm[1] / I.FL(ps))
        stats.append(st)
    return warr, sarr, wcoeffs, scoeffs, disp, stats


def distortion_report(stats, edgecoeffs, lines2d):
    """
    What MC_WAVECAL2D_KNAHACK prints and what TS41CALS_MKWAVECAL overlays on
    the arc in ximgtool, from wavecal2d()'s stats.

    badslines / badwlines are the printed lists of lines with a point
    rejected by the spatial / wavelength fit, built as the IDL builds them
    (orders joined with no separator, the labels as SW holds them);
    bad_spatial / bad_wave are the same per order. The overlay:

      traced     every point kept (SX, SY of each order; the IDL's list
                 starts with order 3 twice), drawn as green squares
      spatial    points the spatial fit rejected (red '+')
      wave_idl   the IDL's 'wavelength-rejected' points (blue triangles):
                 SX[SL] + the bottom edge, for every point of an order with
                 any rejection. SL holds column numbers, so these are not
                 the rejected points; kept as the IDL draws them.
      wave       the points the wavelength fit actually rejected (not drawn
                 by the IDL).
    """
    badslines = "".join(", ".join(s["bad_spatial_lines"]) for s in stats)
    badwlines = "".join(", ".join(s["bad_wave_lines"]) for s in stats)
    cx = [lines2d[0]["sx"]] + [L["sx"] for L in lines2d]
    cy = [lines2d[0]["sy"]] + [L["sy"] for L in lines2d]
    pxs, pys, pxl, pyl, wx, wy = [], [], [], [], [], []
    for i, s in enumerate(stats):
        sx, sy = s["sx"], s["sy_rel"]
        sb = np.where(s["sgb"] == 0)[0]
        lb = np.where(s["wgb"] == 0)[0]
        if sb.size:
            pxs.append(sx[sb])
            pys.append(sy[sb] + I.poly(sx[sb], edgecoeffs[i, 0]))
        if lb.size:
            # a FLOAT subscript truncates, and out-of-range ones are clipped
            k = np.clip(np.trunc(lines2d[i]["sl"]).astype(int), 0, sx.size - 1)
            pxl.append(sx[k])
            pyl.append(sy[k] + I.poly(sx[k], edgecoeffs[i, 0]))
            wx.append(sx[lb])
            wy.append(sy[lb] + I.poly(sx[lb], edgecoeffs[i, 0]))
    cat = lambda v: np.concatenate(v) if v else np.zeros(0)
    return dict(badslines=badslines, badwlines=badwlines,
                badslines_printed=I.splittext(badslines, 80),
                badwlines_printed=I.splittext(badwlines, 80),
                bad_spatial={s["order"]: s["bad_spatial_lines"] for s in stats},
                bad_wave={s["order"]: s["bad_wave_lines"] for s in stats},
                traced=(cat(cx), cat(cy)), spatial=(cat(pxs), cat(pys)),
                wave_idl=(cat(pxl), cat(pyl)), wave=(cat(wx), cat(wy)))


# ------------------------------------------------------------- the whole step

ARCTYPES = ("Sky", "CuHeAr Arc")


def xgetoffset(whome, fhome, guess_home, wobs, fobs, guess_obs):
    """
    XGETOFFSET: the offset the user gives by marking the same line in the
    template (whome, fhome) and in the observed home order (wobs, fobs):
    the observed centre minus the template's, a FLOAT (getoffset.xgetoffset).
    make_wavecal() returns the four arrays as result['xgetoffset'].
    """
    from triplespec_reducer.getoffset import xgetoffset as _xgetoffset
    return _xgetoffset(whome, fhome, wobs, fobs, guess_home, guess_obs)[0]


def make_wavecal(files, flatinfo, arctype="Sky", modeinfo=None, progress=None, offset=None,
                 short_lines="abort"):
    """
    TS41CALS_MKWAVECAL.

    files       the frames to calibrate on, any list (the IDL's Wavecal
                Images field takes any frames): a sky sequence or part of
                one, short exposures, the arcs or some of them.
    flatinfo    the flat to use (the IDL's Full Flat Name): read_flat()'s
                dict, or the path of a flat file.
    arctype     "Sky" (OH lines, OHXD_LineCal) or "CuHeAr Arc" (Arll.dat,
                CuHeArXD_LineCal), the IDL's free choice.
    offset      replaces the automatic cross-correlation (xspextool's
                'Auto X-Correlate' off): the value XGETOFFSET returns, see
                xgetoffset().
    short_lines "abort" (IDL) or "skip"; see findlines2d().

    The result carries the calibration QA: sol1d['qa'] (the _Wavecal
    plot), stats (the _Distortion plot) and report (the printed lists of
    rejected lines and the points the IDL overlays on the arc); write them
    with write_wavecal_qa().
    """
    if arctype not in ARCTYPES:
        raise ValueError(f"arctype must be one of {ARCTYPES}, not {arctype!r}")
    if isinstance(flatinfo, (str, os.PathLike)):
        from triplespec_reducer.flat import read_flat
        flatinfo = read_flat(str(flatinfo))
    files = [str(f) for f in files]
    if not files:
        raise ValueError("no wavecal frames given")
    cuhear = arctype == "CuHeAr Arc"         # IDL: anything else is Sky
    mi = modeinfo or read_modefile()
    say = progress or (lambda *_: None)
    rot = flatinfo["rotation"]
    orders = flatinfo["orders"]
    ec, xr = flatinfo["edgecoeffs"], flatinfo["xranges"]
    nrows, ncols = flatinfo["flat"].shape
    say("Loading the arc frames" if cuhear else "Loading the sky frames")
    data, var, hdrs, _ = read_ts4(files, rotation=rot)
    arc = (I.median(data, even=True, axis=0) if data.shape[0] != 1 else data[0])
    arc = (arc / flatinfo["flat"]).astype(np.float32)
    say("Wavelength calibrating")
    wfake, sfake = simwavecal2d(ncols, nrows, ec, xr, flatinfo["slith_arc"])
    sa = flatinfo["slith_arc"]
    spec = medextspec1d(data[0], var[0], flatinfo["omask"], orders, wfake, sfake,
                        [1], [[sa / 2. - 5, sa / 2 + 5]] * len(orders))
    calname = "CuHeArXD_LineCal.fits" if cuhear else "OHXD_LineCal.fits"
    calinfo = read_linecal(os.path.join(INSTRUMENT_DIR, calname))
    zh = orders.index(calinfo["xcororder"])
    idx = nantrim(calinfo["xcorspec"][0], 2)
    xgo = (calinfo["xcorspec"][0][idx], calinfo["xcorspec"][1][idx], spec[zh][0], spec[zh][1])
    if offset is None:
        offset, xc = xcorr_offset(calinfo, spec[zh], ncols)
    else:
        # XGETOFFSET hands back a FLOAT (pos2 - pos1 of its FLOAT fields)
        offset, xc = float(np.float32(offset)), None
    llname = "Arll.dat" if cuhear else "OHll.dat"
    lines = read_linelist(os.path.join(INSTRUMENT_DIR, llname))
    lt = "Ar" if cuhear else "OH"
    w1 = xdwavecal1d(spec, orders, calinfo["w2pcoeffs"], offset, lines,
                     calinfo["ordertype"], [lt] * lines.size, flatinfo["slitw_pix"],
                     calinfo["homeorder"], mi.dispdeg, mi.ordrdeg, 3.5)
    say("Tracing the lines along the slit")
    plines = rdohll(os.path.join(INSTRUMENT_DIR, llname), res=flatinfo["rp"], thresh=3)
    w2 = findlines2d(arc, ec, xr, orders, plines, calinfo["w2pcoeffs"], offset,
                     w1["wranges"], flatinfo["slitw_pix"], sa, mi.linedeg,
                     mi.linereg / 2., mi.ystep, mi.ysum, short_lines=short_lines)
    say("Fitting the 2D solution")
    wcinfo1d = dict(coeffs=w1["coeffs"], dispdeg=mi.dispdeg, ordrdeg=mi.ordrdeg,
                    homeorder=mi.homeorder)
    warr, sarr, wc, sc, disp, stats = wavecal2d(
        wcinfo1d, w2, ec, flatinfo["omask"], orders, xr, w1["wranges"],
        mi.wxdeg, mi.wydeg, mi.sxdeg, mi.sydeg, ps=mi.ps)
    return dict(wavecal=warr, spatcal=sarr, wcinfo1d=wcinfo1d, wcoeffs=wc,
                scoeffs=sc, disp=disp, orders=orders, rotation=rot, offset=offset,
                xcorr=xc, xgetoffset=xgo, spec1d=spec, sol1d=w1, lines2d=w2, stats=stats,
                report=distortion_report(stats, ec, w2), arc=arc, edgecoeffs=ec, xranges=xr,
                modeinfo=mi, files=files, arctype=arctype)


# ------------------------------------------------------------- QA plots

# IDL colour indices 2, 3, 4, ... (red, green, blue, ...) of MC_MKCT, for the orders
_ORDER_COLORS = ["#eb6834", "#1baf7a", "#2a78d6", "#eda100", "#e87ba4", "#008300",
                 "#7b4fd6", "#5c5c5c"]


def plot_wavecal1d(path, qa, orders):
    """
    The plot MC_XDWAVECAL1D writes to <oname>_Wavecal.eps (8 x 10.5 in):
    the 1D-fit residuals (Angstrom) against order number and against
    column, coloured by order, rejected lines boxed, the mean (dashed) and
    mean +- 1 sigma (dotted) of the lines kept, the RMS in the title.
    """
    from matplotlib.figure import Figure
    so, sx, res, gb = qa["order"], qa["column"], qa["resid"], qa["goodbad"]
    m0, sd = qa["mean"], np.sqrt(qa["var"])
    fig = Figure(figsize=(8, 10.5))
    axs = fig.subplots(2, 1)
    bad = gb == 0
    for ax, xv, xl in ((axs[0], so, "Order Number"), (axs[1], sx, "Column (pixels)")):
        for i, o in enumerate(orders):
            z = so.astype(int) == o
            if z.any():
                ax.plot(xv[z], res[z], "o", ms=5, mfc=_ORDER_COLORS[i % 8], mec="k", mew=0.6,
                        label=f"Order {o}")
        if bad.any():
            ax.plot(xv[bad], res[bad], "s", ms=10, mfc="none", mec="k", label="rejected")
        for v, ls in ((m0, "--"), (m0 - sd, ":"), (m0 + sd, ":")):
            ax.axhline(v, color="k", ls=ls, lw=1)
        ax.set_xlabel(xl)
        ax.set_ylabel("Data-Model (Angstroms)")
    axs[0].set_xlim(so.min() - 0.5, so.max() + 0.5)
    axs[0].set_title(qa["title"])
    axs[0].legend(fontsize=8, ncol=3)
    fig.tight_layout()
    fig.savefig(path)
    return path


def plot_distortion(path, stats):
    """
    The plot MC_WAVECAL2D_KNAHACK writes to <oname>_Distortion.ps, one
    landscape page (10.5 x 8 in) per order: spatial residual (pixels)
    against X and Y, wavelength residual (pixels) against X and Y, points
    the fits rejected in blue, zero in red, the mean (dashed) and +-1 sigma
    (dotted), and the IDL's RMS titles. PDF (several pages) or, for any
    other extension, one image per order named <path stem>_<order>.<ext>.
    """
    from matplotlib.figure import Figure
    pages = []
    for s in stats:
        ps = s.get("ps", 1.0)
        sx, sy = s["sx"], s["sy_rel"]
        ds, dl = s["sres"] / ps, s["lres"]
        sb, lb = s["sgb"] == 0, s["wgb"] == 0
        (sm0, smv), (lm0, lmv) = s["smoment"], s["lmoment"]
        fig = Figure(figsize=(10.5, 8))
        ax = fig.subplots(2, 2)
        panels = [
            (ax[0, 0], sx, ds, sb, "X (pixels)", "Spatial-Model (pixels)", s.get("spatial_title", ""),
             "h", sm0 / ps, np.sqrt(smv) / ps),
            (ax[0, 1], ds, sy, sb, "Spatial-Model (pixels)", "Y (pixels)", "", "v", sm0 / ps,
             np.sqrt(smv) / ps),
            (ax[1, 0], sx, dl, lb, "X (pixels)", "Wavelength-Model (pixels)", s["wave_title"],
             "h", lm0, np.sqrt(lmv)),
            (ax[1, 1], dl, sy, lb, "Wavelength-Model (pixels)", "Y (pixels)", "", "v", lm0,
             np.sqrt(lmv))]
        for a, xv, yv, bad, xl, yl, title, orient, m, sd in panels:
            a.plot(xv, yv, "o", ms=2.5, color="k")
            if bad.any():
                a.plot(xv[bad], yv[bad], "o", ms=2.5, color="#2a78d6")
            line = a.axhline if orient == "h" else a.axvline
            line(0, color="#eb6834", lw=1)
            line(m, color="#eb6834", ls="--", lw=1)
            line(m - sd, color="#eb6834", ls=":", lw=1)
            line(m + sd, color="#eb6834", ls=":", lw=1)
            a.set_xlabel(xl)
            a.set_ylabel(yl)
            if title:
                a.set_title(title)
        ax[0, 0].set_xlim(*s["xrange"])
        ax[1, 0].set_xlim(*s["xrange"])
        fig.suptitle(f"Order {s['order']:02d}", fontsize=16)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        pages.append(fig)
    stem, ext = os.path.splitext(path)
    if ext.lower() == ".pdf":
        from matplotlib.backends.backend_pdf import PdfPages
        with PdfPages(path) as pdf:
            for f in pages:
                pdf.savefig(f)
        return [path]
    out = []
    for s, f in zip(stats, pages):
        p = f"{stem}_{s['order']:02d}{ext}"
        f.savefig(p)
        out.append(p)
    return out


def write_wavecal_qa(r, base, ext="pdf"):
    """
    The QA files TS41CALS_MKWAVECAL leaves next to the wavecal, for the
    output name `base` (path without '.fits'): <base>_Wavecal.<ext>,
    <base>_Distortion.<ext>, and <base>_Distortion.txt with the lists of
    lines whose points were rejected, as MC_WAVECAL2D_KNAHACK prints them.
    Returns the paths written.
    """
    out = [plot_wavecal1d(f"{base}_Wavecal.{ext}", r["sol1d"]["qa"], r["orders"])]
    out += plot_distortion(f"{base}_Distortion.{ext}", r["stats"])
    rep = r["report"]
    txt = f"{base}_Distortion.txt"
    with open(txt, "w") as fh:
        fh.write("\nThe following lines had points identified as bad in the spatial fit:\n\n")
        fh.write("\n".join(rep["badslines_printed"]) + "\n\n")
        fh.write("The following lines had points identified as bad in the wavelength fit:\n\n")
        fh.write("\n".join(rep["badwlines_printed"]) + "\n\n")
    out.append(txt)
    return out


def write_wavecal(path, r):
    """The file MC_WAVECAL2D_KNAHACK writes, read back by MC_RDWAVECAL2D."""
    mi = r["modeinfo"]
    w1 = r["wcinfo1d"]
    hdr = fits.Header()
    # numbers as FXADDPAR/SXADDPAR write them; in MC_WAVECAL2D_KNAHACK's order
    I.fxaddpar(hdr, "WCTYPE", "2D", " Wavelength calibration type")
    I.fxaddpar(hdr, "WDEG", int(w1["dispdeg"]), " X Polynomial degree for 1D wavelength calibration")
    I.fxaddpar(hdr, "ODEG", int(w1["ordrdeg"]), " Y Polynomial degree for 1D wavelength calibration")
    I.fxaddpar(hdr, "HOMEORDR", f"{int(w1['homeorder']):02d}", " Home order for 1D wavelength calibration")
    I.fxaddpar(hdr, "WXDEG", int(mi.wxdeg), ' X Polynomial degree for 2D "wavelength" distortion calibration')
    I.fxaddpar(hdr, "WYDEG", int(mi.wydeg), ' Y Polynomial degree for 2D "wavelength" distortion calibration')
    I.fxaddpar(hdr, "SXDEG", int(mi.sxdeg), ' X Polynomial degree for 2D "spatial" distortion calibration')
    I.fxaddpar(hdr, "SYDEG", int(mi.sydeg), ' Y Polynomial degree for 2D "spatial" distortion calibration')
    I.fxaddpar(hdr, "NORDERS", len(r["orders"]), " Number of orders")
    I.fxaddpar(hdr, "ORDERS", ",".join(str(o) for o in r["orders"]), " Order numbers")
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


def no_wavecal(flatinfo):
    """
    What XSPEXTOOL_LOADIMAGES loads when the Full Wavecal Name field is
    blank: MC_SIMWAVECAL2D's images (column number as 'wavelength', arcsec
    along the slit) and wctype '1D', which sends the extraction to the
    SpeX-type MC_EXTPSSPEC / MC_EXTXSSPEC, with spectra in pixels. Same
    keys as read_wavecal(); there is no dispersion.
    """
    nrows, ncols = flatinfo["flat"].shape
    wavecal, spatcal = simwavecal2d(ncols, nrows, flatinfo["edgecoeffs"],
                                    flatinfo["xranges"], flatinfo["slith_arc"])
    return dict(wavecal=wavecal, spatcal=spatcal, header=None, wctype="1D",
                orders=list(flatinfo["orders"]), disp=None)
