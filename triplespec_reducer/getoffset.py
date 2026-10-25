"""
The manual pixel offset of the wavelength calibration: XGETOFFSET.

When xspextool's 'Auto X-Correlate' is off, KNA_TS4CALS2D shows the
template spectrum of the home order (the line-calibration file's XCORSPEC)
and the observed one, and the user clicks the same line in each. In each
spectrum the line is fitted, over 15 pixels either side of the click, with
MPFITEXPR 'p[0] + p[1]*x + gaussian(x,p[2:4])' started from
FLOAT([0, 0, y at the click, click, 2]); the offset is the observed centre
minus the template centre, in FLOAT (XGETOFFSET's 'Accept').

    line_centre(x, y, click)          XGETOFFSET_ZOOM for one spectrum
    xgetoffset(x1, y1, x2, y2, c1, c2) the offset, as XGETOFFSET returns it
    home_spectra(...)                 the two spectra KNA_TS4CALS2D passes
"""

import math
import os

import numpy as np

from triplespec_reducer import idl as I
from triplespec_reducer.mpfit import mpfit, _is_double

EXPR = "p[0] + p[1]*x + gaussian(x,p[2:4])"


def _mc_tabinv(xarr, x):
    """MC_TABINV: TABINV on the NaN-trimmed array, plus the first index (FLOAT + LONG)."""
    from triplespec_reducer.wavecal import nantrim
    k = nantrim(xarr, 2)
    t = I.tabinv(np.asarray(xarr)[k], x)
    if I.isf32(t):
        return np.float32(t + np.float32(k[0]))
    return np.float64(t + k[0])


def _gaussian(x, p, D):
    """
    astrolib GAUSSIAN without /DOUBLE: the Gaussian is stored in an FLTARR
    and the underflow guard uses the single-precision MACHAR, whatever the
    type of x and p; times p[0] in their type.
    """
    z = ((x - p[1]) / p[2]).astype(D)
    zz = (z * z).astype(D)
    minexp = I.alog(np.float32(np.finfo(np.float32).tiny))
    w = zz < np.float32(-2) * minexp
    g = np.zeros(x.size, dtype=np.float32)
    a = (-zz[w] / D(2)).astype(D)
    if D is np.float32:
        g[w] = I.exp(a)
    else:
        g[w] = np.array([math.exp(v) for v in a], dtype=np.float64)
    return (p[0] * g).astype(D)


def fit_line(xwin, ywin, start):
    """
    MPFITEXPR('p[0] + p[1]*x + gaussian(x,p[2:4])', xwin, ywin, undefined
    error, start, /QUIET, /NOCOVAR): plain residuals, in DOUBLE when the data
    are DOUBLE and in FLOAT otherwise (MPFIT's precision is the residuals').
    Returns the five parameters and MPFIT's status.
    """
    D = np.float64 if _is_double(xwin, ywin) else np.float32
    x = np.asarray(xwin, dtype=D)
    y = np.asarray(ywin, dtype=D)

    def fcn(p):
        f = ((p[0] + (p[1] * x).astype(D)).astype(D) + _gaussian(x, p[2:5], D)).astype(D)
        return (y - f).astype(D)

    p, status = mpfit(fcn, np.asarray(start, dtype=np.float32), dtype=D)
    return p, status


def line_centre(x, y, click):
    """
    XGETOFFSET_ZOOM for one spectrum: `click` is where the user clicked (the
    FLOAT cursor x). The window is x[bot:top] with bot = 0 > (i-15) and top =
    (i+15) < (n-1), i the MC_TABINV index of the click (fractional, truncated
    as a subscript); the start FLOAT([0, 0, y[ROUND(i')], click, 2]), i' the
    index of the click in the window. Returns (centre as the FLOAT pos field
    holds it, the fit, the window's x, the window's y, MPFIT's status).
    """
    x, y = np.asarray(x), np.asarray(y)
    click = np.float32(click)                      # the FLOAT guesspos field
    idx = _mc_tabinv(x, click)
    lo = max(idx - 15, 0)                          # 0 > (idx-15)
    hi = min(idx + 15, x.size - 1)                 # (idx+15) < (n-1)
    bot, top = I.itrunc(lo), I.itrunc(hi)
    xwin, ywin = x[bot:top + 1], y[bot:top + 1]
    i2 = _mc_tabinv(xwin, click)
    start = np.array([0., 0., ywin[I.iround(i2)], click, 2], dtype=np.float32)
    fit, status = fit_line(xwin, ywin, start)
    return np.float32(fit[3]), fit, xwin, ywin, status


def xgetoffset(x1, y1, x2, y2, click1, click2):
    """
    XGETOFFSET: the offset of spectrum 2 (observed) from spectrum 1
    (template), pos2 - pos1 in FLOAT, from one click on the same line in
    each. Returns (offset, details) with each spectrum's line_centre result.
    """
    c1 = line_centre(x1, y1, click1)
    c2 = line_centre(x2, y2, click2)
    offset = np.float32(c2[0] - c1[0])
    return offset, dict(template=c1, observed=c2)


def home_spectra(files, flatinfo, arctype="Sky"):
    """
    The two spectra KNA_TS4CALS2D passes to XGETOFFSET: the template (the
    line-calibration file's XCORSPEC, NaN-trimmed with MC_NANTRIM) and the
    home order of the arc or sky frames extracted as wavecal.make_wavecal
    extracts them. Returns (whome, fhome, wobs, fobs, home order).
    """
    from triplespec_reducer import wavecal as W
    from triplespec_reducer.frames import read_ts4
    from triplespec_reducer.instrument import INSTRUMENT_DIR
    rot = flatinfo["rotation"]
    orders = flatinfo["orders"]
    nrows, ncols = flatinfo["flat"].shape
    data, var, _, _ = read_ts4(files, rotation=rot)
    sa = flatinfo["slith_arc"]
    wfake, sfake = W.simwavecal2d(ncols, nrows, flatinfo["edgecoeffs"], flatinfo["xranges"], sa)
    spec = W.medextspec1d(data[0], var[0], flatinfo["omask"], orders, wfake, sfake,
                          [1], [[sa / 2. - 5, sa / 2 + 5]] * len(orders))
    calname = "CuHeArXD_LineCal.fits" if arctype != "Sky" else "OHXD_LineCal.fits"
    calinfo = W.read_linecal(os.path.join(INSTRUMENT_DIR, calname))
    whome, fhome = calinfo["xcorspec"][0], calinfo["xcorspec"][1]
    k = W.nantrim(whome, 2)
    zh = orders.index(calinfo["xcororder"])
    return whome[k], fhome[k], spec[zh][0], spec[zh][1], calinfo["xcororder"]
