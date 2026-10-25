"""
Numerical building blocks with the behaviour of IDL 7.1's and astrolib's.

The science modules of this package follow the steps of Spextool 4.1 and its
TS4_SOAR module, and their results match the IDL tool's bit for bit only if
the building blocks behave as IDL's do. The places where IDL's behaviour
differs from the obvious numpy spelling are the reason this module exists:

* ``median`` without /EVEN returns the upper of the two middle values of an
  even-length array, and ignores NaNs.
* ``CONVOL`` correlates (it does not flip the kernel) and zeroes the edges.
* ``INTERPOL`` extrapolates linearly; ``INTERPOLATE`` clamps to the edge.
* ``ROUND`` rounds halves away from zero; float subscripts truncate.
* ``GAUSSFIT`` stops CURVEFIT at a relative chi-square change of 1e-3, so a
  fully converged fit is not the same number.

Arrays follow numpy's [row, column] order. Where an IDL routine indexes an
image as [x, y], the caller translates; nothing here transposes silently.
Each function names the IDL routine whose behaviour it matches.
"""

import numpy as np
from scipy.spatial import Delaunay


# ---------------------------------------------------------------- scalars

def isf32(a):
    """True for single-precision data, whatever its byte order (FITS is big-endian)."""
    dt = np.asarray(a).dtype
    return dt.kind == "f" and dt.itemsize == 4


def isf64(a):
    """True for double-precision data, whatever its byte order."""
    dt = np.asarray(a).dtype
    return dt.kind == "f" and dt.itemsize == 8


# ---------------------------------------------------------------- functions
#
# IDL 7.1 evaluates EXP, ALOG, SIN and COS of FLOAT data with the C library's
# single-precision functions (expf, logf, ...), which numpy's own float32
# routines do not reproduce (they differ in a third of the values); numba
# compiles math.exp on float32 to the same libm calls. ATAN, ALOG10 and ^ of
# FLOAT data are computed in double and rounded. Double data: numpy is IDL.
# All checked on 200 000 values each against IDL.

import math as _math

try:
    import numba as _nb

    @_nb.vectorize(["float32(float32)"], cache=True)
    def _expf(v):
        return _math.exp(v)

    @_nb.vectorize(["float32(float32)"], cache=True)
    def _logf(v):
        return _math.log(v)

    @_nb.vectorize(["float32(float32)"], cache=True)
    def _sinf(v):
        return _math.sin(v)

    @_nb.vectorize(["float32(float32)"], cache=True)
    def _cosf(v):
        return _math.cos(v)
    @_nb.vectorize(["float64(float64, float64)"], cache=True)
    def _pow64(v, e):
        return _math.pow(v, e)
except ImportError:                                   # pragma: no cover
    _expf = _logf = _sinf = _cosf = None
    _pow64 = np.power


def _libm1(fn, npfn):
    def f(a):
        if isf32(a):
            a32 = np.asarray(a, dtype=np.float32)
            if fn is None:                            # pragma: no cover
                return npfn(a32.astype(float)).astype(np.float32)
            out = fn(a32)
            return out if np.ndim(a) else np.float32(out)
        return npfn(a)
    return f


exp = _libm1(_expf, np.exp)
alog = _libm1(_logf, np.log)
sin = _libm1(_sinf, np.sin)
cos = _libm1(_cosf, np.cos)


def _via_double(npfn):
    def f(a):
        if isf32(a):
            return np.asarray(npfn(np.asarray(a, dtype=float))).astype(np.float32)
        return npfn(a)
    return f


atan = _via_double(np.arctan)
alog10 = _via_double(np.log10)


def dpow(a, b):
    """
    IDL a^b with a FLOAT or DOUBLE exponent on double data: the C library's
    pow(). numpy's x ** 2.0 is the exact square, which pow() is not always.
    """
    return _pow64(np.asarray(a, dtype=float), np.asarray(b, dtype=float))


def fpow(a, b):
    """IDL a^b for a non-integer b: double pow(), rounded to FLOAT when IDL's result is FLOAT."""
    r = dpow(a, b)
    single = (isf32(a) or np.ndim(a) == 0 and not isinstance(a, float)) and \
        (isf32(b) or np.ndim(b) == 0 and not isinstance(b, float))
    return r.astype(np.float32) if single else r


def sq(x):
    """IDL x^2, which is x*x. numpy's scalar x ** 2 goes through pow() and is
    not always the product, in single or double precision."""
    return x * x


def FL(v):
    """
    An IDL FLOAT literal (2.354, 0.13, ...) as it enters a DOUBLE expression:
    the single-precision value, promoted. Python's 2.354 is a different number.
    """
    return np.float64(np.float32(v))


def iround(x):
    """IDL ROUND: halves go away from zero."""
    x = np.asarray(x, dtype=float)
    out = np.sign(x) * np.floor(np.abs(x) + 0.5)
    return out.astype(np.int64) if out.ndim else int(out)


def itrunc(x):
    """An IDL float subscript or FIX(): truncation toward zero."""
    x = np.asarray(x, dtype=float)
    out = np.trunc(x).astype(np.int64)
    return out if out.ndim else int(out)


# ---------------------------------------------------------------- medians

def median(a, even=False, axis=None):
    """
    IDL MEDIAN(array, /EVEN, DIMENSION=d), the keywords optional.

    NaNs are ignored; infinities are sorted with the other values, as IDL
    does. An even number of values gives the upper middle one unless
    `even`, in which case the two middle values are averaged. `axis` is
    the numpy axis that IDL's DIMENSION collapses.
    """
    P = np.float32 if isf32(a) else np.float64
    a = np.asarray(a, dtype=float)
    if axis is None:
        v = np.sort(a[~np.isnan(a)].ravel())
        n = v.size
        if n == 0:
            return P(np.nan)
        if n % 2 == 1 or not even:
            return P(v[n // 2])
        return P(0.5 * (v[n // 2 - 1] + v[n // 2]))
    s = np.sort(a, axis=axis)                       # NaNs sort to the end
    n = np.sum(~np.isnan(a), axis=axis)
    hi = np.take_along_axis(s, np.expand_dims(n // 2, axis), axis=axis)
    hi = np.squeeze(hi, axis=axis)
    if not even:
        out = hi
    else:
        lo = np.take_along_axis(s, np.expand_dims(np.maximum(n // 2 - 1, 0),
                                                  axis), axis=axis)
        lo = np.squeeze(lo, axis=axis)
        out = np.where(n % 2 == 1, hi, 0.5 * (lo + hi))
    return np.where(n > 0, out, np.nan).astype(P)


def median_filter(img, width):
    """
    IDL MEDIAN(image, width): a width x width box median.

    Pixels closer than width/2 to an edge are returned unchanged, as IDL
    does. The box median of an odd number of values needs no /EVEN rule.
    A 1-D even width w (probed in IDL 7.1) takes x[i-w/2 : i+w/2-1] for
    i = w/2 ... n-1-w/2, NaNs ignored, the upper middle of an even count.
    """
    from numpy.lib.stride_tricks import sliding_window_view
    img = np.asarray(img, dtype=float)
    out = img.copy()
    h = width // 2
    if img.ndim == 1:
        if img.size - 2 * h <= 0:
            return out
        win = sliding_window_view(img, width)[:img.size - 2 * h]
        out[h:img.size - h] = median(win, axis=-1)
        return out
    win = sliding_window_view(img, (width, width))
    win = win.reshape(win.shape[0], win.shape[1], -1)
    out[h:img.shape[0] - h, h:img.shape[1] - h] = median(win, axis=-1)
    return out


# ---------------------------------------------------------------- lookups

def value_locate(xarr, x):
    """IDL VALUE_LOCATE for a monotonic (either direction) vector; a NaN gives -1 (probed)."""
    xarr = np.asarray(xarr, dtype=float)
    x = np.asarray(x, dtype=float)
    if xarr[-1] >= xarr[0]:
        idx = np.searchsorted(xarr, x, side="right") - 1
    else:
        idx = np.searchsorted(-xarr, -x, side="right") - 1
    return np.where(np.isnan(x), -1, idx)


def _is_double(*a):
    """True when any argument would be DOUBLE in IDL (Python floats are)."""
    return any(not isf32(v) for v in a)


def tabinv(xarr, x):
    """
    astrolib TABINV: the effective (fractional) index of x in xarr.

    Computed in double, but stored in a FLOAT array, as TABINV does, unless
    xarr or x is double.
    """
    single = not _is_double(xarr, x)
    xarr = np.asarray(xarr, dtype=float)
    scalar = np.ndim(x) == 0
    x = np.atleast_1d(np.asarray(x, dtype=float))
    npt = xarr.size - 1
    ieff = value_locate(xarr, x).astype(float)
    g = (ieff < npt) & (ieff >= 0)
    if g.any():
        neff = ieff[g].astype(int)
        x0 = xarr[neff]
        ieff[g] = neff + (x[g] - x0) / (xarr[neff + 1] - x0)
    ieff = np.maximum(ieff, 0.0)
    if single:
        ieff = ieff.astype(np.float32)
    return ieff[0] if scalar else ieff


def findidx(xarr, x):
    """Spextool MC_FINDIDX: like TABINV, but NaN outside the table."""
    xarr = np.asarray(xarr, dtype=float)
    scalar = np.ndim(x) == 0
    x = np.atleast_1d(np.asarray(x, dtype=float))
    ieff = tabinv(xarr, x)
    lo, hi = (xarr[0], xarr[-1]) if xarr[-1] >= xarr[0] else (xarr[-1], xarr[0])
    ieff = np.where((x < xarr[0]) | (x > xarr[-1]), np.nan, ieff) \
        if xarr[-1] >= xarr[0] else np.where((x < lo) | (x > hi), np.nan, ieff)
    return np.float64(ieff[0]) if scalar else ieff


def interpolate1d(y, r, missing=None):
    """
    IDL INTERPOLATE(y, r): linear, clamped to the ends unless MISSING.

    IDL 7.1 rounds the subscripts to FLOAT, even double ones, and computes
    y0 + f (y1 - y0) in double, the difference y1 - y0 being taken in the
    precision of y; the result is FLOAT for a float y. A NaN subscript reads
    element 0, or takes MISSING when that is given. Matched to IDL element
    for element, MISSING included.
    """
    single = isf32(y)
    yp = np.asarray(y, dtype=np.float32 if single else float)
    y = yp.astype(float)
    r = np.asarray(r, dtype=float).astype(np.float32).astype(float)
    n = y.size
    bad = ~np.isfinite(r)
    rs = np.where(bad, 0.0, r)
    rc = np.clip(rs, 0, n - 1)
    i0 = np.clip(np.floor(rc).astype(int), 0, max(n - 2, 0))
    f = rc - i0
    i1 = np.minimum(i0 + 1, n - 1)
    out = y[i0] + f * (yp[i1] - yp[i0]).astype(float)
    if missing is not None:
        # IDL 7.1 (1-D): a subscript between n-1 and n still reads the last
        # element; only one below 0 or at/after n takes MISSING.
        out = np.where(bad | (rs < 0) | (rs >= n), missing, out)
    return out.astype(np.float32) if single else out


def linterp(xtab, ytab, xint, missing=None):
    """astrolib LINTERP; FLOAT for a float ytab, as INTERPOLATE returns."""
    xtab = np.asarray(xtab)
    ytab = np.asarray(ytab)
    if not isf32(ytab):
        ytab = ytab.astype(float)
    scalar = np.ndim(xint) == 0
    xi = np.atleast_1d(np.asarray(xint))
    npts = min(xtab.size, ytab.size)
    r = tabinv(xtab[:npts], xi)
    yi = interpolate1d(ytab[:npts], r)
    if missing is not None:
        xmin, xmax = min(xtab[0], xtab[npts - 1]), max(xtab[0], xtab[npts - 1])
        yi = np.where((xi < xmin) | (xi > xmax), missing, yi).astype(yi.dtype)
    return yi[0] if scalar else yi


def interpol(v, x, u):
    """
    IDL INTERPOL(v, x, u), linear, irregular grid: extrapolates past the
    ends. NaN points of v or x are dropped first, and the arithmetic is in
    the inputs' type (FLOAT unless one of them is double), as in IDL.
    """
    P = np.float64 if (isf64(v) or isf64(x) or isf64(u)) else np.float32
    v = np.asarray(v).astype(P)
    x = np.asarray(x).astype(P)
    u = np.asarray(u).astype(P)
    bad = np.isnan(v) | np.isnan(x)
    if bad.any() and (~bad).any():
        v, x = v[~bad], x[~bad]
    m = v.size
    s = np.clip(value_locate(x, u), 0, m - 2)
    diff = (v[s + 1] - v[s]).astype(P)
    return ((u - x[s]) * diff / (x[s + 1] - x[s]) + v[s]).astype(P)


def interpolv(v, x, u):
    """
    Spextool INTERPOLV(v, x, u), linear, irregular grid: INTERPOL's formula
    without its NaN removal, in the inputs' type.
    """
    P = np.float64 if (isf64(v) or isf64(x) or isf64(u)) else np.float32
    v = np.asarray(v).astype(P)
    x = np.asarray(x).astype(P)
    u = np.asarray(u).astype(P)
    s = np.clip(value_locate(x, u), 0, v.size - 2)
    return ((u - x[s]) * (v[s + 1] - v[s]) / (x[s + 1] - x[s]) + v[s]).astype(P)


def _cubic1d(p0, p1, p2, p3, t, a):
    """One cubic-convolution step, grouped as IDL evaluates it."""
    c1 = a * (p0 - p2)
    c2 = a * (2 * p2 - 2 * p0 + p3 - p1) + 3 * (p2 - p1)
    c3 = a * (p0 + p1 - p2 - p3) + 2 * (p1 - p2)
    return p1 + t * (c1 + t * (c2 + t * c3))


def interpolate2d(img, x, y, cubic=None, missing=None):
    """
    IDL INTERPOLATE(img, x, y, CUBIC=c, MISSING=m) at scattered points, the keywords optional.

    `img` is numpy [row, col]; x indexes columns, y rows, as IDL's
    INTERPOLATE(img_idl, x, y) does on the [x, y] array. Points outside the
    image take MISSING when given, otherwise the edge is clamped; a NaN
    coordinate reads 0, or takes MISSING.

    Reproduces IDL 7.1 element for element: the coordinates are rounded to
    FLOAT whatever their type; the cubic is evaluated in the precision of
    the image, as a polynomial in the fraction, along x on each of the four
    rows and then along y.
    """
    img = np.asarray(img)
    P = np.float64 if isf64(img) else np.float32
    img = img.astype(P)
    ny, nx = img.shape
    x = np.asarray(x, dtype=float).astype(np.float32)
    y = np.asarray(y, dtype=float).astype(np.float32)
    shape = np.broadcast(x, y).shape
    x = np.broadcast_to(x, shape).ravel()
    y = np.broadcast_to(y, shape).ravel()
    bad = ~(np.isfinite(x) & np.isfinite(y))
    xs = np.where(np.isfinite(x), x, np.float32(0))
    ys = np.where(np.isfinite(y), y, np.float32(0))
    out_of = bad | (xs < 0) | (xs > nx - 1) | (ys < 0) | (ys > ny - 1)
    xc = np.clip(xs, 0, nx - 1)
    yc = np.clip(ys, 0, ny - 1)
    if cubic is None:
        x0 = np.clip(np.floor(xc).astype(int), 0, nx - 2)
        y0 = np.clip(np.floor(yc).astype(int), 0, ny - 2)
        fx = xc.astype(float) - x0
        fy = yc.astype(float) - y0
        im = img.astype(float)
        r0 = im[y0, x0] + fx * (im[y0, x0 + 1] - im[y0, x0])
        r1 = im[y0 + 1, x0] + fx * (im[y0 + 1, x0 + 1] - im[y0 + 1, x0])
        out = r0 + fy * (r1 - r0)
    else:
        a = P(cubic)
        x0 = np.floor(xc).astype(int)
        y0 = np.floor(yc).astype(int)
        tx = (xc - x0).astype(P)
        ty = (yc - y0).astype(P)
        cols = [np.clip(x0 + i, 0, nx - 1) for i in range(-1, 3)]
        rows = []
        for j in range(-1, 3):
            yy = np.clip(y0 + j, 0, ny - 1)
            rows.append(_cubic1d(*(img[yy, c] for c in cols), tx, a))
        out = _cubic1d(*rows, ty, a)
    out = out.astype(P)
    if missing is not None:
        out = np.where(out_of, P(missing), out)
    return out.reshape(shape)


# ---------------------------------------------------------------- polynomials

def poly(x, c):
    """IDL POLY: FLOAT when x and the coefficients are both float."""
    x = np.asarray(x)
    c = np.asarray(c)
    P = np.float32 if (isf32(x) and isf32(c)) else np.float64
    x = x.astype(P)
    y = np.zeros_like(x)
    for k in range(len(c) - 1, -1, -1):
        y = (y * x + c[k]).astype(P)
    return y


def ipow(x, n):
    """
    IDL x^n for an INTEGER n >= 0: square-and-multiply (x^4 = (x^2)^2,
    x^5 = x x^4), which rounds differently from pow() for n >= 3. A FLOAT or
    FINDGEN exponent goes through pow(), as x ** float(n) does here.
    """
    x = np.asarray(x)
    r = None
    b = x
    n = int(n)
    if n == 0:
        return np.ones_like(x)
    while n:
        if n & 1:
            r = b if r is None else r * b
        n >>= 1
        if n:
            b = b * b
    return r


def poly1d(x, coeffs, coeffs_var=None):
    """Spextool POLY1D: the value and, when a (co)variance is given, its variance."""
    xx = np.asarray(x, dtype=float)
    coeffs = np.asarray(coeffs, dtype=float)
    y = np.zeros_like(xx)
    yvar = np.zeros_like(xx)
    # The loop index is an INT in the IDL, so xx^i is an integer power.
    for i, c in enumerate(coeffs):
        y = y + c * ipow(xx, i)
    if coeffs_var is not None:
        cv = np.asarray(coeffs_var, dtype=float)
        if cv.ndim == 1:
            for i in range(coeffs.size):
                yvar = yvar + cv[i] * sq(ipow(xx, i))
        else:
            for i in range(cv.shape[0]):
                for j in range(cv.shape[0]):
                    yvar = yvar + cv[i, j] * ipow(xx, i) * ipow(xx, j)
    return y, yvar


def poly2d(x, y, xorder, yorder, coeffs):
    """
    Spextool POLY2D: coefficient k = i + (xorder+1)*j multiplies x^i y^j.

    As in the IDL, TOTAL((x^xexp # y^yexp) * coeffs) with FINDGEN exponents:
    the product x^i y^j first, then the coefficient, summed in order.
    """
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    z = np.zeros(np.broadcast(xx, yy).shape)
    k = 0
    for j in range(yorder + 1):
        for i in range(xorder + 1):
            z = z + (dpow(xx, i) * dpow(yy, j)) * coeffs[k]
            k += 1
    return z


# IDL splits a TOTAL of TPOOL_MIN_ELTS or more elements over its thread
# pool: TPOOL_NTHREADS chunks of n // TPOOL_NTHREADS elements, the last one
# taking the remainder, each summed in order in the array's precision, the
# partial sums added in double. The last digit of such a sum therefore
# depends on the machine IDL runs on; these are the values of the machine the
# package was checked on (!CPU.TPOOL_NTHREADS = 64, IDL's defaults).
TPOOL_NTHREADS = 64
TPOOL_MIN_ELTS = 100000


def _tpool_bounds(n):
    c = n // TPOOL_NTHREADS
    b = [(k * c, (k + 1) * c) for k in range(TPOOL_NTHREADS - 1)]
    return b + [((TPOOL_NTHREADS - 1) * c, n)]


def itotal(a):
    """IDL TOTAL of a double array (or TOTAL(/DOUBLE)): in order, or by thread chunks."""
    a = np.asarray(a, dtype=float).ravel()
    if a.size == 0:
        return 0.0
    if a.size < TPOOL_MIN_ELTS:
        return np.float64(np.cumsum(a)[-1])
    parts = np.array([np.cumsum(a[i:j])[-1] for i, j in _tpool_bounds(a.size)])
    return np.float64(np.cumsum(parts)[-1])


def ftotal(a):
    """IDL TOTAL of a float array: a float running sum, or float chunks added in double."""
    a = np.asarray(a, dtype=np.float32).ravel()
    if a.size == 0:
        return np.float32(0)
    if a.size < TPOOL_MIN_ELTS:
        return np.cumsum(a, dtype=np.float32)[-1]
    parts = np.array([np.cumsum(a[i:j], dtype=np.float32)[-1]
                      for i, j in _tpool_bounds(a.size)], dtype=float)
    return np.float32(np.cumsum(parts)[-1])


def invert(a):
    """
    IDL INVERT: returns (inverse, status), in the precision of the input.

    In-place Gauss-Jordan elimination, the pivot of each row being the first
    largest element along that row (columns swapped into place). A zero pivot
    is not an error: IDL sets status 1 and returns the working matrix as it
    stands, swaps not undone, and CURVEFIT carries on with it. Matched to
    IDL 7.1 element for element on singular and regular matrices, in float
    and double; numpy's inverse rounds differently.
    """
    P = np.float64 if isf64(a) else np.float32
    b = np.array(a, dtype=P)
    n = b.shape[0]
    swaps = []
    for k in range(n):
        p = k + int(np.argmax(np.abs(b[k, k:])))
        if b[k, p] == 0:
            return b, 1
        if p != k:
            b[:, [k, p]] = b[:, [p, k]]
            swaps.append((k, p))
        pinv = P(1) / b[k, k]
        col = b[:, k].copy()
        row = (b[k] * pinv).astype(P)
        row[k] = pinv
        for i in range(n):
            if i != k:
                b[i] = b[i] - col[i] * row
                b[i, k] = -col[i] * pinv
        b[k] = row
    for k, p in reversed(swaps):
        b[[k, p]] = b[[p, k]]
    return b, 0


def ludc(a):
    """
    IDL LUDC (Numerical Recipes' ludcmp, Crout with implicit partial
    pivoting) in the matrix's precision. Returns (lu, index).
    """
    P = np.float64 if isf64(a) else np.float32
    a = np.array(a, dtype=P)
    n = a.shape[0]
    vv = np.zeros(n, dtype=P)
    for i in range(n):
        big = P(0)
        for j in range(n):
            big = max(big, abs(a[i, j]))
        vv[i] = P(1) / big
    indx = np.zeros(n, dtype=int)
    for j in range(n):
        for i in range(j):
            s_ = a[i, j]
            for k in range(i):
                s_ = P(s_ - a[i, k] * a[k, j])
            a[i, j] = s_
        big = P(0)
        imax = j
        for i in range(j, n):
            s_ = a[i, j]
            for k in range(j):
                s_ = P(s_ - a[i, k] * a[k, j])
            a[i, j] = s_
            dum = P(vv[i] * abs(s_))
            if dum >= big:
                big = dum
                imax = i
        if j != imax:
            a[[imax, j]] = a[[j, imax]]
            vv[imax] = vv[j]
        indx[j] = imax
        if a[j, j] == 0:
            a[j, j] = P(1e-20)
        if j != n - 1:
            dum = P(P(1) / a[j, j])
            for i in range(j + 1, n):
                a[i, j] = P(a[i, j] * dum)
    return a, indx


def lusol(lu, indx, b):
    """IDL LUSOL (Numerical Recipes' lubksb) in the matrix's precision."""
    P = lu.dtype.type
    b = np.array(b, dtype=P)
    n = b.size
    ii = -1
    for i in range(n):
        ip = indx[i]
        s_ = b[ip]
        b[ip] = b[i]
        if ii >= 0:
            for j in range(ii, i):
                s_ = P(s_ - lu[i, j] * b[j])
        elif s_ != 0:
            ii = i
        b[i] = s_
    for i in range(n - 1, -1, -1):
        s_ = b[i]
        for j in range(i + 1, n):
            s_ = P(s_ - lu[i, j] * b[j])
        b[i] = P(s_ / lu[i, i])
    return b


def poly_fit(x, y, ndegree):
    """
    IDL POLY_FIT (no weights): double normal equations summed in order,
    inverted by INVERT, coefficients b # covar; float unless x or y is
    double. A singular matrix gives NaN, as IDL does with STATUS.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    double = isf64(x) or isf64(y)
    m = ndegree + 1
    covar = np.zeros((m, m))
    b = np.zeros(m)
    z = np.ones(x.size)
    wy = y.astype(float)
    covar[0, 0] = x.size
    b[0] = itotal(wy)
    for p in range(1, 2 * ndegree + 1):
        z = z * x.astype(float)
        if p < m:
            b[p] = itotal(wy * z)
        s = itotal(z)
        for j in range(max(0, p - ndegree), min(ndegree, p) + 1):
            covar[p - j, j] = s
    covar, status = invert(covar)
    if status == 1:
        res = np.full(m, np.nan)
    else:
        res = np.array([itotal(b * covar[j, :]) for j in range(m)])
    return res if double else res.astype(np.float32)


def _gaussj(aa, bb):
    """
    Spextool MC_GAUSSJ: Gauss-Jordan elimination with full pivoting.

    Indexed as in Numerical Recipes (aa(j,k) is aa[j, k] here). The
    normal equations of a cubic in column number are close to singular,
    and a different solver rounds them to a visibly different polynomial:
    numpy's inverse put TS4's order-7 guess 128 pixels from IDL's.
    """
    a = np.array(aa, dtype=float)
    b = np.array(bb, dtype=float).ravel()
    n = a.shape[0]
    ipiv = np.zeros(n, dtype=int)
    indxr = np.zeros(n, dtype=int)
    indxc = np.zeros(n, dtype=int)
    irow = icol = 0
    for i in range(n):
        big = 0.0
        for j in range(n):
            if ipiv[j] != 1:
                for k in range(n):
                    if ipiv[k] == 0:
                        if abs(a[j, k]) >= big:
                            big = abs(a[j, k])
                            irow, icol = j, k
                    elif ipiv[k] > 1:
                        raise np.linalg.LinAlgError("mc_gaussj: singular matrix")
        ipiv[icol] += 1
        if irow != icol:
            a[[irow, icol], :] = a[[icol, irow], :]
            b[[irow, icol]] = b[[icol, irow]]
        indxr[i] = irow
        indxc[i] = icol
        if a[icol, icol] == 0.0:
            raise np.linalg.LinAlgError("mc_gaussj: singular matrix")
        pivinv = 1. / a[icol, icol]
        a[icol, icol] = 1.
        a[icol, :] = a[icol, :] * pivinv
        b[icol] = b[icol] * pivinv
        for ll in range(n):
            if ll != icol:
                dum = a[ll, icol]
                a[ll, icol] = 0.0
                a[ll, :] = a[ll, :] - a[icol, :] * dum
                b[ll] = b[ll] - b[icol] * dum
    for l in range(n - 1, -1, -1):
        if indxr[l] != indxc[l]:
            a[:, [indxr[l], indxc[l]]] = a[:, [indxc[l], indxr[l]]]
    return b, a


def poly_fit1d(x, y, order, yerr=None, goodbad=None, want_covar=False):
    """
    Spextool POLY_FIT1D, the GAUSSJ branch every caller here uses.

    Order 0 is MC_MEANCOMB of the values with no weights (TOTAL(/DOUBLE)/n
    and IDL's one-pass variance, both returned as FLOAT). FINDGEN
    exponents: pow(). As in the IDL, GOODBAD indexes the NaN-filtered values with positions in the full
    arrays, and the errors are not filtered by it.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    yerr = np.ones_like(x) if yerr is None else np.asarray(yerr, dtype=float)
    z = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr)
    xx, yy, ee = x[z], y[z], yerr[z]
    if goodbad is not None:
        g = np.where(np.asarray(goodbad) == 1)[0]
        if g.size:
            g = np.minimum(g, xx.size - 1)              # IDL clips subscript arrays
            xx, yy = xx[g], yy[g]
            ee = ee[:xx.size]
    if order == 0:
        # MC_MEANCOMB returns its mean and variance as FLOAT: [FLOAT, 0.]
        cnt = yy.size
        mean = itotal(yy) / cnt
        var = (itotal(sq(yy)) - cnt * sq(mean)) / (cnt - 1) / cnt if cnt > 1 else np.nan
        coeffs = np.array([mean, 0.0], dtype=np.float32)
        covar = np.array([var, 0.0], dtype=np.float32)
        return (coeffs, covar) if want_covar else coeffs
    n = order + 1
    alpha = np.zeros((n, n))
    beta = np.zeros(n)
    b = yy / ee
    for i in range(n):
        for j in range(i, n):
            at = dpow(xx, i) / ee
            aj = dpow(xx, j) / ee
            alpha[j, i] = itotal(at * aj)
            beta[i] = itotal(at * b)
    alpha = alpha + alpha.T
    alpha[np.diag_indices(n)] /= 2.
    coeffs, covar = _gaussj(alpha, beta)
    return (coeffs, covar) if want_covar else coeffs


def moments(data, robust=None, igoodbad=None):
    """
    Spextool MC_MOMENTS: mean, variance, standard deviation, and goodbad.

    The moments are IDL MOMENT's (variance with n-1, and its corrected
    two-pass formula), in the precision of the data, as IDL computes them.
    With `robust`, points more than that many MADs (1.482 x median absolute
    deviation) from the median are dropped first; the median is taken with
    /EVEN.
    """
    P = np.float32 if isf32(data) else np.float64
    tot = ftotal if P is np.float32 else (lambda v: np.float64(itotal(v)))
    d = np.asarray(data).astype(P).ravel()
    gb = np.ones(d.size, dtype=int) if igoodbad is None else np.array(igoodbad, dtype=int).ravel().copy()
    gb[~np.isfinite(d)] = 2
    z = np.where(gb == 1)[0]
    if z.size <= 1:
        return np.nan, np.nan, np.nan, gb
    dz = d[z]
    if robust is not None:
        med = median(dz, even=True)
        mad = P(1.482) * median(np.abs(dz - med), even=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            good = np.where(np.abs((dz - med) / mad) <= robust)[0]
        if good.size == 0:
            good = np.arange(dz.size)
    else:
        good = np.arange(dz.size)
    v = dz[good]
    v = v[np.isfinite(v)]
    nx = v.size
    mean = P(tot(v) / P(nx))
    if nx > 1:
        resid = (v - mean).astype(P)
        var = P((tot(sq(resid)) - sq(tot(resid)) / P(nx)) / P(nx - 1.0))
    else:
        var = P(np.nan)
    out = gb.copy()
    out[z] = 0
    out[z[good]] = 1
    return mean, var, np.sqrt(var), out


def robustpoly1d(x, y, order, thresh, eps, yerr=None, igoodbad=None,
                 want_covar=False):
    """
    Spextool ROBUSTPOLY1D: an iteratively sigma-clipped polynomial fit.

    Returns the coefficients and OGOODBAD (1 good, 0 clipped, 2 NaN), and
    the covariance of the last fit when asked.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    yerr = np.ones_like(x) if yerr is None else np.asarray(yerr, dtype=float)
    if igoodbad is not None:
        ogoodbad = np.array(igoodbad, dtype=int).copy()
        ogoodbad[~(np.isfinite(y) & np.isfinite(x) & np.isfinite(yerr))] = 2
    else:
        ogoodbad = np.ones(x.size, dtype=int)
        ogoodbad[~(np.isfinite(y) & np.isfinite(x) & np.isfinite(yerr))] = 2
    z = np.where(ogoodbad == 1)[0]
    if z.size == 0 or z.size < order + 1:
        raise ValueError("robustpoly1d: not enough good points")
    xx, yy, ee = x[z], y[z], yerr[z]
    coeff, covar = poly_fit1d(xx, yy, order, yerr=ee, want_covar=True)
    test = yy - poly(xx, coeff)
    mean, _, sig, _ = moments(test)
    good = np.where(np.abs((test - mean) / sig) <= thresh)[0]
    if good.size == 0:
        raise ValueError("robustpoly1d: every point rejected")
    if good.size != z.size:
        for ittr in range(2, 11):
            sig_old = sig
            coeff, covar = poly_fit1d(xx[good], yy[good], order, yerr=ee[good],
                                      want_covar=True)
            test = yy[good] - poly(xx[good], coeff)
            mean, _, sig, _ = moments(test)
            test = yy - poly(xx, coeff)
            if (sig_old - sig) / sig_old < eps:
                break
            good = np.where(np.abs((test - mean) / sig) <= thresh)[0]
            if good.size == 0:
                raise ValueError("robustpoly1d: every point rejected")
    gb = np.zeros(z.size, dtype=int)
    gb[good] = 1
    ogoodbad[z] = gb
    return (coeff, ogoodbad, covar) if want_covar else (coeff, ogoodbad)


def poly_fit2d(x, y, z, xorder, yorder, zerr=None):
    """Spextool POLY_FIT2D."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    zerr = np.ones_like(x) if zerr is None else np.asarray(zerr, dtype=float)
    g = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    xx, yy, zz, ee = x[g], y[g], z[g], zerr[g]
    rows = []
    for j in range(yorder + 1):
        for i in range(xorder + 1):
            rows.append(ipow(xx, i) * ipow(yy, j) / ee)      # ORDER is INTARR
    AT = np.vstack(rows)
    nc = AT.shape[0]
    alpha = np.array([[itotal(AT[p] * AT[q]) for q in range(nc)] for p in range(nc)])
    beta = np.array([itotal(AT[p] * (zz / ee)) for p in range(nc)])
    coeffs, covar = _gaussj(alpha, beta)
    return coeffs


def _fit_rms(xx, yy, zz, xorder, yorder, coeff):
    """POLY_FIT2D's RMS keyword: STDDEV (MOMENT) of the residuals of the fit."""
    return np.float64(np.sqrt(moments(zz - poly2d(xx, yy, xorder, yorder, coeff))[1]))


def robustpoly2d(x, y, z, xorder, yorder, thresh, eps, zerr=None, want_rms=False):
    """
    Spextool ROBUSTPOLY2D: returns coefficients and OGOODBAD (and, with
    want_rms, the RMS keyword: that of the last POLY_FIT2D call, which after
    the tenth iteration is not the final good set).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    zerr = np.ones_like(x) if zerr is None else np.asarray(zerr, dtype=float)
    ogoodbad = np.ones(x.size, dtype=int)
    ogoodbad[~(np.isfinite(x) & np.isfinite(y) & np.isfinite(z))] = 2
    a = np.where(ogoodbad == 1)[0]
    xx, yy, zz, ee = x[a], y[a], z[a], zerr[a]
    coeff = poly_fit2d(xx, yy, zz, xorder, yorder, zerr=ee)
    rms = _fit_rms(xx, yy, zz, xorder, yorder, coeff) if want_rms else None
    test = zz - poly2d(xx, yy, xorder, yorder, coeff)
    mean, _, sig, _ = moments(test)
    good = np.where(np.abs((test - mean) / sig) <= thresh)[0]
    if good.size != a.size:
        for ittr in range(2, 11):
            sig_old = sig
            coeff = poly_fit2d(xx[good], yy[good], zz[good], xorder, yorder,
                               zerr=ee[good])
            if want_rms:
                rms = _fit_rms(xx[good], yy[good], zz[good], xorder, yorder, coeff)
            test = zz[good] - poly2d(xx[good], yy[good], xorder, yorder, coeff)
            mean, _, sig, _ = moments(test)
            test = zz - poly2d(xx, yy, xorder, yorder, coeff)
            if (sig_old - sig) / sig_old < eps:
                break
            good = np.where(np.abs((test - mean) / sig) <= thresh)[0]
    gb = np.zeros(a.size, dtype=int)
    gb[good] = 1
    ogoodbad[a] = gb
    return (coeff, ogoodbad, rms) if want_rms else (coeff, ogoodbad)


# ---------------------------------------------------------------- combining

def findoutliers(data, thresh):
    """Spextool MC_FINDOUTLIERS: 1 good, 0 outlier, 2 NaN; in the data's precision."""
    P = np.float32 if isf32(data) else np.float64
    d = np.asarray(data).astype(P).ravel()
    mask = np.ones(d.size, dtype=int)
    mask[~np.isfinite(d)] = 2
    med = median(d, even=True)
    mad = np.float32(1.482) * median(np.abs(d - med), even=True)   # FLOAT literal
    if mad == 0.0:
        return mask
    with np.errstate(invalid="ignore", divide="ignore"):
        mask[np.abs((d - med) / mad) > thresh] = 0
    return mask


def _seqsum0(a):
    """IDL TOTAL(a, dim, /DOUBLE) over the stacking dimension: in order, per pixel."""
    acc = np.zeros(a.shape[1:])
    for k in range(a.shape[0]):
        acc = acc + a[k].astype(float)
    return acc


def meancomb(data, datavar=None, mask=None, robust=None, rms=False):
    """
    Spextool MC_MEANCOMB, with the stack along numpy axis 0.

    IDL stacks along its last dimension; here data[k] is the k-th member.
    Returns (mean, variance, ogoodbad). The weighted mean uses `datavar`;
    without it the variance is the sample variance, divided by N unless
    `rms` (the standard error of the mean, as IDL's default). Types as in
    the IDL: 1/datavar and data/datavar in the data's precision, summed in
    double in order, the count a FLOAT TOTAL; results FLOAT.
    """
    P = np.float32 if isf32(data) else np.float64
    data = np.asarray(data).astype(P)
    V = None
    if datavar is not None:
        V = np.float32 if isf32(datavar) else np.float64
        datavar = np.asarray(datavar).astype(V)
    imask = np.ones(data.shape, dtype=int) if mask is None else np.array(mask, dtype=int).copy()
    with np.errstate(divide="ignore", invalid="ignore"):
        if data.ndim == 1:
            if robust is not None:
                r = findoutliers(data, robust)
                imask[r == 0] = 0
            z = imask == 1
            if datavar is not None:
                mvar = 1.0 / itotal((V(1.) / datavar[z]))
                mean = itotal(data[z] / datavar[z]) * mvar
            else:
                cnt = int(z.sum())
                mean = itotal(data[z]) / cnt
                mvar = (itotal(sq(data[z].astype(float))) - cnt * sq(mean)) / (cnt - 1)
                if not rms:
                    mvar = mvar / cnt
            return np.float32(mean), np.float32(mvar), imask
        if robust is not None:
            if data.ndim == 2:
                for i in range(data.shape[1]):
                    r = findoutliers(data[:, i], robust)
                    imask[r == 0, i] = 0
            else:
                med = median(data, even=True, axis=0)
                mad = np.float32(1.482) * median(np.abs(data - med), even=True, axis=0)
                imask = (~(np.abs((data - med) / mad) > robust)).astype(int)
        if datavar is not None:
            tmp = np.where(imask == 0, 0, V(1.) / datavar)
            mvar = 1.0 / _seqsum0(tmp)
            tmp = np.where(imask == 0, 0, data / datavar)
            mean = _seqsum0(tmp) * mvar
        else:
            idata = np.where(imask == 0, P(0), data)
            ndat = _seqsum0(imask).astype(np.float32)          # TOTAL of INT: FLOAT
            mean = _seqsum0(idata) / ndat
            mvar = (_seqsum0(sq(idata.astype(float))) - ndat * sq(mean)) / (ndat - np.float32(1))
            if not rms:
                mvar = mvar / ndat
    return mean.astype(np.float32), mvar.astype(np.float32), imask


def medcomb(data, mask=None, mad=False):
    """
    Spextool MC_MEDCOMB along numpy axis 0 (median and its variance), in the
    data's precision (1.4826 is a FLOAT literal); with a mask, zero
    variances become 1, as in the IDL.
    """
    P = np.float32 if isf32(data) else np.float64
    data = np.asarray(data).astype(P)
    n = data.shape[0]
    k = np.float32(1.4826)
    if mask is None:
        med = median(data, even=True, axis=0)
        var = sq(k * median(np.abs(data - med), even=True, axis=0))
        if not mad:
            var = var / np.float32(n)
        return med, var
    m = np.asarray(mask) == 1
    d = np.where(m, data, P(np.nan))
    med = median(d, even=True, axis=0)
    var = sq(k * median(np.abs(d - med), even=True, axis=0))
    if not mad:
        var = var / np.sum(m, axis=0).astype(np.float32)
    var = np.where(var == 0, P(1.), var).astype(var.dtype)
    return med, var


# ---------------------------------------------------------------- fitting

def _gauss_funct(x, a):
    """
    IDL GAUSS_FUNCT and its partial derivatives. The model is computed from
    the double x GAUSSFIT passes; PDER is an FLTARR, each column rounded
    to float as it is stored (column 2 is made from the stored column 1).
    """
    n = len(a)
    x = np.asarray(x)
    if a[2] != 0.0:
        z = (x - a[1]) / a[2]
        ez = np.exp(-sq(z) / 2.)
    else:
        z = np.full(x.size, 100.0)
        ez = z * 0
    f = a[0] * ez
    if n >= 4:
        f = f + a[3]
    if n >= 5:
        f = f + a[4] * x
    if n >= 6:
        f = f + a[5] * sq(x)
    pder = np.zeros((x.size, n), dtype=np.float32)
    pder[:, 0] = ez
    if a[2] != 0.0:
        pder[:, 1] = a[0] * ez * z / a[2]
    pder[:, 2] = pder[:, 1] * z
    if n > 3:
        pder[:, 3] = 1.
    if n > 4:
        pder[:, 4] = x
    if n > 5:
        pder[:, 5] = sq(x)
    return f, pder


def _seqdot(u, v, P):
    """
    One element of an IDL matrix product: each product rounded to the result
    type, then a running sum in that type (cumsum adds strictly in order).
    """
    prod = np.asarray(u).astype(P) * np.asarray(v).astype(P)
    return np.cumsum(prod, dtype=P)[-1] if prod.size else P(0)


def _seqmat(a, b, P):
    """IDL TRANSPOSE(a) # b for column stacks a [n, p], b [n, q]: running sums, P."""
    prod = np.asarray(a).astype(P)[:, :, None] * np.asarray(b).astype(P)[:, None, :]
    return np.cumsum(prod, axis=0, dtype=P)[-1]


def curvefit(x, y, a, funct=_gauss_funct, weights=None, itmax=20, tol=None,
             single=False):
    """
    IDL CURVEFIT (gradient-expansion least squares, Bevington 1969), step by step.

    `single` is CURVEFIT called with a FLOAT parameter vector. The matrix
    products are IDL's #: running sums in the result type (beta in double,
    alpha in float from the FLTARR derivatives); the step (array/c) #
    TRANSPOSE(beta) is summed in double and rounded into the parameters.
    """
    P = np.float32 if single else np.float64
    a = np.array(a, dtype=P)
    nterms = a.size
    ny = y.size
    nfree = ny - nterms
    tol = (FL(1e-3) if single else 1e-3) if tol is None else tol
    flambda = P(0.001)
    w = np.ones(ny, dtype=np.float32) if weights is None else np.asarray(weights)
    y = np.asarray(y, dtype=float)
    for _ in range(itmax):
        yfit, pder = funct(x, a)
        r = (y - yfit) * w
        beta = _seqmat(r[:, None], pder, np.float64)[0]
        PA = np.result_type(pder.dtype, w.dtype).type
        wp = (w[:, None] * np.ones(nterms, dtype=np.float32)[None, :] * pder).astype(PA)
        alpha = _seqmat(pder, wp, PA)
        chisq1 = itotal(w * sq(y - yfit)) / nfree
        if chisq1 < itotal(np.abs(y)) / 1e7 / nfree:
            return a.astype(float), yfit
        cd = np.sqrt(np.diag(alpha)).astype(PA)
        c = np.outer(cd, cd).astype(PA)
        count = 0
        yfit1 = yfit
        while True:
            count += 1
            arr = (alpha / c).astype(PA)
            arr[np.diag_indices(nterms)] = (arr[np.diag_indices(nterms)]
                                            * (PA(1.) + flambda)).astype(PA)
            arr = invert(arr)[0]
            m = (arr / c).astype(PA)
            step = np.array([_seqdot(m[:, j], beta, np.float64) for j in range(nterms)])
            b = (a + step).astype(P)
            yfitb, _ = funct(x, b)
            chisq = itotal(w * sq(y - yfitb)) / nfree
            if not np.isfinite(chisq) or (count > 30 and chisq >= chisq1):
                return a.astype(float), yfit1
            flambda = P(flambda * P(10.))
            if chisq <= chisq1:
                break
        flambda = P(flambda / P(100.))
        a = b
        yfit = yfitb
        if (chisq1 - chisq) / chisq1 <= tol:
            break
    return a.astype(float), yfit


def gaussfit(x, y, nterms=6, estimates=None):
    """
    IDL GAUSSFIT: returns (yfit, a), with IDL's own starting estimates.

    The result is cast to float32 when neither input is double, as GAUSSFIT
    does, so a centroid carries the same precision IDL reports.
    """
    is_double = isf64(x) or isf64(y)
    xd = np.asarray(x, dtype=float)
    yd0 = np.asarray(y, dtype=float)
    n = yd0.size
    a = None
    if estimates is None or estimates[2] == 0:
        if nterms > 3:
            deg = 0 if nterms == 4 else 1
            c = poly_fit(xd, yd0, deg)
            yf = c[deg]
            for k in range(deg - 1, -1, -1):
                yf = c[k] + yf * xd
            yd = yd0 - yf
        else:
            yd = yd0
            c = np.array([0.0])
        imax = int(np.argmax(yd))
        imin = int(np.argmin(yd))
        i0 = imax if abs(yd[imax]) > abs(yd[imin]) else imin
        i0 = min(max(i0, 1), n - 2)
        dy = yd[i0]
        dl = dy / np.e
        i = 0
        while (i0 + i + 1) < n and (i0 - i) > 0 and abs(yd[i0 + i]) > abs(dl) \
                and abs(yd[i0 - i]) > abs(dl):
            i += 1
        a = [yd[i0], xd[i0], abs(xd[i0] - xd[i0 + i])]
        if nterms > 3:
            a.append(c[0])
        if nterms > 4:
            a.append(c[1] if c.size > 1 else 0.0)
        if nterms > 5:
            a.append(0.0)
        a = np.array(a, dtype=float)
    single = False
    if estimates is not None:
        single = not isf64(estimates)
        tmp = np.array(estimates, dtype=float)
        if estimates[2] == 0:
            tmp[2] = a[2]
        a = tmp
    a, yfit = curvefit(xd, yd0, a, single=single)
    if not is_double:
        a = a.astype(np.float32).astype(float)
        yfit = yfit.astype(np.float32)
    return yfit, a


def _mpfitpeak_est(x, y):
    """
    MPFITPEAK_EST: the starting guess MPFITPEAK makes for itself, in the
    precision of the data, with IDL's TOTAL (NaN not skipped, as MPFITPEAK
    is called without /NAN).
    """
    P = np.float64 if (isf64(x) or isf64(y)) else np.float32
    tot = (lambda v: P(itotal(v))) if P is np.float64 else ftotal
    x = np.asarray(x).astype(P)
    y = np.asarray(y).astype(P)
    nx = x.size
    s = np.argsort(x, kind="stable")
    xs, ys = x[s], y[s]
    maxx, minx = xs.max(), xs.min()
    maxy, miny = ys.max(), ys.min()
    dx = (P(0.5) * np.concatenate([[xs[1] - xs[0]], xs[2:] - xs[:-2],
                                   [xs[nx - 1] - xs[nx - 2]]])).astype(P)
    totarea = tot(dx * ys)
    av = P(totarea / (maxx - minx))
    if miny == maxy:
        return np.array([0.0, xs[nx // 2], (xs[nx - 1] - xs[0]) / 2, ys[0]]).astype(P)
    wh1 = np.where(y >= av)[0]
    wh2 = np.where(y <= av)[0]
    sd1 = tot(sq(x[wh1])) / P(wh1.size) - sq(tot(x[wh1]) / P(wh1.size))
    sd2 = tot(sq(x[wh2])) / P(wh2.size) - sq(tot(x[wh2]) / P(wh2.size))
    if sd1 < sd2:
        cent = x[y == maxy][0]
        peak = P(maxy - av)
    else:
        cent = x[y == miny][0]
        peak = P(miny - av)
    peakarea = P(totarea - tot(dx * np.minimum(ys, av)))
    if peak == 0:
        peak = P(0.5 * peakarea)
    width = P(peakarea / (2 * abs(peak)))
    if width == 0 or not np.isfinite(width):
        width = median(dx)
    return np.array([peak, cent, width, av], dtype=P)


def _peak_model(p, x, lorentz=False):
    if lorentz:
        u = sq(((x - p[1]) / p[2]))
        f = p[0] / (u + 1.0)
    else:
        u = sq(((x - p[1]) / p[2]))
        f = p[0] * np.exp(-0.5 * u)
    if p.size > 3:
        f = f + p[3]
    if p.size > 4:
        f = f + p[4] * x
    if p.size > 5:
        f = f + p[5] * sq(x)
    return f


def mpfitpeak(x, y, nterms=4, estimates=None, lorentzian=False, weights=None):
    """
    MPFITPEAK: a Gaussian (or Lorentzian) peak on a constant or a line.

    Without ESTIMATES the starting point is MPFITPEAK_EST's. The fit is
    MPFIT itself (triplespec_reducer.mpfit), in the precision IDL would
    use for these inputs. Returns (yfit, a) with a[2] positive.
    """
    from triplespec_reducer.mpfit import fit_peak
    if estimates is None:
        estimates = _mpfitpeak_est(x, y)
    yfit, a, _ = fit_peak(x, y, estimates, nterms=nterms, lorentzian=lorentzian,
                          weights=weights)
    return yfit, a


def gaussian(x, p):
    """
    astrolib GAUSSIAN(x, p) without /DOUBLE: p = [peak, centre, sigma(, const)].

    The Gaussian is stored in an FLTARR: exp() of the (possibly double)
    argument rounded to float, zero where it would underflow a float; then
    p[0] times it, in the parameters' precision.
    """
    xa = np.asarray(x)
    pa = np.asarray(p)
    if pa.dtype.kind in "iu":
        pa = pa.astype(np.float32)
    P = np.float64 if (isf64(xa) or isf64(pa)) else np.float32
    xa = xa.astype(P)
    z = ((xa - pa[1].astype(P)) / pa[2].astype(P)).astype(P)
    zz = (z * z).astype(P)
    minexp = alog(np.float32(np.finfo(np.float32).tiny))
    g = np.zeros(xa.shape, dtype=np.float32)
    w = zz < -2 * minexp
    g[w] = exp((-zz[w] / 2).astype(P)).astype(np.float32)
    out = pa[0] * g
    return out + pa[3] if pa.size >= 4 else out


# ---------------------------------------------------------------- signals

def convol(a, kernel, scale=None, edge_truncate=False):
    """
    IDL CONVOL(a, kernel, scale, /EDGE_TRUNCATE) for 1-D arrays, /CENTER; scale and the keyword optional.

    R[t] = sum_i a[t + i - K/2] * k[i] / scale: a correlation, the kernel
    is not reversed. Without /EDGE_TRUNCATE the result is zero within K/2
    of either end; with it, the array is extended by repeating its end
    values. As IDL computes it: a running sum over the kernel in the
    precision of the array. Kernel elements equal to zero are skipped, as
    IDL skips them (probed in IDL 7.1: with [0, 1, 0] a NaN next to a pixel
    does not reach it).
    """
    P = np.float32 if isf32(a) else np.float64
    a = np.asarray(a).astype(P)
    k = np.asarray(kernel).astype(P)
    n, nk = a.size, k.size
    h = nk // 2
    out = np.zeros(n, dtype=P)
    if edge_truncate:
        pad = np.concatenate([np.full(h, a[0]), a, np.full(nk - 1 - h, a[-1])]).astype(P)
        acc = np.zeros(n, dtype=P)
        for j in range(nk):
            if k[j] == 0:
                continue
            acc = (acc + pad[j:j + n] * k[j]).astype(P)
        out = acc
    elif n >= nk:
        m = n - nk + 1
        acc = np.zeros(m, dtype=P)
        for j in range(nk):
            if k[j] == 0:
                continue
            acc = (acc + a[j:j + m] * k[j]).astype(P)
        out[h:h + m] = acc
    if scale is not None:
        out = (out / P(scale)).astype(P)
    return out


def psf_gaussian_1d(npixel, fwhm, normalize=True):
    """
    astrolib PSF_GAUSSIAN(NDIMEN=1, NPIXEL=, FWHM=[, /NORMALIZE]): the
    standard deviation computed in double and rounded to float, the PSF
    astrolib GAUSSIAN of FINDGEN - (npixel-1)/2., normalized by its TOTAL.
    """
    st_dev = np.float32(np.float64(fwhm) / (2.0 * np.sqrt(2.0 * np.log(2.0))))
    x = np.arange(npixel, dtype=np.float32) - np.float32((npixel - 1) / 2.)
    psf = gaussian(x, np.array([1, 0, st_dev], dtype=np.float32))
    return (psf / ftotal(psf)).astype(np.float32) if normalize else psf


def c_correlate(x, y, lag):
    """
    IDL C_CORRELATE: in single precision unless x or y is double, the sums
    being IDL TOTALs, exactly as the IDL routine computes them.
    """
    double = isf64(x) or isf64(y)
    P = np.float64 if double else np.float32
    tot = (lambda v: np.float64(itotal(v))) if double else ftotal
    x = np.asarray(x).astype(P)
    y = np.asarray(y).astype(P)
    nx = x.size
    xd = (x - P(tot(x) / P(nx))).astype(P)
    yd = (y - P(tot(y) / P(nx))).astype(P)
    out = np.empty(len(lag), dtype=P)
    for k, L in enumerate(lag):
        L = int(L)
        if L >= 0:
            out[k] = tot(xd[0:nx - L] * yd[L:])
        else:
            out[k] = tot(yd[0:nx + L] * xd[-L:])
    return (out / P(np.sqrt(P(tot(sq(xd)) * tot(sq(yd)))))).astype(P)


def idl_fft(a, inverse=False):
    """IDL FFT: the forward transform carries the 1/N."""
    a = np.asarray(a)
    return np.fft.ifft(a) * a.size if inverse else np.fft.fft(a) / a.size


def spline(x, y, t, sigma=1.0):
    """
    IDL SPLINE(x, y, t, sigma): Cline's spline under tension.

    Step by step in IDL's types: FLOAT data unless x, y or t is
    double; the 0.5d factors make the tridiagonal set-up double; float EXP
    is the C library's; the tension keeps its own type (a double tension
    makes the dependent quantities double).
    """
    dbl = isf64(x) or isf64(y) or isf64(t)
    P = np.float64 if dbl else np.float32
    D = np.float64
    # A plain Python tension is an IDL literal (FLOAT for float data); pass a
    # numpy float64 for a DOUBLE one.
    sig = sigma if isinstance(sigma, np.generic) else P(sigma)
    if sig < P(0.001):
        sig = np.result_type(sig, P).type(0.001)
    xx = np.asarray(x).astype(P).ravel()
    yy = np.asarray(y).astype(P).ravel()
    tt = np.atleast_1d(np.asarray(t).astype(P))
    n = min(xx.size, yy.size)
    xx, yy = xx[:n], yy[:n]
    yp = np.zeros(2 * n, dtype=P)
    delx1 = xx[1] - xx[0]
    dx1 = (yy[1] - yy[0]) / delx1
    nm1 = n - 1
    delx2 = xx[2] - xx[1]
    delx12 = xx[2] - xx[0]
    c1 = -(delx12 + delx1) / delx12 / delx1
    c2 = delx12 / delx1 / delx2
    c3 = -delx1 / delx12 / delx2
    slpp1 = c1 * yy[0] + c2 * yy[1] + c3 * yy[2]
    deln = xx[nm1] - xx[nm1 - 1]
    delnm1 = xx[nm1 - 1] - xx[nm1 - 2]
    delnn = xx[nm1] - xx[nm1 - 2]
    c1 = (delnn + deln) / delnn / deln
    c2 = -delnn / deln / delnm1
    c3 = deln / delnn / delnm1
    slppn = c3 * yy[nm1 - 2] + c2 * yy[nm1 - 1] + c1 * yy[nm1]
    sigmap = sig * nm1 / (xx[nm1] - xx[0])
    dels = sigmap * delx1
    exps = exp(dels)
    sinhs = D(0.5) * (exps - 1. / exps)
    sinhin = 1. / (delx1 * sinhs)
    diag1 = sinhin * (dels * D(0.5) * (exps + 1. / exps) - sinhs)
    diagin = 1. / diag1
    yp[0] = diagin * (dx1 - slpp1)
    spdiag = sinhin * (sinhs - dels)
    yp[n] = diagin * spdiag
    delx2 = xx[1:] - xx[:-1]
    dx2 = (yy[1:] - yy[:-1]) / delx2
    dels = sigmap * delx2
    exps = exp(dels)
    sinhs = D(0.5) * (exps - 1. / exps)
    sinhin = 1. / (delx2 * sinhs)
    diag2 = sinhin * (dels * (D(0.5) * (exps + 1. / exps)) - sinhs)
    diag2 = np.concatenate([[0.0], diag2[:-1] + diag2[1:]])
    dx2nm1 = dx2[nm1 - 1]
    dx2 = np.concatenate([np.zeros(1, P), dx2[1:] - dx2[:-1]])
    spdiag = sinhin * (sinhs - dels)
    for i in range(1, nm1):
        diagin = 1. / (diag2[i] - spdiag[i - 1] * yp[i + n - 1])
        yp[i] = diagin * (dx2[i] - spdiag[i - 1] * yp[i - 1])
        yp[i + n] = diagin * spdiag[i]
    diagin = 1. / (diag1 - spdiag[nm1 - 1] * yp[n + nm1 - 1])
    yp[nm1] = diagin * (slppn - dx2nm1 - spdiag[nm1 - 1] * yp[nm1 - 1])
    for i in range(n - 2, -1, -1):
        yp[i] = yp[i] - yp[i + n] * yp[i + 1]
    m = tt.size
    subs = np.full(m, nm1, dtype=int)
    s_ = xx[nm1] - xx[0]
    sigmap = sig * nm1 / s_
    j = 0
    done = False
    for i in range(1, nm1 + 1):
        while tt[j] < xx[i]:
            subs[j] = i
            j += 1
            if j == m:
                done = True
                break
        if done:
            break
    subs1 = subs - 1
    del1 = tt - xx[subs1]
    del2 = xx[subs] - tt
    dels = xx[subs] - xx[subs1]
    exps1 = exp(sigmap * del1)
    sinhd1 = 0.5 * (exps1 - 1. / exps1)
    exps = exp(sigmap * del2)
    sinhd2 = 0.5 * (exps - 1. / exps)
    exps = exps1 * exps
    sinhs = 0.5 * (exps - 1. / exps)
    spl = (yp[subs] * sinhd1 + yp[subs1] * sinhd2) / sinhs + \
        ((yy[subs] - yp[subs]) * del1 + (yy[subs1] - yp[subs1]) * del2) / dels
    return spl


def qsimp(func, a, b, eps=1e-6, jmax=20):
    """astrolib QSIMP: Simpson's rule refined until the change is below eps."""
    ost = -1.0e30
    os_ = -1.0e30
    st = 0.0
    for j in range(1, jmax + 1):
        if j == 1:
            st = 0.5 * (b - a) * (func(a) + func(b))
            it = 1
        else:
            tnm = it
            dl = (b - a) / tnm
            xs = a + 0.5 * dl + dl * np.arange(it)
            st = 0.5 * (st + (b - a) * np.sum(func(xs)) / tnm)
            it *= 2
        s = (4.0 * st - ost) / 3.0
        if j > 5 and abs(s - os_) < eps * abs(os_):
            return s
        os_ = s
        ost = st
    return s


# ---------------------------------------------------------------- geometry

def triangulate_grid(xs, ys, values, xgrid, ygrid, missing=np.nan, tri=None):
    """
    IDL TRIANGULATE + TRIGRID (linear) onto a regular grid.

    `values` may be a list of arrays sharing one triangulation; returns one
    [len(ygrid), len(xgrid)] array per input. IDL's TRIANGULATE is a
    Delaunay triangulation and TRIGRID interpolates linearly within each
    triangle, which is what this does.
    """
    from scipy.interpolate import LinearNDInterpolator
    pts = np.column_stack([np.asarray(xs, float), np.asarray(ys, float)])
    tri = Delaunay(pts) if tri is None else tri
    gx, gy = np.meshgrid(xgrid, ygrid)
    out = []
    for v in values:
        f = LinearNDInterpolator(tri, np.asarray(v, float), fill_value=missing)
        out.append(f(gx, gy))
    return out, tri


def trigrid_points(xs, ys, values, xo, yo, missing=np.nan, tri=None):
    """IDL GRIDDATA /LINEAR with TRIANGLES (or TRIGRID XOUT/YOUT) at points."""
    from scipy.interpolate import LinearNDInterpolator
    pts = np.column_stack([np.asarray(xs, float), np.asarray(ys, float)])
    tri = Delaunay(pts) if tri is None else tri
    out = []
    for v in values:
        f = LinearNDInterpolator(tri, np.asarray(v, float), fill_value=missing)
        out.append(f(np.asarray(xo, float), np.asarray(yo, float)))
    return out, tri


def trigrid_axes(xmin, xmax, dx, ymin, ymax, dy):
    """
    The output grid of IDL TRIGRID(GS=[dx,dy], LIMITS=[xmin,ymin,xmax,ymax]).

    TRIGRID makes ceil((xmax-xmin)/dx) + 1 points from xmin, the ratio
    taken in single precision: a range that is 5.000000000005 steps long
    gives 6 points, not 7, and one 53.5 steps long gives 55, the last
    beyond xmax (both checked against IDL 7.1).
    """
    # LIMITS is one array: DOUBLE when any limit is (the wavelengths are),
    # so float spatial limits are differenced in double.
    xmin, xmax, ymin, ymax = (float(v) for v in (xmin, xmax, ymin, ymax))
    nx = int(np.ceil(np.float32((xmax - xmin) / dx))) + 1
    ny = int(np.ceil(np.float32((ymax - ymin) / dy))) + 1
    return xmin + np.arange(nx) * dx, ymin + np.arange(ny) * dy


def polyfillaa(px, py, sx, sy):
    """
    JD Smith's POLYFILLAA: pixels a polygon covers, with the area of each.

    Pixel (i, j) spans [i, i+1] x [j, j+1] in polygon coordinates; indices
    are returned as (column, row) pairs with their overlap areas (FLOAT).
    The compiled clipper Spextool builds is what IDL runs: the same as
    extract._clip_area, columns outermost, zero areas dropped.
    """
    from triplespec_reducer.extract import _clip_area
    pxd = np.asarray(px, float)
    pyd = np.asarray(py, float)
    i0 = max(int(np.floor(pxd.min())), 0)
    i1 = min(int(np.floor(pxd.max())), sx - 1)
    j0 = max(int(np.floor(pyd.min())), 0)
    j1 = min(int(np.floor(pyd.max())), sy - 1)
    px32, py32 = pxd.astype(np.float32), pyd.astype(np.float32)
    cols, rows, areas = [], [], []
    for i in range(i0, i1 + 1):
        for j in range(j0, j1 + 1):
            a = _clip_area(px32, py32, i, j)
            if a != 0.0:
                cols.append(i)
                rows.append(j)
                areas.append(a)
    return np.array(cols, int), np.array(rows, int), np.array(areas, np.float32)


# ---------------------------------------------------------------- FITS headers

def _fxtext(value, fmt=None):
    """The text FXADDPAR (or SXADDPAR) writes for a number."""
    if fmt is not None:
        f = fmt.upper()
        if f.startswith("G0."):                       # minimal width, n significant
            return f"{float(value):.{int(f[3:])}G}"
        if f.startswith("F"):
            w, d = (int(t) for t in f[1:].split("."))
            s = f"{float(value):{w}.{d}f}"
            return "*" * w if len(s) > w else s
        if f.startswith("I"):
            w = int(f[1:].split(".")[0])
            s = str(int(np.round(float(value))))
            return "*" * w if len(s) > w else s
        raise ValueError(f"no IDL format {fmt}")
    if isinstance(value, (bool, np.bool_)):
        raise TypeError("logicals are written as strings")
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isf32(value):
        return f"{float(value):#.6G}"                 # FLOAT: STRTRIM(STRING(v)), G13.6
    return f"{float(value):#.12G}"                    # DOUBLE: G19.12


def fxcard(key, value, comment="", fmt=None):
    """
    The 80-character record FXADDPAR/SXADDPAR makes. A number is
    right-justified to column 30 (or starts at 10 when longer), FLOAT values
    (numpy float32) in IDL's default format, DOUBLE (any other real) as
    G19.12, or with `fmt`; a string is quoted from column 10, padded to eight
    characters; 'T'/'F' sit in column 30. Then ' /' at 30, the comment at 32.
    None when FXADDPAR would continue the string on CONTINUE cards.
    """
    from astropy.io import fits
    h = list(f"{key.upper():<8s}= ".ljust(80))
    if isinstance(value, (bool, np.bool_)):
        value = "T" if value else "F"
    if isinstance(value, str):
        if value.upper() in ("T", "F"):
            h[29] = value.upper()
        else:
            val = value.replace("'", "''")
            if len(val) > 68:
                # FXADDPAR_CONTPAR: 67 characters and '&' per record, the
                # rest on CONTINUE records; the comment follows the first
                # record's value and is cut at column 80
                recs = [("".join(h[:10]) + "'" + val[:67] + "&'" + " /" + comment)[:80].ljust(80)]
                rest = val[67:]
                while rest:
                    if len(rest) > 68:
                        recs.append(("CONTINUE  '" + rest[:67] + "&'").ljust(80))
                        rest = rest[67:]
                    else:
                        recs.append(("CONTINUE  '" + rest + "'").ljust(80))
                        rest = ""
                return fits.Card.fromstring("".join(recs))
            if len(val) > 18:
                img = "".join(h[:10]) + "'" + val + "'" + " /" + comment
                return fits.Card.fromstring(img[:80].ljust(80))
            h[10:11 + len(val)] = "'" + val
            q = 11 + max(len(val), 8)
            h[q:q + 1] = "'"
    else:
        v = _fxtext(value, fmt)
        st = max(30 - len(v), 10)
        h[st:st + len(v)] = v
    h[30:32] = " /"
    h[32:32 + len(comment)] = comment
    return fits.Card.fromstring("".join(h[:80]))


def fxaddpar(header, key, value, comment="", fmt=None):
    """FXADDPAR on an astropy header: the record as IDL writes it, replaced in place."""
    if value is not None and not isinstance(value, (str, bool, np.bool_)) \
            and set(_fxtext(value, fmt)) == {"*"}:    # a format overflow: FXPAR reads a string
        value = _fxtext(value, fmt)
    card = None if value is None else fxcard(key, value, comment, fmt)
    if card is None:
        header[key] = (value, comment)
        return
    if key in header:
        i = header.index(key)
        del header[i]
        header.insert(i, card)
    else:
        # a new keyword goes before the first blank, COMMENT or HISTORY
        # record (FXPARPOS BEFORE=''/'COMMENT'/'HISTORY'), else at the end
        pos = [n for n, c in enumerate(header.cards) if c.keyword in ("", "COMMENT", "HISTORY")]
        if pos:
            header.insert(pos[0], card)
        else:
            header.append(card, useblanks=False, bottom=True)
    if len(card.image) > 80 and "LONGSTRN" not in header:
        # FXADDPAR_CONTWARN: LONGSTRN and four COMMENTs before the keyword
        from astropy.io import fits
        at = header.index(key)
        warn = [fxcard("LONGSTRN", "OGIP 1.0", " The OGIP long string convention may be used.")]
        for t in (" This FITS file may contain long string keyword values that are",
                  " continued over multiple keywords.  This convention uses the  '&'",
                  " character at the end of a string which is then continued",
                  " on subsequent keywords whose name = 'CONTINUE'."):
            warn.append(fits.Card.fromstring(("COMMENT " + t).ljust(80)))
        for n, c in enumerate(warn):
            header.insert(at + n, c)


def fxpar(header, key):
    """
    FXPAR: a keyword's value typed as IDL types it. A real is DOUBLE when it
    holds a 'D' or its text up to the '/' (trailing blanks included) is 8 or
    more characters long, FLOAT otherwise; a whole number is LONG.
    """
    card = header.cards[key]
    sv = card.image[9:80].strip()
    if sv.startswith("'") or isinstance(card.value, (str, bool)):
        return card.value
    slash = sv.find("/")
    test = sv[:slash] if slash > 0 else sv
    tok = test.split()
    if not tok:
        return card.value
    if tok[0] in ("T", "F"):
        return int(tok[0] == "T")
    if any(c in test for c in ".ED"):
        if test.find("D") > 0 or len(test) >= 8:
            return np.float64(float(test.replace("D", "E")))
        return np.float32(float(test))
    n = int(test)
    return n if -2 ** 31 <= n <= 2 ** 31 - 1 else np.float64(n)


def fxcomment(header, key):
    """FXPAR's COMMENT=: everything after the '/' that ends the value, as written."""
    sv = header.cards[key].image[9:80].strip()
    start = 0
    if sv.startswith("'"):
        k = 1
        while True:
            k = sv.find("'", k)
            if k < 0 or sv[k + 1:k + 2] != "'":
                break
            k += 2
        start = max(k, 0)
    slash = sv.find("/", start)
    return "" if slash < 0 else sv[slash + 1:]


def splittext(text, length):
    """MC_SPLITTEXT: break at the last blank before `length`, as Spextool wraps history."""
    if len(text) < length:
        return [text]
    out, ntext, left = [], text, len(text)
    while len(ntext) > length:
        pos = ntext.rfind(" ", 0, length)            # STRPOS(..., length-1, /REVERSE_SEARCH)
        if pos == -1:
            pos = length - 1
        out.append(ntext[:pos])
        ntext = ntext[pos + 1:]
        if len(ntext) == left:
            break
        left = len(ntext)
    return out + [ntext]


def addhistlabel(header, label):
    """
    MC_ADDHISTLABEL: a blank line, '##...# label #...##', a blank line. IDL
    writes the label line from column 2, which is not a valid FITS record;
    here it is a standard blank-keyword card holding the same text.
    """
    tot = (80 - (len(label) + 2)) // 2
    # SXADDHIST /BLANK puts ' '+text in an 80-character record: 79 characters kept
    for t in (" ", ("#" * tot + " " + label + " " + "#" * tot)[:79], " "):
        header.append(("", t), useblanks=False, bottom=True)   # SXADDHIST: before END


def addhistory(header, line):
    """SXADDHIST (or FXADDPAR 'HISTORY') of one line: a HISTORY record just before END."""
    header.append(("HISTORY", line), useblanks=False, bottom=True)


def sxaddhist(header, text, length=70):
    """MC_SPLITTEXT then SXADDHIST: the text as HISTORY lines, at the end of the header."""
    for line in splittext(text, length):
        addhistory(header, line)
