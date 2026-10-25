"""
Combining the spectra of one target: XCOMBSPEC.

    every spectrum onto the first one's grid   MC_INTERPSPEC (LEAVENANS)
    optional scale factors from one order      XMC_SCALESPEC 'Median',
                                               'Spectrum' or 'None'
                                               (MC_GETSPECSCALE over the
                                               middle 80% of that order or
                                               a chosen range)
    optional prune of a spectrum or pixels     XMC_MASKSPEC
    optional shape correction                  MC_SPECCOR
    combine, per order and aperture            MC_MEANCOMB / MC_MEDCOMB
    write, header and history                  XCOMBSPEC_WRITEFILE

Once a spectrum is removed from an order, XCOMBSPEC gives the statistics the
spectra in use with the pixel mask of the whole stack; that is reproduced
(meancomb_stackmask, medcomb_stackmask), the intended computation is the
option intended_masks=True.

The statistics are XCOMBSPEC's nine, by name; the default is its own,
"Robust Weighted Mean" with a threshold of 8.
"""

import os

import numpy as np
from astropy.io import fits

from triplespec_reducer import idl as I
from triplespec_reducer.extract import read_spectrum, caaspecflags

STATISTICS = ["Robust Weighted Mean", "Robust Mean (RMS)", "Robust Mean (Std Error)",
              "Weighted Mean", "Mean (RMS)", "Mean (Std Error)", "Median (MAD)",
              "Median (Median Error)", "Sum"]


def interpspec(ix, iy, ox, iyerror=None, leavenans=False):
    """
    MC_INTERPSPEC: linear interpolation that never extrapolates.

    In the inputs' types, as the IDL computes (float spectra after
    MC_CAASPECFLAGS, double output wavelengths), the results stored in
    FLTARRs.
    """
    ix, iy, ox = np.asarray(ix), np.asarray(iy), np.asarray(ox)
    zci = np.isfinite(ix) if leavenans else (np.isfinite(ix) & np.isfinite(iy))
    cix, ciy = ix[zci], iy[zci]
    ciyerr = None if iyerror is None else np.asarray(iyerror)[zci]
    zco = np.isfinite(ox)
    cox = ox[zco]
    # FLTARR outputs: the interpolation is rounded to float as stored.
    coy = np.full(cox.size, np.nan, dtype=np.float32)
    coe = np.full(cox.size, np.nan, dtype=np.float32)
    idx = I.findidx(cix, cox)
    g = np.isfinite(idx)
    if not g.any():
        raise ValueError("mc_interpspec: no overlap")
    idx = idx[g]
    bot = np.floor(idx).astype(int)
    top = np.ceil(idx).astype(int)
    tcox = cox[g]
    tcoy = np.zeros(g.sum(), dtype=np.float32)
    same = bot == top
    tcoy[same] = ciy[bot[same]]
    d = ~same
    m = (ciy[top[d]] - ciy[bot[d]]) / (cix[top[d]] - cix[bot[d]])
    tcoy[d] = ciy[bot[d]] + m * (tcox[d] - cix[bot[d]])
    coy[g] = tcoy
    oy = np.zeros(ox.size, dtype=np.float32)
    oy[zco] = coy
    if ciyerr is None:
        return oy, None
    v = I.sq(ciyerr)
    tv = np.full(g.sum(), np.nan, dtype=np.float32)
    tv[same] = v[bot[same]]
    m1 = I.sq((tcox[d] - cix[bot[d]]) / (cix[top[d]] - cix[bot[d]]))
    v1 = v[bot[d]] + m1 * (v[bot[d]] + v[top[d]])
    m2 = I.sq((cix[top[d]] - tcox[d]) / (cix[top[d]] - cix[bot[d]]))
    v2 = v[top[d]] + m2 * (v[bot[d]] + v[top[d]])
    tv[d] = np.minimum(v1, v2)
    coe[g] = np.sqrt(tv)
    oe = np.zeros(ox.size, dtype=np.float32)
    oe[zco] = coe
    return oy, oe


def interpflagspec(ix, iflag, ox, nbits=8):
    """
    MC_INTERPFLAGSPEC: each bit is interpolated linearly on its own and
    rounded up, so a new pixel carries a bit if either neighbour does.
    """
    ix = np.asarray(ix, dtype=float)
    iy = np.nan_to_num(np.asarray(iflag, dtype=float)).astype(np.uint8)
    ox = np.asarray(ox, dtype=float)
    zci = np.isfinite(ix)
    cix, ciy = ix[zci], iy[zci]
    zco = np.isfinite(ox)
    cox = ox[zco]
    coy = np.zeros(cox.size)
    idx = I.findidx(cix, cox)
    g = np.isfinite(idx)
    idx = idx[g]
    bot = np.floor(idx).astype(int)
    top = np.ceil(idx).astype(int)
    tcox = cox[g]
    for b in range(nbits):
        st = ((ciy >> b) & 1).astype(np.float32)
        t = np.zeros(g.sum(), dtype=np.float32)
        same = bot == top
        t[same] = st[bot[same]]
        d = ~same
        m = (st[top[d]] - st[bot[d]]) / (cix[top[d]] - cix[bot[d]])
        t[d] = st[bot[d]] + m * (tcox[d] - cix[bot[d]])
        coy[np.where(g)[0]] += np.ceil(t) * 2 ** b
    oy = np.zeros(ox.size)
    oy[zco] = coy
    return oy.astype(np.uint8)


def combflagstack(stack):
    """MC_COMBFLAGSTACK: bitwise OR down the stack."""
    stack = np.asarray(stack).astype(np.uint8)
    return np.bitwise_or.reduce(stack, axis=0)


def getspecscale(stack, idx=None):
    """MC_GETSPECSCALE: median of reference/spectrum, for each spectrum."""
    ref = I.medcomb(stack)[0] if idx is None else stack[idx]
    out = np.zeros(stack.shape[0], dtype=np.float32)         # FLTARR
    for i in range(stack.shape[0]):
        with np.errstate(divide="ignore", invalid="ignore"):
            div = ref / stack[i]
        out[i] = I.median(div[np.isfinite(div)], even=True)
    return out


def load(files, combine_apertures=False):
    """
    XCOMBSPEC_LOADSPEC: `spc` holds, per aperture and order, a stack (nspec, 4, npix) on
    the first file's wavelength grid. With `combine_apertures` ('Combine
    apertures: Yes') every aperture of every file goes into one stack and
    one aperture comes out; the integration time then counts each file
    twice, as XCOMBSPEC counts it.

    Each file is interpolated from its pixels where the first file's
    wavelengths are finite, as XCOMBSPEC_LOADSPEC selects them. The
    returned state also keeps the unmodified stacks (`ospc`, from which
    scaling starts), the masks Prune Spectra sets (`spcmask`, `pixmask`,
    all ones) and how far each aperture has been modified (`modspec`:
    1 scaled, 3 shape-corrected), which XCOMBSPEC uses to refuse a second
    scaling or correction.
    """
    first = read_spectrum(files[0])
    orders, naps = first["orders"], first["naps"]
    fs = first["spec"]
    npix = fs.shape[2]
    fnaps = 1 if combine_apertures else naps
    nspec = len(files) * naps if combine_apertures else len(files)
    spc = [[np.full((nspec, 4, npix), np.nan, dtype=np.float32) for _ in orders]
           for _ in range(fnaps)]                         # FLTARR, as XCOMBSPEC
    itot, hdrs = 0, []
    for i, f in enumerate(files):
        d = read_spectrum(f)
        data = d["spec"]
        data = caaspecflags(data)
        t = d["itot"]                         # FXPAR(hdr, 'ITOT'): its text decides the type
        itot = itot + t * 2 if combine_apertures else itot + t
        hdrs.append(d["header"])
        for j in range(len(orders)):
            for k in range(naps):
                n = j * naps + k
                z = np.isfinite(fs[n, 0])
                nf, ne = interpspec(data[n, 0][z], data[n, 1][z], fs[n, 0][z],
                                    iyerror=data[n, 2][z], leavenans=True)
                nb = interpflagspec(data[n, 0][z], data[n, 3][z], fs[n, 0][z])
                a, row = (0, i * naps + k) if combine_apertures else (k, i)
                spc[a][j][row, 0, z] = fs[n, 0][z]
                spc[a][j][row, 1, z] = nf
                spc[a][j][row, 2, z] = ne
                spc[a][j][row, 3, z] = nb
    ospc = [[s_.copy() for s_ in ap] for ap in spc]
    spcmask = [np.ones((len(orders), nspec), dtype=int) for _ in range(fnaps)]
    pixmask = [[np.ones((nspec, npix), dtype=int) for _ in orders] for _ in range(fnaps)]
    return dict(spc=spc, ospc=ospc, orders=orders, naps=fnaps, filenaps=naps, files=list(files),
                itot=itot, headers=hdrs, first=first, combine_apertures=combine_apertures,
                spcmask=spcmask, pixmask=pixmask, modspec=[0] * fnaps)


SCALE_MODES = ["None", "Median", "Spectrum"]


def scale(loaded, order, aperture=1, spcmask=None, srange=None, mode="Median", spectrum=None):
    """
    XCOMBSPEC_SCALESPEC with XMC_SCALESPEC.

    Factors come from `order` and are applied to every order of the
    aperture; each aperture is scaled on its own call, so each can use its
    own order. `mode` is XMC_SCALESPEC's 'Scale to' choice: 'Median' (the
    default: MC_GETSPECSCALE against the median spectrum of the spectra
    not removed from that order), 'Spectrum' (against spectrum `spectrum`,
    0-based, which must not be removed; every spectrum's factor is computed
    and those of the removed ones are set back to 1) or 'None' (all 1).
    `srange` is the wavelength range chosen with 's'; by default the middle
    80% of the order. As XCOMBSPEC does, the scaling starts from the loaded
    spectra and is refused once the aperture has been scaled or corrected.
    Not calling this at all is not clicking 'Scale Spectra': nothing is
    scaled and the history says nothing about scaling. Returns the factors
    and the median RMS of each order after scaling.
    """
    ap = aperture - 1
    if loaded["modspec"][ap] > 0:
        raise ValueError("Cannot perform this operation again. Please reload spectra and start over.")
    if mode not in SCALE_MODES:
        raise ValueError(f"mode must be one of {SCALE_MODES}")
    zo = loaded["orders"].index(order)
    spc = [s_.copy() for s_ in loaded["ospc"][ap]]           # (*state.d.ospec)[ap]
    nf = spc[zo].shape[0]
    mask = (np.asarray(loaded["spcmask"][ap][zo]) if spcmask is None
            else np.asarray(spcmask)).astype(int)
    wave = spc[zo][0, 0]
    flux = spc[zo][:, 1, :]
    # XMC_SCALESPEC_INITCOMMON: a FLOAT range
    wmin, wmax = np.nanmin(wave), np.nanmax(wave)
    if srange is None:
        wr = wmax - wmin
        srange = [wmin + np.float32(0.1) * wr, wmax - np.float32(0.1) * wr]
    else:
        srange = np.sort(np.asarray(srange, dtype=np.float32))
        if not (wmin < srange[0] and srange[1] < wmax):           # MC_CRANGE /KGT /KLT
            raise ValueError("The scale wavelength range must lie inside the order.")
    srange = np.asarray(srange, dtype=np.float32)
    z = (wave > srange[0]) & (wave < srange[1])
    zm = np.where(mask == 1)[0]
    scales = np.ones(nf, dtype=np.float32)
    if mode == "Median":
        scales[zm] = getspecscale(flux[zm][:, z])
    elif mode == "Spectrum":
        if spectrum is None or spectrum not in zm:
            raise ValueError("Choose a spectrum that is not removed from the scale order.")
        scales[zm] = getspecscale(flux[:, z], idx=spectrum)[zm]
    medrms = []
    for j in range(len(loaded["orders"])):
        spc[j][:, 1, :] *= scales[:, None]
        spc[j][:, 2, :] *= scales[:, None]
        if zm.size == 1:                    # REFORM leaves one spectrum: MC_MEANCOMB's 1-D case
            _, var, _ = I.meancomb(spc[j][zm[0], 1, :], rms=True)
        else:
            _, var, _ = I.meancomb(spc[j][zm, 1, :], rms=True)
        medrms.append(I.median(np.sqrt(var), even=True))
    loaded["spc"][ap] = spc                                   # (*state.d.wspec)[ap]
    loaded["modspec"][ap] = 1
    loaded.setdefault("scales", {})[aperture] = scales
    loaded.setdefault("scaleorder", {})[aperture] = order
    loaded.setdefault("medrms", {})[aperture] = medrms
    loaded.setdefault("scalemode", {})[aperture] = mode
    loaded.setdefault("srange", {})[aperture] = srange
    return scales, medrms


def maskspec(wave, nspec, remove=(), regions=()):
    """
    XMC_MASKSPEC: the spectrum mask and pixel mask that one 'Prune Spectra'
    session returns for one order.

    `remove` lists the spectra (0-based) taken out ('Remove'); `regions`
    lists the wavelength ranges masked ('Mask' with 's'), each (w1, w2) for
    every spectrum not removed (the 'All' button) or (w1, w2, spectra).
    A range becomes pixels as XMC_MASKSPEC converts it: TABINV on the
    NaN-trimmed wavelengths, in FLOAT, plus the first good pixel, the two
    indices truncated, both ends included. Returns specmask (1 kept) and
    pixmask (nspec, npix; 1 use, 0 masked).
    """
    from triplespec_reducer.wavecal import nantrim
    wave = np.asarray(wave, dtype=np.float32)
    removespec = np.zeros(nspec, dtype=int)
    removespec[list(remove)] = 1
    pix = np.ones((nspec, wave.size), dtype=np.float32)
    good = nantrim(wave, 2)
    for reg in regions:
        tmp = np.sort(np.asarray(reg[:2], dtype=np.float32))     # state.p.reg is FLOAT
        if len(reg) > 2:
            z = np.asarray(reg[2], dtype=int)
        else:
            z = np.where(removespec == 0)[0]
        if z.size == 0:
            continue
        idx = I.tabinv(wave[good], tmp) + np.float32(good[0])
        lo, hi = I.itrunc(idx[0]), I.itrunc(idx[1])
        pix[np.ix_(z, np.arange(lo, hi + 1))] = np.nan
    specmask = (removespec <= 0).astype(int)
    pixmask = np.where(np.isfinite(pix), pix, 0.0).astype(int)
    return specmask, pixmask


def prune(loaded, order, aperture=1, remove=(), regions=()):
    """
    XCOMBSPEC_PRUNESPEC: one XMC_MASKSPEC session on `order` of one
    aperture. As in XCOMBSPEC, each session starts from no mask and its
    result replaces that order's spectrum mask and pixel mask; the pixel
    mask applies to that order only.
    """
    ap = aperture - 1
    zo = loaded["orders"].index(order)
    wave = loaded["ospc"][ap][zo][0, 0]
    nspec = loaded["ospc"][ap][zo].shape[0]
    specmask, pixmask = maskspec(wave, nspec, remove=remove, regions=regions)
    loaded["spcmask"][ap][zo] = specmask
    loaded["pixmask"][ap][zo] = pixmask
    return specmask, pixmask


def robustsg(x, y, width, thresh, eps, degree=None):
    """
    MC_ROBUSTSG: Savitzky-Golay smoothing that sets outliers aside, in the
    inputs' precision (FINDGEN x, float spectra in XCOMBSPEC).
    """
    from triplespec_reducer.clean import savitzky_golay
    x = np.asarray(x)
    y = np.asarray(y)
    ave = I.moments(x[:-1] - x[1:])[0]
    ww = width / abs(ave)
    gb = np.ones(y.size, dtype=int)
    gb[~(np.isfinite(y) & np.isfinite(x))] = 2
    zg = np.where(gb == 1)[0]
    yy, xx = y[zg], x[zg]
    goodyy = savitzky_golay(yy, ww, degree=degree)
    test = yy - goodyy
    mean, _, sig, _ = I.moments(test)
    good = np.where(np.abs((test - mean) / sig) <= thresh)[0]
    if not (good.size == zg.size or good.size == 0):
        for _ in range(2, 11):
            sig_old = sig
            res = savitzky_golay(yy[good], ww, degree=degree)
            goodyy = I.linterp(xx[good], res, xx)
            test = yy[good] - goodyy[good]
            mean, _, sig, _ = I.moments(test)
            test = yy - goodyy
            good = np.where(np.abs((test - mean) / sig) <= thresh)[0]
            if (sig_old - sig) / sig_old < eps:
                break
    out = np.full(y.size, np.nan, dtype=np.float32)          # FLTARR + NaN
    out[zg] = goodyy
    return np.stack([x, out])


def _clip(idx, n):
    """An IDL subscript array: out-of-range elements are clipped to the ends."""
    return np.clip(np.asarray(idx, dtype=np.int64), 0, n - 1)


def meancomb_stackmask(data, mask, datavar=None, robust=None, rms=False):
    """
    MC_MEANCOMB of the spectra in use (`data`, nz rows) with the pixel mask
    of the full stack (`mask`, nspec >= nz rows), as XCOMBSPEC and MC_SPECCOR
    call it once a spectrum is removed. The IDL then works by subscripts:
    the outlier flags of row r go to mask row r, the masked pixels are found
    as positions in the full mask and zeroed in the smaller array with
    subscripts clipped to its last element, and the unweighted mean divides
    by the count of the full mask's rows. With equal sizes this is MC_MEANCOMB.
    Returns (mean, variance, ogoodbad), ogoodbad with the full mask's rows.
    """
    data = np.asarray(data)
    mask = np.asarray(mask)
    if mask.shape == data.shape:
        return I.meancomb(data, datavar=datavar, mask=mask, robust=robust, rms=rms)
    P = np.float32 if I.isf32(data) else np.float64
    data = data.astype(P)
    imask = mask.copy()
    if robust is not None:
        for i in range(data.shape[1]):
            r = I.findoutliers(data[:, i], robust)
            imask[np.where(r == 0)[0], i] = 0
    zf = np.flatnonzero(imask == 0)                 # WHERE(imask eq 0): full-stack positions
    with np.errstate(divide="ignore", invalid="ignore"):
        if datavar is not None:
            V = np.float32 if I.isf32(datavar) else np.float64
            datavar = np.asarray(datavar).astype(V)
            tmp = (V(1.) / datavar)
            if zf.size:
                tmp.flat[_clip(zf, tmp.size)] = 0
            mvar = 1.0 / I._seqsum0(tmp)
            tmp = data / datavar
            if zf.size:
                tmp.flat[_clip(zf, tmp.size)] = 0
            mean = I._seqsum0(tmp) * mvar
        else:
            idata = data.copy()
            if zf.size:
                idata.flat[_clip(zf, idata.size)] = P(0)
            ndat = I._seqsum0(imask).astype(np.float32)          # over the full stack
            mean = I._seqsum0(idata) / ndat
            mvar = (I._seqsum0(I.sq(idata.astype(float))) - ndat * I.sq(mean)) / (ndat - np.float32(1))
            if not rms:
                mvar = mvar / ndat
    return mean.astype(np.float32), mvar.astype(np.float32), imask


def medcomb_stackmask(data, mask, mad=False):
    """
    MC_MEDCOMB with the pixel mask of the full stack: for each pixel the
    rows the mask keeps are taken from the spectra in use by subscript,
    clipped to the last one (so it can count twice), and the variance is
    divided by the full mask's count. With equal sizes this is MC_MEDCOMB.
    """
    data = np.asarray(data)
    mask = np.asarray(mask)
    if mask.shape == data.shape:
        return I.medcomb(data, mask=mask, mad=mad)
    P = np.float32 if I.isf32(data) else np.float64
    data = data.astype(P)
    nz, ncols = data.shape
    k = np.float32(1.4826)
    med = np.zeros(ncols, dtype=np.float32)                   # FLTARR
    var = np.zeros(ncols, dtype=np.float32)
    for i in range(ncols):
        z = np.where(mask[:, i] == 1)[0]
        ngood = z.size
        if ngood == 0:
            z = np.array([-1])                                # WHERE's -1
        col = data[_clip(z, nz), i]
        med[i] = I.median(col, even=True)
        v = I.sq(k * I.median(np.abs(col - med[i]), even=True))
        if not mad:
            with np.errstate(divide="ignore"):
                v = np.float32(v) / np.float32(ngood)
        var[i] = v
    var[var == 0] = 1.
    return med, var


def correct_shapes(loaded, aperture=1, spcmask=None, fwidth=4, intended_masks=False):
    """
    XCOMBSPEC_CORRECTSPEC: MC_SPECCOR on every order of one aperture.

    As MC_SPECCOR does: the robust mean uses the pixel mask of the full stack
    against the spectra in use (see meancomb_stackmask: with a spectrum
    removed, the mean is scaled by nkept/nspec and the masks shift), the
    corrected stack is a zeroed FLTARR filled only where the mean spectrum
    is finite (other pixels come back 0), and the low-frequency ratio comes
    from IDL's FFT (here numpy's; see docs/validation.md). The errors keep
    their values outside those pixels.

    XCOMBSPEC passes the pixel masks of all apertures; with two or more
    apertures (not combined) that array is 3-D, MC_SPECCOR's argument check
    rejects it and nothing is corrected. Likewise one spectrum alone (a 1-D
    stack). Then this returns False and leaves the spectra as they are.
    `intended_masks=True` (not what IDL does) uses each order's own mask of
    this aperture, the kept spectra's rows only, and corrects in those cases.
    """
    ap = aperture - 1
    if loaded["modspec"][ap] >= 3:
        raise ValueError("Cannot perform this operation again. Please reload spectra and start over.")
    spc = loaded["spc"][ap]
    if not intended_masks and (loaded["naps"] > 1 or spc[0].shape[0] < 2):
        return False
    for j in range(len(loaded["orders"])):
        stack = spc[j][:, 1, :]
        err = spc[j][:, 2, :]
        nf, ncols = stack.shape
        if spcmask is not None:
            specmask = np.asarray(spcmask[j] if np.ndim(spcmask) > 1 else spcmask, np.float32)
        else:
            specmask = np.asarray(loaded["spcmask"][ap][j], np.float32)
        pixmask = np.asarray(loaded["pixmask"][ap][j], np.float32)
        z = np.where(specmask == 1)[0]
        if intended_masks:
            meanspec, _, npm_z = I.meancomb(stack[z], robust=3, mask=pixmask[z])
            npm = np.ones((nf, ncols), dtype=np.float32)
            npm[z] = npm_z
        else:
            meanspec, _, npm = meancomb_stackmask(stack[z], pixmask, robust=3)
        tstack = stack.copy()
        for i in range(nf):
            if not specmask[i]:
                continue
            bad = npm[i] == 0
            tstack[i, bad] = meanspec[bad]
        x = np.arange(ncols, dtype=np.float32)
        for i in range(nf):
            tstack[i] = robustsg(x, tstack[i], 10, 5, 0.1)[1]
        y = (np.concatenate([np.arange(ncols // 2 + 1), np.arange(ncols // 2)[::-1]])
             if ncols % 2 else np.concatenate([np.arange(ncols // 2), np.arange(ncols // 2)[::-1]])
             ).astype(np.float32)
        filt = (np.float32(1.0) / (1 + I.ipow(y / np.float32(fwidth), 10))).astype(np.float32)
        good = np.where(np.isfinite(meanspec))[0]
        ng = good.size
        ref_lo = I.idl_fft(I.idl_fft(meanspec[good]) * filt[:ng], inverse=True)
        cstack = np.zeros((nf, ncols), dtype=np.float32)
        for i in range(nf):
            lo = I.idl_fft(I.idl_fft(tstack[i, good]) * filt[:ng], inverse=True)
            ratio = ref_lo / lo
            cstack[i, good] = (stack[i, good] * ratio).real
            err[i, good] = (err[i, good] * ratio).real
        stack[:, :] = cstack
    loaded.setdefault("corrected", {})[aperture] = True
    loaded["modspec"][ap] = 3
    return True


def combine(loaded, statistic="Robust Weighted Mean", thresh=8.0, spcmask=None,
            pixmask=None, intended_masks=False):
    """
    XCOMBSPEC_COMBINESPEC: one spectrum per order and aperture.

    The masks are those Prune Spectra set (prune()); `spcmask`, per
    aperture and order, and `pixmask`, per aperture and order an
    (nspec, npix) array of 1 (use) and 0 (prune), replace them when given.
    As XCOMBSPEC does, the statistics get the spectra in use but the pixel
    mask of the full stack, so with a spectrum removed an unweighted mean
    is divided by the full count, a median takes the last spectrum twice,
    and the flags are combined from subscript-clipped rows (see
    meancomb_stackmask, medcomb_stackmask). `intended_masks=True` (not what
    IDL does) gives every statistic only the kept spectra's mask rows.
    """
    stat = STATISTICS.index(statistic)
    orders, naps = loaded["orders"], loaded["naps"]
    npix = loaded["spc"][0][0].shape[2]
    out = np.full((len(orders) * naps, 4, npix), np.nan, dtype=np.float32)   # FLTARR
    for i in range(naps):
        for j in range(len(orders)):
            st = loaded["spc"][i][j]
            nf = st.shape[0]
            m = np.asarray(loaded["spcmask"][i][j] if spcmask is None else spcmask[i][j])
            z = np.where(m == 1)[0]
            spec = st[z, 1, :]
            err = st[z, 2, :]
            flag = st[z, 3, :]
            full = np.asarray(loaded["pixmask"][i][j] if pixmask is None else pixmask[i][j]).astype(int)
            mask = full[z] if intended_masks else full
            mc = I.meancomb if intended_masks else meancomb_stackmask
            md = (lambda d, mask, mad=False: I.medcomb(d, mask=mask, mad=mad)) \
                if intended_masks else medcomb_stackmask
            ogb = None
            if stat == 0:
                mean, mvar, ogb = mc(spec, datavar=I.sq(err), mask=mask, robust=thresh)
            elif stat == 1:
                mean, mvar, ogb = mc(spec, rms=True, mask=mask, robust=thresh)
            elif stat == 2:
                mean, mvar, ogb = mc(spec, mask=mask, robust=thresh)
            elif stat == 3:
                mean, mvar, _ = mc(spec, datavar=I.sq(err), mask=mask)
            elif stat == 4:
                mean, mvar, _ = mc(spec, rms=True, mask=mask)
            elif stat == 5:
                mean, mvar, _ = mc(spec, mask=mask)
            elif stat == 6:
                mean, mvar = md(spec, mask, mad=True)
            elif stat == 7:
                mean, mvar = md(spec, mask)
            else:
                # TOTAL(..., 2, /NAN): in order down the stack, NaNs skipped;
                # spec*mask has the smaller array's size (the first mask rows)
                sm = spec * mask[:spec.shape[0]]
                mean = np.array([I.ftotal(c[np.isfinite(c)]) for c in sm.T])
                ev = mask[:spec.shape[0]] * I.sq(err.astype(float))
                mvar = np.array([I.itotal(c[np.isfinite(c)]) for c in ev.T])
            n = j * naps + i
            out[n, 0] = st[0, 0]
            out[n, 1] = mean
            out[n, 2] = np.sqrt(mvar)
            # mask*BYTE(flag), then *ogoodbad: each product the smaller size
            fl = np.nan_to_num(flag).astype(np.uint8)
            tmp = mask[:fl.shape[0]] * fl
            if stat <= 2:
                tmp = tmp * ogb[:fl.shape[0]]
            rows = np.arange(tmp.shape[0]) if intended_masks else _clip(z, tmp.shape[0])
            out[n, 3] = combflagstack(tmp[rows])
    loaded["combined"] = out
    loaded["statistic"] = statistic
    loaded["thresh"] = thresh
    return out


def _idlstring(v):
    """IDL's STRING() of a number with no format: FLOAT %#13.6g, DOUBLE %#16.8g."""
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return f"{float(v):#13.6g}" if I.isf32(v) else f"{float(v):#16.8g}"


def gethdrinfo(header, keywords=None):
    """
    MC_GETHDRINFO: each keyword typed by FXPAR, with its FXPAR comment. With
    no list, every keyword with a value and then HISTORY joined into one
    string; a listed keyword that is missing reads "Doesn't Exist".
    """
    vals, coms = {}, {}
    if keywords is None:
        for c in header.cards:
            k = c.keyword
            if not k or k in ("COMMENT", "HISTORY") or k in vals or c.image[8:9] != "=":
                continue
            vals[k] = I.fxpar(header, k)
            coms[k] = I.fxcomment(header, k)
        vals["HISTORY"] = "".join(str(h) for h in header.get("HISTORY", []))
        coms["HISTORY"] = ""
        return vals, coms
    for k in keywords:
        if k in header:
            vals[k], coms[k] = I.fxpar(header, k), I.fxcomment(header, k)
        else:
            vals[k], coms[k] = "Doesn't Exist", ""
    return vals, coms


def _tenv(d, m, s):
    """astrolib TENV of three DOUBLE parts: negative when any part's text has a '-'."""
    v = abs(d) + abs(m) / 60.0 + abs(s) / 3600.0
    neg = any("-" in _idlstring(np.float64(x)) for x in (d, m, s))
    return -v if neg else v


def _sixty(x):
    """astrolib SIXTY of a DOUBLE: hours (FIX), minutes (FIX), seconds."""
    ss, mm, dd = abs(3600.0 * x), abs(60.0 * x), abs(x)
    r0 = float(int(dd))
    r1 = float(int(mm - 60.0 * r0))
    r2 = ss - 3600.0 * r0 - 60.0 * r1
    if x < 0.0:
        if r0 != 0:
            r0 = -r0
        elif r1 != 0:
            r1 = -r1
        else:
            r2 = -r2
    return r0, r1, r2


def _total(vals):
    """TOTAL of the values in their IDL type: a running sum (FLOAT or DOUBLE)."""
    if all(I.isf32(v) for v in vals):
        return I.ftotal(np.array(vals, dtype=np.float32))
    return I.itotal(np.array(vals, dtype=float))


def _clock(vals):
    """TENV of 'h:m:s' strings, one DOUBLE per header."""
    return np.array([_tenv(*[float(t) for t in str(v).split(":")]) for v in vals])


def avehdrs(infos, time_obs="UT", posangle="ROTANG", ha="HA", airmass="AIRMASS"):
    """
    MC_AVEHDRS: the first header, with the mean time (wrapping midnight),
    position angle (wrapping 360), hour angle and airmass of all of them,
    computed in the types FXPAR gave the values. `infos` is a list of
    gethdrinfo() pairs, one per header.
    """
    n = len(infos)
    vals, coms = dict(infos[0][0]), dict(infos[0][1])
    names = list(vals)

    def col(k):
        # REPLICATE + COPY_STRUCT: every value takes the first header's field type
        t0 = vals[k]
        if isinstance(t0, str):
            return [str(inf[0][k]) for inf in infos]
        if isinstance(t0, (int, np.integer)):
            return [int(np.trunc(float(inf[0][k]))) for inf in infos]
        return [type(t0)(inf[0][k]) for inf in infos]
    if time_obs.upper() in names:
        k = time_obs.upper()
        parts = np.array([[float(t) for t in str(v).split(":")] for v in col(k)])
        if (parts[:, 0] == 23.0).any() and (parts[:, 0] == 0.0).any():
            parts[parts[:, 0] == 0.0] += 24.
            ave = I.itotal([_tenv(*p) for p in parts]) / float(n)
            if ave > 24.:
                ave = ave - 24.
        else:
            ave = I.itotal([_tenv(*p) for p in parts]) / float(n)
        r = _sixty(ave)
        vals[k] = f"{int(r[0]):02d}:{int(r[1]):02d}:" + _idlstring(np.float64(r[2])).strip()
        coms[k] = " Average " + coms[k].strip()
    if posangle.upper() in names and not isinstance(vals[posangle.upper()], str):
        k = posangle.upper()
        pos = col(k)
        fl = I.isf32(pos[0])
        typ = np.float32 if fl else np.float64
        ipos = [int(np.trunc(p)) for p in pos]
        if 359 in ipos and 0 in ipos:
            pos = [typ(p + typ(360.)) if ip == 0 else typ(p) for p, ip in zip(pos, ipos)]
            newp = typ(_total(pos) / (np.float32(n) if fl else float(n)))
            if newp > 360.:
                newp = typ(newp - typ(360.))
        else:
            newp = typ(_total(pos) / (np.float32(n) if fl else float(n)))
        vals[k], coms[k] = newp, " Average " + infos[0][1][k].strip()
        if posangle.upper() != "POSANGLE":
            vals["POSANGLE"], coms["POSANGLE"] = newp, " Average " + infos[0][1][k].strip()
    if ha.upper() in names:
        k = ha.upper()
        ave = I.itotal(_clock(col(k))) / float(n)
        r = _sixty(abs(ave))
        newha = ("-" if ave < 0 else "+") + f"{int(r[0]):02d}:{int(r[1]):02d}:" \
            + _idlstring(np.float64(r[2])).strip()
        vals[k], coms[k] = newha, " Average " + infos[0][1][k].strip()
        if ha.upper() != "HA":
            vals["HA"], coms["HA"] = newha, " Average " + infos[0][1][k].strip()
    if airmass.upper() in names and not isinstance(vals[airmass.upper()], str):
        k = airmass.upper()
        am = col(k)
        if I.isf32(am[0]):
            newa = np.float32(_total(am) / np.float32(n))
        else:
            newa = np.float64(_total(am) / float(n))
        # stored into the first header's field, which keeps its type
        vals[k] = newa if not isinstance(vals[k], int) else int(newa)
        coms[k] = " Average " + coms[k].strip()
    return vals, coms


def history(loaded, basic=False):
    """
    XCOMBSPEC_WRITEFILE's history text: per aperture the files sentence,
    the scale factors (only when that aperture was scaled), the spectra
    removed from each order and the shape correction, then the statistic.
    The removed files are named by FILES[z]; with the apertures combined z
    counts spectra, not files, and the subscripts are clipped, as in IDL.
    """
    files = [os.path.basename(f) for f in loaded["files"]]
    hist = "This FITS header is from the first spectrum in the stack.  " if basic else ""
    for i in range(loaded["naps"]):
        ap = i + 1
        hist += "  The files " + ", ".join(files) + " were combined.  "
        hist += f"  Aperture {ap:02d} Info: "
        if ap in loaded.get("scaleorder", {}):
            hist += (f"The scale factors were determined using order {int(loaded['scaleorder'][ap]):02d}"
                     ".  The scale factors are "
                     + ", ".join(_idlstring(np.float32(v)).strip() for v in loaded["scales"][ap])
                     + ".  The median RMS deviation in each order is "
                     + ", ".join(_idlstring(np.float32(v)).strip() for v in loaded["medrms"][ap])
                     + ".  ")
        for j, o in enumerate(loaded["orders"]):
            z = np.where(np.asarray(loaded["spcmask"][i][j]) == 0)[0]
            if z.size:
                hist += ("  The spectrum(a) from file(s) "
                         + ", ".join(files[k] for k in _clip(z, len(files)))
                         + f" was (ere) removed from order {int(o):02d}.")
        if loaded.get("corrected", {}).get(ap):
            hist += "  The spectral shapes have been corrected.  "
    hist += "The spectra were combined using a " + loaded["statistic"] + ".  "
    if STATISTICS.index(loaded["statistic"]) <= 2:
        hist += "The robust threshold was " + _idlstring(np.float32(loaded["thresh"])).strip() + "."
    return hist


def write(path, loaded, basic=False):
    """
    XCOMBSPEC_WRITEFILE: the first spectrum's header through MC_AVEHDRS,
    every value written back by FXADDPAR in the type FXPAR gave it, the
    earlier history under 'Xspextool History', XCOMBSPEC's own under
    'Xcombspec History' (MC_SPLITTEXT at 70 characters), ITOT (a FLOAT
    field) and the file name. `basic` is xcombspec,/BASIC: the first
    header as it is, not averaged, and a sentence saying so.
    """
    first = loaded["headers"][0]
    # XCOMBSPEC_LOADSPEC replicates the first header for every file
    infos = [gethdrinfo(first)] * len(loaded["headers"])
    if basic:
        vals, coms = dict(infos[0][0]), dict(infos[0][1])
    else:
        vals, coms = avehdrs(infos)
    vals["NAPS"] = type(vals["NAPS"])(loaded["naps"]) if "NAPS" in vals else loaded["naps"]
    hdr = fits.Header()
    for k, v in vals.items():
        if k in ("SIMPLE", "BITPIX", "EXTEND") or k.startswith("NAXIS"):
            continue                         # the data's own, set on writing
        if k == "HISTORY":
            I.addhistlabel(hdr, "Xspextool History")
            for j in range(0, len(v), 70):
                I.addhistory(hdr, v[j:j + 70])              # FXADDPAR 'HISTORY'
            continue
        I.fxaddpar(hdr, k, v, coms[k])
    I.addhistlabel(hdr, "Xcombspec History")
    I.sxaddhist(hdr, history(loaded, basic=basic))
    I.fxaddpar(hdr, "ITOT", np.float32(loaded["itot"]), " Total integration time (sec)")
    I.fxaddpar(hdr, "FILENAME", os.path.basename(path), coms.get("FILENAME", ""))
    fits.PrimaryHDU(loaded["combined"], header=hdr).writeto(path, overwrite=True)
