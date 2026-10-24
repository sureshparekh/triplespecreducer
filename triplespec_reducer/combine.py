"""
Combining the spectra of one target: XCOMBSPEC.

    every spectrum onto the first one's grid   MC_INTERPSPEC (LEAVENANS)
    scale factors from one order               XMC_SCALESPEC 'Median'
                                               (MC_GETSPECSCALE over the
                                               middle 80% of that order)
    optional prune of a spectrum or pixels     XMC_MASKSPEC
    optional shape correction                  MC_SPECCOR
    combine, per order and aperture            MC_MEANCOMB / MC_MEDCOMB

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
                nf, ne = interpspec(data[n, 0], data[n, 1], fs[n, 0][z],
                                    iyerror=data[n, 2], leavenans=True)
                nb = interpflagspec(data[n, 0], data[n, 3], fs[n, 0][z])
                a, row = (0, i * naps + k) if combine_apertures else (k, i)
                spc[a][j][row, 0, z] = fs[n, 0][z]
                spc[a][j][row, 1, z] = nf
                spc[a][j][row, 2, z] = ne
                spc[a][j][row, 3, z] = nb
    return dict(spc=spc, orders=orders, naps=fnaps, files=list(files), itot=itot,
                headers=hdrs, first=first, combine_apertures=combine_apertures)


def scale(loaded, order, aperture=1, spcmask=None, srange=None):
    """
    XCOMBSPEC_SCALESPEC with XMC_SCALESPEC 'Median'.

    Factors come from `order` over the middle 80% of its wavelengths (or
    `srange`) and are applied to every order. Returns the factors and the
    median RMS of each order after scaling.
    """
    spc = loaded["spc"][aperture - 1]
    zo = loaded["orders"].index(order)
    nf = spc[zo].shape[0]
    mask = np.ones(nf, dtype=int) if spcmask is None else np.asarray(spcmask)
    wave = spc[zo][0, 0]
    flux = spc[zo][:, 1, :]
    if srange is None:
        wmin, wmax = np.nanmin(wave), np.nanmax(wave)
        wr = wmax - wmin
        srange = [wmin + 0.1 * wr, wmax - 0.1 * wr]
    z = (wave > srange[0]) & (wave < srange[1])
    zm = np.where(mask == 1)[0]
    scales = np.ones(nf, dtype=np.float32)
    scales[zm] = getspecscale(flux[zm][:, z])
    medrms = []
    for j in range(len(loaded["orders"])):
        spc[j][:, 1, :] *= scales[:, None]
        spc[j][:, 2, :] *= scales[:, None]
        _, var, _ = I.meancomb(spc[j][zm, 1, :], rms=True)
        medrms.append(I.median(np.sqrt(var), even=True))
    loaded.setdefault("scales", {})[aperture] = scales
    loaded.setdefault("scaleorder", {})[aperture] = order
    loaded.setdefault("medrms", {})[aperture] = medrms
    return scales, medrms


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


def correct_shapes(loaded, aperture=1, spcmask=None, fwidth=4):
    """
    XCOMBSPEC_CORRECTSPEC: MC_SPECCOR on every order of one aperture.

    As MC_SPECCOR does: the robust mean uses the pixel mask of the full stack
    against the spectra in use, the corrected stack is a zeroed FLTARR filled
    only where the mean spectrum is finite (other pixels come back 0), and
    the low-frequency ratio comes from IDL's FFT (here numpy's; see
    docs/validation.md). The errors keep their values outside those pixels.
    """
    spc = loaded["spc"][aperture - 1]
    for j in range(len(loaded["orders"])):
        stack = spc[j][:, 1, :]
        err = spc[j][:, 2, :]
        nf, ncols = stack.shape
        specmask = (np.ones(nf, np.float32) if spcmask is None
                    else np.asarray(spcmask[j] if np.ndim(spcmask) > 1 else spcmask, np.float32))
        z = np.where(specmask == 1)[0]
        meanspec, _, npm = I.meancomb(stack[z], robust=3,
                                      mask=np.ones((z.size, ncols), dtype=int))
        tstack = stack.copy()
        for i in range(nf):
            if not specmask[i]:
                continue
            if i < npm.shape[0]:
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


def combine(loaded, statistic="Robust Weighted Mean", thresh=8.0, spcmask=None,
            pixmask=None):
    """
    XCOMBSPEC_COMBINESPEC: one spectrum per order and aperture.

    `spcmask`, per aperture and order, switches whole spectra off; `pixmask`,
    per aperture and order, is an (nfiles, npix) array of 1 (use) and 0 (prune).
    """
    stat = STATISTICS.index(statistic)
    orders, naps = loaded["orders"], loaded["naps"]
    npix = loaded["spc"][0][0].shape[2]
    out = np.full((len(orders) * naps, 4, npix), np.nan, dtype=np.float32)   # FLTARR
    for i in range(naps):
        for j in range(len(orders)):
            st = loaded["spc"][i][j]
            nf = st.shape[0]
            m = np.ones(nf, dtype=int) if spcmask is None else np.asarray(spcmask[i][j])
            z = np.where(m == 1)[0]
            spec = st[z, 1, :]
            err = st[z, 2, :]
            flag = st[z, 3, :]
            mask = np.ones(spec.shape, dtype=int) if pixmask is None else np.asarray(pixmask[i][j])[z]
            ogb = None
            if stat == 0:
                mean, mvar, ogb = I.meancomb(spec, datavar=I.sq(err), mask=mask, robust=thresh)
            elif stat == 1:
                mean, mvar, ogb = I.meancomb(spec, rms=True, mask=mask, robust=thresh)
            elif stat == 2:
                mean, mvar, ogb = I.meancomb(spec, mask=mask, robust=thresh)
            elif stat == 3:
                mean, mvar, _ = I.meancomb(spec, datavar=I.sq(err), mask=mask)
            elif stat == 4:
                mean, mvar, _ = I.meancomb(spec, rms=True, mask=mask)
            elif stat == 5:
                mean, mvar, _ = I.meancomb(spec, mask=mask)
            elif stat == 6:
                mean, mvar = I.medcomb(spec, mask=mask, mad=True)
            elif stat == 7:
                mean, mvar = I.medcomb(spec, mask=mask)
            else:
                # TOTAL(..., 2, /NAN): in order down the stack, NaNs skipped
                sm = spec * mask
                mean = np.array([I.ftotal(c[np.isfinite(c)]) for c in sm.T])
                ev = mask * I.sq(err.astype(float))
                mvar = np.array([I.itotal(c[np.isfinite(c)]) for c in ev.T])
            n = j * naps + i
            out[n, 0] = st[0, 0]
            out[n, 1] = mean
            out[n, 2] = np.sqrt(mvar)
            tmp = (mask * np.nan_to_num(flag).astype(np.uint8))
            if stat <= 2:
                tmp = tmp * ogb
            out[n, 3] = combflagstack(tmp)
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


def write(path, loaded):
    """
    XCOMBSPEC_WRITEFILE: the first spectrum's header through MC_AVEHDRS,
    every value written back by FXADDPAR in the type FXPAR gave it, the
    history, ITOT (a FLOAT field) and the file name.
    """
    first = loaded["headers"][0]
    # XCOMBSPEC_LOADSPEC replicates the first header for every file
    vals, coms = avehdrs([gethdrinfo(first)] * len(loaded["headers"]))
    vals["NAPS"] = type(vals["NAPS"])(loaded["naps"]) if "NAPS" in vals else loaded["naps"]
    hdr = fits.Header()
    for k, v in vals.items():
        if k in ("SIMPLE", "BITPIX", "EXTEND") or k.startswith("NAXIS"):
            continue                         # the data's own, set on writing
        if k == "HISTORY":
            for j in range(0, len(v), 70):
                hdr.add_history(v[j:j + 70])
            continue
        I.fxaddpar(hdr, k, v, coms[k])
    hist = "The files " + ", ".join(os.path.basename(f) for f in loaded["files"]) + " were combined."
    for ap, order in loaded.get("scaleorder", {}).items():
        hist += (f"  Aperture {ap:02d} Info: The scale factors were determined using order "
                 f"{order:02d}.  The scale factors are "
                 + ", ".join(_idlstring(np.float32(s)).strip() for s in loaded["scales"][ap])
                 + ".  The median RMS deviation in each order is "
                 + ", ".join(_idlstring(np.float32(s)).strip() for s in loaded["medrms"][ap]) + ".")
    hist += f"  The spectra were combined using a {loaded['statistic']}."
    if STATISTICS.index(loaded["statistic"]) <= 2:
        hist += f"  The robust threshold was {_idlstring(np.float32(loaded['thresh'])).strip()}."
    hdr.add_history("Xcombspec History: " + hist)
    I.fxaddpar(hdr, "ITOT", np.float32(loaded["itot"]), " Total integration time (sec)")
    I.fxaddpar(hdr, "FILENAME", os.path.basename(path), coms.get("FILENAME", ""))
    fits.PrimaryHDU(loaded["combined"], header=hdr).writeto(path, overwrite=True)
