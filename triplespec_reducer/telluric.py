"""
Telluric correction and flux calibration with an A0 V star: XTELLCOR.

The method is Vacca, Cushing & Rayner (2003). A model of Vega, shifted to
the star's velocity, broadened to the instrument with a kernel, its
hydrogen lines scaled to the star's, reddened to the star's B-V and
scaled to its V magnitude, divided by the observed star gives the
telluric-plus-instrument response; multiplying the object by it corrects
the telluric absorption and calibrates the flux.

    load           XTELLCOR_LOADSPEC      airmass difference warning at 0.1
    kernel         XMC_CONKERN            (Deconvolution: MC_VEGACORR for
                                           the velocity, MC_VEGACONV for
                                           the kernel, on one H line) and
                                           XTELLCOR_CONKERNEL (resampled to
                                           every order); or the IP method
    line scales    XMC_SCALELINES         control points at every H I line
    telluric       MC_MKTELLSPEC          per order, cut regions bridged
    shift          XMC_FINDSHIFT          sub-pixel shift per order
    write          XTELLCOR_WRITEFILE

XMC_CONKERN keeps a copy of the second Vega continuum under the name of the
first (`fc2vega:fcvega` in its state), so the velocity and the normalized
Vega used for the kernel are built on fcvega. That is reproduced.
"""

import os

import numpy as np
from astropy.io import fits
from scipy.io import readsav

from triplespec_reducer import idl as I
from triplespec_reducer._paths import SPEXTOOL_DIR, INSTRUMENT_DIR
from triplespec_reducer.combine import interpspec, interpflagspec, combflagstack
from triplespec_reducer.extract import read_spectrum, caaspecflags
from triplespec_reducer.instrument import read_atran


# XTELLCOR's unit choices and the MC_CHFUNITS code for each; the Vega model,
# and so the telluric spectrum, is in ergs s-1 cm-2 A-1 (code 1).
UNITS = {"ergs s-1 cm-2 A-1": 1, "ergs s-1 cm-2 Hz-1": 3, "W m-2 um-1": 0,
         "W m-2 Hz-1": 2, "Jy": 4}


def chfunits(wave, flux, iwunit, ifunit, ofunit, error=None):
    """
    MC_CHFUNITS: flux density units, wave in um when iwunit is 0. Each
    conversion is written out as the IDL writes it, so the products and
    quotients are taken in the same order.
    """
    if ifunit == ofunit:
        return flux, error
    D, F = np.float64, np.float32          # the IDL's 10D-style and 10.-style literals
    h, c = D(6.6260755E-34), D(2.99792458E+8)
    nw = {0: lambda: wave, 1: lambda: wave / F(1000.), 2: lambda: wave / F(10000.)}[iwunit]()

    def to_cgs(x):          # into ergs s-1 cm-2 A-1
        return {0: lambda: x / D(10.),
                1: lambda: x,
                2: lambda: x / I.sq(nw) * (c / D(1.0E-5)),
                3: lambda: x / I.sq(nw) * (c / D(1.0E-2)),
                4: lambda: x / I.sq(nw) * (c / D(1.0E21)),
                5: lambda: x / nw / D(10.),
                6: lambda: x / (nw * F(10000.))}[ifunit]()

    def to_out(x):
        return {0: lambda: x * D(10.),
                1: lambda: x,
                2: lambda: x * I.sq(nw) * (D(1.0E-5) / c),
                3: lambda: x * I.sq(nw) * (D(1.0E-2) / c),
                4: lambda: x * I.sq(nw) * (D(1.0E21) / c),
                5: lambda: x * D(10.) * nw,
                6: lambda: x * nw * F(10000.),
                7: lambda: x * F(10.) * nw / (h * c * D(1e6))}[ofunit]()

    f = to_out(to_cgs(flux))
    return f, (None if error is None else to_out(to_cgs(error)))


def read_hlines():
    """HI.dat as XTELLCOR reads it: wavelengths (um) and names."""
    w, n = [], []
    with open(os.path.join(SPEXTOOL_DIR, "HI.dat")) as fh:
        for line in fh:
            if "|" not in line:
                continue
            a, b = line.split("|", 1)
            try:
                w.append(float(a))
            except ValueError:
                continue
            n.append(b.strip().replace("!5", "").replace("!7", "").replace("!N", ""))
    return np.array(w), n


def read_vega(obsmode="XD"):
    """lvega5.sav (lvega99.sav for LowRes15): wavelength in um, three fluxes."""
    f = "lvega99.sav" if obsmode == "LowRes15" else "lvega5.sav"
    d = readsav(os.path.join(SPEXTOOL_DIR, f))
    wv = (d["wvin"].astype(np.float32) / np.float32(10000.)).astype(np.float32)
    fv = d["fvin"].astype(np.float32)
    fcv = d["fcvin"].astype(np.float32)
    fc2v = d["fc2vin"].astype(np.float32) if "fc2vin" in d else fcv
    return wv, fv, fcv, fc2v


def load(std_file, obj_file, bmag, vmag):
    """XTELLCOR_LOADSPEC."""
    std = read_spectrum(std_file)
    obj = read_spectrum(obj_file)
    for d in (std, obj):
        d["spec"] = caaspecflags(d["spec"])
    h = std["header"]
    disp = std["disp"]                                   # MC_FXPAR(stdhdr, 'DISP*')
    wv, fv, fcv, fc2v = read_vega(str(h.get("MODENAME", "XD")).strip())
    aw, at = read_atran(obj["rp"])
    # B and V are FLOAT fields (MC_CFLD type 4) in XTELLCOR; the slit width
    # is the OBJECT's (its MC_READSPEC call comes second), kept in FLOAT
    # fields; DAIRMASS is the DOUBLE difference stored in a FLOAT field.
    swp, swa = np.float32(obj["slitw_pix"]), np.float32(obj["slitw_arc"])
    st = dict(std=std, obj=obj, bmag=np.float32(bmag), vmag=np.float32(vmag), std_path=std_file,
              stdorders=std["orders"], objorders=obj["orders"], objnaps=obj["naps"],
              stddisp=disp, slitw_pix=swp, slitw_arc=swa,
              stdfwhm=disp * swp, obsmode=str(h.get("MODENAME", "XD")).strip(),
              stdairmass=std["airmass"], objairmass=obj["airmass"],
              dairmass=np.float32(std["airmass"] - obj["airmass"]),
              wvega=wv, fvega=fv, cfvega=fcv, cf2vega=fc2v, awave=aw, atrans=at,
              vshift=np.float32(0.0),
              shift=np.zeros((len(obj["orders"]), obj["naps"]), dtype=np.float32))   # FLTARR
    st["airmass_warning"] = abs(std["airmass"] - obj["airmass"]) > 0.1      # the DOUBLEs
    st["hlines"], st["hnames"] = read_hlines()
    return st


# ------------------------------------------------------------- the kernel

def vegacorr(wave, flux, wvin, fvin, fcvin):
    """
    MC_VEGACORR: Vega's velocity relative to the star, by cross-correlation.

    In single precision throughout, as the IDL: FLOAT ALOG and EXP for the
    log-wavelength grid, float INTERPOL and C_CORRELATE, MPFITPEAK on the
    INT lags, and the speed of light as the FLOAT literal 2.99792458E5.
    """
    f32 = np.float32
    cspeed = f32(2.99792458E5)
    wvega = np.asarray(wvin, dtype=f32)
    fnvega = (np.asarray(fvin, f32) / np.asarray(fcvin, f32) - f32(1.0)).astype(f32)
    fdata = (np.asarray(flux, f32) - f32(1.0)).astype(f32)
    wave = np.asarray(wave, dtype=f32)
    wmin = max(wvega.min(), wave.min())
    wmax = min(wvega.max(), wave.max())
    ndata = int(np.sum((wave >= wmin) & (wave <= wmax)))
    nvega = int(np.sum((wvega >= wmin) & (wvega <= wmax)))
    num = max(ndata, nvega)
    acoef = f32(f32(num - 1) / (I.alog(wmax) - I.alog(wmin)))
    bcoef = f32(f32(num) - (acoef * I.alog(wmax)))
    xpon = np.arange(num, dtype=f32) + f32(1.0)
    wx = I.exp(((xpon - bcoef) / acoef).astype(f32))
    fxdata = I.interpol(fdata, wave, wx)
    fxvega = I.interpol(fnvega, wvega, wx)
    lag = (np.arange(num) - num // 2).astype(np.int16)          # INDGEN
    res = I.c_correlate(fxvega, fxdata, lag)
    _, p = I.mpfitpeak(lag, res, nterms=3)
    lshift = f32(p[1])
    vshift = f32((cspeed * lshift) / acoef)
    return vshift, f32(vshift / cspeed)


def vegaconv(wave, flux, wvega, fcvega, fnvega, wmin, wmax, twf=10.):
    """
    MC_VEGACONV: the kernel that turns Vega's line into the star's.

    The normalized star and the scaled model line are Fourier transformed,
    divided, filtered with a smooth window at `twf` times the narrower of
    their widths, and transformed back. Returns the kernel, the line
    centre, the equivalent-width scale and the max and rms deviations.
    """
    wave = np.asarray(wave, dtype=np.float32)
    flux = np.asarray(flux, dtype=np.float32)
    wvega = np.asarray(wvega, dtype=np.float32)
    fnvega = np.asarray(fnvega, dtype=np.float32)
    wabsmin = max(float(wave.min()), float(wvega.min()))
    wabsmax = min(float(wave.max()), float(wvega.max()))
    zdata = np.where((wave >= wmin) & (wave <= wmax))[0]
    zmodl = np.where((wvega >= wmin) & (wvega <= wmax))[0]
    if zdata.size % 2 == 1:
        zdata = np.append(zdata, zdata[-1] + 1)
    if zmodl.size % 2 == 1:
        zmodl = np.append(zmodl, zmodl[-1] + 1)
    ndata, nmodl = zdata.size, zmodl.size
    _, fco = I.mpfitpeak(wvega[zmodl], fnvega[zmodl], nterms=4)
    wline = fco[1]
    dstep = (wave[zdata].max() - wave[zdata].min()) / np.float32(ndata - 1)
    mstep = (wvega[zmodl].max() - wvega[zmodl].min()) / np.float32(nmodl - 1)
    ewstar = dstep * np.cumsum((np.float32(1.0) - flux[zdata]), dtype=np.float32)[-1]
    ewvega = mstep * np.cumsum((np.float32(1.0) - fnvega[zmodl]), dtype=np.float32)[-1]
    scale = np.float32(ewstar / ewvega)
    wldata = wave[zdata]
    fldata = (flux[zdata] - np.float32(1.0)).astype(np.float32)
    wlmodl = wvega[zmodl]
    flmodl = (fnvega[zmodl] - np.float32(1.0)).astype(np.float32)
    sflmodl = (flmodl * np.float32(scale)).astype(np.float32)
    rfldata = I.interpol(fldata, wldata, wlmodl).astype(np.float32)
    fftdata = I.idl_fft(rfldata).astype(np.complex64)
    fftmodl = I.idl_fft(sflmodl).astype(np.complex64)
    fftkern = fftdata / fftmodl
    dwmodl = (wlmodl.max() - wlmodl.min()) / np.float32(nmodl - 1)
    freq = np.arange(nmodl // 2, dtype=np.float32) / (np.float32(nmodl) * dwmodl)
    freq = np.concatenate([freq, -1 * freq[1:nmodl // 2][::-1], [0.0]]).astype(np.float32)
    _, datcoef = I.mpfitpeak(freq, (fftdata * fftdata).real.astype(np.float32), nterms=4)
    _, modcoef = I.mpfitpeak(freq, (fftmodl * fftmodl).real.astype(np.float32), nterms=4)
    sfft = min(modcoef[2], datcoef[2])
    f32 = np.float32
    sfft = f32(sfft)
    fwhm = f32(f32(2.0) * sfft * np.sqrt(f32(2.0) * I.alog(f32(2.0))))   # FLOAT
    fmax = f32(f32(twf) * fwhm)
    wind = (f32(1.0) / (f32(1.0) + I.ipow(np.abs(freq / fmax), 10))).astype(f32)
    fftkern = fftkern * wind
    kernel = I.idl_fft(fftkern, inverse=True).real.astype(np.float32)
    nk = kernel.size
    kernel = np.roll(kernel, nk // 2)
    kernel = (kernel / I.ftotal(kernel)).astype(np.float32)
    zv = np.where((wvega >= wabsmin) & (wvega <= wabsmax))[0]
    cvega = I.convol(((fnvega[zv] - f32(1.0)) * scale).astype(f32), kernel, I.ftotal(kernel))
    rcvega = I.interpolv(cvega, wvega[zv], wldata)
    ratio = ((fldata + f32(1.0)) / (rcvega + f32(1.0))).astype(f32)
    maxdev = float(np.max(np.abs(ratio - f32(1.0))))
    rmsdev = float(I.moments(ratio)[2])
    return dict(kernel=kernel, wline=float(wline), scale=scale, maxdev=maxdev,
                rmsdev=rmsdev, plot=dict(w=wldata, data=fldata,
                                         vega=I.interpol(fnvega - 1.0, wvega, wldata),
                                         conv=rcvega, resid=ratio - 1.0))


def normalize_region(wave, flux, normreg, degree=1):
    """XMC_CONKERN_NORMSPEC: divide by a robust polynomial through the windows."""
    gb = np.zeros(wave.size, dtype=int)
    for r in normreg:
        idx = np.sort(I.tabinv(wave, np.asarray(r, dtype=np.float32)))
        gb[I.iround(idx[0]):I.iround(idx[1]) + 1] = 1
    c, _ = I.robustpoly1d(wave, flux, degree, 3, 0.01, igoodbad=gb)
    return (flux / I.poly(wave, c)).astype(np.float32), c


def conkernel(st, order=6, normreg=None, normdegree=1, linereg=None, twf=10.,
              method="Deconvolution", ipcoeffs=None):
    """
    XMC_CONKERN + XTELLCOR_CONKERNEL. `normreg` are continuum windows
    around the line, `linereg` the line window, both in microns, chosen in
    the IDL by pressing n / s on the plot.
    """
    f32 = np.float32
    vdisp = f32(2.96495e-4) if st["obsmode"] == "LowRes15" else f32(1.49774e-5)   # FLOAT literals
    if method == "IP":
        # XTELLCOR_CONKERNEL picks the coefficients by slit width; the
        # instrument's IP_coefficients.dat holds five widths (0.3, 0.5, 0.8,
        # 1.6, 3.0 arcsec; for TS4_SOAR a dummy copy of uSpeX's), and any
        # other width stops XTELLCOR ("CASE statement found no matches").
        widths = [0.3, 0.5, 0.8, 1.6, 3.0]
        sw_ = st["slitw_arc"]
        hit = [k for k, w in enumerate(widths) if sw_ == np.float64(np.float32(w))]
        if not hit:
            raise ValueError(
                f"The IP method has instrument-profile coefficients only for 0.3, 0.5, 0.8, 1.6 "
                f"and 3.0 arcsec slits; XTELLCOR stops for this {sw_:g} arcsec slit. "
                "Use the Deconvolution method.")
        if ipcoeffs is None:
            from triplespec_reducer.instrument import INSTRUMENT_DIR
            ipcoeffs = np.loadtxt(os.path.join(INSTRUMENT_DIR, "IP_coefficients.dat"),
                                  comments="#")[hit[0]]
        p = np.asarray(ipcoeffs, dtype=float)
        kernels = {}
        for i, o in enumerate(st["stdorders"]):
            vkw = st["stdfwhm"][i] / vdisp
            nk = I.iround(10. * vkw)
            if not (nk % 2):
                nk += 1
            kvx = (np.arange(nk, dtype=f32) - f32(int(nk) // 2)) * vdisp
            ksx = kvx / st["stddisp"][i]
            ip = (_erf((ksx + p[1] - p[0]) / p[2]) - _erf((ksx - p[1] - p[0]) / p[2]))
            kernels[o] = ip / I.itotal(ip)
        st.update(kernels=kernels, vshift=f32(0.0), scale=f32(1.0), method="IP",
                  maxdev=np.nan, rmsdev=np.nan)
        return st
    idx = st["stdorders"].index(order)
    sw = st["std"]["spec"][idx, 0]
    sf = st["std"]["spec"][idx, 1]
    keep = _nantrim(sw, 2)
    wave = sw[keep].astype(np.float32)
    flux = sf[keep].astype(np.float32)
    nflux, ncoef = normalize_region(wave, flux, normreg, normdegree)
    z = (wave > linereg[0]) & (wave < linereg[1])
    # The 'fc2vega' XMC_CONKERN carries is fcvega (see the module notes).
    vshift, redshift = vegacorr(wave[z], nflux[z], st["wvega"], st["fvega"], st["cfvega"])
    wvega = (st["wvega"] * np.float32(1.0 + redshift)).astype(np.float32)
    fnvega = (st["fvega"] / st["cfvega"]).astype(np.float32)
    vc = vegaconv(wave, nflux, wvega, st["cfvega"], fnvega, linereg[0], linereg[1], twf=twf)
    kernel = vc["kernel"]
    kernels = resample_kernel(kernel, st["stddisp"], idx, vdisp, st["stdorders"])
    st.update(kernels=kernels, vshift=vshift, scale=vc["scale"], maxdev=vc["maxdev"],
              rmsdev=vc["rmsdev"], method="Deconvolution", kernel_order=order,
              kernel_plot=vc["plot"], normreg=normreg, linereg=linereg,
              normcoeffs=ncoef, wline=vc["wline"])
    return st


def resample_kernel(kernel, stddisp, idx, vdisp, orders):
    """
    XTELLCOR_CONKERNEL after XMC_CONKERN: the kernel, made on Vega's
    velocity grid, resampled to each order's pixels. Types as the IDL:
    GAUSSFIT of the FLOAT kernel, the pixel grids DOUBLE (the dispersions
    are), each kernel normalized by its FLOAT TOTAL.
    """
    kx_vegapix = np.arange(kernel.size, dtype=np.float32)
    _, a = I.gaussfit(kx_vegapix, kernel, nterms=3)
    kx_spexpix = (kx_vegapix - np.float32(a[1])) * vdisp / np.float64(stddisp[idx])   # TOTAL of a DOUBLE scalar
    out = {}
    for i, o in enumerate(orders):
        kx_spexwave = kx_spexpix * np.float64(stddisp[i])
        npix = (kx_spexwave.max() - kx_spexwave.min()) / vdisp
        if np.fmod(npix, 2) == 0:            # IDL: if not npix mod 2 (MOD binds first)
            npix = npix + 1
        kx_vegawave = (np.arange(I.itrunc(npix), dtype=np.float32) - npix / 2.) * vdisp
        nkern = I.linterp(kx_spexwave, kernel, kx_vegawave)
        out[o] = (nkern / I.ftotal(nkern)).astype(np.float32)
    return out


def _erf(x):
    from scipy.special import erf
    return erf(x)


def _nantrim(a, flag):
    from triplespec_reducer.wavecal import nantrim
    return nantrim(a, flag)


# ------------------------------------------------------------- line scales

def initial_points(st):
    """XMC_SCALELINES_INITCOMMON: a control point at each end and each H I line."""
    pts = {}
    for i, o in enumerate(st["stdorders"]):
        w = st["std"]["spec"][i, 0]
        lo, hi = np.nanmin(w), np.nanmax(w)
        arr = [[lo, st["scale"]]]
        for h in st["hlines"][(st["hlines"] > lo) & (st["hlines"] < hi)]:
            arr.append([h, st["scale"]])
        arr.append([hi, st["scale"]])
        arr = np.array(arr)
        pts[o] = arr[np.argsort(arr[:, 0], kind="stable")]
    st["cpoints"] = pts
    st.setdefault("cutreg", {o: [] for o in st["stdorders"]})
    return pts


def scale_spectrum(st, order, tension=None, interactive=False):
    """
    The H I line scales of one order, a spline under tension through the
    control points.

    As XMC_SCALELINES returns them to XTELLCOR (default): an FLTARR over the
    whole wavelength column of the standard, NaN padding included. With
    `interactive`, as XMC_SCALELINES_MAKESCALE makes them for its own
    preview: double, over the NaN-trimmed column. Two control points
    (no H I lines in the order) give 1.0 everywhere.
    """
    i = st["stdorders"].index(order)
    wave = st["std"]["spec"][i, 0]
    if interactive:
        wave = wave[_nantrim(wave, 2)]
    cp = st["cpoints"][order]
    tension = st.get("tension", np.float32(10.)) if tension is None else tension
    if cp.shape[0] == 2:
        return np.ones(wave.size, dtype=np.float32)
    sc = I.spline(cp[:, 0], cp[:, 1], wave, tension)
    return sc if interactive else sc.astype(np.float32)


def mktellspec(std_wave, std_flux, std_error, std_mag, std_bminv, kernel, scales,
               wvega, fvega, cfvega, cf2vega, vshift):
    """MC_MKTELLSPEC: the telluric correction spectrum of one order."""
    wmin, wmax = np.nanmin(std_wave), np.nanmax(std_wave)
    zv = np.where((wvega > wmin) & (wvega < wmax))[0]
    nkern = kernel.size
    nadd = I.iround(nkern / 2.)
    idx = (np.arange(zv.size + 2 * nadd) + (zv[0] - nadd)).astype(int)
    trim = _nantrim(std_wave, 2)
    rscales = I.linterp(std_wave[trim], scales[trim], wvega[idx])
    nfvconv = I.convol((fvega[idx] / cfvega[idx] - 1.0) * rscales, kernel) + 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        cfvconv = I.convol((cfvega[idx] / cf2vega[idx] - 1.0) * rscales, kernel) + 1.0
    cfvconv = cfvconv * cf2vega[idx]
    fvconv = nfvconv * cfvconv
    f32 = np.float32
    rshift = f32(f32(1.0) + (f32(vshift) / f32(2.99792458E5)))
    swvega = (wvega * rshift).astype(f32)
    rfvconv = I.linterp(swvega[idx], fvconv, std_wave)
    vegamag, vegabmv = f32(0.03), f32(0.00)
    magscale = I.fpow(f32(10.0), f32(f32(-0.4) * (f32(std_mag) - vegamag)))
    ebmv = max(f32(f32(std_bminv) - vegabmv), f32(0.0))
    avred = f32(f32(3.10) * ebmv)
    fvred = redden(std_wave, rfvconv, ebmv) * f32(magscale * I.fpow(f32(10.0), f32(f32(0.4) * avred)))
    scvega = fvred
    z = np.isfinite(std_flux)
    wa, fa, ea, fv = std_wave[z], std_flux[z], std_error[z], fvred[z]
    z = fa > 0.0
    wa, fa, ea, fv = wa[z], fa[z], ea[z], fv[z]
    t = fv / fa
    te = np.sqrt(I.sq((fv / I.sq(fa))) * I.sq(ea))
    tellcor, tellerr = interpspec(wa, t, std_wave, iyerror=te)
    return tellcor, tellerr, scvega


def redden(wave, flux, ebmv, rv=None):
    """
    MC_REDDEN: CCM between 0.3 and 3.3 um (Rieke beyond, not reached here).
    The wavelengths go to CCM_UNRED in DOUBLE Angstroms (wave*1.0d4); the
    flux keeps its type.
    """
    f32 = np.float32
    rvccm = f32(3.1) if rv is None else rv
    out = np.array(flux)
    wave = np.asarray(wave)
    zc = (wave > f32(0.3)) & (wave < f32(3.3))
    if zc.any():
        out[zc] = ccm_unred(wave[zc].astype(float) * 1.0e4, out[zc], f32(-1.0) * ebmv, rvccm)
    return out


def ccm_unred(wave, flux, ebv, r_v=None):
    """
    astrolib CCM_UNRED (Cardelli, Clayton & Mathis 1989; O'Donnell 1994),
    in IDL's types: x = 10000./wave in the wavelengths' precision, the
    coefficients FLOAT literals, a and b FLTARRs, A_V and A_lambda FLOAT,
    10.^(0.4 A_lambda) a FLOAT power.
    """
    f32 = np.float32
    r_v = f32(3.1) if r_v is None else r_v
    x = 10000. / np.asarray(wave)
    a = np.zeros(x.size, dtype=f32)
    b = np.zeros(x.size, dtype=f32)
    g = (x > f32(0.3)) & (x < f32(1.1))
    a[g] = f32(0.574) * I.fpow(x[g], f32(1.61))
    b[g] = f32(-0.527) * I.fpow(x[g], f32(1.61))
    g = (x >= f32(1.1)) & (x < f32(3.3))
    y = x[g] - f32(1.82)
    c1 = np.array([1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505], f32)
    c2 = np.array([0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347], f32)
    a[g] = I.poly(y, c1)
    b[g] = I.poly(y, c2)
    a_v = f32(r_v * ebv)
    a_lambda = (a_v * (a + b / f32(r_v))).astype(f32)
    return (np.asarray(flux) * I.fpow(f32(10.), (f32(0.4) * a_lambda).astype(f32))).astype(
        np.result_type(np.asarray(flux).dtype, f32))


def estimate_scale(st, order, ereg):
    """
    XMC_SCALELINES_ESTIMATESCALE: fit a Lorentzian to the line left in the
    telluric spectrum over `ereg` and set the control point(s) there to
    the scale that removes it.
    """
    std_wave, std_flux, _, tel = _interactive_telluric(st, order)
    tel_flux = (np.float32(1.) / tel).astype(np.float32)
    ereg = np.asarray(ereg, dtype=np.float32)
    z = (std_wave > ereg[0]) & (std_wave < ereg[1])
    sel = (st["hlines"] > ereg[0]) & (st["hlines"] < ereg[1])
    if not sel.any():
        return None
    wguess = st["hlines"][sel]
    idx = np.atleast_1d(I.tabinv(std_wave, wguess))
    idlo = I.tabinv(std_wave, ereg[0])
    idhi = I.tabinv(std_wave, ereg[1])
    ti = lambda r: tel_flux[I.itrunc(r)]
    fguess = np.array([ti(r) for r in idx]) - np.float32(0.5) * (ti(idlo) + ti(idhi))
    sguess = (ereg[1] - ereg[0]) / np.float32(5.)
    cguess = I.poly_fit1d(std_wave[z], tel_flux[z], 1)
    mguess = (ti(idlo) - ti(idhi)) / (std_wave[I.itrunc(idlo)] - std_wave[I.itrunc(idhi)])
    # ESTIMATES=[fguess,wguess,sguess,cguess,mguess]: cguess has two terms, so
    # MPFITPEAK (NTERMS=5) starts from [fguess, wguess, sguess, c0, c1] for one
    # line; mguess is never used.
    est = np.concatenate([fguess, wguess, [sguess], cguess, [mguess]]).astype(float)
    _, g = I.mpfitpeak(std_wave[z], tel_flux[z], nterms=5, lorentzian=True, estimates=est[:5])
    st["estimate_fit"] = g
    f32 = np.float32
    wvega = (st["wvega"] * (f32(1) + st["vshift"] / f32(2.99792458E5))).astype(f32)
    wmin, wmax = np.nanmin(std_wave), np.nanmax(std_wave)
    zv = np.where((wvega > wmin) & (wvega < wmax))[0]
    kern = st["kernels"][order]
    nadd = I.iround(kern.size / 2.)
    vi = (np.arange(zv.size + 2 * nadd) + (zv[0] - nadd)).astype(int)
    rvred, vegabmv, vegamag = f32(3.10), f32(0.00), f32(0.03)
    magscale = I.fpow(f32(10.0), f32(f32(-0.4) * (f32(st["vmag"]) - vegamag)))
    ebmv = max(f32(f32(st["bmag"] - st["vmag"]) - vegabmv), f32(0.0))
    avred = f32(rvred * ebmv)
    redfact = I.fpow(f32(10.0), f32(f32(0.4) * avred))
    magfact = np.full(std_wave.size, f32(redfact * magscale), dtype=f32)
    nfvconv = I.convol((st["fvega"][vi] / st["cfvega"][vi] - f32(1.0)).astype(f32), kern)
    rnfvconv = I.linterp(wvega[vi], nfvconv, std_wave)
    rfcvega = I.linterp(wvega, st["cfvega"], std_wave)
    rrfcvega = redden(std_wave, rfcvega, ebmv)
    s = (_poly_smooth(std_flux, 5) / (magfact * I.poly(std_wave, g[3:5]) * rrfcvega)) - f32(1.0)
    s = s / rnfvconv
    cp = st["cpoints"][order]
    zz = (cp[:, 0] > ereg[0]) & (cp[:, 0] < ereg[1])
    if not zz.any():
        return None
    ii = I.tabinv(std_wave, cp[zz, 0])
    cp[zz, 1] = s[I.itrunc(np.atleast_1d(ii))]          # stored in the control-point array's type
    return cp[zz, 1]


def _poly_smooth(y, width):
    """
    astrolib POLY_SMOOTH(y, width): the same filter as Spextool's
    SAVITZKY_GOLAY (degree 2, LUDC/LUSOL) applied with CONVOL /EDGE_TRUNCATE.
    """
    from triplespec_reducer.clean import savitzky_golay
    return savitzky_golay(np.asarray(y), width, degree=2)


def _interactive_telluric(st, order):
    """
    XMC_SCALELINES_MAKETELLURIC: the telluric spectrum XMC_SCALELINES shows
    and 'estimate scale' works on: the NaN-trimmed standard, the double
    interactive scales, cut regions on the flux only.
    """
    i = st["stdorders"].index(order)
    s_ = st["std"]["spec"][i]
    keep = _nantrim(s_[0], 2)
    w, f, e = s_[0][keep], s_[1][keep], s_[2][keep]
    tc, _, _ = mktellspec(w, f, e, st["vmag"], st["bmag"] - st["vmag"], st["kernels"][order],
                          scale_spectrum(st, order, interactive=True), st["wvega"], st["fvega"],
                          st["cfvega"], st["cf2vega"], st["vshift"])
    for (w0, w1) in st.get("cutreg", {}).get(order, []):
        ix = I.iround(I.tabinv(w, np.array([w0, w1], dtype=np.float32)))
        c = I.poly_fit1d(np.array([w[ix[0]], w[ix[1]]]), np.array([tc[ix[0]], tc[ix[1]]]), 1)
        tc[ix[0]:ix[1] + 1] = I.poly(w[ix[0]:ix[1] + 1], c)
    return w, f, e, tc


def telluric_order(st, order):
    """
    XTELLCOR_CONTELLSPEC for one order: MC_MKTELLSPEC with the returned
    scales, then the cut regions bridged by straight lines.

    Two IDL behaviours are kept: the kernel is taken by the position of the
    order among the object's orders, and the bridged errors are written at
    trimmed-array positions of the untrimmed array and then overwritten by
    the unbridged ones (so an error bridge only survives where the untrimmed
    array is NaN).
    """
    i = st["stdorders"].index(order)
    n = st["objorders"].index(order)
    s = st["std"]["spec"][i]
    kern = st["kernels"][st["stdorders"][n]]
    tc, te, sv = mktellspec(s[0], s[1], s[2], st["vmag"], st["bmag"] - st["vmag"],
                            kern, scale_spectrum(st, order), st["wvega"], st["fvega"],
                            st["cfvega"], st["cf2vega"], st["vshift"])
    cuts = st.get("cutreg", {}).get(order, [])
    if cuts:
        nn = _nantrim(s[0], 3)
        tw = s[0][nn]
        ttc, tte = tc[nn].copy(), te[nn].copy()
        for (w0, w1) in cuts:
            ix = I.iround(I.tabinv(tw, np.array([w0, w1], dtype=np.float32)))
            x = np.array([tw[ix[0]], tw[ix[1]]])
            c = I.poly_fit1d(x, np.array([ttc[ix[0]], ttc[ix[1]]]), 1)
            ce = I.poly_fit1d(x, np.array([tte[ix[0]], tte[ix[1]]]), 1)
            ttc[ix[0]:ix[1] + 1] = I.poly(tw[ix[0]:ix[1] + 1], c)
            te[ix[0]:ix[1] + 1] = I.poly(tw[ix[0]:ix[1] + 1], ce)
        tc[nn] = ttc
        te[nn] = tte
    return tc, te, sv


def contellspec(st, units="ergs s-1 cm-2 A-1"):
    """XTELLCOR_CONTELLSPEC: the telluric spectrum for every object order."""
    # copies of the standard (FLOAT), as XTELLCOR_CONTELLSPEC makes them
    tell = np.zeros((len(st["objorders"]), 4, st["std"]["spec"].shape[2]), dtype=st["std"]["spec"].dtype)
    vega = np.zeros_like(tell)
    for n, o in enumerate(st["objorders"]):
        i = st["stdorders"].index(o)
        tc, te, sv = telluric_order(st, o)
        w = st["std"]["spec"][i, 0]
        code = UNITS[units]
        if code != 1:                          # XTELLCOR_CHANGEUNITS
            tc, te = chfunits(w, tc, 0, 1, code, error=te)
            sv, _ = chfunits(w, sv, 0, 1, code)
        tell[n] = st["std"]["spec"][i]
        tell[n, 1], tell[n, 2] = tc, te
        vega[n] = st["std"]["spec"][i]
        vega[n, 1], vega[n, 2] = sv, 1.0
    st["tellspec"], st["vegaspec"] = tell, vega
    st["units"] = units
    return tell


# ------------------------------------------------------------- shift

def findshift(objwave, objflux, telflux, srange, initshift=0.0):
    """
    XMC_FINDSHIFT_AUTOFIND: the shift of the telluric spectrum, in pixels,
    that leaves the least structure in object x telluric over `srange`.
    """
    z = (objwave > srange[0]) & (objwave < srange[1])
    shifts = np.arange(151, dtype=np.float32) / np.float32(50.) + np.float32(initshift - 1.5)
    rms = np.zeros(151, dtype=np.float32)
    n = telflux.size
    x2 = np.arange(n, dtype=np.float32)
    for k, sh in enumerate(shifts):
        y, _ = interpspec(x2 + sh, telflux, x2)
        div = objflux * y
        div = div / I.median(div, even=True)
        rms[k] = I.moments(div[z])[2]
    mi = int(np.argmin(rms))
    if mi - 5 < 0 or mi + 5 > 150:
        # IDL stops here with a subscript error; say why instead.
        raise ValueError("The best shift lies at the edge of the +-1.5 pixel search; "
                         "choose a region with more telluric structure.")
    c = I.poly_fit1d(shifts[mi - 5:mi + 6], rms[mi - 5:mi + 6], 2)
    return float(-c[1] / 2. / c[2]), (shifts, rms)


def getshift(st, order, aperture=1, srange=None):
    """XTELLCOR_GETSHIFT for one object order and aperture."""
    zo = st["objorders"].index(order)
    k = zo * st["objnaps"] + (aperture - 1)
    ow = st["obj"]["spec"][k, 0]
    of = st["obj"]["spec"][k, 1]
    t = _nantrim(ow, 2)
    ow, of = ow[t], of[t]
    zs = st["stdorders"].index(order)
    tw = st["tellspec"][zs, 0]
    tf = st["tellspec"][zs, 1]
    t = _nantrim(tw, 2)
    nf, _ = interpspec(tw[t], tf[t], ow)
    sh, curve = findshift(ow, of, nf, srange, st["shift"][zo, aperture - 1])
    st["shift"][zo, aperture - 1] = sh
    return sh, curve


# ------------------------------------------------------------- output

def corrected(st, propstdunc=True):
    """XTELLCOR_WRITEFILE: the object corrected by the (shifted) telluric spectrum."""
    obj = st["obj"]["spec"]
    out = obj.copy()
    for i, o in enumerate(st["objorders"]):
        if o not in st["stdorders"]:
            continue
        for j in range(st["objnaps"]):
            k = i * st["objnaps"] + j
            t = st["tellspec"][i]
            nflux, nerror = interpspec(t[0], t[1], obj[k, 0], iyerror=t[2])
            nflag = interpflagspec(t[0], t[3], obj[k, 0])
            x = np.arange(obj.shape[2], dtype=np.float32)
            sh = st["shift"][i, j]
            nnflux, nnerror = interpspec(x + sh, nflux, x, iyerror=nerror)
            nnflag = interpflagspec(x + sh, nflag, x)
            out[k, 1] = nnflux * obj[k, 1]
            if propstdunc:
                out[k, 2] = np.sqrt(I.sq(nnflux) * I.sq(obj[k, 2]) + I.sq(obj[k, 1]) * I.sq(nnerror))
            else:
                out[k, 2] = nnflux * obj[k, 2]
            out[k, 3] = combflagstack([np.nan_to_num(obj[k, 3]).astype(np.uint8), nnflag])
    st["corrected"] = out
    return out


_ANG = "!5!sA!r!u!9 %!5!n"
NYTITLE = {                                  # XTELLCOR_CHANGEUNITS
    "ergs s-1 cm-2 A-1": "!5f!D!7k!N!5 (ergs s!E-1!N cm!E-2!N " + _ANG + "!E-1!N",
    "ergs s-1 cm-2 Hz-1": "!5f!D!7m!N!5 (ergs s!E-1!N cm!E-2!N Hz!E-1!N",
    "W m-2 um-1": "!5f!D!7k!N!5 (W m!E-2!N !7l!5m!E-1!N",
    "W m-2 Hz-1": "!5f!D!7m!N!5 (W m!E-2!N Hz!E-1!N",
    "Jy": "!5f!D!7m!N!5 (Jy",
}


def _fstr(v):
    """STRTRIM of a FLOAT, as IDL writes it into a history line."""
    return f"{float(np.float32(v)):#13.6g}".strip()


def write(path, st, tellpath=None, vegapath=None, std_name=None):
    """
    XTELLCOR_WRITEFILE: the corrected spectrum (the object's header kept as
    it is, a few keywords removed, XTELLCOR's added), and optionally the
    telluric and modified-Vega spectra with headers of their own. The
    values come from XTELLCOR's FLOAT state fields.
    """
    f32 = np.float32
    units = st.get("units", "ergs s-1 cm-2 A-1")
    ohdr, shdr = st["obj"]["header"], st["std"]["header"]
    std = I.fxpar(shdr, "OBJECT") if std_name is None else std_name
    bmag, vmag, vshift = f32(st["bmag"]), f32(st["vmag"]), f32(st["vshift"])
    stdairmass = f32(st["stdairmass"])
    dairmass = f32(st["dairmass"])
    deconv = st["method"] == "Deconvolution"
    maxdev, rmsdev = f32(st["maxdev"]), f32(st["rmsdev"])
    obj = os.path.splitext(os.path.basename(path))[0]
    nyt = NYTITLE.get(units, "")

    def dev_history(h):
        return (h + "The Vega model was modified using the deconvolution method.  The maximum "
                "deviation and rms deviation between the modified Vega model and the data over "
                "the kernel wavelength range are " + _fstr(f32(maxdev * f32(100.0))) + "% and "
                + _fstr(f32(rmsdev * f32(100.0))) + "%, respectively.")

    def add_history(h, label, text):
        I.addhistlabel(h, label)
        I.sxaddhist(h, text)

    if tellpath:
        th = fits.Header()
        I.fxaddpar(th, "IRAFNAME", os.path.basename(tellpath))
        I.fxaddpar(th, "ORDERS", I.fxpar(ohdr, "ORDERS"), I.fxcomment(ohdr, "ORDERS"))
        I.fxaddpar(th, "NORDERS", I.fxpar(ohdr, "NORDERS"), I.fxcomment(ohdr, "NORDERS"))
        I.fxaddpar(th, "A0VStd", std, " Tellluric Correction A0 V Standard")
        I.fxaddpar(th, "A0VBmag", bmag, " B-band magnitude")
        I.fxaddpar(th, "A0VVmag", vmag, " V-band magnitude")
        I.fxaddpar(th, "GRAT", I.fxpar(shdr, "MODENAME"), I.fxcomment(shdr, "MODENAME"))
        I.fxaddpar(th, "NAPS", 1, " Number of apertures")
        I.fxaddpar(th, "AIRMASS", stdairmass, "Average airmass")
        I.fxaddpar(th, "DAIRMASS", dairmass, " Average of (std-obj) airmass")
        I.fxaddpar(th, "TELMETH", st["method"].strip(), " Telluric Correction Method")
        I.fxaddpar(th, "Vegadv", vshift, " Vega velocity shift in km s-1")
        hist = ("These telluric correction spectra were constructed from the spectra "
                + os.path.basename(st.get("std_path", "")) + ".  The velocity shift of Vega is "
                + _fstr(vshift) + " km s-1.  ")
        if deconv:
            I.fxaddpar(th, "TCMaxDev", f32(maxdev * 100), " The maxmimum % deviation of Vega-data")
            I.fxaddpar(th, "TCRMSDev", f32(rmsdev * 100), " The RMS % deviation of Vega-model")
            hist = dev_history(hist)
        else:
            hist += " The Vega model was modified using the IP method."
        I.fxaddpar(th, "XUNITS", I.fxpar(shdr, "XUNITS"), I.fxcomment(shdr, "XUNITS"))
        if "XTITLE" in shdr:
            I.fxaddpar(th, "XTITLE", I.fxpar(shdr, "XTITLE"), I.fxcomment(shdr, "XTITLE"))
        else:
            I.fxaddpar(th, "XTITLE", "!7k!5 (" + str(I.fxpar(shdr, "XUNITS")).strip() + ")", "IDL X Title")
        I.fxaddpar(th, "YUNITS", units.replace(" ", "") + "/DNs-1", " Units of the Y axis")
        I.fxaddpar(th, "YTITLE", nyt + " / DN s!U-1!N)", "IDL Y title")
        add_history(th, "Xtellcor History", hist)
        fits.PrimaryHDU(st["tellspec"], header=th).writeto(tellpath, overwrite=True)
    if vegapath:
        vh = fits.Header()
        I.fxaddpar(vh, "IRAFNAME", os.path.basename(vegapath))
        I.fxaddpar(vh, "ORDERS", I.fxpar(ohdr, "ORDERS"), I.fxcomment(ohdr, "ORDERS"))
        I.fxaddpar(vh, "NORDERS", I.fxpar(ohdr, "NORDERS"), I.fxcomment(ohdr, "NORDERS"))
        I.fxaddpar(vh, "OBJECT", "Convolved Vega Spectrum")
        # FXPAR of the standard's GRAT, which TripleSpec spectra do not have: 0
        I.fxaddpar(vh, "GRAT", I.fxpar(shdr, "GRAT") if "GRAT" in shdr else 0,
                   I.fxcomment(shdr, "GRAT") if "GRAT" in shdr else "")
        I.fxaddpar(vh, "NAPS", 1, " Number of apertures")
        hist = ("This Vega spectrum has been scaled to a V magnitude of " + _fstr(vmag)
                + ", shifted by " + _fstr(vshift) + " km s-1, convolved with the kernel and "
                "resampled onto the wavelength grid of " + os.path.basename(st.get("std_path", "")) + ".")
        if deconv:
            I.fxaddpar(vh, "TCMaxDev", f32(maxdev * 100), " The maxmimum % deviation of Vega-data")
            I.fxaddpar(vh, "TCRMSDev", f32(rmsdev * 100), " The RMS % deviation of Vega-data")
            hist = dev_history(hist + "  ")
        else:
            hist += " The Vega model was modified using the IP method."
        I.fxaddpar(vh, "XUNITS", I.fxpar(shdr, "XUNITS"), I.fxcomment(shdr, "XUNITS"))
        I.fxaddpar(vh, "XTITLE", I.fxpar(shdr, "XTITLE") if "XTITLE" in shdr else 0,
                   I.fxcomment(shdr, "XTITLE") if "XTITLE" in shdr else "")
        I.fxaddpar(vh, "YUNITS", units.replace(" ", ""), "Units of the Y axis")
        I.fxaddpar(vh, "YTITLE", nyt + ")", "IDL Y title")
        add_history(vh, "Xtellcor History", hist)
        fits.PrimaryHDU(st["vegaspec"], header=vh).writeto(vegapath, overwrite=True)
    hdr = ohdr.copy()
    for k in ("DIVISOR", "BEAM", "GRSTCNT", "CALMIR", "DIT", "OSF", "QTH_LAMP",
              "INC_LAMP", "IR_SRC", "ARG_SRC", "SHUTTER", "AIMAGE", "SKY"):
        hdr.remove(k, ignore_missing=True, remove_all=True)
    # a literal IRAFNAME keyword (FILENAME stays the object file's)
    I.fxaddpar(hdr, "IRAFNAME", os.path.basename(path),
               I.fxcomment(hdr, "IRAFNAME") if "IRAFNAME" in hdr else "")
    I.fxaddpar(hdr, "YUNITS", units.replace(" ", ""), "Units of the Y axis")
    I.fxaddpar(hdr, "YTITLE", nyt + ")", " IDL Y title")
    I.fxaddpar(hdr, "A0VStd", std, " Tellluric Correction A0 V Standard")
    I.fxaddpar(hdr, "A0VBmag", bmag, " B-band magnitude")
    I.fxaddpar(hdr, "A0VVmag", vmag, " V-band magnitude")
    I.fxaddpar(hdr, "DAIRMASS", dairmass, " Average of (std-obj) airmass")
    I.fxaddpar(hdr, "TELMETH", st["method"].strip(), " Telluric Correction Method")
    I.fxaddpar(hdr, "Vegadv", vshift, " Vega velocity shift in km s-1")
    tname = os.path.basename(tellpath) if tellpath else obj + "_tellspec.fits"
    hist = ("This spectrum was telluric corrected using the telluric correction spectra "
            + tname + ".  The velocity shift of Vega is " + _fstr(vshift) + " km s-1.")
    if deconv:
        I.fxaddpar(hdr, "TCMaxDev", f32(maxdev * 100), " The maxmimum % deviation of Vega-data")
        I.fxaddpar(hdr, "TCRMSDev", f32(rmsdev * 100), " The RMS % deviation of Vega-data")
        hist = dev_history(hist + "  ")
    else:
        hist += " The Vega model was modified using the IP method."
    I.fxaddpar(hdr, "Telfile", tname, " The telluric correction file")
    for j in range(st["objnaps"]):
        hist += (f"  The telluric correction spectra for aperture {j + 1:02d} were shifted by "
                 + ", ".join(_fstr(s_) for s_ in st["shift"][:, j]) + " pixels.")
    add_history(hdr, "Xtellcor History", hist)
    fits.PrimaryHDU(st["corrected"], header=hdr).writeto(path, overwrite=True)


def apply_telluric(tellspec_path, obj_path, shifts=None):
    """
    XTELLCOR_FINISH: correct another object with an existing telluric file.

    Reproduced as written, two differences from XTELLCOR included: the error
    is computed after the flux has been corrected, so its second term uses
    the corrected flux; and an object order with no telluric counterpart
    ends the loop (IDL: goto, cont), leaving the orders after it as they
    were. The pixel shift is applied on the telluric file's own grid.
    """
    with fits.open(tellspec_path) as h:
        tel = h[0].data.astype(float)
        thdr = h[0].header
    obj = read_spectrum(obj_path)
    spec = obj["spec"].copy()
    torders = [int(t) for t in str(thdr["ORDERS"]).split(",")]
    naps = obj["naps"]
    sh = np.zeros((len(obj["orders"]), naps)) if shifts is None else np.asarray(shifts, float)
    for i, o in enumerate(obj["orders"]):
        if o not in torders:
            break
        z = torders.index(o)
        for j in range(naps):
            k = i * naps + j
            nflux, nerror = interpspec(tel[z, 0], tel[z, 1], spec[k, 0], iyerror=tel[z, 2])
            nt, no = tel.shape[2], spec.shape[2]
            if no > nt:
                raise ValueError("XTELLCOR_FINISH cannot apply a telluric spectrum shorter than "
                                 "the object's (IDL stops with a size mismatch).")
            # x runs over the telluric grid; IDL clamps the object-length
            # flux's out-of-range subscripts to its last element.
            pad = lambda a: np.concatenate([a, np.full(nt - no, a[-1])]) if nt > no else a
            x = np.arange(nt, dtype=np.float32)
            nnflux, nnerror = interpspec(x + sh[i, j], pad(nflux), x, iyerror=pad(nerror))
            nnflux, nnerror = nnflux[:no], nnerror[:no]
            spec[k, 1] = nnflux * spec[k, 1]
            spec[k, 2] = np.sqrt(I.sq(nnflux) * I.sq(spec[k, 2]) + I.sq(spec[k, 1]) * I.sq(nnerror))
    return spec, thdr, obj["header"]
