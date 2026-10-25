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
    shift          XMC_FINDSHIFT          sub-pixel shift per order and aperture
    write          XTELLCOR_WRITEFILE     corrected, telluric and (Model) Vega files

and the two other ways XTELLCOR offers:

    any standard   XTELLCOR_BASIC         divide by the star, optionally times
                                          a blackbody of given T scaled to V
    existing file  XTELLCOR_FINISH        a telluric file applied to another object

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

def findshift(objwave, objflux, telflux, srange=None, initshift=0.0):
    """
    XMC_FINDSHIFT_AUTOFIND: the shift of the telluric spectrum, in pixels,
    that leaves the least structure in object x telluric over `srange`.

    As XMC_FINDSHIFT_SETUP does first, the object is divided by the maximum
    of its 4-pixel running median; without `srange` the region is the whole
    order (XMC_FINDSHIFT's default before a region is selected with 's').
    """
    objwave = np.asarray(objwave)
    objflux = np.asarray(objflux)
    smooth = I.median_filter(objflux, 4)
    objflux = objflux / np.asarray(np.nanmax(smooth)).astype(objflux.dtype)
    if srange is None:
        srange = [np.nanmin(objwave), np.nanmax(objwave)]
    srange = np.asarray(srange, dtype=np.float32)          # FLOAT state field
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

def set_shift(st, order, aperture, value):
    """
    The shift of one object order and aperture, as typed in XMC_FINDSHIFT's
    'Shift' field (a FLOAT). XTELLCOR keeps one shift per order and per
    aperture; set every aperture's own.
    """
    st["shift"][st["objorders"].index(order), aperture - 1] = np.float32(value)


def shifts_array(st, shifts):
    """
    `shifts` as XTELLCOR's FLTARR(norders, naps): an array of that shape, or
    a dict {(order, aperture): pixels} (apertures from 1; missing ones 0).
    """
    norders, naps = len(st["objorders"]), st["objnaps"]
    if isinstance(shifts, dict):
        out = np.zeros((norders, naps), dtype=np.float32)
        for (o, ap), v in shifts.items():
            out[st["objorders"].index(o), ap - 1] = np.float32(v)
        return out
    out = np.asarray(shifts, dtype=np.float32)
    if out.shape != (norders, naps):
        raise ValueError(f"shifts must be ({norders}, {naps}): one per order and aperture")
    return out


def corrected(st, propstdunc=True, shifts=None):
    """
    XTELLCOR_WRITEFILE: the object corrected by the (shifted) telluric spectrum.

    The shifts are st["shift"] (norders, naps), or `shifts` (see
    shifts_array), which then replace them. `propstdunc=False` (XTELLCOR's
    'Propagate Std Uncertainty', whose widget TS4_SOAR's XTELLCOR does not
    show; always on there) leaves the standard's error out.
    """
    if shifts is not None:
        st["shift"] = shifts_array(st, shifts)
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


def output_names(path):
    """
    The names XTELLCOR_WRITEFILE gives its three files for output name
    `obj` (path without .fits): obj.fits, obj_tellspec.fits ('Telluric',
    on by default) and obj_modVega.fits ('Model', off by default).
    """
    base = path[:-5] if path.endswith(".fits") else path
    return base + ".fits", base + "_tellspec.fits", base + "_modVega.fits"


def write(path, st, tellpath=None, vegapath=None, std_name=None, model=False):
    """
    XTELLCOR_WRITEFILE: the corrected spectrum (the object's header kept as
    it is, a few keywords removed, XTELLCOR's added), and optionally the
    telluric and modified-Vega spectra with headers of their own. The
    values come from XTELLCOR's FLOAT state fields.

    `model=True` ('Additional output: Model') writes the convolved, scaled
    Vega model to `vegapath`, by default XTELLCOR's name for it
    (output_names). Its IRAFNAME is XTELLCOR's: for a file named
    obj_modVega.fits it reads obj_modvega.fits.
    """
    if model and not vegapath:
        vegapath = output_names(path)[2]
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
        vname = os.path.basename(vegapath)
        if vname.endswith("_modVega.fits"):
            vname = vname[:-len("_modVega.fits")] + "_modvega.fits"     # strtrim(obj+'_modvega.fits')
        I.fxaddpar(vh, "IRAFNAME", vname)
        I.fxaddpar(vh, "ORDERS", I.fxpar(ohdr, "ORDERS"), I.fxcomment(ohdr, "ORDERS"))
        I.fxaddpar(vh, "NORDERS", I.fxpar(ohdr, "NORDERS"), I.fxcomment(ohdr, "NORDERS"))
        I.fxaddpar(vh, "OBJECT", "Convolved Vega Spectrum")
        # FXPAR of the standard's GRAT, which TripleSpec spectra do not have: 0;
        # a keyword FXPAR does not find leaves COMMENT= as it was, so the
        # comments are those the telluric file's block read (MODENAME's)
        cgrat = I.fxcomment(shdr, "GRAT") if "GRAT" in shdr else (
            I.fxcomment(shdr, "MODENAME") if tellpath else "")
        I.fxaddpar(vh, "GRAT", I.fxpar(shdr, "GRAT") if "GRAT" in shdr else 0, cgrat)
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
                   I.fxcomment(shdr, "XTITLE") if "XTITLE" in shdr else (
                       "IDL X Title" if tellpath else ""))
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


# ------------------------------------------------------------- an existing telluric file

def load_finish(tel_file, obj_file):
    """
    XTELLCOR_FINISH_LOADSPEC: a telluric file XTELLCOR wrote and an object,
    both as MC_READSPEC reads them (the object is not given a flag column),
    shifts 0, the ATRAN model for the object's resolving power, and the
    airmass difference (a FLOAT) with XTELLCOR's warning above 0.1.
    """
    tel = read_spectrum(tel_file)
    obj = read_spectrum(obj_file)
    aw, at = read_atran(obj["rp"])
    st = dict(tel=tel, obj=obj, tel_path=tel_file, telorders=tel["orders"],
              objorders=obj["orders"], objnaps=obj["naps"], awave=aw, atrans=at,
              dairmass=np.float32(tel["airmass"] - obj["airmass"]),
              shift=np.zeros((len(obj["orders"]), obj["naps"]), dtype=np.float32))
    st["airmass_warning"] = abs(tel["airmass"] - obj["airmass"]) > 0.1
    return st


def getshift_finish(st, order, aperture=1, srange=None):
    """XTELLCOR_FINISH_GETSHIFT: XMC_FINDSHIFT on the object and the telluric spectrum."""
    zo = st["objorders"].index(order)
    k = zo * st["objnaps"] + (aperture - 1)
    ow, of = st["obj"]["spec"][k, 0], st["obj"]["spec"][k, 1]
    t = _nantrim(ow, 2)
    ow, of = ow[t], of[t]
    z = st["telorders"].index(order)
    tw, tf = st["tel"]["spec"][z, 0], st["tel"]["spec"][z, 1]
    t = _nantrim(tw, 2)
    nf, _ = interpspec(tw[t], tf[t], ow)
    sh, curve = findshift(ow, of, nf, srange, st["shift"][zo, aperture - 1])
    st["shift"][zo, aperture - 1] = sh
    return sh, curve


def tellcor_finish(st, shifts=None):
    """
    XTELLCOR_FINISH_TELLCOR: the object times the shifted telluric spectrum.

    As written there: no flags are combined; the error is computed after the
    flux has been corrected, so its second term uses the corrected flux; an
    object order with no telluric counterpart ends the loop (goto, cont),
    leaving it and the orders after it as they were. The shift runs on the
    telluric file's pixel grid against the object-length spectrum: MC_INTERPSPEC
    keeps the pixels where both are finite, which truncates to the shorter.
    """
    if shifts is not None:
        st["shift"] = shifts_array(st, shifts)
    tel = st["tel"]["spec"]
    spec = st["obj"]["spec"].copy()
    naps = st["objnaps"]
    for i, o in enumerate(st["objorders"]):
        if o not in st["telorders"]:
            break
        z = st["telorders"].index(o)
        for j in range(naps):
            k = i * naps + j
            nflux, nerror = interpspec(tel[z, 0], tel[z, 1], spec[k, 0], iyerror=tel[z, 2])
            nt, no = tel.shape[2], spec.shape[2]
            if no > nt:
                raise ValueError("XTELLCOR_FINISH cannot apply a telluric spectrum shorter than "
                                 "the object's (IDL stops with a size mismatch).")
            x = np.arange(nt, dtype=np.float32)
            m = min(nt, no)
            nnflux, nnerror = interpspec((x + st["shift"][i, j])[:m], nflux[:m], x,
                                         iyerror=nerror[:m])
            nnflux, nnerror = nnflux[:no], nnerror[:no]
            spec[k, 1] = nnflux * spec[k, 1]
            spec[k, 2] = np.sqrt(I.sq(nnflux) * I.sq(spec[k, 2]) + I.sq(spec[k, 1]) * I.sq(nnerror))
    st["corrected"] = spec
    return spec


def _idlstr(v):
    """STRTRIM(v, 2) of a number FXPAR returned: FLOAT %#13.6g, DOUBLE %#16.8g, integers."""
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return (f"{float(v):#13.6g}" if I.isf32(v) else f"{float(v):#16.8g}").strip()


def write_finish(path, st):
    """
    XTELLCOR_FINISH_TELLCOR's file: the object's header, a few keywords
    removed, the telluric file's A0VStd, A0VBmag, A0VVmag, DAIRMASS,
    TELMETH, VEGADV, TELFILE (its IRAFNAME), TCMAXDEV and TCRMSDEV with
    their comments, YUNITS and YTITLE up to the '/' of the telluric file's,
    and the history under 'Xtellcor History'. IRAFNAME is the output name.
    """
    thdr = st["tel"]["header"]
    hdr = st["obj"]["header"].copy()
    for k in ("DIVISOR", "BEAM", "GRSTCNT", "CALMIR", "DIT", "OSF", "QTH_LAMP",
              "INC_LAMP", "IR_SRC", "ARG_SRC", "SHUTTER", "AIMAGE", "SKY"):
        hdr.remove(k, ignore_missing=True, remove_all=True)
    I.fxaddpar(hdr, "IRAFNAME", os.path.basename(path),
               I.fxcomment(hdr, "IRAFNAME") if "IRAFNAME" in hdr else "")
    yunits = [t for t in str(I.fxpar(thdr, "YUNITS")).split("/") if t][0]
    for k in ("A0VStd", "A0VBmag", "A0VVmag", "DAIRMASS"):
        I.fxaddpar(hdr, k, I.fxpar(thdr, k), I.fxcomment(thdr, k))
    telmeth = str(I.fxpar(thdr, "TELMETH")).strip()
    I.fxaddpar(hdr, "TELMETH", telmeth, I.fxcomment(thdr, "TELMETH"))
    vegadv = I.fxpar(thdr, "Vegadv")
    I.fxaddpar(hdr, "VEGADV", vegadv, I.fxcomment(thdr, "Vegadv"))
    telfile = I.fxpar(thdr, "IRAFNAME")
    I.fxaddpar(hdr, "TELFILE", telfile, I.fxcomment(thdr, "IRAFNAME"))
    hist = ("This spectrum was telluric corrected using the telluric correction spectra "
            + str(telfile) + ".  The velocity shift of Vega is " + _idlstr(vegadv) + " km s-1.")
    if telmeth == "Deconvolution":
        maxdev = I.fxpar(thdr, "TCMaxDev")
        I.fxaddpar(hdr, "TCMAXDEV", maxdev, I.fxcomment(thdr, "TCMaxDev"))
        rmsdev = I.fxpar(thdr, "TCRMSDev")
        I.fxaddpar(hdr, "TCRMSDEV", rmsdev, I.fxcomment(thdr, "TCRMSDev"))
        hist += ("  The Vega model was modified using the deconvolution method.  The maximum "
                 "deviation and rms deviation between the modified Vega model and the data over "
                 "the kernel wavelength range are " + _idlstr(maxdev) + "% and "
                 + _idlstr(rmsdev) + "%, respectively.  ")
    else:
        hist += "The Vega model was modified using the IP method.  "
    for j in range(st["objnaps"]):
        hist += (f"The telluric correction spectra for aperture {j + 1:02d} were shifted by "
                 + ", ".join(_fstr(s_) for s_ in st["shift"][:, j]) + " pixels.  ")
    I.fxaddpar(hdr, "YUNITS", yunits, "Units of the Y axis")
    ytitle = [t for t in str(I.fxpar(thdr, "YTITLE")).split("/") if t][0]
    I.fxaddpar(hdr, "YTITLE", ytitle + ")", " IDL Y title")
    I.addhistlabel(hdr, "Xtellcor History")
    I.sxaddhist(hdr, hist)
    fits.PrimaryHDU(st["corrected"], header=hdr).writeto(path, overwrite=True)


def apply_telluric(tellspec_path, obj_path, shifts=None):
    """
    XTELLCOR_FINISH in one call (see load_finish, tellcor_finish): returns
    the corrected spectrum, the telluric file's header and the object's.
    """
    st = load_finish(tellspec_path, obj_path)
    spec = tellcor_finish(st, shifts=shifts)
    return spec, st["tel"]["header"], st["obj"]["header"]


# ------------------------------------------------------------- any standard: XTELLCOR_BASIC

def load_basic(std_file, obj_file):
    """
    XTELLCOR_BASIC_LOADSPEC: a standard of any spectral type and an object,
    both through MC_CAASPECFLAGS, shifts 0, the ATRAN model for the
    object's resolving power and the airmass warning (> 0.1; the object's
    TCS_AM when its AIRMASS is 0).
    """
    std = read_spectrum(std_file)
    obj = read_spectrum(obj_file)
    for d in (std, obj):
        d["spec"] = caaspecflags(d["spec"])
    aw, at = read_atran(obj["rp"])
    st = dict(std=std, obj=obj, std_path=std_file, stdorders=std["orders"],
              objorders=obj["orders"], objnaps=obj["naps"], awave=aw, atrans=at,
              stdairmass=std["airmass"], objairmass=obj["airmass"],
              shift=np.zeros((len(obj["orders"]), obj["naps"]), dtype=np.float32))
    st["airmass_warning"] = abs(std["airmass"] - obj["airmass"]) > 0.1
    return st


def getshift_basic(st, order, aperture=1, srange=None):
    """XTELLCOR_BASIC_GETSHIFT: XMC_FINDSHIFT on the object and 1/standard."""
    zo = st["objorders"].index(order)
    k = zo * st["objnaps"] + (aperture - 1)
    ow, of = st["obj"]["spec"][k, 0], st["obj"]["spec"][k, 1]
    t = _nantrim(ow, 2)
    ow, of = ow[t], of[t]
    z = st["stdorders"].index(order)
    tw = st["std"]["spec"][z, 0]
    tf = np.float32(1.) / st["std"]["spec"][z, 1]
    t = _nantrim(tw, 2)
    nf, _ = interpspec(tw[t], tf[t], ow)
    sh, curve = findshift(ow, of, nf, srange, st["shift"][zo, aperture - 1])
    st["shift"][zo, aperture - 1] = sh
    return sh, curve


def planck(wave, temp):
    """
    astrolib PLANCK: the blackbody (erg s-1 cm-2 A-1) at `wave` (A) and
    `temp` (K), in the types the IDL computes it: the wavelength in cm in
    the input's precision, the exponent in double, the result stored in a
    FLOAT array (for float or integer wavelengths).
    """
    wave = np.asarray(wave)
    single = I.isf32(wave) or wave.dtype.kind in "iu"
    P = np.float32 if single else np.float64
    w = wave.astype(P) / P(1.E8) if single else wave / np.float64(np.float32(1.E8))
    c1, c2 = np.float64(3.7417749e-5), np.float64(1.4387687)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        val = c2 / w / temp
        good = val < np.log(np.finfo(np.float64).max)              # MACHAR(/DOUBLE).XMAX
        bb = np.zeros(w.shape, dtype=P)
        bb[good] = c1 / (I.ipow(w[good], 5) * (np.exp(val[good]) - 1.))
    return (bb * P(1.E-8)).astype(P)


def _basic_units(wave, spec, err, units):
    """XTELLCOR_BASIC_CHANGEUNITS: its own conversions (c a FLOAT literal) and Y title."""
    c = np.float32(2.99792458e+8)
    f = np.float32
    ang = "!5!sA!r!u!9 %!5!n"
    if units == "ergs s-1 cm-2 A-1":
        return spec, err, "!5f!D!7k!N!5 (ergs s!E-1!N cm!E-2!N " + ang + "!E-1!N)"
    if units == "ergs s-1 cm-2 Hz-1":
        k = f(1.0e-2) / c
        return spec * I.sq(wave) * k, err * I.sq(wave) * k, "!5f!D!7m!N!5 (ergs s!E-1!N cm!E-2!N Hz!E-1!N)"
    if units == "W m-2 um-1":
        return spec * f(10.), err * f(10.), "!5f!D!7k!N!5 (W m!E-2!N !7l!5m!E-1!N)"
    if units == "W m-2 Hz-1":
        k = f(1.0e-5) / c
        return spec * I.sq(wave) * k, err * I.sq(wave) * k, "!5f!D!7m!N!5 (W m!E-2!N Hz!E-1!N)"
    if units == "Jy":
        k = f(1.0e21) / c
        return spec * I.sq(wave) * k, err * I.sq(wave) * k, "!5f!D!7m!N!5 (Jy)"
    raise ValueError(f"units must be one of {list(UNITS)}")


def tellcor_basic(st, restore=True, temp=None, vmag=None, units="ergs s-1 cm-2 A-1", shifts=None):
    """
    XTELLCOR_BASIC_TELLCOR: the object divided by the (shifted) standard.

    With `restore` (XTELLCOR_BASIC's default, 'Restore Continuum: Yes') the
    ratio is multiplied by a blackbody of temperature `temp` (K) scaled to
    V magnitude `vmag` (3.46e-9 x 10^(-0.4 (V - 0.03)) erg s-1 cm-2 A-1 at
    5556 A), both FLOAT fields. The units conversion runs in either case,
    so without the restore a unitless ratio is converted too, as in IDL. An
    object order the standard lacks ends that aperture's loop (goto, cont):
    it and the later orders stay as loaded. Flags are not combined.
    """
    if shifts is not None:
        st["shift"] = shifts_array(st, shifts)
    if restore:
        if temp is None or vmag is None:
            raise ValueError("The blackbody needs a temperature and a V magnitude.")
        temp, vmag = np.float32(temp), np.float32(vmag)
    std = st["std"]["spec"]
    obj = st["obj"]["spec"].copy()
    naps = st["objnaps"]
    nytitle = ""
    for i in range(naps):
        for j, o in enumerate(st["objorders"]):
            k = j * naps + i
            if o not in st["stdorders"]:
                break
            z = st["stdorders"].index(o)
            sflux, serror = interpspec(std[z, 0], std[z, 1], obj[k, 0], iyerror=std[z, 2])
            one = np.float32(1.0)
            tflux = one / sflux
            terror = (one / I.sq(sflux)) * serror
            x = np.arange(sflux.size, dtype=np.float32)
            stflux, sterror = interpspec(x + st["shift"][j, i], tflux, x, iyerror=terror)
            cor = obj[k, 1] * stflux
            cerr = np.sqrt(I.sq(stflux) * I.sq(obj[k, 2]) + I.sq(obj[k, 1]) * I.sq(sterror))
            if restore:
                bbv = planck(5556, temp)
                scale = np.float32(3.46e-9) * I.fpow(10, np.float32(-0.4) * (vmag - np.float32(0.03))) / bbv
                bb = planck(obj[k, 0] * np.float32(10000.), temp) * scale
                cor = cor * bb
                cerr = cerr * bb
            cor, cerr, nytitle = _basic_units(obj[k, 0], cor, cerr, units)
            obj[k, 1] = cor
            obj[k, 2] = cerr
    st.update(corrected=obj, restore=restore, temp=temp, vmag=vmag, units=units, nytitle=nytitle)
    return obj


def write_basic(path, st, std_name=None):
    """
    XTELLCOR_BASIC_TELLCOR's file: the object's header with IRAFNAME (the
    output name), YUNITS and YTITLE (units and title; 'None' and '!5Ratio ()'
    without the restore) and the history under 'Xtellcor Basic History'.
    `std_name` is the standard's file name as typed (by default its base name).
    """
    stdfile = os.path.basename(st["std_path"]) if std_name is None else std_name
    hist = "The spectra were divided by the spectra of the standard star " + stdfile.strip() + "."
    if st["restore"]:
        hist += ("  The spectra were also multiplied by a  blackbody with a temperature of "
                 + _fstr(st["temp"]) + " K scaled to a V-band magnitude of "
                 + f"{float(st['vmag']):6.3f}" + ".")
        yunits, ytitle = st["units"], st["nytitle"]
    else:
        yunits, ytitle = "None", "!5Ratio ()"
        hist += " Therefore the spectrum is unitless."
    hdr = st["obj"]["header"].copy()
    I.fxaddpar(hdr, "IRAFNAME", os.path.basename(path).strip(),
               I.fxcomment(hdr, "IRAFNAME") if "IRAFNAME" in hdr else "")
    I.fxaddpar(hdr, "YUNITS", yunits, "Units of the Y axis")
    I.fxaddpar(hdr, "YTITLE", ytitle, "IDL Y Title")
    I.addhistlabel(hdr, "Xtellcor Basic History")
    I.sxaddhist(hdr, hist)
    fits.PrimaryHDU(st["corrected"], header=hdr).writeto(path, overwrite=True)
