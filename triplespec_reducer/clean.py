"""
The last, optional step: XCLEANSPEC.

    smoothing    Savitzky-Golay (the default; SAVITZKY_GOLAY) or a Gaussian
                 (MC_CONVOLVESPEC), the error propagated through either
    'f' fix      a region replaced by the straight line between its ends,
                 flag bit 2 (value 4) set
    'r' remove   a region set to NaN
"""

import numpy as np
from astropy.io import fits

from triplespec_reducer import idl as I


def savitzky_golay(data, width, degree=None, error=None, nl=None, nr=None, order=0):
    """
    SAVITZKY_GOLAY (Spextool's copy of POLY_SMOOTH with an error estimate).

    The filter is built in single precision with LU decomposition, as the
    IDL builds it, and applied with CONVOL /CENTER /EDGE_TRUNCATE. Returns
    the smoothed data, or (data, error) when `error` is given.
    """
    degree = 2 if degree is None else int(degree)
    order = max(min(order, degree - 1), 0)
    width = max(I.itrunc(width), 3)
    if nr is None and nl is None:
        nl = width // 2
        nr = width - nl - 1
    degree = max(degree, 2)
    nj = degree + 1
    nl, nr = max(nl, 0), max(nr, 0)
    nrl = nr + nl + 1
    if nrl <= degree:
        return data if error is None else (data, error)
    ata = np.zeros((nj, nj), dtype=np.float32)
    ata[0, 0] = 1
    iaj = np.add.outer(np.arange(nj), np.arange(nj))
    m1 = (-1.0) ** iaj
    for k in range(1, max(nr, nl) + 1):
        kij = np.float32(k) ** iaj
        if k <= min(nr, nl):
            ata = ata + (kij + kij * m1)
        elif k <= nr:
            ata = ata + kij
        elif k <= nl:
            ata = ata + kij * m1
    ata = ata.astype(np.float32)
    # LUDC and LUSOL in single precision (ATA is symmetric, so /COLUMN
    # changes nothing)
    lu, indx = I.ludc(ata)
    bmat = np.zeros((min(degree, 4), nj), dtype=np.float32)
    for m in range(min(degree - 1, 3) + 1):
        b = np.zeros(nj, dtype=np.float32)
        b[m] = 1
        bmat[m] = I.lusol(lu, indx, b)
    kvec = np.concatenate([(-np.arange(nl) - 1)[::-1], [0], np.arange(nr) + 1]).astype(np.float32)
    kmat = np.ones((nj, nrl), dtype=np.float32)
    for m in range(1, degree + 1):
        kmat[m] = kmat[m - 1] * kvec
    # Kmat # Bmat: a float running sum over the polynomial terms
    coefs = np.array([[I._seqdot(kmat[:, r], bmat[c], np.float32) for c in range(bmat.shape[0])]
                      for r in range(nrl)], dtype=np.float32)    # [nrl, nderiv]
    if nr > nl:
        coefs = np.vstack([np.zeros((nr - nl, coefs.shape[1]), np.float32), coefs])
    elif nl > nr:
        coefs = np.vstack([coefs, np.zeros((nl - nr, coefs.shape[1]), np.float32)])
    fc = coefs[:, order]
    out = I.convol(data, fc, edge_truncate=True)
    if error is None:
        return out
    err = np.sqrt(I.convol(I.sq(np.asarray(error)), I.sq(fc), edge_truncate=True))
    return out, err


def convolvespec(wave, flux, fwhm, error=None, khwidth=4, edge_truncate=False):
    """MC_CONVOLVESPEC: a Gaussian of `fwhm` (in units of `wave`)."""
    wave = np.asarray(wave)
    ndat = wave.size
    conv = (ndat - 1) / (wave[ndat - 1] - wave[0])
    # the FWHM is a FLOAT field; times a float conversion it stays float
    fwhm = np.result_type(wave.dtype, np.float32).type(np.float32(fwhm) * conv)
    npix = 2 * I.iround(fwhm / 2 * khwidth)
    if npix % 2 != 1:
        npix += 1
    kernel = I.psf_gaussian_1d(npix, fwhm)
    flux = np.asarray(flux)
    z = np.isfinite(flux)
    cflux = flux.copy()
    tflux = I.convol(flux[z], kernel, edge_truncate=edge_truncate)
    h = npix // 2
    if error is None:
        if not edge_truncate:
            tflux[:h] = np.nan
            tflux[z.sum() - h:] = np.nan
        cflux[z] = tflux
        return cflux, None
    error = np.asarray(error)
    cerror = error.copy()
    terr = np.sqrt(I.convol(I.sq(error[z]), I.sq(kernel), edge_truncate=edge_truncate))
    if not edge_truncate:
        tflux[:h] = np.nan
        tflux[z.sum() - h:] = np.nan
        terr[:h] = np.nan
        terr[z.sum() - h:] = np.nan
    cflux[z] = tflux
    cerror[z] = terr
    return cflux, cerror


def read(path):
    """
    How XCLEANSPEC reads a spectrum (MC_READSPEC): the stored type, and the
    merged [npix, 4] form as one order. Returns (cube, header, two_d).
    """
    with fits.open(path) as h:
        d = h[0].data
        spec = d.astype(d.dtype.newbyteorder("="))
        hdr = h[0].header
    two_d = spec.ndim == 2
    return (spec[None] if two_d else spec), hdr, two_d


def write(path, cube, hdr, history="", two_d=True):
    """XCLEANSPEC_WRITEFILE: the input header, the smoothing history if any."""
    hdr = hdr.copy()
    if history:
        I.addhistlabel(hdr, "Xcleanspec History")
        I.sxaddhist(hdr, history)
    fits.PrimaryHDU(cube[0] if two_d else cube, header=hdr).writeto(path, overwrite=True)


def smooth(spec, method="Savitzky-Golay", sgwin=None, sgdeg=2, fwhm=None, slitw_pix=2.4):
    """
    XCLEANSPEC_SMOOTHSPEC on a [n, 4, npix] cube. The defaults are the
    panel's: SG width slitw_pix x 1.5 (x 2.5 when the slit is 2 pixels or
    narrower), degree 2; Gaussian FWHM = slitw_pix pixels.
    """
    out = spec.copy()
    sgwin = float(f"{slitw_pix * (1.5 if slitw_pix > 2 else 2.5):.1f}") if sgwin is None else sgwin
    fwhm = slitw_pix if fwhm is None else fwhm
    for i in range(spec.shape[0]):
        z = np.isfinite(spec[i, 1])
        if method == "Savitzky-Golay":
            f, e = savitzky_golay(spec[i, 1][z], sgwin, degree=sgdeg, error=spec[i, 2][z])
            out[i, 1, z] = f
            out[i, 2, z] = e
        else:
            # FINDGEN(stop-start+1), the FWHM a FLOAT field
            x = np.arange(spec.shape[2], dtype=np.float32)
            f, e = convolvespec(x, spec[i, 1][z], np.float32(fwhm), error=spec[i, 2][z])
            out[i, 1, z] = f
            out[i, 2, z] = e
    # XCLEANSPEC_SMOOTHSPEC's history: the FLOAT field values as STRTRIM prints them
    fs = lambda v: f"{float(np.float32(v)):#13.6g}".strip()
    lead = "These spectra have" if spec.shape[0] > 1 else "This spectrum has"
    history = (f"{lead} been convolved with a Savitzky-Golay filter of width {fs(sgwin)} "
               f"pixels and degree {fs(sgdeg)}." if method == "Savitzky-Golay"
               else f"{lead} been convolved with a Gaussian of FWHM={fs(fwhm)} pixels.")
    return out, history


def _mc_tabinv(xarr, x):
    """Spextool MC_TABINV: TABINV on the NaN-trimmed array, offset back."""
    from triplespec_reducer.wavecal import nantrim
    k = nantrim(xarr, 2)
    return I.tabinv(xarr[k], x) + k[0]


def fix_region(spec, idx, w0, w1, flux_too=True):
    """
    'f': a straight line between the two ends of [w0, w1]; flag +4. The
    region comes from the cursor (FLOAT), sorted, through MC_TABINV.

    In XCLEANSPEC the line is drawn into the working copy and kept only
    when 's' is pressed; here every edit is kept, as if 's' followed each.
    """
    out = spec.copy()
    i0, i1 = I.iround(_mc_tabinv(spec[idx, 0], np.sort(np.array([w0, w1], dtype=np.float32))))
    x = spec[idx, 0, [i0, i1]]
    if flux_too:
        c = I.poly_fit1d(x, spec[idx, 1, [i0, i1]], 1)
        out[idx, 1, i0:i1 + 1] = I.poly(spec[idx, 0, i0:i1 + 1], c)
    c = I.poly_fit1d(x, spec[idx, 2, [i0, i1]], 1)
    out[idx, 2, i0:i1 + 1] = I.poly(spec[idx, 0, i0:i1 + 1], c)
    if spec.shape[1] > 3:
        out[idx, 3, i0:i1 + 1] = out[idx, 3, i0:i1 + 1] + 4
    return out


def remove_region(spec, idx, w0, w1, flux_too=True):
    """'r': the region set to NaN (cursor positions FLOAT, through MC_TABINV)."""
    out = spec.copy()
    i0, i1 = sorted(I.iround(_mc_tabinv(spec[idx, 0], np.array([w0, w1], dtype=np.float32))))
    i0, i1 = max(i0, 0), min(i1, spec.shape[2] - 1)
    if flux_too:
        out[idx, 1, i0:i1 + 1] = np.nan
    out[idx, 2, i0:i1 + 1] = np.nan
    return out
