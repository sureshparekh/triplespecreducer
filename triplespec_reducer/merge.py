"""
Joining the orders into one spectrum: XMERGEORDERS.

An anchor order is chosen, then each neighbouring order is added in turn:
scaled (1 unless the user sets a factor or asks for XMERGEORDERS_AUTOSCALE
over a region of the overlap), optionally trimmed at one end, and merged
with MC_MERGESPEC, which takes the error-weighted mean where the two
overlap and leaves a NaN pair at a gap between orders that do not.
"""

import numpy as np
from astropy.io import fits

from triplespec_reducer import idl as I
from triplespec_reducer.combine import interpspec, interpflagspec, combflagstack


def nantrim(array, flag):
    from triplespec_reducer.wavecal import nantrim as _nt
    return _nt(array, flag)


def mergespec(w1, f1, w2, f2, e1, e2, bf1, bf2, sum_=False):
    """MC_MERGESPEC (with errors and flags, as XMERGEORDERS calls it)."""
    # XMERGEORDERS passes DOUBLE(flux) and DOUBLE(error); the wavelengths
    # keep their own type.
    w1, w2 = np.asarray(w1), np.asarray(w2)
    f1, e1 = (np.asarray(a, dtype=float) for a in (f1, e1))
    f2, e2 = (np.asarray(a, dtype=float) for a in (f2, e2))
    bf1 = np.nan_to_num(np.asarray(bf1, dtype=float)).astype(np.uint8)
    bf2 = np.nan_to_num(np.asarray(bf2, dtype=float)).astype(np.uint8)
    i1 = nantrim(f1, 2)
    nw1, nf1, ne1, nbf1 = w1[i1], f1[i1], e1[i1], bf1[i1]
    n1 = nw1.size
    i2 = nantrim(f2, 2)
    nw2, nf2, ne2, nbf2 = w2[i2], f2[i2], e2[i2], bf2[i2]
    n2 = nw2.size
    min1, max1 = nw1.min(), nw1.max()
    min2, max2 = nw2.min(), nw2.max()
    pos = "Left" if min2 < min1 else ("Right" if max2 > max1 else "Inside")

    def combine_overlap(nnf2, nne2):
        if sum_:
            # TOTAL(..., 2) without /NAN: a NaN in either spectrum stays NaN
            return nf1 + nnf2.astype(float), I.sq(ne1) + I.sq(nne2.astype(float))
        stack = np.vstack([nf1, nnf2])
        dv = I.sq(np.vstack([ne1, nne2]))
        msk = (np.isfinite(stack) & np.isfinite(dv)).astype(int)   # MC_MEANCOMB /NAN
        f, v, _ = I.meancomb(stack, datavar=np.where(msk == 1, dv, 1.0), mask=msk)
        return f, v                    # FLOAT, as MC_MEANCOMB returns them

    if pos == "Right":
        if min2 > max1:
            return (np.concatenate([nw1, nw2]),
                    np.concatenate([nf1[:n1 - 1], [np.nan, np.nan], nf2[1:]]),
                    np.concatenate([ne1[:n1 - 1], [np.nan, np.nan], ne2[1:]]),
                    np.concatenate([nbf1[:n1 - 1], [0, 0], nbf2[1:]]).astype(np.uint8),
                    [np.nan, np.nan])
        z = np.where(nw1 < min2)[0]
        lw, lf, le, lb = nw1[z], nf1[z], ne1[z], nbf1[z]
        nnf2, nne2 = interpspec(nw2, nf2, nw1, iyerror=ne2, leavenans=True)
        idx = nantrim(nnf2, 1)
        mf, mv = combine_overlap(nnf2, nne2)
        nnbf2 = combflagstack([nbf1, interpflagspec(nw2, nbf2, nw1)])
        mw, mf, me, mb = nw1[idx], mf[idx], np.sqrt(mv[idx]), nnbf2[idx]
        z = np.where(nw2 >= mw.max())[0]
        rw, rf, re, rb = nw2[z[1:]], nf2[z[1:]], ne2[z[1:]], nbf2[z[1:]]
        return (np.concatenate([lw, mw, rw]), np.concatenate([lf, mf, rf]),
                np.concatenate([le, me, re]), np.concatenate([lb, mb, rb]).astype(np.uint8),
                [min2, mw.max()])
    if pos == "Left":
        if min1 > max2:
            return (np.concatenate([nw2, nw1]),
                    np.concatenate([nf2[:n2 - 1], [np.nan, np.nan], nf1[1:]]),
                    np.concatenate([ne2[:n2 - 1], [np.nan, np.nan], ne1[1:]]),
                    np.concatenate([nbf2[:n2 - 1], [0, 0], nbf1[1:]]).astype(np.uint8),
                    [np.nan, np.nan])
        z = np.where(nw2 <= min1)[0]
        lw, lf, le, lb = nw2[z[:z.size - 1]], nf2[z[:z.size - 1]], ne2[z[:z.size - 1]], nbf2[z[:z.size - 1]]
        nnf2, nne2 = interpspec(nw2, nf2, nw1, iyerror=ne2, leavenans=True)
        idx = nantrim(nnf2, 0)
        mf, mv = combine_overlap(nnf2, nne2)
        nnbf2 = combflagstack([nbf1, interpflagspec(nw2, nbf2, nw1)])
        mw, mf, me, mb = nw1[idx], mf[idx], np.sqrt(mv[idx]), nnbf2[idx]
        z = np.where(w1 > mw.max())[0]
        rw, rf, re, rb = nw1[z], nf1[z], ne1[z], nbf1[z]
        return (np.concatenate([lw, mw, rw]), np.concatenate([lf, mf, rf]),
                np.concatenate([le, me, re]), np.concatenate([lb, mb, rb]).astype(np.uint8),
                [min1, mw.max()])
    z = np.where(w1 <= min2)[0]
    lw, lf, le, lb = nw1[z[:z.size - 1]], nf1[z[:z.size - 1]], ne1[z[:z.size - 1]], nbf1[z[:z.size - 1]]
    z = np.where(nw1 > max2)[0]
    rw, rf, re, rb = nw1[z], nf1[z], ne1[z], nbf1[z]
    z = (nw2 > min1) & (nw2 < max1)
    nnf2, nne2 = interpspec(nw2[z], nf2[z], nw1, iyerror=ne2[z], leavenans=True)
    mf, mv = combine_overlap(nnf2, nne2)
    nnbf2 = combflagstack([nbf1, interpflagspec(nw2[z], nbf2[z], nw1)])
    idx = nantrim(mf, 2)
    mw, mf, me, mb = nw1[idx], mf[idx], np.sqrt(mv[idx]), nnbf2[idx]
    return (np.concatenate([lw, mw, rw]), np.concatenate([lf, mf, rf]),
            np.concatenate([le, me, re]), np.concatenate([lb, mb, rb]).astype(np.uint8),
            [min1, max2])


def autoscale(blue, green, wrange):
    """
    XMERGEORDERS_AUTOSCALE: the error-weighted scale of green to blue in
    wrange, with IDL's TOTALs; returned as the FLOAT scale field holds it.
    """
    wrange = np.sort(np.asarray(wrange, dtype=np.float32))
    zb = (blue[0] > wrange[0]) & (blue[0] < wrange[1])
    zg = (green[0] > wrange[0]) & (green[0] < wrange[1])
    bw, bf, be = blue[0][zb], blue[1][zb], blue[2][zb]
    gw, gf, ge = green[0][zg], green[1][zg], green[2][zg]
    ngf, nge = interpspec(gw, gf, bw, iyerror=ge)
    z = np.isfinite(bf) & np.isfinite(ngf) & np.isfinite(be) & np.isfinite(nge)
    denom = I.sq(be) + I.sq(nge)
    num, den = ngf[z] * bf[z] / denom[z], I.sq(ngf[z]) / denom[z]
    tot = lambda v: I.ftotal(v) if I.isf32(v) else I.itotal(v)
    return np.float32(tot(num) / tot(den))


def merge_orders(spec, orders, naps=1, aperture=1, anchor=None, sequence=None,
                 scales=None, trims=None, weighted=True, sum_=False, regions=None):
    """
    XMERGEORDERS from start to finish.

    `sequence` is the order in which orders are added (by default outward
    from the anchor, lower first, as the add-order list offers them);
    `scales[o]` multiplies order o before it is merged (default 1) and may
    be 'auto', scaled over regions[o] = (w0, w1); `trims[o]` is
    ('Green'|'Blue', 'Left'|'Right', wavelength): which spectrum is cut
    and which side of the wavelength is kept, as XMERGEORDERS_CUTSCALESPEC. Returns
    [4, npix] (wave, flux, error, flag) and the scales used.
    """
    # MC_READSPEC keeps the stored type (FLOAT for a corrected spectrum),
    # MC_CAASPECFLAGS adds the flag column.
    from triplespec_reducer.extract import caaspecflags
    spec = caaspecflags(np.asarray(spec)).copy()
    if not weighted:
        spec[:, 2, :] = 1.
    orders = list(orders)
    anchor = orders[0] if anchor is None else anchor
    ai = orders.index(anchor) * naps + (aperture - 1)
    nn = np.isfinite(spec[ai, 1])
    merged = spec[ai][:, nn]
    done = [anchor]
    used = {o: 1.0 for o in orders}
    if sequence is None:
        sequence = []
        lo = hi = orders.index(anchor)
        while len(sequence) < len(orders) - 1:
            if hi + 1 < len(orders):
                hi += 1
                sequence.append(orders[hi])
            if lo - 1 >= 0:
                lo -= 1
                sequence.append(orders[lo])
    scales = scales or {}
    trims = trims or {}
    for o in sequence:
        k = orders.index(o) * naps + (aperture - 1)
        nn = np.isfinite(spec[k, 1])
        add = spec[k][:, nn]
        blue = merged
        green = add.copy()
        t = trims.get(o)
        if t:
            which, side, w = t
            w = np.float32(w)                     # the FLOAT trim field
            if which == "Green":
                z = add[0] < w if side == "Right" else add[0] > w
                green = add[:, z].copy()
            else:
                z = (merged[0] < w) | ~np.isfinite(merged[0]) if side == "Right" else \
                    (merged[0] > w) | ~np.isfinite(merged[0])
                blue = merged[:, z]
        sc = scales.get(o, 1.0)
        if sc == "auto":
            sc = autoscale(blue, green, (regions or {})[o])
        sc = np.float32(sc)                       # the FLOAT scale field
        green[1] = green[1] * sc
        green[2] = green[2] * abs(sc)
        w, f, e, b, _ = mergespec(blue[0], blue[1], green[0], green[1], blue[2], green[2],
                                  blue[3], green[3], sum_=sum_)
        # [[owave],[oflux],[oerror],[obitflag]] with a DOUBLE flux: DOUBLE
        merged = np.array([w, f, e, b], dtype=np.float64)
        used[o] = float(sc)
        done.append(o)
    if not weighted:
        merged[2] = 1.
    return merged, used, done


def write(path, merged, hdr, used, done, aperture=1, source=""):
    """
    XMERGEORDERS_WRITEFILE: the [npix, 4] merged spectrum (unweighted: the
    errors already set to 1) under the input file's header, with ORDERS,
    NORDERS, NAPS, START and STOP, the output name in IRAFNAME, and the
    history.
    """
    import os
    h = hdr.copy()
    name = os.path.splitext(os.path.basename(path))[0]
    keep = lambda k: I.fxcomment(h, k) if k in h else ""      # no comment given: kept
    I.fxaddpar(h, "IRAFNAME", name, keep("IRAFNAME"))
    I.fxaddpar(h, "ORDERS", "0", keep("ORDERS"))
    I.fxaddpar(h, "NORDERS", 1, keep("NORDERS"))
    I.fxaddpar(h, "NAPS", 1, keep("NAPS"))
    I.fxaddpar(h, "START", 0, keep("START"))
    I.fxaddpar(h, "STOP", int(merged.shape[1] - 1), keep("STOP"))
    m = sorted(done)
    hist = (f"This merged spectrum was derived from aperture {aperture:02d} of file "
            f"{source.strip()}.  The scaled factors for orders "
            + ", ".join(str(o) for o in m) + " are "
            + ", ".join(f"{float(np.float32(used[o])):#13.6g}".strip() for o in m) + " respectively.")
    I.addhistlabel(h, "Xmergeorders History")
    I.sxaddhist(h, hist)
    fits.PrimaryHDU(merged, header=h).writeto(path, overwrite=True)
