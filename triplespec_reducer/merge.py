"""
Joining the orders into one spectrum: XMERGEORDERS.

An anchor order is chosen, then each neighbouring order (addable_orders)
is added in turn: scaled (1 unless the user sets a factor or asks for
XMERGEORDERS_AUTOSCALE over a region of the overlap), optionally cut (the
order being added, the merged spectrum, or both), and merged with
MC_MERGESPEC, which takes the error-weighted mean where the two overlap
and leaves a NaN pair at a gap between orders that do not.
"""

import numpy as np
from astropy.io import fits

from triplespec_reducer import idl as I
from triplespec_reducer.combine import interpspec, interpflagspec, combflagstack


def nantrim(array, flag):
    from triplespec_reducer.wavecal import nantrim as _nt
    return _nt(array, flag)


def mergespec(w1, f1, w2, f2, e1, e2, bf1, bf2, sum_=False):
    """
    MC_MERGESPEC (with errors and flags, as XMERGEORDERS calls it). `sum_`
    is its /SUM, which XMERGEORDERS cannot reach (see merge_orders).
    """
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
    # MC_MERGESPEC passes IYERROR=ne2, the whole error array, with the
    # x and y of the overlap only (nw2[z], nf2[z]); MC_INTERPSPEC takes the
    # errors at the same subscripts as x (iyerror[zci]), so the errors used
    # are the first ones of ne2, not those of the overlap.
    nnf2, nne2 = interpspec(nw2[z], nf2[z], nw1, iyerror=ne2[:int(z.sum())], leavenans=True)
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

    XMERGEORDERS_AUTOSCALE first sets the scale to 1 and calls
    XMERGEORDERS_CUTSCALESPEC with no trim, so `blue` must be the whole
    merged spectrum and `green` the whole (unscaled) order being added;
    merge_orders passes them so.
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


def addable_orders(orders, merged):
    """
    XMERGEORDERS_MODADDORDER: the orders the 'Add Order' list offers next,
    the one below and the one above the merged block in the file's order
    list ([min-1, max+1] of their positions, kept when inside the list).
    Empty when every order is merged ('Merging complete').
    """
    orders = list(orders)
    if len(set(merged)) >= len(orders):
        return []
    pos = [orders.index(o) for o in merged]
    possible = [min(pos) - 1, max(pos) + 1]
    return [orders[k] for k in possible if 0 <= k <= len(orders) - 1]


def default_sequence(orders, anchor):
    """
    The add order this library uses when none is given: outward from the
    anchor, the next position up first, then the next one down, so every
    order added is one addable_orders offers.
    """
    orders = list(orders)
    seq = []
    lo = hi = orders.index(anchor)
    while len(seq) < len(orders) - 1:
        if hi + 1 < len(orders):
            hi += 1
            seq.append(orders[hi])
        if lo - 1 >= 0:
            lo -= 1
            seq.append(orders[lo])
    return seq


def overlap_region(blue, green):
    """
    A region for XMERGEORDERS_AUTOSCALE: the whole wavelength overlap of the
    merged spectrum and the order being added (FLOAT, as the cursor gives
    it). XMERGEORDERS has no default region -- the user selects one with
    's' -- so this is only a choice a user could make there by hand; the
    pipeline uses it when it is asked for (regions[o] = 'overlap').
    """
    lo = max(np.nanmin(blue[0]), np.nanmin(green[0]))
    hi = min(np.nanmax(blue[0]), np.nanmax(green[0]))
    if not lo < hi:
        raise ValueError("the order does not overlap the merged spectrum: no region to scale over")
    return np.float32(lo), np.float32(hi)


TRIM_SPECTRA = ("Green", "Blue")
TRIM_SIDES = ("Right", "Left")


def _trim_list(t):
    """trims[o] as a list of (spectrum, side, wavelength) cuts."""
    if not t:
        return []
    if isinstance(t[0], str):          # a single (spectrum, side, wavelength)
        t = [t]
    out = []
    for which, side, w in t:
        which, side = str(which).capitalize(), str(side).capitalize()
        if which not in TRIM_SPECTRA or side not in TRIM_SIDES:
            raise ValueError(f"trim {which!r}/{side!r}: the spectrum is Green or Blue, the side Right or Left")
        out.append((which, side, w))
    return out


def merge_orders(spec, orders, naps=1, aperture=1, anchor=None, sequence=None,
                 scales=None, trims=None, regions=None, check_adjacent=True,
                 nonidl_unweighted=False, nonidl_sum_flux=False,
                 weighted=None, sum_=None):
    """
    XMERGEORDERS from start to finish, for one aperture.

    `sequence` is the order in which orders are added (by default
    default_sequence: outward from the anchor); with `check_adjacent` each
    one must be an order the 'Add Order' list would offer at that point
    (addable_orders), as XMERGEORDERS allows no other.

    `scales[o]` multiplies order o before it is merged (default 1, the
    'Scale' field) and may be 'auto': XMERGEORDERS_AUTOSCALE over
    `regions[o]`, a (w0, w1) region (the 's' selection) or 'overlap'
    (overlap_region). The scale is found, as XMERGEORDERS finds it, on the
    untrimmed merged spectrum and the untrimmed order. 'auto' without a
    region is refused, as XMERGEORDERS refuses it ("Please select a region
    using 's'.").

    `trims[o]` is one ('Green'|'Blue', 'Right'|'Left', wavelength) cut or a
    list of them, XMERGEORDERS_CUTSCALESPEC: 'Green' cuts the order being
    added, 'Blue' the merged spectrum, 'Right' keeps what lies below the
    wavelength, 'Left' what lies above. A green cut and a blue cut are kept
    together; a later cut of the same spectrum replaces the earlier one (it
    is made from the uncut spectrum again). They act after the scale, as in
    XMERGEORDERS, where a new scale (typed or automatic) cancels the cuts.

    Two choices do not exist in XMERGEORDERS 4.1 and are off by default:
    `nonidl_unweighted` (its 'Weighted Mean' unticked: errors set to 1,
    plain mean) and `nonidl_sum_flux` (its 'Sum Flux': fluxes added, the
    variances added). In XMERGEORDERS both stop with an undefined variable
    (XMERGEORDERS_COMBINESPEC uses BSPEC_E, which exists only when weighted),
    so they produce no merged spectrum there; here they give results with
    no Spextool counterpart. `weighted` and `sum_` are their older names
    (weighted=False is nonidl_unweighted=True, sum_=True is
    nonidl_sum_flux=True).

    Returns [4, npix] (wave, flux, error, flag), the scales used and the
    merged orders in the order they were added.
    """
    if weighted is not None:
        nonidl_unweighted = nonidl_unweighted or not weighted
    if sum_ is not None:
        nonidl_sum_flux = nonidl_sum_flux or bool(sum_)
    unweighted = nonidl_unweighted or nonidl_sum_flux      # 'Sum Flux' unticks 'Weighted Mean'
    # MC_READSPEC keeps the stored type (FLOAT for a corrected spectrum),
    # MC_CAASPECFLAGS adds the flag column.
    from triplespec_reducer.extract import caaspecflags
    spec = caaspecflags(np.asarray(spec)).copy()
    if unweighted:
        spec[:, 2, :] = 1.
    orders = list(orders)
    anchor = orders[0] if anchor is None else anchor
    ai = orders.index(anchor) * naps + (aperture - 1)
    nn = np.isfinite(spec[ai, 1])
    merged = spec[ai][:, nn]
    done = [anchor]
    used = {o: 1.0 for o in orders}
    if sequence is None:
        sequence = default_sequence(orders, anchor)
    scales = scales or {}
    trims = trims or {}
    regions = regions or {}
    for o in sequence:
        if check_adjacent and o not in addable_orders(orders, done):
            raise ValueError(f"order {o} cannot be added after {done}: XMERGEORDERS offers only "
                             f"{addable_orders(orders, done)}")
        k = orders.index(o) * naps + (aperture - 1)
        nn = np.isfinite(spec[k, 1])
        add = spec[k][:, nn]
        sc = scales.get(o, 1.0)
        if isinstance(sc, str):
            if sc != "auto":
                raise ValueError(f"scale {sc!r} for order {o}: a number or 'auto'")
            reg = regions.get(o)
            if reg is None:
                raise ValueError(f"auto scale for order {o} needs a region (regions[{o}] = (w0, w1) "
                                 "or 'overlap'): XMERGEORDERS asks 'Please select a region using 's'.'")
            if isinstance(reg, str):
                if reg != "overlap":
                    raise ValueError(f"region {reg!r} for order {o}: (w0, w1) or 'overlap'")
                reg = overlap_region(merged, add)
            # XMERGEORDERS_AUTOSCALE: scale 1, no trim, so the whole merged
            # spectrum and the whole order
            sc = autoscale(merged, add.copy(), reg)
        sc = np.float32(sc)                       # the FLOAT scale field
        blue = merged
        green = add.copy()
        for which, side, w in _trim_list(trims.get(o)):
            w = np.float32(w)                     # the FLOAT trim field
            if which == "Green":
                z = add[0] < w if side == "Right" else add[0] > w
                green = add[:, z].copy()
            else:
                z = (merged[0] < w) | ~np.isfinite(merged[0]) if side == "Right" else \
                    (merged[0] > w) | ~np.isfinite(merged[0])
                blue = merged[:, z]
        green[1] = green[1] * sc
        green[2] = green[2] * abs(sc)
        w, f, e, b, _ = mergespec(blue[0], blue[1], green[0], green[1], blue[2], green[2],
                                  blue[3], green[3], sum_=nonidl_sum_flux)
        # [[owave],[oflux],[oerror],[obitflag]] with a DOUBLE flux: DOUBLE
        merged = np.array([w, f, e, b], dtype=np.float64)
        used[o] = float(sc)
        done.append(o)
    if unweighted:
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
