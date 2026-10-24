"""
Combining raw images before extraction: XSPEXTOOL_COMBINEIMGS.

Many frames (mode A) or A-B pairs (mode A-B) are averaged into one image,
written to proc/ (or cal/, for a sky frame used in A-Sky mode) as a header
plus three images -- flux, variance, flags -- in the raw orientation, the
file xspextool then reads in A mode.

Two of its options do not work in the IDL and are reproduced as they are:

* 'Scale orders' (MC_SCALESPECSKY) reuses the variable that holds an
  order's pixels for the list of non-finite scale factors; with none, that
  list is IDL's -1 and the subscript stops the run. Asked for, it stops
  here too, with the reason; `scale_orders="intended"` applies the
  per-order scaling the routine was written to make.
* 'Subtract the sky' (MC_SUBSPECSKY) leaves the rest of an order alone
  once a column falls off the array (goto, cont1) and tests rows against
  the column count. Reproduced.

The 'Weighted Mean' statistic passes a variable that does not exist
(DATAVAR=ivcube), so it is a plain mean; reproduced.

One defect is not reproduced, because it leaves nothing to reproduce.
Without 'scale orders' or 'subtract the sky' the IDL reads the frames with
ROTATE=0, and KNA_READTS4FITS stores each [1024, 2048] frame into an array
it allocated as [2048, 1024]: the values are the right ones, in the raw
memory order, under the wrong shape (checked in IDL 7.1). The file written
has every combined pixel value IDL computes but cannot be extracted -- the
A-mode read then mismatches the flat. Here the same values are written in
the correct orientation.
"""

import os

import numpy as np
from astropy.io import fits

from triplespec_reducer import idl as I
from triplespec_reducer.frames import read_ts4, rotate, unrotate

COMBSTATS = ["Robust Weighted Mean", "Robust Mean (RMS)", "Robust Mean (Std Error)",
             "Weighted Mean", "Mean (RMS)", "Mean (Std Error)", "Median (MAD)",
             "Median (Median Error)"]


def abba_to_abab(files):
    """The 'ABBA to ABAB' switch: every second pair swapped, so all are A-B."""
    files = list(files)
    idx = list(range(len(files)))
    for i in range(len(files) // 2):
        if i % 2:
            idx[i * 2], idx[i * 2 + 1] = i * 2 + 1, i * 2
    return [files[k] for k in idx]


def scalespecsky(data, var, omask, orders, intended=False):
    """MC_SCALESPECSKY; see the module notes."""
    for o in orders:
        z = omask == o
        stat = np.array([I.median(d[z], even=True) for d in data])
        med = I.median(stat, even=True)
        scale = med / stat
        bad = ~np.isfinite(scale)
        if not intended:
            if not bad.any():
                raise ValueError("Scale orders stops here in the IDL tool: MC_SCALESPECSKY "
                                 "overwrites the order's pixel list and subscripts with -1. "
                                 "Choose 'intended' to apply the scaling it was meant to make.")
            # IDL's reused Z now holds the positions of the bad factors in
            # the frame list, and those PIXELS of each frame are scaled.
            scale[bad] = 1.0
            zi = np.where(bad)[0]
            for j in range(data.shape[0]):
                fd = data[j].reshape(-1)
                fv = var[j].reshape(-1)
                fd[zi] = fd[zi] * scale[j]
                fv[zi] = fv[zi] * I.sq(scale[j])
            continue
        scale[bad] = 1.0
        for j in range(data.shape[0]):
            data[j][z] = data[j][z] * scale[j]
            var[j][z] = var[j][z] * I.sq(scale[j])
    return data, var


def subspecsky(images, var, edgecoeffs, xranges, ybuffer=3):
    """MC_SUBSPECSKY (no DEG): a robust weighted mean of each column, removed."""
    images = images.astype(np.float32)
    nframes, nrows, ncols = images.shape
    for i in range(nframes):
        for j in range(len(xranges)):
            start, stop = int(xranges[j][0]), int(xranges[j][1])
            x = np.arange(stop - start + 1) + start
            bote = I.iround(I.poly(x.astype(float), edgecoeffs[j, 0]))
            tope = I.iround(I.poly(x.astype(float), edgecoeffs[j, 1]))
            for k in range(x.size):
                if bote[k] > ncols - 1 or tope[k] < 0:
                    break                       # IDL: goto, cont1 -- the next order
                bot = max(0, bote[k]) + ybuffer
                top = min(ncols - 1, tope[k]) - ybuffer
                if top - bot < 2:
                    break
                # FLOAT throughout, with IDL's running TOTALs
                f32 = np.float32
                sf = images[i, bot:top + 1, x[k]].astype(f32)
                sv = var[i, bot:top + 1, x[k]].astype(f32)
                z = I.findoutliers(sf, 3) == 1
                if not z.any():
                    bg = f32(I.ftotal(sf / sv) / I.ftotal(f32(1.) / sv))
                    bgv = f32(f32(1.) / I.ftotal(f32(1.) / sv))
                else:
                    bg = f32(I.ftotal(sf[z] * f32(1.) / sv[z]) / I.ftotal(f32(1.) / sv[z]))
                    bgv = f32(f32(1.) / I.ftotal(f32(1.) / sv[z]))
                bot, top = bot - ybuffer, top + ybuffer
                images[i, bot:top + 1, x[k]] -= np.float32(bg)
                var[i, bot:top + 1, x[k]] += np.float32(bgv)
    return images, var


def combine_images(files, out_path, mode="A", statistic="Robust Mean (Std Error)",
                   thresh=8.0, abba2abab=False, scale_orders=False, subtract_sky=False,
                   flatinfo=None, lincormax=80000.0, progress=None):
    """
    XSPEXTOOL_COMBINEIMGS. `scale_orders` is False, True (as the IDL: stops)
    or 'intended'. Returns the path written.
    """
    say = progress or (lambda *_: None)
    pair = mode == "A-B"
    files = list(files)
    if pair and len(files) % 2:
        raise ValueError("Must be an even number of images for pair subtraction.")
    if pair and abba2abab:
        files = abba_to_abab(files)
    rotation = flatinfo["rotation"] if (scale_orders or subtract_sky) else 0
    say("Loading the images")
    data, var, hdrs, itot = read_ts4(files, rotation=rotation, pair=pair)
    # The linearity flag KNA_READTS4FITS sets: a raw value above LINCORMAX.
    flags = []
    raw = [fits.getdata(f).astype(float) for f in files]
    if pair:
        for k in range(len(files) // 2):
            flags.append(rotate(((raw[2 * k] > lincormax) | (raw[2 * k + 1] > lincormax)).astype(np.uint8), rotation))
    else:
        flags = [rotate((r > lincormax).astype(np.uint8), rotation) for r in raw]
    flags = np.array(flags)
    itot = np.float32(I.ftotal(np.float32(itot)) / np.float32(2.)) if pair else I.ftotal(np.float32(itot))
    if scale_orders:
        say("Scaling the orders")
        data, var = scalespecsky(data, var, flatinfo["omask"], flatinfo["orders"],
                                 intended=(scale_orders == "intended"))
    if subtract_sky:
        say("Subtracting the sky level")
        data, var = subspecsky(data, var, flatinfo["edgecoeffs"], flatinfo["xranges"], 3)
    say("Combining")
    stat = COMBSTATS.index(statistic)
    ogb = None
    if stat == 0:
        mean, mvar, ogb = I.meancomb(data, datavar=var, robust=thresh)
    elif stat == 1:
        mean, mvar, ogb = I.meancomb(data, rms=True, robust=thresh)
    elif stat == 2:
        mean, mvar, ogb = I.meancomb(data, robust=thresh)
    elif stat == 3:
        mean, mvar, _ = I.meancomb(data)          # DATAVAR=ivcube: undefined in the IDL
    elif stat == 4:
        mean, mvar, _ = I.meancomb(data, rms=True)
    elif stat == 5:
        mean, mvar, _ = I.meancomb(data)
    elif stat == 6:
        mean, mvar = I.medcomb(data, mad=True)
    else:
        mean, mvar = I.medcomb(data)
    if ogb is not None:
        flags = flags * ogb.astype(np.uint8)
    bitmask = np.bitwise_or.reduce(flags, axis=0)
    # MC_AVEHDRS over every frame's instrument keywords, written by FXADDPAR
    from triplespec_reducer.combine import avehdrs, gethdrinfo
    from triplespec_reducer.instrument import read_instrfile
    keys = read_instrfile().keywords
    vals, coms = avehdrs([gethdrinfo(h, keys) for h in hdrs])
    vals["FILENAME"] = os.path.basename(out_path)
    hdr = fits.Header()
    for k, v in vals.items():
        I.fxaddpar(hdr, k, v, coms[k])
    I.fxaddpar(hdr, "NIMGS", len(files), " Number of images combined.")
    I.fxaddpar(hdr, "PAIRSUB", int(pair), " Pair subtraction?: 1=yes, 0=no.")
    I.fxaddpar(hdr, "ITOT", np.float32(itot), " Total integration time PER BEAM")
    I.fxaddpar(hdr, "LINMAX", int(lincormax), " Linearity maximum in DN")       # LONG
    hist = (f"This image was created by {statistic} combining the images "
            + ", ".join(os.path.basename(f) for f in files) + ".  ")
    if scale_orders:
        hist += ("Each order within an image was first scaled to the median flux value of all "
                 "of the pixels in that order in all of the images.  ")
    if subtract_sky:
        hist += "This median flux level was subtracted off on a column by column basis within each order."
    hdr.add_history("Xspextool History: " + hist)
    fits.HDUList([fits.PrimaryHDU(header=hdr),
                  fits.ImageHDU(unrotate(np.asarray(mean, np.float32), rotation)),
                  fits.ImageHDU(unrotate(np.asarray(mvar, np.float32), rotation)),
                  fits.ImageHDU(unrotate(bitmask.astype(np.uint8), rotation))]).writeto(out_path, overwrite=True)
    return out_path


def read_combined(path, rotation=3):
    """How XSPEXTOOL_LOADIMAGES reads a combined image (A mode, from proc/)."""
    with fits.open(path) as h:
        hdr = h[0].header
        data = rotate(h[1].data.astype(np.float32), rotation)
        var = rotate(h[2].data.astype(np.float32), rotation)
        flags = rotate(h[3].data.astype(np.uint8), rotation)
    return data, var, flags, hdr, I.fxpar(hdr, "ITOT") if "ITOT" in hdr else 0
