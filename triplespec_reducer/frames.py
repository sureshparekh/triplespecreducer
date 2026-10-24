"""
Reading raw TripleSpec 4.1 frames: KNA_READTS4FITS.

A frame becomes a count rate (DN/s) with a variance, rotated so that the
dispersion runs along columns and order 3 is at the bottom. IDL's
ROTATE(img, 3) on the [x, y] array is numpy's rot90(img, 1) on [row, col];
MC_UNROTATE is the inverse. With PAIR, consecutive frames are differenced:
A/t_A - B/t_B, with the two variances added.
"""

import numpy as np
from astropy.io import fits



def _sq(x):
    """x^2 as IDL computes it, x*x (numpy's scalar ** 2 is not always that)."""
    return x * x

def rotate(img, rotation):
    """IDL ROTATE for the rotations Spextool uses, on a numpy image."""
    if rotation == 0:
        return img
    if rotation == 3:
        return np.rot90(img, 1)
    if rotation == 1:
        return np.rot90(img, -1)
    if rotation == 2:
        return np.rot90(img, 2)
    raise ValueError(f"rotation {rotation} is not used by TS4_SOAR")


def unrotate(img, rotation):
    """MC_UNROTATE."""
    return rotate(img, {0: 0, 1: 3, 2: 2, 3: 1}[rotation])


def _read_one(path):
    with fits.open(path, memmap=False) as hdul:
        hdu = hdul[1] if path.endswith(".fz") else hdul[0]
        return hdu.data, hdu.header


def _rate_and_var(img, hdr):
    """Count rate and its variance from one frame's own header."""
    f32 = np.float32
    itime = f32(hdr["EXPTIME"])
    coadds = f32(hdr["COADDS"])
    ndrs = f32(hdr["FSAMPLE"])
    gain = f32(hdr["GAIN"])
    readnoise = f32(hdr["RON"])
    readtime = f32(hdr["FRAMTIME"])
    rdvar = (f32(2.) * _sq(readnoise)) / ndrs / coadds / _sq(itime) / _sq(gain)
    crtn = (f32(1.0) - readtime * (_sq(ndrs) - f32(1.0)) / f32(3.) / itime / ndrs)
    img = img.astype(np.float32)
    var = np.abs(img) * crtn / (_sq(coadds)) / (_sq(itime)) / gain + rdvar
    rate = img / coadds / itime
    return rate.astype(np.float32), var.astype(np.float32), float(itime * coadds)


def read_ts4(files, rotation=3, pair=False):
    """
    KNA_READTS4FITS.

    Returns (data, var, headers, itot): data and var are stacks [n, rows,
    cols] in the rotated orientation (n = len(files), or half that with
    `pair`), headers the FITS headers and itot the integration times.
    """
    files = list(files)
    if pair and len(files) % 2:
        raise ValueError("A-B needs an even number of frames.")
    hdrs, itot = [], []
    data, var = [], []
    if pair:
        for i in range(len(files) // 2):
            a, ha = _read_one(files[2 * i])
            b, hb = _read_one(files[2 * i + 1])
            ra, va, ta = _rate_and_var(a, ha)
            rb, vb, tb = _rate_and_var(b, hb)
            data.append(rotate(ra - rb, rotation))
            var.append(rotate(va + vb, rotation))
            hdrs += [ha, hb]
            itot += [ta, tb]
    else:
        for f in files:
            img, h = _read_one(f)
            r, v, t = _rate_and_var(img, h)
            data.append(rotate(r, rotation))
            var.append(rotate(v, rotation))
            hdrs.append(h)
            itot.append(t)
    return (np.ascontiguousarray(np.array(data, dtype=np.float32)),
            np.ascontiguousarray(np.array(var, dtype=np.float32)), hdrs, itot)
