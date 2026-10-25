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


# ------------------------------------------------------------- finding files

def find_raw_frames(datadir, suffix="*.fits*"):
    """
    Every raw frame in `datadir` the instrument file's SUFFIX lets
    xspextool see (TS4_SOAR: '*.fits*', so .fits, .fits.gz and .fits.fz
    alike), sorted as FILE_SEARCH returns them.
    """
    import glob
    import os
    return sorted(glob.glob(os.path.join(datadir, "*" + suffix.lstrip("*"))))


def mkfullpath(datapath, files, index=True, nint=4, prefix="", suffix="*.fits*",
               exist=True):
    """
    MC_MKFULLPATH: the paths of frames typed in xspextool's fields.

    With `index` (File Read Mode 'Index'), `files` is an index list as
    MC_FSEXTRACT reads it ('120-131', '7-3', '1,4-6', or a list of ints),
    each becoming <datapath><prefix><index, NINT digits><suffix>; with
    `exist`, the first file FILE_SEARCH finds for that pattern, an error if
    none. Otherwise (File Read Mode 'Filename') `files` is a comma list or a
    list of names inside `datapath`.
    """
    import glob
    import os
    from triplespec_reducer.instrument import fsextract_index, fsextract_filename
    if index:
        idx = fsextract_index(files) if isinstance(files, str) else [int(f) for f in files]
        if any(i > 10 ** nint - 1 for i in idx):
            raise ValueError(f"Image numbers greater than {10 ** nint - 1} not allowed.")
        out = []
        for i in idx:
            name = "".join(os.path.join(datapath, f"{prefix}{i:0{nint}d}{suffix}").split())
            if exist:
                found = sorted(glob.glob(name))
                if not found:
                    raise FileNotFoundError(f"File {name} does not exist.")
                name = found[0]
            out.append(name)
        return out
    names = fsextract_filename(files) if isinstance(files, str) else list(files)
    out = [os.path.join(datapath, n) for n in names]
    if exist:
        for p in out:
            if not glob.glob(p):
                raise FileNotFoundError(f"File {p} does not exist.")
    return out


# ------------------------------------------------------------- processed images

def _native(a):
    return a.astype(a.dtype.newbyteorder("="))


def read_processed(path, rotation=3):
    """
    How XSPEXTOOL_LOADIMAGES reads an image that is already processed (a
    combined image from cal/ or proc/, the Source Images directory not
    raw/): ext 1-3 data, variance and flags, rotated, kept in their stored
    types (MRDFITS), the header of ext 0 and ITOT from it (FXPAR: 0 when
    absent). Returns (data, var, flags, header, itot).
    """
    from triplespec_reducer.idl import fxpar
    with fits.open(path, memmap=False) as h:
        hdr = h[0].header
        data = rotate(_native(h[1].data), rotation)
        var = rotate(_native(h[2].data), rotation)
        flags = rotate(_native(h[3].data), rotation)
    itot = fxpar(hdr, "ITOT") if "ITOT" in hdr else 0
    return data, var, flags, hdr, itot


def is_processed(path):
    """
    The test XSPEXTOOL_LOADIMAGES makes in A-Sky mode: the first file's
    primary header has a non-zero ITOT, so it is a processed image rather
    than a raw frame.
    """
    from triplespec_reducer.idl import fxpar
    hdr = fits.getheader(path, 0)
    return "ITOT" in hdr and fxpar(hdr, "ITOT") != 0


def load_source_image(files, mode, rotation=3, sourcedir="raw/"):
    """
    The object image of XSPEXTOOL_LOADIMAGES before the sky, flat and
    error-propagation steps, for every source directory the IDL allows:

    A       raw/: KNA_READTS4FITS of the frame; cal/ or proc/: the processed
            image (read_processed)
    A-B     always raw, as a pair
    A-Sky   the first file is read as a processed image when its header's
            ITOT is non-zero (is_processed), otherwise as a raw frame

    Returns (data, var, headers, itot) with data and var 2D.
    """
    files = [str(f) for f in files]
    if mode == "A-B":
        d, v, h, t = read_ts4(files, rotation=rotation, pair=True)
        return d[0], v[0], h, t
    processed = (sourcedir != "raw/") if mode == "A" else is_processed(files[0])
    if mode not in ("A", "A-Sky"):
        raise ValueError(f"unknown reduction mode {mode!r}")
    if processed:
        d, v, _, h, t = read_processed(files[0], rotation)
        return d, v, [h], [t]
    d, v, h, t = read_ts4(files[:1], rotation=rotation)
    return d[0], v[0], h, t
