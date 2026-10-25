"""
A spectrum file as text: SPEX2TEXT.

Every header card (as READFITS returns them, the END card included) on a
line of its own after '#', then one line per pixel with every column of
every order and aperture -- wavelength, flux, error, flag of the first
spectrum, then of the next -- each number as IDL's STRING prints it by
default (FLOAT in 13 characters with 6 significant digits, DOUBLE in 16
with 8) and joined by two spaces. The text file is named after the FITS
file, with '.txt' for '.fits'.
"""

import glob
import math
import os

import numpy as np
from astropy.io import fits


def _ieee(v, width, prec, inf_word):
    if math.isnan(v):
        s = "-NaN" if math.copysign(1.0, v) < 0 else "NaN"
    elif math.isinf(v):
        s = ("-" if v < 0 else "") + inf_word
    else:
        s = f"{v:#.{prec}g}"
        # IDL's choice: a value just below 10^prec that rounds up to it is
        # written '1.e+<prec>' (C's %#g writes 1.000...e+<prec>)
        if "e" in s and v != 0 and abs(v) < 10.0 ** prec and \
                math.floor(math.log10(abs(v))) == prec - 1:
            s = ("-" if v < 0 else "") + f"1.e+{prec:02d}"
    return s.rjust(width)


def idl_string(a):
    """IDL's default STRING of each element of `a`, as STRJOIN converts it."""
    a = np.asarray(a)
    k = a.dtype
    if k.kind == "f" and k.itemsize == 4:
        return [_ieee(float(v), 13, 6, "Inf") for v in a.ravel()]
    if k.kind == "f" and k.itemsize == 8:
        return [_ieee(float(v), 16, 8, "Infinity") for v in a.ravel()]
    if k.kind in "iu" and k.itemsize == 2:
        return [f"{int(v):8d}" for v in a.ravel()]
    if k.kind in "iu" and k.itemsize == 4:
        return [f"{int(v):12d}" for v in a.ravel()]
    if k.kind in "iu" and k.itemsize == 8:
        return [f"{int(v):22d}" for v in a.ravel()]
    raise TypeError(f"spex2text: no IDL default format for {k}")


def header_lines(path):
    """READFITS's header: the 80-character cards of the primary header up to END."""
    out = []
    with open(path, "rb") as fh:
        while True:
            block = fh.read(2880)
            if len(block) < 2880:
                break
            for i in range(36):
                card = block[i * 80:(i + 1) * 80].decode("ascii", "replace")
                out.append(card)
                if card.startswith("END") and card[3:].strip() == "":
                    return out
    return out


def text_lines(path):
    """The lines SPEX2TEXT writes for one file (without newlines)."""
    with fits.open(path) as h:
        d = h[0].data
        hdr = h[0].header
    d = d.astype(d.dtype.newbyteorder("="))
    npix = int(hdr["NAXIS1"])
    norders, naps = int(hdr["NORDERS"]), int(hdr["NAPS"])
    cube = d[None] if d.ndim == 2 else d                   # [nspec, ncol, npix]
    narr = cube.shape[1]
    n = narr * naps * norders
    lines = ["#" + c for c in header_lines(path)]
    for j in range(npix):
        row = cube[:, :, j].reshape(-1)
        if row.size != n:
            raise ValueError(f"{path}: {row.size} values per pixel, the header says {n}")
        lines.append("  ".join(idl_string(row)))
    return lines


def spex2text(pattern, outdir=None):
    """
    SPEX2TEXT on every file that matches `pattern` (FILE_SEARCH, a single
    name or a wildcard). Each goes to <basename>.txt in `outdir`: SPEX2TEXT
    writes into IDL's current directory; here the default is the FITS
    file's own directory. Returns the paths written.
    """
    files = sorted(glob.glob(os.path.expanduser(pattern.strip())))
    written = []
    for f in files:
        root = os.path.basename(f)
        if root.endswith(".fits"):
            root = root[:-5]
        out = os.path.join(outdir if outdir is not None else os.path.dirname(os.path.abspath(f)),
                           root + ".txt")
        with open(out, "w") as fh:
            for line in text_lines(f):
                fh.write(line + "\n")
        written.append(out)
    return written
