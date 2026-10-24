"""
The reduction as a chain: each Python stage reads only the files the Python
stage before it wrote (flat -> wavelength calibration -> extraction ->
combination), and every file is compared with IDL's file from IDL's own chain
(validation/idl/ref_flat.pro, ref_wavecal.pro, ref_dg.pro): the data bit for
bit and the header card for card.

The stage-by-stage comparisons feed each stage IDL's file, so they cannot see
what a file's header carries to the next stage: Spextool writes numbers as
text (FXADDPAR) and reads that text back (FXPAR), rounded and typed by it.
This script follows the numbers through that round trip.
"""

import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from triplespec_reducer import (flat as F, wavecal as W, extract as E,   # noqa: E402
                                combine as C)
from triplespec_reducer.instrument import read_atran                      # noqa: E402

RAW = os.path.join(os.path.dirname(HERE), "data", "raw", "2025-11-02")
REF = os.path.join(HERE, "idl_out")
OUT = os.path.join(HERE, "py_chain")
os.makedirs(OUT, exist_ok=True)

# Cards that cannot agree: the software version, the date a file was written,
# and the FITS structure cards astropy writes itself.
SKIP = {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "EXTEND", "COMMENT",
        "HISTORY", "", "DATE", "VERSION",
        # MC_ADDHISTLABEL's banner starts in column 2, not a valid FITS record;
        # written here as a blank-keyword card with the same text
        "#######"}
# Ours, not IDL's: the wavecal's record of the arc type, offset and 1D rms.
OURS = {"TSRTYPE", "TSROFFS", "TSRRMS1D"}
failures = []


def rd(path, ext=0):
    with fits.open(path) as h:
        d = h[ext].data
        return d.astype(d.dtype.newbyteorder("="))


def same_data(name, a, b):
    ok = a.shape == b.shape and a.dtype == b.dtype and \
        bool(np.all((a == b) | (np.isnan(a.astype(float)) & np.isnan(b.astype(float)))))
    print(f"  {name:40s} data {'identical' if ok else 'DIFFERENT'}"
          + ("" if ok else f"  {a.shape}{a.dtype} vs IDL {b.shape}{b.dtype}"))
    if not ok:
        failures.append(name + " data")


def same_header(name, ours, idl):
    a, b = fits.getheader(ours), fits.getheader(idl)
    diff = [k for k in b if k not in SKIP and (k not in a or a.cards[k].image != b.cards[k].image)]
    extra = [k for k in a if k not in b and k not in SKIP | OURS]
    print(f"  {name:40s} header {'identical' if not diff and not extra else 'DIFFERENT'}")
    for k in diff[:8]:
        print(f"      ours {a.cards[k].image.rstrip() if k in a else '(missing)'}\n"
              f"      IDL  {b.cards[k].image.rstrip()}")
    if extra:
        print("      only ours:", extra)
    if diff or extra:
        failures.append(name + " header")


print("flat (ref_flat.pro):")
on = [os.path.join(RAW, f"SPEC_DFLAT_{i:04d}.fits") for i in range(2, 52)]
off = [os.path.join(RAW, f"SPEC_DFLAT_{i:04d}.fits") for i in range(52, 104)]
F.write_flat(os.path.join(OUT, "flat.fits"), F.make_flat(on, off), version="v4.1")
for e, what in ((1, "flat"), (2, "variance"), (3, "flags")):
    same_data(f"flat.fits [{what}]", rd(os.path.join(OUT, "flat.fits"), e), rd(os.path.join(REF, "flat.fits"), e))
same_header("flat.fits", os.path.join(OUT, "flat.fits"), os.path.join(REF, "flat.fits"))
fl = F.read_flat(os.path.join(OUT, "flat.fits"))

print("wavelength calibration on sky (ref_wavecal.pro), from our flat file:")
sky = [os.path.join(RAW, f"SPEC_SO2025D-006-{n:04d}.fits") for n in range(120, 132)]
W.write_wavecal(os.path.join(OUT, "wavecal_NGC134.fits"), W.make_wavecal(sky, fl))
for e, what in ((1, "wavelengths"), (2, "slit positions")):
    same_data(f"wavecal_NGC134.fits [{what}]", rd(os.path.join(OUT, "wavecal_NGC134.fits"), e),
              rd(os.path.join(REF, "wavecal_NGC134.fits"), e))
same_header("wavecal_NGC134.fits", os.path.join(OUT, "wavecal_NGC134.fits"),
            os.path.join(REF, "wavecal_NGC134.fits"))
wc = W.read_wavecal(os.path.join(OUT, "wavecal_NGC134.fits"), rotation=fl["rotation"])

print("extraction (ref_dg.pro), from our flat and wavecal files:")
atran = read_atran(fl["rp"])
for (a, b), kind in (((116, 117), "point"), ((118, 119), "point"),
                     ((120, 121), "extended"), ((126, 127), "extended")):
    files = [os.path.join(RAW, f"SPEC_SO2025D-006-{n:04d}.fits") for n in (a, b)]
    op = [os.path.join(OUT, f"spectra{a:04d}.fits"), os.path.join(OUT, f"spectra{b:04d}.fits")]
    r = E.reduce_pair(files, fl, wc, op, source=kind, atran=atran, xs_positions=[20.1],
                      flatname="flat.fits", wavecalname="wavecal_NGC134.fits")
    for p in r["written"]:
        n = os.path.basename(p)
        same_data(n, rd(p), rd(os.path.join(REF, n)))
        same_header(n, p, os.path.join(REF, n))

print("combination (ref_dg.pro), from our spectra files:")
for names, oname in ((["spectra0116", "spectra0117", "spectra0118", "spectra0119"], "comb_HD225200"),
                     (["spectra0120", "spectra0126"], "comb_NGC134")):
    L = C.load([os.path.join(OUT, f + ".fits") for f in names])
    C.scale(L, 5)
    C.combine(L)
    p = os.path.join(OUT, oname + ".fits")
    C.write(p, L)
    # ref_dg.pro fills wavelength, flux and error; its flag column stays 0
    same_data(oname + ".fits [w, f, e]", rd(p)[:, :3], rd(os.path.join(REF, oname + ".fits"))[:, :3])
    same_header(oname + ".fits", p, os.path.join(REF, oname + ".fits"))

print("\nall identical" if not failures else f"\n{len(failures)} differences: {failures}")
