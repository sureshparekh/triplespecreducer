"""
Calibration branches the main parity runs do not reach, against IDL
(validation/idl/ref_flat_branches.pro, ref_wavecal_branches.pro): a single
lamps-on frame, a single dark, no darks; a single sky frame and a manual
offset. Also the calibration QA (the printed lists of rejected lines and
the RMS titles of the IDL plots), the instrument file and file handling.
"""

import os
import re
import sys

import numpy as np
import pytest
from astropy.io import fits

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
REF = os.path.join(ROOT, "validation", "idl_out")
RAW = os.path.join(ROOT, "data", "raw", "2025-11-02")
have = os.path.isfile(os.path.join(REF, "flatb_one_nodark.fits")) and os.path.isdir(RAW)
needs_ref = pytest.mark.skipif(not have, reason="IDL branch references or raw frames not present")

from triplespec_reducer import flat as F, wavecal as W, frames as FR, instrument as IN   # noqa: E402

SKIP = {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "EXTEND", "COMMENT",
        "", "DATE", "#######", "HISTORY"}


def same(a, b, what):
    a, b = np.asarray(a), np.asarray(b)
    assert a.shape == b.shape, f"{what}: shape {a.shape} vs IDL {b.shape}"
    ok = (a == b) | (np.isnan(a.astype(float)) & np.isnan(b.astype(float)))
    assert ok.all(), f"{what}: {int((~ok).sum())} values differ from IDL"


def flats(n):
    return [os.path.join(RAW, f"SPEC_DFLAT_{i:04d}.fits") for i in n]


@needs_ref
@pytest.mark.parametrize("name,on,off", [
    ("flatb_one_nodark", [6], []),
    ("flatb_five_onedark", range(2, 7), [60]),
    ("flatb_five_nodark", range(2, 7), []),
    ("flatb_one_onedark", [6], [60]),
])
def test_flat_branch(tmp_path, name, on, off):
    r = F.make_flat(flats(on), flats(off) or None)
    same(r["combined_for_orders"], fits.getdata(os.path.join(REF, name + "_med.fits")), "image for order finding")
    out = F.flat_output_path(str(tmp_path), name)
    F.write_flat(out, r)
    with fits.open(out) as a, fits.open(os.path.join(REF, name + ".fits")) as b:
        for e in (1, 2, 3):
            same(a[e].data, b[e].data, f"{name} ext {e}")
        ha, hb = a[0].header, b[0].header
        for k in hb:
            if k not in SKIP:
                assert ha.cards[k].image == hb.cards[k].image, f"{name}: card {k}"
        assert "".join(ha["HISTORY"]) == "".join(hb["HISTORY"])


def _lines_same(r, prefix):
    for i, o in enumerate(r["orders"]):
        ref = fits.getdata(os.path.join(REF, f"{prefix}_lines2d_{o}.fits"))
        L = r["lines2d"][i]
        same(np.vstack([L[k] for k in ("sx", "sy", "ss", "sl", "sw")]), ref, f"{prefix} lines, order {o}")
    ref = W.read_wavecal(os.path.join(REF, prefix + ".fits"))
    same(r["wavecal"], ref["wavecal"], f"{prefix} wavelength image")
    same(r["spatcal"], ref["spatcal"], f"{prefix} slit-position image")


def _idl_lists(log, run):
    """The two printed lists of the run-th MC_WAVECAL2D_KNAHACK call in the log."""
    txt = open(log).read()
    blocks = re.findall(r"spatial fit:\n\n(.*?)\n\nThe following lines had points identified as bad "
                        r"in the wavelength fit:\n\n(.*?)\n\n", txt, re.S)
    s, w = blocks[run]
    return s.split("\n"), w.split("\n")


def _titles(path):
    raw = open(path, "rb").read().decode("latin-1")
    return re.findall(r"RMS = [0-9. ]*(?:pixels|Angstroms)", raw)


@needs_ref
def test_wavecal_single_frame_and_qa(tmp_path):
    flat = os.path.join(REF, "flat.fits")
    r = W.make_wavecal([os.path.join(RAW, "SPEC_SO2025D-006-0120.fits")], flat)
    assert abs(r["offset"] - -0.656294) < 5e-7
    assert r["sol1d"]["nlines"] == [31, 64, 63, 76, 32]
    _lines_same(r, "wcb_single")
    s, w = _idl_lists(os.path.join(REF, "ref_wavecal_branches.log"), 0)
    assert r["report"]["badslines_printed"] == s
    assert r["report"]["badwlines_printed"] == w
    titles = []
    for st in r["stats"]:
        titles += [st["spatial_title"], st["wave_title"]]
    assert titles == _titles(os.path.join(REF, "wcb_single_Distortion.ps"))
    assert [r["sol1d"]["qa"]["title"]] == _titles(os.path.join(REF, "wcb_single_Wavecal.eps"))
    paths = W.write_wavecal_qa(r, str(tmp_path / "wcb_single"))
    assert all(os.path.getsize(p) > 0 for p in paths)


@needs_ref
def test_wavecal_manual_offset():
    files = [os.path.join(RAW, f"SPEC_SO2025D-006-{n:04d}.fits") for n in range(120, 132)]
    # XGETOFFSET returns a FLOAT; the reference used 1.25 - 1.9 in IDL single precision
    r = W.make_wavecal(files, os.path.join(REF, "flat.fits"),
                       offset=np.float32(1.25) - np.float32(1.9))
    assert r["sol1d"]["nlines"] == [31, 64, 63, 77, 32]
    _lines_same(r, "wcb_manual")


@needs_ref
def test_no_wavecal_images():
    fl = F.read_flat(os.path.join(REF, "flat.fits"))
    w = W.no_wavecal(fl)
    ref = fits.getdata(os.path.join(REF, "wcb_simwave.fits"))     # [wavecal, spatcal], rotated frame
    same(w["wavecal"], ref[0], "simulated wavelength image")
    same(w["spatcal"], ref[1], "simulated slit image")
    assert w["wctype"] == "1D"


def test_short_line_aborts_as_idl():
    """A line found at fewer than LINEDEG+1 rows stops the calibration (IDL), or is skipped on request."""
    img = np.zeros((60, 40), dtype=np.float32)
    ec = np.zeros((1, 2, 2))
    ec[0, 0, 0], ec[0, 1, 0] = 10.0, 13.0            # a slit too short to trace anything
    args = (img, ec, np.array([[0, 39]]), [3], np.array([1.0]), np.array([[20.0, 0.0]]), 0.0,
            np.array([[0.5, 1.5]]), 2.4, 38.0, 1, 5.0, 3, 2)
    with pytest.raises(W.ShortLineError, match="short_lines='skip'"):
        W.findlines2d(*args)
    out = W.findlines2d(*args, short_lines="skip")
    assert out[0]["sx"].size == 0 and len(out[0]["skipped"]) == 1


def test_instrument_file_typing():
    info = IN.read_instrfile()
    assert info.suffix == "*.fits*" and info.combodir == "proc/" and info.bdpxmk == "fakebmask.fits"
    assert info.lincormax == 80000 and isinstance(info.lincormax, int)
    assert info.oversamp == "2" and info.irafname == "FILENAME"     # STRTRIM fields stay text
    assert info.switches["AMPCOR"] == [0, 1] and info.switches["FIXBADPIXELS"] == [0, 0]
    d = IN.defaults()
    assert d.ampcor == (False, True) and d.errorprop == (True, True) and d.version == "4.1"
    assert d.badpixmask.endswith("fakebmask.fits") and d.suffix == "*.fits*"
    assert IN.read_badpixmask().dtype == np.uint8


def test_fsextract_and_paths(tmp_path):
    assert IN.fsextract_index("7-3") == [7, 6, 5, 4, 3]
    assert IN.fsextract_index("1,4-6,9") == [1, 4, 5, 6, 9]
    assert IN.fsextract_filename(" a.fits , b .fits") == ["a.fits", "b.fits"]
    for n in ("spc0007.fits", "spc0006.fits.gz", "spc0005.fits.fz", "notes.txt"):
        (tmp_path / n).write_bytes(b"")
    assert [os.path.basename(p) for p in FR.find_raw_frames(str(tmp_path))] == \
        ["spc0005.fits.fz", "spc0006.fits.gz", "spc0007.fits"]
    got = FR.mkfullpath(str(tmp_path) + "/", "7-5", prefix="spc")
    assert [os.path.basename(p) for p in got] == ["spc0007.fits", "spc0006.fits.gz", "spc0005.fits.fz"]
    with pytest.raises(FileNotFoundError):
        FR.mkfullpath(str(tmp_path) + "/", "8", prefix="spc")
