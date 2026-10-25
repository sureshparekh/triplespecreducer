"""
Point-source extraction as the TS4_SOAR xspextool panel makes it, against
IDL (validation/idl/ref_optimal.pro): optimal extraction (the default, as
PS BG Subtraction switches it on), Median Profile, BG subtraction off,
Guess / Fix positions, Import Trace, PS Orders and A mode with PSNAPS
apertures; and the 1D extraction made without a wavecal file
(validation/idl/ref_1d.pro). Every spectrum is required to be
bit-identical to IDL's.
"""

import os
import sys

import numpy as np
import pytest
from astropy.io import fits

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "validation"))
REF = os.path.join(ROOT, "validation", "idl_out")
RAW = os.path.join(ROOT, "data", "raw", "2025-11-02")
pytestmark = pytest.mark.skipif(
    not (os.path.isfile(os.path.join(REF, "op_spec_0.fits")) and os.path.isdir(RAW)),
    reason="IDL optimal-extraction references or raw frames not present")

import compare_optimal as CO                                     # noqa: E402
from triplespec_reducer import flat as F, wavecal as W, extract as E   # noqa: E402
from triplespec_reducer.instrument import read_atran, read_instrfile   # noqa: E402


@pytest.fixture(scope="module")
def cal():
    fl = F.read_flat(os.path.join(REF, "flat.fits"))
    wc = W.read_wavecal(os.path.join(REF, "wavecal_NGC134.fits"), rotation=fl["rotation"])
    return fl, wc, read_atran(fl["rp"])


@pytest.mark.parametrize("case", [c for c in CO.CASES
                                  if os.path.isfile(os.path.join(REF, f"{c}_spec_0.fits"))])
def test_case(cal, tmp_path, case):
    fl, wc, atran = cal
    r = CO.run(case, fl, wc, atran, str(tmp_path))
    nd, nt, dmax = CO.compare(case, r)
    assert nd == 0, f"{case}: {nd} of {nt} values differ from IDL (max {dmax:.3g})"
    if case == "op":
        h = fits.getheader(r["written"][0])
        assert h["PSFRAD"] == 2.2 and "BGSTART" in h
    if case == "nb":
        h = fits.getheader(r["written"][0])
        assert "PSFRAD" not in h and "BGSTART" not in h
    if case == "p1o":
        h = fits.getheader(r["written"][0])
        assert h["XUNITS"] == "pixels" and "WAVECAL" not in h and "DISPO03" not in h
        assert fits.getdata(r["written"][0]).shape[1] == 4          # flags row
    if case == "it":
        ref = fits.getdata(os.path.join(REF, "it_appos.fits"))
        assert np.array_equal(r["positions"], ref, equal_nan=True)


@pytest.mark.parametrize("case", ["fx", "os"])
def test_trace_file_text(cal, tmp_path, case):
    fl, wc, atran = cal
    same, _ = CO.trace_text(case, fl, wc, atran, str(tmp_path))
    assert same, f"{case}: trace file text differs from IDL's"


def test_spatcoeffs(cal):
    fl, wc, atran = cal
    work = fits.getdata(os.path.join(REF, "ps_work.fits")).astype(np.float32)
    _, spc, mspc = E.mkspatprof2d(work, fl["omask"], wc["wavecal"], wc["spatcal"],
                                  fl["slith_arc"] / fl["slith_pix"], 1, *atran, 0.4, coeffs=True)
    for i, o in enumerate(fl["orders"]):
        assert np.array_equal(spc[o][1], fits.getdata(os.path.join(REF, f"op_spc_c_{i}.fits")))
        assert np.array_equal(mspc[o][1], fits.getdata(os.path.join(REF, f"op_mspc_c_{i}.fits")))


def test_input_checks():
    with pytest.raises(ValueError, match="Ap Radius"):
        E.check_ps_inputs(38.0, 2.5, psfradius=2.2)
    with pytest.raises(ValueError, match="Background Start"):
        E.check_ps_inputs(38.0, 1.0, psbg=(0.9, 5.0, 1))
    with pytest.raises(ValueError, match="four apertures"):
        E.check_ps_inputs(38.0, 1.0, naps=5)
    with pytest.raises(ValueError, match="odd"):
        E.check_ps_inputs(38.0, 1.0, sumap=4)
    E.check_ps_inputs(38.0, 2.4, psbgsub=False, optimal=False)     # sum: limit is the slit


def test_instrfile_types():
    info = read_instrfile()
    assert info.tracesigthresh == 0 and isinstance(info.tracesigthresh, int)   # FIX('0.1')
    assert info.pspsfrad == np.float32(2.2) and info.lincormax == 80000
